use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, Write};
use std::path::{Path, PathBuf};
use std::sync::{mpsc, Arc};
use std::thread;

use zstd::Decoder;

use crate::compress;
use crate::decompress;
use crate::records;
use crate::utils::vec2str;

const REQUIRED_ARCHIVE_FILES: [&str; 5] = [
    "filenames_id.txt",
    "positions_kloe.bin",
    "bucket_sizes.txt",
    "id_to_color_id.txt.zst",
    "tigs_kloe.fa",
];

const INDEX_BYTES_PER_KMER_ESTIMATE: usize = 24;
const INDEX_BUDGET_FRACTION_NUMERATOR: usize = 7;
const INDEX_BUDGET_FRACTION_DENOMINATOR: usize = 10;
const PRODUCER_BATCH_SIZE: usize = 2048;

#[derive(Clone, Debug)]
pub struct MergeConfig {
    pub threads: usize,
    pub minimizer_size: usize,
    pub partition_power: u32,
    pub temp_dir: String,
    pub verify_kmers: bool,
    pub skip_sort: bool,
    pub use_unitigs: bool,
    pub use_matchtigs: bool,
    pub use_eulertigs: bool,
}

impl Default for MergeConfig {
    fn default() -> Self {
        Self {
            threads: 1,
            minimizer_size: 12,
            partition_power: 10,
            temp_dir: String::new(),
            verify_kmers: false,
            skip_sort: false,
            use_unitigs: false,
            use_matchtigs: false,
            use_eulertigs: false,
        }
    }
}

#[derive(Debug)]
struct MergeStats {
    passes: usize,
    emitted_records: u64,
    matched_kmers: usize,
}

#[derive(Debug)]
struct ArchiveInfo {
    filenames: Vec<String>,
    positions: Vec<(u64, u64)>,
    cid_to_ids: Vec<Vec<u32>>,
    tigs_path: PathBuf,
    sizes_path: PathBuf,
    estimated_kmers_per_cid: Vec<usize>,
}

impl ArchiveInfo {
    fn load(root: &Path, k: usize) -> io::Result<Self> {
        ensure_archive_complete(root)?;

        let positions_path = root.join("positions_kloe.bin");
        let sizes_path = root.join("bucket_sizes.txt");
        let cid_path = root.join("id_to_color_id.txt.zst");
        let tigs_path = root.join("tigs_kloe.fa");
        let filenames_path = root.join("filenames_id.txt");

        let positions = load_positions(&positions_path)?;
        if positions.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "positions file '{}' is empty",
                    positions_path.to_string_lossy()
                ),
            ));
        }

        let filenames = load_filenames(&filenames_path)?;
        let dataset_to_cids = load_dataset_to_cids(&cid_path)?;
        if dataset_to_cids.len() != filenames.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "filenames count ({}) differs from id_to_color_id count ({}) in '{}'",
                    filenames.len(),
                    dataset_to_cids.len(),
                    root.to_string_lossy()
                ),
            ));
        }

        let mut cid_to_ids = vec![Vec::<u32>::new(); positions.len()];
        for (dataset_idx, cids) in dataset_to_cids.iter().enumerate() {
            let dataset_id = (dataset_idx + 1) as u32;
            for &cid in cids {
                if cid >= cid_to_ids.len() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "cid index {} outside [0, {}) in '{}'",
                            cid,
                            cid_to_ids.len(),
                            cid_path.to_string_lossy()
                        ),
                    ));
                }
                cid_to_ids[cid].push(dataset_id);
            }
        }
        for ids in &mut cid_to_ids {
            ids.sort_unstable();
            ids.dedup();
        }

        let estimated_kmers_per_cid = estimate_kmers_per_cid(&positions, &sizes_path, k)?;

        Ok(Self {
            filenames,
            positions,
            cid_to_ids,
            tigs_path,
            sizes_path,
            estimated_kmers_per_cid,
        })
    }

    fn active_cids(&self) -> Vec<usize> {
        self.cid_to_ids
            .iter()
            .enumerate()
            .filter_map(|(cid, ids)| if ids.is_empty() { None } else { Some(cid) })
            .collect()
    }

    fn for_each_kmer_in_cids<F>(&self, cids: &[usize], k: usize, mut f: F) -> io::Result<()>
    where
        F: FnMut(usize, u64),
    {
        if cids.is_empty() {
            return Ok(());
        }

        let mut sorted_cids = cids.to_vec();
        sorted_cids.sort_unstable();

        let mut tigs_reader = BufReader::new(File::open(&self.tigs_path)?);
        let mut sizes_reader = BufReader::new(File::open(&self.sizes_path)?);

        let cid_limit = self.positions.len().saturating_sub(1);
        for cid in sorted_cids {
            if cid >= cid_limit {
                continue;
            }

            let (tigs_pos, sizes_pos) = self.positions[cid];
            let sizes = read_bucket_sizes_at(&mut sizes_reader, sizes_pos)?;
            if sizes.is_empty() {
                continue;
            }

            tigs_reader.seek(std::io::SeekFrom::Start(tigs_pos))?;
            for size in sizes {
                if size == 0 {
                    continue;
                }

                let read_size = size.div_ceil(4);
                let mut encoded = vec![0u8; read_size];
                tigs_reader.read_exact(&mut encoded)?;

                if size < k {
                    continue;
                }

                let seq = vec2str(&encoded, &size);
                for_each_canonical_kmer(seq.as_bytes(), k, |kmer| f(cid, kmer));
            }
        }
        Ok(())
    }
}

pub fn merge_archives(
    archive_a_dir: &str,
    archive_b_dir: &str,
    output_dir: &str,
    k: usize,
    memory_gb: usize,
) -> io::Result<()> {
    merge_archives_with_config(
        archive_a_dir,
        archive_b_dir,
        output_dir,
        k,
        memory_gb,
        MergeConfig::default(),
    )
}

pub fn merge_archives_with_config(
    archive_a_dir: &str,
    archive_b_dir: &str,
    output_dir: &str,
    k: usize,
    memory_gb: usize,
    cfg: MergeConfig,
) -> io::Result<()> {
    if k == 0 || k > 32 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("merge supports k in [1, 32], got {k}"),
        ));
    }
    if output_dir.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "missing output directory for merge",
        ));
    }

    let archive_a_root = Path::new(archive_a_dir);
    let archive_b_root = Path::new(archive_b_dir);
    let output_root = Path::new(output_dir);

    fs::create_dir_all(output_root)?;
    let output_dir_norm = normalize_output_dir(output_root);

    println!(
        "Loading archive A from {}",
        archive_a_root.to_string_lossy()
    );
    let archive_a = ArchiveInfo::load(archive_a_root, k)?;
    println!(
        "Loading archive B from {}",
        archive_b_root.to_string_lossy()
    );
    let archive_b = ArchiveInfo::load(archive_b_root, k)?;

    let mut merged_filenames =
        Vec::with_capacity(archive_a.filenames.len() + archive_b.filenames.len());
    merged_filenames.extend(archive_a.filenames.iter().cloned());
    merged_filenames.extend(archive_b.filenames.iter().cloned());

    let offset = archive_a.filenames.len() as u32;
    let mut shifted_b_cid_to_ids = Vec::with_capacity(archive_b.cid_to_ids.len());
    for ids in &archive_b.cid_to_ids {
        let mut shifted = Vec::with_capacity(ids.len());
        for &id in ids {
            shifted.push(id.saturating_add(offset));
        }
        shifted_b_cid_to_ids.push(shifted);
    }

    let a_cids = archive_a.active_cids();
    let b_cids = archive_b.active_cids();
    let index_budget_bytes = index_budget_bytes(memory_gb);
    let batches = plan_a_batches(
        &a_cids,
        &archive_a.estimated_kmers_per_cid,
        index_budget_bytes,
    );

    println!(
        "Starting merge with {} A color sets, {} B color sets, {} pass(es), index budget ~{} MB/pass",
        a_cids.len(),
        b_cids.len(),
        batches.len(),
        index_budget_bytes / (1024 * 1024)
    );

    let merge_tmp_root = create_merge_temp_root(output_root, &cfg.temp_dir)?;
    let stage_archive_dir = merge_tmp_root.join("stage_archive");
    let stage_dumps_dir = merge_tmp_root.join("stage_dumps");
    fs::create_dir_all(&stage_archive_dir)?;
    fs::create_dir_all(&stage_dumps_dir)?;

    let stage_output_dir = normalize_output_dir(&stage_archive_dir);
    let total_datasets = merged_filenames.len();
    let (tx, rx) = mpsc::sync_channel::<records::SimplitigBatch>(16);

    let producer = thread::spawn(move || {
        produce_merged_records(
            archive_a,
            archive_b,
            shifted_b_cid_to_ids,
            batches,
            b_cids,
            k,
            tx,
        )
    });

    let id_cid_offsets = compress::sort_by_bucket_streaming(
        &stage_output_dir,
        total_datasets as u32,
        cfg.threads,
        rx,
    );

    let stats = match producer.join() {
        Ok(res) => res?,
        Err(_) => {
            return Err(io::Error::other(
                "merge producer thread panicked while generating records",
            ))
        }
    };

    let stage_names = (0..total_datasets)
        .map(stage_dataset_name)
        .collect::<Vec<_>>();
    compress::write_filenames_id_offsets(&stage_output_dir, &stage_names, &id_cid_offsets)?;

    println!(
        "Merge k-mer stage complete: passes={}, emitted_records={}, matched_kmers={}",
        stats.passes, stats.emitted_records, stats.matched_kmers
    );

    recompact_stage_archive(
        &stage_archive_dir,
        &stage_dumps_dir,
        &merge_tmp_root,
        &output_dir_norm,
        &merged_filenames,
        k,
        memory_gb,
        &cfg,
    )?;

    if let Err(err) = fs::remove_dir_all(&merge_tmp_root) {
        eprintln!(
            "Warning: could not remove merge temp directory '{}': {}",
            merge_tmp_root.to_string_lossy(),
            err
        );
    }

    println!(
        "Merge complete: output archive written to {} ({})",
        output_root.to_string_lossy(),
        merge_mode_name(&cfg)
    );
    Ok(())
}

fn produce_merged_records(
    archive_a: ArchiveInfo,
    archive_b: ArchiveInfo,
    shifted_b_cid_to_ids: Vec<Vec<u32>>,
    a_batches: Vec<Vec<usize>>,
    b_cids: Vec<usize>,
    k: usize,
    sender: mpsc::SyncSender<records::SimplitigBatch>,
) -> io::Result<MergeStats> {
    let mut matched_b_kmers: HashSet<u64> = HashSet::new();
    let mut emitted_records = 0u64;

    for (pass_idx, batch) in a_batches.iter().enumerate() {
        println!(
            "Merge pass {}/{}: indexing {} A color set(s), then scanning all B color sets",
            pass_idx + 1,
            a_batches.len(),
            batch.len()
        );

        let mut indexed_kmers: HashMap<u64, Vec<u32>> = HashMap::new();
        archive_a.for_each_kmer_in_cids(batch, k, |cid, kmer| {
            let a_ids = &archive_a.cid_to_ids[cid];
            match indexed_kmers.entry(kmer) {
                std::collections::hash_map::Entry::Vacant(slot) => {
                    slot.insert(a_ids.clone());
                }
                std::collections::hash_map::Entry::Occupied(mut slot) => {
                    merge_color_ids(slot.get_mut(), a_ids);
                }
            }
        })?;

        archive_b.for_each_kmer_in_cids(&b_cids, k, |b_cid, kmer| {
            if let Some(merged_ids) = indexed_kmers.get_mut(&kmer) {
                merge_color_ids(merged_ids, &shifted_b_cid_to_ids[b_cid]);
                matched_b_kmers.insert(kmer);
            }
        })?;

        emitted_records += emit_indexed_pass_records(indexed_kmers, k, &sender)?;
    }

    println!(
        "Final B-only scan for non-overlapping k-mers across {} color set(s)",
        b_cids.len()
    );
    emitted_records += emit_b_only_records(
        &archive_b,
        &b_cids,
        &shifted_b_cid_to_ids,
        &matched_b_kmers,
        k,
        &sender,
    )?;

    drop(sender);
    Ok(MergeStats {
        passes: a_batches.len(),
        emitted_records,
        matched_kmers: matched_b_kmers.len(),
    })
}

fn recompact_stage_archive(
    stage_archive_dir: &Path,
    stage_dumps_dir: &Path,
    merge_tmp_root: &Path,
    output_dir_norm: &str,
    merged_filenames: &[String],
    k: usize,
    memory_gb: usize,
    cfg: &MergeConfig,
) -> io::Result<()> {
    let stage_archive_dir_norm = normalize_output_dir(stage_archive_dir);
    let stage_dumps_dir_norm = normalize_output_dir(stage_dumps_dir);

    println!(
        "Recompacting merged k-mer stage into {}",
        merge_mode_name(cfg)
    );
    decompress::decompress(
        &String::from("bucket_sizes.txt"),
        &String::from("id_to_color_id.txt.zst"),
        &String::from("tigs_kloe.fa"),
        &String::from("positions_kloe.bin"),
        &String::from("filenames_id.txt"),
        &stage_dumps_dir_norm,
        &String::from(""),
        stage_archive_dir_norm,
    )?;

    let fof_path = merge_tmp_root.join("merge_stage.fof");
    let mut fof_writer = BufWriter::new(File::create(&fof_path)?);
    for idx in 0..merged_filenames.len() {
        let dump_path = stage_dump_path(stage_dumps_dir, idx);
        if !dump_path.exists() {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!(
                    "expected stage dump '{}' not found",
                    dump_path.to_string_lossy()
                ),
            ));
        }
        writeln!(fof_writer, "{}", dump_path.to_string_lossy())?;
    }
    fof_writer.flush()?;

    let ggcat_cfg = compress::GgcatCompressionConfig {
        memory_gb,
        temp_dir: cfg.temp_dir.clone(),
    };
    compress::compress_with_ggcat(
        &output_dir_norm.to_string(),
        &fof_path.to_string_lossy().to_string(),
        cfg.threads.max(1),
        k,
        cfg.minimizer_size,
        cfg.partition_power,
        cfg.verify_kmers,
        cfg.skip_sort,
        cfg.use_unitigs,
        cfg.use_matchtigs,
        cfg.use_eulertigs,
        ggcat_cfg,
    )?;

    let final_offsets =
        load_offsets_from_filenames_id(&Path::new(output_dir_norm).join("filenames_id.txt"))?;
    compress::write_filenames_id_offsets(output_dir_norm, merged_filenames, &final_offsets)?;

    Ok(())
}

fn emit_indexed_pass_records(
    indexed_kmers: HashMap<u64, Vec<u32>>,
    k: usize,
    sender: &mpsc::SyncSender<records::SimplitigBatch>,
) -> io::Result<u64> {
    let mut records = indexed_kmers
        .into_iter()
        .map(|(kmer, ids)| (ids, kmer))
        .collect::<Vec<_>>();
    records.sort_unstable_by(|(ids_a, kmer_a), (ids_b, kmer_b)| {
        ids_a.cmp(ids_b).then_with(|| kmer_a.cmp(kmer_b))
    });

    let mut emitted = 0u64;
    let mut batch = Vec::with_capacity(PRODUCER_BATCH_SIZE);
    let mut current_color: Option<Arc<Vec<u32>>> = None;

    for (ids, kmer) in records {
        let color = match current_color.as_ref() {
            Some(active) if active.as_ref() == &ids => Arc::clone(active),
            _ => {
                let arc = Arc::new(ids);
                current_color = Some(Arc::clone(&arc));
                arc
            }
        };
        batch.push(records::SimplitigRecord {
            color_ids: color,
            seq: decode_kmer_bits(kmer, k),
        });
        emitted += 1;

        if batch.len() >= PRODUCER_BATCH_SIZE {
            let out = std::mem::take(&mut batch);
            sender.send(out).map_err(channel_send_error)?;
            batch = Vec::with_capacity(PRODUCER_BATCH_SIZE);
        }
    }
    if !batch.is_empty() {
        sender.send(batch).map_err(channel_send_error)?;
    }
    Ok(emitted)
}

fn emit_b_only_records(
    archive_b: &ArchiveInfo,
    b_cids: &[usize],
    shifted_b_cid_to_ids: &[Vec<u32>],
    matched_b_kmers: &HashSet<u64>,
    k: usize,
    sender: &mpsc::SyncSender<records::SimplitigBatch>,
) -> io::Result<u64> {
    let mut emitted = 0u64;
    let mut batch = Vec::with_capacity(PRODUCER_BATCH_SIZE);
    let mut current_cid: Option<usize> = None;
    let mut current_color: Option<Arc<Vec<u32>>> = None;
    let mut send_error: Option<io::Error> = None;

    archive_b.for_each_kmer_in_cids(b_cids, k, |b_cid, kmer| {
        if send_error.is_some() || matched_b_kmers.contains(&kmer) {
            return;
        }

        let color = if current_cid == Some(b_cid) {
            Arc::clone(current_color.as_ref().expect("color must be initialized"))
        } else {
            let arc = Arc::new(shifted_b_cid_to_ids[b_cid].clone());
            current_cid = Some(b_cid);
            current_color = Some(Arc::clone(&arc));
            arc
        };

        batch.push(records::SimplitigRecord {
            color_ids: color,
            seq: decode_kmer_bits(kmer, k),
        });
        emitted += 1;

        if batch.len() >= PRODUCER_BATCH_SIZE {
            let out = std::mem::take(&mut batch);
            if let Err(err) = sender.send(out) {
                send_error = Some(channel_send_error(err));
                return;
            }
            batch = Vec::with_capacity(PRODUCER_BATCH_SIZE);
        }
    })?;

    if let Some(err) = send_error {
        return Err(err);
    }
    if !batch.is_empty() {
        sender.send(batch).map_err(channel_send_error)?;
    }
    Ok(emitted)
}

fn plan_a_batches(
    a_cids: &[usize],
    estimated_kmers_per_cid: &[usize],
    budget_bytes: usize,
) -> Vec<Vec<usize>> {
    if a_cids.is_empty() {
        return vec![Vec::new()];
    }

    let budget = budget_bytes.max(16 * 1024 * 1024);
    let mut batches = Vec::new();
    let mut current = Vec::new();
    let mut current_bytes = 0usize;

    for &cid in a_cids {
        let est_kmers = estimated_kmers_per_cid
            .get(cid)
            .copied()
            .unwrap_or(1)
            .max(1);
        let est_bytes = est_kmers.saturating_mul(INDEX_BYTES_PER_KMER_ESTIMATE);

        if !current.is_empty() && current_bytes.saturating_add(est_bytes) > budget {
            batches.push(std::mem::take(&mut current));
            current_bytes = 0;
        }

        current.push(cid);
        current_bytes = current_bytes.saturating_add(est_bytes);
    }

    if !current.is_empty() {
        batches.push(current);
    }
    if batches.is_empty() {
        batches.push(Vec::new());
    }
    batches
}

fn merge_color_ids(target: &mut Vec<u32>, other: &[u32]) {
    if other.is_empty() {
        return;
    }
    if target.is_empty() {
        target.extend_from_slice(other);
        return;
    }

    let mut merged = Vec::with_capacity(target.len() + other.len());
    let mut i = 0usize;
    let mut j = 0usize;
    while i < target.len() && j < other.len() {
        let a = target[i];
        let b = other[j];
        if a == b {
            merged.push(a);
            i += 1;
            j += 1;
        } else if a < b {
            merged.push(a);
            i += 1;
        } else {
            merged.push(b);
            j += 1;
        }
    }
    if i < target.len() {
        merged.extend_from_slice(&target[i..]);
    }
    if j < other.len() {
        merged.extend_from_slice(&other[j..]);
    }
    merged.dedup();
    *target = merged;
}

fn load_positions(path: &Path) -> io::Result<Vec<(u64, u64)>> {
    let mut file = BufReader::new(File::open(path)?);
    let file_size = file.get_ref().metadata()?.len() as usize;
    if file_size % 16 != 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "positions file '{}' has size {} not divisible by 16",
                path.to_string_lossy(),
                file_size
            ),
        ));
    }
    let num_entries = file_size / 16;
    let mut positions = Vec::with_capacity(num_entries);
    let mut buf = [0u8; 16];
    for _ in 0..num_entries {
        file.read_exact(&mut buf)?;
        let tigs_pos = u64::from_le_bytes(buf[..8].try_into().unwrap());
        let sizes_pos = u64::from_le_bytes(buf[8..16].try_into().unwrap());
        positions.push((tigs_pos, sizes_pos));
    }
    Ok(positions)
}

fn load_filenames(path: &Path) -> io::Result<Vec<String>> {
    let reader = BufReader::new(File::open(path)?);
    let mut names = Vec::new();
    for line_result in reader.lines() {
        let line = line_result?;
        if line.trim().is_empty() {
            continue;
        }
        let Some((name, _)) = line.split_once(':') else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "invalid filenames_id row '{}' in '{}'",
                    line,
                    path.to_string_lossy()
                ),
            ));
        };
        names.push(name.to_string());
    }
    Ok(names)
}

fn load_dataset_to_cids(path: &Path) -> io::Result<Vec<Vec<usize>>> {
    let mut reader = BufReader::new(File::open(path)?);
    let mut all = Vec::new();

    loop {
        let mut len_buf = [0u8; 8];
        match reader.read_exact(&mut len_buf) {
            Ok(()) => {}
            Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => break,
            Err(err) => return Err(err),
        }
        let payload_len = usize::from_le_bytes(len_buf);
        if payload_len == 0 {
            break;
        }

        let mut payload = vec![0u8; payload_len];
        reader.read_exact(&mut payload)?;

        let mut decompressed = Vec::new();
        Decoder::new(&payload[..])?.read_to_end(&mut decompressed)?;
        let text = String::from_utf8(decompressed).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "invalid UTF-8 while decoding '{}' entry: {}",
                    path.to_string_lossy(),
                    err
                ),
            )
        })?;

        let mut cids = Vec::new();
        let mut current_cid = 0usize;
        for token in text.split(',') {
            let token = token.trim();
            if token.is_empty() {
                continue;
            }
            let delta = token.parse::<usize>().map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "invalid cid delta '{}' in '{}': {}",
                        token,
                        path.to_string_lossy(),
                        err
                    ),
                )
            })?;
            current_cid = current_cid.checked_add(delta).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "cid delta overflow while decoding '{}'",
                        path.to_string_lossy()
                    ),
                )
            })?;
            cids.push(current_cid);
        }
        cids.sort_unstable();
        cids.dedup();
        all.push(cids);
    }

    Ok(all)
}

fn load_offsets_from_filenames_id(path: &Path) -> io::Result<Vec<usize>> {
    let reader = BufReader::new(File::open(path)?);
    let mut offsets = Vec::new();
    for line_result in reader.lines() {
        let line = line_result?;
        if line.trim().is_empty() {
            continue;
        }
        let Some((_, offset_str)) = line.rsplit_once(':') else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid filenames_id row '{}'", line),
            ));
        };
        let offset = offset_str.trim().parse::<usize>().map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid filenames_id offset '{}' ({})", offset_str, err),
            )
        })?;
        offsets.push(offset);
    }
    Ok(offsets)
}

fn estimate_kmers_per_cid(
    positions: &[(u64, u64)],
    sizes_path: &Path,
    k: usize,
) -> io::Result<Vec<usize>> {
    let cid_limit = positions.len().saturating_sub(1);
    let mut estimates = vec![0usize; positions.len()];
    let mut sizes_reader = BufReader::new(File::open(sizes_path)?);

    for cid in 0..cid_limit {
        let (_, sizes_pos) = positions[cid];
        let sizes = read_bucket_sizes_at(&mut sizes_reader, sizes_pos)?;
        let mut kmers = 0usize;
        for size in sizes {
            if size >= k {
                kmers = kmers.saturating_add(size - k + 1);
            }
        }
        estimates[cid] = kmers;
    }
    Ok(estimates)
}

fn read_bucket_sizes_at<R: Read + Seek>(reader: &mut R, offset: u64) -> io::Result<Vec<usize>> {
    reader.seek(std::io::SeekFrom::Start(offset))?;

    let mut len_buf = [0u8; 8];
    match reader.read_exact(&mut len_buf) {
        Ok(()) => {}
        Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => return Ok(Vec::new()),
        Err(err) => return Err(err),
    }
    let payload_len = usize::from_le_bytes(len_buf);
    if payload_len == 0 {
        return Ok(Vec::new());
    }

    let mut payload = vec![0u8; payload_len];
    reader.read_exact(&mut payload)?;

    let mut decompressed = Vec::new();
    Decoder::new(&payload[..])?.read_to_end(&mut decompressed)?;
    if decompressed.len() % 8 != 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid decompressed bucket size payload ({} bytes)",
                decompressed.len()
            ),
        ));
    }

    let mut prev = 0usize;
    let mut sizes = Vec::with_capacity(decompressed.len() / 8);
    for chunk in decompressed.chunks_exact(8) {
        let delta = usize::from_le_bytes(chunk.try_into().unwrap());
        let size = prev.saturating_add(delta);
        sizes.push(size);
        prev = size;
    }
    Ok(sizes)
}

fn for_each_canonical_kmer<F>(seq: &[u8], k: usize, mut f: F)
where
    F: FnMut(u64),
{
    if seq.len() < k {
        return;
    }

    let mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let rev_shift = 2 * (k - 1);

    let mut fw = 0u64;
    let mut rc = 0u64;
    let mut valid = 0usize;

    for &b in seq {
        let Some(bits) = base_to_bits(b) else {
            fw = 0;
            rc = 0;
            valid = 0;
            continue;
        };

        fw = ((fw << 2) | bits) & mask;
        let comp = 3u64 - bits;
        rc = (rc >> 2) | (comp << rev_shift);
        valid += 1;

        if valid >= k {
            f(fw.min(rc));
        }
    }
}

fn base_to_bits(base: u8) -> Option<u64> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn decode_kmer_bits(kmer: u64, k: usize) -> Vec<u8> {
    let mut seq = vec![b'A'; k];
    for (i, slot) in seq.iter_mut().enumerate() {
        let shift = 2 * (k - 1 - i);
        let base = ((kmer >> shift) & 0b11) as u8;
        *slot = match base {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        };
    }
    seq
}

fn ensure_archive_complete(root: &Path) -> io::Result<()> {
    for filename in REQUIRED_ARCHIVE_FILES {
        let path = root.join(filename);
        if !path.exists() {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!("missing archive file '{}'", path.to_string_lossy()),
            ));
        }
    }
    Ok(())
}

fn normalize_output_dir(root: &Path) -> String {
    let mut s = root.to_string_lossy().to_string();
    if !s.ends_with('/') {
        s.push('/');
    }
    s
}

fn index_budget_bytes(memory_gb: usize) -> usize {
    let total = memory_gb.saturating_mul(1024 * 1024 * 1024);
    total
        .saturating_mul(INDEX_BUDGET_FRACTION_NUMERATOR)
        .checked_div(INDEX_BUDGET_FRACTION_DENOMINATOR)
        .unwrap_or(total)
}

fn channel_send_error<T>(_: mpsc::SendError<T>) -> io::Error {
    io::Error::new(
        io::ErrorKind::BrokenPipe,
        "merge writer receiver closed while sending records",
    )
}

fn create_merge_temp_root(output_root: &Path, configured_temp_dir: &str) -> io::Result<PathBuf> {
    let base = if configured_temp_dir.is_empty() {
        output_root.to_path_buf()
    } else {
        PathBuf::from(configured_temp_dir)
    };
    fs::create_dir_all(&base)?;

    let pid = std::process::id();
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let root = base.join(format!(".kloe-merge-{pid}-{nanos}"));
    fs::create_dir_all(&root)?;
    Ok(root)
}

fn stage_dataset_name(idx: usize) -> String {
    format!("__kloe_merge_ds_{idx:08}.fa")
}

fn stage_dump_path(stage_dumps_dir: &Path, idx: usize) -> PathBuf {
    let stem = format!("__kloe_merge_ds_{idx:08}");
    stage_dumps_dir.join(format!("Dump_{stem}.fa"))
}

fn merge_mode_name(cfg: &MergeConfig) -> &'static str {
    if cfg.use_unitigs {
        "unitigs"
    } else if cfg.use_matchtigs {
        "matchtigs"
    } else if cfg.use_eulertigs {
        "eulertigs"
    } else {
        "simplitigs"
    }
}
