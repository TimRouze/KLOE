use anyhow::{anyhow, bail, Context, Result};
use bio::io::fasta;
use clap::Parser;
use flate2::read::GzDecoder;
use hashbrown::{hash_map::Entry, HashMap};
use num_format::{Locale, ToFormattedString};
use parking_lot::Mutex;
use simd_minimizers::packed_seq::{PackedSeqVec, SeqVec};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::mpsc::{self, SyncSender};
use std::sync::Arc;
use std::thread;
use threadpool::ThreadPool;
use xz2::read::XzDecoder;
use zstd::stream::read::Decoder as ZstdDecoder;
use zstd::stream::write::Encoder as ZstdEncoder;

const BUFFER_TARGET: usize = 128 * 1024;
const ZSTD_LEVEL_FAST4: i32 = -4; // zstd "fast=4" mode for lower memory and higher speed

struct PartitionWriter {
    path: PathBuf,
    encoder: Mutex<Option<ZstdEncoder<'static, File>>>,
}

type SharedEncoders = Arc<Vec<PartitionWriter>>;

struct Stats {
    total_superkmers: AtomicU64,
    total_bases: AtomicU64,
}

impl Stats {
    fn new() -> Self {
        Self {
            total_superkmers: AtomicU64::new(0),
            total_bases: AtomicU64::new(0),
        }
    }

    fn add_batch(&self, superkmers: u64, bases: u64) {
        if superkmers == 0 && bases == 0 {
            return;
        }
        self.total_superkmers
            .fetch_add(superkmers, Ordering::Relaxed);
        self.total_bases.fetch_add(bases, Ordering::Relaxed);
    }
}

fn bitset_words(dataset_count: usize) -> usize {
    (dataset_count + 63) / 64
}

fn ids_from_bitset(bits: &[u64]) -> Vec<usize> {
    let mut ids = Vec::new();
    for (word_idx, &word) in bits.iter().enumerate() {
        if word == 0 {
            continue;
        }
        for bit in 0..64 {
            if (word >> bit) & 1 == 1 {
                ids.push(word_idx * 64 + bit + 1); // stored 0-based, output 1-based
            }
        }
    }
    ids
}
fn ensure_nofile_limit(required: u64) -> Result<()> {
    let (soft, hard) = rlimit::getrlimit(rlimit::Resource::NOFILE)?;
    if soft >= required {
        return Ok(());
    }

    let new_soft = required.min(hard);
    rlimit::setrlimit(rlimit::Resource::NOFILE, new_soft, hard).with_context(|| {
        format!(
            "failed to raise file limit to {} (soft) / {} (hard)",
            new_soft, hard
        )
    })?;
    let (soft_after, _) = rlimit::getrlimit(rlimit::Resource::NOFILE)?;
    if soft_after < required {
        bail!(
            "could not raise open file limit high enough (have {}, need {})",
            soft_after,
            required
        );
    }
    Ok(())
}

fn read_fof(path: &Path) -> Result<Vec<PathBuf>> {
    let file = File::open(path).with_context(|| format!("open input list {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut paths = Vec::new();
    for (idx, line) in reader.lines().enumerate() {
        let line =
            line.with_context(|| format!("read line {} from {}", idx + 1, path.display()))?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        paths.push(PathBuf::from(trimmed));
    }
    Ok(paths)
}

fn open_fasta_reader(path: &Path) -> Result<fasta::Reader<Box<dyn BufRead>>> {
    let file = File::open(path).with_context(|| format!("open input {}", path.display()))?;
    let name = path
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();

    let reader: Box<dyn BufRead> = if name.ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else if name.ends_with(".xz") {
        Box::new(BufReader::new(XzDecoder::new(file)))
    } else if name.ends_with(".zst") || name.ends_with(".zstd") {
        let decoder = ZstdDecoder::new(file)
            .with_context(|| format!("create zstd decoder for {}", path.display()))?;
        Box::new(BufReader::new(decoder))
    } else {
        Box::new(BufReader::new(file))
    };

    Ok(fasta::Reader::from_bufread(reader))
}

fn base_to_bits(b: u8) -> Option<u8> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

fn encode_kmer(seq: &[u8]) -> Option<u64> {
    let mut v = 0u64;
    for &b in seq {
        let bits = base_to_bits(b)?;
        v = (v << 2) | bits as u64;
    }
    Some(v)
}

fn revcomp_bits(kmer: u64, k: usize) -> u64 {
    let mut rc = 0u64;
    let mut val = kmer;
    for _ in 0..k {
        let b = (!val) & 0b11;
        rc = (rc << 2) | b;
        val >>= 2;
    }
    rc
}

fn canonical_bits(kmer: u64, k: usize) -> u64 {
    let rc = revcomp_bits(kmer, k);
    if rc < kmer {
        rc
    } else {
        kmer
    }
}

fn decode_kmer(kmer: u64, k: usize) -> Vec<u8> {
    let mut seq = Vec::with_capacity(k);
    for i in (0..k).rev() {
        let bits = (kmer >> (2 * i)) & 0b11;
        let b = match bits {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        };
        seq.push(b);
    }
    seq
}

fn create_partition_encoders(output_dir: &Path, partitions: u64) -> Result<SharedEncoders> {
    let mut writers = Vec::with_capacity(partitions as usize);
    for i in 0..partitions {
        let file_path = output_dir.join(format!("{i}.fa.zst"));
        let file = File::create(&file_path)
            .with_context(|| format!("create partition file {}", file_path.display()))?;
        let encoder = ZstdEncoder::new(file, ZSTD_LEVEL_FAST4)
            .with_context(|| format!("build zstd encoder for {}", file_path.display()))?;
        writers.push(PartitionWriter {
            path: file_path,
            encoder: Mutex::new(Some(encoder)),
        });
    }
    Ok(Arc::new(writers))
}

fn flush_buffer(
    partition_id: usize,
    buffer: &mut Vec<u8>,
    encoders: &SharedEncoders,
) -> Result<()> {
    if buffer.is_empty() {
        return Ok(());
    }
    let mut guard = encoders[partition_id].encoder.lock();
    if let Some(writer) = guard.as_mut() {
        writer.write_all(buffer)?;
    }
    buffer.clear();
    Ok(())
}

fn flush_all_buffers(
    local_buffers: &mut HashMap<usize, Vec<u8>>,
    encoders: &SharedEncoders,
) -> Result<()> {
    for (partition_id, buffer) in local_buffers.iter_mut() {
        flush_buffer(*partition_id, buffer, encoders)?;
    }
    Ok(())
}

fn build_kmer_map_from_partition(
    partition_path: &Path,
    k: usize,
    dataset_count: usize,
) -> Result<HashMap<u64, Vec<u64>>> {
    let words = bitset_words(dataset_count);
    let mut map: HashMap<u64, Vec<u64>> = HashMap::with_capacity(1024);
    let reader = open_fasta_reader(partition_path)?;
    for record in reader.records() {
        let record = record
            .with_context(|| format!("read partition record in {}", partition_path.display()))?;
        let seq = record.seq();
        if seq.len() < k {
            continue;
        }
        let dataset_id: usize = record
            .id()
            .parse()
            .with_context(|| format!("parse dataset id in {}", partition_path.display()))?;
        if dataset_id == 0 || dataset_id > dataset_count {
            bail!(
                "dataset id {} invalid for dataset count {}",
                dataset_id,
                dataset_count
            );
        }
        let id_idx = dataset_id - 1;
        let word_idx = id_idx / 64;
        let bit_mask = 1u64 << (id_idx % 64);
        for i in 0..=seq.len() - k {
            let kmer_slice = &seq[i..i + k];
            if let Some(bits) = encode_kmer(kmer_slice) {
                let canon = canonical_bits(bits, k);
                let entry = map.entry(canon).or_insert_with(|| vec![0u64; words]);
                entry[word_idx] |= bit_mask;
            }
        }
    }
    Ok(map)
}

fn ids_to_bitset(ids: &[usize], dataset_count: usize) -> Result<Vec<u64>> {
    let words = bitset_words(dataset_count);
    let mut bits = vec![0u64; words];
    for &id in ids {
        if id == 0 || id > dataset_count {
            bail!(
                "dataset id {} invalid for dataset count {}",
                id,
                dataset_count
            );
        }
        let idx = id - 1;
        bits[idx / 64] |= 1u64 << (idx % 64);
    }
    Ok(bits)
}

fn merge_bitsets(into: &mut Vec<u64>, from: &[u64]) {
    for (dst, src) in into.iter_mut().zip(from.iter()) {
        *dst |= *src;
    }
}

fn merge_kmer_maps(
    target: &mut HashMap<u64, Vec<u64>>,
    source: HashMap<u64, Vec<u64>>,
) -> Result<()> {
    for (kmer, ids_bits) in source {
        match target.entry(kmer) {
            Entry::Occupied(mut entry) => {
                let existing = entry.get_mut();
                if existing.len() != ids_bits.len() {
                    bail!(
                        "bitset length mismatch while merging k-mer maps ({} vs {})",
                        existing.len(),
                        ids_bits.len()
                    );
                }
                merge_bitsets(existing, &ids_bits);
            }
            Entry::Vacant(entry) => {
                entry.insert(ids_bits);
            }
        }
    }
    Ok(())
}

fn build_input_kmer_map(
    partition_paths: &[PathBuf],
    k: usize,
    dataset_count: usize,
) -> Result<HashMap<u64, Vec<u64>>> {
    let mut global = HashMap::new();
    for path in partition_paths {
        let map = build_kmer_map_from_partition(path, k, dataset_count)
            .with_context(|| format!("build k-mer map from {}", path.display()))?;
        merge_kmer_maps(&mut global, map)?;
    }
    Ok(global)
}

fn parse_simplitig_ids(header: &str, dataset_count: usize) -> Result<Vec<u64>> {
    let Some(rest) = header.strip_prefix("ids:") else {
        bail!("simplitig header '{}' missing ids: prefix", header);
    };
    let mut parsed = Vec::new();
    for part in rest.split(',') {
        if part.is_empty() {
            continue;
        }
        let id: usize = part
            .parse()
            .with_context(|| format!("parse dataset id '{part}' from header '{header}'"))?;
        parsed.push(id);
    }
    if parsed.is_empty() {
        bail!("no dataset ids found in simplitig header '{}'", header);
    }
    ids_to_bitset(&parsed, dataset_count)
}

fn build_output_kmer_map(
    simplitig_path: &Path,
    k: usize,
    dataset_count: usize,
) -> Result<HashMap<u64, Vec<u64>>> {
    let mut map: HashMap<u64, Vec<u64>> = HashMap::new();
    let reader = open_fasta_reader(simplitig_path)?;
    for record in reader.records() {
        let record = record.with_context(|| {
            format!(
                "read simplitig record from {}",
                simplitig_path.display()
            )
        })?;
        let ids_bits = parse_simplitig_ids(record.id(), dataset_count)?;
        let seq = record.seq();
        if seq.len() < k {
            continue;
        }
        for i in 0..=seq.len() - k {
            let slice = &seq[i..i + k];
            if let Some(bits) = encode_kmer(slice) {
                let canon = canonical_bits(bits, k);
                match map.entry(canon) {
                    Entry::Occupied(mut entry) => {
                        if entry.get().len() != ids_bits.len() {
                            bail!(
                                "bitset length mismatch for k-mer {} in {}",
                                String::from_utf8_lossy(&decode_kmer(canon, k)),
                                simplitig_path.display()
                            );
                        }
                        let existing = entry.get_mut();
                        merge_bitsets(existing, &ids_bits);
                    }
                    Entry::Vacant(entry) => {
                        entry.insert(ids_bits.clone());
                    }
                }
            }
        }
    }
    Ok(map)
}

fn verify_kmer_maps(
    input_map: &HashMap<u64, Vec<u64>>,
    output_map: &HashMap<u64, Vec<u64>>,
    k: usize,
) -> Result<()> {
    let mut missing = 0usize;
    let mut mismatched = 0usize;
    let mut extra = 0usize;
    let mut samples = Vec::new();

    for (kmer, ids_in) in input_map.iter() {
        match output_map.get(kmer) {
            Some(ids_out) => {
                if ids_in != ids_out {
                    mismatched += 1;
                    if samples.len() < 5 {
                        samples.push(format!(
                            "k-mer {} has ids {:?} in input but {:?} in output",
                            String::from_utf8_lossy(&decode_kmer(*kmer, k)),
                            ids_from_bitset(ids_in),
                            ids_from_bitset(ids_out)
                        ));
                    }
                }
            }
            None => {
                missing += 1;
                if samples.len() < 5 {
                    samples.push(format!(
                        "k-mer {} missing from output",
                        String::from_utf8_lossy(&decode_kmer(*kmer, k))
                    ));
                }
            }
        }
    }

    for kmer in output_map.keys() {
        if !input_map.contains_key(kmer) {
            extra += 1;
            if samples.len() < 5 {
                samples.push(format!(
                    "k-mer {} present in output but absent from input",
                    String::from_utf8_lossy(&decode_kmer(*kmer, k))
                ));
            }
        }
    }

    if missing + mismatched + extra > 0 {
        bail!(
            "k-mer verification failed (missing: {}, mismatched datasets: {}, extra: {}). Examples: {}",
            missing,
            mismatched,
            extra,
            samples.join("; ")
        );
    }
    Ok(())
}

fn assemble_simplitigs(
    mut kmer_map: HashMap<u64, Vec<u64>>,
    k: usize,
    tx: &SyncSender<(Vec<u8>, Vec<u64>)>,
) -> Result<()> {
    let mask: u64 = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };

    while let Some((&seed_canon, ids)) = kmer_map.iter().next() {
        let ids = ids.clone();
        kmer_map.remove(&seed_canon);

        let mut forward_bits = seed_canon;
        let mut seq = decode_kmer(forward_bits, k);

        // Extend forward
        loop {
            let mut extended = false;
            for base_bits in [0u64, 1, 2, 3] {
                let next_forward = ((forward_bits << 2) & mask) | base_bits;
                let next_canon = canonical_bits(next_forward, k);
                if let Some(entry) = kmer_map.get(&next_canon) {
                    if entry == &ids {
                        forward_bits = next_forward;
                        let next_base = match base_bits {
                            0 => b'A',
                            1 => b'C',
                            2 => b'G',
                            _ => b'T',
                        };
                        seq.push(next_base);
                        kmer_map.remove(&next_canon);
                        extended = true;
                        break;
                    }
                }
            }
            if !extended {
                break;
            }
        }

        // Extend backward
        let mut left_bits = seed_canon;
        loop {
            let mut extended = false;
            for base_bits in [0u64, 1, 2, 3] {
                let prev_forward = ((base_bits << (2 * (k - 1))) | (left_bits >> 2)) & mask;
                let prev_canon = canonical_bits(prev_forward, k);
                if let Some(entry) = kmer_map.get(&prev_canon) {
                    if entry == &ids {
                        left_bits = prev_forward;
                        let prev_base = match base_bits {
                            0 => b'A',
                            1 => b'C',
                            2 => b'G',
                            _ => b'T',
                        };
                        seq.insert(0, prev_base);
                        kmer_map.remove(&prev_canon);
                        extended = true;
                        break;
                    }
                }
            }
            if !extended {
                break;
            }
        }

        let _ = tx.send((seq, ids));
    }
    Ok(())
}

fn process_file(
    file_path: &Path,
    file_id: usize,
    encoders: &SharedEncoders,
    k: usize,
    _m: usize,
    w: usize,
    partitions: u64,
    stats: &Arc<Stats>,
    dataset_count: usize,
) -> Result<()> {
    let reader = open_fasta_reader(file_path)?;
    let mut local_buffers: HashMap<usize, Vec<u8>> = HashMap::new();
    let mut local_superkmers = 0u64;
    let mut local_bases = 0u64;
    if file_id > dataset_count {
        bail!(
            "file id {} exceeds dataset count {}",
            file_id,
            dataset_count
        );
    }

    for record in reader.records() {
        let record = record.with_context(|| format!("parse record in {}", file_path.display()))?;
        let seq = record.seq();
        if seq.len() < k {
            continue;
        }

        let packed_seq = PackedSeqVec::from_ascii(seq);
        let packed_slice = packed_seq.as_slice();

        let mut superkmer_starts = Vec::new();
        let mut minimizer_positions = simd_minimizers::minimizer_positions(packed_slice, k, w);
        if minimizer_positions.is_empty() {
            continue;
        }

        let minimizer_values: Vec<u64> =
            simd_minimizers::minimizers(k, w).super_kmers(&mut superkmer_starts).run(packed_seq.as_slice(), &mut minimizer_positions).values_u64().collect();

        for (idx, (&start, &min_val)) in superkmer_starts
            .iter()
            .zip(minimizer_values.iter())
            .enumerate()
        {
            let start = start as usize;
            let end = if idx + 1 < superkmer_starts.len() {
                let limit = superkmer_starts[idx + 1] as usize + k - 1;
                limit.min(seq.len())
            } else {
                seq.len()
            };

            if end <= start {
                continue;
            }

            let superkmer_slice = &seq[start..end];
            let partition_id = (min_val.wrapping_mul(0x9e3779b97f4a7c15) % partitions) as usize;
            let buffer = local_buffers
                .entry(partition_id)
                .or_insert_with(|| Vec::with_capacity(BUFFER_TARGET));
            buffer.extend_from_slice(format!(">{}\n", file_id).as_bytes());
            buffer.extend_from_slice(superkmer_slice);
            buffer.push(b'\n');

            local_superkmers += 1;
            local_bases += superkmer_slice.len() as u64;

            if buffer.len() >= BUFFER_TARGET {
                flush_buffer(partition_id, buffer, encoders)?;
            }
        }
    }

    flush_all_buffers(&mut local_buffers, encoders)?;
    stats.add_batch(local_superkmers, local_bases);
    Ok(())
}

fn finalize_encoders(encoders: &SharedEncoders) -> Result<()> {
    for pw in encoders.iter() {
        let encoder = pw.encoder.lock().take();
        if let Some(enc) = encoder {
            enc.finish()
                .with_context(|| format!("finish partition writer {}", pw.path.display()))?;
        }
    }
    Ok(())
}

fn partition_size_stats(encoders: &SharedEncoders) -> Result<(u64, u64, u64)> {
    let mut total = 0u64;
    let mut min_size = u64::MAX;
    let mut max_size = 0u64;
    for writer in encoders.iter() {
        let meta = fs::metadata(&writer.path)
            .with_context(|| format!("stat partition file {}", writer.path.display()))?;
        let len = meta.len();
        total += len;
        min_size = min_size.min(len);
        max_size = max_size.max(len);
    }
    if encoders.is_empty() {
        min_size = 0;
    }
    Ok((total, min_size, max_size))
}

pub fn test_parser(
    k: usize,
    m: usize,
    partition_power: u32,
    output_dir: PathBuf,
    input_fof: PathBuf,
    threads: usize,
) -> Result<()> {
    test_parser_with_verify(
        k,
        m,
        partition_power,
        output_dir,
        input_fof,
        threads,
        false,
    )
}

pub fn test_parser_with_verify(
    k: usize,
    m: usize,
    partition_power: u32,
    output_dir: PathBuf,
    input_fof: PathBuf,
    threads: usize,
    verify_kmers: bool,
) -> Result<()> {
    if k < m {
        bail!(
            "k-mer length ({}) must be >= minimizer length ({})",
            k,
            m
        );
    }
    let partitions = 1u64
        .checked_shl(partition_power)
        .context("partition power too large")?;
    let window = k
        .checked_sub(m)
        .map(|v| v + 1)
        .context("failed to compute window size; ensure k >= m")?;
    let required_limit = partitions + 64;

    ensure_nofile_limit(required_limit)?;
    fs::create_dir_all(&output_dir)?;
    let encoders = create_partition_encoders(&output_dir, partitions)?;
    let file_paths = read_fof(&input_fof)?;
    let dataset_count = file_paths.len();
    let stats = Arc::new(Stats::new());

    println!(
        "Starting superkmer partitioning for {} files (k={}, m={}, partitions={}, threads={})",
        file_paths.len(),
        k,
        m,
        partitions,
        threads
    );

    let pool = ThreadPool::new(threads);
    for (idx, file_path) in file_paths.into_iter().enumerate() {
        let encoders = Arc::clone(&encoders);
        let partitions = partitions;
        let stats = Arc::clone(&stats);
        let dataset_count = dataset_count;
        pool.execute(move || {
            if let Err(err) = process_file(
                &file_path,
                idx + 1,
                &encoders,
                k,
                m,
                window,
                partitions,
                &stats,
                dataset_count,
            ) {
                eprintln!("Failed to process {}: {err}", file_path.display());
            }
        });
    }
    pool.join();
    finalize_encoders(&encoders)?;
    let (total_size, min_size, max_size) = partition_size_stats(&encoders)?;
    let total_superkmers = stats.total_superkmers.load(Ordering::Relaxed);
    let total_bases = stats.total_bases.load(Ordering::Relaxed);
    let total_kmers =
        total_bases.saturating_sub(total_superkmers.saturating_mul((k - 1) as u64));

    let loc = Locale::en;
    println!(
        "Superkmer count: {}",
        total_superkmers.to_formatted_string(&loc)
    );
    println!(
        "Total superkmer length (bp): {}",
        total_bases.to_formatted_string(&loc)
    );
    println!(
        "Total kmers covered: {}",
        total_kmers.to_formatted_string(&loc)
    );
    println!(
        "Partition file sizes (bytes): total={}, min={}, max={}",
        total_size.to_formatted_string(&loc),
        min_size.to_formatted_string(&loc),
        max_size.to_formatted_string(&loc)
    );

    // Second phase: per-partition compaction into simplitigs with identical ID sets.
    println!("Starting per-partition simplitig compaction...");
    let output_simplitigs = output_dir.join("simplitigs.fa.zst");
    let (tx, rx): (SyncSender<(Vec<u8>, Vec<u64>)>, _) = mpsc::sync_channel(16);
    let writer_path = output_simplitigs.clone();
    let writer_handle = thread::spawn(move || -> Result<()> {
        let file = File::create(&writer_path)
            .with_context(|| format!("create simplitig output {}", writer_path.display()))?;
        let mut encoder = ZstdEncoder::new(file, ZSTD_LEVEL_FAST4)
            .with_context(|| format!("build zstd encoder for {}", writer_path.display()))?;
        for (seq, ids_bits) in rx {
            let ids = ids_from_bitset(&ids_bits);
            let header = ids
                .iter()
                .map(|id| id.to_string())
                .collect::<Vec<_>>()
                .join(",");
            writeln!(
                encoder,
                ">ids:{}\n{}",
                header,
                String::from_utf8_lossy(&seq)
            )?;
        }
        encoder
            .finish()
            .with_context(|| format!("finalize simplitig output {}", writer_path.display()))?;
        Ok(())
    });

    let compaction_pool = ThreadPool::new(threads);
    for pw in encoders.iter() {
        let partition_path = pw.path.clone();
        let tx = tx.clone();
        compaction_pool.execute(move || {
            match build_kmer_map_from_partition(&partition_path, k, dataset_count) {
                Ok(map) => {
                    if let Err(e) = assemble_simplitigs(map, k, &tx) {
                        eprintln!(
                            "Failed to assemble simplitigs for {}: {}",
                            partition_path.display(),
                            e
                        );
                    }
                }
                Err(e) => {
                    eprintln!(
                        "Failed to build k-mer map for {}: {}",
                        partition_path.display(),
                        e
                    );
                }
            }
        });
    }
    drop(tx);
    compaction_pool.join();
    writer_handle
        .join()
        .map_err(|_| anyhow!("simplitig writer thread panicked"))??;

    if verify_kmers {
        println!("Verifying k-mer preservation...");
        let partition_paths: Vec<PathBuf> = encoders.iter().map(|pw| pw.path.clone()).collect();
        let input_map = build_input_kmer_map(&partition_paths, k, dataset_count)?;
        let output_map = build_output_kmer_map(&output_simplitigs, k, dataset_count)?;
        verify_kmer_maps(&input_map, &output_map, k)?;
        println!(
            "Verified {} canonical k-mers.",
            input_map.len().to_formatted_string(&Locale::en)
        );
    }

    println!(
        "Simplitig compaction complete. Output: {}",
        output_simplitigs.display()
    );
    println!("Partitioning complete.");
    Ok(())
}

pub fn run_parser(k: usize, m: usize, partition_power: u32, output_dir: PathBuf, input_fof: PathBuf, threads: usize, compaction_threads: usize, verify_kmers: bool) -> Result<()> {
    if k < m {
        bail!(
            "k-mer length ({}) must be >= minimizer length ({})",
            k,
            m
        );
    }
    let partitions = 1u64
        .checked_shl(partition_power)
        .context("partition power too large")?;
    let window = k
        .checked_sub(m)
        .map(|v| v + 1)
        .context("failed to compute window size; ensure k >= m")?;
    let required_limit = partitions + 64;

    ensure_nofile_limit(required_limit)?;
    fs::create_dir_all(&output_dir)?;
    let encoders = create_partition_encoders(&output_dir, partitions)?;
    let file_paths = read_fof(&input_fof)?;
    let dataset_count = file_paths.len();
    let stats = Arc::new(Stats::new());

    println!(
        "Starting superkmer partitioning for {} files (k={}, m={}, partitions={}, threads={})",
        file_paths.len(),
        k,
        m,
        partitions,
        threads
    );

    let pool = ThreadPool::new(threads);
    for (idx, file_path) in file_paths.into_iter().enumerate() {
        let encoders = Arc::clone(&encoders);
        let partitions = partitions;
        let k = k;
        let m = m;
        let stats = Arc::clone(&stats);
        let dataset_count = dataset_count;
        pool.execute(move || {
            if let Err(err) = process_file(
                &file_path,
                idx + 1,
                &encoders,
                k,
                m,
                window,
                partitions,
                &stats,
                dataset_count,
            ) {
                eprintln!("Failed to process {}: {err}", file_path.display());
            }
        });
    }
    pool.join();
    finalize_encoders(&encoders)?;
    let (total_size, min_size, max_size) = partition_size_stats(&encoders)?;
    let total_superkmers = stats.total_superkmers.load(Ordering::Relaxed);
    let total_bases = stats.total_bases.load(Ordering::Relaxed);
    let total_kmers =
        total_bases.saturating_sub(total_superkmers.saturating_mul((k - 1) as u64));

    let loc = Locale::en;
    println!(
        "Superkmer count: {}",
        total_superkmers.to_formatted_string(&loc)
    );
    println!(
        "Total superkmer length (bp): {}",
        total_bases.to_formatted_string(&loc)
    );
    println!(
        "Total kmers covered: {}",
        total_kmers.to_formatted_string(&loc)
    );
    println!(
        "Partition file sizes (bytes): total={}, min={}, max={}",
        total_size.to_formatted_string(&loc),
        min_size.to_formatted_string(&loc),
        max_size.to_formatted_string(&loc)
    );

    // Second phase: per-partition compaction into simplitigs with identical ID sets.
    println!("Starting per-partition simplitig compaction...");
    let output_simplitigs = output_dir.join("simplitigs.fa.zst");
    let (tx, rx): (SyncSender<(Vec<u8>, Vec<u64>)>, _) = mpsc::sync_channel(16);
    let writer_path = output_simplitigs.clone();
    let writer_handle = thread::spawn(move || -> Result<()> {
        let file = File::create(&writer_path)
            .with_context(|| format!("create simplitig output {}", writer_path.display()))?;
        let mut encoder = ZstdEncoder::new(file, ZSTD_LEVEL_FAST4)
            .with_context(|| format!("build zstd encoder for {}", writer_path.display()))?;
        for (seq, ids_bits) in rx {
            let ids = ids_from_bitset(&ids_bits);
            let header = ids
                .iter()
                .map(|id| id.to_string())
                .collect::<Vec<_>>()
                .join(",");
            writeln!(
                encoder,
                ">ids:{}\n{}",
                header,
                String::from_utf8_lossy(&seq)
            )?;
        }
        encoder
            .finish()
            .with_context(|| format!("finalize simplitig output {}", writer_path.display()))?;
        Ok(())
    });

    let compaction_pool = ThreadPool::new(compaction_threads);
    for pw in encoders.iter() {
        let partition_path = pw.path.clone();
        let tx = tx.clone();
        compaction_pool.execute(move || {
            match build_kmer_map_from_partition(&partition_path, k, dataset_count) {
                Ok(map) => {
                    if let Err(e) = assemble_simplitigs(map, k, &tx) {
                        eprintln!(
                            "Failed to assemble simplitigs for {}: {}",
                            partition_path.display(),
                            e
                        );
                    }
                }
                Err(e) => {
                    eprintln!(
                        "Failed to build k-mer map for {}: {}",
                        partition_path.display(),
                        e
                    );
                }
            }
        });
    }
    drop(tx);
    compaction_pool.join();
    writer_handle
        .join()
        .map_err(|_| anyhow!("simplitig writer thread panicked"))??;

    if verify_kmers {
        println!("Verifying k-mer preservation...");
        let partition_paths: Vec<PathBuf> = encoders.iter().map(|pw| pw.path.clone()).collect();
        let input_map = build_input_kmer_map(&partition_paths, k, dataset_count)?;
        let output_map = build_output_kmer_map(&output_simplitigs, k, dataset_count)?;
        verify_kmer_maps(&input_map, &output_map, k)?;
        println!(
            "Verified {} canonical k-mers.",
            input_map.len().to_formatted_string(&Locale::en)
        );
    }

    println!(
        "Simplitig compaction complete. Output: {}",
        output_simplitigs.display()
    );
    println!("Partitioning complete.");
    Ok(())
}

/*
fn main() -> Result<()> {
    let args = Args::parse();
    run(args)
}
*/
