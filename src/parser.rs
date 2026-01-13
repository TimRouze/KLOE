use anyhow::{bail, Context, Result};
use bio::io::fasta;
use chrono::{DateTime, Duration, Utc};
use clap::Parser;
use flate2::read::GzDecoder;
use hashbrown::{hash_map::Entry, HashMap};
use num_format::{Locale, ToFormattedString};
use parking_lot::Mutex;
use rayon::{prelude::*, ThreadPoolBuilder};
use simd_minimizers::packed_seq::{PackedSeqVec, SeqVec};
use std::cmp::Reverse;
use std::collections::{BTreeMap, BinaryHeap};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::mpsc;
use std::sync::Arc;
use std::thread;
use threadpool::ThreadPool;
use xz2::read::XzDecoder;
use zstd::stream::read::Decoder as ZstdDecoder;
use zstd::stream::write::Encoder as ZstdEncoder;

const BUFFER_TARGET: usize = 128 * 1024;
const WRITE_BUFFER_TARGET: usize = 8 * 1024 * 1024;
const ZSTD_LEVEL_FAST: i32 = -4; // zstd "fast=4" mode for high speed

fn remove_intermediate_file_best_effort(path: &Path) {
    if let Err(err) = fs::remove_file(path) {
        if err.kind() != std::io::ErrorKind::NotFound {
            eprintln!(
                "Warning: failed to remove intermediate file {}: {}",
                path.display(),
                err
            );
        }
    }
}

struct PartitionWriter {
    path: PathBuf,
    encoder: Mutex<Option<ZstdEncoder<'static, File>>>,
}

struct PartitionReader {
    reader: Option<Box<dyn BufRead + Send>>,
    path: PathBuf,
    id_width: usize,
    header_buf: Vec<u8>,
    seq_buf: Vec<u8>,
}

impl PartitionReader {
    fn new(path: PathBuf, id_width: usize) -> Result<Self> {
        let file = File::open(&path)
            .with_context(|| format!("open partition simplitigs {}", path.display()))?;
        #[cfg(unix)]
        remove_intermediate_file_best_effort(&path);
        let decoder = ZstdDecoder::new(file)
            .with_context(|| format!("build zstd decoder for {}", path.display()))?;
        Ok(Self {
            reader: Some(Box::new(BufReader::new(decoder))),
            path,
            id_width,
            header_buf: Vec::new(),
            seq_buf: Vec::new(),
        })
    }

    fn next_record(&mut self) -> Result<Option<SimplitigRecord>> {
        let Some(reader) = self.reader.as_mut() else {
            return Ok(None);
        };

        let result = read_simplitig_record(
            reader.as_mut(),
            self.id_width,
            &mut self.header_buf,
            &mut self.seq_buf,
        )
        .with_context(|| format!("read simplitig from {}", self.path.display()));

        match result {
            Ok(Some(record)) => Ok(Some(record)),
            Ok(None) => {
                self.reader.take();
                remove_intermediate_file_best_effort(&self.path);
                Ok(None)
            }
            Err(err) => {
                self.reader.take();
                remove_intermediate_file_best_effort(&self.path);
                Err(err)
            }
        }
    }
}

impl Drop for PartitionReader {
    fn drop(&mut self) {
        self.reader.take();
        remove_intermediate_file_best_effort(&self.path);
    }
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

struct KmerEntry {
    ids: Vec<u64>,
    successor: Option<u8>,
    predecessor: Option<u8>,
    succ_ambig: bool,
    pred_ambig: bool,
    visited: bool,
}

fn format_duration(duration: Duration) -> String {
    let std_duration = duration
        .to_std()
        .unwrap_or_else(|_| std::time::Duration::from_secs(0));
    let hours = std_duration.as_secs() / 3600;
    let minutes = (std_duration.as_secs() % 3600) / 60;
    let seconds = std_duration.as_secs() % 60;
    let millis = std_duration.subsec_millis();
    format!("{hours:02}:{minutes:02}:{seconds:02}.{millis:03}")
}

pub fn log_checkpoint(label: &str, start: DateTime<Utc>) {
    let elapsed = Utc::now().signed_duration_since(start);
    println!("{label} wall time: {}", format_duration(elapsed));
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

fn id_width(dataset_count: usize) -> usize {
    let mut value = dataset_count.max(1);
    let mut width = 1usize;
    while value >= 10 {
        value /= 10;
        width += 1;
    }
    width
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

fn zstd_encoder_mt(path: &Path, level: i32, threads: usize) -> Result<ZstdEncoder<'static, File>> {
    let file = File::create(path).with_context(|| format!("create {}", path.display()))?;
    let mut encoder = ZstdEncoder::new(file, level)
        .with_context(|| format!("build zstd encoder for {}", path.display()))?;
    if threads > 1 {
        let _ = encoder.multithread(threads as u32);
    }
    Ok(encoder)
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

fn bits_to_base(bits: u8) -> u8 {
    match bits {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        _ => b'T',
    }
}

fn complement_bits(bits: u8) -> u8 {
    bits ^ 0b11
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
        let encoder = ZstdEncoder::new(file, ZSTD_LEVEL_FAST)
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

fn merge_successor(entry: &mut KmerEntry, next_bits: u8) {
    match entry.successor {
        None => entry.successor = Some(next_bits),
        Some(existing) if existing == next_bits => {}
        Some(_) => entry.succ_ambig = true,
    }
}

fn merge_predecessor(entry: &mut KmerEntry, prev_bits: u8) {
    match entry.predecessor {
        None => entry.predecessor = Some(prev_bits),
        Some(existing) if existing == prev_bits => {}
        Some(_) => entry.pred_ambig = true,
    }
}

fn build_kmer_map_from_partition(
    partition_path: &Path,
    k: usize,
    dataset_count: usize,
) -> Result<HashMap<u64, KmerEntry>> {
    let words = bitset_words(dataset_count);
    let est_entries = fs::metadata(partition_path)
        .map(|m| (m.len() / (k as u64)).saturating_add(1024) as usize)
        .unwrap_or(1024)
        .max(1024);
    let mut map: HashMap<u64, KmerEntry> = HashMap::with_capacity(est_entries);
    let file = File::open(partition_path)
        .with_context(|| format!("open partition {}", partition_path.display()))?;
    #[cfg(unix)]
    remove_intermediate_file_best_effort(partition_path);
    let decoder = ZstdDecoder::new(file)
        .with_context(|| format!("create zstd decoder for {}", partition_path.display()))?;
    let reader = fasta::Reader::from_bufread(BufReader::new(decoder));
    let mask: u64 = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
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
        let mut fwd: u64 = 0;
        let mut rev: u64 = 0;
        let mut len = 0usize;
        for (idx, &base) in seq.iter().enumerate() {
            let Some(bits) = base_to_bits(base) else {
                len = 0;
                fwd = 0;
                rev = 0;
                continue;
            };
            fwd = ((fwd << 2) | bits as u64) & mask;
            rev = (rev >> 2) | ((complement_bits(bits) as u64) << (2 * (k - 1)));
            len += 1;
            if len < k {
                continue;
            }
            let start = idx + 1 - k;
            let forward_is_canon = fwd <= rev;
            let canon = if forward_is_canon { fwd } else { rev };
            let successor = if forward_is_canon {
                seq.get(idx + 1).and_then(|b| base_to_bits(*b))
            } else if start > 0 {
                seq.get(start - 1)
                    .and_then(|b| base_to_bits(*b))
                    .map(complement_bits)
            } else {
                None
            };
            let predecessor = if forward_is_canon {
                if start > 0 {
                    seq.get(start - 1).and_then(|b| base_to_bits(*b))
                } else {
                    None
                }
            } else {
                seq.get(idx + 1)
                    .and_then(|b| base_to_bits(*b))
                    .map(complement_bits)
            };

            let entry = map.entry(canon).or_insert_with(|| KmerEntry {
                ids: vec![0u64; words],
                successor: None,
                predecessor: None,
                succ_ambig: false,
                pred_ambig: false,
                visited: false,
            });
            entry.ids[word_idx] |= bit_mask;
            if let Some(next_bits) = successor {
                merge_successor(entry, next_bits);
            }
            if let Some(prev_bits) = predecessor {
                merge_predecessor(entry, prev_bits);
            }
        }
    }
    Ok(map)
}

fn insert_kmers_from_sequence(
    map: &mut HashMap<u64, KmerEntry>,
    seq: &[u8],
    k: usize,
    ids_bits: &[u64],
) {
    if seq.len() < k {
        return;
    }
    let mask: u64 = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let mut fwd: u64 = 0;
    let mut rev: u64 = 0;
    let mut len = 0usize;
    for &base in seq {
        let Some(bits) = base_to_bits(base) else {
            len = 0;
            fwd = 0;
            rev = 0;
            continue;
        };
        fwd = ((fwd << 2) | bits as u64) & mask;
        rev = (rev >> 2) | ((complement_bits(bits) as u64) << (2 * (k - 1)));
        len += 1;
        if len < k {
            continue;
        }
        let canon = if fwd <= rev { fwd } else { rev };
        map.entry(canon).or_insert_with(|| KmerEntry {
            ids: ids_bits.to_vec(),
            successor: None,
            predecessor: None,
            succ_ambig: false,
            pred_ambig: false,
            visited: false,
        });
    }
}

fn assemble_simplitigs_bidirected(
    kmer_map: HashMap<u64, KmerEntry>,
    k: usize,
    mut sink: impl FnMut(Vec<u8>, Vec<usize>) -> Result<()>,
) -> Result<()> {
    let mask: u64 = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };

    let mut entries: Vec<(u64, KmerEntry)> = kmer_map.into_iter().collect();
    let mut index: HashMap<u64, usize> = HashMap::with_capacity(entries.len());
    for (idx, (bits, _)) in entries.iter().enumerate() {
        index.insert(*bits, idx);
    }

    for seed_idx in 0..entries.len() {
        if entries[seed_idx].1.visited {
            continue;
        }
        let ids_bits = entries[seed_idx].1.ids.clone();
        entries[seed_idx].1.visited = true;

        let canon_bits = entries[seed_idx].0;
        let mut start_bits = canon_bits;
        let rev_bits = revcomp_bits(canon_bits, k);
        if rev_bits != canon_bits {
            let has_neighbor = |bits: u64| -> bool {
                for base in 0u8..4u8 {
                    let next_bits = ((bits << 2) & mask) | base as u64;
                    let next_canon = canonical_bits(next_bits, k);
                    if let Some(&idx) = index.get(&next_canon) {
                        if !entries[idx].1.visited && entries[idx].1.ids == ids_bits {
                            return true;
                        }
                    }
                    let prev_bits =
                        ((base as u64) << (2 * (k - 1)) | (bits >> 2)) & mask;
                    let prev_canon = canonical_bits(prev_bits, k);
                    if let Some(&idx) = index.get(&prev_canon) {
                        if !entries[idx].1.visited && entries[idx].1.ids == ids_bits {
                            return true;
                        }
                    }
                }
                false
            };
            if !has_neighbor(start_bits) && has_neighbor(rev_bits) {
                start_bits = rev_bits;
            }
        }

        let mut seq_bits = start_bits;
        let mut seq = decode_kmer(seq_bits, k);

        // extend to the right using k-1 overlaps (pick any available edge)
        loop {
            let mut found = None;
            for base in 0u8..4u8 {
                let next_bits = ((seq_bits << 2) & mask) | base as u64;
                let next_canon = canonical_bits(next_bits, k);
                let Some(&next_idx) = index.get(&next_canon) else {
                    continue;
                };
                if entries[next_idx].1.visited || entries[next_idx].1.ids != ids_bits {
                    continue;
                }
                found = Some((next_idx, next_bits, base));
                break;
            }
            let Some((next_idx, next_bits, base)) = found else {
                break;
            };
            seq.push(bits_to_base(base));
            entries[next_idx].1.visited = true;
            seq_bits = next_bits;
        }

        // extend to the left using k-1 overlaps (pick any available edge)
        let mut left_bits = start_bits;
        let mut prefix: Vec<u8> = Vec::new();
        loop {
            let mut found = None;
            for base in 0u8..4u8 {
                let prev_bits =
                    ((base as u64) << (2 * (k - 1)) | (left_bits >> 2)) & mask;
                let prev_canon = canonical_bits(prev_bits, k);
                let Some(&prev_idx) = index.get(&prev_canon) else {
                    continue;
                };
                if entries[prev_idx].1.visited || entries[prev_idx].1.ids != ids_bits {
                    continue;
                }
                found = Some((prev_idx, prev_bits, base));
                break;
            }
            let Some((prev_idx, prev_bits, base)) = found else {
                break;
            };
            prefix.push(bits_to_base(base));
            entries[prev_idx].1.visited = true;
            left_bits = prev_bits;
        }

        if !prefix.is_empty() {
            let mut full = Vec::with_capacity(prefix.len() + seq.len());
            for b in prefix.into_iter().rev() {
                full.push(b);
            }
            full.extend(seq);
            seq = full;
        }

        sink(seq, ids_from_bitset(&ids_bits))?;
    }
    Ok(())
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

fn build_kmer_map_from_inputs(file_paths: &[PathBuf], k: usize) -> Result<HashMap<u64, Vec<u64>>> {
    let dataset_count = file_paths.len();
    let words = bitset_words(dataset_count);
    let mut map: HashMap<u64, Vec<u64>> = HashMap::new();
    for (idx, path) in file_paths.iter().enumerate() {
        let word_idx = idx / 64;
        let bit_mask = 1u64 << (idx % 64);
        let reader =
            open_fasta_reader(path).with_context(|| format!("open input {}", path.display()))?;
        for record in reader.records() {
            let record =
                record.with_context(|| format!("read input record in {}", path.display()))?;
            let seq = record.seq();
            if seq.len() < k {
                continue;
            }
            for i in 0..=seq.len() - k {
                let kmer_slice = &seq[i..i + k];
                if let Some(bits) = encode_kmer(kmer_slice) {
                    let canon = canonical_bits(bits, k);
                    let entry = map.entry(canon).or_insert_with(|| vec![0u64; words]);
                    entry[word_idx] |= bit_mask;
                }
            }
        }
    }
    Ok(map)
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
        let record = record
            .with_context(|| format!("read simplitig record from {}", simplitig_path.display()))?;
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
    let mut missing_samples = Vec::new();
    let mut mismatch_samples = Vec::new();
    let mut extra_samples = Vec::new();

    for (kmer, ids_in) in input_map.iter() {
        match output_map.get(kmer) {
            Some(ids_out) => {
                if ids_in != ids_out {
                    mismatched += 1;
                    if mismatch_samples.len() < 5 {
                        mismatch_samples.push(format!(
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
                if missing_samples.len() < 5 {
                    missing_samples.push(format!(
                        "input k-mer {} missing from output",
                        String::from_utf8_lossy(&decode_kmer(*kmer, k))
                    ));
                }
            }
        }
    }

    for kmer in output_map.keys() {
        if !input_map.contains_key(kmer) {
            extra += 1;
            if extra_samples.len() < 5 {
                extra_samples.push(format!(
                    "k-mer {} present in output but absent from input",
                    String::from_utf8_lossy(&decode_kmer(*kmer, k))
                ));
            }
        }
    }

    if missing + mismatched + extra > 0 {
        let mut problems = Vec::new();
        if missing > 0 {
            let suffix = if missing_samples.is_empty() {
                String::new()
            } else {
                format!(" Examples: {}", missing_samples.join("; "))
            };
            problems.push(format!(
                "{} input k-mers never appeared in the output.{}",
                missing, suffix
            ));
        }
        if mismatched > 0 {
            let suffix = if mismatch_samples.is_empty() {
                String::new()
            } else {
                format!(" Examples: {}", mismatch_samples.join("; "))
            };
            problems.push(format!(
                "{} k-mers had differing dataset ids between input and output.{}",
                mismatched, suffix
            ));
        }
        if extra > 0 {
            let suffix = if extra_samples.is_empty() {
                String::new()
            } else {
                format!(" Examples: {}", extra_samples.join("; "))
            };
            problems.push(format!(
                "{} k-mers were present only in the output.{}",
                extra, suffix
            ));
        }
        bail!("k-mer verification failed: {}", problems.join(" "));
    }
    Ok(())
}

#[derive(Debug)]
struct SimplitigRecord {
    header: Vec<u8>,
    key: Vec<u8>,
    seq: Vec<u8>,
}

#[derive(Debug)]
struct HeapItem {
    run_idx: usize,
    record: SimplitigRecord,
}

impl PartialEq for HeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.record.key == other.record.key && self.record.seq == other.record.seq
    }
}

impl Eq for HeapItem {}

impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        simplitig_cmp(&self.record, &other.record)
    }
}

fn build_padded_key_from_header(header: &[u8], id_width: usize) -> Result<Vec<u8>> {
    let rest = header
        .strip_prefix(b">ids:")
        .with_context(|| format!("invalid simplitig header (expected >ids:): '{:?}'", header))?;
    if rest.is_empty() {
        bail!("no dataset ids found in header '{:?}'", header);
    }
    let mut key = Vec::with_capacity(rest.len() + 16);
    let mut first = true;
    let mut idx = 0usize;
    while idx < rest.len() {
        let start = idx;
        while idx < rest.len() && rest[idx] != b',' {
            let b = rest[idx];
            if !(b'0'..=b'9').contains(&b) {
                bail!("invalid dataset id token in header '{:?}'", header);
            }
            idx += 1;
        }
        let token = &rest[start..idx];
        if token.is_empty() {
            bail!("invalid dataset id token in header '{:?}'", header);
        }
        if !first {
            key.push(b',');
        }
        first = false;
        if token.len() < id_width {
            key.extend(std::iter::repeat_n(b'0', id_width - token.len()));
        }
        key.extend_from_slice(token);
        if idx < rest.len() && rest[idx] == b',' {
            idx += 1;
        }
    }
    Ok(key)
}

fn parse_padded_key_ids(key: &[u8]) -> Result<Vec<usize>> {
    if key.is_empty() {
        bail!("no dataset ids found in simplitig key '{:?}'", key);
    }
    let mut ids = Vec::new();
    let mut idx = 0usize;
    while idx < key.len() {
        let mut value = 0usize;
        let mut saw_digit = false;
        while idx < key.len() && key[idx] != b',' {
            let b = key[idx];
            if !(b'0'..=b'9').contains(&b) {
                bail!("invalid dataset id token in simplitig key '{:?}'", key);
            }
            value = value * 10 + (b - b'0') as usize;
            saw_digit = true;
            idx += 1;
        }
        if !saw_digit {
            bail!("invalid dataset id token in simplitig key '{:?}'", key);
        }
        ids.push(value);
        if idx < key.len() && key[idx] == b',' {
            idx += 1;
        }
    }
    if ids.is_empty() {
        bail!("no dataset ids found in simplitig key '{:?}'", key);
    }
    Ok(ids)
}

fn push_usize_decimal(buf: &mut Vec<u8>, mut value: usize) {
    if value == 0 {
        buf.push(b'0');
        return;
    }
    let mut tmp = [0u8; 20];
    let mut len = 0usize;
    while value > 0 {
        tmp[len] = b'0' + (value % 10) as u8;
        len += 1;
        value /= 10;
    }
    for &b in tmp[..len].iter().rev() {
        buf.push(b);
    }
}

fn push_usize_decimal_padded(buf: &mut Vec<u8>, value: usize, width: usize) {
    let mut tmp = [0u8; 20];
    let mut len = 0usize;
    let mut v = value;
    if v == 0 {
        tmp[0] = b'0';
        len = 1;
    } else {
        while v > 0 {
            tmp[len] = b'0' + (v % 10) as u8;
            len += 1;
            v /= 10;
        }
    }
    if len < width {
        buf.extend(std::iter::repeat_n(b'0', width - len));
    }
    for &b in tmp[..len].iter().rev() {
        buf.push(b);
    }
}

fn build_simplitig_header(ids: &[usize]) -> Vec<u8> {
    let mut header = Vec::with_capacity(16 + ids.len() * 4);
    header.extend_from_slice(b">ids:");
    for (idx, &id) in ids.iter().enumerate() {
        if idx > 0 {
            header.push(b',');
        }
        push_usize_decimal(&mut header, id);
    }
    header
}

fn build_simplitig_key(ids: &[usize], id_width: usize) -> Vec<u8> {
    let mut key = Vec::with_capacity(ids.len() * (id_width + 1));
    for (idx, &id) in ids.iter().enumerate() {
        if idx > 0 {
            key.push(b',');
        }
        push_usize_decimal_padded(&mut key, id, id_width);
    }
    key
}

fn read_line_trimmed<R: BufRead + ?Sized>(
    reader: &mut R,
    buf: &mut Vec<u8>,
) -> std::io::Result<usize> {
    buf.clear();
    let n = reader.read_until(b'\n', buf)?;
    while buf.last() == Some(&b'\n') || buf.last() == Some(&b'\r') {
        buf.pop();
    }
    Ok(n)
}

fn read_simplitig_record<R: BufRead + ?Sized>(
    reader: &mut R,
    id_width: usize,
    header_buf: &mut Vec<u8>,
    seq_buf: &mut Vec<u8>,
) -> Result<Option<SimplitigRecord>> {
    if read_line_trimmed(reader, header_buf)? == 0 {
        return Ok(None);
    }
    if read_line_trimmed(reader, seq_buf)? == 0 {
        bail!(
            "unexpected EOF after header '{}'",
            String::from_utf8_lossy(header_buf)
        );
    }
    let header = header_buf.to_vec();
    let key = build_padded_key_from_header(&header, id_width)?;
    let seq = seq_buf.to_vec();
    Ok(Some(SimplitigRecord { header, key, seq }))
}

fn write_simplitig_record_to_buf(record: &SimplitigRecord, buf: &mut Vec<u8>) {
    buf.extend_from_slice(&record.header);
    buf.push(b'\n');
    buf.extend_from_slice(&record.seq);
    buf.push(b'\n');
}

fn simplitig_cmp(a: &SimplitigRecord, b: &SimplitigRecord) -> std::cmp::Ordering {
    match a.key.cmp(&b.key) {
        std::cmp::Ordering::Equal => a.seq.cmp(&b.seq),
        ord => ord,
    }
}

fn sort_simplitigs(records: &mut [SimplitigRecord]) {
    // Parallel sort only for large batches to avoid oversubscription on small inputs.
    const PAR_SORT_THRESHOLD: usize = 200_000;
    if records.len() >= PAR_SORT_THRESHOLD {
        records.par_sort_unstable_by(simplitig_cmp);
    } else {
        records.sort_unstable_by(simplitig_cmp);
    }
}

fn write_sorted_chunk<W: Write>(chunk: &mut [SimplitigRecord], writer: &mut W) -> Result<()> {
    sort_simplitigs(chunk);
    let mut buffer = Vec::with_capacity(WRITE_BUFFER_TARGET);
    for record in chunk.iter() {
        write_simplitig_record_to_buf(record, &mut buffer);
        if buffer.len() >= WRITE_BUFFER_TARGET {
            writer.write_all(&buffer)?;
            buffer.clear();
        }
    }
    if !buffer.is_empty() {
        writer.write_all(&buffer)?;
    }
    Ok(())
}

#[allow(dead_code)]
fn assemble_simplitigs(
    kmer_map: HashMap<u64, KmerEntry>,
    k: usize,
    mut sink: impl FnMut(Vec<u8>, Vec<usize>) -> Result<()>,
) -> Result<()> {
    let mask: u64 = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };

    let mut entries: Vec<(u64, KmerEntry)> = kmer_map.into_iter().collect();
    let mut index: HashMap<u64, usize> = HashMap::with_capacity(entries.len());
    for (idx, (canon, _)) in entries.iter().enumerate() {
        index.insert(*canon, idx);
    }

    for seed_idx in 0..entries.len() {
        if entries[seed_idx].1.visited {
            continue;
        }
        let ids_bits = entries[seed_idx].1.ids.clone();
        entries[seed_idx].1.visited = true;

        let mut seq = decode_kmer(entries[seed_idx].0, k);

        // extend to the right using successor pointers when unambiguous
        let mut right_idx = seed_idx;
        let mut right_bits = entries[seed_idx].0;
        loop {
            let entry = &entries[right_idx].1;
            let Some(base_bits) = entry.successor else {
                break;
            };
            let next_forward = ((right_bits << 2) & mask) | base_bits as u64;
            let next_canon = canonical_bits(next_forward, k);
            let Some(&next_idx) = index.get(&next_canon) else {
                break;
            };
            if entries[next_idx].1.visited || entries[next_idx].1.ids != ids_bits {
                break;
            }
            seq.push(bits_to_base(base_bits));
            entries[next_idx].1.visited = true;
            right_idx = next_idx;
            right_bits = entries[next_idx].0;
        }

        // extend to the left using predecessor pointers when unambiguous
        let mut left_idx = seed_idx;
        let mut left_bits = entries[seed_idx].0;
        let mut prefix: Vec<u8> = Vec::new();
        loop {
            let entry = &entries[left_idx].1;
            let Some(base_bits) = entry.predecessor else {
                break;
            };
            let prev_forward = ((base_bits as u64) << (2 * (k - 1)) | (left_bits >> 2)) & mask;
            let prev_canon = canonical_bits(prev_forward, k);
            let Some(&next_idx) = index.get(&prev_canon) else {
                break;
            };
            if entries[next_idx].1.visited || entries[next_idx].1.ids != ids_bits {
                break;
            }
            prefix.push(bits_to_base(base_bits));
            entries[next_idx].1.visited = true;
            left_idx = next_idx;
            left_bits = entries[next_idx].0;
        }

        if !prefix.is_empty() {
            let mut full = Vec::with_capacity(prefix.len() + seq.len());
            for b in prefix.into_iter().rev() {
                full.push(b);
            }
            full.extend(seq);
            seq = full;
        }

        sink(seq, ids_from_bitset(&ids_bits))?;
    }
    Ok(())
}

fn process_file(
    file_path: &Path,
    file_id: usize,
    encoders: &SharedEncoders,
    k: usize,
    _m: usize,
    window_kmers: usize,
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
        let mut minimizer_positions = Vec::new();
        let kmers_in_seq = seq.len() - k + 1;
        let effective_window = window_kmers.min(kmers_in_seq).max(1);
        let minimizer_values: Vec<u64> = simd_minimizers::minimizers(k, effective_window)
            .super_kmers(&mut superkmer_starts)
            .run(packed_slice, &mut minimizer_positions)
            .values_u64()
            .collect();
        if superkmer_starts.is_empty() {
            continue;
        }

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

fn write_partition_simplitigs(
    partition_path: &Path,
    k: usize,
    dataset_count: usize,
    output_path: &Path,
    threads: usize,
    sort_records: bool,
) -> Result<()> {
    let map = build_kmer_map_from_partition(partition_path, k, dataset_count)?;
    let mut records: Vec<SimplitigRecord> = Vec::with_capacity(map.len());
    let width = id_width(dataset_count);
    assemble_simplitigs_bidirected(map, k, |seq, ids| {
        let header = build_simplitig_header(&ids);
        let key = build_simplitig_key(&ids, width);
        records.push(SimplitigRecord { header, key, seq });
        Ok(())
    })?;

    let mut encoder = zstd_encoder_mt(output_path, ZSTD_LEVEL_FAST, threads)?;
    if sort_records {
        write_sorted_chunk(&mut records, &mut encoder)?;
    } else {
        let mut buffer = Vec::with_capacity(WRITE_BUFFER_TARGET);
        for record in records.iter() {
            write_simplitig_record_to_buf(record, &mut buffer);
            if buffer.len() >= WRITE_BUFFER_TARGET {
                encoder.write_all(&buffer)?;
                buffer.clear();
            }
        }
        if !buffer.is_empty() {
            encoder.write_all(&buffer)?;
        }
    }
    records.clear();

    encoder
        .finish()
        .with_context(|| format!("finalize simplitig output {}", output_path.display()))?;
    Ok(())
}

struct AssemblyJob {
    idx: usize,
    ids: Vec<usize>,
    map: HashMap<u64, KmerEntry>,
}

fn assemble_group(
    ids: Vec<usize>,
    map: HashMap<u64, KmerEntry>,
    k: usize,
    id_width: usize,
) -> Result<Vec<SimplitigRecord>> {
    let header = build_simplitig_header(&ids);
    let key = build_simplitig_key(&ids, id_width);
    let mut records: Vec<SimplitigRecord> = Vec::new();
    assemble_simplitigs_bidirected(map, k, |seq, _| {
        records.push(SimplitigRecord {
            header: header.clone(),
            key: key.clone(),
            seq,
        });
        Ok(())
    })?;
    sort_simplitigs(&mut records);
    Ok(records)
}

fn assemble_sorted_records_parallel(
    mut next_record: impl FnMut() -> Result<Option<SimplitigRecord>>,
    output_path: &Path,
    encoder_threads: usize,
    dataset_count: usize,
    k: usize,
    assembly_threads: usize,
) -> Result<()> {
    let assembly_threads = assembly_threads.max(1);
    let job_capacity = assembly_threads * 2;
    let (job_tx, job_rx) = mpsc::sync_channel::<AssemblyJob>(job_capacity);
    let (result_tx, result_rx) =
        mpsc::sync_channel::<(usize, Result<Vec<SimplitigRecord>>)>(job_capacity);
    let job_rx = Arc::new(Mutex::new(job_rx));
    let id_width = id_width(dataset_count);

    let mut worker_handles = Vec::new();
    for _ in 0..assembly_threads {
        let job_rx = Arc::clone(&job_rx);
        let result_tx = result_tx.clone();
        let handle = thread::spawn(move || -> Result<()> {
            loop {
                let job = {
                    let rx = job_rx.lock();
                    rx.recv()
                };
                let AssemblyJob { idx, ids, map } = match job {
                    Ok(job) => job,
                    Err(_) => break,
                };
                let result = assemble_group(ids, map, k, id_width);
                if result_tx.send((idx, result)).is_err() {
                    break;
                }
            }
            Ok(())
        });
        worker_handles.push(handle);
    }
    drop(result_tx);

    let output_path = output_path.to_path_buf();
    let writer_handle = thread::spawn(move || -> Result<()> {
        let mut encoder = zstd_encoder_mt(&output_path, ZSTD_LEVEL_FAST, encoder_threads)?;
        let mut out_buf = Vec::with_capacity(WRITE_BUFFER_TARGET);
        let mut pending: BTreeMap<usize, Vec<SimplitigRecord>> = BTreeMap::new();
        let mut next_idx = 0usize;
        for (idx, result) in result_rx {
            let records = result?;
            pending.insert(idx, records);
            while let Some(records) = pending.remove(&next_idx) {
                for record in records {
                    write_simplitig_record_to_buf(&record, &mut out_buf);
                    if out_buf.len() >= WRITE_BUFFER_TARGET {
                        encoder.write_all(&out_buf)?;
                        out_buf.clear();
                    }
                }
                next_idx += 1;
            }
        }
        if !pending.is_empty() {
            bail!("missing simplitig group {}", next_idx);
        }
        if !out_buf.is_empty() {
            encoder.write_all(&out_buf)?;
        }
        encoder
            .finish()
            .with_context(|| format!("finalize simplitig output {}", output_path.display()))?;
        Ok(())
    });

    let mut current_key: Option<Vec<u8>> = None;
    let mut current_ids: Vec<usize> = Vec::new();
    let mut current_bits: Vec<u64> = Vec::new();
    let mut kmer_map: HashMap<u64, KmerEntry> = HashMap::new();
    let mut group_idx = 0usize;

    while let Some(record) = next_record()? {
        let key = record.key;
        let key_changed = match current_key.as_ref() {
            Some(existing) => existing.as_slice() != key.as_slice(),
            None => true,
        };
        if key_changed {
            if current_key.is_some() {
                let job = AssemblyJob {
                    idx: group_idx,
                    ids: std::mem::take(&mut current_ids),
                    map: std::mem::take(&mut kmer_map),
                };
                job_tx.send(job)?;
                group_idx += 1;
            }
            current_ids = parse_padded_key_ids(&key)?;
            current_bits = ids_to_bitset(&current_ids, dataset_count)?;
            current_key = Some(key);
        }
        insert_kmers_from_sequence(&mut kmer_map, &record.seq, k, &current_bits);
    }

    if current_key.is_some() {
        let job = AssemblyJob {
            idx: group_idx,
            ids: std::mem::take(&mut current_ids),
            map: std::mem::take(&mut kmer_map),
        };
        job_tx.send(job)?;
    }
    drop(job_tx);

    let mut first_worker_err: Option<anyhow::Error> = None;
    for handle in worker_handles {
        match handle.join() {
            Ok(Ok(())) => {}
            Ok(Err(e)) => {
                if first_worker_err.is_none() {
                    first_worker_err = Some(e);
                }
            }
            Err(_) => {
                if first_worker_err.is_none() {
                    first_worker_err = Some(anyhow::anyhow!("assembly worker thread panicked"));
                }
            }
        }
    }

    let writer_result = match writer_handle.join() {
        Ok(result) => result,
        Err(_) => Err(anyhow::anyhow!("assembly writer thread panicked")),
    };

    if let Some(err) = first_worker_err {
        return Err(err);
    }
    writer_result
}

fn merge_sorted_partitions(
    partition_paths: &[PathBuf],
    output_path: &Path,
    threads: usize,
    sort_records: bool,
    dataset_count: usize,
    k: usize,
) -> Result<()> {
    if !sort_records {
        concatenate_zstd_frames(partition_paths, output_path, threads)?;
        return Ok(());
    }

    // The k-way merge is inherently sequential, but we can parallelize by first merging subsets of
    // runs into sorted streams in worker threads, then performing a final k-way merge of those
    // streams (far fewer inputs) while writing the final output.
    let threads = threads.max(1);
    if threads > 1 && partition_paths.len() > 8 {
        return parallel_streaming_merge(partition_paths, output_path, threads, dataset_count, k);
    }

    kway_merge_sorted_partitions(partition_paths, output_path, threads, dataset_count, k)
}

fn kway_merge_sorted_partitions(
    partition_paths: &[PathBuf],
    output_path: &Path,
    threads: usize,
    dataset_count: usize,
    k: usize,
) -> Result<()> {
    let width = id_width(dataset_count);
    let mut readers: Vec<PartitionReader> = Vec::new();
    for path in partition_paths {
        readers.push(PartitionReader::new(path.clone(), width)?);
    }

    let mut heap: BinaryHeap<Reverse<HeapItem>> = BinaryHeap::new();
    for (idx, reader) in readers.iter_mut().enumerate() {
        if let Some(record) = reader.next_record()? {
            heap.push(Reverse(HeapItem {
                run_idx: idx,
                record,
            }));
        }
    }

    let assembly_threads = threads.saturating_sub(1).max(1);
    let next_record = || -> Result<Option<SimplitigRecord>> {
        if let Some(Reverse(item)) = heap.pop() {
            let run_idx = item.run_idx;
            if let Some(next) = readers[run_idx].next_record()? {
                heap.push(Reverse(HeapItem {
                    run_idx,
                    record: next,
                }));
            }
            Ok(Some(item.record))
        } else {
            Ok(None)
        }
    };

    assemble_sorted_records_parallel(
        next_record,
        output_path,
        threads.max(1),
        dataset_count,
        k,
        assembly_threads,
    )
}

fn merge_partition_group_to_channel(
    partition_paths: Vec<PathBuf>,
    sender: mpsc::SyncSender<Vec<SimplitigRecord>>,
    id_width: usize,
) -> Result<()> {
    let mut readers: Vec<PartitionReader> = Vec::with_capacity(partition_paths.len());
    for path in partition_paths {
        readers.push(PartitionReader::new(path, id_width)?);
    }

    let mut heap: BinaryHeap<Reverse<HeapItem>> = BinaryHeap::new();
    for (idx, reader) in readers.iter_mut().enumerate() {
        if let Some(record) = reader.next_record()? {
            heap.push(Reverse(HeapItem {
                run_idx: idx,
                record,
            }));
        }
    }

    const BATCH_SIZE: usize = 512;
    let mut batch: Vec<SimplitigRecord> = Vec::with_capacity(BATCH_SIZE);
    while let Some(Reverse(item)) = heap.pop() {
        let run_idx = item.run_idx;
        batch.push(item.record);
        if batch.len() >= BATCH_SIZE {
            // Reverse so the consumer can `pop()` in sorted order (O(1)).
            batch.reverse();
            if sender.send(batch).is_err() {
                // Receiver dropped (e.g. due to an upstream error); stop early.
                return Ok(());
            }
            batch = Vec::with_capacity(BATCH_SIZE);
        }
        if let Some(next) = readers[run_idx].next_record()? {
            heap.push(Reverse(HeapItem {
                run_idx,
                record: next,
            }));
        }
    }
    if !batch.is_empty() {
        batch.reverse();
        let _ = sender.send(batch);
    }
    Ok(())
}

struct BatchedRecordStream {
    receiver: mpsc::Receiver<Vec<SimplitigRecord>>,
    buffer: Vec<SimplitigRecord>,
}

impl BatchedRecordStream {
    fn next_record(&mut self) -> Option<SimplitigRecord> {
        loop {
            if let Some(record) = self.buffer.pop() {
                return Some(record);
            }
            match self.receiver.recv() {
                Ok(batch) => self.buffer = batch,
                Err(_) => return None,
            }
        }
    }
}

fn merge_sorted_streams_to_output(
    receivers: Vec<mpsc::Receiver<Vec<SimplitigRecord>>>,
    output_path: &Path,
    encoder_threads: usize,
    dataset_count: usize,
    k: usize,
    assembly_threads: usize,
) -> Result<()> {
    let mut heap: BinaryHeap<Reverse<HeapItem>> = BinaryHeap::new();

    let mut streams: Vec<BatchedRecordStream> = receivers
        .into_iter()
        .map(|receiver| BatchedRecordStream {
            receiver,
            buffer: Vec::new(),
        })
        .collect();

    for (idx, stream) in streams.iter_mut().enumerate() {
        if let Some(record) = stream.next_record() {
            heap.push(Reverse(HeapItem {
                run_idx: idx,
                record,
            }));
        }
    }

    let next_record = || -> Result<Option<SimplitigRecord>> {
        if let Some(Reverse(item)) = heap.pop() {
            let stream_idx = item.run_idx;
            if let Some(record) = streams[stream_idx].next_record() {
                heap.push(Reverse(HeapItem {
                    run_idx: stream_idx,
                    record,
                }));
            }
            Ok(Some(item.record))
        } else {
            Ok(None)
        }
    };

    assemble_sorted_records_parallel(
        next_record,
        output_path,
        encoder_threads.max(1),
        dataset_count,
        k,
        assembly_threads,
    )
}

fn parallel_streaming_merge(
    partition_paths: &[PathBuf],
    output_path: &Path,
    threads: usize,
    dataset_count: usize,
    k: usize,
) -> Result<()> {
    let threads = threads.max(1);
    if threads <= 1 || partition_paths.len() <= 1 {
        return kway_merge_sorted_partitions(partition_paths, output_path, threads, dataset_count, k);
    }

    // Favor output compression: the final encoder tends to dominate wall time.
    let encoder_threads = (threads * 3 / 4).max(1);
    let worker_threads = threads.saturating_sub(encoder_threads).max(1);
    let group_count = worker_threads.min(partition_paths.len()).max(1);
    if group_count <= 1 {
        return kway_merge_sorted_partitions(partition_paths, output_path, threads, dataset_count, k);
    }

    let chunk_size = (partition_paths.len() + group_count - 1) / group_count;
    let channel_capacity = 1usize;
    let width = id_width(dataset_count);

    let mut receivers: Vec<mpsc::Receiver<Vec<SimplitigRecord>>> = Vec::new();
    let mut handles = Vec::new();

    for chunk in partition_paths.chunks(chunk_size) {
        let (tx, rx) = mpsc::sync_channel::<Vec<SimplitigRecord>>(channel_capacity);
        receivers.push(rx);
        let group: Vec<PathBuf> = chunk.to_vec();
        handles.push(thread::spawn(move || {
            merge_partition_group_to_channel(group, tx, width)
        }));
    }

    let merge_result = merge_sorted_streams_to_output(
        receivers,
        output_path,
        encoder_threads,
        dataset_count,
        k,
        worker_threads.max(1),
    );

    let mut first_worker_err: Option<anyhow::Error> = None;
    for handle in handles {
        match handle.join() {
            Ok(Ok(())) => {}
            Ok(Err(e)) => {
                if first_worker_err.is_none() {
                    first_worker_err = Some(e);
                }
            }
            Err(_) => {
                if first_worker_err.is_none() {
                    first_worker_err = Some(anyhow::anyhow!("merge worker thread panicked"));
                }
            }
        }
    }

    if let Some(err) = first_worker_err {
        return Err(err);
    }
    merge_result
}

fn concatenate_zstd_frames(
    partition_paths: &[PathBuf],
    output_path: &Path,
    threads: usize,
) -> Result<()> {
    if partition_paths.is_empty() {
        let encoder = zstd_encoder_mt(output_path, ZSTD_LEVEL_FAST, threads.max(1))?;
        encoder
            .finish()
            .with_context(|| format!("finalize simplitig output {}", output_path.display()))?;
        return Ok(());
    }

    let threads = threads.max(1);
    if threads <= 1 || partition_paths.len() <= 1 {
        return concatenate_zstd_frames_sequential(partition_paths, output_path);
    }

    #[cfg(any(unix, windows))]
    {
        const COPY_BUFFER: usize = 8 * 1024 * 1024;
        let mut offset = 0u64;
        let mut tasks: Vec<(PathBuf, u64, u64)> = Vec::with_capacity(partition_paths.len());
        for path in partition_paths {
            let len = fs::metadata(path)
                .with_context(|| format!("stat partition {}", path.display()))?
                .len();
            tasks.push((path.clone(), offset, len));
            offset = offset
                .checked_add(len)
                .context("output simplitig file too large")?;
        }

        let output = File::create(output_path)
            .with_context(|| format!("create simplitig output {}", output_path.display()))?;
        output
            .set_len(offset)
            .with_context(|| format!("preallocate {}", output_path.display()))?;
        let output = Arc::new(output);

        let pool = ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .context("build merge thread pool")?;
        pool.install(|| {
            tasks
                .par_iter()
                .try_for_each(|(path, start, len)| -> Result<()> {
                    copy_file_region(path, &output, *start, *len, COPY_BUFFER, output_path)
                })
        })?;
        return Ok(());
    }

    #[cfg(not(any(unix, windows)))]
    {
        concatenate_zstd_frames_sequential(partition_paths, output_path)
    }
}

fn concatenate_zstd_frames_sequential(
    partition_paths: &[PathBuf],
    output_path: &Path,
) -> Result<()> {
    const COPY_BUFFER: usize = 8 * 1024 * 1024;
    let output = File::create(output_path)
        .with_context(|| format!("create simplitig output {}", output_path.display()))?;
    let mut writer = BufWriter::with_capacity(COPY_BUFFER, output);
    let mut buffer = vec![0u8; COPY_BUFFER];
    for path in partition_paths {
        {
            let input =
                File::open(path).with_context(|| format!("open partition {}", path.display()))?;
            #[cfg(unix)]
            remove_intermediate_file_best_effort(path);
            let mut reader = BufReader::with_capacity(COPY_BUFFER, input);
            loop {
                let n = reader
                    .read(&mut buffer)
                    .with_context(|| format!("read partition {}", path.display()))?;
                if n == 0 {
                    break;
                }
                writer
                    .write_all(&buffer[..n])
                    .with_context(|| format!("write simplitig output {}", output_path.display()))?;
            }
        }
        remove_intermediate_file_best_effort(path);
    }
    writer
        .flush()
        .with_context(|| format!("flush simplitig output {}", output_path.display()))?;
    Ok(())
}

#[cfg(any(unix, windows))]
fn copy_file_region(
    input_path: &Path,
    output: &File,
    mut output_offset: u64,
    mut remaining: u64,
    buffer_size: usize,
    output_path: &Path,
) -> Result<()> {
    if remaining == 0 {
        return Ok(());
    }
    let mut input = File::open(input_path)
        .with_context(|| format!("open partition {}", input_path.display()))?;
    #[cfg(unix)]
    remove_intermediate_file_best_effort(input_path);
    let buffer_len = buffer_size.min(remaining as usize).max(1024);
    let mut buffer = vec![0u8; buffer_len];
    while remaining > 0 {
        let want = (remaining.min(buffer.len() as u64)) as usize;
        let n = input
            .read(&mut buffer[..want])
            .with_context(|| format!("read partition {}", input_path.display()))?;
        if n == 0 {
            bail!(
                "unexpected EOF while concatenating {} into {}",
                input_path.display(),
                output_path.display()
            );
        }
        write_all_at(output, &buffer[..n], output_offset).with_context(|| {
            format!(
                "write {} bytes at offset {} into {}",
                n,
                output_offset,
                output_path.display()
            )
        })?;
        output_offset += n as u64;
        remaining -= n as u64;
    }
    drop(input);
    remove_intermediate_file_best_effort(input_path);
    Ok(())
}

#[cfg(any(unix, windows))]
fn write_all_at(file: &File, mut buf: &[u8], mut offset: u64) -> std::io::Result<()> {
    while !buf.is_empty() {
        let written = write_at(file, buf, offset)?;
        if written == 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::WriteZero,
                "failed to write to output",
            ));
        }
        offset += written as u64;
        buf = &buf[written..];
    }
    Ok(())
}

#[cfg(unix)]
fn write_at(file: &File, buf: &[u8], offset: u64) -> std::io::Result<usize> {
    use std::os::unix::fs::FileExt;
    file.write_at(buf, offset)
}

#[cfg(windows)]
fn write_at(file: &File, buf: &[u8], offset: u64) -> std::io::Result<usize> {
    use std::os::windows::fs::FileExt;
    file.seek_write(buf, offset)
}

fn parallel_merge_sorted_partitions(
    partition_paths: Vec<PathBuf>,
    output_path: &Path,
    threads: usize,
    sort_records: bool,
    dataset_count: usize,
    k: usize,
) -> Result<()> {
    if partition_paths.is_empty() {
        let encoder = zstd_encoder_mt(output_path, ZSTD_LEVEL_FAST, threads)?;
        encoder
            .finish()
            .with_context(|| format!("finalize simplitig output {}", output_path.display()))?;
        return Ok(());
    }
    if !sort_records && partition_paths.len() == 1 {
        let only = partition_paths[0].clone();
        match fs::rename(&only, output_path) {
            Ok(_) => return Ok(()),
            Err(_) => {
                fs::copy(&only, output_path).with_context(|| {
                    format!("copy {} to {}", only.display(), output_path.display())
                })?;
                let _ = fs::remove_file(&only);
                return Ok(());
            }
        }
    }

    merge_sorted_partitions(
        &partition_paths,
        output_path,
        threads,
        sort_records,
        dataset_count,
        k,
    )
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
    test_parser_with_verify(k, m, partition_power, output_dir, input_fof, threads, false)
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
    run_parser(
        input_fof,
        output_dir,
        k,
        m,
        partition_power,
        threads,
        threads,
        verify_kmers,
        false,
    )
}

pub fn run_parser(
    input_fof: PathBuf,
    output_dir: PathBuf,
    k: usize,
    m: usize,
    partition_power: u32,
    threads: usize,
    compaction_threads: usize,
    verify_kmers: bool,
    skip_sort: bool,
) -> Result<()> {
    if k < m {
        bail!("k-mer length ({}) must be >= minimizer length ({})", k, m);
    }
    let overall_start = Utc::now();
    let partitions = 1u64
        .checked_shl(partition_power)
        .context("partition power too large")?;
    let window = k
        .checked_sub(m)
        .context("failed to compute window size; ensure k > m")?;
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

    let partitioning_start = Utc::now();
    let pool = ThreadPool::new(threads);
    for (idx, file_path) in file_paths.iter().enumerate() {
        let encoders = Arc::clone(&encoders);
        let partitions = partitions;
        let k = k;
        let m = m;
        let stats = Arc::clone(&stats);
        let dataset_count = dataset_count;
        let file_path = file_path.clone();
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
    let total_kmers = total_bases.saturating_sub(total_superkmers.saturating_mul((k - 1) as u64));

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

    log_checkpoint("Step 1 - superkmer partitioning", partitioning_start);

    // Second phase: per-partition compaction into simplitigs with identical ID sets.
    println!("Starting per-partition simplitig compaction...");
    let compaction_start = Utc::now();
    let output_simplitigs = output_dir.join("simplitigs.fa.zst");
    let partition_outputs = Arc::new(Mutex::new(Vec::new()));
    let compaction_pool = ThreadPoolBuilder::new()
        .num_threads(compaction_threads)
        .build()
        .context("build compaction thread pool")?;
    compaction_pool.install(|| {
        encoders.par_iter().enumerate().for_each(|(idx, pw)| {
            let partition_path = pw.path.clone();
            let partition_outputs = Arc::clone(&partition_outputs);
            let part_output = output_dir.join(format!("simplitigs-part-{idx}.fa.zst"));
            let result = write_partition_simplitigs(
                &partition_path,
                k,
                dataset_count,
                &part_output,
                1,
                !skip_sort,
            );
            match result {
                Ok(_) => partition_outputs.lock().push(part_output.clone()),
                Err(e) => {
                    eprintln!(
                        "Failed to assemble simplitigs for {}: {}",
                        partition_path.display(),
                        e
                    );
                    remove_intermediate_file_best_effort(&part_output);
                }
            }
            remove_intermediate_file_best_effort(&partition_path);
        });
    });
    log_checkpoint("Step 2 - simplitig compaction", compaction_start);
    let mut partition_paths = partition_outputs.lock().clone();
    partition_paths.sort();
    let merge_start = Utc::now();
    let merge_result = parallel_merge_sorted_partitions(
        partition_paths.clone(),
        &output_simplitigs,
        compaction_threads,
        !skip_sort,
        dataset_count,
        k,
    );
    for path in &partition_paths {
        remove_intermediate_file_best_effort(path);
    }
    merge_result?;
    log_checkpoint("Step 3 - simplitig merge/sort", merge_start);

    if verify_kmers {
        println!("Verifying k-mer preservation...");
        let output_map = build_output_kmer_map(&output_simplitigs, k, dataset_count)?;
        let input_map = build_kmer_map_from_inputs(&file_paths, k)?;
        verify_kmer_maps(&input_map, &output_map, k)?;
        println!(
            "Verified {} canonical k-mers; all input k-mers were found in the output.",
            input_map.len().to_formatted_string(&Locale::en)
        );
    } else {
        println!("Step 3 - k-mer verification skipped.");
    }

    println!(
        "Simplitig compaction complete. Output: {}",
        output_simplitigs.display()
    );
    println!(
        "Total wall time: {}",
        format_duration(Utc::now().signed_duration_since(overall_start))
    );
    println!("Partitioning complete.");
    Ok(())
}
