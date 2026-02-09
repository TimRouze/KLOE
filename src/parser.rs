use crate::assembly::*;

use anyhow::{bail, Context, Result};
use bio::io::fasta;
use chrono::{DateTime, Duration, Utc};
use flate2::read::GzDecoder;
use hashbrown::{hash_map::Entry, HashMap};
use num_format::{Locale, ToFormattedString};
use parking_lot::Mutex;
use rayon::{prelude::*, ThreadPoolBuilder};
use simd_minimizers::packed_seq::{PackedSeqVec, SeqVec};
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::mpsc;
use std::sync::Arc;
use std::thread;
use tempfile::{Builder as TempBuilder, TempPath};
use threadpool::ThreadPool;
use xz2::read::XzDecoder;
use zstd::stream::read::Decoder as ZstdDecoder;
use zstd::stream::write::Encoder as ZstdEncoder;

const BUFFER_TARGET: usize = 128 * 1024;
const WRITE_BUFFER_TARGET: usize = 8 * 1024 * 1024;
const ZSTD_LEVEL_FAST: i32 = -4; // zstd "fast=4" mode for high speed
const PAR_SORT_THRESHOLD: usize = 200_000;
const BUCKET_SORT_SPILL_BYTES: usize = 128 * 1024 * 1024;
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
pub struct SimplitigRecord {
    pub header: Vec<u8>,
    pub key: Vec<u8>,
    pub seq: Vec<u8>,
    /// Dataset IDs (0-indexed). Populated at creation time to avoid re-parsing header in compression.
    pub ids: Vec<u16>,
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

/// Parse 0-indexed dataset IDs from a header like b">ids:1,3,5" (1-indexed in header → 0-indexed).
fn parse_ids_from_header(header: &[u8]) -> Vec<u16> {
    let rest = header.strip_prefix(b">ids:").unwrap_or(header);
    if rest.is_empty() {
        return Vec::new();
    }
    let mut ids = Vec::new();
    let mut idx = 0usize;
    while idx < rest.len() {
        let start = idx;
        while idx < rest.len() && rest[idx] != b',' {
            idx += 1;
        }
        if start < idx {
            let mut val = 0u16;
            for &b in &rest[start..idx] {
                val = val * 10 + (b - b'0') as u16;
            }
            ids.push(val - 1); // 1-indexed in header → 0-indexed
        }
        if idx < rest.len() {
            idx += 1; // skip comma
        }
    }
    ids
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
    let ids = parse_ids_from_header(&header);
    let seq = seq_buf.to_vec();
    Ok(Some(SimplitigRecord { header, key, seq, ids }))
}

fn write_simplitig_record_to_buf(record: &SimplitigRecord, buf: &mut Vec<u8>) {
    buf.extend_from_slice(&record.header);
    buf.push(b'\n');
    buf.extend_from_slice(&record.seq);
    buf.push(b'\n');
}

fn simplitig_cmp(a: &SimplitigRecord, b: &SimplitigRecord) -> std::cmp::Ordering {
    match a.key.cmp(&b.key) {
        std::cmp::Ordering::Equal => a
            .seq
            .len()
            .cmp(&b.seq.len())
            .then_with(|| a.seq.cmp(&b.seq)),
        ord => ord,
    }
}

fn sort_simplitigs(records: &mut [SimplitigRecord]) {
    // Parallel sort only for large batches to avoid oversubscription on small inputs.
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
    let file_id_header = format!(">{}\n", file_id).into_bytes();
    if file_id > dataset_count {
        bail!(
            "file id {} exceeds dataset count {}",
            file_id,
            dataset_count
        );
    }

    let mut superkmer_starts = Vec::new();
    let mut minimizer_positions = Vec::new();

    for record in reader.records() {
        let record = record.with_context(|| format!("parse record in {}", file_path.display()))?;
        let seq = record.seq();
        if seq.len() < k {
            continue;
        }

        let packed_seq = PackedSeqVec::from_ascii(seq);
        let packed_slice = packed_seq.as_slice();

        superkmer_starts.clear();
        minimizer_positions.clear();
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
            buffer.extend_from_slice(&file_id_header);
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

fn build_kmer_map_from_partition(
    partition_path: &Path,
    k: usize,
    dataset_count: usize,
) -> Result<(HashMap<u64, KmerEntry>, Vec<u64>, usize)> {
    let mut map = HashMap::new();
    let mut arena = Vec::new();
    let words = build_kmer_map_reuse(partition_path, k, dataset_count, &mut map, &mut arena)?;
    Ok((map, arena, words))
}

fn build_kmer_map_reuse(
    partition_path: &Path,
    k: usize,
    dataset_count: usize,
    map: &mut HashMap<u64, KmerEntry>,
    arena: &mut Vec<u64>,
) -> Result<usize> {
    let words = bitset_words(dataset_count);
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

            let entry = map.entry(canon).or_insert_with(|| {
                let offset = (arena.len() / words) as u32;
                arena.resize(arena.len() + words, 0);
                KmerEntry::new(offset)
            });
            arena[entry.ids_offset as usize * words + word_idx] |= bit_mask;
            if let Some(next_bits) = successor {
                entry.merge_successor(next_bits);
            }
            if let Some(prev_bits) = predecessor {
                entry.merge_predecessor(prev_bits);
            }
        }
    }
    Ok(words)
}

fn write_partition_simplitigs(
    partition_path: &Path,
    k: usize,
    dataset_count: usize,
    output_path: &Path,
    threads: usize,
    sort_records: bool,
    use_unitigs: bool,
    use_matchtigs: bool,
    use_eulertigs: bool,
) -> Result<()> {
    let (mut map, arena, words) =
        build_kmer_map_from_partition(partition_path, k, dataset_count)?;
    let width = id_width(dataset_count);
    let mut records = Vec::new();

    if use_matchtigs {
        let seqs = build_matchtig_sequences_from_kmers(&mut map, &arena, words, k, threads)?;
        // Collect the union of all dataset IDs in this partition
        let mut all_ids_bits = vec![0u64; words];
        for (_, entry) in map.iter() {
            let eids = entry_ids(entry, &arena, words);
            for (dst, src) in all_ids_bits.iter_mut().zip(eids.iter()) {
                *dst |= *src;
            }
        }
        let ids_usize = ids_from_bitset(&all_ids_bits);
        let header = build_simplitig_header(&ids_usize);
        let key = build_simplitig_key(&ids_usize, width);
        let ids: Vec<u16> = ids_usize.iter().map(|&x| x as u16).collect();
        for seq in seqs {
            records.push(SimplitigRecord { header: header.clone(), key: key.clone(), seq, ids: ids.clone() });
        }
    } else if use_eulertigs {
        let seqs = build_eulertig_sequences_from_kmers(&mut map, &arena, words, k)?;
        let mut all_ids_bits = vec![0u64; words];
        for (_, entry) in map.iter() {
            let eids = entry_ids(entry, &arena, words);
            for (dst, src) in all_ids_bits.iter_mut().zip(eids.iter()) {
                *dst |= *src;
            }
        }
        let ids_usize = ids_from_bitset(&all_ids_bits);
        let header = build_simplitig_header(&ids_usize);
        let key = build_simplitig_key(&ids_usize, width);
        let ids: Vec<u16> = ids_usize.iter().map(|&x| x as u16).collect();
        for seq in seqs {
            records.push(SimplitigRecord { header: header.clone(), key: key.clone(), seq, ids: ids.clone() });
        }
    } else if use_unitigs {
        assemble_unitigs_bidirected(&mut map, &arena, words, k, |seq, ids_bits| {
            let ids_usize = ids_from_bitset(ids_bits);
            let header = build_simplitig_header(&ids_usize);
            let key = build_simplitig_key(&ids_usize, width);
            let ids: Vec<u16> = ids_usize.iter().map(|&x| x as u16).collect();
            records.push(SimplitigRecord { header, key, seq, ids });
            Ok(())
        })?;
    } else {
        // Default: simplitigs
        let use_flat = map.len() >= FLAT_TABLE_THRESHOLD;
        if use_flat {
            let mut flat = FlatKmerTable::from_hashmap(&map);
            map.clear();
            assemble_simplitigs_flat(&mut flat, &arena, words, k, |seq, ids_bits| {
                let ids_usize = ids_from_bitset(ids_bits);
                let header = build_simplitig_header(&ids_usize);
                let key = build_simplitig_key(&ids_usize, width);
                let ids: Vec<u16> = ids_usize.iter().map(|&x| x as u16).collect();
                records.push(SimplitigRecord { header, key, seq, ids });
                Ok(())
            })?;
        } else {
            assemble_simplitigs_bidirected(&mut map, &arena, words, k, |seq, ids_bits| {
                let ids_usize = ids_from_bitset(ids_bits);
                let header = build_simplitig_header(&ids_usize);
                let key = build_simplitig_key(&ids_usize, width);
                let ids: Vec<u16> = ids_usize.iter().map(|&x| x as u16).collect();
                records.push(SimplitigRecord { header, key, seq, ids });
                Ok(())
            })?;
        }
    }

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

    encoder
        .finish()
        .with_context(|| format!("finalize simplitig output {}", output_path.display()))?;
    Ok(())
}

struct SeqRunReader {
    reader: BufReader<File>,
    buf: Vec<u8>,
}

impl SeqRunReader {
    fn new(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("open sequence run {}", path.display()))?;
        Ok(Self {
            reader: BufReader::new(file),
            buf: Vec::new(),
        })
    }

    fn next_seq(&mut self) -> Result<Option<Vec<u8>>> {
        if read_line_trimmed(&mut self.reader, &mut self.buf)? == 0 {
            return Ok(None);
        }
        Ok(Some(self.buf.clone()))
    }
}

#[derive(Debug)]
struct SeqHeapItem {
    run_idx: usize,
    len: usize,
    seq: Vec<u8>,
}

impl PartialEq for SeqHeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.len == other.len && self.run_idx == other.run_idx
    }
}

impl Eq for SeqHeapItem {}

impl PartialOrd for SeqHeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for SeqHeapItem {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.len
            .cmp(&other.len)
            .then_with(|| self.run_idx.cmp(&other.run_idx))
    }
}

fn sort_sequences_by_length(seqs: &mut [Vec<u8>]) {
    if seqs.len() >= PAR_SORT_THRESHOLD {
        seqs.par_sort_unstable_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
    } else {
        seqs.sort_unstable_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
    }
}

fn write_sequence_lines<W: Write>(seqs: &[Vec<u8>], writer: &mut W) -> Result<()> {
    let mut buffer = Vec::with_capacity(WRITE_BUFFER_TARGET);
    for seq in seqs {
        buffer.extend_from_slice(seq);
        buffer.push(b'\n');
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

fn spill_sequence_run(
    seqs: &mut Vec<Vec<u8>>,
    run_paths: &mut Vec<TempPath>,
    spill_dir: &Path,
) -> Result<()> {
    if seqs.is_empty() {
        return Ok(());
    }
    sort_sequences_by_length(seqs);
    let mut tmp = TempBuilder::new()
        .prefix("simplitigs-run-")
        .suffix(".seq")
        .tempfile_in(spill_dir)
        .with_context(|| format!("create sequence run in {}", spill_dir.display()))?;
    {
        let mut writer = BufWriter::with_capacity(WRITE_BUFFER_TARGET, tmp.as_file_mut());
        write_sequence_lines(seqs, &mut writer)?;
        writer.flush()?;
    }
    let temp_path = tmp.into_temp_path();
    run_paths.push(temp_path);
    seqs.clear();
    Ok(())
}

/// Sort simplitigs by length and write them out.
fn chain_simplitigs_for_group<W: Write>(
    group_seqs: &mut Vec<Vec<u8>>,
    header: &[u8],
    encoder: &mut W,
    out_buf: &mut Vec<u8>,
    _k: usize,
) -> Result<()> {
    if group_seqs.is_empty() {
        return Ok(());
    }
    sort_sequences_by_length(group_seqs);
    for seq in group_seqs.iter() {
        out_buf.extend_from_slice(header);
        out_buf.push(b'\n');
        out_buf.extend_from_slice(seq);
        out_buf.push(b'\n');
        if out_buf.len() >= WRITE_BUFFER_TARGET {
            encoder.write_all(out_buf)?;
            out_buf.clear();
        }
    }
    Ok(())
}

/// Reads sorted simplitig records from `next_record`, groups by color key,
/// sorts simplitigs within each group by length, and writes to zstd output.
///
/// For huge groups (>256 MB accumulated), spills sorted batches to temp files
/// and performs an external merge-sort to cap memory usage.
fn chain_and_write_sorted_records(
    mut next_record: impl FnMut() -> Result<Option<SimplitigRecord>>,
    output_path: &Path,
    encoder_threads: usize,
    k: usize,
    pre_sorted: bool,
) -> Result<()> {
    let mut encoder = zstd_encoder_mt(output_path, ZSTD_LEVEL_FAST, encoder_threads)?;
    let mut out_buf = Vec::with_capacity(WRITE_BUFFER_TARGET);

    if pre_sorted {
        // Input is already sorted by (key, seq_len, seq) — stream directly,
        // skipping all per-group buffering, sorting, and spill-to-disk overhead.
        let mut current_key: Option<Vec<u8>> = None;
        let mut groups_written = 0usize;
        while let Some(record) = next_record()? {
            let key_changed = match current_key.as_ref() {
                Some(k) => k.as_slice() != record.key.as_slice(),
                None => true,
            };
            if key_changed {
                if current_key.is_some() {
                    groups_written += 1;
                    if groups_written % 100_000 == 0 {
                        eprintln!("  streamed {} color groups so far", groups_written);
                    }
                }
                current_key = Some(record.key);
            }
            out_buf.extend_from_slice(&record.header);
            out_buf.push(b'\n');
            out_buf.extend_from_slice(&record.seq);
            out_buf.push(b'\n');
            if out_buf.len() >= WRITE_BUFFER_TARGET {
                encoder.write_all(&out_buf)?;
                out_buf.clear();
            }
        }
        if !out_buf.is_empty() {
            encoder.write_all(&out_buf)?;
        }
        encoder
            .finish()
            .with_context(|| format!("finalize simplitig output {}", output_path.display()))?;
        if current_key.is_some() { groups_written += 1; }
        eprintln!("  stream+write complete: {} color groups written", groups_written);
        return Ok(());
    }

    // Non-pre-sorted path: buffer each color group, sort by length, spill if needed.
    let spill_dir = output_path
        .parent()
        .unwrap_or_else(|| Path::new("."))
        .to_path_buf();

    let mut current_key: Option<Vec<u8>> = None;
    let mut current_header: Vec<u8> = Vec::new();
    let mut group_seqs: Vec<Vec<u8>> = Vec::new();
    let mut group_bytes: usize = 0;
    let mut run_paths: Vec<TempPath> = Vec::new();
    let mut groups_written = 0usize;

    while let Some(record) = next_record()? {
        let key_changed = match current_key.as_ref() {
            Some(existing) => existing.as_slice() != record.key.as_slice(),
            None => true,
        };
        if key_changed {
            if current_key.is_some() {
                flush_group_sorted(
                    &mut group_seqs,
                    &mut group_bytes,
                    &mut run_paths,
                    &current_header,
                    &mut encoder,
                    &mut out_buf,
                    &spill_dir,
                    k,
                )?;
                groups_written += 1;
                if groups_written % 100_000 == 0 {
                    eprintln!("  sorted {} color groups so far", groups_written);
                }
            }
            current_header = record.header.clone();
            current_key = Some(record.key);
        }
        group_bytes += record.seq.len();
        group_seqs.push(record.seq);
        // Spill to disk if accumulated bytes exceed threshold.
        if group_bytes >= BUCKET_SORT_SPILL_BYTES {
            spill_sequence_run(&mut group_seqs, &mut run_paths, &spill_dir)?;
            group_bytes = 0;
        }
    }

    // Flush last group.
    if current_key.is_some() {
        flush_group_sorted(
            &mut group_seqs,
            &mut group_bytes,
            &mut run_paths,
            &current_header,
            &mut encoder,
            &mut out_buf,
            &spill_dir,
            k,
        )?;
        groups_written += 1;
    }

    if !out_buf.is_empty() {
        encoder.write_all(&out_buf)?;
    }
    encoder
        .finish()
        .with_context(|| format!("finalize simplitig output {}", output_path.display()))?;

    eprintln!("  sort+write complete: {} color groups written", groups_written);
    Ok(())
}

/// Flush a complete color group: sort by length and write.
/// If there are spilled runs on disk, merge-sort them together with any remaining
/// in-memory sequences, then stream the sorted result into the encoder.
fn flush_group_sorted<W: Write>(
    group_seqs: &mut Vec<Vec<u8>>,
    group_bytes: &mut usize,
    run_paths: &mut Vec<TempPath>,
    header: &[u8],
    encoder: &mut W,
    out_buf: &mut Vec<u8>,
    spill_dir: &Path,
    _k: usize,
) -> Result<()> {
    if run_paths.is_empty() {
        // Small group — sort in memory and write directly.
        chain_simplitigs_for_group(group_seqs, header, encoder, out_buf, _k)?;
        group_seqs.clear();
        *group_bytes = 0;
        return Ok(());
    }

    // Spill any remaining in-memory sequences as a final run.
    if !group_seqs.is_empty() {
        spill_sequence_run(group_seqs, run_paths, spill_dir)?;
    }
    *group_bytes = 0;

    // Merge-sort all runs and stream into the encoder.
    let mut readers: Vec<SeqRunReader> = Vec::with_capacity(run_paths.len());
    for path in run_paths.iter() {
        readers.push(SeqRunReader::new(path)?);
    }

    let mut heap: BinaryHeap<Reverse<SeqHeapItem>> = BinaryHeap::new();
    for (idx, reader) in readers.iter_mut().enumerate() {
        if let Some(seq) = reader.next_seq()? {
            let len = seq.len();
            heap.push(Reverse(SeqHeapItem { run_idx: idx, len, seq }));
        }
    }

    while let Some(Reverse(item)) = heap.pop() {
        out_buf.extend_from_slice(header);
        out_buf.push(b'\n');
        out_buf.extend_from_slice(&item.seq);
        out_buf.push(b'\n');
        if out_buf.len() >= WRITE_BUFFER_TARGET {
            encoder.write_all(out_buf)?;
            out_buf.clear();
        }
        if let Some(seq) = readers[item.run_idx].next_seq()? {
            let len = seq.len();
            heap.push(Reverse(SeqHeapItem { run_idx: item.run_idx, len, seq }));
        }
    }

    // Clean up temp files.
    for path in run_paths.drain(..) {
        remove_intermediate_file_best_effort(path.as_ref());
    }

    Ok(())
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

    chain_and_write_sorted_records(
        next_record,
        output_path,
        threads.max(1),
        k,
        true, // k-way merge output is already sorted
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
    k: usize,
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

    chain_and_write_sorted_records(
        next_record,
        output_path,
        encoder_threads.max(1),
        k,
        true, // parallel merge output is already sorted
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

    // ZSTD at level -4 (fast mode) is I/O-bound, not CPU-bound — 1 encoder thread
    // is sufficient. Give remaining threads to merge workers so the parallel merge
    // path actually activates (previously 3/4 went to encoder, leaving 1 worker,
    // which fell through to the single-threaded kway_merge_sorted_partitions).
    let encoder_threads = 1;
    let worker_threads = threads.saturating_sub(encoder_threads).max(1);
    let group_count = worker_threads.min(partition_paths.len()).max(1);
    if group_count <= 1 {
        return kway_merge_sorted_partitions(partition_paths, output_path, threads, dataset_count, k);
    }

    let chunk_size = (partition_paths.len() + group_count - 1) / group_count;
    let channel_capacity = 16usize;
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
        k,
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
            let file = enc.finish()
                .with_context(|| format!("finish partition writer {}", pw.path.display()))?;
            file.sync_all()
                .with_context(|| format!("sync partition file {}", pw.path.display()))?;
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

pub fn run_parser(
    input_fof: PathBuf,
    output_dir: PathBuf,
    k: usize,
    m: usize,
    partition_power: u32,
    threads: usize,
    verify_kmers: bool,
    skip_sort: bool,
    use_unitigs: bool,
    use_matchtigs: bool,
    use_eulertigs: bool,
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
    // Sync output directory metadata to ensure all partition files are visible on NFS
    if let Ok(dir) = File::open(&output_dir) {
        let _ = dir.sync_all();
    }
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
    // Uses explicit worker threads with crossbeam channel for work stealing,
    // reusing HashMap/arena allocations across partitions.
    println!("Starting per-partition simplitig compaction...");
    let compaction_start = Utc::now();
    let output_simplitigs = output_dir.join("simplitigs.fa.zst");
    let partition_outputs = Arc::new(Mutex::new(Vec::new()));
    // Sort partitions largest-first so the biggest ones start early and don't
    // become stragglers at the end of the parallel loop.
    let mut partition_indices: Vec<(usize, u64)> = encoders
        .iter()
        .enumerate()
        .map(|(idx, pw)| {
            let size = fs::metadata(&pw.path).map(|m| m.len()).unwrap_or(0);
            (idx, size)
        })
        .collect();
    partition_indices.sort_by(|a, b| b.1.cmp(&a.1));

    {
        let compaction_pool = ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .context("failed to build compaction thread pool")?;
        compaction_pool.install(|| {
            partition_indices
                .par_iter()
                .for_each(|&(idx, _size)| {
                    let pw = &encoders[idx];
                    let partition_path = pw.path.clone();
                    let part_output =
                        output_dir.join(format!("simplitigs-part-{idx}.fa.zst"));
                    let result = write_partition_simplitigs(
                        &partition_path,
                        k,
                        dataset_count,
                        &part_output,
                        1,
                        !skip_sort,
                        use_unitigs,
                        use_matchtigs,
                        use_eulertigs,
                    );
                    match result {
                        Ok(_) => partition_outputs.lock().push(part_output.clone()),
                        Err(e) => {
                            eprintln!(
                                "Failed to assemble simplitigs for {}: {:#}",
                                partition_path.display(),
                                e
                            );
                            remove_intermediate_file_best_effort(&part_output);
                        }
                    }
                    remove_intermediate_file_best_effort(&partition_path);
                });
        });
    }
    log_checkpoint("Step 2 - simplitig compaction", compaction_start);
    let mut partition_paths = partition_outputs.lock().clone();
    partition_paths.sort();
    let merge_start = Utc::now();
    let merge_result = parallel_merge_sorted_partitions(
        partition_paths.clone(),
        &output_simplitigs,
        threads,
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

/// Streaming variant of run_parser. Runs Steps 1+2 synchronously, then spawns
/// a background thread for the k-way merge (Step 3) that sends sorted records
/// through a channel. Returns the receiver and a join handle for the merge thread.
///
/// The caller should consume all records from the receiver, then join the handle
/// to propagate any merge errors.
pub fn run_parser_streaming(
    input_fof: PathBuf,
    output_dir: PathBuf,
    k: usize,
    m: usize,
    partition_power: u32,
    threads: usize,
    verify_kmers: bool,
    skip_sort: bool,
    use_unitigs: bool,
    use_matchtigs: bool,
    use_eulertigs: bool,
) -> Result<(
    mpsc::Receiver<SimplitigRecord>,
    thread::JoinHandle<Result<()>>,
    usize, // dataset_count
)> {
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

    // Step 1: superkmer partitioning
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
    if let Ok(dir) = File::open(&output_dir) {
        let _ = dir.sync_all();
    }
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

    // Step 2: per-partition compaction
    println!("Starting per-partition simplitig compaction...");
    let compaction_start = Utc::now();
    let partition_outputs = Arc::new(Mutex::new(Vec::new()));
    let mut partition_indices: Vec<(usize, u64)> = encoders
        .iter()
        .enumerate()
        .map(|(idx, pw)| {
            let size = fs::metadata(&pw.path).map(|m| m.len()).unwrap_or(0);
            (idx, size)
        })
        .collect();
    partition_indices.sort_by(|a, b| b.1.cmp(&a.1));

    {
        let compaction_pool = ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .context("failed to build compaction thread pool")?;
        compaction_pool.install(|| {
            partition_indices
                .par_iter()
                .for_each(|&(idx, _size)| {
                    let pw = &encoders[idx];
                    let partition_path = pw.path.clone();
                    let part_output =
                        output_dir.join(format!("simplitigs-part-{idx}.fa.zst"));
                    let result = write_partition_simplitigs(
                        &partition_path,
                        k,
                        dataset_count,
                        &part_output,
                        1,
                        !skip_sort,
                        use_unitigs,
                        use_matchtigs,
                        use_eulertigs,
                    );
                    match result {
                        Ok(_) => partition_outputs.lock().push(part_output.clone()),
                        Err(e) => {
                            eprintln!(
                                "Failed to assemble simplitigs for {}: {:#}",
                                partition_path.display(),
                                e
                            );
                            remove_intermediate_file_best_effort(&part_output);
                        }
                    }
                    remove_intermediate_file_best_effort(&partition_path);
                });
        });
    }
    log_checkpoint("Step 2 - simplitig compaction", compaction_start);

    let mut partition_paths = partition_outputs.lock().clone();
    partition_paths.sort();

    // Step 3: spawn merge thread that sends records via channel
    let merge_start = Utc::now();
    let (record_tx, record_rx) = mpsc::sync_channel::<SimplitigRecord>(8192);
    let sort_records = !skip_sort;

    let merge_handle = thread::spawn(move || -> Result<()> {
        if partition_paths.is_empty() {
            drop(record_tx);
            log_checkpoint("Step 3 - simplitig merge/sort (streaming)", merge_start);
            return Ok(());
        }

        if sort_records {
            // k-way merge sending records to channel
            merge_to_channel(&partition_paths, record_tx, dataset_count, k)?;
        } else {
            // No sorting: just stream records from each partition
            let width = id_width(dataset_count);
            for path in &partition_paths {
                let mut reader = PartitionReader::new(path.clone(), width)?;
                while let Some(record) = reader.next_record()? {
                    if record_tx.send(record).is_err() {
                        break;
                    }
                }
            }
            drop(record_tx);
        }

        for path in &partition_paths {
            remove_intermediate_file_best_effort(path);
        }
        log_checkpoint("Step 3 - simplitig merge/sort (streaming)", merge_start);

        if verify_kmers {
            eprintln!("Warning: --verify-kmers is not supported in streaming mode (skipped).");
        }

        println!(
            "Total parser wall time: {}",
            format_duration(Utc::now().signed_duration_since(overall_start))
        );
        Ok(())
    });

    Ok((record_rx, merge_handle, dataset_count))
}

/// Performs k-way merge of sorted partition files, sending records to a channel.
fn merge_to_channel(
    partition_paths: &[PathBuf],
    sender: mpsc::SyncSender<SimplitigRecord>,
    dataset_count: usize,
    _k: usize,
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

    let mut records_sent = 0u64;
    while let Some(Reverse(item)) = heap.pop() {
        let run_idx = item.run_idx;
        if sender.send(item.record).is_err() {
            break; // receiver dropped
        }
        records_sent += 1;
        if records_sent % 1_000_000 == 0 {
            eprintln!("  merge: sent {} records so far", records_sent);
        }
        if let Some(next) = readers[run_idx].next_record()? {
            heap.push(Reverse(HeapItem {
                run_idx,
                record: next,
            }));
        }
    }
    drop(sender);
    eprintln!("  merge complete: {} records sent", records_sent);
    Ok(())
}
