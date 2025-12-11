use anyhow::{bail, Context, Result};
use bio::io::fasta;
use clap::Parser;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use flate2::read::GzDecoder;
use hashbrown::{hash_map::Entry, HashMap};
use num_format::{Locale, ToFormattedString};
use parking_lot::Mutex;
use simd_minimizers::packed_seq::{PackedSeqVec, SeqVec};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};
use std::sync::Arc;
use std::sync::mpsc;
use tempfile::{Builder as TempFileBuilder, TempPath};
use threadpool::ThreadPool;
use xz2::read::XzDecoder;
use zstd::stream::read::Decoder as ZstdDecoder;
use zstd::stream::write::Encoder as ZstdEncoder;

const BUFFER_TARGET: usize = 128 * 1024;
const ZSTD_LEVEL_FAST4: i32 = -4; // zstd "fast=4" mode for lower memory and higher speed
const SORT_RUN_TARGET_BYTES: usize = 256 * 1024 * 1024; // spill sort runs around 256 MiB
const MERGE_FANIN: usize = 8;
static MERGE_TMP_COUNTER: AtomicUsize = AtomicUsize::new(0);

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

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Parse FASTA files into superkmers grouped by minimizer"
)]
struct Args {
    #[arg(short = 'i', long = "input-fof")]
    input_fof: PathBuf,
    #[arg(short = 'o', long = "output-dir")]
    output_dir: PathBuf,
    #[arg(short = 'k', long = "kmer", default_value_t = 31)]
    k: usize,
    #[arg(short = 'm', long = "minimizer", default_value_t = 9)]
    m: usize,
    /// Number of partitions is 2^partition_power (default 1024 partitions)
    #[arg(short = 'P', long = "partition-power", default_value_t = 10)]
    partition_power: u32,
    #[arg(short = 't', long = "threads", default_value_t = 32)]
    threads: usize,
    /// Concurrency for the compaction phase (per-partition simplitigs)
    #[arg(long = "compaction-threads", default_value_t = num_cpus::get())]
    compaction_threads: usize,
    /// Optionally verify that all canonical k-mers are preserved with the correct dataset IDs
    #[arg(long = "verify-kmers", default_value_t = false)]
    verify_kmers: bool,
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
    ids: Vec<usize>,
    seq: Vec<u8>,
}

fn estimated_record_bytes(ids: &[usize], seq_len: usize) -> usize {
    ids.len() * std::mem::size_of::<usize>() + seq_len
}

fn parse_ids_from_header(header: &str) -> Result<Vec<usize>> {
    let Some(rest) = header.strip_prefix(">ids:") else {
        bail!("invalid simplitig header (expected >ids:): '{header}'");
    };
    let mut ids = Vec::new();
    for part in rest.split(',') {
        if part.is_empty() {
            continue;
        }
        let id: usize = part
            .parse()
            .with_context(|| format!("parse dataset id '{part}' from header '{header}'"))?;
        ids.push(id);
    }
    if ids.is_empty() {
        bail!("no dataset ids found in header '{header}'");
    }
    Ok(ids)
}

fn read_simplitig_record<R: BufRead + ?Sized>(
    reader: &mut R,
) -> Result<Option<SimplitigRecord>> {
    let mut header = String::new();
    if reader.read_line(&mut header)? == 0 {
        return Ok(None);
    }
    let mut seq_line = String::new();
    if reader.read_line(&mut seq_line)? == 0 {
        bail!("unexpected EOF after header '{header}'");
    }
    let header = header.trim_end_matches(['\r', '\n']);
    let seq = seq_line.trim_end_matches(['\r', '\n']).as_bytes().to_vec();
    let ids = parse_ids_from_header(header)?;
    Ok(Some(SimplitigRecord { ids, seq }))
}

fn write_simplitig_record<W: Write>(record: &SimplitigRecord, writer: &mut W) -> Result<()> {
    writer.write_all(b">ids:")?;
    for (idx, id) in record.ids.iter().enumerate() {
        if idx > 0 {
            writer.write_all(b",")?;
        }
        write!(writer, "{id}")?;
    }
    writer.write_all(b"\n")?;
    writer.write_all(&record.seq)?;
    writer.write_all(b"\n")?;
    Ok(())
}

fn simplitig_cmp(a: &SimplitigRecord, b: &SimplitigRecord) -> std::cmp::Ordering {
    match a.ids.cmp(&b.ids) {
        std::cmp::Ordering::Equal => a.seq.cmp(&b.seq),
        ord => ord,
    }
}

fn write_sorted_chunk<W: Write>(chunk: &mut [SimplitigRecord], writer: &mut W) -> Result<()> {
    chunk.sort_by(simplitig_cmp);
    for record in chunk.iter() {
        write_simplitig_record(record, writer)?;
    }
    Ok(())
}

struct RunReader {
    reader: BufReader<File>,
    path: TempPath,
}

impl RunReader {
    fn new(path: TempPath) -> Result<Self> {
        let path_ref: &Path = path.as_ref();
        let file = File::open(path_ref)
            .with_context(|| format!("open temporary run {}", path.display()))?;
        Ok(Self {
            reader: BufReader::new(file),
            path,
        })
    }

    fn next_record(&mut self) -> Result<Option<SimplitigRecord>> {
        read_simplitig_record(&mut self.reader)
            .with_context(|| format!("read from temporary run {}", self.path.display()))
    }
}

#[derive(Debug)]
struct HeapItem {
    run_idx: usize,
    record: SimplitigRecord,
}

impl PartialEq for HeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.record.ids == other.record.ids && self.record.seq == other.record.seq
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

fn merge_runs(
    mut readers: Vec<RunReader>,
    writer: &mut ZstdEncoder<'static, File>,
) -> Result<()> {
    let mut heap: BinaryHeap<Reverse<HeapItem>> = BinaryHeap::new();
    for (idx, reader) in readers.iter_mut().enumerate() {
        if let Some(record) = reader.next_record()? {
            heap.push(Reverse(HeapItem {
                run_idx: idx,
                record,
            }));
        }
    }

    while let Some(Reverse(item)) = heap.pop() {
        write_simplitig_record(&item.record, writer)?;
        if let Some(next) = readers[item.run_idx].next_record()? {
            heap.push(Reverse(HeapItem {
                run_idx: item.run_idx,
                record: next,
            }));
        }
    }
    Ok(())
}

fn flush_chunk_to_run(
    chunk: &mut Vec<SimplitigRecord>,
    chunk_bytes: &mut usize,
    run_dir: &Path,
    run_paths: &mut Vec<TempPath>,
) -> Result<()> {
    if chunk.is_empty() {
        return Ok(());
    }
    let mut tmp = TempFileBuilder::new()
        .prefix("simplitig-run-")
        .tempfile_in(run_dir)
        .with_context(|| format!("create temporary run in {}", run_dir.display()))?;
    {
        let mut writer = BufWriter::new(tmp.as_file_mut());
        write_sorted_chunk(chunk, &mut writer)?;
        writer.flush()?;
    }
    run_paths.push(tmp.into_temp_path());
    chunk.clear();
    *chunk_bytes = 0;
    Ok(())
}

fn assemble_simplitigs(
    mut kmer_map: HashMap<u64, Vec<u64>>,
    k: usize,
    mut sink: impl FnMut(Vec<u8>, Vec<u64>) -> Result<()>,
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

        sink(seq, ids)?;
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
) -> Result<()> {
    let map = build_kmer_map_from_partition(partition_path, k, dataset_count)?;
    let run_dir = output_path
        .parent()
        .map(Path::to_path_buf)
        .unwrap_or_else(|| PathBuf::from("."));
    let mut chunk: Vec<SimplitigRecord> = Vec::new();
    let mut chunk_bytes = 0usize;
    let mut run_paths: Vec<TempPath> = Vec::new();

    assemble_simplitigs(map, k, |seq, ids_bits| {
        let ids = ids_from_bitset(&ids_bits);
        chunk_bytes += estimated_record_bytes(&ids, seq.len());
        chunk.push(SimplitigRecord { ids, seq });
        if chunk_bytes >= SORT_RUN_TARGET_BYTES {
            flush_chunk_to_run(&mut chunk, &mut chunk_bytes, &run_dir, &mut run_paths)?;
        }
        Ok(())
    })?;

    let file = File::create(output_path)
        .with_context(|| format!("create simplitig output {}", output_path.display()))?;
    let mut encoder = ZstdEncoder::new(file, ZSTD_LEVEL_FAST4)
        .with_context(|| format!("build zstd encoder for {}", output_path.display()))?;

    if run_paths.is_empty() {
        write_sorted_chunk(&mut chunk, &mut encoder)?;
    } else {
        if !chunk.is_empty() {
            flush_chunk_to_run(&mut chunk, &mut chunk_bytes, &run_dir, &mut run_paths)?;
        }
        let readers = run_paths
            .into_iter()
            .map(RunReader::new)
            .collect::<Result<Vec<_>>>()?;
        merge_runs(readers, &mut encoder)?;
    }

    encoder
        .finish()
        .with_context(|| format!("finalize simplitig output {}", output_path.display()))?;
    Ok(())
}

struct PartitionReader {
    reader: Box<dyn BufRead + Send>,
    path: PathBuf,
}

impl PartitionReader {
    fn new(path: PathBuf) -> Result<Self> {
        let file = File::open(&path)
            .with_context(|| format!("open partition simplitigs {}", path.display()))?;
        let decoder = ZstdDecoder::new(file)
            .with_context(|| format!("build zstd decoder for {}", path.display()))?;
        Ok(Self {
            reader: Box::new(BufReader::new(decoder)),
            path,
        })
    }

    fn next_record(&mut self) -> Result<Option<SimplitigRecord>> {
        read_simplitig_record(self.reader.as_mut())
            .with_context(|| format!("read simplitig from {}", self.path.display()))
    }
}

fn merge_sorted_partitions(partition_paths: &[PathBuf], output_path: &Path) -> Result<()> {
    let mut readers: Vec<PartitionReader> = Vec::new();
    for path in partition_paths {
        readers.push(PartitionReader::new(path.clone())?);
    }

    let file = File::create(output_path)
        .with_context(|| format!("create simplitig output {}", output_path.display()))?;
    let mut encoder = ZstdEncoder::new(file, ZSTD_LEVEL_FAST4)
        .with_context(|| format!("build zstd encoder for {}", output_path.display()))?;

    let mut heap: BinaryHeap<Reverse<HeapItem>> = BinaryHeap::new();
    for (idx, reader) in readers.iter_mut().enumerate() {
        if let Some(record) = reader.next_record()? {
            heap.push(Reverse(HeapItem {
                run_idx: idx,
                record,
            }));
        } else {
            // remove empty partition outputs immediately
            let _ = fs::remove_file(&reader.path);
        }
    }

    while let Some(Reverse(item)) = heap.pop() {
        write_simplitig_record(&item.record, &mut encoder)?;
        if let Some(next) = readers[item.run_idx].next_record()? {
            heap.push(Reverse(HeapItem {
                run_idx: item.run_idx,
                record: next,
            }));
        } else {
            let path = readers[item.run_idx].path.clone();
            let _ = fs::remove_file(&path);
        }
    }

    encoder
        .finish()
        .with_context(|| format!("finalize simplitig output {}", output_path.display()))?;
    Ok(())
}

fn parallel_merge_sorted_partitions(
    partition_paths: Vec<PathBuf>,
    output_path: &Path,
    threads: usize,
) -> Result<()> {
    if partition_paths.is_empty() {
        let file = File::create(output_path)
            .with_context(|| format!("create simplitig output {}", output_path.display()))?;
        let encoder = ZstdEncoder::new(file, ZSTD_LEVEL_FAST4)
            .with_context(|| format!("build zstd encoder for {}", output_path.display()))?;
        encoder
            .finish()
            .with_context(|| format!("finalize simplitig output {}", output_path.display()))?;
        return Ok(());
    }
    if partition_paths.len() == 1 {
        let only = partition_paths[0].clone();
        match fs::rename(&only, output_path) {
            Ok(_) => return Ok(()),
            Err(_) => {
                fs::copy(&only, output_path).with_context(|| {
                    format!(
                        "copy {} to {}",
                        only.display(),
                        output_path.display()
                    )
                })?;
                let _ = fs::remove_file(&only);
                return Ok(());
            }
        }
    }

    let temp_dir = output_path
        .parent()
        .map(Path::to_path_buf)
        .unwrap_or_else(|| PathBuf::from("."));

    let next_merge_path = |dir: &Path| -> PathBuf {
        loop {
            let id = MERGE_TMP_COUNTER.fetch_add(1, AtomicOrdering::Relaxed);
            let candidate = dir.join(format!("simplitig-merge-{id}.fa.zst"));
            if !candidate.exists() {
                return candidate;
            }
        }
    };

    let mut current = partition_paths;
    while current.len() > 1 {
        let mut next = Vec::new();
        let (tx, rx) = mpsc::channel();
        let batch_threads = threads.max(1).min(current.len());
        let pool = ThreadPool::new(batch_threads);
        for chunk in current.chunks(MERGE_FANIN) {
            let paths = chunk.to_vec();
            let tx = tx.clone();
            let temp_dir = temp_dir.clone();
            pool.execute(move || {
                let result: Result<PathBuf> = (|| {
                    if paths.len() == 1 {
                        return Ok(paths[0].clone());
                    }
                    let out_path = next_merge_path(&temp_dir);
                    merge_sorted_partitions(&paths, &out_path)?;
                    Ok(out_path)
                })();
                let _ = tx.send(result);
            });
        }
        drop(tx);
        for result in rx {
            next.push(result?);
        }
        pool.join();
        current = next;
    }

    let final_path = current.pop().unwrap();
    match fs::rename(&final_path, output_path) {
        Ok(_) => {}
        Err(_) => {
            fs::copy(&final_path, output_path).with_context(|| {
                format!(
                    "copy {} to {}",
                    final_path.display(),
                    output_path.display()
                )
            })?;
            let _ = fs::remove_file(&final_path);
        }
    }
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
    if k < m {
        bail!("k-mer length ({}) must be >= minimizer length ({})", k, m);
    }
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

    let pool = ThreadPool::new(threads);
    for (idx, file_path) in file_paths.iter().enumerate() {
        let encoders = Arc::clone(&encoders);
        let partitions = partitions;
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

    // Second phase: per-partition compaction into simplitigs with identical ID sets.
    println!("Starting per-partition simplitig compaction...");
    let output_simplitigs = output_dir.join("simplitigs.fa.zst");
    let partition_outputs = Arc::new(Mutex::new(Vec::new()));
    let compaction_pool = ThreadPool::new(threads);
    for (idx, pw) in encoders.iter().enumerate() {
        let partition_path = pw.path.clone();
        let partition_outputs = Arc::clone(&partition_outputs);
        let part_output = output_dir.join(format!("simplitigs-part-{idx}.fa.zst"));
        compaction_pool.execute(move || {
            let result = write_partition_simplitigs(
                &partition_path,
                k,
                dataset_count,
                &part_output,
            );
            match result {
                Ok(_) => partition_outputs.lock().push(part_output),
                Err(e) => eprintln!(
                    "Failed to assemble simplitigs for {}: {}",
                    partition_path.display(),
                    e
                ),
            }
        });
    }
    compaction_pool.join();
    let mut partition_paths = partition_outputs.lock().clone();
    partition_paths.sort();
    parallel_merge_sorted_partitions(partition_paths, &output_simplitigs, threads)?;

    if verify_kmers {
        println!("Verifying k-mer preservation...");
        let output_map = build_output_kmer_map(&output_simplitigs, k, dataset_count)?;
        let input_map = build_kmer_map_from_inputs(&file_paths, k)?;
        verify_kmer_maps(&input_map, &output_map, k)?;
        println!(
            "Verified {} canonical k-mers; all input k-mers were found in the output.",
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

pub fn run_parser(
    input_fof: PathBuf,
    output_dir: PathBuf,
    k: usize,
    m: usize,
    partition_power: u32,
    threads: usize,
    compaction_threads: usize,
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
    let partition_outputs = Arc::new(Mutex::new(Vec::new()));
    let compaction_pool = ThreadPool::new(compaction_threads);
    for (idx, pw) in encoders.iter().enumerate() {
        let partition_path = pw.path.clone();
        let partition_outputs = Arc::clone(&partition_outputs);
        let part_output = output_dir
            .join(format!("simplitigs-part-{idx}.fa.zst"));
        compaction_pool.execute(move || {
            let result = write_partition_simplitigs(
                &partition_path,
                k,
                dataset_count,
                &part_output,
            );
            match result {
                Ok(_) => partition_outputs.lock().push(part_output),
                Err(e) => eprintln!(
                    "Failed to assemble simplitigs for {}: {}",
                    partition_path.display(),
                    e
                ),
            }
        });
    }
    compaction_pool.join();
    let mut partition_paths = partition_outputs.lock().clone();
    partition_paths.sort();
    parallel_merge_sorted_partitions(
        partition_paths,
        &output_simplitigs,
        compaction_threads,
    )?;

    if verify_kmers {
        println!("Verifying k-mer preservation...");
        let output_map = build_output_kmer_map(&output_simplitigs, k, dataset_count)?;
        let input_map = build_kmer_map_from_inputs(&file_paths, k)?;
        verify_kmer_maps(&input_map, &output_map, k)?;
        println!(
            "Verified {} canonical k-mers; all input k-mers were found in the output.",
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

/*fn main() -> Result<()> {
    let args = Args::parse();
    run(
        args.input_fof,
        args.output_dir,
        args.k,
        args.m,
        args.partition_power,
        args.threads,
        args.compaction_threads,
        args.verify_kmers,
    )
}*/