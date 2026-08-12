use ggcat_api::{
    install_plain_fasta_output_callback, ColorIndexType, ExtraElaboration, GGCATConfig,
    GGCATInstance, GeneralSequenceBlockData, PlainFastaOutputRecord,
};
use ggcat_colors::colors_manager::ColorMapReader;
use ggcat_colors::storage::deserializer::ColorsDeserializer;
use ggcat_colors::DefaultColorsSerializer;
use rayon::prelude::*;
use rayon::slice::ParallelSliceMut;
use std::cmp::Ordering as CmpOrdering;
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Result, Write};
use std::path::{Path, PathBuf};
use std::sync::{mpsc, Arc, Mutex};
use std::thread;
use std::time::Instant;
use tempfile::Builder as TempBuilder;
use zstd::Encoder;

use crate::records::{SimplitigBatch, SimplitigRecord};
use crate::utils::{Convert, Converter};

const IO_BUFFER_CAPACITY: usize = 16 * 1024 * 1024;
const ENCODED_SEQ_BUFFER_TARGET: usize = 4 * 1024 * 1024;
const ID_CID_SPILL_BUFFER_CAPACITY: usize = 1024 * 1024;
const PAR_SORT_THRESHOLD: usize = 200_000;
const COLOR_RECORD_BATCH_SIZE: usize = 65_536;
const COLOR_CHUNK_TARGET_BYTES: usize = 256 * 1024 * 1024;
const GROUP_SORT_SPILL_BYTES: usize = 128 * 1024 * 1024;
const GROUP_WORKER_QUEUE_DEPTH: usize = 8;
const GROUP_WORK_BATCH_GROUPS: usize = 2048;
const GROUP_WORK_BATCH_SEQ_BYTES: usize = 256 * 1024 * 1024;
const GROUP_IN_MEMORY_MAX_SEQ_BYTES: usize = 2 * 1024 * 1024;
const SPILL_CID_FLUSH_BYTES: usize = 64 * 1024;
const PRODUCER_PARSE_BATCH_ENTRIES: usize = 4_096;
const PRODUCER_PARSE_BATCH_BYTES: usize = 128 * 1024 * 1024;
const SUBSET_QUERY_CHUNK_SIZE: usize = 500_000;
const BUCKET_SIZES_MAGIC: &[u8; 4] = b"KSB2";
const POSITIONS_MAGIC: &[u8; 4] = b"KPS2";
const ID_TO_CID_MAGIC: &[u8; 4] = b"KIC2";
const BUCKET_SIZE_BLOCK_MAX_GROUPS: usize = 8_192;
const BUCKET_SIZE_BLOCK_MAX_UNCOMPRESSED_BYTES: usize = 8 * 1024 * 1024;

pub fn compress(
    output_dir: &String,
    input_fof: &String,
    threads: usize,
    k: usize,
    m: usize,
    partition_power: u32,
    verify_kmers: bool,
    skip_sort: bool,
    use_unitigs: bool,
    use_matchtigs: bool,
    use_eulertigs: bool,
) -> Result<()> {
    compress_with_ggcat(
        output_dir,
        input_fof,
        threads,
        k,
        m,
        partition_power,
        verify_kmers,
        skip_sort,
        use_unitigs,
        use_matchtigs,
        use_eulertigs,
        GgcatCompressionConfig::default(),
    )
}

#[derive(Clone, Debug)]
pub struct GgcatCompressionConfig {
    pub memory_gb: usize,
    pub temp_dir: String,
}

impl Default for GgcatCompressionConfig {
    fn default() -> Self {
        Self {
            memory_gb: 8,
            temp_dir: String::new(),
        }
    }
}

#[derive(Debug)]
struct SortedColorRecord {
    subset: ColorIndexType,
    seq: Vec<u8>,
}

#[derive(Debug)]
struct ChunkHeapItem {
    chunk_index: usize,
    record: SortedColorRecord,
}

impl Eq for ChunkHeapItem {}

impl PartialEq for ChunkHeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.record.subset == other.record.subset && self.chunk_index == other.chunk_index
    }
}

impl Ord for ChunkHeapItem {
    fn cmp(&self, other: &Self) -> CmpOrdering {
        other
            .record
            .subset
            .cmp(&self.record.subset)
            .then_with(|| other.chunk_index.cmp(&self.chunk_index))
    }
}

impl PartialOrd for ChunkHeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<CmpOrdering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug)]
struct SeqHeapItem {
    run_idx: usize,
    seq: Vec<u8>,
}

impl Eq for SeqHeapItem {}

impl PartialEq for SeqHeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.seq == other.seq
    }
}

impl Ord for SeqHeapItem {
    fn cmp(&self, other: &Self) -> CmpOrdering {
        other
            .seq
            .len()
            .cmp(&self.seq.len())
            .then_with(|| other.seq.cmp(&self.seq))
    }
}

impl PartialOrd for SeqHeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<CmpOrdering> {
        Some(self.cmp(other))
    }
}

struct SeqRunReader {
    reader: BufReader<File>,
    buf: Vec<u8>,
}

impl SeqRunReader {
    fn new(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        Ok(Self {
            reader: BufReader::new(file),
            buf: Vec::new(),
        })
    }

    fn next_seq(&mut self) -> Result<Option<Vec<u8>>> {
        self.buf.clear();
        let n = self.reader.read_until(b'\n', &mut self.buf)?;
        if n == 0 {
            return Ok(None);
        }
        while self.buf.last() == Some(&b'\n') || self.buf.last() == Some(&b'\r') {
            self.buf.pop();
        }
        Ok(Some(self.buf.clone()))
    }
}

fn read_input_fof_filenames(input_fof: &str) -> Result<Vec<String>> {
    let input_fof_reader = BufReader::new(File::open(input_fof)?);
    let mut filenames = Vec::new();
    for line_result in input_fof_reader.lines() {
        let line = line_result?;
        let trimmed = line.trim();
        if !trimmed.is_empty() {
            filenames.push(trimmed.to_string());
        }
    }
    Ok(filenames)
}

pub(crate) fn write_filenames_id_offsets(
    output_dir: &str,
    filenames: &[String],
    id_cid_line_sizes: &[usize],
) -> Result<()> {
    if filenames.len() != id_cid_line_sizes.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "filenames count ({}) differs from id->cid offsets count ({})",
                filenames.len(),
                id_cid_line_sizes.len()
            ),
        ));
    }

    let mut fof_id = BufWriter::new(File::create(output_dir.to_owned() + "filenames_id.txt")?);
    for (filename, offset) in filenames.iter().zip(id_cid_line_sizes.iter()) {
        fof_id.write_all(format!("{filename}:{offset}\n").as_bytes())?;
    }
    fof_id.flush()?;
    Ok(())
}

fn to_io_err(context: &str, err: impl std::fmt::Display) -> io::Error {
    io::Error::other(format!("{context}: {err}"))
}

#[derive(Clone, Copy, Debug)]
struct PhaseTiming {
    wall_sec: f64,
    cpu_sec: Option<f64>,
}

#[derive(Debug)]
struct PhaseTimer {
    wall_start: Instant,
    cpu_start: Option<f64>,
}

impl PhaseTimer {
    fn start() -> Self {
        Self {
            wall_start: Instant::now(),
            cpu_start: read_process_cpu_seconds(),
        }
    }

    fn finish(&self) -> PhaseTiming {
        let wall_sec = self.wall_start.elapsed().as_secs_f64();
        let cpu_sec = self
            .cpu_start
            .and_then(|start| read_process_cpu_seconds().map(|now| (now - start).max(0.0)));
        PhaseTiming { wall_sec, cpu_sec }
    }
}

#[derive(Debug, Default)]
struct TimingAccumulator {
    wall_sec: f64,
    cpu_sec: Option<f64>,
}

impl TimingAccumulator {
    fn add(&mut self, timing: PhaseTiming) {
        self.wall_sec += timing.wall_sec;
        self.cpu_sec = match (self.cpu_sec, timing.cpu_sec) {
            (Some(acc), Some(curr)) => Some(acc + curr),
            _ => None,
        };
    }
}

fn read_process_cpu_seconds() -> Option<f64> {
    let stat = fs::read_to_string("/proc/self/stat").ok()?;
    let rparen = stat.rfind(')')?;
    let rest = stat.get((rparen + 2)..)?;
    let fields: Vec<&str> = rest.split_whitespace().collect();
    if fields.len() <= 12 {
        return None;
    }

    let utime_ticks = fields[11].parse::<u64>().ok()? as f64;
    let stime_ticks = fields[12].parse::<u64>().ok()? as f64;
    let clk_tck = unsafe { libc::sysconf(libc::_SC_CLK_TCK) };
    if clk_tck <= 0 {
        return None;
    }

    Some((utime_ticks + stime_ticks) / (clk_tck as f64))
}

fn log_phase_timing(phase: &str, timing: PhaseTiming) {
    if let Some(cpu_sec) = timing.cpu_sec {
        let avg_cores = if timing.wall_sec > 0.0 {
            cpu_sec / timing.wall_sec
        } else {
            0.0
        };
        println!(
            "[phase-timing] phase={} wall_s={:.3} cpu_s={:.3} avg_cores={:.3}",
            phase, timing.wall_sec, cpu_sec, avg_cores
        );
    } else {
        println!(
            "[phase-timing] phase={} wall_s={:.3} cpu_s=NA avg_cores=NA",
            phase, timing.wall_sec
        );
    }
}

fn write_varint_u64(mut value: u64, out: &mut Vec<u8>) {
    while value >= 0x80 {
        out.push((value as u8 & 0x7f) | 0x80);
        value >>= 7;
    }
    out.push(value as u8);
}

fn write_varint_u64_to_writer(mut value: u64, mut out: impl Write) -> io::Result<()> {
    let mut buf = [0u8; 10];
    let mut len = 0usize;
    while value >= 0x80 {
        buf[len] = (value as u8 & 0x7f) | 0x80;
        value >>= 7;
        len += 1;
    }
    buf[len] = value as u8;
    len += 1;
    out.write_all(&buf[..len])
}

#[derive(Debug, Default)]
struct ChunkFlushMetrics {
    calls: usize,
    records: usize,
    timing: TimingAccumulator,
}

#[derive(Debug, Default)]
struct WorkerStats {
    groups: usize,
    timing: TimingAccumulator,
    max_wall_sec: f64,
    max_cpu_sec: Option<f64>,
}

impl WorkerStats {
    fn record_group(&mut self, wall_sec: f64, cpu_sec: Option<f64>) {
        self.groups += 1;
        self.timing.add(PhaseTiming { wall_sec, cpu_sec });
        if wall_sec > self.max_wall_sec {
            self.max_wall_sec = wall_sec;
        }
        self.max_cpu_sec = match (self.max_cpu_sec, cpu_sec) {
            (Some(acc), Some(curr)) => Some(acc.max(curr)),
            (None, Some(curr)) => Some(curr),
            _ => self.max_cpu_sec,
        };
    }
}

#[derive(Debug, Default)]
struct CommitStats {
    groups: usize,
    tigs_bytes: u64,
    sizes_bytes: u64,
    timing: TimingAccumulator,
}

fn ggcat_extra_elaboration(
    use_unitigs: bool,
    use_matchtigs: bool,
    use_eulertigs: bool,
) -> ExtraElaboration {
    if use_unitigs {
        ExtraElaboration::None
    } else if use_matchtigs {
        ExtraElaboration::GreedyMatchtigs
    } else if use_eulertigs {
        ExtraElaboration::FastEulertigs
    } else {
        // For colored-record export we need color-safe output to preserve per-file k-mers.
        ExtraElaboration::None
    }
}

fn write_sorted_color_record(
    writer: &mut BufWriter<File>,
    record: &SortedColorRecord,
) -> Result<()> {
    writer.write_all(&record.subset.to_le_bytes())?;
    writer.write_all(&(record.seq.len() as u32).to_le_bytes())?;
    writer.write_all(&record.seq)?;
    Ok(())
}

fn read_sorted_color_record(reader: &mut BufReader<File>) -> Result<Option<SortedColorRecord>> {
    let mut len_buf = [0u8; 4];
    match reader.read_exact(&mut len_buf) {
        Ok(()) => {}
        Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(err) => return Err(err),
    }
    let subset = u32::from_le_bytes(len_buf);

    reader.read_exact(&mut len_buf)?;
    let seq_len = u32::from_le_bytes(len_buf) as usize;

    let mut seq = vec![0u8; seq_len];
    reader.read_exact(&mut seq)?;

    Ok(Some(SortedColorRecord { subset, seq }))
}

fn read_color_record(
    reader: &mut BufReader<File>,
) -> Result<Option<(Vec<ColorIndexType>, Vec<u8>)>> {
    let mut len_buf = [0u8; 4];
    match reader.read_exact(&mut len_buf) {
        Ok(()) => {}
        Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(err) => return Err(err),
    }
    let subset_count = u32::from_le_bytes(len_buf) as usize;

    reader.read_exact(&mut len_buf)?;
    let seq_len = u32::from_le_bytes(len_buf) as usize;

    let mut subsets = vec![0u32; subset_count];
    for subset in &mut subsets {
        reader.read_exact(&mut len_buf)?;
        *subset = u32::from_le_bytes(len_buf);
    }

    let mut seq = vec![0u8; seq_len];
    reader.read_exact(&mut seq)?;

    Ok(Some((subsets, seq)))
}

fn parse_color_runs_from_header(header: &[u8]) -> Result<Vec<(ColorIndexType, usize)>> {
    let mut runs = Vec::new();
    for token in header.split(|b| *b == b' ') {
        if token.len() < 4 || token[0] != b'C' || token[1] != b':' {
            continue;
        }
        let rest = &token[2..];
        let Some(colon_pos) = rest.iter().position(|b| *b == b':') else {
            continue;
        };
        if colon_pos == 0 || colon_pos + 1 >= rest.len() {
            continue;
        }
        let color_hex = std::str::from_utf8(&rest[..colon_pos]).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "invalid color tag '{}' in header '{}': {err}",
                    String::from_utf8_lossy(token),
                    String::from_utf8_lossy(header)
                ),
            )
        })?;
        let count_text = std::str::from_utf8(&rest[(colon_pos + 1)..]).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "invalid color count tag '{}' in header '{}': {err}",
                    String::from_utf8_lossy(token),
                    String::from_utf8_lossy(header)
                ),
            )
        })?;

        let subset = ColorIndexType::from_str_radix(color_hex, 16).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "invalid color subset '{}' in header '{}': {err}",
                    color_hex,
                    String::from_utf8_lossy(header)
                ),
            )
        })?;
        let count = count_text.parse::<usize>().map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "invalid color run count '{}' in header '{}': {err}",
                    count_text,
                    String::from_utf8_lossy(header)
                ),
            )
        })?;
        if count > 0 {
            runs.push((subset, count));
        }
    }
    Ok(runs)
}

fn flush_sorted_chunk(
    records: &mut Vec<SortedColorRecord>,
    chunk_files: &mut Vec<PathBuf>,
    chunk_dir: &Path,
) -> Result<()> {
    if records.is_empty() {
        return Ok(());
    }

    if records.len() >= PAR_SORT_THRESHOLD {
        records.par_sort_unstable_by(|left, right| left.subset.cmp(&right.subset));
    } else {
        records.sort_unstable_by(|left, right| left.subset.cmp(&right.subset));
    }

    let chunk_path = chunk_dir.join(format!("chunk_{:08}.bin", chunk_files.len()));
    let mut chunk_writer = BufWriter::with_capacity(IO_BUFFER_CAPACITY, File::create(&chunk_path)?);
    for record in records.iter() {
        write_sorted_color_record(&mut chunk_writer, record)?;
    }
    chunk_writer.flush()?;
    chunk_files.push(chunk_path);
    records.clear();
    Ok(())
}

fn flush_sorted_chunk_timed(
    records: &mut Vec<SortedColorRecord>,
    chunk_files: &mut Vec<PathBuf>,
    chunk_dir: &Path,
    metrics: &mut ChunkFlushMetrics,
) -> Result<()> {
    if records.is_empty() {
        return Ok(());
    }
    let records_count = records.len();
    let timer = PhaseTimer::start();
    flush_sorted_chunk(records, chunk_files, chunk_dir)?;
    metrics.calls += 1;
    metrics.records += records_count;
    metrics.timing.add(timer.finish());
    Ok(())
}

fn parse_entry_color_segments(
    header: &[u8],
    seq: &[u8],
    k: usize,
) -> Result<Vec<SortedColorRecord>> {
    let runs = parse_color_runs_from_header(header)?;
    if runs.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "ggcat emitted colored sequence without color runs: '{}'",
                String::from_utf8_lossy(header)
            ),
        ));
    }

    if seq.len() < k {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "ggcat sequence shorter than k (len={}, k={}) for '{}'",
                seq.len(),
                k,
                String::from_utf8_lossy(header)
            ),
        ));
    }

    let seq_kmers = seq.len() - k + 1;
    let mut total_run_kmers = 0usize;
    let mut merged_runs: Vec<(ColorIndexType, usize)> = Vec::with_capacity(runs.len());

    for (subset, count) in runs {
        total_run_kmers = total_run_kmers.checked_add(count).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "overflow while summing GGCAT color-run lengths",
            )
        })?;

        if let Some((last_subset, last_count)) = merged_runs.last_mut() {
            if *last_subset == subset {
                *last_count += count;
                continue;
            }
        }
        merged_runs.push((subset, count));
    }

    if total_run_kmers != seq_kmers {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "GGCAT header color-run lengths mismatch sequence k-mers for '{}': runs={}, expected={}",
                String::from_utf8_lossy(header),
                total_run_kmers,
                seq_kmers
            ),
        ));
    }

    let mut segments = Vec::with_capacity(merged_runs.len());
    let mut run_start = 0usize;
    for (subset, count) in merged_runs {
        let run_end = run_start + count;
        if run_end > seq_kmers {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "GGCAT color run exceeds sequence k-mers for '{}': end={}, total={}",
                    String::from_utf8_lossy(header),
                    run_end,
                    seq_kmers
                ),
            ));
        }
        segments.push(SortedColorRecord {
            subset,
            seq: seq[run_start..(run_end + k - 1)].to_vec(),
        });
        run_start = run_end;
    }

    if run_start != seq_kmers {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "GGCAT color runs do not fully cover sequence k-mers for '{}': covered={}, total={}",
                String::from_utf8_lossy(header),
                run_start,
                seq_kmers
            ),
        ));
    }

    Ok(segments)
}

fn resolve_subsets_to_dataset_ids(
    instance: &GGCATInstance,
    colormap_file: &Path,
    subsets: &HashSet<ColorIndexType>,
    color_index_to_dataset_index: &[usize],
    dataset_count: usize,
) -> Result<HashMap<ColorIndexType, Arc<Vec<u32>>>> {
    if subsets.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "no color subsets were parsed from GGCAT output",
        ));
    }

    let mut subset_list = subsets.iter().copied().collect::<Vec<_>>();
    subset_list.sort_unstable();

    let resolved = Mutex::new(HashMap::<ColorIndexType, Arc<Vec<u32>>>::with_capacity(
        subset_list.len(),
    ));
    let callback_err = Mutex::new(None::<io::Error>);

    for chunk in subset_list.chunks(SUBSET_QUERY_CHUNK_SIZE.max(1)) {
        instance
            .query_colormap(
                colormap_file.to_path_buf(),
                chunk.to_vec(),
                false,
                |subset, colors| {
                    if callback_err
                        .lock()
                        .expect("callback error lock poisoned")
                        .is_some()
                    {
                        return;
                    }

                    let mut dataset_ids = Vec::with_capacity(colors.len());
                    for &color in colors {
                        let color_idx = color as usize;
                        let Some(&dataset_idx) = color_index_to_dataset_index.get(color_idx) else {
                            *callback_err.lock().expect("callback error lock poisoned") =
                                Some(io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    format!("ggcat color index {} out of range", color_idx),
                                ));
                            return;
                        };
                        if dataset_idx >= dataset_count {
                            *callback_err.lock().expect("callback error lock poisoned") =
                                Some(io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    format!("dataset index {} out of range", dataset_idx),
                                ));
                            return;
                        }
                        dataset_ids.push((dataset_idx + 1) as u32);
                    }

                    if dataset_ids.is_empty() {
                        *callback_err.lock().expect("callback error lock poisoned") =
                            Some(io::Error::new(
                                io::ErrorKind::InvalidData,
                                "GGCAT color subset resolved to empty dataset set",
                            ));
                        return;
                    }
                    dataset_ids.sort_unstable();

                    resolved
                        .lock()
                        .expect("resolved-subset lock poisoned")
                        .insert(subset, Arc::new(dataset_ids));
                },
            )
            .map_err(|e| to_io_err("query ggcat colormap", e))?;

        if let Some(err) = callback_err
            .lock()
            .expect("callback error lock poisoned")
            .take()
        {
            return Err(err);
        }
    }

    let map = resolved
        .into_inner()
        .expect("resolved-subset lock poisoned");
    if map.len() != subsets.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "resolved {} color subsets but expected {}",
                map.len(),
                subsets.len()
            ),
        ));
    }
    Ok(map)
}

fn sort_and_spill_sequence_run(
    seqs: &mut Vec<Vec<u8>>,
    run_paths: &mut Vec<PathBuf>,
    spill_dir: &Path,
) -> Result<()> {
    if seqs.is_empty() {
        return Ok(());
    }

    if seqs.len() >= PAR_SORT_THRESHOLD {
        seqs.par_sort_unstable_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
    } else {
        seqs.sort_unstable_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
    }

    let mut tmp = TempBuilder::new()
        .prefix("kloe-group-run-")
        .suffix(".seq")
        .tempfile_in(spill_dir)?;
    {
        let mut writer = BufWriter::with_capacity(IO_BUFFER_CAPACITY, tmp.as_file_mut());
        for seq in seqs.iter() {
            writer.write_all(seq)?;
            writer.write_all(b"\n")?;
        }
        writer.flush()?;
    }

    let run_path = tmp
        .into_temp_path()
        .keep()
        .map_err(|err| io::Error::other(format!("persist group spill file: {err}")))?;
    run_paths.push(run_path);
    seqs.clear();
    Ok(())
}

fn emit_sorted_group_sequences(
    seqs: &mut Vec<Vec<u8>>,
    run_paths: &mut Vec<PathBuf>,
    spill_dir: &Path,
    mut on_seq: impl FnMut(&[u8]) -> Result<()>,
) -> Result<()> {
    if run_paths.is_empty() {
        if seqs.len() >= PAR_SORT_THRESHOLD {
            seqs.par_sort_unstable_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
        } else {
            seqs.sort_unstable_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
        }
        for seq in seqs.iter() {
            on_seq(seq)?;
        }
        seqs.clear();
        return Ok(());
    }

    if !seqs.is_empty() {
        sort_and_spill_sequence_run(seqs, run_paths, spill_dir)?;
    }

    let mut readers = Vec::with_capacity(run_paths.len());
    for run in run_paths.iter() {
        readers.push(SeqRunReader::new(run)?);
    }

    let mut heap = BinaryHeap::new();
    for (run_idx, reader) in readers.iter_mut().enumerate() {
        if let Some(seq) = reader.next_seq()? {
            heap.push(SeqHeapItem { run_idx, seq });
        }
    }

    while let Some(item) = heap.pop() {
        on_seq(&item.seq)?;
        if let Some(next) = readers[item.run_idx].next_seq()? {
            heap.push(SeqHeapItem {
                run_idx: item.run_idx,
                seq: next,
            });
        }
    }

    for run_path in run_paths.drain(..) {
        let _ = fs::remove_file(run_path);
    }
    Ok(())
}

struct StreamWriterState {
    omni_file: BufWriter<File>,
    size_file: BufWriter<File>,
    spill_writers: Vec<BufWriter<File>>,
    spill_cid_buffers: Vec<Vec<u8>>,
    size_block_groups: usize,
    size_block_offsets: Vec<u32>,
    size_block_uncompressed: Vec<u8>,
    pos_nb_unitig: Vec<(u64, u64)>,
    prev_tigs_size: u64,
    prev_bucket_pos: u64,
    cid: usize,
    encoded_seq_buffer: Vec<u8>,
}

impl StreamWriterState {
    fn new(
        unitigs_file_path: String,
        output_dir: &str,
        nb_files: usize,
    ) -> Result<(Self, Vec<PathBuf>, PathBuf)> {
        let mut spill_paths = Vec::with_capacity(nb_files);
        let mut spill_writers = Vec::with_capacity(nb_files);
        let mut spill_cid_buffers = Vec::with_capacity(nb_files);
        let spill_dir = create_id_cid_spill_dir(output_dir)?;
        for id in 0..nb_files {
            let path = spill_dir.join(format!("id_{id}.cids.bin"));
            let writer =
                BufWriter::with_capacity(ID_CID_SPILL_BUFFER_CAPACITY, File::create(&path)?);
            spill_paths.push(path);
            spill_writers.push(writer);
            spill_cid_buffers.push(Vec::with_capacity(SPILL_CID_FLUSH_BYTES));
        }

        Ok((
            Self {
                omni_file: BufWriter::with_capacity(
                    IO_BUFFER_CAPACITY,
                    File::create(unitigs_file_path)?,
                ),
                size_file: {
                    let mut out = BufWriter::with_capacity(
                        IO_BUFFER_CAPACITY,
                        File::create(output_dir.to_owned() + "bucket_sizes.txt")?,
                    );
                    out.write_all(BUCKET_SIZES_MAGIC)?;
                    out
                },
                spill_writers,
                spill_cid_buffers,
                size_block_groups: 0,
                size_block_offsets: vec![0],
                size_block_uncompressed: Vec::new(),
                pos_nb_unitig: vec![(0, 0)],
                prev_tigs_size: 0,
                prev_bucket_pos: 0,
                cid: 0,
                encoded_seq_buffer: Vec::with_capacity(ENCODED_SEQ_BUFFER_TARGET),
            },
            spill_paths,
            spill_dir,
        ))
    }

    fn flush_size_block(&mut self) -> Result<()> {
        if self.size_block_groups == 0 {
            return Ok(());
        }

        let mut compressed = Vec::new();
        {
            let mut encoder = Encoder::new(&mut compressed, 1)?;
            encoder.write_all(&self.size_block_uncompressed)?;
            encoder.finish()?;
        }

        self.size_file
            .write_all(&(self.size_block_groups as u32).to_le_bytes())?;
        self.size_file
            .write_all(&(compressed.len() as u64).to_le_bytes())?;
        for &offset in &self.size_block_offsets {
            self.size_file.write_all(&offset.to_le_bytes())?;
        }
        self.size_file.write_all(&compressed)?;

        self.size_block_groups = 0;
        self.size_block_offsets.clear();
        self.size_block_offsets.push(0);
        self.size_block_uncompressed.clear();
        Ok(())
    }

    fn append_bucket_group_sizes_payload(&mut self, payload: &[u8]) -> Result<()> {
        self.size_block_uncompressed.extend_from_slice(payload);
        let next_offset = u32::try_from(self.size_block_uncompressed.len()).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "bucket sizes block exceeds 4GiB uncompressed payload",
            )
        })?;
        self.size_block_offsets.push(next_offset);
        self.size_block_groups += 1;

        if self.size_block_groups >= BUCKET_SIZE_BLOCK_MAX_GROUPS
            || self.size_block_uncompressed.len() >= BUCKET_SIZE_BLOCK_MAX_UNCOMPRESSED_BYTES
        {
            self.flush_size_block()?;
        }
        Ok(())
    }

    fn append_cid_for_dataset(&mut self, dataset_id_zero_based: usize, cid: usize) -> Result<()> {
        let buf = &mut self.spill_cid_buffers[dataset_id_zero_based];
        buf.extend_from_slice(&(cid as u64).to_le_bytes());
        if buf.len() >= SPILL_CID_FLUSH_BYTES {
            self.spill_writers[dataset_id_zero_based].write_all(buf)?;
            buf.clear();
        }
        Ok(())
    }

    fn append_group_cids(&mut self, dataset_ids_zero_based: &[usize]) -> Result<()> {
        for &id in dataset_ids_zero_based {
            self.append_cid_for_dataset(id, self.cid)?;
        }
        self.cid += 1;
        Ok(())
    }

    fn write_group(
        &mut self,
        dataset_ids_zero_based: &[usize],
        seqs: &mut Vec<Vec<u8>>,
        run_paths: &mut Vec<PathBuf>,
        spill_dir: &Path,
    ) -> Result<()> {
        if seqs.is_empty() && run_paths.is_empty() {
            return Ok(());
        }

        let mut group_sizes_buffer: Vec<u8> = Vec::new();

        let mut prev_size: usize = 0;
        emit_sorted_group_sequences(seqs, run_paths, spill_dir, |seq| {
            let encoded_seq = <Converter as Convert<&[u8]>>::str2num(seq);
            self.prev_tigs_size += encoded_seq.len() as u64;
            self.encoded_seq_buffer.extend_from_slice(&encoded_seq);
            if self.encoded_seq_buffer.len() >= ENCODED_SEQ_BUFFER_TARGET {
                self.omni_file.write_all(&self.encoded_seq_buffer)?;
                self.encoded_seq_buffer.clear();
            }

            let size = seq.len();
            let delta = size.checked_sub(prev_size).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "group sizes are not nondecreasing",
                )
            })?;
            write_varint_u64(delta as u64, &mut group_sizes_buffer);
            prev_size = size;
            Ok(())
        })?;
        self.append_bucket_group_sizes_payload(&group_sizes_buffer)?;
        self.prev_bucket_pos += 1;
        self.pos_nb_unitig
            .push((self.prev_tigs_size, self.prev_bucket_pos));

        self.append_group_cids(dataset_ids_zero_based)?;

        Ok(())
    }

    fn write_precomputed_group_files(
        &mut self,
        dataset_ids_zero_based: &[usize],
        encoded_tigs_path: &Path,
        encoded_tigs_len: u64,
        group_sizes_path: &Path,
        group_sizes_len: u64,
    ) -> Result<()> {
        let mut encoded_reader =
            BufReader::with_capacity(IO_BUFFER_CAPACITY, File::open(encoded_tigs_path)?);
        io::copy(&mut encoded_reader, &mut self.omni_file)?;
        self.prev_tigs_size += encoded_tigs_len;

        self.prev_bucket_pos += 1;
        self.pos_nb_unitig
            .push((self.prev_tigs_size, self.prev_bucket_pos));
        let mut group_sizes_payload = Vec::with_capacity(group_sizes_len as usize);
        let mut sizes_reader =
            BufReader::with_capacity(IO_BUFFER_CAPACITY, File::open(group_sizes_path)?);
        sizes_reader.read_to_end(&mut group_sizes_payload)?;
        self.append_bucket_group_sizes_payload(&group_sizes_payload)?;

        self.append_group_cids(dataset_ids_zero_based)?;
        Ok(())
    }

    fn write_precomputed_group_buffers(
        &mut self,
        dataset_ids_zero_based: &[usize],
        encoded_tigs: &[u8],
        group_sizes: &[u8],
    ) -> Result<()> {
        self.omni_file.write_all(encoded_tigs)?;
        self.prev_tigs_size += encoded_tigs.len() as u64;

        self.prev_bucket_pos += 1;
        self.pos_nb_unitig
            .push((self.prev_tigs_size, self.prev_bucket_pos));
        self.append_bucket_group_sizes_payload(group_sizes)?;

        self.append_group_cids(dataset_ids_zero_based)?;
        Ok(())
    }

    fn finalize(mut self) -> Result<Vec<(u64, u64)>> {
        if !self.encoded_seq_buffer.is_empty() {
            self.omni_file.write_all(&self.encoded_seq_buffer)?;
            self.encoded_seq_buffer.clear();
        }
        for (idx, buf) in self.spill_cid_buffers.iter_mut().enumerate() {
            if !buf.is_empty() {
                self.spill_writers[idx].write_all(buf)?;
                buf.clear();
            }
        }
        for writer in &mut self.spill_writers {
            writer.flush()?;
        }
        self.flush_size_block()?;
        self.omni_file.flush()?;
        self.size_file.flush()?;

        println!(
            "Completed compression: total tigs={}, total sizes={}",
            self.prev_tigs_size, self.prev_bucket_pos
        );
        Ok(self.pos_nb_unitig)
    }
}

#[derive(Debug)]
struct GroupTask {
    cid: usize,
    dataset_ids_zero_based: Vec<usize>,
    seqs: Vec<Vec<u8>>,
    run_paths: Vec<PathBuf>,
}

#[derive(Debug)]
struct GroupWorkItem {
    tasks: Vec<GroupTask>,
    spill_dir: PathBuf,
}

#[derive(Debug)]
enum GroupWriteData {
    InMemory {
        encoded_tigs: Vec<u8>,
        group_sizes: Vec<u8>,
    },
    Spilled {
        encoded_tigs_path: PathBuf,
        encoded_tigs_len: u64,
        group_sizes_path: PathBuf,
        group_sizes_len: u64,
    },
}

#[derive(Debug)]
struct GroupWriteResult {
    cid: usize,
    dataset_ids_zero_based: Vec<usize>,
    data: GroupWriteData,
    worker_wall_sec: f64,
    worker_cpu_sec: Option<f64>,
}

#[derive(Debug)]
struct GroupWorkResultBatch {
    results: Vec<GroupWriteResult>,
}

fn map_color_ids_to_zero_based(color_ids: &[u32], max_datasets: usize) -> io::Result<Vec<usize>> {
    let ids = color_ids
        .iter()
        .map(|&id| {
            if id == 0 {
                Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "dataset id 0 is invalid",
                ))
            } else {
                Ok(id as usize - 1)
            }
        })
        .collect::<io::Result<Vec<usize>>>()?;

    for &id in &ids {
        if id >= max_datasets {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("dataset id {} outside [0, {})", id, max_datasets),
            ));
        }
    }
    Ok(ids)
}

fn process_single_group_task(
    mut task: GroupTask,
    spill_dir: &Path,
) -> io::Result<GroupWriteResult> {
    let group_timer = PhaseTimer::start();
    let seq_bytes = task.seqs.iter().map(Vec::len).sum::<usize>();
    let use_in_memory = task.run_paths.is_empty() && seq_bytes <= GROUP_IN_MEMORY_MAX_SEQ_BYTES;

    let data = if use_in_memory {
        let mut encoded_tigs = Vec::with_capacity(seq_bytes / 4 + 1024);
        let mut group_sizes = Vec::new();
        let mut prev_size: usize = 0;
        emit_sorted_group_sequences(&mut task.seqs, &mut task.run_paths, spill_dir, |seq| {
            let encoded_seq = <Converter as Convert<&[u8]>>::str2num(seq);
            encoded_tigs.extend_from_slice(&encoded_seq);
            let size = seq.len();
            let delta = size.checked_sub(prev_size).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "group sizes are not nondecreasing",
                )
            })?;
            write_varint_u64(delta as u64, &mut group_sizes);
            prev_size = size;
            Ok(())
        })?;
        GroupWriteData::InMemory {
            encoded_tigs,
            group_sizes,
        }
    } else {
        let encoded_tigs_path = spill_dir.join(format!("group_{:012}.tigs.tmp", task.cid));
        let group_sizes_path = spill_dir.join(format!("group_{:012}.sizes.tmp", task.cid));

        let encoded_file = File::create(&encoded_tigs_path)?;
        let mut encoded_writer = BufWriter::with_capacity(IO_BUFFER_CAPACITY, encoded_file);
        let mut sizes_writer =
            BufWriter::with_capacity(IO_BUFFER_CAPACITY, File::create(&group_sizes_path)?);
        let mut prev_size: usize = 0;

        let emit_res =
            emit_sorted_group_sequences(&mut task.seqs, &mut task.run_paths, spill_dir, |seq| {
                let encoded_seq = <Converter as Convert<&[u8]>>::str2num(seq);
                encoded_writer.write_all(&encoded_seq)?;

                let size = seq.len();
                let delta = size.checked_sub(prev_size).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "group sizes are not nondecreasing",
                    )
                })?;
                write_varint_u64_to_writer(delta as u64, &mut sizes_writer)?;
                prev_size = size;
                Ok(())
            });

        if emit_res.is_err() {
            for run_path in task.run_paths.drain(..) {
                let _ = fs::remove_file(run_path);
            }
            let _ = fs::remove_file(&encoded_tigs_path);
            let _ = fs::remove_file(&group_sizes_path);
        }
        emit_res?;
        sizes_writer.flush()?;
        encoded_writer.flush()?;
        drop(sizes_writer);
        drop(encoded_writer);

        let encoded_tigs_len = fs::metadata(&encoded_tigs_path)?.len();
        let group_sizes_len = fs::metadata(&group_sizes_path)?.len();
        GroupWriteData::Spilled {
            encoded_tigs_path,
            encoded_tigs_len,
            group_sizes_path,
            group_sizes_len,
        }
    };

    let group_timing = group_timer.finish();

    Ok(GroupWriteResult {
        cid: task.cid,
        dataset_ids_zero_based: task.dataset_ids_zero_based,
        data,
        worker_wall_sec: group_timing.wall_sec,
        worker_cpu_sec: group_timing.cpu_sec,
    })
}

fn process_group_work_item(work: GroupWorkItem) -> io::Result<GroupWorkResultBatch> {
    let mut results = Vec::with_capacity(work.tasks.len());
    for task in work.tasks {
        results.push(process_single_group_task(task, &work.spill_dir)?);
    }
    Ok(GroupWorkResultBatch { results })
}

fn absorb_group_result(
    received: io::Result<GroupWorkResultBatch>,
    pending: &mut HashMap<usize, GroupWriteResult>,
    writer: &mut StreamWriterState,
    worker_stats: &mut WorkerStats,
    commit_stats: &mut CommitStats,
) -> io::Result<()> {
    let batch = received?;
    for result in batch.results {
        worker_stats.record_group(result.worker_wall_sec, result.worker_cpu_sec);
        if pending.insert(result.cid, result).is_some() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "duplicate group result for cid",
            ));
        }
    }

    while let Some(next) = pending.remove(&writer.cid) {
        let commit_timer = PhaseTimer::start();
        match &next.data {
            GroupWriteData::InMemory {
                encoded_tigs,
                group_sizes,
            } => {
                writer.write_precomputed_group_buffers(
                    &next.dataset_ids_zero_based,
                    encoded_tigs,
                    group_sizes,
                )?;
                commit_stats.tigs_bytes += encoded_tigs.len() as u64;
                commit_stats.sizes_bytes += group_sizes.len() as u64;
            }
            GroupWriteData::Spilled {
                encoded_tigs_path,
                encoded_tigs_len,
                group_sizes_path,
                group_sizes_len,
            } => {
                writer.write_precomputed_group_files(
                    &next.dataset_ids_zero_based,
                    encoded_tigs_path,
                    *encoded_tigs_len,
                    group_sizes_path,
                    *group_sizes_len,
                )?;
                commit_stats.tigs_bytes += *encoded_tigs_len;
                commit_stats.sizes_bytes += *group_sizes_len;
                let _ = fs::remove_file(encoded_tigs_path);
                let _ = fs::remove_file(group_sizes_path);
            }
        }
        commit_stats.groups += 1;
        commit_stats.timing.add(commit_timer.finish());
    }
    Ok(())
}

fn drain_available_group_results(
    result_rx: &mpsc::Receiver<io::Result<GroupWorkResultBatch>>,
    pending: &mut HashMap<usize, GroupWriteResult>,
    writer: &mut StreamWriterState,
    worker_stats: &mut WorkerStats,
    commit_stats: &mut CommitStats,
) -> io::Result<()> {
    loop {
        match result_rx.try_recv() {
            Ok(received) => {
                absorb_group_result(received, pending, writer, worker_stats, commit_stats)?
            }
            Err(mpsc::TryRecvError::Empty) => break,
            Err(mpsc::TryRecvError::Disconnected) => break,
        }
    }
    Ok(())
}

fn write_compressed_from_stream(
    unitigs_file_path: String,
    output_dir: &String,
    nb_files: u32,
    worker_threads: usize,
    record_rx: mpsc::Receiver<SimplitigBatch>,
) -> Result<(Vec<(u64, u64)>, Vec<PathBuf>, PathBuf)> {
    let total_timer = PhaseTimer::start();
    let (mut writer, spill_paths, spill_dir) =
        StreamWriterState::new(unitigs_file_path, output_dir, nb_files as usize)?;

    let group_spill_dir = Path::new(output_dir).join(format!(
        ".kloe-group-sort-spill-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0)
    ));
    fs::create_dir_all(&group_spill_dir)?;

    let setup_workers_timer = PhaseTimer::start();
    let worker_count = worker_threads.max(1);
    let (result_tx, result_rx) = mpsc::channel::<io::Result<GroupWorkResultBatch>>();
    let mut worker_inputs = Vec::with_capacity(worker_count);
    let mut worker_handles = Vec::with_capacity(worker_count);

    for _ in 0..worker_count {
        let (work_tx, work_rx) =
            mpsc::sync_channel::<Option<GroupWorkItem>>(GROUP_WORKER_QUEUE_DEPTH);
        worker_inputs.push(work_tx);
        let worker_result_tx = result_tx.clone();
        worker_handles.push(thread::spawn(move || {
            while let Ok(msg) = work_rx.recv() {
                match msg {
                    Some(work) => {
                        let send_res = worker_result_tx.send(process_group_work_item(work));
                        if send_res.is_err() {
                            break;
                        }
                    }
                    None => break,
                }
            }
        }));
    }
    drop(result_tx);
    let setup_workers_timing = setup_workers_timer.finish();

    let mut pending_results: HashMap<usize, GroupWriteResult> = HashMap::new();
    let mut worker_stats = WorkerStats::default();
    let mut commit_stats = CommitStats::default();
    let mut current_color_ids: Option<Arc<Vec<u32>>> = None;
    let mut current_ids_zero_based: Option<Vec<usize>> = None;
    let mut group_seqs: Vec<Vec<u8>> = Vec::new();
    let mut group_run_paths: Vec<PathBuf> = Vec::new();
    let mut pending_group_tasks: Vec<GroupTask> = Vec::with_capacity(GROUP_WORK_BATCH_GROUPS);
    let mut pending_group_task_seq_bytes: usize = 0;
    let mut group_bytes: usize = 0;
    let mut group_total_seq_bytes: usize = 0;
    let mut submitted_groups: usize = 0;
    let mut rr_index = 0usize;
    let dispatch_timer = PhaseTimer::start();

    for batch in record_rx {
        for record in batch {
            let key_changed = match current_color_ids.as_ref() {
                None => true,
                Some(ids) => {
                    !(Arc::ptr_eq(ids, &record.color_ids)
                        || ids.as_ref() == record.color_ids.as_ref())
                }
            };

            if key_changed {
                if let Some(ids) = current_ids_zero_based.take() {
                    if !group_seqs.is_empty() || !group_run_paths.is_empty() {
                        pending_group_tasks.push(GroupTask {
                            cid: submitted_groups,
                            dataset_ids_zero_based: ids,
                            seqs: std::mem::take(&mut group_seqs),
                            run_paths: std::mem::take(&mut group_run_paths),
                        });
                        pending_group_task_seq_bytes =
                            pending_group_task_seq_bytes.saturating_add(group_total_seq_bytes);
                        submitted_groups += 1;
                        if pending_group_tasks.len() >= GROUP_WORK_BATCH_GROUPS
                            || pending_group_task_seq_bytes >= GROUP_WORK_BATCH_SEQ_BYTES
                        {
                            let work = GroupWorkItem {
                                tasks: std::mem::take(&mut pending_group_tasks),
                                spill_dir: group_spill_dir.clone(),
                            };
                            pending_group_task_seq_bytes = 0;
                            worker_inputs[rr_index % worker_inputs.len()]
                                .send(Some(work))
                                .map_err(|err| {
                                    io::Error::other(format!("group worker channel closed: {err}"))
                                })?;
                            rr_index += 1;
                            drain_available_group_results(
                                &result_rx,
                                &mut pending_results,
                                &mut writer,
                                &mut worker_stats,
                                &mut commit_stats,
                            )?;
                        }
                    }
                }

                current_ids_zero_based = Some(map_color_ids_to_zero_based(
                    &record.color_ids,
                    writer.spill_writers.len(),
                )?);
                current_color_ids = Some(Arc::clone(&record.color_ids));
                group_bytes = 0;
                group_total_seq_bytes = 0;
            }

            group_bytes += record.seq.len();
            group_total_seq_bytes += record.seq.len();
            group_seqs.push(record.seq);
            if group_bytes >= GROUP_SORT_SPILL_BYTES {
                sort_and_spill_sequence_run(
                    &mut group_seqs,
                    &mut group_run_paths,
                    &group_spill_dir,
                )?;
                group_bytes = 0;
            }
        }
    }

    if let Some(ids) = current_ids_zero_based.take() {
        if !group_seqs.is_empty() || !group_run_paths.is_empty() {
            pending_group_tasks.push(GroupTask {
                cid: submitted_groups,
                dataset_ids_zero_based: ids,
                seqs: std::mem::take(&mut group_seqs),
                run_paths: std::mem::take(&mut group_run_paths),
            });
            submitted_groups += 1;
        }
    }

    if !pending_group_tasks.is_empty() {
        let work = GroupWorkItem {
            tasks: std::mem::take(&mut pending_group_tasks),
            spill_dir: group_spill_dir.clone(),
        };
        worker_inputs[rr_index % worker_inputs.len()]
            .send(Some(work))
            .map_err(|err| io::Error::other(format!("group worker channel closed: {err}")))?;
    }
    let dispatch_timing = dispatch_timer.finish();

    for input in worker_inputs {
        let _ = input.send(None);
    }

    let wait_results_timer = PhaseTimer::start();
    while writer.cid < submitted_groups {
        let received = result_rx.recv().map_err(|err| {
            io::Error::new(
                io::ErrorKind::BrokenPipe,
                format!("group worker results closed early: {err}"),
            )
        })?;
        absorb_group_result(
            received,
            &mut pending_results,
            &mut writer,
            &mut worker_stats,
            &mut commit_stats,
        )?;
    }
    let wait_results_timing = wait_results_timer.finish();

    let join_workers_timer = PhaseTimer::start();
    for handle in worker_handles {
        if handle.join().is_err() {
            return Err(io::Error::other("group worker panicked"));
        }
    }
    let join_workers_timing = join_workers_timer.finish();

    let pos_nb_unitig = writer.finalize()?;
    let _ = fs::remove_dir_all(&group_spill_dir);

    log_phase_timing("post_ggcat_stream.setup_workers", setup_workers_timing);
    log_phase_timing("post_ggcat_stream.dispatch_records", dispatch_timing);
    log_phase_timing(
        "post_ggcat_stream.worker_groups_total",
        PhaseTiming {
            wall_sec: worker_stats.timing.wall_sec,
            cpu_sec: worker_stats.timing.cpu_sec,
        },
    );
    log_phase_timing(
        "post_ggcat_stream.writer_commit_groups",
        PhaseTiming {
            wall_sec: commit_stats.timing.wall_sec,
            cpu_sec: commit_stats.timing.cpu_sec,
        },
    );
    log_phase_timing("post_ggcat_stream.wait_results", wait_results_timing);
    log_phase_timing("post_ggcat_stream.join_workers", join_workers_timing);
    log_phase_timing("post_ggcat_stream.total", total_timer.finish());
    println!(
        "[phase-stats] phase=post_ggcat_stream workers={} groups_submitted={} groups_processed={} groups_committed={} commit_tigs_bytes={} commit_sizes_bytes={} worker_max_group_wall_s={:.3} worker_max_group_cpu_s={}",
        worker_count,
        submitted_groups,
        worker_stats.groups,
        commit_stats.groups,
        commit_stats.tigs_bytes,
        commit_stats.sizes_bytes,
        worker_stats.max_wall_sec,
        worker_stats
            .max_cpu_sec
            .map(|v| format!("{v:.3}"))
            .unwrap_or_else(|| "NA".to_string())
    );

    Ok((pos_nb_unitig, spill_paths, spill_dir))
}

fn write_positions(pos_nb_unitigs: Vec<(u64, u64)>, filepath: String) -> Result<()> {
    let mut pos_file = BufWriter::with_capacity(IO_BUFFER_CAPACITY, File::create(filepath)?);
    pos_file.write_all(POSITIONS_MAGIC)?;
    write_varint_u64_to_writer(pos_nb_unitigs.len() as u64, &mut pos_file)?;

    let mut prev_tigs = 0u64;
    let mut prev_sizes = 0u64;
    for (tigs_pos, sizes_pos) in &pos_nb_unitigs {
        if *tigs_pos < prev_tigs || *sizes_pos < prev_sizes {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "positions are not monotonic",
            ));
        }
        write_varint_u64_to_writer(*tigs_pos - prev_tigs, &mut pos_file)?;
        write_varint_u64_to_writer(*sizes_pos - prev_sizes, &mut pos_file)?;
        prev_tigs = *tigs_pos;
        prev_sizes = *sizes_pos;
    }
    pos_file.flush()?;
    Ok(())
}

fn create_id_cid_spill_dir(output_dir: &str) -> Result<PathBuf> {
    let pid = std::process::id();
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let path = Path::new(output_dir).join(format!(".kloe-id-cid-spill-{pid}-{nanos}"));
    fs::create_dir_all(&path)?;
    Ok(path)
}

fn build_id_cid_payload_from_spill(path: &Path, cid_file_path: &str) -> std::io::Result<Vec<u8>> {
    let mut reader = BufReader::with_capacity(IO_BUFFER_CAPACITY, File::open(path)?);
    let mut raw_payload = Vec::new();
    let mut prev_cid = 0usize;
    loop {
        let mut raw = [0u8; 8];
        match std::io::Read::read_exact(&mut reader, &mut raw) {
            Ok(()) => {
                let raw_cid = u64::from_le_bytes(raw);
                let cid = usize::try_from(raw_cid).map_err(|_| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "color id {} cannot be represented as usize while writing {}",
                            raw_cid, cid_file_path
                        ),
                    )
                })?;
                if cid < prev_cid {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "color ids are not sorted ({} before {}) while writing {}",
                            prev_cid, cid, cid_file_path
                        ),
                    ));
                }
                let delta = cid - prev_cid;
                write_varint_u64(delta as u64, &mut raw_payload);
                prev_cid = cid;
            }
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => break,
            Err(e) => return Err(e),
        }
    }

    let mut payload = Vec::new();
    {
        let mut cid_encoder = Encoder::new(&mut payload, 1)?;
        cid_encoder.write_all(&raw_payload)?;
        cid_encoder.finish()?;
    }
    Ok(payload)
}

fn write_id_to_color_id_from_spills(
    cid_file_path: String,
    spill_paths: &[PathBuf],
    spill_dir: &Path,
    worker_threads: usize,
) -> std::io::Result<Vec<usize>> {
    let total_timer = PhaseTimer::start();
    let mut cid_file = BufWriter::with_capacity(IO_BUFFER_CAPACITY, File::create(&cid_file_path)?);
    cid_file.write_all(ID_TO_CID_MAGIC)?;
    let mut id_cid_line_sizes = Vec::with_capacity(spill_paths.len());
    let mut tot_size = ID_TO_CID_MAGIC.len();
    let mut per_dataset_timing = TimingAccumulator::default();
    let mut total_payload_bytes = 0usize;
    let mut payloads: Vec<Option<Vec<u8>>> = vec![None; spill_paths.len()];

    if spill_paths.len() <= 1 || worker_threads <= 1 {
        for (idx, path) in spill_paths.iter().enumerate() {
            let dataset_timer = PhaseTimer::start();
            let payload = build_id_cid_payload_from_spill(path, &cid_file_path)?;
            let _ = fs::remove_file(path);
            per_dataset_timing.add(dataset_timer.finish());
            payloads[idx] = Some(payload);
        }
    } else {
        let max_threads = worker_threads.max(1).min(spill_paths.len());
        let build_payloads = || -> io::Result<Vec<(usize, Vec<u8>, PhaseTiming)>> {
            spill_paths
                .par_iter()
                .enumerate()
                .map(|(idx, path)| -> io::Result<(usize, Vec<u8>, PhaseTiming)> {
                    let dataset_timer = PhaseTimer::start();
                    let payload = build_id_cid_payload_from_spill(path, &cid_file_path)?;
                    let _ = fs::remove_file(path);
                    Ok((idx, payload, dataset_timer.finish()))
                })
                .collect()
        };

        let results = match rayon::ThreadPoolBuilder::new()
            .num_threads(max_threads)
            .build()
        {
            Ok(pool) => pool.install(build_payloads),
            Err(_) => build_payloads(),
        }?;

        for (idx, payload, timing) in results {
            per_dataset_timing.add(timing);
            payloads[idx] = Some(payload);
        }
    }

    for payload in payloads {
        let payload = payload
            .ok_or_else(|| io::Error::other("missing id->cid payload after spill processing"))?;
        id_cid_line_sizes.push(tot_size);
        tot_size += 8 + payload.len();
        total_payload_bytes += payload.len();
        cid_file.write_all(&(payload.len() as u64).to_le_bytes())?;
        cid_file.write_all(&payload)?;
    }

    cid_file.write_all(&(0_u64).to_le_bytes())?;
    cid_file.flush()?;
    let _ = fs::remove_dir(spill_dir);
    log_phase_timing(
        "post_ggcat.id_to_cid.per_dataset_total",
        PhaseTiming {
            wall_sec: per_dataset_timing.wall_sec,
            cpu_sec: per_dataset_timing.cpu_sec,
        },
    );
    log_phase_timing("post_ggcat.id_to_cid.total", total_timer.finish());
    println!(
        "[phase-stats] phase=post_ggcat.id_to_cid datasets={} payload_bytes={}",
        spill_paths.len(),
        total_payload_bytes
    );
    Ok(id_cid_line_sizes)
}

pub(crate) fn sort_by_bucket_streaming(
    output_dir: &String,
    nb_files: u32,
    worker_threads: usize,
    record_rx: mpsc::Receiver<SimplitigBatch>,
) -> Vec<usize> {
    let total_timer = PhaseTimer::start();
    println!("Starting writing compressed sequences (streaming).");
    let write_stream_timer = PhaseTimer::start();
    let triple = match write_compressed_from_stream(
        output_dir.clone() + "tigs_kloe.fa",
        output_dir,
        nb_files,
        worker_threads,
        record_rx,
    ) {
        Ok(res_pair) => res_pair,
        Err(e) => panic!("Error writing compressed unitigs: {e:?}"),
    };
    let write_stream_timing = write_stream_timer.finish();
    println!(
        "Writing compressed sequences wall time: {:.3}s",
        write_stream_timing.wall_sec
    );
    log_phase_timing("post_ggcat.write_stream", write_stream_timing);
    let pos_nb_unitig = triple.0;
    let spill_paths = triple.1;
    let spill_dir = triple.2;

    let position_timer = PhaseTimer::start();
    println!("Starting to write positions");
    if let Err(e) = write_positions(
        pos_nb_unitig,
        String::from(output_dir.clone() + "positions_kloe.bin"),
    ) {
        panic!("Error writting positions: {e:?}");
    }
    let position_timing = position_timer.finish();
    println!(
        "Write positions wall time: {:.3}s",
        position_timing.wall_sec
    );
    log_phase_timing("post_ggcat.write_positions", position_timing);
    let id_timer = PhaseTimer::start();
    let write_id_cid = match write_id_to_color_id_from_spills(
        output_dir.clone() + "id_to_color_id.txt.zst",
        &spill_paths,
        &spill_dir,
        worker_threads,
    ) {
        Ok(id_cid_line_sizes) => id_cid_line_sizes,
        Err(e) => panic!("error writting id to color id list: {e:?}"),
    };
    let id_timing = id_timer.finish();
    println!("Write id to cid wall time: {:.3}s", id_timing.wall_sec);
    log_phase_timing("post_ggcat.write_id_to_cid", id_timing);
    let total_timing = total_timer.finish();
    println!("Compression took: {:.3}s", total_timing.wall_sec);
    log_phase_timing("post_ggcat.total", total_timing);
    write_id_cid
}

fn stream_sorted_records_from_chunks(
    chunk_paths: &[PathBuf],
    mut emit: impl FnMut(SortedColorRecord) -> Result<()>,
) -> Result<()> {
    let mut readers = Vec::with_capacity(chunk_paths.len());
    for path in chunk_paths {
        readers.push(BufReader::with_capacity(
            IO_BUFFER_CAPACITY,
            File::open(path)?,
        ));
    }

    let mut heap = BinaryHeap::new();
    for (idx, reader) in readers.iter_mut().enumerate() {
        if let Some(record) = read_sorted_color_record(reader)? {
            heap.push(ChunkHeapItem {
                chunk_index: idx,
                record,
            });
        }
    }

    while let Some(item) = heap.pop() {
        emit(item.record)?;
        if let Some(next) = read_sorted_color_record(&mut readers[item.chunk_index])? {
            heap.push(ChunkHeapItem {
                chunk_index: item.chunk_index,
                record: next,
            });
        }
    }

    Ok(())
}

struct FastaBlockReader {
    receiver: mpsc::Receiver<Option<Vec<u8>>>,
    current: io::Cursor<Vec<u8>>,
}

impl FastaBlockReader {
    fn new(receiver: mpsc::Receiver<Option<Vec<u8>>>) -> Self {
        Self {
            receiver,
            current: io::Cursor::new(Vec::new()),
        }
    }
}

impl Read for FastaBlockReader {
    fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
        loop {
            let read = self.current.read(output)?;
            if read != 0 {
                return Ok(read);
            }
            match self.receiver.recv() {
                Ok(Some(block)) => self.current = io::Cursor::new(block),
                Ok(None) | Err(_) => return Ok(0),
            }
        }
    }
}

fn read_structured_varint(data: &[u8], position: &mut usize) -> Result<u64> {
    let mut value = 0u64;
    let mut shift = 0u32;
    loop {
        let byte = *data.get(*position).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "truncated GGCAT structured color payload",
            )
        })?;
        *position += 1;
        if shift >= 64 && byte & 0x7f != 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "overflow in GGCAT structured color payload",
            ));
        }
        value |= u64::from(byte & 0x7f) << shift;
        if byte & 0x80 == 0 {
            return Ok(value);
        }
        shift += 7;
    }
}

fn parse_structured_color_segments(
    record: PlainFastaOutputRecord,
    k: usize,
) -> Result<Vec<SortedColorRecord>> {
    let PlainFastaOutputRecord {
        sequence,
        color_data,
    } = record;
    if sequence.len() < k {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "GGCAT structured sequence shorter than k (len={}, k={})",
                sequence.len(),
                k
            ),
        ));
    }

    let mut position = 0usize;
    let runs_count = usize::try_from(read_structured_varint(&color_data, &mut position)?)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "too many GGCAT color runs"))?;
    if runs_count == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "GGCAT emitted structured sequence without color runs",
        ));
    }

    let mut total_run_kmers = 0usize;
    let mut merged_runs: Vec<(ColorIndexType, usize)> = Vec::with_capacity(runs_count);
    for _ in 0..runs_count {
        let subset = ColorIndexType::try_from(read_structured_varint(
            &color_data,
            &mut position,
        )?)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "GGCAT subset id overflow"))?;
        let count = usize::try_from(read_structured_varint(&color_data, &mut position)?)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "GGCAT run length overflow"))?;
        total_run_kmers = total_run_kmers.checked_add(count).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "overflow while summing GGCAT structured color-run lengths",
            )
        })?;
        if let Some((last_subset, last_count)) = merged_runs.last_mut() {
            if *last_subset == subset {
                *last_count += count;
                continue;
            }
        }
        if count > 0 {
            merged_runs.push((subset, count));
        }
    }
    if position != color_data.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "trailing bytes in GGCAT structured color payload",
        ));
    }

    let seq_kmers = sequence.len() - k + 1;
    if total_run_kmers != seq_kmers {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "GGCAT structured color-run lengths mismatch: runs={}, expected={}",
                total_run_kmers, seq_kmers
            ),
        ));
    }

    let mut segments = Vec::with_capacity(merged_runs.len());
    let mut run_start = 0usize;
    for (subset, count) in merged_runs {
        let run_end = run_start + count;
        segments.push(SortedColorRecord {
            subset,
            seq: sequence[run_start..(run_end + k - 1)].to_vec(),
        });
        run_start = run_end;
    }
    Ok(segments)
}

fn push_grouped_record(
    groups_by_subset: &mut Vec<Vec<Vec<u8>>>,
    unique_subsets: &mut HashSet<ColorIndexType>,
    record: SortedColorRecord,
) -> Result<()> {
    let subset_index = usize::try_from(record.subset).map_err(|_| {
        io::Error::new(io::ErrorKind::InvalidData, "GGCAT subset id overflow")
    })?;
    if subset_index >= groups_by_subset.len() {
        groups_by_subset.resize_with(subset_index + 1, Vec::new);
    }
    unique_subsets.insert(record.subset);
    groups_by_subset[subset_index].push(record.seq);
    Ok(())
}

struct ParsedColorRecords {
    groups_by_subset: Vec<Vec<Vec<u8>>>,
    input_sequences_count: usize,
    emitted_segments_count: usize,
    unique_subsets: HashSet<ColorIndexType>,
    parse_timing: PhaseTiming,
}

fn parse_streamed_color_records(
    receiver: mpsc::Receiver<Option<Vec<PlainFastaOutputRecord>>>,
    k: usize,
) -> Result<ParsedColorRecords> {
    let parse_timer = PhaseTimer::start();
    let mut groups_by_subset = Vec::<Vec<Vec<u8>>>::new();
    let mut input_sequences_count = 0usize;
    let mut emitted_segments_count = 0usize;
    let mut unique_subsets = HashSet::<ColorIndexType>::new();

    while let Ok(message) = receiver.recv() {
        let Some(batch) = message else {
            break;
        };
        input_sequences_count += batch.len();
        let parsed_batches = batch
            .into_par_iter()
            .map(|record| parse_structured_color_segments(record, k))
            .collect::<Result<Vec<_>>>()?;

        for segments in parsed_batches {
            emitted_segments_count += segments.len();
            for record in segments {
                push_grouped_record(&mut groups_by_subset, &mut unique_subsets, record)?;
            }
        }
    }

    Ok(ParsedColorRecords {
        groups_by_subset,
        input_sequences_count,
        emitted_segments_count,
        unique_subsets,
        parse_timing: parse_timer.finish(),
    })
}

fn produce_ggcat_records(
    filenames: Vec<String>,
    threads: usize,
    k: usize,
    m: usize,
    use_unitigs: bool,
    use_matchtigs: bool,
    use_eulertigs: bool,
    ggcat_cfg: GgcatCompressionConfig,
    sender: mpsc::SyncSender<SimplitigBatch>,
) -> Result<()> {
    let dataset_count = filenames.len();
    if dataset_count == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "input file list is empty",
        ));
    }

    let ggcat_temp_root = if ggcat_cfg.temp_dir.is_empty() {
        std::env::temp_dir().join(format!("kloe-ggcat-shared-{}", std::process::id()))
    } else {
        PathBuf::from(&ggcat_cfg.temp_dir)
    };
    fs::create_dir_all(&ggcat_temp_root)?;

    let pid = std::process::id();
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);

    let raw_records_file = ggcat_temp_root.join(format!("kloe-color-records-{pid}-{nanos}.tmp"));
    let chunk_dir = ggcat_temp_root.join(format!("kloe-color-chunks-{pid}-{nanos}"));
    fs::create_dir_all(&chunk_dir)?;

    let instance = GGCATInstance::create(GGCATConfig {
        temp_dir: Some(ggcat_temp_root.clone()),
        memory: ggcat_cfg.memory_gb.max(1) as f64,
        prefer_memory: true,
        total_threads_count: threads.max(1),
        intermediate_compression_level: None,
        stats_file: None,
        messages_callback: None,
    })
    .map_err(|e| to_io_err("create ggcat instance", e))?;

    let input_streams = filenames
        .iter()
        .map(|file| {
            let path = PathBuf::from(file);
            let resolved = fs::canonicalize(&path).unwrap_or(path);
            GeneralSequenceBlockData::FASTA((resolved, None))
        })
        .collect::<Vec<_>>();

    let color_names = (1..=dataset_count)
        .map(|idx| idx.to_string())
        .collect::<Vec<_>>();

    let total_threads = threads.max(1);
    let parser_threads = if total_threads >= 4 {
        (total_threads / 4).max(1)
    } else {
        1
    };
    let ggcat_threads = total_threads.saturating_sub(parser_threads).max(1);
    println!(
        "Running embedded ggcat streaming path (k={}, m={}, ggcat_threads={}, parser_threads={}, memory={}GB)",
        k,
        m,
        ggcat_threads,
        parser_threads,
        ggcat_cfg.memory_gb.max(1)
    );

    let (fasta_block_tx, fasta_block_rx) =
        mpsc::sync_channel::<Option<Vec<PlainFastaOutputRecord>>>(16);
    let fasta_block_end_tx = fasta_block_tx.clone();
    let parser_handle = thread::spawn(move || {
        let parser_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(parser_threads)
            .build()
            .map_err(|err| io::Error::other(format!("build GGCAT parser pool: {err}")))?;
        parser_pool.install(|| parse_streamed_color_records(fasta_block_rx, k))
    });
    let callback_guard = install_plain_fasta_output_callback(
        &raw_records_file,
        Arc::new(move |block| {
            let _ = fasta_block_tx.send(Some(block));
        }),
    )?;

    let ggcat_build_timer = PhaseTimer::start();
    let build_result = instance.build_graph(
        input_streams,
        raw_records_file,
        Some(&color_names),
        k,
        ggcat_threads,
        false,
        Some(m),
        true,
        1,
        ggcat_extra_elaboration(use_unitigs, use_matchtigs, use_eulertigs),
        None,
    );
    log_phase_timing("ggcat.build_graph_api_call", ggcat_build_timer.finish());
    let _ = fasta_block_end_tx.send(None);
    drop(callback_guard);

    let parsed_records = parser_handle
        .join()
        .map_err(|_| io::Error::other("streamed GGCAT parser panicked"))??;
    let records_output = build_result.map_err(|e| to_io_err("ggcat build_graph", e))?;

    let colormap_file = GGCATInstance::get_colormap_file(&records_output);
    let colors_deserializer =
        ColorsDeserializer::<DefaultColorsSerializer>::new(&colormap_file, true)
            .map_err(|e| to_io_err("open ggcat colormap", e))?;

    let colors_count = colors_deserializer.colors_count();
    if colors_count != dataset_count {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "ggcat returned {colors_count} colors but KLOE loaded {dataset_count} datasets"
            ),
        ));
    }

    let dumped_colors = GGCATInstance::dump_colors(&colormap_file)
        .map_err(|e| to_io_err("dump ggcat colors", e))?
        .collect::<Vec<_>>();
    if dumped_colors.len() != colors_count {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "ggcat dumped {} color names but reports {} colors",
                dumped_colors.len(),
                colors_count
            ),
        ));
    }
    let mut color_index_to_dataset_index = Vec::with_capacity(colors_count);
    for color_name in dumped_colors {
        let parsed = color_name.parse::<usize>().map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "ggcat color name '{}' is not a dataset index; expected numeric names",
                    color_name
                ),
            )
        })?;
        if parsed == 0 || parsed > dataset_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "ggcat color name '{}' maps outside dataset range [1, {}]",
                    color_name, dataset_count
                ),
            ));
        }
        color_index_to_dataset_index.push(parsed - 1);
    }

    let keep_ggcat = std::env::var_os("KLOE_KEEP_GGCAT").is_some();
    let process_result = (|| -> Result<()> {
        let producer_total_timer = PhaseTimer::start();
        let ParsedColorRecords {
            groups_by_subset,
            input_sequences_count,
            emitted_segments_count,
            unique_subsets,
            parse_timing,
        } = parsed_records;

        let resolve_subsets_timer = PhaseTimer::start();
        let subset_to_dataset_ids = resolve_subsets_to_dataset_ids(
            &instance,
            &colormap_file,
            &unique_subsets,
            &color_index_to_dataset_index,
            dataset_count,
        )?;
        let resolve_subsets_timing = resolve_subsets_timer.finish();

        let mut batch: SimplitigBatch = Vec::with_capacity(COLOR_RECORD_BATCH_SIZE);
        let mut grouped_records_count = 0usize;
        let mut emitted_batches_count = 0usize;
        let mut emitted_records_count = 0usize;
        let grouped_emit_timer = PhaseTimer::start();

        for (subset_index, sequences) in groups_by_subset.into_iter().enumerate() {
            if sequences.is_empty() {
                continue;
            }
            let subset = ColorIndexType::try_from(subset_index).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "GGCAT subset id overflow")
            })?;
            let color_ids = subset_to_dataset_ids.get(&subset).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("missing resolved dataset ids for subset {}", subset),
                )
            })?;

            for seq in sequences {
                grouped_records_count += 1;
                batch.push(SimplitigRecord {
                    color_ids: Arc::clone(color_ids),
                    seq,
                });

                if batch.len() >= COLOR_RECORD_BATCH_SIZE {
                    let out = std::mem::take(&mut batch);
                    emitted_batches_count += 1;
                    emitted_records_count += out.len();
                    sender.send(out).map_err(|e| {
                        io::Error::new(
                            io::ErrorKind::BrokenPipe,
                            format!("record receiver dropped: {e}"),
                        )
                    })?;
                    batch = Vec::with_capacity(COLOR_RECORD_BATCH_SIZE);
                }
            };
        }

        if !batch.is_empty() {
            emitted_batches_count += 1;
            emitted_records_count += batch.len();
            sender.send(batch).map_err(|e| {
                io::Error::new(
                    io::ErrorKind::BrokenPipe,
                    format!("record receiver dropped: {e}"),
                )
            })?;
        }

        let grouped_emit_timing = grouped_emit_timer.finish();
        log_phase_timing("producer.parse_structured_direct_group", parse_timing);
        log_phase_timing(
            "producer.resolve_subsets_dataset_ids",
            resolve_subsets_timing,
        );
        log_phase_timing("producer.emit_grouped_records", grouped_emit_timing);
        log_phase_timing("producer.total_post_ggcat", producer_total_timer.finish());
        println!(
            "[phase-stats] phase=producer input_sequences={} emitted_segments={} subsets_unique={} grouped_records={} emitted_batches={} emitted_records={}",
            input_sequences_count,
            emitted_segments_count,
            unique_subsets.len(),
            grouped_records_count,
            emitted_batches_count,
            emitted_records_count
        );

        drop(sender);
        Ok(())
    })();

    if keep_ggcat {
        eprintln!(
            "DEBUG_KEEP_GGCAT graph={} colormap={} chunks_dir={}",
            records_output.display(),
            colormap_file.display(),
            chunk_dir.display()
        );
    } else {
        let _ = fs::remove_file(&records_output);
        let _ = fs::remove_file(&colormap_file);
        let _ = fs::remove_dir_all(&chunk_dir);
    }

    process_result
}

pub fn compress_with_ggcat(
    output_dir: &String,
    input_fof: &String,
    threads: usize,
    k: usize,
    m: usize,
    partition_power: u32,
    verify_kmers: bool,
    skip_sort: bool,
    use_unitigs: bool,
    use_matchtigs: bool,
    use_eulertigs: bool,
    ggcat_cfg: GgcatCompressionConfig,
) -> Result<()> {
    let _unused_partition_power = partition_power;
    if verify_kmers {
        eprintln!("Warning: --verify-kmers is currently ignored in the ggcat hard-switch backend.");
    }
    if skip_sort {
        eprintln!(
            "Warning: --skip-sort is ignored; ggcat build-colored-fasta always sorts by colors."
        );
    }

    let filenames = read_input_fof_filenames(input_fof)?;
    let dataset_count = filenames.len();
    if dataset_count == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "input file list is empty",
        ));
    }

    println!("Compression backend: embedded ggcat build-colored-fasta");

    let (record_tx, record_rx) = mpsc::sync_channel::<SimplitigBatch>(16);
    let producer_filenames = filenames.clone();
    let producer_cfg = ggcat_cfg.clone();

    let producer = thread::spawn(move || {
        produce_ggcat_records(
            producer_filenames,
            threads,
            k,
            m,
            use_unitigs,
            use_matchtigs,
            use_eulertigs,
            producer_cfg,
            record_tx,
        )
    });

    let sort_start = Instant::now();
    let id_cid_line_sizes =
        sort_by_bucket_streaming(output_dir, dataset_count as u32, threads, record_rx);
    println!(
        "Streamed ggcat colored records into KLOE archive in {:.3}s",
        sort_start.elapsed().as_secs_f64()
    );

    match producer.join() {
        Ok(Ok(())) => {}
        Ok(Err(err)) => {
            return Err(io::Error::other(format!("ggcat producer failed: {err}")));
        }
        Err(_) => {
            return Err(io::Error::other("ggcat producer panicked"));
        }
    }

    write_filenames_id_offsets(output_dir, &filenames, &id_cid_line_sizes)?;
    Ok(())
}
