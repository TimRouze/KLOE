use ggcat_api::{
    ColorIndexType, ExtraElaboration, GGCATConfig, GGCATInstance, GeneralSequenceBlockData,
    register_channel_output, unregister_channel_output,
};
use ggcat_colors::colors_manager::ColorMapReader;
use ggcat_colors::storage::deserializer::ColorsDeserializer;
use ggcat_colors::DefaultColorsSerializer;
use parking_lot::Mutex as FastMutex;
use rayon::prelude::*;
use rayon::slice::ParallelSliceMut;
use std::cmp::Ordering as CmpOrdering;
use std::collections::BinaryHeap;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Result, Seek, SeekFrom, Write};
use std::os::unix::fs::FileExt;
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
const PAR_SORT_THRESHOLD: usize = 200_000;
const COLOR_RECORD_BATCH_SIZE: usize = 65_536;
const COLOR_CHUNK_TARGET_BYTES: usize = 256 * 1024 * 1024;
const GROUP_SORT_SPILL_BYTES: usize = 128 * 1024 * 1024;
const GROUP_WORK_BATCH_GROUPS: usize = 2048;
const GROUP_IN_MEMORY_MAX_SEQ_BYTES: usize = 2 * 1024 * 1024;
const MIN_MEMORY_BUDGET_BYTES: usize = 256 * 1024 * 1024;
const MIN_GROUP_BATCH_BYTES: usize = 16 * 1024 * 1024;
const MAX_GROUP_BATCH_BYTES: usize = 512 * 1024 * 1024;
const MIN_SORT_CHUNK_BYTES: usize = 32 * 1024 * 1024;
const DATASETS_PER_CID_PARTITION: usize = 512;
const MAX_CID_PARTITIONS: usize = 256;
const CID_PARTITION_BUFFER_BYTES: usize = 256 * 1024;
const CID_TRANSPOSE_BLOCK_BYTES: usize = 64 * 1024;
const COLOR_RUN_BUFFER_BYTES: usize = 64 * 1024;
const COLOR_MERGE_FAN_IN: usize = 64;
const BUCKET_SIZES_MAGIC: &[u8; 4] = b"KSB2";
const POSITIONS_MAGIC: &[u8; 4] = b"KPS2";
const ID_TO_CID_MAGIC: &[u8; 4] = b"KIC2";
pub(crate) const CID_TO_DATASET_MAGIC: &[u8; 4] = b"KCD2";
pub(crate) const CID_TO_DATASET_FILE: &str = "cid_to_dataset_id.bin";
const BUCKET_SIZE_BLOCK_MAX_GROUPS: usize = 8_192;
const BUCKET_SIZE_BLOCK_MAX_UNCOMPRESSED_BYTES: usize = 8 * 1024 * 1024;

#[derive(Debug)]
struct CidDatasetSidecarBlock {
    first_group: usize,
    group_count: usize,
    offsets_file_offset: u64,
    data_offset: u64,
    data_len: usize,
}

#[derive(Debug)]
struct CidDatasetSidecarCache {
    file: File,
    block_index: usize,
    offsets: Vec<u32>,
    data: Vec<u8>,
    decoded_ids: Vec<u32>,
}

#[derive(Debug)]
pub(crate) struct CidDatasetSidecar {
    blocks: Vec<CidDatasetSidecarBlock>,
    group_count: usize,
    cache: Mutex<CidDatasetSidecarCache>,
}

struct CidDatasetSidecarReader<'a> {
    sidecar: &'a CidDatasetSidecar,
    file: File,
    block_index: usize,
    offsets: Vec<u32>,
    data: Vec<u8>,
}

impl CidDatasetSidecar {
    pub(crate) fn open(path: &Path) -> io::Result<Self> {
        let mut file = File::open(path)?;
        let mut magic = [0u8; 4];
        file.read_exact(&mut magic)?;
        if &magic != CID_TO_DATASET_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid CID-to-dataset sidecar '{}': bad magic", path.display()),
            ));
        }
        let mut blocks = Vec::new();
        let mut group_count = 0usize;
        loop {
            let mut count_bytes = [0u8; 4];
            match file.read_exact(&mut count_bytes) {
                Ok(()) => {}
                Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => break,
                Err(err) => return Err(err),
            }
            let block_groups = u32::from_le_bytes(count_bytes) as usize;
            if block_groups == 0 {
                break;
            }
            let mut len_bytes = [0u8; 8];
            file.read_exact(&mut len_bytes)?;
            let data_len = usize::try_from(u64::from_le_bytes(len_bytes)).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "CID sidecar block is too large")
            })?;
            let offsets_file_offset = file.stream_position()?;
            let offsets_bytes = (block_groups + 1)
                .checked_mul(std::mem::size_of::<u32>())
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "CID sidecar offsets overflow")
                })?;
            let data_offset = offsets_file_offset
                .checked_add(offsets_bytes as u64)
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "CID sidecar offset overflow")
                })?;
            let next_block = data_offset.checked_add(data_len as u64).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "CID sidecar offset overflow")
            })?;
            blocks.push(CidDatasetSidecarBlock {
                first_group: group_count,
                group_count: block_groups,
                offsets_file_offset,
                data_offset,
                data_len,
            });
            group_count = group_count.checked_add(block_groups).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "CID sidecar group count overflow")
            })?;
            file.seek(SeekFrom::Start(next_block))?;
        }
        Ok(Self {
            blocks,
            group_count,
            cache: Mutex::new(CidDatasetSidecarCache {
                file: File::open(path)?,
                block_index: usize::MAX,
                offsets: Vec::new(),
                data: Vec::new(),
                decoded_ids: Vec::new(),
            }),
        })
    }

    pub(crate) fn len(&self) -> usize {
        self.group_count
    }

    fn reader(&self) -> io::Result<CidDatasetSidecarReader<'_>> {
        let file = self
            .cache
            .lock()
            .expect("CID sidecar cache lock poisoned")
            .file
            .try_clone()?;
        Ok(CidDatasetSidecarReader {
            sidecar: self,
            file,
            block_index: usize::MAX,
            offsets: Vec::new(),
            data: Vec::new(),
        })
    }

    fn merge_group_into(
        &self,
        group: usize,
        target: &mut Vec<u32>,
        dataset_offset: u32,
    ) -> io::Result<()> {
        let block_index = self
            .blocks
            .partition_point(|block| block.first_group <= group)
            .checked_sub(1)
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "CID sidecar group is out of range")
            })?;
        let block = self.blocks.get(block_index).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "CID sidecar group is out of range")
        })?;
        let local = group - block.first_group;
        if local >= block.group_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "CID sidecar group is outside its block",
            ));
        }
        let mut cache = self.cache.lock().expect("CID sidecar cache lock poisoned");
        if cache.block_index != block_index {
            cache.file.seek(SeekFrom::Start(block.offsets_file_offset))?;
            let mut offsets = vec![0u32; block.group_count + 1];
            let mut bytes = [0u8; 4];
            for offset in &mut offsets {
                cache.file.read_exact(&mut bytes)?;
                *offset = u32::from_le_bytes(bytes);
            }
            let mut compressed = vec![0u8; block.data_len];
            cache.file.seek(SeekFrom::Start(block.data_offset))?;
            cache.file.read_exact(&mut compressed)?;
            let data = zstd::decode_all(compressed.as_slice())?;
            if offsets.windows(2).any(|range| range[0] > range[1])
                || offsets.last().copied().map(|v| v as usize) != Some(data.len())
            {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CID sidecar block offsets are invalid for its payload",
                ));
            }
            cache.offsets = offsets;
            cache.data = data;
            cache.block_index = block_index;
        }
        let start = cache.offsets[local] as usize;
        let end = cache.offsets[local + 1] as usize;
        let mut ids = std::mem::take(&mut cache.decoded_ids);
        ids.clear();
        let encoded = &cache.data[start..end];
        let mut cursor = 0usize;
        let mut previous = 0u64;
        while cursor < encoded.len() {
            let delta = read_varint_field(encoded, &mut cursor, "CID dataset delta")?;
            previous = previous.checked_add(delta).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "CID sidecar dataset ID overflow")
            })?;
            let one_based = previous.checked_add(1).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "CID sidecar dataset ID overflow")
            })?;
            let id = u32::try_from(one_based).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "CID sidecar dataset ID overflow")
            })?;
            ids.push(id);
        }
        if ids.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "CID sidecar group has no dataset membership",
            ));
        }
        let result = merge_sorted_dataset_ids_with_offset(target, &ids, dataset_offset);
        cache.decoded_ids = ids;
        result
    }
}

impl CidDatasetSidecarReader<'_> {
    fn load_group_into(
        &mut self,
        group: usize,
        target: &mut Vec<u32>,
        dataset_offset: u32,
    ) -> io::Result<()> {
        let block_index = self
            .sidecar
            .blocks
            .partition_point(|block| block.first_group <= group)
            .checked_sub(1)
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "CID sidecar group is out of range")
            })?;
        let block = self.sidecar.blocks.get(block_index).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "CID sidecar group is out of range")
        })?;
        let local = group - block.first_group;
        if local >= block.group_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "CID sidecar group is outside its block",
            ));
        }
        if self.block_index != block_index {
            self.offsets.resize(block.group_count + 1, 0);
            let mut offset_bytes = vec![0u8; self.offsets.len() * std::mem::size_of::<u32>()];
            self.file
                .read_exact_at(&mut offset_bytes, block.offsets_file_offset)?;
            for (offset, bytes) in self.offsets.iter_mut().zip(offset_bytes.chunks_exact(4)) {
                *offset = u32::from_le_bytes(bytes.try_into().expect("four-byte offset"));
            }
            let mut compressed = vec![0u8; block.data_len];
            self.file.read_exact_at(&mut compressed, block.data_offset)?;
            self.data = zstd::decode_all(compressed.as_slice())?;
            if self.offsets.windows(2).any(|range| range[0] > range[1])
                || self.offsets.last().copied().map(|v| v as usize) != Some(self.data.len())
            {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CID sidecar block offsets are invalid for its payload",
                ));
            }
            self.block_index = block_index;
        }

        target.clear();
        let start = self.offsets[local] as usize;
        let end = self.offsets[local + 1] as usize;
        let encoded = &self.data[start..end];
        let mut cursor = 0usize;
        let mut previous = 0u64;
        while cursor < encoded.len() {
            let delta = read_varint_field(encoded, &mut cursor, "CID dataset delta")?;
            previous = previous.checked_add(delta).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "CID sidecar dataset ID overflow")
            })?;
            let one_based = previous
                .checked_add(1)
                .and_then(|id| id.checked_add(dataset_offset as u64))
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "CID sidecar dataset ID overflow")
                })?;
            target.push(u32::try_from(one_based).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "CID sidecar dataset ID overflow")
            })?);
        }
        if target.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "CID sidecar group has no dataset membership",
            ));
        }
        Ok(())
    }
}

pub(crate) struct CidDatasetSidecarWriter {
    file: BufWriter<File>,
    block_groups: usize,
    block_offsets: Vec<u32>,
    block_data: Vec<u8>,
}

impl CidDatasetSidecarWriter {
    pub(crate) fn create(path: &Path) -> io::Result<Self> {
        let mut file = BufWriter::with_capacity(IO_BUFFER_CAPACITY, File::create(path)?);
        file.write_all(CID_TO_DATASET_MAGIC)?;
        Ok(Self {
            file,
            block_groups: 0,
            block_offsets: vec![0],
            block_data: Vec::new(),
        })
    }

    pub(crate) fn append_zero_based(&mut self, dataset_ids: &[usize]) -> io::Result<()> {
        if dataset_ids.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "cannot store an empty CID dataset set",
            ));
        }
        let mut previous = 0usize;
        for (index, &dataset) in dataset_ids.iter().enumerate() {
            if index > 0 && dataset <= previous {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CID dataset IDs are not strictly increasing",
                ));
            }
            write_varint_u64((dataset - previous) as u64, &mut self.block_data);
            previous = dataset;
        }
        self.finish_group()
    }

    pub(crate) fn append_one_based(&mut self, dataset_ids: &[u32]) -> io::Result<()> {
        if dataset_ids.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "cannot store an empty CID dataset set",
            ));
        }
        let mut previous = 0u32;
        let mut first = true;
        for &dataset in dataset_ids {
            let zero_based = dataset.checked_sub(1).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "dataset ID zero is invalid")
            })?;
            if !first && zero_based <= previous {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CID dataset IDs are not strictly increasing",
                ));
            }
            write_varint_u64((zero_based - previous) as u64, &mut self.block_data);
            previous = zero_based;
            first = false;
        }
        self.finish_group()
    }

    fn finish_group(&mut self) -> io::Result<()> {
        self.block_offsets
            .push(u32::try_from(self.block_data.len()).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "CID sidecar block exceeds 4GiB")
            })?);
        self.block_groups += 1;
        if self.block_groups >= BUCKET_SIZE_BLOCK_MAX_GROUPS
            || self.block_data.len() >= BUCKET_SIZE_BLOCK_MAX_UNCOMPRESSED_BYTES
        {
            self.flush_block()?;
        }
        Ok(())
    }

    fn flush_block(&mut self) -> io::Result<()> {
        if self.block_groups == 0 {
            return Ok(());
        }
        let compressed = zstd::encode_all(self.block_data.as_slice(), 1)?;
        self.file
            .write_all(&(self.block_groups as u32).to_le_bytes())?;
        self.file.write_all(&(compressed.len() as u64).to_le_bytes())?;
        for &offset in &self.block_offsets {
            self.file.write_all(&offset.to_le_bytes())?;
        }
        self.file.write_all(&compressed)?;
        self.block_groups = 0;
        self.block_offsets.clear();
        self.block_offsets.push(0);
        self.block_data.clear();
        Ok(())
    }

    pub(crate) fn finish(mut self) -> io::Result<()> {
        self.flush_block()?;
        self.file.write_all(&0u32.to_le_bytes())?;
        self.file.flush()
    }
}

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

#[derive(Clone, Copy, Debug)]
struct CompressionMemoryBudget {
    total_bytes: usize,
}

impl CompressionMemoryBudget {
    fn from_gb(memory_gb: usize) -> Self {
        Self {
            total_bytes: memory_gb
                .max(1)
                .saturating_mul(1024 * 1024 * 1024)
                .max(MIN_MEMORY_BUDGET_BYTES),
        }
    }

    fn group_batch_bytes(self) -> usize {
        (self.total_bytes / 8).clamp(MIN_GROUP_BATCH_BYTES, MAX_GROUP_BATCH_BYTES)
    }

    fn color_chunk_bytes(self) -> usize {
        (self.total_bytes / 16).clamp(MIN_SORT_CHUNK_BYTES, COLOR_CHUNK_TARGET_BYTES)
    }

    fn subset_window_sequence_bytes(self) -> usize {
        (self.total_bytes / 16).clamp(8 * 1024 * 1024, 1024 * 1024 * 1024)
    }

    fn subset_window_count(self, dataset_count: usize) -> usize {
        let worst_case_subset_bytes = dataset_count
            .saturating_mul(std::mem::size_of::<u32>())
            .saturating_add(96)
            .max(1);
        (self.total_bytes / 8 / worst_case_subset_bytes).max(1)
    }

    fn ggcat_memory_gb(self) -> f64 {
        ((self.total_bytes as f64) / (1024.0 * 1024.0 * 1024.0) * 0.75).max(0.25)
    }
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
    fn record_batch(&mut self, groups: usize, timing: PhaseTiming) {
        self.groups += groups;
        self.timing.add(timing);
        let wall_sec = timing.wall_sec;
        let cpu_sec = timing.cpu_sec;
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

fn flush_sorted_chunk(
    records: &mut Vec<SortedColorRecord>,
    chunk_files: &mut Vec<PathBuf>,
    chunk_dir: &Path,
    sort_pool: &rayon::ThreadPool,
) -> Result<()> {
    if records.is_empty() {
        return Ok(());
    }

    if records.len() >= PAR_SORT_THRESHOLD {
        // GGCAT can have every worker in the global Rayon pool blocked while its
        // bounded structured-output channel is full.  Sorting on that same pool
        // would then deadlock the consumer which is needed to drain the channel.
        sort_pool.install(|| {
            records.par_sort_unstable_by(|left, right| left.subset.cmp(&right.subset));
        });
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
    sort_pool: &rayon::ThreadPool,
    metrics: &mut ChunkFlushMetrics,
) -> Result<()> {
    if records.is_empty() {
        return Ok(());
    }
    let records_count = records.len();
    let timer = PhaseTimer::start();
    flush_sorted_chunk(records, chunk_files, chunk_dir, sort_pool)?;
    metrics.calls += 1;
    metrics.records += records_count;
    metrics.timing.add(timer.finish());
    Ok(())
}

fn split_color_run_segments(
    runs: Vec<(ColorIndexType, usize)>,
    seq: &[u8],
    k: usize,
    source: &str,
) -> Result<Vec<SortedColorRecord>> {
    if runs.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("ggcat emitted colored sequence without color runs: {source}"),
        ));
    }

    if seq.len() < k {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("ggcat sequence shorter than k (len={}, k={k}) for {source}", seq.len()),
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
                "GGCAT color-run lengths mismatch sequence k-mers for {source}: runs={}, expected={}",
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
                    "GGCAT color run exceeds sequence k-mers for {source}: end={}, total={}",
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
                "GGCAT color runs do not fully cover sequence k-mers for {source}: covered={}, total={}",
                run_start,
                seq_kmers
            ),
        ));
    }

    Ok(segments)
}

fn read_u64_field(block: &[u8], cursor: &mut usize, field: &str) -> Result<u64> {
    let end = cursor
        .checked_add(8)
        .ok_or_else(|| io::Error::other("structured-output cursor overflow"))?;
    let bytes = block.get(*cursor..end).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("truncated GGCAT structured-output {field}"),
        )
    })?;
    *cursor = end;
    Ok(u64::from_le_bytes(bytes.try_into().unwrap()))
}

fn read_varint_field(data: &[u8], cursor: &mut usize, field: &str) -> Result<u64> {
    let mut value = 0u64;
    for shift in (0..=63).step_by(7) {
        let byte = *data.get(*cursor).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("truncated GGCAT structured-output {field}"),
            )
        })?;
        *cursor += 1;
        if shift == 63 && byte > 1 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("overflow in GGCAT structured-output {field}"),
            ));
        }
        value |= ((byte & 0x7f) as u64) << shift;
        if byte & 0x80 == 0 {
            return Ok(value);
        }
    }
    Err(io::Error::new(
        io::ErrorKind::InvalidData,
        format!("unterminated GGCAT structured-output {field}"),
    ))
}

fn visit_structured_ggcat_block(
    block: &[u8],
    k: usize,
    mut emit: impl FnMut(SortedColorRecord) -> Result<()>,
) -> Result<(usize, usize)> {
    let mut cursor = 0usize;
    let mut sequences = 0usize;
    let mut segments = 0usize;
    while cursor < block.len() {
        let sequence_index = read_u64_field(block, &mut cursor, "sequence index")?;
        let sequence_len = usize::try_from(read_u64_field(
            block,
            &mut cursor,
            "sequence length",
        )?)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "sequence length overflow"))?;
        let color_len = usize::try_from(read_u64_field(block, &mut cursor, "color length")?)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "color length overflow"))?;
        let links_len = usize::try_from(read_u64_field(block, &mut cursor, "links length")?)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "links length overflow"))?;
        let record_len = sequence_len
            .checked_add(color_len)
            .and_then(|len| len.checked_add(links_len))
            .ok_or_else(|| io::Error::other("structured-output record length overflow"))?;
        let record_end = cursor
            .checked_add(record_len)
            .ok_or_else(|| io::Error::other("structured-output cursor overflow"))?;
        if record_end > block.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("truncated GGCAT structured-output sequence {sequence_index}"),
            ));
        }
        let seq = &block[cursor..cursor + sequence_len];
        cursor += sequence_len;
        let color_data = &block[cursor..cursor + color_len];
        cursor += color_len;
        cursor += links_len;

        let mut color_cursor = 0usize;
        let colors_count = usize::try_from(read_varint_field(
            color_data,
            &mut color_cursor,
            "color count",
        )?)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "color count overflow"))?;
        let mut runs = Vec::with_capacity(colors_count);
        for _ in 0..colors_count {
            let subset = ColorIndexType::try_from(read_varint_field(
                color_data,
                &mut color_cursor,
                "color subset",
            )?)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "color subset overflow"))?;
            let count = usize::try_from(read_varint_field(
                color_data,
                &mut color_cursor,
                "color run length",
            )?)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "color run overflow"))?;
            if count > 0 {
                runs.push((subset, count));
            }
        }
        if color_cursor != color_data.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "unused bytes in GGCAT structured-output colors for sequence {sequence_index}"
                ),
            ));
        }
        let parsed = split_color_run_segments(
            runs,
            seq,
            k,
            &format!("structured sequence {sequence_index}"),
        )?;
        sequences += 1;
        segments += parsed.len();
        for record in parsed {
            emit(record)?;
        }
    }
    Ok((sequences, segments))
}

fn merge_sorted_dataset_ids_with_offset(
    target: &mut Vec<u32>,
    other: &[u32],
    offset: u32,
) -> io::Result<()> {
    if offset != 0 && other.last().copied().unwrap_or(0) > u32::MAX - offset {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "merged dataset id overflow",
        ));
    }
    if target.is_empty() {
        target.reserve(other.len());
        target.extend(other.iter().map(|&id| id + offset));
        return Ok(());
    }
    if other.is_empty() {
        return Ok(());
    }
    let mut merged = Vec::with_capacity(target.len().saturating_add(other.len()));
    let (mut left, mut right) = (0usize, 0usize);
    while left < target.len() && right < other.len() {
        let right_id = other[right] + offset;
        match target[left].cmp(&right_id) {
            CmpOrdering::Less => {
                merged.push(target[left]);
                left += 1;
            }
            CmpOrdering::Greater => {
                merged.push(right_id);
                right += 1;
            }
            CmpOrdering::Equal => {
                merged.push(target[left]);
                left += 1;
                right += 1;
            }
        }
    }
    merged.extend_from_slice(&target[left..]);
    merged.extend(other[right..].iter().map(|&id| id + offset));
    *target = merged;
    Ok(())
}

fn merge_sorted_dataset_ids(target: &mut Vec<u32>, other: &[u32]) {
    merge_sorted_dataset_ids_with_offset(target, other, 0)
        .expect("zero-offset dataset IDs cannot overflow");
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct SourceDatasetSpan {
    storage: usize,
    start: usize,
    end: usize,
    offset: u32,
}

#[derive(Debug)]
pub(crate) struct SourceDatasetMap {
    storages: Vec<Arc<Vec<u32>>>,
    spans: Vec<SourceDatasetSpan>,
    dense_ranges: Vec<DenseSourceDatasetRange>,
    disk_ranges: Vec<DiskSourceDatasetRange>,
    source_count: usize,
}

#[derive(Debug)]
struct DenseSourceDatasetRange {
    source_start: usize,
    source_end: usize,
    storage: usize,
    offsets: Arc<Vec<usize>>,
    offset: u32,
}

#[derive(Debug)]
struct DiskSourceDatasetRange {
    source_start: usize,
    source_end: usize,
    sidecar: Arc<CidDatasetSidecar>,
    offset: u32,
}

struct SourceDatasetReader<'a> {
    map: &'a SourceDatasetMap,
    disk_readers: Vec<CidDatasetSidecarReader<'a>>,
}

impl SourceDatasetMap {
    pub(crate) fn new() -> Self {
        Self {
            storages: Vec::new(),
            spans: Vec::new(),
            dense_ranges: Vec::new(),
            disk_ranges: Vec::new(),
            source_count: 0,
        }
    }

    pub(crate) fn add_storage(&mut self, values: Arc<Vec<u32>>) -> usize {
        let index = self.storages.len();
        self.storages.push(values);
        index
    }

    pub(crate) fn push_span(
        &mut self,
        storage: usize,
        start: usize,
        end: usize,
        offset: u32,
    ) -> io::Result<()> {
        let values = self.storages.get(storage).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "source dataset storage is missing")
        })?;
        if start > end || end > values.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "source dataset span is outside its storage",
            ));
        }
        self.spans.push(SourceDatasetSpan {
            storage,
            start,
            end,
            offset,
        });
        self.source_count = self.source_count.checked_add(1).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "too many GGCAT sources")
        })?;
        Ok(())
    }

    pub(crate) fn add_dense_storage(
        &mut self,
        values: Arc<Vec<u32>>,
        offsets: Arc<Vec<usize>>,
        offset: u32,
    ) -> io::Result<()> {
        if offsets.is_empty() || offsets[0] != 0 || offsets.last().copied() != Some(values.len()) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "dense source offsets do not cover their dataset storage",
            ));
        }
        if offsets.windows(2).any(|range| range[0] >= range[1]) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "archive contains a source group without any dataset membership",
            ));
        }
        let count = offsets.len() - 1;
        let source_start = self.source_count;
        let source_end = source_start.checked_add(count).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "too many GGCAT sources")
        })?;
        let storage = self.add_storage(values);
        self.dense_ranges.push(DenseSourceDatasetRange {
            source_start,
            source_end,
            storage,
            offsets,
            offset,
        });
        self.source_count = source_end;
        Ok(())
    }

    pub(crate) fn add_disk_storage(
        &mut self,
        sidecar: Arc<CidDatasetSidecar>,
        offset: u32,
    ) -> io::Result<()> {
        let source_start = self.source_count;
        let source_end = source_start.checked_add(sidecar.len()).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "too many GGCAT sources")
        })?;
        self.disk_ranges.push(DiskSourceDatasetRange {
            source_start,
            source_end,
            sidecar,
            offset,
        });
        self.source_count = source_end;
        Ok(())
    }

    fn from_singletons(dataset_count: usize) -> io::Result<Self> {
        let mut values = Vec::with_capacity(dataset_count);
        for dataset in 1..=dataset_count {
            values.push(u32::try_from(dataset).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidInput, "too many input datasets")
            })?);
        }
        let mut map = Self::new();
        let storage = map.add_storage(Arc::new(values));
        for dataset in 0..dataset_count {
            map.push_span(storage, dataset, dataset + 1, 0)?;
        }
        Ok(map)
    }

    pub(crate) fn len(&self) -> usize {
        self.source_count
    }

    fn reader(&self) -> io::Result<SourceDatasetReader<'_>> {
        let mut disk_readers = Vec::with_capacity(self.disk_ranges.len());
        for range in &self.disk_ranges {
            disk_readers.push(range.sidecar.reader()?);
        }
        Ok(SourceDatasetReader {
            map: self,
            disk_readers,
        })
    }

    fn merge_source_into(&self, source: usize, target: &mut Vec<u32>) -> io::Result<()> {
        if let Some(span) = self.spans.get(source) {
            let values = self.storages.get(span.storage).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "source dataset storage is missing")
            })?;
            let values = values.get(span.start..span.end).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "source dataset span is invalid")
            })?;
            return merge_sorted_dataset_ids_with_offset(target, values, span.offset);
        }
        if let Some(range) = self
            .disk_ranges
            .iter()
            .find(|range| source >= range.source_start && source < range.source_end)
        {
            return range.sidecar.merge_group_into(
                source - range.source_start,
                target,
                range.offset,
            );
        }
        let range = self
            .dense_ranges
            .iter()
            .find(|range| source >= range.source_start && source < range.source_end)
            .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("ggcat source index {} is out of range", source),
            )
        })?;
        let local = source - range.source_start;
        let start = range.offsets[local];
        let end = range.offsets[local + 1];
        let values = self.storages.get(range.storage).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "source dataset storage is missing")
        })?;
        let values = values.get(start..end).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "source dataset span is invalid")
        })?;
        merge_sorted_dataset_ids_with_offset(target, values, range.offset)
    }
}

impl SourceDatasetReader<'_> {
    fn load_source_into(&mut self, source: usize, target: &mut Vec<u32>) -> io::Result<()> {
        if let Some(span) = self.map.spans.get(source) {
            let values = self.map.storages.get(span.storage).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "source dataset storage is missing")
            })?;
            let values = values.get(span.start..span.end).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "source dataset span is invalid")
            })?;
            target.clear();
            target.reserve(values.len());
            for &value in values {
                target.push(value.checked_add(span.offset).ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "dataset ID overflow")
                })?);
            }
            return Ok(());
        }
        if let Some((index, range)) = self
            .map
            .disk_ranges
            .iter()
            .enumerate()
            .find(|(_, range)| source >= range.source_start && source < range.source_end)
        {
            return self.disk_readers[index].load_group_into(
                source - range.source_start,
                target,
                range.offset,
            );
        }
        let range = self
            .map
            .dense_ranges
            .iter()
            .find(|range| source >= range.source_start && source < range.source_end)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("ggcat source index {} is out of range", source),
                )
            })?;
        let local = source - range.source_start;
        let start = range.offsets[local];
        let end = range.offsets[local + 1];
        let values = self.map.storages.get(range.storage).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "source dataset storage is missing")
        })?;
        let values = values.get(start..end).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "source dataset span is invalid")
        })?;
        target.clear();
        target.reserve(values.len());
        for &value in values {
            target.push(value.checked_add(range.offset).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "dataset ID overflow")
            })?);
        }
        Ok(())
    }
}

struct ResolvedSubsets {
    subsets: Vec<ColorIndexType>,
    dataset_ids: Vec<Arc<Vec<u32>>>,
}

fn resolve_subsets_to_dataset_ids(
    instance: &GGCATInstance,
    colormap_file: &Path,
    subsets: &[ColorIndexType],
    color_index_to_source_index: Option<&[usize]>,
    source_dataset_ids: &SourceDatasetMap,
    dataset_count: usize,
) -> Result<ResolvedSubsets> {
    if subsets.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "no color subsets were parsed from GGCAT output",
        ));
    }

    let mut subset_list = subsets.to_vec();
    subset_list.sort_unstable();
    subset_list.dedup();

    let callback_err = Mutex::new(None::<io::Error>);

    for chunk in subset_list.chunks(subset_list.len().max(1)) {
        let subset_sources = Mutex::new(Vec::<(ColorIndexType, Vec<usize>)>::with_capacity(
            chunk.len(),
        ));
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

                    let mut source_indices = Vec::with_capacity(colors.len());
                    for &color in colors {
                        let color_idx = color as usize;
                        let source_idx = match color_index_to_source_index {
                            Some(mapping) => match mapping.get(color_idx) {
                                Some(&source_idx) => source_idx,
                                None => {
                                    *callback_err.lock().expect("callback error lock poisoned") =
                                        Some(io::Error::new(
                                            io::ErrorKind::InvalidData,
                                            format!(
                                                "ggcat color index {} is missing from its source mapping",
                                                color_idx
                                            ),
                                        ));
                                    return;
                                }
                            },
                            None => color_idx,
                        };
                        if source_idx >= source_dataset_ids.len() {
                            *callback_err.lock().expect("callback error lock poisoned") =
                                Some(io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    format!("ggcat color index {} out of range", color_idx),
                                ));
                            return;
                        }
                        source_indices.push(source_idx);
                    }
                    if source_indices.is_empty() {
                        *callback_err.lock().expect("callback error lock poisoned") =
                            Some(io::Error::new(
                                io::ErrorKind::InvalidData,
                                "GGCAT color subset resolved to empty dataset set",
                            ));
                        return;
                    }
                    subset_sources
                        .lock()
                        .expect("subset-source lock poisoned")
                        .push((subset, source_indices));
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

        let subset_sources = subset_sources
            .into_inner()
            .expect("subset-source lock poisoned");
        if subset_sources.len() != chunk.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "queried {} color subsets but GGCAT returned {}",
                    chunk.len(),
                    subset_sources.len()
                ),
            ));
        }
        let mut subset_sources = subset_sources;
        subset_sources.sort_unstable_by_key(|(subset, _)| *subset);
        if subset_sources
            .iter()
            .map(|(subset, _)| *subset)
            .ne(chunk.iter().copied())
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "GGCAT returned unexpected color subset IDs",
            ));
        }
        let request_count: usize = subset_sources
            .iter()
            .map(|(_, sources)| sources.len())
            .sum();
        let mut requests = Vec::<u64>::with_capacity(request_count);
        for (subset_index, (_, sources)) in subset_sources.into_iter().enumerate() {
            let subset_index = u32::try_from(subset_index).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "too many subsets in one window")
            })?;
            for source in sources {
                let source = u32::try_from(source).map_err(|_| {
                    io::Error::new(io::ErrorKind::InvalidData, "source index exceeds u32")
                })?;
                requests.push(((source as u64) << 32) | subset_index as u64);
            }
        }
        requests.par_sort_unstable();
        requests.dedup();
        let partial: Vec<FastMutex<Vec<u32>>> = (0..chunk.len())
            .map(|_| FastMutex::new(Vec::new()))
            .collect();
        let target_tasks = rayon::current_num_threads().max(1).saturating_mul(4);
        let request_chunk_size = requests.len().div_ceil(target_tasks).max(1);
        requests
            .par_chunks(request_chunk_size)
            .try_for_each(|request_chunk| -> Result<()> {
                let mut reader = source_dataset_ids.reader()?;
                let mut current_source = None;
                let mut current_dataset_ids = Vec::new();
                for &request in request_chunk {
                    let source = (request >> 32) as usize;
                    let subset_index = request as u32 as usize;
                    if current_source != Some(source) {
                        reader.load_source_into(source, &mut current_dataset_ids)?;
                        current_source = Some(source);
                    }
                    let mut dataset_ids = partial[subset_index].lock();
                    dataset_ids.extend_from_slice(&current_dataset_ids);
                }
                Ok(())
            })?;

        let chunk_dataset_ids = partial
            .into_par_iter()
            .map(|dataset_ids| -> Result<Arc<Vec<u32>>> {
                let mut dataset_ids = dataset_ids.into_inner();
                dataset_ids.sort_unstable();
                dataset_ids.dedup();
                if dataset_ids.is_empty() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "GGCAT color subset resolved to empty dataset set",
                    ));
                }
                if let Some(&invalid) = dataset_ids
                    .iter()
                    .find(|&&dataset_id| dataset_id == 0 || dataset_id as usize > dataset_count)
                {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "source color maps to dataset id {} outside [1, {}]",
                            invalid, dataset_count
                        ),
                    ));
                }
                Ok(Arc::new(dataset_ids))
            })
            .collect::<Result<Vec<_>>>()?;
        debug_assert_eq!(chunk_dataset_ids.len(), chunk.len());
        return Ok(ResolvedSubsets {
            subsets: chunk.to_vec(),
            dataset_ids: chunk_dataset_ids,
        });
    }
    unreachable!("non-empty subset list is processed as one chunk")
}

#[allow(clippy::too_many_arguments)]
fn emit_resolved_record_window(
    records: &mut Vec<SortedColorRecord>,
    instance: &GGCATInstance,
    colormap_file: &Path,
    color_index_to_source_index: Option<&[usize]>,
    source_dataset_ids: &SourceDatasetMap,
    dataset_count: usize,
    sender: &mpsc::SyncSender<SimplitigBatch>,
    output_batch: &mut SimplitigBatch,
    output_batch_bytes: &mut usize,
    output_batch_byte_limit: usize,
    emitted_batches_count: &mut usize,
    emitted_records_count: &mut usize,
) -> Result<()> {
    if records.is_empty() {
        return Ok(());
    }
    let mut subsets = Vec::new();
    let mut previous = None;
    for record in records.iter() {
        if previous != Some(record.subset) {
            subsets.push(record.subset);
            previous = Some(record.subset);
        }
    }
    let resolved = resolve_subsets_to_dataset_ids(
        instance,
        colormap_file,
        &subsets,
        color_index_to_source_index,
        source_dataset_ids,
        dataset_count,
    )?;
    let mut resolved_index = 0usize;
    for record in records.drain(..) {
        while resolved
            .subsets
            .get(resolved_index)
            .is_some_and(|&subset| subset < record.subset)
        {
            resolved_index += 1;
        }
        if resolved.subsets.get(resolved_index).copied() != Some(record.subset) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("missing resolved dataset ids for subset {}", record.subset),
            ));
        }
        let color_ids = resolved.dataset_ids.get(resolved_index).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("missing resolved dataset ids for subset {}", record.subset),
            )
        })?;
        *output_batch_bytes = output_batch_bytes.saturating_add(record.seq.len());
        output_batch.push(SimplitigRecord {
            color_ids: Arc::clone(color_ids),
            seq: record.seq,
        });
        if output_batch.len() >= COLOR_RECORD_BATCH_SIZE
            || *output_batch_bytes >= output_batch_byte_limit
        {
            let out = std::mem::take(output_batch);
            *emitted_batches_count += 1;
            *emitted_records_count += out.len();
            sender.send(out).map_err(|err| {
                io::Error::new(
                    io::ErrorKind::BrokenPipe,
                    format!("record receiver dropped: {err}"),
                )
            })?;
            *output_batch_bytes = 0;
            *output_batch = Vec::with_capacity(COLOR_RECORD_BATCH_SIZE);
        }
    }
    if !output_batch.is_empty() {
        let out = std::mem::take(output_batch);
        *emitted_batches_count += 1;
        *emitted_records_count += out.len();
        sender.send(out).map_err(|err| {
            io::Error::new(
                io::ErrorKind::BrokenPipe,
                format!("record receiver dropped: {err}"),
            )
        })?;
        *output_batch_bytes = 0;
        *output_batch = Vec::with_capacity(COLOR_RECORD_BATCH_SIZE);
    }
    Ok(())
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

struct DatasetCidSpill {
    paths: Vec<PathBuf>,
    writers: Vec<BufWriter<File>>,
    last_cids: Vec<u64>,
    directory: PathBuf,
    datasets_per_partition: usize,
    dataset_count: usize,
}

struct DatasetCidSpillFiles {
    paths: Vec<PathBuf>,
    directory: PathBuf,
    datasets_per_partition: usize,
    dataset_count: usize,
}

impl DatasetCidSpill {
    fn new(output_dir: &str, dataset_count: usize) -> Result<Self> {
        let directory = create_id_cid_spill_dir(output_dir)?;
        let partition_count = dataset_count
            .div_ceil(DATASETS_PER_CID_PARTITION)
            .clamp(1, MAX_CID_PARTITIONS);
        let datasets_per_partition = dataset_count.div_ceil(partition_count).max(1);
        let mut paths = Vec::with_capacity(partition_count);
        let mut writers = Vec::with_capacity(partition_count);
        let last_cids = vec![0u64; partition_count];
        for partition in 0..partition_count {
            let path = directory.join(format!("partition_{partition:04}.pairs.bin"));
            writers.push(BufWriter::with_capacity(
                CID_PARTITION_BUFFER_BYTES,
                File::create(&path)?,
            ));
            paths.push(path);
        }
        println!(
            "Dataset-to-CID spill: datasets={}, partitions={}, datasets_per_partition={}",
            dataset_count, partition_count, datasets_per_partition
        );
        Ok(Self {
            paths,
            writers,
            last_cids,
            directory,
            datasets_per_partition,
            dataset_count,
        })
    }

    fn append_group(&mut self, dataset_ids: &[usize], cid: usize) -> Result<()> {
        let cid = u64::try_from(cid).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "CID cannot be represented as u64",
            )
        })?;
        let mut start = 0usize;
        while start < dataset_ids.len() {
            let first_dataset = dataset_ids[start];
            if first_dataset >= self.dataset_count {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "dataset id {} outside [0, {}) while spilling CID {}",
                        first_dataset, self.dataset_count, cid
                    ),
                ));
            }
            let partition =
                (first_dataset / self.datasets_per_partition).min(self.writers.len() - 1);
            let mut end = start + 1;
            while end < dataset_ids.len()
                && dataset_ids[end] / self.datasets_per_partition == partition
            {
                end += 1;
            }

            let cid_delta = cid.checked_sub(self.last_cids[partition]).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CID order regressed while writing dataset transpose spill",
                )
            })?;
            let writer = &mut self.writers[partition];
            write_varint_u64_to_writer(cid_delta, &mut *writer)?;
            write_varint_u64_to_writer((end - start) as u64, &mut *writer)?;
            let partition_start = partition.saturating_mul(self.datasets_per_partition);
            let mut previous_dataset = partition_start;
            for &dataset in &dataset_ids[start..end] {
                if dataset >= self.dataset_count || dataset < previous_dataset {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "dataset IDs are not sorted within a color set",
                    ));
                }
                write_varint_u64_to_writer((dataset - previous_dataset) as u64, &mut *writer)?;
                previous_dataset = dataset;
            }
            self.last_cids[partition] = cid;
            start = end;
        }
        Ok(())
    }

    fn finalize(mut self) -> Result<DatasetCidSpillFiles> {
        for writer in &mut self.writers {
            writer.flush()?;
        }
        drop(self.writers);
        Ok(DatasetCidSpillFiles {
            paths: self.paths,
            directory: self.directory,
            datasets_per_partition: self.datasets_per_partition,
            dataset_count: self.dataset_count,
        })
    }
}

struct PositionSpill {
    path: PathBuf,
    writer: BufWriter<File>,
    entries: u64,
    previous_tigs: u64,
    previous_sizes: u64,
}

impl PositionSpill {
    fn new(directory: &Path) -> Result<Self> {
        let path = directory.join("positions.deltas.tmp");
        let writer = BufWriter::with_capacity(IO_BUFFER_CAPACITY, File::create(&path)?);
        let mut spill = Self {
            path,
            writer,
            entries: 0,
            previous_tigs: 0,
            previous_sizes: 0,
        };
        spill.append(0, 0)?;
        Ok(spill)
    }

    fn append(&mut self, tigs: u64, sizes: u64) -> Result<()> {
        if tigs < self.previous_tigs || sizes < self.previous_sizes {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "positions are not monotonic",
            ));
        }
        write_varint_u64_to_writer(tigs - self.previous_tigs, &mut self.writer)?;
        write_varint_u64_to_writer(sizes - self.previous_sizes, &mut self.writer)?;
        self.previous_tigs = tigs;
        self.previous_sizes = sizes;
        self.entries += 1;
        Ok(())
    }

    fn finalize(mut self) -> Result<(PathBuf, u64)> {
        self.writer.flush()?;
        Ok((self.path, self.entries))
    }
}

struct StreamFinalize {
    position_path: PathBuf,
    position_entries: u64,
    cid_spill: DatasetCidSpillFiles,
}

struct StreamWriterState {
    omni_file: BufWriter<File>,
    size_file: BufWriter<File>,
    cid_dataset_writer: Option<CidDatasetSidecarWriter>,
    cid_spill: DatasetCidSpill,
    position_spill: PositionSpill,
    size_block_groups: usize,
    size_block_offsets: Vec<u32>,
    size_block_uncompressed: Vec<u8>,
    prev_tigs_size: u64,
    prev_bucket_pos: u64,
    cid: usize,
    encoded_seq_buffer: Vec<u8>,
}

impl StreamWriterState {
    fn new(unitigs_file_path: String, output_dir: &str, nb_files: usize) -> Result<Self> {
        let cid_spill = DatasetCidSpill::new(output_dir, nb_files)?;
        let position_spill = PositionSpill::new(&cid_spill.directory)?;
        let cid_dataset_writer = CidDatasetSidecarWriter::create(
            &Path::new(output_dir).join(CID_TO_DATASET_FILE),
        )?;
        Ok(Self {
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
            cid_dataset_writer: Some(cid_dataset_writer),
            cid_spill,
            position_spill,
            size_block_groups: 0,
            size_block_offsets: vec![0],
            size_block_uncompressed: Vec::new(),
            prev_tigs_size: 0,
            prev_bucket_pos: 0,
            cid: 0,
            encoded_seq_buffer: Vec::with_capacity(ENCODED_SEQ_BUFFER_TARGET),
        })
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

    fn append_group_cids(&mut self, dataset_ids_zero_based: &[usize]) -> Result<()> {
        self.cid_dataset_writer
            .as_mut()
            .expect("CID dataset sidecar writer must be initialized")
            .append_zero_based(dataset_ids_zero_based)?;
        self.cid_spill
            .append_group(dataset_ids_zero_based, self.cid)?;
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
        self.position_spill
            .append(self.prev_tigs_size, self.prev_bucket_pos)?;

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
        self.position_spill
            .append(self.prev_tigs_size, self.prev_bucket_pos)?;
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
        self.position_spill
            .append(self.prev_tigs_size, self.prev_bucket_pos)?;
        self.append_bucket_group_sizes_payload(group_sizes)?;

        self.append_group_cids(dataset_ids_zero_based)?;
        Ok(())
    }

    fn finalize(mut self) -> Result<StreamFinalize> {
        if !self.encoded_seq_buffer.is_empty() {
            self.omni_file.write_all(&self.encoded_seq_buffer)?;
            self.encoded_seq_buffer.clear();
        }
        self.flush_size_block()?;
        self.omni_file.flush()?;
        self.size_file.flush()?;
        self.cid_dataset_writer
            .take()
            .expect("CID dataset sidecar writer must be initialized")
            .finish()?;

        println!(
            "Completed compression: total tigs={}, total sizes={}",
            self.prev_tigs_size, self.prev_bucket_pos
        );
        let (position_path, position_entries) = self.position_spill.finalize()?;
        let cid_spill = self.cid_spill.finalize()?;
        Ok(StreamFinalize {
            position_path,
            position_entries,
            cid_spill,
        })
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

    Ok(GroupWriteResult {
        cid: task.cid,
        dataset_ids_zero_based: task.dataset_ids_zero_based,
        data,
    })
}

fn commit_group_result(
    next: GroupWriteResult,
    writer: &mut StreamWriterState,
    commit_stats: &mut CommitStats,
) -> Result<()> {
    if next.cid != writer.cid {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "out-of-order group result: expected CID {}, received {}",
                writer.cid, next.cid
            ),
        ));
    }
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
    Ok(())
}

fn process_and_commit_group_batch(
    tasks: &mut Vec<GroupTask>,
    resident_bytes: &mut usize,
    spill_dir: &Path,
    pool: &rayon::ThreadPool,
    writer: &mut StreamWriterState,
    worker_stats: &mut WorkerStats,
    commit_stats: &mut CommitStats,
) -> Result<()> {
    if tasks.is_empty() {
        return Ok(());
    }
    let tasks = std::mem::take(tasks);
    *resident_bytes = 0;
    let groups_count = tasks.len();
    let worker_timer = PhaseTimer::start();
    let results = pool.install(|| {
        tasks
            .into_par_iter()
            .map(|task| process_single_group_task(task, spill_dir))
            .collect::<io::Result<Vec<_>>>()
    })?;
    worker_stats.record_batch(groups_count, worker_timer.finish());
    let commit_timer = PhaseTimer::start();
    for result in results {
        commit_group_result(result, writer, commit_stats)?;
    }
    commit_stats.timing.add(commit_timer.finish());
    Ok(())
}

fn group_task_resident_bytes(task: &GroupTask) -> usize {
    task.dataset_ids_zero_based
        .capacity()
        .saturating_mul(std::mem::size_of::<usize>())
        .saturating_add(
            task.seqs
                .iter()
                .map(|seq| {
                    seq.capacity()
                        .saturating_add(std::mem::size_of::<Vec<u8>>())
                })
                .sum::<usize>(),
        )
        .saturating_add(
            task.run_paths
                .capacity()
                .saturating_mul(std::mem::size_of::<PathBuf>()),
        )
}

fn write_compressed_from_stream(
    unitigs_file_path: String,
    output_dir: &String,
    nb_files: u32,
    worker_threads: usize,
    memory_budget: CompressionMemoryBudget,
    record_rx: mpsc::Receiver<SimplitigBatch>,
) -> Result<StreamFinalize> {
    let total_timer = PhaseTimer::start();
    let mut writer = StreamWriterState::new(unitigs_file_path, output_dir, nb_files as usize)?;

    let group_spill_dir = Path::new(output_dir).join(format!(
        ".kloe-group-sort-spill-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0)
    ));
    fs::create_dir_all(&group_spill_dir)?;

    let worker_count = worker_threads.max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(worker_count)
        .build()
        .map_err(|err| io::Error::other(format!("create bounded group worker pool: {err}")))?;
    let mut worker_stats = WorkerStats::default();
    let mut commit_stats = CommitStats::default();
    let mut current_color_ids: Option<Arc<Vec<u32>>> = None;
    let mut current_ids_zero_based: Option<Vec<usize>> = None;
    let mut group_seqs: Vec<Vec<u8>> = Vec::new();
    let mut group_run_paths: Vec<PathBuf> = Vec::new();
    let mut pending_group_tasks: Vec<GroupTask> = Vec::with_capacity(GROUP_WORK_BATCH_GROUPS);
    let mut pending_group_resident_bytes: usize = 0;
    let mut group_bytes: usize = 0;
    let mut submitted_groups: usize = 0;
    let group_batch_budget = memory_budget.group_batch_bytes();
    let group_spill_threshold = GROUP_SORT_SPILL_BYTES.min((group_batch_budget / 2).max(1));
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
                        let task = GroupTask {
                            cid: submitted_groups,
                            dataset_ids_zero_based: ids,
                            seqs: std::mem::take(&mut group_seqs),
                            run_paths: std::mem::take(&mut group_run_paths),
                        };
                        pending_group_resident_bytes = pending_group_resident_bytes
                            .saturating_add(group_task_resident_bytes(&task));
                        pending_group_tasks.push(task);
                        submitted_groups += 1;
                        if pending_group_tasks.len() >= GROUP_WORK_BATCH_GROUPS
                            || pending_group_resident_bytes >= group_batch_budget
                        {
                            process_and_commit_group_batch(
                                &mut pending_group_tasks,
                                &mut pending_group_resident_bytes,
                                &group_spill_dir,
                                &pool,
                                &mut writer,
                                &mut worker_stats,
                                &mut commit_stats,
                            )?;
                        }
                    }
                }

                current_ids_zero_based = Some(map_color_ids_to_zero_based(
                    &record.color_ids,
                    nb_files as usize,
                )?);
                current_color_ids = Some(Arc::clone(&record.color_ids));
                group_bytes = 0;
            }

            group_bytes += record.seq.len();
            group_seqs.push(record.seq);
            if group_bytes >= group_spill_threshold {
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
            let task = GroupTask {
                cid: submitted_groups,
                dataset_ids_zero_based: ids,
                seqs: std::mem::take(&mut group_seqs),
                run_paths: std::mem::take(&mut group_run_paths),
            };
            pending_group_resident_bytes =
                pending_group_resident_bytes.saturating_add(group_task_resident_bytes(&task));
            pending_group_tasks.push(task);
            submitted_groups += 1;
        }
    }

    process_and_commit_group_batch(
        &mut pending_group_tasks,
        &mut pending_group_resident_bytes,
        &group_spill_dir,
        &pool,
        &mut writer,
        &mut worker_stats,
        &mut commit_stats,
    )?;
    let dispatch_timing = dispatch_timer.finish();
    let finalized = writer.finalize()?;
    let _ = fs::remove_dir_all(&group_spill_dir);

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
    log_phase_timing("post_ggcat_stream.total", total_timer.finish());
    println!(
        "[phase-stats] phase=post_ggcat_stream workers={} groups_submitted={} groups_processed={} groups_committed={} commit_tigs_bytes={} commit_sizes_bytes={} worker_max_batch_wall_s={:.3} worker_max_batch_cpu_s={}",
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

    Ok(finalized)
}

fn write_positions_from_spill(
    position_spill_path: &Path,
    entries: u64,
    filepath: String,
) -> Result<()> {
    let mut pos_file = BufWriter::with_capacity(IO_BUFFER_CAPACITY, File::create(filepath)?);
    pos_file.write_all(POSITIONS_MAGIC)?;
    write_varint_u64_to_writer(entries, &mut pos_file)?;
    let mut spill_reader =
        BufReader::with_capacity(IO_BUFFER_CAPACITY, File::open(position_spill_path)?);
    io::copy(&mut spill_reader, &mut pos_file)?;
    pos_file.flush()?;
    fs::remove_file(position_spill_path)?;
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

fn read_optional_varint_u64(reader: &mut impl Read) -> Result<Option<u64>> {
    let mut value = 0u64;
    let mut shift = 0u32;
    let mut saw_byte = false;
    loop {
        let mut byte = [0u8; 1];
        match reader.read_exact(&mut byte) {
            Ok(()) => saw_byte = true,
            Err(err) if err.kind() == io::ErrorKind::UnexpectedEof && !saw_byte => {
                return Ok(None)
            }
            Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "truncated varint in dataset-to-CID spill record",
                ))
            }
            Err(err) => return Err(err),
        }
        value |= ((byte[0] & 0x7f) as u64) << shift;
        if byte[0] & 0x80 == 0 {
            return Ok(Some(value));
        }
        shift += 7;
        if shift >= 64 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "CID delta varint overflow in dataset-to-CID spill record",
            ));
        }
    }
}

fn read_required_varint_u64(reader: &mut impl Read) -> Result<u64> {
    read_optional_varint_u64(reader)?.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "missing varint in dataset-to-CID spill record",
        )
    })
}

#[derive(Debug, Default)]
struct DatasetCidBlocks {
    last_cid: u64,
    seen: bool,
    encoded_deltas: Vec<u8>,
    blocks: Vec<(u64, u32)>,
}

fn flush_cid_transpose_block(
    dataset: &mut DatasetCidBlocks,
    writer: &mut BufWriter<File>,
    spool_offset: &mut u64,
) -> Result<()> {
    if dataset.encoded_deltas.is_empty() {
        return Ok(());
    }
    let block_len = u32::try_from(dataset.encoded_deltas.len()).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "dataset-to-CID transpose block is larger than u32",
        )
    })?;
    dataset.blocks.push((*spool_offset, block_len));
    writer.write_all(&dataset.encoded_deltas)?;
    *spool_offset = spool_offset
        .checked_add(block_len as u64)
        .ok_or_else(|| io::Error::other("dataset-to-CID transpose offset overflow"))?;
    dataset.encoded_deltas.clear();
    Ok(())
}

fn transpose_cid_partition(
    path: &Path,
    partition_index: usize,
    spill_dir: &Path,
    datasets_per_partition: usize,
    dataset_count: usize,
    payload_writer: &mut DatasetPayloadWriter,
) -> Result<(u64, u64)> {
    let dataset_start = partition_index.saturating_mul(datasets_per_partition);
    let dataset_end = dataset_start
        .saturating_add(datasets_per_partition)
        .min(dataset_count);
    let mut datasets = (dataset_start..dataset_end)
        .map(|_| DatasetCidBlocks::default())
        .collect::<Vec<_>>();
    let spool_path = spill_dir.join(format!("partition_{partition_index:04}.transpose.bin"));
    let mut spool_writer =
        BufWriter::with_capacity(IO_BUFFER_CAPACITY, File::create(&spool_path)?);
    let mut spool_offset = 0u64;
    let mut records = 0u64;

    let mut reader = BufReader::with_capacity(IO_BUFFER_CAPACITY, File::open(path)?);
    let mut cid = 0u64;
    while let Some(cid_delta) = read_optional_varint_u64(&mut reader)? {
        cid = cid.checked_add(cid_delta).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "CID overflow in transpose spill")
        })?;
        let count = usize::try_from(read_required_varint_u64(&mut reader)?).map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, "dataset count cannot fit usize")
        })?;
        if count == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "empty CID group in dataset transpose spill",
            ));
        }
        let mut dataset = dataset_start;
        for _ in 0..count {
            let dataset_delta = usize::try_from(read_required_varint_u64(&mut reader)?).map_err(
                |_| io::Error::new(io::ErrorKind::InvalidData, "dataset delta cannot fit usize"),
            )?;
            dataset = dataset.checked_add(dataset_delta).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "dataset index overflow")
            })?;
            if dataset < dataset_start || dataset >= dataset_end {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "dataset {} is outside partition {} range [{}, {})",
                        dataset, partition_index, dataset_start, dataset_end
                    ),
                ));
            }
            let state = &mut datasets[dataset - dataset_start];
            if state.seen && cid < state.last_cid {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "CIDs for dataset {} are not monotonic: {} before {}",
                        dataset, state.last_cid, cid
                    ),
                ));
            }
            let delta = if state.seen { cid - state.last_cid } else { cid };
            write_varint_u64(delta, &mut state.encoded_deltas);
            state.last_cid = cid;
            state.seen = true;
            records += 1;
            if state.encoded_deltas.len() >= CID_TRANSPOSE_BLOCK_BYTES {
                flush_cid_transpose_block(state, &mut spool_writer, &mut spool_offset)?;
            }
        }
    }

    for state in &mut datasets {
        flush_cid_transpose_block(state, &mut spool_writer, &mut spool_offset)?;
    }
    spool_writer.flush()?;
    drop(spool_writer);
    fs::remove_file(path)?;

    let mut spool_reader = File::open(&spool_path)?;
    let mut block_buffer = Vec::new();
    for (local_dataset, state) in datasets.into_iter().enumerate() {
        if !state.seen {
            continue;
        }
        let dataset = dataset_start + local_dataset;
        for (offset, len) in state.blocks {
            spool_reader.seek(SeekFrom::Start(offset))?;
            block_buffer.resize(len as usize, 0);
            spool_reader.read_exact(&mut block_buffer)?;
            payload_writer.append_encoded_deltas(dataset, &block_buffer)?;
        }
    }
    drop(spool_reader);
    fs::remove_file(&spool_path)?;
    Ok((records, spool_offset))
}

struct DatasetPayloadWriter {
    output: BufWriter<File>,
    offsets: Vec<usize>,
    dataset_count: usize,
    total_size: usize,
    current_dataset: Option<usize>,
    payload_path: PathBuf,
    payload_encoder: Option<Encoder<'static, File>>,
    payload_bytes: u64,
}

impl DatasetPayloadWriter {
    fn new(path: &str, spill_dir: &Path, dataset_count: usize) -> Result<Self> {
        let mut output = BufWriter::with_capacity(IO_BUFFER_CAPACITY, File::create(path)?);
        output.write_all(ID_TO_CID_MAGIC)?;
        Ok(Self {
            output,
            offsets: Vec::with_capacity(dataset_count),
            dataset_count,
            total_size: ID_TO_CID_MAGIC.len(),
            current_dataset: None,
            payload_path: spill_dir.join("dataset_payload.zst.tmp"),
            payload_encoder: None,
            payload_bytes: 0,
        })
    }

    fn start_dataset(&mut self, dataset: usize) -> Result<()> {
        if self.current_dataset.is_some() {
            return Err(io::Error::other(
                "cannot start dataset payload before finishing the previous payload",
            ));
        }
        let payload_file = File::options()
            .read(true)
            .write(true)
            .create(true)
            .truncate(true)
            .open(&self.payload_path)?;
        let encoder = Encoder::new(payload_file, 1)?;
        self.current_dataset = Some(dataset);
        self.payload_encoder = Some(encoder);
        Ok(())
    }

    fn finish_dataset(&mut self) -> Result<()> {
        let dataset = self
            .current_dataset
            .take()
            .ok_or_else(|| io::Error::other("no active dataset payload"))?;
        let encoder = self
            .payload_encoder
            .take()
            .ok_or_else(|| io::Error::other("missing dataset payload encoder"))?;
        let mut payload_file = encoder.finish()?;
        let payload_len = payload_file.seek(SeekFrom::End(0))?;
        payload_file.seek(SeekFrom::Start(0))?;
        self.offsets.push(self.total_size);
        self.output.write_all(&payload_len.to_le_bytes())?;
        io::copy(&mut payload_file, &mut self.output)?;
        self.total_size = self
            .total_size
            .checked_add(8)
            .and_then(|size| size.checked_add(payload_len as usize))
            .ok_or_else(|| io::Error::other("id-to-CID output size overflow"))?;
        self.payload_bytes = self.payload_bytes.saturating_add(payload_len);
        drop(payload_file);
        debug_assert_eq!(self.offsets.len(), dataset + 1);
        Ok(())
    }

    fn advance_to(&mut self, dataset: usize) -> Result<()> {
        if let Some(current) = self.current_dataset {
            if dataset < current {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "dataset spill order regressed from {} to {}",
                        current, dataset
                    ),
                ));
            }
            if dataset == current {
                return Ok(());
            }
            self.finish_dataset()?;
        }
        while self.offsets.len() < dataset {
            let missing = self.offsets.len();
            self.start_dataset(missing)?;
            self.finish_dataset()?;
        }
        self.start_dataset(dataset)
    }

    fn append_encoded_deltas(&mut self, dataset: usize, encoded_deltas: &[u8]) -> Result<()> {
        if dataset >= self.dataset_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "dataset {} outside [0, {}) in transposed dataset-to-CID spill",
                    dataset, self.dataset_count
                ),
            ));
        }
        self.advance_to(dataset)?;
        self.payload_encoder
            .as_mut()
            .ok_or_else(|| io::Error::other("missing dataset payload encoder"))?
            .write_all(encoded_deltas)
    }

    fn finish(mut self, dataset_count: usize) -> Result<(Vec<usize>, u64)> {
        if self.current_dataset.is_some() {
            self.finish_dataset()?;
        }
        while self.offsets.len() < dataset_count {
            let dataset = self.offsets.len();
            self.start_dataset(dataset)?;
            self.finish_dataset()?;
        }
        if self.offsets.len() != dataset_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "dataset-to-CID spill contains an out-of-range dataset",
            ));
        }
        self.output.write_all(&0u64.to_le_bytes())?;
        self.output.flush()?;
        if self.payload_path.exists() {
            fs::remove_file(&self.payload_path)?;
        }
        Ok((self.offsets, self.payload_bytes))
    }
}

fn write_id_to_color_id_from_partitions(
    cid_file_path: String,
    spill: DatasetCidSpillFiles,
    _memory_budget: CompressionMemoryBudget,
) -> std::io::Result<Vec<usize>> {
    let total_timer = PhaseTimer::start();
    let mut payload_writer =
        DatasetPayloadWriter::new(&cid_file_path, &spill.directory, spill.dataset_count)?;
    let mut transposed_records = 0u64;
    let mut transpose_bytes = 0u64;
    for (partition_index, path) in spill.paths.iter().enumerate() {
        let (records, bytes) = transpose_cid_partition(
            path,
            partition_index,
            &spill.directory,
            spill.datasets_per_partition,
            spill.dataset_count,
            &mut payload_writer,
        )?;
        transposed_records = transposed_records.saturating_add(records);
        transpose_bytes = transpose_bytes.saturating_add(bytes);
    }
    let (offsets, total_payload_bytes) = payload_writer.finish(spill.dataset_count)?;
    let _ = fs::remove_dir(&spill.directory);
    log_phase_timing("post_ggcat.id_to_cid.total", total_timer.finish());
    println!(
        "[phase-stats] phase=post_ggcat.id_to_cid datasets={} partitions={} datasets_per_partition={} transposed_records={} transpose_bytes={} payload_bytes={}",
        spill.dataset_count,
        spill.paths.len(),
        spill.datasets_per_partition,
        transposed_records,
        transpose_bytes,
        total_payload_bytes
    );
    Ok(offsets)
}

pub(crate) fn sort_by_bucket_streaming(
    output_dir: &String,
    nb_files: u32,
    worker_threads: usize,
    memory_gb: usize,
    record_rx: mpsc::Receiver<SimplitigBatch>,
) -> Vec<usize> {
    let total_timer = PhaseTimer::start();
    let memory_budget = CompressionMemoryBudget::from_gb(memory_gb);
    println!("Starting writing compressed sequences (streaming).");
    let write_stream_timer = PhaseTimer::start();
    let triple = match write_compressed_from_stream(
        output_dir.clone() + "tigs_kloe.fa",
        output_dir,
        nb_files,
        worker_threads,
        memory_budget,
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
    let position_timer = PhaseTimer::start();
    println!("Starting to write positions");
    if let Err(e) = write_positions_from_spill(
        &triple.position_path,
        triple.position_entries,
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
    let write_id_cid = match write_id_to_color_id_from_partitions(
        output_dir.clone() + "id_to_color_id.txt.zst",
        triple.cid_spill,
        memory_budget,
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
            COLOR_RUN_BUFFER_BYTES,
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

fn reduce_sorted_color_chunks(mut chunks: Vec<PathBuf>, chunk_dir: &Path) -> Result<Vec<PathBuf>> {
    let mut merge_pass = 0usize;
    while chunks.len() > COLOR_MERGE_FAN_IN {
        let mut next_chunks = Vec::with_capacity(chunks.len().div_ceil(COLOR_MERGE_FAN_IN));
        for (group_index, group) in chunks.chunks(COLOR_MERGE_FAN_IN).enumerate() {
            if group.len() == 1 {
                next_chunks.push(group[0].clone());
                continue;
            }
            let path = chunk_dir.join(format!("color_merge_{merge_pass:03}_{group_index:06}.bin"));
            let mut writer = BufWriter::with_capacity(IO_BUFFER_CAPACITY, File::create(&path)?);
            stream_sorted_records_from_chunks(group, |record| {
                write_sorted_color_record(&mut writer, &record)
            })?;
            writer.flush()?;
            for old_path in group {
                fs::remove_file(old_path)?;
            }
            next_chunks.push(path);
        }
        chunks = next_chunks;
        merge_pass += 1;
    }
    Ok(chunks)
}

struct StructuredGgcatCapture {
    chunk_files: Vec<PathBuf>,
    chunk_flush_metrics: ChunkFlushMetrics,
    input_sequences_count: usize,
    emitted_segments_count: usize,
    initial_chunk_count: usize,
}

fn capture_structured_ggcat_output(
    receiver: mpsc::Receiver<Vec<u8>>,
    chunk_dir: PathBuf,
    memory_budget: CompressionMemoryBudget,
    k: usize,
    worker_threads: usize,
) -> Result<StructuredGgcatCapture> {
    // This pool must remain independent from GGCAT's pool: the producer may be
    // blocked waiting for this consumer while a chunk is being sorted.
    let sort_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(worker_threads.max(1))
        .thread_name(|index| format!("kloe-color-sort-{index}"))
        .build()
        .map_err(|err| io::Error::other(format!("create color-sort worker pool: {err}")))?;
    let mut records = Vec::new();
    let mut chunk_files = Vec::new();
    let mut chunk_bytes = 0usize;
    let mut chunk_flush_metrics = ChunkFlushMetrics::default();
    let mut input_sequences_count = 0usize;
    let mut emitted_segments_count = 0usize;

    for block in receiver {
        let (sequences, segments) = visit_structured_ggcat_block(&block, k, |record| {
            chunk_bytes = chunk_bytes
                .saturating_add(std::mem::size_of::<ColorIndexType>() + record.seq.len());
            records.push(record);
            if chunk_bytes >= memory_budget.color_chunk_bytes() {
                flush_sorted_chunk_timed(
                    &mut records,
                    &mut chunk_files,
                    &chunk_dir,
                    &sort_pool,
                    &mut chunk_flush_metrics,
                )?;
                chunk_bytes = 0;
            }
            Ok(())
        })?;
        input_sequences_count = input_sequences_count.saturating_add(sequences);
        emitted_segments_count = emitted_segments_count.saturating_add(segments);
    }

    flush_sorted_chunk_timed(
        &mut records,
        &mut chunk_files,
        &chunk_dir,
        &sort_pool,
        &mut chunk_flush_metrics,
    )?;
    let initial_chunk_count = chunk_files.len();
    let chunk_files = reduce_sorted_color_chunks(chunk_files, &chunk_dir)?;
    Ok(StructuredGgcatCapture {
        chunk_files,
        chunk_flush_metrics,
        input_sequences_count,
        emitted_segments_count,
        initial_chunk_count,
    })
}

fn produce_ggcat_records(
    input_streams: Vec<GeneralSequenceBlockData>,
    source_dataset_ids: SourceDatasetMap,
    dataset_count: usize,
    threads: usize,
    k: usize,
    m: usize,
    use_unitigs: bool,
    use_matchtigs: bool,
    use_eulertigs: bool,
    ggcat_cfg: GgcatCompressionConfig,
    sender: mpsc::SyncSender<SimplitigBatch>,
) -> Result<()> {
    let source_count = source_dataset_ids.len();
    let memory_budget = CompressionMemoryBudget::from_gb(ggcat_cfg.memory_gb);
    if source_count == 0 || input_streams.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "ggcat source list is empty",
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

    let structured_output_key =
        ggcat_temp_root.join(format!("kloe-color-stream-{pid}-{nanos}.tmp"));
    let chunk_dir = ggcat_temp_root.join(format!("kloe-color-chunks-{pid}-{nanos}"));
    fs::create_dir_all(&chunk_dir)?;

    let instance = GGCATInstance::create(GGCATConfig {
        temp_dir: Some(ggcat_temp_root.clone()),
        memory: memory_budget.ggcat_memory_gb(),
        prefer_memory: false,
        total_threads_count: threads.max(1),
        intermediate_compression_level: None,
        stats_file: None,
        messages_callback: None,
    })
    .map_err(|e| to_io_err("create ggcat instance", e))?;

    println!(
        "Running embedded ggcat structured-output path (k={}, m={}, threads={}, sources={}, datasets={}, total_memory={}GB, ggcat_memory={:.2}GB, disk_backed=true)",
        k,
        m,
        threads.max(1),
        source_count,
        dataset_count,
        ggcat_cfg.memory_gb.max(1),
        memory_budget.ggcat_memory_gb()
    );

    let (structured_tx, structured_rx) = mpsc::sync_channel::<Vec<u8>>(2);
    if !register_channel_output(structured_output_key.clone(), structured_tx) {
        return Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            format!(
                "a GGCAT structured output is already registered for {}",
                structured_output_key.display()
            ),
        ));
    }
    let capture_chunk_dir = chunk_dir.clone();
    let capture_thread = thread::spawn(move || {
        capture_structured_ggcat_output(
            structured_rx,
            capture_chunk_dir,
            memory_budget,
            k,
            threads,
        )
    });

    let ggcat_build_timer = PhaseTimer::start();
    let build_result = instance
        .build_graph(
            input_streams,
            structured_output_key.clone(),
            None,
            Some(source_count),
            k,
            threads.max(1),
            false,
            Some(m),
            true,
            1,
            ggcat_extra_elaboration(use_unitigs, use_matchtigs, use_eulertigs),
            None,
        );
    unregister_channel_output(&structured_output_key);
    let capture_result = capture_thread
        .join()
        .map_err(|_| io::Error::other("GGCAT structured-output capture thread panicked"))?;
    let records_output = build_result.map_err(|e| to_io_err("ggcat build_graph", e))?;
    let captured = capture_result?;
    log_phase_timing("ggcat.build_graph_api_call", ggcat_build_timer.finish());

    let colormap_file = GGCATInstance::get_colormap_file(&records_output);
    let colors_deserializer =
        ColorsDeserializer::<DefaultColorsSerializer>::new(&colormap_file, false)
            .map_err(|e| to_io_err("open ggcat colormap", e))?;

    let colors_count = colors_deserializer.colors_count();
    if colors_count != source_count {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "ggcat returned {colors_count} colors but KLOE registered {source_count} sources"
            ),
        ));
    }

    let keep_ggcat = std::env::var_os("KLOE_KEEP_GGCAT").is_some();
    let process_result = (|| -> Result<()> {
        let producer_total_timer = PhaseTimer::start();
        let StructuredGgcatCapture {
            chunk_files,
            chunk_flush_metrics,
            input_sequences_count,
            emitted_segments_count,
            initial_chunk_count,
        } = captured;

        let resolve_and_emit_timer = PhaseTimer::start();
        let mut batch: SimplitigBatch = Vec::with_capacity(COLOR_RECORD_BATCH_SIZE);
        let mut batch_bytes = 0usize;
        let mut record_window = Vec::<SortedColorRecord>::new();
        let mut record_window_bytes = 0usize;
        let mut record_window_subsets = 0usize;
        let mut window_last_subset = None;
        let mut global_last_subset = None;
        let mut unique_subsets_count = 0usize;
        let mut merged_records_count = 0usize;
        let mut emitted_batches_count = 0usize;
        let mut emitted_records_count = 0usize;
        let window_sequence_limit = memory_budget.subset_window_sequence_bytes();
        let window_subset_limit = memory_budget.subset_window_count(dataset_count);

        stream_sorted_records_from_chunks(&chunk_files, |record| {
            merged_records_count += 1;
            if global_last_subset != Some(record.subset) {
                unique_subsets_count += 1;
                global_last_subset = Some(record.subset);
            }
            let is_new_window_subset = window_last_subset != Some(record.subset);
            if !record_window.is_empty()
                && (record_window_bytes >= window_sequence_limit
                    || (is_new_window_subset && record_window_subsets >= window_subset_limit))
            {
                emit_resolved_record_window(
                    &mut record_window,
                    &instance,
                    &colormap_file,
                    None,
                    &source_dataset_ids,
                    dataset_count,
                    &sender,
                    &mut batch,
                    &mut batch_bytes,
                    window_sequence_limit,
                    &mut emitted_batches_count,
                    &mut emitted_records_count,
                )?;
                record_window_bytes = 0;
                record_window_subsets = 0;
                window_last_subset = None;
            }
            if window_last_subset != Some(record.subset) {
                record_window_subsets += 1;
                window_last_subset = Some(record.subset);
            }
            record_window_bytes = record_window_bytes.saturating_add(record.seq.len());
            record_window.push(record);
            Ok(())
        })?;

        emit_resolved_record_window(
            &mut record_window,
            &instance,
            &colormap_file,
            None,
            &source_dataset_ids,
            dataset_count,
            &sender,
            &mut batch,
            &mut batch_bytes,
            window_sequence_limit,
            &mut emitted_batches_count,
            &mut emitted_records_count,
        )?;

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

        let resolve_and_emit_timing = resolve_and_emit_timer.finish();
        log_phase_timing(
            "ggcat.structured_capture_chunk_flush_total",
            PhaseTiming {
                wall_sec: chunk_flush_metrics.timing.wall_sec,
                cpu_sec: chunk_flush_metrics.timing.cpu_sec,
            },
        );
        log_phase_timing(
            "producer.resolve_subsets_and_emit_records",
            resolve_and_emit_timing,
        );
        log_phase_timing("producer.total_post_ggcat", producer_total_timer.finish());
        println!(
            "[phase-stats] phase=producer input_sequences={} emitted_segments={} subsets_unique={} chunks={} chunk_flush_calls={} chunk_flush_records={} merged_records={} emitted_batches={} emitted_records={}",
            input_sequences_count,
            emitted_segments_count,
            unique_subsets_count,
            initial_chunk_count,
            chunk_flush_metrics.calls,
            chunk_flush_metrics.records,
            merged_records_count,
            emitted_batches_count,
            emitted_records_count
        );

        drop(sender);
        Ok(())
    })();

    if keep_ggcat {
        eprintln!(
            "DEBUG_KEEP_GGCAT structured_output=no_file colormap={} chunks_dir={}",
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

#[allow(clippy::too_many_arguments)]
pub(crate) fn compress_ggcat_sources(
    output_dir: &str,
    filenames: Vec<String>,
    input_streams: Vec<GeneralSequenceBlockData>,
    source_dataset_ids: SourceDatasetMap,
    threads: usize,
    k: usize,
    m: usize,
    use_unitigs: bool,
    use_matchtigs: bool,
    use_eulertigs: bool,
    ggcat_cfg: GgcatCompressionConfig,
) -> Result<()> {
    let dataset_count = filenames.len();
    if dataset_count == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "input file list is empty",
        ));
    }

    println!("Compression backend: embedded ggcat structured output");

    let (record_tx, record_rx) = mpsc::sync_channel::<SimplitigBatch>(2);
    let producer_cfg = ggcat_cfg.clone();

    let producer = thread::spawn(move || {
        produce_ggcat_records(
            input_streams,
            source_dataset_ids,
            dataset_count,
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
    let id_cid_line_sizes = sort_by_bucket_streaming(
        &output_dir.to_string(),
        dataset_count as u32,
        threads,
        ggcat_cfg.memory_gb,
        record_rx,
    );
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
            "Warning: --skip-sort is ignored; KLOE groups GGCAT structured records by colors."
        );
    }

    let filenames = read_input_fof_filenames(input_fof)?;
    let input_streams = filenames
        .iter()
        .map(|file| {
            let path = PathBuf::from(file);
            let resolved = fs::canonicalize(&path).unwrap_or(path);
            GeneralSequenceBlockData::FASTA((resolved, None))
        })
        .collect::<Vec<_>>();
    let source_dataset_ids = SourceDatasetMap::from_singletons(filenames.len())?;

    compress_ggcat_sources(
        output_dir,
        filenames,
        input_streams,
        source_dataset_ids,
        threads,
        k,
        m,
        use_unitigs,
        use_matchtigs,
        use_eulertigs,
        ggcat_cfg,
    )
}

#[cfg(test)]
mod structured_output_tests {
    use super::*;
    use ggcat_colors::storage::run_length::RunLengthColorsSerializer;
    use ggcat_colors::storage::serializer::ColorsSerializer;

    #[test]
    fn source_color_union_is_linear_sorted_and_deduplicated() {
        let mut ids = vec![1, 3, 7, 9];
        merge_sorted_dataset_ids(&mut ids, &[2, 3, 8, 9, 10]);
        assert_eq!(ids, vec![1, 2, 3, 7, 8, 9, 10]);

        merge_sorted_dataset_ids_with_offset(&mut ids, &[1, 2, 5], 10).unwrap();
        assert_eq!(ids, vec![1, 2, 3, 7, 8, 9, 10, 11, 12, 15]);
    }

    #[test]
    fn cid_dataset_sidecar_roundtrips_and_applies_archive_offset() {
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join(CID_TO_DATASET_FILE);
        let mut writer = CidDatasetSidecarWriter::create(&path).unwrap();
        writer.append_one_based(&[1, 3, 9]).unwrap();
        writer.append_zero_based(&[1, 4]).unwrap();
        writer.finish().unwrap();

        let sidecar = Arc::new(CidDatasetSidecar::open(&path).unwrap());
        assert_eq!(sidecar.len(), 2);
        let mut map = SourceDatasetMap::new();
        map.add_disk_storage(sidecar, 10).unwrap();

        let mut ids = vec![2, 11, 20];
        map.merge_source_into(0, &mut ids).unwrap();
        assert_eq!(ids, vec![2, 11, 13, 19, 20]);
        let mut ids = Vec::new();
        map.merge_source_into(1, &mut ids).unwrap();
        assert_eq!(ids, vec![12, 15]);
    }

    #[test]
    fn structured_color_runs_are_split_without_fasta_headers() {
        let sequence = b"AACCGG";
        let mut colors = Vec::new();
        write_varint_u64(2, &mut colors);
        write_varint_u64(5, &mut colors);
        write_varint_u64(2, &mut colors);
        write_varint_u64(7, &mut colors);
        write_varint_u64(2, &mut colors);

        let mut block = Vec::new();
        block.extend_from_slice(&11u64.to_le_bytes());
        block.extend_from_slice(&(sequence.len() as u64).to_le_bytes());
        block.extend_from_slice(&(colors.len() as u64).to_le_bytes());
        block.extend_from_slice(&0u64.to_le_bytes());
        block.extend_from_slice(sequence);
        block.extend_from_slice(&colors);

        let mut records = Vec::new();
        let counts = visit_structured_ggcat_block(&block, 3, |record| {
            records.push(record);
            Ok(())
        })
        .unwrap();

        assert_eq!(counts, (1, 2));
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].subset, 5);
        assert_eq!(records[0].seq, b"AACC");
        assert_eq!(records[1].subset, 7);
        assert_eq!(records[1].seq, b"CCGG");
    }

    #[test]
    fn parallel_checkpoint_pipeline_preserves_assigned_color_ids() {
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("colors.dat");
        let serializer = Arc::new(
            ColorsSerializer::<RunLengthColorsSerializer>::new(
                &path,
                &[],
                Some(128),
                4,
                false,
            )
            .unwrap(),
        );

        let handles = (0..8)
            .map(|thread_id| {
                let serializer = Arc::clone(&serializer);
                thread::spawn(move || {
                    let mut assigned = Vec::new();
                    for index in 0..4_000u32 {
                        let color = (thread_id * 17 + index) % 128;
                        let subset_id = serializer.serialize_colors(&[color]);
                        assigned.push((subset_id, color));
                    }
                    assigned
                })
            })
            .collect::<Vec<_>>();

        let mut assigned = handles
            .into_iter()
            .flat_map(|handle| handle.join().unwrap())
            .collect::<Vec<_>>();
        assigned.sort_unstable_by_key(|(subset_id, _)| *subset_id);
        Arc::try_unwrap(serializer).ok().unwrap().finalize();

        let mut deserializer =
            ColorsDeserializer::<RunLengthColorsSerializer>::new(&path, false).unwrap();
        assert_eq!(deserializer.colors_count(), 128);
        let mut decoded = Vec::new();
        for (subset_id, expected_color) in assigned {
            decoded.clear();
            deserializer.get_color_mappings(subset_id, &mut decoded);
            assert_eq!(decoded, [expected_color]);
        }
    }
}
