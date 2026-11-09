use crate::buckets::writers::{finalize_bucket_file, initialize_bucket_file, THREADS_BUSY_WRITING};
use crate::buckets::{CheckpointData, LockFreeBucket};
use crate::memory_data_size::MemoryDataSize;
use crate::memory_fs::file::flush::GlobalFlush;
use crate::memory_fs::file::internal::MemoryFileMode;
use crate::memory_fs::file::reader::FileRangeReference;
use crate::memory_fs::file::writer::FileWriter;
use crate::utils::memory_size_to_log2;
use crate::DEFAULT_BINCODE_CONFIG;
use bincode::Encode;
use lz4::{BlockMode, BlockSize, ContentChecksum};
use mt_debug_counters::counter::AtomicCounterGuardSum;
use parking_lot::Mutex;
use std::collections::VecDeque;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::{Arc, LazyLock, Weak};

use super::BucketHeader;

pub const COMPRESSED_BUCKET_MAGIC: &[u8; 16] = b"CPLZ4_INTR_BKT_M";

#[derive(Clone)]
pub struct CompressedCheckpointSize(u8);
impl CompressedCheckpointSize {
    pub const fn new_from_size(size: MemoryDataSize) -> Self {
        Self(memory_size_to_log2(size))
    }
    pub const fn new_from_log2(val: u8) -> Self {
        Self(val)
    }
}

fn create_lz4_stream<W: Write>(writer: W, level: CompressionLevelInfo) -> lz4::Encoder<W> {
    let (queue_occupation, queue_size) = GlobalFlush::global_queue_occupation();

    let level = if queue_size < 2 * queue_occupation {
        level.slow_disk
    } else {
        level.fast_disk
    };

    lz4::EncoderBuilder::new()
        .level(level)
        .checksum(ContentChecksum::NoChecksum)
        .block_mode(BlockMode::Linked)
        .block_size(BlockSize::Max64KB)
        .build(writer)
        .unwrap()
}

struct CompressedBinaryWriterInternal {
    writer: Option<lz4::Encoder<FileWriter>>,
    idle_file: Option<FileWriter>,
    checkpoint_max_size: u64,
    checkpoints: Vec<CheckpointData>,
    current_chunk_size: u64,
    level: CompressionLevelInfo,
    data_format_info: Vec<u8>,

    checkpoint_data: Option<Vec<u8>>,
}

pub struct CompressedBinaryWriter {
    inner: Arc<Mutex<CompressedBinaryWriterInternal>>,
    path: PathBuf,
}
unsafe impl Send for CompressedBinaryWriter {}

// A large collection can have 8192 buckets. Keeping an LZ4 context for
// every bucket costs memory even before any sequence is downloaded. Retain
// only a bounded working set; eviction finishes an ordinary indexed frame.
static COMPRESSORS: LazyLock<Mutex<VecDeque<Weak<Mutex<CompressedBinaryWriterInternal>>>>> =
    LazyLock::new(|| Mutex::new(VecDeque::new()));

fn compressor_cache_capacity() -> usize {
    // Reserve at most 1/16 of the configured temporary-memory budget for
    // retained LZ4 contexts (conservative 512 KiB per context). Large-memory
    // builds can keep all 8192 contexts and avoid unnecessary frame restarts.
    (crate::memory_fs::allocator::CHUNKS_ALLOCATOR
        .get_total_memory()
        .as_bytes()
        / 16
        / (512 * 1024))
        .clamp(16, 8192)
}

impl CompressedBinaryWriterInternal {
    fn file(&self) -> &FileWriter {
        match &self.writer {
            Some(writer) => writer.writer(),
            None => self.idle_file.as_ref().unwrap(),
        }
    }

    fn finish_block(&mut self) {
        if let Some(writer) = self.writer.take() {
            let (file, result) = writer.finish();
            result.unwrap();
            file.flush_async();
            self.idle_file = Some(file);
        }
        self.current_chunk_size = 0;
    }

    fn activate(owner: &Arc<Mutex<Self>>, inner: &mut Self) {
        if inner.writer.is_some() {
            return;
        }
        let mut cache = COMPRESSORS.lock();
        // Never wait on another bucket lock while holding our own. Busy
        // writers may temporarily add at most one context per active caller.
        let mut remaining = cache.len();
        while cache.len() >= compressor_cache_capacity() && remaining > 0 {
            remaining -= 1;
            if let Some(candidate) = cache.pop_front().unwrap().upgrade() {
                if let Some(mut old) = candidate.try_lock() {
                    // Metadata can describe an entire checkpoint (e.g. its
                    // sequence count); never split such a logical block.
                    if old.checkpoint_data.is_none() && old.checkpoint_max_size < (1 << 62) {
                        old.finish_block();
                    }
                } else {
                    cache.push_back(Arc::downgrade(&candidate));
                }
            }
        }
        let file = inner.idle_file.take().unwrap();
        inner.checkpoints.push(CheckpointData {
            offset: file.len(),
            data: inner.checkpoint_data.clone(),
        });
        inner.writer = Some(create_lz4_stream(file, inner.level));
        if inner.checkpoint_data.is_none() && inner.checkpoint_max_size < (1 << 62) {
            cache.push_back(Arc::downgrade(owner));
        }
    }

    fn create_new_block(&mut self, passtrough_range: Option<FileRangeReference>) {
        self.finish_block();
        if let Some(range) = passtrough_range {
            let file = self.idle_file.as_ref().unwrap();
            self.checkpoints.push(CheckpointData {
                offset: file.len(),
                data: self.checkpoint_data.clone(),
            });
            unsafe {
                range.copy_to_unsync(file);
            }
        }
    }
}

impl CompressedBinaryWriter {
    pub const CHECKPOINT_SIZE_UNLIMITED: CompressedCheckpointSize =
        CompressedCheckpointSize::new_from_log2(62);
}

#[derive(Copy, Clone)]
pub struct CompressionLevelInfo {
    pub fast_disk: u32,
    pub slow_disk: u32,
}

impl LockFreeBucket for CompressedBinaryWriter {
    type InitData = (
        MemoryFileMode,
        CompressedCheckpointSize,
        CompressionLevelInfo,
    );

    fn new_serialized_data_format(
        path_prefix: &Path,
        (file_mode, checkpoint_max_size, compression_level): &(
            MemoryFileMode,
            CompressedCheckpointSize,
            CompressionLevelInfo,
        ),
        index: usize,
        data_format_info: &[u8],
    ) -> Self {
        assert!(
            data_format_info.len() <= BucketHeader::MAX_DATA_FORMAT_INFO_SIZE,
            "Serialized data format info is too big, this is a bug"
        );

        let path = path_prefix.parent().unwrap().join(format!(
            "{}.{}",
            path_prefix.file_name().unwrap().to_str().unwrap(),
            index
        ));

        let mut file = FileWriter::create(&path, *file_mode);

        initialize_bucket_file(&mut file);
        file.flush_async();

        Self {
            inner: Arc::new(Mutex::new(CompressedBinaryWriterInternal {
                writer: None,
                idle_file: Some(file),
                checkpoint_max_size: (1 << checkpoint_max_size.0),
                checkpoints: Vec::new(),
                current_chunk_size: 0,
                level: *compression_level,
                data_format_info: data_format_info.to_vec(),
                checkpoint_data: None,
            })),
            path,
        }
    }

    fn set_checkpoint_data<T: Encode>(
        &self,
        data: Option<&T>,
        passtrough_range: Option<FileRangeReference>,
    ) {
        let mut inner = self.inner.lock();
        inner.checkpoint_data =
            data.map(|data| bincode::encode_to_vec(data, DEFAULT_BINCODE_CONFIG).unwrap());
        // Always create a new block on checkpoint data change
        inner.create_new_block(passtrough_range);
    }

    fn get_bucket_size(&self) -> u64 {
        self.inner.lock().file().len()
    }

    fn write_data(&self, bytes: &[u8]) -> u64 {
        let stat_raii = AtomicCounterGuardSum::new(&THREADS_BUSY_WRITING, 1);
        //
        let mut inner = self.inner.lock();

        CompressedBinaryWriterInternal::activate(&self.inner, &mut inner);
        inner.writer.as_mut().unwrap().write_all(bytes).unwrap();
        inner.current_chunk_size += bytes.len() as u64;
        if inner.current_chunk_size > inner.checkpoint_max_size {
            inner.create_new_block(None);
        }

        drop(stat_raii);
        inner.file().len()
    }

    fn get_path(&self) -> PathBuf {
        self.path.clone()
    }
    fn finalize(self) {
        let mut inner = self.inner.lock();
        // Keep an empty frame for empty buckets, as in the original format.
        if inner.checkpoints.is_empty() {
            CompressedBinaryWriterInternal::activate(&self.inner, &mut inner);
        }
        inner.finish_block();
        let file = inner.idle_file.take().unwrap();
        let checkpoints = std::mem::take(&mut inner.checkpoints);
        finalize_bucket_file(
            file,
            COMPRESSED_BUCKET_MAGIC,
            checkpoints,
            &inner.data_format_info,
        );
    }
}
