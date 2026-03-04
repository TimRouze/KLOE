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
use replace_with::replace_with_or_abort;
use std::io::Write;
use std::path::{Path, PathBuf};

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
    writer: lz4::Encoder<FileWriter>,
    checkpoint_max_size: u64,
    checkpoints: Vec<CheckpointData>,
    current_chunk_size: u64,
    level: CompressionLevelInfo,
    data_format_info: Vec<u8>,

    checkpoint_data: Option<Vec<u8>>,
}

pub struct CompressedBinaryWriter {
    inner: Mutex<CompressedBinaryWriterInternal>,
    path: PathBuf,
}
unsafe impl Send for CompressedBinaryWriter {}

impl CompressedBinaryWriterInternal {
    fn create_new_block(&mut self, passtrough_range: Option<FileRangeReference>) {
        replace_with_or_abort(&mut self.writer, |writer| {
            let (file_buf, res) = writer.finish();
            res.unwrap();

            let checkpoint_pos = if let Some(passtrough_range) = passtrough_range {
                // Add an optional passtrough block
                self.checkpoints.push(CheckpointData {
                    offset: file_buf.len() as u64,
                    data: self.checkpoint_data.clone(),
                });

                unsafe {
                    passtrough_range.copy_to_unsync(&file_buf);
                }

                file_buf.len()
            } else {
                file_buf.len()
            };

            self.checkpoints.push(CheckpointData {
                offset: checkpoint_pos as u64,
                data: self.checkpoint_data.clone(),
            });

            create_lz4_stream(file_buf, self.level)
        });
        self.current_chunk_size = 0;
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

        let first_checkpoint = initialize_bucket_file(&mut file);

        let writer = create_lz4_stream(file, *compression_level);

        Self {
            inner: Mutex::new(CompressedBinaryWriterInternal {
                writer,
                checkpoint_max_size: (1 << checkpoint_max_size.0),
                checkpoints: vec![CheckpointData {
                    offset: first_checkpoint,
                    data: None,
                }],
                current_chunk_size: 0,
                level: *compression_level,
                data_format_info: data_format_info.to_vec(),
                checkpoint_data: None,
            }),
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
        self.inner.lock().writer.writer().len()
    }

    fn write_data(&self, bytes: &[u8]) -> u64 {
        let stat_raii = AtomicCounterGuardSum::new(&THREADS_BUSY_WRITING, 1);
        //
        let mut inner = self.inner.lock();

        inner.writer.write_all(bytes).unwrap();
        inner.current_chunk_size += bytes.len() as u64;
        if inner.current_chunk_size > inner.checkpoint_max_size {
            inner.create_new_block(None);
        }

        drop(stat_raii);
        inner.writer.writer().len()
    }

    fn get_path(&self) -> PathBuf {
        self.path.clone()
    }
    fn finalize(self) {
        let mut inner = self.inner.into_inner();

        inner.writer.flush().unwrap();
        let (file, res) = inner.writer.finish();
        res.unwrap();

        finalize_bucket_file(
            file,
            COMPRESSED_BUCKET_MAGIC,
            inner.checkpoints,
            &inner.data_format_info,
        );
    }
}
