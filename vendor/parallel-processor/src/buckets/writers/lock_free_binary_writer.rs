use crate::buckets::writers::{
    finalize_bucket_file, initialize_bucket_file, BucketHeader, THREADS_BUSY_WRITING,
};
use crate::buckets::{CheckpointData, LockFreeBucket};
use crate::memory_data_size::MemoryDataSize;
use crate::memory_fs::file::internal::MemoryFileMode;
use crate::memory_fs::file::reader::FileRangeReference;
use crate::memory_fs::file::writer::FileWriter;
use crate::utils::memory_size_to_log2;
use crate::DEFAULT_BINCODE_CONFIG;
use bincode::Encode;
use mt_debug_counters::counter::AtomicCounterGuardSum;
use parking_lot::Mutex;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

pub const LOCK_FREE_BUCKET_MAGIC: &[u8; 16] = b"PLAIN_INTR_BKT_M";

#[derive(Clone)]
pub struct LockFreeCheckpointSize(u8);
impl LockFreeCheckpointSize {
    pub const fn new_from_size(size: MemoryDataSize) -> Self {
        Self(memory_size_to_log2(size))
    }
    pub const fn new_from_log2(val: u8) -> Self {
        Self(val)
    }
}

pub struct LockFreeBinaryWriter {
    writer: FileWriter,
    checkpoint_max_size_log2: u8,
    checkpoints: Mutex<Vec<CheckpointData>>,
    checkpoint_data: Mutex<Option<Vec<u8>>>,

    file_size: AtomicU64,
    data_format_info: Vec<u8>,
}
unsafe impl Send for LockFreeBinaryWriter {}

impl LockFreeBinaryWriter {
    pub const CHECKPOINT_SIZE_UNLIMITED: LockFreeCheckpointSize =
        LockFreeCheckpointSize::new_from_log2(62);
}

impl LockFreeBucket for LockFreeBinaryWriter {
    type InitData = (MemoryFileMode, LockFreeCheckpointSize);

    fn new_serialized_data_format(
        path_prefix: &Path,
        (file_mode, checkpoint_max_size): &(MemoryFileMode, LockFreeCheckpointSize),
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

        let mut writer = FileWriter::create(path, *file_mode);

        let first_checkpoint = initialize_bucket_file(&mut writer);

        Self {
            writer,
            checkpoint_max_size_log2: checkpoint_max_size.0,
            checkpoints: Mutex::new(vec![CheckpointData {
                offset: first_checkpoint,
                data: None,
            }]),
            file_size: AtomicU64::new(0),
            data_format_info: data_format_info.to_vec(),
            checkpoint_data: Mutex::new(None),
        }
    }

    fn set_checkpoint_data<T: Encode>(
        &self,
        data: Option<&T>,
        passtrough_range: Option<FileRangeReference>,
    ) {
        let data = data.map(|data| bincode::encode_to_vec(data, DEFAULT_BINCODE_CONFIG).unwrap());
        *self.checkpoint_data.lock() = data.clone();
        // Always create a new block on checkpoint data change

        if let Some(passtrough_range) = passtrough_range {
            let position = self.writer.write_all_parallel(&[], 1);
            self.checkpoints.lock().push(CheckpointData {
                offset: position as u64,
                data: data.clone(),
            });
            unsafe {
                // TODO: Be sure that the file has exclusive access
                passtrough_range.copy_to_unsync(&self.writer);
            }
        }

        let position = self.writer.write_all_parallel(&[], 1);
        self.checkpoints.lock().push(CheckpointData {
            offset: position as u64,
            data,
        });
    }

    fn get_bucket_size(&self) -> u64 {
        self.file_size.load(Ordering::Relaxed)
    }

    fn write_data(&self, bytes: &[u8]) -> u64 {
        let stat_raii = AtomicCounterGuardSum::new(&THREADS_BUSY_WRITING, 1);

        // let _lock = self.checkpoints.lock();
        let position = self.writer.write_all_parallel(bytes, 1);

        let old_size = self
            .file_size
            .fetch_add(bytes.len() as u64, Ordering::Relaxed);
        if old_size >> self.checkpoint_max_size_log2
            != (old_size + bytes.len() as u64) >> self.checkpoint_max_size_log2
        {
            self.checkpoints.lock().push(CheckpointData {
                offset: position as u64,
                data: self.checkpoint_data.lock().clone(),
            });
        }

        drop(stat_raii);

        old_size + bytes.len() as u64
    }

    fn get_path(&self) -> PathBuf {
        self.writer.get_path()
    }
    fn finalize(self) {
        finalize_bucket_file(
            self.writer,
            LOCK_FREE_BUCKET_MAGIC,
            {
                let mut checkpoints = self.checkpoints.into_inner();
                checkpoints.sort();
                checkpoints.dedup();
                checkpoints
            },
            &self.data_format_info,
        );
    }
}
