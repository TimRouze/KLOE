use rustc_hash::FxHashMap;
use std::fs::File;

use std::path::{Path, PathBuf};

use std::sync::Arc;

use crate::memory_fs::file::flush::*;

use parking_lot::Mutex;

use crate::memory_data_size::MemoryDataSize;
use crate::memory_fs::allocator::{MemoryAllocationLimits, CHUNKS_ALLOCATOR};
use crate::memory_fs::file::internal::{MemoryFileInternal, SWAPPABLE_FILES};

pub const O_DIRECT: i32 = 0x4000;

#[macro_use]
pub mod allocator;
pub mod file;
pub mod flushable_buffer;
pub(crate) mod stats;

static FILES_FLUSH_HASH_MAP: Mutex<Option<FxHashMap<PathBuf, Vec<Arc<(PathBuf, Mutex<File>)>>>>> =
    Mutex::new(None);

pub struct MemoryFs;

#[derive(Copy, Clone, Debug)]
pub enum RemoveFileMode {
    Keep,
    Remove { remove_fs: bool },
}

impl MemoryFs {
    pub fn init(
        memory_size: MemoryDataSize,
        flush_queue_size: usize,
        threads_count: usize,
        min_chunks_count: usize,
        alloc_limits: Option<MemoryAllocationLimits>,
    ) {
        let chunk_size = (memory_size / (min_chunks_count as f64)).as_bytes();

        let mut suggested_chunk_size_log = 1;

        while (1 << (suggested_chunk_size_log + 1)) <= chunk_size {
            suggested_chunk_size_log += 1;
        }

        CHUNKS_ALLOCATOR.initialize(
            memory_size,
            suggested_chunk_size_log,
            min_chunks_count,
            alloc_limits,
        );
        *FILES_FLUSH_HASH_MAP.lock() = Some(FxHashMap::with_capacity_and_hasher(
            8192,
            Default::default(),
        ));
        GlobalFlush::init(flush_queue_size, threads_count);
    }

    pub fn remove_file(file: impl AsRef<Path>, remove_mode: RemoveFileMode) -> Result<(), ()> {
        match remove_mode {
            RemoveFileMode::Keep => {
                // Do nothing
                Ok(())
            }
            RemoveFileMode::Remove { remove_fs } => {
                if MemoryFileInternal::delete(file, remove_fs) {
                    Ok(())
                } else {
                    Err(())
                }
            }
        }
    }

    pub fn get_file_size(file: impl AsRef<Path>) -> Option<usize> {
        MemoryFileInternal::retrieve_reference(&file)
            .map(|f| f.read().len())
            .or_else(|| std::fs::metadata(&file).map(|m| m.len() as usize).ok())
    }

    pub fn ensure_flushed(file: impl AsRef<Path>) {
        FILES_FLUSH_HASH_MAP
            .lock()
            .as_mut()
            .unwrap()
            .remove(&file.as_ref().to_path_buf());
    }

    pub fn remove_directory(dir: impl AsRef<Path>, remove_fs: bool) -> bool {
        MemoryFileInternal::delete_directory(dir, remove_fs)
    }

    pub fn flush_to_disk(flush_all: bool) {
        if flush_all {
            MemoryFileInternal::flush_all_to_disk();
        }
        GlobalFlush::flush_to_disk();
    }

    pub fn free_memory() {
        CHUNKS_ALLOCATOR.giveback_free_memory()
    }

    pub fn terminate() {
        GlobalFlush::terminate();
        CHUNKS_ALLOCATOR.deinitialize();
    }

    pub fn get_stats() -> stats::MemoryFsStats {
        stats::get_stats()
    }

    pub fn stats_reset() {
        stats::reset();
    }

    pub fn reduce_pressure() -> bool {
        // crate::log_info!("Reducing pressure!");
        let (current, max_size) = GlobalFlush::global_queue_occupation();
        if current * 3 < max_size {
            let mut map_lock = SWAPPABLE_FILES.lock();
            if let Some(file) = map_lock.get_next() {
                drop(map_lock);
                let mut file = file.write();
                file.change_to_disk_only();
                file.flush_chunks(usize::MAX);
                return true;
            }
        }

        return !GlobalFlush::is_queue_empty();
    }
}

#[cfg(test)]
mod tests {
    use crate::memory_data_size::MemoryDataSize;
    use crate::memory_fs::file::flush::GlobalFlush;
    use crate::memory_fs::file::internal::MemoryFileMode;
    use crate::memory_fs::file::reader::FileReader;
    use crate::memory_fs::file::writer::FileWriter;
    use crate::memory_fs::MemoryFs;
    use rayon::prelude::*;
    use std::io::{Read, Seek, SeekFrom, Write};

    #[test]
    #[ignore]
    pub fn memory_fs_test() {
        MemoryFs::init(MemoryDataSize::from_mebioctets(100 * 1024), 1024, 3, 0);
        let data = (0..3337).map(|x| (x % 256) as u8).collect::<Vec<u8>>();

        (0..400).into_par_iter().for_each(|i: u32| {
            crate::log_info!("Writing file {}", i);
            let mut file = FileWriter::create(
                format!("/home/andrea/genome-assembly/test1234/{}.tmp", i),
                MemoryFileMode::PreferMemory { swap_priority: 3 },
            );
            for _ in 0..(1024 * 64) {
                file.write(data.as_slice()).unwrap();
            }
            drop(file);
            let mut file2 = FileReader::open(
                format!("/home/andrea/genome-assembly/test1234/{}.tmp", i),
                None,
            )
            .unwrap();

            file2.seek(SeekFrom::Start(17 + 3337 * 12374)).unwrap();
            let mut buffer = [0; 4];
            file2.read_exact(&mut buffer).unwrap();
            assert_eq!(&buffer, &data[17..21]);
        });

        GlobalFlush::flush_to_disk();

        (0..400).into_par_iter().for_each(|i: u32| {
            crate::log_info!("Reading file {}", i);
            let mut datar = vec![0; 3337];
            let mut file = FileReader::open(
                format!("/home/andrea/genome-assembly/test1234/{}.tmp", i),
                None,
            )
            .unwrap();
            for _ in 0..(1024 * 64) {
                file.read_exact(datar.as_mut_slice()).unwrap();
                assert_eq!(datar, data);
            }
            assert_eq!(file.read(datar.as_mut_slice()).unwrap(), 0);
            crate::log_info!("Read file {}", i);
        });

        MemoryFs::terminate();
    }
}
