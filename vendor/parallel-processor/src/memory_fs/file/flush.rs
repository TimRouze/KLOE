use crate::memory_fs::allocator::AllocatedChunk;
use crate::memory_fs::file::internal::FileChunk;
use crate::memory_fs::flushable_buffer::{FileFlushMode, FlushableItem};
use crate::memory_fs::stats;
use crossbeam::channel::*;
use mt_debug_counters::counter::{AtomicCounter, AtomicCounterGuardSum, MaxMode, SumMode};
use parking_lot::lock_api::{RawMutex, RawRwLock};
use parking_lot::{Mutex, RwLock};
use std::cmp::max;
use std::io::{Seek, SeekFrom, Write};
use std::ops::DerefMut;
use std::sync::{Arc, atomic::{AtomicUsize, Ordering}};
use std::thread::JoinHandle;
use std::time::Duration;

use super::handle::FileHandle;

static mut GLOBAL_FLUSH_QUEUE: Option<Sender<FlushableItem>> = None;
static FLUSH_THREADS: Mutex<Vec<JoinHandle<()>>> =
    Mutex::const_new(parking_lot::RawMutex::INIT, vec![]);

static TAKE_FROM_QUEUE_MUTEX: Mutex<()> = Mutex::const_new(parking_lot::RawMutex::INIT, ());
static WRITING_CHECK: RwLock<()> = RwLock::const_new(parking_lot::RawRwLock::INIT, ());

// Includes queued items, senders blocked on the bounded queue, and items already
// taken by a writer. An empty queue alone does not mean memory cannot be freed.
static PENDING_WRITES: AtomicUsize = AtomicUsize::new(0);

static COUNTER_WRITING_APPEND: AtomicCounter<SumMode> =
    declare_counter_i64!("threads_file_append_count", SumMode, false);

static COUNTER_WRITE_AT: AtomicCounter<SumMode> =
    declare_counter_i64!("threads_write_at_count", SumMode, false);

static COUNTER_DISK_FLUSHES: AtomicCounter<SumMode> =
    declare_counter_i64!("disk_flushes", SumMode, false);

static COUNTER_BYTES_WRITTEN: AtomicCounter<SumMode> =
    declare_counter_i64!("bytes_written_count", SumMode, false);

static GLOBAL_QUEUE_MAX_SIZE_NOW: AtomicCounter<MaxMode> =
    declare_counter_i64!("global_queue_max_size_now", MaxMode, true);

pub struct GlobalFlush;

#[allow(static_mut_refs)]
impl GlobalFlush {
    pub fn global_queue_occupation() -> (usize, usize) {
        unsafe {
            (
                GLOBAL_FLUSH_QUEUE.as_ref().unwrap().len(),
                GLOBAL_FLUSH_QUEUE
                    .as_ref()
                    .unwrap()
                    .capacity()
                    .unwrap_or(usize::MAX),
            )
        }
    }

    pub fn is_queue_empty() -> bool {
        unsafe { GLOBAL_FLUSH_QUEUE.as_ref().unwrap().len() == 0 }
    }

    pub fn has_pending_writes() -> bool {
        PENDING_WRITES.load(Ordering::Acquire) != 0
    }

    pub fn add_item_to_flush_queue(item: FlushableItem) {
        PENDING_WRITES.fetch_add(1, Ordering::AcqRel);
        unsafe { GLOBAL_FLUSH_QUEUE.as_ref().unwrap().send(item).unwrap() }
    }

    fn flush_thread(flush_channel_receiver: Receiver<FlushableItem>) {
        let mut queue_take_lock = TAKE_FROM_QUEUE_MUTEX.lock();

        while let Ok(file_to_flush) = flush_channel_receiver.recv() {
            // Here it's held an exclusive lock to ensure that sequential writes to a file are not processed concurrently
            let mut file_lock = file_to_flush.underlying_file.lock();
            let mut file = file_lock.get_file();
            let _writing_check = WRITING_CHECK.read();
            drop(queue_take_lock);

            match file_to_flush.mode {
                FileFlushMode::Append {
                    chunk: mut file_chunk,
                } => {
                    if let FileChunk::OnMemory { chunk } = file_chunk.deref_mut() {
                        let _stat = AtomicCounterGuardSum::new(&COUNTER_WRITING_APPEND, 1);
                        GLOBAL_QUEUE_MAX_SIZE_NOW.max(flush_channel_receiver.len() as i64);
                        COUNTER_DISK_FLUSHES.inc();

                        let offset = file.stream_position().unwrap();

                        file.write_all(chunk.get())
                            .map_err(|e| {
                                format!(
                                    "Error while writing file: '{}' => {}",
                                    file_lock.get_path().display(),
                                    e
                                )
                            })
                            .unwrap();

                        stats::add_disk_usage(chunk.len() as u64);

                        let len = chunk.len();
                        *file_chunk = FileChunk::OnDisk { offset, len };
                        COUNTER_BYTES_WRITTEN.inc_by(len as i64);
                    }
                }
                FileFlushMode::WriteAt { buffer, offset } => {
                    let _writing_check = WRITING_CHECK.read();
                    GLOBAL_QUEUE_MAX_SIZE_NOW.max(flush_channel_receiver.len() as i64);

                    COUNTER_DISK_FLUSHES.inc();
                    let _stat = AtomicCounterGuardSum::new(&COUNTER_WRITE_AT, 1);
                    file.seek(SeekFrom::Start(offset)).unwrap();
                    file.write_all(buffer.get()).unwrap();
                    COUNTER_BYTES_WRITTEN.inc_by(buffer.get().len() as i64);
                }
            }

            drop(_writing_check);
            drop(file);
            drop(file_lock);
            PENDING_WRITES.fetch_sub(1, Ordering::AcqRel);
            // Try lock the queue again
            queue_take_lock = TAKE_FROM_QUEUE_MUTEX.lock();
        }
    }

    pub fn is_initialized() -> bool {
        return unsafe { GLOBAL_FLUSH_QUEUE.is_some() };
    }

    pub fn init(flush_queue_size: usize, threads_count: usize) {
        let (flush_channel_sender, flush_channel_receiver) = bounded(flush_queue_size);

        unsafe {
            GLOBAL_FLUSH_QUEUE = Some(flush_channel_sender);
        }

        for _ in 0..max(1, threads_count) {
            let flush_channel_receiver = flush_channel_receiver.clone();
            FLUSH_THREADS.lock().push(
                std::thread::Builder::new()
                    .name(String::from("flushing-thread"))
                    .spawn(move || Self::flush_thread(flush_channel_receiver))
                    .unwrap(),
            );
        }
    }

    pub fn schedule_disk_write(file: Arc<Mutex<FileHandle>>, buffer: AllocatedChunk, offset: u64) {
        Self::add_item_to_flush_queue(FlushableItem {
            underlying_file: file,
            mode: FileFlushMode::WriteAt { buffer, offset },
        })
    }

    pub fn flush_to_disk() {
        while Self::has_pending_writes() {
            std::thread::sleep(Duration::from_millis(50));
        }
        // Ensure that no writers are still writing!
        drop(WRITING_CHECK.write());
    }

    pub fn terminate() {
        Self::flush_to_disk();
        unsafe {
            GLOBAL_FLUSH_QUEUE.take();
        }
        let mut threads = FLUSH_THREADS.lock();
        for thread in threads.drain(..) {
            thread.join().unwrap();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::memory_data_size::MemoryDataSize;
    use crate::memory_fs::{MemoryFs, RemoveFileMode};
    use crate::memory_fs::allocator::CHUNKS_ALLOCATOR;
    use crate::memory_fs::file::{internal::MemoryFileMode, writer::FileWriter};
    use std::io::Read;
    use std::sync::mpsc;
    use std::time::{Instant, SystemTime, UNIX_EPOCH};

    fn scratch_path(label: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!("kloe-{label}-{}-{}",
            std::process::id(), SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos()))
    }

    #[test]
    fn memory_pressure_waits_for_dequeued_writer() {
        let path = scratch_path("slow-flush");
        MemoryFs::init(MemoryDataSize::from_mebioctets(8), 2, 1, 32, None);
        let reserved = CHUNKS_ALLOCATOR.get_total_memory().as_bytes();
        let mut chunks: Vec<_> = (0..32)
            .map(|_| CHUNKS_ALLOCATOR.request_chunk(chunk_usage!(TemporarySpace)))
            .collect();
        let writing = chunks.pop().unwrap();
        writing.write_bytes_noextend(b"slow writer").unwrap();
        let handle = Arc::new(Mutex::new(FileHandle::new(path.clone())));
        let held_file = handle.lock();
        // Exhaust the pool, then stop the writer on its file lock after it
        // dequeues the final item. Releasing this write can free a buffer.
        GlobalFlush::add_item_to_flush_queue(FlushableItem {
            underlying_file: handle.clone(),
            mode: FileFlushMode::WriteAt { buffer: writing, offset: 0 },
        });
        let deadline = Instant::now() + Duration::from_secs(5);
        while !GlobalFlush::is_queue_empty() && Instant::now() < deadline {
            std::thread::sleep(Duration::from_millis(1));
        }
        assert!(GlobalFlush::is_queue_empty());
        assert!(GlobalFlush::has_pending_writes());
        assert!(MemoryFs::reduce_pressure(), "active writes must prevent emergency growth");

        let (allocated, allocation) = mpsc::channel();
        let requester = std::thread::spawn(move || {
            allocated.send(CHUNKS_ALLOCATOR.request_chunk(chunk_usage!(TemporarySpace))).unwrap();
        });
        // The former allocator started an extra 256 MiB allocation after
        // eleven 25 ms waits, despite the outstanding disk write.
        assert!(allocation.recv_timeout(Duration::from_millis(400)).is_err(),
            "allocator expanded instead of waiting for the outstanding write");
        assert_eq!(CHUNKS_ALLOCATOR.get_total_memory().as_bytes(), reserved);

        let (done, completion) = mpsc::channel();
        let waiter = std::thread::spawn(move || {
            GlobalFlush::flush_to_disk();
            done.send(()).unwrap();
        });
        assert!(completion.recv_timeout(Duration::from_millis(100)).is_err(),
            "flush returned while a dequeued writer still held its data");
        drop(held_file);
        drop(allocation.recv_timeout(Duration::from_secs(5)).unwrap());
        requester.join().unwrap();
        completion.recv_timeout(Duration::from_secs(5)).unwrap();
        waiter.join().unwrap();
        assert!(!GlobalFlush::has_pending_writes());
        drop(chunks);
        MemoryFs::terminate();
        drop(handle);
        assert_eq!(std::fs::read(&path).unwrap(), b"slow writer");
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn memory_pressure_disk_stream_reuses_fixed_pool() {
        let path = scratch_path("bounded-disk-stream");
        MemoryFs::init(MemoryDataSize::from_mebioctets(8), 2, 1, 32, None);
        let reserved = CHUNKS_ALLOCATOR.get_total_memory().as_bytes();
        let block: Vec<u8> = (0..64 * 1024).map(|i| (i % 251) as u8).collect();
        let mut writer = FileWriter::create(&path, MemoryFileMode::DiskOnly);
        // Stream 64 MiB through an 8 MiB pool and a two-item disk queue.
        for _ in 0..1024 {
            writer.write_all(&block).unwrap();
        }
        drop(writer);
        MemoryFs::flush_to_disk(false);
        assert_eq!(CHUNKS_ALLOCATOR.get_total_memory().as_bytes(), reserved);
        let mut input = std::fs::File::open(&path).unwrap();
        let mut restored = vec![0; block.len()];
        for _ in 0..1024 {
            input.read_exact(&mut restored).unwrap();
            assert_eq!(restored, block);
        }
        assert_eq!(input.read(&mut [0]).unwrap(), 0);
        drop(input);
        MemoryFs::remove_file(&path, RemoveFileMode::Remove { remove_fs: true }).unwrap();
        MemoryFs::terminate();
    }
}
