use arc_swap::ArcSwap;
use bincode::{Decode, Encode};
use parking_lot::Mutex;
use std::num::NonZeroU64;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::memory_fs::file::reader::FileRangeReference;
use crate::DEFAULT_BINCODE_CONFIG;

pub mod bucket_writer;
pub mod concurrent;
pub mod readers;
pub mod single;
pub mod writers;

/// This enum serves as a way to specifying the behavior of the
/// bucket portions created after setting the checkpoint data.
/// If set on passtrough there is the option to directly read binary data and copy it somewhere else
#[derive(Encode, Decode, Copy, Clone, Debug, PartialEq, Eq)]
pub enum CheckpointStrategy {
    Decompress,
    Passtrough,
}

#[derive(Encode, Decode, Clone, Debug, PartialEq, Eq)]
pub(crate) struct CheckpointData {
    offset: u64,
    data: Option<Vec<u8>>,
}

impl PartialOrd for CheckpointData {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.offset.partial_cmp(&other.offset)
    }
}

impl Ord for CheckpointData {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.offset.cmp(&other.offset)
    }
}

pub trait LockFreeBucket: Sized {
    type InitData: Clone;

    fn new_serialized_data_format(
        path: &Path,
        data: &Self::InitData,
        index: usize,
        data_format: &[u8],
    ) -> Self;

    fn new<T: Encode>(path: &Path, data: &Self::InitData, index: usize, data_format: &T) -> Self {
        Self::new_serialized_data_format(
            path,
            data,
            index,
            &bincode::encode_to_vec(data_format, DEFAULT_BINCODE_CONFIG).unwrap(),
        )
    }

    fn set_checkpoint_data<T: Encode>(
        &self,
        data: Option<&T>,
        passtrough_range: Option<FileRangeReference>,
    );

    fn get_bucket_size(&self) -> u64;
    fn write_data(&self, bytes: &[u8]) -> u64;
    fn get_path(&self) -> PathBuf;
    fn finalize(self);
}

#[derive(Encode, Decode, Debug, Clone)]
pub struct MultiChunkBucket {
    pub index: usize,
    pub chunks: Vec<PathBuf>,
    pub extra_bucket_data: Option<ExtraBucketData>,
}

impl MultiChunkBucket {
    pub fn into_single(mut self) -> SingleBucket {
        assert!(self.chunks.len() == 1);
        SingleBucket {
            index: self.index,
            path: self.chunks.pop().unwrap(),
            extra_bucket_data: self.extra_bucket_data,
        }
    }
}

#[derive(Encode, Decode)]
pub struct SingleBucket {
    pub index: usize,
    pub path: PathBuf,
    pub extra_bucket_data: Option<ExtraBucketData>,
}

impl SingleBucket {
    pub fn to_multi_chunk(self) -> MultiChunkBucket {
        MultiChunkBucket {
            index: self.index,
            chunks: vec![self.path],
            extra_bucket_data: self.extra_bucket_data,
        }
    }
}

#[derive(Debug, Clone)]
pub enum ChunkingStatus {
    SameChunk,
    NewChunk,
}

#[derive(Encode, Decode, Copy, Clone, Debug, PartialEq, Eq)]
pub struct BucketsCount {
    pub normal_buckets_count: usize,
    pub normal_buckets_count_log: usize,
    pub total_buckets_count: usize,
    pub extra_buckets_count: ExtraBuckets,
}

impl BucketsCount {
    pub const ONE: Self = Self::new(0, ExtraBuckets::None);

    pub fn from_power_of_two(size: usize, extra_buckets: ExtraBuckets) -> Self {
        assert_eq!(size, size.next_power_of_two());
        Self::new(size.ilog2() as usize, extra_buckets)
    }

    pub const fn new(size_log: usize, extra_buckets: ExtraBuckets) -> Self {
        let normal_buckets_count = 1 << size_log;
        let extra_buckets_count = match extra_buckets {
            ExtraBuckets::None => 0,
            ExtraBuckets::Extra { count, .. } => count,
        };

        Self {
            normal_buckets_count,
            normal_buckets_count_log: size_log,
            total_buckets_count: normal_buckets_count + extra_buckets_count,
            extra_buckets_count: extra_buckets,
        }
    }

    pub fn get_extra_buckets_count(&self) -> usize {
        match self.extra_buckets_count {
            ExtraBuckets::None => 0,
            ExtraBuckets::Extra { count, .. } => count,
        }
    }
}

pub struct MultiThreadBuckets<B: LockFreeBucket> {
    active_buckets: Vec<ArcSwap<B>>,
    stored_buckets: Vec<Mutex<MultiChunkBucket>>,
    chunking_size_threshold: Option<NonZeroU64>,
    bucket_count_lock: Mutex<usize>,
    base_path: Option<PathBuf>,
    init_data: Option<B::InitData>,
    serialized_format_info: Vec<u8>,
    size: BucketsCount,
}

#[derive(Encode, Decode, Copy, Clone, Debug, PartialEq, Eq)]
pub struct ExtraBucketData(pub usize);

#[derive(Encode, Decode, Copy, Clone, Debug, PartialEq, Eq)]
pub enum ExtraBuckets {
    None,
    Extra { count: usize, data: ExtraBucketData },
}

impl<B: LockFreeBucket> MultiThreadBuckets<B> {
    pub const EMPTY: Self = Self {
        active_buckets: vec![],
        stored_buckets: vec![],
        chunking_size_threshold: None,
        bucket_count_lock: Mutex::new(0),
        base_path: None,
        init_data: None,
        serialized_format_info: vec![],
        size: BucketsCount {
            normal_buckets_count: 0,
            normal_buckets_count_log: 0,
            total_buckets_count: 0,
            extra_buckets_count: ExtraBuckets::None,
        },
    };

    pub fn get_buckets_count(&self) -> &BucketsCount {
        &self.size
    }

    pub fn create_matching_multichunks(&self) -> Vec<Mutex<MultiChunkBucket>> {
        self.stored_buckets
            .iter()
            .map(|s| {
                let s = s.lock();
                Mutex::new(MultiChunkBucket {
                    index: s.index,
                    chunks: vec![],
                    extra_bucket_data: s.extra_bucket_data.clone(),
                })
            })
            .collect()
    }

    pub fn new(
        size: BucketsCount,
        path: PathBuf,
        chunking_size_threshold: Option<u64>,
        init_data: &B::InitData,
        format_info: &impl Encode,
    ) -> MultiThreadBuckets<B> {
        let mut buckets = Vec::with_capacity(size.total_buckets_count);

        for i in 0..size.total_buckets_count {
            buckets.push(ArcSwap::from_pointee(B::new(
                &path,
                init_data,
                i,
                format_info,
            )));
        }
        MultiThreadBuckets {
            active_buckets: buckets,
            stored_buckets: (0..size.total_buckets_count)
                .map(|index| {
                    Mutex::new(MultiChunkBucket {
                        index,
                        chunks: vec![],
                        extra_bucket_data: match size.extra_buckets_count {
                            ExtraBuckets::None => None,
                            ExtraBuckets::Extra { data, .. } => {
                                if index >= size.normal_buckets_count {
                                    Some(data)
                                } else {
                                    None
                                }
                            }
                        },
                    })
                })
                .collect(),
            chunking_size_threshold: chunking_size_threshold.map(NonZeroU64::new).flatten(),
            bucket_count_lock: Mutex::new(size.total_buckets_count),
            base_path: Some(path),
            init_data: Some(init_data.clone()),
            serialized_format_info: bincode::encode_to_vec(format_info, DEFAULT_BINCODE_CONFIG)
                .unwrap(),
            size,
        }
    }

    pub fn get_stored_buckets(&self) -> &Vec<Mutex<MultiChunkBucket>> {
        &self.stored_buckets
    }

    pub fn into_buckets(mut self) -> impl Iterator<Item = B> {
        assert!(
            self.stored_buckets
                .iter_mut()
                .all(|bucket| bucket.get_mut().chunks.is_empty())
                && self.chunking_size_threshold.is_none()
        );
        let buckets = std::mem::take(&mut self.active_buckets);
        buckets
            .into_iter()
            .map(|bucket| Arc::into_inner(bucket.into_inner()).unwrap())
    }

    pub fn get_path(&self, bucket: u16) -> PathBuf {
        self.active_buckets[bucket as usize].load().get_path()
    }

    pub fn add_data(&self, index: u16, data: &[u8]) -> ChunkingStatus {
        let bucket_guard = self.active_buckets[index as usize].load();
        let last_bucket_size = bucket_guard.write_data(data);

        drop(bucket_guard);

        // If the disk usage limit is set, check if the disk usage is not exceeded
        if let Some(chunk_threshold) = self.chunking_size_threshold {
            if last_bucket_size >= chunk_threshold.get() {
                let mut buckets_count = self.bucket_count_lock.lock();

                if self.active_buckets[index as usize].load().get_bucket_size()
                    < chunk_threshold.get()
                {
                    // Do not add a new chunk if a previous writer already swapped it
                    return ChunkingStatus::SameChunk;
                }

                // Take the largest bucket and add it to the stored buckets
                let mut stored_bucket = self.active_buckets[index as usize].swap(Arc::new(
                    B::new_serialized_data_format(
                        &self.base_path.as_deref().unwrap(),
                        &self.init_data.as_ref().unwrap(),
                        *buckets_count,
                        &self.serialized_format_info,
                    ),
                ));

                let stored_bucket = loop {
                    // Wait for the bucket to end all the pending writes before finalizing it
                    match Arc::try_unwrap(stored_bucket) {
                        Ok(bucket) => break bucket,
                        Err(waiting_arc) => {
                            stored_bucket = waiting_arc;
                            std::hint::spin_loop();
                        }
                    }
                };

                // Add the bucket to the stored buckets and clear its active usage
                let bucket_path = stored_bucket.get_path();
                stored_bucket.finalize();
                self.stored_buckets[index as usize]
                    .lock()
                    .chunks
                    .push(bucket_path);

                *buckets_count += 1;
                ChunkingStatus::NewChunk
            } else {
                ChunkingStatus::SameChunk
            }
        } else {
            ChunkingStatus::SameChunk
        }
    }

    pub fn finalize_single(self: Arc<Self>) -> Vec<SingleBucket> {
        assert!(self.chunking_size_threshold.is_none());
        let buckets = self.finalize();
        buckets
            .into_iter()
            .map(|mut bucket| {
                assert!(bucket.chunks.len() == 1);
                SingleBucket {
                    index: bucket.index,
                    path: bucket.chunks.pop().unwrap(),
                    extra_bucket_data: bucket.extra_bucket_data,
                }
            })
            .collect()
    }

    pub fn finalize(self: Arc<Self>) -> Vec<MultiChunkBucket> {
        let mut self_ = Arc::try_unwrap(self)
            .unwrap_or_else(|_| panic!("Cannot take full ownership of multi thread buckets!"));

        self_
            .active_buckets
            .drain(..)
            .zip(self_.stored_buckets.drain(..))
            .map(|(bucket, stored)| {
                let mut stored = stored.into_inner();
                let bucket = Arc::into_inner(bucket.into_inner()).unwrap();
                stored.chunks.push(bucket.get_path());
                bucket.finalize();
                stored
            })
            .collect()
    }
}

impl<B: LockFreeBucket> Drop for MultiThreadBuckets<B> {
    fn drop(&mut self) {
        self.active_buckets.drain(..).for_each(|bucket| {
            let bucket = Arc::into_inner(bucket.into_inner()).unwrap();
            bucket.finalize();
        });
    }
}

unsafe impl<B: LockFreeBucket> Send for MultiThreadBuckets<B> {}

unsafe impl<B: LockFreeBucket> Sync for MultiThreadBuckets<B> {}
