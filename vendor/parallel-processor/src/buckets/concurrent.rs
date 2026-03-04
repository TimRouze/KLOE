use crate::buckets::bucket_writer::BucketItemSerializer;
use crate::buckets::{BucketsCount, LockFreeBucket, MultiThreadBuckets};
use crate::memory_data_size::MemoryDataSize;
use crate::utils::panic_on_drop::PanicOnDrop;
use std::sync::Arc;

use super::ChunkingStatus;

pub struct BucketsThreadBuffer {
    buffers: Vec<Vec<u8>>,
}

impl BucketsThreadBuffer {
    pub const EMPTY: Self = Self { buffers: vec![] };

    pub fn new(max_buffer_size: MemoryDataSize, buckets_count: &BucketsCount) -> Self {
        let mut buffers = Vec::with_capacity(buckets_count.total_buckets_count);
        let capacity = max_buffer_size.as_bytes();
        for _ in 0..buckets_count.total_buckets_count {
            buffers.push(Vec::with_capacity(capacity));
        }

        Self { buffers }
    }
}

pub struct BucketsThreadDispatcher<B: LockFreeBucket, S: BucketItemSerializer> {
    mtb: Arc<MultiThreadBuckets<B>>,
    thread_data: BucketsThreadBuffer,
    drop_panic: PanicOnDrop,
    serializers: Vec<S>,
}

impl<B: LockFreeBucket, S: BucketItemSerializer> BucketsThreadDispatcher<B, S> {
    pub fn new(
        mtb: &Arc<MultiThreadBuckets<B>>,
        thread_data: BucketsThreadBuffer,
        deserializer_init_data: S::InitData,
    ) -> Self
    where
        S::InitData: Copy,
    {
        assert_eq!(mtb.active_buckets.len(), thread_data.buffers.len());
        Self {
            mtb: mtb.clone(),
            thread_data,
            drop_panic: PanicOnDrop::new("buckets thread dispatcher not finalized"),
            serializers: (0..mtb.active_buckets.len())
                .map(|_| S::new(deserializer_init_data))
                .collect(),
        }
    }

    #[inline(always)]
    pub fn get_buckets_count(&self) -> &BucketsCount {
        self.mtb.get_buckets_count()
    }

    #[inline]
    pub fn add_element_extended(
        &mut self,
        bucket: u16,
        extra_data: &S::ExtraData,
        extra_data_buffer: &S::ExtraDataBuffer,
        element: &S::InputElementType<'_>,
    ) -> ChunkingStatus {
        let bucket_buf = &mut self.thread_data.buffers[bucket as usize];
        let mut chunking_status = ChunkingStatus::SameChunk;
        if self.serializers[bucket as usize].get_size(element, extra_data) + bucket_buf.len()
            > bucket_buf.capacity()
            && bucket_buf.len() > 0
        {
            chunking_status = self.mtb.add_data(bucket, bucket_buf.as_slice());
            bucket_buf.clear();
            self.serializers[bucket as usize].reset();
        }
        self.serializers[bucket as usize].write_to(
            element,
            bucket_buf,
            extra_data,
            extra_data_buffer,
        );
        chunking_status
    }

    #[inline]
    pub fn add_element(
        &mut self,
        bucket: u16,
        extra_data: &S::ExtraData,
        element: &S::InputElementType<'_>,
    ) where
        S: BucketItemSerializer<ExtraDataBuffer = ()>,
    {
        self.add_element_extended(bucket, extra_data, &(), element);
    }

    pub fn finalize(mut self) -> (BucketsThreadBuffer, Arc<MultiThreadBuckets<B>>) {
        for (index, vec) in self.thread_data.buffers.iter_mut().enumerate() {
            if vec.len() == 0 {
                continue;
            }
            self.mtb.add_data(index as u16, vec.as_slice());
            vec.clear();
        }
        self.drop_panic.disengage();
        (self.thread_data, self.mtb)
    }
}
