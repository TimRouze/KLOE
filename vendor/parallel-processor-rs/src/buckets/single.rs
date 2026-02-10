use crate::buckets::bucket_writer::BucketItemSerializer;
use crate::buckets::{LockFreeBucket, MultiThreadBuckets};
use crate::memory_data_size::MemoryDataSize;
use std::path::PathBuf;

pub struct SingleBucketThreadDispatcher<'a, B: LockFreeBucket, S: BucketItemSerializer> {
    buckets: &'a MultiThreadBuckets<B>,
    bucket_index: u16,
    buffer: Vec<u8>,
    serializer: S,
}

impl<'a, B: LockFreeBucket, S: BucketItemSerializer> SingleBucketThreadDispatcher<'a, B, S> {
    pub fn new(
        buffer_size: MemoryDataSize,
        bucket_index: u16,
        buckets: &'a MultiThreadBuckets<B>,
    ) -> Self {
        let buffer = Vec::with_capacity(buffer_size.as_bytes());

        Self {
            buckets,
            bucket_index,
            buffer,
            serializer: S::new(),
        }
    }

    pub fn get_bucket_index(&self) -> u16 {
        self.bucket_index
    }

    pub fn get_path(&self) -> PathBuf {
        self.buckets.get_path(self.bucket_index)
    }

    fn flush_buffer(&mut self) {
        if self.buffer.len() == 0 {
            return;
        }

        self.buckets.add_data(self.bucket_index, &self.buffer);
        self.buffer.clear();
    }

    pub fn add_element_extended(
        &mut self,
        extra_data: &S::ExtraData,
        extra_buffer: &S::ExtraDataBuffer,
        element: &S::InputElementType<'_>,
    ) {
        if self.serializer.get_size(element, extra_data) + self.buffer.len()
            > self.buffer.capacity()
        {
            self.flush_buffer();
            self.serializer.reset();
        }
        self.serializer
            .write_to(element, &mut self.buffer, extra_data, extra_buffer);
    }

    pub fn add_element(&mut self, extra_data: &S::ExtraData, element: &S::InputElementType<'_>)
    where
        S: BucketItemSerializer<ExtraDataBuffer = ()>,
    {
        self.add_element_extended(extra_data, &(), element);
    }

    pub fn finalize(self) {}
}

impl<'a, B: LockFreeBucket, S: BucketItemSerializer> Drop
    for SingleBucketThreadDispatcher<'a, B, S>
{
    fn drop(&mut self) {
        self.flush_buffer();
    }
}
