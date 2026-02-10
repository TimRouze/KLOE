use crate::buckets::bucket_writer::BucketItemSerializer;
use std::path::PathBuf;

pub mod async_binary_reader;
pub mod compressed_binary_reader;
pub mod generic_binary_reader;
pub mod lock_free_binary_reader;
pub mod unbuffered_compressed_binary_reader;

pub trait BucketReader {
    fn decode_all_bucket_items<
        S: BucketItemSerializer,
        F: for<'a> FnMut(S::ReadType<'a>, &mut S::ExtraDataBuffer),
    >(
        self,
        buffer: S::ReadBuffer,
        extra_buffer: &mut S::ExtraDataBuffer,
        func: F,
    );

    fn get_name(&self) -> PathBuf;
}
