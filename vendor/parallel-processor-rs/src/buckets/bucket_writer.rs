use std::{
    io::{Read, Write},
    marker::PhantomData,
};

pub trait BucketItemSerializer {
    type InputElementType<'a>: ?Sized;
    type ExtraData;
    type ReadBuffer;
    type ExtraDataBuffer;
    type ReadType<'a>;

    /// Creates a new instance
    fn new() -> Self;
    /// Reset on non continuous data
    fn reset(&mut self);

    fn write_to(
        &mut self,
        element: &Self::InputElementType<'_>,
        bucket: &mut Vec<u8>,
        extra_data: &Self::ExtraData,
        extra_read_buffer: &Self::ExtraDataBuffer,
    );
    fn read_from<'a, S: Read>(
        &mut self,
        stream: S,
        read_buffer: &'a mut Self::ReadBuffer,
        extra_read_buffer: &mut Self::ExtraDataBuffer,
    ) -> Option<Self::ReadType<'a>>;

    fn get_size(&self, element: &Self::InputElementType<'_>, extra: &Self::ExtraData) -> usize;
}

pub struct BytesArraySerializer<const SIZE: usize>(PhantomData<[(); SIZE]>);
impl<const SIZE: usize> BucketItemSerializer for BytesArraySerializer<SIZE> {
    type InputElementType<'a> = [u8; SIZE];
    type ExtraData = ();
    type ExtraDataBuffer = ();
    type ReadBuffer = [u8; SIZE];
    type ReadType<'a> = &'a [u8; SIZE];

    #[inline(always)]
    fn new() -> Self {
        Self(PhantomData)
    }

    #[inline(always)]
    fn reset(&mut self) {}

    #[inline(always)]
    fn write_to(
        &mut self,
        element: &Self::InputElementType<'_>,
        bucket: &mut Vec<u8>,
        _: &Self::ExtraData,
        _: &Self::ExtraDataBuffer,
    ) {
        bucket.write(element).unwrap();
    }

    fn read_from<'a, S: Read>(
        &mut self,
        mut stream: S,
        read_buffer: &'a mut Self::ReadBuffer,
        _: &mut Self::ExtraDataBuffer,
    ) -> Option<Self::ReadType<'a>> {
        stream.read_exact(read_buffer).ok()?;
        Some(read_buffer)
    }

    #[inline(always)]
    fn get_size(&self, element: &Self::InputElementType<'_>, _: &()) -> usize {
        element.len()
    }
}

pub struct BytesSliceSerializer;
impl BucketItemSerializer for BytesSliceSerializer {
    type InputElementType<'a> = [u8];
    type ExtraData = ();
    type ExtraDataBuffer = ();
    type ReadBuffer = ();
    type ReadType<'a> = ();

    #[inline(always)]
    fn new() -> Self {
        Self
    }

    #[inline(always)]
    fn reset(&mut self) {}

    #[inline(always)]
    fn write_to(
        &mut self,
        element: &Self::InputElementType<'_>,
        bucket: &mut Vec<u8>,
        _extra_data: &Self::ExtraData,
        _: &Self::ExtraDataBuffer,
    ) {
        bucket.write(element).unwrap();
    }

    fn read_from<'a, S: Read>(
        &mut self,
        _stream: S,
        _read_buffer: &'a mut Self::ReadBuffer,
        _: &mut Self::ExtraDataBuffer,
    ) -> Option<Self::ReadType<'a>> {
        unimplemented!("Cannot read slices of unknown size!")
    }

    #[inline(always)]
    fn get_size(&self, element: &Self::InputElementType<'_>, _: &()) -> usize {
        element.len()
    }
}
