use crossbeam::channel::{bounded, unbounded, Receiver, Sender};
use parking_lot::Mutex;
use std::{collections::VecDeque, io::Read, sync::Arc, thread::JoinHandle};

use crate::{
    buckets::{
        bucket_writer::BucketItemSerializer,
        readers::{
            binary_reader::{BinaryChunkReader, BinaryReaderChunk, ChunkDecoder, DecoderType},
            compressed_decoder::CompressedStreamDecoder,
            lock_free_decoder::LockFreeStreamDecoder,
        },
    },
    execution_manager::objects_pool::{ObjectsPool, PoolObject, PoolObjectTrait},
};

// Passtrough {
//     file_range: FileRangeReference,
//     checkpoint_data: Option<Vec<u8>>,
// },

struct AsyncBuffer {
    buffer: Vec<u8>,
}

impl PoolObjectTrait for AsyncBuffer {
    type InitData = usize;

    fn allocate_new(init_data: &Self::InitData) -> Self {
        Self {
            buffer: Vec::with_capacity(*init_data),
        }
    }

    fn reset(&mut self) {
        self.buffer.clear();
    }
}

enum AsyncReaderBuffer {
    Decompressed { data: PoolObject<AsyncBuffer> },
    Closed,
}

struct AsyncReaderJob {
    chunks: Vec<BinaryReaderChunk>,
}

pub struct AsyncReaderThread {
    output_data: Receiver<AsyncReaderBuffer>,
    jobs: Sender<AsyncReaderJob>,
    _thread: JoinHandle<()>,
}

impl AsyncReaderThread {
    #[inline]
    fn process_chunk<D: ChunkDecoder>(
        chunk: BinaryReaderChunk,
        buffers_pool: &ObjectsPool<AsyncBuffer>,
        buffer: &mut PoolObject<AsyncBuffer>,
        data_sender: &Sender<AsyncReaderBuffer>,
    ) {
        let mut reader = BinaryChunkReader::<D>::new(chunk);

        'chunk_read: loop {
            let buffer_start = buffer.buffer.len();
            unsafe {
                let capacity = buffer.buffer.capacity();
                buffer.buffer.set_len(capacity);
            }

            let bytes_read = reader.read(&mut buffer.buffer[buffer_start..]).unwrap();
            buffer.buffer.truncate(buffer_start + bytes_read);

            if bytes_read == 0 {
                // Chunk finished, go to next one
                break 'chunk_read;
            }

            // Buffer full, send it
            if buffer.buffer.len() == buffer.buffer.capacity() {
                let buffer = std::mem::replace(buffer, buffers_pool.alloc_object());
                data_sender
                    .send(AsyncReaderBuffer::Decompressed { data: buffer })
                    .unwrap();
            }
        }
    }

    pub fn new(buffers_size: usize, buffers_count: usize) -> Arc<Self> {
        let (data_sender, output_data) = bounded(buffers_count);
        let (jobs, jobs_receiver) = unbounded::<AsyncReaderJob>();

        let thread = std::thread::spawn(move || {
            let buffers_pool: ObjectsPool<AsyncBuffer> =
                ObjectsPool::new(buffers_count + 1, buffers_size);
            while let Ok(job) = jobs_receiver.recv() {
                let mut buffer = buffers_pool.alloc_object();

                for chunk in job.chunks {
                    match chunk.decoder_type {
                        DecoderType::LockFree => Self::process_chunk::<LockFreeStreamDecoder>(
                            chunk,
                            &buffers_pool,
                            &mut buffer,
                            &data_sender,
                        ),
                        DecoderType::Compressed => Self::process_chunk::<CompressedStreamDecoder>(
                            chunk,
                            &buffers_pool,
                            &mut buffer,
                            &data_sender,
                        ),
                    }
                }

                if buffer.buffer.len() > 0 {
                    data_sender
                        .send(AsyncReaderBuffer::Decompressed { data: buffer })
                        .unwrap();
                }

                // Notify the end of the job
                data_sender.send(AsyncReaderBuffer::Closed).unwrap();
            }
        });

        Arc::new(Self {
            output_data,
            jobs,
            _thread: thread,
        })
    }
}

struct AsyncStreamReader {
    reader_thread: Arc<AsyncReaderThread>,
    buffer: PoolObject<AsyncBuffer>,
    position: usize,
    finished: bool,
}

impl Read for AsyncStreamReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let mut total_read = 0;
        while total_read < buf.len() {
            if self.position >= self.buffer.buffer.len() {
                if self.finished {
                    return Ok(0);
                }

                match self.reader_thread.output_data.recv().unwrap() {
                    AsyncReaderBuffer::Decompressed { data } => {
                        self.buffer = data;
                        self.position = 0;
                    }
                    AsyncReaderBuffer::Closed => {
                        self.finished = true;
                        break;
                    }
                }
            }

            let remaining = buf.len() - total_read;
            let available = self.buffer.buffer.len() - self.position;
            let to_copy = remaining.min(available);

            buf[total_read..total_read + to_copy]
                .copy_from_slice(&self.buffer.buffer[self.position..self.position + to_copy]);

            self.position += to_copy;
            total_read += to_copy;
        }

        Ok(total_read)
    }
}

pub struct TypedStreamReaderBuffers<S: BucketItemSerializer> {
    pub buffer: S::ReadBuffer,
    pub extra_buffer: S::ExtraDataBuffer,
}

pub struct TypedStreamReader;

impl TypedStreamReader {
    #[inline(always)]
    fn decode_single_chunk<S: BucketItemSerializer, D: ChunkDecoder>(
        buffer: &mut S::ReadBuffer,
        extra_buffer: &mut S::ExtraDataBuffer,
        deserializer: &mut S,
        chunk: BinaryReaderChunk,
        mut items_callback: impl FnMut(S::ReadType<'_>, &mut S::ExtraDataBuffer),
    ) {
        match chunk.decoder_type {
            DecoderType::LockFree => {
                let mut reader = BinaryChunkReader::<LockFreeStreamDecoder>::new(chunk);
                while let Some(item) = deserializer.read_from(&mut reader, buffer, extra_buffer) {
                    items_callback(item, extra_buffer)
                }
            }
            DecoderType::Compressed => {
                let mut reader = BinaryChunkReader::<CompressedStreamDecoder>::new(chunk);
                while let Some(item) = deserializer.read_from(&mut reader, buffer, extra_buffer) {
                    items_callback(item, extra_buffer)
                }
            }
        }
    }

    pub fn get_items_parallel<S: BucketItemSerializer>(
        deserializer_init_data: S::InitData,
        parallel_chunks_list: &Mutex<VecDeque<BinaryReaderChunk>>,
        mut items_callback: impl FnMut(S::ReadType<'_>, &mut S::ExtraDataBuffer),
    ) -> TypedStreamReaderBuffers<S> {
        let mut deserializer = S::new(deserializer_init_data);
        let mut buffer = S::ReadBuffer::default();
        let mut extra_buffer = S::ExtraDataBuffer::default();

        while let Some(chunk) = {
            let chunk = parallel_chunks_list.lock().pop_front();
            chunk
        } {
            match chunk.decoder_type {
                DecoderType::LockFree => {
                    Self::decode_single_chunk::<S, LockFreeStreamDecoder>(
                        &mut buffer,
                        &mut extra_buffer,
                        &mut deserializer,
                        chunk,
                        &mut items_callback,
                    );
                }
                DecoderType::Compressed => {
                    Self::decode_single_chunk::<S, CompressedStreamDecoder>(
                        &mut buffer,
                        &mut extra_buffer,
                        &mut deserializer,
                        chunk,
                        &mut items_callback,
                    );
                }
            }
        }

        TypedStreamReaderBuffers {
            buffer,
            extra_buffer,
        }
    }

    pub fn get_items<S: BucketItemSerializer>(
        reader_thread: Option<Arc<AsyncReaderThread>>,
        deserializer_init_data: S::InitData,
        chunks_list: Vec<BinaryReaderChunk>,
        mut items_callback: impl FnMut(S::ReadType<'_>, &mut S::ExtraDataBuffer),
    ) -> TypedStreamReaderBuffers<S> {
        let mut deserializer = S::new(deserializer_init_data);
        let mut buffer = S::ReadBuffer::default();
        let mut extra_buffer = S::ExtraDataBuffer::default();

        if let Some(reader_thread) = reader_thread {
            // Send the async job
            reader_thread
                .jobs
                .send(AsyncReaderJob {
                    chunks: chunks_list,
                })
                .unwrap();

            let first_buffer = match reader_thread.output_data.recv().unwrap() {
                AsyncReaderBuffer::Decompressed { data } => data,
                AsyncReaderBuffer::Closed => {
                    return TypedStreamReaderBuffers {
                        buffer,
                        extra_buffer,
                    }
                }
            };

            let mut reader = AsyncStreamReader {
                reader_thread,
                buffer: first_buffer,
                position: 0,
                finished: false,
            };

            while let Some(item) =
                deserializer.read_from(&mut reader, &mut buffer, &mut extra_buffer)
            {
                items_callback(item, &mut extra_buffer)
            }
        } else {
            for chunk in chunks_list {
                match chunk.decoder_type {
                    DecoderType::LockFree => {
                        Self::decode_single_chunk::<S, LockFreeStreamDecoder>(
                            &mut buffer,
                            &mut extra_buffer,
                            &mut deserializer,
                            chunk,
                            &mut items_callback,
                        );
                    }
                    DecoderType::Compressed => {
                        Self::decode_single_chunk::<S, CompressedStreamDecoder>(
                            &mut buffer,
                            &mut extra_buffer,
                            &mut deserializer,
                            chunk,
                            &mut items_callback,
                        );
                    }
                }
            }
        }
        TypedStreamReaderBuffers {
            buffer,
            extra_buffer,
        }
    }
}
