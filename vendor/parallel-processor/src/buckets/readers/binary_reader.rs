use crate::buckets::readers::compressed_decoder::CompressedStreamDecoder;
use crate::buckets::readers::lock_free_decoder::LockFreeStreamDecoder;
use crate::buckets::writers::{BucketCheckpoints, BucketHeader};
use crate::memory_fs::file::reader::{FileRangeReference, FileReader};
use crate::memory_fs::{MemoryFs, RemoveFileMode};
use crate::DEFAULT_BINCODE_CONFIG;
use bincode::Decode;
use desse::Desse;
use desse::DesseSized;
use parking_lot::Mutex;
use std::collections::VecDeque;
use std::fmt::Debug;
use std::io::{Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::Arc;

pub trait ChunkDecoder {
    const MAGIC_HEADER: &'static [u8; 16];
    type ReaderType: Read;
    fn decode_stream(reader: FileReader, size: u64) -> Self::ReaderType;
    fn dispose_stream(stream: Self::ReaderType) -> FileReader;
}

pub(crate) struct RemoveFileGuard {
    path: PathBuf,
    remove_mode: RemoveFileMode,
}

impl Drop for RemoveFileGuard {
    fn drop(&mut self) {
        MemoryFs::remove_file(&self.path, self.remove_mode).unwrap();
    }
}

#[derive(Clone)]
pub struct BinaryReaderChunk {
    pub(crate) reader: FileReader,
    pub(crate) length: u64,
    pub(crate) remove_guard: Arc<RemoveFileGuard>,
    pub(crate) extra_data: Option<Vec<u8>>,
    pub(crate) decoder_type: DecoderType,
}

impl Debug for BinaryReaderChunk {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BinaryReaderChunk")
            .field("length", &self.length)
            .field("extra_data", &self.extra_data)
            .field("decoder_type", &self.decoder_type)
            .field("file_name", &self.reader.get_file_path().display())
            .finish()
    }
}

impl BinaryReaderChunk {
    pub fn get_extra_data<E: Decode<()>>(&self) -> Option<E> {
        self.extra_data.as_ref().map(|data| {
            bincode::decode_from_slice(&data, DEFAULT_BINCODE_CONFIG)
                .unwrap()
                .0
        })
    }

    pub fn get_decoder_type(&self) -> DecoderType {
        self.decoder_type
    }

    pub fn get_unique_file_id(&self) -> usize {
        self.reader.get_unique_file_id()
    }

    pub fn get_length(&self) -> u64 {
        self.length
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum DecoderType {
    LockFree,
    Compressed,
}

pub struct ChunkedBinaryReaderIndex {
    file_path: PathBuf,
    chunks: Vec<BinaryReaderChunk>,
    format_data_info: Vec<u8>,
    file_size: u64,
}

pub enum DecodeItemsStatus<T> {
    Decompressed,
    Passtrough {
        file_range: FileRangeReference,
        data: Option<T>,
    },
}

impl ChunkedBinaryReaderIndex {
    fn get_chunk_size(
        checkpoints: &BucketCheckpoints,
        last_byte_position: u64,
        index: usize,
    ) -> u64 {
        if checkpoints.index.len() > (index + 1) as usize {
            checkpoints.index[(index + 1) as usize].offset
                - checkpoints.index[index as usize].offset
        } else {
            last_byte_position - checkpoints.index[index as usize].offset
        }
    }

    pub fn get_file_size(&self) -> u64 {
        self.file_size
    }

    pub fn from_file(
        name: impl AsRef<Path>,
        remove_file: RemoveFileMode,
        prefetch_amount: Option<usize>,
    ) -> Self {
        let mut file = FileReader::open(&name, prefetch_amount)
            .unwrap_or_else(|| panic!("Cannot open file {}", name.as_ref().display()));

        let mut header_buffer = [0; BucketHeader::SIZE];
        file.read_exact(&mut header_buffer)
            .unwrap_or_else(|_| panic!("File {} is corrupted", name.as_ref().display()));

        let header: BucketHeader = BucketHeader::deserialize_from(&header_buffer);

        file.seek(SeekFrom::Start(header.index_offset)).unwrap();
        let index: BucketCheckpoints =
            bincode::decode_from_std_read(&mut file, DEFAULT_BINCODE_CONFIG).unwrap();

        let remove_guard = Arc::new(RemoveFileGuard {
            path: name.as_ref().to_path_buf(),
            remove_mode: remove_file,
        });

        let mut chunks = Vec::with_capacity(index.index.len());

        let decoder_type = match &header.magic {
            LockFreeStreamDecoder::MAGIC_HEADER => DecoderType::LockFree,
            CompressedStreamDecoder::MAGIC_HEADER => DecoderType::Compressed,
            _ => panic!(
                "Invalid decode bucket header magic: {:?}. This is a bug",
                header.magic
            ),
        };

        if index.index.len() > 0 {
            file.seek(SeekFrom::Start(index.index[0].offset)).unwrap();

            for (chunk_idx, chunk) in index.index.iter().enumerate() {
                let length = Self::get_chunk_size(&index, header.index_offset, chunk_idx);

                chunks.push(BinaryReaderChunk {
                    reader: file.clone(),
                    length,
                    remove_guard: remove_guard.clone(),
                    extra_data: chunk.data.clone(),
                    decoder_type,
                });

                file.seek(SeekFrom::Current(length as i64)).unwrap();
            }
        }

        Self {
            file_path: name.as_ref().to_path_buf(),
            chunks,
            format_data_info: header.data_format_info.to_vec(),
            file_size: MemoryFs::get_file_size(name).unwrap() as u64,
        }
    }

    pub fn get_data_format_info<T: Decode<()>>(&self) -> T {
        bincode::decode_from_slice(&self.format_data_info, DEFAULT_BINCODE_CONFIG)
            .unwrap()
            .0
    }

    pub fn get_path(&self) -> &Path {
        &self.file_path
    }

    pub fn into_chunks(self) -> Vec<BinaryReaderChunk> {
        self.chunks
    }

    pub fn into_parallel_chunks(self) -> Mutex<VecDeque<BinaryReaderChunk>> {
        Mutex::new(self.chunks.into())
    }
}

pub struct BinaryChunkReader<D: ChunkDecoder> {
    reader: D::ReaderType,
    _remove_guard: Arc<RemoveFileGuard>,
}

impl<D: ChunkDecoder> BinaryChunkReader<D> {
    pub fn new(chunk: BinaryReaderChunk) -> Self {
        Self {
            reader: D::decode_stream(chunk.reader, chunk.length),
            _remove_guard: chunk.remove_guard,
        }
    }
}

impl<D: ChunkDecoder> Read for BinaryChunkReader<D> {
    #[inline(always)]
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        self.reader.read(buf)
    }
}
