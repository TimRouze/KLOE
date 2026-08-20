use std::fs::File;
use std::io::{self, BufWriter, Read, Seek, SeekFrom, Write};
use std::os::unix::fs::FileExt;
use std::path::Path;

const PACKED_TIGS_MAGIC: &[u8; 4] = b"KTG3";
const PACKED_TIGS_BLOCK_BYTES: usize = 8 * 1024 * 1024;
const OUTPUT_BUFFER_BYTES: usize = 16 * 1024 * 1024;

pub(crate) struct PackedTigsWriter {
    file: BufWriter<File>,
    block: Vec<u8>,
    block_bytes: usize,
}

impl PackedTigsWriter {
    pub(crate) fn create(path: impl AsRef<Path>) -> io::Result<Self> {
        Self::create_with_block_size(path, PACKED_TIGS_BLOCK_BYTES)
    }

    fn create_with_block_size(path: impl AsRef<Path>, block_bytes: usize) -> io::Result<Self> {
        let mut file = BufWriter::with_capacity(OUTPUT_BUFFER_BYTES, File::create(path)?);
        file.write_all(PACKED_TIGS_MAGIC)?;
        Ok(Self {
            file,
            block: Vec::with_capacity(block_bytes),
            block_bytes,
        })
    }

    fn flush_block(&mut self) -> io::Result<()> {
        if self.block.is_empty() {
            return Ok(());
        }
        let compressed = zstd::encode_all(self.block.as_slice(), 1)?;
        let raw_len = u32::try_from(self.block.len()).map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, "packed-tig block exceeds 4 GiB")
        })?;
        let compressed_len = u32::try_from(compressed.len()).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "compressed packed-tig block exceeds 4 GiB",
            )
        })?;
        self.file.write_all(&raw_len.to_le_bytes())?;
        self.file.write_all(&compressed_len.to_le_bytes())?;
        self.file.write_all(&compressed)?;
        self.block.clear();
        Ok(())
    }

    pub(crate) fn finish(mut self) -> io::Result<()> {
        self.flush_block()?;
        self.file.write_all(&0u32.to_le_bytes())?;
        self.file.flush()
    }
}

impl Write for PackedTigsWriter {
    fn write(&mut self, mut input: &[u8]) -> io::Result<usize> {
        let input_len = input.len();
        while !input.is_empty() {
            let available = self.block_bytes - self.block.len();
            let count = available.min(input.len());
            self.block.extend_from_slice(&input[..count]);
            input = &input[count..];
            if self.block.len() == self.block_bytes {
                self.flush_block()?;
            }
        }
        Ok(input_len)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.flush_block()?;
        self.file.flush()
    }
}

#[derive(Debug)]
struct PackedTigsBlock {
    logical_start: u64,
    logical_len: usize,
    data_offset: u64,
    compressed_len: usize,
}

enum PackedTigsEncoding {
    Raw,
    BlockCompressed(Vec<PackedTigsBlock>),
}

pub(crate) struct PackedTigsReader {
    file: File,
    encoding: PackedTigsEncoding,
    logical_len: u64,
    cached_block: usize,
    cached_data: Vec<u8>,
}

impl PackedTigsReader {
    pub(crate) fn open(path: impl AsRef<Path>) -> io::Result<Self> {
        let path = path.as_ref();
        let mut file = File::open(path)?;
        let physical_len = file.metadata()?.len();
        let mut magic = [0u8; 4];
        if file.read_exact(&mut magic).is_err() || &magic != PACKED_TIGS_MAGIC {
            return Ok(Self {
                file: File::open(path)?,
                encoding: PackedTigsEncoding::Raw,
                logical_len: physical_len,
                cached_block: usize::MAX,
                cached_data: Vec::new(),
            });
        }

        let mut blocks = Vec::new();
        let mut logical_len = 0u64;
        loop {
            let mut raw_len_bytes = [0u8; 4];
            file.read_exact(&mut raw_len_bytes)?;
            let raw_len = u32::from_le_bytes(raw_len_bytes) as usize;
            if raw_len == 0 {
                break;
            }
            let mut compressed_len_bytes = [0u8; 4];
            file.read_exact(&mut compressed_len_bytes)?;
            let compressed_len = u32::from_le_bytes(compressed_len_bytes) as usize;
            if compressed_len == 0 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "packed-tig block has an empty compressed payload",
                ));
            }
            let data_offset = file.stream_position()?;
            blocks.push(PackedTigsBlock {
                logical_start: logical_len,
                logical_len: raw_len,
                data_offset,
                compressed_len,
            });
            logical_len = logical_len.checked_add(raw_len as u64).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "packed-tig logical size overflow",
                )
            })?;
            file.seek(SeekFrom::Current(compressed_len as i64))?;
        }
        Ok(Self {
            file: File::open(path)?,
            encoding: PackedTigsEncoding::BlockCompressed(blocks),
            logical_len,
            cached_block: usize::MAX,
            cached_data: Vec::new(),
        })
    }

    pub(crate) fn logical_len(&self) -> u64 {
        self.logical_len
    }

    pub(crate) fn read_exact_at(
        &mut self,
        mut logical_offset: u64,
        mut output: &mut [u8],
    ) -> io::Result<()> {
        let requested_end = logical_offset
            .checked_add(output.len() as u64)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "tig read overflow"))?;
        if requested_end > self.logical_len {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "packed-tig read extends past the logical payload",
            ));
        }
        match &self.encoding {
            PackedTigsEncoding::Raw => self.file.read_exact_at(output, logical_offset),
            PackedTigsEncoding::BlockCompressed(_) => {
                while !output.is_empty() {
                    let block_index = self.block_index(logical_offset)?;
                    self.load_block(block_index)?;
                    let blocks = match &self.encoding {
                        PackedTigsEncoding::BlockCompressed(blocks) => blocks,
                        PackedTigsEncoding::Raw => unreachable!(),
                    };
                    let block = &blocks[block_index];
                    let local = (logical_offset - block.logical_start) as usize;
                    let count = output.len().min(block.logical_len - local);
                    output[..count].copy_from_slice(&self.cached_data[local..local + count]);
                    logical_offset += count as u64;
                    output = &mut output[count..];
                }
                Ok(())
            }
        }
    }

    fn block_index(&self, logical_offset: u64) -> io::Result<usize> {
        let blocks = match &self.encoding {
            PackedTigsEncoding::BlockCompressed(blocks) => blocks,
            PackedTigsEncoding::Raw => unreachable!(),
        };
        let index = blocks
            .partition_point(|block| block.logical_start <= logical_offset)
            .checked_sub(1)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "packed-tig offset is not indexed",
                )
            })?;
        let block = &blocks[index];
        if logical_offset >= block.logical_start + block.logical_len as u64 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "packed-tig offset falls between indexed blocks",
            ));
        }
        Ok(index)
    }

    fn load_block(&mut self, block_index: usize) -> io::Result<()> {
        if self.cached_block == block_index {
            return Ok(());
        }
        let block = match &self.encoding {
            PackedTigsEncoding::BlockCompressed(blocks) => &blocks[block_index],
            PackedTigsEncoding::Raw => unreachable!(),
        };
        let mut compressed = vec![0u8; block.compressed_len];
        self.file
            .read_exact_at(&mut compressed, block.data_offset)?;
        self.cached_data = zstd::decode_all(compressed.as_slice())?;
        if self.cached_data.len() != block.logical_len {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "packed-tig block decompressed to the wrong size",
            ));
        }
        self.cached_block = block_index;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn block_compressed_roundtrip_supports_cross_block_reads() {
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("tigs.bin");
        let input = (0..251u16).map(|value| value as u8).collect::<Vec<_>>();
        let mut writer = PackedTigsWriter::create_with_block_size(&path, 31).unwrap();
        writer.write_all(&input[..67]).unwrap();
        writer.write_all(&input[67..]).unwrap();
        writer.finish().unwrap();

        let mut reader = PackedTigsReader::open(&path).unwrap();
        assert_eq!(reader.logical_len(), input.len() as u64);
        let mut output = vec![0u8; 143];
        reader.read_exact_at(19, &mut output).unwrap();
        assert_eq!(output, input[19..162]);
    }

    #[test]
    fn raw_archives_remain_readable() {
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("legacy-tigs.bin");
        let input = (0..97u8).collect::<Vec<_>>();
        std::fs::write(&path, &input).unwrap();

        let mut reader = PackedTigsReader::open(&path).unwrap();
        assert_eq!(reader.logical_len(), input.len() as u64);
        let mut output = vec![0u8; 41];
        reader.read_exact_at(13, &mut output).unwrap();
        assert_eq!(output, input[13..54]);
    }
}
