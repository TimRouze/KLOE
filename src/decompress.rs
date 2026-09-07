use core::panic;
use std::cmp::Reverse;
use std::collections::{BTreeMap, BinaryHeap, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Result, Seek, Write};
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::thread;

use crossbeam::channel::{bounded, Receiver};
use flate2::{write::GzEncoder, Compression};
use ggcat_api::{ExtraElaboration, GGCATConfig, GGCATInstance, GeneralSequenceBlockData};
use xz2::write::XzEncoder;
use zstd::stream::write::Encoder as ZstdEncoder;
use zstd::Decoder;

use crate::compress::{
    stable_filename_hash, CidDatasetSidecar, ABUNDANCE_BASE_FILE, ABUNDANCE_BASE_INDEX_FILE,
    BUCKET_SIZES_INDEX_FILE, BUCKET_SIZES_INDEX_MAGIC, CID_TO_DATASET_FILE, DATASET_POSTING_BITMAP,
    DATASET_POSTING_DELTA, DATASET_POSTING_MAJORITY_XOR, DATASET_TO_CID_FILE,
    DATASET_TO_CID_FLAG_ABUNDANCE, DATASET_TO_CID_FLAG_ADAPTIVE, DATASET_TO_CID_FLAG_MAJORITY,
    DATASET_TO_CID_MAGIC, FILENAME_INDEX_FILE, FILENAME_INDEX_MAGIC, POSITIONS_INDEX_FILE,
    POSITIONS_INDEX_MAGIC, TIGS_INDEX_FILE,
};
use crate::packed_tigs::PackedTigsReader;
use crate::utils::vec2str;

const BUCKET_SIZES_MAGIC: &[u8; 4] = b"KSB2";
const BUCKET_SIZES_COMPACT_MAGIC: &[u8; 4] = b"KSB3";
const POSITIONS_MAGIC: &[u8; 4] = b"KPS2";
const POSITIONS_COMPACT_MAGIC: &[u8; 4] = b"KPS3";
const ID_TO_CID_MAGIC: &[u8; 4] = b"KIC2";
const DECOMPRESS_BATCH_BASES: usize = 8 * 1024 * 1024;
const FILENAME_INDEX_HEADER_BYTES: u64 = 20;
const DATASET_BOUNDARY_BYTES: u64 = 8;
const DATASET_NAME_RECORD_BYTES: u64 = 12;
const FILENAME_HASH_SLOT_BYTES: u64 = 12;
const POSITION_INDEX_HEADER_BYTES: u64 = 24;
const SIZE_INDEX_HEADER_BYTES: u64 = 28;
const SIZE_INDEX_RECORD_BYTES: u64 = 32;
const PACKED_BASES: [[u8; 4]; 256] = {
    let mut table = [[b'A'; 4]; 256];
    let mut byte = 0usize;
    while byte < table.len() {
        let mut base = 0usize;
        while base < 4 {
            table[byte][base] = [b'A', b'C', b'G', b'T'][(byte >> (2 * base)) & 0b11];
            base += 1;
        }
        byte += 1;
    }
    table
};

//   =========================================================================================== DECOMPRESSION ==============================================================================

#[derive(Clone, Copy)]
struct FilenameIndexRecord {
    hash: u64,
    dataset_id: u32,
}

struct FilenameIndex {
    index: File,
    names: File,
    count: u64,
    capacity: u64,
    names_offset: u64,
    hash_offset: u64,
}

impl FilenameIndex {
    fn open(index_path: &Path, names_path: &Path) -> Result<Self> {
        let mut index = File::open(index_path)?;
        let mut magic = [0u8; 4];
        index.read_exact(&mut magic)?;
        if &magic != FILENAME_INDEX_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid filename index",
            ));
        }
        let mut count = [0u8; 8];
        index.read_exact(&mut count)?;
        let count = u64::from_le_bytes(count);
        let mut capacity = [0u8; 8];
        index.read_exact(&mut capacity)?;
        let capacity = u64::from_le_bytes(capacity);
        if capacity == 0 || !capacity.is_power_of_two() || count > capacity {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid filename hash table capacity",
            ));
        }
        let names_offset = FILENAME_INDEX_HEADER_BYTES
            .checked_add(
                count
                    .saturating_add(1)
                    .saturating_mul(DATASET_BOUNDARY_BYTES),
            )
            .ok_or_else(|| io::Error::other("dataset directory size overflow"))?;
        let hash_offset = names_offset
            .checked_add(count.saturating_mul(DATASET_NAME_RECORD_BYTES))
            .ok_or_else(|| io::Error::other("dataset directory size overflow"))?;
        let expected = hash_offset
            .checked_add(capacity.saturating_mul(FILENAME_HASH_SLOT_BYTES))
            .ok_or_else(|| io::Error::other("filename index size overflow"))?;
        if index.metadata()?.len() != expected {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid filename index length",
            ));
        }
        Ok(Self {
            index,
            names: File::open(names_path)?,
            count,
            capacity,
            names_offset,
            hash_offset,
        })
    }

    fn record(&self, slot: u64) -> Result<FilenameIndexRecord> {
        let mut bytes = [0u8; FILENAME_HASH_SLOT_BYTES as usize];
        self.index.read_exact_at(
            &mut bytes,
            self.hash_offset + slot * FILENAME_HASH_SLOT_BYTES,
        )?;
        Ok(FilenameIndexRecord {
            hash: u64::from_le_bytes(bytes[..8].try_into().unwrap()),
            dataset_id: u32::from_le_bytes(bytes[8..12].try_into().unwrap()),
        })
    }

    fn name_record(&self, dataset_id: u32) -> Result<(u64, u32)> {
        if u64::from(dataset_id) >= self.count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "filename hash slot contains an invalid dataset ID",
            ));
        }
        let mut bytes = [0u8; DATASET_NAME_RECORD_BYTES as usize];
        self.index.read_exact_at(
            &mut bytes,
            self.names_offset + u64::from(dataset_id) * DATASET_NAME_RECORD_BYTES,
        )?;
        Ok((
            u64::from_le_bytes(bytes[..8].try_into().unwrap()),
            u32::from_le_bytes(bytes[8..12].try_into().unwrap()),
        ))
    }

    fn lookup(&self, name: &str) -> Result<Option<u32>> {
        let hash = stable_filename_hash(name.as_bytes());
        let mut slot = hash & (self.capacity - 1);
        for _ in 0..self.capacity {
            let record = self.record(slot)?;
            if record.dataset_id == u32::MAX {
                return Ok(None);
            }
            if record.hash == hash {
                let (name_offset, name_len) = self.name_record(record.dataset_id)?;
                if name_len as usize == name.len() {
                    let mut stored = vec![0u8; name_len as usize];
                    self.names.read_exact_at(&mut stored, name_offset)?;
                    if stored == name.as_bytes() {
                        return Ok(Some(record.dataset_id));
                    }
                }
            }
            slot = (slot + 1) & (self.capacity - 1);
        }
        Ok(None)
    }
}

struct DatasetPostingIndex {
    index: File,
    data_path: PathBuf,
    dataset_count: u64,
    abundance_log_base: Option<f64>,
    adaptive: bool,
    cid_count: Option<u64>,
    majority_frame: Option<(u64, u64)>,
}

type PostingDecoder = Decoder<'static, BufReader<std::io::Take<File>>>;

impl DatasetPostingIndex {
    fn open(index_path: &Path, data_path: &Path) -> Result<Self> {
        let mut index = File::open(index_path)?;
        let mut magic = [0u8; 4];
        index.read_exact(&mut magic)?;
        if &magic != FILENAME_INDEX_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid dataset directory",
            ));
        }
        let mut count = [0u8; 8];
        index.read_exact(&mut count)?;
        let dataset_count = u64::from_le_bytes(count);
        let mut capacity = [0u8; 8];
        index.read_exact(&mut capacity)?;
        let capacity = u64::from_le_bytes(capacity);
        let expected = FILENAME_INDEX_HEADER_BYTES
            .saturating_add(dataset_count.saturating_add(1).saturating_mul(8))
            .saturating_add(dataset_count.saturating_mul(DATASET_NAME_RECORD_BYTES))
            .saturating_add(capacity.saturating_mul(FILENAME_HASH_SLOT_BYTES));
        if capacity == 0 || !capacity.is_power_of_two() || index.metadata()?.len() != expected {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid dataset directory length",
            ));
        }

        let mut data = File::open(data_path)?;
        data.read_exact(&mut magic)?;
        if &magic != DATASET_TO_CID_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid dataset-to-CID payload",
            ));
        }
        let mut flags = [0u8; 1];
        data.read_exact(&mut flags)?;
        if flags[0]
            & !(DATASET_TO_CID_FLAG_ABUNDANCE
                | DATASET_TO_CID_FLAG_ADAPTIVE
                | DATASET_TO_CID_FLAG_MAJORITY)
            != 0
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "unsupported dataset-to-CID flags",
            ));
        }
        let abundance_log_base = if flags[0] & DATASET_TO_CID_FLAG_ABUNDANCE != 0 {
            let mut base = [0u8; 8];
            data.read_exact(&mut base)?;
            let base = f64::from_le_bytes(base);
            if !base.is_finite() || base <= 1.0 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid abundance log base",
                ));
            }
            Some(base)
        } else {
            None
        };
        let adaptive = flags[0] & DATASET_TO_CID_FLAG_ADAPTIVE != 0;
        let cid_count = if adaptive {
            let mut count = [0u8; 8];
            data.read_exact(&mut count)?;
            Some(u64::from_le_bytes(count))
        } else {
            None
        };
        let majority_frame = if flags[0] & DATASET_TO_CID_FLAG_MAJORITY != 0 {
            if !adaptive || abundance_log_base.is_some() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid majority dataset posting flags",
                ));
            }
            let mut len = [0u8; 8];
            data.read_exact(&mut len)?;
            let len = u64::from_le_bytes(len);
            Some((data.stream_position()?, len))
        } else {
            None
        };
        Ok(Self {
            index,
            data_path: data_path.to_path_buf(),
            dataset_count,
            abundance_log_base,
            adaptive,
            cid_count,
            majority_frame,
        })
    }

    fn posting_reader(&self, dataset_id: u32) -> Result<DatasetPostingReader> {
        if u64::from(dataset_id) >= self.dataset_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "dataset ID outside posting index",
            ));
        }
        let mut offset_bytes = [0u8; 16];
        self.index.read_exact_at(
            &mut offset_bytes,
            FILENAME_INDEX_HEADER_BYTES + u64::from(dataset_id) * 8,
        )?;
        let offset = u64::from_le_bytes(offset_bytes[..8].try_into().unwrap());
        let end = u64::from_le_bytes(offset_bytes[8..].try_into().unwrap());
        let compressed_len = end.checked_sub(offset).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "dataset posting offsets are not monotonic",
            )
        })?;
        let mut data = File::open(&self.data_path)?;
        data.seek(std::io::SeekFrom::Start(offset))?;
        let (codec, payload_offset, payload_len) = if self.adaptive {
            if compressed_len < 1 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "truncated adaptive dataset posting frame",
                ));
            }
            let mut codec = [0u8; 1];
            data.read_exact(&mut codec)?;
            (codec[0], offset + 1, compressed_len - 1)
        } else {
            (DATASET_POSTING_DELTA, offset, compressed_len)
        };
        if codec == DATASET_POSTING_BITMAP || codec == DATASET_POSTING_MAJORITY_XOR {
            if self.abundance_log_base.is_some() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "bitmap postings cannot contain abundance values",
                ));
            }
            let decoder = Decoder::new(data.take(payload_len))?.single_frame();
            let reference_decoder = if codec == DATASET_POSTING_MAJORITY_XOR {
                let (reference_offset, reference_len) = self.majority_frame.ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "majority-XOR posting has no reference bitmap",
                    )
                })?;
                let mut reference = File::open(&self.data_path)?;
                reference.seek(std::io::SeekFrom::Start(reference_offset))?;
                Some(Decoder::new(reference.take(reference_len))?.single_frame())
            } else {
                None
            };
            return Ok(DatasetPostingReader {
                encoding: PostingEncoding::Bitmap {
                    decoder,
                    reference_decoder,
                    current_byte: 0,
                    bit: 8,
                    next_cid: 0,
                    cid_count: self.cid_count.ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidData, "bitmap CID universe is missing")
                    })?,
                },
                abundance_decoder: None,
            });
        }
        if codec != DATASET_POSTING_DELTA {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "unsupported dataset posting codec",
            ));
        }
        let (decoder, abundance_decoder) = if self.abundance_log_base.is_some() {
            if payload_len < 8 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "truncated abundance posting frame",
                ));
            }
            let mut raw_len = [0u8; 8];
            data.read_exact(&mut raw_len)?;
            let cid_len = u64::from_le_bytes(raw_len);
            let abundance_len = payload_len
                .checked_sub(8)
                .and_then(|len| len.checked_sub(cid_len))
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "invalid abundance posting lengths",
                    )
                })?;
            let decoder = Decoder::new(data.take(cid_len))?.single_frame();
            let mut abundance_file = File::open(&self.data_path)?;
            abundance_file.seek(std::io::SeekFrom::Start(payload_offset + 8 + cid_len))?;
            let abundance_decoder =
                Decoder::new(abundance_file.take(abundance_len))?.single_frame();
            (decoder, Some(abundance_decoder))
        } else {
            (Decoder::new(data.take(payload_len))?.single_frame(), None)
        };
        Ok(DatasetPostingReader {
            encoding: PostingEncoding::Delta { decoder, cid: 0 },
            abundance_decoder,
        })
    }
}

enum PostingEncoding {
    Delta {
        decoder: PostingDecoder,
        cid: u64,
    },
    Bitmap {
        decoder: PostingDecoder,
        reference_decoder: Option<PostingDecoder>,
        current_byte: u8,
        bit: u8,
        next_cid: u64,
        cid_count: u64,
    },
}

struct DatasetPostingReader {
    encoding: PostingEncoding,
    abundance_decoder: Option<PostingDecoder>,
}

impl DatasetPostingReader {
    fn next(&mut self) -> Result<Option<(u64, Option<u8>)>> {
        let cid = match &mut self.encoding {
            PostingEncoding::Delta { decoder, cid } => {
                let Some(delta) = read_varint_u64_from_reader(decoder)? else {
                    return Ok(None);
                };
                *cid = cid.checked_add(delta).ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "dataset CID overflow")
                })?;
                *cid
            }
            PostingEncoding::Bitmap {
                decoder,
                reference_decoder,
                current_byte,
                bit,
                next_cid,
                cid_count,
            } => loop {
                if *next_cid >= *cid_count {
                    return Ok(None);
                }
                if *bit == 8 {
                    let mut byte = [0u8; 1];
                    decoder.read_exact(&mut byte)?;
                    if let Some(reference) = reference_decoder.as_mut() {
                        let mut base = [0u8; 1];
                        reference.read_exact(&mut base)?;
                        byte[0] ^= base[0];
                    }
                    *current_byte = byte[0];
                    *bit = 0;
                }
                let cid = *next_cid;
                let present = (*current_byte & (1u8 << *bit)) != 0;
                *bit += 1;
                *next_cid += 1;
                if present {
                    break cid;
                }
            },
        };
        let abundance = if let Some(decoder) = self.abundance_decoder.as_mut() {
            let mut deviation = [0u8; 1];
            decoder.read_exact(&mut deviation)?;
            Some(deviation[0])
        } else {
            None
        };
        Ok(Some((cid, abundance)))
    }
}

struct PositionLookup {
    positions: File,
    index: File,
    data_start: u64,
    entries: u64,
    stride: u64,
    checkpoints: u64,
    cached_checkpoint: u64,
    cached_positions: Vec<u64>,
}

impl PositionLookup {
    fn open(positions_path: &Path, index_path: &Path) -> Result<Self> {
        let mut positions = File::open(positions_path)?;
        let mut magic = [0u8; 4];
        positions.read_exact(&mut magic)?;
        if &magic != POSITIONS_COMPACT_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "indexed positions require KPS3",
            ));
        }
        let entries = read_varint_u64_from_reader(&mut positions)?.ok_or_else(|| {
            io::Error::new(io::ErrorKind::UnexpectedEof, "missing positions count")
        })?;
        let data_start = positions.stream_position()?;

        let mut index = File::open(index_path)?;
        index.read_exact(&mut magic)?;
        if &magic != POSITIONS_INDEX_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid positions index",
            ));
        }
        let mut stride = [0u8; 4];
        index.read_exact(&mut stride)?;
        let stride = u32::from_le_bytes(stride) as u64;
        let mut value = [0u8; 8];
        index.read_exact(&mut value)?;
        if u64::from_le_bytes(value) != entries || stride == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "positions index does not match payload",
            ));
        }
        index.read_exact(&mut value)?;
        let checkpoints = u64::from_le_bytes(value);
        if index.metadata()?.len() != POSITION_INDEX_HEADER_BYTES + checkpoints.saturating_mul(16) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid positions index length",
            ));
        }
        Ok(Self {
            positions,
            index,
            data_start,
            entries,
            stride,
            checkpoints,
            cached_checkpoint: u64::MAX,
            cached_positions: Vec::new(),
        })
    }

    fn load_checkpoint(&mut self, checkpoint: u64) -> Result<()> {
        if self.cached_checkpoint == checkpoint {
            return Ok(());
        }
        if checkpoint >= self.checkpoints {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "position checkpoint outside index",
            ));
        }
        let mut record = [0u8; 16];
        self.index
            .read_exact_at(&mut record, POSITION_INDEX_HEADER_BYTES + checkpoint * 16)?;
        let relative_after = u64::from_le_bytes(record[..8].try_into().unwrap());
        let mut position = u64::from_le_bytes(record[8..].try_into().unwrap());
        let first = checkpoint * self.stride;
        let end = (first + self.stride).min(self.entries);
        self.cached_positions.clear();
        self.cached_positions.push(position);
        self.positions.seek(std::io::SeekFrom::Start(
            self.data_start
                .checked_add(relative_after)
                .ok_or_else(|| io::Error::other("position checkpoint offset overflow"))?,
        ))?;
        for _ in first + 1..end {
            let delta = read_varint_u64_from_reader(&mut self.positions)?.ok_or_else(|| {
                io::Error::new(io::ErrorKind::UnexpectedEof, "truncated indexed positions")
            })?;
            position = position.checked_add(delta).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "indexed position overflow")
            })?;
            self.cached_positions.push(position);
        }
        self.cached_checkpoint = checkpoint;
        Ok(())
    }

    fn position(&mut self, entry: u64) -> Result<u64> {
        if entry >= self.entries {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "position entry outside index",
            ));
        }
        let checkpoint = entry / self.stride;
        self.load_checkpoint(checkpoint)?;
        Ok(self.cached_positions[(entry % self.stride) as usize])
    }

    fn range(&mut self, cid: u64) -> Result<(u64, u64)> {
        let start = self.position(cid)?;
        let end = self.position(cid + 1)?;
        Ok((start, end))
    }
}

#[derive(Clone, Copy)]
struct SizeBlockDescriptor {
    first_group: u64,
    group_count: u32,
    data_offset: u64,
    data_len: u64,
}

struct SizeLookup {
    data: File,
    index: File,
    total_groups: u64,
    block_count: u64,
    stride: u64,
    routing_offset: u64,
    cached_block: u64,
    cached_ranges: Vec<(usize, usize)>,
    cached_data: Vec<u8>,
}

impl SizeLookup {
    fn open(data_path: &Path, index_path: &Path) -> Result<Self> {
        let mut index = File::open(index_path)?;
        let mut magic = [0u8; 4];
        index.read_exact(&mut magic)?;
        if &magic != BUCKET_SIZES_INDEX_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid bucket-size index",
            ));
        }
        let mut value = [0u8; 8];
        index.read_exact(&mut value)?;
        let total_groups = u64::from_le_bytes(value);
        index.read_exact(&mut value)?;
        let block_count = u64::from_le_bytes(value);
        let mut stride = [0u8; 4];
        index.read_exact(&mut stride)?;
        let stride = u32::from_le_bytes(stride) as u64;
        let mut reserved = [0u8; 4];
        index.read_exact(&mut reserved)?;
        if stride == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "zero bucket-size index stride",
            ));
        }
        let routing_offset = SIZE_INDEX_HEADER_BYTES
            .checked_add(block_count.saturating_mul(SIZE_INDEX_RECORD_BYTES))
            .ok_or_else(|| io::Error::other("bucket-size routing offset overflow"))?;
        index.read_exact_at(&mut value, routing_offset)?;
        let routing_count = u64::from_le_bytes(value);
        if routing_count != total_groups.div_ceil(stride)
            || index.metadata()?.len() != routing_offset + 8 + routing_count.saturating_mul(4)
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid bucket-size routing table",
            ));
        }
        Ok(Self {
            data: File::open(data_path)?,
            index,
            total_groups,
            block_count,
            stride,
            routing_offset,
            cached_block: u64::MAX,
            cached_ranges: Vec::new(),
            cached_data: Vec::new(),
        })
    }

    fn descriptor(&self, block: u64) -> Result<SizeBlockDescriptor> {
        if block >= self.block_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "size block outside index",
            ));
        }
        let mut bytes = [0u8; SIZE_INDEX_RECORD_BYTES as usize];
        self.index.read_exact_at(
            &mut bytes,
            SIZE_INDEX_HEADER_BYTES + block * SIZE_INDEX_RECORD_BYTES,
        )?;
        Ok(SizeBlockDescriptor {
            first_group: u64::from_le_bytes(bytes[..8].try_into().unwrap()),
            group_count: u32::from_le_bytes(bytes[8..12].try_into().unwrap()),
            data_offset: u64::from_le_bytes(bytes[16..24].try_into().unwrap()),
            data_len: u64::from_le_bytes(bytes[24..32].try_into().unwrap()),
        })
    }

    fn block_for_group(&self, group: u64) -> Result<u64> {
        if group >= self.total_groups {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "size group outside index",
            ));
        }
        let route = group / self.stride;
        let mut raw = [0u8; 4];
        self.index
            .read_exact_at(&mut raw, self.routing_offset + 8 + route * 4)?;
        let mut block = u32::from_le_bytes(raw) as u64;
        loop {
            let descriptor = self.descriptor(block)?;
            if group >= descriptor.first_group
                && group < descriptor.first_group + u64::from(descriptor.group_count)
            {
                return Ok(block);
            }
            block += 1;
        }
    }

    fn load_block(&mut self, block: u64) -> Result<()> {
        if self.cached_block == block {
            return Ok(());
        }
        let descriptor = self.descriptor(block)?;
        let mut compressed = vec![
            0u8;
            usize::try_from(descriptor.data_len).map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "compressed size block is too large",
                )
            })?
        ];
        self.data
            .read_exact_at(&mut compressed, descriptor.data_offset)?;
        self.cached_data = zstd::decode_all(compressed.as_slice())?;
        self.cached_ranges.clear();
        let mut cursor = std::io::Cursor::new(self.cached_data.as_slice());
        for _ in 0..descriptor.group_count {
            let length = read_varint_u64_from_reader(&mut cursor)?.ok_or_else(|| {
                io::Error::new(io::ErrorKind::UnexpectedEof, "missing size-group length")
            })?;
            let start = cursor.position() as usize;
            let end = start
                .checked_add(usize::try_from(length).map_err(|_| {
                    io::Error::new(io::ErrorKind::InvalidData, "size-group length is too large")
                })?)
                .ok_or_else(|| io::Error::other("size-group offset overflow"))?;
            if end > self.cached_data.len() {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "truncated size group",
                ));
            }
            self.cached_ranges.push((start, end));
            cursor.set_position(end as u64);
        }
        if cursor.position() as usize != self.cached_data.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "trailing size-block data",
            ));
        }
        self.cached_block = block;
        Ok(())
    }

    fn group_into(&mut self, group: u64, sizes: &mut Vec<usize>) -> Result<()> {
        let block = self.block_for_group(group)?;
        self.load_block(block)?;
        let descriptor = self.descriptor(block)?;
        let local = usize::try_from(group - descriptor.first_group).unwrap();
        let (start, end) = self.cached_ranges[local];
        decode_varint_deltas_to_sizes_into(
            &self.cached_data[start..end],
            "indexed size group",
            sizes,
        )
    }
}

/// High-level decompression entry point.
pub fn decompress(
    size_filename: &String,
    color_id_filename: &String,
    tigs_filename: &String,
    positions_filename: &String,
    filename_id: &String,
    out_dir: &String,
    wanted_files_path: &String,
    input_dir: String,
) -> std::io::Result<()> {
    decompress_with_options(
        size_filename,
        color_id_filename,
        tigs_filename,
        positions_filename,
        filename_id,
        out_dir,
        wanted_files_path,
        input_dir,
        GgcatRebuildConfig::default(),
    )
}

#[derive(Clone, Debug)]
pub struct GgcatRebuildConfig {
    pub enabled: bool,
    pub threads: usize,
    pub memory_gb: usize,
    pub k: usize,
    pub temp_dir: String,
    pub use_unitigs: bool,
    pub use_matchtigs: bool,
    pub use_eulertigs: bool,
    pub restore_abundance: bool,
    pub output_compression: OutputCompression,
    pub color_set_operation: Option<ColorSetOperation>,
}

impl Default for GgcatRebuildConfig {
    fn default() -> Self {
        Self {
            enabled: false,
            threads: 1,
            memory_gb: 8,
            k: 31,
            temp_dir: String::new(),
            use_unitigs: false,
            use_matchtigs: false,
            use_eulertigs: false,
            restore_abundance: false,
            output_compression: OutputCompression::Zstd,
            color_set_operation: None,
        }
    }
}

#[derive(clap::ValueEnum, Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum OutputCompression {
    #[value(alias = "uncompressed", alias = "fa")]
    Fasta,
    #[value(alias = "gzip")]
    Gz,
    #[default]
    #[value(alias = "zst")]
    Zstd,
    Xz,
}

impl std::fmt::Display for OutputCompression {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Fasta => "fasta",
            Self::Gz => "gz",
            Self::Zstd => "zstd",
            Self::Xz => "xz",
        })
    }
}

impl OutputCompression {
    fn fasta_suffix(self) -> &'static str {
        match self {
            Self::Fasta => ".fa",
            Self::Gz => ".fa.gz",
            Self::Zstd => ".fa.zst",
            Self::Xz => ".fa.xz",
        }
    }
}

#[derive(clap::ValueEnum, Clone, Copy, Debug, Eq, PartialEq)]
pub enum ColorSetOperation {
    Union,
    Intersection,
}

impl std::fmt::Display for ColorSetOperation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Union => "union",
            Self::Intersection => "intersection",
        })
    }
}

enum FastaOutputWriter {
    Fasta(BufWriter<File>),
    Gz(GzEncoder<BufWriter<File>>),
    Zstd(ZstdEncoder<'static, BufWriter<File>>),
    Xz(XzEncoder<BufWriter<File>>),
}

impl FastaOutputWriter {
    fn create(path: &Path, compression: OutputCompression) -> Result<Self> {
        let writer = BufWriter::new(File::create(path)?);
        match compression {
            OutputCompression::Fasta => Ok(Self::Fasta(writer)),
            OutputCompression::Gz => Ok(Self::Gz(GzEncoder::new(writer, Compression::default()))),
            OutputCompression::Zstd => Ok(Self::Zstd(ZstdEncoder::new(writer, 0)?)),
            OutputCompression::Xz => Ok(Self::Xz(XzEncoder::new(writer, 6))),
        }
    }

    fn finish(self) -> Result<()> {
        match self {
            Self::Fasta(mut writer) => writer.flush(),
            Self::Gz(writer) => writer.finish().and_then(|mut writer| writer.flush()),
            Self::Zstd(writer) => writer.finish().and_then(|mut writer| writer.flush()),
            Self::Xz(writer) => writer.finish().and_then(|mut writer| writer.flush()),
        }
    }
}

impl Write for FastaOutputWriter {
    fn write(&mut self, buf: &[u8]) -> Result<usize> {
        match self {
            Self::Fasta(writer) => writer.write(buf),
            Self::Gz(writer) => writer.write(buf),
            Self::Zstd(writer) => writer.write(buf),
            Self::Xz(writer) => writer.write(buf),
        }
    }

    fn flush(&mut self) -> Result<()> {
        match self {
            Self::Fasta(writer) => writer.flush(),
            Self::Gz(writer) => writer.flush(),
            Self::Zstd(writer) => writer.flush(),
            Self::Xz(writer) => writer.flush(),
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum RebuildMode {
    Simplitig,
    Unitig,
    Matchtig,
    Eulertig,
}

fn resolve_rebuild_mode(
    use_unitigs: bool,
    use_matchtigs: bool,
    use_eulertigs: bool,
) -> RebuildMode {
    if use_unitigs {
        RebuildMode::Unitig
    } else if use_matchtigs {
        RebuildMode::Matchtig
    } else if use_eulertigs {
        RebuildMode::Eulertig
    } else {
        RebuildMode::Simplitig
    }
}

fn ggcat_extra_elaboration(mode: RebuildMode) -> ExtraElaboration {
    match mode {
        RebuildMode::Unitig => ExtraElaboration::None,
        RebuildMode::Matchtig => ExtraElaboration::GreedyMatchtigs,
        RebuildMode::Eulertig => ExtraElaboration::FastEulertigs,
        RebuildMode::Simplitig => ExtraElaboration::FastSimplitigs,
    }
}

fn rebuild_mode_name(mode: RebuildMode) -> &'static str {
    match mode {
        RebuildMode::Simplitig => "simplitigs",
        RebuildMode::Unitig => "unitigs",
        RebuildMode::Matchtig => "matchtigs",
        RebuildMode::Eulertig => "eulertigs",
    }
}

fn normalize_out_dir(out_dir: &str) -> PathBuf {
    if out_dir.is_empty() {
        PathBuf::from(".")
    } else {
        PathBuf::from(out_dir)
    }
}

fn ensure_output_dir(out_dir: &str) -> io::Result<PathBuf> {
    let path = normalize_out_dir(out_dir);
    fs::create_dir_all(&path)?;
    Ok(path)
}

fn dump_output_path(out_dir: &str, source_path: &str, compression: OutputCompression) -> PathBuf {
    let out_dir_path = normalize_out_dir(out_dir);
    let trunc_filename = Path::new(source_path).file_stem().unwrap_or_default();
    out_dir_path.join(format!(
        "Dump_{}{}",
        trunc_filename.to_str().unwrap_or("unknown"),
        compression.fasta_suffix()
    ))
}

fn color_set_output_path(
    out_dir: &str,
    operation: ColorSetOperation,
    compression: OutputCompression,
) -> PathBuf {
    normalize_out_dir(out_dir).join(format!("Dump_{}{}", operation, compression.fasta_suffix()))
}

fn collect_dump_fastas(
    out_dir: &str,
    compression: OutputCompression,
) -> std::io::Result<Vec<PathBuf>> {
    let out_dir_path = normalize_out_dir(out_dir);
    let suffix = compression.fasta_suffix();
    let mut dump_files = Vec::new();
    for entry_result in fs::read_dir(&out_dir_path)? {
        let entry = entry_result?;
        let path = entry.path();
        let is_dump = path
            .file_name()
            .and_then(|name| name.to_str())
            .map(|name| name.starts_with("Dump_") && name.ends_with(suffix))
            .unwrap_or(false);
        if is_dump {
            dump_files.push(path);
        }
    }
    dump_files.sort();
    Ok(dump_files)
}

fn run_ggcat_rebuild(out_dir: &str, cfg: &GgcatRebuildConfig) -> std::io::Result<()> {
    if cfg.restore_abundance {
        eprintln!(
            "Warning: decompressed Dump_* files retain abundance headers; the combined GGCAT rebuild output does not preserve per-dataset abundance."
        );
    }
    let mode = resolve_rebuild_mode(cfg.use_unitigs, cfg.use_matchtigs, cfg.use_eulertigs);
    let dump_fastas = collect_dump_fastas(out_dir, cfg.output_compression)?;
    if dump_fastas.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            "no decompressed Dump_* files found to rebuild with ggcat",
        ));
    }

    let out_dir_path = ensure_output_dir(out_dir)?;

    let rebuilt_output = out_dir_path.join(format!(
        "rebuilt_{}{}",
        rebuild_mode_name(mode),
        cfg.output_compression.fasta_suffix()
    ));
    let temp_dir = if cfg.temp_dir.is_empty() {
        out_dir_path.join("ggcat_rebuild_tmp")
    } else {
        PathBuf::from(&cfg.temp_dir)
    };
    fs::create_dir_all(&temp_dir)?;

    // GGCAT supports plain, gzip, and zstd FASTA streams, but not xz. Keep xz
    // available to users by staging plain FASTA inputs and recompressing the
    // rebuilt graph. The temporary directory is removed automatically.
    let xz_bridge = if cfg.output_compression == OutputCompression::Xz {
        Some(
            tempfile::Builder::new()
                .prefix("kloe-xz-rebuild-")
                .tempdir_in(&temp_dir)?,
        )
    } else {
        None
    };
    let rebuild_inputs = if let Some(bridge) = xz_bridge.as_ref() {
        let mut inputs = Vec::with_capacity(dump_fastas.len());
        for (index, dump) in dump_fastas.iter().enumerate() {
            let staged = bridge.path().join(format!("input_{index}.fa"));
            let mut reader = xz2::read::XzDecoder::new(BufReader::new(File::open(dump)?));
            let mut writer = BufWriter::new(File::create(&staged)?);
            io::copy(&mut reader, &mut writer)?;
            writer.flush()?;
            inputs.push(staged);
        }
        inputs
    } else {
        dump_fastas
    };
    let ggcat_output = xz_bridge
        .as_ref()
        .map(|bridge| {
            bridge
                .path()
                .join(format!("rebuilt_{}.fa", rebuild_mode_name(mode)))
        })
        .unwrap_or_else(|| rebuilt_output.clone());

    let ggcat_memory_gb = cfg.memory_gb.max(1);
    let instance = GGCATInstance::create(GGCATConfig {
        temp_dir: Some(temp_dir.clone()),
        memory: ggcat_memory_gb as f64,
        prefer_memory: false,
        total_threads_count: cfg.threads.max(1),
        intermediate_compression_level: None,
        stats_file: None,
        messages_callback: None,
    })
    .map_err(|err| io::Error::other(format!("create ggcat instance: {err}")))?;

    let streams = rebuild_inputs
        .iter()
        .map(|dump| {
            let resolved = fs::canonicalize(dump).unwrap_or_else(|_| dump.clone());
            GeneralSequenceBlockData::FASTA((resolved, None))
        })
        .collect::<Vec<_>>();

    println!(
        "Running embedded ggcat rebuild: output={}, mode={}, k={}, threads={}, memory={}GB",
        rebuilt_output.display(),
        rebuild_mode_name(mode),
        cfg.k,
        cfg.threads.max(1),
        ggcat_memory_gb
    );
    let graph_path = instance
        .build_graph(
            streams,
            ggcat_output,
            None,
            None,
            cfg.k,
            cfg.threads.max(1),
            false,
            None,
            false,
            1,
            ggcat_extra_elaboration(mode),
            None,
        )
        .map_err(|err| io::Error::other(format!("run ggcat rebuild: {err}")))?;

    if xz_bridge.is_some() {
        let mut reader = BufReader::new(File::open(&graph_path)?);
        let mut writer = FastaOutputWriter::create(&rebuilt_output, OutputCompression::Xz)?;
        io::copy(&mut reader, &mut writer)?;
        writer.finish()?;
    }

    println!(
        "ggcat rebuild complete: {} (mode={})",
        rebuilt_output.display(),
        rebuild_mode_name(mode)
    );
    Ok(())
}

pub fn decompress_with_options(
    size_filename: &String,
    color_id_filename: &String,
    tigs_filename: &String,
    positions_filename: &String,
    filename_id: &String,
    out_dir: &String,
    wanted_files_path: &String,
    input_dir: String,
    ggcat_cfg: GgcatRebuildConfig,
) -> std::io::Result<()> {
    println!("Writing decompressed data in {out_dir}");
    ensure_output_dir(out_dir)?;

    let archive_root = Path::new(&input_dir);
    let dataset_data_path = archive_root.join(DATASET_TO_CID_FILE);
    let dataset_index_path = archive_root.join(FILENAME_INDEX_FILE);
    let filename_index_path = archive_root.join(FILENAME_INDEX_FILE);
    let positions_index_path = archive_root.join(POSITIONS_INDEX_FILE);
    let sizes_index_path = archive_root.join(BUCKET_SIZES_INDEX_FILE);
    let tigs_index_path = archive_root.join(TIGS_INDEX_FILE);
    if let Some(operation) = ggcat_cfg.color_set_operation {
        if wanted_files_path.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "a color-set union or intersection requires a non-empty wanted-files list",
            ));
        }
        if ggcat_cfg.restore_abundance {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "abundance restoration is not defined for a combined color-set output",
            ));
        }
        if dataset_data_path.is_file()
            && dataset_index_path.is_file()
            && filename_index_path.is_file()
            && positions_index_path.is_file()
            && sizes_index_path.is_file()
            && tigs_index_path.is_file()
        {
            decompress_indexed_color_set(
                archive_root,
                &(input_dir.clone() + filename_id),
                &(input_dir.clone() + positions_filename),
                &(input_dir.clone() + size_filename),
                &(input_dir.clone() + tigs_filename),
                wanted_files_path,
                out_dir,
                operation,
                ggcat_cfg.output_compression,
            )?;
        } else {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "color-set union/intersection requires an indexed KLOE archive",
            ));
        }
        if ggcat_cfg.enabled {
            run_ggcat_rebuild(out_dir, &ggcat_cfg)?;
        }
        return Ok(());
    }
    if !wanted_files_path.is_empty()
        && dataset_data_path.is_file()
        && dataset_index_path.is_file()
        && filename_index_path.is_file()
        && positions_index_path.is_file()
        && sizes_index_path.is_file()
        && tigs_index_path.is_file()
    {
        decompress_indexed_targeted(
            archive_root,
            &(input_dir.clone() + filename_id),
            &(input_dir.clone() + positions_filename),
            &(input_dir.clone() + size_filename),
            &(input_dir.clone() + tigs_filename),
            wanted_files_path,
            out_dir,
            ggcat_cfg.restore_abundance,
            ggcat_cfg.output_compression,
        )?;
        if ggcat_cfg.enabled {
            run_ggcat_rebuild(out_dir, &ggcat_cfg)?;
        }
        return Ok(());
    }
    if wanted_files_path.is_empty()
        && dataset_data_path.is_file()
        && dataset_index_path.is_file()
        && positions_index_path.is_file()
        && sizes_index_path.is_file()
        && tigs_index_path.is_file()
    {
        let filenames = load_archive_filenames(&(input_dir.clone() + filename_id))?;
        for batch in filenames.chunks(128) {
            decompress_indexed_selected_batch(
                archive_root,
                &(input_dir.clone() + positions_filename),
                &(input_dir.clone() + size_filename),
                &(input_dir.clone() + tigs_filename),
                batch,
                out_dir,
                ggcat_cfg.restore_abundance,
                ggcat_cfg.output_compression,
            )?;
        }
        if ggcat_cfg.enabled {
            run_ggcat_rebuild(out_dir, &ggcat_cfg)?;
        }
        return Ok(());
    }

    let filenames_path = input_dir.clone() + filename_id;
    let filenames = load_archive_filenames(&filenames_path)?;
    let sidecar_path = Path::new(&input_dir).join(CID_TO_DATASET_FILE);
    if sidecar_path.is_file() {
        let wanted_filenames = if wanted_files_path.is_empty() {
            filenames.clone()
        } else {
            select_wanted_filenames(&filenames, wanted_files_path)?
        };
        println!(
            "Using compact CID-to-dataset index for {} decompression",
            if wanted_files_path.is_empty() {
                "complete"
            } else {
                "targeted"
            }
        );
        decompress_sidecar(
            &sidecar_path,
            &(input_dir.clone() + positions_filename),
            &(input_dir.to_owned() + tigs_filename),
            &(input_dir.to_owned() + size_filename),
            out_dir,
            &filenames,
            if wanted_files_path.is_empty() {
                None
            } else {
                Some(&wanted_filenames)
            },
            ggcat_cfg.threads,
            ggcat_cfg.restore_abundance,
            ggcat_cfg.output_compression,
        )?;
        if ggcat_cfg.enabled {
            run_ggcat_rebuild(out_dir, &ggcat_cfg)?;
        }
        return Ok(());
    }

    if ggcat_cfg.restore_abundance {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--abundance was requested but this legacy archive has no abundance metadata",
        ));
    }

    if wanted_files_path != "" {
        let input_file = File::open(input_dir.clone() + filename_id).unwrap();
        let input_reader = BufReader::new(input_file);
        let mut filenames_id_map = HashMap::new();
        let mut file_id: u32 = 0;
        for line_result in input_reader.lines() {
            let line = line_result?;
            if let Some((path, size)) = line.split_once(':') {
                filenames_id_map.insert(path.to_string(), (file_id, size.parse::<u64>().unwrap()));
                println!("File: {}, ID: {}, offset: {}", path, file_id, size);
            }
            file_id += 1;
        }

        let cid_map_out_filenames = match get_cid_to_id_targeted(
            &(input_dir.clone() + &color_id_filename),
            &filenames_id_map,
            &wanted_files_path,
        ) {
            Ok(map) => map,
            Err(e) => panic!("error gathering color set ids: {e:?}"),
        };
        let cid_to_id_map = cid_map_out_filenames.0;
        let wanted_filenames = cid_map_out_filenames.1;

        println!("Query file given, decompressing only subpart of archive....");
        decompress_wanted(
            &wanted_filenames,
            &(input_dir.clone() + positions_filename),
            cid_to_id_map,
            &(input_dir.to_owned() + tigs_filename),
            &(input_dir.to_owned() + size_filename),
            out_dir,
            ggcat_cfg.output_compression,
        )?;
    } else {
        println!("No query file given, decompressing entire archive....");
        let cid_to_id_map = match get_cid_to_id(&(input_dir.clone() + &color_id_filename)) {
            Ok(map) => map,
            Err(e) => panic!("Error getting cid to id map {e:?}"),
        };
        let input_file = File::open(input_dir.clone() + filename_id).unwrap();
        let input_reader = BufReader::new(input_file);
        let mut filenames_id = Vec::new();
        let mut file_id: u32 = 0;
        for line_result in input_reader.lines() {
            let line = line_result?;
            if let Some((path, _)) = line.split_once(":") {
                filenames_id.push((path.to_owned(), file_id));
            }
            file_id += 1;
        }
        println!("{}", input_dir.to_owned() + size_filename);
        decompress_all(
            &(input_dir.to_owned() + size_filename),
            &(input_dir.clone() + positions_filename),
            &(input_dir.to_owned() + tigs_filename),
            out_dir,
            filenames_id,
            cid_to_id_map,
            ggcat_cfg.output_compression,
        )?;
    }
    if ggcat_cfg.enabled {
        run_ggcat_rebuild(out_dir, &ggcat_cfg)?;
    }
    Ok(())
}

fn indexed_wanted_datasets(
    filename_index_path: &Path,
    filename_table_path: &Path,
    wanted_files_path: &str,
) -> Result<Vec<(String, u32)>> {
    let filename_index = FilenameIndex::open(filename_index_path, filename_table_path)?;
    let wanted_reader = BufReader::new(File::open(wanted_files_path)?);
    let mut selected = Vec::new();
    let mut seen = HashSet::new();
    for line in wanted_reader.lines() {
        let name = line?;
        match filename_index.lookup(&name)? {
            Some(dataset_id) if seen.insert(dataset_id) => selected.push((name, dataset_id)),
            Some(_) => {}
            None => println!(
                "FILE {} NOT FOUND IN ARCHIVE, CHECK SPELLING OR ACTUAL PRESENCE IN ARCHIVE",
                name
            ),
        }
    }
    Ok(selected)
}

fn indexed_required_datasets(
    filename_index_path: &Path,
    filename_table_path: &Path,
    wanted_files_path: &str,
) -> Result<Vec<(String, u32)>> {
    let filename_index = FilenameIndex::open(filename_index_path, filename_table_path)?;
    let wanted_reader = BufReader::new(File::open(wanted_files_path)?);
    let mut selected = Vec::new();
    let mut seen = HashSet::new();
    for line in wanted_reader.lines() {
        let name = line?;
        let name = name.trim();
        if name.is_empty() {
            continue;
        }
        let dataset_id = filename_index.lookup(name)?.ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("selected color '{name}' was not found in the archive"),
            )
        })?;
        if seen.insert(dataset_id) {
            selected.push((name.to_owned(), dataset_id));
        }
    }
    if selected.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "the wanted-files list contains no colors",
        ));
    }
    Ok(selected)
}

#[allow(clippy::too_many_arguments)]
fn decompress_indexed_color_set(
    archive_root: &Path,
    filename_table_path: &str,
    positions_path: &str,
    sizes_path: &str,
    tigs_path: &str,
    wanted_files_path: &str,
    out_dir: &str,
    operation: ColorSetOperation,
    output_compression: OutputCompression,
) -> Result<()> {
    let selected = indexed_required_datasets(
        &archive_root.join(FILENAME_INDEX_FILE),
        Path::new(filename_table_path),
        wanted_files_path,
    )?;
    let postings = DatasetPostingIndex::open(
        &archive_root.join(FILENAME_INDEX_FILE),
        &archive_root.join(DATASET_TO_CID_FILE),
    )?;
    let mut readers = selected
        .iter()
        .map(|(_, dataset_id)| postings.posting_reader(*dataset_id))
        .collect::<Result<Vec<_>>>()?;
    let mut heap = BinaryHeap::<Reverse<(u64, usize, u8)>>::new();
    for (reader_index, reader) in readers.iter_mut().enumerate() {
        if let Some((cid, abundance)) = reader.next()? {
            heap.push(Reverse((cid, reader_index, abundance.unwrap_or(0))));
        }
    }

    let mut positions = PositionLookup::open(
        Path::new(positions_path),
        &archive_root.join(POSITIONS_INDEX_FILE),
    )?;
    let mut sizes = SizeLookup::open(
        Path::new(sizes_path),
        &archive_root.join(BUCKET_SIZES_INDEX_FILE),
    )?;
    let mut tigs =
        PackedTigsReader::open_indexed(Path::new(tigs_path), archive_root.join(TIGS_INDEX_FILE))?;
    let output_path = color_set_output_path(out_dir, operation, output_compression);
    let mut output = FastaOutputWriter::create(&output_path, output_compression)?;
    let mut lengths = Vec::<usize>::new();
    let mut selected_cids = 0u64;
    let mut total_tigs = 0u64;

    while let Some(Reverse((cid, reader_index, _))) = heap.pop() {
        let mut membership_count = 1usize;
        if let Some((next_cid, abundance)) = readers[reader_index].next()? {
            heap.push(Reverse((next_cid, reader_index, abundance.unwrap_or(0))));
        }
        while heap
            .peek()
            .is_some_and(|Reverse((next_cid, _, _))| *next_cid == cid)
        {
            let Reverse((_, next_reader, _)) = heap.pop().unwrap();
            membership_count += 1;
            if let Some((next_cid, abundance)) = readers[next_reader].next()? {
                heap.push(Reverse((next_cid, next_reader, abundance.unwrap_or(0))));
            }
        }
        let include = match operation {
            ColorSetOperation::Union => true,
            ColorSetOperation::Intersection => membership_count == selected.len(),
        };
        if !include {
            continue;
        }

        sizes.group_into(cid, &mut lengths)?;
        let (packed_start, packed_end) = positions.range(cid)?;
        let expected_packed = lengths
            .iter()
            .try_fold(0u64, |total, length| {
                total.checked_add(length.div_ceil(4) as u64)
            })
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "CID size overflow"))?;
        if packed_start.checked_add(expected_packed) != Some(packed_end) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("CID {cid} lengths do not match its packed interval"),
            ));
        }

        let mut packed_position = packed_start;
        let mut first = 0usize;
        while first < lengths.len() {
            let mut end = first;
            let mut output_bytes = 0usize;
            let mut packed_bytes = 0usize;
            while end < lengths.len() {
                let record_bytes = lengths[end].saturating_add(3);
                if end > first && output_bytes.saturating_add(record_bytes) > DECOMPRESS_BATCH_BASES
                {
                    break;
                }
                output_bytes = output_bytes.saturating_add(record_bytes);
                packed_bytes = packed_bytes
                    .checked_add(lengths[end].div_ceil(4))
                    .ok_or_else(|| io::Error::other("color-set packed batch overflow"))?;
                end += 1;
            }
            let mut encoded = vec![0u8; packed_bytes];
            tigs.read_exact_at(packed_position, &mut encoded)?;
            packed_position += packed_bytes as u64;
            let decoded = decode_fasta_batch(DecodeJob {
                ordinal: 0,
                sizes: lengths[first..end].to_vec(),
                encoded,
                dataset_ids: Arc::new(vec![0]),
                abundance_codes: None,
                abundance_log_base: None,
            })?;
            output.write_all(decoded.fasta.as_slice())?;
            first = end;
        }
        if packed_position != packed_end {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "color-set CID ended at the wrong packed position",
            ));
        }
        selected_cids += 1;
        total_tigs = total_tigs.saturating_add(lengths.len() as u64);
    }
    output.finish()?;
    println!(
        "Color-set {} complete: colors={}, CIDs={}, tigs={}, output={}",
        operation,
        selected.len(),
        selected_cids,
        total_tigs,
        output_path.display()
    );
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn decompress_indexed_targeted(
    archive_root: &Path,
    filename_table_path: &str,
    positions_path: &str,
    sizes_path: &str,
    tigs_path: &str,
    wanted_files_path: &str,
    out_dir: &str,
    restore_abundance: bool,
    output_compression: OutputCompression,
) -> Result<()> {
    let selected = indexed_wanted_datasets(
        &archive_root.join(FILENAME_INDEX_FILE),
        Path::new(filename_table_path),
        wanted_files_path,
    )?;
    if selected.is_empty() {
        println!("No requested datasets were found in the archive");
        return Ok(());
    }

    for batch in selected.chunks(128) {
        decompress_indexed_selected_batch(
            archive_root,
            positions_path,
            sizes_path,
            tigs_path,
            batch,
            out_dir,
            restore_abundance,
            output_compression,
        )?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn decompress_indexed_selected_batch(
    archive_root: &Path,
    positions_path: &str,
    sizes_path: &str,
    tigs_path: &str,
    selected: &[(String, u32)],
    out_dir: &str,
    restore_abundance: bool,
    output_compression: OutputCompression,
) -> Result<()> {
    let postings = DatasetPostingIndex::open(
        &archive_root.join(FILENAME_INDEX_FILE),
        &archive_root.join(DATASET_TO_CID_FILE),
    )?;
    if restore_abundance && postings.abundance_log_base.is_none() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--abundance was requested but this archive has no abundance metadata",
        ));
    }
    let abundance_log_base = restore_abundance
        .then_some(postings.abundance_log_base)
        .flatten();

    let mut readers = selected
        .iter()
        .map(|(_, dataset_id)| postings.posting_reader(*dataset_id))
        .collect::<Result<Vec<_>>>()?;
    let mut heap = BinaryHeap::<Reverse<(u64, usize, u8)>>::new();
    for (reader_index, reader) in readers.iter_mut().enumerate() {
        if let Some((cid, abundance)) = reader.next()? {
            heap.push(Reverse((cid, reader_index, abundance.unwrap_or(0))));
        }
    }

    let mut positions = PositionLookup::open(
        Path::new(positions_path),
        &archive_root.join(POSITIONS_INDEX_FILE),
    )?;
    let mut sizes = SizeLookup::open(
        Path::new(sizes_path),
        &archive_root.join(BUCKET_SIZES_INDEX_FILE),
    )?;
    let mut tigs =
        PackedTigsReader::open_indexed(Path::new(tigs_path), archive_root.join(TIGS_INDEX_FILE))?;
    let mut abundance_bases = if abundance_log_base.is_some() {
        let reader = PackedTigsReader::open_indexed(
            archive_root.join(ABUNDANCE_BASE_FILE),
            archive_root.join(ABUNDANCE_BASE_INDEX_FILE),
        )?;
        if reader.logical_len() != positions.entries.saturating_sub(1) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "abundance base stream does not match archive CID count",
            ));
        }
        Some(reader)
    } else {
        None
    };
    let mut outputs = selected
        .iter()
        .map(|(name, _)| {
            FastaOutputWriter::create(
                &dump_output_path(out_dir, name, output_compression),
                output_compression,
            )
        })
        .collect::<Result<Vec<_>>>()?;

    let mut lengths = Vec::<usize>::new();
    let mut selected_cids = 0u64;
    let mut total_tigs = 0u64;
    while let Some(Reverse((cid, reader_index, abundance))) = heap.pop() {
        let mut targets = vec![(reader_index, abundance)];
        if let Some((next_cid, next_abundance)) = readers[reader_index].next()? {
            heap.push(Reverse((
                next_cid,
                reader_index,
                next_abundance.unwrap_or(0),
            )));
        }
        while heap
            .peek()
            .is_some_and(|Reverse((next_cid, _, _))| *next_cid == cid)
        {
            let Reverse((_, next_reader, next_code)) = heap.pop().unwrap();
            targets.push((next_reader, next_code));
            if let Some((next_cid, next_abundance)) = readers[next_reader].next()? {
                heap.push(Reverse((
                    next_cid,
                    next_reader,
                    next_abundance.unwrap_or(0),
                )));
            }
        }

        if let Some(reader) = abundance_bases.as_mut() {
            let mut base = [0u8; 1];
            reader.read_exact_at(cid, &mut base)?;
            for (_, code) in &mut targets {
                *code = base[0].wrapping_add(*code);
            }
        }

        sizes.group_into(cid, &mut lengths)?;
        let (packed_start, packed_end) = positions.range(cid)?;
        let expected_packed = lengths
            .iter()
            .try_fold(0u64, |total, length| {
                total.checked_add(length.div_ceil(4) as u64)
            })
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "CID packed length overflow")
            })?;
        if packed_start.checked_add(expected_packed) != Some(packed_end) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("CID {cid} lengths do not match its packed interval"),
            ));
        }

        let target_ids = Arc::new(
            targets
                .iter()
                .map(|(target, _)| *target as u32)
                .collect::<Vec<_>>(),
        );
        let target_codes = abundance_log_base
            .map(|_| Arc::new(targets.iter().map(|(_, code)| *code).collect::<Vec<_>>()));
        let mut packed_position = packed_start;
        let mut first = 0usize;
        while first < lengths.len() {
            let mut end = first;
            let mut output_bytes = 0usize;
            let mut packed_bytes = 0usize;
            while end < lengths.len() {
                let record_bytes = lengths[end].saturating_add(3);
                if end > first && output_bytes.saturating_add(record_bytes) > DECOMPRESS_BATCH_BASES
                {
                    break;
                }
                output_bytes = output_bytes.saturating_add(record_bytes);
                packed_bytes = packed_bytes
                    .checked_add(lengths[end].div_ceil(4))
                    .ok_or_else(|| io::Error::other("targeted packed batch overflow"))?;
                end += 1;
            }
            let mut encoded = vec![0u8; packed_bytes];
            tigs.read_exact_at(packed_position, &mut encoded)?;
            packed_position += packed_bytes as u64;
            let decoded = decode_fasta_batch(DecodeJob {
                ordinal: 0,
                sizes: lengths[first..end].to_vec(),
                encoded,
                dataset_ids: Arc::clone(&target_ids),
                abundance_codes: target_codes.as_ref().map(Arc::clone),
                abundance_log_base,
            })?;
            if let Some(sequences) = decoded.abundance_sequences.as_ref() {
                let mut rendered = HashMap::<u8, Vec<u8>>::new();
                for &(target, code) in &targets {
                    let fasta = if let Some(existing) = rendered.get(&code) {
                        existing
                    } else {
                        let abundance = crate::compress::decode_abundance(
                            code,
                            abundance_log_base.expect("abundance checked above"),
                        );
                        let fasta = build_abundance_fasta(sequences, abundance)?;
                        rendered.entry(code).or_insert(fasta)
                    };
                    outputs[target].write_all(fasta)?;
                }
            } else {
                for &(target, _) in &targets {
                    outputs[target].write_all(decoded.fasta.as_slice())?;
                }
            }
            first = end;
        }
        if packed_position != packed_end {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "targeted CID ended at the wrong packed position",
            ));
        }
        selected_cids += 1;
        total_tigs = total_tigs.saturating_add(lengths.len() as u64);
    }
    for output in outputs {
        output.finish()?;
    }
    println!(
        "Output-sensitive targeted decompression complete: datasets={}, CIDs={}, tigs={}",
        selected.len(),
        selected_cids,
        total_tigs
    );
    Ok(())
}

fn load_archive_filenames(filename_id_path: &str) -> Result<Vec<(String, u32)>> {
    let input_reader = BufReader::new(File::open(filename_id_path)?);
    let mut filenames = Vec::new();
    for (file_id, line_result) in input_reader.lines().enumerate() {
        let line = line_result?;
        let path = line.rsplit_once(':').map(|(path, _)| path).unwrap_or(&line);
        let file_id = u32::try_from(file_id).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "archive has more than u32::MAX datasets",
            )
        })?;
        filenames.push((path.to_owned(), file_id));
    }
    Ok(filenames)
}

fn select_wanted_filenames(
    filenames: &[(String, u32)],
    wanted_files_path: &str,
) -> Result<Vec<(String, u32)>> {
    let by_name: HashMap<&str, u32> = filenames
        .iter()
        .map(|(path, file_id)| (path.as_str(), *file_id))
        .collect();
    let wanted_reader = BufReader::new(File::open(wanted_files_path)?);
    let mut wanted = Vec::new();
    for line_result in wanted_reader.lines() {
        let path = line_result?;
        if let Some(&file_id) = by_name.get(path.as_str()) {
            wanted.push((path, file_id));
        } else {
            println!(
                "FILE {} NOT FOUND IN ARCHIVE, CHECK SPELLING OR ACTUAL PRESENCE IN ARCHIVE",
                path
            );
        }
    }
    Ok(wanted)
}

fn decompress_sidecar(
    sidecar_filename: &Path,
    positions_filename: &str,
    tigs_filename: &str,
    size_filename: &str,
    out_dir: &str,
    filenames: &[(String, u32)],
    wanted_filenames: Option<&[(String, u32)]>,
    threads: usize,
    restore_abundance: bool,
    output_compression: OutputCompression,
) -> Result<()> {
    let sidecar = CidDatasetSidecar::open(sidecar_filename)?;
    if restore_abundance && sidecar.abundance_log_base().is_none() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--abundance was requested but this archive has no abundance metadata",
        ));
    }
    if restore_abundance {
        println!("Restoring logarithmically discretized abundance values in FASTA headers");
    }
    let tig_positions = preload_tig_positions(positions_filename)?;
    if tig_positions.len() != sidecar.len().saturating_add(1) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "CID sidecar has {} groups but positions contain {} entries",
                sidecar.len(),
                tig_positions.len()
            ),
        ));
    }

    let wanted_ids = wanted_filenames.map(|wanted| {
        wanted
            .iter()
            .map(|(_, file_id)| *file_id)
            .collect::<HashSet<_>>()
    });
    if let Some(wanted) = wanted_filenames {
        for (path, _) in wanted {
            let output_path = dump_output_path(out_dir, path, output_compression);
            if output_path.exists() {
                fs::remove_file(output_path)?;
            }
        }
    }
    let selected_datasets = wanted_ids
        .as_ref()
        .map_or(filenames.len(), HashSet::len)
        .max(1);
    let thread_budget = threads.max(1);
    let writer_count = selected_datasets.min((thread_budget / 2).max(1)).min(32);
    let decoder_count = thread_budget.saturating_sub(writer_count).max(1).min(16);
    let output_paths = Arc::new(
        filenames
            .iter()
            .map(|(path, _)| dump_output_path(out_dir, path, output_compression))
            .collect::<Vec<_>>(),
    );
    let mut path_shards = HashMap::<PathBuf, usize>::new();
    let writer_routes = Arc::new(
        output_paths
            .iter()
            .map(|path| {
                let next_shard = path_shards.len() % writer_count;
                *path_shards.entry(path.clone()).or_insert(next_shard)
            })
            .collect::<Vec<_>>(),
    );
    println!(
        "Parallel decompression pipeline: {} decoders, {} writer shards, {} MiB batches",
        decoder_count,
        writer_count,
        DECOMPRESS_BATCH_BASES / (1024 * 1024)
    );

    let (job_tx, job_rx) = bounded::<DecodeJob>(decoder_count * 2);
    let (decoded_tx, decoded_rx) = bounded::<Result<DecodedBatch>>(decoder_count * 2);
    let producer_groups = sidecar.len();
    let producer_tigs = tigs_filename.to_owned();
    let producer_sizes = size_filename.to_owned();
    let producer_filenames = filenames.len();
    let producer = thread::spawn(move || {
        produce_decode_jobs(
            sidecar,
            tig_positions,
            &producer_tigs,
            &producer_sizes,
            producer_filenames,
            wanted_ids,
            restore_abundance,
            job_tx,
        )
    });

    let mut decoders = Vec::with_capacity(decoder_count);
    for _ in 0..decoder_count {
        let receiver = job_rx.clone();
        let sender = decoded_tx.clone();
        decoders.push(thread::spawn(move || {
            while let Ok(job) = receiver.recv() {
                let result = decode_fasta_batch(job);
                let failed = result.is_err();
                if sender.send(result).is_err() || failed {
                    while receiver.recv().is_ok() {}
                    break;
                }
            }
        }));
    }
    drop(job_rx);
    drop(decoded_tx);

    let mut writer_senders = Vec::with_capacity(writer_count);
    let mut writers = Vec::with_capacity(writer_count);
    for shard in 0..writer_count {
        let (sender, receiver) = bounded::<WriteBatch>(2);
        writer_senders.push(sender);
        let paths = Arc::clone(&output_paths);
        let routes = Arc::clone(&writer_routes);
        writers.push(thread::spawn(move || {
            write_output_shard(receiver, paths, routes, shard, output_compression)
        }));
    }

    let mut pending = BTreeMap::<u64, DecodedBatch>::new();
    let mut next_ordinal = 0u64;
    let mut pipeline_error = None;
    while let Ok(decoded) = decoded_rx.recv() {
        match decoded {
            Ok(batch) => {
                if pipeline_error.is_some() {
                    continue;
                }
                pending.insert(batch.ordinal, batch);
                while let Some(batch) = pending.remove(&next_ordinal) {
                    if pipeline_error.is_none() {
                        if let Err(err) =
                            dispatch_write_batch(batch, &writer_senders, &writer_routes)
                        {
                            pipeline_error = Some(err);
                        }
                    }
                    next_ordinal += 1;
                }
            }
            Err(err) if pipeline_error.is_none() => pipeline_error = Some(err),
            Err(_) => {}
        }
    }
    drop(writer_senders);

    let producer_stats = producer
        .join()
        .map_err(|_| io::Error::other("decompression producer thread panicked"))?;
    if let Err(err) = producer_stats.as_ref() {
        if pipeline_error.is_none() {
            pipeline_error = Some(io::Error::new(err.kind(), err.to_string()));
        }
    }
    for decoder in decoders {
        if decoder.join().is_err() && pipeline_error.is_none() {
            pipeline_error = Some(io::Error::other("decompression decoder thread panicked"));
        }
    }
    for writer in writers {
        match writer.join() {
            Ok(Ok(())) => {}
            Ok(Err(err)) if pipeline_error.is_none() => pipeline_error = Some(err),
            Err(_) if pipeline_error.is_none() => {
                pipeline_error = Some(io::Error::other("decompression writer thread panicked"));
            }
            _ => {}
        }
    }
    if let Some(err) = pipeline_error {
        return Err(err);
    }
    let stats = producer_stats?;
    if stats.archive_groups != producer_groups {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "decompression producer did not visit every archive group",
        ));
    }
    println!(
        "Decompression complete: {} unitigs written across {} selected CIDs",
        stats.total_unitigs, stats.selected_cids
    );
    Ok(())
}

struct DecodeJob {
    ordinal: u64,
    sizes: Vec<usize>,
    encoded: Vec<u8>,
    dataset_ids: Arc<Vec<u32>>,
    abundance_codes: Option<Arc<Vec<u8>>>,
    abundance_log_base: Option<f64>,
}

struct DecodedBatch {
    ordinal: u64,
    unitigs: usize,
    fasta: Arc<Vec<u8>>,
    abundance_sequences: Option<Arc<DecodedSequences>>,
    dataset_ids: Arc<Vec<u32>>,
    abundance_codes: Option<Arc<Vec<u8>>>,
    abundance_log_base: Option<f64>,
}

struct DecodedSequences {
    bases: Vec<u8>,
    offsets: Vec<u32>,
}

struct WriteBatch {
    dataset_ids: Vec<u32>,
    fasta: Arc<Vec<u8>>,
}

struct DecompressionStats {
    selected_cids: usize,
    total_unitigs: u64,
    archive_groups: usize,
}

#[allow(clippy::too_many_arguments)]
fn produce_decode_jobs(
    sidecar: CidDatasetSidecar,
    tig_positions: Vec<u64>,
    tigs_filename: &str,
    size_filename: &str,
    filename_count: usize,
    wanted_ids: Option<HashSet<u32>>,
    restore_abundance: bool,
    sender: crossbeam::channel::Sender<DecodeJob>,
) -> Result<DecompressionStats> {
    let mut sizes_reader = SequentialCompactSizesReader::open(size_filename)?;
    let mut tigs_file = PackedTigsReader::open(tigs_filename)?;
    let mut selected_ids = Vec::<u32>::new();
    let mut sizes = Vec::<usize>::new();
    let mut selected_cids = 0usize;
    let mut total_unitigs = 0u64;
    let mut ordinal = 0u64;
    let archive_groups = sidecar.len();

    let abundance_log_base = restore_abundance
        .then(|| sidecar.abundance_log_base())
        .flatten();
    sidecar.visit_groups_with_abundance(|cid, dataset_ids, group_abundance_codes| {
        if !sizes_reader.next_group_into(&mut sizes)? {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!("bucket sizes ended before CID {}", cid),
            ));
        }
        selected_ids.clear();
        let mut selected_abundance_codes = Vec::new();
        for (membership_index, &one_based) in dataset_ids.iter().enumerate() {
            let file_id = one_based.checked_sub(1).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "sidecar contains dataset ID zero",
                )
            })?;
            if file_id as usize >= filename_count {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("sidecar dataset ID {} is out of range", one_based),
                ));
            }
            if wanted_ids
                .as_ref()
                .is_none_or(|wanted| wanted.contains(&file_id))
            {
                selected_ids.push(file_id);
                if abundance_log_base.is_some() {
                    let codes = group_abundance_codes.ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            "archive abundance metadata is incomplete",
                        )
                    })?;
                    selected_abundance_codes.push(codes[membership_index]);
                }
            }
        }
        if selected_ids.is_empty() {
            return Ok(());
        }

        selected_cids += 1;
        if selected_cids % 1000 == 1 || cid + 1 == archive_groups {
            println!(
                "Queuing selected CID {} (archive CID {}/{}, {} unitigs queued so far)",
                selected_cids,
                cid + 1,
                archive_groups,
                total_unitigs
            );
        }
        let dataset_ids = Arc::new(selected_ids.clone());
        let abundance_codes =
            abundance_log_base.map(|_| Arc::new(selected_abundance_codes.clone()));
        let mut packed_position = tig_positions[cid];
        let expected_end = tig_positions[cid + 1];
        let mut start = 0usize;
        while start < sizes.len() {
            let mut end = start;
            let mut output_bytes = 0usize;
            let mut encoded_bytes = 0usize;
            while end < sizes.len() {
                let record_bytes = sizes[end].checked_add(3).ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "FASTA record size overflow")
                })?;
                if end > start && output_bytes.saturating_add(record_bytes) > DECOMPRESS_BATCH_BASES
                {
                    break;
                }
                output_bytes = output_bytes.checked_add(record_bytes).ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "FASTA batch size overflow")
                })?;
                encoded_bytes = encoded_bytes
                    .checked_add(sizes[end].div_ceil(4))
                    .ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidData, "packed-tig batch overflow")
                    })?;
                end += 1;
            }
            let mut encoded = vec![0u8; encoded_bytes];
            tigs_file.read_exact_at(packed_position, &mut encoded)?;
            packed_position = packed_position
                .checked_add(encoded_bytes as u64)
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "tig offset overflow"))?;
            sender
                .send(DecodeJob {
                    ordinal,
                    sizes: sizes[start..end].to_vec(),
                    encoded,
                    dataset_ids: Arc::clone(&dataset_ids),
                    abundance_codes: abundance_codes.as_ref().map(Arc::clone),
                    abundance_log_base,
                })
                .map_err(|_| io::Error::new(io::ErrorKind::BrokenPipe, "decoders stopped"))?;
            ordinal = ordinal.checked_add(1).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "batch count overflow")
            })?;
            total_unitigs = total_unitigs
                .checked_add((end - start) as u64)
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "unitig count overflow")
                })?;
            start = end;
        }
        if packed_position != expected_end {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "CID {} sizes describe packed range ending at {}, expected {}",
                    cid, packed_position, expected_end
                ),
            ));
        }
        Ok(())
    })?;
    if sizes_reader.next_group_into(&mut sizes)? {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "bucket sizes contain more groups than the CID sidecar",
        ));
    }
    Ok(DecompressionStats {
        selected_cids,
        total_unitigs,
        archive_groups,
    })
}

fn decode_fasta_batch(job: DecodeJob) -> Result<DecodedBatch> {
    let output_bytes = job.sizes.iter().try_fold(0usize, |total, &size| {
        size.checked_add(3)
            .and_then(|record| total.checked_add(record))
    });
    let output_bytes = output_bytes.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "decoded FASTA batch size overflow",
        )
    })?;
    let build_fasta = |abundance: Option<u64>| -> Result<Vec<u8>> {
        let header_extra = abundance.map_or(0, |value| value.to_string().len() + 8);
        let mut fasta = Vec::with_capacity(
            output_bytes.saturating_add(header_extra.saturating_mul(job.sizes.len())),
        );
        let mut cursor = 0usize;
        for &size in &job.sizes {
            let encoded_len = size.div_ceil(4);
            let end = cursor.checked_add(encoded_len).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "packed-tig cursor overflow")
            })?;
            let encoded = job.encoded.get(cursor..end).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "truncated packed-tig decode batch",
                )
            })?;
            if let Some(value) = abundance {
                write!(fasta, "> ka:f:{value}\n")?;
            } else {
                fasta.extend_from_slice(b">\n");
            }
            let mut remaining = size;
            for &packed in encoded {
                let count = remaining.min(4);
                fasta.extend_from_slice(&PACKED_BASES[packed as usize][..count]);
                remaining -= count;
            }
            fasta.push(b'\n');
            cursor = end;
        }
        if cursor != job.encoded.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "packed-tig decode batch contains trailing bytes",
            ));
        }
        Ok(fasta)
    };

    let (fasta, abundance_sequences) = match (job.abundance_codes.as_ref(), job.abundance_log_base)
    {
        (None, None) => (Arc::new(build_fasta(None)?), None),
        (Some(codes), Some(_)) => {
            if codes.len() != job.dataset_ids.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "decoded abundance count differs from dataset membership count",
                ));
            }
            let total_bases = job
                .sizes
                .iter()
                .try_fold(0usize, |total, size| total.checked_add(*size))
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "decoded base count overflow")
                })?;
            let mut bases = Vec::with_capacity(total_bases);
            let mut offsets = Vec::with_capacity(job.sizes.len() + 1);
            offsets.push(0);
            let mut cursor = 0usize;
            for &size in &job.sizes {
                let encoded_len = size.div_ceil(4);
                let end = cursor.checked_add(encoded_len).ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "packed-tig cursor overflow")
                })?;
                let encoded = job.encoded.get(cursor..end).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::UnexpectedEof,
                        "truncated packed-tig decode batch",
                    )
                })?;
                let mut remaining = size;
                for &packed in encoded {
                    let count = remaining.min(4);
                    bases.extend_from_slice(&PACKED_BASES[packed as usize][..count]);
                    remaining -= count;
                }
                offsets.push(u32::try_from(bases.len()).map_err(|_| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "decoded abundance batch exceeds 4 GiB",
                    )
                })?);
                cursor = end;
            }
            if cursor != job.encoded.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "packed-tig decode batch contains trailing bytes",
                ));
            }
            (
                Arc::new(Vec::new()),
                Some(Arc::new(DecodedSequences { bases, offsets })),
            )
        }
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "archive abundance metadata is incomplete",
            ))
        }
    };
    Ok(DecodedBatch {
        ordinal: job.ordinal,
        unitigs: job.sizes.len(),
        fasta,
        abundance_sequences,
        dataset_ids: job.dataset_ids,
        abundance_codes: job.abundance_codes,
        abundance_log_base: job.abundance_log_base,
    })
}

fn build_abundance_fasta(sequences: &DecodedSequences, abundance: u64) -> Result<Vec<u8>> {
    let header = format!("> ka:f:{abundance}\n");
    let mut fasta = Vec::with_capacity(
        sequences
            .bases
            .len()
            .saturating_add(header.len().saturating_mul(sequences.offsets.len())),
    );
    for range in sequences.offsets.windows(2) {
        let start = range[0] as usize;
        let end = range[1] as usize;
        let sequence = sequences.bases.get(start..end).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "decoded abundance sequence offsets are invalid",
            )
        })?;
        fasta.extend_from_slice(header.as_bytes());
        fasta.extend_from_slice(sequence);
        fasta.push(b'\n');
    }
    Ok(fasta)
}

fn dispatch_write_batch(
    batch: DecodedBatch,
    writer_senders: &[crossbeam::channel::Sender<WriteBatch>],
    writer_routes: &[usize],
) -> Result<()> {
    if batch.unitigs == 0 {
        return Ok(());
    }
    let mut shard_ids = (0..writer_senders.len())
        .map(|_| Vec::<(u32, Option<u8>)>::new())
        .collect::<Vec<_>>();
    for (index, &file_id) in batch.dataset_ids.iter().enumerate() {
        let shard = *writer_routes.get(file_id as usize).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "output dataset ID is out of range",
            )
        })?;
        let abundance = batch.abundance_codes.as_ref().map(|codes| codes[index]);
        shard_ids[shard].push((file_id, abundance));
    }
    if let Some(sequences) = batch.abundance_sequences.as_ref() {
        let log_base = batch.abundance_log_base.ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "archive abundance log base is missing",
            )
        })?;
        let mut codes = batch
            .abundance_codes
            .as_ref()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "archive abundance codes are missing",
                )
            })?
            .as_ref()
            .clone();
        codes.sort_unstable();
        codes.dedup();
        for code in codes {
            let abundance = crate::compress::decode_abundance(code, log_base);
            let fasta = Arc::new(build_abundance_fasta(sequences, abundance)?);
            for (shard, memberships) in shard_ids.iter().enumerate() {
                let dataset_ids = memberships
                    .iter()
                    .filter_map(|&(file_id, membership_code)| {
                        (membership_code == Some(code)).then_some(file_id)
                    })
                    .collect::<Vec<_>>();
                if !dataset_ids.is_empty() {
                    writer_senders[shard]
                        .send(WriteBatch {
                            dataset_ids,
                            fasta: Arc::clone(&fasta),
                        })
                        .map_err(|_| {
                            io::Error::new(io::ErrorKind::BrokenPipe, "writer shard stopped")
                        })?;
                }
            }
        }
    } else {
        for (shard, memberships) in shard_ids.into_iter().enumerate() {
            if memberships.is_empty() {
                continue;
            }
            writer_senders[shard]
                .send(WriteBatch {
                    dataset_ids: memberships
                        .into_iter()
                        .map(|(file_id, _)| file_id)
                        .collect(),
                    fasta: Arc::clone(&batch.fasta),
                })
                .map_err(|_| io::Error::new(io::ErrorKind::BrokenPipe, "writer shard stopped"))?;
        }
    }
    Ok(())
}

fn write_output_shard(
    receiver: Receiver<WriteBatch>,
    output_paths: Arc<Vec<PathBuf>>,
    writer_routes: Arc<Vec<usize>>,
    shard: usize,
    output_compression: OutputCompression,
) -> Result<()> {
    let mut writers = HashMap::<u32, FastaOutputWriter>::new();
    while let Ok(batch) = receiver.recv() {
        for file_id in batch.dataset_ids {
            if writer_routes.get(file_id as usize).copied() != Some(shard) {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "dataset was routed to the wrong writer shard",
                ));
            }
            let path = output_paths.get(file_id as usize).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "output dataset ID is out of range",
                )
            })?;
            let writer = match writers.entry(file_id) {
                std::collections::hash_map::Entry::Occupied(entry) => entry.into_mut(),
                std::collections::hash_map::Entry::Vacant(entry) => {
                    entry.insert(FastaOutputWriter::create(path, output_compression)?)
                }
            };
            writer.write_all(batch.fasta.as_slice())?;
        }
    }
    for (_, writer) in writers {
        writer.finish()?;
    }
    Ok(())
}

#[cfg(test)]
mod parallel_decompression_tests {
    use super::*;

    fn decode_test_output(path: &Path, compression: OutputCompression) -> Vec<u8> {
        let file = File::open(path).unwrap();
        let mut decoded = Vec::new();
        match compression {
            OutputCompression::Fasta => BufReader::new(file).read_to_end(&mut decoded).unwrap(),
            OutputCompression::Gz => flate2::read::GzDecoder::new(file)
                .read_to_end(&mut decoded)
                .unwrap(),
            OutputCompression::Zstd => zstd::Decoder::new(file)
                .unwrap()
                .read_to_end(&mut decoded)
                .unwrap(),
            OutputCompression::Xz => xz2::read::XzDecoder::new(file)
                .read_to_end(&mut decoded)
                .unwrap(),
        };
        decoded
    }

    #[test]
    fn output_writers_roundtrip_all_formats() {
        let temp = tempfile::tempdir().unwrap();
        let fasta = b">record\nACGTACGT\n";
        for compression in [
            OutputCompression::Fasta,
            OutputCompression::Gz,
            OutputCompression::Zstd,
            OutputCompression::Xz,
        ] {
            let path = temp
                .path()
                .join(format!("output{}", compression.fasta_suffix()));
            let mut writer = FastaOutputWriter::create(&path, compression).unwrap();
            writer.write_all(fasta).unwrap();
            writer.finish().unwrap();
            assert_eq!(decode_test_output(&path, compression), fasta);
        }
    }

    #[test]
    fn default_output_is_zstd() {
        assert_eq!(
            GgcatRebuildConfig::default().output_compression,
            OutputCompression::Zstd
        );
        assert_eq!(OutputCompression::Zstd.fasta_suffix(), ".fa.zst");
    }

    #[test]
    fn packed_batch_decodes_directly_to_ordered_fasta() {
        let decoded = decode_fasta_batch(DecodeJob {
            ordinal: 7,
            sizes: vec![5, 3],
            encoded: vec![0xe4, 0x00, 0x1b],
            dataset_ids: Arc::new(vec![0, 2]),
            abundance_codes: None,
            abundance_log_base: None,
        })
        .unwrap();

        assert_eq!(decoded.ordinal, 7);
        assert_eq!(decoded.unitigs, 2);
        assert_eq!(decoded.fasta.as_slice(), b">\nACGTA\n>\nTGC\n");
        assert_eq!(decoded.dataset_ids.as_slice(), &[0, 2]);
    }

    #[test]
    fn abundance_batch_decodes_bases_once_and_restores_header() {
        let code = crate::compress::encode_abundance(17, 1.05);
        let decoded = decode_fasta_batch(DecodeJob {
            ordinal: 9,
            sizes: vec![5, 3],
            encoded: vec![0xe4, 0x00, 0x1b],
            dataset_ids: Arc::new(vec![0]),
            abundance_codes: Some(Arc::new(vec![code])),
            abundance_log_base: Some(1.05),
        })
        .unwrap();
        assert!(decoded.fasta.is_empty());
        let sequences = decoded.abundance_sequences.unwrap();
        let abundance = crate::compress::decode_abundance(code, 1.05);
        let fasta = build_abundance_fasta(&sequences, abundance).unwrap();
        assert_eq!(fasta, b"> ka:f:17\nACGTA\n> ka:f:17\nTGC\n");
    }
}

struct SequentialCompactSizesReader {
    file: BufReader<File>,
    length_prefixed: bool,
    ranges: Vec<(usize, usize)>,
    data: Vec<u8>,
    next_range: usize,
    finished: bool,
}

impl SequentialCompactSizesReader {
    fn open(path: &str) -> Result<Self> {
        let mut file = BufReader::new(File::open(path)?);
        let mut magic = [0u8; 4];
        file.read_exact(&mut magic)?;
        let length_prefixed = if &magic == BUCKET_SIZES_MAGIC {
            false
        } else if &magic == BUCKET_SIZES_COMPACT_MAGIC {
            true
        } else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("unsupported compact bucket-size encoding in '{}'", path),
            ));
        };
        Ok(Self {
            file,
            length_prefixed,
            ranges: Vec::new(),
            data: Vec::new(),
            next_range: 0,
            finished: false,
        })
    }

    fn next_group_into(&mut self, sizes: &mut Vec<usize>) -> Result<bool> {
        if self.next_range == self.ranges.len() && !self.load_next_block()? {
            return Ok(false);
        }
        let (start, end) = self.ranges[self.next_range];
        self.next_range += 1;
        decode_varint_deltas_to_sizes_into(&self.data[start..end], "compact bucket group", sizes)?;
        Ok(true)
    }

    fn load_next_block(&mut self) -> Result<bool> {
        if self.finished {
            return Ok(false);
        }
        let mut groups_buf = [0u8; 4];
        match self.file.read_exact(&mut groups_buf) {
            Ok(()) => {}
            Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => {
                self.finished = true;
                return Ok(false);
            }
            Err(err) => return Err(err),
        }
        let groups = u32::from_le_bytes(groups_buf) as usize;
        if groups == 0 {
            self.finished = true;
            return Ok(false);
        }
        let mut compressed_len_buf = [0u8; 8];
        self.file.read_exact(&mut compressed_len_buf)?;
        let compressed_len =
            usize::try_from(u64::from_le_bytes(compressed_len_buf)).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "bucket-size block is too large")
            })?;

        let mut offsets = Vec::new();
        if !self.length_prefixed {
            offsets.resize(groups + 1, 0u32);
            let mut bytes = [0u8; 4];
            for offset in &mut offsets {
                self.file.read_exact(&mut bytes)?;
                *offset = u32::from_le_bytes(bytes);
            }
        }
        let mut compressed = vec![0u8; compressed_len];
        self.file.read_exact(&mut compressed)?;
        self.data.clear();
        Decoder::new(compressed.as_slice())?.read_to_end(&mut self.data)?;
        self.ranges.clear();
        self.ranges.reserve(groups);
        self.next_range = 0;

        if self.length_prefixed {
            let mut cursor = std::io::Cursor::new(self.data.as_slice());
            for _ in 0..groups {
                let group_len = read_varint_u64_from_reader(&mut cursor)?.ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::UnexpectedEof,
                        "missing compact bucket group length",
                    )
                })?;
                let start = cursor.position() as usize;
                let end = start
                    .checked_add(usize::try_from(group_len).map_err(|_| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            "compact bucket group is too large",
                        )
                    })?)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            "compact bucket group offset overflow",
                        )
                    })?;
                if end > self.data.len() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "compact bucket group extends past its block",
                    ));
                }
                self.ranges.push((start, end));
                cursor.set_position(end as u64);
            }
            if cursor.position() as usize != self.data.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "compact bucket block has trailing data",
                ));
            }
        } else {
            if offsets.windows(2).any(|range| range[0] > range[1])
                || offsets.last().copied().map(|value| value as usize) != Some(self.data.len())
            {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "compact bucket block offsets are invalid",
                ));
            }
            self.ranges.extend(
                offsets
                    .windows(2)
                    .map(|range| (range[0] as usize, range[1] as usize)),
            );
        }
        Ok(true)
    }
}

fn preload_tig_positions(positions_filename: &str) -> Result<Vec<u64>> {
    let mut file = BufReader::new(File::open(positions_filename)?);
    let mut magic = [0u8; 4];
    file.read_exact(&mut magic)?;
    if &magic == POSITIONS_MAGIC || &magic == POSITIONS_COMPACT_MAGIC {
        let implicit_sizes = &magic == POSITIONS_COMPACT_MAGIC;
        let entries = read_varint_u64_from_reader(&mut file)?.ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "missing positions entry count",
            )
        })?;
        let entries = usize::try_from(entries).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "positions entry count is too large",
            )
        })?;
        let mut positions = Vec::with_capacity(entries);
        let mut position = 0u64;
        for _ in 0..entries {
            let delta = read_varint_u64_from_reader(&mut file)?.ok_or_else(|| {
                io::Error::new(io::ErrorKind::UnexpectedEof, "truncated tig position")
            })?;
            position = position.checked_add(delta).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "tig position overflow")
            })?;
            positions.push(position);
            if !implicit_sizes {
                read_varint_u64_from_reader(&mut file)?.ok_or_else(|| {
                    io::Error::new(io::ErrorKind::UnexpectedEof, "truncated size position")
                })?;
            }
        }
        return Ok(positions);
    }

    file.seek(std::io::SeekFrom::Start(0))?;
    let file_size = file.get_ref().metadata()?.len() as usize;
    if file_size % 16 != 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "legacy positions size is not divisible by 16",
        ));
    }
    let mut positions = Vec::with_capacity(file_size / 16);
    let mut entry = [0u8; 16];
    while file.read_exact(&mut entry).is_ok() {
        positions.push(u64::from_le_bytes(entry[..8].try_into().unwrap()));
    }
    Ok(positions)
}

/// Preload all positions from the raw binary positions file.
///
/// The positions file contains N entries of exactly 16 bytes each:
/// [u64 tigs_pos LE][u64 sizes_pos LE].
/// Returns a Vec where index i corresponds to CID i.
fn read_varint_u64_from_reader(reader: &mut impl Read) -> io::Result<Option<u64>> {
    let mut shift = 0u32;
    let mut value = 0u64;
    let mut saw_byte = false;
    loop {
        let mut b = [0u8; 1];
        match reader.read_exact(&mut b) {
            Ok(()) => {}
            Err(err) if err.kind() == io::ErrorKind::UnexpectedEof && !saw_byte => return Ok(None),
            Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "truncated varint value",
                ))
            }
            Err(err) => return Err(err),
        }
        saw_byte = true;
        let byte = b[0];
        value |= ((byte & 0x7f) as u64) << shift;
        if byte & 0x80 == 0 {
            return Ok(Some(value));
        }
        shift += 7;
        if shift >= 64 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "varint overflow while decoding u64",
            ));
        }
    }
}

fn decode_varint_delta_cids(payload: &[u8], context: &str) -> Result<Vec<usize>> {
    let mut cursor = std::io::Cursor::new(payload);
    let mut cids = Vec::new();
    let mut current = 0usize;
    loop {
        let Some(delta_u64) = read_varint_u64_from_reader(&mut cursor)? else {
            break;
        };
        let delta = usize::try_from(delta_u64).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "cid delta {} cannot fit usize while decoding '{}'",
                    delta_u64, context
                ),
            )
        })?;
        current = current.checked_add(delta).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("cid delta overflow while decoding '{}'", context),
            )
        })?;
        cids.push(current);
    }
    Ok(cids)
}

fn decode_varint_deltas_to_sizes(payload: &[u8], context: &str) -> Result<Vec<usize>> {
    let mut sizes = Vec::new();
    decode_varint_deltas_to_sizes_into(payload, context, &mut sizes)?;
    Ok(sizes)
}

fn decode_varint_deltas_to_sizes_into(
    payload: &[u8],
    context: &str,
    sizes: &mut Vec<usize>,
) -> Result<()> {
    let mut cursor = std::io::Cursor::new(payload);
    sizes.clear();
    let mut prev = 0usize;
    loop {
        let Some(delta_u64) = read_varint_u64_from_reader(&mut cursor)? else {
            break;
        };
        let delta = usize::try_from(delta_u64).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "size delta {} cannot fit usize while decoding '{}'",
                    delta_u64, context
                ),
            )
        })?;
        let size = prev.checked_add(delta).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("size delta overflow while decoding '{}'", context),
            )
        })?;
        sizes.push(size);
        prev = size;
    }
    Ok(())
}

fn preload_positions(positions_filename: &str) -> Result<Vec<(u64, u64)>> {
    let mut file = BufReader::new(File::open(positions_filename)?);
    let mut magic = [0u8; 4];
    match file.read_exact(&mut magic) {
        Ok(()) => {}
        Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => return Ok(Vec::new()),
        Err(err) => return Err(err),
    }

    if &magic == POSITIONS_MAGIC || &magic == POSITIONS_COMPACT_MAGIC {
        let implicit_sizes = &magic == POSITIONS_COMPACT_MAGIC;
        let entries = read_varint_u64_from_reader(&mut file)?.ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "missing positions entry count in compact positions file",
            )
        })?;
        let entries = usize::try_from(entries).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "positions entry count cannot fit usize",
            )
        })?;
        let mut positions = Vec::with_capacity(entries);
        let mut tigs_pos = 0u64;
        let mut sizes_pos = 0u64;
        for entry in 0..entries {
            let dt = read_varint_u64_from_reader(&mut file)?.ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "truncated tigs delta in positions file",
                )
            })?;
            let ds = if implicit_sizes {
                u64::from(entry > 0)
            } else {
                read_varint_u64_from_reader(&mut file)?.ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::UnexpectedEof,
                        "truncated sizes delta in positions file",
                    )
                })?
            };
            tigs_pos = tigs_pos.checked_add(dt).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "tigs position delta overflow in positions file",
                )
            })?;
            sizes_pos = sizes_pos.checked_add(ds).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "sizes position delta overflow in positions file",
                )
            })?;
            positions.push((tigs_pos, sizes_pos));
        }
        println!(
            "Preloaded {} compact position entries from {}",
            positions.len(),
            positions_filename
        );
        return Ok(positions);
    }

    file.seek(std::io::SeekFrom::Start(0))?;
    let file_size = file.get_ref().metadata()?.len() as usize;
    let num_entries = file_size / 16;
    let mut positions = Vec::with_capacity(num_entries);
    let mut buf = [0u8; 16];
    for _ in 0..num_entries {
        file.read_exact(&mut buf)?;
        let tigs_pos = u64::from_le_bytes(buf[..8].try_into().unwrap());
        let sizes_pos = u64::from_le_bytes(buf[8..16].try_into().unwrap());
        positions.push((tigs_pos, sizes_pos));
    }
    println!(
        "Preloaded {} legacy position entries from {}",
        num_entries, positions_filename
    );
    Ok(positions)
}

/// Decode comma-separated delta-encoded CIDs into absolute CID indexes.
///
/// Example: "3,2,0,5" => [3,5,5,10]
fn decode_delta_cids(text: &str, context: &str) -> Result<Vec<usize>> {
    let mut cids = Vec::new();
    let mut current = 0usize;
    for token in text.split(',') {
        let token = token.trim();
        if token.is_empty() {
            continue;
        }
        let delta = token.parse::<usize>().map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid cid delta '{}' in '{}': {}", token, context, err),
            )
        })?;
        current = current.checked_add(delta).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("cid delta overflow while decoding '{}'", context),
            )
        })?;
        cids.push(current);
    }
    Ok(cids)
}

fn decode_cid_entry_payload(
    payload: &[u8],
    context: &str,
    binary_varints: bool,
) -> Result<Vec<usize>> {
    if binary_varints {
        decode_varint_delta_cids(payload, context)
    } else {
        let text = String::from_utf8(payload.to_vec()).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid UTF-8 while decoding '{}': {}", context, err),
            )
        })?;
        decode_delta_cids(&text, context)
    }
}

/// Preload all bucket sizes from the sizes file into a Vec indexed by byte offset.
///
/// The sizes file format: repeated [u64 compressed_len][compressed_data...].
/// Returns a HashMap from byte_offset -> Vec<usize> of actual sizes.
fn preload_sizes(size_filename: &str) -> Result<HashMap<u64, Vec<usize>>> {
    let mut file = BufReader::new(File::open(size_filename)?);
    let mut sizes_map = HashMap::new();
    let mut magic = [0u8; 4];
    match file.read_exact(&mut magic) {
        Ok(()) => {}
        Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => return Ok(HashMap::new()),
        Err(err) => return Err(err),
    }

    if &magic == BUCKET_SIZES_MAGIC || &magic == BUCKET_SIZES_COMPACT_MAGIC {
        let length_prefixed = &magic == BUCKET_SIZES_COMPACT_MAGIC;
        let mut global_group_index = 0u64;
        loop {
            let mut groups_buf = [0u8; 4];
            match file.read_exact(&mut groups_buf) {
                Ok(()) => {}
                Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => break,
                Err(err) => return Err(err),
            }
            let groups = u32::from_le_bytes(groups_buf) as usize;
            if groups == 0 {
                break;
            }
            let mut compressed_len_buf = [0u8; 8];
            file.read_exact(&mut compressed_len_buf)?;
            let compressed_len = u64::from_le_bytes(compressed_len_buf) as usize;

            let mut offsets = Vec::new();
            if !length_prefixed {
                offsets.resize(groups + 1, 0u32);
                for off in &mut offsets {
                    let mut buf = [0u8; 4];
                    file.read_exact(&mut buf)?;
                    *off = u32::from_le_bytes(buf);
                }
            }

            let mut compressed = vec![0u8; compressed_len];
            file.read_exact(&mut compressed)?;

            let mut decompressed = Vec::new();
            Decoder::new(&compressed[..])?.read_to_end(&mut decompressed)?;

            let mut compact_cursor = std::io::Cursor::new(decompressed.as_slice());
            for i in 0..groups {
                let (start, end) = if length_prefixed {
                    let group_len =
                        read_varint_u64_from_reader(&mut compact_cursor)?.ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::UnexpectedEof,
                                "missing compact bucket group length",
                            )
                        })?;
                    let start = compact_cursor.position() as usize;
                    let end = start
                        .checked_add(usize::try_from(group_len).map_err(|_| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                "compact bucket group is too large",
                            )
                        })?)
                        .ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                "compact bucket group offset overflow",
                            )
                        })?;
                    compact_cursor.set_position(end as u64);
                    (start, end)
                } else {
                    (offsets[i] as usize, offsets[i + 1] as usize)
                };
                if end < start || end > decompressed.len() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "invalid bucket block offsets while preloading sizes",
                    ));
                }
                let key = global_group_index + i as u64;
                let context = format!("{}#group{}", size_filename, key);
                let sizes = decode_varint_deltas_to_sizes(&decompressed[start..end], &context)?;
                sizes_map.insert(key, sizes);
            }
            if length_prefixed && compact_cursor.position() as usize != decompressed.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "compact bucket block has trailing group data",
                ));
            }
            global_group_index += groups as u64;
        }

        println!(
            "Preloaded {} compact size buckets from {}",
            sizes_map.len(),
            size_filename
        );
        return Ok(sizes_map);
    }

    file.seek(std::io::SeekFrom::Start(0))?;
    let mut offset: u64 = 0;
    loop {
        let mut len_buf = [0u8; 8];
        match file.read_exact(&mut len_buf) {
            Ok(()) => {}
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => break,
            Err(e) => return Err(e),
        }
        let compressed_size = u64::from_le_bytes(len_buf) as usize;
        if compressed_size == 0 {
            break;
        }
        let mut compressed_buffer = vec![0; compressed_size];
        file.read_exact(&mut compressed_buffer)?;

        let mut decoder = Decoder::new(&compressed_buffer[..])?;
        let mut decompressed = Vec::new();
        decoder.read_to_end(&mut decompressed)?;

        let mut sizes = Vec::new();
        let mut prev = 0usize;
        for chunk in decompressed.chunks_exact(8) {
            let delta = usize::from_le_bytes(chunk.try_into().unwrap());
            let actual_size = prev
                .checked_add(delta)
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "size overflow"))?;
            sizes.push(actual_size);
            prev = actual_size;
        }

        sizes_map.insert(offset, sizes);
        offset += 8 + compressed_size as u64;
    }
    println!(
        "Preloaded {} legacy size buckets from {}",
        sizes_map.len(),
        size_filename
    );
    Ok(sizes_map)
}

/// Decompress the entire archive to individual FASTA shards.
///
/// Uses preloaded positions and sizes to avoid per-CID file reopens.
/// Keeps output file handles open in a HashMap to avoid per-unitig open/close.
fn decompress_all(
    size_filename: &String,
    positions_filename: &String,
    tigs_filename: &String,
    out_dir: &String,
    filenames: Vec<(String, u32)>,
    cid_to_id_map: HashMap<usize, Vec<u32>>,
    output_compression: OutputCompression,
) -> Result<()> {
    let mut tigs_file =
        BufReader::new(File::open(&tigs_filename).expect("Error opening tigs file"));

    let total_cids = cid_to_id_map.len();
    println!("nb cid: {}", total_cids);

    // Preload all positions and sizes into memory
    let all_positions = preload_positions(positions_filename).expect("Failed to preload positions");
    let all_sizes = preload_sizes(size_filename).expect("Failed to preload sizes");

    // Keep output file handles open
    let mut writers: HashMap<u32, FastaOutputWriter> = HashMap::new();

    let mut sorted_cids: Vec<_> = cid_to_id_map.keys().cloned().collect();
    sorted_cids.sort();

    let mut total_unitigs = 0_u64;
    for (idx, cid) in sorted_cids.iter().enumerate() {
        let file_ids = cid_to_id_map.get(cid).unwrap();

        // Progress logging
        if idx % 1000 == 0 || idx == total_cids - 1 {
            println!(
                "Processing CID {}/{} ({} unitigs written so far)",
                idx + 1,
                total_cids,
                total_unitigs
            );
        }

        // Look up position from preloaded data
        let pos_index = *cid;
        if pos_index >= all_positions.len() {
            eprintln!(
                "Error: position index {} out of range for CID {}",
                pos_index, cid
            );
            continue;
        }
        let (tigs_pos, sizes_pos) = all_positions[pos_index];

        // Look up sizes from preloaded data
        let sizes = match all_sizes.get(&sizes_pos) {
            Some(s) => s,
            None => {
                eprintln!(
                    "Error: no sizes found at offset {} for CID {}",
                    sizes_pos, cid
                );
                continue;
            }
        };

        if sizes.is_empty() {
            continue;
        }

        tigs_file
            .seek(std::io::SeekFrom::Start(tigs_pos))
            .expect("Failed to seek in tigs file");

        for size in sizes {
            if *size < 31 {
                eprintln!("Warning: unitig size {} is less than k-mer size", size);
                continue;
            }

            let read_size = size.div_ceil(4);
            let mut tig_buffer = vec![0; read_size];
            tigs_file
                .read_exact(&mut tig_buffer)
                .expect("Failed to read tig");

            let tig = vec2str(&tig_buffer, size);
            for file_id in file_ids {
                if let std::collections::hash_map::Entry::Vacant(entry) = writers.entry(*file_id) {
                    let curr_filename = &filenames[*file_id as usize];
                    let output_path =
                        dump_output_path(out_dir, &curr_filename.0, output_compression);
                    entry.insert(FastaOutputWriter::create(&output_path, output_compression)?);
                }
                let writer = writers.get_mut(file_id).expect("writer inserted above");
                writeln!(writer, ">")?;
                writeln!(writer, "{}", tig)?;
            }
            total_unitigs += 1;
        }
    }

    // Flush all writers at the end
    for (_, writer) in writers {
        writer.finish()?;
    }
    println!(
        "Decompression complete: {} unitigs written across {} CIDs",
        total_unitigs, total_cids
    );
    Ok(())
}

/// Read full id->color_id file and build color id -> list of file ids.
fn get_cid_to_id(color_id_filename: &String) -> Result<HashMap<usize, Vec<u32>>> {
    let mut color_id_file = BufReader::new(
        File::open(color_id_filename)
            .expect("Error opening color id file, are you sure you gave the right path?"),
    );
    let mut cid_ids_map = HashMap::new();
    let mut counter: u32 = 0;
    println!("Reading CID to ID file: {}", color_id_filename);

    let mut magic = [0u8; 4];
    let binary_varints = match color_id_file.read_exact(&mut magic) {
        Ok(()) if &magic == ID_TO_CID_MAGIC => true,
        Ok(()) => {
            color_id_file.seek(std::io::SeekFrom::Start(0))?;
            false
        }
        Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => return Ok(HashMap::new()),
        Err(err) => return Err(err),
    };

    let mut buffer_size = [0; 8];
    color_id_file.read_exact(&mut buffer_size)?;
    let mut size_read = u64::from_le_bytes(buffer_size) as usize;
    while size_read != 0 {
        let mut buffer = vec![0; size_read];
        color_id_file.read_exact(&mut buffer)?;
        let mut decompressed_data = Vec::new();
        {
            let mut decoder = Decoder::new(&buffer[..])?;
            decoder.read_to_end(&mut decompressed_data)?;
        }
        let context = format!("{}#{}", color_id_filename, counter);
        let decoded_cids = decode_cid_entry_payload(&decompressed_data, &context, binary_varints)?;
        for cid in decoded_cids {
            cid_ids_map
                .entry(cid)
                .and_modify(|list: &mut Vec<_>| list.push(counter))
                .or_insert(Vec::from([counter]));
        }
        color_id_file.read_exact(&mut buffer_size)?;
        size_read = u64::from_le_bytes(buffer_size) as usize;
        counter += 1;
    }

    Ok(cid_ids_map)
}

/// Build cid -> file id mapping for a targeted subset of files.
fn get_cid_to_id_targeted(
    color_id_filename: &String,
    filenames_id_map: &HashMap<String, (u32, u64)>,
    wanted_files_path: &String,
) -> std::io::Result<(HashMap<usize, Vec<u32>>, Vec<(String, u32)>)> {
    let mut color_id_file = BufReader::new(
        File::open(color_id_filename)
            .expect("Error opening color id file, are you sure you gave the right path?"),
    );
    let mut cid_ids_map = HashMap::new();
    let mut wanted_filenames = Vec::new();

    let mut magic = [0u8; 4];
    let binary_varints = match color_id_file.read_exact(&mut magic) {
        Ok(()) if &magic == ID_TO_CID_MAGIC => true,
        Ok(()) => {
            color_id_file.seek(std::io::SeekFrom::Start(0))?;
            false
        }
        Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => {
            return Ok((HashMap::new(), Vec::new()))
        }
        Err(err) => return Err(err),
    };

    let wanted_file = File::open(wanted_files_path)?;
    let wanted_reader = BufReader::new(wanted_file);
    for line_result in wanted_reader.lines() {
        let line = line_result?;
        if filenames_id_map.contains_key(&line) {
            let entry = filenames_id_map.get(&line).unwrap();
            color_id_file.seek(std::io::SeekFrom::Start(entry.1))?;
            let mut buffer_size = [0; 8];
            color_id_file.read_exact(&mut buffer_size)?;
            let size_read = u64::from_le_bytes(buffer_size) as usize;
            let mut buffer = vec![0; size_read];
            color_id_file.read_exact(&mut buffer)?;
            let mut decompressed_data = Vec::new();
            {
                let mut decoder = Decoder::new(&buffer[..])?;
                decoder.read_to_end(&mut decompressed_data)?;
            }
            let context = format!("{}@{}", color_id_filename, entry.1);
            let decoded_cids =
                decode_cid_entry_payload(&decompressed_data, &context, binary_varints)?;
            for cid in decoded_cids {
                cid_ids_map
                    .entry(cid)
                    .and_modify(|list: &mut Vec<_>| list.push(entry.0))
                    .or_insert(Vec::from([entry.0]));
            }
            wanted_filenames.push((line.clone(), entry.0));
        } else {
            println!(
                "FILE {} NOT FOUND IN ARCHIVE, CHECK SPELLING OR ACTUAL PRESENCE IN ARCHIVE",
                line
            );
        }
    }
    Ok((cid_ids_map, wanted_filenames))
}

/// Decompress only the wanted files from the archive.
///
/// Uses preloaded positions and sizes. Keeps output file handles open.
fn decompress_wanted(
    wanted_files: &Vec<(String, u32)>,
    positions_filename: &String,
    cid_to_id_map: HashMap<usize, Vec<u32>>,
    tigs_filename: &String,
    size_filename: &String,
    out_dir: &String,
    output_compression: OutputCompression,
) -> Result<()> {
    let total_cids = cid_to_id_map.len();
    println!("NB COLOR TO DECOMPRESS: {}", total_cids);
    println!("Wanted files:");
    for elem in wanted_files {
        println!("{} : {}", elem.0, elem.1);
    }

    let mut tigs_file =
        BufReader::new(File::open(&tigs_filename).expect("Error opening tigs file"));

    // Preload all positions and sizes into memory
    let all_positions = preload_positions(positions_filename).expect("Failed to preload positions");
    let all_sizes = preload_sizes(size_filename).expect("Failed to preload sizes");

    // Build a set of wanted file_ids for quick lookup
    let wanted_ids: std::collections::HashSet<u32> = wanted_files.iter().map(|w| w.1).collect();

    // Keep output file handles open
    let mut writers: HashMap<u32, FastaOutputWriter> = HashMap::new();
    // Pre-open all wanted output files
    for wanted_file in wanted_files {
        let output_path = dump_output_path(out_dir, &wanted_file.0, output_compression);
        let writer = FastaOutputWriter::create(&output_path, output_compression)?;
        writers.insert(wanted_file.1, writer);
    }

    let mut sorted_cids: Vec<_> = cid_to_id_map.keys().cloned().collect();
    sorted_cids.sort();

    let mut total_unitigs = 0_u64;
    for (idx, cid) in sorted_cids.iter().enumerate() {
        let file_ids = cid_to_id_map.get(cid).unwrap();

        // Progress logging
        if idx % 1000 == 0 || idx == total_cids - 1 {
            println!(
                "Processing CID {}/{} ({} unitigs written so far)",
                idx + 1,
                total_cids,
                total_unitigs
            );
        }

        // Look up position from preloaded data
        let pos_index = *cid;
        if pos_index >= all_positions.len() {
            eprintln!(
                "Error: position index {} out of range for CID {}",
                pos_index, cid
            );
            continue;
        }
        let (tigs_pos, sizes_pos) = all_positions[pos_index];

        // Look up sizes from preloaded data
        let sizes = match all_sizes.get(&sizes_pos) {
            Some(s) => s,
            None => {
                eprintln!(
                    "Error: no sizes found at offset {} for CID {}",
                    sizes_pos, cid
                );
                continue;
            }
        };

        if sizes.is_empty() {
            continue;
        }
        tigs_file
            .seek(std::io::SeekFrom::Start(tigs_pos))
            .expect("Failed to seek in tigs file");

        for size in sizes {
            if *size < 31 {
                eprintln!("Warning: unitig size {} is less than k-mer size", size);
                continue;
            }

            let read_size = size.div_ceil(4);
            let mut tig_buffer = vec![0; read_size];
            tigs_file
                .read_exact(&mut tig_buffer)
                .expect("Failed to read tig");

            let tig = vec2str(&tig_buffer, size);

            for file_id in file_ids {
                if wanted_ids.contains(file_id) {
                    if let Some(writer) = writers.get_mut(file_id) {
                        writeln!(writer, ">")?;
                        writeln!(writer, "{}", tig)?;
                    }
                }
            }
            total_unitigs += 1;
        }
    }

    // Flush all writers at the end
    for (_, writer) in writers {
        writer.finish()?;
    }
    println!(
        "Decompression complete: {} unitigs written across {} CIDs",
        total_unitigs, total_cids
    );
    Ok(())
}
