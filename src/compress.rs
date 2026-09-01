use ggcat_api::{
    register_channel_output, unregister_channel_output, ColorIndexType, DnaSequence,
    DynamicSequencesStream, ExtraElaboration, GGCATConfig, GGCATInstance, GeneralSequenceBlockData,
    SequenceInfo, SequencesReader,
};
use ggcat_colors::colors_manager::ColorMapReader;
use ggcat_colors::storage::deserializer::ColorsDeserializer;
use ggcat_colors::DefaultColorsSerializer;
use parking_lot::Mutex as FastMutex;
use rayon::prelude::*;
use rayon::slice::ParallelSliceMut;
use std::cmp::Ordering as CmpOrdering;
use std::collections::BinaryHeap;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Result, Seek, SeekFrom, Write};
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};
use std::sync::{mpsc, Arc, Mutex};
use std::thread;
use std::time::Instant;
use tempfile::Builder as TempBuilder;
use zstd::Encoder;

use crate::packed_tigs::{PackedTigsReader, PackedTigsWriter};
use crate::records::{SimplitigBatch, SimplitigRecord};
use crate::utils::{Convert, Converter};

const IO_BUFFER_CAPACITY: usize = 16 * 1024 * 1024;
const ENCODED_SEQ_BUFFER_TARGET: usize = 4 * 1024 * 1024;
const PAR_SORT_THRESHOLD: usize = 200_000;
const COLOR_RECORD_BATCH_SIZE: usize = 65_536;
const COLOR_CHUNK_TARGET_BYTES: usize = 256 * 1024 * 1024;
const GROUP_SORT_SPILL_BYTES: usize = 128 * 1024 * 1024;
const GROUP_WORK_BATCH_GROUPS: usize = 2048;
const GROUP_IN_MEMORY_MAX_SEQ_BYTES: usize = 2 * 1024 * 1024;
const MIN_MEMORY_BUDGET_BYTES: usize = 256 * 1024 * 1024;
const MIN_GROUP_BATCH_BYTES: usize = 16 * 1024 * 1024;
const MAX_GROUP_BATCH_BYTES: usize = 512 * 1024 * 1024;
const MIN_SORT_CHUNK_BYTES: usize = 32 * 1024 * 1024;
const DATASETS_PER_CID_PARTITION: usize = 512;
const MAX_CID_PARTITIONS: usize = 256;
const CID_PARTITION_BUFFER_BYTES: usize = 256 * 1024;
const CID_TRANSPOSE_BLOCK_BYTES: usize = 64 * 1024;
const COLOR_RUN_BUFFER_BYTES: usize = 64 * 1024;
const COLOR_MERGE_FAN_IN: usize = 64;
const BUCKET_SIZES_MAGIC: &[u8; 4] = b"KSB2";
const BUCKET_SIZES_COMPACT_MAGIC: &[u8; 4] = b"KSB3";
const POSITIONS_MAGIC: &[u8; 4] = b"KPS2";
const POSITIONS_COMPACT_MAGIC: &[u8; 4] = b"KPS3";
const ID_TO_CID_MAGIC: &[u8; 4] = b"KIC2";
pub(crate) const DATASET_TO_CID_MAGIC: &[u8; 4] = b"KDI4";
pub(crate) const DATASET_TO_CID_FLAG_ABUNDANCE: u8 = 1;
pub(crate) const DATASET_TO_CID_FLAG_ADAPTIVE: u8 = 2;
pub(crate) const DATASET_TO_CID_FLAG_MAJORITY: u8 = 4;
pub(crate) const DATASET_POSTING_DELTA: u8 = 0;
pub(crate) const DATASET_POSTING_BITMAP: u8 = 1;
pub(crate) const DATASET_POSTING_MAJORITY_XOR: u8 = 2;
const MAX_DATASET_POSTING_BITMAP_BYTES: u64 = 64 * 1024 * 1024;
pub(crate) const DATASET_TO_CID_FILE: &str = "dataset_to_cid.bin";
pub(crate) const FILENAME_INDEX_MAGIC: &[u8; 4] = b"KFN4";
pub(crate) const FILENAME_INDEX_FILE: &str = "filenames.idx";
pub(crate) const POSITIONS_INDEX_MAGIC: &[u8; 4] = b"KPX4";
pub(crate) const POSITIONS_INDEX_FILE: &str = "positions_kloe.idx";
pub(crate) const BUCKET_SIZES_INDEX_MAGIC: &[u8; 4] = b"KSX4";
pub(crate) const BUCKET_SIZES_INDEX_FILE: &str = "bucket_sizes.idx";
pub(crate) const TIGS_INDEX_FILE: &str = "tigs_kloe.idx";
pub(crate) const ABUNDANCE_BASE_FILE: &str = "abundance_base.bin";
pub(crate) const ABUNDANCE_BASE_INDEX_FILE: &str = "abundance_base.idx";
pub(crate) const ARCHIVE_MANIFEST_FILE: &str = "manifest.kloe";
const ARCHIVE_MANIFEST_MAGIC: &[u8; 4] = b"KLM4";
pub(crate) const CID_TO_DATASET_MAGIC: &[u8; 4] = b"KCD2";
const CID_TO_DATASET_COMPACT_MAGIC: &[u8; 4] = b"KCD3";
const CID_TO_DATASET_ABUNDANCE_MAGIC: &[u8; 4] = b"KCD4";
const CID_TO_DATASET_COLUMNAR_ABUNDANCE_MAGIC: &[u8; 4] = b"KCD5";
const CID_TO_DATASET_HYBRID_ABUNDANCE_MAGIC: &[u8; 4] = b"KCD6";
const CID_TO_DATASET_SPLIT_ABUNDANCE_MAGIC: &[u8; 4] = b"KCD7";
const CID_TO_DATASET_ADAPTIVE_ABUNDANCE_MAGIC: &[u8; 4] = b"KCD8";
const ABUNDANCE_CODEC_SHIFT: u32 = 62;
const ABUNDANCE_LENGTH_MASK: u64 = (1u64 << ABUNDANCE_CODEC_SHIFT) - 1;
const ABUNDANCE_CODEC_RAW: u8 = 0;
const ABUNDANCE_CODEC_GROUP_DELTA: u8 = 1;
const ABUNDANCE_CODEC_GROUP_FIRST: u8 = 2;
pub(crate) const CID_TO_DATASET_FILE: &str = "cid_to_dataset_id.bin";
const ABUNDANCE_CODES_PER_DATASET: usize = 256;
const BUCKET_SIZE_BLOCK_MAX_GROUPS: usize = 8_192;
const BUCKET_SIZE_BLOCK_MAX_UNCOMPRESSED_BYTES: usize = 8 * 1024 * 1024;
const POSITION_INDEX_STRIDE: usize = 256;
const SIZE_INDEX_STRIDE: usize = 256;

#[derive(Debug)]
struct CidDatasetSidecarBlock {
    first_group: usize,
    group_count: usize,
    offsets_file_offset: u64,
    data_offset: u64,
    data_len: usize,
    abundance_data_offset: Option<u64>,
    abundance_data_len: usize,
    abundance_codec: u8,
}

#[derive(Debug)]
struct CidDatasetSidecarCache {
    file: File,
    block_index: usize,
    ranges: Vec<(usize, usize)>,
    data: Vec<u8>,
    abundance_ranges: Vec<(usize, usize)>,
    abundance_codes: Vec<u8>,
    bitmap_modes: Vec<bool>,
    decoded_ids: Vec<u32>,
}

#[derive(Debug)]
pub(crate) struct CidDatasetSidecar {
    blocks: Vec<CidDatasetSidecarBlock>,
    group_count: usize,
    length_prefixed: bool,
    columnar_abundance: bool,
    hybrid_abundance: bool,
    split_abundance: bool,
    dataset_count: Option<usize>,
    abundance_log_base: Option<f64>,
    cache: Mutex<CidDatasetSidecarCache>,
}

struct CidDatasetSidecarReader<'a> {
    sidecar: &'a CidDatasetSidecar,
    file: File,
    block_index: usize,
    ranges: Vec<(usize, usize)>,
    data: Vec<u8>,
    abundance_ranges: Vec<(usize, usize)>,
    abundance_codes: Vec<u8>,
    bitmap_modes: Vec<bool>,
}

impl CidDatasetSidecar {
    pub(crate) fn open(path: &Path) -> io::Result<Self> {
        let mut file = File::open(path)?;
        let mut magic = [0u8; 4];
        file.read_exact(&mut magic)?;
        let (
            length_prefixed,
            columnar_abundance,
            hybrid_abundance,
            split_abundance,
            adaptive_abundance,
            abundance_log_base,
        ) = if &magic == CID_TO_DATASET_MAGIC {
            (false, false, false, false, false, None)
        } else if &magic == CID_TO_DATASET_COMPACT_MAGIC {
            (true, false, false, false, false, None)
        } else if &magic == CID_TO_DATASET_ABUNDANCE_MAGIC
            || &magic == CID_TO_DATASET_COLUMNAR_ABUNDANCE_MAGIC
            || &magic == CID_TO_DATASET_HYBRID_ABUNDANCE_MAGIC
            || &magic == CID_TO_DATASET_SPLIT_ABUNDANCE_MAGIC
            || &magic == CID_TO_DATASET_ADAPTIVE_ABUNDANCE_MAGIC
        {
            let mut base = [0u8; 8];
            file.read_exact(&mut base)?;
            let base = f64::from_le_bytes(base);
            if !base.is_finite() || base <= 1.0 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "invalid CID abundance sidecar '{}': log base must be finite and > 1",
                        path.display()
                    ),
                ));
            }
            (
                true,
                &magic != CID_TO_DATASET_ABUNDANCE_MAGIC,
                &magic == CID_TO_DATASET_HYBRID_ABUNDANCE_MAGIC
                    || &magic == CID_TO_DATASET_SPLIT_ABUNDANCE_MAGIC
                    || &magic == CID_TO_DATASET_ADAPTIVE_ABUNDANCE_MAGIC,
                &magic == CID_TO_DATASET_SPLIT_ABUNDANCE_MAGIC
                    || &magic == CID_TO_DATASET_ADAPTIVE_ABUNDANCE_MAGIC,
                &magic == CID_TO_DATASET_ADAPTIVE_ABUNDANCE_MAGIC,
                Some(base),
            )
        } else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "invalid CID-to-dataset sidecar '{}': bad magic",
                    path.display()
                ),
            ));
        };
        let dataset_count = if hybrid_abundance {
            let mut count = [0u8; 8];
            file.read_exact(&mut count)?;
            let count = usize::try_from(u64::from_le_bytes(count)).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "CID dataset count is too large")
            })?;
            if count == 0 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CID hybrid abundance sidecar has no datasets",
                ));
            }
            Some(count)
        } else {
            None
        };
        let mut blocks = Vec::new();
        let mut group_count = 0usize;
        loop {
            let mut count_bytes = [0u8; 4];
            match file.read_exact(&mut count_bytes) {
                Ok(()) => {}
                Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => break,
                Err(err) => return Err(err),
            }
            let block_groups = u32::from_le_bytes(count_bytes) as usize;
            if block_groups == 0 {
                break;
            }
            let mut len_bytes = [0u8; 8];
            file.read_exact(&mut len_bytes)?;
            let data_len = usize::try_from(u64::from_le_bytes(len_bytes)).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "CID sidecar block is too large")
            })?;
            let (abundance_data_len, abundance_codec) = if split_abundance {
                file.read_exact(&mut len_bytes)?;
                let encoded_len = u64::from_le_bytes(len_bytes);
                let codec = if adaptive_abundance {
                    (encoded_len >> ABUNDANCE_CODEC_SHIFT) as u8
                } else {
                    ABUNDANCE_CODEC_RAW
                };
                if codec > ABUNDANCE_CODEC_GROUP_FIRST {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "CID abundance block uses an unsupported codec",
                    ));
                }
                let data_len = if adaptive_abundance {
                    encoded_len & ABUNDANCE_LENGTH_MASK
                } else {
                    encoded_len
                };
                let data_len = usize::try_from(data_len).map_err(|_| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "CID abundance block is too large",
                    )
                })?;
                (data_len, codec)
            } else {
                (0, ABUNDANCE_CODEC_RAW)
            };
            let offsets_file_offset = file.stream_position()?;
            let offsets_bytes = if length_prefixed {
                0
            } else {
                (block_groups + 1)
                    .checked_mul(std::mem::size_of::<u32>())
                    .ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidData, "CID sidecar offsets overflow")
                    })?
            };
            let data_offset = offsets_file_offset
                .checked_add(offsets_bytes as u64)
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "CID sidecar offset overflow")
                })?;
            let abundance_data_offset = if split_abundance {
                Some(data_offset.checked_add(data_len as u64).ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "CID sidecar offset overflow")
                })?)
            } else {
                None
            };
            let next_block = abundance_data_offset
                .unwrap_or(data_offset)
                .checked_add(if split_abundance {
                    abundance_data_len as u64
                } else {
                    data_len as u64
                })
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "CID sidecar offset overflow")
                })?;
            blocks.push(CidDatasetSidecarBlock {
                first_group: group_count,
                group_count: block_groups,
                offsets_file_offset,
                data_offset,
                data_len,
                abundance_data_offset,
                abundance_data_len,
                abundance_codec,
            });
            group_count = group_count.checked_add(block_groups).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CID sidecar group count overflow",
                )
            })?;
            file.seek(SeekFrom::Start(next_block))?;
        }
        Ok(Self {
            blocks,
            group_count,
            length_prefixed,
            columnar_abundance,
            hybrid_abundance,
            split_abundance,
            dataset_count,
            abundance_log_base,
            cache: Mutex::new(CidDatasetSidecarCache {
                file: File::open(path)?,
                block_index: usize::MAX,
                ranges: Vec::new(),
                data: Vec::new(),
                abundance_ranges: Vec::new(),
                abundance_codes: Vec::new(),
                bitmap_modes: Vec::new(),
                decoded_ids: Vec::new(),
            }),
        })
    }

    pub(crate) fn len(&self) -> usize {
        self.group_count
    }

    pub(crate) fn abundance_log_base(&self) -> Option<f64> {
        self.abundance_log_base
    }

    pub(crate) fn visit_groups(
        &self,
        mut visit: impl FnMut(usize, &[u32]) -> io::Result<()>,
    ) -> io::Result<()> {
        let mut reader = self.reader()?;
        let mut dataset_ids = Vec::new();
        for group in 0..self.group_count {
            reader.load_group_into(group, &mut dataset_ids, 0)?;
            visit(group, &dataset_ids)?;
        }
        Ok(())
    }

    pub(crate) fn visit_groups_with_abundance(
        &self,
        mut visit: impl FnMut(usize, &[u32], Option<&[u8]>) -> io::Result<()>,
    ) -> io::Result<()> {
        let mut reader = self.reader()?;
        let mut dataset_ids = Vec::new();
        let mut abundance_codes = Vec::new();
        for group in 0..self.group_count {
            reader.load_group_with_abundance_into(
                group,
                &mut dataset_ids,
                &mut abundance_codes,
                0,
            )?;
            visit(
                group,
                &dataset_ids,
                self.abundance_log_base.map(|_| abundance_codes.as_slice()),
            )?;
        }
        Ok(())
    }

    fn reader(&self) -> io::Result<CidDatasetSidecarReader<'_>> {
        let file = self
            .cache
            .lock()
            .expect("CID sidecar cache lock poisoned")
            .file
            .try_clone()?;
        Ok(CidDatasetSidecarReader {
            sidecar: self,
            file,
            block_index: usize::MAX,
            ranges: Vec::new(),
            data: Vec::new(),
            abundance_ranges: Vec::new(),
            abundance_codes: Vec::new(),
            bitmap_modes: Vec::new(),
        })
    }

    fn merge_group_into(
        &self,
        group: usize,
        target: &mut Vec<u32>,
        dataset_offset: u32,
    ) -> io::Result<()> {
        let block_index = self
            .blocks
            .partition_point(|block| block.first_group <= group)
            .checked_sub(1)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CID sidecar group is out of range",
                )
            })?;
        let block = self.blocks.get(block_index).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "CID sidecar group is out of range",
            )
        })?;
        let local = group - block.first_group;
        if local >= block.group_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "CID sidecar group is outside its block",
            ));
        }
        let mut cache = self.cache.lock().expect("CID sidecar cache lock poisoned");
        if cache.block_index != block_index {
            let CidDatasetSidecarCache {
                file,
                ranges,
                data,
                abundance_ranges,
                abundance_codes,
                bitmap_modes,
                ..
            } = &mut *cache;
            load_cid_sidecar_block(
                file,
                block,
                self.length_prefixed,
                self.columnar_abundance,
                self.hybrid_abundance,
                self.split_abundance,
                self.dataset_count,
                ranges,
                data,
                abundance_ranges,
                abundance_codes,
                bitmap_modes,
            )?;
            cache.block_index = block_index;
        }
        let (start, end) = cache.ranges[local];
        let bitmap_mode = cache.bitmap_modes.get(local).copied().unwrap_or(false);
        let mut ids = std::mem::take(&mut cache.decoded_ids);
        ids.clear();
        let encoded = &cache.data[start..end];
        if bitmap_mode {
            let dataset_count = self.dataset_count.expect("hybrid sidecar dataset count");
            for (byte_index, &byte) in encoded.iter().enumerate() {
                for bit in 0..8 {
                    if byte & (1 << bit) != 0 {
                        let dataset = byte_index * 8 + bit;
                        if dataset >= dataset_count {
                            return Err(io::Error::new(
                                io::ErrorKind::InvalidData,
                                "CID bitmap contains an out-of-range dataset",
                            ));
                        }
                        ids.push(u32::try_from(dataset + 1).map_err(|_| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                "CID sidecar dataset ID overflow",
                            )
                        })?);
                    }
                }
            }
        } else {
            let mut cursor = 0usize;
            let mut previous = 0u64;
            while cursor < encoded.len() {
                let delta = read_varint_field(encoded, &mut cursor, "CID dataset delta")?;
                previous = previous.checked_add(delta).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "CID sidecar dataset ID overflow",
                    )
                })?;
                let one_based = previous.checked_add(1).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "CID sidecar dataset ID overflow",
                    )
                })?;
                let id = u32::try_from(one_based).map_err(|_| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "CID sidecar dataset ID overflow",
                    )
                })?;
                ids.push(id);
                if self.abundance_log_base.is_some() && !self.columnar_abundance {
                    cursor = cursor.checked_add(1).ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidData, "CID abundance cursor overflow")
                    })?;
                    if cursor > encoded.len() {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "CID abundance code is missing after dataset ID",
                        ));
                    }
                }
            }
        }
        if ids.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "CID sidecar group has no dataset membership",
            ));
        }
        let result = merge_sorted_dataset_ids_with_offset(target, &ids, dataset_offset);
        cache.decoded_ids = ids;
        result
    }
}

impl CidDatasetSidecarReader<'_> {
    fn load_group_into(
        &mut self,
        group: usize,
        target: &mut Vec<u32>,
        dataset_offset: u32,
    ) -> io::Result<()> {
        let mut abundance_codes = Vec::new();
        self.load_group_with_abundance_into(group, target, &mut abundance_codes, dataset_offset)
    }

    fn load_group_with_abundance_into(
        &mut self,
        group: usize,
        target: &mut Vec<u32>,
        abundance_codes: &mut Vec<u8>,
        dataset_offset: u32,
    ) -> io::Result<()> {
        let block_index = self
            .sidecar
            .blocks
            .get(self.block_index)
            .filter(|block| {
                group >= block.first_group && group < block.first_group + block.group_count
            })
            .map(|_| self.block_index)
            .or_else(|| {
                self.sidecar
                    .blocks
                    .partition_point(|block| block.first_group <= group)
                    .checked_sub(1)
            })
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CID sidecar group is out of range",
                )
            })?;
        let block = self.sidecar.blocks.get(block_index).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "CID sidecar group is out of range",
            )
        })?;
        let local = group - block.first_group;
        if local >= block.group_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "CID sidecar group is outside its block",
            ));
        }
        if self.block_index != block_index {
            load_cid_sidecar_block(
                &self.file,
                block,
                self.sidecar.length_prefixed,
                self.sidecar.columnar_abundance,
                self.sidecar.hybrid_abundance,
                self.sidecar.split_abundance,
                self.sidecar.dataset_count,
                &mut self.ranges,
                &mut self.data,
                &mut self.abundance_ranges,
                &mut self.abundance_codes,
                &mut self.bitmap_modes,
            )?;
            self.block_index = block_index;
        }

        target.clear();
        abundance_codes.clear();
        if self.sidecar.columnar_abundance {
            let (start, end) = self.abundance_ranges[local];
            abundance_codes.extend_from_slice(&self.abundance_codes[start..end]);
        }
        let (start, end) = self.ranges[local];
        let encoded = &self.data[start..end];
        if self.bitmap_modes.get(local).copied().unwrap_or(false) {
            let dataset_count = self
                .sidecar
                .dataset_count
                .expect("hybrid sidecar dataset count");
            for (byte_index, &byte) in encoded.iter().enumerate() {
                for bit in 0..8 {
                    if byte & (1 << bit) != 0 {
                        let dataset = byte_index * 8 + bit;
                        if dataset >= dataset_count {
                            return Err(io::Error::new(
                                io::ErrorKind::InvalidData,
                                "CID bitmap contains an out-of-range dataset",
                            ));
                        }
                        let one_based = dataset
                            .checked_add(1)
                            .and_then(|id| id.checked_add(dataset_offset as usize))
                            .ok_or_else(|| {
                                io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    "CID sidecar dataset ID overflow",
                                )
                            })?;
                        target.push(u32::try_from(one_based).map_err(|_| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                "CID sidecar dataset ID overflow",
                            )
                        })?);
                    }
                }
            }
        } else {
            let mut cursor = 0usize;
            let mut previous = 0u64;
            while cursor < encoded.len() {
                let delta = read_varint_field(encoded, &mut cursor, "CID dataset delta")?;
                previous = previous.checked_add(delta).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "CID sidecar dataset ID overflow",
                    )
                })?;
                let one_based = previous
                    .checked_add(1)
                    .and_then(|id| id.checked_add(dataset_offset as u64))
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            "CID sidecar dataset ID overflow",
                        )
                    })?;
                target.push(u32::try_from(one_based).map_err(|_| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "CID sidecar dataset ID overflow",
                    )
                })?);
                if self.sidecar.abundance_log_base.is_some() && !self.sidecar.columnar_abundance {
                    let code = *encoded.get(cursor).ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            "CID abundance code is missing after dataset ID",
                        )
                    })?;
                    cursor += 1;
                    abundance_codes.push(code);
                }
            }
        }
        if target.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "CID sidecar group has no dataset membership",
            ));
        }
        if self.sidecar.abundance_log_base.is_some() && abundance_codes.len() != target.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "CID abundance count differs from dataset membership count",
            ));
        }
        Ok(())
    }
}

fn load_cid_sidecar_block(
    file: &File,
    block: &CidDatasetSidecarBlock,
    length_prefixed: bool,
    columnar_abundance: bool,
    hybrid_abundance: bool,
    split_abundance: bool,
    dataset_count: Option<usize>,
    ranges: &mut Vec<(usize, usize)>,
    data: &mut Vec<u8>,
    abundance_ranges: &mut Vec<(usize, usize)>,
    abundance_codes: &mut Vec<u8>,
    bitmap_modes: &mut Vec<bool>,
) -> io::Result<()> {
    let mut compressed = vec![0u8; block.data_len];
    file.read_exact_at(&mut compressed, block.data_offset)?;
    let decoded = zstd::decode_all(compressed.as_slice())?;
    let decoded_abundance = if split_abundance {
        let mut compressed = vec![0u8; block.abundance_data_len];
        file.read_exact_at(
            &mut compressed,
            block
                .abundance_data_offset
                .expect("split abundance block offset"),
        )?;
        Some(zstd::decode_all(compressed.as_slice())?)
    } else {
        None
    };
    ranges.clear();
    ranges.reserve(block.group_count);
    abundance_ranges.clear();
    abundance_codes.clear();
    bitmap_modes.clear();

    if columnar_abundance {
        let mut cursor = 0usize;
        let mut counts = Vec::with_capacity(block.group_count);
        let mut total_memberships = 0usize;
        for _ in 0..block.group_count {
            let descriptor = read_varint_field(
                &decoded,
                &mut cursor,
                "CID abundance group membership count",
            )?;
            let bitmap_mode = hybrid_abundance && descriptor & 1 != 0;
            let encoded_count = if hybrid_abundance {
                descriptor >> 1
            } else {
                descriptor
            };
            let count = usize::try_from(encoded_count).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "CID sidecar group is too large")
            })?;
            if count == 0 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CID sidecar group has no dataset membership",
                ));
            }
            total_memberships = total_memberships.checked_add(count).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CID abundance membership count overflow",
                )
            })?;
            counts.push(count);
            bitmap_modes.push(bitmap_mode);
        }
        let ids_start = cursor;
        let bitmap_bytes = dataset_count.map(|count| count.div_ceil(8)).unwrap_or(0);
        let mut membership_cursor = 0usize;
        for (&count, &bitmap_mode) in counts.iter().zip(bitmap_modes.iter()) {
            let group_start = cursor - ids_start;
            if bitmap_mode {
                let end = cursor.checked_add(bitmap_bytes).ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "CID bitmap range overflow")
                })?;
                let bitmap = decoded.get(cursor..end).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "CID bitmap extends past its block payload",
                    )
                })?;
                if bitmap
                    .iter()
                    .map(|byte| byte.count_ones() as usize)
                    .sum::<usize>()
                    != count
                {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "CID bitmap membership count is invalid",
                    ));
                }
                cursor = end;
            } else {
                for _ in 0..count {
                    read_varint_field(&decoded, &mut cursor, "CID dataset delta")?;
                }
            }
            ranges.push((group_start, cursor - ids_start));
            abundance_ranges.push((membership_cursor, membership_cursor + count));
            membership_cursor += count;
        }
        let codes_start = cursor;
        data.clear();
        data.extend_from_slice(&decoded[ids_start..codes_start]);
        if let Some(decoded_abundance) = decoded_abundance.as_deref() {
            if codes_start != decoded.len() || decoded_abundance.len() != total_memberships {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CID split abundance block has an invalid stream length",
                ));
            }
            match block.abundance_codec {
                ABUNDANCE_CODEC_RAW => abundance_codes.extend_from_slice(decoded_abundance),
                ABUNDANCE_CODEC_GROUP_DELTA => {
                    let mut code_cursor = 0usize;
                    for &count in &counts {
                        let mut previous = 0u8;
                        for membership in 0..count {
                            let encoded = decoded_abundance[code_cursor];
                            let code = if membership == 0 {
                                encoded
                            } else {
                                previous.wrapping_add(encoded)
                            };
                            abundance_codes.push(code);
                            previous = code;
                            code_cursor += 1;
                        }
                    }
                    debug_assert_eq!(code_cursor, total_memberships);
                }
                ABUNDANCE_CODEC_GROUP_FIRST => {
                    let mut code_cursor = 0usize;
                    for &count in &counts {
                        let first = decoded_abundance[code_cursor];
                        abundance_codes.push(first);
                        code_cursor += 1;
                        for _ in 1..count {
                            abundance_codes
                                .push(first.wrapping_add(decoded_abundance[code_cursor]));
                            code_cursor += 1;
                        }
                    }
                    debug_assert_eq!(code_cursor, total_memberships);
                }
                _ => unreachable!("abundance codec validated while opening sidecar"),
            }
        } else {
            let expected_end = codes_start.checked_add(total_memberships).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "CID abundance payload overflow")
            })?;
            if expected_end != decoded.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CID columnar abundance block has an invalid code stream length",
                ));
            }
            abundance_codes.extend_from_slice(&decoded[codes_start..]);
        }
        debug_assert_eq!(membership_cursor, abundance_codes.len());
    } else if length_prefixed {
        *data = decoded;
        let mut cursor = 0usize;
        for _ in 0..block.group_count {
            let group_len = usize::try_from(read_varint_field(
                data,
                &mut cursor,
                "CID dataset group length",
            )?)
            .map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "CID sidecar group is too large")
            })?;
            let end = cursor.checked_add(group_len).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CID sidecar group offset overflow",
                )
            })?;
            if end > data.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CID sidecar group extends past its block payload",
                ));
            }
            ranges.push((cursor, end));
            cursor = end;
        }
        if cursor != data.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "CID sidecar block has trailing group data",
            ));
        }
    } else {
        *data = decoded;
        let mut offset_bytes = vec![0u8; (block.group_count + 1) * std::mem::size_of::<u32>()];
        file.read_exact_at(&mut offset_bytes, block.offsets_file_offset)?;
        let offsets = offset_bytes
            .chunks_exact(4)
            .map(|bytes| u32::from_le_bytes(bytes.try_into().expect("four-byte offset")) as usize)
            .collect::<Vec<_>>();
        if offsets.windows(2).any(|range| range[0] > range[1])
            || offsets.last().copied() != Some(data.len())
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "CID sidecar block offsets are invalid for its payload",
            ));
        }
        ranges.extend(offsets.windows(2).map(|range| (range[0], range[1])));
    }
    Ok(())
}

pub(crate) struct CidDatasetSidecarWriter {
    file: BufWriter<File>,
    abundance_log_base: Option<f64>,
    dataset_count: Option<usize>,
    block_groups: usize,
    block_offsets: Vec<u32>,
    block_data: Vec<u8>,
    block_membership_counts: Vec<u32>,
    block_bitmap_modes: Vec<bool>,
    block_abundance_codes: Vec<u8>,
}

impl CidDatasetSidecarWriter {
    pub(crate) fn create(path: &Path) -> io::Result<Self> {
        Self::create_with_abundance(path, None, None)
    }

    pub(crate) fn create_with_abundance(
        path: &Path,
        abundance_log_base: Option<f64>,
        dataset_count: Option<usize>,
    ) -> io::Result<Self> {
        let mut file = BufWriter::with_capacity(IO_BUFFER_CAPACITY, File::create(path)?);
        if let Some(base) = abundance_log_base {
            if !base.is_finite() || base <= 1.0 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "abundance log base must be finite and > 1",
                ));
            }
            let dataset_count = dataset_count.filter(|&count| count > 0).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "abundance sidecar requires a nonzero dataset count",
                )
            })?;
            file.write_all(CID_TO_DATASET_ADAPTIVE_ABUNDANCE_MAGIC)?;
            file.write_all(&base.to_le_bytes())?;
            file.write_all(&(dataset_count as u64).to_le_bytes())?;
        } else {
            file.write_all(CID_TO_DATASET_COMPACT_MAGIC)?;
        }
        Ok(Self {
            file,
            abundance_log_base,
            dataset_count,
            block_groups: 0,
            block_offsets: vec![0],
            block_data: Vec::new(),
            block_membership_counts: Vec::new(),
            block_bitmap_modes: Vec::new(),
            block_abundance_codes: Vec::new(),
        })
    }

    pub(crate) fn append_zero_based(&mut self, dataset_ids: &[usize]) -> io::Result<()> {
        self.append_zero_based_with_abundance(dataset_ids, None)
    }

    pub(crate) fn append_zero_based_with_abundance(
        &mut self,
        dataset_ids: &[usize],
        abundance_codes: Option<&[u8]>,
    ) -> io::Result<()> {
        if dataset_ids.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "cannot store an empty CID dataset set",
            ));
        }
        match (self.abundance_log_base, abundance_codes) {
            (Some(_), Some(codes)) if codes.len() == dataset_ids.len() => {}
            (Some(_), Some(_)) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CID abundance count differs from dataset membership count",
                ))
            }
            (Some(_), None) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "abundance-enabled archive is missing CID abundance codes",
                ))
            }
            (None, Some(_)) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "cannot write abundance codes to an abundance-disabled archive",
                ))
            }
            (None, None) => {}
        }
        let group_start = self.block_data.len();
        let mut previous = 0usize;
        for (index, &dataset) in dataset_ids.iter().enumerate() {
            if index > 0 && dataset <= previous {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CID dataset IDs are not strictly increasing",
                ));
            }
            write_varint_u64((dataset - previous) as u64, &mut self.block_data);
            previous = dataset;
        }
        if let Some(dataset_count) = self.dataset_count {
            if dataset_ids
                .last()
                .is_some_and(|&dataset| dataset >= dataset_count)
            {
                self.block_data.truncate(group_start);
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CID sidecar dataset ID exceeds its declared dataset count",
                ));
            }
            let bitmap_bytes = dataset_count.div_ceil(8);
            let delta_bytes = self.block_data.len() - group_start;
            if bitmap_bytes < delta_bytes {
                self.block_data.truncate(group_start);
                self.block_data.resize(group_start + bitmap_bytes, 0);
                for &dataset in dataset_ids {
                    self.block_data[group_start + dataset / 8] |= 1 << (dataset % 8);
                }
                self.block_bitmap_modes.push(true);
            } else {
                self.block_bitmap_modes.push(false);
            }
        }
        if let Some(codes) = abundance_codes {
            self.block_abundance_codes.extend_from_slice(codes);
        }
        self.finish_group(dataset_ids.len())
    }

    pub(crate) fn append_one_based(&mut self, dataset_ids: &[u32]) -> io::Result<()> {
        if self.abundance_log_base.is_some() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "abundance-enabled sidecar requires abundance codes",
            ));
        }
        if dataset_ids.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "cannot store an empty CID dataset set",
            ));
        }
        let mut previous = 0u32;
        let mut first = true;
        for &dataset in dataset_ids {
            let zero_based = dataset.checked_sub(1).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "dataset ID zero is invalid")
            })?;
            if !first && zero_based <= previous {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CID dataset IDs are not strictly increasing",
                ));
            }
            write_varint_u64((zero_based - previous) as u64, &mut self.block_data);
            previous = zero_based;
            first = false;
        }
        self.finish_group(dataset_ids.len())
    }

    fn finish_group(&mut self, membership_count: usize) -> io::Result<()> {
        self.block_offsets
            .push(u32::try_from(self.block_data.len()).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "CID sidecar block exceeds 4GiB")
            })?);
        if self.abundance_log_base.is_some() {
            self.block_membership_counts
                .push(u32::try_from(membership_count).map_err(|_| {
                    io::Error::new(io::ErrorKind::InvalidData, "CID sidecar group is too large")
                })?);
            debug_assert_eq!(
                self.block_bitmap_modes.len(),
                self.block_membership_counts.len()
            );
        }
        self.block_groups += 1;
        if self.block_groups >= BUCKET_SIZE_BLOCK_MAX_GROUPS
            || self
                .block_data
                .len()
                .saturating_add(self.block_abundance_codes.len())
                >= BUCKET_SIZE_BLOCK_MAX_UNCOMPRESSED_BYTES
        {
            self.flush_block()?;
        }
        Ok(())
    }

    fn flush_block(&mut self) -> io::Result<()> {
        if self.block_groups == 0 {
            return Ok(());
        }
        let mut encoded = Vec::with_capacity(
            self.block_data
                .len()
                .saturating_add(self.block_abundance_codes.len())
                .saturating_add(self.block_groups.saturating_mul(2)),
        );
        let abundance_enabled = self.abundance_log_base.is_some();
        if abundance_enabled {
            debug_assert_eq!(self.block_membership_counts.len(), self.block_groups);
            for (&count, &bitmap_mode) in self
                .block_membership_counts
                .iter()
                .zip(self.block_bitmap_modes.iter())
            {
                write_varint_u64(((count as u64) << 1) | u64::from(bitmap_mode), &mut encoded);
            }
            encoded.extend_from_slice(&self.block_data);
        } else {
            for offsets in self.block_offsets.windows(2) {
                let start = offsets[0] as usize;
                let end = offsets[1] as usize;
                write_varint_u64((end - start) as u64, &mut encoded);
                encoded.extend_from_slice(&self.block_data[start..end]);
            }
        }
        let compressed = zstd::encode_all(encoded.as_slice(), 1)?;
        let (abundance_compressed, abundance_codec) = if abundance_enabled {
            let raw = zstd::encode_all(self.block_abundance_codes.as_slice(), 1)?;
            let mut transformed_codes = Vec::with_capacity(self.block_abundance_codes.len());
            let mut cursor = 0usize;
            for &count in &self.block_membership_counts {
                let mut previous = 0u8;
                for membership in 0..count as usize {
                    let code = self.block_abundance_codes[cursor];
                    transformed_codes.push(if membership == 0 {
                        code
                    } else {
                        code.wrapping_sub(previous)
                    });
                    previous = code;
                    cursor += 1;
                }
            }
            debug_assert_eq!(cursor, self.block_abundance_codes.len());
            let delta = zstd::encode_all(transformed_codes.as_slice(), 1)?;

            transformed_codes.clear();
            cursor = 0;
            for &count in &self.block_membership_counts {
                let first = self.block_abundance_codes[cursor];
                transformed_codes.push(first);
                cursor += 1;
                for _ in 1..count {
                    transformed_codes.push(self.block_abundance_codes[cursor].wrapping_sub(first));
                    cursor += 1;
                }
            }
            debug_assert_eq!(cursor, self.block_abundance_codes.len());
            let first = zstd::encode_all(transformed_codes.as_slice(), 1)?;

            let candidates = [
                (raw, ABUNDANCE_CODEC_RAW),
                (delta, ABUNDANCE_CODEC_GROUP_DELTA),
                (first, ABUNDANCE_CODEC_GROUP_FIRST),
            ];
            let (compressed, codec) = candidates
                .into_iter()
                .min_by_key(|(compressed, _)| compressed.len())
                .expect("three abundance codecs");
            (Some(compressed), codec)
        } else {
            (None, ABUNDANCE_CODEC_RAW)
        };
        self.file
            .write_all(&(self.block_groups as u32).to_le_bytes())?;
        self.file
            .write_all(&(compressed.len() as u64).to_le_bytes())?;
        if let Some(abundance_compressed) = abundance_compressed.as_ref() {
            let encoded_len = u64::try_from(abundance_compressed.len()).map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CID abundance block is too large",
                )
            })?;
            if encoded_len > ABUNDANCE_LENGTH_MASK {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CID abundance block exceeds its format limit",
                ));
            }
            self.file.write_all(
                &(encoded_len | ((abundance_codec as u64) << ABUNDANCE_CODEC_SHIFT)).to_le_bytes(),
            )?;
        }
        self.file.write_all(&compressed)?;
        if let Some(abundance_compressed) = abundance_compressed {
            self.file.write_all(&abundance_compressed)?;
        }
        self.block_groups = 0;
        self.block_offsets.clear();
        self.block_offsets.push(0);
        self.block_data.clear();
        self.block_membership_counts.clear();
        self.block_bitmap_modes.clear();
        self.block_abundance_codes.clear();
        Ok(())
    }

    pub(crate) fn finish(mut self) -> io::Result<()> {
        self.flush_block()?;
        self.file.write_all(&0u32.to_le_bytes())?;
        self.file.flush()
    }
}

pub fn compress(
    output_dir: &String,
    input_fof: &String,
    threads: usize,
    k: usize,
    m: usize,
    partition_power: u32,
    verify_kmers: bool,
    skip_sort: bool,
    use_unitigs: bool,
    use_matchtigs: bool,
    use_eulertigs: bool,
) -> Result<()> {
    compress_with_ggcat(
        output_dir,
        input_fof,
        threads,
        k,
        m,
        partition_power,
        verify_kmers,
        skip_sort,
        use_unitigs,
        use_matchtigs,
        use_eulertigs,
        GgcatCompressionConfig::default(),
    )
}

#[derive(Clone, Debug)]
pub struct GgcatCompressionConfig {
    pub memory_gb: usize,
    pub temp_dir: String,
    pub abundance_log_base: Option<f64>,
}

#[derive(Clone, Copy, Debug)]
struct CompressionMemoryBudget {
    total_bytes: usize,
}

impl CompressionMemoryBudget {
    fn from_gb(memory_gb: usize) -> Self {
        Self {
            total_bytes: memory_gb
                .max(1)
                .saturating_mul(1024 * 1024 * 1024)
                .max(MIN_MEMORY_BUDGET_BYTES),
        }
    }

    fn group_batch_bytes(self) -> usize {
        (self.total_bytes / 8).clamp(MIN_GROUP_BATCH_BYTES, MAX_GROUP_BATCH_BYTES)
    }

    fn color_chunk_bytes(self) -> usize {
        (self.total_bytes / 16).clamp(MIN_SORT_CHUNK_BYTES, COLOR_CHUNK_TARGET_BYTES)
    }

    fn subset_window_sequence_bytes(self) -> usize {
        (self.total_bytes / 16).clamp(8 * 1024 * 1024, 1024 * 1024 * 1024)
    }

    fn subset_window_count(self, dataset_count: usize) -> usize {
        let worst_case_subset_bytes = dataset_count
            .saturating_mul(std::mem::size_of::<u32>())
            .saturating_add(96)
            .max(1);
        (self.total_bytes / 8 / worst_case_subset_bytes).max(1)
    }

    fn ggcat_memory_gb(self) -> f64 {
        ((self.total_bytes as f64) / (1024.0 * 1024.0 * 1024.0) * 0.75).max(0.25)
    }
}

impl Default for GgcatCompressionConfig {
    fn default() -> Self {
        Self {
            memory_gb: 8,
            temp_dir: String::new(),
            abundance_log_base: None,
        }
    }
}

#[derive(Debug)]
struct SortedColorRecord {
    subset: ColorIndexType,
    seq: Vec<u8>,
}

#[derive(Debug)]
struct ChunkHeapItem {
    chunk_index: usize,
    record: SortedColorRecord,
}

impl Eq for ChunkHeapItem {}

impl PartialEq for ChunkHeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.record.subset == other.record.subset && self.chunk_index == other.chunk_index
    }
}

impl Ord for ChunkHeapItem {
    fn cmp(&self, other: &Self) -> CmpOrdering {
        other
            .record
            .subset
            .cmp(&self.record.subset)
            .then_with(|| other.chunk_index.cmp(&self.chunk_index))
    }
}

impl PartialOrd for ChunkHeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<CmpOrdering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug)]
struct SeqHeapItem {
    run_idx: usize,
    seq: Vec<u8>,
}

impl Eq for SeqHeapItem {}

impl PartialEq for SeqHeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.seq == other.seq
    }
}

impl Ord for SeqHeapItem {
    fn cmp(&self, other: &Self) -> CmpOrdering {
        other
            .seq
            .len()
            .cmp(&self.seq.len())
            .then_with(|| other.seq.cmp(&self.seq))
    }
}

impl PartialOrd for SeqHeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<CmpOrdering> {
        Some(self.cmp(other))
    }
}

struct SeqRunReader {
    reader: BufReader<File>,
    buf: Vec<u8>,
}

impl SeqRunReader {
    fn new(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        Ok(Self {
            reader: BufReader::new(file),
            buf: Vec::new(),
        })
    }

    fn next_seq(&mut self) -> Result<Option<Vec<u8>>> {
        self.buf.clear();
        let n = self.reader.read_until(b'\n', &mut self.buf)?;
        if n == 0 {
            return Ok(None);
        }
        while self.buf.last() == Some(&b'\n') || self.buf.last() == Some(&b'\r') {
            self.buf.pop();
        }
        Ok(Some(self.buf.clone()))
    }
}

fn read_input_fof_filenames(input_fof: &str) -> Result<Vec<String>> {
    let input_fof_reader = BufReader::new(File::open(input_fof)?);
    let mut filenames = Vec::new();
    for line_result in input_fof_reader.lines() {
        let line = line_result?;
        let trimmed = line.trim();
        if !trimmed.is_empty() {
            filenames.push(trimmed.to_string());
        }
    }
    Ok(filenames)
}

fn parse_header_abundance(header: &[u8]) -> Option<u64> {
    let header = std::str::from_utf8(header).ok()?;
    for token in header.split_ascii_whitespace() {
        for prefix in ["ka:f:", "ka:i:", "km:f:", "km:i:", "abundance="] {
            if let Some(value) = token.strip_prefix(prefix) {
                let parsed = value.parse::<f64>().ok()?;
                if parsed.is_finite() && parsed >= 0.0 {
                    return Some(parsed.round() as u64);
                }
            }
        }
    }
    None
}

pub(crate) fn encode_abundance(value: u64, log_base: f64) -> u8 {
    if value == 0 {
        return 0;
    }
    let code = 1.0 + (value as f64).ln() / log_base.ln();
    code.round().clamp(1.0, u8::MAX as f64) as u8
}

pub(crate) fn decode_abundance(code: u8, log_base: f64) -> u64 {
    if code == 0 {
        0
    } else {
        log_base.powi(code as i32 - 1).round().max(1.0) as u64
    }
}

struct AbundanceSequencesStream {
    files: Vec<PathBuf>,
    log_base: f64,
}

impl AbundanceSequencesStream {
    fn new(files: &[String], log_base: f64) -> Self {
        Self {
            files: files
                .iter()
                .map(|file| {
                    let path = PathBuf::from(file);
                    fs::canonicalize(&path).unwrap_or(path)
                })
                .collect(),
            log_base,
        }
    }
}

impl DynamicSequencesStream for AbundanceSequencesStream {
    fn read_block(
        &self,
        block: usize,
        _copy_ident_data: bool,
        partial_read_copyback: Option<usize>,
        callback: &mut dyn FnMut(DnaSequence, SequenceInfo),
    ) {
        let path = self
            .files
            .get(block)
            .unwrap_or_else(|| panic!("invalid abundance input block {block}"));
        let color_base = block
            .checked_mul(ABUNDANCE_CODES_PER_DATASET)
            .unwrap_or_else(|| panic!("abundance color index overflow for block {block}"));
        let mut reader = SequencesReader::new();
        reader.process_file_extended(
            path,
            |sequence| {
                let abundance = parse_header_abundance(sequence.ident_data).unwrap_or_else(|| {
                    panic!(
                        "input record in '{}' has no supported abundance field (expected ka:f:, ka:i:, km:f:, or km:i:)",
                        path.display()
                    )
                });
                let code = encode_abundance(abundance, self.log_base);
                let color = u32::try_from(color_base + code as usize)
                    .unwrap_or_else(|_| panic!("abundance color index exceeds u32"));
                callback(sequence, SequenceInfo { color: Some(color) });
            },
            partial_read_copyback,
            true,
            false,
        );
    }

    fn estimated_base_count(&self, block: usize) -> u64 {
        self.files
            .get(block)
            .and_then(|path| fs::metadata(path).ok())
            .map_or(1, |metadata| metadata.len().max(1))
    }
}

pub(crate) fn write_filenames_id_offsets(
    output_dir: &str,
    filenames: &[String],
    id_cid_line_sizes: &[usize],
) -> Result<()> {
    if id_cid_line_sizes.len() != filenames.len()
        && id_cid_line_sizes.len() != filenames.len().saturating_add(1)
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "filenames count ({}) is incompatible with dataset boundary count ({})",
                filenames.len(),
                id_cid_line_sizes.len()
            ),
        ));
    }

    #[derive(Clone, Copy)]
    struct FilenameRecord {
        hash: u64,
        dataset_id: u32,
        name_offset: u64,
        name_len: u32,
    }

    let mut fof_id = BufWriter::new(File::create(output_dir.to_owned() + "filenames_id.txt")?);
    let mut file_offset = 0u64;
    let mut records = Vec::with_capacity(filenames.len());
    let indexed_v4 = id_cid_line_sizes.len() == filenames.len().saturating_add(1);
    for (dataset_id, filename) in filenames.iter().enumerate() {
        let offset = if indexed_v4 {
            0
        } else {
            id_cid_line_sizes[dataset_id]
        };
        let line = format!("{filename}:{offset}\n");
        fof_id.write_all(line.as_bytes())?;
        if dataset_id >= u32::MAX as usize {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "too many datasets for filename index",
            ));
        }
        records.push(FilenameRecord {
            hash: stable_filename_hash(filename.as_bytes()),
            dataset_id: dataset_id as u32,
            name_offset: file_offset,
            name_len: u32::try_from(filename.len()).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "dataset filename is too long")
            })?,
        });
        file_offset = file_offset
            .checked_add(line.len() as u64)
            .ok_or_else(|| io::Error::other("filename table offset overflow"))?;
    }
    fof_id.flush()?;

    // Structural legacy archives do not have dataset-major postings.
    if !indexed_v4 {
        return Ok(());
    }

    // Keep the load factor at or below 75%. The resulting open-addressed table
    // is directly queryable on disk and avoids reading or rebuilding a map at
    // decompression time.
    let minimum_capacity = records.len().saturating_mul(4).div_ceil(3).max(1);
    let capacity = minimum_capacity
        .checked_next_power_of_two()
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "filename hash table capacity overflow",
            )
        })?;
    let mut slots = vec![None; capacity];
    for record in &records {
        let mut slot = record.hash as usize & (capacity - 1);
        loop {
            if slots[slot].is_none() {
                slots[slot] = Some((record.hash, record.dataset_id));
                break;
            }
            slot = (slot + 1) & (capacity - 1);
        }
    }

    let mut index = BufWriter::with_capacity(
        IO_BUFFER_CAPACITY,
        File::create(output_dir.to_owned() + FILENAME_INDEX_FILE)?,
    );
    index.write_all(FILENAME_INDEX_MAGIC)?;
    index.write_all(&(filenames.len() as u64).to_le_bytes())?;
    index.write_all(&(capacity as u64).to_le_bytes())?;
    for &boundary in id_cid_line_sizes {
        index.write_all(&(boundary as u64).to_le_bytes())?;
    }
    for record in &records {
        index.write_all(&record.name_offset.to_le_bytes())?;
        index.write_all(&record.name_len.to_le_bytes())?;
    }
    for slot in slots {
        if let Some((hash, dataset_id)) = slot {
            index.write_all(&hash.to_le_bytes())?;
            index.write_all(&dataset_id.to_le_bytes())?;
        } else {
            index.write_all(&0u64.to_le_bytes())?;
            index.write_all(&u32::MAX.to_le_bytes())?;
        }
    }
    index.flush()?;
    Ok(())
}

pub(crate) fn stable_filename_hash(bytes: &[u8]) -> u64 {
    let mut hash = 0xcbf29ce484222325u64;
    for &byte in bytes {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct MembershipEdge {
    cid: u64,
    dataset: u32,
    abundance: u8,
}

fn write_membership_edge(writer: &mut impl Write, edge: MembershipEdge) -> Result<()> {
    writer.write_all(&edge.cid.to_le_bytes())?;
    writer.write_all(&edge.dataset.to_le_bytes())?;
    writer.write_all(&[edge.abundance, 0, 0, 0])
}

fn read_membership_edge(reader: &mut impl Read) -> Result<Option<MembershipEdge>> {
    let mut bytes = [0u8; 16];
    match reader.read_exact(&mut bytes) {
        Ok(()) => Ok(Some(MembershipEdge {
            cid: u64::from_le_bytes(bytes[..8].try_into().unwrap()),
            dataset: u32::from_le_bytes(bytes[8..12].try_into().unwrap()),
            abundance: bytes[12],
        })),
        Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => Ok(None),
        Err(err) => Err(err),
    }
}

fn flush_membership_run(
    edges: &mut Vec<MembershipEdge>,
    runs: &mut Vec<PathBuf>,
    directory: &Path,
) -> Result<()> {
    if edges.is_empty() {
        return Ok(());
    }
    edges.par_sort_unstable_by_key(|edge| (edge.cid, edge.dataset));
    let path = directory.join(format!("membership_{:06}.run", runs.len()));
    let mut writer = BufWriter::with_capacity(IO_BUFFER_CAPACITY, File::create(&path)?);
    for edge in edges.drain(..) {
        write_membership_edge(&mut writer, edge)?;
    }
    writer.flush()?;
    runs.push(path);
    Ok(())
}

fn merge_membership_runs(
    paths: &[PathBuf],
    output: &Path,
    mut visit: Option<&mut dyn FnMut(MembershipEdge) -> Result<()>>,
) -> Result<()> {
    let mut readers = paths
        .iter()
        .map(|path| BufReader::with_capacity(COLOR_RUN_BUFFER_BYTES, File::open(path).unwrap()))
        .collect::<Vec<_>>();
    let mut heap = BinaryHeap::<std::cmp::Reverse<(u64, u32, u8, usize)>>::new();
    for (index, reader) in readers.iter_mut().enumerate() {
        if let Some(edge) = read_membership_edge(reader)? {
            heap.push(std::cmp::Reverse((
                edge.cid,
                edge.dataset,
                edge.abundance,
                index,
            )));
        }
    }
    let mut writer = if visit.is_none() {
        Some(BufWriter::with_capacity(
            IO_BUFFER_CAPACITY,
            File::create(output)?,
        ))
    } else {
        None
    };
    while let Some(std::cmp::Reverse((cid, dataset, abundance, index))) = heap.pop() {
        let edge = MembershipEdge {
            cid,
            dataset,
            abundance,
        };
        if let Some(visit) = visit.as_deref_mut() {
            visit(edge)?;
        } else {
            write_membership_edge(writer.as_mut().unwrap(), edge)?;
        }
        if let Some(next) = read_membership_edge(&mut readers[index])? {
            heap.push(std::cmp::Reverse((
                next.cid,
                next.dataset,
                next.abundance,
                index,
            )));
        }
    }
    if let Some(writer) = writer.as_mut() {
        writer.flush()?;
    }
    Ok(())
}

pub(crate) fn build_temporary_cid_sidecar_from_dataset_index(
    data_path: &Path,
    index_path: &Path,
    output_path: &Path,
    work_directory: &Path,
    expected_groups: usize,
) -> Result<Arc<CidDatasetSidecar>> {
    let mut index = File::open(index_path)?;
    let mut magic = [0u8; 4];
    index.read_exact(&mut magic)?;
    if &magic != FILENAME_INDEX_MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid dataset directory",
        ));
    }
    let mut value = [0u8; 8];
    index.read_exact(&mut value)?;
    let dataset_count = usize::try_from(u64::from_le_bytes(value))
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "dataset count is too large"))?;
    index.read_exact(&mut value)?; // hash-table capacity
    let mut offsets = vec![0u8; dataset_count.saturating_add(1).saturating_mul(8)];
    index.read_exact(&mut offsets)?;

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
        data.read_exact(&mut value)?;
        Some(f64::from_le_bytes(value))
    } else {
        None
    };
    let adaptive = flags[0] & DATASET_TO_CID_FLAG_ADAPTIVE != 0;
    let cid_count = if adaptive {
        data.read_exact(&mut value)?;
        Some(u64::from_le_bytes(value))
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
        data.read_exact(&mut value)?;
        Some((data.stream_position()?, u64::from_le_bytes(value)))
    } else {
        None
    };

    let run_directory = work_directory.join(format!(
        ".kloe-membership-transpose-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|duration| duration.as_nanos())
            .unwrap_or(0)
    ));
    fs::create_dir_all(&run_directory)?;
    let edge_capacity = (64 * 1024 * 1024 / std::mem::size_of::<MembershipEdge>()).max(1);
    let mut edges = Vec::with_capacity(edge_capacity);
    let mut runs = Vec::new();
    for dataset in 0..dataset_count {
        let offset = u64::from_le_bytes(offsets[dataset * 8..dataset * 8 + 8].try_into().unwrap());
        let end = u64::from_le_bytes(
            offsets[dataset * 8 + 8..dataset * 8 + 16]
                .try_into()
                .unwrap(),
        );
        let compressed_len = end.checked_sub(offset).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "dataset posting offsets are not monotonic",
            )
        })?;
        let mut cid_file = File::open(data_path)?;
        cid_file.seek(SeekFrom::Start(offset))?;
        let (codec, payload_len) = if adaptive {
            if compressed_len < 1 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "truncated adaptive dataset posting frame",
                ));
            }
            let mut codec = [0u8; 1];
            cid_file.read_exact(&mut codec)?;
            (codec[0], compressed_len - 1)
        } else {
            (DATASET_POSTING_DELTA, compressed_len)
        };
        if codec == DATASET_POSTING_BITMAP || codec == DATASET_POSTING_MAJORITY_XOR {
            if abundance_log_base.is_some() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "bitmap postings cannot contain abundance values",
                ));
            }
            let universe = cid_count.ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "bitmap CID universe is missing")
            })?;
            if universe != expected_groups as u64 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "bitmap CID universe does not match archive",
                ));
            }
            let mut decoder = zstd::Decoder::new(cid_file.take(payload_len))?.single_frame();
            let mut reference_decoder = if codec == DATASET_POSTING_MAJORITY_XOR {
                let (reference_offset, reference_len) = majority_frame.ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "majority-XOR posting has no reference bitmap",
                    )
                })?;
                let mut reference = File::open(data_path)?;
                reference.seek(SeekFrom::Start(reference_offset))?;
                Some(zstd::Decoder::new(reference.take(reference_len))?.single_frame())
            } else {
                None
            };
            let mut cid = 0u64;
            while cid < universe {
                let mut byte = [0u8; 1];
                decoder.read_exact(&mut byte)?;
                if let Some(reference) = reference_decoder.as_mut() {
                    let mut base = [0u8; 1];
                    reference.read_exact(&mut base)?;
                    byte[0] ^= base[0];
                }
                for bit in 0..8 {
                    if cid >= universe {
                        break;
                    }
                    if byte[0] & (1u8 << bit) != 0 {
                        edges.push(MembershipEdge {
                            cid,
                            dataset: u32::try_from(dataset + 1).map_err(|_| {
                                io::Error::new(io::ErrorKind::InvalidData, "dataset ID overflow")
                            })?,
                            abundance: 0,
                        });
                        if edges.len() == edge_capacity {
                            flush_membership_run(&mut edges, &mut runs, &run_directory)?;
                        }
                    }
                    cid += 1;
                }
            }
            continue;
        }
        if codec != DATASET_POSTING_DELTA {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "unsupported dataset posting codec",
            ));
        }
        let (mut decoder, mut abundance_decoder) = if abundance_log_base.is_some() {
            if payload_len < 8 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "truncated abundance posting frame",
                ));
            }
            cid_file.read_exact(&mut value)?;
            let cid_len = u64::from_le_bytes(value);
            let abundance_len = payload_len
                .checked_sub(8)
                .and_then(|len| len.checked_sub(cid_len))
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "invalid abundance posting lengths",
                    )
                })?;
            let decoder = zstd::Decoder::new(cid_file.take(cid_len))?.single_frame();
            let mut abundance_file = File::open(data_path)?;
            abundance_file.seek(SeekFrom::Start(offset + u64::from(adaptive) + 8 + cid_len))?;
            let abundance_decoder =
                zstd::Decoder::new(abundance_file.take(abundance_len))?.single_frame();
            (decoder, Some(abundance_decoder))
        } else {
            (
                zstd::Decoder::new(cid_file.take(payload_len))?.single_frame(),
                None,
            )
        };
        let mut cid = 0u64;
        while let Some(delta) = read_optional_varint_u64(&mut decoder)? {
            cid = cid.checked_add(delta).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "dataset posting CID overflow")
            })?;
            if cid >= expected_groups as u64 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "dataset posting CID outside archive",
                ));
            }
            let abundance = if abundance_log_base.is_some() {
                let mut residual = [0u8; 1];
                abundance_decoder
                    .as_mut()
                    .expect("abundance decoder initialized")
                    .read_exact(&mut residual)?;
                residual[0]
            } else {
                0
            };
            edges.push(MembershipEdge {
                cid,
                dataset: u32::try_from(dataset + 1).map_err(|_| {
                    io::Error::new(io::ErrorKind::InvalidData, "dataset ID overflow")
                })?,
                abundance,
            });
            if edges.len() == edge_capacity {
                flush_membership_run(&mut edges, &mut runs, &run_directory)?;
            }
        }
    }
    flush_membership_run(&mut edges, &mut runs, &run_directory)?;
    if runs.is_empty() && expected_groups != 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "dataset index has no memberships",
        ));
    }
    let mut pass = 0usize;
    while runs.len() > COLOR_MERGE_FAN_IN {
        let mut next = Vec::new();
        for (group, paths) in runs.chunks(COLOR_MERGE_FAN_IN).enumerate() {
            let path = run_directory.join(format!("merge_{pass:03}_{group:06}.run"));
            merge_membership_runs(paths, &path, None)?;
            for old in paths {
                fs::remove_file(old)?;
            }
            next.push(path);
        }
        runs = next;
        pass += 1;
    }

    let mut writer = CidDatasetSidecarWriter::create_with_abundance(
        output_path,
        abundance_log_base,
        abundance_log_base.map(|_| dataset_count),
    )?;
    let mut current_cid = 0u64;
    let mut current_ids = Vec::<u32>::new();
    let mut current_codes = Vec::<u8>::new();
    let mut saw_edge = false;
    let mut abundance_bases = if abundance_log_base.is_some() {
        let root = data_path.parent().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "dataset index has no archive directory",
            )
        })?;
        let reader = PackedTigsReader::open_indexed(
            root.join(ABUNDANCE_BASE_FILE),
            root.join(ABUNDANCE_BASE_INDEX_FILE),
        )?;
        if reader.logical_len() != expected_groups as u64 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "abundance base stream does not match archive CID count",
            ));
        }
        Some(reader)
    } else {
        None
    };
    let mut current_base = 0u8;
    let mut visit = |edge: MembershipEdge| -> Result<()> {
        if saw_edge && edge.cid != current_cid {
            if edge.cid != current_cid + 1 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "dataset index omits a CID",
                ));
            }
            writer.append_zero_based_with_abundance(
                &current_ids
                    .iter()
                    .map(|id| *id as usize - 1)
                    .collect::<Vec<_>>(),
                abundance_log_base.map(|_| current_codes.as_slice()),
            )?;
            current_ids.clear();
            current_codes.clear();
            current_cid = edge.cid;
            if let Some(reader) = abundance_bases.as_mut() {
                reader.read_exact_at(edge.cid, std::slice::from_mut(&mut current_base))?;
            }
        } else if !saw_edge {
            if edge.cid != 0 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "dataset index does not begin at CID zero",
                ));
            }
            current_cid = edge.cid;
            saw_edge = true;
            if let Some(reader) = abundance_bases.as_mut() {
                reader.read_exact_at(edge.cid, std::slice::from_mut(&mut current_base))?;
            }
        }
        if current_ids.last().is_some_and(|last| *last >= edge.dataset) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "duplicate or unsorted dataset membership",
            ));
        }
        current_ids.push(edge.dataset);
        if abundance_log_base.is_some() {
            current_codes.push(current_base.wrapping_add(edge.abundance));
        }
        Ok(())
    };
    let sink = run_directory.join("unused");
    merge_membership_runs(&runs, &sink, Some(&mut visit))?;
    drop(visit);
    if saw_edge {
        writer.append_zero_based_with_abundance(
            &current_ids
                .iter()
                .map(|id| *id as usize - 1)
                .collect::<Vec<_>>(),
            abundance_log_base.map(|_| current_codes.as_slice()),
        )?;
    }
    writer.finish()?;
    for run in runs {
        fs::remove_file(run)?;
    }
    let _ = fs::remove_dir(&run_directory);
    let sidecar = Arc::new(CidDatasetSidecar::open(output_path)?);
    if sidecar.len() != expected_groups {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "temporary CID transpose has wrong group count",
        ));
    }
    fs::remove_file(output_path)?;
    Ok(sidecar)
}

fn to_io_err(context: &str, err: impl std::fmt::Display) -> io::Error {
    io::Error::other(format!("{context}: {err}"))
}

#[derive(Clone, Copy, Debug)]
struct PhaseTiming {
    wall_sec: f64,
    cpu_sec: Option<f64>,
}

#[derive(Debug)]
struct PhaseTimer {
    wall_start: Instant,
    cpu_start: Option<f64>,
}

impl PhaseTimer {
    fn start() -> Self {
        Self {
            wall_start: Instant::now(),
            cpu_start: read_process_cpu_seconds(),
        }
    }

    fn finish(&self) -> PhaseTiming {
        let wall_sec = self.wall_start.elapsed().as_secs_f64();
        let cpu_sec = self
            .cpu_start
            .and_then(|start| read_process_cpu_seconds().map(|now| (now - start).max(0.0)));
        PhaseTiming { wall_sec, cpu_sec }
    }
}

#[derive(Debug, Default)]
struct TimingAccumulator {
    wall_sec: f64,
    cpu_sec: Option<f64>,
}

impl TimingAccumulator {
    fn add(&mut self, timing: PhaseTiming) {
        self.wall_sec += timing.wall_sec;
        self.cpu_sec = match (self.cpu_sec, timing.cpu_sec) {
            (Some(acc), Some(curr)) => Some(acc + curr),
            _ => None,
        };
    }
}

fn read_process_cpu_seconds() -> Option<f64> {
    let stat = fs::read_to_string("/proc/self/stat").ok()?;
    let rparen = stat.rfind(')')?;
    let rest = stat.get((rparen + 2)..)?;
    let fields: Vec<&str> = rest.split_whitespace().collect();
    if fields.len() <= 12 {
        return None;
    }

    let utime_ticks = fields[11].parse::<u64>().ok()? as f64;
    let stime_ticks = fields[12].parse::<u64>().ok()? as f64;
    let clk_tck = unsafe { libc::sysconf(libc::_SC_CLK_TCK) };
    if clk_tck <= 0 {
        return None;
    }

    Some((utime_ticks + stime_ticks) / (clk_tck as f64))
}

fn log_phase_timing(phase: &str, timing: PhaseTiming) {
    if let Some(cpu_sec) = timing.cpu_sec {
        let avg_cores = if timing.wall_sec > 0.0 {
            cpu_sec / timing.wall_sec
        } else {
            0.0
        };
        println!(
            "[phase-timing] phase={} wall_s={:.3} cpu_s={:.3} avg_cores={:.3}",
            phase, timing.wall_sec, cpu_sec, avg_cores
        );
    } else {
        println!(
            "[phase-timing] phase={} wall_s={:.3} cpu_s=NA avg_cores=NA",
            phase, timing.wall_sec
        );
    }
}

fn write_varint_u64(mut value: u64, out: &mut Vec<u8>) {
    while value >= 0x80 {
        out.push((value as u8 & 0x7f) | 0x80);
        value >>= 7;
    }
    out.push(value as u8);
}

fn write_varint_u64_to_writer(mut value: u64, mut out: impl Write) -> io::Result<()> {
    let mut buf = [0u8; 10];
    let mut len = 0usize;
    while value >= 0x80 {
        buf[len] = (value as u8 & 0x7f) | 0x80;
        value >>= 7;
        len += 1;
    }
    buf[len] = value as u8;
    len += 1;
    out.write_all(&buf[..len])
}

#[derive(Debug, Default)]
struct ChunkFlushMetrics {
    calls: usize,
    records: usize,
    timing: TimingAccumulator,
}

#[derive(Debug, Default)]
struct WorkerStats {
    groups: usize,
    timing: TimingAccumulator,
    max_wall_sec: f64,
    max_cpu_sec: Option<f64>,
}

impl WorkerStats {
    fn record_batch(&mut self, groups: usize, timing: PhaseTiming) {
        self.groups += groups;
        self.timing.add(timing);
        let wall_sec = timing.wall_sec;
        let cpu_sec = timing.cpu_sec;
        if wall_sec > self.max_wall_sec {
            self.max_wall_sec = wall_sec;
        }
        self.max_cpu_sec = match (self.max_cpu_sec, cpu_sec) {
            (Some(acc), Some(curr)) => Some(acc.max(curr)),
            (None, Some(curr)) => Some(curr),
            _ => self.max_cpu_sec,
        };
    }
}

#[derive(Debug, Default)]
struct CommitStats {
    groups: usize,
    tigs_bytes: u64,
    sizes_bytes: u64,
    timing: TimingAccumulator,
}

fn ggcat_extra_elaboration(
    use_unitigs: bool,
    use_matchtigs: bool,
    use_eulertigs: bool,
) -> ExtraElaboration {
    if use_unitigs {
        ExtraElaboration::None
    } else if use_matchtigs {
        ExtraElaboration::GreedyMatchtigs
    } else if use_eulertigs {
        ExtraElaboration::FastEulertigs
    } else {
        // For colored-record export we need color-safe output to preserve per-file k-mers.
        ExtraElaboration::None
    }
}

fn write_sorted_color_record(
    writer: &mut BufWriter<File>,
    record: &SortedColorRecord,
) -> Result<()> {
    writer.write_all(&record.subset.to_le_bytes())?;
    writer.write_all(&(record.seq.len() as u32).to_le_bytes())?;
    writer.write_all(&record.seq)?;
    Ok(())
}

fn read_sorted_color_record(reader: &mut BufReader<File>) -> Result<Option<SortedColorRecord>> {
    let mut len_buf = [0u8; 4];
    match reader.read_exact(&mut len_buf) {
        Ok(()) => {}
        Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(err) => return Err(err),
    }
    let subset = u32::from_le_bytes(len_buf);

    reader.read_exact(&mut len_buf)?;
    let seq_len = u32::from_le_bytes(len_buf) as usize;

    let mut seq = vec![0u8; seq_len];
    reader.read_exact(&mut seq)?;

    Ok(Some(SortedColorRecord { subset, seq }))
}

fn read_color_record(
    reader: &mut BufReader<File>,
) -> Result<Option<(Vec<ColorIndexType>, Vec<u8>)>> {
    let mut len_buf = [0u8; 4];
    match reader.read_exact(&mut len_buf) {
        Ok(()) => {}
        Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(err) => return Err(err),
    }
    let subset_count = u32::from_le_bytes(len_buf) as usize;

    reader.read_exact(&mut len_buf)?;
    let seq_len = u32::from_le_bytes(len_buf) as usize;

    let mut subsets = vec![0u32; subset_count];
    for subset in &mut subsets {
        reader.read_exact(&mut len_buf)?;
        *subset = u32::from_le_bytes(len_buf);
    }

    let mut seq = vec![0u8; seq_len];
    reader.read_exact(&mut seq)?;

    Ok(Some((subsets, seq)))
}

fn flush_sorted_chunk(
    records: &mut Vec<SortedColorRecord>,
    chunk_files: &mut Vec<PathBuf>,
    chunk_dir: &Path,
    sort_pool: &rayon::ThreadPool,
) -> Result<()> {
    if records.is_empty() {
        return Ok(());
    }

    if records.len() >= PAR_SORT_THRESHOLD {
        // GGCAT can have every worker in the global Rayon pool blocked while its
        // bounded structured-output channel is full.  Sorting on that same pool
        // would then deadlock the consumer which is needed to drain the channel.
        sort_pool.install(|| {
            records.par_sort_unstable_by(|left, right| left.subset.cmp(&right.subset));
        });
    } else {
        records.sort_unstable_by(|left, right| left.subset.cmp(&right.subset));
    }

    let chunk_path = chunk_dir.join(format!("chunk_{:08}.bin", chunk_files.len()));
    let mut chunk_writer = BufWriter::with_capacity(IO_BUFFER_CAPACITY, File::create(&chunk_path)?);
    for record in records.iter() {
        write_sorted_color_record(&mut chunk_writer, record)?;
    }
    chunk_writer.flush()?;
    chunk_files.push(chunk_path);
    records.clear();
    Ok(())
}

fn flush_sorted_chunk_timed(
    records: &mut Vec<SortedColorRecord>,
    chunk_files: &mut Vec<PathBuf>,
    chunk_dir: &Path,
    sort_pool: &rayon::ThreadPool,
    metrics: &mut ChunkFlushMetrics,
) -> Result<()> {
    if records.is_empty() {
        return Ok(());
    }
    let records_count = records.len();
    let timer = PhaseTimer::start();
    flush_sorted_chunk(records, chunk_files, chunk_dir, sort_pool)?;
    metrics.calls += 1;
    metrics.records += records_count;
    metrics.timing.add(timer.finish());
    Ok(())
}

fn split_color_run_segments(
    runs: Vec<(ColorIndexType, usize)>,
    seq: &[u8],
    k: usize,
    source: &str,
) -> Result<Vec<SortedColorRecord>> {
    if runs.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("ggcat emitted colored sequence without color runs: {source}"),
        ));
    }

    if seq.len() < k {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "ggcat sequence shorter than k (len={}, k={k}) for {source}",
                seq.len()
            ),
        ));
    }

    let seq_kmers = seq.len() - k + 1;
    let mut total_run_kmers = 0usize;
    let mut merged_runs: Vec<(ColorIndexType, usize)> = Vec::with_capacity(runs.len());

    for (subset, count) in runs {
        total_run_kmers = total_run_kmers.checked_add(count).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "overflow while summing GGCAT color-run lengths",
            )
        })?;

        if let Some((last_subset, last_count)) = merged_runs.last_mut() {
            if *last_subset == subset {
                *last_count += count;
                continue;
            }
        }
        merged_runs.push((subset, count));
    }

    if total_run_kmers != seq_kmers {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "GGCAT color-run lengths mismatch sequence k-mers for {source}: runs={}, expected={}",
                total_run_kmers,
                seq_kmers
            ),
        ));
    }

    let mut segments = Vec::with_capacity(merged_runs.len());
    let mut run_start = 0usize;
    for (subset, count) in merged_runs {
        let run_end = run_start + count;
        if run_end > seq_kmers {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "GGCAT color run exceeds sequence k-mers for {source}: end={}, total={}",
                    run_end, seq_kmers
                ),
            ));
        }
        segments.push(SortedColorRecord {
            subset,
            seq: seq[run_start..(run_end + k - 1)].to_vec(),
        });
        run_start = run_end;
    }

    if run_start != seq_kmers {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "GGCAT color runs do not fully cover sequence k-mers for {source}: covered={}, total={}",
                run_start,
                seq_kmers
            ),
        ));
    }

    Ok(segments)
}

fn read_u64_field(block: &[u8], cursor: &mut usize, field: &str) -> Result<u64> {
    let end = cursor
        .checked_add(8)
        .ok_or_else(|| io::Error::other("structured-output cursor overflow"))?;
    let bytes = block.get(*cursor..end).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("truncated GGCAT structured-output {field}"),
        )
    })?;
    *cursor = end;
    Ok(u64::from_le_bytes(bytes.try_into().unwrap()))
}

fn read_varint_field(data: &[u8], cursor: &mut usize, field: &str) -> Result<u64> {
    let mut value = 0u64;
    for shift in (0..=63).step_by(7) {
        let byte = *data.get(*cursor).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("truncated GGCAT structured-output {field}"),
            )
        })?;
        *cursor += 1;
        if shift == 63 && byte > 1 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("overflow in GGCAT structured-output {field}"),
            ));
        }
        value |= ((byte & 0x7f) as u64) << shift;
        if byte & 0x80 == 0 {
            return Ok(value);
        }
    }
    Err(io::Error::new(
        io::ErrorKind::InvalidData,
        format!("unterminated GGCAT structured-output {field}"),
    ))
}

fn visit_structured_ggcat_block(
    block: &[u8],
    k: usize,
    mut emit: impl FnMut(SortedColorRecord) -> Result<()>,
) -> Result<(usize, usize)> {
    let mut cursor = 0usize;
    let mut sequences = 0usize;
    let mut segments = 0usize;
    while cursor < block.len() {
        let sequence_index = read_u64_field(block, &mut cursor, "sequence index")?;
        let sequence_len = usize::try_from(read_u64_field(block, &mut cursor, "sequence length")?)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "sequence length overflow"))?;
        let color_len = usize::try_from(read_u64_field(block, &mut cursor, "color length")?)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "color length overflow"))?;
        let links_len = usize::try_from(read_u64_field(block, &mut cursor, "links length")?)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "links length overflow"))?;
        let record_len = sequence_len
            .checked_add(color_len)
            .and_then(|len| len.checked_add(links_len))
            .ok_or_else(|| io::Error::other("structured-output record length overflow"))?;
        let record_end = cursor
            .checked_add(record_len)
            .ok_or_else(|| io::Error::other("structured-output cursor overflow"))?;
        if record_end > block.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("truncated GGCAT structured-output sequence {sequence_index}"),
            ));
        }
        let seq = &block[cursor..cursor + sequence_len];
        cursor += sequence_len;
        let color_data = &block[cursor..cursor + color_len];
        cursor += color_len;
        cursor += links_len;

        let mut color_cursor = 0usize;
        let colors_count = usize::try_from(read_varint_field(
            color_data,
            &mut color_cursor,
            "color count",
        )?)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "color count overflow"))?;
        let mut runs = Vec::with_capacity(colors_count);
        for _ in 0..colors_count {
            let subset = ColorIndexType::try_from(read_varint_field(
                color_data,
                &mut color_cursor,
                "color subset",
            )?)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "color subset overflow"))?;
            let count = usize::try_from(read_varint_field(
                color_data,
                &mut color_cursor,
                "color run length",
            )?)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "color run overflow"))?;
            if count > 0 {
                runs.push((subset, count));
            }
        }
        if color_cursor != color_data.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "unused bytes in GGCAT structured-output colors for sequence {sequence_index}"
                ),
            ));
        }
        let parsed = split_color_run_segments(
            runs,
            seq,
            k,
            &format!("structured sequence {sequence_index}"),
        )?;
        sequences += 1;
        segments += parsed.len();
        for record in parsed {
            emit(record)?;
        }
    }
    Ok((sequences, segments))
}

fn merge_sorted_dataset_ids_with_offset(
    target: &mut Vec<u32>,
    other: &[u32],
    offset: u32,
) -> io::Result<()> {
    if offset != 0 && other.last().copied().unwrap_or(0) > u32::MAX - offset {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "merged dataset id overflow",
        ));
    }
    if target.is_empty() {
        target.reserve(other.len());
        target.extend(other.iter().map(|&id| id + offset));
        return Ok(());
    }
    if other.is_empty() {
        return Ok(());
    }
    let mut merged = Vec::with_capacity(target.len().saturating_add(other.len()));
    let (mut left, mut right) = (0usize, 0usize);
    while left < target.len() && right < other.len() {
        let right_id = other[right] + offset;
        match target[left].cmp(&right_id) {
            CmpOrdering::Less => {
                merged.push(target[left]);
                left += 1;
            }
            CmpOrdering::Greater => {
                merged.push(right_id);
                right += 1;
            }
            CmpOrdering::Equal => {
                merged.push(target[left]);
                left += 1;
                right += 1;
            }
        }
    }
    merged.extend_from_slice(&target[left..]);
    merged.extend(other[right..].iter().map(|&id| id + offset));
    *target = merged;
    Ok(())
}

fn merge_sorted_dataset_ids(target: &mut Vec<u32>, other: &[u32]) {
    merge_sorted_dataset_ids_with_offset(target, other, 0)
        .expect("zero-offset dataset IDs cannot overflow");
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct SourceDatasetSpan {
    storage: usize,
    start: usize,
    end: usize,
    offset: u32,
}

#[derive(Debug)]
pub(crate) struct SourceDatasetMap {
    storages: Vec<Arc<Vec<u32>>>,
    spans: Vec<SourceDatasetSpan>,
    dense_ranges: Vec<DenseSourceDatasetRange>,
    disk_ranges: Vec<DiskSourceDatasetRange>,
    source_count: usize,
    virtual_abundance_datasets: Option<usize>,
    abundance_log_base: Option<f64>,
}

#[derive(Debug)]
struct DenseSourceDatasetRange {
    source_start: usize,
    source_end: usize,
    storage: usize,
    offsets: Arc<Vec<usize>>,
    offset: u32,
}

#[derive(Debug)]
struct DiskSourceDatasetRange {
    source_start: usize,
    source_end: usize,
    sidecar: Arc<CidDatasetSidecar>,
    offset: u32,
}

struct SourceDatasetReader<'a> {
    map: &'a SourceDatasetMap,
    disk_readers: Vec<CidDatasetSidecarReader<'a>>,
}

impl SourceDatasetMap {
    pub(crate) fn new() -> Self {
        Self {
            storages: Vec::new(),
            spans: Vec::new(),
            dense_ranges: Vec::new(),
            disk_ranges: Vec::new(),
            source_count: 0,
            virtual_abundance_datasets: None,
            abundance_log_base: None,
        }
    }

    fn from_abundance_datasets(dataset_count: usize, log_base: f64) -> io::Result<Self> {
        if !log_base.is_finite() || log_base <= 1.0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "abundance log base must be finite and > 1",
            ));
        }
        let source_count = dataset_count
            .checked_mul(ABUNDANCE_CODES_PER_DATASET)
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidInput, "too many input datasets")
            })?;
        if source_count > u32::MAX as usize {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "dataset count is too large for 8-bit abundance colors",
            ));
        }
        Ok(Self {
            storages: Vec::new(),
            spans: Vec::new(),
            dense_ranges: Vec::new(),
            disk_ranges: Vec::new(),
            source_count,
            virtual_abundance_datasets: Some(dataset_count),
            abundance_log_base: Some(log_base),
        })
    }

    pub(crate) fn add_storage(&mut self, values: Arc<Vec<u32>>) -> usize {
        let index = self.storages.len();
        self.storages.push(values);
        index
    }

    pub(crate) fn push_span(
        &mut self,
        storage: usize,
        start: usize,
        end: usize,
        offset: u32,
    ) -> io::Result<()> {
        let values = self.storages.get(storage).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "source dataset storage is missing",
            )
        })?;
        if start > end || end > values.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "source dataset span is outside its storage",
            ));
        }
        self.spans.push(SourceDatasetSpan {
            storage,
            start,
            end,
            offset,
        });
        self.source_count = self
            .source_count
            .checked_add(1)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "too many GGCAT sources"))?;
        Ok(())
    }

    pub(crate) fn add_dense_storage(
        &mut self,
        values: Arc<Vec<u32>>,
        offsets: Arc<Vec<usize>>,
        offset: u32,
    ) -> io::Result<()> {
        if offsets.is_empty() || offsets[0] != 0 || offsets.last().copied() != Some(values.len()) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "dense source offsets do not cover their dataset storage",
            ));
        }
        if offsets.windows(2).any(|range| range[0] >= range[1]) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "archive contains a source group without any dataset membership",
            ));
        }
        let count = offsets.len() - 1;
        let source_start = self.source_count;
        let source_end = source_start
            .checked_add(count)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "too many GGCAT sources"))?;
        let storage = self.add_storage(values);
        self.dense_ranges.push(DenseSourceDatasetRange {
            source_start,
            source_end,
            storage,
            offsets,
            offset,
        });
        self.source_count = source_end;
        Ok(())
    }

    pub(crate) fn add_disk_storage(
        &mut self,
        sidecar: Arc<CidDatasetSidecar>,
        offset: u32,
    ) -> io::Result<()> {
        if let Some(base) = sidecar.abundance_log_base() {
            match self.abundance_log_base {
                Some(existing) if existing.to_bits() != base.to_bits() => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "cannot combine archives with different abundance log bases",
                    ))
                }
                None if self.source_count != 0 => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "cannot combine abundance and non-abundance sources",
                    ))
                }
                None => self.abundance_log_base = Some(base),
                Some(_) => {}
            }
        } else if self.abundance_log_base.is_some() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "cannot combine abundance and non-abundance sources",
            ));
        }
        let source_start = self.source_count;
        let source_end = source_start
            .checked_add(sidecar.len())
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "too many GGCAT sources"))?;
        self.disk_ranges.push(DiskSourceDatasetRange {
            source_start,
            source_end,
            sidecar,
            offset,
        });
        self.source_count = source_end;
        Ok(())
    }

    fn from_singletons(dataset_count: usize) -> io::Result<Self> {
        let mut values = Vec::with_capacity(dataset_count);
        for dataset in 1..=dataset_count {
            values.push(u32::try_from(dataset).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidInput, "too many input datasets")
            })?);
        }
        let mut map = Self::new();
        let storage = map.add_storage(Arc::new(values));
        for dataset in 0..dataset_count {
            map.push_span(storage, dataset, dataset + 1, 0)?;
        }
        Ok(map)
    }

    pub(crate) fn len(&self) -> usize {
        self.source_count
    }

    pub(crate) fn abundance_log_base(&self) -> Option<f64> {
        self.abundance_log_base
    }

    fn reader(&self) -> io::Result<SourceDatasetReader<'_>> {
        let mut disk_readers = Vec::with_capacity(self.disk_ranges.len());
        for range in &self.disk_ranges {
            disk_readers.push(range.sidecar.reader()?);
        }
        Ok(SourceDatasetReader {
            map: self,
            disk_readers,
        })
    }

    fn merge_source_into(&self, source: usize, target: &mut Vec<u32>) -> io::Result<()> {
        if let Some(span) = self.spans.get(source) {
            let values = self.storages.get(span.storage).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "source dataset storage is missing",
                )
            })?;
            let values = values.get(span.start..span.end).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "source dataset span is invalid")
            })?;
            return merge_sorted_dataset_ids_with_offset(target, values, span.offset);
        }
        if let Some(range) = self
            .disk_ranges
            .iter()
            .find(|range| source >= range.source_start && source < range.source_end)
        {
            return range.sidecar.merge_group_into(
                source - range.source_start,
                target,
                range.offset,
            );
        }
        let range = self
            .dense_ranges
            .iter()
            .find(|range| source >= range.source_start && source < range.source_end)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("ggcat source index {} is out of range", source),
                )
            })?;
        let local = source - range.source_start;
        let start = range.offsets[local];
        let end = range.offsets[local + 1];
        let values = self.storages.get(range.storage).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "source dataset storage is missing",
            )
        })?;
        let values = values.get(start..end).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "source dataset span is invalid")
        })?;
        merge_sorted_dataset_ids_with_offset(target, values, range.offset)
    }
}

impl SourceDatasetReader<'_> {
    fn load_source_with_abundance_into(
        &mut self,
        source: usize,
        target: &mut Vec<u32>,
        abundance_codes: &mut Vec<u8>,
    ) -> io::Result<()> {
        abundance_codes.clear();
        if let Some(dataset_count) = self.map.virtual_abundance_datasets {
            if source >= self.map.source_count {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("ggcat source index {} is out of range", source),
                ));
            }
            let dataset = source / ABUNDANCE_CODES_PER_DATASET;
            if dataset >= dataset_count {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "virtual abundance color maps outside the dataset list",
                ));
            }
            target.clear();
            target.push(u32::try_from(dataset + 1).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "dataset ID exceeds u32")
            })?);
            abundance_codes.push((source % ABUNDANCE_CODES_PER_DATASET) as u8);
            return Ok(());
        }
        if let Some(span) = self.map.spans.get(source) {
            let values = self.map.storages.get(span.storage).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "source dataset storage is missing",
                )
            })?;
            let values = values.get(span.start..span.end).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "source dataset span is invalid")
            })?;
            target.clear();
            target.reserve(values.len());
            for &value in values {
                target.push(value.checked_add(span.offset).ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "dataset ID overflow")
                })?);
            }
            return Ok(());
        }
        if let Some((index, range)) = self
            .map
            .disk_ranges
            .iter()
            .enumerate()
            .find(|(_, range)| source >= range.source_start && source < range.source_end)
        {
            return self.disk_readers[index].load_group_with_abundance_into(
                source - range.source_start,
                target,
                abundance_codes,
                range.offset,
            );
        }
        let range = self
            .map
            .dense_ranges
            .iter()
            .find(|range| source >= range.source_start && source < range.source_end)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("ggcat source index {} is out of range", source),
                )
            })?;
        let local = source - range.source_start;
        let start = range.offsets[local];
        let end = range.offsets[local + 1];
        let values = self.map.storages.get(range.storage).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "source dataset storage is missing",
            )
        })?;
        let values = values.get(start..end).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "source dataset span is invalid")
        })?;
        target.clear();
        target.reserve(values.len());
        for &value in values {
            target.push(value.checked_add(range.offset).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "dataset ID overflow")
            })?);
        }
        Ok(())
    }
}

struct ResolvedSubsets {
    subsets: Vec<ColorIndexType>,
    dataset_ids: Vec<Arc<Vec<u32>>>,
    abundance_codes: Option<Vec<Arc<Vec<u8>>>>,
}

type AbundanceMemberships = (Vec<Arc<Vec<u32>>>, Vec<Arc<Vec<u8>>>);

fn decode_virtual_abundance_sources(
    mut sources: Vec<usize>,
    source_limit: usize,
) -> Result<(Arc<Vec<u32>>, Arc<Vec<u8>>)> {
    sources.sort_unstable();
    sources.dedup();
    if sources.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "GGCAT color subset resolved to empty dataset set",
        ));
    }
    let mut dataset_ids = Vec::with_capacity(sources.len());
    let mut abundance_codes = Vec::with_capacity(sources.len());
    let mut previous_dataset = None;
    for source in sources {
        if source >= source_limit {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("GGCAT abundance color index {source} is out of range"),
            ));
        }
        let dataset = source / ABUNDANCE_CODES_PER_DATASET;
        if previous_dataset == Some(dataset) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "one dataset assigns conflicting abundance codes to the same k-mer",
            ));
        }
        dataset_ids.push(
            u32::try_from(dataset + 1).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "dataset ID exceeds u32")
            })?,
        );
        abundance_codes.push((source % ABUNDANCE_CODES_PER_DATASET) as u8);
        previous_dataset = Some(dataset);
    }
    Ok((Arc::new(dataset_ids), Arc::new(abundance_codes)))
}

fn resolve_virtual_abundance_subsets(
    subset_sources: Vec<(ColorIndexType, Vec<usize>)>,
    dataset_count: usize,
) -> Result<AbundanceMemberships> {
    let source_limit = dataset_count
        .checked_mul(ABUNDANCE_CODES_PER_DATASET)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "too many input datasets"))?;
    let memberships = subset_sources
        .into_par_iter()
        .map(|(_, sources)| decode_virtual_abundance_sources(sources, source_limit))
        .collect::<Result<Vec<_>>>()?;
    let mut dataset_ids = Vec::with_capacity(memberships.len());
    let mut abundance_codes = Vec::with_capacity(memberships.len());
    for (ids, codes) in memberships {
        dataset_ids.push(ids);
        abundance_codes.push(codes);
    }
    Ok((dataset_ids, abundance_codes))
}

fn resolve_subsets_to_dataset_ids(
    instance: &GGCATInstance,
    colormap_file: &Path,
    subsets: &[ColorIndexType],
    color_index_to_source_index: Option<&[usize]>,
    source_dataset_ids: &SourceDatasetMap,
    dataset_count: usize,
) -> Result<ResolvedSubsets> {
    if subsets.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "no color subsets were parsed from GGCAT output",
        ));
    }

    let mut subset_list = subsets.to_vec();
    subset_list.sort_unstable();
    subset_list.dedup();

    let callback_err = Mutex::new(None::<io::Error>);

    for chunk in subset_list.chunks(subset_list.len().max(1)) {
        let subset_sources = Mutex::new(Vec::<(ColorIndexType, Vec<usize>)>::with_capacity(
            chunk.len(),
        ));
        instance
            .query_colormap(
                colormap_file.to_path_buf(),
                chunk.to_vec(),
                false,
                |subset, colors| {
                    if callback_err
                        .lock()
                        .expect("callback error lock poisoned")
                        .is_some()
                    {
                        return;
                    }

                    let mut source_indices = Vec::with_capacity(colors.len());
                    for &color in colors {
                        let color_idx = color as usize;
                        let source_idx = match color_index_to_source_index {
                            Some(mapping) => match mapping.get(color_idx) {
                                Some(&source_idx) => source_idx,
                                None => {
                                    *callback_err.lock().expect("callback error lock poisoned") =
                                        Some(io::Error::new(
                                            io::ErrorKind::InvalidData,
                                            format!(
                                                "ggcat color index {} is missing from its source mapping",
                                                color_idx
                                            ),
                                        ));
                                    return;
                                }
                            },
                            None => color_idx,
                        };
                        if source_idx >= source_dataset_ids.len() {
                            *callback_err.lock().expect("callback error lock poisoned") =
                                Some(io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    format!("ggcat color index {} out of range", color_idx),
                                ));
                            return;
                        }
                        source_indices.push(source_idx);
                    }
                    if source_indices.is_empty() {
                        *callback_err.lock().expect("callback error lock poisoned") =
                            Some(io::Error::new(
                                io::ErrorKind::InvalidData,
                                "GGCAT color subset resolved to empty dataset set",
                            ));
                        return;
                    }
                    subset_sources
                        .lock()
                        .expect("subset-source lock poisoned")
                        .push((subset, source_indices));
                },
            )
            .map_err(|e| to_io_err("query ggcat colormap", e))?;

        if let Some(err) = callback_err
            .lock()
            .expect("callback error lock poisoned")
            .take()
        {
            return Err(err);
        }

        let subset_sources = subset_sources
            .into_inner()
            .expect("subset-source lock poisoned");
        if subset_sources.len() != chunk.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "queried {} color subsets but GGCAT returned {}",
                    chunk.len(),
                    subset_sources.len()
                ),
            ));
        }
        let mut subset_sources = subset_sources;
        subset_sources.sort_unstable_by_key(|(subset, _)| *subset);
        if subset_sources
            .iter()
            .map(|(subset, _)| *subset)
            .ne(chunk.iter().copied())
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "GGCAT returned unexpected color subset IDs",
            ));
        }
        if source_dataset_ids.virtual_abundance_datasets.is_some() {
            let (dataset_ids, abundance_codes) =
                resolve_virtual_abundance_subsets(subset_sources, dataset_count)?;
            return Ok(ResolvedSubsets {
                subsets: chunk.to_vec(),
                dataset_ids,
                abundance_codes: Some(abundance_codes),
            });
        }
        let request_count: usize = subset_sources
            .iter()
            .map(|(_, sources)| sources.len())
            .sum();
        let mut requests = Vec::<u64>::with_capacity(request_count);
        for (subset_index, (_, sources)) in subset_sources.into_iter().enumerate() {
            let subset_index = u32::try_from(subset_index).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "too many subsets in one window")
            })?;
            for source in sources {
                let source = u32::try_from(source).map_err(|_| {
                    io::Error::new(io::ErrorKind::InvalidData, "source index exceeds u32")
                })?;
                requests.push(((source as u64) << 32) | subset_index as u64);
            }
        }
        requests.par_sort_unstable();
        requests.dedup();
        let abundance_enabled = source_dataset_ids.abundance_log_base().is_some();
        let partial: Vec<FastMutex<Vec<(u32, u8)>>> = (0..chunk.len())
            .map(|_| FastMutex::new(Vec::new()))
            .collect();
        let target_tasks = rayon::current_num_threads().max(1).saturating_mul(4);
        let request_chunk_size = requests.len().div_ceil(target_tasks).max(1);
        requests
            .par_chunks(request_chunk_size)
            .try_for_each(|request_chunk| -> Result<()> {
                let mut reader = source_dataset_ids.reader()?;
                let mut current_source = None;
                let mut current_dataset_ids = Vec::new();
                let mut current_abundance_codes = Vec::new();
                for &request in request_chunk {
                    let source = (request >> 32) as usize;
                    let subset_index = request as u32 as usize;
                    if current_source != Some(source) {
                        reader.load_source_with_abundance_into(
                            source,
                            &mut current_dataset_ids,
                            &mut current_abundance_codes,
                        )?;
                        if abundance_enabled
                            && current_abundance_codes.len() != current_dataset_ids.len()
                        {
                            return Err(io::Error::new(
                                io::ErrorKind::InvalidData,
                                "abundance source does not have one code per dataset",
                            ));
                        }
                        current_source = Some(source);
                    }
                    let mut memberships = partial[subset_index].lock();
                    if abundance_enabled {
                        memberships.extend(
                            current_dataset_ids
                                .iter()
                                .copied()
                                .zip(current_abundance_codes.iter().copied()),
                        );
                    } else {
                        memberships.extend(
                            current_dataset_ids
                                .iter()
                                .copied()
                                .map(|dataset| (dataset, 0)),
                        );
                    }
                }
                Ok(())
            })?;

        let chunk_memberships = partial
            .into_par_iter()
            .map(|memberships| -> Result<(Arc<Vec<u32>>, Arc<Vec<u8>>)> {
                let mut memberships = memberships.into_inner();
                memberships.sort_unstable();
                memberships.dedup();
                if memberships.is_empty() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "GGCAT color subset resolved to empty dataset set",
                    ));
                }
                if abundance_enabled
                    && memberships
                        .windows(2)
                        .any(|pair| pair[0].0 == pair[1].0 && pair[0].1 != pair[1].1)
                {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "one dataset assigns conflicting abundance codes to the same k-mer",
                    ));
                }
                memberships.dedup_by_key(|membership| membership.0);
                if let Some(&(invalid, _)) = memberships.iter().find(|&&(dataset_id, _)| {
                    dataset_id == 0 || dataset_id as usize > dataset_count
                }) {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "source color maps to dataset id {} outside [1, {}]",
                            invalid, dataset_count
                        ),
                    ));
                }
                let mut dataset_ids = Vec::with_capacity(memberships.len());
                let mut abundance_codes = Vec::with_capacity(memberships.len());
                for (dataset_id, code) in memberships {
                    dataset_ids.push(dataset_id);
                    abundance_codes.push(code);
                }
                Ok((Arc::new(dataset_ids), Arc::new(abundance_codes)))
            })
            .collect::<Result<Vec<_>>>()?;
        debug_assert_eq!(chunk_memberships.len(), chunk.len());
        let mut chunk_dataset_ids = Vec::with_capacity(chunk_memberships.len());
        let mut chunk_abundance_codes =
            abundance_enabled.then(|| Vec::with_capacity(chunk_memberships.len()));
        for (dataset_ids, abundance_codes) in chunk_memberships {
            chunk_dataset_ids.push(dataset_ids);
            if let Some(codes) = chunk_abundance_codes.as_mut() {
                codes.push(abundance_codes);
            }
        }
        return Ok(ResolvedSubsets {
            subsets: chunk.to_vec(),
            dataset_ids: chunk_dataset_ids,
            abundance_codes: chunk_abundance_codes,
        });
    }
    unreachable!("non-empty subset list is processed as one chunk")
}

#[allow(clippy::too_many_arguments)]
fn emit_resolved_record_window(
    records: &mut Vec<SortedColorRecord>,
    instance: &GGCATInstance,
    colormap_file: &Path,
    color_index_to_source_index: Option<&[usize]>,
    source_dataset_ids: &SourceDatasetMap,
    dataset_count: usize,
    sender: &mpsc::SyncSender<SimplitigBatch>,
    output_batch: &mut SimplitigBatch,
    output_batch_bytes: &mut usize,
    output_batch_byte_limit: usize,
    emitted_batches_count: &mut usize,
    emitted_records_count: &mut usize,
) -> Result<()> {
    if records.is_empty() {
        return Ok(());
    }
    let mut subsets = Vec::new();
    let mut previous = None;
    for record in records.iter() {
        if previous != Some(record.subset) {
            subsets.push(record.subset);
            previous = Some(record.subset);
        }
    }
    let resolved = resolve_subsets_to_dataset_ids(
        instance,
        colormap_file,
        &subsets,
        color_index_to_source_index,
        source_dataset_ids,
        dataset_count,
    )?;
    let mut resolved_index = 0usize;
    for record in records.drain(..) {
        while resolved
            .subsets
            .get(resolved_index)
            .is_some_and(|&subset| subset < record.subset)
        {
            resolved_index += 1;
        }
        if resolved.subsets.get(resolved_index).copied() != Some(record.subset) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("missing resolved dataset ids for subset {}", record.subset),
            ));
        }
        let color_ids = resolved.dataset_ids.get(resolved_index).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("missing resolved dataset ids for subset {}", record.subset),
            )
        })?;
        let abundance_codes = resolved
            .abundance_codes
            .as_ref()
            .map(|codes| Arc::clone(&codes[resolved_index]));
        *output_batch_bytes = output_batch_bytes.saturating_add(record.seq.len());
        output_batch.push(SimplitigRecord {
            color_ids: Arc::clone(color_ids),
            abundance_codes,
            seq: record.seq,
        });
        if output_batch.len() >= COLOR_RECORD_BATCH_SIZE
            || *output_batch_bytes >= output_batch_byte_limit
        {
            let out = std::mem::take(output_batch);
            *emitted_batches_count += 1;
            *emitted_records_count += out.len();
            sender.send(out).map_err(|err| {
                io::Error::new(
                    io::ErrorKind::BrokenPipe,
                    format!("record receiver dropped: {err}"),
                )
            })?;
            *output_batch_bytes = 0;
            *output_batch = Vec::with_capacity(COLOR_RECORD_BATCH_SIZE);
        }
    }
    if !output_batch.is_empty() {
        let out = std::mem::take(output_batch);
        *emitted_batches_count += 1;
        *emitted_records_count += out.len();
        sender.send(out).map_err(|err| {
            io::Error::new(
                io::ErrorKind::BrokenPipe,
                format!("record receiver dropped: {err}"),
            )
        })?;
        *output_batch_bytes = 0;
        *output_batch = Vec::with_capacity(COLOR_RECORD_BATCH_SIZE);
    }
    Ok(())
}

fn sort_and_spill_sequence_run(
    seqs: &mut Vec<Vec<u8>>,
    run_paths: &mut Vec<PathBuf>,
    spill_dir: &Path,
) -> Result<()> {
    if seqs.is_empty() {
        return Ok(());
    }

    if seqs.len() >= PAR_SORT_THRESHOLD {
        seqs.par_sort_unstable_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
    } else {
        seqs.sort_unstable_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
    }

    let mut tmp = TempBuilder::new()
        .prefix("kloe-group-run-")
        .suffix(".seq")
        .tempfile_in(spill_dir)?;
    {
        let mut writer = BufWriter::with_capacity(IO_BUFFER_CAPACITY, tmp.as_file_mut());
        for seq in seqs.iter() {
            writer.write_all(seq)?;
            writer.write_all(b"\n")?;
        }
        writer.flush()?;
    }

    let run_path = tmp
        .into_temp_path()
        .keep()
        .map_err(|err| io::Error::other(format!("persist group spill file: {err}")))?;
    run_paths.push(run_path);
    seqs.clear();
    Ok(())
}

fn emit_sorted_group_sequences(
    seqs: &mut Vec<Vec<u8>>,
    run_paths: &mut Vec<PathBuf>,
    spill_dir: &Path,
    mut on_seq: impl FnMut(&[u8]) -> Result<()>,
) -> Result<()> {
    if run_paths.is_empty() {
        if seqs.len() >= PAR_SORT_THRESHOLD {
            seqs.par_sort_unstable_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
        } else {
            seqs.sort_unstable_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
        }
        for seq in seqs.iter() {
            on_seq(seq)?;
        }
        seqs.clear();
        return Ok(());
    }

    if !seqs.is_empty() {
        sort_and_spill_sequence_run(seqs, run_paths, spill_dir)?;
    }

    let mut readers = Vec::with_capacity(run_paths.len());
    for run in run_paths.iter() {
        readers.push(SeqRunReader::new(run)?);
    }

    let mut heap = BinaryHeap::new();
    for (run_idx, reader) in readers.iter_mut().enumerate() {
        if let Some(seq) = reader.next_seq()? {
            heap.push(SeqHeapItem { run_idx, seq });
        }
    }

    while let Some(item) = heap.pop() {
        on_seq(&item.seq)?;
        if let Some(next) = readers[item.run_idx].next_seq()? {
            heap.push(SeqHeapItem {
                run_idx: item.run_idx,
                seq: next,
            });
        }
    }

    for run_path in run_paths.drain(..) {
        let _ = fs::remove_file(run_path);
    }
    Ok(())
}

struct DatasetCidSpill {
    paths: Vec<PathBuf>,
    writers: Vec<BufWriter<File>>,
    last_cids: Vec<u64>,
    directory: PathBuf,
    datasets_per_partition: usize,
    dataset_count: usize,
    cid_count: usize,
    majority_path: PathBuf,
    majority_writer: BufWriter<File>,
    majority_byte: u8,
    majority_bits: u8,
    abundance_log_base: Option<f64>,
}

struct DatasetCidSpillFiles {
    paths: Vec<PathBuf>,
    directory: PathBuf,
    datasets_per_partition: usize,
    dataset_count: usize,
    cid_count: usize,
    majority_path: PathBuf,
    abundance_log_base: Option<f64>,
}

impl DatasetCidSpill {
    fn new(
        output_dir: &str,
        dataset_count: usize,
        abundance_log_base: Option<f64>,
    ) -> Result<Self> {
        let directory = create_id_cid_spill_dir(output_dir)?;
        let partition_count = dataset_count
            .div_ceil(DATASETS_PER_CID_PARTITION)
            .clamp(1, MAX_CID_PARTITIONS);
        let datasets_per_partition = dataset_count.div_ceil(partition_count).max(1);
        let mut paths = Vec::with_capacity(partition_count);
        let mut writers = Vec::with_capacity(partition_count);
        let last_cids = vec![0u64; partition_count];
        for partition in 0..partition_count {
            let path = directory.join(format!("partition_{partition:04}.pairs.bin"));
            writers.push(BufWriter::with_capacity(
                CID_PARTITION_BUFFER_BYTES,
                File::create(&path)?,
            ));
            paths.push(path);
        }
        let majority_path = directory.join("majority_membership.bin");
        let majority_writer =
            BufWriter::with_capacity(CID_PARTITION_BUFFER_BYTES, File::create(&majority_path)?);
        println!(
            "Dataset-to-CID spill: datasets={}, partitions={}, datasets_per_partition={}",
            dataset_count, partition_count, datasets_per_partition
        );
        Ok(Self {
            paths,
            writers,
            last_cids,
            directory,
            datasets_per_partition,
            dataset_count,
            cid_count: 0,
            majority_path,
            majority_writer,
            majority_byte: 0,
            majority_bits: 0,
            abundance_log_base,
        })
    }

    fn append_group(
        &mut self,
        dataset_ids: &[usize],
        abundance_codes: Option<&[u8]>,
        cid: usize,
    ) -> Result<()> {
        match (self.abundance_log_base, abundance_codes) {
            (Some(_), Some(codes)) if codes.len() == dataset_ids.len() => {}
            (Some(_), _) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "dataset-to-CID abundance codes do not match memberships",
                ))
            }
            (None, Some(_)) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "abundance codes supplied for a non-abundance dataset-to-CID index",
                ))
            }
            (None, None) => {}
        }
        let cid = u64::try_from(cid).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "CID cannot be represented as u64",
            )
        })?;
        if cid != self.cid_count as u64 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "CID order is not contiguous while building majority membership",
            ));
        }
        if dataset_ids.len().saturating_mul(2) >= self.dataset_count {
            self.majority_byte |= 1u8 << self.majority_bits;
        }
        self.majority_bits += 1;
        if self.majority_bits == 8 {
            self.majority_writer.write_all(&[self.majority_byte])?;
            self.majority_byte = 0;
            self.majority_bits = 0;
        }
        self.cid_count += 1;
        let mut start = 0usize;
        while start < dataset_ids.len() {
            let first_dataset = dataset_ids[start];
            if first_dataset >= self.dataset_count {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "dataset id {} outside [0, {}) while spilling CID {}",
                        first_dataset, self.dataset_count, cid
                    ),
                ));
            }
            let partition =
                (first_dataset / self.datasets_per_partition).min(self.writers.len() - 1);
            let mut end = start + 1;
            while end < dataset_ids.len()
                && dataset_ids[end] / self.datasets_per_partition == partition
            {
                end += 1;
            }

            let cid_delta = cid.checked_sub(self.last_cids[partition]).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CID order regressed while writing dataset transpose spill",
                )
            })?;
            let writer = &mut self.writers[partition];
            write_varint_u64_to_writer(cid_delta, &mut *writer)?;
            write_varint_u64_to_writer((end - start) as u64, &mut *writer)?;
            let partition_start = partition.saturating_mul(self.datasets_per_partition);
            let mut previous_dataset = partition_start;
            for (local_index, &dataset) in dataset_ids[start..end].iter().enumerate() {
                if dataset >= self.dataset_count || dataset < previous_dataset {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "dataset IDs are not sorted within a color set",
                    ));
                }
                write_varint_u64_to_writer((dataset - previous_dataset) as u64, &mut *writer)?;
                if let Some(codes) = abundance_codes {
                    // The CID's first abundance is stored once in the indexed
                    // base stream. Dataset postings keep only deviations from
                    // that base, which are strongly concentrated around zero.
                    writer.write_all(&[codes[start + local_index].wrapping_sub(codes[0])])?;
                }
                previous_dataset = dataset;
            }
            self.last_cids[partition] = cid;
            start = end;
        }
        Ok(())
    }

    fn finalize(mut self) -> Result<DatasetCidSpillFiles> {
        for writer in &mut self.writers {
            writer.flush()?;
        }
        drop(self.writers);
        if self.majority_bits != 0 {
            self.majority_writer.write_all(&[self.majority_byte])?;
        }
        self.majority_writer.flush()?;
        drop(self.majority_writer);
        Ok(DatasetCidSpillFiles {
            paths: self.paths,
            directory: self.directory,
            datasets_per_partition: self.datasets_per_partition,
            dataset_count: self.dataset_count,
            cid_count: self.cid_count,
            majority_path: self.majority_path,
            abundance_log_base: self.abundance_log_base,
        })
    }
}

struct PositionSpill {
    path: PathBuf,
    writer: BufWriter<File>,
    entries: u64,
    previous_tigs: u64,
    previous_sizes: u64,
}

impl PositionSpill {
    fn new(directory: &Path) -> Result<Self> {
        let path = directory.join("positions.deltas.tmp");
        let writer = BufWriter::with_capacity(IO_BUFFER_CAPACITY, File::create(&path)?);
        let mut spill = Self {
            path,
            writer,
            entries: 0,
            previous_tigs: 0,
            previous_sizes: 0,
        };
        spill.append(0, 0)?;
        Ok(spill)
    }

    fn append(&mut self, tigs: u64, sizes: u64) -> Result<()> {
        if tigs < self.previous_tigs || sizes < self.previous_sizes {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "positions are not monotonic",
            ));
        }
        write_varint_u64_to_writer(tigs - self.previous_tigs, &mut self.writer)?;
        self.previous_tigs = tigs;
        self.previous_sizes = sizes;
        self.entries += 1;
        Ok(())
    }

    fn finalize(mut self) -> Result<(PathBuf, u64)> {
        self.writer.flush()?;
        Ok((self.path, self.entries))
    }
}

struct StreamFinalize {
    position_path: PathBuf,
    position_entries: u64,
    spill_directory: PathBuf,
    dataset_cid_spill: DatasetCidSpillFiles,
}

struct StreamWriterState {
    tigs_file: PackedTigsWriter,
    abundance_bases: Option<PackedTigsWriter>,
    size_file: BufWriter<File>,
    size_index_path: PathBuf,
    size_file_bytes: u64,
    size_groups_written: u64,
    size_blocks: Vec<(u64, u32, u64, u64)>,
    dataset_cid_spill: Option<DatasetCidSpill>,
    spill_directory: PathBuf,
    position_spill: PositionSpill,
    size_block_groups: usize,
    size_block_offsets: Vec<u32>,
    size_block_uncompressed: Vec<u8>,
    prev_tigs_size: u64,
    prev_bucket_pos: u64,
    cid: usize,
    encoded_seq_buffer: Vec<u8>,
}

impl StreamWriterState {
    fn new(
        unitigs_file_path: String,
        output_dir: &str,
        nb_files: usize,
        abundance_log_base: Option<f64>,
    ) -> Result<Self> {
        let spill_directory = create_id_cid_spill_dir(output_dir)?;
        let position_spill = PositionSpill::new(&spill_directory)?;
        let dataset_cid_spill = DatasetCidSpill::new(output_dir, nb_files, abundance_log_base)?;
        Ok(Self {
            tigs_file: PackedTigsWriter::create_indexed(
                unitigs_file_path,
                Path::new(output_dir).join(TIGS_INDEX_FILE),
            )?,
            abundance_bases: if abundance_log_base.is_some() {
                Some(PackedTigsWriter::create_indexed(
                    Path::new(output_dir).join(ABUNDANCE_BASE_FILE),
                    Path::new(output_dir).join(ABUNDANCE_BASE_INDEX_FILE),
                )?)
            } else {
                None
            },
            size_file: {
                let mut out = BufWriter::with_capacity(
                    IO_BUFFER_CAPACITY,
                    File::create(output_dir.to_owned() + "bucket_sizes.txt")?,
                );
                out.write_all(BUCKET_SIZES_COMPACT_MAGIC)?;
                out
            },
            size_index_path: Path::new(output_dir).join(BUCKET_SIZES_INDEX_FILE),
            size_file_bytes: BUCKET_SIZES_COMPACT_MAGIC.len() as u64,
            size_groups_written: 0,
            size_blocks: Vec::new(),
            dataset_cid_spill: Some(dataset_cid_spill),
            spill_directory,
            position_spill,
            size_block_groups: 0,
            size_block_offsets: vec![0],
            size_block_uncompressed: Vec::new(),
            prev_tigs_size: 0,
            prev_bucket_pos: 0,
            cid: 0,
            encoded_seq_buffer: Vec::with_capacity(ENCODED_SEQ_BUFFER_TARGET),
        })
    }

    fn flush_size_block(&mut self) -> Result<()> {
        if self.size_block_groups == 0 {
            return Ok(());
        }

        let mut length_prefixed = Vec::with_capacity(
            self.size_block_uncompressed
                .len()
                .saturating_add(self.size_block_groups.saturating_mul(2)),
        );
        for offsets in self.size_block_offsets.windows(2) {
            let start = offsets[0] as usize;
            let end = offsets[1] as usize;
            write_varint_u64((end - start) as u64, &mut length_prefixed);
            length_prefixed.extend_from_slice(&self.size_block_uncompressed[start..end]);
        }
        let mut compressed = Vec::new();
        {
            let mut encoder = Encoder::new(&mut compressed, 1)?;
            encoder.write_all(&length_prefixed)?;
            encoder.finish()?;
        }

        let group_count = u32::try_from(self.size_block_groups).map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, "too many groups in size block")
        })?;
        let compressed_len = u64::try_from(compressed.len()).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "compressed size block is too large",
            )
        })?;
        let data_offset = self
            .size_file_bytes
            .checked_add(12)
            .ok_or_else(|| io::Error::other("size file offset overflow"))?;
        self.size_blocks.push((
            self.size_groups_written,
            group_count,
            data_offset,
            compressed_len,
        ));
        self.size_file
            .write_all(&(self.size_block_groups as u32).to_le_bytes())?;
        self.size_file
            .write_all(&(compressed.len() as u64).to_le_bytes())?;
        self.size_file.write_all(&compressed)?;
        self.size_file_bytes = data_offset
            .checked_add(compressed_len)
            .ok_or_else(|| io::Error::other("size file offset overflow"))?;
        self.size_groups_written = self
            .size_groups_written
            .checked_add(u64::from(group_count))
            .ok_or_else(|| io::Error::other("size group count overflow"))?;

        self.size_block_groups = 0;
        self.size_block_offsets.clear();
        self.size_block_offsets.push(0);
        self.size_block_uncompressed.clear();
        Ok(())
    }

    fn append_bucket_group_sizes_payload(&mut self, payload: &[u8]) -> Result<()> {
        self.size_block_uncompressed.extend_from_slice(payload);
        let next_offset = u32::try_from(self.size_block_uncompressed.len()).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "bucket sizes block exceeds 4GiB uncompressed payload",
            )
        })?;
        self.size_block_offsets.push(next_offset);
        self.size_block_groups += 1;

        if self.size_block_groups >= BUCKET_SIZE_BLOCK_MAX_GROUPS
            || self.size_block_uncompressed.len() >= BUCKET_SIZE_BLOCK_MAX_UNCOMPRESSED_BYTES
        {
            self.flush_size_block()?;
        }
        Ok(())
    }

    fn append_group_cids(
        &mut self,
        dataset_ids_zero_based: &[usize],
        abundance_codes: Option<&[u8]>,
    ) -> Result<()> {
        match (&mut self.abundance_bases, abundance_codes) {
            (Some(writer), Some(codes)) if !codes.is_empty() => writer.write_all(&[codes[0]])?,
            (Some(_), _) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "abundance color set has no base code",
                ))
            }
            (None, Some(_)) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "unexpected abundance codes in non-abundance archive",
                ))
            }
            (None, None) => {}
        }
        self.dataset_cid_spill
            .as_mut()
            .expect("dataset-to-CID spill must be initialized")
            .append_group(dataset_ids_zero_based, abundance_codes, self.cid)?;
        self.cid += 1;
        Ok(())
    }

    fn write_group(
        &mut self,
        dataset_ids_zero_based: &[usize],
        abundance_codes: Option<&[u8]>,
        seqs: &mut Vec<Vec<u8>>,
        run_paths: &mut Vec<PathBuf>,
        spill_dir: &Path,
    ) -> Result<()> {
        if seqs.is_empty() && run_paths.is_empty() {
            return Ok(());
        }

        let mut group_sizes_buffer: Vec<u8> = Vec::new();

        let mut prev_size: usize = 0;
        emit_sorted_group_sequences(seqs, run_paths, spill_dir, |seq| {
            let encoded_seq = <Converter as Convert<&[u8]>>::str2num(seq);
            self.prev_tigs_size += encoded_seq.len() as u64;
            self.encoded_seq_buffer.extend_from_slice(&encoded_seq);
            if self.encoded_seq_buffer.len() >= ENCODED_SEQ_BUFFER_TARGET {
                self.tigs_file.write_all(&self.encoded_seq_buffer)?;
                self.encoded_seq_buffer.clear();
            }

            let size = seq.len();
            let delta = size.checked_sub(prev_size).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "group sizes are not nondecreasing",
                )
            })?;
            write_varint_u64(delta as u64, &mut group_sizes_buffer);
            prev_size = size;
            Ok(())
        })?;
        self.append_bucket_group_sizes_payload(&group_sizes_buffer)?;
        self.prev_bucket_pos += 1;
        self.position_spill
            .append(self.prev_tigs_size, self.prev_bucket_pos)?;

        self.append_group_cids(dataset_ids_zero_based, abundance_codes)?;

        Ok(())
    }

    fn write_precomputed_group_files(
        &mut self,
        dataset_ids_zero_based: &[usize],
        abundance_codes: Option<&[u8]>,
        encoded_tigs_path: &Path,
        encoded_tigs_len: u64,
        group_sizes_path: &Path,
        group_sizes_len: u64,
    ) -> Result<()> {
        let mut encoded_reader =
            BufReader::with_capacity(IO_BUFFER_CAPACITY, File::open(encoded_tigs_path)?);
        io::copy(&mut encoded_reader, &mut self.tigs_file)?;
        self.prev_tigs_size += encoded_tigs_len;

        self.prev_bucket_pos += 1;
        self.position_spill
            .append(self.prev_tigs_size, self.prev_bucket_pos)?;
        let mut group_sizes_payload = Vec::with_capacity(group_sizes_len as usize);
        let mut sizes_reader =
            BufReader::with_capacity(IO_BUFFER_CAPACITY, File::open(group_sizes_path)?);
        sizes_reader.read_to_end(&mut group_sizes_payload)?;
        self.append_bucket_group_sizes_payload(&group_sizes_payload)?;

        self.append_group_cids(dataset_ids_zero_based, abundance_codes)?;
        Ok(())
    }

    fn write_precomputed_group_buffers(
        &mut self,
        dataset_ids_zero_based: &[usize],
        abundance_codes: Option<&[u8]>,
        encoded_tigs: &[u8],
        group_sizes: &[u8],
    ) -> Result<()> {
        self.tigs_file.write_all(encoded_tigs)?;
        self.prev_tigs_size += encoded_tigs.len() as u64;

        self.prev_bucket_pos += 1;
        self.position_spill
            .append(self.prev_tigs_size, self.prev_bucket_pos)?;
        self.append_bucket_group_sizes_payload(group_sizes)?;

        self.append_group_cids(dataset_ids_zero_based, abundance_codes)?;
        Ok(())
    }

    fn finalize(mut self) -> Result<StreamFinalize> {
        if !self.encoded_seq_buffer.is_empty() {
            self.tigs_file.write_all(&self.encoded_seq_buffer)?;
            self.encoded_seq_buffer.clear();
        }
        self.flush_size_block()?;
        self.tigs_file.finish()?;
        if let Some(abundance_bases) = self.abundance_bases.take() {
            abundance_bases.finish()?;
        }
        self.size_file.flush()?;
        let mut size_index =
            BufWriter::with_capacity(IO_BUFFER_CAPACITY, File::create(&self.size_index_path)?);
        size_index.write_all(BUCKET_SIZES_INDEX_MAGIC)?;
        size_index.write_all(&self.size_groups_written.to_le_bytes())?;
        size_index.write_all(&(self.size_blocks.len() as u64).to_le_bytes())?;
        size_index.write_all(&(SIZE_INDEX_STRIDE as u32).to_le_bytes())?;
        size_index.write_all(&0u32.to_le_bytes())?;
        for (first_group, group_count, data_offset, data_len) in &self.size_blocks {
            size_index.write_all(&first_group.to_le_bytes())?;
            size_index.write_all(&group_count.to_le_bytes())?;
            size_index.write_all(&0u32.to_le_bytes())?;
            size_index.write_all(&data_offset.to_le_bytes())?;
            size_index.write_all(&data_len.to_le_bytes())?;
        }
        let routing_count = self.size_groups_written.div_ceil(SIZE_INDEX_STRIDE as u64);
        size_index.write_all(&routing_count.to_le_bytes())?;
        let mut block_index = 0usize;
        for route in 0..routing_count {
            let group = route * SIZE_INDEX_STRIDE as u64;
            while block_index + 1 < self.size_blocks.len()
                && self.size_blocks[block_index].0 + u64::from(self.size_blocks[block_index].1)
                    <= group
            {
                block_index += 1;
            }
            size_index.write_all(&(block_index as u32).to_le_bytes())?;
        }
        size_index.flush()?;
        let dataset_cid_spill = self
            .dataset_cid_spill
            .take()
            .expect("dataset-to-CID spill must be initialized")
            .finalize()?;

        println!(
            "Completed compression: total tigs={}, total sizes={}",
            self.prev_tigs_size, self.prev_bucket_pos
        );
        let (position_path, position_entries) = self.position_spill.finalize()?;
        Ok(StreamFinalize {
            position_path,
            position_entries,
            spill_directory: self.spill_directory,
            dataset_cid_spill,
        })
    }
}

#[derive(Debug)]
struct GroupTask {
    cid: usize,
    dataset_ids_zero_based: Vec<usize>,
    abundance_codes: Option<Vec<u8>>,
    seqs: Vec<Vec<u8>>,
    run_paths: Vec<PathBuf>,
}

#[derive(Debug)]
enum GroupWriteData {
    InMemory {
        encoded_tigs: Vec<u8>,
        group_sizes: Vec<u8>,
    },
    Spilled {
        encoded_tigs_path: PathBuf,
        encoded_tigs_len: u64,
        group_sizes_path: PathBuf,
        group_sizes_len: u64,
    },
}

#[derive(Debug)]
struct GroupWriteResult {
    cid: usize,
    dataset_ids_zero_based: Vec<usize>,
    abundance_codes: Option<Vec<u8>>,
    data: GroupWriteData,
}

fn map_color_ids_to_zero_based(color_ids: &[u32], max_datasets: usize) -> io::Result<Vec<usize>> {
    let ids = color_ids
        .iter()
        .map(|&id| {
            if id == 0 {
                Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "dataset id 0 is invalid",
                ))
            } else {
                Ok(id as usize - 1)
            }
        })
        .collect::<io::Result<Vec<usize>>>()?;

    for &id in &ids {
        if id >= max_datasets {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("dataset id {} outside [0, {})", id, max_datasets),
            ));
        }
    }
    Ok(ids)
}

fn process_single_group_task(
    mut task: GroupTask,
    spill_dir: &Path,
) -> io::Result<GroupWriteResult> {
    let seq_bytes = task.seqs.iter().map(Vec::len).sum::<usize>();
    let use_in_memory = task.run_paths.is_empty() && seq_bytes <= GROUP_IN_MEMORY_MAX_SEQ_BYTES;

    let data = if use_in_memory {
        let mut encoded_tigs = Vec::with_capacity(seq_bytes / 4 + 1024);
        let mut group_sizes = Vec::new();
        let mut prev_size: usize = 0;
        emit_sorted_group_sequences(&mut task.seqs, &mut task.run_paths, spill_dir, |seq| {
            let encoded_seq = <Converter as Convert<&[u8]>>::str2num(seq);
            encoded_tigs.extend_from_slice(&encoded_seq);
            let size = seq.len();
            let delta = size.checked_sub(prev_size).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "group sizes are not nondecreasing",
                )
            })?;
            write_varint_u64(delta as u64, &mut group_sizes);
            prev_size = size;
            Ok(())
        })?;
        GroupWriteData::InMemory {
            encoded_tigs,
            group_sizes,
        }
    } else {
        let encoded_tigs_path = spill_dir.join(format!("group_{:012}.tigs.tmp", task.cid));
        let group_sizes_path = spill_dir.join(format!("group_{:012}.sizes.tmp", task.cid));

        let encoded_file = File::create(&encoded_tigs_path)?;
        let mut encoded_writer = BufWriter::with_capacity(IO_BUFFER_CAPACITY, encoded_file);
        let mut sizes_writer =
            BufWriter::with_capacity(IO_BUFFER_CAPACITY, File::create(&group_sizes_path)?);
        let mut prev_size: usize = 0;

        let emit_res =
            emit_sorted_group_sequences(&mut task.seqs, &mut task.run_paths, spill_dir, |seq| {
                let encoded_seq = <Converter as Convert<&[u8]>>::str2num(seq);
                encoded_writer.write_all(&encoded_seq)?;

                let size = seq.len();
                let delta = size.checked_sub(prev_size).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "group sizes are not nondecreasing",
                    )
                })?;
                write_varint_u64_to_writer(delta as u64, &mut sizes_writer)?;
                prev_size = size;
                Ok(())
            });

        if emit_res.is_err() {
            for run_path in task.run_paths.drain(..) {
                let _ = fs::remove_file(run_path);
            }
            let _ = fs::remove_file(&encoded_tigs_path);
            let _ = fs::remove_file(&group_sizes_path);
        }
        emit_res?;
        sizes_writer.flush()?;
        encoded_writer.flush()?;
        drop(sizes_writer);
        drop(encoded_writer);

        let encoded_tigs_len = fs::metadata(&encoded_tigs_path)?.len();
        let group_sizes_len = fs::metadata(&group_sizes_path)?.len();
        GroupWriteData::Spilled {
            encoded_tigs_path,
            encoded_tigs_len,
            group_sizes_path,
            group_sizes_len,
        }
    };

    Ok(GroupWriteResult {
        cid: task.cid,
        dataset_ids_zero_based: task.dataset_ids_zero_based,
        abundance_codes: task.abundance_codes,
        data,
    })
}

fn commit_group_result(
    next: GroupWriteResult,
    writer: &mut StreamWriterState,
    commit_stats: &mut CommitStats,
) -> Result<()> {
    if next.cid != writer.cid {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "out-of-order group result: expected CID {}, received {}",
                writer.cid, next.cid
            ),
        ));
    }
    match &next.data {
        GroupWriteData::InMemory {
            encoded_tigs,
            group_sizes,
        } => {
            writer.write_precomputed_group_buffers(
                &next.dataset_ids_zero_based,
                next.abundance_codes.as_deref(),
                encoded_tigs,
                group_sizes,
            )?;
            commit_stats.tigs_bytes += encoded_tigs.len() as u64;
            commit_stats.sizes_bytes += group_sizes.len() as u64;
        }
        GroupWriteData::Spilled {
            encoded_tigs_path,
            encoded_tigs_len,
            group_sizes_path,
            group_sizes_len,
        } => {
            writer.write_precomputed_group_files(
                &next.dataset_ids_zero_based,
                next.abundance_codes.as_deref(),
                encoded_tigs_path,
                *encoded_tigs_len,
                group_sizes_path,
                *group_sizes_len,
            )?;
            commit_stats.tigs_bytes += *encoded_tigs_len;
            commit_stats.sizes_bytes += *group_sizes_len;
            let _ = fs::remove_file(encoded_tigs_path);
            let _ = fs::remove_file(group_sizes_path);
        }
    }
    commit_stats.groups += 1;
    Ok(())
}

fn process_and_commit_group_batch(
    tasks: &mut Vec<GroupTask>,
    resident_bytes: &mut usize,
    spill_dir: &Path,
    pool: &rayon::ThreadPool,
    writer: &mut StreamWriterState,
    worker_stats: &mut WorkerStats,
    commit_stats: &mut CommitStats,
) -> Result<()> {
    if tasks.is_empty() {
        return Ok(());
    }
    let tasks = std::mem::take(tasks);
    *resident_bytes = 0;
    let groups_count = tasks.len();
    let worker_timer = PhaseTimer::start();
    let results = pool.install(|| {
        tasks
            .into_par_iter()
            .map(|task| process_single_group_task(task, spill_dir))
            .collect::<io::Result<Vec<_>>>()
    })?;
    worker_stats.record_batch(groups_count, worker_timer.finish());
    let commit_timer = PhaseTimer::start();
    for result in results {
        commit_group_result(result, writer, commit_stats)?;
    }
    commit_stats.timing.add(commit_timer.finish());
    Ok(())
}

fn group_task_resident_bytes(task: &GroupTask) -> usize {
    task.dataset_ids_zero_based
        .capacity()
        .saturating_mul(std::mem::size_of::<usize>())
        .saturating_add(
            task.seqs
                .iter()
                .map(|seq| {
                    seq.capacity()
                        .saturating_add(std::mem::size_of::<Vec<u8>>())
                })
                .sum::<usize>(),
        )
        .saturating_add(
            task.run_paths
                .capacity()
                .saturating_mul(std::mem::size_of::<PathBuf>()),
        )
        .saturating_add(
            task.abundance_codes
                .as_ref()
                .map_or(0, |codes| codes.capacity()),
        )
}

fn write_compressed_from_stream(
    unitigs_file_path: String,
    output_dir: &String,
    nb_files: u32,
    worker_threads: usize,
    memory_budget: CompressionMemoryBudget,
    abundance_log_base: Option<f64>,
    record_rx: mpsc::Receiver<SimplitigBatch>,
) -> Result<StreamFinalize> {
    let total_timer = PhaseTimer::start();
    let mut writer = StreamWriterState::new(
        unitigs_file_path,
        output_dir,
        nb_files as usize,
        abundance_log_base,
    )?;

    let group_spill_dir = Path::new(output_dir).join(format!(
        ".kloe-group-sort-spill-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0)
    ));
    fs::create_dir_all(&group_spill_dir)?;

    let worker_count = worker_threads.max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(worker_count)
        .build()
        .map_err(|err| io::Error::other(format!("create bounded group worker pool: {err}")))?;
    let mut worker_stats = WorkerStats::default();
    let mut commit_stats = CommitStats::default();
    let mut current_color_ids: Option<Arc<Vec<u32>>> = None;
    let mut current_abundance_codes: Option<Arc<Vec<u8>>> = None;
    let mut current_ids_zero_based: Option<Vec<usize>> = None;
    let mut group_seqs: Vec<Vec<u8>> = Vec::new();
    let mut group_run_paths: Vec<PathBuf> = Vec::new();
    let mut pending_group_tasks: Vec<GroupTask> = Vec::with_capacity(GROUP_WORK_BATCH_GROUPS);
    let mut pending_group_resident_bytes: usize = 0;
    let mut group_bytes: usize = 0;
    let mut submitted_groups: usize = 0;
    let group_batch_budget = memory_budget.group_batch_bytes();
    let group_spill_threshold = GROUP_SORT_SPILL_BYTES.min((group_batch_budget / 2).max(1));
    let dispatch_timer = PhaseTimer::start();

    for batch in record_rx {
        for record in batch {
            let key_changed = match current_color_ids.as_ref() {
                None => true,
                Some(ids) => {
                    let same_ids = Arc::ptr_eq(ids, &record.color_ids)
                        || ids.as_ref() == record.color_ids.as_ref();
                    let same_abundance = match (
                        current_abundance_codes.as_ref(),
                        record.abundance_codes.as_ref(),
                    ) {
                        (None, None) => true,
                        (Some(left), Some(right)) => {
                            Arc::ptr_eq(left, right) || left.as_ref() == right.as_ref()
                        }
                        _ => false,
                    };
                    !(same_ids && same_abundance)
                }
            };

            if key_changed {
                if let Some(ids) = current_ids_zero_based.take() {
                    if !group_seqs.is_empty() || !group_run_paths.is_empty() {
                        let task = GroupTask {
                            cid: submitted_groups,
                            dataset_ids_zero_based: ids,
                            abundance_codes: current_abundance_codes
                                .take()
                                .map(|codes| codes.as_ref().clone()),
                            seqs: std::mem::take(&mut group_seqs),
                            run_paths: std::mem::take(&mut group_run_paths),
                        };
                        pending_group_resident_bytes = pending_group_resident_bytes
                            .saturating_add(group_task_resident_bytes(&task));
                        pending_group_tasks.push(task);
                        submitted_groups += 1;
                        if pending_group_tasks.len() >= GROUP_WORK_BATCH_GROUPS
                            || pending_group_resident_bytes >= group_batch_budget
                        {
                            process_and_commit_group_batch(
                                &mut pending_group_tasks,
                                &mut pending_group_resident_bytes,
                                &group_spill_dir,
                                &pool,
                                &mut writer,
                                &mut worker_stats,
                                &mut commit_stats,
                            )?;
                        }
                    }
                }

                current_ids_zero_based = Some(map_color_ids_to_zero_based(
                    &record.color_ids,
                    nb_files as usize,
                )?);
                current_color_ids = Some(Arc::clone(&record.color_ids));
                current_abundance_codes = record.abundance_codes.as_ref().map(Arc::clone);
                group_bytes = 0;
            }

            group_bytes += record.seq.len();
            group_seqs.push(record.seq);
            if group_bytes >= group_spill_threshold {
                sort_and_spill_sequence_run(
                    &mut group_seqs,
                    &mut group_run_paths,
                    &group_spill_dir,
                )?;
                group_bytes = 0;
            }
        }
    }

    if let Some(ids) = current_ids_zero_based.take() {
        if !group_seqs.is_empty() || !group_run_paths.is_empty() {
            let task = GroupTask {
                cid: submitted_groups,
                dataset_ids_zero_based: ids,
                abundance_codes: current_abundance_codes
                    .take()
                    .map(|codes| codes.as_ref().clone()),
                seqs: std::mem::take(&mut group_seqs),
                run_paths: std::mem::take(&mut group_run_paths),
            };
            pending_group_resident_bytes =
                pending_group_resident_bytes.saturating_add(group_task_resident_bytes(&task));
            pending_group_tasks.push(task);
            submitted_groups += 1;
        }
    }

    process_and_commit_group_batch(
        &mut pending_group_tasks,
        &mut pending_group_resident_bytes,
        &group_spill_dir,
        &pool,
        &mut writer,
        &mut worker_stats,
        &mut commit_stats,
    )?;
    let dispatch_timing = dispatch_timer.finish();
    let finalized = writer.finalize()?;
    let _ = fs::remove_dir_all(&group_spill_dir);

    log_phase_timing("post_ggcat_stream.dispatch_records", dispatch_timing);
    log_phase_timing(
        "post_ggcat_stream.worker_groups_total",
        PhaseTiming {
            wall_sec: worker_stats.timing.wall_sec,
            cpu_sec: worker_stats.timing.cpu_sec,
        },
    );
    log_phase_timing(
        "post_ggcat_stream.writer_commit_groups",
        PhaseTiming {
            wall_sec: commit_stats.timing.wall_sec,
            cpu_sec: commit_stats.timing.cpu_sec,
        },
    );
    log_phase_timing("post_ggcat_stream.total", total_timer.finish());
    println!(
        "[phase-stats] phase=post_ggcat_stream workers={} groups_submitted={} groups_processed={} groups_committed={} commit_tigs_bytes={} commit_sizes_bytes={} worker_max_batch_wall_s={:.3} worker_max_batch_cpu_s={}",
        worker_count,
        submitted_groups,
        worker_stats.groups,
        commit_stats.groups,
        commit_stats.tigs_bytes,
        commit_stats.sizes_bytes,
        worker_stats.max_wall_sec,
        worker_stats
            .max_cpu_sec
            .map(|v| format!("{v:.3}"))
            .unwrap_or_else(|| "NA".to_string())
    );

    Ok(finalized)
}

fn write_positions_from_spill(
    position_spill_path: &Path,
    entries: u64,
    filepath: String,
    index_filepath: String,
) -> Result<()> {
    let mut pos_file = BufWriter::with_capacity(IO_BUFFER_CAPACITY, File::create(filepath)?);
    pos_file.write_all(POSITIONS_COMPACT_MAGIC)?;
    write_varint_u64_to_writer(entries, &mut pos_file)?;
    let mut spill_reader =
        BufReader::with_capacity(IO_BUFFER_CAPACITY, File::open(position_spill_path)?);
    let mut checkpoints = Vec::<(u64, u64)>::new();
    let mut absolute = 0u64;
    let mut seen = 0u64;
    while let Some(delta) = read_optional_varint_u64(&mut spill_reader)? {
        absolute = absolute
            .checked_add(delta)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "position index overflow"))?;
        write_varint_u64_to_writer(delta, &mut pos_file)?;
        if seen as usize % POSITION_INDEX_STRIDE == 0 {
            checkpoints.push((spill_reader.stream_position()?, absolute));
        }
        seen += 1;
    }
    if seen != entries {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("position spill has {seen} entries, expected {entries}"),
        ));
    }
    pos_file.flush()?;

    let mut index = BufWriter::with_capacity(IO_BUFFER_CAPACITY, File::create(index_filepath)?);
    index.write_all(POSITIONS_INDEX_MAGIC)?;
    index.write_all(&(POSITION_INDEX_STRIDE as u32).to_le_bytes())?;
    index.write_all(&entries.to_le_bytes())?;
    index.write_all(&(checkpoints.len() as u64).to_le_bytes())?;
    for (relative_after, position) in checkpoints {
        index.write_all(&relative_after.to_le_bytes())?;
        index.write_all(&position.to_le_bytes())?;
    }
    index.flush()?;
    fs::remove_file(position_spill_path)?;
    Ok(())
}

fn create_id_cid_spill_dir(output_dir: &str) -> Result<PathBuf> {
    let pid = std::process::id();
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let path = Path::new(output_dir).join(format!(".kloe-id-cid-spill-{pid}-{nanos}"));
    fs::create_dir_all(&path)?;
    Ok(path)
}

fn read_optional_varint_u64(reader: &mut impl Read) -> Result<Option<u64>> {
    let mut value = 0u64;
    let mut shift = 0u32;
    let mut saw_byte = false;
    loop {
        let mut byte = [0u8; 1];
        match reader.read_exact(&mut byte) {
            Ok(()) => saw_byte = true,
            Err(err) if err.kind() == io::ErrorKind::UnexpectedEof && !saw_byte => return Ok(None),
            Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "truncated varint in dataset-to-CID spill record",
                ))
            }
            Err(err) => return Err(err),
        }
        value |= ((byte[0] & 0x7f) as u64) << shift;
        if byte[0] & 0x80 == 0 {
            return Ok(Some(value));
        }
        shift += 7;
        if shift >= 64 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "CID delta varint overflow in dataset-to-CID spill record",
            ));
        }
    }
}

fn read_required_varint_u64(reader: &mut impl Read) -> Result<u64> {
    read_optional_varint_u64(reader)?.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "missing varint in dataset-to-CID spill record",
        )
    })
}

#[derive(Debug, Default)]
struct DatasetCidBlocks {
    last_cid: u64,
    seen: bool,
    encoded_deltas: Vec<u8>,
    blocks: Vec<(u64, u32)>,
}

fn flush_cid_transpose_block(
    dataset: &mut DatasetCidBlocks,
    writer: &mut BufWriter<File>,
    spool_offset: &mut u64,
) -> Result<()> {
    if dataset.encoded_deltas.is_empty() {
        return Ok(());
    }
    let block_len = u32::try_from(dataset.encoded_deltas.len()).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "dataset-to-CID transpose block is larger than u32",
        )
    })?;
    dataset.blocks.push((*spool_offset, block_len));
    writer.write_all(&dataset.encoded_deltas)?;
    *spool_offset = spool_offset
        .checked_add(block_len as u64)
        .ok_or_else(|| io::Error::other("dataset-to-CID transpose offset overflow"))?;
    dataset.encoded_deltas.clear();
    Ok(())
}

fn transpose_cid_partition(
    path: &Path,
    partition_index: usize,
    spill_dir: &Path,
    datasets_per_partition: usize,
    dataset_count: usize,
    abundance_enabled: bool,
    payload_writer: &mut DatasetPayloadWriter,
) -> Result<(u64, u64)> {
    let dataset_start = partition_index.saturating_mul(datasets_per_partition);
    let dataset_end = dataset_start
        .saturating_add(datasets_per_partition)
        .min(dataset_count);
    let mut datasets = (dataset_start..dataset_end)
        .map(|_| DatasetCidBlocks::default())
        .collect::<Vec<_>>();
    let spool_path = spill_dir.join(format!("partition_{partition_index:04}.transpose.bin"));
    let mut spool_writer = BufWriter::with_capacity(IO_BUFFER_CAPACITY, File::create(&spool_path)?);
    let mut spool_offset = 0u64;
    let mut records = 0u64;

    let mut reader = BufReader::with_capacity(IO_BUFFER_CAPACITY, File::open(path)?);
    let mut cid = 0u64;
    while let Some(cid_delta) = read_optional_varint_u64(&mut reader)? {
        cid = cid.checked_add(cid_delta).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "CID overflow in transpose spill",
            )
        })?;
        let count = usize::try_from(read_required_varint_u64(&mut reader)?).map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, "dataset count cannot fit usize")
        })?;
        if count == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "empty CID group in dataset transpose spill",
            ));
        }
        let mut dataset = dataset_start;
        for _ in 0..count {
            let dataset_delta =
                usize::try_from(read_required_varint_u64(&mut reader)?).map_err(|_| {
                    io::Error::new(io::ErrorKind::InvalidData, "dataset delta cannot fit usize")
                })?;
            dataset = dataset.checked_add(dataset_delta).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "dataset index overflow")
            })?;
            if dataset < dataset_start || dataset >= dataset_end {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "dataset {} is outside partition {} range [{}, {})",
                        dataset, partition_index, dataset_start, dataset_end
                    ),
                ));
            }
            let state = &mut datasets[dataset - dataset_start];
            let abundance = if abundance_enabled {
                let mut code = [0u8; 1];
                reader.read_exact(&mut code)?;
                Some(code[0])
            } else {
                None
            };
            if state.seen && cid < state.last_cid {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "CIDs for dataset {} are not monotonic: {} before {}",
                        dataset, state.last_cid, cid
                    ),
                ));
            }
            let delta = if state.seen {
                cid - state.last_cid
            } else {
                cid
            };
            write_varint_u64(delta, &mut state.encoded_deltas);
            if let Some(code) = abundance {
                state.encoded_deltas.push(code);
            }
            state.last_cid = cid;
            state.seen = true;
            records += 1;
            if state.encoded_deltas.len() >= CID_TRANSPOSE_BLOCK_BYTES {
                flush_cid_transpose_block(state, &mut spool_writer, &mut spool_offset)?;
            }
        }
    }

    for state in &mut datasets {
        flush_cid_transpose_block(state, &mut spool_writer, &mut spool_offset)?;
    }
    spool_writer.flush()?;
    drop(spool_writer);
    fs::remove_file(path)?;

    let mut spool_reader = File::open(&spool_path)?;
    let mut block_buffer = Vec::new();
    for (local_dataset, state) in datasets.into_iter().enumerate() {
        if !state.seen {
            continue;
        }
        let dataset = dataset_start + local_dataset;
        for (offset, len) in state.blocks {
            spool_reader.seek(SeekFrom::Start(offset))?;
            block_buffer.resize(len as usize, 0);
            spool_reader.read_exact(&mut block_buffer)?;
            payload_writer.append_encoded_deltas(dataset, &block_buffer)?;
        }
    }
    drop(spool_reader);
    fs::remove_file(&spool_path)?;
    Ok((records, spool_offset))
}

struct DatasetPayloadWriter {
    output: BufWriter<File>,
    offsets: Vec<usize>,
    dataset_count: usize,
    total_size: usize,
    current_dataset: Option<usize>,
    cid_payload_path: PathBuf,
    abundance_payload_path: PathBuf,
    cid_encoder: Option<Encoder<'static, File>>,
    abundance_encoder: Option<Encoder<'static, File>>,
    abundance_enabled: bool,
    cid_count: u64,
    current_postings: u64,
    bitmap_payload_path: PathBuf,
    xor_payload_path: PathBuf,
    majority_bitmap_path: Option<PathBuf>,
    delta_datasets: u64,
    bitmap_datasets: u64,
    majority_xor_datasets: u64,
    payload_bytes: u64,
}

impl DatasetPayloadWriter {
    fn new(
        path: &str,
        spill_dir: &Path,
        dataset_count: usize,
        cid_count: usize,
        majority_bitmap_path: &Path,
        abundance_log_base: Option<f64>,
    ) -> Result<Self> {
        let mut output = BufWriter::with_capacity(IO_BUFFER_CAPACITY, File::create(path)?);
        output.write_all(DATASET_TO_CID_MAGIC)?;
        let flags = DATASET_TO_CID_FLAG_ADAPTIVE
            | u8::from(abundance_log_base.is_some()) * DATASET_TO_CID_FLAG_ABUNDANCE
            | u8::from(abundance_log_base.is_none()) * DATASET_TO_CID_FLAG_MAJORITY;
        output.write_all(&[flags])?;
        if let Some(base) = abundance_log_base {
            output.write_all(&base.to_le_bytes())?;
        }
        let cid_count = u64::try_from(cid_count)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "CID count is too large"))?;
        output.write_all(&cid_count.to_le_bytes())?;
        let mut header_size = DATASET_TO_CID_MAGIC.len()
            + 1
            + usize::from(abundance_log_base.is_some()) * std::mem::size_of::<f64>()
            + std::mem::size_of::<u64>();
        let majority_bitmap_path = if abundance_log_base.is_none() {
            let compressed_path = spill_dir.join("majority_membership.zst.tmp");
            let compressed_file = File::options()
                .read(true)
                .write(true)
                .create(true)
                .truncate(true)
                .open(&compressed_path)?;
            let mut encoder = Encoder::new(compressed_file, 1)?;
            let mut majority = File::open(majority_bitmap_path)?;
            io::copy(&mut majority, &mut encoder)?;
            let mut compressed_file = encoder.finish()?;
            let compressed_len = compressed_file.seek(SeekFrom::End(0))?;
            compressed_file.seek(SeekFrom::Start(0))?;
            output.write_all(&compressed_len.to_le_bytes())?;
            io::copy(&mut compressed_file, &mut output)?;
            header_size = header_size
                .checked_add(8 + compressed_len as usize)
                .ok_or_else(|| io::Error::other("dataset posting header size overflow"))?;
            fs::remove_file(compressed_path)?;
            Some(majority_bitmap_path.to_path_buf())
        } else {
            None
        };
        Ok(Self {
            output,
            offsets: Vec::with_capacity(dataset_count),
            dataset_count,
            total_size: header_size,
            current_dataset: None,
            cid_payload_path: spill_dir.join("dataset_cids.zst.tmp"),
            abundance_payload_path: spill_dir.join("dataset_abundance.zst.tmp"),
            cid_encoder: None,
            abundance_encoder: None,
            abundance_enabled: abundance_log_base.is_some(),
            cid_count,
            current_postings: 0,
            bitmap_payload_path: spill_dir.join("dataset_bitmap.zst.tmp"),
            xor_payload_path: spill_dir.join("dataset_majority_xor.zst.tmp"),
            majority_bitmap_path,
            delta_datasets: 0,
            bitmap_datasets: 0,
            majority_xor_datasets: 0,
            payload_bytes: 0,
        })
    }

    fn start_dataset(&mut self, dataset: usize) -> Result<()> {
        if self.current_dataset.is_some() {
            return Err(io::Error::other(
                "cannot start dataset payload before finishing the previous payload",
            ));
        }
        let payload_file = File::options()
            .read(true)
            .write(true)
            .create(true)
            .truncate(true)
            .open(&self.cid_payload_path)?;
        let encoder = Encoder::new(payload_file, 1)?;
        let abundance_encoder = if self.abundance_enabled {
            let file = File::options()
                .read(true)
                .write(true)
                .create(true)
                .truncate(true)
                .open(&self.abundance_payload_path)?;
            Some(Encoder::new(file, 1)?)
        } else {
            None
        };
        self.current_dataset = Some(dataset);
        self.current_postings = 0;
        self.cid_encoder = Some(encoder);
        self.abundance_encoder = abundance_encoder;
        Ok(())
    }

    fn finish_dataset(&mut self) -> Result<()> {
        let dataset = self
            .current_dataset
            .take()
            .ok_or_else(|| io::Error::other("no active dataset payload"))?;
        let encoder = self
            .cid_encoder
            .take()
            .ok_or_else(|| io::Error::other("missing dataset payload encoder"))?;
        let mut cid_file = encoder.finish()?;
        let cid_len = cid_file.seek(SeekFrom::End(0))?;
        cid_file.seek(SeekFrom::Start(0))?;
        self.offsets.push(self.total_size);
        if let Some(encoder) = self.abundance_encoder.take() {
            let mut abundance_file = encoder.finish()?;
            let abundance_len = abundance_file.seek(SeekFrom::End(0))?;
            abundance_file.seek(SeekFrom::Start(0))?;
            self.output.write_all(&[DATASET_POSTING_DELTA])?;
            self.output.write_all(&cid_len.to_le_bytes())?;
            io::copy(&mut cid_file, &mut self.output)?;
            io::copy(&mut abundance_file, &mut self.output)?;
            self.total_size = self
                .total_size
                .checked_add(1 + 8)
                .and_then(|size| size.checked_add(cid_len as usize))
                .and_then(|size| size.checked_add(abundance_len as usize))
                .ok_or_else(|| io::Error::other("id-to-CID output size overflow"))?;
            self.payload_bytes = self
                .payload_bytes
                .saturating_add(cid_len)
                .saturating_add(abundance_len);
            self.delta_datasets += 1;
        } else if let Some((codec, mut bitmap_file, bitmap_len)) =
            self.build_smaller_bitmap_payload(cid_len)?
        {
            self.output.write_all(&[codec])?;
            io::copy(&mut bitmap_file, &mut self.output)?;
            self.total_size = self
                .total_size
                .checked_add(1 + bitmap_len as usize)
                .ok_or_else(|| io::Error::other("id-to-CID output size overflow"))?;
            self.payload_bytes = self.payload_bytes.saturating_add(bitmap_len);
            if codec == DATASET_POSTING_BITMAP {
                self.bitmap_datasets += 1;
            } else {
                self.majority_xor_datasets += 1;
            }
        } else {
            self.output.write_all(&[DATASET_POSTING_DELTA])?;
            io::copy(&mut cid_file, &mut self.output)?;
            self.total_size = self
                .total_size
                .checked_add(1 + cid_len as usize)
                .ok_or_else(|| io::Error::other("id-to-CID output size overflow"))?;
            self.payload_bytes = self.payload_bytes.saturating_add(cid_len);
            self.delta_datasets += 1;
        }
        debug_assert_eq!(self.offsets.len(), dataset + 1);
        Ok(())
    }

    fn build_smaller_bitmap_payload(&self, delta_len: u64) -> Result<Option<(u8, File, u64)>> {
        let bitmap_bytes = self.cid_count.div_ceil(8);
        // A full bitmap scan is allowed only when it is bounded by eight times
        // the number of returned CIDs. This preserves O(output) query time.
        if self.cid_count == 0
            || self.current_postings.saturating_mul(8) < self.cid_count
            || bitmap_bytes > MAX_DATASET_POSTING_BITMAP_BYTES
        {
            return Ok(None);
        }
        let mut bitmap = vec![0u8; bitmap_bytes as usize];
        // Reopen rather than clone: duplicated Unix file descriptors share a
        // cursor, which would otherwise consume the delta payload selected as
        // the fallback representation.
        let source = File::open(&self.cid_payload_path)?;
        let mut decoder = zstd::Decoder::new(source)?.single_frame();
        let mut cid = 0u64;
        while let Some(delta) = read_optional_varint_u64(&mut decoder)? {
            cid = cid.checked_add(delta).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "dataset posting CID overflow")
            })?;
            if cid >= self.cid_count {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "dataset posting CID outside bitmap universe",
                ));
            }
            bitmap[(cid / 8) as usize] |= 1u8 << (cid % 8);
        }
        let bitmap_file = File::options()
            .read(true)
            .write(true)
            .create(true)
            .truncate(true)
            .open(&self.bitmap_payload_path)?;
        let mut encoder = Encoder::new(bitmap_file, 1)?;
        encoder.write_all(&bitmap)?;
        let mut bitmap_file = encoder.finish()?;
        let bitmap_len = bitmap_file.seek(SeekFrom::End(0))?;
        bitmap_file.seek(SeekFrom::Start(0))?;
        let mut best_codec = DATASET_POSTING_BITMAP;
        let mut best_len = bitmap_len;
        if let Some(majority_path) = self.majority_bitmap_path.as_ref() {
            let mut majority = File::open(majority_path)?;
            let mut reference = vec![0u8; 1024 * 1024];
            let mut offset = 0usize;
            while offset < bitmap.len() {
                let chunk_len = (bitmap.len() - offset).min(reference.len());
                majority.read_exact(&mut reference[..chunk_len])?;
                for (value, base) in bitmap[offset..offset + chunk_len]
                    .iter_mut()
                    .zip(&reference[..chunk_len])
                {
                    *value ^= *base;
                }
                offset += chunk_len;
            }
            let xor_file = File::options()
                .read(true)
                .write(true)
                .create(true)
                .truncate(true)
                .open(&self.xor_payload_path)?;
            let mut encoder = Encoder::new(xor_file, 1)?;
            encoder.write_all(&bitmap)?;
            let xor_file = encoder.finish()?;
            let xor_len = xor_file.metadata()?.len();
            if xor_len < best_len {
                best_codec = DATASET_POSTING_MAJORITY_XOR;
                best_len = xor_len;
            }
        }
        if best_len >= delta_len {
            return Ok(None);
        }
        let selected_path = if best_codec == DATASET_POSTING_BITMAP {
            &self.bitmap_payload_path
        } else {
            &self.xor_payload_path
        };
        let mut selected = File::open(selected_path)?;
        selected.seek(SeekFrom::Start(0))?;
        Ok(Some((best_codec, selected, best_len)))
    }

    fn advance_to(&mut self, dataset: usize) -> Result<()> {
        if let Some(current) = self.current_dataset {
            if dataset < current {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "dataset spill order regressed from {} to {}",
                        current, dataset
                    ),
                ));
            }
            if dataset == current {
                return Ok(());
            }
            self.finish_dataset()?;
        }
        while self.offsets.len() < dataset {
            let missing = self.offsets.len();
            self.start_dataset(missing)?;
            self.finish_dataset()?;
        }
        self.start_dataset(dataset)
    }

    fn append_encoded_deltas(&mut self, dataset: usize, encoded_deltas: &[u8]) -> Result<()> {
        if dataset >= self.dataset_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "dataset {} outside [0, {}) in transposed dataset-to-CID spill",
                    dataset, self.dataset_count
                ),
            ));
        }
        self.advance_to(dataset)?;
        if self.abundance_enabled {
            let mut cursor = std::io::Cursor::new(encoded_deltas);
            while let Some(delta) = read_optional_varint_u64(&mut cursor)? {
                self.current_postings = self.current_postings.saturating_add(1);
                write_varint_u64_to_writer(
                    delta,
                    self.cid_encoder
                        .as_mut()
                        .ok_or_else(|| io::Error::other("missing dataset CID encoder"))?,
                )?;
                let mut deviation = [0u8; 1];
                cursor.read_exact(&mut deviation)?;
                self.abundance_encoder
                    .as_mut()
                    .ok_or_else(|| io::Error::other("missing dataset abundance encoder"))?
                    .write_all(&deviation)?;
            }
            Ok(())
        } else {
            let mut cursor = std::io::Cursor::new(encoded_deltas);
            while read_optional_varint_u64(&mut cursor)?.is_some() {
                self.current_postings = self.current_postings.saturating_add(1);
            }
            self.cid_encoder
                .as_mut()
                .ok_or_else(|| io::Error::other("missing dataset CID encoder"))?
                .write_all(encoded_deltas)
        }
    }

    fn finish(mut self, dataset_count: usize) -> Result<(Vec<usize>, u64, u64, u64, u64)> {
        if self.current_dataset.is_some() {
            self.finish_dataset()?;
        }
        while self.offsets.len() < dataset_count {
            let dataset = self.offsets.len();
            self.start_dataset(dataset)?;
            self.finish_dataset()?;
        }
        if self.offsets.len() != dataset_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "dataset-to-CID spill contains an out-of-range dataset",
            ));
        }
        // A terminal boundary lets readers derive each compressed frame length
        // from two adjacent index entries. This avoids one redundant u64 per
        // dataset in the payload.
        self.offsets.push(self.total_size);
        self.output.flush()?;
        if self.cid_payload_path.exists() {
            fs::remove_file(&self.cid_payload_path)?;
        }
        if self.abundance_payload_path.exists() {
            fs::remove_file(&self.abundance_payload_path)?;
        }
        if self.bitmap_payload_path.exists() {
            fs::remove_file(&self.bitmap_payload_path)?;
        }
        if self.xor_payload_path.exists() {
            fs::remove_file(&self.xor_payload_path)?;
        }
        Ok((
            self.offsets,
            self.payload_bytes,
            self.delta_datasets,
            self.bitmap_datasets,
            self.majority_xor_datasets,
        ))
    }
}

fn write_dataset_to_cid_from_partitions(
    cid_file_path: String,
    spill: DatasetCidSpillFiles,
    _memory_budget: CompressionMemoryBudget,
) -> std::io::Result<Vec<usize>> {
    let total_timer = PhaseTimer::start();
    let mut payload_writer = DatasetPayloadWriter::new(
        &cid_file_path,
        &spill.directory,
        spill.dataset_count,
        spill.cid_count,
        &spill.majority_path,
        spill.abundance_log_base,
    )?;
    let mut transposed_records = 0u64;
    let mut transpose_bytes = 0u64;
    for (partition_index, path) in spill.paths.iter().enumerate() {
        let (records, bytes) = transpose_cid_partition(
            path,
            partition_index,
            &spill.directory,
            spill.datasets_per_partition,
            spill.dataset_count,
            spill.abundance_log_base.is_some(),
            &mut payload_writer,
        )?;
        transposed_records = transposed_records.saturating_add(records);
        transpose_bytes = transpose_bytes.saturating_add(bytes);
    }
    let (offsets, total_payload_bytes, delta_datasets, bitmap_datasets, majority_xor_datasets) =
        payload_writer.finish(spill.dataset_count)?;
    if spill.majority_path.exists() {
        fs::remove_file(&spill.majority_path)?;
    }
    let _ = fs::remove_dir(&spill.directory);
    log_phase_timing("post_ggcat.dataset_to_cid.total", total_timer.finish());
    println!(
        "[phase-stats] phase=post_ggcat.dataset_to_cid datasets={} partitions={} datasets_per_partition={} transposed_records={} transpose_bytes={} payload_bytes={} delta_datasets={} bitmap_datasets={} majority_xor_datasets={}",
        spill.dataset_count,
        spill.paths.len(),
        spill.datasets_per_partition,
        transposed_records,
        transpose_bytes,
        total_payload_bytes,
        delta_datasets,
        bitmap_datasets,
        majority_xor_datasets
    );
    Ok(offsets)
}

pub(crate) fn sort_by_bucket_streaming(
    output_dir: &String,
    nb_files: u32,
    worker_threads: usize,
    memory_gb: usize,
    abundance_log_base: Option<f64>,
    record_rx: mpsc::Receiver<SimplitigBatch>,
) -> (u64, Vec<usize>) {
    let total_timer = PhaseTimer::start();
    let memory_budget = CompressionMemoryBudget::from_gb(memory_gb);
    println!("Starting writing compressed sequences (streaming).");
    let write_stream_timer = PhaseTimer::start();
    let triple = match write_compressed_from_stream(
        output_dir.clone() + "tigs_kloe.fa",
        output_dir,
        nb_files,
        worker_threads,
        memory_budget,
        abundance_log_base,
        record_rx,
    ) {
        Ok(res_pair) => res_pair,
        Err(e) => panic!("Error writing compressed unitigs: {e:?}"),
    };
    let write_stream_timing = write_stream_timer.finish();
    println!(
        "Writing compressed sequences wall time: {:.3}s",
        write_stream_timing.wall_sec
    );
    log_phase_timing("post_ggcat.write_stream", write_stream_timing);
    let position_timer = PhaseTimer::start();
    println!("Starting to write positions");
    let position_entries = triple.position_entries;
    if let Err(e) = write_positions_from_spill(
        &triple.position_path,
        triple.position_entries,
        String::from(output_dir.clone() + "positions_kloe.bin"),
        String::from(output_dir.clone() + POSITIONS_INDEX_FILE),
    ) {
        panic!("Error writting positions: {e:?}");
    }
    let position_timing = position_timer.finish();
    println!(
        "Write positions wall time: {:.3}s",
        position_timing.wall_sec
    );
    log_phase_timing("post_ggcat.write_positions", position_timing);
    let transpose_timer = PhaseTimer::start();
    let dataset_offsets = match write_dataset_to_cid_from_partitions(
        output_dir.clone() + DATASET_TO_CID_FILE,
        triple.dataset_cid_spill,
        memory_budget,
    ) {
        Ok(offsets) => offsets,
        Err(e) => panic!("Error writing dataset-to-CID postings: {e:?}"),
    };
    log_phase_timing("post_ggcat.dataset_to_cid", transpose_timer.finish());
    let _ = fs::remove_dir(&triple.spill_directory);
    let total_timing = total_timer.finish();
    println!("Compression took: {:.3}s", total_timing.wall_sec);
    log_phase_timing("post_ggcat.total", total_timing);
    (position_entries, dataset_offsets)
}

fn stream_sorted_records_from_chunks(
    chunk_paths: &[PathBuf],
    mut emit: impl FnMut(SortedColorRecord) -> Result<()>,
) -> Result<()> {
    let mut readers = Vec::with_capacity(chunk_paths.len());
    for path in chunk_paths {
        readers.push(BufReader::with_capacity(
            COLOR_RUN_BUFFER_BYTES,
            File::open(path)?,
        ));
    }

    let mut heap = BinaryHeap::new();
    for (idx, reader) in readers.iter_mut().enumerate() {
        if let Some(record) = read_sorted_color_record(reader)? {
            heap.push(ChunkHeapItem {
                chunk_index: idx,
                record,
            });
        }
    }

    while let Some(item) = heap.pop() {
        emit(item.record)?;
        if let Some(next) = read_sorted_color_record(&mut readers[item.chunk_index])? {
            heap.push(ChunkHeapItem {
                chunk_index: item.chunk_index,
                record: next,
            });
        }
    }

    Ok(())
}

fn reduce_sorted_color_chunks(mut chunks: Vec<PathBuf>, chunk_dir: &Path) -> Result<Vec<PathBuf>> {
    let mut merge_pass = 0usize;
    while chunks.len() > COLOR_MERGE_FAN_IN {
        let mut next_chunks = Vec::with_capacity(chunks.len().div_ceil(COLOR_MERGE_FAN_IN));
        for (group_index, group) in chunks.chunks(COLOR_MERGE_FAN_IN).enumerate() {
            if group.len() == 1 {
                next_chunks.push(group[0].clone());
                continue;
            }
            let path = chunk_dir.join(format!("color_merge_{merge_pass:03}_{group_index:06}.bin"));
            let mut writer = BufWriter::with_capacity(IO_BUFFER_CAPACITY, File::create(&path)?);
            stream_sorted_records_from_chunks(group, |record| {
                write_sorted_color_record(&mut writer, &record)
            })?;
            writer.flush()?;
            for old_path in group {
                fs::remove_file(old_path)?;
            }
            next_chunks.push(path);
        }
        chunks = next_chunks;
        merge_pass += 1;
    }
    Ok(chunks)
}

struct StructuredGgcatCapture {
    chunk_files: Vec<PathBuf>,
    chunk_flush_metrics: ChunkFlushMetrics,
    input_sequences_count: usize,
    emitted_segments_count: usize,
    initial_chunk_count: usize,
}

fn capture_structured_ggcat_output(
    receiver: mpsc::Receiver<Vec<u8>>,
    chunk_dir: PathBuf,
    memory_budget: CompressionMemoryBudget,
    k: usize,
    worker_threads: usize,
) -> Result<StructuredGgcatCapture> {
    // Sort/write one half-budget chunk while GGCAT fills the other.  A rendezvous
    // channel permits exactly one chunk in each stage, so pipelining does not
    // increase the previous peak record-buffer memory.
    let chunk_target = memory_budget.color_chunk_bytes().div_ceil(2).max(1);
    let sort_threads = (worker_threads / 4).max(1);
    let (sort_tx, sort_rx) = mpsc::sync_channel::<Vec<SortedColorRecord>>(0);
    let sort_chunk_dir = chunk_dir.clone();
    let sort_thread = thread::spawn(move || -> Result<(Vec<PathBuf>, ChunkFlushMetrics)> {
        let sort_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(sort_threads)
            .thread_name(|index| format!("kloe-color-sort-{index}"))
            .build()
            .map_err(|err| io::Error::other(format!("create color-sort worker pool: {err}")))?;
        let mut chunk_files = Vec::new();
        let mut chunk_flush_metrics = ChunkFlushMetrics::default();
        for mut chunk in sort_rx {
            flush_sorted_chunk_timed(
                &mut chunk,
                &mut chunk_files,
                &sort_chunk_dir,
                &sort_pool,
                &mut chunk_flush_metrics,
            )?;
        }
        Ok((chunk_files, chunk_flush_metrics))
    });

    let mut records = Vec::new();
    let mut chunk_bytes = 0usize;
    let mut input_sequences_count = 0usize;
    let mut emitted_segments_count = 0usize;

    let capture_result = (|| -> Result<()> {
        for block in receiver {
            let (sequences, segments) = visit_structured_ggcat_block(&block, k, |record| {
                chunk_bytes = chunk_bytes
                    .saturating_add(std::mem::size_of::<ColorIndexType>() + record.seq.len());
                records.push(record);
                if chunk_bytes >= chunk_target {
                    sort_tx.send(std::mem::take(&mut records)).map_err(|_| {
                        io::Error::other("structured color-sort worker disconnected")
                    })?;
                    chunk_bytes = 0;
                }
                Ok(())
            })?;
            input_sequences_count = input_sequences_count.saturating_add(sequences);
            emitted_segments_count = emitted_segments_count.saturating_add(segments);
        }
        if !records.is_empty() {
            sort_tx
                .send(std::mem::take(&mut records))
                .map_err(|_| io::Error::other("structured color-sort worker disconnected"))?;
        }
        Ok(())
    })();
    drop(sort_tx);
    let sort_result = sort_thread
        .join()
        .map_err(|_| io::Error::other("structured color-sort worker panicked"))?;
    capture_result?;
    let (chunk_files, chunk_flush_metrics) = sort_result?;
    let initial_chunk_count = chunk_files.len();
    let chunk_files = reduce_sorted_color_chunks(chunk_files, &chunk_dir)?;
    Ok(StructuredGgcatCapture {
        chunk_files,
        chunk_flush_metrics,
        input_sequences_count,
        emitted_segments_count,
        initial_chunk_count,
    })
}

fn produce_ggcat_records(
    input_streams: Vec<GeneralSequenceBlockData>,
    source_dataset_ids: SourceDatasetMap,
    dataset_count: usize,
    threads: usize,
    k: usize,
    m: usize,
    use_unitigs: bool,
    use_matchtigs: bool,
    use_eulertigs: bool,
    ggcat_cfg: GgcatCompressionConfig,
    sender: mpsc::SyncSender<SimplitigBatch>,
) -> Result<()> {
    let source_count = source_dataset_ids.len();
    let memory_budget = CompressionMemoryBudget::from_gb(ggcat_cfg.memory_gb);
    if source_count == 0 || input_streams.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "ggcat source list is empty",
        ));
    }
    let ggcat_temp_root = if ggcat_cfg.temp_dir.is_empty() {
        std::env::temp_dir().join(format!("kloe-ggcat-shared-{}", std::process::id()))
    } else {
        PathBuf::from(&ggcat_cfg.temp_dir)
    };
    fs::create_dir_all(&ggcat_temp_root)?;

    let pid = std::process::id();
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);

    let structured_output_key =
        ggcat_temp_root.join(format!("kloe-color-stream-{pid}-{nanos}.tmp"));
    let chunk_dir = ggcat_temp_root.join(format!("kloe-color-chunks-{pid}-{nanos}"));
    fs::create_dir_all(&chunk_dir)?;

    let instance = GGCATInstance::create(GGCATConfig {
        temp_dir: Some(ggcat_temp_root.clone()),
        memory: memory_budget.ggcat_memory_gb(),
        prefer_memory: false,
        total_threads_count: threads.max(1),
        intermediate_compression_level: None,
        stats_file: None,
        messages_callback: None,
    })
    .map_err(|e| to_io_err("create ggcat instance", e))?;

    println!(
        "Running embedded ggcat structured-output path (k={}, m={}, threads={}, sources={}, datasets={}, total_memory={}GB, ggcat_memory={:.2}GB, disk_backed=true)",
        k,
        m,
        threads.max(1),
        source_count,
        dataset_count,
        ggcat_cfg.memory_gb.max(1),
        memory_budget.ggcat_memory_gb()
    );

    let (structured_tx, structured_rx) = mpsc::sync_channel::<Vec<u8>>(2);
    if !register_channel_output(structured_output_key.clone(), structured_tx) {
        return Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            format!(
                "a GGCAT structured output is already registered for {}",
                structured_output_key.display()
            ),
        ));
    }
    let capture_chunk_dir = chunk_dir.clone();
    let capture_thread = thread::spawn(move || {
        capture_structured_ggcat_output(structured_rx, capture_chunk_dir, memory_budget, k, threads)
    });

    let ggcat_build_timer = PhaseTimer::start();
    let build_result = instance.build_graph(
        input_streams,
        structured_output_key.clone(),
        None,
        Some(source_count),
        k,
        threads.max(1),
        false,
        Some(m),
        true,
        1,
        ggcat_extra_elaboration(use_unitigs, use_matchtigs, use_eulertigs),
        None,
    );
    unregister_channel_output(&structured_output_key);
    let capture_result = capture_thread
        .join()
        .map_err(|_| io::Error::other("GGCAT structured-output capture thread panicked"))?;
    let records_output = build_result.map_err(|e| to_io_err("ggcat build_graph", e))?;
    let captured = capture_result?;
    log_phase_timing("ggcat.build_graph_api_call", ggcat_build_timer.finish());

    let colormap_file = GGCATInstance::get_colormap_file(&records_output);
    let colors_deserializer =
        ColorsDeserializer::<DefaultColorsSerializer>::new(&colormap_file, false)
            .map_err(|e| to_io_err("open ggcat colormap", e))?;

    let colors_count = colors_deserializer.colors_count();
    if colors_count != source_count {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "ggcat returned {colors_count} colors but KLOE registered {source_count} sources"
            ),
        ));
    }

    let keep_ggcat = std::env::var_os("KLOE_KEEP_GGCAT").is_some();
    let process_result = (|| -> Result<()> {
        let producer_total_timer = PhaseTimer::start();
        let StructuredGgcatCapture {
            chunk_files,
            chunk_flush_metrics,
            input_sequences_count,
            emitted_segments_count,
            initial_chunk_count,
        } = captured;

        let resolve_and_emit_timer = PhaseTimer::start();
        let mut batch: SimplitigBatch = Vec::with_capacity(COLOR_RECORD_BATCH_SIZE);
        let mut batch_bytes = 0usize;
        let mut record_window = Vec::<SortedColorRecord>::new();
        let mut record_window_bytes = 0usize;
        let mut record_window_subsets = 0usize;
        let mut window_last_subset = None;
        let mut global_last_subset = None;
        let mut unique_subsets_count = 0usize;
        let mut merged_records_count = 0usize;
        let mut emitted_batches_count = 0usize;
        let mut emitted_records_count = 0usize;
        let window_sequence_limit = memory_budget.subset_window_sequence_bytes();
        let window_subset_limit = memory_budget.subset_window_count(dataset_count);

        stream_sorted_records_from_chunks(&chunk_files, |record| {
            merged_records_count += 1;
            if global_last_subset != Some(record.subset) {
                unique_subsets_count += 1;
                global_last_subset = Some(record.subset);
            }
            let is_new_window_subset = window_last_subset != Some(record.subset);
            if !record_window.is_empty()
                && (record_window_bytes >= window_sequence_limit
                    || (is_new_window_subset && record_window_subsets >= window_subset_limit))
            {
                emit_resolved_record_window(
                    &mut record_window,
                    &instance,
                    &colormap_file,
                    None,
                    &source_dataset_ids,
                    dataset_count,
                    &sender,
                    &mut batch,
                    &mut batch_bytes,
                    window_sequence_limit,
                    &mut emitted_batches_count,
                    &mut emitted_records_count,
                )?;
                record_window_bytes = 0;
                record_window_subsets = 0;
                window_last_subset = None;
            }
            if window_last_subset != Some(record.subset) {
                record_window_subsets += 1;
                window_last_subset = Some(record.subset);
            }
            record_window_bytes = record_window_bytes.saturating_add(record.seq.len());
            record_window.push(record);
            Ok(())
        })?;

        emit_resolved_record_window(
            &mut record_window,
            &instance,
            &colormap_file,
            None,
            &source_dataset_ids,
            dataset_count,
            &sender,
            &mut batch,
            &mut batch_bytes,
            window_sequence_limit,
            &mut emitted_batches_count,
            &mut emitted_records_count,
        )?;

        if !batch.is_empty() {
            emitted_batches_count += 1;
            emitted_records_count += batch.len();
            sender.send(batch).map_err(|e| {
                io::Error::new(
                    io::ErrorKind::BrokenPipe,
                    format!("record receiver dropped: {e}"),
                )
            })?;
        }

        let resolve_and_emit_timing = resolve_and_emit_timer.finish();
        log_phase_timing(
            "ggcat.structured_capture_chunk_flush_total",
            PhaseTiming {
                wall_sec: chunk_flush_metrics.timing.wall_sec,
                cpu_sec: chunk_flush_metrics.timing.cpu_sec,
            },
        );
        log_phase_timing(
            "producer.resolve_subsets_and_emit_records",
            resolve_and_emit_timing,
        );
        log_phase_timing("producer.total_post_ggcat", producer_total_timer.finish());
        println!(
            "[phase-stats] phase=producer input_sequences={} emitted_segments={} subsets_unique={} chunks={} chunk_flush_calls={} chunk_flush_records={} merged_records={} emitted_batches={} emitted_records={}",
            input_sequences_count,
            emitted_segments_count,
            unique_subsets_count,
            initial_chunk_count,
            chunk_flush_metrics.calls,
            chunk_flush_metrics.records,
            merged_records_count,
            emitted_batches_count,
            emitted_records_count
        );

        drop(sender);
        Ok(())
    })();

    if keep_ggcat {
        eprintln!(
            "DEBUG_KEEP_GGCAT structured_output=no_file colormap={} chunks_dir={}",
            colormap_file.display(),
            chunk_dir.display()
        );
    } else {
        let _ = fs::remove_file(&records_output);
        let _ = fs::remove_file(&colormap_file);
        let _ = fs::remove_dir_all(&chunk_dir);
    }

    process_result
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn compress_ggcat_sources(
    output_dir: &str,
    filenames: Vec<String>,
    input_streams: Vec<GeneralSequenceBlockData>,
    source_dataset_ids: SourceDatasetMap,
    threads: usize,
    k: usize,
    m: usize,
    use_unitigs: bool,
    use_matchtigs: bool,
    use_eulertigs: bool,
    ggcat_cfg: GgcatCompressionConfig,
) -> Result<()> {
    let dataset_count = filenames.len();
    if dataset_count == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "input file list is empty",
        ));
    }

    println!("Compression backend: embedded ggcat structured output");

    let (record_tx, record_rx) = mpsc::sync_channel::<SimplitigBatch>(2);
    let producer_cfg = ggcat_cfg.clone();
    let abundance_log_base = source_dataset_ids.abundance_log_base();

    let producer = thread::spawn(move || {
        produce_ggcat_records(
            input_streams,
            source_dataset_ids,
            dataset_count,
            threads,
            k,
            m,
            use_unitigs,
            use_matchtigs,
            use_eulertigs,
            producer_cfg,
            record_tx,
        )
    });

    let sort_start = Instant::now();
    let (position_entries, dataset_offsets) = sort_by_bucket_streaming(
        &output_dir.to_string(),
        dataset_count as u32,
        threads,
        ggcat_cfg.memory_gb,
        abundance_log_base,
        record_rx,
    );
    println!(
        "Streamed ggcat colored records into KLOE archive in {:.3}s",
        sort_start.elapsed().as_secs_f64()
    );

    match producer.join() {
        Ok(Ok(())) => {}
        Ok(Err(err)) => {
            return Err(io::Error::other(format!("ggcat producer failed: {err}")));
        }
        Err(_) => {
            return Err(io::Error::other("ggcat producer panicked"));
        }
    }

    let cid_count = position_entries.checked_sub(1).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "archive positions omit the initial boundary",
        )
    })?;
    let mut manifest = BufWriter::new(File::create(
        Path::new(output_dir).join(ARCHIVE_MANIFEST_FILE),
    )?);
    manifest.write_all(ARCHIVE_MANIFEST_MAGIC)?;
    manifest.write_all(&4u32.to_le_bytes())?;
    manifest.write_all(&(k as u32).to_le_bytes())?;
    manifest.write_all(&(m as u32).to_le_bytes())?;
    manifest.write_all(&(dataset_count as u64).to_le_bytes())?;
    manifest.write_all(&cid_count.to_le_bytes())?;
    manifest.write_all(&u32::from(abundance_log_base.is_some()).to_le_bytes())?;
    manifest.write_all(&abundance_log_base.unwrap_or(0.0).to_le_bytes())?;
    manifest.flush()?;

    // Dataset-to-CID is the sole persistent membership relation in v4.  Merge
    // reconstructs the inverse as a bounded temporary external transpose.
    write_filenames_id_offsets(output_dir, &filenames, &dataset_offsets)?;
    Ok(())
}

pub fn compress_with_ggcat(
    output_dir: &String,
    input_fof: &String,
    threads: usize,
    k: usize,
    m: usize,
    partition_power: u32,
    verify_kmers: bool,
    skip_sort: bool,
    use_unitigs: bool,
    use_matchtigs: bool,
    use_eulertigs: bool,
    ggcat_cfg: GgcatCompressionConfig,
) -> Result<()> {
    let _unused_partition_power = partition_power;
    if verify_kmers {
        eprintln!("Warning: --verify-kmers is currently ignored in the ggcat hard-switch backend.");
    }
    if skip_sort {
        eprintln!(
            "Warning: --skip-sort is ignored; KLOE groups GGCAT structured records by colors."
        );
    }

    let filenames = read_input_fof_filenames(input_fof)?;
    let (input_streams, source_dataset_ids) = if let Some(log_base) = ggcat_cfg.abundance_log_base {
        let stream = Arc::new(AbundanceSequencesStream::new(&filenames, log_base));
        let stream: Arc<dyn DynamicSequencesStream> = stream;
        let input_streams = (0..filenames.len())
            .map(|block| GeneralSequenceBlockData::Dynamic((Arc::clone(&stream), block)))
            .collect::<Vec<_>>();
        println!("Abundance preservation enabled: 8-bit logarithmic codes, base={log_base}");
        (
            input_streams,
            SourceDatasetMap::from_abundance_datasets(filenames.len(), log_base)?,
        )
    } else {
        let input_streams = filenames
            .iter()
            .map(|file| {
                let path = PathBuf::from(file);
                let resolved = fs::canonicalize(&path).unwrap_or(path);
                GeneralSequenceBlockData::FASTA((resolved, None))
            })
            .collect::<Vec<_>>();
        (
            input_streams,
            SourceDatasetMap::from_singletons(filenames.len())?,
        )
    };

    compress_ggcat_sources(
        output_dir,
        filenames,
        input_streams,
        source_dataset_ids,
        threads,
        k,
        m,
        use_unitigs,
        use_matchtigs,
        use_eulertigs,
        ggcat_cfg,
    )
}

#[cfg(test)]
mod structured_output_tests {
    use super::*;
    use ggcat_colors::storage::run_length::RunLengthColorsSerializer;
    use ggcat_colors::storage::serializer::ColorsSerializer;

    #[test]
    fn source_color_union_is_linear_sorted_and_deduplicated() {
        let mut ids = vec![1, 3, 7, 9];
        merge_sorted_dataset_ids(&mut ids, &[2, 3, 8, 9, 10]);
        assert_eq!(ids, vec![1, 2, 3, 7, 8, 9, 10]);

        merge_sorted_dataset_ids_with_offset(&mut ids, &[1, 2, 5], 10).unwrap();
        assert_eq!(ids, vec![1, 2, 3, 7, 8, 9, 10, 11, 12, 15]);
    }

    #[test]
    fn cid_dataset_sidecar_roundtrips_and_applies_archive_offset() {
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join(CID_TO_DATASET_FILE);
        let mut writer = CidDatasetSidecarWriter::create(&path).unwrap();
        writer.append_one_based(&[1, 3, 9]).unwrap();
        writer.append_zero_based(&[1, 4]).unwrap();
        writer.finish().unwrap();

        let sidecar = Arc::new(CidDatasetSidecar::open(&path).unwrap());
        assert_eq!(sidecar.len(), 2);
        let mut map = SourceDatasetMap::new();
        map.add_disk_storage(sidecar, 10).unwrap();

        let mut ids = vec![2, 11, 20];
        map.merge_source_into(0, &mut ids).unwrap();
        assert_eq!(ids, vec![2, 11, 13, 19, 20]);
        let mut ids = Vec::new();
        map.merge_source_into(1, &mut ids).unwrap();
        assert_eq!(ids, vec![12, 15]);
    }

    #[test]
    fn cid_dataset_sidecar_reads_legacy_offset_blocks() {
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("legacy-cid-to-dataset.bin");
        let payload = vec![0, 2, 6, 1, 3];
        let compressed = zstd::encode_all(payload.as_slice(), 1).unwrap();
        let mut writer = BufWriter::new(File::create(&path).unwrap());
        writer.write_all(CID_TO_DATASET_MAGIC).unwrap();
        writer.write_all(&2u32.to_le_bytes()).unwrap();
        writer
            .write_all(&(compressed.len() as u64).to_le_bytes())
            .unwrap();
        for offset in [0u32, 3, 5] {
            writer.write_all(&offset.to_le_bytes()).unwrap();
        }
        writer.write_all(&compressed).unwrap();
        writer.write_all(&0u32.to_le_bytes()).unwrap();
        writer.flush().unwrap();

        let sidecar = CidDatasetSidecar::open(&path).unwrap();
        let mut groups = Vec::new();
        sidecar
            .visit_groups(|_, ids| {
                groups.push(ids.to_vec());
                Ok(())
            })
            .unwrap();
        assert_eq!(groups, vec![vec![1, 3, 9], vec![2, 5]]);
    }

    #[test]
    fn abundance_sidecar_roundtrips_codes_aligned_with_datasets() {
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("abundance-cid-to-dataset.bin");
        let mut writer =
            CidDatasetSidecarWriter::create_with_abundance(&path, Some(1.05), Some(9)).unwrap();
        writer
            .append_zero_based_with_abundance(&[0, 4, 8], Some(&[7, 19, 31]))
            .unwrap();
        writer
            .append_zero_based_with_abundance(&[2], Some(&[42]))
            .unwrap();
        writer
            .append_zero_based_with_abundance(
                &[0, 1, 2, 3, 4, 5, 6, 7, 8],
                Some(&[11, 11, 11, 11, 11, 11, 11, 11, 11]),
            )
            .unwrap();
        writer
            .append_zero_based_with_abundance(
                &[0, 1, 2, 3, 4, 5, 6, 7, 8],
                Some(&[12, 12, 12, 12, 12, 12, 12, 12, 12]),
            )
            .unwrap();
        writer.finish().unwrap();

        let sidecar = CidDatasetSidecar::open(&path).unwrap();
        assert_eq!(sidecar.abundance_log_base(), Some(1.05));
        let mut groups = Vec::new();
        sidecar
            .visit_groups_with_abundance(|_, ids, codes| {
                groups.push((ids.to_vec(), codes.unwrap().to_vec()));
                Ok(())
            })
            .unwrap();
        assert_eq!(
            groups,
            vec![
                (vec![1, 5, 9], vec![7, 19, 31]),
                (vec![3], vec![42]),
                (
                    vec![1, 2, 3, 4, 5, 6, 7, 8, 9],
                    vec![11, 11, 11, 11, 11, 11, 11, 11, 11]
                ),
                (
                    vec![1, 2, 3, 4, 5, 6, 7, 8, 9],
                    vec![12, 12, 12, 12, 12, 12, 12, 12, 12]
                )
            ]
        );
        let mut merged = Vec::new();
        sidecar.merge_group_into(3, &mut merged, 10).unwrap();
        assert_eq!(merged, vec![11, 12, 13, 14, 15, 16, 17, 18, 19]);
    }

    #[test]
    fn abundance_sidecar_reads_legacy_interleaved_blocks() {
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("legacy-abundance-cid-to-dataset.bin");
        let mut payload = Vec::new();
        write_varint_u64(6, &mut payload);
        payload.extend_from_slice(&[0, 7, 4, 19, 4, 31]);
        write_varint_u64(2, &mut payload);
        payload.extend_from_slice(&[2, 42]);
        let compressed = zstd::encode_all(payload.as_slice(), 1).unwrap();
        let mut writer = BufWriter::new(File::create(&path).unwrap());
        writer.write_all(CID_TO_DATASET_ABUNDANCE_MAGIC).unwrap();
        writer.write_all(&1.05f64.to_le_bytes()).unwrap();
        writer.write_all(&2u32.to_le_bytes()).unwrap();
        writer
            .write_all(&(compressed.len() as u64).to_le_bytes())
            .unwrap();
        writer.write_all(&compressed).unwrap();
        writer.write_all(&0u32.to_le_bytes()).unwrap();
        writer.flush().unwrap();

        let sidecar = CidDatasetSidecar::open(&path).unwrap();
        let mut groups = Vec::new();
        sidecar
            .visit_groups_with_abundance(|_, ids, codes| {
                groups.push((ids.to_vec(), codes.unwrap().to_vec()));
                Ok(())
            })
            .unwrap();
        assert_eq!(
            groups,
            vec![(vec![1, 5, 9], vec![7, 19, 31]), (vec![3], vec![42])]
        );
    }

    #[test]
    fn abundance_sidecar_reads_group_first_residual_blocks() {
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("first-residual-abundance.bin");
        let membership = zstd::encode_all([6, 0, 4, 4].as_slice(), 1).unwrap();
        let abundance = zstd::encode_all([7, 12, 24].as_slice(), 1).unwrap();
        let abundance_len = (abundance.len() as u64)
            | ((ABUNDANCE_CODEC_GROUP_FIRST as u64) << ABUNDANCE_CODEC_SHIFT);
        let mut writer = BufWriter::new(File::create(&path).unwrap());
        writer
            .write_all(CID_TO_DATASET_ADAPTIVE_ABUNDANCE_MAGIC)
            .unwrap();
        writer.write_all(&1.05f64.to_le_bytes()).unwrap();
        writer.write_all(&9u64.to_le_bytes()).unwrap();
        writer.write_all(&1u32.to_le_bytes()).unwrap();
        writer
            .write_all(&(membership.len() as u64).to_le_bytes())
            .unwrap();
        writer.write_all(&abundance_len.to_le_bytes()).unwrap();
        writer.write_all(&membership).unwrap();
        writer.write_all(&abundance).unwrap();
        writer.write_all(&0u32.to_le_bytes()).unwrap();
        writer.flush().unwrap();

        let sidecar = CidDatasetSidecar::open(&path).unwrap();
        let mut groups = Vec::new();
        sidecar
            .visit_groups_with_abundance(|_, ids, codes| {
                groups.push((ids.to_vec(), codes.unwrap().to_vec()));
                Ok(())
            })
            .unwrap();
        assert_eq!(groups, vec![(vec![1, 5, 9], vec![7, 19, 31])]);
    }

    #[test]
    fn abundance_headers_and_log_codes_are_supported() {
        assert_eq!(parse_header_abundance(b">42 LN:i:100 ka:f:17"), Some(17));
        assert_eq!(parse_header_abundance(b">42 km:i:23"), Some(23));
        for abundance in [1, 3, 17, 1_000, 50_000] {
            let decoded = decode_abundance(encode_abundance(abundance, 1.05), 1.05);
            let relative_error = decoded.abs_diff(abundance) as f64 / abundance as f64;
            assert!(relative_error <= 0.026, "{abundance} decoded as {decoded}");
        }
    }

    #[test]
    fn virtual_abundance_subsets_decode_without_generic_transposition() {
        let (dataset_ids, abundance_codes) = decode_virtual_abundance_sources(
            vec![2 * ABUNDANCE_CODES_PER_DATASET + 31, 9, 9],
            3 * ABUNDANCE_CODES_PER_DATASET,
        )
        .unwrap();
        assert_eq!(dataset_ids.as_slice(), &[1, 3]);
        assert_eq!(abundance_codes.as_slice(), &[9, 31]);

        assert!(decode_virtual_abundance_sources(
            vec![
                ABUNDANCE_CODES_PER_DATASET + 7,
                ABUNDANCE_CODES_PER_DATASET + 8,
            ],
            3 * ABUNDANCE_CODES_PER_DATASET,
        )
        .is_err());
    }

    #[test]
    fn structured_color_runs_are_split_without_fasta_headers() {
        let sequence = b"AACCGG";
        let mut colors = Vec::new();
        write_varint_u64(2, &mut colors);
        write_varint_u64(5, &mut colors);
        write_varint_u64(2, &mut colors);
        write_varint_u64(7, &mut colors);
        write_varint_u64(2, &mut colors);

        let mut block = Vec::new();
        block.extend_from_slice(&11u64.to_le_bytes());
        block.extend_from_slice(&(sequence.len() as u64).to_le_bytes());
        block.extend_from_slice(&(colors.len() as u64).to_le_bytes());
        block.extend_from_slice(&0u64.to_le_bytes());
        block.extend_from_slice(sequence);
        block.extend_from_slice(&colors);

        let mut records = Vec::new();
        let counts = visit_structured_ggcat_block(&block, 3, |record| {
            records.push(record);
            Ok(())
        })
        .unwrap();

        assert_eq!(counts, (1, 2));
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].subset, 5);
        assert_eq!(records[0].seq, b"AACC");
        assert_eq!(records[1].subset, 7);
        assert_eq!(records[1].seq, b"CCGG");
    }

    #[test]
    fn parallel_checkpoint_pipeline_preserves_assigned_color_ids() {
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("colors.dat");
        let serializer = Arc::new(
            ColorsSerializer::<RunLengthColorsSerializer>::new(&path, &[], Some(128), 4, false)
                .unwrap(),
        );

        let handles = (0..8)
            .map(|thread_id| {
                let serializer = Arc::clone(&serializer);
                thread::spawn(move || {
                    let mut assigned = Vec::new();
                    for index in 0..4_000u32 {
                        let color = (thread_id * 17 + index) % 128;
                        let subset_id = serializer.serialize_colors(&[color]);
                        assigned.push((subset_id, color));
                    }
                    assigned
                })
            })
            .collect::<Vec<_>>();

        let mut assigned = handles
            .into_iter()
            .flat_map(|handle| handle.join().unwrap())
            .collect::<Vec<_>>();
        assigned.sort_unstable_by_key(|(subset_id, _)| *subset_id);
        Arc::try_unwrap(serializer).ok().unwrap().finalize();

        let mut deserializer =
            ColorsDeserializer::<RunLengthColorsSerializer>::new(&path, false).unwrap();
        assert_eq!(deserializer.colors_count(), 128);
        let mut decoded = Vec::new();
        for (subset_id, expected_color) in assigned {
            decoded.clear();
            deserializer.get_color_mappings(subset_id, &mut decoded);
            assert_eq!(decoded, [expected_color]);
        }
    }
}
