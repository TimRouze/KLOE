use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Read, Seek};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use ggcat_api::{
    ColorIndexType, DnaSequence, DnaSequencesFileType, DynamicSequencesStream,
    GeneralSequenceBlockData, SequenceInfo,
};
use zstd::Decoder;

use crate::compress;

const BUCKET_SIZES_MAGIC: &[u8; 4] = b"KSB2";
const POSITIONS_MAGIC: &[u8; 4] = b"KPS2";
const ID_TO_CID_MAGIC: &[u8; 4] = b"KIC2";

const REQUIRED_ARCHIVE_FILES: [&str; 5] = [
    "filenames_id.txt",
    "positions_kloe.bin",
    "bucket_sizes.txt",
    "id_to_color_id.txt.zst",
    "tigs_kloe.fa",
];

const DECODE_BASES: [[u8; 4]; 256] = {
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

fn decode_packed_sequence(encoded: &[u8], size: usize, sequence: &mut Vec<u8>) {
    sequence.resize(size, b'A');
    for (&packed, output) in encoded.iter().zip(sequence.chunks_mut(4)) {
        let output_len = output.len();
        output.copy_from_slice(&DECODE_BASES[packed as usize][..output_len]);
    }
}

#[derive(Clone, Debug)]
pub struct MergeConfig {
    pub threads: usize,
    pub minimizer_size: usize,
    pub partition_power: u32,
    pub temp_dir: String,
    pub verify_kmers: bool,
    pub skip_sort: bool,
    pub use_unitigs: bool,
    pub use_matchtigs: bool,
    pub use_eulertigs: bool,
}

impl Default for MergeConfig {
    fn default() -> Self {
        Self {
            threads: 1,
            minimizer_size: 12,
            partition_power: 10,
            temp_dir: String::new(),
            verify_kmers: false,
            skip_sort: false,
            use_unitigs: false,
            use_matchtigs: false,
            use_eulertigs: false,
        }
    }
}

#[derive(Debug)]
struct ArchiveInfo {
    filenames: Vec<String>,
    positions: Arc<ArchivePositions>,
    cid_to_ids: CidDatasetIds,
    tigs_path: PathBuf,
    sizes_path: PathBuf,
    compact_sizes: Option<Arc<CompactSizesIndex>>,
}

#[derive(Debug)]
struct ArchivePositions {
    tigs: Vec<u64>,
    legacy_sizes: Option<Vec<u64>>,
}

impl ArchivePositions {
    fn len(&self) -> usize {
        self.tigs.len()
    }

    fn is_empty(&self) -> bool {
        self.tigs.is_empty()
    }

    fn tigs(&self, cid: usize) -> Option<u64> {
        self.tigs.get(cid).copied()
    }

    fn sizes(&self, cid: usize) -> Option<u64> {
        match &self.legacy_sizes {
            Some(sizes) => sizes.get(cid).copied(),
            None => u64::try_from(cid).ok(),
        }
    }
}

#[derive(Debug)]
struct CidDatasetIds {
    offsets: Vec<usize>,
    values: Arc<Vec<u32>>,
}

impl CidDatasetIds {
    fn span(&self, cid: usize) -> Option<(usize, usize)> {
        Some((*self.offsets.get(cid)?, *self.offsets.get(cid + 1)?))
    }

    fn active_cids(&self) -> impl Iterator<Item = usize> + '_ {
        self.offsets
            .windows(2)
            .enumerate()
            .filter_map(|(cid, range)| (range[0] != range[1]).then_some(cid))
    }
}

#[derive(Debug)]
struct CompactSizesBlock {
    first_group: usize,
    group_count: usize,
    offsets: Vec<u32>,
    data_offset: u64,
    data_len: usize,
}

#[derive(Debug)]
struct CompactSizesIndex {
    path: PathBuf,
    blocks: Vec<CompactSizesBlock>,
}

#[derive(Clone)]
struct ArchiveSequencesStream {
    positions: Arc<ArchivePositions>,
    tigs_path: PathBuf,
    sizes_path: PathBuf,
    compact_sizes: Option<Arc<CompactSizesIndex>>,
    blocks: Vec<ArchiveSequenceBlock>,
}

#[derive(Clone)]
struct ArchiveSequenceBlock {
    cids: Vec<(usize, ColorIndexType)>,
    estimated_bases: u64,
}

impl ArchiveSequencesStream {
    fn read_block_inner(
        &self,
        block: usize,
        callback: &mut dyn FnMut(DnaSequence, SequenceInfo),
    ) -> io::Result<()> {
        let block_data = self.blocks.get(block).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("archive input block {} is outside the configured block range", block),
            )
        })?;
        let mut tigs_reader = BufReader::with_capacity(1024 * 1024, File::open(&self.tigs_path)?);
        let mut legacy_sizes_reader = if self.compact_sizes.is_none() {
            Some(BufReader::with_capacity(
                1024 * 1024,
                File::open(&self.sizes_path)?,
            ))
        } else {
            None
        };
        let mut compact_sizes_reader = if self.compact_sizes.is_some() {
            Some(File::open(&self.sizes_path)?)
        } else {
            None
        };
        let mut cached_block_index = usize::MAX;
        let mut cached_block_data = Vec::new();
        let mut encoded = Vec::new();
        let mut sequence = Vec::new();
        let mut sizes = Vec::new();
        let mut current_tigs_pos = None;

        for &(cid, color) in &block_data.cids {
            if cid >= self.positions.len().saturating_sub(1) {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("archive CID {} is outside the sequence range", cid),
                ));
            }
            let tigs_pos = self.positions.tigs(cid).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "archive CID has no tig position")
            })?;
            let sizes_pos = self.positions.sizes(cid).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "archive CID has no size position")
            })?;
            if let Some(compact) = self.compact_sizes.as_deref() {
                let group = usize::try_from(sizes_pos).map_err(|_| {
                    io::Error::new(io::ErrorKind::InvalidData, "size group cannot fit usize")
                })?;
                let block_index = compact.block_index_for_group(group)?;
                if cached_block_index != block_index {
                    cached_block_data = compact.load_block_from(
                        block_index,
                        compact_sizes_reader
                            .as_mut()
                            .expect("compact sizes reader must be initialized"),
                    )?;
                    cached_block_index = block_index;
                }
                compact.decode_group_into(block_index, group, &cached_block_data, &mut sizes)?;
            } else {
                read_bucket_sizes_at_into(
                    legacy_sizes_reader
                        .as_mut()
                        .expect("legacy sizes reader must be initialized"),
                    sizes_pos,
                    &mut sizes,
                )?;
            }

            if current_tigs_pos != Some(tigs_pos) {
                tigs_reader.seek(std::io::SeekFrom::Start(tigs_pos))?;
                current_tigs_pos = Some(tigs_pos);
            }
            for &size in &sizes {
                if size == 0 {
                    continue;
                }
                encoded.resize(size.div_ceil(4), 0);
                tigs_reader.read_exact(&mut encoded)?;
                current_tigs_pos = Some(
                    current_tigs_pos
                        .unwrap_or(tigs_pos)
                        .saturating_add(encoded.len() as u64),
                );
                decode_packed_sequence(&encoded, size, &mut sequence);
                callback(
                    DnaSequence {
                        ident_data: &[],
                        seq: &sequence,
                        format: DnaSequencesFileType::FASTA,
                    },
                    SequenceInfo { color: Some(color) },
                );
            }
        }
        Ok(())
    }

    fn estimated_block_bases(&self, block: usize) -> io::Result<u64> {
        let block = self.blocks.get(block).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "invalid archive input block")
        })?;
        Ok(block.estimated_bases.max(1))
    }
}

impl DynamicSequencesStream for ArchiveSequencesStream {
    fn read_block(
        &self,
        block: usize,
        _copy_ident_data: bool,
        _partial_read_copyback: Option<usize>,
        callback: &mut dyn FnMut(DnaSequence, SequenceInfo),
    ) {
        self.read_block_inner(block, callback)
            .unwrap_or_else(|err| panic!("cannot stream archive block {block} into ggcat: {err}"));
    }

    fn estimated_base_count(&self, block: usize) -> u64 {
        self.estimated_block_bases(block).unwrap_or(1).max(1)
    }
}
impl ArchiveInfo {
    fn load(root: &Path, _k: usize) -> io::Result<Self> {
        ensure_archive_complete(root)?;

        let positions_path = root.join("positions_kloe.bin");
        let sizes_path = root.join("bucket_sizes.txt");
        let cid_path = root.join("id_to_color_id.txt.zst");
        let tigs_path = root.join("tigs_kloe.fa");
        let filenames_path = root.join("filenames_id.txt");

        let compact_sizes = if is_compact_bucket_sizes_file(&sizes_path)? {
            Some(Arc::new(load_compact_bucket_sizes_index(&sizes_path)?))
        } else {
            None
        };
        let positions = Arc::new(load_positions(&positions_path, compact_sizes.is_some())?);
        if positions.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "positions file '{}' is empty",
                    positions_path.to_string_lossy()
                ),
            ));
        }

        let filenames = load_filenames(&filenames_path)?;
        let cid_to_ids = load_cid_to_ids(
            &cid_path,
            filenames.len(),
            positions.len().saturating_sub(1),
        )?;

        Ok(Self {
            filenames,
            positions,
            cid_to_ids,
            tigs_path,
            sizes_path,
            compact_sizes,
        })
    }

    fn sequence_stream(
        &self,
        colored_cids: Vec<(usize, ColorIndexType)>,
        target_blocks: usize,
    ) -> ArchiveSequencesStream {
        let total_bytes = colored_cids.iter().fold(0u64, |total, &(cid, _)| {
            let start = self.positions.tigs(cid).unwrap_or(0);
            let end = self.positions.tigs(cid + 1).unwrap_or(start);
            total.saturating_add(end.saturating_sub(start))
        });
        let target_bytes = total_bytes
            .div_ceil(target_blocks.max(1) as u64)
            .max(1);
        let mut blocks = Vec::with_capacity(target_blocks.min(colored_cids.len()).max(1));
        let mut current_cids = Vec::new();
        let mut current_bytes = 0u64;
        for colored_cid in colored_cids {
            let start = self.positions.tigs(colored_cid.0).unwrap_or(0);
            let end = self.positions.tigs(colored_cid.0 + 1).unwrap_or(start);
            let cid_bytes = end.saturating_sub(start);
            if !current_cids.is_empty()
                && current_bytes >= target_bytes
                && blocks.len() + 1 < target_blocks.max(1)
            {
                blocks.push(ArchiveSequenceBlock {
                    cids: std::mem::take(&mut current_cids),
                    estimated_bases: current_bytes.saturating_mul(4).max(1),
                });
                current_bytes = 0;
            }
            current_cids.push(colored_cid);
            current_bytes = current_bytes.saturating_add(cid_bytes);
        }
        if !current_cids.is_empty() {
            blocks.push(ArchiveSequenceBlock {
                cids: current_cids,
                estimated_bases: current_bytes.saturating_mul(4).max(1),
            });
        }
        ArchiveSequencesStream {
            positions: Arc::clone(&self.positions),
            tigs_path: self.tigs_path.clone(),
            sizes_path: self.sizes_path.clone(),
            compact_sizes: self.compact_sizes.clone(),
            blocks,
        }
    }

}

pub fn merge_archives(
    archive_a_dir: &str,
    archive_b_dir: &str,
    output_dir: &str,
    k: usize,
    memory_gb: usize,
) -> io::Result<()> {
    merge_archives_with_config(
        archive_a_dir,
        archive_b_dir,
        output_dir,
        k,
        memory_gb,
        MergeConfig::default(),
    )
}

pub fn merge_archives_with_config(
    archive_a_dir: &str,
    archive_b_dir: &str,
    output_dir: &str,
    k: usize,
    memory_gb: usize,
    cfg: MergeConfig,
) -> io::Result<()> {
    if k == 0 || k > 32 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("merge supports k in [1, 32], got {k}"),
        ));
    }
    if output_dir.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "missing output directory for merge",
        ));
    }

    let archive_a_root = Path::new(archive_a_dir);
    let archive_b_root = Path::new(archive_b_dir);
    let output_root = Path::new(output_dir);

    fs::create_dir_all(output_root)?;
    let output_dir_norm = normalize_output_dir(output_root);

    println!(
        "Loading archive A from {}",
        archive_a_root.to_string_lossy()
    );
    let archive_a = ArchiveInfo::load(archive_a_root, k)?;
    println!(
        "Loading archive B from {}",
        archive_b_root.to_string_lossy()
    );
    let archive_b = ArchiveInfo::load(archive_b_root, k)?;

    let mut merged_filenames =
        Vec::with_capacity(archive_a.filenames.len() + archive_b.filenames.len());
    merged_filenames.extend(archive_a.filenames.iter().cloned());
    merged_filenames.extend(archive_b.filenames.iter().cloned());

    let offset = u32::try_from(archive_a.filenames.len()).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "archive A has more than u32::MAX datasets",
        )
    })?;
    u32::try_from(merged_filenames.len()).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "merged archive has more than u32::MAX datasets",
        )
    })?;

    let source_count = archive_a.cid_to_ids.active_cids().count()
        + archive_b.cid_to_ids.active_cids().count();
    let mut source_dataset_ids = compress::SourceDatasetMap::new();
    let a_storage = source_dataset_ids.add_storage(Arc::clone(&archive_a.cid_to_ids.values));
    let b_storage = source_dataset_ids.add_storage(Arc::clone(&archive_b.cid_to_ids.values));
    let mut a_colored_cids = Vec::with_capacity(archive_a.cid_to_ids.active_cids().count());
    let mut b_colored_cids = Vec::with_capacity(archive_b.cid_to_ids.active_cids().count());

    for cid in archive_a.cid_to_ids.active_cids() {
        let color = ColorIndexType::try_from(source_dataset_ids.len()).map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, "too many merge source colors")
        })?;
        a_colored_cids.push((cid, color));
        let (start, end) = archive_a
            .cid_to_ids
            .span(cid)
            .expect("active archive CID must have a dataset span");
        source_dataset_ids.push_span(a_storage, start, end, 0)?;
    }
    for cid in archive_b.cid_to_ids.active_cids() {
        let color = ColorIndexType::try_from(source_dataset_ids.len()).map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, "too many merge source colors")
        })?;
        b_colored_cids.push((cid, color));
        let (start, end) = archive_b
            .cid_to_ids
            .span(cid)
            .expect("active archive CID must have a dataset span");
        source_dataset_ids.push_span(b_storage, start, end, offset)?;
    }
    debug_assert_eq!(source_dataset_ids.len(), source_count);

    let blocks_per_archive = cfg.threads.max(1).div_ceil(2);
    let a_stream = Arc::new(archive_a.sequence_stream(a_colored_cids, blocks_per_archive));
    let b_stream = Arc::new(archive_b.sequence_stream(b_colored_cids, blocks_per_archive));
    let a_block_count = a_stream.blocks.len();
    let b_block_count = b_stream.blocks.len();
    let a_stream: Arc<dyn DynamicSequencesStream> = a_stream;
    let b_stream: Arc<dyn DynamicSequencesStream> = b_stream;
    drop(archive_a);
    drop(archive_b);
    let mut input_streams = Vec::with_capacity(a_block_count + b_block_count);
    for block in 0..a_block_count {
        input_streams.push(GeneralSequenceBlockData::Dynamic((
            Arc::clone(&a_stream),
            block,
        )));
    }
    for block in 0..b_block_count {
        input_streams.push(GeneralSequenceBlockData::Dynamic((
            Arc::clone(&b_stream),
            block,
        )));
    }

    println!(
        "Starting direct merge with {} archive color sets in {} input blocks and {} output datasets; packed SPSS input is streamed directly into ggcat",
        source_dataset_ids.len(),
        input_streams.len(),
        merged_filenames.len()
    );
    if cfg.verify_kmers {
        eprintln!("Warning: --verify-kmers is ignored by the embedded ggcat merge backend.");
    }
    if cfg.skip_sort {
        eprintln!("Warning: --skip-sort is ignored; merged GGCAT records are grouped by colors.");
    }

    compress::compress_ggcat_sources(
        &output_dir_norm,
        merged_filenames,
        input_streams,
        source_dataset_ids,
        cfg.threads.max(1),
        k,
        cfg.minimizer_size,
        cfg.use_unitigs,
        cfg.use_matchtigs,
        cfg.use_eulertigs,
        compress::GgcatCompressionConfig {
            memory_gb,
            temp_dir: cfg.temp_dir.clone(),
        },
    )?;

    println!(
        "Merge complete: output archive written to {} ({})",
        output_root.to_string_lossy(),
        merge_mode_name(&cfg)
    );
    Ok(())
}

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

fn decode_varint_delta_cids(payload: &[u8], context: &str) -> io::Result<Vec<usize>> {
    let mut cursor = std::io::Cursor::new(payload);
    let mut out = Vec::new();
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
        out.push(current);
    }
    Ok(out)
}

fn decode_varint_deltas_to_sizes_into(
    payload: &[u8],
    context: &str,
    sizes: &mut Vec<usize>,
) -> io::Result<()> {
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

fn is_compact_bucket_sizes_file(path: &Path) -> io::Result<bool> {
    let mut reader = BufReader::new(File::open(path)?);
    let mut magic = [0u8; 4];
    match reader.read_exact(&mut magic) {
        Ok(()) => Ok(&magic == BUCKET_SIZES_MAGIC),
        Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => Ok(false),
        Err(err) => Err(err),
    }
}

impl CompactSizesIndex {
    fn block_index_for_group(&self, group: usize) -> io::Result<usize> {
        let index = self
            .blocks
            .partition_point(|block| block.first_group + block.group_count <= group);
        match self.blocks.get(index) {
            Some(block) if group >= block.first_group => Ok(index),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("compact size group {} is not indexed", group),
            )),
        }
    }

    fn load_block_from(&self, block_index: usize, file: &mut File) -> io::Result<Vec<u8>> {
        let block = self.blocks.get(block_index).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "compact size block is not indexed")
        })?;
        file.seek(std::io::SeekFrom::Start(block.data_offset))?;
        let mut compressed = vec![0u8; block.data_len];
        file.read_exact(&mut compressed)?;
        let mut decompressed = Vec::new();
        Decoder::new(&compressed[..])?.read_to_end(&mut decompressed)?;
        let expected_min = block.offsets.last().copied().unwrap_or(0) as usize;
        if decompressed.len() < expected_min {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("truncated compact size block in '{}'", self.path.display()),
            ));
        }
        Ok(decompressed)
    }

    fn decode_group_into(
        &self,
        block_index: usize,
        group: usize,
        decompressed: &[u8],
        sizes: &mut Vec<usize>,
    ) -> io::Result<()> {
        let block = self.blocks.get(block_index).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "compact size block is not indexed")
        })?;
        let local = group.checked_sub(block.first_group).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "compact size group precedes block")
        })?;
        if local >= block.group_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "compact size group exceeds block",
            ));
        }
        let start = block.offsets[local] as usize;
        let end = block.offsets[local + 1] as usize;
        if end < start || end > decompressed.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid compact size offsets in '{}'", self.path.display()),
            ));
        }
        decode_varint_deltas_to_sizes_into(
            &decompressed[start..end],
            &format!("{}#group{}", self.path.display(), group),
            sizes,
        )
    }

}

fn load_compact_bucket_sizes_index(path: &Path) -> io::Result<CompactSizesIndex> {
    let mut reader = BufReader::new(File::open(path)?);
    let mut magic = [0u8; 4];
    reader.read_exact(&mut magic)?;
    if &magic != BUCKET_SIZES_MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("expected compact bucket sizes magic in '{}'", path.display()),
        ));
    }

    let mut blocks = Vec::new();
    let mut first_group = 0usize;
    loop {
        let mut groups_buf = [0u8; 4];
        match reader.read_exact(&mut groups_buf) {
            Ok(()) => {}
            Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => break,
            Err(err) => return Err(err),
        }
        let group_count = u32::from_le_bytes(groups_buf) as usize;
        if group_count == 0 {
            break;
        }

        let mut compressed_len_buf = [0u8; 8];
        reader.read_exact(&mut compressed_len_buf)?;
        let data_len = usize::try_from(u64::from_le_bytes(compressed_len_buf)).map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, "compact size block is too large")
        })?;
        let mut offsets = vec![0u32; group_count + 1];
        for offset in &mut offsets {
            let mut bytes = [0u8; 4];
            reader.read_exact(&mut bytes)?;
            *offset = u32::from_le_bytes(bytes);
        }
        let data_offset = reader.stream_position()?;
        let next_offset = data_offset.checked_add(data_len as u64).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "compact size file offset overflow")
        })?;
        reader.seek(std::io::SeekFrom::Start(next_offset))?;
        blocks.push(CompactSizesBlock {
            first_group,
            group_count,
            offsets,
            data_offset,
            data_len,
        });
        first_group = first_group.saturating_add(group_count);
    }

    Ok(CompactSizesIndex {
        path: path.to_path_buf(),
        blocks,
    })
}
fn load_positions(path: &Path, compact_sizes: bool) -> io::Result<ArchivePositions> {
    let mut file = BufReader::new(File::open(path)?);
    let mut magic = [0u8; 4];
    match file.read_exact(&mut magic) {
        Ok(()) => {}
        Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => {
            return Ok(ArchivePositions {
                tigs: Vec::new(),
                legacy_sizes: (!compact_sizes).then(Vec::new),
            })
        }
        Err(err) => return Err(err),
    }

    if &magic == POSITIONS_MAGIC {
        let count_u64 = read_varint_u64_from_reader(&mut file)?.ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!(
                    "missing compact positions count in '{}'",
                    path.to_string_lossy()
                ),
            )
        })?;
        let count = usize::try_from(count_u64).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "positions count {} cannot fit usize in '{}'",
                    count_u64,
                    path.to_string_lossy()
                ),
            )
        })?;
        let mut tigs = Vec::with_capacity(count);
        let mut legacy_sizes = (!compact_sizes).then(|| Vec::with_capacity(count));
        let mut tigs_pos = 0u64;
        let mut sizes_pos = 0u64;
        for _ in 0..count {
            let dt = read_varint_u64_from_reader(&mut file)?.ok_or_else(|| {
                io::Error::new(io::ErrorKind::UnexpectedEof, "truncated tigs delta")
            })?;
            let ds = read_varint_u64_from_reader(&mut file)?.ok_or_else(|| {
                io::Error::new(io::ErrorKind::UnexpectedEof, "truncated sizes delta")
            })?;
            tigs_pos = tigs_pos.checked_add(dt).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "tigs position overflow")
            })?;
            sizes_pos = sizes_pos.checked_add(ds).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "sizes position overflow")
            })?;
            tigs.push(tigs_pos);
            if let Some(positions) = &mut legacy_sizes {
                positions.push(sizes_pos);
            }
        }
        return Ok(ArchivePositions { tigs, legacy_sizes });
    }

    file.seek(std::io::SeekFrom::Start(0))?;
    let file_size = file.get_ref().metadata()?.len() as usize;
    if file_size % 16 != 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "positions file '{}' has size {} not divisible by 16",
                path.to_string_lossy(),
                file_size
            ),
        ));
    }
    let num_entries = file_size / 16;
    let mut tigs = Vec::with_capacity(num_entries);
    let mut legacy_sizes = (!compact_sizes).then(|| Vec::with_capacity(num_entries));
    let mut buf = [0u8; 16];
    for _ in 0..num_entries {
        file.read_exact(&mut buf)?;
        let tigs_pos = u64::from_le_bytes(buf[..8].try_into().unwrap());
        let sizes_pos = u64::from_le_bytes(buf[8..16].try_into().unwrap());
        tigs.push(tigs_pos);
        if let Some(positions) = &mut legacy_sizes {
            positions.push(sizes_pos);
        }
    }
    Ok(ArchivePositions { tigs, legacy_sizes })
}

fn load_filenames(path: &Path) -> io::Result<Vec<String>> {
    let reader = BufReader::new(File::open(path)?);
    let mut names = Vec::new();
    for line_result in reader.lines() {
        let line = line_result?;
        if line.trim().is_empty() {
            continue;
        }
        let Some((name, _)) = line.split_once(':') else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "invalid filenames_id row '{}' in '{}'",
                    line,
                    path.to_string_lossy()
                ),
            ));
        };
        names.push(name.to_string());
    }
    Ok(names)
}

fn for_each_dataset_cids(
    path: &Path,
    expected_dataset_count: usize,
    mut visit: impl FnMut(usize, &[usize]) -> io::Result<()>,
) -> io::Result<()> {
    let mut reader = BufReader::new(File::open(path)?);
    let mut dataset_idx = 0usize;
    let mut magic = [0u8; 4];
    let binary_varints = match reader.read_exact(&mut magic) {
        Ok(()) if &magic == ID_TO_CID_MAGIC => true,
        Ok(()) => {
            reader.seek(std::io::SeekFrom::Start(0))?;
            false
        }
        Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => false,
        Err(err) => return Err(err),
    };

    loop {
        let mut len_buf = [0u8; 8];
        match reader.read_exact(&mut len_buf) {
            Ok(()) => {}
            Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => break,
            Err(err) => return Err(err),
        }
        let payload_len = usize::try_from(u64::from_le_bytes(len_buf)).map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, "dataset CID payload is too large")
        })?;
        if payload_len == 0 {
            break;
        }
        if dataset_idx >= expected_dataset_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "id_to_color_id contains more than {} datasets in '{}'",
                    expected_dataset_count,
                    path.to_string_lossy()
                ),
            ));
        }

        let mut payload = vec![0u8; payload_len];
        reader.read_exact(&mut payload)?;
        let mut decompressed = Vec::new();
        Decoder::new(&payload[..])?.read_to_end(&mut decompressed)?;
        let context = format!("{}#{}", path.to_string_lossy(), dataset_idx);
        let mut cids = if binary_varints {
            decode_varint_delta_cids(&decompressed, &context)?
        } else {
            let text = String::from_utf8(decompressed).map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid UTF-8 while decoding '{}': {}", context, err),
                )
            })?;
            let mut out = Vec::new();
            let mut current_cid = 0usize;
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
                current_cid = current_cid.checked_add(delta).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("cid delta overflow while decoding '{}'", context),
                    )
                })?;
                out.push(current_cid);
            }
            out
        };
        cids.dedup();
        visit(dataset_idx, &cids)?;
        dataset_idx += 1;
    }

    if dataset_idx != expected_dataset_count {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "filenames count ({}) differs from id_to_color_id count ({}) in '{}'",
                expected_dataset_count,
                dataset_idx,
                path.to_string_lossy()
            ),
        ));
    }
    Ok(())
}

fn load_cid_to_ids(
    path: &Path,
    expected_dataset_count: usize,
    cid_count: usize,
) -> io::Result<CidDatasetIds> {
    let mut counts = vec![0usize; cid_count];
    for_each_dataset_cids(path, expected_dataset_count, |_, cids| {
        for &cid in cids {
            let count = counts.get_mut(cid).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "cid index {} outside [0, {}) in '{}'",
                        cid,
                        cid_count,
                        path.to_string_lossy()
                    ),
                )
            })?;
            *count = count.checked_add(1).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "CID membership count overflow")
            })?;
        }
        Ok(())
    })?;

    let mut offsets = Vec::with_capacity(cid_count + 1);
    offsets.push(0usize);
    for &count in &counts {
        let next = offsets.last().copied().unwrap().checked_add(count).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "CID membership storage overflow")
        })?;
        offsets.push(next);
    }
    let mut values = vec![0u32; offsets.last().copied().unwrap_or(0)];
    counts.fill(0);
    for_each_dataset_cids(path, expected_dataset_count, |dataset_idx, cids| {
        let dataset_id = u32::try_from(dataset_idx + 1).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "merge supports at most u32::MAX datasets",
            )
        })?;
        for &cid in cids {
            let written = counts.get_mut(cid).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "CID disappeared between passes")
            })?;
            let index = offsets[cid].checked_add(*written).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "CID membership offset overflow")
            })?;
            let slot = values.get_mut(index).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "CID membership count changed")
            })?;
            *slot = dataset_id;
            *written += 1;
        }
        Ok(())
    })?;
    if counts
        .iter()
        .zip(offsets.windows(2))
        .any(|(&written, span)| written != span[1] - span[0])
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "CID membership counts changed while loading the archive",
        ));
    }

    Ok(CidDatasetIds {
        offsets,
        values: Arc::new(values),
    })
}

fn read_bucket_sizes_at_into<R: Read + Seek>(
    reader: &mut R,
    offset: u64,
    sizes: &mut Vec<usize>,
) -> io::Result<()> {
    reader.seek(std::io::SeekFrom::Start(offset))?;
    sizes.clear();
    let mut len_buf = [0u8; 8];
    match reader.read_exact(&mut len_buf) {
        Ok(()) => {}
        Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => return Ok(()),
        Err(err) => return Err(err),
    }
    let payload_len = u64::from_le_bytes(len_buf) as usize;
    if payload_len == 0 {
        return Ok(());
    }

    let mut payload = vec![0u8; payload_len];
    reader.read_exact(&mut payload)?;
    let mut decompressed = Vec::new();
    Decoder::new(&payload[..])?.read_to_end(&mut decompressed)?;
    if decompressed.len() % 8 != 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid decompressed bucket size payload ({} bytes)",
                decompressed.len()
            ),
        ));
    }

    let mut previous = 0usize;
    sizes.reserve(decompressed.len() / 8);
    for chunk in decompressed.chunks_exact(8) {
        let delta = u64::from_le_bytes(chunk.try_into().unwrap()) as usize;
        let size = previous.saturating_add(delta);
        sizes.push(size);
        previous = size;
    }
    Ok(())
}

fn ensure_archive_complete(root: &Path) -> io::Result<()> {
    for filename in REQUIRED_ARCHIVE_FILES {
        let path = root.join(filename);
        if !path.exists() {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!("missing archive file '{}'", path.to_string_lossy()),
            ));
        }
    }
    Ok(())
}

fn normalize_output_dir(root: &Path) -> String {
    let mut s = root.to_string_lossy().to_string();
    if !s.ends_with('/') {
        s.push('/');
    }
    s
}

fn merge_mode_name(cfg: &MergeConfig) -> &'static str {
    if cfg.use_unitigs {
        "unitigs"
    } else if cfg.use_matchtigs {
        "matchtigs"
    } else if cfg.use_eulertigs {
        "eulertigs"
    } else {
        "simplitigs"
    }
}
