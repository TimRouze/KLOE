use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use ggcat_api::{
    ColorIndexType, DnaSequence, DnaSequencesFileType, DynamicSequencesStream,
    GeneralSequenceBlockData, SequenceInfo,
};
use zstd::Decoder;

use crate::compress;
use crate::packed_tigs::PackedTigsReader;

const BUCKET_SIZES_MAGIC: &[u8; 4] = b"KSB2";
const BUCKET_SIZES_COMPACT_MAGIC: &[u8; 4] = b"KSB3";
const POSITIONS_MAGIC: &[u8; 4] = b"KPS2";
const POSITIONS_COMPACT_MAGIC: &[u8; 4] = b"KPS3";
const ID_TO_CID_MAGIC: &[u8; 4] = b"KIC2";

const REQUIRED_ARCHIVE_FILES: [&str; 4] = [
    "filenames_id.txt",
    "positions_kloe.bin",
    "bucket_sizes.txt",
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
    pub structural: bool,
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
            structural: false,
        }
    }
}

#[derive(Debug)]
struct ArchiveInfo {
    filenames: Vec<String>,
    positions: Arc<ArchivePositions>,
    cid_to_ids: ArchiveDatasetIds,
    tigs_path: PathBuf,
    sizes_path: PathBuf,
    compact_sizes: Option<Arc<CompactSizesIndex>>,
}

#[derive(Debug)]
struct CompactArchiveLayout {
    root: PathBuf,
    filenames: Vec<String>,
    position_entries: u64,
    group_count: u64,
    packed_tigs_bytes: u64,
}

fn require_magic(path: &Path, expected: &[u8; 4]) -> io::Result<BufReader<File>> {
    let mut reader = BufReader::new(File::open(path)?);
    let mut magic = [0u8; 4];
    reader.read_exact(&mut magic)?;
    if &magic != expected {
        return Err(io::Error::new(
            io::ErrorKind::Unsupported,
            format!(
                "'{}' uses a legacy archive encoding and cannot be structurally joined",
                path.display()
            ),
        ));
    }
    Ok(reader)
}

fn inspect_compact_positions(path: &Path) -> io::Result<(u64, u64, u64)> {
    let mut reader = BufReader::new(File::open(path)?);
    let mut magic = [0u8; 4];
    reader.read_exact(&mut magic)?;
    if &magic != POSITIONS_MAGIC && &magic != POSITIONS_COMPACT_MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::Unsupported,
            format!("'{}' uses an unsupported positions encoding", path.display()),
        ));
    }
    let implicit_sizes = &magic == POSITIONS_COMPACT_MAGIC;
    let entries = read_varint_u64_from_reader(&mut reader)?.ok_or_else(|| {
        io::Error::new(io::ErrorKind::UnexpectedEof, "missing compact positions count")
    })?;
    if entries == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("compact positions file '{}' has no entries", path.display()),
        ));
    }
    let mut tigs_position = 0u64;
    let mut sizes_position = 0u64;
    for entry in 0..entries {
        let tigs_delta = read_varint_u64_from_reader(&mut reader)?.ok_or_else(|| {
            io::Error::new(io::ErrorKind::UnexpectedEof, "truncated compact tig positions")
        })?;
        let sizes_delta = if implicit_sizes {
            u64::from(entry > 0)
        } else {
            read_varint_u64_from_reader(&mut reader)?.ok_or_else(|| {
                io::Error::new(io::ErrorKind::UnexpectedEof, "truncated compact size positions")
            })?
        };
        tigs_position = tigs_position.checked_add(tigs_delta).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "compact tig position overflow")
        })?;
        sizes_position = sizes_position.checked_add(sizes_delta).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "compact size position overflow")
        })?;
    }
    if read_varint_u64_from_reader(&mut reader)?.is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("compact positions file '{}' has trailing entries", path.display()),
        ));
    }
    let groups = entries - 1;
    if sizes_position != groups {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "compact positions in '{}' describe {} size groups, expected {}",
                path.display(),
                sizes_position,
                groups
            ),
        ));
    }
    Ok((entries, groups, tigs_position))
}

fn inspect_compact_bucket_sizes(path: &Path) -> io::Result<u64> {
    let mut reader = BufReader::new(File::open(path)?);
    let mut magic = [0u8; 4];
    reader.read_exact(&mut magic)?;
    if &magic != BUCKET_SIZES_MAGIC && &magic != BUCKET_SIZES_COMPACT_MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::Unsupported,
            format!("'{}' uses an unsupported bucket-size encoding", path.display()),
        ));
    }
    let length_prefixed = &magic == BUCKET_SIZES_COMPACT_MAGIC;
    let mut groups = 0u64;
    loop {
        let mut count_buf = [0u8; 4];
        match reader.read_exact(&mut count_buf) {
            Ok(()) => {}
            Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => break,
            Err(err) => return Err(err),
        }
        let count = u32::from_le_bytes(count_buf) as u64;
        if count == 0 {
            break;
        }
        let mut len_buf = [0u8; 8];
        reader.read_exact(&mut len_buf)?;
        let compressed_len = u64::from_le_bytes(len_buf);
        let offsets_bytes = if length_prefixed {
            0
        } else {
            (count + 1).checked_mul(4).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "bucket-size offset overflow")
            })?
        };
        let skip = offsets_bytes.checked_add(compressed_len).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "bucket-size block overflow")
        })?;
        reader.seek(SeekFrom::Current(i64::try_from(skip).map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, "bucket-size block is too large")
        })?))?;
        groups = groups.checked_add(count).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "bucket-size group count overflow")
        })?;
    }
    Ok(groups)
}

impl CompactArchiveLayout {
    fn inspect(root: &Path) -> io::Result<Self> {
        ensure_archive_complete(root)?;
        let positions_path = root.join("positions_kloe.bin");
        let sizes_path = root.join("bucket_sizes.txt");
        let tigs_path = root.join("tigs_kloe.fa");
        let legacy_index = root.join("id_to_color_id.txt.zst");
        if !legacy_index.is_file() {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                format!(
                    "compact archive '{}' omits the redundant dataset-to-CID index; use the default deduplicating merge",
                    root.display()
                ),
            ));
        }
        let _ = require_magic(&legacy_index, ID_TO_CID_MAGIC)?;
        let (position_entries, group_count, packed_tigs_bytes) =
            inspect_compact_positions(&positions_path)?;
        let sizes_groups = inspect_compact_bucket_sizes(&sizes_path)?;
        if sizes_groups != group_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "archive '{}' has {} position groups but {} bucket-size groups",
                    root.display(),
                    group_count,
                    sizes_groups
                ),
            ));
        }
        let actual_tigs_bytes = fs::metadata(&tigs_path)?.len();
        if actual_tigs_bytes != packed_tigs_bytes {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "archive '{}' positions end at byte {}, but packed tigs contain {} bytes",
                    root.display(),
                    packed_tigs_bytes,
                    actual_tigs_bytes
                ),
            ));
        }
        Ok(Self {
            root: root.to_path_buf(),
            filenames: load_filenames(&root.join("filenames_id.txt"))?,
            position_entries,
            group_count,
            packed_tigs_bytes,
        })
    }
}

fn write_varint_u64_to_writer(mut value: u64, writer: &mut impl Write) -> io::Result<()> {
    loop {
        let mut byte = (value & 0x7f) as u8;
        value >>= 7;
        if value != 0 {
            byte |= 0x80;
        }
        writer.write_all(&[byte])?;
        if value == 0 {
            return Ok(());
        }
    }
}

fn copy_file_into(path: &Path, output: &mut BufWriter<File>) -> io::Result<u64> {
    let mut input = BufReader::with_capacity(16 * 1024 * 1024, File::open(path)?);
    io::copy(&mut input, output)
}

fn append_compact_file_body(
    path: &Path,
    magic: &[u8; 4],
    output: &mut BufWriter<File>,
) -> io::Result<u64> {
    let mut input = require_magic(path, magic)?;
    io::copy(&mut input, output)
}

fn append_position_deltas(
    path: &Path,
    skip_initial_entry: bool,
    output: &mut BufWriter<File>,
) -> io::Result<()> {
    let mut input = require_magic(path, POSITIONS_MAGIC)?;
    let entries = read_varint_u64_from_reader(&mut input)?.ok_or_else(|| {
        io::Error::new(io::ErrorKind::UnexpectedEof, "missing compact positions count")
    })?;
    if entries == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "compact positions contain no initial entry",
        ));
    }
    if skip_initial_entry {
        let initial_tigs = read_varint_u64_from_reader(&mut input)?.ok_or_else(|| {
            io::Error::new(io::ErrorKind::UnexpectedEof, "missing initial tig position")
        })?;
        let initial_sizes = read_varint_u64_from_reader(&mut input)?.ok_or_else(|| {
            io::Error::new(io::ErrorKind::UnexpectedEof, "missing initial size position")
        })?;
        if initial_tigs != 0 || initial_sizes != 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "compact positions do not begin at zero",
            ));
        }
    }
    io::copy(&mut input, output)?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn append_shifted_dataset_cids(
    path: &Path,
    expected_datasets: usize,
    source_group_count: u64,
    cid_offset: u64,
    output: &mut BufWriter<File>,
    output_position: &mut u64,
    output_offsets: &mut Vec<usize>,
    scratch: &mut File,
) -> io::Result<()> {
    let mut input = require_magic(path, ID_TO_CID_MAGIC)?;
    for dataset in 0..expected_datasets {
        let mut len_buf = [0u8; 8];
        input.read_exact(&mut len_buf).map_err(|err| {
            io::Error::new(
                err.kind(),
                format!("missing dataset {} CID payload in '{}': {err}", dataset, path.display()),
            )
        })?;
        let payload_len = u64::from_le_bytes(len_buf);
        if payload_len == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("premature CID terminator at dataset {} in '{}'", dataset, path.display()),
            ));
        }

        scratch.set_len(0)?;
        scratch.seek(SeekFrom::Start(0))?;
        {
            let limited = (&mut input).take(payload_len);
            let mut decoder = Decoder::new(limited)?;
            let mut encoder = zstd::Encoder::new(&mut *scratch, 1)?;
            let mut current_cid = 0u64;
            let mut first = true;
            while let Some(delta) = read_varint_u64_from_reader(&mut decoder)? {
                current_cid = current_cid.checked_add(delta).ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "dataset CID overflow")
                })?;
                if current_cid >= source_group_count {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "dataset {} CID {} is outside [0, {}) in '{}'",
                            dataset,
                            current_cid,
                            source_group_count,
                            path.display()
                        ),
                    ));
                }
                let output_delta = if first {
                    current_cid.checked_add(cid_offset).ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidData, "shifted CID overflow")
                    })?
                } else {
                    delta
                };
                write_varint_u64_to_writer(output_delta, &mut encoder)?;
                first = false;
            }
            encoder.finish()?;
        }

        let output_payload_len = scratch.seek(SeekFrom::End(0))?;
        scratch.seek(SeekFrom::Start(0))?;
        output_offsets.push(usize::try_from(*output_position).map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, "id-to-CID offset cannot fit usize")
        })?);
        output.write_all(&output_payload_len.to_le_bytes())?;
        io::copy(scratch, output)?;
        *output_position = output_position
            .checked_add(8)
            .and_then(|position| position.checked_add(output_payload_len))
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "id-to-CID size overflow"))?;
    }

    let mut terminator = [0u8; 8];
    input.read_exact(&mut terminator)?;
    if u64::from_le_bytes(terminator) != 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("'{}' contains more datasets than filenames", path.display()),
        ));
    }
    Ok(())
}

fn append_unshifted_dataset_cids(
    path: &Path,
    expected_datasets: usize,
    output: &mut BufWriter<File>,
    output_position: &mut u64,
    output_offsets: &mut Vec<usize>,
) -> io::Result<()> {
    let mut input = require_magic(path, ID_TO_CID_MAGIC)?;
    for dataset in 0..expected_datasets {
        let mut len_buf = [0u8; 8];
        input.read_exact(&mut len_buf).map_err(|err| {
            io::Error::new(
                err.kind(),
                format!("missing dataset {} CID payload in '{}': {err}", dataset, path.display()),
            )
        })?;
        let payload_len = u64::from_le_bytes(len_buf);
        if payload_len == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("premature CID terminator at dataset {} in '{}'", dataset, path.display()),
            ));
        }
        output_offsets.push(usize::try_from(*output_position).map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, "id-to-CID offset cannot fit usize")
        })?);
        output.write_all(&len_buf)?;
        let copied = io::copy(&mut (&mut input).take(payload_len), output)?;
        if copied != payload_len {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!("truncated dataset {} CID payload in '{}'", dataset, path.display()),
            ));
        }
        *output_position = (*output_position)
            .checked_add(8)
            .and_then(|position| position.checked_add(payload_len))
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "id-to-CID size overflow"))?;
    }
    let mut terminator = [0u8; 8];
    input.read_exact(&mut terminator)?;
    if u64::from_le_bytes(terminator) != 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("'{}' contains more datasets than filenames", path.display()),
        ));
    }
    Ok(())
}

fn merge_compact_archives_structurally(
    archive_a_root: &Path,
    archive_b_root: &Path,
    output_root: &Path,
    temp_dir: &str,
) -> io::Result<()> {
    let timer = Instant::now();
    let archive_a = CompactArchiveLayout::inspect(archive_a_root)?;
    let archive_b = CompactArchiveLayout::inspect(archive_b_root)?;
    let output_groups = archive_a
        .group_count
        .checked_add(archive_b.group_count)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "merged group count overflow"))?;
    let output_entries = archive_a
        .position_entries
        .checked_add(archive_b.position_entries)
        .and_then(|entries| entries.checked_sub(1))
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "merged position count overflow"))?;

    fs::create_dir_all(output_root)?;
    let output_dir = normalize_output_dir(output_root);

    let mut tigs_output = BufWriter::with_capacity(
        16 * 1024 * 1024,
        File::create(output_root.join("tigs_kloe.fa"))?,
    );
    let copied_a = copy_file_into(&archive_a.root.join("tigs_kloe.fa"), &mut tigs_output)?;
    let copied_b = copy_file_into(&archive_b.root.join("tigs_kloe.fa"), &mut tigs_output)?;
    tigs_output.flush()?;
    if copied_a != archive_a.packed_tigs_bytes || copied_b != archive_b.packed_tigs_bytes {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "packed tig file changed while it was being merged",
        ));
    }

    let mut sizes_output = BufWriter::with_capacity(
        16 * 1024 * 1024,
        File::create(output_root.join("bucket_sizes.txt"))?,
    );
    sizes_output.write_all(BUCKET_SIZES_MAGIC)?;
    append_compact_file_body(
        &archive_a.root.join("bucket_sizes.txt"),
        BUCKET_SIZES_MAGIC,
        &mut sizes_output,
    )?;
    append_compact_file_body(
        &archive_b.root.join("bucket_sizes.txt"),
        BUCKET_SIZES_MAGIC,
        &mut sizes_output,
    )?;
    sizes_output.flush()?;

    let mut positions_output = BufWriter::with_capacity(
        16 * 1024 * 1024,
        File::create(output_root.join("positions_kloe.bin"))?,
    );
    positions_output.write_all(POSITIONS_MAGIC)?;
    write_varint_u64_to_writer(output_entries, &mut positions_output)?;
    append_position_deltas(
        &archive_a.root.join("positions_kloe.bin"),
        false,
        &mut positions_output,
    )?;
    append_position_deltas(
        &archive_b.root.join("positions_kloe.bin"),
        true,
        &mut positions_output,
    )?;
    positions_output.flush()?;

    let scratch_root = if temp_dir.is_empty() {
        output_root
    } else {
        Path::new(temp_dir)
    };
    fs::create_dir_all(scratch_root)?;
    let mut scratch = tempfile::tempfile_in(scratch_root)?;
    let mut cid_output = BufWriter::with_capacity(
        16 * 1024 * 1024,
        File::create(output_root.join("id_to_color_id.txt.zst"))?,
    );
    cid_output.write_all(ID_TO_CID_MAGIC)?;
    let mut cid_output_position = ID_TO_CID_MAGIC.len() as u64;
    let mut dataset_offsets = Vec::with_capacity(
        archive_a
            .filenames
            .len()
            .saturating_add(archive_b.filenames.len()),
    );
    append_unshifted_dataset_cids(
        &archive_a.root.join("id_to_color_id.txt.zst"),
        archive_a.filenames.len(),
        &mut cid_output,
        &mut cid_output_position,
        &mut dataset_offsets,
    )?;
    append_shifted_dataset_cids(
        &archive_b.root.join("id_to_color_id.txt.zst"),
        archive_b.filenames.len(),
        archive_b.group_count,
        archive_a.group_count,
        &mut cid_output,
        &mut cid_output_position,
        &mut dataset_offsets,
        &mut scratch,
    )?;
    cid_output.write_all(&0u64.to_le_bytes())?;
    cid_output.flush()?;

    let mut filenames = archive_a.filenames;
    filenames.extend(archive_b.filenames);
    compress::write_filenames_id_offsets(&output_dir, &filenames, &dataset_offsets)?;

    println!(
        "Structural archive merge complete: datasets={}, groups={}, packed_tigs_bytes={}, elapsed={:.3}s",
        filenames.len(),
        output_groups,
        copied_a.saturating_add(copied_b),
        timer.elapsed().as_secs_f64()
    );
    Ok(())
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
    offsets: Arc<Vec<usize>>,
    values: Arc<Vec<u32>>,
}

#[derive(Debug)]
enum ArchiveDatasetIds {
    Memory(CidDatasetIds),
    Disk(Arc<compress::CidDatasetSidecar>),
}

impl ArchiveDatasetIds {
    fn len(&self) -> usize {
        match self {
            Self::Memory(ids) => ids.len(),
            Self::Disk(sidecar) => sidecar.len(),
        }
    }

    fn add_to_source_map(
        &self,
        map: &mut compress::SourceDatasetMap,
        dataset_offset: u32,
    ) -> io::Result<()> {
        match self {
            Self::Memory(ids) => map.add_dense_storage(
                Arc::clone(&ids.values),
                Arc::clone(&ids.offsets),
                dataset_offset,
            ),
            Self::Disk(sidecar) => map.add_disk_storage(Arc::clone(sidecar), dataset_offset),
        }
    }

    fn ensure_disk_backed(&mut self, temp_path: &Path) -> io::Result<()> {
        let Self::Memory(ids) = self else {
            return Ok(());
        };
        let mut writer = compress::CidDatasetSidecarWriter::create(temp_path)?;
        for range in ids.offsets.windows(2) {
            writer.append_one_based(&ids.values[range[0]..range[1]])?;
        }
        writer.finish()?;
        let sidecar = Arc::new(compress::CidDatasetSidecar::open(temp_path)?);
        if sidecar.len() != ids.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "temporary CID sidecar group count mismatch",
            ));
        }
        *self = Self::Disk(sidecar);
        fs::remove_file(temp_path)?;
        Ok(())
    }
}

impl CidDatasetIds {
    fn len(&self) -> usize {
        self.offsets.len().saturating_sub(1)
    }
}

#[derive(Debug)]
struct CompactSizesBlock {
    first_group: usize,
    group_count: usize,
    offsets_file_offset: u64,
    data_offset: u64,
    data_len: usize,
}

#[derive(Debug)]
struct CompactSizesIndex {
    path: PathBuf,
    blocks: Vec<CompactSizesBlock>,
    length_prefixed: bool,
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
    cid_start: usize,
    cid_end: usize,
    first_color: ColorIndexType,
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
        let mut tigs_reader = PackedTigsReader::open(&self.tigs_path)?;
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
        let mut cached_block_ranges = Vec::new();
        let mut cached_block_data = Vec::new();
        let mut encoded = Vec::new();
        let mut sequence = Vec::new();
        let mut sizes = Vec::new();

        for cid in block_data.cid_start..block_data.cid_end {
            let color_offset = cid - block_data.cid_start;
            let color = ColorIndexType::try_from(block_data.first_color as usize + color_offset)
                .expect("validated archive color range must fit ColorIndexType");
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
                    (cached_block_ranges, cached_block_data) = compact.load_block_from(
                        block_index,
                        compact_sizes_reader
                            .as_mut()
                            .expect("compact sizes reader must be initialized"),
                    )?;
                    cached_block_index = block_index;
                }
                compact.decode_group_into(
                    block_index,
                    group,
                    &cached_block_ranges,
                    &cached_block_data,
                    &mut sizes,
                )?;
            } else {
                read_bucket_sizes_at_into(
                    legacy_sizes_reader
                        .as_mut()
                        .expect("legacy sizes reader must be initialized"),
                    sizes_pos,
                    &mut sizes,
                )?;
            }

            let mut current_tigs_pos = tigs_pos;
            for &size in &sizes {
                if size == 0 {
                    continue;
                }
                encoded.resize(size.div_ceil(4), 0);
                tigs_reader.read_exact_at(current_tigs_pos, &mut encoded)?;
                current_tigs_pos = current_tigs_pos
                    .checked_add(encoded.len() as u64)
                    .ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidData, "packed-tig offset overflow")
                    })?;
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
        let cid_sidecar_path = root.join(compress::CID_TO_DATASET_FILE);
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
        let expected_tigs_len = positions.tigs(positions.len() - 1).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "positions omit final tig offset")
        })?;
        let actual_tigs_len = PackedTigsReader::open(&tigs_path)?.logical_len();
        if actual_tigs_len != expected_tigs_len {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "packed-tig logical size mismatch in '{}': positions={}, tigs={}",
                    root.display(),
                    expected_tigs_len,
                    actual_tigs_len
                ),
            ));
        }

        let filenames = load_filenames(&filenames_path)?;
        let expected_groups = positions.len().saturating_sub(1);
        let cid_to_ids = if cid_sidecar_path.is_file() {
            let sidecar = Arc::new(compress::CidDatasetSidecar::open(&cid_sidecar_path)?);
            if sidecar.len() != expected_groups {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "CID sidecar '{}' contains {} groups, expected {}",
                        cid_sidecar_path.display(),
                        sidecar.len(),
                        expected_groups
                    ),
                ));
            }
            println!(
                "Using disk-backed CID memberships from {}",
                cid_sidecar_path.display()
            );
            ArchiveDatasetIds::Disk(sidecar)
        } else {
            eprintln!(
                "Archive '{}' predates disk-backed CID memberships; using the compatibility in-memory transpose",
                root.display()
            );
            ArchiveDatasetIds::Memory(load_cid_to_ids(
                &cid_path,
                filenames.len(),
                expected_groups,
            )?)
        };

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
        first_color: usize,
        target_blocks: usize,
    ) -> io::Result<ArchiveSequencesStream> {
        let cid_count = self.cid_to_ids.len();
        let last_color = first_color.checked_add(cid_count).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "merge source color count overflow")
        })?;
        ColorIndexType::try_from(last_color.saturating_sub(1)).map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, "too many merge source colors")
        })?;
        let mut blocks = Vec::with_capacity(target_blocks.min(cid_count).max(1));
        if cid_count > 0 {
            let first_byte = self.positions.tigs(0).unwrap_or(0);
            let end_byte = self.positions.tigs(cid_count).unwrap_or(first_byte);
            let total_bytes = end_byte.saturating_sub(first_byte);
            let block_count = target_blocks.max(1).min(cid_count);
            let mut current_start = 0usize;
            for block_number in 1..block_count {
                let target_byte = first_byte.saturating_add(
                    ((total_bytes as u128 * block_number as u128) / block_count as u128) as u64,
                );
                let mut current_end = self
                    .positions
                    .tigs
                    .partition_point(|&position| position < target_byte)
                    .clamp(current_start + 1, cid_count);
                if cid_count - current_end < block_count - block_number {
                    current_end = cid_count - (block_count - block_number);
                }
                let start_byte = self.positions.tigs(current_start).unwrap_or(first_byte);
                let block_end_byte = self.positions.tigs(current_end).unwrap_or(start_byte);
                blocks.push(ArchiveSequenceBlock {
                    cid_start: current_start,
                    cid_end: current_end,
                    first_color: ColorIndexType::try_from(first_color + current_start)
                        .expect("validated archive color range must fit ColorIndexType"),
                    estimated_bases: block_end_byte
                        .saturating_sub(start_byte)
                        .saturating_mul(4)
                        .max(1),
                });
                current_start = current_end;
            }
            let start_byte = self.positions.tigs(current_start).unwrap_or(first_byte);
            blocks.push(ArchiveSequenceBlock {
                cid_start: current_start,
                cid_end: cid_count,
                first_color: ColorIndexType::try_from(first_color + current_start)
                    .expect("validated archive color range must fit ColorIndexType"),
                estimated_bases: end_byte
                    .saturating_sub(start_byte)
                    .saturating_mul(4)
                    .max(1),
            });
        }
        Ok(ArchiveSequencesStream {
            positions: Arc::clone(&self.positions),
            tigs_path: self.tigs_path.clone(),
            sizes_path: self.sizes_path.clone(),
            compact_sizes: self.compact_sizes.clone(),
            blocks,
        })
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

    let requested_tig_rebuild = cfg.use_unitigs || cfg.use_matchtigs || cfg.use_eulertigs;
    if cfg.structural {
        if requested_tig_rebuild {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--structural-merge cannot change the SPSS mode; remove the tig output flag",
            ));
        }
        println!(
            "Structurally joining packed KLOE archives (no graph rebuild or cross-archive deduplication)"
        );
        if cfg.verify_kmers {
            eprintln!(
                "Warning: --verify-kmers is not needed for the byte-preserving structural merge."
            );
        }
        if cfg.skip_sort {
            eprintln!("Warning: --skip-sort has no effect on structural archive merge.");
        }
        merge_compact_archives_structurally(
            archive_a_root,
            archive_b_root,
            output_root,
            &cfg.temp_dir,
        )?;
        println!(
            "Merge complete: output archive written to {} (structural SPSS)",
            output_root.to_string_lossy()
        );
        return Ok(());
    }

    fs::create_dir_all(output_root)?;
    let output_dir_norm = normalize_output_dir(output_root);

    println!(
        "Loading archive A from {}",
        archive_a_root.to_string_lossy()
    );
    let mut archive_a = ArchiveInfo::load(archive_a_root, k)?;
    let legacy_a_sidecar = output_root.join(format!(
        ".kloe-merge-source-a-{}.bin",
        std::process::id()
    ));
    archive_a
        .cid_to_ids
        .ensure_disk_backed(&legacy_a_sidecar)?;
    println!(
        "Loading archive B from {}",
        archive_b_root.to_string_lossy()
    );
    let mut archive_b = ArchiveInfo::load(archive_b_root, k)?;
    let legacy_b_sidecar = output_root.join(format!(
        ".kloe-merge-source-b-{}.bin",
        std::process::id()
    ));
    archive_b
        .cid_to_ids
        .ensure_disk_backed(&legacy_b_sidecar)?;

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

    let a_source_count = archive_a.cid_to_ids.len();
    let b_source_count = archive_b.cid_to_ids.len();
    let source_count = a_source_count.checked_add(b_source_count).ok_or_else(|| {
        io::Error::new(io::ErrorKind::InvalidData, "merge source color count overflow")
    })?;
    let mut source_dataset_ids = compress::SourceDatasetMap::new();
    archive_a.cid_to_ids.add_to_source_map(&mut source_dataset_ids, 0)?;
    archive_b
        .cid_to_ids
        .add_to_source_map(&mut source_dataset_ids, offset)?;
    debug_assert_eq!(source_dataset_ids.len(), source_count);

    let blocks_per_archive = cfg.threads.max(1).div_ceil(2);
    let a_stream = Arc::new(archive_a.sequence_stream(0, blocks_per_archive)?);
    let b_stream = Arc::new(archive_b.sequence_stream(a_source_count, blocks_per_archive)?);
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
        Ok(()) => Ok(&magic == BUCKET_SIZES_MAGIC || &magic == BUCKET_SIZES_COMPACT_MAGIC),
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

    fn load_block_from(
        &self,
        block_index: usize,
        file: &mut File,
    ) -> io::Result<(Vec<(usize, usize)>, Vec<u8>)> {
        let block = self.blocks.get(block_index).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "compact size block is not indexed")
        })?;
        let mut offsets = Vec::new();
        if !self.length_prefixed {
            file.seek(std::io::SeekFrom::Start(block.offsets_file_offset))?;
            offsets.resize(block.group_count + 1, 0u32);
            let mut bytes = [0u8; 4];
            for offset in &mut offsets {
                file.read_exact(&mut bytes)?;
                *offset = u32::from_le_bytes(bytes);
            }
        }
        file.seek(std::io::SeekFrom::Start(block.data_offset))?;
        let mut compressed = vec![0u8; block.data_len];
        file.read_exact(&mut compressed)?;
        let mut decompressed = Vec::new();
        Decoder::new(&compressed[..])?.read_to_end(&mut decompressed)?;
        let mut ranges = Vec::with_capacity(block.group_count);
        if self.length_prefixed {
            let mut cursor = std::io::Cursor::new(decompressed.as_slice());
            for _ in 0..block.group_count {
                let group_len = read_varint_u64_from_reader(&mut cursor)?.ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::UnexpectedEof,
                        "missing compact size group length",
                    )
                })?;
                let start = cursor.position() as usize;
                let end = start
                    .checked_add(usize::try_from(group_len).map_err(|_| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            "compact size group is too large",
                        )
                    })?)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            "compact size group offset overflow",
                        )
                    })?;
                if end > decompressed.len() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("truncated compact size block in '{}'", self.path.display()),
                    ));
                }
                ranges.push((start, end));
                cursor.set_position(end as u64);
            }
            if cursor.position() as usize != decompressed.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("trailing compact size data in '{}'", self.path.display()),
                ));
            }
        } else {
            let expected_min = offsets.last().copied().unwrap_or(0) as usize;
            if decompressed.len() < expected_min
                || offsets.windows(2).any(|range| range[0] > range[1])
            {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("truncated compact size block in '{}'", self.path.display()),
                ));
            }
            ranges.extend(
                offsets
                    .windows(2)
                    .map(|range| (range[0] as usize, range[1] as usize)),
            );
        }
        Ok((ranges, decompressed))
    }

    fn decode_group_into(
        &self,
        block_index: usize,
        group: usize,
        ranges: &[(usize, usize)],
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
        let (start, end) = *ranges.get(local).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "missing compact size group range")
        })?;
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
    if &magic != BUCKET_SIZES_MAGIC && &magic != BUCKET_SIZES_COMPACT_MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("expected compact bucket sizes magic in '{}'", path.display()),
        ));
    }
    let length_prefixed = &magic == BUCKET_SIZES_COMPACT_MAGIC;

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
        let offsets_file_offset = reader.stream_position()?;
        let offsets_bytes = if length_prefixed {
            0
        } else {
            (group_count + 1)
                .checked_mul(std::mem::size_of::<u32>())
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "compact size offsets overflow")
                })?
        };
        reader.seek(std::io::SeekFrom::Current(offsets_bytes as i64))?;
        let data_offset = reader.stream_position()?;
        let next_offset = data_offset.checked_add(data_len as u64).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "compact size file offset overflow")
        })?;
        reader.seek(std::io::SeekFrom::Start(next_offset))?;
        blocks.push(CompactSizesBlock {
            first_group,
            group_count,
            offsets_file_offset,
            data_offset,
            data_len,
        });
        first_group = first_group.saturating_add(group_count);
    }

    Ok(CompactSizesIndex {
        path: path.to_path_buf(),
        blocks,
        length_prefixed,
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

    if &magic == POSITIONS_MAGIC || &magic == POSITIONS_COMPACT_MAGIC {
        let implicit_sizes = &magic == POSITIONS_COMPACT_MAGIC;
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
        for entry in 0..count {
            let dt = read_varint_u64_from_reader(&mut file)?.ok_or_else(|| {
                io::Error::new(io::ErrorKind::UnexpectedEof, "truncated tigs delta")
            })?;
            let ds = if implicit_sizes {
                u64::from(entry > 0)
            } else {
                read_varint_u64_from_reader(&mut file)?.ok_or_else(|| {
                    io::Error::new(io::ErrorKind::UnexpectedEof, "truncated sizes delta")
                })?
            };
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
        offsets: Arc::new(offsets),
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
    let legacy_index = root.join("id_to_color_id.txt.zst");
    let compact_index = root.join(compress::CID_TO_DATASET_FILE);
    if !legacy_index.is_file() && !compact_index.is_file() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!(
                "archive '{}' has neither a CID-to-dataset nor dataset-to-CID membership index",
                root.display()
            ),
        ));
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
