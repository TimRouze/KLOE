use ggcat_api::{
    ColorIndexType, ExtraElaboration, GGCATConfig, GGCATInstance, GeneralSequenceBlockData,
};
use ggcat_colors::colors_manager::ColorMapReader;
use ggcat_colors::storage::deserializer::ColorsDeserializer;
use ggcat_colors::DefaultColorsSerializer;
use rayon::slice::ParallelSliceMut;
use std::cmp::Ordering as CmpOrdering;
use std::collections::{BinaryHeap, HashMap};
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Result, Write};
use std::path::{Path, PathBuf};
use std::sync::{mpsc, Arc};
use std::thread;
use std::time::Instant;
use tempfile::{Builder as TempBuilder, TempPath};
use zstd::Encoder;

use crate::records::{SimplitigBatch, SimplitigRecord};
use crate::utils::{Convert, Converter};

const IO_BUFFER_CAPACITY: usize = 8 * 1024 * 1024;
const ENCODED_SEQ_BUFFER_TARGET: usize = 4 * 1024 * 1024;
const ID_CID_SPILL_BUFFER_CAPACITY: usize = 256 * 1024;
const PAR_SORT_THRESHOLD: usize = 200_000;
const COLOR_RECORD_BATCH_SIZE: usize = 4_096;
// Larger chunk targets reduce temporary chunk count and k-way merge overhead
// in the post-ggcat external sort stage.
const COLOR_CHUNK_TARGET_BYTES: usize = 512 * 1024 * 1024;
const GROUP_SORT_SPILL_BYTES: usize = 128 * 1024 * 1024;
const BUCKET_SIZES_FORMAT_MAGIC: [u8; 8] = *b"KLSZV001";
const BUCKET_SIZE_COMPRESSION_LEVEL: i32 = 6;

#[inline]
fn write_uvarint<W: Write>(writer: &mut W, mut value: usize) -> io::Result<()> {
    let mut buf = [0u8; 10];
    let mut idx = 0usize;
    while value >= 0x80 {
        buf[idx] = ((value as u8) & 0x7f) | 0x80;
        value >>= 7;
        idx += 1;
    }
    buf[idx] = value as u8;
    writer.write_all(&buf[..=idx])
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
}

impl Default for GgcatCompressionConfig {
    fn default() -> Self {
        Self {
            memory_gb: 8,
            temp_dir: String::new(),
        }
    }
}

#[derive(Debug)]
struct SortedColorRecord {
    key: Vec<u8>,
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
        self.record.key == other.record.key && self.record.seq == other.record.seq
    }
}

impl Ord for ChunkHeapItem {
    fn cmp(&self, other: &Self) -> CmpOrdering {
        other
            .record
            .key
            .cmp(&self.record.key)
            .then_with(|| other.record.seq.cmp(&self.record.seq))
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

pub(crate) fn write_filenames_id_offsets(
    output_dir: &str,
    filenames: &[String],
    id_cid_line_sizes: &[usize],
) -> Result<()> {
    if filenames.len() != id_cid_line_sizes.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "filenames count ({}) differs from id->cid offsets count ({})",
                filenames.len(),
                id_cid_line_sizes.len()
            ),
        ));
    }

    let mut fof_id = BufWriter::new(File::create(output_dir.to_owned() + "filenames_id.txt")?);
    for (filename, offset) in filenames.iter().zip(id_cid_line_sizes.iter()) {
        fof_id.write_all(format!("{filename}:{offset}\n").as_bytes())?;
    }
    fof_id.flush()?;
    Ok(())
}

fn to_io_err(context: &str, err: impl std::fmt::Display) -> io::Error {
    io::Error::other(format!("{context}: {err}"))
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
        // Hard-switch default mode to simplitigs when no explicit tig flag is provided.
        ExtraElaboration::FastSimplitigs
    }
}

fn write_sorted_color_record(
    writer: &mut BufWriter<File>,
    record: &SortedColorRecord,
) -> Result<()> {
    writer.write_all(&(record.key.len() as u32).to_le_bytes())?;
    writer.write_all(&(record.seq.len() as u32).to_le_bytes())?;
    writer.write_all(&record.key)?;
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
    let key_len = u32::from_le_bytes(len_buf) as usize;

    reader.read_exact(&mut len_buf)?;
    let seq_len = u32::from_le_bytes(len_buf) as usize;

    let mut key = vec![0u8; key_len];
    reader.read_exact(&mut key)?;

    let mut seq = vec![0u8; seq_len];
    reader.read_exact(&mut seq)?;

    Ok(Some(SortedColorRecord { key, seq }))
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
) -> Result<()> {
    if records.is_empty() {
        return Ok(());
    }

    if records.len() >= PAR_SORT_THRESHOLD {
        records.par_sort_unstable_by(|left, right| {
            left.key
                .cmp(&right.key)
                .then_with(|| left.seq.cmp(&right.seq))
        });
    } else {
        records.sort_unstable_by(|left, right| {
            left.key
                .cmp(&right.key)
                .then_with(|| left.seq.cmp(&right.seq))
        });
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

fn bitset_to_dataset_ids(bitset: &[u8], dataset_count: usize) -> Result<Vec<u32>> {
    if bitset.len() != dataset_count {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "bitset width ({}) differs from dataset count ({dataset_count})",
                bitset.len()
            ),
        ));
    }

    let mut ids = Vec::new();
    for (idx, bit) in bitset.iter().enumerate() {
        match bit {
            b'1' => ids.push((idx + 1) as u32),
            b'0' => {}
            other => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid bitset character '{}'", *other as char),
                ));
            }
        }
    }

    if ids.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "encountered empty color bitset",
        ));
    }

    Ok(ids)
}

fn sort_and_spill_sequence_run(
    seqs: &mut Vec<Vec<u8>>,
    run_paths: &mut Vec<TempPath>,
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

    run_paths.push(tmp.into_temp_path());
    seqs.clear();
    Ok(())
}

fn emit_sorted_group_sequences(
    seqs: &mut Vec<Vec<u8>>,
    run_paths: &mut Vec<TempPath>,
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

    run_paths.clear();
    Ok(())
}

struct StreamWriterState {
    omni_file: BufWriter<File>,
    size_file: BufWriter<File>,
    spill_writers: Vec<BufWriter<File>>,
    pos_nb_unitig: Vec<(u64, u64)>,
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
    ) -> Result<(Self, Vec<PathBuf>, PathBuf)> {
        let mut spill_paths = Vec::with_capacity(nb_files);
        let mut spill_writers = Vec::with_capacity(nb_files);
        let spill_dir = create_id_cid_spill_dir(output_dir)?;
        for id in 0..nb_files {
            let path = spill_dir.join(format!("id_{id}.cids.bin"));
            let writer =
                BufWriter::with_capacity(ID_CID_SPILL_BUFFER_CAPACITY, File::create(&path)?);
            spill_paths.push(path);
            spill_writers.push(writer);
        }

        let mut size_file = BufWriter::with_capacity(
            IO_BUFFER_CAPACITY,
            File::create(output_dir.to_owned() + "bucket_sizes.txt")?,
        );
        size_file.write_all(&BUCKET_SIZES_FORMAT_MAGIC)?;

        Ok((
            Self {
                omni_file: BufWriter::with_capacity(
                    IO_BUFFER_CAPACITY,
                    File::create(unitigs_file_path)?,
                ),
                size_file,
                spill_writers,
                pos_nb_unitig: vec![(0, BUCKET_SIZES_FORMAT_MAGIC.len() as u64)],
                prev_tigs_size: 0,
                prev_bucket_pos: BUCKET_SIZES_FORMAT_MAGIC.len() as u64,
                cid: 0,
                encoded_seq_buffer: Vec::with_capacity(ENCODED_SEQ_BUFFER_TARGET),
            },
            spill_paths,
            spill_dir,
        ))
    }

    fn write_group(
        &mut self,
        dataset_ids_zero_based: &[usize],
        seqs: &mut Vec<Vec<u8>>,
        run_paths: &mut Vec<TempPath>,
        spill_dir: &Path,
    ) -> Result<()> {
        if seqs.is_empty() && run_paths.is_empty() {
            return Ok(());
        }

        let mut group_sizes_buffer: Vec<u8> = Vec::new();
        let mut group_encoder =
            Encoder::new(&mut group_sizes_buffer, BUCKET_SIZE_COMPRESSION_LEVEL)?;

        let mut prev_size: usize = 0;
        emit_sorted_group_sequences(seqs, run_paths, spill_dir, |seq| {
            let encoded_seq = <Converter as Convert<&[u8]>>::str2num(seq);
            self.prev_tigs_size += encoded_seq.len() as u64;
            self.encoded_seq_buffer.extend_from_slice(&encoded_seq);
            if self.encoded_seq_buffer.len() >= ENCODED_SEQ_BUFFER_TARGET {
                self.omni_file.write_all(&self.encoded_seq_buffer)?;
                self.encoded_seq_buffer.clear();
            }

            let size = seq.len();
            let delta = size.checked_sub(prev_size).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "group sizes are not nondecreasing",
                )
            })?;
            write_uvarint(&mut group_encoder, delta)?;
            prev_size = size;
            Ok(())
        })?;

        group_encoder.finish()?;

        self.prev_bucket_pos += (8 + group_sizes_buffer.len()) as u64;
        self.pos_nb_unitig
            .push((self.prev_tigs_size, self.prev_bucket_pos));
        self.size_file
            .write_all(&(group_sizes_buffer.len() as u64).to_le_bytes())?;
        self.size_file.write_all(&group_sizes_buffer)?;

        for &id in dataset_ids_zero_based {
            self.spill_writers[id].write_all(&(self.cid as u64).to_le_bytes())?;
        }
        self.cid += 1;

        Ok(())
    }

    fn finalize(mut self) -> Result<Vec<(u64, u64)>> {
        if !self.encoded_seq_buffer.is_empty() {
            self.omni_file.write_all(&self.encoded_seq_buffer)?;
            self.encoded_seq_buffer.clear();
        }
        for writer in &mut self.spill_writers {
            writer.flush()?;
        }
        self.omni_file.flush()?;
        self.size_file.flush()?;

        println!(
            "Completed compression: total tigs={}, total sizes={}",
            self.prev_tigs_size, self.prev_bucket_pos
        );
        Ok(self.pos_nb_unitig)
    }
}

fn write_compressed_from_stream(
    unitigs_file_path: String,
    output_dir: &String,
    nb_files: u32,
    record_rx: mpsc::Receiver<SimplitigBatch>,
) -> Result<(Vec<(u64, u64)>, Vec<PathBuf>, PathBuf)> {
    let (mut writer, spill_paths, spill_dir) =
        StreamWriterState::new(unitigs_file_path, output_dir, nb_files as usize)?;

    let group_spill_dir = Path::new(output_dir).join(format!(
        ".kloe-group-sort-spill-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0)
    ));
    fs::create_dir_all(&group_spill_dir)?;

    let mut current_color_ids: Option<Arc<Vec<u32>>> = None;
    let mut current_ids_zero_based: Option<Vec<usize>> = None;
    let mut group_seqs: Vec<Vec<u8>> = Vec::new();
    let mut group_run_paths: Vec<TempPath> = Vec::new();
    let mut group_bytes: usize = 0;

    for batch in record_rx {
        for record in batch {
            let key_changed = current_color_ids
                .as_ref()
                .is_none_or(|ids| ids.as_ref() != record.color_ids.as_ref());

            if key_changed {
                if let Some(ref ids) = current_ids_zero_based {
                    writer.write_group(
                        ids,
                        &mut group_seqs,
                        &mut group_run_paths,
                        &group_spill_dir,
                    )?;
                }

                let ids = record
                    .color_ids
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
                    if id >= writer.spill_writers.len() {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!(
                                "dataset id {} outside [0, {})",
                                id,
                                writer.spill_writers.len()
                            ),
                        ));
                    }
                }

                current_ids_zero_based = Some(ids);
                current_color_ids = Some(Arc::clone(&record.color_ids));
                group_bytes = 0;
            }

            group_bytes += record.seq.len();
            group_seqs.push(record.seq);
            if group_bytes >= GROUP_SORT_SPILL_BYTES {
                sort_and_spill_sequence_run(
                    &mut group_seqs,
                    &mut group_run_paths,
                    &group_spill_dir,
                )?;
                group_bytes = 0;
            }
        }
    }

    if let Some(ref ids) = current_ids_zero_based {
        writer.write_group(ids, &mut group_seqs, &mut group_run_paths, &group_spill_dir)?;
    }

    let pos_nb_unitig = writer.finalize()?;
    let _ = fs::remove_dir_all(&group_spill_dir);

    Ok((pos_nb_unitig, spill_paths, spill_dir))
}

fn write_positions(pos_nb_unitigs: Vec<(u64, u64)>, filepath: String) -> Result<()> {
    let mut pos_file = BufWriter::new(File::create(filepath)?);
    for elem in &pos_nb_unitigs {
        pos_file.write_all(&elem.0.to_le_bytes())?;
        pos_file.write_all(&elem.1.to_le_bytes())?;
    }
    pos_file.flush()?;
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

fn write_id_to_color_id_from_spills(
    cid_file_path: String,
    spill_paths: &[PathBuf],
    spill_dir: &Path,
) -> std::io::Result<Vec<usize>> {
    let mut cid_file = BufWriter::with_capacity(IO_BUFFER_CAPACITY, File::create(&cid_file_path)?);
    let mut id_cid_line_sizes = Vec::with_capacity(spill_paths.len());
    let mut tot_size = 0usize;

    for path in spill_paths {
        let mut reader = BufReader::with_capacity(IO_BUFFER_CAPACITY, File::open(path)?);
        let mut payload = Vec::new();
        {
            let mut cid_encoder = Encoder::new(&mut payload, 1)?;
            let mut first = true;
            let mut prev_cid = 0usize;
            loop {
                let mut raw = [0u8; 8];
                match std::io::Read::read_exact(&mut reader, &mut raw) {
                    Ok(()) => {
                        let raw_cid = u64::from_le_bytes(raw);
                        let cid = usize::try_from(raw_cid).map_err(|_| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!(
                                    "color id {} cannot be represented as usize while writing {}",
                                    raw_cid, cid_file_path
                                ),
                            )
                        })?;
                        if cid < prev_cid {
                            return Err(io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!(
                                    "color ids are not sorted ({} before {}) while writing {}",
                                    prev_cid, cid, cid_file_path
                                ),
                            ));
                        }
                        let delta = cid - prev_cid;
                        if !first {
                            cid_encoder.write_all(b",")?;
                        }
                        first = false;
                        write!(&mut cid_encoder, "{}", delta)?;
                        prev_cid = cid;
                    }
                    Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => break,
                    Err(e) => return Err(e),
                }
            }
            cid_encoder.finish()?;
        }

        id_cid_line_sizes.push(tot_size);
        tot_size += 8 + payload.len();
        cid_file.write_all(&(payload.len() as u64).to_le_bytes())?;
        cid_file.write_all(&payload)?;

        let _ = fs::remove_file(path);
    }

    cid_file.write_all(&(0_u64).to_le_bytes())?;
    cid_file.flush()?;
    let _ = fs::remove_dir(spill_dir);
    Ok(id_cid_line_sizes)
}

pub(crate) fn sort_by_bucket_streaming(
    output_dir: &String,
    nb_files: u32,
    record_rx: mpsc::Receiver<SimplitigBatch>,
) -> Vec<usize> {
    let write_time = Instant::now();
    println!("Starting writing compressed sequences (streaming).");
    let triple = match write_compressed_from_stream(
        output_dir.clone() + "tigs_kloe.fa",
        output_dir,
        nb_files,
        record_rx,
    ) {
        Ok(res_pair) => res_pair,
        Err(e) => panic!("Error writing compressed unitigs: {e:?}"),
    };
    println!(
        "Writing compressed sequences wall time: {:.3}s",
        write_time.elapsed().as_secs_f64()
    );
    let position_time = Instant::now();
    let pos_nb_unitig = triple.0;
    let spill_paths = triple.1;
    let spill_dir = triple.2;

    println!("Starting to write positions");
    if let Err(e) = write_positions(
        pos_nb_unitig,
        String::from(output_dir.clone() + "positions_kloe.bin"),
    ) {
        panic!("Error writting positions: {e:?}");
    }
    println!(
        "Write positions wall time: {:.3}s",
        position_time.elapsed().as_secs_f64()
    );
    let id_time = Instant::now();
    let write_id_cid = match write_id_to_color_id_from_spills(
        output_dir.clone() + "id_to_color_id.txt.zst",
        &spill_paths,
        &spill_dir,
    ) {
        Ok(id_cid_line_sizes) => id_cid_line_sizes,
        Err(e) => panic!("error writting id to color id list: {e:?}"),
    };
    println!(
        "Write id to cid wall time: {:.3}s",
        id_time.elapsed().as_secs_f64()
    );
    println!(
        "Compression took: {:.3}s",
        write_time.elapsed().as_secs_f64()
    );
    write_id_cid
}

fn stream_sorted_records_from_chunks(
    chunk_paths: &[PathBuf],
    mut emit: impl FnMut(SortedColorRecord) -> Result<()>,
) -> Result<()> {
    let mut readers = Vec::with_capacity(chunk_paths.len());
    for path in chunk_paths {
        readers.push(BufReader::with_capacity(
            IO_BUFFER_CAPACITY,
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

fn produce_ggcat_records(
    filenames: Vec<String>,
    threads: usize,
    k: usize,
    m: usize,
    use_unitigs: bool,
    use_matchtigs: bool,
    use_eulertigs: bool,
    ggcat_cfg: GgcatCompressionConfig,
    sender: mpsc::SyncSender<SimplitigBatch>,
) -> Result<()> {
    let dataset_count = filenames.len();
    if dataset_count == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "input file list is empty",
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

    let raw_records_file = ggcat_temp_root.join(format!("kloe-color-records-{pid}-{nanos}.tmp"));
    let chunk_dir = ggcat_temp_root.join(format!("kloe-color-chunks-{pid}-{nanos}"));
    fs::create_dir_all(&chunk_dir)?;

    let instance = GGCATInstance::create(GGCATConfig {
        temp_dir: Some(ggcat_temp_root.clone()),
        memory: ggcat_cfg.memory_gb.max(1) as f64,
        prefer_memory: true,
        total_threads_count: threads.max(1),
        intermediate_compression_level: None,
        stats_file: None,
        messages_callback: None,
    })
    .map_err(|e| to_io_err("create ggcat instance", e))?;

    let input_streams = filenames
        .iter()
        .map(|file| {
            let path = PathBuf::from(file);
            let resolved = fs::canonicalize(&path).unwrap_or(path);
            GeneralSequenceBlockData::FASTA((resolved, None))
        })
        .collect::<Vec<_>>();

    let color_names = (1..=dataset_count)
        .map(|idx| idx.to_string())
        .collect::<Vec<_>>();

    println!(
        "Running embedded ggcat build-colored-fasta path (k={}, m={}, threads={}, memory={}GB)",
        k,
        m,
        threads.max(1),
        ggcat_cfg.memory_gb.max(1)
    );

    let records_output = instance
        .build_graph_color_records(
            input_streams,
            raw_records_file,
            Some(&color_names),
            k,
            threads.max(1),
            false,
            Some(m),
            1,
            ggcat_extra_elaboration(use_unitigs, use_matchtigs, use_eulertigs),
        )
        .map_err(|e| to_io_err("ggcat build_graph_color_records", e))?;

    let colormap_file = GGCATInstance::get_colormap_file(&records_output);
    let mut colors_deserializer =
        ColorsDeserializer::<DefaultColorsSerializer>::new(&colormap_file, true)
            .map_err(|e| to_io_err("open ggcat colormap", e))?;

    let colors_count = colors_deserializer.colors_count();
    if colors_count != dataset_count {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "ggcat returned {colors_count} colors but KLOE loaded {dataset_count} datasets"
            ),
        ));
    }

    let dumped_colors = GGCATInstance::dump_colors(&colormap_file)
        .map_err(|e| to_io_err("dump ggcat colors", e))?
        .collect::<Vec<_>>();
    if dumped_colors.len() != colors_count {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "ggcat dumped {} color names but reports {} colors",
                dumped_colors.len(),
                colors_count
            ),
        ));
    }
    let mut color_index_to_dataset_index = Vec::with_capacity(colors_count);
    for color_name in dumped_colors {
        let parsed = color_name.parse::<usize>().map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "ggcat color name '{}' is not a dataset index; expected numeric names",
                    color_name
                ),
            )
        })?;
        if parsed == 0 || parsed > dataset_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "ggcat color name '{}' maps outside dataset range [1, {}]",
                    color_name, dataset_count
                ),
            ));
        }
        color_index_to_dataset_index.push(parsed - 1);
    }

    let mut subset_cache: HashMap<ColorIndexType, Vec<ColorIndexType>> = HashMap::new();
    let mut bitset = vec![b'0'; colors_count];
    let mut records = Vec::new();
    let mut chunk_files = Vec::new();
    let mut chunk_bytes = 0usize;

    let mut reader = BufReader::with_capacity(IO_BUFFER_CAPACITY, File::open(&records_output)?);
    while let Some((subsets, seq)) = read_color_record(&mut reader)? {
        bitset.fill(b'0');

        if subsets.is_empty() {
            if colors_count == 1 {
                bitset[0] = b'1';
            } else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "ggcat emitted a sequence without color subsets",
                ));
            }
        } else {
            for subset in subsets {
                let mapped = subset_cache.entry(subset).or_insert_with(|| {
                    let mut colors = Vec::new();
                    colors_deserializer.get_color_mappings(subset, &mut colors);
                    colors
                });
                for &color in mapped.iter() {
                    let color_idx = color as usize;
                    let Some(&dataset_idx) = color_index_to_dataset_index.get(color_idx) else {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("ggcat color index {} out of range", color_idx),
                        ));
                    };
                    if dataset_idx >= bitset.len() {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("dataset index {} out of range", dataset_idx),
                        ));
                    }
                    bitset[dataset_idx] = b'1';
                }
            }
        }

        let record = SortedColorRecord {
            key: bitset.clone(),
            seq,
        };
        chunk_bytes += record.key.len() + record.seq.len();
        records.push(record);

        if chunk_bytes >= COLOR_CHUNK_TARGET_BYTES {
            flush_sorted_chunk(&mut records, &mut chunk_files, &chunk_dir)?;
            chunk_bytes = 0;
        }
    }

    flush_sorted_chunk(&mut records, &mut chunk_files, &chunk_dir)?;

    let mut current_key: Option<Vec<u8>> = None;
    let mut current_color_ids: Option<Arc<Vec<u32>>> = None;
    let mut batch: SimplitigBatch = Vec::with_capacity(COLOR_RECORD_BATCH_SIZE);

    stream_sorted_records_from_chunks(&chunk_files, |record| {
        let color_ids = match current_key.as_ref() {
            Some(active_key) if active_key.as_slice() == record.key.as_slice() => {
                Arc::clone(current_color_ids.as_ref().expect("color ids must be set"))
            }
            _ => {
                let ids = Arc::new(bitset_to_dataset_ids(&record.key, dataset_count)?);
                current_key = Some(record.key.clone());
                current_color_ids = Some(Arc::clone(&ids));
                ids
            }
        };

        batch.push(SimplitigRecord {
            color_ids,
            seq: record.seq,
        });

        if batch.len() >= COLOR_RECORD_BATCH_SIZE {
            let out = std::mem::take(&mut batch);
            sender.send(out).map_err(|e| {
                io::Error::new(
                    io::ErrorKind::BrokenPipe,
                    format!("record receiver dropped: {e}"),
                )
            })?;
            batch = Vec::with_capacity(COLOR_RECORD_BATCH_SIZE);
        }

        Ok(())
    })?;

    if !batch.is_empty() {
        sender.send(batch).map_err(|e| {
            io::Error::new(
                io::ErrorKind::BrokenPipe,
                format!("record receiver dropped: {e}"),
            )
        })?;
    }

    drop(sender);

    let _ = fs::remove_file(&records_output);
    let _ = fs::remove_file(&colormap_file);
    for chunk in chunk_files {
        let _ = fs::remove_file(chunk);
    }
    let _ = fs::remove_dir_all(&chunk_dir);

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
            "Warning: --skip-sort is ignored; ggcat build-colored-fasta always sorts by colors."
        );
    }

    let filenames = read_input_fof_filenames(input_fof)?;
    let dataset_count = filenames.len();
    if dataset_count == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "input file list is empty",
        ));
    }

    println!("Compression backend: embedded ggcat build-colored-fasta");

    let (record_tx, record_rx) = mpsc::sync_channel::<SimplitigBatch>(16);
    let producer_filenames = filenames.clone();
    let producer_cfg = ggcat_cfg.clone();

    let producer = thread::spawn(move || {
        produce_ggcat_records(
            producer_filenames,
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
    let id_cid_line_sizes = sort_by_bucket_streaming(output_dir, dataset_count as u32, record_rx);
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

    write_filenames_id_offsets(output_dir, &filenames, &id_cid_line_sizes)?;
    Ok(())
}
