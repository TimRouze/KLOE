use chrono::Utc;
use std::fmt::Write as FmtWrite;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Result, Write};
use std::path::{Path, PathBuf};
use std::sync::mpsc;

use zstd::Encoder;

use crate::parser;
use crate::utils::{Convert, Converter};

const IO_BUFFER_CAPACITY: usize = 8 * 1024 * 1024;
const ENCODED_SEQ_BUFFER_TARGET: usize = 4 * 1024 * 1024;
const ID_CID_SPILL_BUFFER_CAPACITY: usize = 256 * 1024;

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
    let _overall_start = Utc::now();
    let filenames = read_input_fof_filenames(input_fof)?;
    println!("Compression backend: native KLOE parser (ggcat only for oversized partitions)");
    let ggcat_memory_gb = (ggcat_cfg.memory_gb / 2).max(1);
    let parser_ggcat_cfg = parser::PartitionGgcatConfig {
        memory_gb: ggcat_memory_gb,
        native_memory_budget_gb: ggcat_cfg.memory_gb as u64,
        temp_dir: ggcat_cfg.temp_dir.clone(),
        ..parser::PartitionGgcatConfig::default()
    };
    let (record_rx, merge_handle, dataset_count) = parser::run_parser_streaming(
        std::path::PathBuf::from(input_fof),
        std::path::PathBuf::from(output_dir),
        k,
        m,
        partition_power,
        threads,
        verify_kmers,
        skip_sort,
        use_unitigs,
        use_matchtigs,
        use_eulertigs,
        parser_ggcat_cfg,
    )
    .expect("run_parser_streaming failed");

    if dataset_count != filenames.len() {
        eprintln!(
            "Warning: parser dataset count ({}) differs from input file count ({})",
            dataset_count,
            filenames.len()
        );
    }

    println!("Compressing sequences from stream by color bucket");
    let sort_time = Utc::now();
    let id_cid_line_sizes = sort_by_bucket_streaming(output_dir, dataset_count as u32, record_rx);
    parser::log_checkpoint("Compression took:", sort_time);

    match merge_handle.join() {
        Ok(Ok(())) => {}
        Ok(Err(e)) => eprintln!("Merge thread error: {:#}", e),
        Err(_) => eprintln!("Merge thread panicked"),
    }

    write_filenames_id_offsets(output_dir, &filenames, &id_cid_line_sizes)?;
    Ok(())
}

/// Main compression function, runs each compression steps.
/// - read color-set information,
/// - compress and write unitigs grouped by color,
/// - write positions of each color bucket and tigs bucket,
/// - generate mapping from input file ID to color-bucket positions.
///
/// PARAM
/// - `output_dir`: base output directory where Fulgor outputs are located.
/// - `nb_files`: number of input files (used to allocate structures).
///
/// RETURNS
/// - Vector of sizes (cursor positions) per input file used to annotate filenames,
///   the position of cid list in id to cid file for each id (used later during decompression).
pub fn sort_by_bucket(output_dir: &String, nb_files: u32) -> Vec<usize> {
    let write_time = Utc::now();
    // PROCESS AND COMPRESS UNITIGS
    println!("Starting writing compressed sequences.");
    let pair = match write_compressed(output_dir.clone() + "tigs_kloe.fa", output_dir, nb_files) {
        Ok(res_pair) => res_pair,
        Err(e) => panic!("Error writing compressed unitigs: {e:?}"),
    };
    parser::log_checkpoint("Writing compressed sequences wall time:", write_time);
    let position_time = Utc::now();
    let pos_nb_unitig = pair.0;
    let id_to_color_vec = pair.1;

    println!("Starting to write positions");
    let cursor_positions = match write_positions(
        pos_nb_unitig,
        String::from(output_dir.clone() + "positions_kloe.bin"),
    ) {
        Ok(vec) => vec,
        Err(e) => panic!("Error writting positions: {e:?}"),
    };
    parser::log_checkpoint("Write positions wall time:", position_time);
    let id_time = Utc::now();
    // WRITE FILE ID TO COLOR ID FILE
    let write_id_cid = match write_id_to_color_id(
        output_dir.clone() + "id_to_color_id.txt.zst",
        id_to_color_vec,
        cursor_positions,
    ) {
        Ok(id_cid_line_sizes) => id_cid_line_sizes,
        Err(e) => panic!("error writting id to color id list: {e:?}"),
    };
    parser::log_checkpoint("Write id to cid wall time:", id_time);
    parser::log_checkpoint("Compression took:", write_time);
    write_id_cid
}

/// Streaming variant: consumes records from a channel instead of reading from a file.
pub(crate) fn sort_by_bucket_streaming(
    output_dir: &String,
    nb_files: u32,
    record_rx: mpsc::Receiver<parser::SimplitigBatch>,
) -> Vec<usize> {
    let write_time = Utc::now();
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
    parser::log_checkpoint("Writing compressed sequences wall time:", write_time);
    let position_time = Utc::now();
    let pos_nb_unitig = triple.0;
    let spill_paths = triple.1;
    let spill_dir = triple.2;

    println!("Starting to write positions");
    let cursor_positions = match write_positions(
        pos_nb_unitig,
        String::from(output_dir.clone() + "positions_kloe.bin"),
    ) {
        Ok(vec) => vec,
        Err(e) => panic!("Error writting positions: {e:?}"),
    };
    parser::log_checkpoint("Write positions wall time:", position_time);
    let id_time = Utc::now();
    let write_id_cid = match write_id_to_color_id_from_spills(
        output_dir.clone() + "id_to_color_id.txt.zst",
        &spill_paths,
        &cursor_positions,
        &spill_dir,
    ) {
        Ok(id_cid_line_sizes) => id_cid_line_sizes,
        Err(e) => panic!("error writting id to color id list: {e:?}"),
    };
    parser::log_checkpoint("Write id to cid wall time:", id_time);
    parser::log_checkpoint("Compression took:", write_time);
    write_id_cid
}

/// Consume sorted SimplitigRecords from a channel and compress them.
/// Same logic as write_compressed but reads from mpsc::Receiver instead of FASTA file.
fn write_compressed_from_stream(
    unitigs_file_path: String,
    output_dir: &String,
    nb_files: u32,
    record_rx: mpsc::Receiver<parser::SimplitigBatch>,
) -> Result<(Vec<(u64, u64)>, Vec<PathBuf>, PathBuf)> {
    let mut omni_file =
        BufWriter::with_capacity(IO_BUFFER_CAPACITY, File::create(unitigs_file_path)?);
    let mut size_file = BufWriter::with_capacity(
        IO_BUFFER_CAPACITY,
        File::create(output_dir.clone() + "bucket_sizes.txt")?,
    );
    let spill_dir = create_id_cid_spill_dir(output_dir)?;
    let mut spill_paths = Vec::with_capacity(nb_files as usize);
    let mut spill_writers = Vec::with_capacity(nb_files as usize);
    for id in 0..nb_files as usize {
        let path = spill_dir.join(format!("id_{id}.cids.bin"));
        let writer = BufWriter::with_capacity(ID_CID_SPILL_BUFFER_CAPACITY, File::create(&path)?);
        spill_paths.push(path);
        spill_writers.push(writer);
    }

    let mut pos_nb_unitig: Vec<(u64, u64)> = vec![(0, 0)];

    let mut prev_tigs_size: u64 = 0;
    let mut prev_bucket_pos: u64 = 0;
    let mut cid = 0_usize;

    let mut current_color_ids: Option<std::sync::Arc<Vec<u32>>> = None;
    let mut current_ids: Option<Vec<usize>> = None;
    let mut prev_size: usize = 0;
    let mut group_sizes_buffer: Vec<u8> = Vec::new();
    let mut group_encoder: Option<Encoder<&mut Vec<u8>>> = None;
    let mut encoded_seq_buffer: Vec<u8> = Vec::with_capacity(ENCODED_SEQ_BUFFER_TARGET);
    let mut first_group = true;

    for batch in record_rx {
        for record in batch {
            let key_changed = current_color_ids
                .as_ref()
                .map_or(true, |ids| ids.as_ref() != record.color_ids.as_ref());
            if key_changed {
                if let Some(ref current) = current_ids {
                    if let Some(encoder) = group_encoder.take() {
                        encoder.finish()?;
                    }
                    prev_bucket_pos += (8 + group_sizes_buffer.len()) as u64;
                    pos_nb_unitig.push((prev_tigs_size, prev_bucket_pos));
                    size_file.write_all(&(group_sizes_buffer.len() as u64).to_le_bytes())?;
                    size_file.write_all(&group_sizes_buffer)?;

                    for &id in current {
                        spill_writers[id].write_all(&(cid as u64).to_le_bytes())?;
                    }
                    cid += 1;
                    prev_size = 0;
                    group_sizes_buffer.clear();
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
                    if id >= spill_writers.len() {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("dataset id {} outside [0, {})", id, spill_writers.len()),
                        ));
                    }
                }
                current_ids = Some(ids);
                current_color_ids = Some(record.color_ids.clone());
                let level = if first_group { 4 } else { 1 };
                group_encoder = Some(Encoder::new(&mut group_sizes_buffer, level)?);
                first_group = false;
            }

            let encoded_seq = <Converter as Convert<&[u8]>>::str2num(&record.seq);
            prev_tigs_size += encoded_seq.len() as u64;
            encoded_seq_buffer.extend_from_slice(&encoded_seq);
            if encoded_seq_buffer.len() >= ENCODED_SEQ_BUFFER_TARGET {
                omni_file.write_all(&encoded_seq_buffer)?;
                encoded_seq_buffer.clear();
            }

            let size = record.seq.len();
            let delta = size.saturating_sub(prev_size);
            if let Some(ref mut encoder) = group_encoder {
                encoder.write_all(&delta.to_le_bytes())?;
            }
            prev_size = size;
        }
    }

    if !encoded_seq_buffer.is_empty() {
        omni_file.write_all(&encoded_seq_buffer)?;
        encoded_seq_buffer.clear();
    }

    if let Some(current) = current_ids {
        if let Some(encoder) = group_encoder.take() {
            encoder.finish()?;
        }
        prev_bucket_pos += (8 + group_sizes_buffer.len()) as u64;
        pos_nb_unitig.push((prev_tigs_size, prev_bucket_pos));
        size_file.write_all(&(group_sizes_buffer.len() as u64).to_le_bytes())?;
        size_file.write_all(&group_sizes_buffer)?;

        for id in current {
            spill_writers[id].write_all(&(cid as u64).to_le_bytes())?;
        }
    }

    for writer in &mut spill_writers {
        writer.flush()?;
    }
    omni_file.flush()?;
    size_file.flush()?;

    println!(
        "Completed compression: total tigs={}, total sizes={}",
        prev_tigs_size, prev_bucket_pos
    );
    Ok((pos_nb_unitig, spill_paths, spill_dir))
}

/// Read unitigs from file, group them by color bucket then compress.
/// (Legacy path, kept for sort_by_bucket which is used by the non-streaming code path)
fn write_compressed(
    unitigs_file_path: String,
    output_dir: &String,
    nb_files: u32,
) -> Result<(Vec<(u64, u64)>, Vec<Vec<usize>>)> {
    let mut omni_file = BufWriter::new(File::create(unitigs_file_path)?);
    let mut size_file = BufWriter::new(File::create(output_dir.clone() + "bucket_sizes.txt")?);

    let unitigs_file = File::open(output_dir.clone() + "simplitigs.fa.zst")?;
    let decoder = zstd::Decoder::new(unitigs_file)?;
    let reader: Box<dyn BufRead> = Box::new(BufReader::new(decoder));
    let fa_reader = bio::io::fasta::Reader::from_bufread(reader);

    let mut id_to_color_vec: Vec<Vec<usize>> = vec![Vec::new(); nb_files as usize];
    let mut pos_nb_unitig: Vec<(u64, u64)> = vec![(0, 0)];

    let mut prev_tigs_size: u64 = 0;
    let mut prev_bucket_pos: u64 = 0;
    let mut cid = 0_usize;

    let mut current_ids: Option<Vec<usize>> = None;
    let mut prev_size: usize = 0;
    let mut group_sizes_buffer: Vec<u8> = Vec::new();
    let mut group_encoder: Option<Encoder<&mut Vec<u8>>> = None;
    let mut buffer_encoded_seq: Vec<Vec<u8>> = Vec::new();

    for record_result in fa_reader.records() {
        let record = record_result?;
        let header = record.id();
        let ids_str = header.strip_prefix("ids:").unwrap();
        let ids: Vec<usize> = ids_str
            .split(',')
            .filter(|s| !s.is_empty())
            .map(|s| s.parse::<usize>().unwrap() - 1)
            .collect();

        let seq = record.seq();
        let encoded_seq = <Converter as Convert<&[u8]>>::str2num(seq);
        let size = seq.len();

        if let Some(ref current) = current_ids {
            if &ids != current {
                for elem in &buffer_encoded_seq {
                    omni_file.write_all(elem)?;
                    prev_tigs_size += elem.len() as u64;
                }
                buffer_encoded_seq.clear();

                if let Some(encoder) = group_encoder.take() {
                    encoder.finish()?;
                }
                prev_bucket_pos += (8 + group_sizes_buffer.len()) as u64;
                pos_nb_unitig.push((prev_tigs_size, prev_bucket_pos));
                size_file.write_all(&(group_sizes_buffer.len() as u64).to_le_bytes())?;
                size_file.write_all(&group_sizes_buffer)?;

                for id in current {
                    id_to_color_vec[*id].push(cid);
                }
                cid += 1;
                prev_size = 0;
                group_sizes_buffer.clear();
                current_ids = Some(ids);
                group_encoder = Some(Encoder::new(&mut group_sizes_buffer, 1)?);
            }
        } else {
            current_ids = Some(ids);
            group_encoder = Some(Encoder::new(&mut group_sizes_buffer, 4)?);
        }

        buffer_encoded_seq.push(encoded_seq);

        if buffer_encoded_seq.len() >= 1000 {
            for elem in &buffer_encoded_seq {
                omni_file.write_all(elem)?;
                prev_tigs_size += elem.len() as u64;
            }
            buffer_encoded_seq.clear();
        }

        let delta = size - prev_size;
        if let Some(ref mut encoder) = group_encoder {
            encoder.write_all(&delta.to_le_bytes())?;
        }
        prev_size = size;
    }

    if let Some(current) = current_ids {
        for elem in &buffer_encoded_seq {
            omni_file.write_all(elem)?;
            prev_tigs_size += elem.len() as u64;
        }
        buffer_encoded_seq.clear();

        if let Some(encoder) = group_encoder.take() {
            encoder.finish()?;
        }
        prev_bucket_pos += (8 + group_sizes_buffer.len()) as u64;
        pos_nb_unitig.push((prev_tigs_size, prev_bucket_pos));
        size_file.write_all(&(group_sizes_buffer.len() as u64).to_le_bytes())?;
        size_file.write_all(&group_sizes_buffer)?;

        for id in current {
            id_to_color_vec[id].push(cid);
        }
    }

    omni_file.flush()?;
    size_file.flush()?;

    println!(
        "Completed compression: total tigs={}, total sizes={}",
        prev_tigs_size, prev_bucket_pos
    );
    Ok((pos_nb_unitig, id_to_color_vec))
}

/// Write position pairs (tigs cursor, sizes cursor) on disk as raw bytes.
///
/// Each entry is exactly 16 bytes: [u64 tigs_pos LE][u64 sizes_pos LE].
/// No compression — the entire file for 10M buckets is only ~160 MB.
/// Entry i is at byte offset i*16.
///
/// PARAM
/// - `pos_nb_unitigs`: vector of (tigs_cursor, sizes_cursor) pairs.
/// - `filepath`: path to position file (should end in .bin, not .zst).
///
/// RETURNS
/// - Vec<usize> of starting byte offsets in positions file (simply i*16).
fn write_positions(pos_nb_unitigs: Vec<(u64, u64)>, filepath: String) -> Result<Vec<usize>> {
    let mut pos_file = BufWriter::new(File::create(filepath).expect("unable to create file"));
    let mut vec_cursor_position = Vec::with_capacity(pos_nb_unitigs.len());
    for (i, elem) in pos_nb_unitigs.iter().enumerate() {
        vec_cursor_position.push(i * 16);
        pos_file.write_all(&elem.0.to_le_bytes())?;
        pos_file.write_all(&elem.1.to_le_bytes())?;
    }
    pos_file.flush()?;
    Ok(vec_cursor_position)
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
    cursor_positions: &[usize],
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
            loop {
                let mut raw = [0u8; 8];
                match std::io::Read::read_exact(&mut reader, &mut raw) {
                    Ok(()) => {
                        let cid = u64::from_le_bytes(raw) as usize;
                        if cid >= cursor_positions.len() {
                            return Err(io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!(
                                    "color id {} outside [0, {}) while writing {}",
                                    cid,
                                    cursor_positions.len(),
                                    cid_file_path
                                ),
                            ));
                        }
                        if !first {
                            cid_encoder.write_all(b",")?;
                        }
                        first = false;
                        write!(&mut cid_encoder, "{}", cursor_positions[cid])?;
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

fn write_id_to_color_id(
    cid_file_path: String,
    id_to_color_vec: Vec<Vec<usize>>,
    cursor_positions: Vec<usize>,
) -> std::io::Result<Vec<usize>> {
    let mut cid_file = BufWriter::new(File::create(&cid_file_path)?);
    let mut id_cid_line_sizes = Vec::with_capacity(id_to_color_vec.len());
    let mut tot_size = 0;

    for elem in id_to_color_vec {
        let mut to_write = String::new();

        for (i, e) in elem.iter().enumerate() {
            let pos = cursor_positions[*e];
            if i > 0 {
                to_write.push(',');
            }
            write!(&mut to_write, "{}", pos).unwrap();
        }

        let mut buffer = Vec::new();
        {
            let mut cid_encoder = Encoder::new(&mut buffer, 1)?;
            cid_encoder.write_all(to_write.as_bytes())?;
            cid_encoder.finish()?;
        }

        id_cid_line_sizes.push(tot_size);
        tot_size += 8 + buffer.len();
        cid_file.write_all(&(buffer.len() as u64).to_le_bytes())?;
        cid_file.write_all(&buffer)?;
    }

    cid_file.write_all(&(0_u64).to_le_bytes())?;
    Ok(id_cid_line_sizes)
}
