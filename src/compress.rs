
use core::panic;
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Result, Write};
use std::path::PathBuf;
use std::sync::mpsc;
use chrono::Utc;

use zstd::Encoder;

use crate::utils::{Converter, Convert};
use crate::parser;

pub fn compress(output_dir: &String, input_fof: &String, threads: usize, k: usize, m: usize, partition_power: u32, verify_kmers: bool, skip_sort: bool, use_unitigs: bool, use_matchtigs: bool, use_eulertigs: bool, only_step1: bool, only_step2: bool, only_step3: bool) -> Result<()>{
    let _overall_start = Utc::now();

    // Start parser in streaming mode: Steps 1+2 run synchronously,
    // Step 3 (merge) runs in a background thread sending records via channel.
    let (record_rx, merge_handle, dataset_count) = parser::run_parser_streaming(
        PathBuf::from(input_fof),
        PathBuf::from(output_dir),
        k, m, partition_power, threads, verify_kmers, skip_sort,
        use_unitigs, use_matchtigs, use_eulertigs,
        only_step1, only_step2, only_step3,
    ).expect("run_parser_streaming failed");

    if only_step1 {
        // Step 1 only — partition files are on disk, nothing more to do.
        return Ok(());
    }

    println!("Simplitigs created, processing sequences (streaming)");
    let mut input_fof_reader = BufReader::new(File::open(input_fof).expect("unable to open fof"));
    let mut filename = String::new();
    let mut filenames = Vec::new();
    while input_fof_reader.read_line(&mut filename)? != 0{
        filename.pop();
        filenames.push(filename.clone());
        filename.clear();
    }

    println!("Compressing sequences from stream by color bucket");
    let sort_time = Utc::now();
    let id_cid_line_sizes = sort_by_bucket_streaming(output_dir, dataset_count as u32, record_rx);
    parser::log_checkpoint("Compression took:", sort_time);

    // Wait for merge thread to finish
    match merge_handle.join() {
        Ok(Ok(())) => {}
        Ok(Err(e)) => eprintln!("Merge thread error: {:#}", e),
        Err(_) => eprintln!("Merge thread panicked"),
    }

    let mut fof_id = BufWriter::new(File::create(output_dir.clone() + "filenames_id.txt").expect("Failed to create fof file"));
    let mut file_cpt: usize = 0;
    for filename in filenames{
        println!("a{}a", filename);
        fof_id.write_all((filename + ":" + id_cid_line_sizes.get(file_cpt).unwrap().to_string().as_str() + "\n").as_bytes())?;
        file_cpt += 1;
    }
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
pub fn sort_by_bucket(output_dir: &String, nb_files: u32) -> Vec<usize>{
    let write_time = Utc::now();
    // PROCESS AND COMPRESS UNITIGS
    println!("Starting writing compressed sequences.");
    let pair = match write_compressed(output_dir.clone()+"tigs_kloe.fa", output_dir, nb_files){
        Ok(res_pair) => res_pair,
        Err(e) => panic!("Error writing compressed unitigs: {e:?}"),
    };
    parser::log_checkpoint("Writing compressed sequences wall time:", write_time);
    let position_time = Utc::now();
    let pos_nb_unitig = pair.0;
    let id_to_color_vec = pair.1;

    println!("Starting to write positions");
    let cursor_positions = match write_positions(pos_nb_unitig, String::from(output_dir.clone()+"positions_kloe.bin")){
        Ok(vec) => vec,
        Err(e) => panic!("Error writting positions: {e:?}"),
    };
    parser::log_checkpoint("Write positions wall time:", position_time);
    let id_time = Utc::now();
    // WRITE FILE ID TO COLOR ID FILE
    let write_id_cid = match write_id_to_color_id_test(output_dir.clone()+"id_to_color_id.txt.zst", id_to_color_vec, cursor_positions){
        Ok(id_cid_line_sizes) => id_cid_line_sizes,
        Err(e) => panic!("error writting id to color id list: {e:?}"),
    };
    parser::log_checkpoint("Write id to cid wall time:", id_time);
    parser::log_checkpoint("Compression took:", write_time);
    write_id_cid
}

/// Streaming variant: consumes records from a channel instead of reading from a file.
fn sort_by_bucket_streaming(output_dir: &String, nb_files: u32, record_rx: mpsc::Receiver<parser::SimplitigRecord>) -> Vec<usize>{
    let write_time = Utc::now();
    println!("Starting writing compressed sequences (streaming).");
    let pair = match write_compressed_from_stream(output_dir.clone()+"tigs_kloe.fa", output_dir, nb_files, record_rx){
        Ok(res_pair) => res_pair,
        Err(e) => panic!("Error writing compressed unitigs: {e:?}"),
    };
    parser::log_checkpoint("Writing compressed sequences wall time:", write_time);
    let position_time = Utc::now();
    let pos_nb_unitig = pair.0;
    let id_to_color_vec = pair.1;

    println!("Starting to write positions");
    let cursor_positions = match write_positions(pos_nb_unitig, String::from(output_dir.clone()+"positions_kloe.bin")){
        Ok(vec) => vec,
        Err(e) => panic!("Error writting positions: {e:?}"),
    };
    parser::log_checkpoint("Write positions wall time:", position_time);
    let id_time = Utc::now();
    let write_id_cid = match write_id_to_color_id_test(output_dir.clone()+"id_to_color_id.txt.zst", id_to_color_vec, cursor_positions){
        Ok(id_cid_line_sizes) => id_cid_line_sizes,
        Err(e) => panic!("error writting id to color id list: {e:?}"),
    };
    parser::log_checkpoint("Write id to cid wall time:", id_time);
    parser::log_checkpoint("Compression took:", write_time);
    write_id_cid
}

/// Consume sorted SimplitigRecords from a channel and compress them.
/// Same logic as write_compressed but reads from mpsc::Receiver instead of FASTA file.
/// Uses bitset-based group detection (u64 compare) instead of header string parsing.
fn write_compressed_from_stream(
    unitigs_file_path: String,
    output_dir: &String,
    nb_files: u32,
    record_rx: mpsc::Receiver<parser::SimplitigRecord>,
) -> Result<(Vec<(u64, u64)>, Vec<Vec<usize>>)> {
    let mut omni_file = BufWriter::new(File::create(unitigs_file_path)?);
    let mut size_file = BufWriter::new(File::create(output_dir.clone() + "bucket_sizes.txt")?);

    let mut id_to_color_vec: Vec<Vec<usize>> = vec![Vec::new(); nb_files as usize];
    let mut pos_nb_unitig: Vec<(u64, u64)> = vec![(0, 0)];

    let mut prev_tigs_size: u64 = 0;
    let mut prev_bucket_pos: u64 = 0;
    let mut cid = 0_usize;

    let mut current_bitset: u64 = u64::MAX; // sentinel: no group yet
    let mut current_words: Vec<u64> = Vec::new();
    let mut has_group = false;
    let mut prev_size: usize = 0;
    let mut group_sizes_buffer: Vec<u8> = Vec::new();
    let mut group_encoder: Option<Encoder<&mut Vec<u8>>> = None;
    let mut buffer_encoded_seq: Vec<Vec<u8>> = Vec::new();

    for record in record_rx {
        let seq = &record.seq;
        let encoded_seq = <Converter as Convert<&[u8]>>::str2num(seq);
        let size = seq.len();

        // Detect group change via bitset comparison (u64 for ≤64 datasets)
        let group_changed = if !has_group {
            true
        } else if record.ids_words.is_empty() && current_words.is_empty() {
            record.ids_bitset != current_bitset
        } else {
            record.ids_bitset != current_bitset || record.ids_words != current_words
        };

        if group_changed {
            if has_group {
                // Flush previous group
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

                // Iterate bits directly into id_to_color_vec (no intermediate Vec)
                if current_words.is_empty() {
                    let mut word = current_bitset;
                    while word != 0 {
                        let bit = word.trailing_zeros() as usize;
                        id_to_color_vec[bit].push(cid);
                        word &= word - 1;
                    }
                } else {
                    for (word_idx, &word) in current_words.iter().enumerate() {
                        let mut w = word;
                        while w != 0 {
                            let bit = w.trailing_zeros() as usize;
                            id_to_color_vec[word_idx * 64 + bit].push(cid);
                            w &= w - 1;
                        }
                    }
                }
                cid += 1;
                prev_size = 0;
                group_sizes_buffer.clear();
                group_encoder = Some(Encoder::new(&mut group_sizes_buffer, 1)?);
            } else {
                group_encoder = Some(Encoder::new(&mut group_sizes_buffer, 4)?);
            }
            current_bitset = record.ids_bitset;
            current_words = record.ids_words.clone();
            has_group = true;
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

    if has_group {
        for elem in &buffer_encoded_seq{
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

        // Final group: iterate bits into id_to_color_vec
        if current_words.is_empty() {
            let mut word = current_bitset;
            while word != 0 {
                let bit = word.trailing_zeros() as usize;
                id_to_color_vec[bit].push(cid);
                word &= word - 1;
            }
        } else {
            for (word_idx, &word) in current_words.iter().enumerate() {
                let mut w = word;
                while w != 0 {
                    let bit = w.trailing_zeros() as usize;
                    id_to_color_vec[word_idx * 64 + bit].push(cid);
                    w &= w - 1;
                }
            }
        }
    }

    omni_file.flush()?;
    size_file.flush()?;

    println!("Completed compression: total tigs={}, total sizes={}", prev_tigs_size, prev_bucket_pos);
    Ok((pos_nb_unitig, id_to_color_vec))
}


/// Read unitigs from file, group them by color bucket then compress.
/// (Legacy path, kept for sort_by_bucket which is used by the non-streaming code path)
fn write_compressed(unitigs_file_path: String, output_dir: &String, nb_files: u32) -> Result<(Vec<(u64, u64)>, Vec<Vec<usize>>)>{
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
        for elem in &buffer_encoded_seq{
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

    println!("Completed compression: total tigs={}, total sizes={}", prev_tigs_size, prev_bucket_pos);
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
fn write_positions(pos_nb_unitigs: Vec<(u64, u64)>, filepath: String) -> Result<Vec<usize>>{
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


fn write_id_to_color_id_test(cid_file_path: String, id_to_color_vec: Vec<Vec<usize>>, cursor_positions: Vec<usize>) -> std::io::Result<Vec<usize>>{
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
