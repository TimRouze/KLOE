
use core::panic;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Result, Seek, Write};
use std::path::{Path, PathBuf};
use std::u32;
use bio::io::fasta;
use chrono::{DateTime, Duration, Utc};

use num_traits::ToPrimitive;
use zstd::stream::read::Decoder;
use zstd::Encoder;

use crate::utils::{Converter, Convert, vec2str};
use crate::parser;

pub fn compress(output_dir: &String, input_fof: &String, threads: usize, k: usize, m: usize, partition_power: u32, compaction_threads: usize) -> Result<()>{
    let overall_start = Utc::now();
    parser::run_parser(PathBuf::from(input_fof), PathBuf::from(output_dir), k, m, 10_u32, threads, compaction_threads, false, false);
    let parsing_time = Utc::now();
    println!("Simplitigs created, processing sequences");
    let mut input_fof_reader = BufReader::new(File::open(input_fof).expect("unable to open fof"));
    let mut filename = String::new();
    let mut filenames = Vec::new();
    while input_fof_reader.read_line(&mut filename)? != 0{
        filename.pop();
        filenames.push(filename.clone());
        filename.clear();

    }
    println!("Sorting sequences by color bucket");
    let sort_time = Utc::now();
    parser::log_checkpoint("Wall time:", parsing_time);
    let id_cid_line_sizes = sort_by_bucket(&output_dir, filenames.len() as u32);
    parser::log_checkpoint("Sorting took:", sort_time);
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
    // let pair = match process_simplitigs(
    //                         &(output_dir.clone() + "simplitigs.fa.zst"),
    //                         &(output_dir.clone()+"tigs_kloe.fa"),
    //                         &(output_dir.clone() + "bucket_sizes.txt")){
    //     Ok(res_pair) => res_pair,
    //     Err(e) => panic!("Error writing compressed unitigs: {e:?}"),
    // };
    let pair = match write_compressed_test(output_dir.clone()+"tigs_kloe.fa", output_dir, nb_files){
        Ok(res_pair) => res_pair,
        Err(e) => panic!("Error writing compressed unitigs: {e:?}"),
    };
    parser::log_checkpoint("Writing compressed sequences wall time:", write_time);
    let position_time = Utc::now();
    let pos_nb_unitig = pair.0;
    let id_to_color_vec = pair.1;
    // POS NB UNITIGS: 
        // STARTING POSITIONS IN THE TIGS FILE FOR EACH COLOR BUCKET (ACTUALLY THE SIZE OF EACH COLOR BUCKET) + NB UNITIG IN EACH COLOR BUCKET
        // NB UNITIGS GIVES THE NUMBER OF SIZES == 64B * NB UNITIGS = POS OF COLOR BUCKET IN THE POS FILE

    println!("Starting to write positions");
    let cursor_positions = match write_positions(pos_nb_unitig, String::from(output_dir.clone()+"positions_kloe.txt.zst")){
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


fn write_compressed_test(unitigs_file_path: String, output_dir: &String, nb_files: u32) -> Result<(Vec<(u32, u32)>, Vec<Vec<usize>>)>{
    let mut omni_file = BufWriter::new(File::create(unitigs_file_path)?);
    let mut size_file = BufWriter::new(File::create(output_dir.clone() + "bucket_sizes.txt")?);

    let unitigs_file = File::open(output_dir.clone() + "simplitigs.fa.zst")?;
    let decoder = zstd::Decoder::new(unitigs_file)?;
    let reader: Box<dyn BufRead> = Box::new(BufReader::new(decoder));
    let fa_reader = fasta::Reader::from_bufread(reader);

    let mut id_to_color_vec: Vec<Vec<usize>> = vec![Vec::new(); nb_files as usize];
    let mut pos_nb_unitig: Vec<(u32, u32)> = vec![(0, 0)];
    
    let mut prev_tigs_size: u32 = 0;
    let mut prev_bucket_pos: u32 = 0;
    let mut cid = 0_usize;

    let mut current_ids: Option<Vec<usize>> = None;
    let mut prev_size: usize = 0;  // For delta encoding within a group
    let mut total_tigs_bytes: u32 = 0;  // Cumulative tigs bytes for current group
    let mut group_sizes_buffer: Vec<u8> = Vec::new();  // Zstd buffer for current group

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
        let encoded = <Converter as Convert<&[u8]>>::str2num(seq);
        let size = seq.len();

        // Check if group changed
        if let Some(ref current) = current_ids {
            if &ids != current {
                // Flush previous group to disk
                prev_bucket_pos += (8 + group_sizes_buffer.len()) as u32;
                pos_nb_unitig.push((prev_tigs_size, prev_bucket_pos));
                size_file.write_all(&(group_sizes_buffer.len() as u64).to_le_bytes())?;
                size_file.write_all(&group_sizes_buffer)?;

                // Reset for new group
                for id in current {
                    id_to_color_vec[*id].push(cid);
                }
                cid += 1;
                prev_size = 0;
                total_tigs_bytes = 0;
                group_sizes_buffer.clear();
                current_ids = Some(ids);
            }
        } else {
            current_ids = Some(ids);
        }

        // Write encoded sequence directly to file
        omni_file.write_all(&encoded)?;
        total_tigs_bytes += encoded.len() as u32;
        prev_tigs_size += encoded.len() as u32;

        // Delta encode size and write to buffer
        let delta = size - prev_size;
        {
            let mut encoder = Encoder::new(&mut group_sizes_buffer, 4)?;
            encoder.write_all(&delta.to_le_bytes())?;
            encoder.finish()?;
        }
        prev_size = size;
    }

    // Process final group
    if let Some(current) = current_ids {
        prev_bucket_pos += (8 + group_sizes_buffer.len()) as u32;
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


/*
WRITE COMPRESSED READS SIMPLITIGS DUMP FILE.
FOR EACH COLOR BUCKET
    - GET THE UNITIGS
    - SORT THEM BY SIZE
    - COMPRESS THEM
    - KEEP THE COMPRESSED SIZES
    - WRITE THE UNITIGS IN KLOE'S TIGS FILE
    - WRITE SORTED SIZES IN SIZES FILEZ
*/
/// Read unitigs, group them by color bucket then, for each bucket
/// - write tigs (2bit encoded) in tigs file,
/// - write unitigs sizes (delta encoded + zstd compressed).
/// Finally, return the recorded positions (tigs-size cursor pairs).
///
/// PARAM
/// - `unitigs_file_path`: path where the encoded tigs output will be written.
/// - `output_dir`: base directory where the Fulgor unitigs dump is located.
///
/// RETURNS
/// - Vec<(u32,u32)> containing pairs (tigs bucket pos in tigs file, tigs sizes positions in the sizes file) for each bucket.
fn write_compressed(unitigs_file_path: String, output_dir: &String, nb_files: u32) -> Result<(Vec<(u32, u32)>, Vec<Vec<usize>>)>{
    //let mut debug_tigs =  BufWriter::new(File::create("mini_unitigs.txt").expect("debug file error"));
    //let mut debug_sizes =  BufWriter::new(File::create("mini_sizes.txt").expect("debug file error"));
    //let mut debug_sizes_clear = BufWriter::new(File::create("mini_sizes_clear.txt").expect("debug file error"));

    let mut cpt_debug_bucket = 0;

    let mut omni_file = BufWriter::new(File::create(unitigs_file_path).expect("unable to create file"));

    let mut size_file = BufWriter::new(File::create(output_dir.clone()+"bucket_sizes.txt").expect("unable to create file"));
    //let mut sizes_encoder = Encoder::new(BufWriter::new(size_file), 12).expect("Failed to create zstd encoder");

    let mut nb_lines_dna: usize = 0;
    let mut nb_kmer = 0;
    let unitigs_file = File::open(output_dir.clone() + "simplitigs.fa.zst")?;
    let decoder = zstd::Decoder::new(unitigs_file)?;
    let reader: Box<dyn BufRead> = Box::new(BufReader::new(decoder));
    let fa_reader = fasta::Reader::from_bufread(reader);
    let mut curr_id = 0;
    //let mut to_write = Vec::new();
    let mut line = String::new();
    //let mut size_read = reader.read_line(&mut line);
    let mut curr_id_vec = Vec::new();

    let mut id_to_color_vec: Vec<Vec<_>> = vec![Vec::new(); nb_files as usize];
    //let mut ids = Vec::new();
    let mut cid = 0_usize;

    let mut pos_nb_unitig: Vec<(u32,u32)> = Vec::new();
    let mut unitigs_sizes_list: Vec<(Vec<u8>, usize)> = Vec::new();
    println!("READING SIMPLITIGS DUMP");

    let mut prev_tigs_size: usize = 0;
    let mut prev_bucket_pos: usize = 0;

    pos_nb_unitig.push((prev_tigs_size as u32, prev_bucket_pos as u32));
    for record in fa_reader.records(){//.read_line(&mut line).unwrap() != 0{
        let record = record?;
        let seq = record.seq();
        //println!("{}", record.id());
        let header = record.id();

        let ids_str = header.strip_prefix("ids:").unwrap();
        let ids: Vec<usize> = ids_str.split(',')
            .filter(|s| !s.is_empty())
            .map(|s| s.parse::<usize>().unwrap() -1) // Convert to 0-based
            .collect();

        /*for id in &ids{
            print!("{id}");
        }
        println!("");*/
        if curr_id_vec.is_empty(){
            curr_id_vec = ids.clone();
        }
        
        //let color_id: usize = line.split(":").collect::<Vec<_>>()[1].split(",").collect::<Vec<_>>().parse().unwrap();
        if ids != curr_id_vec{

            //println!("CHANGING COLOR, DUMPING CURRENT COLOR");
            unitigs_sizes_list.sort_by(|a, b| a.1.cmp(&b.1));
            // 2BIT / NUC (A = 00, T = 10, C = 01, G = 11)
            let mut prev: usize = 0;
            //let mut size_bucket_pos: usize = 0;
            let mut total_size = 0;

            let mut vec_sizes = Vec::new();
                
            for pair in &unitigs_sizes_list{

                omni_file.write_all(&pair.0)?;
                total_size += pair.0.len();
                let delta_encoded: usize = pair.1 - prev;
                vec_sizes.push(delta_encoded);
                //sizes_encoder.write_all(&delta_encoded.to_le_bytes());
                //size_bucket_pos += 8;
                prev = pair.1;
            }
            prev_tigs_size += total_size;
            unitigs_sizes_list.clear();


            let mut buffer = Vec::new();
            {
                let mut sizes_encoder = Encoder::new(&mut buffer, 1).expect("Failed to create zstd encoder");
                for elem in vec_sizes{
                    sizes_encoder.write_all(&elem.to_le_bytes())?;
                }
                sizes_encoder.finish()?;
            }
            prev_bucket_pos += 8 + buffer.len();

            pos_nb_unitig.push((prev_tigs_size as u32, prev_bucket_pos as u32));
            size_file.write_all(&buffer.len().to_le_bytes())?;
            size_file.write_all(&buffer)?;
            cpt_debug_bucket += 1;
            for elem in &curr_id_vec{
                id_to_color_vec[*elem].push(cid);
                //println!("Writting color id {cid} in {elem}");
            }
            curr_id_vec = ids.clone();
            //println!("{cid}");
            //let mut input = String::new();
            //std::io::stdin().read_line(&mut input).expect("error: unable to read user input");

            cid += 1;
            
        }
        nb_lines_dna += 1;
        nb_kmer += seq.len()-30;
        let converter = Converter;
        unitigs_sizes_list.push((<Converter as Convert<&[u8]>>::str2num(seq), seq.len()));

    }
    if !unitigs_sizes_list.is_empty(){
        let mut prev = 0;


        unitigs_sizes_list.sort_by(|a, b| a.1.cmp(&b.1));
        //println!("nb unitigs: {}", unitigs_sizes_list.len());

        //let mut size_bucket_pos: usize = 0;
        let mut total_size = 0;

        let mut vec_sizes = Vec::new();
        //let mut size_file_test = BufWriter::new(File::create(output_dir.clone()+"sizes_test.txt").expect("unable to create file"));
    
        for pair in &unitigs_sizes_list{
            omni_file.write_all(&pair.0);
            //println!("SIZE: {}", pair.1);
            //println!("PREV: {prev}");
            //println!("SIZE TIG ENCODED: {}", pair.0.len());
            total_size += pair.0.len();
            vec_sizes.push(pair.1 - prev);
            prev = pair.1;
            //size_bucket_pos += 8;
        }
        prev_tigs_size += total_size;
        let mut buffer = Vec::new();
        {
            let mut sizes_encoder = Encoder::new(&mut buffer, 1).expect("Failed to create zstd encoder");
            //size_file_test.write_all(String::from(vec_sizes.len().to_string()).as_bytes())?;
            for elem in vec_sizes{
                //println!("a{}a", elem);
                sizes_encoder.write_all(&elem.to_le_bytes())?;
                //size_file_test.write_all(String::from(elem.to_string()).as_bytes())?;
            }
            sizes_encoder.finish()?;
        }
        prev_bucket_pos += 8 + buffer.len();
        pos_nb_unitig.push((prev_tigs_size as u32, prev_bucket_pos as u32));
        size_file.write_all(&buffer.len().to_le_bytes())?;
        //println!("buffer len: {}", buffer.len());
        size_file.write_all(&buffer)?;

        unitigs_sizes_list.clear();
        for elem in curr_id_vec{
            id_to_color_vec[elem].push(cid);
            //println!("Writting color id {cid} in {elem}");
        }

        //println!("{cid}");
        //let mut input = String::new();
        //std::io::stdin().read_line(&mut input).expect("error: unable to read user input");
    }
    println!("SIZE TOTALE TIGS: {}\nSIZE TOTALE SIZES: {}", prev_tigs_size, prev_bucket_pos);
    println!("I HAVE SEEN {} LINES WITH DNA", nb_lines_dna);
    println!("I HAVE SEEN {} K-MERS", nb_kmer);

    let mut i = 0;
    /*for elem in id_to_color_vec.clone(){
        println!("CID = {i}");
        for id in elem{
            print!("{id}");
        }
        print!("\n");
        i += 1;
    }*/

    //println!("NB BUCKETS {cpt_debug_bucket}");
    /*for elem in &pos_nb_unitig{
        println!("SIZE: {}", elem.1);
    }*/
    //debug_sizes.flush();
    //debug_tigs.flush();

    omni_file.flush()?;
    size_file.flush()?;

    //pos_nb_unitig.push((prev_tigs_size as u32, prev_bucket_pos as u32));
    // STARTING POSITIONS IN THE TIGS FILE FOR EACH COLOR BUCKET (ACTUALLY THE SIZE OF EACH COLOR BUCKET) + NB UNITIG IN EACH COLOR BUCKET
    // NB UNITIGS GIVES THE NUMBER OF SIZES == 64B * NB UNITIGS = POS OF COLOR BUCKET IN THE POS FILE
    Ok((pos_nb_unitig, id_to_color_vec))
}


/// Write position pairs (tigs cursor, sizes cursor) on disk.
///
/// Positions are written as little endian bytes (u32) + zstd compressed
///
/// PARAM
/// - `pos_nb_unitigs`: vector of (tigs_cursor, sizes_cursor) pairs.
/// - `filepath`: path to position file.
///
/// RETURNS
/// - Vec<usize> of starting positions in positions file.
fn write_positions(pos_nb_unitigs: Vec<(u32, u32)>, filepath: String) -> Result<Vec<usize>>{
    let mut pos_file = BufWriter::new(File::create(filepath).expect("unable to create file"));
    let mut vec_cursor_position = Vec::new();
    let mut cursor_pos = 0_usize;
    for elem in pos_nb_unitigs{
        let mut buffer = Vec::new();
        {
            let mut pos_encoder = Encoder::new(&mut buffer, 1).expect("Failed to create zstd encoder");
            pos_encoder.write_all(&elem.0.to_le_bytes())?;
            pos_encoder.write_all(&elem.1.to_le_bytes())?;
            pos_encoder.finish()?;
        }
        vec_cursor_position.push(cursor_pos);
        pos_file.write_all(&(buffer.len() as u32).to_le_bytes())?;
        cursor_pos += 4 + buffer.len();
        pos_file.write_all(&buffer)?;
    }
    Ok(vec_cursor_position)
}


fn write_id_to_color_id_test(cid_file_path: String, id_to_color_vec: Vec<Vec<usize>>, cursor_positions: Vec<usize>) -> std::io::Result<Vec<usize>>{
    let mut cid_file = BufWriter::new(File::create(&cid_file_path)?);
    let mut id_cid_line_sizes = Vec::with_capacity(id_to_color_vec.len());
    let mut tot_size = 0;

    for elem in id_to_color_vec {
        // Pre-allocate string with estimated size
        let mut to_write = String::new();
        
        for (i, e) in elem.iter().enumerate() {
            let pos = cursor_positions[*e];
            if i > 0 {
                to_write.push(',');
            }
            // Use fmt::write instead of string concat
            use std::fmt::Write;
            write!(&mut to_write, "{}", pos).unwrap();
        }

        let mut buffer = Vec::new();
        {
            let mut cid_encoder = Encoder::new(&mut buffer, 12)?;
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

//TODO SORT CID AND WRITE BIT VECTOR (1 COLOR IS PRESENT, 0 OTHERWISE)
// THEN COMPRESS FURTHER BY WRITING DIFFERENCES BETWEEN VECS
//CHECK NOTES FOR WHAT TO DO ??


/// Writes cursors for cid buckets in position file.
/// One line = list of positions = cids, line N° = file id.
///
/// The mapping is written to `cid_file_path` in a sequence of records:
/// [u64 size_of_record][zstd_compressed_ascii_comma_separated_positions] ...
/// and terminated by a zero-size u64. The function also returns a vector of
/// file offsets (byte cursors) for each input file's record in the written file.
///
/// PARAM
/// - `cid_file_path`: destination filename for id->color positions.
/// - `id_to_color_vec`: mapping from file id -> list of color IDs (as produced by get_colors).
/// - `cursor_positions`: mapping from color id -> (byte) cursor position in positions file.
///
/// RETURNS
/// - Vec<usize> of offsets (byte cursors) pointing to the start of each written record.
fn write_id_to_color_id(cid_file_path: String, id_to_color_vec: Vec<Vec<usize>>, cursor_positions: Vec<usize>) -> std::io::Result<Vec<usize>>{
    println!("WRITTING ID TO COLOR ID FILE: {}", cid_file_path);
    let mut cid_file = BufWriter::new(File::create(&cid_file_path).expect("unable to create file"));
    let mut id_cid_line_sizes = Vec::new();
    let mut tot_size = 0;
    for elem in id_to_color_vec{
        let mut to_write = String::new();
        let mut i: u16 = 0;
        //let mut prev: u16 = 0;
        for e in elem{
            //to_write += &(e.to_u16().unwrap()-prev).to_string();
            //cid_encoder.write_all(&(e.to_u16().unwrap()-prev).to_le_bytes())?;
            let pos = cursor_positions.get(e).unwrap();
            //println!("{pos}");
            //println!("{e}");    
            if i != 0{
                to_write = to_write + "," + &pos.to_string();// - prev).to_string();
            }else {
                to_write = to_write + &pos.to_string();// - prev).to_string();
                i += 1;
            }
            //prev = e.to_u16().unwrap();
        }
        let mut buffer = Vec::new();
        {
            let mut cid_encoder = Encoder::new(&mut buffer, 12).expect("Failed to create zstd encoder");
            cid_encoder.write_all(to_write.as_bytes())?;
            cid_encoder.finish()?;
        }
        id_cid_line_sizes.push(tot_size);
        tot_size += 8 + buffer.len();
        cid_file.write_all(&buffer.len().to_le_bytes())?;
        cid_file.write_all(&buffer)?;
        //println!("{to_write}");
    }
    cid_file.write_all(&(0_u64).to_le_bytes())?;
    Ok(id_cid_line_sizes)
}