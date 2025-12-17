use core::panic;
use std::collections::HashMap;
use std::process::Command;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Result, Seek, Write};
use std::path::Path;
use std::u32;

use num_traits::ToPrimitive;
use zstd::{Decoder, Encoder};

use crate::utils::{Converter, Convert, vec2str};

const SIZE_PAIR_POS: usize = 17;

//   =========================================================================================== DECOMPRESSION ==============================================================================

/// High-level decompression entry point.
///
/// Depending on whether `wanted_files_path` is provided, the function either
/// decompresses the entire archive or only the subset of files listed in the
/// `wanted_files_path`. It reads auxiliary files (positions, id->cid mapping,
/// tigs and sizes) and invokes the corresponding decompression routine.
///
/// PARAM
/// - `size_filename`: name of the sizes file in the input directory.
/// - `color_id_filename`: name of the id->color mapping file.
/// - `tigs_filename`: name of the tigs (encoded unitigs) file.
/// - `positions_filename`: name of the positions file.
/// - `filename_id`: filename that maps original file paths to ids (created at build time).
/// - `out_dir`: output directory where decompressed FASTA shards will be written.
/// - `wanted_files_path`: optional path to a file listing input files to extract (one per line).
/// - `input_dir`: directory where compressed inputs reside (contains sizes/tigs/positions).
///
/// RETURNS
/// - Result<()> indicating success or IO error.
pub fn decompress(size_filename: &String, color_id_filename: &String, tigs_filename: &String, positions_filename: &String, filename_id: &String, out_dir: &String, wanted_files_path: &String, input_dir: String) -> std::io::Result<()>{
    println!("Writting decompressed data in {out_dir}");

    

    if wanted_files_path != ""{
        let input_file = File::open(input_dir.clone() + filename_id).unwrap();
        let input_reader = BufReader::new(input_file);
        let mut filenames_id_map = HashMap::new();
        let mut file_id: u32 = 0;
        for line_result in input_reader.lines(){
            let line = line_result?;
            if let Some((path, size)) = line.split_once(':'){
                filenames_id_map.insert(path.to_string(), (file_id, size.parse::<u64>().unwrap()));
                println!("File: {}, ID: {}, LINE SIZE IN ID TO CIDS FILE: {}, LINE SIZE IN ID TO CIDS FILE AS U64: {}", path.to_string(), file_id, size, size.parse::<u64>().unwrap());
            }
            file_id += 1;
        }

        let cid_map_out_filenames = match get_cid_to_id_targeted(&(input_dir.clone() + &color_id_filename), &filenames_id_map, &wanted_files_path){
            Ok(map) => map,
            Err(e) => panic!("error gathering color set ids: {e:?}"),
        };
        let cid_to_id_map = cid_map_out_filenames.0;
        let wanted_filenames = cid_map_out_filenames.1;
        
        println!("Query file given, decompressing only subpart of archive....");
        decompress_wanted(&wanted_filenames, &(input_dir.clone()+positions_filename), cid_to_id_map, &(input_dir.to_owned()+tigs_filename), &(input_dir.to_owned()+size_filename), out_dir);
    }else{
        println!("No query file given, decompressing entire archive....");
        let cid_to_id_map = match get_cid_to_id(&(input_dir.clone() + &color_id_filename)){
            Ok(map) => map,
            Err(e) => panic!("Error getting cid to id map {e:?}"),
        };
        for elem in &cid_to_id_map{
            println!("CID: {}", elem.0);
        }
        let input_file = File::open(input_dir.clone() + filename_id).unwrap();
        let input_reader = BufReader::new(input_file);
        let mut filenames_id = Vec::new();
        let mut file_id: u32 = 0;
        //ID IS POS IN ID TO CID FILE
        for line_result in input_reader.lines(){
            let line = line_result?;
            if let Some((path, _)) = line.split_once(":"){
                ("File: {}, ID: {}", path, file_id);
                filenames_id.push((path.to_owned(), file_id));
            }
            file_id += 1;
        }
        println!("{}", input_dir.to_owned()+size_filename);
        decompress_all(&(input_dir.to_owned()+size_filename), &(input_dir.clone()+positions_filename), &(input_dir.to_owned()+tigs_filename), out_dir, filenames_id, cid_to_id_map);
    }
    Ok(())
}

/// Read and return all position pairs from the positions file.
///
/// The positions file stores many individually zstd-compressed records.
/// This function reads and decompresses `cid_to_id_map.len()` records and
/// returns a vector of (tigs_cursor, sizes_cursor) pairs.
///
/// PARAM
/// - `position_filename`: path to the positions file (zstd records).
/// - `cid_to_id_map`: used only to know how many records to read.
///
/// RETURNS
/// - Vec<(u32,u32)> list of (tigs_byte_offset, sizes_byte_offset).
fn get_positions(position_filename: &String, cid_to_id_map: &HashMap<usize, Vec<u32>>) -> Result<Vec<(u32, u32)>>{
    let position_file = File::open(position_filename)?;
    let mut positions_reader = BufReader::new(position_file);
    let mut positions: Vec<(u32, u32)> = Vec::new();
    for _i in 0..cid_to_id_map.len(){

        let mut size_buffer = [0; 4];
        positions_reader.read_exact(&mut size_buffer)?;
        let compressed_size = u32::from_le_bytes(size_buffer) as usize;

        let mut current_buffer = vec![0;compressed_size];
        positions_reader.read_exact(&mut current_buffer)?;
        let mut decompressed_positions = Vec::new();

        {
            let mut decoder_positions = Decoder::new(&current_buffer[..])?;
            decoder_positions.read_to_end(&mut decompressed_positions)?;
            println!("a{}a", decompressed_positions.len());
            let pos_size = u32::from_le_bytes(decompressed_positions[..4].try_into().unwrap());
            let pos_tigs = u32::from_le_bytes(decompressed_positions[4..].try_into().unwrap());
            println!("SIZE: {}", pos_size);
            println!("TIGS {}", pos_tigs);
            positions.push((pos_tigs, pos_size));
        }
    }
    Ok(positions)
}

/// Read individual position records for a targeted set of color ids.
///
/// This function seeks into the positions file to the exact byte offsets
/// corresponding to color ids present in `cid_to_id_map` and returns the
/// decompressed pairs.
///
/// PARAM
/// - `position_filename`: path to the positions file.
/// - `cid_to_id_map`: keys are color ids required; the function reads only those.
///
/// RETURNS
/// - Vec<(u32,u32)> for the requested color ids.
fn get_specific_pos(position_filename: &String, cid_to_id_map: &HashMap<usize, Vec<u32>>) -> Result<Vec<(u32, u32)>>{
    let position_file = File::open(position_filename)?;
    let mut positions_reader = BufReader::new(position_file);
    let mut positions: Vec<(u32, u32)> = Vec::new();
    let mut buffer_position_tigs = [0; 4];
    let mut buffer_position_size = [0; 4];
    for elem in cid_to_id_map.keys(){
        let cursor_pos : u64 = (elem+SIZE_PAIR_POS) as u64;
        positions_reader.seek(std::io::SeekFrom::Start(cursor_pos as u64))?;
        positions_reader.read_exact(&mut buffer_position_tigs)?;
        positions_reader.read_exact(&mut buffer_position_size)?;
        let mut decompressed_pos_tigs = Vec::new();
        let mut decompressed_pos_size = Vec::new();
        {
            let mut decoder_pos_tigs = Decoder::new(&buffer_position_tigs[..])?;
            decoder_pos_tigs.read_to_end(&mut decompressed_pos_tigs)?;
            let mut decoder_pos_size = Decoder::new(&buffer_position_size[..])?;
            decoder_pos_size.read_to_end(&mut decompressed_pos_size)?;
            let pos_size = u32::from_le_bytes(decompressed_pos_size.try_into().unwrap());
            let pos_tigs = u32::from_le_bytes(decompressed_pos_tigs.try_into().unwrap());
            println!("SIZE: {}", pos_size);
            println!("TIGS {}", pos_tigs);
            positions.push((pos_tigs, pos_size));
        }
    }
    Ok(positions)
}

/// Decompress the entire archive to individual FASTA shards.
///
/// Iterates over all color buckets and writes each unitig to every file that
/// belongs to the color-set. The function uses the `positions` vector to seek
/// in the size/tigs files appropriately.
///
/// PARAM
/// - `size_filename`: path to decompressed sizes file.
/// - `positions`: vector of (tigs_cursor,sizes_cursor) pairs for all buckets.
/// - `tigs_filename`: path to the tigs (encoded) file.
/// - `out_dir`: output directory for resulting FASTA fragments.
/// - `filenames`: Vec of (original_path, file_id) used to name outputs.
/// - `cid_to_id_map`: mapping color_id -> list of file ids that contain it.
fn decompress_all(size_filename: &String, positions_filename: &String, tigs_filename: &String, out_dir: &String, filenames: Vec<(String, u32)>, cid_to_id_map: HashMap<usize, Vec<u32>>) {
    let mut tigs_file = BufReader::new(
        File::open(&tigs_filename).expect("Error opening tigs file")
    );
    
    println!("nb cid: {}", cid_to_id_map.len());
    
    let mut sorted_cids: Vec<_> = cid_to_id_map.keys().cloned().collect();
    sorted_cids.sort();
    
    for (idx, cid) in sorted_cids.iter().enumerate() {
        let file_ids = cid_to_id_map.get(cid).unwrap();
        for elem in file_ids.clone(){
            println!("{elem}");
        }
        println!("Processing CID: {} (position in positions file)", cid);
        
        let (tigs_pos, sizes_pos) = match read_position_at_cid(positions_filename, *cid) {
            Ok(pos) => pos,
            Err(e) => {
                eprintln!("Error reading position at CID {}: {:?}", cid, e);
                continue;
            }
        };
        
        let (_next_tigs_pos, next_sizes_pos) = if idx + 1 < sorted_cids.len() {
            match read_position_at_cid(positions_filename, sorted_cids[idx + 1]) {
                Ok(pos) => pos,
                Err(_) => {
                    get_file_end_positions(tigs_filename, size_filename).unwrap()
                }
            }
        } else {
            get_file_end_positions(tigs_filename, size_filename).unwrap()
        };
        println!("{}", next_sizes_pos);
        let sizes = match read_bucket_sizes_at_position(size_filename, sizes_pos) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("Error reading sizes for CID {}: {:?}", cid, e);
                continue;
            }
        };
        
        if sizes.is_empty() {
            println!("No unitigs in bucket at CID {}", cid);
            continue;
        }
        
        tigs_file.seek(std::io::SeekFrom::Start(tigs_pos as u64)).expect("Failed to seek in tigs file");
        
        for size in sizes {
            if size < 31 {
                eprintln!("Warning: unitig size {} is less than k-mer size", size);
                continue;
            }
            
            let read_size = size.div_ceil(4);
            let mut tig_buffer = vec![0; read_size];
            tigs_file.read_exact(&mut tig_buffer).expect("Failed to read tig");
            
            let tig = vec2str(&tig_buffer, &size);
            for file_id in file_ids {
                println!("{file_id}");
                let curr_filename = &filenames[*file_id as usize];
                let trunc_filename = Path::new(&curr_filename.0).file_stem().unwrap();
                
                let output_path = format!("{}Dump_{}.fa", out_dir, trunc_filename.to_str().unwrap());
                
                let mut out_file = BufWriter::new(File::options().append(true).create(true).open(output_path).expect("Unable to create file"));
                writeln!(out_file, ">").unwrap();
                writeln!(out_file, "{}", tig).unwrap();
                out_file.flush().unwrap();
            }
        }
    }
}


/// Read full id->color_id file and build color id -> list of file ids.
///
/// The function expects the id->color file format written by write_id_to_color_id:
/// size(u64) + zstd_compressed_csv_positions repeated, terminated by 0_u64.
///
/// PARAM
/// - `color_id_filename`: path to the id->color file.
///
/// RETURNS
/// - HashMap<color_id, Vec<file_id>>
fn get_cid_to_id(color_id_filename: &String) -> Result<HashMap<usize, Vec<u32>>>{
    let mut color_id_file = BufReader::new(File::open(color_id_filename).expect("Error opening color id file, are you sur you gave the right path ?"));
    let mut cid_ids_map = HashMap::new();
    let mut counter: u32 = 0;
    println!("READING CID TO ID FILE: {}", color_id_filename);
    let mut buffer_size = [0; 8];
    color_id_file.read_exact(&mut buffer_size)?;
    let mut size_read = usize::from_le_bytes(buffer_size);
    while size_read != 0 {
        let mut buffer = vec![0; size_read];
        color_id_file.read_exact(&mut buffer)?;
        let mut decompressed_data = Vec::new();
        {
            let mut decoder = Decoder::new(&buffer[..])?;
            decoder.read_to_end(&mut decompressed_data)?;
        }
        let str_tmp = String::from_utf8(decompressed_data).expect("Error reading cids");
        let temp_cids = str_tmp.split(',').collect::<Vec<_>>();
        for cid in temp_cids{
            if cid != "" {
                cid_ids_map.entry(cid.parse::<usize>().unwrap())
                    .and_modify(|list: &mut Vec<_>| list.push(counter))
                    .or_insert(Vec::from([counter]));
            }
        }
        color_id_file.read_exact(&mut buffer_size)?;
        size_read = usize::from_le_bytes(buffer_size);
        counter += 1;
    }
    
    Ok(cid_ids_map)
}

/// Build cid -> file id mapping for a targeted subset of files.
///
/// This function reads the id->color file at the offsets provided by
/// `filenames_id_map` for the requested file paths listed in `wanted_files_path`.
/// It returns:
/// - a HashMap from color id -> list of file ids found in the requested set,
/// - a Vec of (filename, file_id) for the requested files (used by downstream code).
///
/// PARAM
/// - `color_id_filename`: path to id->color file.
/// - `filenames_id_map`: map original path -> (file_id, byte_offset_in_id_file).
/// - `wanted_files_path`: path containing one filename per line to extract.
///
/// RETURNS
/// - Result<(HashMap<usize, Vec<u32>>, Vec<(String,u32)>)>
fn get_cid_to_id_targeted(color_id_filename: &String, filenames_id_map: &HashMap<String, (u32, u64)>, wanted_files_path: &String) -> std::io::Result<(HashMap<usize, Vec<u32>>, Vec<(String, u32)>)>{    
    // FILE ID TO COLOR IDS
    let mut color_id_file = BufReader::new(File::open(color_id_filename).expect("Error opening color id file, are you sur you gave the right path ?"));

    // CID TO ID HASHMAP
    let mut cid_ids_map = HashMap::new();

    // WANTED FILENAMES VEC
    let mut wanted_filenames = Vec::new();

    // QUERIED FILES TO DECOMPRESS
    let wanted_file = File::open(wanted_files_path)?;
    let wanted_reader = BufReader::new(wanted_file);
    for line_result in wanted_reader.lines(){
        let line = line_result?;
        if filenames_id_map.contains_key(&line){
            println!("COUCOUC");
            let entry = filenames_id_map.get(&line).unwrap();
            println!("{}", entry.1);
            println!("{}", entry.0);
            color_id_file.seek(std::io::SeekFrom::Start(entry.1))?;
            let mut buffer_size = [0; 8];
            color_id_file.read_exact(&mut buffer_size)?;
            let size_read = usize::from_le_bytes(buffer_size);
            println!("{size_read}");
            let mut buffer = vec![0; size_read];
            color_id_file.read_exact(&mut buffer)?;
            let mut decompressed_data = Vec::new();
            {
                let mut decoder = Decoder::new(&buffer[..])?;
                decoder.read_to_end(&mut decompressed_data)?;
            }
            let str_tmp = String::from_utf8(decompressed_data).expect("Error reading cids");
            let temp_cids = str_tmp.split(',').collect::<Vec<_>>();
            for cid in temp_cids{
                if cid != "" {
                    cid_ids_map.entry(cid.parse::<usize>().unwrap())
                        .and_modify(|list: &mut Vec<_>| list.push(entry.0))
                        .or_insert(Vec::from([entry.0]));
                }
            }
            println!("REEEEE");
            wanted_filenames.push((line.clone(), entry.0));
        }else {
            println!("FILE {} NOT FOUND IN ARCHIVE, CHECK SPELLING OR ACTUAL PRESENCE IN ARCHIVE", line);
        }
    }
    Ok((cid_ids_map, wanted_filenames))
}

/// Convenience helper to fully decompress the bucket sizes file to a temporary text file.
///
/// This function is a small utility used in debugging and expects the input
/// `size_filename` to be zstd-compressed. It writes decompressed sizes to
/// `tmp_sizes.txt`.
fn decompress_sizes(size_filename: &String){
    let size_file = BufReader::new(File::open(size_filename).expect("Error opening size file, are you sure path is good?"));
    let mut size_decoder = Decoder::new(size_file).expect("Failed to decode color file");
    let mut tmp_size = BufWriter::new(File::create("tmp_sizes.txt").expect("Unable to create temporary size file"));
    io::copy(&mut size_decoder, &mut tmp_size);
}


fn decompress_wanted(wanted_files: &Vec<(String, u32)>, positions_filename: &String, cid_to_id_map: HashMap<usize, Vec<u32>>, tigs_filename: &String, size_filename: &String, out_dir: &String) {
    
    println!("NB COLOR TO DECOMPRESS: {}", cid_to_id_map.len());
    println!("Wanted files:");
    for elem in wanted_files {
        println!("{} : {}", elem.0, elem.1);
    }
    
    let mut tigs_file = BufReader::new(
        File::open(&tigs_filename).expect("Error opening tigs file")
    );
    
    let mut sorted_cids: Vec<_> = cid_to_id_map.keys().cloned().collect();
    sorted_cids.sort();
    
    for (_idx, cid) in sorted_cids.iter().enumerate() {
        let file_ids = cid_to_id_map.get(cid).unwrap();
        
        println!("Processing CID: {} (byte position in positions file)", cid);
        
        let (tigs_pos, sizes_pos) = match read_position_at_cid(positions_filename, *cid) {
            Ok(pos) => pos,
            Err(e) => {
                eprintln!("Error reading position at CID {}: {:?}", cid, e);
                continue;
            }
        };
        
        let sizes = match read_bucket_sizes_at_position(size_filename, sizes_pos) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("Error reading sizes for CID {}: {:?}", cid, e);
                continue;
            }
        };
        
        if sizes.is_empty() {
            println!("No unitigs at CID {}", cid);
            continue;
        }
        tigs_file.seek(std::io::SeekFrom::Start(tigs_pos as u64)).expect("Failed to seek in tigs file");
        
        for size in sizes {
            if size < 31 {
                eprintln!("Warning: unitig size {} is less than k-mer size", size);
                continue;
            }
            
            let read_size = size.div_ceil(4);
            let mut tig_buffer = vec![0; read_size];
            tigs_file.read_exact(&mut tig_buffer).expect("Failed to read tig");
            
            let tig = vec2str(&tig_buffer, &size);
            
            for wanted_file in wanted_files {
                if file_ids.contains(&wanted_file.1) {
                    let trunc_filename = Path::new(&wanted_file.0).file_stem().unwrap();
                    
                    let output_path = format!(
                        "{}Dump_{}.fa",
                        out_dir,
                        trunc_filename.to_str().unwrap()
                    );
                    
                    let mut out_file = BufWriter::new(
                        File::options().append(true).create(true).open(output_path).expect("Unable to create file")
                    );
                    
                    writeln!(out_file, ">").unwrap();
                    writeln!(out_file, "{}", tig).unwrap();
                    out_file.flush().unwrap();
                }
            }
        }
    }
}


fn read_position_at_cid(positions_filename: &str, cid: usize) -> Result<(u32, u32)> {
    let mut positions_file = BufReader::new(File::open(positions_filename)?);
    
    positions_file.seek(std::io::SeekFrom::Start(cid as u64))?;
    
    let mut size_buffer = [0; 4];
    positions_file.read_exact(&mut size_buffer)?;
    let compressed_size = u32::from_le_bytes(size_buffer) as usize;
    
    let mut compressed_buffer = vec![0; compressed_size];
    positions_file.read_exact(&mut compressed_buffer)?;
    
    let mut decompressed = Vec::new();
    let mut decoder = Decoder::new(&compressed_buffer[..])?;
    decoder.read_to_end(&mut decompressed)?;
    
    let pos_tigs = u32::from_le_bytes(decompressed[..4].try_into().unwrap());
    let pos_sizes = u32::from_le_bytes(decompressed[4..8].try_into().unwrap());
    
    Ok((pos_tigs, pos_sizes))
}


fn read_bucket_sizes_at_position(size_file_path: &str, start_pos: u32) -> Result<Vec<usize>> {
    let mut size_file = BufReader::new(File::open(size_file_path)?);
    size_file.seek(std::io::SeekFrom::Start(start_pos as u64))?;
    
    let mut size_buffer = [0; 8];
    size_file.read_exact(&mut size_buffer)?;
    let compressed_size = usize::from_le_bytes(size_buffer);
    println!("compressed size: {compressed_size}");
    println!("start pos:{start_pos}");
    let mut compressed_buffer = vec![0; compressed_size];
    size_file.read_exact(&mut compressed_buffer)?;
    

    let mut decoder = Decoder::new(&compressed_buffer[..])?;
    let mut decompressed = Vec::new();
    decoder.read_to_end(&mut decompressed)?;
    println!("aa");
    let mut sizes = Vec::new();
    let mut prev = 0;
    for chunk in decompressed.chunks_exact(8) {
        let delta = usize::from_le_bytes(chunk.try_into().unwrap());
        let actual_size = delta + prev;
        println!("actual:{actual_size}");
        println!("encoded:{delta}");
        println!("prev:{prev}");
        sizes.push(actual_size);
        prev = actual_size;
    }
    
    Ok(sizes)
}

fn get_file_end_positions(tigs_filename: &str, sizes_filename: &str) -> Result<(u32, u32)> {
    let tigs_file = File::open(tigs_filename)?;
    let tigs_size = tigs_file.metadata()?.len() as u32;
    
    let sizes_file = File::open(sizes_filename)?;
    let sizes_size = sizes_file.metadata()?.len() as u32;
    
    Ok((tigs_size, sizes_size))
}
