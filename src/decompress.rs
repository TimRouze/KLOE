use core::panic;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Result, Seek, Write};
use std::path::Path;

use zstd::Decoder;

use crate::utils::vec2str;


//   =========================================================================================== DECOMPRESSION ==============================================================================

/// High-level decompression entry point.
pub fn decompress(size_filename: &String, color_id_filename: &String, tigs_filename: &String, positions_filename: &String, filename_id: &String, out_dir: &String, wanted_files_path: &String, input_dir: String) -> std::io::Result<()>{
    println!("Writing decompressed data in {out_dir}");

    if wanted_files_path != ""{
        let input_file = File::open(input_dir.clone() + filename_id).unwrap();
        let input_reader = BufReader::new(input_file);
        let mut filenames_id_map = HashMap::new();
        let mut file_id: u32 = 0;
        for line_result in input_reader.lines(){
            let line = line_result?;
            if let Some((path, size)) = line.split_once(':'){
                filenames_id_map.insert(path.to_string(), (file_id, size.parse::<u64>().unwrap()));
                println!("File: {}, ID: {}, offset: {}", path, file_id, size);
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
        let input_file = File::open(input_dir.clone() + filename_id).unwrap();
        let input_reader = BufReader::new(input_file);
        let mut filenames_id = Vec::new();
        let mut file_id: u32 = 0;
        for line_result in input_reader.lines(){
            let line = line_result?;
            if let Some((path, _)) = line.split_once(":"){
                filenames_id.push((path.to_owned(), file_id));
            }
            file_id += 1;
        }
        println!("{}", input_dir.to_owned()+size_filename);
        decompress_all(&(input_dir.to_owned()+size_filename), &(input_dir.clone()+positions_filename), &(input_dir.to_owned()+tigs_filename), out_dir, filenames_id, cid_to_id_map);
    }
    Ok(())
}

/// Preload all positions from the raw binary positions file.
///
/// The positions file contains N entries of exactly 16 bytes each:
/// [u64 tigs_pos LE][u64 sizes_pos LE].
/// Returns a Vec where index i corresponds to the position pair at byte offset i*16.
fn preload_positions(positions_filename: &str) -> Result<Vec<(u64, u64)>> {
    let mut file = BufReader::new(File::open(positions_filename)?);
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
    println!("Preloaded {} position entries from {}", num_entries, positions_filename);
    Ok(positions)
}

/// Build a lookup from byte-offset (as stored in id_to_color_id) to position index.
///
/// The id_to_color_id file stores byte offsets into the positions file as CID keys.
/// With the raw format, byte_offset = index * 16, so index = byte_offset / 16.
fn byte_offset_to_position_index(byte_offset: usize) -> usize {
    byte_offset / 16
}

/// Preload all bucket sizes from the sizes file into a Vec indexed by byte offset.
///
/// The sizes file format: repeated [u64 compressed_len][compressed_data...].
/// Returns a HashMap from byte_offset -> Vec<usize> of actual sizes.
fn preload_sizes(size_filename: &str) -> Result<HashMap<u64, Vec<usize>>> {
    let mut file = BufReader::new(File::open(size_filename)?);
    let mut sizes_map = HashMap::new();
    let mut offset: u64 = 0;
    loop {
        let mut len_buf = [0u8; 8];
        match file.read_exact(&mut len_buf) {
            Ok(()) => {},
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => break,
            Err(e) => return Err(e),
        }
        let compressed_size = usize::from_le_bytes(len_buf);
        if compressed_size == 0 {
            break;
        }
        let mut compressed_buffer = vec![0; compressed_size];
        file.read_exact(&mut compressed_buffer)?;

        let mut decoder = Decoder::new(&compressed_buffer[..])?;
        let mut decompressed = Vec::new();
        decoder.read_to_end(&mut decompressed)?;

        let mut sizes = Vec::new();
        let mut prev = 0;
        for chunk in decompressed.chunks_exact(8) {
            let delta = usize::from_le_bytes(chunk.try_into().unwrap());
            let actual_size = delta + prev;
            sizes.push(actual_size);
            prev = actual_size;
        }

        sizes_map.insert(offset, sizes);
        offset += 8 + compressed_size as u64;
    }
    println!("Preloaded {} size buckets from {}", sizes_map.len(), size_filename);
    Ok(sizes_map)
}

/// Decompress the entire archive to individual FASTA shards.
///
/// Uses preloaded positions and sizes to avoid per-CID file reopens.
/// Keeps output file handles open in a HashMap to avoid per-unitig open/close.
fn decompress_all(size_filename: &String, positions_filename: &String, tigs_filename: &String, out_dir: &String, filenames: Vec<(String, u32)>, cid_to_id_map: HashMap<usize, Vec<u32>>) {
    let mut tigs_file = BufReader::new(
        File::open(&tigs_filename).expect("Error opening tigs file")
    );

    let total_cids = cid_to_id_map.len();
    println!("nb cid: {}", total_cids);

    // Preload all positions and sizes into memory
    let all_positions = preload_positions(positions_filename).expect("Failed to preload positions");
    let all_sizes = preload_sizes(size_filename).expect("Failed to preload sizes");

    // Keep output file handles open
    let mut writers: HashMap<u32, BufWriter<File>> = HashMap::new();

    let mut sorted_cids: Vec<_> = cid_to_id_map.keys().cloned().collect();
    sorted_cids.sort();

    let mut total_unitigs = 0_u64;
    for (idx, cid) in sorted_cids.iter().enumerate() {
        let file_ids = cid_to_id_map.get(cid).unwrap();

        // Progress logging
        if idx % 1000 == 0 || idx == total_cids - 1 {
            println!("Processing CID {}/{} ({} unitigs written so far)", idx + 1, total_cids, total_unitigs);
        }

        // Look up position from preloaded data
        let pos_index = byte_offset_to_position_index(*cid);
        if pos_index >= all_positions.len() {
            eprintln!("Error: position index {} out of range for CID {}", pos_index, cid);
            continue;
        }
        let (tigs_pos, sizes_pos) = all_positions[pos_index];

        // Look up sizes from preloaded data
        let sizes = match all_sizes.get(&sizes_pos) {
            Some(s) => s,
            None => {
                eprintln!("Error: no sizes found at offset {} for CID {}", sizes_pos, cid);
                continue;
            }
        };

        if sizes.is_empty() {
            continue;
        }

        tigs_file.seek(std::io::SeekFrom::Start(tigs_pos)).expect("Failed to seek in tigs file");

        for size in sizes {
            if *size < 31 {
                eprintln!("Warning: unitig size {} is less than k-mer size", size);
                continue;
            }

            let read_size = size.div_ceil(4);
            let mut tig_buffer = vec![0; read_size];
            tigs_file.read_exact(&mut tig_buffer).expect("Failed to read tig");

            let tig = vec2str(&tig_buffer, size);
            for file_id in file_ids {
                let writer = writers.entry(*file_id).or_insert_with(|| {
                    let curr_filename = &filenames[*file_id as usize];
                    let trunc_filename = Path::new(&curr_filename.0).file_stem().unwrap();
                    let output_path = format!("{}Dump_{}.fa", out_dir, trunc_filename.to_str().unwrap());
                    BufWriter::new(File::options().append(true).create(true).open(output_path).expect("Unable to create file"))
                });
                writeln!(writer, ">").unwrap();
                writeln!(writer, "{}", tig).unwrap();
            }
            total_unitigs += 1;
        }
    }

    // Flush all writers at the end
    for (_, mut writer) in writers {
        writer.flush().unwrap();
    }
    println!("Decompression complete: {} unitigs written across {} CIDs", total_unitigs, total_cids);
}


/// Read full id->color_id file and build color id -> list of file ids.
fn get_cid_to_id(color_id_filename: &String) -> Result<HashMap<usize, Vec<u32>>>{
    let mut color_id_file = BufReader::new(File::open(color_id_filename).expect("Error opening color id file, are you sure you gave the right path?"));
    let mut cid_ids_map = HashMap::new();
    let mut counter: u32 = 0;
    println!("Reading CID to ID file: {}", color_id_filename);
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
fn get_cid_to_id_targeted(color_id_filename: &String, filenames_id_map: &HashMap<String, (u32, u64)>, wanted_files_path: &String) -> std::io::Result<(HashMap<usize, Vec<u32>>, Vec<(String, u32)>)>{
    let mut color_id_file = BufReader::new(File::open(color_id_filename).expect("Error opening color id file, are you sure you gave the right path?"));
    let mut cid_ids_map = HashMap::new();
    let mut wanted_filenames = Vec::new();

    let wanted_file = File::open(wanted_files_path)?;
    let wanted_reader = BufReader::new(wanted_file);
    for line_result in wanted_reader.lines(){
        let line = line_result?;
        if filenames_id_map.contains_key(&line){
            let entry = filenames_id_map.get(&line).unwrap();
            color_id_file.seek(std::io::SeekFrom::Start(entry.1))?;
            let mut buffer_size = [0; 8];
            color_id_file.read_exact(&mut buffer_size)?;
            let size_read = usize::from_le_bytes(buffer_size);
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
            wanted_filenames.push((line.clone(), entry.0));
        }else {
            println!("FILE {} NOT FOUND IN ARCHIVE, CHECK SPELLING OR ACTUAL PRESENCE IN ARCHIVE", line);
        }
    }
    Ok((cid_ids_map, wanted_filenames))
}

/// Decompress only the wanted files from the archive.
///
/// Uses preloaded positions and sizes. Keeps output file handles open.
fn decompress_wanted(wanted_files: &Vec<(String, u32)>, positions_filename: &String, cid_to_id_map: HashMap<usize, Vec<u32>>, tigs_filename: &String, size_filename: &String, out_dir: &String) {

    let total_cids = cid_to_id_map.len();
    println!("NB COLOR TO DECOMPRESS: {}", total_cids);
    println!("Wanted files:");
    for elem in wanted_files {
        println!("{} : {}", elem.0, elem.1);
    }

    let mut tigs_file = BufReader::new(
        File::open(&tigs_filename).expect("Error opening tigs file")
    );

    // Preload all positions and sizes into memory
    let all_positions = preload_positions(positions_filename).expect("Failed to preload positions");
    let all_sizes = preload_sizes(size_filename).expect("Failed to preload sizes");

    // Build a set of wanted file_ids for quick lookup
    let wanted_ids: std::collections::HashSet<u32> = wanted_files.iter().map(|w| w.1).collect();

    // Keep output file handles open
    let mut writers: HashMap<u32, BufWriter<File>> = HashMap::new();
    // Pre-open all wanted output files
    for wanted_file in wanted_files {
        let trunc_filename = Path::new(&wanted_file.0).file_stem().unwrap();
        let output_path = format!("{}Dump_{}.fa", out_dir, trunc_filename.to_str().unwrap());
        let writer = BufWriter::new(
            File::options().append(true).create(true).open(output_path).expect("Unable to create file")
        );
        writers.insert(wanted_file.1, writer);
    }

    let mut sorted_cids: Vec<_> = cid_to_id_map.keys().cloned().collect();
    sorted_cids.sort();

    let mut total_unitigs = 0_u64;
    for (idx, cid) in sorted_cids.iter().enumerate() {
        let file_ids = cid_to_id_map.get(cid).unwrap();

        // Progress logging
        if idx % 1000 == 0 || idx == total_cids - 1 {
            println!("Processing CID {}/{} ({} unitigs written so far)", idx + 1, total_cids, total_unitigs);
        }

        // Look up position from preloaded data
        let pos_index = byte_offset_to_position_index(*cid);
        if pos_index >= all_positions.len() {
            eprintln!("Error: position index {} out of range for CID {}", pos_index, cid);
            continue;
        }
        let (tigs_pos, sizes_pos) = all_positions[pos_index];

        // Look up sizes from preloaded data
        let sizes = match all_sizes.get(&sizes_pos) {
            Some(s) => s,
            None => {
                eprintln!("Error: no sizes found at offset {} for CID {}", sizes_pos, cid);
                continue;
            }
        };

        if sizes.is_empty() {
            continue;
        }
        tigs_file.seek(std::io::SeekFrom::Start(tigs_pos)).expect("Failed to seek in tigs file");

        for size in sizes {
            if *size < 31 {
                eprintln!("Warning: unitig size {} is less than k-mer size", size);
                continue;
            }

            let read_size = size.div_ceil(4);
            let mut tig_buffer = vec![0; read_size];
            tigs_file.read_exact(&mut tig_buffer).expect("Failed to read tig");

            let tig = vec2str(&tig_buffer, size);

            for file_id in file_ids {
                if wanted_ids.contains(file_id) {
                    if let Some(writer) = writers.get_mut(file_id) {
                        writeln!(writer, ">").unwrap();
                        writeln!(writer, "{}", tig).unwrap();
                    }
                }
            }
            total_unitigs += 1;
        }
    }

    // Flush all writers at the end
    for (_, mut writer) in writers {
        writer.flush().unwrap();
    }
    println!("Decompression complete: {} unitigs written across {} CIDs", total_unitigs, total_cids);
}
