use clap::error::Result;
use rayon::str::Chars;
use seq_io::fasta::{Reader, Record};
use std::ops::Div;
use std::path::{Path, PathBuf};
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, Read, Seek, Write};
use::rayon::prelude::*;
use zstd::bulk::Decompressor;
use crate::utils::vec2str;

/// DECOMPRESS
/// Omnicolor: file containing omnicolored kmers
/// Multicolor: file containing multicolored kmers
/// Input names: fof containing genomes to extract. SAME ORDER AS FOF USED TO CREATE COMPRESSED FILES
/// 
/// DECOMPRESS will dump simplitigs from the requested files. 
/// Creating one output file per genome requested.
/// Omnicolored kmers are dumped in every out file.
/// Multicolored kmers are split according to their colors
pub fn decompress(omnicolor: &str, multicolor: &str, input_names: &str, out_dir: PathBuf, wanted_files_path: &str){
    //TODO DELETE OUTFILES IF THEY ALREADY EXISTS BEFORE DECOMPRESSING.
    //TODO CLEAN CODE
    
    let input_fof: File;   
    if wanted_files_path != ""{
        input_fof = File::open(wanted_files_path).unwrap();
    }else{
        //TODO ASK CONFIRMATION TO DECOMPRESS EVERYTHING.
        input_fof = File::open(input_names).unwrap();
    }
    let reader = BufReader::new(input_fof);
    let filenames: Vec<_> = reader.lines().collect::<Result<_, _>>().unwrap();
    for filename in filenames.iter(){
        let path = out_dir.join(String::from("Dump_")+ Path::new(filename).file_name().unwrap().to_str().unwrap());
        if path.is_file(){
            std::fs::remove_file(path).expect("Unable to remove file");
        }
    }
    decompress_multicolor("multicolor_bucket_size.txt", "filename_to_color.txt", wanted_files_path, multicolor, &out_dir);
   
    (0..filenames.len()).into_par_iter().for_each(|file_number|{
        let filename = filenames.get(file_number).unwrap();
        println!("{}", filename);
        let mut curr_file = File::open(omnicolor).unwrap();
        let mut curr_reader = BufReader::new(&curr_file);
        //let ( reader, _compression) = niffler::get_reader(Box::new(File::open(omnicolor).unwrap())).unwrap();
        //let mut counter = 0;
        let mut cursor = 0;

        let metadata = curr_file.metadata().unwrap();
        let mut file_size: usize = metadata.len() as usize;
        //println!("FILE SIZE = {}", file_size);
        while cursor < file_size{
            let mut size_buf = [0; 4];
            cursor += 4;
            curr_reader.read_exact(&mut size_buf).expect("Error reading simplitig size in temp file");
            let size_to_read: u32 = u32::from_le_bytes(size_buf).div_ceil(4);
            let size_simplitig: u32 = u32::from_le_bytes(size_buf);
            println!("READING CURSOR = {}", cursor);
            cursor += size_to_read as usize;
            println!("SIZE: {}", size_simplitig);
            let mut simplitig = vec![0; size_to_read as usize];
            println!("Reading {} Bytes", simplitig.len());
            curr_reader.read_exact(&mut simplitig).expect("Error reading simplitig");
            let to_write = vec2str(&simplitig, &(size_simplitig as usize));
            //counter += to_write.len()-31+1;
            let mut content = to_write;
            println!("CURR PATH: {}", filename);
            write_output(&content, &filename, &out_dir);
        }
        //println!("During decompression, I have seen {} k-mers", counter);
        //dump_file.finish().expect("Error writing decompressed data");
    }); 
}

fn decompress_multicolor(color_size_path: &str, filename_color_path: &str, wanted_files_path: &str, multicolor: &str, out_dir: &PathBuf){

    //Filename to color: PATH/TO/FILE.fa:0
    //Number = position of file in color array (e.g. 011001).
    let filename_color_file = File::open(out_dir.clone().join(filename_color_path)).unwrap();
    let filename_color_reader = BufReader::new(filename_color_file);
    let mut filename_to_color: Vec<_> = filename_color_reader.lines().collect::<Result<_, _>>().unwrap();
    //Color to encoded size of bucket for this color.
    //Color is the array (e.g. 011001) The size is the size in bytes to be read for this specific bucket.
    let color_size_file = File::open(out_dir.clone().join(color_size_path)).unwrap();
    let color_size_reader = BufReader::new(color_size_file);
    let color_to_pos: Vec<_> = color_size_reader.lines().collect::<Result<_, _>>().unwrap();
    let mut filenames = Vec::new();
    let mut positions_in_color = Vec::new();
    if wanted_files_path == ""{
        println!("No query provided, decompressing everything.");
        for line in filename_to_color.iter(){
            filenames.push(String::from(line.split(':').collect::<Vec<_>>()[0]));
            positions_in_color.push(line.split(':').collect::<Vec<_>>()[1].parse::<usize>().unwrap());
        }
    }else{
        (filenames, positions_in_color) = filter_filenames_multicolor(wanted_files_path, &mut filename_to_color);
    }
    decompress_needed(&color_to_pos, multicolor, out_dir, &filenames, &positions_in_color);
    
}

fn read_at_pos(multicolor_reader: &mut BufReader<File>, size: usize, prev_cursor: &mut u64) -> String{
    //let mut decompressor = Decompressor::new().unwrap();
    multicolor_reader.seek(std::io::SeekFrom::Start(*prev_cursor)).expect("Unable to seek to specified position");
    let mut buffer = vec![0; size.div_ceil(4)];
    multicolor_reader.read_exact(&mut buffer).unwrap();
    println!("Cursor is at: {}", prev_cursor);
    println!("Size to read: {}", size);
    println!("Buffer size: {}", buffer.len());
    println!("I have read {} bytes this time", multicolor_reader.stream_position().unwrap()- *prev_cursor);
    if vec2str(&buffer.to_vec(), &(size)).contains("AAAAAAAAAA"){
        println!("CURSOR AFTER READ: {}", *prev_cursor + size.div_ceil(4) as u64);
        println!("SIMPLITIG: {}", vec2str(&buffer.to_vec(), &size));
        let mut input = String::new();
        io::stdin().read_line(&mut input).expect("error: unable to read user input");
    }
    //let mut input = String::new();
    //std::io::stdin().read_line(&mut input).expect("error: unable to read user input");
    //let decompressed = decompressor.decompress(&buffer, size*100).unwrap();
    //println!("Decompressed: {}", vec2str(&buffer, &size));
    *prev_cursor += size.div_ceil(4) as u64;
    vec2str(&buffer, &size)
}

fn increment_cursor(cursor: &mut u64, end_pos: &String, to_read: bool, sizes: &Vec<&str>){
    println!("CURSOR BEFORE INCREMENT: {}", cursor);
    let mut already_read = 0;
    if !to_read{
        println!("SKIPPING COLOR");
        *cursor += end_pos.parse::<u64>().unwrap();
    }else{
        /*for size in sizes.iter(){
            already_read += size.parse::<u64>().unwrap().div_ceil(4);
        }*/
        *cursor += end_pos.parse::<u64>().unwrap() - *cursor;
    }
    println!("end_pos: {}", end_pos);
    println!("parsed end_pos: {}", end_pos.parse::<u64>().unwrap());
    println!("already read = {}", already_read);
    println!("CURSOR AFTER INCREMENT: {}", cursor);
}

fn filter_filenames_multicolor(wanted_files_path: &str, filename_to_color: &mut Vec<String>) -> (Vec<String>, Vec<usize>){
    let wanted_files = File::open(wanted_files_path).unwrap();
    let wanted_files_reader = BufReader::new(wanted_files);
    let wanted_files_list: Vec<_> = wanted_files_reader.lines().collect::<Result<_, _>>().unwrap();
    println!("Decompressing {} files...", wanted_files_list.len());
    //FILTER FILES AND COLORS OF INTEREST.
    let mut filenames = Vec::new();
    let mut positions_in_color = Vec::new();
    for elem in filename_to_color.iter(){
        let filename = elem.split(':').collect::<Vec<_>>()[0];
        for line in wanted_files_list.iter(){
            if line == filename{
                println!("{} is in wanted files, added to decompression list.", line);
                println!("Pos in color is: {}", elem.split(':').collect::<Vec<_>>()[1].parse::<usize>().unwrap());
                filenames.push(String::from(filename));
                positions_in_color.push(elem.split(':').collect::<Vec<_>>()[1].parse::<usize>().unwrap());
            }
        }
    }
    //ERROR HANDLING
    println!("NB file found: {}", filename_to_color.len());
    if filenames.len() != wanted_files_list.len(){
        println!("ERROR: I found a different number of files than the number asked.");
        println!("Maybe you asked files that where not in the initial input ?");
        println!("Files in wanted list:");
        for elem in wanted_files_list{
            println!("{}", elem);
        }
        println!("Files found:");
        for elem in &filenames{
            println!("{}", elem)
        }
        //TODO Ask if should continue decompressing.
    }
    (filenames, positions_in_color)
}

fn decompress_needed(color_to_pos: &Vec<String>, multicolor: &str, out_dir: &PathBuf, filenames: &Vec<String>, positions_in_color: &Vec<usize>){
    let mut prev_cursor: u64 = 0;
    let mut content: String = String::new();
    let color_to_sizes = organise_interface_data(color_to_pos);
    //OPEN MULTICOLOR FILE
    let multicolor_file = File::open(multicolor).unwrap();
    let mut multicolor_reader = BufReader::new(multicolor_file);

    //let mut seen_color = false;
    let mut to_read = false;
    let mut pos = 0;
    for (color, sizes, end_cursor_pos) in color_to_sizes.iter(){
        for i in positions_in_color.iter(){
            if color.chars().nth(*i).unwrap() == '1'{
                to_read = true;
                pos = *i as usize;
            }
        }
        if to_read{
            if color == "01101100"{
                println!("color: {}", color);
                println!("pos: {}", pos);
            
            }
            for size in sizes.iter(){
                content = read_at_pos(&mut multicolor_reader, size.parse::<usize>().unwrap(), &mut prev_cursor);
                if color == "01101100"{
                    println!("Reading {} bytes", size);
                    println!("Cursor then: {}", prev_cursor-size.parse::<u64>().unwrap());
                    println!("cursor now: {}", prev_cursor);
                    println!("current size: {}", size);
                    println!("Decompressing color: {}", color);
                    let mut input = String::new();
                    io::stdin().read_line(&mut input).expect("error: unable to read user input");
                }
                
                let pos = positions_in_color.iter().position(|&r| r == pos).unwrap();
                write_output(&content, filenames.get(pos).unwrap(), out_dir);
            }
        }
        println!("Color is: {}", color);
        increment_cursor(&mut prev_cursor, end_cursor_pos, to_read, sizes);
        to_read = false;
    }
}

//TRAITE LES DONNES DU FICHIER D'INTERFACE ET LES ORGANISE POUR LES RENDRE UTILISABLE EFFICACEMENT.
fn organise_interface_data(color_to_pos: &Vec<String>) -> Vec<(String, Vec<&str>, String)>{
    let mut color_to_sizes: Vec<(String, Vec<&str>, String)> = Vec::new();
    for color_size in color_to_pos{
        //01000000,32984:11349,24210,13602,4395,11693,4809,9139,...
        let first_part = color_size.split(':').collect::<Vec<_>>()[0];
        let sizes = color_size.split(':').collect::<Vec<_>>()[1];
        let color = first_part.split(',').collect::<Vec<_>>()[0];
        let end_cursor_pos = first_part.split(',').collect::<Vec<_>>()[1];
        println!("COLOR: {}", color);
        if sizes.contains(','){
            let sizes_vec = sizes.split(',').collect::<Vec<_>>();
            color_to_sizes.push((String::from(color), sizes_vec, String::from(end_cursor_pos)));
        }else{
            let mut vec = Vec::new();
            vec.push(sizes);
            color_to_sizes.push((String::from(color), vec, String::from(end_cursor_pos)));
        }
    }
    color_to_sizes
}

fn write_output(content: &String, filename: &str, out_dir: &PathBuf){
    let path = out_dir.join(String::from("Dump_")+ Path::new(filename).file_name().unwrap().to_str().unwrap());
    /*else{
        let mut file = File::options().write(true).read(true).create_new(true).open(path);
    }*/
    let mut out_file = File::options().write(true).read(true).create(true).open(&path).expect("Unable to create file");
    println!("OUTFILENAME: {}", path.display());
    out_file.seek(std::io::SeekFrom::End(0)).expect("unable to seek to end of file.");
    out_file.write_all(&">\n".as_bytes()).unwrap();
    out_file.write_all(content.as_bytes()).unwrap();
    out_file.write_all(&"\n".as_bytes()).unwrap();
}