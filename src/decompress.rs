use clap::error::Result;
use std::path::{Path, PathBuf};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, Write, BufWriter, Stdin};
use zstd::stream::read::Decoder;
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
    //let filename_to_nb_kmer = HashMap::new();

    for filename in filenames.iter(){
        let stem_filename = Path::new(filename).file_stem().unwrap();
        let path = out_dir.join(String::from("Dump_")+ stem_filename.to_str().unwrap());
        println!("{}", path.to_str().unwrap());
        if path.is_file(){
            std::fs::remove_file(path).expect("Unable to remove file");
        }
    }
    let full_path = Path::new(omnicolor);

    let mut color_size_path = String::from("multicolor_bucket_size.txt.zst");
    let mut filename_color_path = String::from("filename_to_color.txt");
    println!("{}", full_path.display());
    if let Some(parent_path) = full_path.parent() {
        println!("a{}a", parent_path.display());
        if parent_path.to_str().unwrap() != ""{
            color_size_path = String::from(parent_path.to_str().unwrap())+"/multicolor_bucket_size.txt.zst";
            filename_color_path = String::from(parent_path.to_str().unwrap())+"/filename_to_color.txt";
        }
        println!("{}", color_size_path);
        println!("{}", filename_color_path);
    }
    decompress_multicolor(&color_size_path, &filename_color_path, wanted_files_path, multicolor, &out_dir);
   
    //(0..filenames.len()).into_par_iter().for_each(|file_number|{
        //let filename = filenames.get(file_number).unwrap();
        //println!("{}", filename);
    let omni_file = File::open(omnicolor).unwrap();
    let mut omni_reader = BufReader::new(&omni_file);
    //let ( reader, _compression) = niffler::get_reader(Box::new(File::open(omnicolor).unwrap())).unwrap();
    let mut cursor = 0;
    let metadata = omni_file.metadata().unwrap();
    let file_size: usize = metadata.len() as usize;
    let mut counter_kmer = 0;
    println!("FILE SIZE = {}", file_size);
    while cursor < file_size{
        let mut size_buf = [0; 4];
        cursor += 4;
        omni_reader.read_exact(&mut size_buf).expect("Error reading simplitig size in temp file");
        let size_to_read: u32 = u32::from_le_bytes(size_buf).div_ceil(4);
        let size_simplitig: u32 = u32::from_le_bytes(size_buf);
        //println!("READING CURSOR = {}", cursor);
        cursor += size_to_read as usize;
        //println!("SIZE: {}", size_simplitig);
        let mut simplitig = vec![0; size_to_read as usize];
        //println!("Reading {} Bytes", simplitig.len());
        /*if cursor > 900000{
            println!("SIZE: {}", size_simplitig);
            println!("Reading {} Bytes", simplitig.len());
            let to_write = vec2str(&simplitig, &(size_simplitig as usize));
            //println!("simplitig = {}", to_write);
            println!("Cursor = {}", cursor);
            println!("SIZE READ: {}", size_to_read);
            let mut input = String::new();
            std::io::stdin().read_line(&mut input).expect("error: unable to read user input"); 
        }*/
        omni_reader.read_exact(&mut simplitig).expect("Error reading simplitig");
        let to_write = vec2str(&simplitig, &(size_simplitig as usize));
        //println!("SIMPLITIG: {}", to_write);
        //let mut input = String::new();
        //std::io::stdin().read_line(&mut input).expect("error: unable to read user input"); 
        //println!("simplitig = {}", to_write);
        let content = to_write;
        counter_kmer += content.len()-30;
        //println!("CURR PATH: {}", filename);
        write_out_omni(&content, &filenames, &out_dir);
    }
    println!("NB KMER SEEN IN OMNI {}", counter_kmer);
        //dump_file.finish().expect("Error writing decompressed data");
    //}); 
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
    let color_size_decoder = Decoder::new(color_size_reader).unwrap();
    let decoder_reader = BufReader::new(color_size_decoder);
    let color_to_pos: Vec<_> = decoder_reader.lines().collect::<Result<Vec<String>, _>>().unwrap();
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
    //let mut input = String::new();
    //std::io::stdin().read_line(&mut input).expect("error: unable to read user input");
    //let decompressed = decompressor.decompress(&buffer, size*100).unwrap();
    //println!("Decompressed: {}", vec2str(&buffer, &size));
    *prev_cursor += size.div_ceil(4) as u64;
    vec2str(&buffer, &size)
}

fn increment_cursor(cursor: &mut u64, end_pos: &String, to_read: bool){
    //println!("CURSOR BEFORE INCREMENT: {}", cursor);
    if !to_read{
        *cursor = end_pos.parse::<u64>().unwrap();
    }else{
        *cursor += end_pos.parse::<u64>().unwrap() - *cursor;
    }
    //println!("cursor after increment: {}", cursor);
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
                filenames.push(String::from(filename));
                positions_in_color.push(elem.split(':').collect::<Vec<_>>()[1].parse::<usize>().unwrap());
            }
        }
    }
    //ERROR HANDLING
    println!("There are {} files in the compressed archive", filename_to_color.len());
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
    let color_to_sizes = organise_interface_data(color_to_pos);
    println!("NB FILES {}", filenames.len());
    //OPEN MULTICOLOR FILE
    let multicolor_file = File::open(multicolor).unwrap();
    let mut multicolor_reader = BufReader::new(multicolor_file);
    let mut to_read = false;
    let mut counter_kmer = 0;
    for (color, sizes, end_cursor_pos) in color_to_sizes.iter(){
        for i in positions_in_color.iter(){
            if color.chars().nth(*i).unwrap() == '1'{
                to_read = true;
                break;
            }
        }
        if to_read{
            for size in sizes.iter(){
                let content = read_at_pos(&mut multicolor_reader, size.parse::<usize>().unwrap(), &mut prev_cursor);
                for elem in positions_in_color{
                    /*println!("ELEM: {}", elem);
                    if *elem == 0{
                        println!("COLOR: {}", color)
                    }*/
                    if color.chars().nth(*elem).unwrap() == '1'{
                        //println!("{}", color);
                        /*if *elem == 0{
                            println!("POS IN COLOR: {}", elem);
                            println!("IS AT POS: {}", positions_in_color.iter().position(|pos| pos == elem).unwrap());
                            let mut input = String::new();
                            std::io::stdin().read_line(&mut input).expect("error: unable to read user input");
                        }*/
                        counter_kmer += content.len()-30;
                        write_output(&content, filenames.get(positions_in_color.iter().position(|pos| pos == elem).unwrap()).unwrap(), out_dir);
                    }
                    //let mut input = String::new();
                    //std::io::stdin().read_line(&mut input).expect("error: unable to read user input");
                    
                }
            }
        }
        //println!("filename: {}", filenames.get(pos).unwrap());
        //println!("NB KMER = {}", counter);
        increment_cursor(&mut prev_cursor, end_cursor_pos, to_read);
        //let mut input = String::new();
        //std::io::stdin().read_line(&mut input).expect("error: unable to read user input");
        to_read = false;
    }

    println!("IN MULT, I HAVE READ {} KMERS", counter_kmer);
}

//TRAITE LES DONNES DU FICHIER D'INTERFACE ET LES ORGANISE POUR LES RENDRE UTILISABLE EFFICACEMENT.
fn organise_interface_data(color_to_pos: &Vec<String>) -> Vec<(String, Vec<&str>, String)>{
    let mut color_to_sizes: Vec<(String, Vec<&str>, String)> = Vec::new();
    for color_size in color_to_pos{
        let first_part = color_size.split(':').collect::<Vec<_>>()[0];
        let sizes = color_size.split(':').collect::<Vec<_>>()[1];
        let color = first_part.split(',').collect::<Vec<_>>()[0];
        let end_cursor_pos = first_part.split(',').collect::<Vec<_>>()[1];
        //println!("COLOR: {}", color);
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

fn write_out_omni(content: &String, filenames: &Vec<String>, out_dir: &PathBuf){
    for filename in filenames{
        write_output(content, filename, out_dir);
    }
}

fn write_output(content: &String, filename: &str, out_dir: &PathBuf){
    let path = out_dir.join(String::from("Dump_")+ Path::new(filename).file_stem().unwrap().to_str().unwrap());
    /*else{
        let mut file = File::options().write(true).read(true).create_new(true).open(path);
    }*/
    //println!("BONJOURENT");
    let mut out_file = BufWriter::new(File::options().write(true).read(true).create(true).open(&path).expect("Unable to create file"));
    //println!("OUTFILENAME: {}", path.display());
    out_file.seek(std::io::SeekFrom::End(0)).expect("unable to seek to end of file.");
    out_file.write_all(&">\n".as_bytes()).unwrap();
    out_file.write_all(content.as_bytes()).unwrap();
    out_file.write_all(&"\n".as_bytes()).unwrap();    
}
