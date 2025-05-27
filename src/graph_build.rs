
use std::collections::HashMap;
use std::process::Command;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, Write, Result};
use std::path::Path;

use rayon::prelude::*;
use zstd::{Decoder, Encoder};

pub fn build_graphs(output_dir: &String, input_fof: &String, threads: &usize, tmp_dir: &String){

    println!("Building Fulgor index");

    let output = Command::new("sh")
        .arg("-c")
        .arg("./../fulgor/build/fulgor build --force -k 31 -m 19 -l ".to_owned() + input_fof + " -o " + output_dir + "fulgor_index_unitigs -t " + &threads.to_string() + " -d " + tmp_dir)
        .output()
        .expect("failed to execute process");
    println!("{}", String::from_utf8(output.stdout).unwrap());
    println!("{}", String::from_utf8(output.stderr).unwrap());
    println!("Fulgor index created, dumping");

    let mut input_fof = BufReader::new(File::open(input_fof).expect("unable to create file"));
    let mut input_list = String::new();
    let mut fof_id = BufWriter::new(File::create(output_dir.clone() + "filenames_id.txt").expect("Failed to create fof file"));
    input_fof.read_to_string(&mut input_list);
    let mut file_cpt: u32 = 0;
    for filename in input_list.lines(){
        println!("a{}a", filename);
        fof_id.write_all((filename.to_owned() + ":" + &file_cpt.to_string() + "\n").as_bytes());
        file_cpt += 1;
    }
    fof_id.flush();
    println!("./../fulgor/build/fulgor dump -i {} fulgor_index_unitigs.fur", output_dir);
    let output = Command::new("sh")
        .arg("-c")
        .arg("./../fulgor/build/fulgor dump -i ".to_owned() + output_dir + "fulgor_index_unitigs.fur")
        .output()
        .expect("failed to execute process");
    println!("{}", String::from_utf8(output.stdout).unwrap());
    println!("{}", String::from_utf8(output.stderr).unwrap());
    sort_by_bucket(&output_dir, file_cpt);
}

pub fn sort_by_bucket(output_dir: &String, nb_files: u32){
    let size_file = BufWriter::new(File::create(output_dir.clone()+"bucket_sizes.txt.zst").expect("unable to create file"));
    let mut zstd_sizes = Encoder::new(size_file, 12).expect("Failed to create zstd encoder");
    let id_to_color_vec: Vec<Vec<_>>;
    //let mut multi_file_vec: Vec<BufWriter<File>> = Vec::new();

    
    // GET COLORING INFORMATION FOR DECOMPRESSION PURPOSES
    id_to_color_vec = get_colors(String::from(output_dir.clone()+"fulgor_index_unitigs.color_sets.txt"), nb_files);
    // COMPRESS UNITIGS
    let sizes = write_compressed(output_dir.clone()+"tigs_kloe.fa", output_dir);
    // WRITE BUCKET SIZES
    for elem in sizes{
        zstd_sizes.write_all(&elem.to_le_bytes());
        //println!("{}", elem);
    }
    zstd_sizes.finish();
    // WRITE FILE ID TO COLOR ID FILE
    write_sizes(output_dir.clone()+"id_to_color_id.txt.zst", id_to_color_vec);
    
}

fn get_colors(color_file_path: String, nb_files: u32) -> Vec<Vec<usize>>{
    let color_file = File::open(color_file_path).unwrap();
    let color_reader = BufReader::new(color_file);
    let color_list: Vec<_> = color_reader.lines().collect::<Result<_>>().unwrap();
    let mut id_to_color_vec: Vec<Vec<_>> = vec![Vec::new(); nb_files as usize];

    println!("GATHERING COLORS");

    for line in color_list.iter(){
        let split_space_line = line.split(" ").collect::<Vec<_>>();
        let id: usize = split_space_line[0].split("=").collect::<Vec<_>>()[1].parse().unwrap();
        //println!("CURR COLOR ID: {}", id);
        let color_set_vec = &split_space_line[2..];
        for color in color_set_vec{
            //println!("{}", color);
            id_to_color_vec[color.parse::<usize>().unwrap()].push(id.clone());
        }
    }
    id_to_color_vec
}

fn write_compressed(unitigs_file_path: String, output_dir: &String) -> Vec<usize>{
    let mut omni_file = BufWriter::new(File::create(unitigs_file_path).expect("unable to create file"));
    let mut cpt = 0;
    let mut nb_kmer = 0;
    let unitigs_file = File::open(output_dir.clone() + "fulgor_index_unitigs.unitigs.fa").unwrap();
    let mut reader = BufReader::new(unitigs_file);
    let mut curr_id = 0;
    let mut to_write = String::new();
    let mut sizes = Vec::new();
    let mut line = String::new();
    let mut size_read = reader.read_line(&mut line);
    while size_read.unwrap() != 0{
        line.pop();
        if line.starts_with(">"){
            let color_id: usize = line.split(" ").collect::<Vec<_>>()[2].split("=").collect::<Vec<_>>()[1].parse().unwrap();
            if curr_id != color_id{
                curr_id = color_id;
                
                sizes.push(write_in_file(&mut omni_file, to_write.clone()));
                //omni_file.write_all(to_write.as_bytes());
                to_write.clear();
            }
        }else{
            cpt += 1;
            nb_kmer += line.len()-30;
            if to_write.is_empty(){
                to_write.push_str(&line);
            }else{
                to_write.push_str(&("\n".to_owned() + &line));
            }
        }
        line.clear();
        size_read = reader.read_line(&mut line);
    }
    /*while let Some(line) = line_iter.next(){
        let line = line.expect("Error reading record").clone();
        if line.starts_with(">"){
            let color_id: usize = line.split(" ").collect::<Vec<_>>()[2].split("=").collect::<Vec<_>>()[1].parse().unwrap();
            if curr_id != color_id{
                curr_id = color_id;
                
                sizes.push(write_in_file(&mut omni_file, to_write.clone()));
                //omni_file.write_all(to_write.as_bytes());
                to_write.clear();
            }
        }else{
            cpt += 1;
            nb_kmer += line.len()-30;
            if to_write.is_empty(){
                to_write.push_str(&line);
            }else{
                to_write.push_str(&("\n".to_owned() + &line));
            }
        }
    }*/
    if !to_write.is_empty(){
        sizes.push(write_in_file(&mut omni_file, to_write.clone()));
        to_write.clear();
    }
    println!("I HAVE SEEN {} LINES WITH DNA", cpt);
    println!("I HAVE SEEN {} K-MERS", nb_kmer);

    omni_file.flush();
    sizes
}

fn write_in_file(file: &mut BufWriter<File>, data: String) -> usize{
    let mut compressed = Vec::new();
    let mut zstd_encoder = Encoder::new(&mut compressed, 4).expect("Failed to create zstd encoder");
    zstd_encoder.write_all(data.as_bytes());
    zstd_encoder.finish();
    let curr_size = compressed.len();
    file.write_all(&compressed);
    curr_size
}

// ID TO CID -> COMPRESS FURTHER BY WRITING DIFFERENCES BETWEEN CID VECS

fn write_sizes(cid_file_path: String, id_to_color_vec: Vec<Vec<usize>>){
    let cid_file = BufWriter::new(File::create(cid_file_path).expect("unable to create file"));
    let mut cid_encoder = Encoder::new(cid_file, 12).expect("Failed to create zstd encoder");
    for elem in id_to_color_vec{
        let mut to_write = String::new();
        let mut i = 0;
        for e in elem{
            if i != 0{
                to_write = to_write + "," + &e.to_string();
            }else {
                to_write = to_write + &e.to_string();
                i += 1;
            }
        }
        to_write = to_write + "\n";
        cid_encoder.write_all(to_write.as_bytes());
    }
    cid_encoder.finish();
}

pub fn init_decompress(size_filename: String, color_id_filename: String, tigs_filename: String, out_dir: &String, wanted_files_path: &str, input_dir: &String) -> std::io::Result<()>{
    //INPUT DIR + LES FILENAMES NECESSAIRES
    //FIRST: FILENAME_ID
    //SECOND: ID_COLOR_ID
    //THIRD: BUCKET_SIZES
    //FOURTH: TIGS_KLOE
    let input_file = File::open(input_dir.to_owned()+"filenames_id.txt").unwrap();
    let reader = BufReader::new(input_file);
    let mut filenames = Vec::new();
    for line_result in reader.lines() {
        let line = line_result?;
        if let Some((path, _useless)) = line.split_once(':') {
            filenames.push(path.to_string());
        }
    }
    println!("OUTDIR: {}", out_dir);
    // ===================================== READ IDS TO COLOR IDS AND CREATE CID TO IDS MAP ============================================================== //
    let cid_to_id_map = get_cid_to_id(&(out_dir.clone() + &color_id_filename));
    
    // =================================================== READ COLOR IDS TO COLOR BUCKET SIZES ============================================================= //
    let sizes = get_cid_to_bucket_size(&(out_dir.clone() + &size_filename), cid_to_id_map.len());

    // ============================================================== FILTER NECESSARY IDS =================================================================== //

    if wanted_files_path != ""{
        decompress_wanted(wanted_files_path, filenames, cid_to_id_map, tigs_filename, &sizes, out_dir);
    }else{
        decompress_all(&sizes, cid_to_id_map, tigs_filename, out_dir, filenames);
    }
    Ok(())
    //TODO
    // build unitig graph for each output using ggcat

}

fn get_cid_to_id(color_id_filename: &String) -> HashMap<String, Vec<u32>>{
    let color_id_file = BufReader::new(File::open(color_id_filename).expect("Error opening color id file, are you sur you gave the right path ?"));
    let mut color_id_decoder = Decoder::new(color_id_file).expect("Failed to decode color file");
    let mut colors = String::new();
    color_id_decoder.read_to_string(&mut colors);
    let lines = colors.split("\n").collect::<Vec<_>>();
    let mut ids = HashMap::new();
    let mut counter: u32 = 0;
    for line in lines{
        let temp_cids = line.split(",").collect::<Vec<_>>();
        for cid in temp_cids{
            if cid != ""{
                ids.entry(String::from(cid))
                    .and_modify(|list: &mut Vec<_>| list.push(counter))
                    .or_insert(Vec::from([counter]));
            }
        }
        counter += 1;
    }
    ids
}

fn get_cid_to_bucket_size(size_filename: &String, nb_color: usize) -> Vec<usize>{
    let size_file = BufReader::new(File::open(size_filename).expect("Error opening size file, are you sure path is good?"));
    let mut size_decoder = Decoder::new(size_file).expect("Failed to decode color file");
    let mut sizes = Vec::new();
    for _i in 0..nb_color{
        let mut size = [0; 8];
        size_decoder.read_exact(&mut size);
        sizes.push(usize::from_le_bytes(size));
    }
    sizes
}

fn decompress_wanted(wanted_files_path: &str, input_filenames: Vec<String>, cid_to_id_map: HashMap<String, Vec<u32>>, tigs_filename: String, sizes: &Vec<usize>, out_dir: &String){
    let mut wanted_ids = get_to_decompress(wanted_files_path, input_filenames, &cid_to_id_map);
    //let mut tigs_decoder = Decoder::new(tigs_file).expect("Failed to decode tigs file");

    // CREATE/OPEN FILENAMES
    let mut dump_files = Vec::new();
    let mut input_fof = File::open(wanted_files_path).unwrap();
    let mut reader = BufReader::new(input_fof);
    let mut filenames = String::new();
    reader.read_to_string(&mut filenames);

    for line in filenames.lines(){
        let stem_filename = Path::new(line).file_stem().unwrap();
        println!("FILENAME: {}", (out_dir.to_owned() + "Dump_" + stem_filename.to_str().unwrap()));
        let mut out_file = BufWriter::new(File::create(out_dir.to_owned() + "Dump_" + stem_filename.to_str().unwrap()).expect("Unable to create file"));
        out_file.write_all(b">\n");
        dump_files.push(out_file);
    }


    //let mut concat_file = BufWriter::new(File::create(String::from("concat_wanted") + ".fa").expect("Unable to create file"));
    //concat_file.write_all(b">\n");
    decompress_needed(sizes, wanted_ids, tigs_filename, &mut dump_files);

    
    //concat_file.flush();
}

fn decompress_needed(sizes: &Vec<usize>, wanted_ids: HashMap<String, Vec<usize>>, tigs_filename: String, dump_files: &mut Vec<BufWriter<File>>){
    let mut counter: usize = 0;
    let mut cursor: usize = 0;
    let mut tigs_file = BufReader::new(File::open(&tigs_filename).expect("Error opening tigs file, are you sure path is good?"));
    
    for size in sizes{
        if let Some(out_files) = wanted_ids.get(&counter.to_string()){
            println!("COLOR BUCKET N°: {}", counter);
            println!("SIZE READ: {}", size);
            println!("OUT FILES:");
            for elem in out_files{
                println!("{}", elem);
            }
            //println!("SIZE: {}", String::from_utf8(clear_color_bucket.clone()).unwrap().len());
            //println!("SEQUENCES: {}", String::from_utf8(clear_color_bucket).unwrap());
            let mut bucket_buffer = vec![0; *size];
            tigs_file.seek(std::io::SeekFrom::Start(cursor.try_into().unwrap()));
            tigs_file.read_exact(&mut bucket_buffer);
            let mut tigs_decoder = Decoder::new(&bucket_buffer[..]).expect("Failed to create zstd decoder");
            let mut clear_color_bucket = Vec::new();
            tigs_decoder.read_to_end(&mut clear_color_bucket);
            let output = String::from_utf8(clear_color_bucket).unwrap();
            for out_file in dump_files.iter_mut(){
                for line in output.lines(){
                    out_file.write_all(line.as_bytes());
                    out_file.write_all(b"\n>\n");
                }
            }
        }          

        cursor += size;
        counter += 1;
    }
    for file in dump_files.iter_mut(){
        file.flush();
    }
}

fn get_to_decompress(wanted_files_path: &str, input_filenames: Vec<String>, cids_to_files: &HashMap<String, Vec<u32>>) -> HashMap<String, Vec<usize>>{
    let input_fof = File::open(wanted_files_path).unwrap();
    let mut reader = BufReader::new(input_fof);
    let mut filenames = String::new();
    reader.read_to_string(&mut filenames);
    
    
    let mut wanted_ids = Vec::new();
    println!("{}", filenames);
    for filename in input_filenames{
        println!("a{}a", filename);
        let curr_filename = filename.split(":").collect::<Vec<_>>()[0].parse::<String>().unwrap();
        if filenames.contains(&curr_filename){
            println!("b{}b", filename);
            let id = filename.split(":").collect::<Vec<_>>()[1].parse::<usize>().unwrap();
            wanted_ids.push(id);
        }
    }

    let mut to_decompress = HashMap::new();

    for (cid, files) in cids_to_files.iter(){
        for id in &wanted_ids{
            if files.contains(&(*id as u32)){
                to_decompress.entry(cid.clone())
                .and_modify(|list: &mut Vec<_>| list.push(*id))
                .or_insert(Vec::from([*id]));
            }
        }   
    }
    to_decompress
}

fn decompress_all(sizes: &Vec<usize>, cid_to_ids: HashMap<String, Vec<u32>>, tigs_filename: String, out_dir: &String, filenames: Vec<String>){
    let mut nb_kmers = 0;
    let mut concat_file = BufWriter::new(File::create(String::from("concat_all") + ".fa").expect("Unable to create file"));
    let mut tigs_file = BufReader::new(File::open(&tigs_filename).expect("Error opening tigs file, are you sure path is good?"));
    //let mut tigs_decoder = Decoder::new(tigs_file).expect("Failed to decode tigs file");
    let mut counter = 0;
    concat_file.write_all(b">\n");
    for size in sizes{
        let mut bucket_buffer = vec![0; *size];
        tigs_file.read_exact(&mut bucket_buffer);
        let mut tigs_decoder = Decoder::new(&bucket_buffer[..]).expect("Failed to create zstd decoder");
        let mut clear_color_bucket = Vec::new();
        tigs_decoder.read_to_end(&mut clear_color_bucket);
        println!("COLOR BUCKET N°: {}", counter);
        println!("SIZE READ: {}", size);
        //println!("SIZE: {}", String::from_utf8(clear_color_bucket.clone()).unwrap().len());
        //println!("SEQUENCES: {}", String::from_utf8(clear_color_bucket).unwrap());
        println!("COUNTER: a{}a", counter);
        let curr_output = cid_to_ids.get(&counter.to_string()).unwrap();
        for elem in curr_output{
            let curr_filename = filenames.get(*elem as usize).unwrap();
            println!("{}", curr_filename);
            let trunc_filename = Path::new(curr_filename).file_stem().unwrap();
            let mut out_file = BufWriter::new(File::options().append(true).create(true).open(out_dir.to_owned() + "Dump_" + trunc_filename.to_str().unwrap()).expect("Unable to create file"));
            out_file.write_all(&clear_color_bucket);
            out_file.flush();
        }
        let output = String::from_utf8(clear_color_bucket).unwrap();
        for line in output.lines(){
            nb_kmers += line.len() - 30;
            concat_file.write_all(&line.as_bytes());
            concat_file.write_all(b"\n>\n");
        }
        

        counter += 1;
    }
    concat_file.flush();
    println!("I HAVE SEEN {} K-MERS", nb_kmers);
    // TODO ENCODE OMNICOLOR AS NORMAL COLOR 

}
