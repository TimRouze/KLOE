
use std::collections::HashMap;
use std::process::Command;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, Write, Result};
use std::path::Path;

use num_traits::ToPrimitive;
use zstd::{Decoder, Encoder};

use crate::utils::{str2num, vec2str};
//TODO DEBUG MODE -> WRITE SUBPART OF FILES (UNITIGS, COLOR ID TO SIZES ETC ETC) UNCOMPRESSED + READABLE TO HAVE A REF

//FILES: 
    // ID TO CID
    // POS FOR CID
    // SIZES
        // POUR CHAQUE BUCKET: POSITION DE DEPART DE LA COULEUR + LES TAILLES
    // TIGS


// TODO REDUCE + CHANGE NAME OF WRITE_COMPRESSED

//TODO ID TO CID DELTA ENCODE
// TODO SAVE TOTAL BUCKET SIZE IN CID TO POS TO ENABLE "RANDOM ACCESS"
    // COMPRESS UNITIGS SIZES BY COLOR BLOCK AND SAVE COMPRESSED SIZE OF SIZES COLOR BLOCK.
    // STORE POSITIONS IN SIZE FILE -> DELTA ENCODING
// TODO MONOCHROMATIGS SORTED PAR TAILLE AU SEIN D'UN BUCKET DE COULEUR
    // ENABLES DELTA ENCODING OF SIZES 
    // POUR SORT, BUCKET PAR BUCKET, EN 2BIT /NUC NORMALEMENT ENVIRON < 4GB
pub fn build_graphs(output_dir: &String, input_fof: &String, threads: &usize, tmp_dir: &String, memory: &usize){

    let output;
    println!("Building Fulgor index");
    if tmp_dir == ""{
        output = Command::new("sh")
            .arg("-c")
            .arg("\\time ./../fulgor/build/fulgor build --force -k 31 -m 19 -l ".to_owned() + input_fof + " -o " + output_dir + "fulgor_index_unitigs -t " + &threads.to_string() + " -g " + &memory.to_string())
            .output()
            .expect("failed to execute process");
    }else{
        output = Command::new("sh")
            .arg("-c")
            .arg("\\time ./../fulgor/build/fulgor build --force -k 31 -m 19 -l ".to_owned() + input_fof + " -o " + output_dir + "fulgor_index_unitigs -t " + &threads.to_string() + " -d " + tmp_dir + " -g " + &memory.to_string())
            .output()
            .expect("failed to execute process");
    }
    println!("Command was: {}", "./../fulgor/build/fulgor build --force -k 31 -m 19 -l ".to_owned() + input_fof + " -o " + output_dir + "fulgor_index_unitigs -t " + &threads.to_string() + " -d " + tmp_dir + " -g " + &memory.to_string());
    println!("{}", String::from_utf8(output.stdout).unwrap());
    println!("{}", String::from_utf8(output.stderr).unwrap());
    println!("Fulgor index created, dumping");

    let mut input_fof_reader = BufReader::new(File::open(input_fof).expect("unable to create file"));
    let mut input_list = String::new();
    let mut fof_id = BufWriter::new(File::create(output_dir.clone() + "filenames_id.txt").expect("Failed to create fof file"));
    input_fof_reader.read_to_string(&mut input_list);
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
        .arg("\\time ./../fulgor/build/fulgor dump -i ".to_owned() + output_dir + "fulgor_index_unitigs.fur")
        .output()
        .expect("failed to execute process");
    println!("{}", String::from_utf8(output.stdout).unwrap());
    println!("{}", String::from_utf8(output.stderr).unwrap());
    sort_by_bucket(&output_dir, file_cpt);
}

pub fn sort_by_bucket(output_dir: &String, nb_files: u32){
    let id_to_color_vec: Vec<Vec<_>>;
    //let mut multi_file_vec: Vec<BufWriter<File>> = Vec::new();

    
    // GET COLORING INFORMATION FOR DECOMPRESSION PURPOSES
    id_to_color_vec = get_colors(String::from(output_dir.clone()+"fulgor_index_unitigs.color_sets.txt"), nb_files);
    // PROCESS AND COMPRESS UNITIGS
    let mut pos_nb_unitig = write_compressed(output_dir.clone()+"tigs_kloe.fa", output_dir);
    // POS NB UNITIGS: 
        // STARTING POSITIONS IN THE TIGS FILE FOR EACH COLOR BUCKET (ACTUALLY THE SIZE OF EACH COLOR BUCKET) + NB UNITIG IN EACH COLOR BUCKET
        // NB UNITIGS GIVES THE NUMBER OF SIZES == 64B * NB UNITIGS = POS OF COLOR BUCKET IN THE POS FILE

    write_positions(pos_nb_unitig, String::from(output_dir.clone()+"positions_kloe.txt.zst"));

    // WRITE FILE ID TO COLOR ID FILE
    let write_id_cid = write_id_to_color_id(output_dir.clone()+"id_to_color_id.txt.zst", id_to_color_vec);
    match write_id_cid {
        Ok(()) => println!("Writting id to color id list: Finished"),
        Err(e) => println!("error writting id to color id list: {e:?}"),
    }
    
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

/*
WRITE COMPRESSED READS FULGOR'S DUMP FILE CONTAINING UNITIGS.
FOR EACH COLOR BUCKET
    - GET THE UNITIGS
    - SORT THEM BY SIZE
    - COMPRESS THEM
    - KEEP THE COMPRESSED SIZES
    - WRITE THE UNITIGS IN KLOE'S TIGS FILE
    - WRITE SORTED SIZES IN SIZES FILEZ
*/
fn write_compressed(unitigs_file_path: String, output_dir: &String) -> Vec<(usize, usize)>{
    let mut omni_file = BufWriter::new(File::create(unitigs_file_path).expect("unable to create file"));

    let size_file = BufWriter::new(File::create(output_dir.clone()+"bucket_sizes.txt.zst").expect("unable to create file"));
    let mut sizes_encoder = Encoder::new(size_file, 12).expect("Failed to create zstd encoder");
    
    let mut nb_lines_dna: usize = 0;
    let mut nb_kmer = 0;
    let unitigs_file = File::open(output_dir.clone() + "fulgor_index_unitigs.unitigs.fa").unwrap();
    let mut reader = BufReader::new(unitigs_file);
    let mut curr_id = 0;
    //let mut to_write = Vec::new();
    let mut line = String::new();
    //let mut size_read = reader.read_line(&mut line);

    let mut pos_nb_unitig: Vec<(usize,usize)> = Vec::new();
    let mut unitigs_sizes_list: Vec<(Vec<u8>, usize)> = Vec::new();
    println!("READING FULGOR DUMP");
    while reader.read_line(&mut line).unwrap() != 0{
    //while size_read.unwrap() != 0{
        line.pop();
        // println!("{}", line);
        // let mut input = String::new();
        // std::io::stdin().read_line(&mut input).expect("error: unable to read user input");
        if line.starts_with(">"){
            let color_id: usize = line.split(" ").collect::<Vec<_>>()[2].split("=").collect::<Vec<_>>()[1].parse().unwrap();
            if curr_id != color_id{
                //prev = 0;
                println!("CHANGING COLOR, DUMPING CURRENT COLOR");
                unitigs_sizes_list.sort_by(|a, b| a.1.cmp(&b.1));
                // println!("================================ COLOR ID {} ================================================", curr_id);
                // for e in &unitigs_list{ 
                //     println!("Size: {}", e.len());
                // }
                // 2BIT / NUC (A = 00, T = 10, C = 01, G = 11)
                let mut prev: usize = 0_usize;
                let mut total_size = 0;
                let mut nb_unitig: usize = 0;
                for pair in &unitigs_sizes_list{

                    omni_file.write_all(&pair.0);
                    total_size += pair.0.len();
                    let delta_encoded: usize = pair.1 - prev;
                    sizes_encoder.write_all(&delta_encoded.to_le_bytes());
                    nb_unitig += 1;
                    prev = pair.1;
                    //println!("{curr_id}");
                }
                pos_nb_unitig.push((total_size, nb_unitig));
                unitigs_sizes_list.clear();

                //NB UNITIG + TOTAL SIZES 

                curr_id = color_id;
            }
        }else{
            nb_lines_dna += 1;
            nb_kmer += line.len()-30;
            // TODO DELTA ENCODING MUST BE DONE AFTER SORTING, NOT HERE BUT WHEN CHANGING COLOR BUCKET. 
            //      FIND A WAY TO DO SO
            //sizes_encoder.write_all(&(line.len() - prev).to_le_bytes());
            //println!("Size: {}\nDelta encoded size: {}\nPrev: {}", line.len(), (line.len() - prev), prev);
            //println!("{curr_id}");
            //prev = line.len();
            /*println!("TIG CLEAR: {}", &line);
            println!("TIG TO VEC TO STR: {}", vec2str(&str2num(&line), &line.len()));
            println!("SIZE BEFORE: {}\nSIZE AFTER: {}", &line.len(), str2num(&line).len());
            let mut input = String::new();
            std::io::stdin().read_line(&mut input).expect("error: unable to read user input");*/
            unitigs_sizes_list.push((str2num(&line), line.len()));
        }
        line.clear();
        //size_read = reader.read_line(&mut line);
    }
    if !unitigs_sizes_list.is_empty(){
        let mut prev = 0;
        let mut total_size = 0;
        let mut nb_unitig: usize = 0;
        for pair in &unitigs_sizes_list{
            omni_file.write_all(&pair.0);
            //let size = unitig.len();
            total_size += pair.0.len();
            sizes_encoder.write_all(&(pair.1 - prev).to_le_bytes());
            prev = pair.1;
            nb_unitig += 1;
        }
        pos_nb_unitig.push((total_size, nb_unitig));
        unitigs_sizes_list.clear();
    }

    println!("I HAVE SEEN {} LINES WITH DNA", nb_lines_dna);
    println!("I HAVE SEEN {} K-MERS", nb_kmer);

    omni_file.flush();
    sizes_encoder.finish();

    // STARTING POSITIONS IN THE TIGS FILE FOR EACH COLOR BUCKET (ACTUALLY THE SIZE OF EACH COLOR BUCKET) + NB UNITIG IN EACH COLOR BUCKET
    // NB UNITIGS GIVES THE NUMBER OF SIZES == 64B * NB UNITIGS = POS OF COLOR BUCKET IN THE POS FILE
    pos_nb_unitig
}

fn write_positions(pos_nb_unitigs: Vec<(usize, usize)>, filepath: String){
    let pos_file = BufWriter::new(File::create(filepath).expect("unable to create file"));
    let mut pos_encoder = Encoder::new(pos_file, 12).expect("Failed to create zstd encoder");
    let mut prev: usize = 0;
    for elem in pos_nb_unitigs{
        pos_encoder.write_all(&prev.to_le_bytes());
        prev = elem.0;
        let pos_in_sizes: usize = 64_usize*elem.1;
        pos_encoder.write_all(&pos_in_sizes.to_le_bytes());
    }
    pos_encoder.write_all(&prev.to_le_bytes());
    pos_encoder.finish();
}


// ID TO CID -> COMPRESS FURTHER BY WRITING DIFFERENCES BETWEEN CID VECS

fn write_id_to_color_id(cid_file_path: String, id_to_color_vec: Vec<Vec<usize>>) -> std::io::Result<()>{
    println!("WRITTING ID TO COLOR ID FILE: {}", cid_file_path);
    let cid_file = BufWriter::new(File::create(&cid_file_path).expect("unable to create file"));
    let mut cid_encoder = Encoder::new(cid_file, 12).expect("Failed to create zstd encoder");
    for elem in id_to_color_vec{
        let mut to_write = String::new();
        let mut i: u16 = 0;
        let mut prev: u16 = 0;
        for e in elem{
            //to_write += &(e.to_u16().unwrap()-prev).to_string();
            //cid_encoder.write_all(&(e.to_u16().unwrap()-prev).to_le_bytes())?;
            
            if i != 0{
                to_write = to_write + "," + &(e.to_u16().unwrap()).to_string();// - prev).to_string();
            }else {
                to_write = to_write + &(e.to_u16().unwrap()).to_string();// - prev).to_string();
                i += 1;
            }
            prev = e.to_u16().unwrap();
        }
        to_write = to_write + "\n";
        cid_encoder.write_all(to_write.as_bytes())?;
    }
    cid_encoder.finish()?;
    Ok(())
}




//   =========================================================================================== DECOMPRESSION ==============================================================================

pub fn decompress(size_filename: &String, color_id_filename: &String, tigs_filename: &String, positions_filename: &String, filename_id: &String, out_dir: &String, wanted_files_path: &str, input_dir: String) -> std::io::Result<()>{

    let input_file = File::open(input_dir.clone() + filename_id).unwrap();
    let input_reader = BufReader::new(input_file);
    let mut filenames = Vec::new();
    for line_result in input_reader.lines(){
        let line = line_result?;
        if let Some((path, _)) = line.split_once(':'){
            filenames.push(path.to_string());
            println!("File: {}", path.to_string());
        }
    }
    println!("Writting decompressed data in {out_dir}");

    let position_file = File::open(input_dir.clone()+positions_filename).unwrap();
    let positions_reader = BufReader::new(position_file);
    let mut positions_decoder = Decoder::new(positions_reader).expect("Failed to decode color file");
    let mut positions: Vec<(usize, usize)> = Vec::new();

    let cid_to_id_map = get_cid_to_id(&(out_dir.clone() + &color_id_filename));
    let mut buffer_position_tigs = [0; 8];
    for _i in 0_usize..cid_to_id_map.len(){
    let mut buffer_position_size = [0; 8];
        positions_decoder.read_exact(&mut buffer_position_tigs);
        positions_decoder.read_exact(&mut buffer_position_size);
        positions.push((usize::from_le_bytes(buffer_position_tigs), usize::from_le_bytes(buffer_position_size)));
    }
    positions_decoder.read_exact(&mut buffer_position_tigs);
    positions.push((usize::from_le_bytes(buffer_position_tigs), 0_usize));

    if wanted_files_path != ""{
        //decompress_wanted(wanted_files_path, filenames, cid_to_id_map, tigs_filename, &size_filename, out_dir);
    }else{
        decompress_all(&(input_dir.to_owned()+size_filename), &positions, &(input_dir.to_owned()+tigs_filename), out_dir, filenames, cid_to_id_map);
    }

    Ok(())
}

fn decompress_all(size_filename: &String, positions: &Vec<(usize, usize)>, tigs_filename: &String, out_dir: &String, filenames: Vec<String>, cid_to_id_map: HashMap<String, Vec<u32>>){

    let mut tigs_file = BufReader::new(File::open(&tigs_filename).expect("Error opening tigs file, are you sure path is good?"));

    let size_file = BufReader::new(File::open(size_filename).expect("Error opening size file, are you sure path is good?"));
    let mut size_decoder = Decoder::new(size_file).expect("Failed to decode color file");
    let mut size_buf = [0; 8];
    let mut cursor: usize = 0;
    let mut bucket_cursor = 0;
    println!("nb pos: {}", positions.len());
    println!("nb cid: {}", cid_to_id_map.len());
    for i in 0_usize..cid_to_id_map.len(){
    //for i in 0_usize..filenames.len(){
        let mut prev: usize = 0;
        println!("BUCKET SIZE: {}", positions.get(i+1).unwrap().0);
        println!("I = {}", i);
        println!("CURSOR: {}", bucket_cursor);
        //let mut input = String::new();
        //std::io::stdin().read_line(&mut input).expect("error: unable to read user input");
        while bucket_cursor < positions.get(i+1).unwrap().0{
            size_decoder.read_exact(&mut size_buf);
            let size = usize::from_le_bytes(size_buf)+prev;
            prev = size;
            cursor += size.div_ceil(4);
            bucket_cursor += size.div_ceil(4);
            let read_size = size.div_ceil(4);
            let mut bucket_buffer = vec![0; read_size];
            tigs_file.read_exact(&mut bucket_buffer);
            let tig = vec2str(&bucket_buffer, &size);
            let curr_output_files = cid_to_id_map.get(&i.to_string()).unwrap();
            // if tig.contains("TAAACCAACGTATTCGATAAGACCGTCAACA")  || tig.contains("TGTTGACGGTCTTATCGAATACGTTGGTTTA"){
            //     println!("SIZE READ: {}", usize::from_le_bytes(size_buf));
            //     println!("UNITIG SIZE (STRING) {}", size);
            //     println!("READING SIZE: {}", size.div_ceil(4));
            //     println!("POS FIN CURR BUCKET: {}", positions.get(i+1).unwrap().0);
            //     println!("CURR POS IN BUCKET: {}", bucket_cursor);
            //     println!("I = {}", i);
            //     let mut input = String::new();
            //     std::io::stdin().read_line(&mut input).expect("error: unable to read user input");
            //     println!("TIG: {}", tig);
            //     println!("NB FILE: {}", curr_output_files.len());
            // }
            // WRITE TIG IN EVERY FILE PRESENT IN CURRENT COLOR SET
            for elem in curr_output_files{
                let curr_filename = filenames.get(*elem as usize).unwrap();
                //println!("FILE NUMBER: {}", elem);
                //println!("a{}a", curr_filename);
                let trunc_filename = Path::new(curr_filename).file_stem().unwrap();
                let mut out_file = BufWriter::new(File::options().append(true).create(true).open(out_dir.to_owned() + "Dump_" + trunc_filename.to_str().unwrap()+".fa").expect("Unable to create file"));
                out_file.write_all(">\n".as_bytes());
                out_file.write_all(&tig.as_bytes());
                out_file.write_all("\n".as_bytes());
                out_file.flush();
            }
        }
        bucket_cursor = 0;
    }
}

fn get_cid_to_id(color_id_filename: &String) -> HashMap<String, Vec<u32>>{
    let color_id_file = BufReader::new(File::open(color_id_filename).expect("Error opening color id file, are you sur you gave the right path ?"));
    let mut color_id_decoder = Decoder::new(color_id_file).expect("Failed to decode color file");
    let mut colors = String::new();
    color_id_decoder.read_to_string(&mut colors);
    let lines = colors.split("\n").collect::<Vec<_>>();
    let mut ids = HashMap::new();
    let mut counter: u32 = 0;
    println!("READING CID TO ID FILE: {}", color_id_filename);
    println!("LINES: {}", colors);
    for line in lines{
        println!("LINE: {}", line);
        let mut input = String::new();
        std::io::stdin().read_line(&mut input).expect("error: unable to read user input");
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
/*
fn filter_cid_to_id(wanted_files_path: &String, cid_to_id_map: HashMap<String, Vec<u32>>, filenames: Vec<String>){
    let wanted_file = File::open(wanted_files_path).unwrap();
    let wanted_reader = BufReader::new(wanted_file);
    let mut wanted_filenames = Vec::new();
    let mut wanted_cids_map = HashMap::new();
    for line_result in wanted_reader.lines(){
        let line = line_result?;
        if filenames.contains(&line){
            let mut cpt: u32 = 0;
            let wanted_id;
            for file in filenames{
                if file == line{
                    wanted_id = cpt;
                    cpt += 1;
                }
            }
            for i in 0..cid_to_id_map.len(){
                if cid_to_id_map.get(&i.to_string()).unwrap().contains(&wanted_id){
                    wanted_cids_map.insert(i, cid_to_id_map.iter().get(&i.to_string()).unwrap());
                }
            }
            wanted_filenames.push(line);
            println!("File: {}", line);
        }else{
            println!("File: {} is not in the archive, are you sure accession is correct ?", line);
        }
        
    }

}

fn decompress_wanted(wanted_files_path: &String, filenames: Vec<String>, cid_to_id_map: HashMap<String, Vec<u32>>, tigs_filename: &String, size_filename: &String, out_dir: &String){

}*/
/*
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
    //let sizes = get_cid_to_bucket_size(&(out_dir.clone() + &size_filename), cid_to_id_map.len());

    // ============================================================== FILTER NECESSARY IDS =================================================================== //

    if wanted_files_path != ""{
        decompress_wanted(wanted_files_path, filenames, cid_to_id_map, tigs_filename, &size_filename, out_dir);
    }else{
        decompress_all(&size_filename, cid_to_id_map, tigs_filename, out_dir, filenames);
    }
    Ok(())
    //TODO
    // build unitig graph for each output using ggcat

}


// TODO READ BUCKET BY BUCKET, WHEN 0 IS READ, CHECK IF BUCKET SHOULD BE READ IF YES GO IF NOT FLUSH DATA AND GO TO NEXT BUCKET
// DO THAT WITH VARIABLE BYTES -> LIB BYTESVARINT ? 
// HANDLE EOF, NB OF SIZES UNKNOWN
// ELSE SAVE NB UNITIGS IN OTHER FILE SOMEWHERE. 
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

fn decompress_wanted(wanted_files_path: &str, input_filenames: Vec<String>, cid_to_id_map: HashMap<String, Vec<u32>>, tigs_filename: String, size_filename: &String, out_dir: &String){
    let wanted_ids = get_to_decompress(wanted_files_path, input_filenames, &cid_to_id_map);
    //let mut tigs_decoder = Decoder::new(tigs_file).expect("Failed to decode tigs file");

    // CREATE/OPEN FILENAMES
    let mut dump_files = Vec::new();
    let input_fof = File::open(wanted_files_path).unwrap();
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
    decompress_needed(size_filename, wanted_ids, tigs_filename, &mut dump_files);

    
    //concat_file.flush();
}

fn decompress_needed(size_filename: &String, wanted_ids: HashMap<String, Vec<usize>>, tigs_filename: String, dump_files: &mut Vec<BufWriter<File>>){
    let mut counter: usize = 0;
    let mut cursor: usize = 0;
    let mut tigs_file = BufReader::new(File::open(&tigs_filename).expect("Error opening tigs file, are you sure path is good?"));
    // ADD LOOP ON EVERY SIZE PER BUCKET UP UNITL EITHER A 0 OR TOTAL SIZE THEN AND ONLY THEN, INCREMENT COUNTER BY ONE TO CHANGE ID.
    
    let size_file = BufReader::new(File::open(size_filename).expect("Error opening size file, are you sure path is good?"));
    let mut size_decoder = Decoder::new(size_file).expect("Failed to decode color file");
    let mut size_buf = [0; 8];
    let mut size: usize;

    for id in wanted_ids{
        size_decoder.read_exact(&mut size_buf);
        size = usize::from_le_bytes(size_buf);
        while size != 0_usize{
            if id.0 == counter.to_string(){
                let mut bucket_buffer = vec![0; size];
                tigs_file.seek(std::io::SeekFrom::Start(cursor.try_into().unwrap()));
                tigs_file.read_exact(&mut bucket_buffer);
                let mut tigs_decoder = Decoder::new(&bucket_buffer[..]).expect("Failed to create zstd decoder");
                let mut clear_color_bucket = Vec::new();
                tigs_decoder.read_to_end(&mut clear_color_bucket);
                let output = String::from_utf8(clear_color_bucket).unwrap();
                for out_file in dump_files.iter_mut(){
                    out_file.write_all(output.as_bytes());
                    out_file.write_all(b"\n>\n");
                }
            }
            cursor += size;
            size_decoder.read_exact(&mut size_buf);
        }
        counter += 1;
    }
    for file in dump_files.iter_mut(){
        file.flush();
    }
    //OLD
    /*
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
        if size == &0_usize {
            counter += 1;
        }
    }
    for file in dump_files.iter_mut(){
        file.flush();
    }
    */
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
// TODO FOR DECOMPRESSION (BOTH TARGETED AND ALL) GET TOTAL NUMBER OF UNITIGS AND STOP READING SIZES WHEN TOTAL NUMBER HAS BEEN READ
fn decompress_all(size_filename: &String, cid_to_ids: HashMap<String, Vec<u32>>, tigs_filename: String, out_dir: &String, filenames: Vec<String>){
    //let mut nb_kmers = 0;
    //let mut concat_file = BufWriter::new(File::create(String::from("concat_all") + ".fa").expect("Unable to create file"));
    let mut tigs_file = BufReader::new(File::open(&tigs_filename).expect("Error opening tigs file, are you sure path is good?"));
    //let mut tigs_decoder = Decoder::new(tigs_file).expect("Failed to decode tigs file");
    let mut counter = 0;
    //let mut cursor: usize = 0;
    let size_file = BufReader::new(File::open(size_filename).expect("Error opening size file, are you sure path is good?"));
    let mut size_decoder = Decoder::new(size_file).expect("Failed to decode color file");
    let mut size_buf = [0; 8];
    let mut size: usize;
    //concat_file.write_all(b">\n");

    //TODO READ SIZE BY SIZE DECOMPRESS UNITIG BY UNITIG
    // TODO WIP
    size_decoder.read_exact(&mut size_buf);
    size = usize::from_le_bytes(size_buf);
    let mut bucket_buffer = vec![0; size];
    let mut read_unitig_res = tigs_file.read_exact(&mut bucket_buffer);
    while read_unitig_res.is_ok(){
        while size != 0_usize{

            println!("size = {}", size);
            //tigs_file.seek(std::io::SeekFrom::Start(cursor.try_into().unwrap()));
            tigs_file.read_exact(&mut bucket_buffer);
            let mut tigs_decoder = Decoder::new(&bucket_buffer[..]).expect("Failed to create zstd decoder");
            let mut clear_color_bucket = Vec::new();
            tigs_decoder.read_to_end(&mut clear_color_bucket);
            println!("TIG = {}", String::from_utf8(clear_color_bucket.clone()).unwrap());
            //let output = String::from_utf8(clear_color_bucket).unwrap();
            let curr_output = cid_to_ids.get(&counter.to_string()).unwrap();
            for elem in curr_output{
                let curr_filename = filenames.get(*elem as usize).unwrap();
                println!("a{}a", curr_filename);
                let trunc_filename = Path::new(curr_filename).file_stem().unwrap();
                let mut out_file = BufWriter::new(File::options().append(true).create(true).open(out_dir.to_owned() + "Dump_" + trunc_filename.to_str().unwrap()).expect("Unable to create file"));
                out_file.write_all(&clear_color_bucket);
                out_file.flush();
            }
            //cursor += size;
            size_decoder.read_exact(&mut size_buf);
            size = usize::from_le_bytes(size_buf);
            bucket_buffer = vec![0; size];
            read_unitig_res = tigs_file.read_exact(&mut bucket_buffer);
        }
        println!("COUNTER = {}", counter);
        counter += 1;
    }
    //CONCAT FILE IS FOR DEBUGGING, CHECK IF NB KMER IS CONSISTENT
    /*
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
     */
    //concat_file.flush();
    //println!("I HAVE SEEN {} K-MERS", nb_kmers);
    // TODO ENCODE OMNICOLOR AS NORMAL COLOR 

} */
