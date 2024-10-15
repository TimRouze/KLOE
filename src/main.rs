#![allow(dead_code)]

mod utils;
mod kmer;
mod decompress;
mod stats;
use clap::Parser;
use stats::compute_stats;
use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use std::sync::atomic::{AtomicBool, Ordering};
use std::fs::{remove_file, File};
use std::path::PathBuf;
use std::io::{self, BufRead, Read, Seek, SeekFrom, Stdin, Write, BufReader};
use std::env;
use::rayon::prelude::*;
use bitvec::prelude::*;
use kmer::{Kmer, RawKmer};
use std::cell::Cell;
use needletail::{parse_fastx_file, Sequence};
use indexmap::IndexMap;

use crate::utils::{num2str, str2num, };

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    ///output: Option<String>,
    /// Decompression "compress" to compress input, "decompress" to decompress input, "stats" to get kmer stats
    decompress: Option<String>,
    ///multicolor file (only for decompression)
    #[arg(long, default_value_t=String::from(""))]
    multicolor_file: String,
    /// Omnicolor file (only for decompression)
    #[arg(long, default_value_t=String::from(""))]
    omnicolor_file: String,
    /// Number of threads (defaults to all available threads)
    #[arg(short, long, default_value_t = 1)]
    threads: usize,
    ///Output directory
    #[arg(short, long, default_value_t = String::from(""))]
    out_dir: String,
    ///Number of expected k-mers
    #[arg(short = 'N', long, default_value_t = 1_000_000)]
    nb_elem: usize,
    ///List of files to decompress
    #[arg(short = 'Q', long, default_value_t = String::from(""))]
    wanted_files: String,
}
pub mod constants {
    include!("constants.rs");
}
use constants::{ARRAY_SIZE, NB_FILES, INPUT_FOF, K, KT};
pub type COLORPAIR = (bitvec::prelude::BitArray<[u8; ARRAY_SIZE]>, Cell<bool>);

//TODO QUERY
    //TODO INTERFACE FILE: FILENAME TO COLOR
    //TODO INTERFACE FILE: COLOR TO SIMPLITIG SIZE
    //TODO MAP COLOR -> SIMPLITIGS TOTAL SIZE
//TODO SET COMPRESSION RATIO AS ARG
fn main() {
    let args = Args::parse();
    let nb_elem = args.nb_elem;
    println!("FILENAME: {}", INPUT_FOF);
    let output_dir = args.out_dir;
    env::set_var("RAYON_NUM_THREADS", args.threads.to_string());
    let file = File::open(INPUT_FOF).unwrap();
    let reader = io::BufReader::new(file);
    let kmer_map_mutex: Arc<Mutex<HashMap<KT, COLORPAIR>>> = Arc::new(Mutex::new(HashMap::with_capacity(nb_elem)));
    let filenames: Vec<_> = reader.lines().collect::<Result<_, _>>().unwrap();
    let wanted_path = args.wanted_files;
    if let Some(do_decompress) = args.decompress{
        if do_decompress == "decompress"{
            let multi_file = args.multicolor_file;
            let omni_file = args.omnicolor_file;
            if multi_file != "" && omni_file != ""{
                decompress::decompress(&omni_file, &multi_file, INPUT_FOF, PathBuf::from(output_dir), &wanted_path);
            }else{
                println!("Error, multicolor and/or omnicolor file(s) are mandatory");
            }
        }else if do_decompress == "compress"{
            let mut filename_color = File::create("filename_to_color.txt").unwrap();
            for i in 0..NB_FILES{
                writeln!(filename_color, "{}:{}", filenames.get(i).unwrap(), i).unwrap();
            }
            (0..NB_FILES).into_par_iter().for_each(|file_number|{
                let mut kmer_map_mutex = Arc::clone(&kmer_map_mutex);
                println!("{}", filenames.get(file_number).unwrap());
                read_fasta(&filenames.get(file_number).unwrap(), &mut kmer_map_mutex, file_number, &nb_elem);
            });
            let mut kmer_map = Arc::try_unwrap(kmer_map_mutex).expect("Failed to unwrap Arc").into_inner().expect("Failed to get Mutex");
            println!("NB KMER = {}", kmer_map.len());
            kmer_map.shrink_to_fit();
            let mut color_simplitig = compute_colored_simplitigs(&mut kmer_map, &output_dir);
            sort_simplitigs(&mut color_simplitig, &output_dir)
        }else if do_decompress == "stats"{
            let multi_file = args.multicolor_file;
            let omni_file = args.omnicolor_file;
            let k = K;
            if multi_file != "" && omni_file != ""{
                compute_stats(&omni_file, &multi_file, &output_dir, &k);
            }else{
                println!("Error, multicolor and/or omnicolor file(s) are mandatory");
            }
        }
    }else {
        println!("Wrong positional arguments given. Values are 'compress' or 'decompress'");
        println!("Ex: if compression: I=my/fof.txt cargo r -r -- compress -o out_dir/ -t 12");
        println!("Ex: if decompression: I=my/fof.txt cargo r -r -- decompress -omnicolor-file out_dir/omnicolor.fa.zstd --multicolor-file out_dir/multicolor.fa.zstd -t 12");
    }
}
/*
READS FASTA FILES.
CREATES K-MER -> COLOR MAP.
*/
fn read_fasta(filename: &str, kmer_map_mutex: &Arc<Mutex<HashMap<KT, COLORPAIR>>>, file_number: usize, nb_elem: &usize){
    let mut fa_reader = parse_fastx_file(filename).expect("Error while opening file");
    let mut counter = 0;
    let mut counter_insert = 0;
    let mut counter_modify = 0;
    //let mut to_add = HashSet::new();
    while let Some(record) = fa_reader.next(){
        let seqrec = record.expect("Error reading record");
        let norm_seq = seqrec.normalize(false);
        let norm_rc = norm_seq.reverse_complement();
        let canon_kmers = norm_seq.canonical_kmers(K as u8, &norm_rc);
        canon_kmers.for_each(|kmer|{
            //println!("{}", std::str::from_utf8(&kmer.1).unwrap());
            
            let curr_kmer: RawKmer<K, KT> = RawKmer::from_nucs(kmer.1);
            counter += 1;
            let mut kmer_map = kmer_map_mutex.lock().unwrap();
            //let file_nb = *file_number.read();
            kmer_map.entry(curr_kmer.to_int()).and_modify(|pair|{
                counter_modify += 1;
                pair.0.set(file_number, true);
            }).or_insert_with(|| {
                let mut bv: BitArr!(for NB_FILES, in u8) = BitArray::<_>::ZERO;
                bv.set(file_number, true);
                counter_insert += 1;
                (bv, Cell::new(false))
            });
            //TODO FIND A WAY TO RESIZE WITH DASHMAP
            /*if kmer_map.capacity() <= (20/100)*nb_elem{
                kmer_map.try_reserve((30/100)*nb_elem);
            }*/
        });
    }
    println!("I have inserted {} k-mers", counter_insert);
    println!("I have seen {} already inserted k-mers", counter_modify);
    println!("I have seen {} k-mers", counter);
}


/*
CONSTRUCTS SIMPLITIGS FROM K-MERS GATHERED IN "READ_FASTA".
SIMPLITIGS ARE MONOCHROMATIC.
OMNICOLORED (SEEN IN EVERY INPUT FILE) SIMPLITIGS ARE WRITTEN IN THE SAME FILE (omnicolor.fa.zstd)
MULTICOLORED SIMPLITIGS ARE WRITTEN IN ANOTHER FILE (multicolor.fa.zstd).
SIMPLITIGS ARE CONSTRUCTED USING PROPHASM'S GREEDY ALGORITHM.
*/
fn compute_colored_simplitigs(kmer_map:  &mut HashMap<KT, COLORPAIR>, output_dir: &String) -> HashMap<bitvec::prelude::BitArray<[u8; ARRAY_SIZE]>, (usize, Vec<usize>)>{
    println!("I will reconstruct simplitigs from {} kmers", kmer_map.len());
    let mut multi_f = File::create(output_dir.clone()+"temp_multicolor.fa").expect("Unable to create file");
    //let mut omni_f = Encoder::new(File::create(output_dir.clone()+"omnicolor.kloe").expect("Unable to create file"), 0).unwrap();
    let mut omni_f = File::options().write(true).read(true).open(output_dir.clone()+"omnicolor.kloe").expect("unable to create file");
    let mut kmer_iterator = kmer_map.iter();
    let mut color_simplitig_size: HashMap<bitvec::prelude::BitArray<[u8; ARRAY_SIZE]>, (usize, Vec<usize>)> = HashMap::new();
    let mut omni: BitArr!(for NB_FILES, in u8) = BitArray::<_>::ZERO;
    let mut is_omni = false;
    for i in 0..NB_FILES{
        omni.set(i, true);
    }
    while let Some((key, curr_cell)) = kmer_iterator.next(){
        if !curr_cell.1.get(){
            if curr_cell.0.eq(&omni){
                is_omni = true;
            }
            curr_cell.1.set(true);
            let mut forward = true;
            let mut backward = true;
            let mut simplitig = num2str(*key);
            while forward{
                let curr_kmer = RawKmer::<K, KT>::from_nucs(&simplitig[(simplitig.len()-K)..].as_bytes());
                forward = extend_forward(&curr_kmer, &kmer_map, &mut simplitig, &curr_cell.0);
            }
            while backward{
                let curr_kmer = RawKmer::<K, KT>::from_nucs(&simplitig[0..K].as_bytes());
                backward = extend_backward(&curr_kmer, &kmer_map, &mut simplitig, &curr_cell.0);
            }
            let simplitig_size = simplitig.len();
            color_simplitig_size.entry(curr_cell.0)
                        .and_modify(|total_size| {
                            total_size.0 += simplitig_size.div_ceil(4);
                            total_size.1.push(simplitig_size);
                        })
                        .or_insert((simplitig_size.div_ceil(4), vec![simplitig_size]));
            write_simplitig(&str2num(&simplitig), &is_omni, &mut multi_f, &mut omni_f, &curr_cell.0, &simplitig.len());
            is_omni = false;
        }
    }
    color_simplitig_size
}
/*
FORWARD EXTENSION FOR SIMPLITIG CREATION.
CHECKS IF SUCCESSORS ARE THE SAME COLOR AS CURRENT K-MER.
*/
fn extend_forward(curr_kmer: &RawKmer<K, KT>, kmer_map:  &HashMap<KT, COLORPAIR>, simplitig: &mut String, color: &BitArray<[u8;ARRAY_SIZE]>) -> bool{
    for succs in curr_kmer.successors(){
        //forward = false;
        if kmer_map.contains_key(&succs.canonical().to_int()){
            let succ_pair = kmer_map.get(&succs.canonical().to_int()).unwrap();
            if succ_pair.0.eq(color) & !succ_pair.1.get() {
                simplitig.push(*succs.to_nucs().last().unwrap() as char);
                succ_pair.1.set(true);
                return true;
            }
            drop(succ_pair);
        }
    }
    false
}

/*
BACKWARD EXTENSION FOR SIMPLITIG CREATION.
CHECKS IF PREDECESSORS ARE THE SAME COLOR AS CURRENT K-MER.
INSERTS FIRST NUCLEOTIDE OF PREDECESSOR (CHECKED MULTIPLE TIMES)
*/
fn extend_backward(curr_kmer: &RawKmer<K, KT>, kmer_map:  &HashMap<KT, COLORPAIR>, simplitig: &mut String, color: &BitArray<[u8;ARRAY_SIZE]>) -> bool{
    for preds in curr_kmer.predecessors(){
        //backward = false;
        if kmer_map.contains_key(&preds.canonical().to_int()){
            let pred_pair = kmer_map.get(&preds.canonical().to_int()).unwrap();
            if pred_pair.0.eq(color) & !pred_pair.1.get() {
                simplitig.insert(0,*preds.to_nucs().first().unwrap() as char);
                pred_pair.1.set(true);
                return true;
            }
            drop(pred_pair);
        }
    }
    false
}

fn compute_cursor_start_pos(color_simplitig: &HashMap<bitvec::prelude::BitArray<[u8; ARRAY_SIZE]>, (usize, Vec<usize>)>) -> IndexMap<bitvec::prelude::BitArray<[u8; ARRAY_SIZE]>, u64>{
    let mut color_to_cursor = IndexMap::new();
    let mut prev_cursor: u64 = 0;
    for (key, value) in color_simplitig.iter(){
        if key.count_ones() != NB_FILES{
            color_to_cursor.entry(*key).or_insert({
                let proxy = prev_cursor;
                prev_cursor += value.0 as u64;
                proxy
            });
        }
    }
    color_to_cursor
}
// save simplitig as vec of u8, save taille reelle, pour obtenir la taille u8: taille reelle/4 arrondi sup.
// COMME SUIT: 00110001 taille_a_lire U8U8U8U8U8
// THEN ON LES LIT UN A UN POUR LES TRIER, LA TAILLE DOIT ETRE STOCKEE QQPART (COLOR_SIMPLITIG_SIZE ?)
// ADDITION DES TAILLES/4 ARRONDI SUP POUR AVOIR LA TAILLE DU BUCKET ET AINSI CALCULER LES CURSEURS
// REECRITURE TRIÉ PUIS COMPRESSION
fn sort_simplitigs(color_simplitig: &mut HashMap<bitvec::prelude::BitArray<[u8; ARRAY_SIZE]>, (usize, Vec<usize>)>, output_dir: &String ) {
    let path = output_dir.clone()+"multicolor.kloe";
    //let compressed_path = output_dir.clone()+"multicolor.zstd";
    let mut out_mult_file = File::options().write(true).read(true).create(true).open(path.clone()).expect("Unable to create file");
    //let mut color_set = HashSet::new();
    let mut color_cursor_pos = compute_cursor_start_pos(color_simplitig);
    let temp_multicolor_file = File::open(output_dir.clone()+"temp_multicolor.fa").unwrap();
    let mut temp_multicolor_reader = BufReader::new(&temp_multicolor_file);
    const color_size: usize = NB_FILES.div_ceil(8);
    let mut cursor: usize = 0;
    let metadata = temp_multicolor_file.metadata().unwrap();
    let file_size: usize = metadata.len() as usize;
    //println!("FILE SIZE = {}", file_size);
    while cursor < file_size{
        cursor+= color_size;
        let mut id = [0; color_size];
        temp_multicolor_reader.read_exact(&mut id).expect("Error reading color in temp file");
        let mut size_buf = [0; 4];
        cursor += 4;
        temp_multicolor_reader.read_exact(&mut size_buf).expect("Error reading simplitig size in temp file");
        let size_to_read: u32 = u32::from_le_bytes(size_buf).div_ceil(4);
        let size_simplitig: u32 = u32::from_le_bytes(size_buf);
        //println!("READING CURSOR = {}", cursor);
        cursor += size_to_read as usize;
        //println!("COLOR: {}\nSIZE: {}", &id.to_vec()[0], size_simplitig);
        let mut simplitig = vec![0; size_to_read as usize];
        //println!("Reading {} Bytes", simplitig.len());
        temp_multicolor_reader.read_exact(&mut simplitig).expect("Error reading simplitig");
        //println!("SIMPLITIG: {}", vec2str(&simplitig.to_vec(), &(size_simplitig as usize)));
        let _ = write_sorted(&mut out_mult_file, simplitig, &mut color_cursor_pos, &id.to_vec()[0]);
        //println!("READING CURSOR = {}", cursor);
        //let mut input = String::new();
        //io::stdin().read_line(&mut input).expect("error: unable to read user input");
    }
    remove_file(output_dir.clone()+"temp_multicolor.fa");
    //zstd_compress_file(&out_mult_file, compressed_path);
    write_interface_file(color_simplitig, color_cursor_pos, output_dir);
}

/* fn zstd_compress_file(in_file: &File, output_path: &Path) -> io::Result<()> {
    let in_file = File::open(input_path)?;
    let mut reader = BufReader::new(input_file);

    let output_file = File::create(output_path)?;
    let mut writer = BufWriter::new(output_file);

    copy_encode(&mut reader, &mut writer, 9)?;
    Ok(())
} */

fn write_sorted(out_mult_file: &mut File, content: Vec<u8>, color_cursor_pos: &mut IndexMap<bitvec::prelude::BitArray<[u8; ARRAY_SIZE]>, u64>, color: &u8) -> io::Result<()> {
    //let mut encoder = Compressor::new(12).unwrap();
    //let compressed = encoder.compress(content.as_bytes()).unwrap();
    let mut color_vec: BitArr!(for NB_FILES, in u8) = BitArray::<_>::ZERO;
    for i in 0..NB_FILES{
        let bit = (color >> i) & 1;
        if bit == 1{
            color_vec.set(i,true);
        }
    }
    let position = color_cursor_pos.get(&color_vec).unwrap();
    out_mult_file.seek(SeekFrom::Start(*position))?;
    //let mut encoder = zstd::Encoder::new(out_mult_file, 9);
    out_mult_file.write_all(&content)?;
    //*position += compressed.len() as u64;
    //encoder.finish()?;
    color_cursor_pos.entry(color_vec).and_modify(|cursor| {*cursor += content.len() as u64});
    Ok(())
}



fn write_interface_file(color_simplitig: &HashMap<bitvec::prelude::BitArray<[u8; ARRAY_SIZE]>, (usize, Vec<usize>)>, color_cursor_pos: IndexMap<bitvec::prelude::BitArray<[u8; ARRAY_SIZE]>, u64>, output_dir: &String){
    let mut file = File::create(output_dir.clone()+"multicolor_bucket_size.txt").expect("Unable to create interface file");
    let mut nb_seen = 0;
    for (key, cursor_end) in color_cursor_pos{
        let mut color = String::new();
        for e in key.iter(){
            if *e{
                color.push('1');
                nb_seen += 1;
            }else{
                color.push('0');
            }
        }
        let mut bucket_sizes = String::new();
        let curr_values = color_simplitig.get(&key).unwrap();
        for e in curr_values.1.iter(){
            bucket_sizes += &(e.to_string().to_owned()+",");
        }
        if nb_seen < NB_FILES{
            bucket_sizes.pop();
            writeln!(file, "{},{}:{}", color, cursor_end, bucket_sizes).unwrap();
        }
        nb_seen = 0;
    }
}

fn write_simplitig(simplitig: &Vec<u8>, is_omni: &bool, multi_f: &mut File, omni_f: &mut File, color: &BitArray<[u8;ARRAY_SIZE]>, size_str : &usize){
    let size:u32 = *size_str as u32;
    if *is_omni{
        omni_f.write_all(&size.to_le_bytes()).unwrap();  
        omni_f.write_all(&simplitig).expect("Unable to write data");
    }else{
        let mut id = Vec::new();
        let mut cpt = 0;
        let mut curr_int: u8 = 0;
        let mut str_color = String::new();
        for e in color.iter(){
            if *e{
                curr_int +=1 << cpt;
                str_color.push('1');
            }else {
                str_color.push('0');
            }
            cpt += 1;
            if cpt == 8{
                cpt = 0;
                id.push(curr_int);
                curr_int = 0;
            }
        }
        //println!("");
        //println!("str color: {}", str_color);
        //println!("color size = {}", id.len());
        //println!("simplitig size = {}", simplitig.len());
        //println!("Simplitig after convert: {}", vec2str(&simplitig, size_str));
        //println!("SIZE: {}", simplitig.len());
        //println!("Color size {}", id.len());
        //let mut input = String::new();
        //io::stdin().read_line(&mut input).expect("error: unable to read user input");
        multi_f.write_all(&id).expect("Unable to write color");
        multi_f.write_all(&size.to_le_bytes()).expect("Unable to write simplitig size");  
        multi_f.write_all(&simplitig).expect("Unable to write data");
    }
}

#[test]
fn test_bitvec_eq_omni(){
    let mut omni: BitArr!(for 8, in u8) = BitArray::<_>::ZERO;
    for i in 0..8{
        omni.set(i, true);
    }
    assert_eq!(omni.all(), true);
}

#[test]
fn test_bitvec_eq_non_omni(){
    let mut omni: BitArr!(for 8, in u8) = BitArray::<_>::ZERO;
    for i in 0..7{
        omni.set(i, true);
    }
    assert_eq!(omni.all(), false);
}

#[test]
fn test_rev_comp_str(){
    let kmer = String::from("AGCTTTTCATTCTGACTGCAACGGGCAATAT");
    let revcomp = rev_comp_str(&kmer);
    assert_eq!(revcomp, "ATATTGCCCGTTGCAGTCAGAATGAAAAGCT");
}
/*
#[test]
fn test_process_fof_parallel() {
    // Replace "path/to/your/test_file.txt" with the path to your test file
    let filename = "../Data/fof_test.txt";
    let modimizer = false; // Adjust as needed
    let nb_files = 2; // Adjust as needed
    let size = 1_000_000_000; // Adjust as needed
    env::set_var("RAYON_NUM_THREADS", "1");
    // Run the function under test
    let result = process_fof_parallel(filename, modimizer, nb_files, size);

    // Assert that the function returns successfully
    assert!(result.is_ok());

    // Extract the result vector from the mutex
    let hist_mutex = result.unwrap();
    let hist = hist_mutex.lock().unwrap();

    // Add additional assertions based on the expected behavior of your function
    assert_eq!(hist.len(), nb_files + 1);
    assert_eq!(hist[2], 97);
    assert_eq!(hist[0], 0);
    // Add more assertions as needed
}*/