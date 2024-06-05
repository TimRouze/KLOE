#![allow(dead_code)]

mod kmer;
mod decompress;
mod stats;
use clap::Parser;
use stats::compute_stats;
use std::sync::{Arc, Mutex};
use std::fs::File;
use std::path::{Path, PathBuf};
use std::io::{Write, self, BufRead, Stdin};
use std::env;
use::rayon::prelude::*;
use zstd::stream::write::Encoder;
use bitvec::prelude::*;
use kmer::{Kmer, RawKmer};
use std::cell::Cell;
use hashbrown::HashMap;
use needletail::{parse_fastx_file, Sequence};


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
}
pub mod constants {
    include!("constants.rs");
}
use constants::{ARRAY_SIZE, NB_FILES, INPUT_FOF};
const K: usize = 31;
pub type KT = u64;
pub type COLORPAIR = (bitvec::prelude::BitArray<[u8; ARRAY_SIZE]>, Cell<bool>);


//TODO SORT COLORS
//TODO QUERY

fn main() {
    let args = Args::parse();
    let nb_elem = args.nb_elem;
    println!("FILENAME: {}", INPUT_FOF);
    let output_dir = args.out_dir;
    env::set_var("RAYON_NUM_THREADS", args.threads.to_string());
    let file = File::open(INPUT_FOF).unwrap();
    let reader = io::BufReader::new(file);
    let kmer_map_mutex: Arc<Mutex<HashMap<u64, COLORPAIR>>> = Arc::new(Mutex::new(HashMap::with_capacity(nb_elem)));
    let filenames: Vec<_> = reader.lines().collect::<Result<_, _>>().unwrap();
    if let Some(do_decompress) = args.decompress{
        if do_decompress == "decompress"{
            let multi_file = args.multicolor_file;
            let omni_file = args.omnicolor_file;
            if multi_file != "" && omni_file != ""{
                decompress::decompress(&omni_file, &multi_file, INPUT_FOF, PathBuf::from(output_dir));
            }else{
                println!("Error, multicolor and/or omnicolor file(s) are mandatory");
            }
        }else if do_decompress == "compress"{
            (0..NB_FILES).into_par_iter().for_each(|file_number|{
                let kmer_map_mutex = Arc::clone(&kmer_map_mutex);
                println!("{}", filenames.get(file_number).unwrap());
                read_fasta(&filenames.get(file_number).unwrap(), &kmer_map_mutex, file_number, &nb_elem);
            });
            let mut kmer_map = Arc::try_unwrap(kmer_map_mutex).expect("Failed to unwrap Arc").into_inner().expect("Failed to get Mutex");
            println!("NB KMER = {}", kmer_map.len());
            kmer_map.shrink_to_fit();
            compute_colored_simplitigs(&mut kmer_map, &output_dir);
        }else if do_decompress == "stats"{
            let multi_file = args.multicolor_file;
            let omni_file = args.omnicolor_file;
            let k = K;
            if multi_file != "" && omni_file != ""{
                compute_stats(&omni_file, &multi_file, "kmer_per_color_stats.txt", &k);
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

fn read_fasta(filename: &str, kmer_map_mutex: &Arc<Mutex<HashMap<u64, COLORPAIR>>>, file_number: usize, nb_elem: &usize){
    let mut fa_reader = parse_fastx_file(filename).expect("Error while opening file");
    //let mut counter = 0;
    let mut counter_insert = 0;
    let mut counter_modify = 0;
    //let mut to_add = HashSet::new();
    while let Some(record) = fa_reader.next(){
        let seqrec = record.expect("Error reading record");
        let norm_seq = seqrec.normalize(false);
        let norm_rc = norm_seq.reverse_complement();
        let mut canon_kmers = norm_seq.canonical_kmers(31, &norm_rc);
        canon_kmers.for_each(|kmer|{
            //println!("{}", std::str::from_utf8(&kmer.1).unwrap());
            
            let curr_kmer: RawKmer<K, u64> = RawKmer::from_nucs(kmer.1);
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
            if kmer_map.len() >= (80/100)*nb_elem{
                kmer_map.reserve((30/100)*nb_elem);
            }
            drop(kmer_map);
        });
    }
    println!("I have inserted {} k-mers", counter_insert);
    println!("I have seen {} already inserted k-mers", counter_modify);
    //println!("I have seen {} k-mers", counter);
}

fn compute_colored_simplitigs(kmer_map:  &mut HashMap<u64, COLORPAIR>, output_dir: &String){
    println!("I will reconstruct simplitigs from {} kmers", kmer_map.len());
    let mut multi_f = Encoder::new(File::create(output_dir.clone()+"multicolor.fa.zstd").expect("Unable to create file"), 0).unwrap();
    let mut omni_f = Encoder::new(File::create(output_dir.clone()+"omnicolor.fa.zstd").expect("Unable to create file"), 0).unwrap();
    let mut iterator = kmer_map.iter();
    let mut color_nbkmer: HashMap<bitvec::prelude::BitArray<[u8; ARRAY_SIZE]>, usize> = HashMap::new();
    let mut omni: BitArr!(for NB_FILES, in u8) = BitArray::<_>::ZERO;
    let mut is_omni = false;
    for i in 0..NB_FILES{
        omni.set(i, true);
    }
    while let Some((key, curr_cell)) = iterator.next(){
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
            color_nbkmer.entry(curr_cell.0)
                        .and_modify(|nb_kmer| *nb_kmer += simplitig.len()-K+1)
                        .or_insert(simplitig.len()-K+1);
            write_simplitig(&simplitig, &is_omni, &mut multi_f, &mut omni_f, &curr_cell.0);
            is_omni = false;
        }
    }
    multi_f.finish().expect("Error finishing writing in file");
    omni_f.finish().expect("Error finishing writing in file");
    println!("Color to nb kmer map has {} lines", color_nbkmer.len());
    let _ = write_hashmap_to_file(&color_nbkmer);
    
}

fn extend_forward(curr_kmer: &RawKmer<K, u64>, kmer_map:  &HashMap<u64, COLORPAIR>, simplitig: &mut String, color: &BitArray<[u8;1]>) -> bool{
    for succs in curr_kmer.successors(){
        //forward = false;
        if kmer_map.contains_key(&succs.canonical().to_int()){
            let succ_pair = kmer_map.get(&succs.canonical().to_int()).unwrap();
            if succ_pair.0.eq(color) & !succ_pair.1.get() {
                simplitig.push(*succs.to_nucs().last().unwrap() as char);
                succ_pair.1.set(true);
                return true;
            }
        }
    }
    false
}

fn extend_backward(curr_kmer: &RawKmer<K, u64>, kmer_map:  &HashMap<u64, COLORPAIR>, simplitig: &mut String, color: &BitArray<[u8;1]>) -> bool{
    for preds in curr_kmer.predecessors(){
        //backward = false;
        if kmer_map.contains_key(&preds.canonical().to_int()){
            let pred_pair = kmer_map.get(&preds.canonical().to_int()).unwrap();
            if pred_pair.0.eq(color) & !pred_pair.1.get() {
                simplitig.insert(0,*preds.to_nucs().last().unwrap() as char);
                pred_pair.1.set(true);
                return true;
            }
        }
    }
    false
}

fn write_simplitig(simplitig: &String, is_omni: &bool, multi_f: &mut Encoder<'static, File>, omni_f: &mut Encoder<'static, File>, color: &BitArray<[u8;1]>){
    let mut header = String::from(">");
    if *is_omni{
        header.push('\n');
        omni_f.write_all((header).as_bytes()).unwrap();                                                                                                                                                       
        omni_f.write_all((simplitig).as_bytes()).expect("Unable to write data");
        omni_f.write_all(b"\n").unwrap();
    }else{
        for e in color.iter(){
            if *e{
                header.push('1');
            }else{
                header.push('0');
            }
        }
        header.push('\n');
        multi_f.write_all((header).as_bytes()).unwrap();                                                                                                                                                       
        multi_f.write_all((simplitig).as_bytes()).expect("Unable to write data");
        multi_f.write_all(b"\n").unwrap();
    }
}

fn write_hashmap_to_file(map: &HashMap<bitvec::prelude::BitArray<[u8; ARRAY_SIZE]>, usize>) -> io::Result<()> {
    let mut file = File::create("color_kmer_stats.txt")?;
    let mut counter = 0;
    println!("Writing color stats in file...");
    for (_key, value) in map {
        /*let mut key_str = String::new();
        for e in key.iter(){
            if *e{
                key_str.push('1');
            }else{
                key_str.push('0');
            }
            count_files += 1;
            if count_files >= NB_FILES{
                break;
            }
        }*/
        counter += 1;
        //writeln!(file, "{}: {}", key_str, value)?;
        writeln!(file, "{}", value)?;
    }
    println!("Wrote {} lines in stat file.", counter);
    Ok(())
}

fn extract_filename(path: &str) -> Option<&str> {
    // Split the path by '/'
    let parts: Vec<&str> = path.split('/').collect();
    // Get the last element of the path
    if let Some(filename) = parts.last() {
        // Split the filename by '.' and take the first part
        if let Some(idx) = filename.find('.') {
            Some(&filename[..idx])
        } else {
            Some(filename)
        }
    } else {
        None
    }
}

fn num2str(mut k_mer: u64) -> String{
    let mut res = String::from("");
    let mut nuc: u64;
    for _i in 0..K{
        nuc = k_mer%4;
        if nuc == 0{
            res.push('A');
        }else if nuc == 1{
            res.push('C');
        }else if nuc == 2{//bebou
            res.push('G');
        }else if nuc == 3{
            res.push('T');
        }
        k_mer >>= 2;
    }
    res.chars().rev().collect()
}

fn nuc2int(b: &u8) -> Option<u64> {
    match b {
        b'A' | b'C' | b'G' | b'T' => Some(((b / 2) % 4) as u64),
        _ => None,
    }
}

fn rev_comp_str(seq: &str) -> String{
    let mut res = String::new();
    for nuc in seq.chars(){
        if nuc == 'A' {
            res = format!("T{}", res);
        }else if nuc == 'C' {
            res = format!("G{}", res);
        }else if nuc == 'T' {
            res = format!("A{}", res);
        }else{
            res = format!("C{}", res);
        }
    }
    res
}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
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