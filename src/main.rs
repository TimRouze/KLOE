#![allow(dead_code)]

mod kmer;
use clap::Parser;
use flate2::Compression;
use seq_io::fasta::{Reader, Record};
use std::collections::{BTreeSet, HashSet};
use std::sync::{Arc, Mutex, RwLock};
use std::fs::File;
use std::path::Path;
use std::io::{Write, self, BufRead, Stdin};
use std::cmp::min;//bebou
use std::env;
use std::time::Instant;
use::rayon::prelude::*;
use flate2::write::GzEncoder;
use bitvec::prelude::*;
use kmer::{Kmer, RawKmer, Base};
use std::cell::Cell;
use hashbrown::HashMap;


#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    ///output: Option<String>,
    /// Number of threads (defaults to all available threads)
    #[arg(short, long, default_value_t = 1)]
    threads: usize,
    ///Output directory
    #[arg(short, long, default_value_t = String::from(""))]
    out_dir: String,
    /// Number of hashes used in Bloom filters
    #[arg(short = 'H', long, default_value_t = 1)]
    hashes: usize,
    /// Seed used for hash functions
    #[arg(short, long, default_value_t = 101010)]
    seed: u64,
    ///Number of expected elements
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

// TODO:
    // REPLACE BIT VECTOR BY U64 SETTING BITS TO 1 OR 0
fn main() {
    let args = Args::parse();
    let nb_elem = args.nb_elem;
    println!("FILENAME: {}", INPUT_FOF);
    let output_dir = args.out_dir;
    env::set_var("RAYON_NUM_THREADS", args.threads.to_string());
    let file = File::open(INPUT_FOF).unwrap();
    let reader = io::BufReader::new(file);
    let kmer_map_mutex: Arc<Mutex<HashMap<u64, COLORPAIR>>> = Arc::new(Mutex::new(HashMap::with_capacity(nb_elem)));
    //let set_mutex: Arc<Mutex<HashSet<u64>>> = Arc::new(Mutex::new(HashSet::new()));
    let file_number_mutex: RwLock<usize> = RwLock::new(0 as usize);
    reader.lines()
            .par_bridge()
            .for_each(|line| {
            let filename = line.unwrap().clone();
            println!("{}", filename);
            read_fasta(&filename, &kmer_map_mutex, &file_number_mutex, &nb_elem);
            let mut file_number = file_number_mutex.write().unwrap();
            *file_number += 1;
            drop(file_number);
            });
    let mut kmer_map = Arc::try_unwrap(kmer_map_mutex).expect("Failed to unwrap Arc").into_inner().expect("Failed to get Mutex");
    println!("NB KMER = {}", kmer_map.len());
    kmer_map.shrink_to_fit();
    let omnicolored_kmers = get_omnicolor(&mut kmer_map);
    kmers_assemble(&RwLock::new(omnicolored_kmers.clone()), &output_dir);
    handle_other_colors(&mut kmer_map, &output_dir);
}

fn read_fasta(filename: &str, kmer_map_mutex: &Arc<Mutex<HashMap<u64, COLORPAIR>>>, file_number: &RwLock<usize>, nb_elem: &usize){
    let ( reader, _compression) = niffler::get_reader(Box::new(File::open(filename).unwrap())).unwrap();
    let mut fa_reader = Reader::new(reader);
    let mut kmer = RawKmer::<K, KT>::new();
    //let mut to_add = HashSet::new();
    while let Some(record) = fa_reader.next(){
        let record = record.expect("Error reading record");
        for (i, nuc) in record.seq().iter().filter_map(KT::from_nuc).enumerate(){
            if i < K - 1{
                kmer = kmer.extend(nuc);
            }else{
                kmer = kmer.append(nuc);
                let canon = kmer.canonical().to_int();
                //to_add.insert(canon);
                let mut kmer_map = kmer_map_mutex.lock().unwrap();
                kmer_map.entry(canon).and_modify(|pair|{
                    pair.0.set(*file_number.read().unwrap(), true);
                }).or_insert({let bv: BitArr!(for NB_FILES, in u8) = BitArray::<_>::ZERO;
                    (bv, Cell::new(false))
                });
                if kmer_map.len() >= (80/100)*nb_elem{
                    kmer_map.reserve((30/100)*nb_elem);
                }
                drop(kmer_map);
                /*if to_add.len() >= 10000{
                    //let mut kmer_map = kmer_map_mutex.lock().unwrap();
                    let mut kmer_set = set_mutex.lock().unwrap();
                    for e in to_add.iter(){
                        /*kmer_map.entry(*e).and_modify(|pair|{
                            pair.0.set(*file_number.read().unwrap(), true);
                        }).or_insert({let mut bv = bitvec![u8, Msb0;];
                            bv.extend(vec![false; nb_files].iter());
                            bv.set(*file_number.read().unwrap(), true);
                            (bv, Cell::new(false))
                        });*/
                        kmer_set.insert(*e);
                    }
                    //drop(kmer_map);
                    drop(kmer_set);
                    to_add.clear();
                }*/
            }
            
        }
    }
}
fn handle_other_colors(kmer_map:  &mut HashMap<u64, COLORPAIR>, output_dir: &String){
    println!("I will reconstruct simplitigs from {} kmers from different colors", kmer_map.len());
    let mut f = GzEncoder::new(File::create(output_dir.clone()+"multicolor.fa.gz").expect("Unable to create file"), Compression::default());
    //let mut to_remove = Vec::new();
    let mut iterator = kmer_map.iter();
    while let Some((key, curr_cell)) = iterator.next(){//get_first_key_value(kmer_map){
        if !curr_cell.1.get(){
            curr_cell.1.set(true);
            //println!("{}", key);
            //to_remove.push(*key);
            let mut forward = true;
            let mut backward = true;
            let mut simplitig = num2str(*key);
            while forward{
                let curr_kmer = RawKmer::<K, KT>::from_nucs(&simplitig[(simplitig.len()-K)..].as_bytes());
                for succs in curr_kmer.successors(){
                    forward = false;
                    if kmer_map.contains_key(&succs.canonical().to_int()){
                        let succ_pair = kmer_map.get(&succs.canonical().to_int()).unwrap();
                        if succ_pair.0.eq(&curr_cell.0) & !succ_pair.1.get() {
                            //println!("{}", succ_pair.1.get());
                            simplitig.push(*succs.to_nucs().last().unwrap() as char);
                            succ_pair.1.set(true);
                            //kmer_map.entry(succs.canonical().to_int()).and_modify(|cel| cel.1.set(true));
                            //println!("{}", succ_pair.1.get());
                            //let mut ret = String::new();
                            //io::stdin().read_line(&mut ret).expect("Error reading input");
                            //to_remove.push(succs.canonical().to_int());
                            forward = true;
                            break;
                        }
                        //println!("{}", succ_pair.1.get());
                    }
                }
            }
            while backward{
                let curr_kmer = RawKmer::<K, KT>::from_nucs(&simplitig[0..K].as_bytes());
                for preds in curr_kmer.predecessors(){
                    backward = false;
                    if kmer_map.contains_key(&preds.canonical().to_int()){
                        let pred_pair = kmer_map.get(&preds.canonical().to_int()).unwrap();
                        if pred_pair.0.eq(&curr_cell.0) & !pred_pair.1.get() {
                            simplitig.insert(0,*preds.to_nucs().last().unwrap() as char);
                            pred_pair.1.set(true);
                            //to_remove.push(preds.canonical().to_int());
                            backward = true;
                            break;
                        }
                        //println!("{}", pred_pair.1.get());
                    }
                }
            }
            //println!("writing simplitig");
            let mut header = String::from(">");
            for e in curr_cell.0.iter(){
                if *e{
                    header.push('1');
                }else{
                    header.push('0');
                }
            }
            header.push('\n');
            f.write((header).as_bytes()).unwrap();                                                                                                                                                       
            f.write((simplitig).as_bytes()).expect("Unable to write data");
            f.write(b"\n").unwrap();
            /*for e in to_remove.iter(){
                kmer_map.remove(e);
            }                
            to_remove.clear();*/
        }
    }
    
}

fn get_omnicolor(kmer_map:  &mut HashMap<u64, COLORPAIR>) -> BTreeSet<u64>{
    println!("getting omnicolored kmers");
    let mut omnicolored_kmer: BTreeSet<u64> = BTreeSet::new();
    let mut to_remove = Vec::new();
    let mut omni: BitArr!(for 8, in u8) = BitArray::<_>::ZERO;
    for i in 0..NB_FILES{
        omni.set(i, true);
    }
    let mut iterator = kmer_map.iter();
    while let Some(pair) = iterator.next(){//get_first_key_value(kmer_map) {
        //println!("{}", pair.0);
        if pair.1.0.eq(&omni){
            //println!("coucou");
            to_remove.push(*pair.0);
            omnicolored_kmer.insert(*pair.0);
        }
    }
    for e in to_remove.iter(){
        kmer_map.remove(e);
    }
    omnicolored_kmer
}

fn kmers_assemble(kmer_set_mutex: &RwLock<BTreeSet<u64>>, output_dir: &String){
    let mut f = GzEncoder::new(File::create(output_dir.clone()+"omnicolor.fa.gz").expect("Unable to create file"), Compression::default());
    println!("I will reconstruct simplitigs from {} kmers", kmer_set_mutex.read().unwrap().len());
    // let mut simplitigs = Vec::new();
    while let Some(seed) =  {
        let mut kmer_set = kmer_set_mutex.write().unwrap();
        kmer_set.pop_last()
    }
    {
        // println!("kmer = {}", num2str(seed));
        let mut simplitig = String::new();
        // let before_max_simplitig = Instant::now();
        let extended_simplitig = compute_max_from_kmer(&kmer_set_mutex, &seed);
        // println!("Assembly of 1 simplitig took: {}ms total.", before_max_simplitig.elapsed().as_micros());
        simplitig.push_str(&extended_simplitig);
        // simplitigs.push(simplitig);
        f.write(b">\n").unwrap();                                                                                                                                                       
        f.write_all((*simplitig).as_bytes()).expect("Unable to write data");
        f.write(b"\n").unwrap();
    }
    // simplitigs
}

fn compute_max_from_kmer(available_kmers: &RwLock<BTreeSet<u64>>, seed: &u64) -> String {
    let mut simplitig = num2str(*seed);
    // let before_forward_canon = Instant::now();
    // println!("Simplitig = {}", simplitig);
    /*let extended_simplitig = */extend_forward(available_kmers, &mut simplitig);
    // println!("Extended simplitig = {}", extended_simplitig);
    // println!("Forward extension canonical left took: {}ms total.", before_forward_canon.elapsed().as_micros());
    let mut reversed_simplitig = rev_comp_str(&simplitig);
    // let before_forward_revcomp = Instant::now();
    simplitig = extend_forward(available_kmers, &mut reversed_simplitig);
    // println!("Final simplitig = {}", simplitig);
    // println!("Forward extension reverse way took: {}ms total.", before_forward_revcomp.elapsed().as_micros());
    simplitig
}

fn extend_forward(available_kmers_mutex: &RwLock<BTreeSet<u64>>, simplitig: &mut String) -> String {
    let mut extend = true;
    let mut query = RawKmer::<K, KT>::from_nucs(&simplitig[(simplitig.len()-K)..].as_bytes());
    while extend {
        // println!("Query kmer = {}", String::from_utf8(query.to_nucs().to_ascii_uppercase()).unwrap());
        // let mut query= RawKmer::<K, KT>::from_nucs(&simplitig[(simplitig.len()-K)..].as_bytes()).canonical();
        extend = false;
        for succs in query.successors(){
            let mut available_kmers = available_kmers_mutex.write().unwrap();
            // println!("Candidate successor = {}", String::from_utf8(succs.to_nucs().to_ascii_uppercase()).unwrap());
            if available_kmers.contains(&succs.canonical().to_int()){
                // println!("Successor = {}", String::from_utf8(succs.to_nucs().to_ascii_uppercase()).unwrap());
                simplitig.push(*succs.to_nucs().last().unwrap() as char);
                available_kmers.remove(&succs.canonical().to_int());
                extend = true;
                break;
            }
        }
        query = RawKmer::<K, KT>::from_nucs(&simplitig[(simplitig.len()-K)..].as_bytes()).canonical();
    }
    simplitig.to_string()
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