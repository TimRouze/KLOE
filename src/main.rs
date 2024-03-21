#![allow(dead_code)]

mod kmer;
use clap::Parser;
use flate2::Compression;
use seq_io::fasta::{Reader, Record};
use std::collections::{BTreeSet, BinaryHeap, HashMap, HashSet};
use std::sync::{Arc, Mutex, RwLock};
use std::fs::File;
use std::path::Path;
use std::io::{Write, self, BufRead, stdin};
use std::cmp::{min, max, Reverse};//bebou
use std::env;
use std::time::Instant;
use::rayon::prelude::*;
use flate2::write::GzEncoder;
use bit_vec::BitVec;

use kmer::{Kmer, RawKmer, Base};

use crate::kmer::RevComp;

const K: usize = 31;
pub type KT = u64;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input file (.fasta, .fa)
    input: String,
    /// Output file (defaults to <input>.csv)
    ///#[arg(short, long)]
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
}
// TODO:
    // CREATE FOF INPUT FILE WITH LINK TO COMPRESSED FILES CONTAINING KMERS
        // ASSOCIATE INPUT FILENAME WITH OUTPUT FILENAMES OF COLORS CONTAINING THIS SET
        // HASH TABLE
    // QUERY FUNCTION TO QUERY FILENAME FROM COMPRESSED SET
fn main() {
    let args = Args::parse();
    let input_fof = args.input.as_str();
    println!("FILENAME: {}", input_fof);
    let output_dir = args.out_dir;
    env::set_var("RAYON_NUM_THREADS", args.threads.to_string());
    if let Ok(lines_temp) = read_lines(input_fof){
        let nb_files = lines_temp.count();
        let file = File::open(input_fof).unwrap();
        let reader = io::BufReader::new(file);
        let kmer_map_mutex: Arc<Mutex<HashMap<u64, BitVec>>> = Arc::new(Mutex::new(HashMap::new()));
        let file_number_mutex: RwLock<usize> = RwLock::new(0 as usize);
        reader.lines()
              .par_bridge()
              .for_each(|line| {
                let filename = line.unwrap().clone();
                println!("{}", filename);
                let curr_vec = read_fasta(&filename, &kmer_map_mutex, nb_files, &file_number_mutex);
                let mut file_number = file_number_mutex.write().unwrap();
                *file_number += 1;
                drop(file_number);
              });
        let mut kmer_map = Arc::try_unwrap(kmer_map_mutex).expect("Failed to unwrap Arc").into_inner().expect("Failed to get Mutex");
        let omnicolored_kmers = get_omnicolor(&mut kmer_map, nb_files);
        let to_write_omnicolor = kmers_assemble(&RwLock::new(omnicolored_kmers.clone()));
        let mut f = GzEncoder::new(File::create(output_dir.clone()+"omnicolor.fa.gz").expect("Unable to create file"), Compression::default());
        println!("I write {} simplitigs from {} kmers when starting.", to_write_omnicolor.len(), omnicolored_kmers.len());                                                                                               
        for i in &to_write_omnicolor{           
            f.write(b">\n").unwrap();                                                                                                                                                       
            f.write_all((*i).as_bytes()).expect("Unable to write data");
            f.write(b"\n").unwrap();                                                                                                                         
        } 
        handle_other_colors(&mut kmer_map);
        /*for (color, kmer_map) in kmer_colors.iter_mut(){
            println!("Reconstructing simplitigs for {} color", color);
            let before_kmers_assemble = Instant::now();
            let to_write = kmers_assemble(&Arc::new(RwLock::new(kmer_map.clone())));
            println!("Assembly took: {}ms total.", before_kmers_assemble.elapsed().as_micros());
            let mut path = color.to_string();
            path.push_str("simplitigs.fa.gz");
            let mut f = GzEncoder::new(File::create(output_dir.clone()+&path).expect("Unable to create file"), Compression::default());
            println!("I write {} simplitigs from {} kmers when starting.", to_write.len(), kmer_map.len());                                                                                               
            for i in &to_write{           
                f.write(b">\n").unwrap();                                                                                                                                                       
                f.write_all((*i).as_bytes()).expect("Unable to write data");
                f.write(b"\n").unwrap();                                                                                                                         
            }     
        }*/
    }
}

fn read_fasta(filename: &str, kmer_map_mutex: &Arc<Mutex<HashMap<u64, BitVec>>>, nb_files: usize, file_number: &RwLock<usize>){
    let ( reader, _compression) = niffler::get_reader(Box::new(File::open(filename).unwrap())).unwrap();
    let mut fa_reader = Reader::new(reader);
    let mut kmer = RawKmer::<K, KT>::new();
    while let Some(record) = fa_reader.next(){
        let record = record.expect("Error reading record");
        for (i, nuc) in record.seq().iter().filter_map(KT::from_nuc).enumerate(){
            if i < K - 1{
                kmer = kmer.extend(nuc);
            }else{
                kmer = kmer.append(nuc);
                let canon = kmer.canonical().to_int();
                let mut kmer_map = kmer_map_mutex.lock().unwrap();
                let mut modified = false;
                kmer_map.entry(canon).and_modify(|bit_vec|{
                    bit_vec.set(*file_number.read().unwrap(), true);
                    modified = true;
                }).or_insert(BitVec::from_elem(nb_files, false));
                if !modified{
                    kmer_map.entry(canon).and_modify(|bit_vec| bit_vec.set(*file_number.read().unwrap(), true));
                }
                drop(kmer_map);
            }
        }
    }
}
fn handle_other_colors(kmer_map:  &mut HashMap<u64, BitVec>){
    println!("I will reconstruct simplitigs from {} kmers from different colors", kmer_map.len());
    let nb_kmers = kmer_map.len();
    let mut simplitigs = Vec::new();
    let mut to_remove = Vec::new();
    while let Some(pair) = get_first_key_value(kmer_map){
        let key = pair.0;
        let curr_color = pair.1;
        to_remove.push(*key);
        let mut forward = true;
        let mut backward = true;
        let mut simplitig = num2str(*key);
        while forward{
            let mut curr_kmer = RawKmer::<K, KT>::from_nucs(&simplitig[(simplitig.len()-K)..].as_bytes());
            for succs in curr_kmer.successors(){
                forward = false;
                if kmer_map.contains_key(&succs.canonical().to_int()){
                    let bit_vector = kmer_map.get(&succs.canonical().to_int()).unwrap();
                    if bit_vector.eq(&curr_color) {
                        simplitig.push(*succs.to_nucs().last().unwrap() as char);
                        to_remove.push(succs.canonical().to_int());
                        forward = true;
                        break;
                    }
                }
            }
        }
        while backward{
            let mut curr_kmer = RawKmer::<K, KT>::from_nucs(&simplitig[0..K].as_bytes());
            for preds in curr_kmer.predecessors(){
                backward = false;
                if kmer_map.contains_key(&preds.canonical().to_int()){
                    let bit_vector = kmer_map.get(&preds.canonical().to_int()).unwrap();
                    if bit_vector.eq(&curr_color) {
                        simplitig.insert(0,*preds.to_nucs().last().unwrap() as char);
                        to_remove.push(preds.canonical().to_int());
                        backward = true;
                        break;
                    }
                }
            }
        }
        for e in to_remove.iter(){
            kmer_map.remove(e);
        }
        simplitigs.push(simplitig);
        to_remove.clear();
    }
    println!("I have {} simplitigs from {} kmers when starting.", simplitigs.len(), nb_kmers);   
}
fn get_first_key_value(kmer_map:  &HashMap<u64, BitVec>) -> Option<(&u64, &BitVec)>{
    kmer_map.iter().next()
}

fn get_omnicolor(kmer_map:  &mut HashMap<u64, BitVec>, nb_files: usize) -> BTreeSet<u64>{
    let mut omnicolored_kmer: BTreeSet<u64> = BTreeSet::new();
    let omni_vec = vec![true; nb_files];
    for(key, val) in kmer_map.clone().iter(){
        if val.eq_vec(&omni_vec){
            omnicolored_kmer.insert(*key);
            kmer_map.remove(key);
        }
    }
    omnicolored_kmer
}

fn kmers_assemble(kmer_set_mutex: &RwLock<BTreeSet<u64>>) -> Vec<String>{
    println!("I will reconstruct simplitigs from {} kmers", kmer_set_mutex.read().unwrap().len());
    let mut simplitigs = Vec::new();
    while let Some(seed) =  {
        let mut kmer_set = kmer_set_mutex.write().unwrap();
        kmer_set.pop_last()
    }
    {
        // println!("kmer = {}", num2str(seed));
        let mut simplitig = String::new();
        let before_max_simplitig = Instant::now();
        let extended_simplitig = compute_max_from_kmer(&kmer_set_mutex, &seed);
        // println!("Assembly of 1 simplitig took: {}ms total.", before_max_simplitig.elapsed().as_micros());
        simplitig.push_str(&extended_simplitig);
        simplitigs.push(simplitig);
    }
    simplitigs
}

fn compute_max_from_kmer(available_kmers: &RwLock<BTreeSet<u64>>, seed: &u64) -> String {
    let mut simplitig = num2str(*seed);
    // let before_forward_canon = Instant::now();
    // println!("Simplitig = {}", simplitig);
    let extended_simplitig = extend_forward(available_kmers, &mut simplitig);
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
    let mask : u64 = ((1_u64) << (2*K))-1;
    while extend {
        let mut query = RawKmer::<K, KT>::from_nucs(&simplitig[(simplitig.len()-K)..].as_bytes());
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
            query = RawKmer::<K, KT>::from_nucs(&simplitig[(simplitig.len()-K)..].as_bytes()).canonical();
        }
    }
    simplitig.to_string()
}

/*fn process_fof_parallel(filename: &str, nb_files: usize, size: usize) -> io::Result<bool>{
    let file = File::open(filename)?;
    let header = ">\n";
    let new_line = "\n";
    let reader = io::BufReader::new(file);
    let kmers_mutex: Arc<Mutex<HashMap<String,HashSet<u64>>>> = Arc::new(Mutex::new(HashMap::new()));
    let colors_mutex: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new()));
    // Process lines in parallel using rayon
    reader
        .lines()
        .par_bridge()
        .for_each(|line| {
            let filename = line.unwrap();
            println!("{}", filename);
            handle_fasta(filename, &kmers_mutex, &colors_mutex);
        });
        //TODO RECONSTRUCT UNITIGS
    let colors = Arc::try_unwrap(colors_mutex).expect("Failed to unwrap Arc...").into_inner().expect("Failed to access Mutex...");
    let mut kmer_map = Arc::try_unwrap(kmers_mutex).expect("Failed to unwrap Arc...").into_inner().expect("Failed to access Mutex...");
    colors.into_iter().for_each(|color|{
        println!("{}", color);
        let mut name = color.clone().to_string();
        name.push_str(".fa");
        let mut file = File::create(name).unwrap();
        file.write_all(b">\n").unwrap();
        kmer_map.retain(|key, value| {
            if value == &color{
                file.write_all(num2str(*key).as_bytes()).unwrap();
                file.write_all(b"\n").unwrap();
            }
            !(value == &color)
        });
    });
    Ok(true)
}
fn handle_fasta(filename: String, kmer_mutex: &Arc<Mutex<HashMap<String,HashSet<u64>>>>, color_mutex: &Arc<Mutex<Vec<String>>>){
    let filename_as_vec: Vec<&str> = filename.split('.').collect();
    let temp: Vec<&str> = filename_as_vec[0].split('/').collect();
    let mut id: String = temp.last().unwrap().to_string();
    id.retain(|c| !r#"(),"_"#.contains(c));
    {
    let mut colors = color_mutex.lock().unwrap();
    colors.push(id.to_string());//1);
    }
    /* let mut s=String::new();
    stdin().read_line(&mut s).expect("Did not enter a correct string"); */
    let ( reader, _compression) = niffler::get_reader(Box::new(File::open(&filename).unwrap())).unwrap();
    let mut fa_reader = Reader::new(reader);
    while let Some(record) = fa_reader.next(){
        let record = record.expect("Error reading record");
        for s in record.seq_lines(){
            let seq = String::from_utf8_lossy(s);
            if seq.len() >= 31{
                for _i in 0..(seq.len()-K){
                    let k_mer = str2num(&seq[_i.._i+K]);
                    let canon = canon(k_mer, rev_comp(k_mer));
                    let mut kmer_map = kmer_mutex.lock().unwrap();
                    let mut kmers = kmer_map.entry(id).or_default();
                    if !kmers.contains(&canon){
                        kmers.insert(canon);
                        let mut colors = color_mutex.lock().unwrap();
                        colors.into_iter().for_each(|color|{
                            
                        });
                        if !colors.contains(&format!("{}_{}", min(&color_id, &id), max(&color_id, &id))){//&3){//
                            colors.push(format!("{}_{}", min(&color_id, &id), max(&color_id, &id)));//3);
                        }
                        /* if color_id.len() > 13{
                            println!("New id = {}\ncurr id = {}", color_id, &id);
                            let mut s=String::new();
                            stdin().read_line(&mut s).expect("Did not enter a correct string");
                        } */
                        /* println!("New id = {}\ncurr id = {}", color_id, id);
                        let mut s=String::new();
                        stdin().read_line(&mut s).expect("Did not enter a correct string"); */
                    }
                }
            }
        }
    }
}*/

/*
fn write_output(hist: Vec<u64>, nb_files: usize) -> Result<(), Box<dyn Error>>{

    let mut wtr = Writer::from_path("out.csv")?;
    let header: Vec<u16> = (1..(nb_files+1) as u16).collect();
    wtr.serialize(header)?;
    wtr.serialize(&hist[1..(nb_files+1)])?;
    wtr.flush()?;
    Ok(())
}
*/
//bebou

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

fn canon(k_mer1: u64, k_mer2:u64) -> u64{
    min(k_mer1, k_mer2)
}

fn str2num(k_mer: &str) -> u64{
    let mut res : u64 = 0;
    for character in k_mer.chars(){
        res <<=2;
        res += (character as u64/2)%4;
    }
    res
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

fn rev_comp(k_mer : u64) -> u64 {
    let mut res = k_mer.reverse_bits();
    res = (res >> 1 & 0x5555_5555_5555_5555) | (res & 0x5555_5555_5555_5555) << 1;
    res ^= 0xAAAA_AAAA_AAAA_AAAA;
    res >> (2 * (32 - K))
}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}
#[test]
fn test_num2str(){
    let kmer = "AGCTTTTCATTCTGACTGCAACGGGCAATAT";
    let num = str2num(kmer);
    assert_eq!(num2str(num), kmer);
}
#[test]
fn test_rc_64() {
    let kmer = "AGCTTTTCATTCTGACTGCAACGGGCAATAT";
    let revcomp = rev_comp(str2num(kmer));
    let res = num2str(revcomp);
    assert_eq!(res, "ATATTGCCCGTTGCAGTCAGAATGAAAAGCT");
}

#[test]
fn test_canon(){
    let seq = "AGCTTTTCATTCTGACTGCAACGGGCAATAT";
    let kmer = str2num(seq);
    let revcomp = rev_comp(kmer);
    let k_mer_canon = canon(kmer, revcomp);
    assert_eq!(num2str(k_mer_canon), "ATATTGCCCGTTGCAGTCAGAATGAAAAGCT");
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