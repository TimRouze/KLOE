#![allow(dead_code)]

use clap::Parser;
use flate2::Compression;
use seq_io::fasta::{Reader, Record};
use std::collections::{HashMap, BinaryHeap};
use std::hash::Hash;
use std::sync::{Arc, Mutex, RwLock};
use std::fs::File;
use std::path::Path;
use std::io::{Write, self, BufRead, stdin};
use std::cmp::{min, max, Reverse};//bebou
use std::env;
use std::time::Instant;
use::rayon::prelude::*;
use flate2::write::GzEncoder;


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
    /// Memory (in GB) allocated to Bloom filters (defaults to 4GB)
    #[arg(short, long, default_value_t = 4)]
    memory: usize,
    /// Number of hashes used in Bloom filters
    #[arg(short = 'H', long, default_value_t = 1)]
    hashes: usize,
    /// Seed used for hash functions
    #[arg(short, long, default_value_t = 101010)]
    seed: u64,
    /// Modimizer ?
    #[arg(short = 'M', long)]
    modimizer: bool,
}
fn main() {
    let args = Args::parse();
    let input_fof = args.input.as_str();
    let size =  args.memory * 1_000_000_000;
    println!("FILENAME: {}", input_fof);
    env::set_var("RAYON_NUM_THREADS", args.threads.to_string());
    if let Ok(lines_temp) = read_lines(input_fof){
        let nb_files = lines_temp.count();
        let file = File::open(input_fof).unwrap();
        let reader = io::BufReader::new(file);
        let vec_of_kmer_vec_mutex: Arc<Mutex<Vec<Vec<u64>>>> = Arc::new(Mutex::new(Vec::new()));
        let color_names_mutex = Arc::new(Mutex::new(Vec::new()));
        reader.lines()
              .par_bridge()
              .for_each(|line| {
                let filename = line.unwrap();
                println!("{}", filename);
                let temp: Vec<&str> = filename.split(&['/', '.']).collect();
                let mut color_names = color_names_mutex.lock().unwrap();
                color_names.push(temp[temp.len()-2].to_string());
                let curr_vec = read_fasta(&filename);
                let mut vec_of_kmer_vec = vec_of_kmer_vec_mutex.lock().unwrap();
                vec_of_kmer_vec.push(curr_vec);
              });
        let vec_of_kmer_vecs = Arc::try_unwrap(vec_of_kmer_vec_mutex).expect("Failed to unwrap Arc").into_inner().expect("Failed to get Mutex");
        let color_names = Arc::try_unwrap(color_names_mutex).expect("Failed to unwrap Arc").into_inner().expect("Failed to get Mutex");
        let mut kmer_colors = process_vectors(&vec_of_kmer_vecs, &color_names);
        for (color, kmer_map) in kmer_colors.iter_mut(){
            println!("Reconstructing simplitigs for {} color", color);
            let before_kmers_assemble = Instant::now();
            let to_write = kmers_assemble(&Arc::new(RwLock::new(kmer_map.clone())));
            println!("Assembly took: {}ms total.", before_kmers_assemble.elapsed().as_micros());
            let mut path = color.to_string();
            path.push_str("simplitigs.fa.gz");
            let mut f = GzEncoder::new(File::create(path).expect("Unable to create file"), Compression::default());
            println!("I write {} simplitigs from {} kmers when starting.", to_write.len(), kmer_map.len());                                                                                               
            for i in &to_write{           
                f.write(b">\n").unwrap();                                                                                                                                                       
                f.write_all((*i).as_bytes()).expect("Unable to write data");
                f.write(b"\n").unwrap();                                                                                                                         
            }     
        }// TO DO DEBUG ARC RWLOCK DE MAP DE KMERS
    }
}

fn read_fasta(filename: &str) -> Vec<u64>{
    let ( reader, _compression) = niffler::get_reader(Box::new(File::open(filename).unwrap())).unwrap();
    let mut fa_reader = Reader::new(reader);
    let mut kmer_vec = Vec::new();
    while let Some(record) = fa_reader.next(){
        let record = record.expect("Error reading record");
        for s in record.seq_lines(){
            let seq = String::from_utf8_lossy(s);
            if seq.len() >= 31{
                for _i in 0..(seq.len()-K){
                    let k_mer = str2num(&seq[_i.._i+K]);
                    let canon = canon(k_mer, rev_comp(k_mer));
                    kmer_vec.push(canon);
                }
            }
        }
    }
    kmer_vec.par_sort();
    kmer_vec.dedup();
    kmer_vec
}

fn process_vectors(vec_of_kmer_vecs:  &[Vec<u64>], color_names: &[String]) -> HashMap<String, HashMap<u64, bool>>{
    println!("NB lists = {}", vec_of_kmer_vecs.len());
    let mut indices = vec![0; vec_of_kmer_vecs.len()];
    let mut kmer_heap = BinaryHeap::new();
    let mut kmer_color: HashMap<String, HashMap<u64, bool>> = HashMap::new();
    let mut run = true;
    let mut counter = 0;
    vec_of_kmer_vecs.iter().for_each(|kmer_vec|{
        kmer_heap.push(Reverse(kmer_vec[0]));
    });
    while run{
        if let Some(min_kmer) = kmer_heap.pop(){
            let mut new_color: String = String::from("");
            for i in 0..vec_of_kmer_vecs.len(){
                if indices[i] < vec_of_kmer_vecs[i].len() && min_kmer.0 == vec_of_kmer_vecs[i][indices[i]]{
                    new_color += &(color_names[i].to_string()+"_");
                    indices[i] += 1;
                    counter += 1;
                    if indices[i] < vec_of_kmer_vecs[i].len(){
                        kmer_heap.push(Reverse(vec_of_kmer_vecs[i][indices[i]]));
                    }
                }
            }
            if !new_color.is_empty(){
                kmer_color
                .entry(new_color)
                .and_modify(|kmer_map|{ kmer_map.insert(min_kmer.0, false);})
                .or_insert(HashMap::from([(min_kmer.0, false)]));
            }
        }else{
            run = false;
        }
    }
    println!("NB kmer in 1 = {}", vec_of_kmer_vecs[0].len());
    println!("NB kmer in 2 = {}", vec_of_kmer_vecs[1].len());
    println!("I have seen {} matches!", counter);
    kmer_color
}

fn kmers_assemble(kmers_mutex: &Arc<RwLock<HashMap<u64, bool>>>) -> Vec<String>{
    println!("I will reconstruct simplitigs from {} kmers", kmers_mutex.read().unwrap().len());
    let mut simplitigs = Vec::new();
    loop {
        let seed = {
            let kmers = kmers_mutex.write().unwrap();
            kmers.iter().find(|&(_, &seen)| !seen).map(|(&seed, _)| seed).clone()
        };
        if let Some(seed) = seed{
            let mut kmers = kmers_mutex.write().unwrap();
            if let Some(value) = kmers.get_mut(&seed) {
                *value = true;
            }
            drop(kmers);
            let mut simplitig = String::new();
            let before_max_simplitig = Instant::now();
            let extended_simplitig = compute_max_from_kmer(kmers_mutex, &seed);
            println!("Assembly of 1 simplitig took: {}ms total.", before_max_simplitig.elapsed().as_micros());
            simplitig.push_str(&extended_simplitig);
            simplitigs.push(simplitig);

        }else {
            break;
        }
    }
    simplitigs
}

fn compute_max_from_kmer(available_kmers_mutex: &Arc<RwLock<HashMap<u64, bool>>>, seed: &u64) -> String {
    let mut simplitig = num2str(*seed);
    let before_forward_canon = Instant::now();
    let extended_simplitig = extend_forward(available_kmers_mutex, &mut simplitig);
    println!("Forward extension canonical left took: {}ms total.", before_forward_canon.elapsed().as_micros());
    let mut reversed_simplitig = rev_comp_str(&extended_simplitig);
    let before_forward_revcomp = Instant::now();
    simplitig = extend_forward(available_kmers_mutex, &mut reversed_simplitig);
    println!("Forward extension reverse way took: {}ms total.", before_forward_revcomp.elapsed().as_micros());
    simplitig
}

fn extend_forward(available_kmers_mutex: &Arc<RwLock<HashMap<u64,bool>>>, simplitig: &mut String) -> String {
    let mut extend = true;
    let mask : u64 = ((1_u64) << (2*K))-1;
    while extend {
        let mut query = str2num(&simplitig[(simplitig.len()-K)..]);
        extend = false; // Set extend to false by default
        for x in [b'A', b'C', b'G', b'T'].iter() {
            query <<= 2;
            query += nuc2int(x).unwrap();
            query &= mask;
            query = canon(query, rev_comp(query));
            let mut available_kmers = available_kmers_mutex.write().unwrap();
            if let Some(false) = available_kmers.get(&query){
                simplitig.push(*x as char);
                available_kmers.entry(query).and_modify(|seen| *seen = true);
                extend = true; // Set extend to true if a match is found
                break;
            }
            query = str2num(&simplitig[(simplitig.len()-K)..]);
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
            res.push('T');
        }else if nuc == 3{
            res.push('G');
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