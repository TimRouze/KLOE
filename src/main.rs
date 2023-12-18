#![allow(dead_code)]

use clap::Parser;
use seq_io::fasta::{Reader, Record};
use std::collections::HashMap;
use std::collections::btree_map::Keys;
use std::sync::{Arc, Mutex, RwLock};
use std::fs::File;
use std::path::Path;
use std::io::{Write, self, BufRead, stdin};
use std::cmp::{min, max};//bebou
use std::env;
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

    env::set_var("RAYON_NUM_THREADS", args.threads.to_string());
    if let Ok(lines_temp) = read_lines(input_fof){
        let nb_files = lines_temp.count();
        match process_fof_parallel(input_fof, nb_files, size) {
            Ok(hist_mutex) => {
                println!("All {} files have been read...\nWriting output...", nb_files);
                //let hist = Arc::try_unwrap(hist_mutex).expect("Failed to unnwrap Arc").into_inner().expect("Failed to get Mutex");
                //write_output(hist, nb_files).unwrap();
            }
            Err(err) => eprintln!("Error reading or processing file: {}", err),
        }
    }
//    var |-ma-variable-est-un-kebab-|
//    var ma_variable_est_un_serpent
//    var maVariableEstUnChameaubebou
}

fn process_fof_parallel(filename: &str, nb_files: usize, size: usize) -> io::Result<bool>{
    let file = File::open(filename)?;
    let header = ">\n";
    let new_line = "\n";
    let reader = io::BufReader::new(file);
    let kmers_mutex: Arc<Mutex<HashMap<u64,String>>> = Arc::new(Mutex::new(HashMap::new()));
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
    let colors = Arc::try_unwrap(colors_mutex).expect("Failed to unwrap Arc...").into_inner().expect("Failed to access Mutex...");
    let mut kmer_map = Arc::try_unwrap(kmers_mutex).expect("Failed to unwrap Arc...").into_inner().expect("Failed to access Mutex...");
    colors.into_iter().for_each(|color|{
        println!("{}", color);
        let mut name = color.clone();
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
//  TODO ADD MERGE PROCESS OR HANDLE SIMETRICAL KEYS
fn handle_fasta(filename: String, kmer_mutex: &Arc<Mutex<HashMap<u64,String>>>, color_mutex: &Arc<Mutex<Vec<String>>>){
    let filename_as_vec: Vec<&str> = filename.split('.').collect();
    let temp: Vec<&str> = filename_as_vec[0].split('/').collect();
    let mut id: String = temp.last().unwrap().to_string();
    id.retain(|c| !r#"(),"_"#.contains(c));
    {
    let mut colors = color_mutex.lock().unwrap();
    colors.push(id.to_string());
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
                    if let Some(color_id) = kmer_map.insert(canon, id.to_string()){
                        if !color_id.contains(&id){
                            //println!("prev id = {}", color_id);
                            let new_id = format!("{}_{}", color_id, id);
                            kmer_map.insert(canon, new_id);
                            let mut colors = color_mutex.lock().unwrap();
                            if !colors.contains(&format!("{}_{}", min(&color_id, &id), max(&color_id, &id))){
                                colors.push(format!("{}_{}", min(&color_id, &id), max(&color_id, &id)));
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
    }
}
/* fn handle_fasta(filename: String, mutex: &Arc<RwLock<HashMap<String,Vec<u64>>>>){
    let filename_as_vec: Vec<&str> = filename.split('.').collect();
    let temp: Vec<&str> = filename_as_vec[0].split('/').collect();
    let id = temp.last().unwrap();
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
                    let map_reader = mutex.read().unwrap();
                    //let vector = kmer_map.entry(id.to_string()).or_insert(vec![]);
                    //vector.push(canon);
                    if map_reader.keys().len() > 1{
                        let mut tmp = id.to_string();//String::new();
                        let clone_dunno_why = map_reader.clone();
                        clone_dunno_why.into_iter().for_each(|(key, vector)|{
                            /* println!("{}\n{}\n{}", key, num2str(canon), vector.contains(&canon));
                            let mut s=String::new();
                            stdin().read_line(&mut s).expect("Did not enter a correct string"); */
                            if vector.contains(&canon){
                                tmp += &key;
                                let mut map_writer = mutex.write().unwrap();
                                let curr_vec = map_writer.entry(key.to_string()).or_default();
                                curr_vec.swap_remove(curr_vec.iter().position(|x| *x == canon).expect("kmer not found"));
                            }
                        });
                        let mut map_writer = mutex.write().unwrap();
                        let kmer_vec = map_writer.entry(tmp).or_insert(vec![]);
                        kmer_vec.push(canon);
                        /* let display = kmer_map.clone();
                        display.into_iter().for_each(|(key,val)|{
                            println!("Key = {}", key);
                            println!("{}", val.len());
                            let mut s=String::new();
                            stdin().read_line(&mut s).expect("Did not enter a correct string");
                        }); */
                    }else{
                        let mut map_writer = mutex.write().unwrap();
                        let vector = map_writer.entry(id.to_string()).or_insert(vec![]);
                        vector.push(canon);
                    }
                }
            }
        }
    }
} */
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