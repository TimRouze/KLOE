#![allow(dead_code)]
// WHILE NEEDLETAIL DOES NOT HANDLE ZST.
extern crate needletail;

mod utils;
mod kmer;
mod decompress;
mod stats;
use clap::Parser;
use num_traits::ToPrimitive;
use stats::compute_stats;
use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use std::fs::{remove_file, File};
use std::ffi::OsStr;
use std::path::{PathBuf, Path};
use std::io::{self, BufRead, Read, Seek, SeekFrom, Stdin, Write, BufReader, BufWriter};
use std::env;
use std::time::Instant;
use::rayon::prelude::*;
use bitvec::prelude::*;
use kmer::{Kmer, RawKmer, Base};
use std::cell::Cell;
use seq_io::fasta::Reader;
use indexmap::IndexMap;
use zstd::stream::write::Encoder;
use zstd::stream::read::Decoder;

use crate::utils::{num2str, str2num, find_min};

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
const BLOCK_SIZE: usize = 1 << (12 - 3);
const SHARD_AMOUNT: usize = 1024;
const M: u8 = 7;

fn main() {
    let args = Args::parse();
    let nb_elem = args.nb_elem;
    println!("FILENAME: {}", INPUT_FOF);
    let output_dir = args.out_dir;
    env::set_var("RAYON_NUM_THREADS", args.threads.to_string());
    let file = File::open(INPUT_FOF).unwrap();
    let reader = io::BufReader::new(file);
    //let mut omni_kmer_map: Arc<Mutex<HashMap<KT, Cell<bool>>>> = Arc::new(Mutex::new(HashMap::with_capacity(nb_elem/2)));
    //let mut multi_kmer_map: Arc<Mutex<HashMap<KT, COLORPAIR>>> = Arc::new(Mutex::new(HashMap::with_capacity(nb_elem/2)));
    let mut omni_vec_o_maps: Vec<Arc<Mutex<HashMap<KT, Cell<bool>>>>> = Vec::new();
    let mut multi_vec_o_maps: Vec<Arc<Mutex<HashMap<KT, COLORPAIR>>>> = Vec::new();
    for _i in 0..SHARD_AMOUNT{
        omni_vec_o_maps.push(Arc::new(Mutex::new(HashMap::with_capacity(nb_elem/SHARD_AMOUNT))));
        multi_vec_o_maps.push(Arc::new(Mutex::new(HashMap::with_capacity(nb_elem/SHARD_AMOUNT))));
    }
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
            let now = Instant::now();
            let mut tot_nb_kmer: u64 = 0;
            (0..NB_FILES).for_each(|file_number|{
                println!("{}", filenames.get(file_number).unwrap());
                let nb_kmer = read_fasta(&filenames.get(file_number).unwrap(), &mut multi_vec_o_maps, &mut omni_vec_o_maps, file_number, &nb_elem);
                tot_nb_kmer += Arc::try_unwrap(nb_kmer).unwrap().into_inner().unwrap();
            });
            let elapsed = now.elapsed();
            let mut filled_omni_vec = Vec::new();
            let mut filled_multi_vec = Vec::new();
            for _i in 0..SHARD_AMOUNT{
                let omni_elem = omni_vec_o_maps.remove(0);
                filled_omni_vec.push(Arc::try_unwrap(omni_elem).expect("Failed to get multi map mutex").into_inner().expect("Failed to access inside multi map mutex"));
                let multi_elem = multi_vec_o_maps.remove(0);
                filled_multi_vec.push(Arc::try_unwrap(multi_elem).expect("Failed to get multi map mutex").into_inner().expect("Failed to access inside multi map mutex"));
            }
            //let mut filled_omni_kmer_map = Arc::try_unwrap(omni_kmer_map).expect("Failed to unwrap Arc").into_inner().expect("Failed to get Mutex");
            //let mut filled_multi_kmer_map = Arc::try_unwrap(multi_kmer_map).expect("Failed to unwrap Arc").into_inner().expect("Failed to get Mutex");
            println!("Filling map took: {:.2?} seconds.", elapsed);
            println!("COMPUTED NB KMER: {}", tot_nb_kmer);
            //println!("NB KMER = {}", filled_omni_kmer_map.len()+filled_multi_kmer_map.len());
            let now = Instant::now();
            //let mut color_simplitig = HashMap::new();
            /*rayon::scope(|s| {
                s.spawn(|_| compute_omnicolored_simplitigs(&mut filled_omni_vec, &output_dir));

                s.spawn(|_| {
                    //let color_simplitig_mutex = compute_multicolored_simplitigs(&mut filled_multi_kmer_map, &output_dir);
                    //color_simplitig = Arc::try_unwrap(color_simplitig_mutex).unwrap().into_inner().unwrap();
                    color_simplitig = compute_multicolored_simplitigs(&mut filled_multi_vec, &output_dir);
                });
            });*/
            compute_omnicolored_simplitigs(&mut filled_omni_vec, &output_dir);
            let mut color_simplitig = compute_multicolored_simplitigs(&mut filled_multi_vec, &output_dir);
            let elapsed = now.elapsed();
            println!("Construction of simplitigs took: {:.2?} seconds.", elapsed);
            let now = Instant::now();
            sort_simplitigs(&mut color_simplitig, &output_dir);
            let elapsed = now.elapsed();
            println!("Sorting Simplitigs took: {:.2?} seconds.", elapsed);
        }else if do_decompress == "stats"{
            let multi_file = args.multicolor_file;
            let omni_file = args.omnicolor_file;
            let k = K;
            if multi_file != "" && omni_file != ""{
                compute_stats(&omni_file, &multi_file, &output_dir, &k);
            }else{
                println!("Error, multicolor and/or omnicolor file(s) are mandatory");
            }
        }/*else if do_decompress == "test"{
            let kmer_map_test_mutex: Arc<Mutex<HashMap<KT, u32>>> = Arc::new(Mutex::new(HashMap::with_capacity(nb_elem)));
            let now = Instant::now();
            (0..NB_FILES).into_par_iter().for_each(|file_number|{
                let mut kmer_map_mutex = Arc::clone(&kmer_map_test_mutex);
                println!("{}", filenames.get(file_number).unwrap());
                test_fill_map(&filenames.get(file_number).unwrap(), &mut kmer_map_mutex, file_number, &nb_elem);
            });
            let elapsed = now.elapsed();
            println!("Filling map took: {:.2?} seconds.", elapsed);
        }*/
    }else {
        println!("Wrong positional arguments given. Values are 'compress' or 'decompress'");
        println!("Ex: if compression: I=my/fof.txt cargo r -r -- compress -o out_dir/ -t 12");
        println!("Ex: if decompression: I=my/fof.txt cargo r -r -- decompress -omnicolor-file out_dir/omnicolor.fa.zstd --multicolor-file out_dir/multicolor.fa.zstd -t 12");
    }
}

/*fn test_fill_map(filename: &str, kmer_map_mutex: &Arc<Mutex<HashMap<KT, u32>>>, file_number: usize, nb_elem: &usize){
    let mut fa_reader = parse_fastx_file(filename).expect("Error while opening file");
    let mut counter : u64 = 0;
    let mut counter_insert : u64 = 0;
    let mut counter_modify : u64 = 0;
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
            kmer_map.entry(curr_kmer.to_int()).and_modify(|cpt|{
                counter_modify += 1;
                *cpt = *cpt+1;
            }).or_insert_with(|| {
                counter_insert += 1;
                let mut cpt : u32 = 0;
                cpt
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
}*/


fn get_extension_from_filename(filename: &str) -> Option<&str> {
    Path::new(filename)
        .extension()
        .and_then(OsStr::to_str)
}

/*
READS FASTA FILES.
CREATES K-MER -> COLOR MAP.
*/
fn read_fasta(filename: &str, multi_kmer_map: &mut Vec<Arc<Mutex<HashMap<KT, COLORPAIR>>>>, omni_kmer_map: &mut 
    Vec<Arc<Mutex<HashMap<KT, Cell<bool>>>>>, file_number: usize, nb_elem: &usize) -> Arc<Mutex<u64>>{
    //let mut fa_reader = parse_fastx_file(filename).expect("Error while opening file");
    //let mut kmer_buffer = Vec::new();

    let extension = get_extension_from_filename(filename).unwrap();

    let file = File::open(filename).unwrap();
    let mut reader = BufReader::new(file);
    let records: Vec<_>;
    println!("{}", extension);
    if extension == "zst"{
        let decoder = Decoder::new(reader).unwrap();
        let decompressed_reader = BufReader::new(decoder);
        let mut fa_reader = Reader::new(decompressed_reader);
        records = fa_reader.records().collect::<Result<Vec<_>, _>>().expect("Failed to access to record.");
    }else if extension == "fa" || extension == "fna" || extension == "fasta" {
        let mut fa_reader = Reader::new(reader);
        records = fa_reader.records().collect::<Result<Vec<_>, _>>().expect("Failed to access to record.");
    }else{
        println!("Error, file extension not handled.");
        std::process::exit(1);
    }
    let kmer_counter_mutex: Arc<Mutex<u64>> = Arc::new(Mutex::new(0));
    println!("FILE NUMBER: {}", file_number);
    // Process each record in parallel
    records.par_iter().for_each(|record| {
        let mut counter : u64 = 1;
        let mut curr_kmer= RawKmer::<K, KT>::from_nucs(&record.seq[..K]);
        let mut canon = curr_kmer.canonical().to_int();
        handle_kmer(omni_kmer_map, multi_kmer_map, file_number, canon);
        for (_i, nuc) in record.seq[K..].iter().filter_map(KT::from_nuc).enumerate() {
            curr_kmer = curr_kmer.append(nuc);
            canon = curr_kmer.canonical().to_int();
            counter +=1;
            handle_kmer(omni_kmer_map, multi_kmer_map, file_number, canon);
        }

        *kmer_counter_mutex.lock().unwrap() += counter;
    });
    if file_number != 0{
        check_omni(omni_kmer_map, multi_kmer_map, file_number);
    }
    println!("I have seen {} total k-mers", kmer_counter_mutex.lock().unwrap());
    kmer_counter_mutex
}

fn handle_kmer(omni_kmer_map: &Vec<Arc<Mutex<HashMap<KT, Cell<bool>>>>>, multi_kmer_map: &Vec<Arc<Mutex<HashMap<KT, COLORPAIR>>>>, file_number: usize, canon: KT){
    let pos = canon.to_usize().unwrap()%SHARD_AMOUNT;
    //counter += 1;
    let mut omni_kmer_map_lock = omni_kmer_map.get(pos).unwrap().lock().unwrap();
    if file_number == 0{
        omni_kmer_map_lock.insert(canon, Cell::new(false));
    }else if let Some(elem) = omni_kmer_map_lock.get_mut(&canon) {
        elem.set(true);
    }else{
        let mut multi_kmer_map_lock = multi_kmer_map.get(pos).unwrap().lock().unwrap();
        multi_kmer_map_lock.entry(canon).and_modify(|pair|{
            pair.0.set(file_number, true);
        }).or_insert_with(|| {
            let mut bv: BitArr!(for NB_FILES, in u8) = BitArray::<_>::ZERO;
            bv.set(file_number, true);
            (bv, Cell::new(false))
        });
        drop(multi_kmer_map_lock);
    }
    drop(omni_kmer_map_lock);
}

fn check_omni(omni_kmer_map: &Vec<Arc<Mutex<HashMap<KT, Cell<bool>>>>>, multi_kmer_map: &mut Vec<Arc<Mutex<HashMap<KT, COLORPAIR>>>>, file_number: usize){
    let kmer_to_remove_mutex: Arc<Mutex<u64>> = Arc::new(Mutex::new(0));
    let kmer_counter_mutex: Arc<Mutex<u64>> = Arc::new(Mutex::new(0));
    (0..SHARD_AMOUNT).into_par_iter().for_each(|i|{
        let mut counter: u64 = 0;
        let mut to_remove: u64 = 0;
        let mut omni_kmer_map_lock = omni_kmer_map.get(i).unwrap().lock().unwrap();
        let mut omni_iterator = omni_kmer_map_lock.iter();
        while let Some(entry) = omni_iterator.next(){
            let curr_k = num2str(*entry.0);
            counter+=1;
            if !(entry.1.get()){
                to_remove += 1;
                let mut multi_kmer_map_lock = multi_kmer_map.get(i).unwrap().lock().unwrap();
                multi_kmer_map_lock.entry(*entry.0).and_modify(|pair|{
                    pair.0.set(file_number, false);
                }).or_insert_with(|| {
                    let mut bv: BitArr!(for NB_FILES, in u8) = BitArray::<_>::ZERO;
                    let mut i = 0;
                    while i < file_number{
                        bv.set(i, true);
                        i += 1;
                    }
                    (bv, Cell::new(false))
                });
                drop(multi_kmer_map_lock);
            }
        }
        *kmer_counter_mutex.lock().unwrap() += counter;
        *kmer_to_remove_mutex.lock().unwrap() += to_remove;
        omni_kmer_map_lock.retain(|_k,v| v.get());
        for elem in omni_kmer_map_lock.iter(){
            elem.1.set(false);
        }
    });
    println!("I have seen {} k-mers total in the omnicolor", kmer_counter_mutex.lock().unwrap());
    println!("I have seen {} k-mers to be removed from omnicolor", kmer_to_remove_mutex.lock().unwrap());
    
}

/*
CONSTRUCTS SIMPLITIGS FROM K-MERS GATHERED IN "READ_FASTA".
SIMPLITIGS ARE MONOCHROMATIC.
OMNICOLORED (SEEN IN EVERY INPUT FILE) SIMPLITIGS ARE WRITTEN IN THE SAME FILE (omnicolor.kloe)
SIMPLITIGS ARE CONSTRUCTED USING PROPHASM'S GREEDY ALGORITHM.
*/
fn compute_omnicolored_simplitigs(omni_kmer_map:  &mut Vec<HashMap<KT, Cell<bool>>>, output_dir: &String){
    let mut nb_kmer = 0;
    for i in 0..SHARD_AMOUNT{
        nb_kmer += omni_kmer_map.get(i).unwrap().len();
        let kmer_map = omni_kmer_map.get(i).unwrap();
        let mut kmer_iterator = kmer_map.iter();
        while let Some(entry) = kmer_iterator.next(){
            if entry.1.get(){
                entry.1.set(false);
            }
        }
    }
    println!("I will reconstruct omnicolored simplitigs from {} kmers", nb_kmer);

    let omni_simpli = Instant::now();
    //let mut omni_f = Encoder::new(File::create(output_dir.clone()+"omnicolor.kloe").expect("Unable to create file"), 0).unwrap();
    let mut omni_f = BufWriter::new(File::create(output_dir.clone()+"omnicolor.kloe").expect("unable to create file"));
    let mut cpt = 0;
    for i in 0..SHARD_AMOUNT{
        let kmer_map = omni_kmer_map.get(i).unwrap();
        let mut kmer_iterator = kmer_map.iter();
        while let Some(entry) = kmer_iterator.next(){
            if !entry.1.get(){
                entry.1.set(true);
                let mut forward = true;
                let mut backward = true;
                let mut simplitig = num2str(entry.0.clone());
                while forward{
                    let curr_kmer = RawKmer::<K, KT>::from_nucs(&simplitig[(simplitig.len()-K)..].as_bytes());
                    forward = extend_omnicolor_forward(&curr_kmer, omni_kmer_map, &mut simplitig);
                }
                while backward{
                    let curr_kmer = RawKmer::<K, KT>::from_nucs(&simplitig[0..K].as_bytes());
                    backward = extend_omnicolor_backward(&curr_kmer, omni_kmer_map, &mut simplitig);
                }
                let simplitig_size:u32 = simplitig.len().to_u32().unwrap();
                //println!("CURR KMER = {}", num2str(*entry.0));
                /*println!("SIMPLITIG SIZE = {}", simplitig_size);
                if usize::BITS == 64{
                    println!("OUI");
                }else {
                    println!("NON");
                }
                //println!("SIMPLITIG = {}", &simplitig);
                let mut input = String::new();
                io::stdin().read_line(&mut input).expect("error: unable to read user input"); */
                cpt += simplitig.len() - K + 1;
                omni_f.write_all(&simplitig_size.to_le_bytes()).expect("Unable to write simplitig size");  
                omni_f.write_all(&str2num(&simplitig)).expect("Unable to write data");
            }
        }
    }    
    omni_f.flush();         
    println!("NB KMERS: {}", cpt);
    let fin_omni_simpli = omni_simpli.elapsed();
    println!("Construction of omnicolored simplitigs took: {:.2?} seconds.", fin_omni_simpli);
}

/*
FORWARD EXTENSION FOR SIMPLITIG CREATION.
CHECKS IF SUCCESSORS ARE THE SAME COLOR AS CURRENT K-MER.
*/
fn extend_omnicolor_forward(curr_kmer: &RawKmer<K, KT>, kmer_map:  &Vec<HashMap<KT, Cell<bool>>>, simplitig: &mut String) -> bool{
    for succs in curr_kmer.successors(){
        let canon = succs.canonical();
        //let mini = find_min(succs.canonical());
        //let pos = mini.to_usize().unwrap()%SHARD_AMOUNT;
        //let canon_succs = succs.clone().canonical();
        let pos = canon.to_int().to_usize().unwrap()%SHARD_AMOUNT;
        let kmer_map = kmer_map.get(pos).unwrap();
        if kmer_map.contains_key(&canon.to_int()){
            if !kmer_map.get(&canon.to_int()).unwrap().get() {
                simplitig.push(*succs.to_nucs().last().unwrap() as char);
                kmer_map.get(&canon.to_int()).unwrap().set(true);
                return true;
            }
        }
    }
    false
}

/*
BACKWARD EXTENSION FOR SIMPLITIG CREATION.
CHECKS IF PREDECESSORS ARE THE SAME COLOR AS CURRENT K-MER.
INSERTS FIRST NUCLEOTIDE OF PREDECESSOR (CHECKED MULTIPLE TIMES)
*/
fn extend_omnicolor_backward(curr_kmer: &RawKmer<K, KT>, kmer_map:  &Vec<HashMap<KT, Cell<bool>>>, simplitig: &mut String) -> bool{
    for preds in curr_kmer.predecessors(){
        let canon = preds.canonical();
        //let mini = find_min(preds.canonical());
        //let pos = mini.to_usize().unwrap()%SHARD_AMOUNT;
        let pos = canon.to_int().to_usize().unwrap()%SHARD_AMOUNT;
        let kmer_map = kmer_map.get(pos).unwrap();
        if kmer_map.contains_key(&canon.to_int()){
            if !kmer_map.get(&canon.to_int()).unwrap().get() {
                simplitig.insert(0,*preds.to_nucs().first().unwrap() as char);
                kmer_map.get(&canon.to_int()).unwrap().set(true);
                return true;
            }
        }
    }
    false
}


/*
CONSTRUCTS SIMPLITIGS FROM K-MERS GATHERED IN "READ_FASTA".
SIMPLITIGS ARE MONOCHROMATIC.
OMNICOLORED (SEEN IN EVERY INPUT FILE) SIMPLITIGS ARE WRITTEN IN THE SAME FILE (omnicolor.fa.zstd)
MULTICOLORED SIMPLITIGS ARE WRITTEN IN ANOTHER FILE (multicolor.fa.zstd).
SIMPLITIGS ARE CONSTRUCTED USING PROPHASM'S GREEDY ALGORITHM.
*/
fn compute_multicolored_simplitigs(multi_kmer_map:  &mut Vec<HashMap<KT, COLORPAIR>>, output_dir: &String) -> HashMap<bitvec::prelude::BitArray<[u8; ARRAY_SIZE]>, (usize, Vec<usize>)>{//Arc<Mutex<HashMap<bitvec::prelude::BitArray<[u8; ARRAY_SIZE]>, (usize, Vec<usize>)>>>{
    let mut nb_kmer = 0;
    for i in 0..SHARD_AMOUNT{
        nb_kmer += multi_kmer_map.get(i).unwrap().len();
    }
    println!("I will reconstruct multicolored simplitigs from {} kmers", nb_kmer);

    let multi_simpli = Instant::now();
    let mut multi_f = BufWriter::new(File::create(output_dir.clone()+"temp_multicolor.fa").expect("Unable to create file"));
    //let mut color_simplitig_size: Arc<Mutex<HashMap<bitvec::prelude::BitArray<[u8; ARRAY_SIZE]>, (usize, Vec<usize>)>>> = Arc::new(Mutex::new(HashMap::new()));
    let mut color_simplitig_size: HashMap<bitvec::prelude::BitArray<[u8; ARRAY_SIZE]>, (usize, Vec<usize>)> = HashMap::new();
    for i in 0..SHARD_AMOUNT{
        let kmer_map = multi_kmer_map.get(i).unwrap();
        let mut kmer_iterator = kmer_map.iter();
        while let Some((key, curr_cell)) = kmer_iterator.next(){
            if !curr_cell.1.get(){
                curr_cell.1.set(true);
                let mut forward = true;
                let mut backward = true;
                let mut simplitig = num2str(*key);
                while forward{
                    let curr_kmer = RawKmer::<K, KT>::from_nucs(&simplitig[(simplitig.len()-K)..].as_bytes());
                    forward = extend_forward(&curr_kmer, &multi_kmer_map, &mut simplitig, &curr_cell.0, i);
                }
                while backward{
                    let curr_kmer = RawKmer::<K, KT>::from_nucs(&simplitig[0..K].as_bytes());
                    backward = extend_backward(&curr_kmer, &multi_kmer_map, &mut simplitig, &curr_cell.0, i);
                }
                let simplitig_size = simplitig.len();
                color_simplitig_size.entry(curr_cell.0)
                            .and_modify(|total_size| {
                                total_size.0 += simplitig_size.div_ceil(4);
                                total_size.1.push(simplitig_size);
                            })
                            .or_insert((simplitig_size.div_ceil(4), vec![simplitig_size]));
                write_multicolored_simplitig(&str2num(&simplitig), &mut multi_f,&curr_cell.0, &simplitig.len());
            }
        }
    }
    multi_f.flush();
    let fin_multi_simpli = multi_simpli.elapsed();
    println!("Construction of multicolored simplitigs took: {:.2?} seconds.", fin_multi_simpli);
    color_simplitig_size
}
/*
FORWARD EXTENSION FOR SIMPLITIG CREATION.
CHECKS IF SUCCESSORS ARE THE SAME COLOR AS CURRENT K-MER.
*/
fn extend_forward(curr_kmer: &RawKmer<K, KT>, vec_kmer_map:  &Vec<HashMap<KT, COLORPAIR>>, simplitig: &mut String, color: &BitArray<[u8;ARRAY_SIZE]>, pos_curr_kmer: usize) -> bool{
    for succs in curr_kmer.successors(){
        let canon = succs.canonical().to_int();
        if let Some(succ_pair) = vec_kmer_map.get(pos_curr_kmer).unwrap().get(&canon){
            if succ_pair.0.eq(color) & !succ_pair.1.get() {
                simplitig.push(*succs.to_nucs().last().unwrap() as char);
                succ_pair.1.set(true);
                return true;
            }
        }else{
            //let mini = find_min(succs.canonical());
            //let pos = mini.to_usize().unwrap()%SHARD_AMOUNT;
            let pos = canon.to_usize().unwrap()%SHARD_AMOUNT;
            let kmer_map = vec_kmer_map.get(pos).unwrap();
            if let Some(succ_pair) = kmer_map.get(&canon){
                if succ_pair.0.eq(color) & !succ_pair.1.get() {
                    simplitig.push(*succs.to_nucs().last().unwrap() as char);
                    succ_pair.1.set(true);
                    return true;
                }
            }
        }
    }
    false
}

/*
BACKWARD EXTENSION FOR SIMPLITIG CREATION.
CHECKS IF PREDECESSORS ARE THE SAME COLOR AS CURRENT K-MER.
INSERTS FIRST NUCLEOTIDE OF PREDECESSOR (CHECKED MULTIPLE TIMES)
*/
fn extend_backward(curr_kmer: &RawKmer<K, KT>, vec_kmer_map:  &Vec<HashMap<KT, COLORPAIR>>, simplitig: &mut String, color: &BitArray<[u8;ARRAY_SIZE]>, pos_curr_kmer: usize) -> bool{
    for preds in curr_kmer.predecessors(){
        let canon = preds.canonical().to_int();
        if let Some(pred_pair) = vec_kmer_map.get(pos_curr_kmer).unwrap().get(&canon){
            if pred_pair.0.eq(color) & !pred_pair.1.get() {
                simplitig.insert(0, *preds.to_nucs().first().unwrap() as char);
                pred_pair.1.set(true);
                return true;
            }
        }else{
            //let mini = find_min(preds.canonical());
            //let pos = mini.to_usize().unwrap()%SHARD_AMOUNT;
            let pos = canon.to_usize().unwrap()%SHARD_AMOUNT;
            let kmer_map = vec_kmer_map.get(pos).unwrap();
            if let Some(pred_pair) = kmer_map.get(&canon){
                if pred_pair.0.eq(color) & !pred_pair.1.get() {
                    simplitig.insert(0,*preds.to_nucs().first().unwrap() as char);
                    pred_pair.1.set(true);
                    return true;
                }
            }
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
    let mut out_mult_file = BufWriter::new(File::options().write(true).read(true).create(true).open(path.clone()).expect("Unable to create file"));
    //let mut color_set = HashSet::new();
    let mut color_cursor_pos = compute_cursor_start_pos(color_simplitig);
    let temp_multicolor_file = File::open(output_dir.clone()+"temp_multicolor.fa").unwrap();
    let mut temp_multicolor_reader = BufReader::new(&temp_multicolor_file);
    const COLOR_SIZE: usize = NB_FILES.div_ceil(8);
    let mut cursor: usize = 0;
    let metadata = temp_multicolor_file.metadata().unwrap();
    let file_size: usize = metadata.len() as usize;
    println!("FILE SIZE = {}", file_size);
    while cursor < file_size{
        cursor+= COLOR_SIZE;
        let mut id = [0; ARRAY_SIZE];
        temp_multicolor_reader.read_exact(&mut id).expect("Error reading color in temp file");
        let mut size_buf = [0; 4];
        cursor += 4;
        temp_multicolor_reader.read_exact(&mut size_buf).expect("Error reading simplitig size in temp file");
        let size_to_read: u32 = u32::from_le_bytes(size_buf).div_ceil(4);
        //let size_simplitig: u32 = u32::from_le_bytes(size_buf);
        //println!("READING CURSOR = {}", cursor);
        cursor += size_to_read as usize;
        //println!("COLOR: {}\nSIZE: {}", &id.to_vec()[0], size_simplitig);
        let mut simplitig = vec![0; size_to_read as usize];
        //println!("Reading {} Bytes", simplitig.len());
        temp_multicolor_reader.read_exact(&mut simplitig).expect("Error reading simplitig");
        //println!("SIMPLITIG: {}", vec2str(&simplitig.to_vec(), &(size_simplitig as usize)));
        let _ = write_sorted(&mut out_mult_file, simplitig, &mut color_cursor_pos, &id.to_vec());
        //println!("READING CURSOR = {}", cursor);
        //let mut input = String::new();
        //io::stdin().read_line(&mut input).expect("error: unable to read user input");
    }
    out_mult_file.flush();
    let _ = remove_file(output_dir.clone()+"temp_multicolor.fa");
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

fn write_sorted(out_mult_file: &mut BufWriter<File>, content: Vec<u8>, color_cursor_pos: &mut IndexMap<bitvec::prelude::BitArray<[u8; ARRAY_SIZE]>, u64>, color: &Vec<u8>) -> io::Result<()> {
    //let mut encoder = Compressor::new(12).unwrap();
    //let compressed = encoder.compress(content.as_bytes()).unwrap();
    let mut color_vec: BitArr!(for NB_FILES, in u8) = BitArray::<_>::ZERO;
    let mut str_color= String::new();
    for i in 0..NB_FILES{
        let part = i/8;
        let bit = (color.get(part).unwrap() >> i) & 1;
        if bit == 1{
            color_vec.set(i,true);
            str_color.push('1');
        }else {
            str_color.push('0');
        }
    }
    //println!("COLOR = {}", str_color);
    
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
    //let mut file = Encoder::new(File::create(output_dir.clone()+"multicolor_bucket_size.txt.zst").expect("Unable to create interface file"));
    let mut nb_seen = 0;
    let mut color_file = Encoder::new(File::create(output_dir.clone()+"multicolor_bucket_size.txt.zst").expect("Unable to create file"), 12).unwrap();
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
            let buf : String = String::from(color + "," + &(cursor_end.to_string()) + ":" + &bucket_sizes + "\n");
            let _ = color_file.write_all(buf.as_bytes());
            
            //writeln!(file, "{},{}:{}", color, cursor_end, bucket_sizes).unwrap();
        }
        nb_seen = 0;
    }
    color_file.flush();
    let _ = color_file.finish();
}





fn write_multicolored_simplitig(simplitig: &Vec<u8>, multi_f: &mut BufWriter<File>, color: &BitArray<[u8;ARRAY_SIZE]>, size_str : &usize){
    let size:u32 = *size_str as u32;
    let mut id = Vec::new();
    let mut cpt = 0;
    let mut curr_int: u8 = 0;
    //let mut str_color = String::new();
    for e in color.iter(){
        if *e{
            curr_int +=1 << cpt;
            //str_color.push('1');
        }/*else {
            str_color.push('0');
        }*/
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