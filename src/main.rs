#![allow(dead_code)]

mod utils;
mod decompress;
mod parser;
mod compress;
use clap::Parser;
use std::env;
use std::path::{Path, PathBuf};

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Result, Seek, Write};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    ///output: Option<String>,
    /// Decompression "compress" to compress input, "decompress" to decompress input, "stats" to get kmer stats
    decompress: Option<String>,
    ///Input file list (Compression)
    #[arg(short, long, default_value_t=String::from(""))]
    input_list: String,
    /// Number of threads (defaults to 1)
    #[arg(short, long, default_value_t = num_cpus::get())]
    threads: usize,
    ///Output directory
    #[arg(short, long, default_value_t = String::from(""))]
    out_dir: String,
    ///input directory for decompression
    #[arg(short, long, default_value_t = String::from(""))]
    compressed_dir: String,
    ///List of files to decompress
    #[arg(short = 'Q', long, default_value_t = String::from(""))]
    wanted_files: String,
    ///Temporary directory for unitigs parsing
    #[arg(short = 'd', long, default_value_t = String::from(""))]
    temp_dir: String,
    ///Max memory (RAM) for fulgor, default = 8GB
    #[arg(short = 'r', long, default_value_t = 8)]
    memory: usize,
    ///K value, default = 31
    #[arg(short, long, default_value_t = 31)]
    k_size: usize,
    ///Minimizer size (< k), default = 7
    #[arg(short, long, default_value_t = 7)]
    minimizer_size: usize,
    /// Optionally verify that all canonical k-mers are preserved with the correct dataset IDs
    #[arg(long = "verify-kmers", default_value_t = false)]
    verify_kmers: bool,
    /// Number of partitions is 2^partition_power (default 1024 partitions)
    #[arg(short = 'P', long = "partition-power", default_value_t = 10)]
    partition_power: u32,
    /// Skip sorting within partitions and during final merge (output will not be globally sorted)
    #[arg(long = "skip-sort", default_value_t = false)]
    skip_sort: bool,
    /// Unitigs
    #[arg(long = "unitigs", default_value_t = false)]
    unitigs: bool,
    /// Matchtigs
    #[arg(long = "matchtigs", default_value_t = false)]
    matchtigs: bool,
    /// Eulertigs
    #[arg(long = "eulertigs", default_value_t = false)]
    eulertigs: bool,
    
}   
const BLOCK_SIZE: usize = 1 << (12 - 3);
const SHARD_AMOUNT: usize = 1024;
const M: u8 = 7;

fn main() {
    let args = Args::parse();
    
    let output_dir = args.out_dir;
    let input_dir = args.compressed_dir;
    //env::set_var("RAYON_NUM_THREADS", args.threads.to_string());
    let input_fof = args.input_list;
    let threads = args.threads;
    let temp_dir = args.temp_dir;
    let memory = args.memory;
    let k = args.k_size;
    let m = args.minimizer_size;
    //TODO HANDLE ERRORS FOR COMP AND DECOMP
    let wanted_path = args.wanted_files;
    if let Some(do_decompress) = args.decompress{
        if do_decompress == "decompress"{
            println!("Checking archive integrity...");
            is_compressed_dir_complete(input_dir.clone());
            let _ = decompress::decompress(&String::from("bucket_sizes.txt"), &String::from("id_to_color_id.txt.zst"), &String::from("tigs_kloe.fa"), &String::from("positions_kloe.txt.zst"), &String::from("filenames_id.txt"), &output_dir, &wanted_path, input_dir);
            //let _ = graph_build::init_decompress(String::from("bucket_sizes.txt.zst"), String::from("id_to_color_id.txt.zst"), unitigs_file, &output_dir, &wanted_path, &input_dir);
        }else if do_decompress == "compress"{
            /*let _ = compress::compress(
                &output_dir, 
                &input_fof, 
                threads, 
                &temp_dir, 
                k, 
                m, 
                partition_power, 
                compaction_threads
            );*/
            //parser::run_parser(k, m, 10_u32, PathBuf::from(output_dir), PathBuf::from(input_fof), threads, compaction_threads, false);
            let _ = compress::compress(&output_dir, &input_fof, threads, k, m, args.partition_power, args.unitigs, args.matchtigs, args.eulertigs);
            //let _ = graph_build::build_graphs(&output_dir, &input_fof, &threads, &temp_dir, &memory);
        }
    }else {
        //parser::run_parser(PathBuf::from(input_fof), PathBuf::from(output_dir), k, m, 10_u32, threads, compaction_threads, false, false);

        let mut input_fof_reader = BufReader::new(File::open(input_fof).expect("unable to open fof"));
        let mut filename = String::new();
        let mut filenames = Vec::new();
        while input_fof_reader.read_line(&mut filename).unwrap() != 0{
            filename.pop();
            filenames.push(filename.clone());
            filename.clear();

        }
        let mut sequence_type;
        if args.unitigs{
            sequence_type = String::from("unitigs");
        }else if args.matchtigs{
            sequence_type = String::from("matchtigs");
        }else if args.eulertigs{
            sequence_type = String::from("eulertigs");
        }else{
            sequence_type = String::from("simplitigs");
        }
        let id_cid_line_sizes = compress::sort_by_bucket(&output_dir, 256, sequence_type);
        let mut fof_id = BufWriter::new(File::create(output_dir.clone() + "filenames_id.txt").expect("Failed to create fof file"));
        let mut file_cpt: usize = 0;
        for filename in filenames{
            println!("a{}a", filename);
            fof_id.write_all((filename + ":" + id_cid_line_sizes.get(file_cpt).unwrap().to_string().as_str() + "\n").as_bytes()).unwrap();
            file_cpt += 1;
        }

        println!("Wrong positional arguments given. Values are 'compress' or 'decompress'");
        println!("Ex: if compression: I=my/fof.txt cargo r -r -- compress -f my_file_of_file.txt -o out_dir/ -t 12");
        println!("Ex: if decompression: I=my/fof.txt cargo r -r -- decompress -f my_file_of_file.txt --omnicolor-file out_dir/omnicolor.fa.zstd --multicolor-file out_dir/multicolor.fa.zstd -t 12");
    }
}

fn is_compressed_dir_complete(input_dir: String){
    if !Path::new(&format!("{input_dir}/filenames_id.txt")).exists(){
        panic!("file not found: {input_dir}/filenames_id.txt");
    }else if !Path::new(&format!("{input_dir}/positions_kloe.txt.zst")).exists(){
        panic!("Positions file not found");
    }else if !Path::new(&format!("{input_dir}/bucket_sizes.txt")).exists(){
        panic!("Tigs sizes file not found");
    }else if !Path::new(&format!("{input_dir}/id_to_color_id.txt.zst")).exists(){
        panic!("id to color id file not found");
    }else if !Path::new(&format!("{input_dir}/tigs_kloe.fa")).exists(){
        panic!("Tigs file not found");
    }else{
        println!("Archive complete, starting decompression...");
    }
}