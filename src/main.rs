#![allow(dead_code)]

mod utils;
mod kmer;
mod decompress;
mod stats;
mod graph_build;
mod parser;
use clap::Parser;
use std::env;
use std::path::{Path, PathBuf};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    ///output: Option<String>,
    /// Decompression "compress" to compress input, "decompress" to decompress input, "stats" to get kmer stats
    decompress: Option<String>,
    ///Input file list (Compression)
    #[arg(short, long, default_value_t=String::from(""))]
    input_list: String,
    /// Number of threads (defaults to all available threads)
    #[arg(short, long, default_value_t = 1)]
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
    ///Temporary directory for graph construction
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

}
pub mod constants {
    include!("constants.rs");
}
const BLOCK_SIZE: usize = 1 << (12 - 3);
const SHARD_AMOUNT: usize = 1024;
const M: u8 = 7;

fn main() {
    let args = Args::parse();
    
    let output_dir = args.out_dir;
    let input_dir = args.compressed_dir;
    env::set_var("RAYON_NUM_THREADS", args.threads.to_string());
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
            let _ = graph_build::decompress(&String::from("bucket_sizes.txt"), &String::from("id_to_color_id.txt.zst"), &String::from("tigs_kloe.fa"), &String::from("positions_kloe.txt.zst"), &String::from("filenames_id.txt"), &output_dir, &wanted_path, input_dir);
            //let _ = graph_build::init_decompress(String::from("bucket_sizes.txt.zst"), String::from("id_to_color_id.txt.zst"), unitigs_file, &output_dir, &wanted_path, &input_dir);
        }else if do_decompress == "compress"{
            let _ = graph_build::build_graphs(&output_dir, &input_fof, &threads, &temp_dir, &memory);
        }/*else if do_decompress == "stats"{
            let unitigs_file = args.unitigs_file;
            let k = K;
            if unitigs_file != ""{
                let _ = compute_stats( &unitigs_file, &output_dir, &k);
            }else{
                println!("Error, multicolor and/or omnicolor file(s) are mandatory");
            }
        }*/
    }else {
        parser::test_parser(k, m, 10_u32, PathBuf::from(output_dir), PathBuf::from(&input_fof), threads);
        
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