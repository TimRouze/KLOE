use seq_io::fasta::{Reader, Record};
use std::path::{Path, PathBuf};
use std::fs::File;
use std::io::{BufReader, Write, BufRead};
use zstd::stream::write::Encoder;
use::rayon::prelude::*;

/// DECOMPRESS
/// Omnicolor: file containing omnicolored kmers
/// Multicolor: file containing multicolored kmers
/// Input names: fof containing genomes to extract. SAME ORDER AS FOF USED TO CREATE COMPRESSED FILES
/// 
/// DECOMPRESS will dump simplitigs from the requested files. 
/// Creating one output file per genome requested.
/// Omnicolored kmers are dumped in every out file.
/// Multicolored kmers are split according to their colors
pub fn decompress(omnicolor: &str, multicolor: &str, input_names: &str, out_dir: PathBuf){
    let input_fof = File::open(input_names).unwrap();
    let reader = BufReader::new(input_fof);
    let filenames: Vec<_> = reader.lines().collect::<Result<_, _>>().unwrap();
    (0..filenames.len()).into_par_iter().for_each(|file_number|{
        let filename = filenames.get(file_number).unwrap();
        println!("{}", filename);
        let path = Path::new(&filename).file_stem().unwrap().to_str().unwrap();
        let without_path = Path::new(path).file_name().unwrap();
        let dump_curr_file = out_dir.join(format!("Dump_{}.zstd", without_path.to_str().unwrap()));
        let mut dump_file = Encoder::new(File::create(dump_curr_file).expect("Unable to create file"), 0).unwrap();
        let ( reader, _compression) = niffler::get_reader(Box::new(File::open(multicolor).unwrap())).unwrap();
        let mut fa_reader = Reader::new(reader);
        while let Some(data) = fa_reader.next(){
            let line = data.unwrap();
            let header: &str = line.id().unwrap();
            if header.chars().nth(file_number).unwrap() == '1'{
                let mut header_fa = String::from(">");
                header_fa.push_str(header);
                header_fa.push('\n');
                let to_write = std::str::from_utf8(line.seq()).unwrap();
                dump_file.write_all((header_fa).as_bytes()).unwrap();
                dump_file.write_all(to_write.as_bytes()).unwrap(); 
                dump_file.write_all(b"\n").unwrap();
            }
        }
        let ( reader, _compression) = niffler::get_reader(Box::new(File::open(omnicolor).unwrap())).unwrap();
        let mut fa_reader = Reader::new(reader);
        //let mut counter = 0;
        while let Some(data) = fa_reader.next(){
            let line = data.unwrap();
            let to_write = std::str::from_utf8(line.seq()).unwrap();
            let mut header_fa = String::from(">");
            //counter += to_write.len()-31+1;
            header_fa.push_str("OMNI\n");
            dump_file.write_all((header_fa).as_bytes()).unwrap();
            dump_file.write_all(to_write.as_bytes()).unwrap();
            dump_file.write_all(b"\n").unwrap();
        }
        //println!("During decompression, I have seen {} k-mers", counter);
        dump_file.finish().expect("Error writing decompressed data");
    }); 
}
