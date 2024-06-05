use std::collections::HashMap;
use std::fs::File;
use std::io::{Write};
use seq_io::fasta::{Reader, Record};


pub fn compute_stats(omnicolor: &str, multicolor: &str, out_file: &str, k: &usize){
    let mut stat_file = File::create(out_file).expect("Unable to create file");
    let ( reader, _compression) = niffler::get_reader(Box::new(File::open(multicolor).unwrap())).unwrap();
    let mut fa_reader = Reader::new(reader);
    let mut color_nb_kmer_map: HashMap<String, usize> = HashMap::new();
    while let Some(data) = fa_reader.next(){
        let line = data.unwrap();
        let id = String::from(line.id().unwrap());
        color_nb_kmer_map.entry(id)
                        .and_modify(|nb_kmer| *nb_kmer += line.seq().len()-k+1)
                        .or_insert(line.seq().len()-k+1);
    }
    let ( reader, _compression) = niffler::get_reader(Box::new(File::open(omnicolor).unwrap())).unwrap();
    let mut fa_reader = Reader::new(reader);
    //let mut counter = 0;
    while let Some(data) = fa_reader.next(){
        let line = data.unwrap();
        let id = String::from(line.id().unwrap());
        color_nb_kmer_map.entry(id)
                        .and_modify(|nb_kmer| *nb_kmer += line.seq().len()-k+1)
                        .or_insert(line.seq().len()-k+1);
    }
    let mut counter = 0;
    for (key, value) in color_nb_kmer_map {
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
        writeln!(stat_file, "{}", value);
    }
    println!("Wrote {} lines in stat file.", counter);
}