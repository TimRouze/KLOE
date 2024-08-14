use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use seq_io::fasta::{Reader, Record};

/*
COMPUTES STATS FROM COMPRESSED OMNICOLOR AND MULTICOLOR FILES.

1) NB KMER PER COLOR
2) NB SIMPLITIG PER EXISTING SIMPLITIG SIZE

*/
pub fn compute_stats(omnicolor: &str, multicolor: &str, out_dir: &str, k: &usize){
    let mut name = String::from(out_dir)+"kmer_per_color.txt";
    let mut kmer_per_color_file = File::create(name).expect("Unable to create file");
    name = String::from(out_dir)+"simplitig_per_simplitig size.txt";
    let mut simplitig_per_simplitig_size_file = File::create(name).expect("Unable to create file");
    let ( reader, _compression) = niffler::get_reader(Box::new(File::open(multicolor).unwrap())).unwrap();
    let mut fa_reader = Reader::new(reader);
    let mut color_nb_kmer_map: HashMap<String, usize> = HashMap::new();
    let mut size_nb_simplitig_map: HashMap<usize, usize> = HashMap::new();
    while let Some(data) = fa_reader.next(){
        let line = data.unwrap();
        let id = String::from(line.id().unwrap());
        color_nb_kmer_map.entry(id)
                        .and_modify(|nb_kmer| *nb_kmer += line.seq().len()-k+1)
                        .or_insert(line.seq().len()-k+1);
        size_nb_simplitig_map.entry(line.seq().len())
                        .and_modify(|nb_kmer| *nb_kmer += 1)
                        .or_insert(1);
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
        size_nb_simplitig_map.entry(line.seq().len())
                        .and_modify(|nb_kmer| *nb_kmer += 1)
                        .or_insert(1);
    }
    let mut counter = 0;
    let mut counter_2 = 0;
    for (_key, value) in color_nb_kmer_map {
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
        writeln!(kmer_per_color_file, "{}", value).expect("Unable to write in stat file");
    }
    for(key, value) in size_nb_simplitig_map{
        counter_2 += 1;
        writeln!(simplitig_per_simplitig_size_file, "{}: {}", key, value).expect("Unable to write in stat file");
    }
    println!("Wrote {} lines in kmer per color file.", counter);
    println!("Wrote {} lines in kmer per simplitig size file.", counter_2)
}