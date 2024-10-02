use bitvec::prelude::BitArray;
use std::path::Path;
use std::io;
use std::io::BufRead;
use std::fs::File;
use std::cell::Cell;

pub mod constants {
    include!("constants.rs");
}
use constants::{ARRAY_SIZE, K, KT};
pub type COLORPAIR = (bitvec::prelude::BitArray<[u8; ARRAY_SIZE]>, Cell<bool>);


pub fn extract_filename(path: &str) -> Option<&str> {
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

pub fn char_array_to_bitarray(seq: &[u8]) -> bitvec::prelude::BitArray<[u8; ARRAY_SIZE]>{
    let mut bit_array = BitArray::<[u8; ARRAY_SIZE]>::ZERO;
    let mut cpt = 0;
    let seq_str = std::str::from_utf8(seq).unwrap();
    seq_str.chars().enumerate().for_each(|c|{
        //println!("{}", c.1);
        if c.1 != '0'{
            bit_array.set(cpt, true);
        }else{
            bit_array.set(cpt, false);
        }
        cpt += 1;
    });
    bit_array
}

pub fn vec2str(seq: &Vec<u8>, size: &usize) -> String{
    let mut res = String::from("");
    let mask = 3;
    for elem in seq.iter(){
        res += nuc2str(&(elem&mask));
        res += nuc2str(&((elem >> 2)&mask));
        res += nuc2str(&((elem >> 4)&mask));
        res += nuc2str(&((elem >> 6)&mask));
    }
    let _ = res.drain(size..);
    res
}

pub fn num2str(mut k_mer: KT) -> String{
    let mut res = String::from("");
    let mut nuc: KT;
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

pub fn nuc2str(nuc: &u8) -> &str{

    if nuc%4 == 0{
        "A"
    }else if nuc%4 == 1{
        "C"
    }else if nuc%4 == 2{//bebou
        "G"
    }else{
        "T"
    }
}

pub fn str2num(sequence: &String) -> Vec<u8>{
    let mut res = Vec::new();
    let mut tmp_res: u8 = 0;
    let mut i = 0;
    let mut shift = 0;
    let mut char_list = sequence.chars();
    //println!("{}", sequence.len());
    while let Some(nuc) = char_list.next(){
        tmp_res += nuc2int(&(nuc as u8)).unwrap() << shift;
        shift += 2;
        i += 1;
        //tmp_res += nuc2int(&(curr_byte.chars().nth(1).unwrap() as u8)).unwrap() << 2;
        //tmp_res += nuc2int(&(curr_byte.chars().nth(2).unwrap() as u8)).unwrap() << 4;
        //tmp_res += nuc2int(&(curr_byte.chars().nth(3).unwrap() as u8)).unwrap() << 6;
        if i%4 == 0{
            res.push(tmp_res);
            tmp_res = 0;
            shift = 0;
        }
    }
    if shift != 0{
        res.push(tmp_res);
    }
    //println!("SIZE AS BYTES: {}", res.len());
    res
}

pub fn nuc2int(b: &u8) -> Option<u8> {
    match b {
        b'A' | b'C' | b'T' | b'G' => Some((b / 3-1) % 4),
        _ => None,
    }
}

pub fn rev_comp_str(seq: &str) -> String{
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

pub fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}
