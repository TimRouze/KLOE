pub trait Convert<T> {
    fn str2num(input: T) -> Vec<u8>;
}
pub struct Converter;

impl Convert<&String> for Converter {
    fn str2num(sequence: &String) -> Vec<u8>{
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
            if i%4 == 0{
                res.push(tmp_res);
                tmp_res = 0;
                shift = 0;
            }
        }
        if shift != 0{
            res.push(tmp_res);
        }
        res
    }
}

impl Convert<&[u8]> for Converter {
    fn str2num(sequence: &[u8]) -> Vec<u8>{
        let mut res = Vec::new();
        let mut tmp_res: u8 = 0;
        let mut i = 0;
        let mut shift = 0;
        //println!("{}", sequence.len());
        for nuc in sequence.iter(){
            tmp_res += nuc2int(nuc).unwrap() << shift;
            shift += 2;
            i += 1;
            if i%4 == 0{
                res.push(tmp_res);
                tmp_res = 0;
                shift = 0;
            }
        }
        if shift != 0{
            res.push(tmp_res);
        }
        res
    }
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


pub fn nuc2int(b: &u8) -> Option<u8> {
    match b {
        b'A' | b'C' | b'T' | b'G' => Some((b / 3-1) % 4),
        _ => None,
    }
}

