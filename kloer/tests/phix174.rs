use std::fs;

use kloer::{compress, decompress, CompressOptions};

fn fasta_sequences(text: &str) -> Vec<Vec<u8>> {
    let mut sequences = Vec::new();
    let mut sequence = Vec::new();
    for line in text.lines() {
        if line.starts_with('>') {
            if !sequence.is_empty() {
                sequences.push(std::mem::take(&mut sequence));
            }
        } else {
            sequence.extend_from_slice(line.trim().as_bytes());
        }
    }
    if !sequence.is_empty() {
        sequences.push(sequence);
    }
    sequences
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => *other,
        })
        .collect()
}

#[test]
fn phix174_reference_relative_roundtrip() {
    let reference_path = std::path::Path::new("tests/data/phix174.fa");
    let reference_text = fs::read_to_string(reference_path).unwrap();
    let reference = fasta_sequences(&reference_text).remove(0);
    assert_eq!(reference.len(), 5_386);

    let exact = reference[300..1_300].to_vec();
    let mut divergent = reference[1_700..2_700].to_vec();
    divergent[37] = if divergent[37] == b'A' { b'C' } else { b'A' };
    divergent[503] = if divergent[503] == b'G' { b'T' } else { b'G' };
    divergent.splice(650..650, b"ACTG".iter().copied());
    divergent.drain(802..805);
    let reverse = reverse_complement(&reference[3_200..4_200]);
    let expected = vec![exact, divergent, reverse];

    let temp = tempfile::tempdir().unwrap();
    let unitigs_path = temp.path().join("unitigs.fa");
    let archive_path = temp.path().join("unitigs.kloer");
    let decoded_path = temp.path().join("decoded.fa");
    let mut input = String::new();
    for (index, sequence) in expected.iter().enumerate() {
        input.push_str(&format!(
            ">unitig-{index}\n{}\n",
            std::str::from_utf8(sequence).unwrap()
        ));
    }
    fs::write(&unitigs_path, &input).unwrap();

    let compression = compress(
        reference_path,
        &unitigs_path,
        &archive_path,
        &CompressOptions {
            threads: 2,
            ..CompressOptions::default()
        },
    )
    .unwrap();
    assert_eq!(compression.records, 3);
    assert_eq!(compression.mapped, 3);
    assert_eq!(compression.literal, 0);
    let archive_bytes = fs::metadata(&archive_path).unwrap().len();
    eprintln!(
        "phiX174 KLOER: input={} bytes archive={} bytes ratio={:.3}",
        input.len(),
        archive_bytes,
        archive_bytes as f64 / input.len() as f64
    );
    assert!(archive_bytes < input.len() as u64 / 4);

    let decompression = decompress(reference_path, &archive_path, &decoded_path).unwrap();
    assert_eq!(decompression.records, 3);
    let decoded = fasta_sequences(&fs::read_to_string(decoded_path).unwrap());
    assert_eq!(decoded, expected);
}
