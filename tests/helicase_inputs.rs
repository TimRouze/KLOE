use ggcat_api::{read_fasta_with_helicase, SequencesReader};
use std::io::{self, BufReader, Read, Write};

fn parse(
    input: &[u8],
    capacity: usize,
    overlap: Option<usize>,
    chunk: usize,
) -> io::Result<(u64, Vec<(Vec<u8>, Vec<u8>)>)> {
    let mut records = Vec::new();
    let bases = read_fasta_with_helicase(
        BufReader::with_capacity(capacity, input),
        overlap,
        true,
        false,
        chunk,
        &mut |s| {
            records.push((s.ident_data.to_vec(), s.seq.to_vec()));
            Ok(())
        },
    )?;
    Ok((bases, records))
}

#[test]
fn helicase_batches_preserve_wrapping_headers_and_ambiguous_bases() {
    let mut input = Vec::new();
    let mut expected = Vec::new();
    for i in 0..5000 {
        let header = format!(">sample_{i} ka:f:12.5 extra=hello");
        writeln!(input, "{header}\r\nacGTNrY.-\r\nTgCaNNNNACGT").unwrap();
        expected.push((header.into_bytes(), b"ACGTNNNNNTGCANNNNACGT".to_vec()));
    }
    for capacity in [7, 4093, 65536] {
        let (bases, got) = parse(&input, capacity, Some(30), 1024 * 1024).unwrap();
        assert_eq!(got.len(), expected.len());
        for (record, expected_record) in got.iter().zip(&expected) {
            assert_eq!(record, expected_record);
        }
        assert_eq!(
            bases,
            expected.iter().map(|r| r.1.len() as u64).sum::<u64>()
        );
    }
}

#[test]
fn oversized_sequences_stay_chunked_with_exact_overlap() {
    let sequence = "ACGTNNtgca".repeat(100_003);
    let input = format!(">long km:f:10\n{sequence}\n>next km:f:20\nTGCATGCA");
    let (bases, records) = parse(input.as_bytes(), 97, Some(30), 1024).unwrap();
    assert_eq!(bases, sequence.len() as u64 + 8);
    let mut restored = Vec::new();
    for (i, (header, chunk)) in records[..records.len() - 1].iter().enumerate() {
        assert_eq!(header, b">long km:f:10");
        assert!(chunk.len() <= 1024);
        restored.extend_from_slice(&chunk[if i == 0 { 0 } else { 30 }..]);
    }
    assert_eq!(restored, sequence.to_ascii_uppercase().as_bytes());
    assert_eq!(records.last().unwrap().1, b"TGCATGCA");
    let (_, whole) = parse(input.as_bytes(), 4096, None, 1024).unwrap();
    assert_eq!(whole.len(), 2);
    assert_eq!(whole[0].1, restored);
}

#[test]
fn comments_empty_sequences_and_embedded_cr_remain_compatible() {
    let input = b";comment\n>empty\n>a ka:i:4\r\nac\rgt\n;ignored\nNN\n>b ka:i:5\nTTTT";
    let (_, records) = parse(input, 3, Some(2), 1024).unwrap();
    assert_eq!(
        records,
        vec![
            (b">a ka:i:4".to_vec(), b"ACGTNN".to_vec()),
            (b">b ka:i:5".to_vec(), b"TTTT".to_vec())
        ]
    );
    let (_, records) = parse(b">empty\n>a\nACGT\n>empty2\n", 5, None, 1024).unwrap();
    assert_eq!(records, vec![(b">a".to_vec(), b"ACGT".to_vec())]);
}

#[test]
fn parser_propagates_read_and_callback_errors() {
    struct Broken {
        data: io::Cursor<Vec<u8>>,
    }
    impl Read for Broken {
        fn read(&mut self, b: &mut [u8]) -> io::Result<usize> {
            if self.data.position() == self.data.get_ref().len() as u64 {
                Err(io::Error::other("download interrupted"))
            } else {
                self.data.read(b)
            }
        }
    }
    let error = read_fasta_with_helicase(
        BufReader::new(Broken {
            data: io::Cursor::new(b">a\nACGT\n".to_vec()),
        }),
        Some(2),
        true,
        false,
        1024,
        &mut |_| Ok(()),
    )
    .unwrap_err();
    assert_eq!(error.to_string(), "download interrupted");
    let error = read_fasta_with_helicase(
        &b">a\nACGT\n>b\nTGCA\n"[..],
        Some(2),
        true,
        false,
        1024,
        &mut |_| Err(io::Error::other("callback failed")),
    )
    .unwrap_err();
    assert_eq!(error.to_string(), "callback failed");
    assert!(parse(b"ACGT without a header", 3, Some(2), 1024).is_err());
}

#[test]
fn local_fasta_compression_formats_and_concatenated_frames_agree() {
    let dir = tempfile::tempdir().unwrap();
    let records = [
        b">a ka:i:4\nacgtNN\n".as_slice(),
        b">b ka:i:5\nTGCATGCA\n".as_slice(),
    ];
    let plain = records.concat();
    let mut gzip = Vec::new();
    let mut zstd = Vec::new();
    for r in records {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        encoder.write_all(r).unwrap();
        gzip.extend(encoder.finish().unwrap());
        zstd.extend(zstd::stream::encode_all(r, 1).unwrap());
    }
    for (suffix, data) in [("fa", plain), ("fa.gz", gzip), ("fa.zstd", zstd)] {
        let path = dir.path().join(format!("input.{suffix}"));
        std::fs::write(&path, data).unwrap();
        let mut got = Vec::new();
        SequencesReader::new().process_file_extended(
            &path,
            |s| got.push((s.ident_data.to_vec(), s.seq.to_vec())),
            Some(30),
            true,
            false,
        );
        assert_eq!(
            got,
            vec![
                (b">a ka:i:4".to_vec(), b"ACGTNN".to_vec()),
                (b">b ka:i:5".to_vec(), b"TGCATGCA".to_vec())
            ]
        );
    }
}
