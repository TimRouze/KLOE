use crate::lines_reader::LinesReader;
use config::DEFAULT_OUTPUT_BUFFER_SIZE;
use nightly_quirks::branch_pred::unlikely;
use std::path::Path;

const IDENT_STATE: usize = 0;
const SEQ_STATE: usize = 1;
const QUAL_STATE: usize = 2;

#[derive(Copy, Clone)]
pub enum DnaSequencesFileType {
    FASTA,
    FASTQ,
    GFA,
    BINARY,
}

#[derive(Copy, Clone)]
pub struct DnaSequence<'a> {
    pub ident_data: &'a [u8],
    pub seq: &'a [u8],
    pub format: DnaSequencesFileType,
}

const SEQ_LETTERS_MAPPING: [u8; 256] = {
    let mut lookup = [b'N'; 256];
    lookup[b'A' as usize] = b'A';
    lookup[b'C' as usize] = b'C';
    lookup[b'G' as usize] = b'G';
    lookup[b'T' as usize] = b'T';
    lookup[b'a' as usize] = b'A';
    lookup[b'c' as usize] = b'C';
    lookup[b'g' as usize] = b'G';
    lookup[b't' as usize] = b'T';
    lookup
};

pub struct SequencesReader {
    lines_reader: LinesReader,
}

impl SequencesReader {
    pub fn new() -> Self {
        Self {
            lines_reader: LinesReader::new(),
        }
    }

    fn normalize_sequence(seq: &mut [u8]) {
        for el in seq.iter_mut() {
            *el = SEQ_LETTERS_MAPPING[*el as usize];
        }
    }

    pub fn process_file_extended<F: FnMut(DnaSequence)>(
        &mut self,
        source: impl AsRef<Path>,
        func: F,
        line_split_copyback: Option<usize>,
        copy_ident: bool,
        remove_file: bool,
    ) {
        const FASTQ_EXTS: &[&str] = &["fq", "fastq"];
        const FASTA_EXTS: &[&str] = &["fa", "fasta", "fna", "ffn"];

        let mut file_type = None;
        let mut tmp = source.as_ref().file_name().unwrap().to_str().unwrap();
        let mut path: &Path = tmp.as_ref();

        while let Some(ext) = path.extension() {
            if FASTQ_EXTS.contains(&ext.to_str().unwrap()) {
                file_type = Some(DnaSequencesFileType::FASTQ);
                break;
            }
            if FASTA_EXTS.contains(&ext.to_str().unwrap()) {
                file_type = Some(DnaSequencesFileType::FASTA);
                break;
            }
            tmp = &tmp[0..tmp.len() - ext.len() - 1];
            path = tmp.as_ref()
        }

        match file_type {
            None => panic!(
                "Cannot recognize file type of '{}'",
                source.as_ref().display()
            ),
            Some(ftype) => match ftype {
                DnaSequencesFileType::FASTA => {
                    self.process_fasta(source, func, line_split_copyback, copy_ident, remove_file);
                }
                DnaSequencesFileType::FASTQ => {
                    self.process_fastq(source, func, remove_file);
                }
                DnaSequencesFileType::GFA => {
                    todo!()
                }
                DnaSequencesFileType::BINARY => {
                    todo!()
                }
            },
        }
    }

    fn process_fasta(
        &mut self,
        source: impl AsRef<Path>,
        mut func: impl FnMut(DnaSequence),
        line_split_copyback: Option<usize>,
        copy_ident: bool,
        remove_file: bool,
    ) {
        let path = source.as_ref();
        let file = std::fs::File::open(path)
            .unwrap_or_else(|e| panic!("Cannot open {}: {e}", path.display()));
        let reader: Box<dyn std::io::Read> = match path.extension().and_then(|s| s.to_str()) {
            Some("gz") => Box::new(flate2::read::MultiGzDecoder::new(file)),
            Some("zst" | "zstd") => Box::new(zstd::stream::read::Decoder::new(file).unwrap()),
            Some("lz4") => Box::new(lz4::Decoder::new(file).unwrap()),
            _ => Box::new(file),
        };
        crate::fasta_reader::read_fasta_with_helicase(
            std::io::BufReader::with_capacity(64 * 1024, reader),
            line_split_copyback,
            copy_ident,
            false,
            DEFAULT_OUTPUT_BUFFER_SIZE,
            &mut |sequence| { func(sequence); Ok(()) },
        ).unwrap_or_else(|e| panic!("Cannot parse {}: {e}", path.display()));
        if remove_file {
            std::fs::remove_file(path).unwrap();
        }
    }

    fn process_fastq(
        &mut self,
        source: impl AsRef<Path>,
        mut func: impl FnMut(DnaSequence),
        // get_quality: bool,
        remove_file: bool,
    ) {
        let mut state = IDENT_STATE;
        let mut skipped_plus = false;

        let mut intermediate = [Vec::new(), Vec::new(), Vec::new()];

        self.lines_reader.process_lines(
            source,
            |line: &[u8], partial, finished| {
                if unlikely(finished) {
                    return;
                }

                if state == QUAL_STATE {
                    if !skipped_plus {
                        if !partial {
                            skipped_plus = true;
                        }
                        return;
                    }

                    // if get_quality {
                    //     intermediate[state].extend_from_slice(line);
                    // }

                    if !partial {
                        Self::normalize_sequence(&mut intermediate[SEQ_STATE]);
                        func(DnaSequence {
                            ident_data: &intermediate[IDENT_STATE],
                            seq: &intermediate[SEQ_STATE],
                            // qual: if get_quality {
                            //     Some(&intermediate[QUAL_STATE])
                            // } else {
                            //     None
                            // },
                            format: DnaSequencesFileType::FASTQ,
                        });

                        intermediate[IDENT_STATE].clear();
                        intermediate[SEQ_STATE].clear();
                        intermediate[QUAL_STATE].clear();

                        skipped_plus = false;
                    }
                } else {
                    intermediate[state].extend_from_slice(line);
                }

                if !partial {
                    state = (state + 1) % 3;
                }
            },
            remove_file,
        );
    }
}
