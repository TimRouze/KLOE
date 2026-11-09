use ggcat_logging::UnrecoverableErrorLogging;

use crate::sequences_reader::{DnaSequence, SequencesReader};
use crate::sequences_stream::{GenericSequencesStream, SequenceInfo};
use std::path::PathBuf;

pub struct FastaFileSequencesStream {
    sequences_reader: SequencesReader,
}

/// Compatibility estimate used to plan buckets for file and remote FASTA inputs.
/// This is a planning heuristic, not a measured count of decoded bases.
pub fn estimated_fasta_bases_from_size(length: u64, compressed: bool) -> u64 {
    if compressed {
        (length as f64 * 0.5) as u64
    } else {
        length
    }
}

impl FastaFileSequencesStream {
    pub fn get_estimated_bases_count(file: &PathBuf) -> anyhow::Result<u64> {
        let length = std::fs::metadata(file)
            .log_unrecoverable_error_with_data("Error while opening file", file.display())?
            .len();

        let compressed = file
            .extension()
            .map(|x| x == "gz" || x == "lz4" || x == "zst" || x == "zstd")
            .unwrap_or(false);
        Ok(estimated_fasta_bases_from_size(length, compressed))
    }
}

impl GenericSequencesStream for FastaFileSequencesStream {
    type SequenceBlockData = (PathBuf, Option<u32>);

    fn new() -> Self {
        Self {
            sequences_reader: SequencesReader::new(),
        }
    }

    fn read_block(
        &mut self,
        block: &Self::SequenceBlockData,
        copy_ident_data: bool,
        partial_read_copyback: Option<usize>,
        mut callback: impl FnMut(DnaSequence, SequenceInfo),
    ) {
        self.sequences_reader.process_file_extended(
            &block.0,
            |x| callback(x, SequenceInfo { color: block.1 }),
            partial_read_copyback,
            copy_ident_data,
            false,
        );
    }
}
