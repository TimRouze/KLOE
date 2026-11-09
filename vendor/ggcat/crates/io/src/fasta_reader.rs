//! Helicase parsing of bounded, record-aligned FASTA batches.
//!
//! Passing an arbitrarily long record to Helicase's reader API would buffer the
//! entire record. Keep the established overlapping-chunk path for those inputs,
//! and for FASTA comments / embedded CRs accepted by KLOE's existing reader.
use crate::sequences_reader::{DnaSequence, DnaSequencesFileType};
use helicase::input::FromSlice;
use helicase::parser::FastaParser;
use helicase::{Config, HelicaseParser, ParserOptions};
use std::io::{self, BufRead, BufReader, Cursor, Read};

const BATCH_BYTES: usize = 256 * 1024;
const MAX_HEADER_BYTES: usize = 1024 * 1024;
const WITH_HEADERS: Config = ParserOptions::default().config();
const WITHOUT_HEADERS: Config = ParserOptions::default().ignore_headers().config();

/// Parse decompressed FASTA while preserving headers and KLOE normalization.
/// With `copyback`, oversized sequences are split with that many overlapping
/// bases. Without it, complete sequences are returned, as for the old local
/// reader. Header-only scans never buffer a complete sequence.
pub fn read_fasta_with_helicase(
    mut reader: impl BufRead,
    copyback: Option<usize>,
    copy_ident: bool,
    headers_only: bool,
    chunk_bases: usize,
    callback: &mut dyn FnMut(DnaSequence) -> io::Result<()>,
) -> io::Result<u64> {
    if headers_only {
        return read_fasta_compat(reader, copyback, copy_ident, true, chunk_bases, callback);
    }
    // A parsed batch must not exceed the caller's requested sequence chunk.
    let batch_limit = BATCH_BYTES.min(chunk_bases.max(1));
    let mut batch = Vec::with_capacity(batch_limit);
    let mut bases = 0;
    loop {
        let mut eof = false;
        while batch.len() < batch_limit {
            let input = reader.fill_buf()?;
            if input.is_empty() {
                eof = true;
                break;
            }
            let take = input.len().min(batch_limit - batch.len());
            let start = batch.len();
            batch.extend_from_slice(&input[..take]);
            reader.consume(take);
            // Deliver complete records before waiting for more network data.
            // Only scan the newly read bytes, including the preceding byte
            // when checking line starts, so tiny reads stay linear-time.
            if memchr::memchr_iter(b'>', &batch[start..])
                .any(|i| start + i > 0 && batch[start + i - 1] == b'\n')
            {
                break;
            }
        }
        if batch.is_empty() {
            return Ok(bases);
        }
        let boundary = if eof {
            batch.len()
        } else {
            memchr::memchr_iter(b'>', &batch)
                .filter(|&i| i > 0 && batch[i - 1] == b'\n')
                .last()
                .unwrap_or(0)
        };
        // Inspect only complete records: a CRLF pair in the trailing partial
        // record may straddle the next read, which is not a compatibility case.
        let complete = &batch[..boundary];
        let needs_compat = complete.first().is_some_and(|&b| b != b'>')
            || memchr::memchr_iter(b';', complete).any(|i| i == 0 || complete[i - 1] == b'\n')
            || memchr::memchr_iter(b'>', complete).any(|i| i > 0 && complete[i - 1] != b'\n')
            || memchr::memchr_iter(b'\r', complete).any(|i| complete.get(i + 1) != Some(&b'\n'));
        if needs_compat || boundary == 0 {
            let input = BufReader::new(Cursor::new(batch).chain(reader));
            return read_fasta_compat(input, copyback, copy_ident, false, chunk_bases, callback)
                .map(|count| bases + count);
        }
        bases += if copy_ident {
            parse_batch::<WITH_HEADERS>(&batch[..boundary], true, callback)?
        } else {
            parse_batch::<WITHOUT_HEADERS>(&batch[..boundary], false, callback)?
        };
        batch.copy_within(boundary.., 0);
        batch.truncate(batch.len() - boundary);
        if eof {
            return Ok(bases);
        }
    }
}

fn parse_batch<const CONFIG: Config>(
    input: &[u8],
    copy_ident: bool,
    callback: &mut dyn FnMut(DnaSequence) -> io::Result<()>,
) -> io::Result<u64> {
    // Helicase detects a line-ending width for the record. Normalize mixed
    // CRLF/LF input first so a short LF ending cannot truncate a nucleotide.
    let lf_input;
    let input = if memchr::memchr(b'\r', input).is_some() {
        lf_input = input
            .iter()
            .copied()
            .filter(|&b| b != b'\r')
            .collect::<Vec<_>>();
        lf_input.as_slice()
    } else {
        input
    };
    let mut parser = FastaParser::<CONFIG, _>::from_slice(input)
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error.to_string()))?;
    let mut header = Vec::new();
    let mut normalized = Vec::new();
    let mut bases = 0;
    while let Some(event) = parser.next() {
        event.map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error.to_string()))?;
        let sequence = parser.get_dna_string();
        if sequence.is_empty() {
            continue;
        }
        bases += sequence.len() as u64;
        if copy_ident {
            header.clear();
            header.push(b'>');
            header.extend_from_slice(parser.get_header());
            if header.len() > MAX_HEADER_BYTES {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "FASTA header exceeds 1 MiB",
                ));
            }
        }
        normalized.clear();
        normalized.extend(sequence.iter().map(|&b| DNA_MAPPING[b as usize]));
        callback(DnaSequence {
            ident_data: &header,
            seq: &normalized,
            format: DnaSequencesFileType::FASTA,
        })?;
    }
    Ok(bases)
}
const DNA_MAPPING: [u8; 256] = {
    let mut map = [b'N'; 256];
    map[b'A' as usize] = b'A';
    map[b'a' as usize] = b'A';
    map[b'C' as usize] = b'C';
    map[b'c' as usize] = b'C';
    map[b'G' as usize] = b'G';
    map[b'g' as usize] = b'G';
    map[b'T' as usize] = b'T';
    map[b't' as usize] = b'T';
    map
};

fn read_fasta_compat(
    mut reader: impl BufRead,
    copyback: Option<usize>,
    copy_ident: bool,
    headers_only: bool,
    chunk_bases: usize,
    callback: &mut dyn FnMut(DnaSequence) -> io::Result<()>,
) -> io::Result<u64> {
    let overlap = copyback.unwrap_or(0);
    let chunk_size = if copyback.is_some() {
        chunk_bases.max(1).max(overlap.saturating_mul(2))
    } else {
        usize::MAX
    };
    let mut header = Vec::with_capacity(if copy_ident { 256 } else { 0 });
    let mut sequence = Vec::with_capacity(if headers_only {
        0
    } else {
        chunk_size.min(BATCH_BYTES)
    });
    let mut line_start = true;
    let mut in_header = false;
    let mut comment = false;
    let mut have_header = false;
    let mut record_bases = 0u64;
    let mut bases = 0u64;
    loop {
        let buffer = reader.fill_buf()?;
        if buffer.is_empty() {
            break;
        }
        let newline = memchr::memchr(b'\n', buffer);
        let end = newline.unwrap_or(buffer.len());
        let mut fragment = &buffer[..end];
        // Remove the usual CRLF terminator without allocating. Embedded CRs are rare.
        if fragment.last() == Some(&b'\r') {
            fragment = &fragment[..fragment.len() - 1];
        }
        let cleaned;
        if memchr::memchr(b'\r', fragment).is_some() {
            cleaned = fragment
                .iter()
                .copied()
                .filter(|&b| b != b'\r')
                .collect::<Vec<_>>();
            fragment = &cleaned;
        }
        if !fragment.is_empty() {
            if line_start && fragment[0] == b'>' {
                if (headers_only && record_bases > 0) || !sequence.is_empty() {
                    callback(DnaSequence {
                        ident_data: &header,
                        seq: &sequence,
                        format: DnaSequencesFileType::FASTA,
                    })?;
                }
                sequence.clear();
                header.clear();
                record_bases = 0;
                have_header = true;
                in_header = true;
            } else if line_start && fragment[0] == b';' {
                comment = true;
            }
            line_start = false;
            if !comment {
                if in_header {
                    if copy_ident {
                        if header.len().saturating_add(fragment.len()) > MAX_HEADER_BYTES {
                            return Err(io::Error::new(
                                io::ErrorKind::InvalidData,
                                "FASTA header exceeds 1 MiB",
                            ));
                        }
                        header.extend_from_slice(fragment);
                    }
                } else {
                    if !have_header {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "sequence data before first FASTA header",
                        ));
                    }
                    bases += fragment.len() as u64;
                    record_bases += fragment.len() as u64;
                    if !headers_only {
                        while !fragment.is_empty() {
                            let count = fragment.len().min(chunk_size - sequence.len());
                            sequence
                                .extend(fragment[..count].iter().map(|&b| DNA_MAPPING[b as usize]));
                            fragment = &fragment[count..];
                            if copyback.is_some() && sequence.len() == chunk_size {
                                callback(DnaSequence {
                                    ident_data: &header,
                                    seq: &sequence,
                                    format: DnaSequencesFileType::FASTA,
                                })?;
                                let start = sequence.len() - overlap;
                                sequence.copy_within(start.., 0);
                                sequence.truncate(overlap);
                            }
                        }
                    }
                }
            }
        }
        let consumed = end + usize::from(newline.is_some());
        reader.consume(consumed);
        if newline.is_some() {
            line_start = true;
            in_header = false;
            comment = false;
        }
    }
    if (headers_only && record_bases > 0) || !sequence.is_empty() {
        callback(DnaSequence {
            ident_data: &header,
            seq: &sequence,
            format: DnaSequencesFileType::FASTA,
        })?;
    }
    Ok(bases)
}
