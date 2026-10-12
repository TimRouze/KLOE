use std::collections::HashMap;
use std::fs::{self, File};
use std::hash::{Hash, Hasher};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

use anyhow::{bail, ensure, Context, Result};
use rammap::{dp_global, encode_nt4, Aligner, DpScoring, Preset, Strand};
use rayon::prelude::*;

const MAGIC: &[u8; 8] = b"KLOER001";
const MULTI_MAGIC: &[u8; 8] = b"KLOER002";
const DEDUP_MULTI_MAGIC: &[u8; 8] = b"KLOER003";
const COMPACT_DEDUP_MULTI_MAGIC: &[u8; 8] = b"KLOER004";
const COMPACT_HAMMING_MULTI_MAGIC: &[u8; 8] = b"KLOER005";
const ENTRY_NEW: u8 = 0;
const ENTRY_REFERENCE: u8 = 1;
const RECORD_MAPPED: u8 = 1;
const RECORD_REVERSE: u8 = 2;
const RECORD_HAMMING: u8 = 4;
const OP_COPY: u8 = 0;
const OP_SUBSTITUTE: u8 = 1;
const OP_INSERT: u8 = 2;
const OP_DELETE: u8 = 3;

#[derive(Clone, Debug)]
struct FastaRecord {
    name: String,
    sequence: Vec<u8>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct SequenceDigest([u8; 32]);

// BLAKE3 output is already uniformly distributed. Hashing all 32 bytes again
// with HashMap's keyed hasher adds work without improving bucket quality; use
// one 64-bit lane for table placement while retaining all 256 bits for equality.
impl Hash for SequenceDigest {
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_u64(u64::from_le_bytes(
            self.0[..8].try_into().expect("eight-byte digest prefix"),
        ));
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum EditOp {
    Copy(usize),
    Substitute(Vec<u8>),
    Insert(Vec<u8>),
    Delete(usize),
}

#[derive(Clone, Debug)]
enum EncodedRecord {
    Mapped {
        target_id: usize,
        target_start: usize,
        target_len: usize,
        reverse: bool,
        query_len: usize,
        ops: Vec<EditOp>,
    },
    Literal(Vec<u8>),
}

#[derive(Clone, Copy, Debug)]
enum RecordPlan {
    New { id: usize, local_index: usize },
    Reference { id: usize },
}

#[derive(Clone, Copy, Debug, Default)]
struct RecordContribution {
    mapped: bool,
    bases: u64,
    divergence_bases: u64,
}

const EXACT_SEED_LENGTH: usize = 21;

#[derive(Debug)]
struct ExactReferenceIndex {
    /// `(two-bit seed, target_id << 32 | start)` sorted by seed.
    seeds: Vec<(u64, u64)>,
}

impl ExactReferenceIndex {
    fn new(references: &[FastaRecord]) -> Result<Self> {
        let total_bases = references
            .iter()
            .map(|reference| reference.sequence.len())
            .sum::<usize>();
        let mut seeds = Vec::with_capacity(total_bases);
        let mask = (1u64 << (2 * EXACT_SEED_LENGTH)) - 1;
        for (target_id, reference) in references.iter().enumerate() {
            let target_id = u32::try_from(target_id).context("too many reference sequences")?;
            let mut key = 0u64;
            let mut valid_run = 0usize;
            for (index, base) in reference.sequence.iter().enumerate() {
                let Some(code) = two_bit_code(*base) else {
                    key = 0;
                    valid_run = 0;
                    continue;
                };
                key = ((key << 2) | code) & mask;
                valid_run += 1;
                if valid_run >= EXACT_SEED_LENGTH {
                    let start = index + 1 - EXACT_SEED_LENGTH;
                    let start = u32::try_from(start).context("reference sequence exceeds 4 Gbp")?;
                    let location = (u64::from(target_id) << 32) | u64::from(start);
                    seeds.push((key, location));
                }
            }
        }
        seeds.par_sort_unstable_by_key(|entry| entry.0);
        Ok(Self { seeds })
    }

    /// Find a same-length reference occurrence supported by an exact seed and
    /// minimize its Hamming distance. Most assembly unitigs differ from their
    /// reference locus only by sparse substitutions; resolving those here
    /// avoids invoking the substantially more expensive general mapper.
    fn find_hamming(
        &self,
        query: &[u8],
        references: &[FastaRecord],
        max_divergence: f64,
    ) -> Option<(usize, usize, usize)> {
        if query.len() < EXACT_SEED_LENGTH {
            return None;
        }
        let maximum_mismatches = (query.len() as f64 * max_divergence).floor() as usize;
        // At <=1% divergence an exact seed overwhelmingly identifies the
        // intended assembly locus. Accept it immediately instead of probing
        // two more seeds (and usually the reverse complement) for a marginally
        // better placement that does not improve reconstruction.
        let confident_mismatches = query.len() / 100;
        let last_offset = query.len() - EXACT_SEED_LENGTH;
        let offsets = [0, last_offset / 2, last_offset];
        let mut best = None::<(usize, usize, usize)>;
        let mut examined = 0usize;
        const MAX_CANDIDATES: usize = 8_192;

        for (offset_index, &query_offset) in offsets.iter().enumerate() {
            if offsets[..offset_index].contains(&query_offset) {
                continue;
            }
            let key = encode_seed(&query[query_offset..query_offset + EXACT_SEED_LENGTH])?;
            let start = self.seeds.partition_point(|entry| entry.0 < key);
            let end = self.seeds.partition_point(|entry| entry.0 <= key);
            for &(_, location) in &self.seeds[start..end] {
                examined += 1;
                if examined > MAX_CANDIDATES {
                    return best;
                }
                let target_id = (location >> 32) as usize;
                let seed_start = (location as u32) as usize;
                let Some(target_start) = seed_start.checked_sub(query_offset) else {
                    continue;
                };
                let Some(target_end) = target_start.checked_add(query.len()) else {
                    continue;
                };
                let Some(target) = references
                    .get(target_id)
                    .and_then(|reference| reference.sequence.get(target_start..target_end))
                else {
                    continue;
                };
                let cutoff = best
                    .map(|(_, _, mismatches)| mismatches.saturating_sub(1))
                    .unwrap_or(maximum_mismatches);
                let mut mismatches = 0usize;
                for (&query_base, &target_base) in query.iter().zip(target) {
                    mismatches += usize::from(query_base != target_base);
                    if mismatches > cutoff {
                        break;
                    }
                }
                if mismatches <= cutoff {
                    best = Some((target_id, target_start, mismatches));
                    if mismatches <= confident_mismatches {
                        return best;
                    }
                }
            }
        }
        best
    }
}

#[derive(Clone, Debug)]
pub struct CompressOptions {
    pub threads: usize,
    pub max_divergence: f64,
    pub min_coverage: f64,
    pub bandwidth: i32,
    pub compression_level: i32,
}

impl Default for CompressOptions {
    fn default() -> Self {
        Self {
            threads: std::thread::available_parallelism().map_or(1, usize::from),
            max_divergence: 0.20,
            min_coverage: 0.80,
            bandwidth: 512,
            compression_level: 3,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct CompressionStats {
    pub records: usize,
    pub mapped: usize,
    pub literal: usize,
    pub input_bases: u64,
    pub mapped_bases: u64,
    pub literal_bases: u64,
    pub divergence_bases: u64,
    pub unique_records: usize,
    pub unique_mapped: usize,
    pub unique_literal: usize,
    pub unique_bases: u64,
    pub unique_mapped_bases: u64,
    pub unique_literal_bases: u64,
}

#[derive(Clone, Debug, Default)]
pub struct DecompressionStats {
    pub records: usize,
    pub bases: u64,
}

pub fn compress(
    reference_path: &Path,
    input_path: &Path,
    output_path: &Path,
    options: &CompressOptions,
) -> Result<CompressionStats> {
    validate_compress_options(options)?;

    let references = read_fasta(reference_path)
        .with_context(|| format!("read reference '{}'", reference_path.display()))?;
    ensure!(!references.is_empty(), "reference FASTA is empty");
    let unitigs = read_fasta_sequences(input_path)
        .with_context(|| format!("read unitigs '{}'", input_path.display()))?;

    let exact_index = ExactReferenceIndex::new(&references)?;
    let aligner = build_aligner(reference_path)?;
    let pool = build_mapping_pool(options.threads)?;
    let encoded = pool.install(|| {
        unitigs
            .par_iter()
            .map(|record| encode_record(record, &references, &exact_index, &aligner, options))
            .collect::<Vec<_>>()
    });

    let output = BufWriter::new(
        File::create(output_path)
            .with_context(|| format!("create archive '{}'", output_path.display()))?,
    );
    let mut output = output;
    output.write_all(MAGIC)?;
    let mut encoder = zstd::Encoder::new(output, options.compression_level)?;
    let fingerprint = reference_fingerprint(&references);
    encoder.write_all(fingerprint.as_bytes())?;
    encoder.write_all(sequence_fingerprint(&unitigs).as_bytes())?;
    write_varint(references.len() as u64, &mut encoder)?;
    for reference in &references {
        write_bytes(reference.name.as_bytes(), &mut encoder)?;
        write_varint(reference.sequence.len() as u64, &mut encoder)?;
    }
    write_varint(encoded.len() as u64, &mut encoder)?;

    let mut stats = CompressionStats::default();
    for record in encoded {
        update_compression_stats(&mut stats, &record);
        update_unique_compression_stats(&mut stats, &record);
        write_record(&record, &mut encoder)?;
    }
    let mut output = encoder.finish()?;
    output.flush()?;
    Ok(stats)
}

/// Compress every FASTA listed in `input_list_path` into one reference-relative
/// archive. Each non-empty, non-comment line names one dataset FASTA. Relative
/// paths are resolved from the list file's parent directory.
pub fn compress_list(
    reference_path: &Path,
    input_list_path: &Path,
    output_path: &Path,
    options: &CompressOptions,
) -> Result<CompressionStats> {
    validate_compress_options(options)?;
    let input_paths = read_input_list(input_list_path)?;
    ensure!(
        !input_paths.is_empty(),
        "input list contains no FASTA files"
    );

    let references = read_fasta(reference_path)
        .with_context(|| format!("read reference '{}'", reference_path.display()))?;
    ensure!(!references.is_empty(), "reference FASTA is empty");
    let exact_index = ExactReferenceIndex::new(&references)?;
    let aligner = build_aligner(reference_path)?;
    let pool = build_mapping_pool(options.threads)?;

    let output = BufWriter::new(
        File::create(output_path)
            .with_context(|| format!("create archive '{}'", output_path.display()))?,
    );
    let mut output = output;
    output.write_all(COMPACT_HAMMING_MULTI_MAGIC)?;
    let mut encoder = zstd::Encoder::new(output, options.compression_level)?;
    // The encoded record stream for a bacterial cohort is small enough for a
    // shared 128 MiB history. This captures repeated literal/edit records across
    // datasets without retaining all unitigs in KLOER's own memory.
    encoder.window_log(27)?;
    encoder.long_distance_matching(true)?;
    write_reference_header(&references, &mut encoder)?;
    write_varint(input_paths.len() as u64, &mut encoder)?;

    let mut stats = CompressionStats::default();
    let mut sequence_ids = HashMap::<SequenceDigest, usize>::new();
    let mut contributions = Vec::<RecordContribution>::new();
    let mut completed_datasets = 0usize;
    // A human unitig FASTA can occupy several GiB after decompression. Mapping
    // one dataset at a time keeps the resident input and encoded-record waves
    // bounded without reducing parallelism inside a dataset.
    const DATASET_BATCH_SIZE: usize = 1;
    for input_batch in input_paths.chunks(DATASET_BATCH_SIZE) {
        let datasets = input_batch
            .iter()
            .map(|input_path| {
                read_fasta_sequences(input_path)
                    .with_context(|| format!("read unitigs '{}'", input_path.display()))
            })
            .collect::<Result<Vec<_>>>()?;
        let mut plans_by_dataset = Vec::with_capacity(datasets.len());
        let mut fingerprints_by_dataset = Vec::with_capacity(datasets.len());
        let mut unique_inputs = Vec::<(usize, &FastaRecord)>::new();
        for unitigs in &datasets {
            let mut plans = Vec::with_capacity(unitigs.len());
            let mut fingerprint = blake3::Hasher::new();
            for record in unitigs {
                update_sequence_fingerprint(&mut fingerprint, &record.sequence);
                let hash = SequenceDigest(*blake3::hash(&record.sequence).as_bytes());
                if let Some(&id) = sequence_ids.get(&hash) {
                    plans.push(RecordPlan::Reference { id });
                } else {
                    let id = contributions.len();
                    let local_index = unique_inputs.len();
                    sequence_ids.insert(hash, id);
                    contributions.push(RecordContribution::default());
                    unique_inputs.push((id, record));
                    plans.push(RecordPlan::New { id, local_index });
                }
            }
            plans_by_dataset.push(plans);
            fingerprints_by_dataset.push(fingerprint.finalize());
        }
        let encoded_unique = pool.install(|| {
            unique_inputs
                .par_iter()
                .map(|(_, record)| {
                    encode_record(record, &references, &exact_index, &aligner, options)
                })
                .collect::<Vec<_>>()
        });
        for ((id, _), record) in unique_inputs.iter().zip(&encoded_unique) {
            contributions[*id] = record_contribution(record);
            update_unique_compression_stats(&mut stats, record);
        }

        for ((input_path, plans), fingerprint) in input_batch
            .iter()
            .zip(&plans_by_dataset)
            .zip(&fingerprints_by_dataset)
        {
            let dataset_name = input_path
                .file_name()
                .and_then(|name| name.to_str())
                .context("dataset filename is not valid UTF-8")?;
            write_bytes(dataset_name.as_bytes(), &mut encoder)?;
            encoder.write_all(fingerprint.as_bytes())?;
            write_varint(plans.len() as u64, &mut encoder)?;
            let mut previous_sequence_id = 0usize;
            for plan in plans {
                let id = match *plan {
                    RecordPlan::New { id, local_index } => {
                        write_varint(0, &mut encoder)?;
                        write_record_compact(&encoded_unique[local_index], &mut encoder)?;
                        id
                    }
                    RecordPlan::Reference { id } => {
                        write_sequence_reference(id, previous_sequence_id, &mut encoder)?;
                        id
                    }
                };
                previous_sequence_id = id;
                apply_record_contribution(&mut stats, contributions[id]);
            }
        }
        completed_datasets += input_batch.len();
        eprintln!(
            "KLOER progress: datasets={completed_datasets}/{} records={} input_bases={} unique_records={}",
            input_paths.len(),
            stats.records,
            stats.input_bases,
            stats.unique_records,
        );
    }
    let mut output = encoder.finish()?;
    output.flush()?;
    Ok(stats)
}

fn validate_compress_options(options: &CompressOptions) -> Result<()> {
    ensure!(options.threads > 0, "thread count must be positive");
    ensure!(
        (0.0..=1.0).contains(&options.max_divergence),
        "max divergence must be in [0, 1]"
    );
    ensure!(
        (0.0..=1.0).contains(&options.min_coverage),
        "minimum coverage must be in [0, 1]"
    );
    Ok(())
}

fn build_aligner(reference_path: &Path) -> Result<Aligner> {
    Aligner::from_fasta(
        reference_path
            .to_str()
            .context("reference path is not valid UTF-8")?,
        Preset::Asm20,
    )
    .with_context(|| format!("build rammap index for '{}'", reference_path.display()))
}

fn build_mapping_pool(threads: usize) -> Result<rayon::ThreadPool> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .context("create KLOER mapping thread pool")
}

fn write_reference_header(references: &[FastaRecord], writer: &mut impl Write) -> Result<()> {
    writer.write_all(reference_fingerprint(references).as_bytes())?;
    write_varint(references.len() as u64, writer)?;
    for reference in references {
        write_bytes(reference.name.as_bytes(), writer)?;
        write_varint(reference.sequence.len() as u64, writer)?;
    }
    Ok(())
}

fn update_compression_stats(stats: &mut CompressionStats, record: &EncodedRecord) {
    apply_record_contribution(stats, record_contribution(record));
}

fn record_contribution(record: &EncodedRecord) -> RecordContribution {
    match record {
        EncodedRecord::Mapped { query_len, ops, .. } => {
            let divergence_bases = ops
                .iter()
                .map(|op| match op {
                    EditOp::Substitute(sequence) | EditOp::Insert(sequence) => {
                        sequence.len() as u64
                    }
                    EditOp::Delete(length) => *length as u64,
                    EditOp::Copy(_) => 0,
                })
                .sum::<u64>();
            RecordContribution {
                mapped: true,
                bases: *query_len as u64,
                divergence_bases,
            }
        }
        EncodedRecord::Literal(sequence) => RecordContribution {
            mapped: false,
            bases: sequence.len() as u64,
            divergence_bases: 0,
        },
    }
}

fn apply_record_contribution(stats: &mut CompressionStats, contribution: RecordContribution) {
    stats.records += 1;
    stats.input_bases += contribution.bases;
    if contribution.mapped {
        stats.mapped += 1;
        stats.mapped_bases += contribution.bases;
        stats.divergence_bases += contribution.divergence_bases;
    } else {
        stats.literal += 1;
        stats.literal_bases += contribution.bases;
    }
}

fn update_unique_compression_stats(stats: &mut CompressionStats, record: &EncodedRecord) {
    let contribution = record_contribution(record);
    stats.unique_records += 1;
    stats.unique_bases += contribution.bases;
    if contribution.mapped {
        stats.unique_mapped += 1;
        stats.unique_mapped_bases += contribution.bases;
    } else {
        stats.unique_literal += 1;
        stats.unique_literal_bases += contribution.bases;
    }
}

pub fn decompress(
    reference_path: &Path,
    input_path: &Path,
    output_path: &Path,
) -> Result<DecompressionStats> {
    let references = read_fasta(reference_path)
        .with_context(|| format!("read reference '{}'", reference_path.display()))?;
    let mut input = BufReader::new(
        File::open(input_path)
            .with_context(|| format!("open archive '{}'", input_path.display()))?,
    );
    let mut magic = [0u8; MAGIC.len()];
    input.read_exact(&mut magic)?;
    ensure!(&magic == MAGIC, "input is not a KLOER archive");
    let mut decoder = zstd::Decoder::new(input)?;
    let mut stored_fingerprint = [0u8; 32];
    decoder.read_exact(&mut stored_fingerprint)?;
    ensure!(
        stored_fingerprint == *reference_fingerprint(&references).as_bytes(),
        "reference fingerprint differs from the reference used for compression"
    );
    let mut expected_sequences_fingerprint = [0u8; 32];
    decoder.read_exact(&mut expected_sequences_fingerprint)?;
    let reference_count = read_varint(&mut decoder)? as usize;
    ensure!(
        reference_count == references.len(),
        "reference sequence count differs from archive"
    );
    for reference in &references {
        let name = read_bytes(&mut decoder)?;
        let length = read_varint(&mut decoder)? as usize;
        ensure!(
            name == reference.name.as_bytes(),
            "reference name differs from archive"
        );
        ensure!(
            length == reference.sequence.len(),
            "reference length differs from archive"
        );
    }
    let record_count = read_varint(&mut decoder)? as usize;
    let mut output = BufWriter::new(
        File::create(output_path)
            .with_context(|| format!("create FASTA '{}'", output_path.display()))?,
    );
    let mut stats = DecompressionStats::default();
    let mut sequences_hasher = blake3::Hasher::new();
    for ordinal in 0..record_count {
        let record = read_record(&mut decoder)?;
        let sequence = decode_record(&record, &references)?;
        update_sequence_fingerprint(&mut sequences_hasher, &sequence);
        writeln!(output, ">{ordinal}")?;
        output.write_all(&sequence)?;
        output.write_all(b"\n")?;
        stats.records += 1;
        stats.bases += sequence.len() as u64;
    }
    ensure!(
        *sequences_hasher.finalize().as_bytes() == expected_sequences_fingerprint,
        "decompressed sequence checksum differs from archive"
    );
    output.flush()?;
    Ok(stats)
}

/// Decompress a multi-dataset KLOER archive into one FASTA per original input
/// dataset. Output names are ordinal-prefixed to remain unique and path-safe.
pub fn decompress_list(
    reference_path: &Path,
    input_path: &Path,
    output_dir: &Path,
) -> Result<DecompressionStats> {
    let references = read_fasta(reference_path)
        .with_context(|| format!("read reference '{}'", reference_path.display()))?;
    let mut input = BufReader::new(
        File::open(input_path)
            .with_context(|| format!("open archive '{}'", input_path.display()))?,
    );
    let mut magic = [0u8; MULTI_MAGIC.len()];
    input.read_exact(&mut magic)?;
    let compact_deduplicated =
        &magic == COMPACT_DEDUP_MULTI_MAGIC || &magic == COMPACT_HAMMING_MULTI_MAGIC;
    let deduplicated = &magic == DEDUP_MULTI_MAGIC || compact_deduplicated;
    ensure!(
        &magic == MULTI_MAGIC || deduplicated,
        "input is not a multi-dataset KLOER archive"
    );
    let mut decoder = zstd::Decoder::new(input)?;
    read_and_check_reference_header(&mut decoder, &references)?;
    let dataset_count = read_usize(&mut decoder, "dataset count")?;
    ensure!(
        dataset_count <= 1_000_000,
        "archive dataset count is unreasonable"
    );
    fs::create_dir_all(output_dir)
        .with_context(|| format!("create output directory '{}'", output_dir.display()))?;

    let mut stats = DecompressionStats::default();
    let mut unique_sequences = Vec::<Vec<u8>>::new();
    for dataset_ordinal in 0..dataset_count {
        let stored_name = String::from_utf8(read_bytes(&mut decoder)?)
            .context("dataset name in archive is not valid UTF-8")?;
        let mut expected_fingerprint = [0u8; 32];
        decoder.read_exact(&mut expected_fingerprint)?;
        let record_count = read_usize(&mut decoder, "dataset record count")?;
        ensure!(
            record_count <= 1_000_000_000,
            "archive dataset record count is unreasonable"
        );
        let output_name = format!(
            "{dataset_ordinal:06}_{}",
            sanitize_dataset_name(&stored_name)
        );
        let output_path = output_dir.join(output_name);
        let mut output = BufWriter::new(
            File::create(&output_path)
                .with_context(|| format!("create FASTA '{}'", output_path.display()))?,
        );
        let mut hasher = blake3::Hasher::new();
        let mut previous_sequence_id = 0usize;
        for record_ordinal in 0..record_count {
            if deduplicated {
                let id = if compact_deduplicated {
                    let token = read_varint(&mut decoder)?;
                    if token == 0 {
                        let record = read_record(&mut decoder)?;
                        unique_sequences.push(decode_record(&record, &references)?);
                        unique_sequences.len() - 1
                    } else {
                        read_sequence_reference(token, previous_sequence_id)?
                    }
                } else {
                    let mut entry_kind = [0u8; 1];
                    decoder.read_exact(&mut entry_kind)?;
                    match entry_kind[0] {
                        ENTRY_NEW => {
                            let record = read_record(&mut decoder)?;
                            unique_sequences.push(decode_record(&record, &references)?);
                            unique_sequences.len() - 1
                        }
                        ENTRY_REFERENCE => read_usize(&mut decoder, "sequence reference ID")?,
                        other => bail!("invalid sequence entry kind {other}"),
                    }
                };
                previous_sequence_id = id;
                let sequence = unique_sequences
                    .get(id)
                    .context("sequence reference ID is invalid")?;
                write_decoded_record(
                    &mut output,
                    &mut hasher,
                    sequence,
                    record_ordinal,
                    &mut stats,
                )?;
            } else {
                let record = read_record(&mut decoder)?;
                let sequence = decode_record(&record, &references)?;
                write_decoded_record(
                    &mut output,
                    &mut hasher,
                    &sequence,
                    record_ordinal,
                    &mut stats,
                )?;
            }
        }
        ensure!(
            *hasher.finalize().as_bytes() == expected_fingerprint,
            "decompressed sequence checksum differs for dataset '{stored_name}'"
        );
        output.flush()?;
    }
    Ok(stats)
}

fn write_decoded_record(
    output: &mut impl Write,
    hasher: &mut blake3::Hasher,
    sequence: &[u8],
    ordinal: usize,
    stats: &mut DecompressionStats,
) -> io::Result<()> {
    update_sequence_fingerprint(hasher, sequence);
    writeln!(output, ">{ordinal}")?;
    output.write_all(sequence)?;
    output.write_all(b"\n")?;
    stats.records += 1;
    stats.bases += sequence.len() as u64;
    Ok(())
}

fn read_and_check_reference_header(
    reader: &mut impl Read,
    references: &[FastaRecord],
) -> Result<()> {
    let mut stored_fingerprint = [0u8; 32];
    reader.read_exact(&mut stored_fingerprint)?;
    ensure!(
        stored_fingerprint == *reference_fingerprint(references).as_bytes(),
        "reference fingerprint differs from the reference used for compression"
    );
    let reference_count = read_usize(reader, "reference sequence count")?;
    ensure!(
        reference_count == references.len(),
        "reference sequence count differs from archive"
    );
    for reference in references {
        let name = read_bytes(reader)?;
        let length = read_usize(reader, "reference sequence length")?;
        ensure!(
            name == reference.name.as_bytes(),
            "reference name differs from archive"
        );
        ensure!(
            length == reference.sequence.len(),
            "reference length differs from archive"
        );
    }
    Ok(())
}

fn sanitize_dataset_name(name: &str) -> String {
    let sanitized = name
        .chars()
        .map(|character| {
            if character.is_ascii_alphanumeric() || matches!(character, '.' | '-' | '_') {
                character
            } else {
                '_'
            }
        })
        .collect::<String>();
    if sanitized.is_empty() || sanitized == "." || sanitized == ".." {
        "dataset.fa".to_string()
    } else {
        sanitized
    }
}

fn encode_record(
    record: &FastaRecord,
    references: &[FastaRecord],
    exact_index: &ExactReferenceIndex,
    aligner: &Aligner,
    options: &CompressOptions,
) -> EncodedRecord {
    let forward_match =
        exact_index.find_hamming(&record.sequence, references, options.max_divergence);
    if let Some((target_id, target_start, mismatches)) = forward_match {
        if mismatches <= record.sequence.len() / 100 {
            return encode_hamming_match(
                &record.sequence,
                target_id,
                target_start,
                false,
                references,
            );
        }
    }
    let reverse_sequence = reverse_complement(&record.sequence);
    let reverse_match =
        exact_index.find_hamming(&reverse_sequence, references, options.max_divergence);
    let hamming_match = match (forward_match, reverse_match) {
        (Some(forward), Some(reverse)) if reverse.2 < forward.2 => Some((reverse, true)),
        (Some(forward), _) => Some((forward, false)),
        (None, Some(reverse)) => Some((reverse, true)),
        (None, None) => None,
    };
    if let Some(((target_id, target_start, _), reverse)) = hamming_match {
        return encode_hamming_match(
            if reverse {
                &reverse_sequence
            } else {
                &record.sequence
            },
            target_id,
            target_start,
            reverse,
            references,
        );
    }
    let mappings = aligner.map_seq(&record.name, &record.sequence).mappings;
    let mapping = mappings
        .iter()
        .filter(|mapping| mapping.is_primary)
        .max_by_key(|mapping| mapping.score)
        .or_else(|| mappings.iter().max_by_key(|mapping| mapping.score));
    let Some(mapping) = mapping else {
        return EncodedRecord::Literal(record.sequence.clone());
    };
    let coverage = mapping.query_end.saturating_sub(mapping.query_start) as f64
        / record.sequence.len().max(1) as f64;
    if coverage < options.min_coverage || mapping.divergence > options.max_divergence {
        return EncodedRecord::Literal(record.sequence.clone());
    }
    let Some(reference) = references.get(mapping.target_id) else {
        return EncodedRecord::Literal(record.sequence.clone());
    };
    let Some(target) = reference
        .sequence
        .get(mapping.target_start..mapping.target_end)
    else {
        return EncodedRecord::Literal(record.sequence.clone());
    };
    let reverse = mapping.strand == Strand::Reverse;
    let oriented = if reverse {
        reverse_sequence
    } else {
        record.sequence.clone()
    };
    // rammap already performs base-level alignment and normally returns the
    // exact CIGAR. Reuse it instead of repeating the same dynamic program.
    // Partial/local CIGARs fall back to a full global alignment below.
    let ops = mapping
        .cigar
        .as_deref()
        .and_then(|cigar| edit_ops_from_cigar(cigar, &oriented, target).ok())
        .or_else(|| {
            let alignment = dp_global(
                &encode_nt4(&oriented),
                &encode_nt4(target),
                &DpScoring::default(),
                options.bandwidth,
            );
            edit_ops_from_cigar(&alignment.cigar, &oriented, target).ok()
        });
    let Some(ops) = ops else {
        return EncodedRecord::Literal(record.sequence.clone());
    };
    let candidate = EncodedRecord::Mapped {
        target_id: mapping.target_id,
        target_start: mapping.target_start,
        target_len: mapping.target_end - mapping.target_start,
        reverse,
        query_len: record.sequence.len(),
        ops,
    };
    match decode_record(&candidate, references) {
        Ok(decoded) if decoded == record.sequence => candidate,
        _ => EncodedRecord::Literal(record.sequence.clone()),
    }
}

fn encode_hamming_match(
    oriented_query: &[u8],
    target_id: usize,
    target_start: usize,
    reverse: bool,
    references: &[FastaRecord],
) -> EncodedRecord {
    let target = &references[target_id].sequence[target_start..target_start + oriented_query.len()];
    EncodedRecord::Mapped {
        target_id,
        target_start,
        target_len: oriented_query.len(),
        reverse,
        query_len: oriented_query.len(),
        ops: edit_ops_from_equal_length(oriented_query, target),
    }
}

fn edit_ops_from_equal_length(query: &[u8], target: &[u8]) -> Vec<EditOp> {
    debug_assert_eq!(query.len(), target.len());
    let mut ops = Vec::new();
    for (&query_base, &target_base) in query.iter().zip(target) {
        if query_base == target_base {
            push_copy(&mut ops, 1);
        } else {
            push_literal(&mut ops, OP_SUBSTITUTE, &[query_base]);
        }
    }
    ops
}

fn edit_ops_from_cigar(cigar: &str, query: &[u8], target: &[u8]) -> Result<Vec<EditOp>> {
    let mut query_cursor = 0usize;
    let mut target_cursor = 0usize;
    let mut ops = Vec::new();
    for (length, operation) in parse_cigar(cigar)? {
        match operation {
            'M' | '=' | 'X' => {
                ensure!(
                    query_cursor + length <= query.len() && target_cursor + length <= target.len(),
                    "CIGAR match exceeds sequence"
                );
                for offset in 0..length {
                    if query[query_cursor + offset] == target[target_cursor + offset] {
                        push_copy(&mut ops, 1);
                    } else {
                        push_literal(&mut ops, OP_SUBSTITUTE, &[query[query_cursor + offset]]);
                    }
                }
                query_cursor += length;
                target_cursor += length;
            }
            'I' => {
                ensure!(
                    query_cursor + length <= query.len(),
                    "CIGAR insertion exceeds query"
                );
                push_literal(
                    &mut ops,
                    OP_INSERT,
                    &query[query_cursor..query_cursor + length],
                );
                query_cursor += length;
            }
            'D' => {
                ensure!(
                    target_cursor + length <= target.len(),
                    "CIGAR deletion exceeds target"
                );
                push_delete(&mut ops, length);
                target_cursor += length;
            }
            _ => bail!("unsupported global-alignment CIGAR operation '{operation}'"),
        }
    }
    ensure!(
        query_cursor == query.len(),
        "CIGAR does not consume the full query"
    );
    ensure!(
        target_cursor == target.len(),
        "CIGAR does not consume the full target"
    );
    Ok(ops)
}

fn parse_cigar(cigar: &str) -> Result<Vec<(usize, char)>> {
    let mut result = Vec::new();
    let mut length = 0usize;
    for character in cigar.chars() {
        if let Some(digit) = character.to_digit(10) {
            length = length
                .checked_mul(10)
                .and_then(|value| value.checked_add(digit as usize))
                .context("CIGAR length overflow")?;
        } else {
            ensure!(length > 0, "CIGAR operation has zero or missing length");
            result.push((length, character));
            length = 0;
        }
    }
    ensure!(length == 0, "CIGAR ends with a length but no operation");
    Ok(result)
}

fn push_copy(ops: &mut Vec<EditOp>, length: usize) {
    if length == 0 {
        return;
    }
    if let Some(EditOp::Copy(previous)) = ops.last_mut() {
        *previous += length;
    } else {
        ops.push(EditOp::Copy(length));
    }
}

fn push_delete(ops: &mut Vec<EditOp>, length: usize) {
    if let Some(EditOp::Delete(previous)) = ops.last_mut() {
        *previous += length;
    } else {
        ops.push(EditOp::Delete(length));
    }
}

fn push_literal(ops: &mut Vec<EditOp>, kind: u8, sequence: &[u8]) {
    let matching = match (kind, ops.last_mut()) {
        (OP_SUBSTITUTE, Some(EditOp::Substitute(previous))) => Some(previous),
        (OP_INSERT, Some(EditOp::Insert(previous))) => Some(previous),
        _ => None,
    };
    if let Some(previous) = matching {
        previous.extend_from_slice(sequence);
    } else if kind == OP_SUBSTITUTE {
        ops.push(EditOp::Substitute(sequence.to_vec()));
    } else {
        ops.push(EditOp::Insert(sequence.to_vec()));
    }
}

fn decode_record(record: &EncodedRecord, references: &[FastaRecord]) -> Result<Vec<u8>> {
    match record {
        EncodedRecord::Literal(sequence) => Ok(sequence.clone()),
        EncodedRecord::Mapped {
            target_id,
            target_start,
            target_len,
            reverse,
            query_len,
            ops,
        } => {
            let reference = references
                .get(*target_id)
                .context("mapped target ID is invalid")?;
            let target_end = target_start
                .checked_add(*target_len)
                .context("mapped target coordinate overflow")?;
            let target = reference
                .sequence
                .get(*target_start..target_end)
                .context("mapped target coordinates are invalid")?;
            let mut target_cursor = 0usize;
            let mut sequence = Vec::with_capacity(*query_len);
            for op in ops {
                match op {
                    EditOp::Copy(length) => {
                        let end = target_cursor
                            .checked_add(*length)
                            .context("copy cursor overflow")?;
                        sequence.extend_from_slice(
                            target
                                .get(target_cursor..end)
                                .context("copy exceeds target span")?,
                        );
                        target_cursor = end;
                    }
                    EditOp::Substitute(bases) => {
                        let end = target_cursor
                            .checked_add(bases.len())
                            .context("substitution cursor overflow")?;
                        ensure!(end <= target.len(), "substitution exceeds target span");
                        sequence.extend_from_slice(bases);
                        target_cursor = end;
                    }
                    EditOp::Insert(bases) => sequence.extend_from_slice(bases),
                    EditOp::Delete(length) => {
                        target_cursor = target_cursor
                            .checked_add(*length)
                            .context("deletion cursor overflow")?;
                        ensure!(
                            target_cursor <= target.len(),
                            "deletion exceeds target span"
                        );
                    }
                }
            }
            ensure!(
                target_cursor == target.len(),
                "edit script does not consume target span"
            );
            ensure!(
                sequence.len() == *query_len,
                "edit script reconstructs wrong query length"
            );
            if *reverse {
                Ok(reverse_complement(&sequence))
            } else {
                Ok(sequence)
            }
        }
    }
}

fn write_record_compact(record: &EncodedRecord, writer: &mut impl Write) -> io::Result<()> {
    let EncodedRecord::Mapped {
        target_id,
        target_start,
        target_len,
        reverse,
        query_len,
        ops,
    } = record
    else {
        return write_record(record, writer);
    };
    let Some(mismatch_count) = hamming_mismatch_count(ops, *target_len, *query_len) else {
        return write_record(record, writer);
    };

    writer
        .write_all(&[RECORD_MAPPED | RECORD_HAMMING | if *reverse { RECORD_REVERSE } else { 0 }])?;
    write_varint(*target_id as u64, writer)?;
    write_varint(*target_start as u64, writer)?;
    write_varint(*query_len as u64, writer)?;
    write_varint(mismatch_count as u64, writer)?;

    let mut position = 0usize;
    let mut previous_mismatch = 0usize;
    for op in ops {
        match op {
            EditOp::Copy(length) => position += length,
            EditOp::Substitute(sequence) => {
                for &base in sequence {
                    let code = two_bit_code(base).expect("validated two-bit substitution");
                    let delta = position - previous_mismatch;
                    let token = ((delta as u64) << 2) | code;
                    write_varint(token, writer)?;
                    previous_mismatch = position;
                    position += 1;
                }
            }
            EditOp::Insert(_) | EditOp::Delete(_) => unreachable!("validated Hamming script"),
        }
    }
    Ok(())
}

fn hamming_mismatch_count(ops: &[EditOp], target_len: usize, query_len: usize) -> Option<usize> {
    if target_len != query_len {
        return None;
    }
    let mut consumed = 0usize;
    let mut mismatches = 0usize;
    for op in ops {
        match op {
            EditOp::Copy(length) => consumed = consumed.checked_add(*length)?,
            EditOp::Substitute(sequence) if is_two_bit_dna(sequence) => {
                consumed = consumed.checked_add(sequence.len())?;
                mismatches = mismatches.checked_add(sequence.len())?;
            }
            EditOp::Substitute(_) | EditOp::Insert(_) | EditOp::Delete(_) => return None,
        }
    }
    (consumed == target_len).then_some(mismatches)
}

fn write_record(record: &EncodedRecord, writer: &mut impl Write) -> io::Result<()> {
    match record {
        EncodedRecord::Literal(sequence) => {
            writer.write_all(&[0])?;
            write_dna(sequence, writer)
        }
        EncodedRecord::Mapped {
            target_id,
            target_start,
            target_len,
            reverse,
            query_len,
            ops,
        } => {
            writer.write_all(&[RECORD_MAPPED | if *reverse { RECORD_REVERSE } else { 0 }])?;
            write_varint(*target_id as u64, writer)?;
            write_varint(*target_start as u64, writer)?;
            write_varint(*target_len as u64, writer)?;
            write_varint(*query_len as u64, writer)?;
            write_varint(ops.len() as u64, writer)?;
            for op in ops {
                match op {
                    EditOp::Copy(length) => {
                        writer.write_all(&[OP_COPY])?;
                        write_varint(*length as u64, writer)?;
                    }
                    EditOp::Substitute(sequence) => {
                        writer.write_all(&[OP_SUBSTITUTE])?;
                        write_dna(sequence, writer)?;
                    }
                    EditOp::Insert(sequence) => {
                        writer.write_all(&[OP_INSERT])?;
                        write_dna(sequence, writer)?;
                    }
                    EditOp::Delete(length) => {
                        writer.write_all(&[OP_DELETE])?;
                        write_varint(*length as u64, writer)?;
                    }
                }
            }
            Ok(())
        }
    }
}

fn read_record(reader: &mut impl Read) -> Result<EncodedRecord> {
    let mut flags = [0u8; 1];
    reader.read_exact(&mut flags)?;
    if flags[0] & RECORD_MAPPED == 0 {
        return Ok(EncodedRecord::Literal(read_dna(reader)?));
    }
    let target_id = read_usize(reader, "target ID")?;
    let target_start = read_usize(reader, "target start")?;
    if flags[0] & RECORD_HAMMING != 0 {
        let query_len = read_usize(reader, "query length")?;
        let mismatch_count = read_usize(reader, "mismatch count")?;
        ensure!(
            mismatch_count <= query_len,
            "mismatch count exceeds query length"
        );
        let mut ops = Vec::with_capacity(mismatch_count.saturating_mul(2).saturating_add(1));
        let mut previous_mismatch = 0usize;
        let mut consumed = 0usize;
        for mismatch_index in 0..mismatch_count {
            let token = read_varint(reader)?;
            let delta = usize::try_from(token >> 2).context("mismatch delta exceeds usize")?;
            let position = previous_mismatch
                .checked_add(delta)
                .context("mismatch position overflow")?;
            ensure!(
                position < query_len,
                "mismatch position exceeds query length"
            );
            ensure!(
                mismatch_index == 0 || position > previous_mismatch,
                "mismatch positions are not strictly increasing"
            );
            ensure!(position >= consumed, "mismatch positions overlap");
            push_copy(&mut ops, position - consumed);
            push_literal(&mut ops, OP_SUBSTITUTE, &[b"ACGT"[(token & 3) as usize]]);
            consumed = position + 1;
            previous_mismatch = position;
        }
        push_copy(&mut ops, query_len - consumed);
        return Ok(EncodedRecord::Mapped {
            target_id,
            target_start,
            target_len: query_len,
            reverse: flags[0] & RECORD_REVERSE != 0,
            query_len,
            ops,
        });
    }
    let target_len = read_usize(reader, "target length")?;
    let query_len = read_usize(reader, "query length")?;
    let op_count = read_usize(reader, "edit operation count")?;
    ensure!(
        op_count <= query_len.saturating_add(target_len).saturating_add(1),
        "edit operation count is inconsistent with record lengths"
    );
    let mut ops = Vec::with_capacity(op_count);
    for _ in 0..op_count {
        let mut kind = [0u8; 1];
        reader.read_exact(&mut kind)?;
        let op = match kind[0] {
            OP_COPY => EditOp::Copy(read_usize(reader, "copy length")?),
            OP_SUBSTITUTE => EditOp::Substitute(read_dna(reader)?),
            OP_INSERT => EditOp::Insert(read_dna(reader)?),
            OP_DELETE => EditOp::Delete(read_usize(reader, "deletion length")?),
            other => bail!("invalid edit operation {other}"),
        };
        ops.push(op);
    }
    Ok(EncodedRecord::Mapped {
        target_id,
        target_start,
        target_len,
        reverse: flags[0] & RECORD_REVERSE != 0,
        query_len,
        ops,
    })
}

fn read_input_list(path: &Path) -> Result<Vec<std::path::PathBuf>> {
    let reader = BufReader::new(
        File::open(path).with_context(|| format!("open input list '{}'", path.display()))?,
    );
    let parent = path.parent().unwrap_or_else(|| Path::new("."));
    let mut inputs = Vec::new();
    for (line_number, line) in reader.lines().enumerate() {
        let line = line?;
        let entry = line.trim();
        if entry.is_empty() || entry.starts_with('#') {
            continue;
        }
        let entry_path = Path::new(entry);
        let resolved = if entry_path.is_absolute() {
            entry_path.to_path_buf()
        } else {
            parent.join(entry_path)
        };
        ensure!(
            resolved.is_file(),
            "input list line {} is not a file: '{}'",
            line_number + 1,
            resolved.display()
        );
        inputs.push(resolved);
    }
    Ok(inputs)
}

fn read_fasta(path: &Path) -> Result<Vec<FastaRecord>> {
    read_fasta_impl(path, true)
}

fn read_fasta_sequences(path: &Path) -> Result<Vec<FastaRecord>> {
    read_fasta_impl(path, false)
}

fn read_fasta_impl(path: &Path, retain_names: bool) -> Result<Vec<FastaRecord>> {
    let file = File::open(path).with_context(|| format!("open FASTA '{}'", path.display()))?;
    let mut reader: Box<dyn BufRead> =
        if path.extension().is_some_and(|extension| extension == "zst") {
            Box::new(BufReader::new(zstd::Decoder::new(file).with_context(
                || format!("open Zstandard FASTA '{}'", path.display()),
            )?))
        } else {
            Box::new(BufReader::new(file))
        };
    let mut records = Vec::new();
    let mut current_name = None::<String>;
    let mut current_sequence = Vec::new();
    let mut line = Vec::new();
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        while matches!(line.last(), Some(b'\n' | b'\r')) {
            line.pop();
        }
        if let Some(header) = line.strip_prefix(b">") {
            if let Some(name) = current_name.take() {
                ensure!(
                    !current_sequence.is_empty(),
                    "FASTA record '{name}' is empty"
                );
                records.push(FastaRecord {
                    name,
                    sequence: std::mem::take(&mut current_sequence),
                });
            }
            let name = if retain_names {
                let name = header
                    .split(|byte| byte.is_ascii_whitespace())
                    .next()
                    .unwrap_or_default();
                String::from_utf8(name.to_vec()).context("FASTA header is not valid UTF-8")?
            } else {
                String::new()
            };
            current_name = Some(name);
        } else if line.iter().any(|byte| !byte.is_ascii_whitespace()) {
            ensure!(
                current_name.is_some(),
                "FASTA sequence occurs before first header"
            );
            current_sequence.extend(
                line.iter()
                    .copied()
                    .filter(|byte| !byte.is_ascii_whitespace())
                    .map(|byte| byte.to_ascii_uppercase()),
            );
        }
    }
    if let Some(name) = current_name {
        ensure!(
            !current_sequence.is_empty(),
            "FASTA record '{name}' is empty"
        );
        records.push(FastaRecord {
            name,
            sequence: current_sequence,
        });
    }
    Ok(records)
}

fn reference_fingerprint(references: &[FastaRecord]) -> blake3::Hash {
    let mut hasher = blake3::Hasher::new();
    for reference in references {
        hasher.update(&(reference.name.len() as u64).to_le_bytes());
        hasher.update(reference.name.as_bytes());
        hasher.update(&(reference.sequence.len() as u64).to_le_bytes());
        hasher.update(&reference.sequence);
    }
    hasher.finalize()
}

fn sequence_fingerprint(records: &[FastaRecord]) -> blake3::Hash {
    let mut hasher = blake3::Hasher::new();
    for record in records {
        update_sequence_fingerprint(&mut hasher, &record.sequence);
    }
    hasher.finalize()
}

fn update_sequence_fingerprint(hasher: &mut blake3::Hasher, sequence: &[u8]) {
    hasher.update(&(sequence.len() as u64).to_le_bytes());
    hasher.update(sequence);
}

fn two_bit_code(base: u8) -> Option<u64> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn encode_seed(sequence: &[u8]) -> Option<u64> {
    let mut key = 0u64;
    for &base in sequence {
        key = (key << 2) | two_bit_code(base)?;
    }
    Some(key)
}

fn complement_base(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| complement_base(*base))
        .collect()
}

fn write_bytes(bytes: &[u8], writer: &mut impl Write) -> io::Result<()> {
    write_varint(bytes.len() as u64, writer)?;
    writer.write_all(bytes)
}

fn write_dna(sequence: &[u8], writer: &mut impl Write) -> io::Result<()> {
    write_varint(sequence.len() as u64, writer)?;
    let raw = !is_two_bit_dna(sequence);
    writer.write_all(&[u8::from(raw)])?;
    write_dna_payload(sequence, raw, writer)
}

fn is_two_bit_dna(sequence: &[u8]) -> bool {
    sequence
        .iter()
        .all(|base| matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
}

fn write_dna_payload(sequence: &[u8], raw: bool, writer: &mut impl Write) -> io::Result<()> {
    if raw {
        writer.write_all(sequence)
    } else {
        let mut packed = vec![0u8; sequence.len().div_ceil(4)];
        for (index, base) in sequence.iter().enumerate() {
            let code = match base.to_ascii_uppercase() {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => unreachable!(),
            };
            packed[index / 4] |= code << (2 * (index % 4));
        }
        writer.write_all(&packed)
    }
}

fn read_dna(reader: &mut impl Read) -> Result<Vec<u8>> {
    let length = read_usize(reader, "DNA sequence length")?;
    ensure!(
        length <= 1 << 30,
        "archive DNA sequence is unreasonably large"
    );
    let mut encoding = [0u8; 1];
    reader.read_exact(&mut encoding)?;
    match encoding[0] {
        0 => read_dna_payload(length, false, reader),
        1 => read_dna_payload(length, true, reader),
        other => bail!("invalid DNA encoding {other}"),
    }
}

fn read_dna_payload(length: usize, raw: bool, reader: &mut impl Read) -> Result<Vec<u8>> {
    if raw {
        let mut sequence = vec![0u8; length];
        reader.read_exact(&mut sequence)?;
        Ok(sequence)
    } else {
        let mut packed = vec![0u8; length.div_ceil(4)];
        reader.read_exact(&mut packed)?;
        let mut sequence = Vec::with_capacity(length);
        for index in 0..length {
            let code = (packed[index / 4] >> (2 * (index % 4))) & 0b11;
            sequence.push(b"ACGT"[code as usize]);
        }
        Ok(sequence)
    }
}

fn read_bytes(reader: &mut impl Read) -> Result<Vec<u8>> {
    let length = read_usize(reader, "byte string length")?;
    ensure!(
        length <= 1 << 20,
        "archive byte string is unreasonably large"
    );
    let mut bytes = vec![0u8; length];
    reader.read_exact(&mut bytes)?;
    Ok(bytes)
}

fn read_usize(reader: &mut impl Read, field: &str) -> Result<usize> {
    usize::try_from(read_varint(reader)?).with_context(|| format!("{field} exceeds usize"))
}

/// Compact duplicate-sequence references as a tagged delta from the previous
/// logical record ID. Token zero is reserved for a newly defined sequence;
/// reference tokens are odd so that small positive and negative deltas usually
/// occupy a single byte.
fn write_sequence_reference(
    id: usize,
    previous_id: usize,
    writer: &mut impl Write,
) -> io::Result<()> {
    let id = i64::try_from(id)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "sequence ID exceeds i64"))?;
    let previous_id = i64::try_from(previous_id).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "previous sequence ID exceeds i64",
        )
    })?;
    let delta = id
        .checked_sub(previous_id)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "sequence ID delta overflow"))?;
    let zigzag = if delta >= 0 {
        (delta as u64) << 1
    } else {
        ((-1 - delta) as u64) << 1 | 1
    };
    let token = zigzag
        .checked_mul(2)
        .and_then(|value| value.checked_add(1))
        .ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "sequence reference overflow")
        })?;
    write_varint(token, writer)
}

fn read_sequence_reference(token: u64, previous_id: usize) -> Result<usize> {
    ensure!(token & 1 == 1, "invalid compact sequence reference token");
    let zigzag = token >> 1;
    let magnitude = i64::try_from(zigzag >> 1).context("sequence ID delta exceeds i64")?;
    let delta = if zigzag & 1 == 0 {
        magnitude
    } else {
        -magnitude - 1
    };
    let previous_id = i64::try_from(previous_id).context("previous sequence ID exceeds i64")?;
    let id = previous_id
        .checked_add(delta)
        .context("sequence reference ID overflow")?;
    ensure!(id >= 0, "sequence reference ID is negative");
    usize::try_from(id).context("sequence reference ID exceeds usize")
}

fn write_varint(mut value: u64, writer: &mut impl Write) -> io::Result<()> {
    loop {
        let mut byte = (value & 0x7f) as u8;
        value >>= 7;
        if value != 0 {
            byte |= 0x80;
        }
        writer.write_all(&[byte])?;
        if value == 0 {
            return Ok(());
        }
    }
}

fn read_varint(reader: &mut impl Read) -> Result<u64> {
    let mut value = 0u64;
    let mut shift = 0u32;
    loop {
        let mut byte = [0u8; 1];
        reader.read_exact(&mut byte)?;
        ensure!(shift < 64, "varint overflow");
        value |= u64::from(byte[0] & 0x7f) << shift;
        if byte[0] & 0x80 == 0 {
            return Ok(value);
        }
        shift += 7;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_zstandard_compressed_fasta() {
        let temp = tempfile::tempdir().unwrap();
        let input_path = temp.path().join("unitigs.fa.zst");
        let mut encoder = zstd::Encoder::new(File::create(&input_path).unwrap(), 1).unwrap();
        encoder.write_all(b">first\nACGT\n>second\nTTAA\n").unwrap();
        encoder.finish().unwrap();

        let records = read_fasta_sequences(&input_path).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].sequence, b"ACGT");
        assert_eq!(records[1].sequence, b"TTAA");
    }

    #[test]
    fn compact_sequence_reference_roundtrips_signed_deltas() {
        for (previous, id) in [(0, 0), (10, 11), (10, 9), (50_000, 2), (2, 50_000)] {
            let mut encoded = Vec::new();
            write_sequence_reference(id, previous, &mut encoded).unwrap();
            let token = read_varint(&mut encoded.as_slice()).unwrap();
            assert_ne!(token, 0);
            assert_eq!(read_sequence_reference(token, previous).unwrap(), id);
        }
    }

    #[test]
    fn edit_script_roundtrips_substitutions_and_indels() {
        let target = b"AACCGGTTAACCGGTT";
        let query = b"AATCGGTTGGAACGTT";
        let alignment = dp_global(
            &encode_nt4(query),
            &encode_nt4(target),
            &DpScoring::default(),
            -1,
        );
        let ops = edit_ops_from_cigar(&alignment.cigar, query, target).unwrap();
        let references = vec![FastaRecord {
            name: "ref".to_string(),
            sequence: target.to_vec(),
        }];
        let record = EncodedRecord::Mapped {
            target_id: 0,
            target_start: 0,
            target_len: target.len(),
            reverse: false,
            query_len: query.len(),
            ops,
        };
        assert_eq!(decode_record(&record, &references).unwrap(), query);
    }

    #[test]
    fn archive_roundtrip_is_exact_and_checks_reference() {
        let temp = tempfile::tempdir().unwrap();
        let reference_path = temp.path().join("reference.fa");
        let input_path = temp.path().join("unitigs.fa");
        let archive_path = temp.path().join("unitigs.kloer");
        let output_path = temp.path().join("decoded.fa");
        let mut state = 0x9e37_79b9_u32;
        let reference = (0..12_000)
            .map(|_| {
                state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                b"ACGT"[(state >> 30) as usize]
            })
            .collect::<Vec<_>>();
        std::fs::write(
            &reference_path,
            format!(">chr1\n{}\n", String::from_utf8(reference.clone()).unwrap()),
        )
        .unwrap();
        let exact = &reference[1_000..2_000];
        let mut changed = reference[3_000..4_000].to_vec();
        changed[11] = if changed[11] == b'A' { b'C' } else { b'A' };
        changed.splice(331..331, b"GGA".iter().copied());
        changed.drain(758..760);
        let reverse = reverse_complement(&reference[5_000..6_000]);
        let literal = vec![b'A'; 1_000];
        let input = format!(
            ">exact\n{}\n>changed\n{}\n>reverse\n{}\n>literal\n{}\n",
            String::from_utf8(exact.to_vec()).unwrap(),
            String::from_utf8(changed.clone()).unwrap(),
            String::from_utf8(reverse.clone()).unwrap(),
            String::from_utf8(literal.clone()).unwrap()
        );
        std::fs::write(&input_path, input).unwrap();

        let stats = compress(
            &reference_path,
            &input_path,
            &archive_path,
            &CompressOptions {
                threads: 2,
                min_coverage: 0.5,
                ..CompressOptions::default()
            },
        )
        .unwrap();
        assert_eq!(stats.records, 4);
        assert!(stats.mapped >= 3, "mapped only {} records", stats.mapped);
        let decoded = decompress(&reference_path, &archive_path, &output_path).unwrap();
        assert_eq!(decoded.records, 4);
        let decoded_records = read_fasta(&output_path).unwrap();
        assert_eq!(decoded_records[0].sequence, exact);
        assert_eq!(decoded_records[1].sequence, changed);
        assert_eq!(decoded_records[2].sequence, reverse);
        assert_eq!(decoded_records[3].sequence, literal);

        let mut altered_reference = reference.clone();
        altered_reference[0] = if altered_reference[0] == b'A' {
            b'C'
        } else {
            b'A'
        };
        std::fs::write(
            &reference_path,
            format!(">chr1\n{}\n", String::from_utf8(altered_reference).unwrap()),
        )
        .unwrap();
        assert!(decompress(&reference_path, &archive_path, &output_path).is_err());
    }

    #[test]
    fn multi_dataset_archive_roundtrip_preserves_boundaries() {
        let temp = tempfile::tempdir().unwrap();
        let reference_path = temp.path().join("reference.fa");
        let first_path = temp.path().join("first.fa");
        let second_path = temp.path().join("second.fa");
        let list_path = temp.path().join("inputs.fof");
        let archive_path = temp.path().join("cohort.kloer");
        let output_dir = temp.path().join("decoded");
        let mut state = 0x51f1_5e5d_u32;
        let reference = (0..16_000)
            .map(|_| {
                state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                b"ACGT"[(state >> 30) as usize]
            })
            .collect::<Vec<_>>();
        std::fs::write(
            &reference_path,
            format!(">chr1\n{}\n", String::from_utf8(reference.clone()).unwrap()),
        )
        .unwrap();
        let first = reference[1_000..3_000].to_vec();
        let mut second = reference[8_000..10_000].to_vec();
        second[123] = if second[123] == b'A' { b'C' } else { b'A' };
        second.splice(777..777, b"GATTACA".iter().copied());
        std::fs::write(
            &first_path,
            format!(">first\n{}\n", String::from_utf8(first.clone()).unwrap()),
        )
        .unwrap();
        std::fs::write(
            &second_path,
            format!(
                ">second\n{}\n>duplicate-first\n{}\n",
                String::from_utf8(second.clone()).unwrap(),
                String::from_utf8(first.clone()).unwrap()
            ),
        )
        .unwrap();
        std::fs::write(&list_path, "first.fa\nsecond.fa\n").unwrap();

        let stats = compress_list(
            &reference_path,
            &list_path,
            &archive_path,
            &CompressOptions {
                threads: 2,
                ..CompressOptions::default()
            },
        )
        .unwrap();
        assert_eq!(stats.records, 3);
        assert_eq!(stats.unique_records, 2);
        assert_eq!(stats.mapped, 3);
        let decoded = decompress_list(&reference_path, &archive_path, &output_dir).unwrap();
        assert_eq!(decoded.records, 3);
        assert_eq!(
            read_fasta(&output_dir.join("000000_first.fa")).unwrap()[0].sequence,
            first
        );
        let decoded_second = read_fasta(&output_dir.join("000001_second.fa")).unwrap();
        assert_eq!(decoded_second[0].sequence, second);
        assert_eq!(decoded_second[1].sequence, first);
    }
}
