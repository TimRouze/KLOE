//! Logan HTTP inputs: embedded HTTP client -> Zstandard decoder -> bounded FASTA parser.
//! Neither the compressed input nor its decoded FASTA is written to disk.
use ggcat_api::{DnaSequence, DnaSequencesFileType, DynamicSequencesStream, SequenceInfo};
use rayon::prelude::*;
use std::fs;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Condvar, Mutex};
use std::time::Duration;

pub(crate) const DEFAULT_BASE_URL: &str = "https://s3.amazonaws.com/logan-pub";
const SEQUENCE_CHUNK: usize = 1024 * 1024;
const MAX_HEADER_BYTES: usize = 1024 * 1024;

#[derive(Debug)]
struct DownloadGate {
    active: Mutex<usize>,
    available: Condvar,
    limit: usize,
}

struct DownloadPermit<'a>(&'a DownloadGate);

impl Drop for DownloadPermit<'_> {
    fn drop(&mut self) {
        *self.0.active.lock().unwrap() -= 1;
        self.0.available.notify_one();
    }
}

#[derive(Clone, Debug)]
pub(crate) struct LoganConfig {
    base_url: String,
    gate: Arc<DownloadGate>,
    failure: Arc<Mutex<Option<String>>>,
    failed: Arc<AtomicBool>,
    client: reqwest::blocking::Client,
}

impl LoganConfig {
    pub(crate) fn new(base_url: &str, downloads: usize) -> io::Result<Self> {
        if downloads == 0
            || !(base_url.starts_with("https://") || base_url.starts_with("http://"))
            || base_url
                .chars()
                .any(|c| c.is_whitespace() || c == '?' || c == '#')
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Logan requires a HTTP(S) base URL and a positive download limit",
            ));
        }
        let client = reqwest::blocking::Client::builder()
            .connect_timeout(Duration::from_secs(30))
            // The blocking client's default timeout is too short for large datasets.
            .timeout(None)
            .retry(reqwest::retry::never())
            .redirect(reqwest::redirect::Policy::custom(|attempt| {
                if !matches!(attempt.url().scheme(), "http" | "https") {
                    attempt.error("Logan redirects must use HTTP(S)")
                } else if attempt.previous().len() >= 10 {
                    attempt.error("too many Logan redirects")
                } else {
                    attempt.follow()
                }
            }))
            .no_gzip()
            .no_brotli()
            .no_deflate()
            .no_zstd()
            .user_agent(concat!("kloe/", env!("CARGO_PKG_VERSION")))
            .build()
            .map_err(io::Error::other)?;
        Ok(Self {
            base_url: base_url.trim_end_matches('/').to_owned(),
            gate: Arc::new(DownloadGate {
                active: Mutex::new(0),
                available: Condvar::new(),
                limit: downloads,
            }),
            failure: Arc::new(Mutex::new(None)),
            failed: Arc::new(AtomicBool::new(false)),
            client,
        })
    }

    pub(crate) fn source(&self, accession: &str) -> io::Result<LoganSource> {
        validate_accession(accession)?;
        Ok(LoganSource {
            accession: accession.to_owned(),
            url: format!("{}/u/{accession}/{accession}.unitigs.fa.zst", self.base_url),
            config: self.clone(),
            estimated_bases: Arc::new(AtomicU64::new(1)),
            fingerprint: Some(Arc::new(Mutex::new(None))),
        })
    }

    pub(crate) fn check(&self) -> io::Result<()> {
        if !self.failed.load(Ordering::Acquire) {
            return Ok(());
        }
        match self.failure.lock().unwrap().as_ref() {
            Some(error) => Err(io::Error::other(error.clone())),
            None => Ok(()),
        }
    }

    pub(crate) fn record_error(&self, error: impl ToString) {
        // Use the same lock order as acquire, so a waiter cannot miss cancellation.
        let _active = self.gate.active.lock().unwrap();
        self.failure
            .lock()
            .unwrap()
            .get_or_insert_with(|| error.to_string());
        self.failed.store(true, Ordering::Release);
        self.gate.available.notify_all();
    }

    pub(crate) fn parallel_pool(&self, count: usize) -> io::Result<rayon::ThreadPool> {
        // GGCAT owns the global Rayon pool; keep HTTP/abundance preparation local.
        rayon::ThreadPoolBuilder::new()
            .num_threads(self.gate.limit.min(count).max(1))
            .build()
            .map_err(io::Error::other)
    }

    fn acquire(&self) -> io::Result<DownloadPermit<'_>> {
        let mut active = self.gate.active.lock().unwrap();
        while *active >= self.gate.limit {
            self.check()?;
            active = self.gate.available.wait(active).unwrap();
        }
        self.check()?;
        *active += 1;
        Ok(DownloadPermit(&self.gate))
    }
}

fn validate_accession(accession: &str) -> io::Result<()> {
    let valid = accession.len() > 3
        && ["SRR", "ERR", "DRR"]
            .iter()
            .any(|prefix| accession.starts_with(prefix))
        && accession.as_bytes()[3..].iter().all(u8::is_ascii_digit);
    if !valid {
        return Err(io::Error::new(io::ErrorKind::InvalidInput,
            format!("invalid Logan run accession '{accession}'; expected SRR, ERR, or DRR followed by digits")));
    }
    Ok(())
}

pub(crate) fn read_accessions(path: &Path) -> io::Result<Vec<String>> {
    let text = fs::read_to_string(path)?;
    let mut accessions = Vec::new();
    let mut seen = std::collections::HashSet::new();
    for (line, accession) in text.lines().enumerate() {
        let accession = accession.trim();
        if accession.is_empty() || accession.starts_with('#') {
            continue;
        }
        validate_accession(accession).map_err(|error| {
            io::Error::new(
                error.kind(),
                format!("{}:{}: {error}", path.display(), line + 1),
            )
        })?;
        if !seen.insert(accession.to_owned()) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "{}:{}: duplicate accession {accession}",
                    path.display(),
                    line + 1
                ),
            ));
        }
        accessions.push(accession.to_owned());
    }
    if accessions.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Logan accession list is empty",
        ));
    }
    Ok(accessions)
}

#[derive(Clone)]
pub(crate) struct LoganSource {
    accession: String,
    url: String,
    config: LoganConfig,
    estimated_bases: Arc<AtomicU64>,
    fingerprint: Option<Arc<Mutex<Option<blake3::Hash>>>>,
}

struct FingerprintReader<R> {
    inner: R,
    hash: Option<blake3::Hasher>,
    config: LoganConfig,
}

impl<R: Read> Read for FingerprintReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.config.check()?;
        let n = self.inner.read(buf)?;
        if let Some(hash) = &mut self.hash {
            hash.update(&buf[..n]);
        }
        Ok(n)
    }
}

impl LoganSource {
    pub(crate) fn estimated_bases(&self) -> u64 {
        self.estimated_bases.load(Ordering::Relaxed)
    }

    fn inspect_size(&self) -> io::Result<()> {
        let _permit = self.config.acquire()?;
        let response = self
            .config
            .client
            .head(&self.url)
            .timeout(Duration::from_secs(60))
            .send()
            .and_then(reqwest::blocking::Response::error_for_status)
            .map_err(|error| {
                io::Error::other(format!(
                    "{}: Logan metadata request failed: {error}",
                    self.accession
                ))
            })?;
        let bytes = response
            .headers()
            .get(reqwest::header::CONTENT_LENGTH)
            .and_then(|value| value.to_str().ok())
            .and_then(|value| value.parse::<u64>().ok());
        // A scheduling/bucketing estimate only; no body download in this request.
        self.estimated_bases.store(
            bytes.map_or(1_000_000_000, |n| n.saturating_mul(4).max(1)),
            Ordering::Relaxed,
        );
        Ok(())
    }

    pub(crate) fn record_error(&self, error: impl ToString) {
        self.config
            .record_error(format!("{}: {}", self.accession, error.to_string()));
    }

    pub(crate) fn read(
        &self,
        copyback: Option<usize>,
        copy_ident: bool,
        callback: &mut dyn FnMut(DnaSequence),
    ) -> io::Result<()> {
        self.read_impl(copyback, copy_ident, false, &mut |sequence| {
            callback(sequence);
            Ok(())
        })
    }

    pub(crate) fn read_headers(
        &self,
        callback: &mut dyn FnMut(DnaSequence) -> io::Result<()>,
    ) -> io::Result<()> {
        self.read_impl(None, true, true, callback)
    }

    fn read_impl(
        &self,
        copyback: Option<usize>,
        copy_ident: bool,
        headers_only: bool,
        callback: &mut dyn FnMut(DnaSequence) -> io::Result<()>,
    ) -> io::Result<()> {
        let _permit = self.config.acquire()?;
        eprintln!("Streaming Logan unitigs: {}", self.accession);
        // A partially consumed body is never retried/replayed into the graph.
        let response = self
            .config
            .client
            .get(&self.url)
            .send()
            .and_then(reqwest::blocking::Response::error_for_status)
            .map_err(|error| {
                io::Error::other(format!(
                    "{}: Logan download failed: {error}",
                    self.accession
                ))
            })?;
        let parsed = (|| {
            let input = FingerprintReader {
                inner: response,
                hash: self.fingerprint.as_ref().map(|_| blake3::Hasher::new()),
                config: self.config.clone(),
            };
            let decoder = zstd::stream::read::Decoder::new(input)?;
            let mut reader = BufReader::with_capacity(64 * 1024, decoder);
            let bases =
                process_fasta_stream(&mut reader, copyback, copy_ident, headers_only, callback)?;
            Ok::<_, io::Error>((
                bases,
                reader
                    .get_ref()
                    .get_ref()
                    .get_ref()
                    .hash
                    .as_ref()
                    .map(blake3::Hasher::finalize),
            ))
        })();
        let (bases, fingerprint) = parsed.map_err(|error| {
            io::Error::new(
                error.kind(),
                format!("{}: Logan download failed: {error}", self.accession),
            )
        })?;
        if let Some(shared) = &self.fingerprint {
            let fingerprint = fingerprint.expect("fingerprinting enabled");
            let mut previous = shared.lock().unwrap();
            if previous.is_some_and(|old| old != fingerprint) {
                return Err(io::Error::other(format!(
                    "{}: Logan object changed between streaming passes",
                    self.accession
                )));
            }
            *previous = Some(fingerprint);
        }
        self.estimated_bases.store(bases.max(1), Ordering::Relaxed);
        Ok(())
    }
}

// Parse whole line fragments from fixed-size buffers. The header pass counts
// bases without copying/normalizing DNA; graph input retains only one chunk.
#[cfg(test)]
fn read_fasta_stream(
    reader: impl BufRead,
    copyback: Option<usize>,
    copy_ident: bool,
    callback: &mut dyn FnMut(DnaSequence),
) -> io::Result<u64> {
    process_fasta_stream(reader, copyback, copy_ident, false, &mut |sequence| {
        callback(sequence);
        Ok(())
    })
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

fn process_fasta_stream(
    mut reader: impl BufRead,
    copyback: Option<usize>,
    copy_ident: bool,
    headers_only: bool,
    callback: &mut dyn FnMut(DnaSequence) -> io::Result<()>,
) -> io::Result<u64> {
    let overlap = copyback.unwrap_or(0);
    let chunk_size = SEQUENCE_CHUNK.max(overlap.saturating_mul(2));
    let mut header = Vec::with_capacity(if copy_ident { 256 } else { 0 });
    let mut sequence = Vec::with_capacity(if headers_only { 0 } else { chunk_size });
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
                            if sequence.len() == chunk_size {
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
    if bases == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Logan stream contains no FASTA sequences",
        ));
    }
    Ok(bases)
}

pub(crate) struct LoganSequencesStream {
    sources: Vec<LoganSource>,
    config: LoganConfig,
}

impl LoganSequencesStream {
    pub(crate) fn new(accessions: &[String], config: &LoganConfig) -> io::Result<Self> {
        let sources: Vec<_> = accessions
            .iter()
            .map(|a| {
                config.source(a).map(|mut source| {
                    source.fingerprint = None;
                    source
                })
            })
            .collect::<io::Result<_>>()?;
        // GGCAT initializes the global Rayon pool itself. Keep HTTP work local.
        let pool = config.parallel_pool(sources.len())?;
        pool.install(|| sources.par_iter().try_for_each(LoganSource::inspect_size))?;
        Ok(Self {
            sources,
            config: config.clone(),
        })
    }
}

impl DynamicSequencesStream for LoganSequencesStream {
    fn estimated_base_count(&self, block: usize) -> u64 {
        self.sources[block].estimated_bases()
    }
    fn read_block(
        &self,
        block: usize,
        copy_ident: bool,
        copyback: Option<usize>,
        callback: &mut dyn FnMut(DnaSequence, SequenceInfo),
    ) {
        if let Err(error) = self.sources[block].read(copyback, copy_ident, &mut |sequence| {
            callback(sequence, SequenceInfo { color: None })
        }) {
            self.config.record_error(error);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_pass_fingerprints_detect_sequence_changes_with_identical_headers() {
        use std::io::Write;
        use std::net::TcpListener;
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let config =
            LoganConfig::new(&format!("http://{}", listener.local_addr().unwrap()), 1).unwrap();
        let source = config.source("SRR1").unwrap();
        let server = std::thread::spawn(move || {
            for fasta in [">a ka:f:7\nACGTACGT\n", ">a ka:f:7\nATGTACGT\n"] {
                let (mut stream, _) = listener.accept().unwrap();
                stream
                    .set_read_timeout(Some(Duration::from_secs(2)))
                    .unwrap();
                let mut reader = BufReader::new(stream.try_clone().unwrap());
                loop {
                    let mut line = String::new();
                    reader.read_line(&mut line).unwrap();
                    if line == "\r\n" || line.is_empty() {
                        break;
                    }
                }
                let body = zstd::stream::encode_all(fasta.as_bytes(), 1).unwrap();
                write!(
                    stream,
                    "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                    body.len()
                )
                .unwrap();
                stream.write_all(&body).unwrap();
            }
        });
        source.read_headers(&mut |_| Ok(())).unwrap();
        let error = source.read(None, true, &mut |_| {}).unwrap_err();
        assert!(error
            .to_string()
            .contains("changed between streaming passes"));
        server.join().unwrap();
    }

    #[test]
    fn header_scan_matches_full_parser_across_buffer_boundaries() {
        let input = b";intro\r\n>a ka:f:7\r\nacgt\r\n;skip this\nNn\r\n>empty\n>b km:f:12\nTgC\rA";
        for capacity in 1..=input.len() {
            let mut full = Vec::new();
            let bases = read_fasta_stream(
                BufReader::with_capacity(capacity, &input[..]),
                Some(2),
                true,
                &mut |s| full.push((s.ident_data.to_vec(), s.seq.to_vec())),
            )
            .unwrap();
            assert_eq!(bases, 10);
            assert_eq!(
                full,
                vec![
                    (b">a ka:f:7".to_vec(), b"ACGTNN".to_vec()),
                    (b">b km:f:12".to_vec(), b"TGCA".to_vec())
                ],
                "buffer={capacity}"
            );
            let mut headers = Vec::new();
            let scanned = process_fasta_stream(
                BufReader::with_capacity(capacity, &input[..]),
                None,
                true,
                true,
                &mut |s| {
                    assert!(s.seq.is_empty());
                    headers.push(s.ident_data.to_vec());
                    Ok(())
                },
            )
            .unwrap();
            assert_eq!(scanned, bases);
            assert_eq!(headers, full.into_iter().map(|s| s.0).collect::<Vec<_>>());
        }
    }

    #[test]
    fn header_scan_propagates_callback_errors_immediately() {
        let input = b">a ka:f:7\nACGT\n>b ka:f:8\nACGT\n";
        let mut reader = BufReader::with_capacity(4, &input[..]);
        let mut calls = 0;
        let error = process_fasta_stream(&mut reader, None, true, true, &mut |_| {
            calls += 1;
            Err(io::Error::other("dictionary full"))
        })
        .unwrap_err();
        assert_eq!(calls, 1);
        assert_eq!(error.to_string(), "dictionary full");
        assert!(!reader.buffer().is_empty() || !reader.get_ref().is_empty());
    }

    #[test]
    fn cancelling_a_download_wakes_every_queued_reader() {
        let config = LoganConfig::new("http://127.0.0.1", 1).unwrap();
        let _permit = config.acquire().unwrap();
        let (started_tx, started_rx) = std::sync::mpsc::channel();
        let (done_tx, done_rx) = std::sync::mpsc::channel();
        let handles: Vec<_> = (0..8)
            .map(|_| {
                let config = config.clone();
                let started = started_tx.clone();
                let done = done_tx.clone();
                std::thread::spawn(move || {
                    started.send(()).unwrap();
                    done.send(config.acquire().is_err()).unwrap();
                })
            })
            .collect();
        for _ in 0..8 {
            started_rx.recv_timeout(Duration::from_secs(2)).unwrap();
        }
        config.record_error("cancelled download");
        // Cancellation must wake everyone even while the occupied slot stays held.
        for _ in 0..8 {
            assert!(done_rx.recv_timeout(Duration::from_secs(2)).unwrap());
        }
        for handle in handles {
            handle.join().unwrap();
        }
    }

    #[test]
    fn accession_list_rejects_invalid_duplicate_and_empty_inputs() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("accessions.txt");
        fs::write(&path, "# inputs\n SRR1 \n\nERR2\nDRR3\n").unwrap();
        assert_eq!(read_accessions(&path).unwrap(), ["SRR1", "ERR2", "DRR3"]);
        for text in ["SRR1\nSRR1", "../SRR1", "SRR", "GCA_123", "# no inputs"] {
            fs::write(&path, text).unwrap();
            assert!(read_accessions(&path).is_err(), "{text}");
        }
    }

    #[test]
    fn fasta_parser_preserves_wrapped_records_and_headers() {
        let mut records = Vec::new();
        read_fasta_stream(
            &b">a ka:f:7\r\nacgt\n;comment\nNN\n>b km:f:12\nTGC\nA"[..],
            Some(2),
            true,
            &mut |s| records.push((s.ident_data.to_vec(), s.seq.to_vec())),
        )
        .unwrap();
        assert_eq!(
            records,
            vec![
                (b">a ka:f:7".to_vec(), b"ACGTNN".to_vec()),
                (b">b km:f:12".to_vec(), b"TGCA".to_vec())
            ]
        );
        assert!(read_fasta_stream(&b"not FASTA"[..], None, false, &mut |_| {}).is_err());
    }

    #[test]
    fn long_lines_are_split_without_losing_or_inventing_kmers() {
        let input = format!(">long ka:f:10\n{}", "ACGT".repeat(SEQUENCE_CHUNK / 4 + 4));
        let mut chunks = Vec::new();
        read_fasta_stream(input.as_bytes(), Some(30), true, &mut |s| {
            chunks.push(s.seq.to_vec())
        })
        .unwrap();
        assert!(chunks.iter().all(|s| s.len() <= SEQUENCE_CHUNK));
        let mut restored = chunks[0].clone();
        for chunk in &chunks[1..] {
            restored.extend_from_slice(&chunk[30..]);
        }
        assert_eq!(
            restored,
            input.as_bytes().split(|&b| b == b'\n').nth(1).unwrap()
        );
    }
}
