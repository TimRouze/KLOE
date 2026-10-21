//! Logan HTTP inputs: embedded HTTP client -> Zstandard decoder -> bounded FASTA parser.
//! Neither the compressed input nor its decoded FASTA is written to disk.
use ggcat_api::{DnaSequence, DnaSequencesFileType, DynamicSequencesStream, SequenceInfo};
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::fs;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Condvar, Mutex};
use std::time::{Duration, Instant, SystemTime};

pub(crate) const DEFAULT_BASE_URL: &str = "https://s3.amazonaws.com/logan-pub";
const SEQUENCE_CHUNK: usize = 1024 * 1024;
const MAX_HEADER_BYTES: usize = 1024 * 1024;
const UNKNOWN_COMPRESSED_BYTES: u64 = u64::MAX;

fn retry_after(headers: &reqwest::header::HeaderMap, now: SystemTime) -> Option<Duration> {
    let value = headers
        .get(reqwest::header::RETRY_AFTER)?
        .to_str()
        .ok()?
        .trim();
    if let Ok(seconds) = value.parse::<u64>() {
        return Some(Duration::from_secs(seconds));
    }
    // HTTP-date permits IMF-fixdate and the two obsolete formats used by older mirrors.
    let date = chrono::DateTime::parse_from_rfc2822(value)
        .map(|d| d.with_timezone(&chrono::Utc))
        .ok()
        .or_else(|| {
            ["%A, %d-%b-%y %H:%M:%S GMT", "%a %b %e %H:%M:%S %Y"]
                .iter()
                .find_map(|format| {
                    chrono::NaiveDateTime::parse_from_str(value, format)
                        .ok()
                        .map(|date| date.and_utc())
                })
        })?;
    Some(
        SystemTime::from(date)
            .duration_since(now)
            .unwrap_or(Duration::ZERO),
    )
}

fn send_with_retry_hint(
    request: reqwest::blocking::RequestBuilder,
) -> Result<reqwest::blocking::Response, (reqwest::Error, Option<Duration>)> {
    let response = request.send().map_err(|error| (error, None))?;
    let delay = retry_after(response.headers(), SystemTime::now());
    response.error_for_status().map_err(|error| (error, delay))
}

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
    strict: bool,
    retries: u8,
    skipped: Arc<Mutex<BTreeMap<String, String>>>,
    threshold: Option<u64>,
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
            // Blocking Response::read applies this limit to each read operation,
            // not to the total lifetime of a large streaming download.
            .timeout(Duration::from_secs(120))
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
            strict: false,
            retries: 2,
            skipped: Arc::new(Mutex::new(BTreeMap::new())),
            threshold: None,
        })
    }

    pub(crate) fn with_policy(mut self, strict: bool, retries: u8) -> Self {
        self.strict = strict;
        self.retries = retries;
        self
    }

    pub(crate) fn with_threshold(mut self, threshold: Option<u64>) -> Self {
        self.threshold = threshold;
        self
    }

    pub(crate) fn threshold(&self) -> Option<u64> {
        self.threshold
    }

    pub(crate) fn is_skipped(&self, accession: &str) -> bool {
        self.skipped.lock().unwrap().contains_key(accession)
    }

    pub(crate) fn skip(&self, accession: &str, error: impl ToString) {
        let error = error.to_string();
        if self.strict {
            self.record_error(format!("{accession}: {error}"));
        }
        let mut skipped = self.skipped.lock().unwrap();
        if !skipped.contains_key(accession) {
            eprintln!("Warning: skipping Logan accession {accession}: {error}");
            skipped.insert(accession.to_owned(), error);
        }
    }

    pub(crate) fn write_report(&self, output: &Path) -> io::Result<()> {
        fs::create_dir_all(output)?;
        let mut file = io::BufWriter::new(fs::File::create(output.join("logan-skipped.tsv"))?);
        writeln!(file, "accession\treason")?;
        for (accession, reason) in self.skipped.lock().unwrap().iter() {
            writeln!(
                file,
                "{}\t{}",
                accession.replace(['\t', '\r', '\n'], " "),
                reason.replace(['\t', '\r', '\n'], " ")
            )?;
        }
        file.flush()
    }

    pub(crate) fn active_indices(&self, accessions: &[String]) -> Vec<usize> {
        let skipped = self.skipped.lock().unwrap();
        accessions
            .iter()
            .enumerate()
            .filter_map(|(i, a)| (!skipped.contains_key(a)).then_some(i))
            .collect()
    }

    pub(crate) fn read_accessions(&self, path: &Path) -> io::Result<Vec<String>> {
        read_accessions_impl(path, Some(self))
    }

    fn request(&self, url: &str, head: bool) -> io::Result<(reqwest::blocking::Response, u8)> {
        for attempt in 0..=self.retries {
            self.check()?;
            let request = if head {
                self.client.head(url).timeout(Duration::from_secs(60))
            } else {
                self.client.get(url)
            };
            match send_with_retry_hint(request) {
                Ok(response) => {
                    if !head && response.status() != reqwest::StatusCode::OK {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "expected a complete HTTP 200 Logan response",
                        ));
                    }
                    return Ok((response, attempt));
                }
                Err((error, delay)) => {
                    let transient = error.status().map_or(
                        error.is_connect() || error.is_timeout() || error.is_request(),
                        |status| {
                            status.is_server_error()
                                || status.as_u16() == 429
                                || status.as_u16() == 408
                        },
                    );
                    if !transient || attempt == self.retries {
                        return Err(io::Error::other(error));
                    }
                    eprintln!(
                        "Warning: retrying Logan {} request (attempt {} of {}): {error}",
                        if head { "metadata" } else { "download" },
                        u16::from(attempt) + 2,
                        u16::from(self.retries) + 1
                    );
                    self.wait_before_retry(
                        delay.unwrap_or_else(|| {
                            Duration::from_millis(250 * (1u64 << attempt.min(4)))
                        }),
                    )?;
                }
            }
        }
        unreachable!()
    }

    fn wait_before_retry(&self, delay: Duration) -> io::Result<()> {
        let deadline = Instant::now().checked_add(delay).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "Logan Retry-After delay is too large",
            )
        })?;
        let mut active = self.gate.active.lock().unwrap();
        loop {
            self.check()?;
            let remaining = deadline.saturating_duration_since(Instant::now());
            if remaining.is_zero() {
                return Ok(());
            }
            active = self
                .gate
                .available
                .wait_timeout(active, remaining)
                .unwrap()
                .0;
        }
    }

    pub(crate) fn source(&self, accession: &str) -> io::Result<LoganSource> {
        validate_accession(accession)?;
        Ok(LoganSource {
            accession: accession.to_owned(),
            url: format!("{}/u/{accession}/{accession}.unitigs.fa.zst", self.base_url),
            config: self.clone(),
            compressed_bytes: Arc::new(AtomicU64::new(UNKNOWN_COMPRESSED_BYTES)),
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
    read_accessions_impl(path, None)
}

fn read_accessions_impl(path: &Path, config: Option<&LoganConfig>) -> io::Result<Vec<String>> {
    let text = fs::read_to_string(path)?;
    let mut accessions = Vec::new();
    let mut seen = std::collections::HashSet::new();
    for (line, accession) in text.lines().enumerate() {
        let accession = accession.trim();
        if accession.is_empty() || accession.starts_with('#') {
            continue;
        }
        if let Err(error) = validate_accession(accession).map_err(|error| {
            io::Error::new(
                error.kind(),
                format!("{}:{}: {error}", path.display(), line + 1),
            )
        }) {
            if let Some(config) = config {
                config.skip(accession, &error);
                config.check()?;
                continue;
            }
            return Err(error);
        }
        if !seen.insert(accession.to_owned()) {
            if config.is_some_and(|config| !config.strict) {
                eprintln!("Warning: ignoring duplicate Logan accession {accession}");
                continue;
            }
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
    // UNKNOWN_COMPRESSED_BYTES means unknown; zero is a known empty object.
    // Never replace the file-size hint with decoded bases:
    // local input planning uses file bytes, including in abundance mode.
    compressed_bytes: Arc<AtomicU64>,
    fingerprint: Option<Arc<Mutex<Option<blake3::Hash>>>>,
}

struct FingerprintReader<R> {
    inner: R,
    hash: Option<blake3::Hasher>,
    config: LoganConfig,
}

/// Resume at the exact compressed byte offset, never replaying decoded records.
/// A strong object validator and a matching 206 Content-Range are required.
struct ResumableBody {
    response: Option<reqwest::blocking::Response>,
    config: LoganConfig,
    url: String,
    offset: u64,
    total: Option<u64>,
    validator: Option<reqwest::header::HeaderValue>,
    retries_left: u8,
}

impl ResumableBody {
    fn new(
        response: reqwest::blocking::Response,
        config: &LoganConfig,
        url: &str,
        used: u8,
    ) -> Self {
        let total = response.content_length();
        let validator = response
            .headers()
            .get(reqwest::header::ETAG)
            .filter(|value| !value.as_bytes().starts_with(b"W/"))
            .cloned();
        Self {
            response: Some(response),
            config: config.clone(),
            url: url.to_owned(),
            offset: 0,
            total,
            validator,
            retries_left: config.retries.saturating_sub(used),
        }
    }

    fn resume(&mut self, mut error: io::Error) -> io::Result<()> {
        // Close the old body before another request, keeping the concurrency
        // limit valid even when a stalled connection has not closed itself.
        self.response.take();
        let mut delay = Duration::from_millis(250);
        while self.retries_left > 0 {
            self.config.check()?;
            if self.offset > 0 && self.validator.is_none() {
                return Err(error);
            }
            self.retries_left -= 1;
            eprintln!(
                "Warning: retrying interrupted Logan download at byte {}",
                self.offset
            );
            self.config.wait_before_retry(delay)?;
            let mut request = self.config.client.get(&self.url);
            if self.offset > 0 {
                request = request
                    .header(reqwest::header::RANGE, format!("bytes={}-", self.offset))
                    .header(reqwest::header::IF_RANGE, self.validator.clone().unwrap());
            }
            match send_with_retry_hint(request) {
                Ok(response) => {
                    if self.offset > 0 {
                        let range = response
                            .headers()
                            .get(reqwest::header::CONTENT_RANGE)
                            .and_then(|value| value.to_str().ok())
                            .and_then(|value| value.strip_prefix("bytes "))
                            .and_then(|value| value.split_once('/'));
                        let valid_range = range
                            .and_then(|(span, total)| {
                                let (start, end) = span.split_once('-')?;
                                Some((
                                    start.parse::<u64>().ok()?,
                                    end.parse::<u64>().ok()?,
                                    total.parse::<u64>().ok()?,
                                ))
                            })
                            .is_some_and(|(start, end, total)| {
                                start == self.offset
                                    && start <= end
                                    && end.checked_add(1) == Some(total)
                                    && self.total.is_none_or(|expected| total == expected)
                                    && response
                                        .content_length()
                                        .is_none_or(|len| len == total - start)
                            });
                        if response.status() != reqwest::StatusCode::PARTIAL_CONTENT
                            || !valid_range
                            || response.headers().get(reqwest::header::ETAG)
                                != self.validator.as_ref()
                        {
                            return Err(io::Error::new(io::ErrorKind::InvalidData,
                                "Logan server cannot safely resume the same object; excluding the dataset"));
                        }
                        self.total = range.and_then(|(_, total)| total.parse().ok());
                    } else {
                        if response.status() != reqwest::StatusCode::OK {
                            return Err(io::Error::new(
                                io::ErrorKind::InvalidData,
                                "expected a complete HTTP 200 Logan response",
                            ));
                        }
                        self.total = response.content_length();
                        self.validator = response
                            .headers()
                            .get(reqwest::header::ETAG)
                            .filter(|v| !v.as_bytes().starts_with(b"W/"))
                            .cloned();
                    }
                    self.response = Some(response);
                    return Ok(());
                }
                Err((failure, retry_delay)) => {
                    delay = retry_delay.unwrap_or(Duration::from_millis(250));
                    let transient = failure
                        .status()
                        .is_none_or(|s| s.is_server_error() || matches!(s.as_u16(), 408 | 429));
                    error = io::Error::other(failure);
                    if !transient {
                        return Err(error);
                    }
                }
            }
        }
        Err(error)
    }
}

impl Read for ResumableBody {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        if buffer.is_empty() {
            return Ok(0);
        }
        loop {
            self.config.check()?;
            match self
                .response
                .as_mut()
                .ok_or_else(|| io::Error::other("interrupted Logan stream"))?
                .read(buffer)
            {
                Ok(0) if self.total.is_some_and(|total| self.offset < total) => {
                    self.resume(io::Error::new(
                        io::ErrorKind::UnexpectedEof,
                        "truncated Logan response",
                    ))?;
                }
                Ok(n) => {
                    self.offset += n as u64;
                    return Ok(n);
                }
                Err(error) => self.resume(error)?,
            }
        }
    }
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
        self.planning_estimate(false)
    }

    pub(crate) fn estimated_annotated_bases(&self) -> u64 {
        self.planning_estimate(true)
    }

    fn planning_estimate(&self, annotated: bool) -> u64 {
        if self.is_skipped() {
            return 0;
        }
        match self.compressed_bytes.load(Ordering::Relaxed) {
            UNKNOWN_COMPRESSED_BYTES => 1_000_000_000, // HEAD unavailable; still attempt GET.
            bytes if annotated => bytes.max(1), // Same as local abundance/filter readers.
            bytes => ggcat_api::estimated_fasta_bases_from_size(bytes, true),
        }
    }

    fn inspect_size(&self) -> io::Result<()> {
        let _permit = self.config.acquire()?;
        let (response, _) = match self.config.request(&self.url, true) {
            Ok(response) => response,
            Err(error) => {
                self.config.check()?;
                let status = error
                    .get_ref()
                    .and_then(|e| e.downcast_ref::<reqwest::Error>())
                    .and_then(reqwest::Error::status);
                if status.is_some_and(|s| matches!(s.as_u16(), 404 | 410)) {
                    return Err(error);
                }
                // HEAD is only a size hint. GET decides whether this dataset is usable.
                eprintln!(
                    "Warning: {}: metadata unavailable ({error}); trying the download",
                    self.accession
                );
                return Ok(());
            }
        };
        let bytes = response
            .headers()
            .get(reqwest::header::CONTENT_LENGTH)
            .and_then(|value| value.to_str().ok())
            .and_then(|value| value.parse::<u64>().ok());
        // A scheduling/bucketing estimate only; no body download in this request.
        self.compressed_bytes.store(bytes.unwrap_or(UNKNOWN_COMPRESSED_BYTES), Ordering::Relaxed);
        Ok(())
    }

    pub(crate) fn record_error(&self, error: impl ToString) {
        self.config.skip(&self.accession, error);
    }

    pub(crate) fn is_skipped(&self) -> bool {
        self.config.is_skipped(&self.accession)
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
        if self.is_skipped() {
            return Ok(());
        }
        let _permit = self.config.acquire()?;
        eprintln!("Streaming Logan unitigs: {}", self.accession);
        // A partially consumed body is never retried/replayed into the graph.
        let (response, retries_used) = self.config.request(&self.url, false).map_err(|error| {
            io::Error::other(format!(
                "{}: Logan download failed: {error}",
                self.accession
            ))
        })?;
        let parsed = (|| {
            let input = FingerprintReader {
                inner: ResumableBody::new(response, &self.config, &self.url, retries_used),
                hash: self.fingerprint.as_ref().map(|_| blake3::Hasher::new()),
                config: self.config.clone(),
            };
            let decoder = zstd::stream::read::Decoder::new(input)?;
            let mut reader = BufReader::with_capacity(64 * 1024, decoder);
            process_fasta_stream(
                &mut reader,
                copyback,
                copy_ident || self.config.threshold.is_some(),
                headers_only,
                &mut |sequence| {
                    if headers_only
                        || crate::abundance::passes_threshold(
                            sequence.ident_data,
                            self.config.threshold,
                        )?
                    {
                        callback(sequence)?;
                    }
                    Ok(())
                },
            )?;
            let input = reader.get_ref().get_ref().get_ref();
            Ok::<_, io::Error>((
                input.inner.offset,
                input.hash.as_ref().map(blake3::Hasher::finalize),
            ))
        })();
        let (compressed_bytes, fingerprint) = parsed.map_err(|error| {
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
        self.compressed_bytes.store(compressed_bytes, Ordering::Relaxed);
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
        pool.install(|| {
            sources.par_iter().try_for_each(|source| {
                if let Err(error) = source.inspect_size() {
                    source.record_error(error);
                }
                config.check()
            })
        })?;
        Ok(Self {
            sources,
            config: config.clone(),
        })
    }
}

impl DynamicSequencesStream for LoganSequencesStream {
    fn estimated_base_count(&self, block: usize) -> u64 {
        if self.config.threshold.is_some() {
            self.sources[block].estimated_annotated_bases()
        } else {
            self.sources[block].estimated_bases()
        }
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
            self.sources[block].record_error(error);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn logan_planning_matches_local_input_before_and_after_download() {
        use ggcat_api::GeneralSequenceBlockData;
        use std::net::TcpListener;

        for threshold in [None, Some(10)] {
            let fasta = format!(">a ka:i:20\n{}\n", "ACGT".repeat(4096));
            let body = zstd::stream::encode_all(fasta.as_bytes(), 1).unwrap();
            let dir = tempfile::tempdir().unwrap();
            let local_path = dir.path().join("SRR1.unitigs.fa.zst");
            fs::write(&local_path, &body).unwrap();
            let expected = if let Some(threshold) = threshold {
                crate::abundance::FilteredSequencesStream::prepare(
                    &[local_path.display().to_string()], threshold,
                ).unwrap().estimated_base_count(0)
            } else {
                GeneralSequenceBlockData::FASTA((local_path, None))
                    .estimated_bases_count().unwrap()
            };
            // A highly compressible fixture exposes decoded-size substitutions.
            assert!(fasta.len() > body.len() * 10);
            let listener = TcpListener::bind("127.0.0.1:0").unwrap();
            let config = LoganConfig::new(
                &format!("http://{}", listener.local_addr().unwrap()), 1,
            ).unwrap().with_threshold(threshold);
            let server = std::thread::spawn(move || {
                for method in ["HEAD", "GET"] {
                    let (mut stream, _) = listener.accept().unwrap();
                    assert!(http_request(&stream).starts_with(method));
                    write!(stream, "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n", body.len()).unwrap();
                    if method == "GET" { stream.write_all(&body).unwrap(); }
                }
            });
            let input = LoganSequencesStream::new(&["SRR1".to_owned()], &config).unwrap();
            assert_eq!(input.estimated_base_count(0), expected);
            let mut bases = 0;
            input.read_block(0, false, Some(30), &mut |sequence, _| {
                bases += sequence.seq.len();
            });
            server.join().unwrap();
            config.check().unwrap();
            assert!(!input.sources[0].is_skipped());
            assert_eq!(bases, 4 * 4096);
            assert_eq!(input.estimated_base_count(0), expected);
        }
    }

    #[test]
    fn logan_planning_excludes_failed_inputs_and_handles_unknown_sizes() {
        let config = LoganConfig::new("http://127.0.0.1:1", 1).unwrap();
        let source = config.source("SRR1").unwrap();
        assert_eq!(source.estimated_bases(), 1_000_000_000);
        assert_eq!(source.estimated_annotated_bases(), 1_000_000_000);
        for bytes in [0, 1, 101, 1024 * 1024 * 1024 * 1024] {
            source.compressed_bytes.store(bytes, Ordering::Relaxed);
            assert_eq!(source.estimated_bases(), (bytes as f64 * 0.5) as u64);
            assert_eq!(source.estimated_annotated_bases(), bytes.max(1));
        }
        source.record_error("fixture excludes an accession after a partial read");
        assert_eq!(source.estimated_bases(), 0);
        assert_eq!(source.estimated_annotated_bases(), 0);
    }

    #[test]
    fn retry_after_accepts_seconds_and_http_dates() {
        let date = chrono::DateTime::parse_from_rfc2822("Sun, 06 Nov 1994 08:49:37 GMT").unwrap();
        let now = SystemTime::from(date) - Duration::from_secs(10);
        for (value, expected) in [
            ("3", Some(3)),
            ("0", Some(0)),
            ("Sun, 06 Nov 1994 08:49:37 GMT", Some(10)),
            ("Sunday, 06-Nov-94 08:49:37 GMT", Some(10)),
            ("Sun Nov  6 08:49:37 1994", Some(10)),
            ("Sun, 06 Nov 1994 08:49:00 GMT", Some(0)),
            ("-1", None),
            ("invalid", None),
        ] {
            let mut headers = reqwest::header::HeaderMap::new();
            headers.insert(reqwest::header::RETRY_AFTER, value.parse().unwrap());
            assert_eq!(
                retry_after(&headers, now),
                expected.map(Duration::from_secs),
                "{value}"
            );
        }
    }

    #[test]
    fn strict_cancellation_wakes_a_server_requested_retry_wait() {
        let config = LoganConfig::new("http://127.0.0.1:1", 1).unwrap();
        let waiter = config.clone();
        let (tx, rx) = std::sync::mpsc::channel();
        let worker = std::thread::spawn(move || {
            tx.send(waiter.wait_before_retry(Duration::from_secs(60)))
                .unwrap();
        });
        std::thread::sleep(Duration::from_millis(20));
        config.record_error("cancelled test request");
        assert!(rx.recv_timeout(Duration::from_secs(2)).unwrap().is_err());
        worker.join().unwrap();
    }

    #[test]
    fn unsupported_or_failed_head_is_advisory_even_in_strict_mode() {
        use std::net::TcpListener;
        for status in [
            "405 Method Not Allowed",
            "501 Not Implemented",
            "503 Service Unavailable",
        ] {
            let listener = TcpListener::bind("127.0.0.1:0").unwrap();
            let config = LoganConfig::new(&format!("http://{}", listener.local_addr().unwrap()), 1)
                .unwrap()
                .with_policy(true, 0);
            let source = config.source("SRR1").unwrap();
            let server = std::thread::spawn(move || {
                let (mut stream, _) = listener.accept().unwrap();
                assert!(http_request(&stream).starts_with("HEAD "));
                write!(
                    stream,
                    "HTTP/1.1 {status}\r\nContent-Length: 0\r\nConnection: close\r\n\r\n"
                )
                .unwrap();
                drop(stream);
                let (mut stream, _) = listener.accept().unwrap();
                assert!(http_request(&stream).starts_with("GET "));
                let body = zstd::stream::encode_all(&b">a\nACGTACGT\n"[..], 1).unwrap();
                write!(
                    stream,
                    "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                    body.len()
                )
                .unwrap();
                stream.write_all(&body).unwrap();
            });
            source.inspect_size().unwrap();
            let mut sequences = Vec::new();
            source
                .read(None, true, &mut |s| sequences.push(s.seq.to_vec()))
                .unwrap();
            assert_eq!(sequences, [b"ACGTACGT".to_vec()]);
            assert!(!source.is_skipped());
            config.check().unwrap();
            server.join().unwrap();
        }
    }

    #[test]
    fn retry_after_is_honored_on_initial_and_resumed_downloads() {
        use std::net::TcpListener;
        for resume in [false, true] {
            let listener = TcpListener::bind("127.0.0.1:0").unwrap();
            let config =
                LoganConfig::new(&format!("http://{}", listener.local_addr().unwrap()), 1).unwrap();
            let source = config.source("SRR1").unwrap();
            let server = std::thread::spawn(move || {
                if resume {
                    let (mut stream, _) = listener.accept().unwrap();
                    http_request(&stream);
                    stream.write_all(b"HTTP/1.1 200 OK\r\nContent-Length: 10\r\nETag: \"same\"\r\nConnection: close\r\n\r\n0123").unwrap();
                }
                let (mut stream, _) = listener.accept().unwrap();
                let request = http_request(&stream).to_ascii_lowercase();
                if resume {
                    assert!(request.contains("range: bytes=4-"));
                }
                stream.write_all(b"HTTP/1.1 429 Too Many Requests\r\nRetry-After: 1\r\nContent-Length: 0\r\nConnection: close\r\n\r\n").unwrap();
                drop(stream);
                let (mut stream, _) = listener.accept().unwrap();
                let request = http_request(&stream).to_ascii_lowercase();
                if resume {
                    assert!(request.contains("range: bytes=4-"));
                    stream.write_all(b"HTTP/1.1 206 Partial Content\r\nContent-Length: 6\r\nContent-Range: bytes 4-9/10\r\nETag: \"same\"\r\nConnection: close\r\n\r\n456789").unwrap();
                } else {
                    stream.write_all(b"HTTP/1.1 200 OK\r\nContent-Length: 10\r\nConnection: close\r\n\r\n0123456789").unwrap();
                }
            });
            let start = Instant::now();
            let (response, used) = config.request(&source.url, false).unwrap();
            let mut body = ResumableBody::new(response, &config, &source.url, used);
            let mut bytes = Vec::new();
            body.read_to_end(&mut bytes).unwrap();
            assert_eq!(bytes, b"0123456789");
            assert!(
                start.elapsed() >= Duration::from_secs(1),
                "must honor Retry-After"
            );
            server.join().unwrap();
        }
    }

    fn http_request(stream: &std::net::TcpStream) -> String {
        stream
            .set_read_timeout(Some(Duration::from_secs(2)))
            .unwrap();
        let mut reader = BufReader::new(stream.try_clone().unwrap());
        let mut request = String::new();
        loop {
            let mut line = String::new();
            reader.read_line(&mut line).unwrap();
            if line.is_empty() || line == "\r\n" {
                break;
            }
            request.push_str(&line);
        }
        request
    }

    #[test]
    fn transient_http_failure_is_retried_but_missing_accession_is_not() {
        use std::net::TcpListener;
        for (status, attempts, succeeds) in [
            ("503 Service Unavailable", 2, true),
            ("503 Service Unavailable", 3, false),
            ("404 Not Found", 1, false),
        ] {
            let listener = TcpListener::bind("127.0.0.1:0").unwrap();
            let config =
                LoganConfig::new(&format!("http://{}", listener.local_addr().unwrap()), 1).unwrap();
            let source = config.source("SRR1").unwrap();
            let server = std::thread::spawn(move || {
                for attempt in 0..attempts {
                    let (mut stream, _) = listener.accept().unwrap();
                    let request = http_request(&stream);
                    assert!(request.starts_with("GET /u/SRR1/"));
                    if attempt == 0 || !succeeds {
                        write!(
                            stream,
                            "HTTP/1.1 {status}\r\nContent-Length: 0\r\nConnection: close\r\n\r\n"
                        )
                        .unwrap();
                    } else {
                        let body =
                            zstd::stream::encode_all(&b">r ka:f:7\nACGTACGT\n"[..], 1).unwrap();
                        write!(
                            stream,
                            "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                            body.len()
                        )
                        .unwrap();
                        stream.write_all(&body).unwrap();
                    }
                }
            });
            let mut records = 0;
            assert_eq!(
                source.read(None, true, &mut |_| records += 1).is_ok(),
                succeeds
            );
            assert_eq!(records, usize::from(succeeds));
            server.join().unwrap();
        }
    }

    #[test]
    fn interrupted_body_resumes_without_replay_and_rejects_changed_objects() {
        use std::net::TcpListener;
        for changed in [false, true] {
            let first = b">a ka:f:7\nACGTACGT\n>b ka:f:8\nTGCATGCA\n";
            let mut body = zstd::stream::encode_all(&first[..], 1).unwrap();
            let cut = body.len();
            body.extend(zstd::stream::encode_all(&b">c ka:f:9\nGATTACAG\n"[..], 1).unwrap());
            let listener = TcpListener::bind("127.0.0.1:0").unwrap();
            let config =
                LoganConfig::new(&format!("http://{}", listener.local_addr().unwrap()), 1).unwrap();
            let source = config.source("SRR1").unwrap();
            let server = std::thread::spawn(move || {
                let (mut stream, _) = listener.accept().unwrap();
                http_request(&stream);
                write!(stream, "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nETag: \"object-a\"\r\nConnection: close\r\n\r\n", body.len()).unwrap();
                stream.write_all(&body[..cut]).unwrap();
                drop(stream);
                let (mut stream, _) = listener.accept().unwrap();
                let request = http_request(&stream).to_ascii_lowercase();
                assert!(
                    request.contains(&format!("range: bytes={cut}-")),
                    "{request}"
                );
                assert!(request.contains("if-range: \"object-a\""), "{request}");
                let etag = if changed { "object-b" } else { "object-a" };
                write!(stream, "HTTP/1.1 206 Partial Content\r\nContent-Length: {}\r\nContent-Range: bytes {}-{}/{}\r\nETag: \"{etag}\"\r\nConnection: close\r\n\r\n",
                    body.len() - cut, cut, body.len() - 1, body.len()).unwrap();
                let _ = stream.write_all(&body[cut..]);
            });
            let mut records = Vec::new();
            let result = source.read(None, true, &mut |s| records.push(s.seq.to_vec()));
            if changed {
                assert!(result.unwrap_err().to_string().contains("safely resume"));
                assert!(
                    !records.is_empty(),
                    "the failure must occur after a record was emitted"
                );
            } else {
                result.unwrap();
                assert_eq!(
                    records,
                    [
                        b"ACGTACGT".to_vec(),
                        b"TGCATGCA".to_vec(),
                        b"GATTACAG".to_vec()
                    ]
                );
            }
            server.join().unwrap();
        }
    }

    #[test]
    fn download_limit_is_respected_and_one_failure_does_not_cancel_queued_downloads() {
        use std::net::TcpListener;
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        listener.set_nonblocking(true).unwrap();
        let config = LoganConfig::new(&format!("http://{}", listener.local_addr().unwrap()), 2)
            .unwrap()
            .with_policy(false, 0);
        let (done_tx, done_rx) = std::sync::mpsc::channel();
        let handles: Vec<_> = (1..=6)
            .map(|i| {
                let config = config.clone();
                let done = done_tx.clone();
                std::thread::spawn(move || {
                    let source = config.source(&format!("SRR{i}")).unwrap();
                    let result = source.read(None, false, &mut |_| {});
                    if let Err(error) = &result {
                        source.record_error(error);
                    }
                    done.send(result.is_ok()).unwrap();
                })
            })
            .collect();
        let deadline = std::time::Instant::now() + Duration::from_secs(5);
        for _ in 0..3 {
            let mut batch = Vec::new();
            while batch.len() < 2 {
                match listener.accept() {
                    Ok((stream, _)) => {
                        let request = http_request(&stream);
                        batch.push((stream, request));
                    }
                    Err(error) if error.kind() == io::ErrorKind::WouldBlock => {
                        assert!(
                            std::time::Instant::now() < deadline,
                            "queued downloads did not continue"
                        );
                        std::thread::sleep(Duration::from_millis(5));
                    }
                    Err(error) => panic!("{error}"),
                }
            }
            assert_eq!(*config.gate.active.lock().unwrap(), 2);
            assert_eq!(
                listener.accept().unwrap_err().kind(),
                io::ErrorKind::WouldBlock
            );
            for (mut stream, request) in batch {
                if request.starts_with("GET /u/SRR1/") {
                    write!(
                        stream,
                        "HTTP/1.1 404 Not Found\r\nContent-Length: 0\r\nConnection: close\r\n\r\n"
                    )
                    .unwrap();
                } else {
                    let body = zstd::stream::encode_all(&b">r ka:f:7\nACGTACGT\n"[..], 1).unwrap();
                    write!(
                        stream,
                        "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                        body.len()
                    )
                    .unwrap();
                    stream.write_all(&body).unwrap();
                }
            }
        }
        let successes = (0..6)
            .filter(|_| done_rx.recv_timeout(Duration::from_secs(2)).unwrap())
            .count();
        for handle in handles {
            handle.join().unwrap();
        }
        assert_eq!(successes, 5);
        assert!(config.is_skipped("SRR1"));
        config.check().unwrap();
    }

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
