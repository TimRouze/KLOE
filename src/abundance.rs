//! Exact integer abundance support. A base of 1 is an internal mode marker,
//! accepted on disk only together with the exact-abundance format flag.
use ggcat_api::{DnaSequence, DynamicSequencesStream, SequenceInfo, SequencesReader};
use std::collections::BTreeMap;
use std::io::{self, Read};
use std::path::PathBuf;
use std::sync::Arc;

pub(crate) const EXACT_BASE: f64 = 1.0;
pub(crate) const EXACT_FLAG: u8 = 8;

pub(crate) fn width(base: Option<f64>) -> usize {
    if base == Some(EXACT_BASE) {
        8
    } else {
        1
    }
}

pub(crate) fn read_code(reader: &mut impl Read, width: usize) -> io::Result<u64> {
    let mut bytes = [0u8; 8];
    reader.read_exact(&mut bytes[..width])?;
    Ok(u64::from_le_bytes(bytes))
}

pub(crate) fn restore_code(base: u64, residual: u64, log_base: f64) -> u64 {
    let code = base.wrapping_add(residual);
    if log_base == EXACT_BASE {
        code
    } else {
        code & 255
    }
}

pub(crate) fn decode(code: u64, base: f64) -> u64 {
    if base == EXACT_BASE {
        code
    } else {
        crate::compress::decode_abundance(code as u8, base)
    }
}

/// Parse integral decimal/scientific notation without a floating-point roundtrip.
fn parse_integer(value: &str) -> Option<u64> {
    let value = value.strip_prefix('+').unwrap_or(value);
    let (mantissa, exponent) = match value.split_once(['e', 'E']) {
        Some((m, e)) => (m, e.parse::<i32>().ok()?),
        None => (value, 0),
    };
    let mut digits = String::new();
    let mut fraction = None;
    for byte in mantissa.bytes() {
        match byte {
            b'0'..=b'9' => {
                digits.push(byte as char);
                if let Some(count) = fraction.as_mut() {
                    *count += 1i64;
                }
            }
            b'.' if fraction.is_none() => fraction = Some(0i64),
            _ => return None,
        }
    }
    if digits.is_empty() {
        return None;
    }
    if digits.bytes().all(|b| b == b'0') {
        return Some(0);
    }
    let scale = i64::from(exponent) - fraction.unwrap_or(0);
    if scale < 0 {
        let remove = usize::try_from(-scale).ok()?;
        let keep = digits.len().checked_sub(remove)?;
        if !digits[keep..].bytes().all(|b| b == b'0') {
            return None;
        }
        digits.truncate(keep);
    }
    let mut result = digits.parse::<u64>().ok()?;
    if scale > 0 {
        if scale > 19 {
            return None;
        }
        for _ in 0..scale {
            result = result.checked_mul(10)?;
        }
    }
    Some(result)
}

pub(crate) fn parse_header(header: &[u8]) -> io::Result<u64> {
    if let Ok(header) = std::str::from_utf8(header) {
        for token in header.split_ascii_whitespace() {
            for prefix in ["ka:f:", "ka:i:", "km:f:", "km:i:", "abundance="] {
                if let Some(value) = token.strip_prefix(prefix) {
                    return parse_integer(value).ok_or_else(|| io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("exact abundance requires an integer in 0..=u64::MAX, got '{value}'"),
                    ));
                }
            }
        }
    }
    Err(io::Error::new(
        io::ErrorKind::InvalidInput,
        "missing supported abundance field",
    ))
}

pub(crate) type ExactSources = Arc<Vec<(u32, u64)>>;

pub(crate) struct ExactSequencesStream {
    files: Vec<PathBuf>,
    colors: Vec<BTreeMap<u64, u32>>,
}

impl ExactSequencesStream {
    pub(crate) fn prepare(files: &[String], memory_gb: usize) -> io::Result<(Self, ExactSources)> {
        let mut paths = Vec::with_capacity(files.len());
        let mut colors = Vec::with_capacity(files.len());
        let mut sources = Vec::new();
        // Conservatively budget the dictionaries and reverse mapping; reject
        // excessive cardinality instead of exceeding the requested RAM budget.
        let max_sources = memory_gb.max(1).saturating_mul(1024 * 1024 * 1024) / 8 / 96;
        for (dataset, file) in files.iter().enumerate() {
            let dataset = u32::try_from(dataset + 1)
                .map_err(|_| io::Error::other("too many exact-abundance datasets"))?;
            let path = std::fs::canonicalize(file)?;
            std::fs::File::open(&path)?;
            let mut values = BTreeMap::new();
            let mut error = None;
            SequencesReader::new().process_file_extended(&path, |sequence| {
                if error.is_some() { return; }
                match parse_header(sequence.ident_data) {
                    Ok(value) => {
                        if !values.contains_key(&value) {
                            if sources.len() >= max_sources || sources.len() >= u32::MAX as usize {
                                error = Some(io::Error::other("exact abundance dictionary exceeds its memory/color budget"));
                                return;
                            }
                            values.insert(value, sources.len() as u32);
                            sources.push((dataset, value));
                        }
                    }
                    Err(err) => error = Some(io::Error::new(err.kind(), format!("{}: {err}", path.display()))),
                }
            }, None, true, false);
            if let Some(err) = error {
                return Err(err);
            }
            paths.push(path);
            colors.push(values);
        }
        Ok((
            Self {
                files: paths,
                colors,
            },
            Arc::new(sources),
        ))
    }
}

impl DynamicSequencesStream for ExactSequencesStream {
    fn estimated_base_count(&self, block: usize) -> u64 {
        std::fs::metadata(&self.files[block]).map_or(1, |meta| meta.len().max(1))
    }
    fn read_block(
        &self,
        block: usize,
        _copy_ident_data: bool,
        partial_read_copyback: Option<usize>,
        callback: &mut dyn FnMut(DnaSequence, SequenceInfo),
    ) {
        SequencesReader::new().process_file_extended(
            &self.files[block],
            |sequence| {
                let value = parse_header(sequence.ident_data)
                    .expect("validated exact abundance input changed");
                let color = self.colors[block][&value];
                callback(sequence, SequenceInfo { color: Some(color) });
            },
            partial_read_copyback,
            true,
            false,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exact_abundance_rejects_fractional_input_before_construction() {
        let work = tempfile::tempdir().unwrap();
        let path = work.path().join("fraction.fa");
        std::fs::write(&path, ">r ka:f:17.3\nACGTACGT\n").unwrap();
        let result = ExactSequencesStream::prepare(&[path.display().to_string()], 1);
        let err = result.err().expect("fractional abundance must be rejected");
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
        assert!(err.to_string().contains("17.3"));
    }

    #[test]
    fn exact_abundance_residuals_roundtrip_at_integer_limits() {
        for base in [0u64, 1, 255, 256, 50_000, u64::MAX] {
            for value in [0u64, 1, 255, 256, 50_000, u64::MAX] {
                let residual = value.wrapping_sub(base);
                let restored = restore_code(base, residual, EXACT_BASE);
                assert_eq!(decode(restored, EXACT_BASE), value);
                let encoded = residual.to_le_bytes();
                assert_eq!(read_code(&mut encoded.as_slice(), 8).unwrap(), residual);
            }
        }
        assert_eq!(restore_code(250, 10, 1.05), 4);
    }

    #[test]
    fn exact_abundance_integer_parser() {
        for (value, expected) in [
            ("17", 17),
            ("17.000", 17),
            ("5e4", 50_000),
            ("1.23e2", 123),
            ("9007199254740993", 9_007_199_254_740_993),
            ("18446744073709551615", u64::MAX),
            ("18446744073709551615.0", u64::MAX),
        ] {
            assert_eq!(parse_integer(value), Some(expected), "{value}");
        }
        for value in [
            "17.3",
            "-1",
            "NaN",
            "inf",
            "18446744073709551616",
            "1e100",
            "1e-100",
            "",
        ] {
            assert_eq!(parse_integer(value), None, "{value}");
        }
    }
}
