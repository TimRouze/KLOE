//! Exact integer abundance support. A base of 1 is an internal mode marker,
//! accepted on disk only together with the exact-abundance format flag.
use ggcat_api::{DnaSequence, DynamicSequencesStream, SequenceInfo, SequencesReader};
use std::collections::BTreeMap;
use std::io::{self, Read};
use std::path::PathBuf;
use std::sync::Arc;

pub(crate) const EXACT_BASE: f64 = 1.0;
pub(crate) const EXACT_FLAG: u8 = 8;
pub(crate) const WIDTH8_FLAG: u8 = 16;
pub(crate) const WIDTH16_FLAG: u8 = 32;
pub(crate) const WIDTH32_FLAG: u8 = 64;
pub(crate) const WIDTH_FLAGS: u8 = WIDTH8_FLAG | WIDTH16_FLAG | WIDTH32_FLAG;
pub(crate) const DENSE_FLAG: u8 = 128;

/// The code width is explicit on disk; legacy quantized archives remain bytes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Encoding {
    pub base: f64,
    pub bytes: usize,
    /// Densely number logarithmic representatives; retain the same decoded values.
    pub dense: bool,
    cutoff_code: u64,
    cutoff_value: u64,
}

impl Encoding {
    pub(crate) const EXACT: Self = Self {
        base: EXACT_BASE,
        bytes: 8,
        dense: false,
        cutoff_code: 0,
        cutoff_value: 0,
    };

    pub(crate) const U16: Self = Self {
        base: EXACT_BASE,
        bytes: 2,
        dense: false,
        cutoff_code: 0,
        cutoff_value: 0,
    };

    pub(crate) fn quantized(base: f64, bits: u8) -> Self {
        Self {
            base,
            bytes: usize::from(bits).div_ceil(8).next_power_of_two(),
            dense: false,
            cutoff_code: 0,
            cutoff_value: 0,
        }
    }

    pub(crate) fn flags(self) -> u8 {
        if self.base == EXACT_BASE {
            return EXACT_FLAG | match self.bytes {
                1 => WIDTH8_FLAG,
                2 => WIDTH16_FLAG,
                4 => WIDTH32_FLAG,
                8 => 0,
                _ => unreachable!(),
            };
        }
        let width_flags = match self.bytes {
            1 => 0,
            2 => WIDTH16_FLAG,
            4 => WIDTH32_FLAG,
            8 => EXACT_FLAG,
            _ => unreachable!(),
        };
        width_flags | if self.dense { DENSE_FLAG } else { 0 }
    }

    pub(crate) fn from_flags(base: f64, flags: u8) -> io::Result<Self> {
        let kind = flags & (EXACT_FLAG | WIDTH_FLAGS);
        let dense = flags & DENSE_FLAG != 0;
        if dense && kind & EXACT_FLAG != 0 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "conflicting abundance representation flags"));
        }
        let bytes = match kind {
            0 => 1,
            WIDTH16_FLAG => 2,
            WIDTH32_FLAG => 4,
            EXACT_FLAG if base == EXACT_BASE => 8,
            kind if kind == (EXACT_FLAG | WIDTH8_FLAG) && base == EXACT_BASE => 1,
            kind if kind == (EXACT_FLAG | WIDTH16_FLAG) && base == EXACT_BASE => 2,
            kind if kind == (EXACT_FLAG | WIDTH32_FLAG) && base == EXACT_BASE => 4,
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid abundance encoding flags",
                ))
            }
        };
        if kind & EXACT_FLAG == 0 && (!base.is_finite() || base <= 1.0) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid abundance log base",
            ));
        }
        let (cutoff_code, cutoff_value) = if dense { Self::dense_cutoffs(base) } else { (0, 0) };
        Ok(Self { base, bytes, dense, cutoff_code, cutoff_value })
    }

    pub(crate) fn mask(self) -> u64 {
        u64::MAX >> (64 - 8 * self.bytes)
    }

    /// Select lossless integers whenever the requested logical code space fits
    /// the configured ceiling; otherwise span 1..maximum logarithmically.
    pub(crate) fn automatic(bits: u8, maximum: u64) -> io::Result<Self> {
        if !(2..=32).contains(&bits) || maximum == 0 {
            return Err(io::Error::new(io::ErrorKind::InvalidInput,
                "abundance requires 2..=32 bits and a positive maximum"));
        }
        let base = if maximum <= (1u64 << bits) - 1 {
            EXACT_BASE
        } else if maximum == 50_000 {
            crate::compress::default_abundance_log_base(bits)
        } else {
            (maximum as f64).powf(1.0 / ((1u64 << bits) - 2) as f64)
        };
        let mut encoding = Self::quantized(base, bits);
        if base != EXACT_BASE {
            let (code, value) = Self::dense_cutoffs(base);
            if value < code {
                encoding.dense = true;
                encoding.cutoff_code = code;
                encoding.cutoff_value = value;
            }
        }
        Ok(encoding)
    }

    fn dense_cutoffs(base: f64) -> (u64, u64) {
        if base >= 2.0 { return (1, 1); }
        // Before this code, successive real representatives differ by < 1:
        // their rounded values cover consecutive integers with repetitions.
        // From this code onward the gap is >= 1, so every code is distinct.
        let mut code = ((1.0 / (base - 1.0)).ln() / base.ln()).ceil() as u64 + 1;
        while base.powf((code - 1) as f64) * (base - 1.0) < 1.0 { code += 1; }
        while code > 1 && base.powf((code - 2) as f64) * (base - 1.0) >= 1.0 { code -= 1; }
        (code, crate::compress::decode_abundance(code, base))
    }

    pub(crate) fn pack_log_code(self, code: u64) -> u64 {
        if !self.dense { return code; }
        if code <= self.cutoff_code {
            crate::compress::decode_abundance(code, self.base)
        } else {
            self.cutoff_value + (code - self.cutoff_code)
        }
    }
}

pub(crate) fn width(base: Option<Encoding>) -> usize {
    base.map_or(1, |encoding| encoding.bytes)
}

pub(crate) fn read_code(reader: &mut impl Read, width: usize) -> io::Result<u64> {
    let mut bytes = [0u8; 8];
    reader.read_exact(&mut bytes[..width])?;
    Ok(u64::from_le_bytes(bytes))
}

pub(crate) fn restore_code(base: u64, residual: u64, encoding: Encoding) -> u64 {
    let residual = if encoding.dense {
        (residual >> 1) ^ 0u64.wrapping_sub(residual & 1)
    } else { residual };
    base.wrapping_add(residual) & encoding.mask()
}

pub(crate) fn residual_code(base: u64, code: u64, encoding: Encoding) -> u64 {
    let residual = code.wrapping_sub(base) & encoding.mask();
    if encoding.dense {
        let sign = residual >> (encoding.bytes * 8 - 1);
        ((residual << 1) ^ 0u64.wrapping_sub(sign)) & encoding.mask()
    } else { residual }
}

pub(crate) fn decode(code: u64, encoding: Encoding) -> u64 {
    if encoding.base == EXACT_BASE {
        code
    } else if encoding.dense {
        if code <= encoding.cutoff_value { code }
        else { crate::compress::decode_abundance(encoding.cutoff_code + (code - encoding.cutoff_value), encoding.base) }
    } else {
        crate::compress::decode_abundance(code, encoding.base)
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
    quantization: Option<(f64, u8)>,
    maximum: Option<u64>,
    storage: Option<Encoding>,
}

impl ExactSequencesStream {
    pub(crate) fn prepare(files: &[String], memory_gb: usize) -> io::Result<(Self, ExactSources)> {
        Self::prepare_impl(files, memory_gb, None, None)
    }

    pub(crate) fn prepare_quantized(
        files: &[String],
        memory_gb: usize,
        base: f64,
        bits: u8,
    ) -> io::Result<(Self, ExactSources)> {
        if !base.is_finite() || base <= 1.0 {
            return Err(io::Error::new(io::ErrorKind::InvalidInput, "logarithmic abundance base must be finite and > 1"));
        }
        Encoding::from_flags(base, Encoding::quantized(base, bits).flags())?;
        Self::prepare_impl(files, memory_gb, Some((base, bits)), None)
    }

    /// Round fractional input as in quantized mode, then store the integer directly.
    /// Overflow is rejected rather than silently clipped.
    pub(crate) fn prepare_u16(files: &[String], memory_gb: usize) -> io::Result<(Self, ExactSources)> {
        Self::prepare_impl(files, memory_gb, Some((EXACT_BASE, 16)), None)
    }

    pub(crate) fn prepare_automatic_values(
        files: &[String], memory_gb: usize, bits: u8, maximum: u64,
    ) -> io::Result<(Self, ExactSources)> {
        let encoding = Encoding::automatic(bits, maximum)?;
        if encoding.base != EXACT_BASE && !encoding.dense {
            return Err(io::Error::new(io::ErrorKind::InvalidInput, "representatives do not fit the code container"));
        }
        Self::prepare_impl(files, memory_gb, Some((encoding.base, bits)), Some(maximum))
    }

    fn code_with_max(header: &[u8], quantization: Option<(f64, u8)>, maximum: Option<u64>, storage: Option<Encoding>) -> io::Result<u64> {
        if let Some(maximum) = maximum {
            // Keep the logarithmic approximation but remove unused/duplicate
            // integer representatives from its code numbering.
            return crate::compress::parse_header_abundance(header)
                .map(|value| {
                    let value = value.min(maximum);
                    match quantization {
                        Some((base, bits)) if base != EXACT_BASE => {
                            let code = crate::compress::encode_abundance_with_bits(value, base, bits);
                            storage.expect("automatic encoding prepared above").pack_log_code(code)
                        }
                        _ => value,
                    }
                })
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing supported abundance field"));
        }
        Self::code(header, quantization)
    }

    fn code(header: &[u8], quantization: Option<(f64, u8)>) -> io::Result<u64> {
        match quantization {
            None => parse_header(header),
            Some((base, 16)) if base == EXACT_BASE => crate::compress::parse_header_abundance(header)
                .filter(|&value| value <= u16::MAX as u64)
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput,
                    "u16 abundance requires a supported nonnegative value that rounds into 0..=65535")),
            Some((base, bits)) => crate::compress::parse_header_abundance(header)
                .map(|value| crate::compress::encode_abundance_with_bits(value, base, bits))
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "missing supported abundance field",
                    )
                }),
        }
    }

    fn prepare_impl(
        files: &[String],
        memory_gb: usize,
        quantization: Option<(f64, u8)>,
        maximum: Option<u64>,
    ) -> io::Result<(Self, ExactSources)> {
        let storage = quantization.zip(maximum)
            .map(|((_, bits), maximum)| Encoding::automatic(bits, maximum)).transpose()?;
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
                match Self::code_with_max(sequence.ident_data, quantization, maximum, storage) {
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
                quantization,
                maximum,
                storage,
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
                let value = Self::code_with_max(sequence.ident_data, self.quantization, self.maximum, self.storage)
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
    fn dense_log_codes_preserve_reconstruction_and_bit_budget() {
        for bits in 2..=15 {
            let encoding = Encoding::automatic(bits, 50_000).unwrap();
            assert_eq!(Encoding::from_flags(encoding.base, encoding.flags()).unwrap(), encoding);
            let mut previous = 0;
            for code in 0..=(1u64 << bits)-1 {
                let packed = encoding.pack_log_code(code);
                assert!(packed <= code && packed >= previous);
                assert_eq!(decode(packed, encoding), crate::compress::decode_abundance(code, encoding.base),
                    "bits={bits}, code={code}");
                previous = packed;
            }
            for value in 0..=50_000 {
                let code = crate::compress::encode_abundance_with_bits(value, encoding.base, bits);
                let packed = encoding.pack_log_code(code);
                assert_eq!(decode(packed, encoding), crate::compress::decode_abundance(code, encoding.base));
            }
        }
        for (bits, maximum) in [(8, 256), (10, 5000), (24, 50_000_000), (32, u64::MAX)] {
            let encoding = Encoding::automatic(bits, maximum).unwrap();
            assert_eq!(Encoding::from_flags(encoding.base, encoding.flags()).unwrap(), encoding);
            for value in [0, 1, 2, 3, maximum/2, maximum] {
                let code = crate::compress::encode_abundance_with_bits(value, encoding.base, bits);
                let packed = encoding.pack_log_code(code);
                assert!(packed <= code);
                assert_eq!(decode(packed, encoding), crate::compress::decode_abundance(code, encoding.base));
            }
        }
        assert!(Encoding::from_flags(EXACT_BASE, EXACT_FLAG | DENSE_FLAG).is_err());
    }

    #[test]
    fn dense_signed_differences_roundtrip_at_each_storage_width() {
        for (bits, maximum) in [(8, 50_000), (14, 50_000), (24, 50_000_000)] {
            let encoding = Encoding::automatic(bits, maximum).unwrap();
            assert!(encoding.dense);
            let mask = encoding.mask();
            for base in [0, 1, mask/2, mask] {
                for code in [0, 1, 2, mask/2, mask/2+1, mask-1, mask] {
                    let residual = residual_code(base, code, encoding);
                    assert!(residual <= mask);
                    assert_eq!(restore_code(base, residual, encoding), code);
                }
            }
            assert_eq!(residual_code(3, 2, encoding), 1);
            assert_eq!(residual_code(2, 3, encoding), 2);
        }
    }

    #[test]
    fn automatic_encoding_switches_at_integer_capacity() {
        for bits in 2..=32 {
            let capacity = (1u64 << bits) - 1;
            let exact = Encoding::automatic(bits, capacity).unwrap();
            assert_eq!(exact.base, EXACT_BASE);
            assert_eq!(Encoding::from_flags(exact.base, exact.flags()).unwrap(), exact);
            for value in [0, 1, capacity / 2, capacity] {
                let header = format!(">r ka:f:{value}");
                let code = ExactSequencesStream::code_with_max(header.as_bytes(),
                    Some((EXACT_BASE, bits)), Some(capacity), None).unwrap();
                let residual = code.wrapping_sub(capacity).to_le_bytes();
                let restored = restore_code(capacity, read_code(&mut &residual[..exact.bytes], exact.bytes).unwrap(), exact);
                assert_eq!(decode(restored, exact), value);
            }
            assert!(Encoding::automatic(bits, capacity + 1).unwrap().base > 1.0);
        }
        assert!(Encoding::automatic(15, 50_000).unwrap().base > 1.0);
        assert_eq!(Encoding::automatic(16, 50_000).unwrap(), Encoding::U16);
        assert_eq!(Encoding::automatic(8, 50_000).unwrap().base,
            crate::compress::DEFAULT_ABUNDANCE_LOG_BASE);
        assert_eq!(Encoding::automatic(4, 500).unwrap().base, 500f64.powf(1.0 / 14.0));
        assert_eq!(ExactSequencesStream::code_with_max(b">r ka:f:50001",
            Some((EXACT_BASE, 16)), Some(50_000), None).unwrap(), 50_000);
        assert!(Encoding::automatic(8, 0).is_err());
        assert!(Encoding::automatic(33, 50_000).is_err());
    }

    #[test]
    fn u16_abundance_preserves_all_integers_and_rejects_overflow() {
        let encoding = Encoding::U16;
        assert_eq!(Encoding::from_flags(encoding.base, encoding.flags()).unwrap(), encoding);
        for value in 0..=u16::MAX as u64 {
            let header = format!(">r ka:f:{value}");
            let code = ExactSequencesStream::code(header.as_bytes(), Some((EXACT_BASE, 16))).unwrap();
            assert_eq!(code, value);
            let residual = code.wrapping_sub(50_000).to_le_bytes();
            let restored = restore_code(50_000, read_code(&mut &residual[..2], 2).unwrap(), encoding);
            assert_eq!(decode(restored, encoding), value);
        }
        for (text, expected) in [("0.4", 0), ("6142.75", 6143), ("65535.49", 65535)] {
            let header = format!(">r ka:f:{text}");
            assert_eq!(ExactSequencesStream::code(header.as_bytes(), Some((EXACT_BASE, 16))).unwrap(), expected);
        }
        for text in ["65535.5", "65536", "100000", "NaN", "-1"] {
            let header = format!(">r ka:f:{text}");
            assert!(ExactSequencesStream::code(header.as_bytes(), Some((EXACT_BASE, 16))).is_err());
        }
        assert!(Encoding::from_flags(EXACT_BASE, WIDTH16_FLAG).is_err());
        assert!(Encoding::from_flags(EXACT_BASE, EXACT_FLAG | WIDTH16_FLAG | WIDTH32_FLAG).is_err());
    }

    #[test]
    fn wide_abundance_codes_wrap_at_their_declared_width() {
        for bits in [2, 4, 8, 16, 32] {
            let encoding =
                Encoding::quantized(crate::compress::default_abundance_log_base(bits), bits);
            assert_eq!(
                Encoding::from_flags(encoding.base, encoding.flags()).unwrap(),
                encoding
            );
            let limit = encoding.mask();
            for base in [0, 1, limit / 2, limit] {
                for code in [0, 1, limit / 2, limit] {
                    let bytes = code.wrapping_sub(base).to_le_bytes();
                    let residual =
                        read_code(&mut &bytes[..encoding.bytes], encoding.bytes).unwrap();
                    assert_eq!(restore_code(base, residual, encoding), code);
                }
            }
        }
        for flags in [
            WIDTH16_FLAG | WIDTH32_FLAG,
            EXACT_FLAG | WIDTH16_FLAG,
            EXACT_FLAG | WIDTH32_FLAG,
        ] {
            assert!(Encoding::from_flags(1.05, flags).is_err());
        }
    }
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
                let restored = restore_code(base, residual, Encoding::EXACT);
                assert_eq!(decode(restored, Encoding::EXACT), value);
                let encoded = residual.to_le_bytes();
                assert_eq!(read_code(&mut encoded.as_slice(), 8).unwrap(), residual);
            }
        }
        assert_eq!(restore_code(250, 10, Encoding::quantized(1.05, 8)), 4);
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
