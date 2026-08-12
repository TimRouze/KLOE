use crate::concurrent::structured_sequences::{IdentSequenceWriter, StructuredSequenceBackend};
use crate::concurrent::temp_reads::extra_data::{
    SequenceExtraData, SequenceExtraDataConsecutiveCompression,
};
use config::{DEFAULT_OUTPUT_BUFFER_SIZE, DEFAULT_PER_CPU_BUFFER_SIZE};
use dynamic_dispatch::dynamic_dispatch;
use flate2::Compression;
use flate2::write::GzEncoder;
use lz4::{BlockMode, BlockSize, ContentChecksum};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::marker::PhantomData;
use std::path::{Path, PathBuf};
use std::sync::{Arc, LazyLock, Mutex};

use super::stream_finish::SequencesWriterWrapper;

#[cfg(feature = "support_kmer_counters")]
use super::SequenceAbundance;
use super::{StructuredSequenceBackendInit, StructuredSequenceBackendWrapper};

pub struct FastaWriterWrapper;
#[derive(Debug)]
pub struct PlainFastaOutputRecord {
    pub sequence: Vec<u8>,
    pub color_data: Vec<u8>,
}

pub type PlainFastaOutputCallback = Arc<dyn Fn(Vec<PlainFastaOutputRecord>) + Send + Sync>;

static PLAIN_FASTA_OUTPUT_CALLBACKS: LazyLock<Mutex<HashMap<PathBuf, PlainFastaOutputCallback>>> =
    LazyLock::new(|| Mutex::new(HashMap::new()));

pub struct PlainFastaOutputCallbackGuard {
    path: PathBuf,
}

pub fn install_plain_fasta_output_callback(
    path: impl AsRef<Path>,
    callback: PlainFastaOutputCallback,
) -> io::Result<PlainFastaOutputCallbackGuard> {
    let path = path.as_ref().to_path_buf();
    let mut callbacks = PLAIN_FASTA_OUTPUT_CALLBACKS
        .lock()
        .map_err(|_| io::Error::other("plain FASTA callback lock is poisoned"))?;
    if callbacks.contains_key(&path) {
        return Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            format!(
                "a plain FASTA output callback is already installed for {}",
                path.display()
            ),
        ));
    }
    callbacks.insert(path.clone(), callback);
    Ok(PlainFastaOutputCallbackGuard { path })
}

impl Drop for PlainFastaOutputCallbackGuard {
    fn drop(&mut self) {
        if let Ok(mut callbacks) = PLAIN_FASTA_OUTPUT_CALLBACKS.lock() {
            callbacks.remove(&self.path);
        }
    }
}

#[dynamic_dispatch]
impl StructuredSequenceBackendWrapper for FastaWriterWrapper {
    type Backend<
        ColorInfo: IdentSequenceWriter + SequenceExtraDataConsecutiveCompression,
        LinksInfo: IdentSequenceWriter + SequenceExtraData,
    > = FastaWriter<ColorInfo, LinksInfo>;
}

pub struct FastaWriter<ColorInfo: IdentSequenceWriter, LinksInfo: IdentSequenceWriter> {
    writer: Box<dyn Write>,
    output_callback: Option<PlainFastaOutputCallback>,
    path: PathBuf,
    _phantom: PhantomData<(ColorInfo, LinksInfo)>,
}

pub enum FastaSequenceTempBuffer {
    Text(Vec<u8>),
    Structured(Vec<PlainFastaOutputRecord>),
}

unsafe impl<ColorInfo: IdentSequenceWriter, LinksInfo: IdentSequenceWriter> Send
    for FastaWriter<ColorInfo, LinksInfo>
{
}

unsafe impl<ColorInfo: IdentSequenceWriter, LinksInfo: IdentSequenceWriter> Sync
    for FastaWriter<ColorInfo, LinksInfo>
{
}

impl<ColorInfo: IdentSequenceWriter, LinksInfo: IdentSequenceWriter> StructuredSequenceBackendInit
    for FastaWriter<ColorInfo, LinksInfo>
{
    fn new_compressed_gzip(path: impl AsRef<Path>, level: u32) -> Self {
        let compress_stream = GzEncoder::new(
            BufWriter::with_capacity(DEFAULT_OUTPUT_BUFFER_SIZE, File::create(&path).unwrap()),
            Compression::new(level),
        );

        FastaWriter {
            writer: Box::new(SequencesWriterWrapper::new(BufWriter::with_capacity(
                DEFAULT_OUTPUT_BUFFER_SIZE,
                compress_stream,
            ))),
            path: path.as_ref().to_path_buf(),
            output_callback: None,
            _phantom: PhantomData,
        }
    }

    fn new_compressed_lz4(path: impl AsRef<Path>, level: u32) -> Self {
        let compress_stream = lz4::EncoderBuilder::new()
            .level(level)
            .checksum(ContentChecksum::NoChecksum)
            .block_mode(BlockMode::Linked)
            .block_size(BlockSize::Max1MB)
            .build(BufWriter::with_capacity(
                DEFAULT_OUTPUT_BUFFER_SIZE,
                File::create(&path).unwrap(),
            ))
            .unwrap();

        FastaWriter {
            writer: Box::new(SequencesWriterWrapper::new(BufWriter::with_capacity(
                DEFAULT_OUTPUT_BUFFER_SIZE,
                compress_stream,
            ))),
            path: path.as_ref().to_path_buf(),
            output_callback: None,
            _phantom: PhantomData,
        }
    }

    fn new_plain(path: impl AsRef<Path>) -> Self {
        let output_callback = PLAIN_FASTA_OUTPUT_CALLBACKS
            .lock()
            .expect("plain FASTA callback lock is poisoned")
            .get(path.as_ref())
            .cloned();
        let writer: Box<dyn Write> = if output_callback.is_some() {
            Box::new(io::sink())
        } else {
            Box::new(SequencesWriterWrapper::new(BufWriter::with_capacity(
                DEFAULT_OUTPUT_BUFFER_SIZE,
                File::create(&path).unwrap(),
            )))
        };
        FastaWriter {
            writer,
            path: path.as_ref().to_path_buf(),
            output_callback,
            _phantom: PhantomData,
        }
    }
}

impl<ColorInfo: IdentSequenceWriter, LinksInfo: IdentSequenceWriter>
    StructuredSequenceBackend<ColorInfo, LinksInfo> for FastaWriter<ColorInfo, LinksInfo>
{
    type SequenceTempBuffer = FastaSequenceTempBuffer;

    fn alloc_temp_buffer(&self, _: usize) -> Self::SequenceTempBuffer {
        if self.output_callback.is_some() {
            FastaSequenceTempBuffer::Structured(Vec::new())
        } else {
            FastaSequenceTempBuffer::Text(Vec::with_capacity(
                DEFAULT_PER_CPU_BUFFER_SIZE.as_bytes(),
            ))
        }
    }

    fn write_sequence(
        _k: usize,
        buffer: &mut Self::SequenceTempBuffer,
        sequence_index: u64,
        sequence: &[u8],

        color_info: ColorInfo,
        links_info: LinksInfo,
        extra_buffers: &(ColorInfo::TempBuffer, LinksInfo::TempBuffer),

        #[cfg(feature = "support_kmer_counters")] abundance: SequenceAbundance,
    ) {
        match buffer {
            FastaSequenceTempBuffer::Text(buffer) => {
                #[cfg(feature = "support_kmer_counters")]
                write!(
                    buffer,
                    ">{} LN:i:{} KC:i:{} km:f:{:.1}",
                    sequence_index,
                    sequence.len(),
                    abundance.sum,
                    abundance.sum as f64 / (sequence.len() - _k + 1) as f64
                )
                .unwrap();

                #[cfg(not(feature = "support_kmer_counters"))]
                write!(buffer, ">{} LN:i:{}", sequence_index, sequence.len(),).unwrap();

                color_info.write_as_ident(buffer, &extra_buffers.0);
                links_info.write_as_ident(buffer, &extra_buffers.1);
                buffer.extend_from_slice(b"\n");
                buffer.extend_from_slice(sequence);
                buffer.extend_from_slice(b"\n");
            }
            FastaSequenceTempBuffer::Structured(records) => {
                let mut color_data = Vec::with_capacity(color_info.max_size());
                color_info.encode_extended(
                    &extra_buffers.0,
                    &mut color_data,
                    Default::default(),
                );
                records.push(PlainFastaOutputRecord {
                    sequence: sequence.to_vec(),
                    color_data,
                });
            }
        }
    }

    fn get_path(&self) -> PathBuf {
        self.path.clone()
    }

    fn flush_temp_buffer(&mut self, buffer: &mut Self::SequenceTempBuffer) {
        match buffer {
            FastaSequenceTempBuffer::Text(buffer) => {
                self.writer.write_all(buffer).unwrap();
                buffer.clear();
            }
            FastaSequenceTempBuffer::Structured(records) => {
                if !records.is_empty() {
                    let records = std::mem::take(records);
                    self.output_callback
                        .as_ref()
                        .expect("structured FASTA buffer requires an output callback")(records);
                }
            }
        }
    }

    fn finalize(self) {}
}

impl<ColorInfo: IdentSequenceWriter, LinksInfo: IdentSequenceWriter> Drop
    for FastaWriter<ColorInfo, LinksInfo>
{
    fn drop(&mut self) {
        self.writer.flush().unwrap();
    }
}
