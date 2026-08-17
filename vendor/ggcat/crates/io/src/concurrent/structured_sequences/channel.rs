use super::{
    IdentSequenceWriter, StructuredSequenceBackend, StructuredSequenceBackendInit,
    StructuredSequenceBackendWrapper,
};
use crate::concurrent::temp_reads::extra_data::{
    SequenceExtraData, SequenceExtraDataConsecutiveCompression,
};
use config::DEFAULT_PER_CPU_BUFFER_SIZE;
use dynamic_dispatch::dynamic_dispatch;
use parking_lot::Mutex;
use std::collections::HashMap;
use std::mem;
use std::path::{Path, PathBuf};
use std::sync::LazyLock;
use std::sync::mpsc::SyncSender;

#[cfg(feature = "support_kmer_counters")]
use super::SequenceAbundance;

static CHANNEL_OUTPUTS: LazyLock<Mutex<HashMap<PathBuf, SyncSender<Vec<u8>>>>> =
    LazyLock::new(|| Mutex::new(HashMap::new()));

/// Registers a bounded output channel for a subsequent assembler invocation.
/// The path is only an opaque key; this backend never creates the file.
pub fn register_channel_output(path: PathBuf, sender: SyncSender<Vec<u8>>) -> bool {
    let mut outputs = CHANNEL_OUTPUTS.lock();
    if outputs.contains_key(&path) {
        false
    } else {
        outputs.insert(path, sender);
        true
    }
}

pub fn unregister_channel_output(path: &Path) {
    CHANNEL_OUTPUTS.lock().remove(path);
}

pub fn has_channel_output(path: &Path) -> bool {
    CHANNEL_OUTPUTS.lock().contains_key(path)
}

pub struct ChannelWriterWrapper;

#[dynamic_dispatch]
impl StructuredSequenceBackendWrapper for ChannelWriterWrapper {
    type Backend<
        ColorInfo: IdentSequenceWriter + SequenceExtraDataConsecutiveCompression,
        LinksInfo: IdentSequenceWriter + SequenceExtraData,
    > = ChannelWriter<ColorInfo, LinksInfo>;
}

pub struct ChannelBuffer {
    data: Vec<u8>,
    extra: Vec<u8>,
}

pub struct ChannelWriter<ColorInfo: IdentSequenceWriter, LinksInfo: IdentSequenceWriter> {
    sender: SyncSender<Vec<u8>>,
    path: PathBuf,
    _phantom: std::marker::PhantomData<(ColorInfo, LinksInfo)>,
}

impl<ColorInfo: IdentSequenceWriter, LinksInfo: IdentSequenceWriter> StructuredSequenceBackendInit
    for ChannelWriter<ColorInfo, LinksInfo>
{
    fn new_plain(path: impl AsRef<Path>) -> Self {
        let path = path.as_ref().to_path_buf();
        let sender = CHANNEL_OUTPUTS
            .lock()
            .remove(&path)
            .unwrap_or_else(|| panic!("no channel output registered for {}", path.display()));
        Self {
            sender,
            path,
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<ColorInfo: IdentSequenceWriter, LinksInfo: IdentSequenceWriter>
    StructuredSequenceBackend<ColorInfo, LinksInfo> for ChannelWriter<ColorInfo, LinksInfo>
{
    type SequenceTempBuffer = ChannelBuffer;

    fn alloc_temp_buffer(_: usize) -> Self::SequenceTempBuffer {
        ChannelBuffer {
            data: Vec::with_capacity(DEFAULT_PER_CPU_BUFFER_SIZE.as_bytes()),
            extra: Vec::new(),
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
        #[cfg(feature = "support_kmer_counters")] _abundance: SequenceAbundance,
    ) {
        buffer.extra.clear();
        <ColorInfo as SequenceExtraDataConsecutiveCompression>::encode_extended(
            &color_info,
            &extra_buffers.0,
            &mut buffer.extra,
            Default::default(),
        );
        let color_len = buffer.extra.len() as u64;

        let links_start = buffer.extra.len();
        <LinksInfo as SequenceExtraDataConsecutiveCompression>::encode_extended(
            &links_info,
            &extra_buffers.1,
            &mut buffer.extra,
            Default::default(),
        );
        let links_len = (buffer.extra.len() - links_start) as u64;

        buffer.data.extend_from_slice(&sequence_index.to_le_bytes());
        buffer
            .data
            .extend_from_slice(&(sequence.len() as u64).to_le_bytes());
        buffer.data.extend_from_slice(&color_len.to_le_bytes());
        buffer.data.extend_from_slice(&links_len.to_le_bytes());
        buffer.data.extend_from_slice(sequence);
        buffer.data.extend_from_slice(&buffer.extra);
    }

    fn get_path(&self) -> PathBuf {
        self.path.clone()
    }

    fn flush_temp_buffer(&mut self, buffer: &mut Self::SequenceTempBuffer) {
        if buffer.data.is_empty() {
            return;
        }
        let block = mem::replace(
            &mut buffer.data,
            Vec::with_capacity(DEFAULT_PER_CPU_BUFFER_SIZE.as_bytes()),
        );
        self.sender
            .send(block)
            .expect("KLOE structured-output receiver disconnected");
    }

    fn finalize(self) {}
}
