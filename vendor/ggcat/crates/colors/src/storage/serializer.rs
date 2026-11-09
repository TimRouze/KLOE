use crate::storage::ColorsSerializerTrait;
use config::{COLORS_SINGLE_BATCH_SIZE, ColorIndexType};
use crossbeam::channel::Sender;
use desse::{Desse, DesseSized};
use ggcat_logging::UnrecoverableErrorLogging;
use parallel_processor::execution_manager::objects_pool::{ObjectsPool, PoolObject};

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Seek, Write};
use std::marker::PhantomData;
use std::path::Path;
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::Arc;
use std::thread::JoinHandle;

pub const COLORMAP_STORAGE_VERSION: u64 = 1;

struct PreserializedColors<P: parallel_processor::execution_manager::objects_pool::PoolObjectTrait> {
    index: ColorIndexType,
    colors: PoolObject<P>,
}

#[derive(Debug, Desse, DesseSized, Default)]
pub(crate) struct ColorsFileHeader {
    pub magic: [u8; 16],
    pub version: u64,
    pub index_offset: u64,
    pub colors_count: u64,
    pub subsets_count: u64,
    pub total_size: u64,
    pub total_uncompressed_size: u64,
}

#[derive(Clone, Copy, Serialize, Deserialize, Ord, PartialOrd, Eq, PartialEq)]
pub struct ColorsIndexEntry {
    pub start_index: ColorIndexType,
    pub file_offset: u64,
}

#[derive(Serialize, Deserialize)]
pub(crate) struct ColorsIndexMap {
    pub pairs: Vec<ColorsIndexEntry>,
    pub subsets_count: u64,
}

pub struct ColorsSerializer<SI: ColorsSerializerTrait> {
    colors_subset_count: AtomicU32,
    buffers_pool: ObjectsPool<SI::PreSerializer>,
    colors_sender: Sender<PreserializedColors<SI::PreSerializer>>,

    serializer: Arc<SI>,
    checkpoint_thread: JoinHandle<SI::CheckpointTracker>,
    ordered_writer_thread: JoinHandle<()>,
    writing_threads: Vec<JoinHandle<()>>,
    print_stats: bool,
    _phantom: PhantomData<SI>,
}

impl<SI: ColorsSerializerTrait> ColorsSerializer<SI> {
    pub fn new(
        file: impl AsRef<Path>,
        color_names: &[String],
        implicit_colors_count: Option<u64>,
        threads_count: usize,
        print_stats: bool,
    ) -> anyhow::Result<Self> {
        let mut colormap_file = File::create(file.as_ref()).log_unrecoverable_error_with_data(
            "Cannot create colormap file",
            file.as_ref().display(),
        )?;

        colormap_file
            .write_all(&ColorsFileHeader::default().serialize()[..])
            .log_unrecoverable_error_with_data(
                "Cannot write colormap header",
                file.as_ref().display(),
            )?;

        colormap_file = {
            let mut color_names_stream = lz4::EncoderBuilder::new()
                .level(4)
                .build(colormap_file)
                .unwrap();
            bincode::serialize_into(&mut color_names_stream, color_names)
                .log_unrecoverable_error_with_data(
                    "Cannot serialize color names",
                    file.as_ref().display(),
                )?;

            let (cf, res) = color_names_stream.finish();
            res.log_unrecoverable_error_with_data(
                "Cannot finish color names stream",
                file.as_ref().display(),
            )?;
            cf
        };

        let file_offset = colormap_file
            .stream_position()
            .log_unrecoverable_error_with_data(
                "Cannot seek colormap file",
                file.as_ref().display(),
            )?;

        let color_processor = ColorsFlushProcessing {
            colormap_file: BufWriter::new(colormap_file),
            colormap_index: ColorsIndexMap {
                pairs: vec![],
                subsets_count: 0,
            },
            offset: file_offset,
            uncompressed_size: 0,
        };

        let colors_count = implicit_colors_count.unwrap_or(color_names.len() as u64);
        anyhow::ensure!(
            implicit_colors_count.is_none() || color_names.is_empty(),
            "implicit color count requires an empty color-name list"
        );

        let (colors_sender, receiver) =
            crossbeam::channel::bounded::<PreserializedColors<SI::PreSerializer>>(128);
        let buffers_pool = ObjectsPool::new(128, ());

        let (serializer, checkpoint_tracker) = SI::new(
            color_processor,
            COLORS_SINGLE_BATCH_SIZE as usize,
            colors_count,
        );

        let serializer = Arc::new(serializer);

        // Checkpoints are assigned round-robin to dedicated compression workers. Each
        // worker retains at most one compressed result, and the writer receives from
        // workers in assignment order. This preserves color-ID order without shared
        // condition variables, predecessor wait chains, or an unbounded reorder queue.
        let workers_count = threads_count.clamp(1, 8);
        let mut checkpoint_senders = Vec::with_capacity(workers_count);
        let mut result_receivers = Vec::with_capacity(workers_count);
        let mut recycle_senders = Vec::with_capacity(workers_count);
        let mut writing_threads = Vec::with_capacity(workers_count);
        for thread in 0..workers_count {
            let serializer = serializer.clone();
            let (checkpoint_sender, checkpoint_receiver) =
                crossbeam::channel::bounded::<SI::CheckpointWriter>(1);
            let (result_sender, result_receiver) =
                crossbeam::channel::bounded::<SI::CompressedCheckpointWriter>(1);
            let (recycle_sender, recycle_receiver) =
                crossbeam::channel::bounded::<SI::CompressedCheckpointBuffer>(1);
            checkpoint_senders.push(checkpoint_sender);
            result_receivers.push(result_receiver);
            recycle_senders.push(recycle_sender);

            writing_threads.push(
                std::thread::Builder::new()
                    .name(format!("cmap-write-{}", thread))
                    .spawn(move || {
                        let mut compressed_buffer = SI::CompressedCheckpointBuffer::default();
                        while let Ok(checkpoint) = checkpoint_receiver.recv() {
                            let compressed =
                                serializer.compress_checkpoint(checkpoint, compressed_buffer);
                            if result_sender.send(compressed).is_err() {
                                return;
                            }
                            compressed_buffer = match recycle_receiver.recv() {
                                Ok(buffer) => buffer,
                                Err(_) => return,
                            };
                        }
                    })
                    .unwrap(),
            );
        }

        let ordered_serializer = serializer.clone();
        let ordered_writer_thread = std::thread::Builder::new()
            .name("cmap-ordered-writer".to_string())
            .spawn(move || {
                let mut checkpoint_index = 0usize;
                loop {
                    let worker = checkpoint_index % workers_count;
                    let checkpoint = match result_receivers[worker].recv() {
                        Ok(checkpoint) => checkpoint,
                        Err(_) => break,
                    };
                    let recycled = ordered_serializer.commit_checkpoint(checkpoint);
                    if recycle_senders[worker].send(recycled).is_err() {
                        break;
                    }
                    checkpoint_index += 1;
                }
            })
            .unwrap();

        let checkpoint_thread = std::thread::Builder::new()
            .name("cmap-checkpoints".to_string())
            .spawn(move || {
                let mut checkpoint_tracker = checkpoint_tracker;
                let mut checkpoint_buffer = SI::CheckpointBuffer::default();
                let mut checkpoint_index = 0usize;
                let mut next_subset = 0u32;
                let mut pending = BTreeMap::new();
                while let Ok(colors) = receiver.recv() {
                    pending.insert(colors.index, colors.colors);
                    while let Some(colors) = pending.remove(&next_subset) {
                        if let Some(checkpoint) = SI::write_color_subset(
                            &mut checkpoint_tracker,
                            &mut checkpoint_buffer,
                            &colors,
                        ) {
                            let worker = checkpoint_index % workers_count;
                            checkpoint_senders[worker].send(checkpoint).unwrap();
                            checkpoint_index += 1;
                        }
                        next_subset = next_subset.checked_add(1).unwrap();
                    }
                }
                assert!(pending.is_empty(), "missing preserialized color subset");
                if let Some(checkpoint) =
                    SI::take_final_checkpoint(&mut checkpoint_tracker, checkpoint_buffer)
                {
                    let worker = checkpoint_index % workers_count;
                    checkpoint_senders[worker].send(checkpoint).unwrap();
                }
                drop(checkpoint_senders);
                checkpoint_tracker
            })
            .unwrap();

        Ok(Self {
            colors_subset_count: AtomicU32::new(0),
            buffers_pool,
            colors_sender,
            serializer,
            checkpoint_thread,
            ordered_writer_thread,
            writing_threads,
            print_stats,
            _phantom: PhantomData,
        })
    }

    #[inline(always)]
    pub fn serialize_colors(&self, colors: &[ColorIndexType]) -> ColorIndexType {
        let mut colors_buffer = self.buffers_pool.alloc_object();
        SI::preserialize_colors(&mut colors_buffer, colors);
        let new_color = self.colors_subset_count.fetch_add(1, Ordering::Relaxed);
        assert_ne!(new_color, ColorIndexType::MAX, "too many color subsets");
        self.colors_sender
            .send(PreserializedColors {
                index: new_color,
                colors: colors_buffer,
            })
            .unwrap();
        new_color
    }

    pub fn finalize(self) {
        // Drop the sender to free the writing threads
        drop(self.colors_sender);

        let tracker = self.checkpoint_thread.join().unwrap();
        self.ordered_writer_thread.join().unwrap();
        for thread in self.writing_threads {
            thread.join().unwrap();
        }

        let serializer = Arc::try_unwrap(self.serializer).unwrap_or_else(|_| unreachable!());

        if self.print_stats {
            SI::print_stats(&tracker);
        }

        serializer.finalize(tracker);
    }
}

pub struct ColorsFlushProcessing {
    pub(crate) colormap_file: BufWriter<File>,
    pub(crate) colormap_index: ColorsIndexMap,
    pub(crate) offset: u64,
    pub(crate) uncompressed_size: u64,
}

impl ColorsFlushProcessing {
    pub fn compress_chunk(data: &[u8], out_data: &mut Vec<u8>) {
        let mut encoder = lz4::EncoderBuilder::new().level(4).build(out_data).unwrap();
        encoder.write_all(data).unwrap();

        let (_, res) = encoder.finish();
        res.unwrap();
    }

    pub fn write_compressed_chunk(
        &mut self,
        start_index: ColorIndexType,
        uncompressed_size: usize,
        compressed_data: &[u8],
    ) {
        self.uncompressed_size += uncompressed_size as u64;
        let file_offset = self.offset;
        self.colormap_file.write_all(compressed_data).unwrap();
        self.offset += compressed_data.len() as u64;

        self.colormap_index.pairs.push(ColorsIndexEntry {
            start_index,
            file_offset,
        });
    }
}
