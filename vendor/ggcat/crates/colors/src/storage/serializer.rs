use crate::storage::ColorsSerializerTrait;
use config::{COLORS_SINGLE_BATCH_SIZE, ColorIndexType};
use crossbeam::channel::Sender;
use desse::{Desse, DesseSized};
use ggcat_logging::UnrecoverableErrorLogging;
use parallel_processor::execution_manager::objects_pool::{ObjectsPool, PoolObject};

use parking_lot::Mutex;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufWriter, Seek, Write};
use std::marker::PhantomData;
use std::path::Path;
use std::sync::Arc;
use std::thread::JoinHandle;

pub const COLORMAP_STORAGE_VERSION: u64 = 1;

struct OrderedCheckpoint<C> {
    checkpoint: C,
    wait_for_previous: crossbeam::channel::Receiver<()>,
    release_next: crossbeam::channel::Sender<()>,
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
    // This lock intentionally covers both ID assignment and channel enqueue: concurrent
    // producers must enqueue subsets in exactly the same order as their assigned IDs.
    colors_subset_count: Mutex<ColorIndexType>,
    buffers_pool: ObjectsPool<SI::PreSerializer>,
    colors_sender: Sender<PoolObject<SI::PreSerializer>>,

    serializer: Arc<SI>,
    checkpoint_thread: JoinHandle<SI::CheckpointTracker>,
    writing_threads: Vec<JoinHandle<()>>,
    print_stats: bool,
    _phantom: PhantomData<SI>,
}

impl<SI: ColorsSerializerTrait> ColorsSerializer<SI> {
    pub fn new(
        file: impl AsRef<Path>,
        color_names: &[String],
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

        let colors_count = color_names.len() as u64;

        let (colors_sender, receiver) =
            crossbeam::channel::bounded::<PoolObject<SI::PreSerializer>>(128);
        let buffers_pool = ObjectsPool::new(128, ());

        let (serializer, checkpoint_tracker) = SI::new(
            color_processor,
            COLORS_SINGLE_BATCH_SIZE as usize,
            colors_count,
        );

        let serializer = Arc::new(serializer);

        // Subsets must enter checkpoints in color-ID order. A single inexpensive
        // checkpoint builder provides that ordering, while completed checkpoints are
        // compressed concurrently by the writer pool below.
        let (checkpoint_sender, checkpoint_receiver) =
            crossbeam::channel::bounded::<OrderedCheckpoint<SI::CheckpointWriter>>(
                (threads_count.max(1) * 2).max(2),
            );

        let mut writing_threads = Vec::with_capacity(threads_count.max(1));
        for thread in 0..threads_count.max(1) {
            let serializer = serializer.clone();
            let checkpoint_receiver = checkpoint_receiver.clone();

            writing_threads.push(
                std::thread::Builder::new()
                    .name(format!("cmap-write-{}", thread))
                    .spawn(move || {
                        let mut compressed_checkpoint_buffer =
                            SI::CompressedCheckpointBuffer::default();
                        while let Ok(checkpoint) = checkpoint_receiver.recv() {
                            serializer.flush_checkpoint(
                                checkpoint.checkpoint,
                                &mut compressed_checkpoint_buffer,
                                checkpoint.wait_for_previous,
                                checkpoint.release_next,
                            );
                        }
                    })
                    .unwrap(),
            );
        }

        let checkpoint_thread = std::thread::Builder::new()
            .name("cmap-checkpoints".to_string())
            .spawn(move || {
                let mut checkpoint_tracker = checkpoint_tracker;
                let mut checkpoint_buffer = SI::CheckpointBuffer::default();
                let (first_release, first_wait) = crossbeam::channel::bounded(1);
                first_release.send(()).unwrap();
                let mut wait_for_previous = Some(first_wait);
                while let Ok(colors) = receiver.recv() {
                    if let Some(checkpoint) = SI::write_color_subset(
                        &mut checkpoint_tracker,
                        &mut checkpoint_buffer,
                        &colors,
                    ) {
                        let (release_next, next_wait) = crossbeam::channel::bounded(1);
                        checkpoint_sender
                            .send(OrderedCheckpoint {
                                checkpoint,
                                wait_for_previous: wait_for_previous.take().unwrap(),
                                release_next,
                            })
                            .unwrap();
                        wait_for_previous = Some(next_wait);
                    }
                }
                if let Some(checkpoint) =
                    SI::take_final_checkpoint(&mut checkpoint_tracker, checkpoint_buffer)
                {
                    let (release_next, _next_wait) = crossbeam::channel::bounded(1);
                    checkpoint_sender
                        .send(OrderedCheckpoint {
                            checkpoint,
                            wait_for_previous: wait_for_previous.take().unwrap(),
                            release_next,
                        })
                        .unwrap();
                }
                drop(checkpoint_sender);
                checkpoint_tracker
            })
            .unwrap();

        Ok(Self {
            colors_subset_count: Mutex::new(0),
            buffers_pool,
            colors_sender,
            serializer,
            checkpoint_thread,
            writing_threads,
            print_stats,
            _phantom: PhantomData,
        })
    }

    #[inline(always)]
    pub fn serialize_colors(&self, colors: &[ColorIndexType]) -> ColorIndexType {
        let mut colors_buffer = self.buffers_pool.alloc_object();
        SI::preserialize_colors(&mut colors_buffer, colors);
        let mut colors_subset_count = self.colors_subset_count.lock();
        let new_color = *colors_subset_count;
        *colors_subset_count += 1;
        self.colors_sender.send(colors_buffer).unwrap();
        new_color
    }

    pub fn finalize(self) {
        // Drop the sender to free the writing threads
        drop(self.colors_sender);

        let tracker = self.checkpoint_thread.join().unwrap();
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
