use crate::{memory_fs::file::writer::FileWriter, DEFAULT_BINCODE_CONFIG};
use bincode::{Decode, Encode};
use desse::{Desse, DesseSized};
use mt_debug_counters::counter::{AtomicCounter, SumMode};
use std::io::Write;

use super::CheckpointData;

pub mod compressed_binary_writer;
pub mod lock_free_binary_writer;

pub(crate) static THREADS_BUSY_WRITING: AtomicCounter<SumMode> =
    declare_counter_i64!("threads_busy_writing", SumMode, false);

#[derive(Debug, Desse, DesseSized, Default)]
pub(crate) struct BucketHeader {
    pub magic: [u8; 16],
    pub index_offset: u64,
    pub data_format_info: [u8; Self::MAX_DATA_FORMAT_INFO_SIZE],
}

impl BucketHeader {
    pub const MAX_DATA_FORMAT_INFO_SIZE: usize = 32;
}

#[derive(Encode, Decode)]
pub(crate) struct BucketCheckpoints {
    pub index: Vec<CheckpointData>,
}

pub(crate) fn initialize_bucket_file(file: &mut FileWriter) -> u64 {
    // Write empty header
    file.write_all(&BucketHeader::default().serialize()[..])
        .unwrap();

    file.len()
}

pub(crate) fn finalize_bucket_file(
    mut file: FileWriter,
    magic: &[u8; 16],
    checkpoints: Vec<CheckpointData>,
    format_info: &[u8],
) {
    file.flush().unwrap();
    let index_position = file.len() as u64;
    bincode::encode_into_writer(
        &BucketCheckpoints { index: checkpoints },
        &mut file,
        DEFAULT_BINCODE_CONFIG,
    )
    .unwrap();

    let data_format_info = {
        let mut array = [0; 32];
        array[0..format_info.len()].copy_from_slice(format_info);
        array
    };

    file.write_at_start(
        &BucketHeader {
            magic: *magic,
            index_offset: index_position,
            data_format_info,
        }
        .serialize()[..],
    )
    .unwrap();
    file.flush_async();
}
