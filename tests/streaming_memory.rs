//! Reproduce the 8192-bucket startup without downloading a large collection.
#[path = "../vendor/parallel-processor/src/memory_budget.rs"]
mod memory_budget;

use parallel_processor::buckets::writers::compressed_binary_writer::{
    CompressedBinaryWriter, CompressedCheckpointSize, CompressionLevelInfo,
};
use parallel_processor::buckets::{
    readers::binary_reader::{BinaryChunkReader, ChunkedBinaryReaderIndex},
    readers::compressed_decoder::CompressedStreamDecoder,
    LockFreeBucket,
};
use parallel_processor::memory_data_size::MemoryDataSize;
use parallel_processor::memory_fs::file::internal::MemoryFileMode;
use parallel_processor::memory_fs::{MemoryFs, RemoveFileMode};
use std::io::Read;

fn payload(round: u8) -> Vec<u8> {
    let mut state = 0x123456789abcdefu64 ^ u64::from(round);
    (0..9000)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state as u8
        })
        .collect()
}

#[test]
fn effective_budget_respects_slurm_allocation() {
    if std::env::var_os("KLOE_MEMORY_BUDGET_TEST_CHILD").is_some() {
        let budget = memory_budget::effective_memory_budget(700 * 1024 * 1024 * 1024);
        assert!(budget > 0 && budget <= (1024 * 1024 * 1024 / 5 * 4));
        return;
    }
    let result = std::process::Command::new(std::env::current_exe().unwrap())
        .args(["--exact", "effective_budget_respects_slurm_allocation"])
        .env("KLOE_MEMORY_BUDGET_TEST_CHILD", "1")
        .env("SLURM_MEM_PER_NODE", "1024")
        .env_remove("SLURM_MEM_PER_CPU")
        .status()
        .unwrap();
    assert!(result.success());
}

#[test]
fn thousands_of_buckets_round_trip_with_eviction() {
    let dir = tempfile::tempdir().unwrap();
    eprintln!("Bucket test directory: {}", dir.path().display());
    MemoryFs::init(MemoryDataSize::from_mebioctets(64), 8, 2, 8192, None);
    let buckets: Vec<_> = (0..8192)
        .map(|i| {
            CompressedBinaryWriter::new(
                &dir.path().join("bucket"),
                &(
                    MemoryFileMode::DiskOnly,
                    CompressedCheckpointSize::new_from_log2(20),
                    CompressionLevelInfo {
                        fast_disk: 0,
                        slow_disk: 0,
                    },
                ),
                i,
                &0u8,
            )
        })
        .collect();
    // Leave one bucket empty; revisit all the others after their compressor
    // has been evicted, and cross physical buffer boundaries.
    for round in 0..3u8 {
        let data = payload(round);
        std::thread::scope(|scope| {
            for part in buckets[..8191].chunks(1024) {
                let data = &data;
                scope.spawn(move || {
                    for bucket in part {
                        bucket.write_data(data);
                    }
                });
            }
        });
    }
    // Tagged checkpoints must not be split by cache eviction. Their metadata
    // can contain the sequence count for the entire logical block.
    buckets[0].set_checkpoint_data(Some(&42u32), None);
    buckets[0].write_data(b"tagged-");
    for bucket in &buckets[1..512] {
        bucket.write_data(b"extra");
    }
    buckets[0].write_data(b"checkpoint");
    let paths: Vec<_> = buckets.iter().map(|b| b.get_path()).collect();
    for bucket in buckets {
        bucket.finalize();
    }
    MemoryFs::flush_to_disk(true);

    // Copy compressed frames without decoding, as bucket compaction does.
    // Checkpoints must still delimit each original frame after eviction.
    let copied = CompressedBinaryWriter::new(
        &dir.path().join("copied"),
        &(
            MemoryFileMode::DiskOnly,
            CompressedCheckpointSize::new_from_log2(20),
            CompressionLevelInfo {
                fast_disk: 0,
                slow_disk: 0,
            },
        ),
        0,
        &0u8,
    );
    let source =
        parallel_processor::memory_fs::file::reader::FileReader::open(&paths[2], None).unwrap();
    let source_index = ChunkedBinaryReaderIndex::from_file(&paths[2], RemoveFileMode::Keep, None);
    // The on-disk bucket header is magic[16], index_offset[8], format_info[32].
    let mut offset = 56;
    let mut expected_copied = Vec::new();
    for chunk in source_index.into_chunks() {
        let length = chunk.get_length();
        copied.set_checkpoint_data(
            Some(&7u32),
            Some(source.get_range_reference(offset..offset + length)),
        );
        BinaryChunkReader::<CompressedStreamDecoder>::new(chunk)
            .read_to_end(&mut expected_copied)
            .unwrap();
        offset += length;
    }
    drop(source);
    copied.write_data(b"tail");
    expected_copied.extend_from_slice(b"tail");
    let copied_path = copied.get_path();
    copied.finalize();
    MemoryFs::flush_to_disk(true);
    let mut actual_copied = Vec::new();
    for chunk in
        ChunkedBinaryReaderIndex::from_file(&copied_path, RemoveFileMode::Keep, None).into_chunks()
    {
        assert_eq!(chunk.get_extra_data::<u32>(), Some(7));
        BinaryChunkReader::<CompressedStreamDecoder>::new(chunk)
            .read_to_end(&mut actual_copied)
            .unwrap();
    }
    assert_eq!(actual_copied, expected_copied);

    for (i, path) in paths.iter().enumerate() {
        let index = ChunkedBinaryReaderIndex::from_file(path, RemoveFileMode::Keep, None);
        let mut decoded = Vec::new();
        let mut tagged = 0;
        for chunk in index.into_chunks() {
            if chunk.get_extra_data::<u32>() == Some(42) {
                tagged += 1;
            }
            BinaryChunkReader::<CompressedStreamDecoder>::new(chunk)
                .read_to_end(&mut decoded)
                .unwrap();
        }
        let mut expected = Vec::new();
        if i != 8191 {
            for round in 0..3u8 {
                expected.extend(payload(round));
            }
        }
        if i == 0 {
            expected.extend_from_slice(b"tagged-checkpoint");
        }
        if (1..512).contains(&i) {
            expected.extend_from_slice(b"extra");
        }
        assert_eq!(decoded, expected, "bucket {i}");
        assert_eq!(tagged, usize::from(i == 0));
    }
    assert_eq!(
        parallel_processor::memory_fs::allocator::CHUNKS_ALLOCATOR
            .get_total_memory()
            .as_bytes(),
        64 * 1024 * 1024,
        "Disk bucket buffers must not trigger emergency allocator growth"
    );
    MemoryFs::remove_directory(dir.path(), true);
    MemoryFs::terminate();
}
