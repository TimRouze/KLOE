//! Small correctness regressions for the CLI, archive publication and merge metadata.
use crate::{archive, compress, decompress, merge, tests, Args};
use clap::Parser;
use std::collections::BTreeMap;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};

const DNA_A: &str = "ACGTTAGCCATGATCG";
const DNA_B: &str = "GGTACCGATCGTACCA";

fn build(files: &[PathBuf], output: &Path, k: usize, abundance: bool) -> io::Result<()> {
    let fof = output.with_extension("fof");
    fs::write(
        &fof,
        files
            .iter()
            .map(|p| format!("{}\n", p.display()))
            .collect::<String>(),
    )?;
    compress::compress_with_ggcat(
        &output.display().to_string(),
        &fof.display().to_string(),
        1,
        k,
        5,
        10,
        false,
        false,
        false,
        false,
        false,
        compress::GgcatCompressionConfig {
            memory_gb: 1,
            temp_dir: output.parent().unwrap().join("tmp").display().to_string(),
            abundance_bits: 16,
            abundance_max: abundance.then_some(50_000),
            ..Default::default()
        },
    )
}

fn restore(
    input: &Path,
    output: &Path,
    wanted: &str,
    compression: decompress::OutputCompression,
    abundance: bool,
) {
    decompress::decompress_with_options(
        &"bucket_sizes.txt".into(),
        &"id_to_color_id.txt.zst".into(),
        &"tigs_kloe.fa".into(),
        &"positions_kloe.bin".into(),
        &"filenames_id.txt".into(),
        &output.display().to_string(),
        &wanted.to_owned(),
        input.display().to_string(),
        decompress::GgcatRebuildConfig {
            enabled: false,
            restore_abundance: abundance,
            output_compression: compression,
            ..Default::default()
        },
    )
    .unwrap();
}

fn snapshot(root: &Path) -> BTreeMap<String, Vec<u8>> {
    fs::read_dir(root)
        .unwrap()
        .map(|entry| entry.unwrap())
        .filter(|entry| entry.path().is_file())
        .map(|entry| {
            (
                entry.file_name().to_string_lossy().into_owned(),
                fs::read(entry.path()).unwrap(),
            )
        })
        .collect()
}

fn merge_into(a: &Path, b: &Path, output: &Path, k: usize, structural: bool) -> io::Result<()> {
    merge::merge_archives_with_config(
        &a.display().to_string(),
        &b.display().to_string(),
        &output.display().to_string(),
        k,
        1,
        merge::MergeConfig {
            minimizer_size: 5,
            structural,
            temp_dir: output.parent().unwrap().join("tmp").display().to_string(),
            ..Default::default()
        },
    )
}

#[test]
fn cli_rejects_missing_or_unknown_modes_and_reports_real_defaults() {
    assert!(Args::try_parse_from(["kloe"]).is_err());
    assert!(Args::try_parse_from(["kloe", "unknown", "-r", "1"]).is_err());
    assert!(Args::try_parse_from(["kloe", "merge", "-k", "0"]).is_err());
    for mode in ["compress", "decompress", "merge"] {
        let args = Args::try_parse_from(["kloe", mode]).unwrap();
        assert_eq!(args.mode, mode);
        assert_eq!(args.threads, num_cpus::get());
        assert_eq!(args.minimizer_size, 12);
        assert_eq!(args.k_size, None);
        assert_eq!(args.abundance_bits, 8);
        assert_eq!(args.abundance_max, 50_000);
    }
}

#[test]
fn an_output_reservation_blocks_other_writers_and_is_released_on_drop() {
    let work = tempfile::tempdir().unwrap();
    let root = work.path().join("out");
    let lock = archive::OutputLock::acquire(&root).unwrap();
    assert!(archive::OutputLock::acquire(&root).is_err());
    #[cfg(unix)]
    {
        let alias = work.path().join("alias");
        std::os::unix::fs::symlink(&root, &alias).unwrap();
        assert!(archive::OutputLock::acquire(&alias).is_err());
    }
    drop(lock);
    assert!(!root.join(".kloe-write.lock").exists());
    assert!(archive::OutputLock::acquire(&root).is_ok());
}

#[test]
fn slashless_output_roundtrips_and_existing_archives_are_never_modified() {
    let _guard = tests::ggcat_test_lock();
    let work = tempfile::tempdir().unwrap();
    let source = work.path().join("sample.fa");
    fs::write(&source, format!(">a\n{DNA_A}\n")).unwrap();
    let output = work.path().join("archive");
    build(&[source.clone()], &output, 7, false).unwrap();
    for name in [
        "manifest.kloe",
        "filenames_id.txt",
        "filenames.idx",
        "tigs_kloe.fa",
        "tigs_kloe.idx",
        "positions_kloe.bin",
        "positions_kloe.idx",
        "bucket_sizes.txt",
        "bucket_sizes.idx",
        "dataset_to_cid.bin",
    ] {
        assert!(output.join(name).is_file(), "missing {name}");
        assert!(!work.path().join(format!("archive{name}")).exists());
    }
    let restored = work.path().join("restored");
    restore(
        &output,
        &restored,
        "",
        decompress::OutputCompression::Fasta,
        false,
    );
    tests::assert_kmer_equivalent(&[DNA_A.to_owned()], &restored.join("Dump_0_sample.fa"), 7);
    let before = snapshot(&output);
    assert_eq!(
        build(&[source], &output, 7, false).unwrap_err().kind(),
        io::ErrorKind::AlreadyExists
    );
    assert_eq!(snapshot(&output), before);
    // The Logan failure-report path must also respect the same guard, before any HTTP access.
    let error = compress::compress_with_ggcat(
        &output.display().to_string(),
        &"nonexistent-accessions.txt".into(),
        1,
        7,
        5,
        10,
        false,
        false,
        false,
        false,
        false,
        compress::GgcatCompressionConfig {
            logan: Some(crate::logan::LoganConfig::new("http://127.0.0.1:1", 1).unwrap()),
            ..Default::default()
        },
    )
    .unwrap_err();
    assert_eq!(error.kind(), io::ErrorKind::AlreadyExists);
    assert_eq!(snapshot(&output), before);
    assert!(!output.join("logan-skipped.tsv").exists());
}

#[test]
fn write_failures_never_publish_a_completion_manifest() {
    let _guard = tests::ggcat_test_lock();
    let work = tempfile::tempdir().unwrap();
    let source = work.path().join("sample.fa");
    fs::write(&source, format!(">a\n{DNA_A}\n")).unwrap();
    // Fail once during DNA writing, once while writing the final filename index.
    for blocker in [compress::TIGS_INDEX_FILE, compress::FILENAME_INDEX_FILE] {
        let output = work.path().join(blocker.replace('.', "_"));
        fs::create_dir_all(output.join(blocker)).unwrap();
        assert!(build(&[source.clone()], &output, 7, false).is_err());
        assert!(!output.join(compress::ARCHIVE_MANIFEST_FILE).exists());
        if blocker == compress::FILENAME_INDEX_FILE {
            assert!(
                output.join("filenames_id.txt").is_file(),
                "must reach final index writing"
            );
        }
    }
}

#[test]
fn stored_tigs_with_identical_basenames_have_distinct_stable_outputs() {
    let _guard = tests::ggcat_test_lock();
    let work = tempfile::tempdir().unwrap();
    let sources: Vec<_> = ["a", "b"]
        .iter()
        .map(|dir| {
            let root = work.path().join(dir);
            fs::create_dir_all(&root).unwrap();
            root.join("sample.fa")
        })
        .collect();
    fs::write(&sources[0], format!(">a ka:f:7\n{DNA_A}\n")).unwrap();
    fs::write(&sources[1], format!(">b ka:f:9\n{DNA_B}\n")).unwrap();
    let wanted = work.path().join("wanted.txt");
    fs::write(&wanted, format!("{}\n", sources[1].display())).unwrap();
    for abundance in [false, true] {
        let archive = work.path().join(format!("archive_{abundance}"));
        build(&sources, &archive, 7, abundance).unwrap();
        for (compression, suffix) in [
            (decompress::OutputCompression::Fasta, ".fa"),
            (decompress::OutputCompression::Zstd, ".fa.zst"),
        ] {
            for targeted in [false, true] {
                let out = work
                    .path()
                    .join(format!("out_{abundance}_{compression:?}_{targeted}"));
                restore(
                    &archive,
                    &out,
                    if targeted {
                        wanted.to_str().unwrap()
                    } else {
                        ""
                    },
                    compression,
                    abundance,
                );
                for (id, dna) in [DNA_A, DNA_B].iter().enumerate() {
                    let path = out.join(format!("Dump_{id}_sample{suffix}"));
                    if targeted && id == 0 {
                        assert!(!path.exists());
                    } else {
                        tests::assert_kmer_equivalent(&[dna.to_string()], &path, 7);
                    }
                }
            }
        }
    }
}

#[test]
fn merge_infers_k_and_rejects_mismatched_metadata_overwrites_and_aliases() {
    let _guard = tests::ggcat_test_lock();
    let work = tempfile::tempdir().unwrap();
    let a_source = work.path().join("a.fa");
    let b_source = work.path().join("b.fa");
    fs::write(&a_source, format!(">a\n{DNA_A}\n")).unwrap();
    fs::write(&b_source, format!(">b\n{DNA_B}\n")).unwrap();
    let a = work.path().join("a");
    let b = work.path().join("b");
    let different = work.path().join("different");
    build(&[a_source], &a, 7, false).unwrap();
    build(&[b_source.clone()], &b, 7, false).unwrap();
    build(&[b_source], &different, 9, false).unwrap();
    let merged = work.path().join("merged");
    merge_into(&a, &b, &merged, 0, false).unwrap();
    assert_eq!(archive::recorded_k(&merged).unwrap(), Some(7));
    let restored = work.path().join("restored");
    restore(
        &merged,
        &restored,
        "",
        decompress::OutputCompression::Fasta,
        false,
    );
    tests::assert_kmer_equivalent(&[DNA_A.to_owned()], &restored.join("Dump_0_a.fa"), 7);
    tests::assert_kmer_equivalent(&[DNA_B.to_owned()], &restored.join("Dump_1_b.fa"), 7);
    for (other, k) in [(&b, 31), (&different, 0), (&different, 7)] {
        let rejected = work.path().join("rejected");
        assert!(merge_into(&a, other, &rejected, k, false)
            .unwrap_err()
            .to_string()
            .contains("same k"));
        assert!(!rejected.exists());
    }
    let before = snapshot(&merged);
    assert!(merge_into(&a, &b, &merged, 0, false).is_err());
    assert_eq!(snapshot(&merged), before);
    let before_a = snapshot(&a);
    for structural in [false, true] {
        assert!(merge_into(&a, &b, &a, 0, structural).is_err());
        #[cfg(unix)]
        {
            let alias = work.path().join(format!("alias_{structural}"));
            std::os::unix::fs::symlink(&a, &alias).unwrap();
            assert!(merge_into(&a, &b, &alias, 0, structural).is_err());
        }
    }
    assert_eq!(snapshot(&a), before_a);
    fs::remove_file(a.join(compress::ARCHIVE_MANIFEST_FILE)).unwrap();
    let legacy = work.path().join("legacy_merge");
    assert!(merge_into(&a, &b, &legacy, 0, false)
        .unwrap_err()
        .to_string()
        .contains("original k"));
    assert!(!legacy.exists());
    merge_into(&a, &b, &legacy, 7, false).unwrap();
    fs::write(a.join(compress::ARCHIVE_MANIFEST_FILE), [0u8; 44]).unwrap();
    let corrupt = work.path().join("corrupt_merge");
    assert!(merge_into(&a, &b, &corrupt, 7, false).is_err());
    assert!(!corrupt.exists());
}
