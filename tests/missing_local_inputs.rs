use std::collections::BTreeSet;
use std::fs;
use std::path::Path;
use std::process::{Command, Output};

fn run(root: &Path, args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_kloe"))
        .current_dir(root)
        .args(args)
        .output()
        .unwrap()
}

fn assert_success(output: &Output) {
    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn dna(mut seed: u64) -> String {
    (0..160)
        .map(|_| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            b"ACGT"[(seed & 3) as usize] as char
        })
        .collect()
}

fn kmers(sequences: &[String]) -> BTreeSet<Vec<u8>> {
    sequences
        .iter()
        .flat_map(|s| {
            s.as_bytes().windows(17).map(|kmer| {
                let reverse: Vec<_> = kmer
                    .iter()
                    .rev()
                    .map(|b| match b {
                        b'A' => b'T',
                        b'C' => b'G',
                        b'G' => b'C',
                        b'T' => b'A',
                        _ => panic!(),
                    })
                    .collect();
                kmer.to_vec().min(reverse)
            })
        })
        .collect()
}

fn sequences(fasta: &str) -> Vec<String> {
    fasta
        .split('>')
        .skip(1)
        .map(|record| record.lines().skip(1).collect())
        .collect()
}

#[test]
fn missing_paths_are_excluded_before_ids_and_abundance_processing() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let first = dna(12345);
    let second = dna(67890);
    let low = dna(98765);
    fs::write(
        root.join("first.fa"),
        format!(">first ka:i:20\n{first}\n>low ka:i:3\n{low}\n"),
    )
    .unwrap();
    fs::write(
        root.join("second.fa.zst"),
        zstd::encode_all(format!(">second ka:i:40\n{second}\n").as_bytes(), 1).unwrap(),
    )
    .unwrap();
    fs::write(
        root.join("inputs.fof"),
        "\nmissing-before.fa\nfirst.fa\nmissing-middle.fa\n\nsecond.fa.zst\nmissing-after.fa\n",
    )
    .unwrap();
    fs::write(root.join("selected.txt"), "second.fa.zst\n").unwrap();

    for (mode, flags, abundance, filtered) in [
        ("plain", vec![], false, false),
        ("filtered", vec!["--abundance-threshold", "3"], false, true),
        (
            "log",
            vec!["--abundance", "--abundance-bits", "8"],
            true,
            false,
        ),
        (
            "exact",
            vec!["--abundance", "--abundance-bits", "16"],
            true,
            false,
        ),
        (
            "log-filtered",
            vec![
                "--abundance",
                "--abundance-bits",
                "8",
                "--abundance-threshold",
                "3",
            ],
            true,
            true,
        ),
        (
            "exact-filtered",
            vec![
                "--abundance",
                "--abundance-bits",
                "16",
                "--abundance-threshold",
                "3",
            ],
            true,
            true,
        ),
    ] {
        let archive = format!("archive-{mode}");
        let output_dir = format!("restored-{mode}");
        let mut args = vec![
            "compress",
            "-i",
            "inputs.fof",
            "-o",
            &archive,
            "-d",
            "scratch",
            "-k",
            "17",
            "-m",
            "7",
            "-t",
            "1",
            "-r",
            "1",
        ];
        args.extend(flags);
        let output = run(root, &args);
        assert_success(&output);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(
            stderr
                .matches("Warning: skipping missing input file")
                .count(),
            3,
            "{stderr}"
        );
        for (path, line) in [
            ("missing-before.fa", 2),
            ("missing-middle.fa", 4),
            ("missing-after.fa", 7),
        ] {
            assert!(
                stderr
                    .lines()
                    .any(|text| text.contains(path) && text.contains(&format!("line {line}"))),
                "{stderr}"
            );
        }
        assert!(root.join(&archive).join("manifest.kloe").is_file());
        let names = fs::read_to_string(root.join(&archive).join("filenames_id.txt")).unwrap();
        let names: Vec<_> = names
            .lines()
            .map(|line| line.rsplit_once(':').unwrap().0)
            .collect();
        assert_eq!(names, ["first.fa", "second.fa.zst"]);

        let mut args = vec![
            "decompress",
            "-c",
            &archive,
            "-o",
            &output_dir,
            "--stored-tigs",
            "--output-compression",
            "fasta",
            "-t",
            "1",
            "-r",
            "1",
        ];
        if abundance {
            args.push("--abundance");
        }
        assert_success(&run(root, &args));
        let mut outputs: Vec<_> = fs::read_dir(root.join(&output_dir))
            .unwrap()
            .map(|entry| entry.unwrap().path())
            .collect();
        outputs.sort();
        assert_eq!(outputs.len(), 2);
        for (id, path) in outputs.iter().enumerate() {
            assert!(path
                .file_name()
                .unwrap()
                .to_str()
                .unwrap()
                .starts_with(&format!("Dump_{id}_")));
            let fasta = fs::read_to_string(path).unwrap();
            let expected = if id == 1 {
                vec![second.clone()]
            } else if filtered {
                vec![first.clone()]
            } else {
                vec![first.clone(), low.clone()]
            };
            assert_eq!(
                kmers(&sequences(&fasta)),
                kmers(&expected),
                "{mode}: dataset {id}"
            );
            if abundance {
                let low_kmers = kmers(&[low.clone()]);
                for record in fasta.split('>').skip(1) {
                    let header = record.lines().next().unwrap();
                    let value: f64 = header
                        .split_whitespace()
                        .find_map(|token| token.strip_prefix("ka:f:"))
                        .expect("abundance header")
                        .parse()
                        .unwrap();
                    let sequence: String = record.lines().skip(1).collect();
                    for kmer in kmers(&[sequence]) {
                        let expected = if id == 1 {
                            40.0
                        } else if low_kmers.contains(&kmer) {
                            3.0
                        } else {
                            20.0
                        };
                        let tolerance = if mode.starts_with("exact") { 0.0 } else { 1.0 };
                        assert!(
                            (value - expected).abs() <= tolerance,
                            "{mode}: dataset {id}: abundance {value}, expected {expected}"
                        );
                    }
                }
            }
        }

        // The second retained file must also be selectable by its original name.
        let targeted = format!("targeted-{mode}");
        assert_success(&run(
            root,
            &[
                "decompress",
                "-c",
                &archive,
                "-o",
                &targeted,
                "-Q",
                "selected.txt",
                "--stored-tigs",
                "--output-compression",
                "fasta",
                "-t",
                "1",
                "-r",
                "1",
            ],
        ));
        let outputs: Vec<_> = fs::read_dir(root.join(targeted))
            .unwrap()
            .map(|e| e.unwrap().path())
            .collect();
        assert_eq!(outputs.len(), 1);
        assert_eq!(
            kmers(&sequences(&fs::read_to_string(&outputs[0]).unwrap())),
            kmers(&[second.clone()])
        );
    }
}

#[test]
fn no_available_inputs_fail_without_a_complete_archive() {
    let temp = tempfile::tempdir().unwrap();
    for (name, contents) in [
        ("missing", "absent.fa\nalso-absent.fa\n"),
        ("empty", "\n  \n"),
    ] {
        let list = format!("{name}.fof");
        fs::write(temp.path().join(&list), contents).unwrap();
        let output = run(
            temp.path(),
            &["compress", "-i", &list, "-o", name, "-r", "1"],
        );
        assert!(!output.status.success());
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(stderr.contains("no available input files"), "{stderr}");
        assert!(!temp.path().join(name).join("manifest.kloe").exists());
        assert!(!temp.path().join(name).join(".kloe-write.lock").exists());
    }
}

#[test]
fn missing_list_and_other_input_errors_are_not_silently_skipped() {
    let temp = tempfile::tempdir().unwrap();
    let output = run(
        temp.path(),
        &[
            "compress",
            "-i",
            "absent.fof",
            "-o",
            "missing-list",
            "-r",
            "1",
        ],
    );
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("cannot open input list 'absent.fof'"));

    // ENOTDIR is a malformed path, not a missing file to ignore.
    fs::write(temp.path().join("not-a-directory"), "data").unwrap();
    fs::write(temp.path().join("bad.fof"), "not-a-directory/input.fa\n").unwrap();
    let output = run(
        temp.path(),
        &["compress", "-i", "bad.fof", "-o", "bad-input", "-r", "1"],
    );
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("cannot open input file 'not-a-directory/input.fa'"),
        "{stderr}"
    );
    assert!(!stderr.contains("skipping missing input file"), "{stderr}");
    assert!(!temp.path().join("bad-input/manifest.kloe").exists());
}
