use std::collections::BTreeSet;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

const K: usize = 17;

fn run(args: &[&str]) {
    let output = Command::new(env!("CARGO_BIN_EXE_kloe"))
        .args(args)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "args={args:?}\n{}\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn read_fasta(path: &Path) -> Vec<String> {
    let (mut reader, format) =
        niffler::get_reader(Box::new(fs::File::open(path).unwrap())).unwrap();
    let expected_format = match path.extension().and_then(|ext| ext.to_str()) {
        Some("zst") => niffler::compression::Format::Zstd,
        Some("gz") => niffler::compression::Format::Gzip,
        Some("xz") => niffler::compression::Format::Lzma,
        Some("fa") => niffler::compression::Format::No,
        other => panic!("unexpected FASTA extension: {other:?}"),
    };
    assert_eq!(
        format,
        expected_format,
        "wrong compression format in {}",
        path.display()
    );
    let mut text = String::new();
    reader.read_to_string(&mut text).unwrap();
    let mut sequences = Vec::new();
    for record in text.split('>').skip(1) {
        sequences.push(record.lines().skip(1).collect::<String>());
    }
    sequences
}

fn rc(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|base| match base {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => panic!(),
        })
        .collect()
}

fn kmers(sequences: &[String]) -> BTreeSet<String> {
    sequences
        .iter()
        .flat_map(|sequence| {
            sequence.as_bytes().windows(K).map(|window| {
                let forward = String::from_utf8(window.to_vec()).unwrap();
                forward.clone().min(rc(&forward))
            })
        })
        .collect()
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

#[test]
fn cli_rebuilds_each_dataset_and_explicit_sets_in_all_output_modes() {
    let temp = tempfile::tempdir().unwrap();
    let archive = temp.path().join("archive");
    let scratch = temp.path().join("scratch");
    fs::create_dir(&archive).unwrap();
    let a = dna(12345);
    let b = a[40..120].to_owned();
    let c = dna(67890);
    assert!(kmers(&[a.clone()]).is_disjoint(&kmers(&[c.clone()])));
    // Same basename deliberately exercises independent, collision-free output.
    let sources: Vec<PathBuf> = (0..3)
        .map(|id| {
            let dir = temp.path().join(format!("input{id}"));
            fs::create_dir(&dir).unwrap();
            dir.join("sample.fa")
        })
        .collect();
    for (source, sequence) in sources.iter().zip([&a, &b, &c]) {
        fs::write(source, format!(">input\n{sequence}\n")).unwrap();
    }
    let fof = temp.path().join("input.fof");
    fs::write(
        &fof,
        sources
            .iter()
            .map(|p| format!("{}\n", p.display()))
            .collect::<String>(),
    )
    .unwrap();
    run(&[
        "compress",
        "-i",
        fof.to_str().unwrap(),
        "-o",
        &(archive.display().to_string() + "/"),
        "-d",
        scratch.to_str().unwrap(),
        "-k",
        "17",
        "-m",
        "7",
        "-t",
        "1",
        "-r",
        "1",
    ]);

    // No trailing slash and no -k: reconstruction must use the manifest's k=17.
    let base = [
        "decompress",
        "-c",
        archive.to_str().unwrap(),
        "-r",
        "1",
        "-t",
        "1",
    ];
    for (flag, mode) in [
        (None, "unitigs"),
        (Some("--simplitig"), "simplitigs"),
        (Some("--matchtig"), "matchtigs"),
        (Some("--eulertig"), "eulertigs"),
    ] {
        for (compression, suffix) in [
            ("zstd", ".fa.zst"),
            ("fasta", ".fa"),
            ("gz", ".fa.gz"),
            ("xz", ".fa.xz"),
        ] {
            let out = temp.path().join(format!("{mode}-{compression}"));
            let mut args = base.to_vec();
            args.extend(["-o", out.to_str().unwrap()]);
            if let Some(flag) = flag {
                args.push(flag);
            }
            // Exercise both defaults together for the default unitig/zstd case.
            if compression != "zstd" {
                args.extend(["--output-compression", compression]);
            }
            run(&args);
            assert_eq!(fs::read_dir(&out).unwrap().count(), 3);
            for (id, expected) in [&a, &b, &c].into_iter().enumerate() {
                let sequences = read_fasta(&out.join(format!("Dump_{id}_sample.{mode}{suffix}")));
                assert_eq!(
                    kmers(&sequences),
                    kmers(&[expected.clone()]),
                    "{mode}/{compression} dataset {id}"
                );
                if mode == "unitigs" {
                    // A's stored graph is split by B's color boundaries. Its
                    // independent graph is a single maximal, nonbranching path.
                    assert_eq!(sequences.len(), 1);
                    assert!(sequences[0] == *expected || sequences[0] == rc(expected));
                }
            }
        }
    }

    let selected = temp.path().join("selected.txt");
    fs::write(&selected, format!("{}\n", sources[1].display())).unwrap();
    let targeted = temp.path().join("targeted");
    let mut args = base.to_vec();
    args.extend([
        "-Q",
        selected.to_str().unwrap(),
        "-o",
        targeted.to_str().unwrap(),
        "--unitig",
    ]);
    run(&args);
    assert_eq!(fs::read_dir(&targeted).unwrap().count(), 1);
    assert_eq!(
        kmers(&read_fasta(&targeted.join("Dump_1_sample.unitigs.fa.zst"))),
        kmers(&[b.clone()])
    );

    fs::write(
        &selected,
        format!("{}\n{}\n", sources[0].display(), sources[1].display()),
    )
    .unwrap();
    let sets = temp.path().join("sets");
    fs::create_dir(&sets).unwrap();
    let stale = sets.join("Dump_stale.fa.zst");
    let stale_bytes = zstd::encode_all(format!(">stale\n{c}\n").as_bytes(), 1).unwrap();
    fs::write(&stale, &stale_bytes).unwrap();
    for (operation, expected) in [("union", &a), ("intersection", &b)] {
        let mut args = base.to_vec();
        args.extend([
            "-Q",
            selected.to_str().unwrap(),
            "-o",
            sets.to_str().unwrap(),
            "--color-set-operation",
            operation,
        ]);
        run(&args);
        let records = read_fasta(&sets.join(format!("Dump_{operation}.unitigs.fa.zst")));
        assert_eq!(kmers(&records), kmers(&[expected.clone()]));
        assert_eq!(records.len(), 1);
    }
    assert_eq!(fs::read(&stale).unwrap(), stale_bytes);
    assert_eq!(fs::read_dir(&sets).unwrap().count(), 3);

    fs::write(
        &selected,
        format!("{}\n{}\n", sources[0].display(), sources[2].display()),
    )
    .unwrap();
    let mut args = base.to_vec();
    args.extend([
        "-Q",
        selected.to_str().unwrap(),
        "-o",
        sets.to_str().unwrap(),
        "--color-set-operation",
        "intersection",
    ]);
    run(&args);
    assert!(read_fasta(&sets.join("Dump_intersection.unitigs.fa.zst")).is_empty());
}
