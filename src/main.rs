#![allow(dead_code)]

mod compress;
mod decompress;
mod merge;
mod parser;
mod utils;
use clap::Parser;
use std::path::Path;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    ///output: Option<String>,
    /// Command mode: "compress", "decompress", or "merge"
    decompress: Option<String>,
    ///Input file list (Compression)
    #[arg(short, long, default_value_t=String::from(""))]
    input_list: String,
    /// Number of threads (defaults to 1)
    #[arg(short, long, default_value_t = num_cpus::get())]
    threads: usize,
    ///Output directory
    #[arg(short, long, default_value_t = String::from(""))]
    out_dir: String,
    ///input directory for decompression
    #[arg(short, long, default_value_t = String::from(""))]
    compressed_dir: String,
    ///second input archive directory for merge mode
    #[arg(long, default_value_t = String::from(""))]
    compressed_dir_b: String,
    ///List of files to decompress
    #[arg(short = 'Q', long, default_value_t = String::from(""))]
    wanted_files: String,
    ///Temporary directory for unitigs parsing
    #[arg(short = 'd', long, default_value_t = String::from(""))]
    temp_dir: String,
    ///Required memory budget in GB for compression/decompression workflows
    #[arg(short = 'r', long)]
    memory: Option<usize>,
    ///K value, default = 31
    #[arg(short, long, default_value_t = 31)]
    k_size: usize,
    ///Minimizer size (< k), default = 7
    #[arg(short, long, default_value_t = 7)]
    minimizer_size: usize,
    /// Optionally verify that all canonical k-mers are preserved with the correct dataset IDs
    #[arg(long = "verify-kmers", default_value_t = false)]
    verify_kmers: bool,
    /// Number of partitions is 2^partition_power (default 1024 partitions)
    #[arg(short = 'P', long = "partition-power", default_value_t = 10)]
    partition_power: u32,
    /// Skip sorting within partitions and during final merge (output will not be globally sorted)
    #[arg(long = "skip-sort", default_value_t = false)]
    skip_sort: bool,
    /// Produce monochromatic unitigs instead of simplitigs
    #[arg(long = "unitig", default_value_t = false)]
    unitig: bool,
    /// Produce monochromatic matchtigs instead of simplitigs
    #[arg(long = "matchtig", default_value_t = false)]
    matchtig: bool,
    /// Produce monochromatic eulertigs instead of simplitigs
    #[arg(long = "eulertig", default_value_t = false)]
    eulertig: bool,
    /// During decompression, rebuild tigs from Dump_*.fa using ggcat
    #[arg(long = "ggcat-rebuild", default_value_t = false)]
    ggcat_rebuild: bool,
}
fn main() {
    let args = Args::parse();

    let output_dir = args.out_dir;
    let input_dir = args.compressed_dir;
    let input_dir_b = args.compressed_dir_b;
    //env::set_var("RAYON_NUM_THREADS", args.threads.to_string());
    let input_fof = args.input_list;
    let threads = args.threads;
    let temp_dir = args.temp_dir;
    let k = args.k_size;
    let m = args.minimizer_size;
    //TODO HANDLE ERRORS FOR COMP AND DECOMP
    let wanted_path = args.wanted_files;
    let use_unitigs = args.unitig;
    let use_matchtigs = args.matchtig;
    let use_eulertigs = args.eulertig;
    let ggcat_rebuild = args.ggcat_rebuild;
    let tig_flags_set = [use_unitigs, use_matchtigs, use_eulertigs]
        .iter()
        .filter(|&&f| f)
        .count();
    if tig_flags_set > 1 {
        eprintln!("Error: only one of --unitig, --matchtig, --eulertig can be set at a time.");
        std::process::exit(1);
    }
    if let Some(do_decompress) = args.decompress {
        let Some(memory) = args.memory else {
            eprintln!("Error: please provide a memory budget with -r/--memory <GB>.");
            std::process::exit(2);
        };
        if do_decompress == "decompress" {
            println!("Checking archive integrity...");
            is_compressed_dir_complete(input_dir.clone());
            let rebuild_cfg = decompress::GgcatRebuildConfig {
                enabled: ggcat_rebuild,
                threads,
                memory_gb: memory,
                k,
                temp_dir: temp_dir.clone(),
                use_unitigs,
                use_matchtigs,
                use_eulertigs,
            };
            if let Err(err) = decompress::decompress_with_options(
                &String::from("bucket_sizes.txt"),
                &String::from("id_to_color_id.txt.zst"),
                &String::from("tigs_kloe.fa"),
                &String::from("positions_kloe.bin"),
                &String::from("filenames_id.txt"),
                &output_dir,
                &wanted_path,
                input_dir,
                rebuild_cfg,
            ) {
                eprintln!("decompression failed: {err}");
                std::process::exit(1);
            }
            //let _ = graph_build::init_decompress(String::from("bucket_sizes.txt.zst"), String::from("id_to_color_id.txt.zst"), unitigs_file, &output_dir, &wanted_path, &input_dir);
        } else if do_decompress == "compress" {
            /*let _ = compress::compress(
                &output_dir,
                &input_fof,
                threads,
                &temp_dir,
                k,
                m,
                partition_power,
                compaction_threads
            );*/
            //parser::run_parser(k, m, 10_u32, PathBuf::from(output_dir), PathBuf::from(input_fof), threads, compaction_threads, false);
            let ggcat_cfg = compress::GgcatCompressionConfig {
                memory_gb: memory,
                temp_dir,
            };
            if let Err(err) = compress::compress_with_ggcat(
                &output_dir,
                &input_fof,
                threads,
                k,
                m,
                args.partition_power,
                args.verify_kmers,
                args.skip_sort,
                use_unitigs,
                use_matchtigs,
                use_eulertigs,
                ggcat_cfg,
            ) {
                eprintln!("compression failed: {err}");
                std::process::exit(1);
            }
            //let _ = graph_build::build_graphs(&output_dir, &input_fof, &threads, &temp_dir, &memory);
        } else if do_decompress == "merge" {
            if input_dir.is_empty() || input_dir_b.is_empty() {
                eprintln!(
                    "merge requires two archives: use -c/--compressed-dir for archive A and --compressed-dir-b for archive B"
                );
                std::process::exit(2);
            }
            let merge_cfg = merge::MergeConfig {
                threads,
                minimizer_size: m,
                partition_power: args.partition_power,
                temp_dir: temp_dir.clone(),
                verify_kmers: args.verify_kmers,
                skip_sort: args.skip_sort,
                use_unitigs,
                use_matchtigs,
                use_eulertigs,
            };
            if let Err(err) = merge::merge_archives_with_config(
                &input_dir,
                &input_dir_b,
                &output_dir,
                k,
                memory,
                merge_cfg,
            ) {
                eprintln!("merge failed: {err}");
                std::process::exit(1);
            }
        }
    } else {
        println!(
            "Wrong positional arguments given. Values are 'compress', 'decompress', or 'merge'"
        );
        println!("Ex: if compression: I=my/fof.txt cargo r -r -- compress -f my_file_of_file.txt -o out_dir/ -t 12");
        println!("Ex: if decompression: I=my/fof.txt cargo r -r -- decompress -f my_file_of_file.txt --omnicolor-file out_dir/omnicolor.fa.zstd --multicolor-file out_dir/multicolor.fa.zstd -t 12");
        println!("Ex: if merge: cargo r -r -- merge -c archive_a --compressed-dir-b archive_b -o merged_archive -r 16");
    }
}

fn is_compressed_dir_complete(input_dir: String) {
    if !Path::new(&format!("{input_dir}/filenames_id.txt")).exists() {
        panic!("file not found: {input_dir}/filenames_id.txt");
    } else if !Path::new(&format!("{input_dir}/positions_kloe.bin")).exists() {
        panic!("Positions file not found");
    } else if !Path::new(&format!("{input_dir}/bucket_sizes.txt")).exists() {
        panic!("Tigs sizes file not found");
    } else if !Path::new(&format!("{input_dir}/id_to_color_id.txt.zst")).exists() {
        panic!("id to color id file not found");
    } else if !Path::new(&format!("{input_dir}/tigs_kloe.fa")).exists() {
        panic!("Tigs file not found");
    } else {
        println!("Archive complete, starting decompression...");
    }
}

#[cfg(test)]
mod tests {
    use super::{compress, decompress, merge};
    use std::collections::BTreeSet;
    use std::fs::{self, File};
    use std::io::Write;
    use std::path::{Path, PathBuf};
    use tempfile::TempDir;

    const K: usize = 31;
    const MODES: [&str; 4] = ["simplitigs", "unitigs", "matchtigs", "eulertigs"];
    const SEQ_SHARED_12: &str = "ACGTACGTTGCAACGTACGTTGCAACGTACGTACGTA";
    const SEQ_ONLY_1: &str = "TTGCAACGTACGTTTGCAACGTACGTTTGCAACGTAC";
    const SEQ_ONLY_2: &str = "GGATCCGGATCCGGATCCGGATCCGGATCCGGATCCA";
    const SEQ_SHARED_23: &str = "CCGTAACCGTAACCGTAACCGTAACCGTAACCGTAAC";
    const SEQ_ONLY_3: &str = "ATGCTTATGCTTATGCTTATGCTTATGCTTATGCTTA";

    fn write_fasta(path: &Path, seqs: &[&str]) {
        let mut fasta = String::new();
        for (i, seq) in seqs.iter().enumerate() {
            fasta.push_str(&format!(">seq{}\n{}\n", i + 1, seq));
        }
        fs::write(path, fasta).expect("write source fasta");
    }

    fn parse_fasta_sequences(path: &Path) -> Vec<String> {
        let content = fs::read_to_string(path).expect("read fasta");
        let mut seqs = Vec::new();
        let mut current = String::new();
        for line in content.lines() {
            if line.starts_with('>') {
                if !current.is_empty() {
                    seqs.push(current.clone());
                    current.clear();
                }
            } else if !line.trim().is_empty() {
                current.push_str(line.trim());
            }
        }
        if !current.is_empty() {
            seqs.push(current);
        }
        seqs
    }

    fn revcomp(seq: &str) -> String {
        seq.chars()
            .rev()
            .map(|b| match b {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                'T' => 'A',
                _ => 'N',
            })
            .collect()
    }

    fn canonical_kmer(kmer: &str) -> String {
        let rc = revcomp(kmer);
        if rc.as_str() < kmer {
            rc
        } else {
            kmer.to_string()
        }
    }

    fn kmer_set_from_seqs(seqs: &[String], k: usize) -> BTreeSet<String> {
        let mut kmers = BTreeSet::new();
        for seq in seqs {
            assert!(
                seq.len() >= k,
                "sequence shorter than k: len={}, k={}",
                seq.len(),
                k
            );
            for i in 0..=seq.len() - k {
                let kmer = &seq[i..i + k];
                kmers.insert(canonical_kmer(kmer));
            }
        }
        kmers
    }

    fn assert_kmer_equivalent(expected: &[String], dumped_fasta: &Path, k: usize) {
        let expected_kmers = kmer_set_from_seqs(expected, k);
        let dumped_kmers = kmer_set_from_seqs(&parse_fasta_sequences(dumped_fasta), k);
        assert_eq!(dumped_kmers, expected_kmers, "k-mer content mismatch");
    }

    fn kmer_set_from_fasta(path: &Path, k: usize) -> BTreeSet<String> {
        let seqs = parse_fasta_sequences(path);
        kmer_set_from_seqs(&seqs, k)
    }

    fn assert_kmer_maps_equal(left: &Path, right: &Path, k: usize, context: &str) {
        let left_map = kmer_set_from_fasta(left, k);
        let right_map = kmer_set_from_fasta(right, k);
        assert_eq!(left_map, right_map, "{}", context);
    }

    fn mode_archive_path(archive_dir: &Path, mode: &str) -> PathBuf {
        archive_dir.join(format!("{mode}.fa.zst"))
    }

    fn make_dump_path(out_dir: &Path, source_path: &Path) -> PathBuf {
        let stem = source_path
            .file_stem()
            .expect("missing file stem")
            .to_str()
            .expect("invalid utf8 stem");
        out_dir.join(format!("Dump_{stem}.fa"))
    }

    fn count_dump_files(out_dir: &Path) -> usize {
        fs::read_dir(out_dir)
            .expect("read output directory")
            .filter_map(|entry| entry.ok())
            .map(|entry| entry.path())
            .filter(|p| {
                p.file_name()
                    .and_then(|n| n.to_str())
                    .map(|name| name.starts_with("Dump_") && name.ends_with(".fa"))
                    .unwrap_or(false)
            })
            .count()
    }

    fn multi_case_expected() -> (Vec<String>, Vec<String>, Vec<String>) {
        (
            vec![SEQ_SHARED_12.to_string(), SEQ_ONLY_1.to_string()],
            vec![
                SEQ_SHARED_12.to_string(),
                SEQ_ONLY_2.to_string(),
                SEQ_SHARED_23.to_string(),
            ],
            vec![SEQ_SHARED_23.to_string(), SEQ_ONLY_3.to_string()],
        )
    }

    fn multi_case_records() -> Vec<(String, String)> {
        vec![
            ("1,2".to_string(), SEQ_SHARED_12.to_string()),
            ("1".to_string(), SEQ_ONLY_1.to_string()),
            ("2".to_string(), SEQ_ONLY_2.to_string()),
            ("2,3".to_string(), SEQ_SHARED_23.to_string()),
            ("3".to_string(), SEQ_ONLY_3.to_string()),
        ]
    }

    fn multi_case_records_permuted() -> Vec<(String, String)> {
        let records = multi_case_records();
        vec![
            records[3].clone(),
            records[1].clone(),
            records[4].clone(),
            records[0].clone(),
            records[2].clone(),
        ]
    }

    fn run_compression(
        input_files: &[PathBuf],
        output_dir: &Path,
        k: usize,
    ) -> Result<(), Box<dyn std::error::Error>> {
        // Create a file-of-files list
        let fof_path = output_dir.join("input.fof");
        let mut fof = File::create(&fof_path)?;
        for file in input_files {
            writeln!(fof, "{}", file.display())?;
        }
        drop(fof);

        // Run the streaming compression
        let output_dir_str = output_dir.display().to_string() + "/";
        compress::compress(
            &output_dir_str,
            &fof_path.display().to_string(),
            1, // threads
            k,
            7,     // minimizer size
            10,    // partition power
            false, // verify_kmers
            false, // skip_sort
            false, // use_unitigs
            false, // use_matchtigs
            false, // use_eulertigs
        )?;
        Ok(())
    }

    fn write_filenames_id(path: &Path, files: &[PathBuf], offsets: &[usize]) {
        assert_eq!(
            files.len(),
            offsets.len(),
            "files and offsets length mismatch"
        );
        let mut out = File::create(path).expect("create filenames_id.txt");
        for i in 0..files.len() {
            writeln!(out, "{}:{}", files[i].display(), offsets[i]).expect("write filenames_id row");
        }
    }

    struct MultiFixture {
        _workdir: TempDir,
        archive_dir: PathBuf,
        full_out_dir: PathBuf,
        targeted_out_dir: PathBuf,
        file1: PathBuf,
        file3: PathBuf,
        dump1_full: PathBuf,
        dump2_full: PathBuf,
        dump3_full: PathBuf,
        dump1_targeted: PathBuf,
        dump2_targeted: PathBuf,
        dump3_targeted: PathBuf,
        expected1: Vec<String>,
        expected2: Vec<String>,
        expected3: Vec<String>,
    }

    struct LargeFixture {
        _workdir: TempDir,
        archive_dir: PathBuf,
        full_out_dir: PathBuf,
        targeted_out_dir: PathBuf,
        source_files: Vec<PathBuf>,
        expected_seqs: Vec<String>,
        dump_full: Vec<PathBuf>,
        dump_targeted: Vec<PathBuf>,
        targeted_indices: Vec<usize>,
    }

    fn setup_multi_fixture(_mode: &str, _records: Vec<(String, String)>) -> MultiFixture {
        let workdir = tempfile::tempdir().expect("create temp workdir");
        let archive_dir = workdir.path().join("archive");
        let full_out_dir = workdir.path().join("out_full");
        let targeted_out_dir = workdir.path().join("out_targeted");
        fs::create_dir_all(&archive_dir).expect("create archive directory");
        fs::create_dir_all(&full_out_dir).expect("create full output directory");
        fs::create_dir_all(&targeted_out_dir).expect("create targeted output directory");

        let file1 = workdir.path().join("sample1.fa");
        let file2 = workdir.path().join("sample2.fa");
        let file3 = workdir.path().join("sample3.fa");

        let (expected1, expected2, expected3) = multi_case_expected();
        let exp1_refs: Vec<&str> = expected1.iter().map(String::as_str).collect();
        let exp2_refs: Vec<&str> = expected2.iter().map(String::as_str).collect();
        let exp3_refs: Vec<&str> = expected3.iter().map(String::as_str).collect();
        write_fasta(&file1, &exp1_refs);
        write_fasta(&file2, &exp2_refs);
        write_fasta(&file3, &exp3_refs);

        // Run streaming compression instead of write_records_archive
        run_compression(
            &[file1.clone(), file2.clone(), file3.clone()],
            &archive_dir,
            K,
        )
        .expect("compression failed");

        MultiFixture {
            _workdir: workdir,
            archive_dir,
            full_out_dir: full_out_dir.clone(),
            targeted_out_dir: targeted_out_dir.clone(),
            file1: file1.clone(),
            file3: file3.clone(),
            dump1_full: make_dump_path(&full_out_dir, &file1),
            dump2_full: make_dump_path(&full_out_dir, &file2),
            dump3_full: make_dump_path(&full_out_dir, &file3),
            dump1_targeted: make_dump_path(&targeted_out_dir, &file1),
            dump2_targeted: make_dump_path(&targeted_out_dir, &file2),
            dump3_targeted: make_dump_path(&targeted_out_dir, &file3),
            expected1,
            expected2,
            expected3,
        }
    }

    fn generate_len100_sequence(index: usize) -> String {
        let bases = ['A', 'C', 'G', 'T'];
        let mut seq = String::with_capacity(100);

        let mut x = index;
        let mut prefix = ['A'; 20];
        for i in (0..20).rev() {
            prefix[i] = bases[x & 0b11];
            x >>= 2;
        }
        for c in prefix {
            seq.push(c);
        }

        let mut state = (index as u64)
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        for _ in 20..100 {
            state = state
                .wrapping_mul(2862933555777941757)
                .wrapping_add(3037000493);
            seq.push(match state & 0b11 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'T',
            });
        }

        seq
    }

    fn setup_large_fixture(_mode: &str, n_files: usize) -> LargeFixture {
        assert!(n_files >= 100, "stress test expects at least 100 files");

        let workdir = tempfile::tempdir().expect("create temp workdir");
        let archive_dir = workdir.path().join("archive");
        let full_out_dir = workdir.path().join("out_full");
        let targeted_out_dir = workdir.path().join("out_targeted");
        fs::create_dir_all(&archive_dir).expect("create archive directory");
        fs::create_dir_all(&full_out_dir).expect("create full output directory");
        fs::create_dir_all(&targeted_out_dir).expect("create targeted output directory");

        let mut source_files = Vec::with_capacity(n_files);
        let mut expected_seqs = Vec::with_capacity(n_files);
        for i in 0..n_files {
            let source = workdir.path().join(format!("sample_{i:05}.fa"));
            let seq = generate_len100_sequence(i);
            expected_seqs.push(seq.clone());
            write_fasta(&source, &[seq.as_str()]);
            source_files.push(source);
        }

        // Run streaming compression (mode parameter is unused but kept for test consistency)
        run_compression(&source_files, &archive_dir, K).expect("compression failed");

        let dump_full: Vec<PathBuf> = source_files
            .iter()
            .map(|p| make_dump_path(&full_out_dir, p))
            .collect();
        let dump_targeted: Vec<PathBuf> = source_files
            .iter()
            .map(|p| make_dump_path(&targeted_out_dir, p))
            .collect();

        let targeted_indices = vec![0, 1, 17, 23, 50, 79, 99];

        LargeFixture {
            _workdir: workdir,
            archive_dir,
            full_out_dir,
            targeted_out_dir,
            source_files,
            expected_seqs,
            dump_full,
            dump_targeted,
            targeted_indices,
        }
    }

    fn reverse_complement(seq: &str) -> String {
        seq.chars()
            .rev()
            .map(|base| match base {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                'T' => 'A',
                _ => panic!("invalid nucleotide in test sequence: {}", base),
            })
            .collect()
    }

    fn hamming_distance(a: &str, b: &str) -> usize {
        assert_eq!(a.len(), b.len(), "hamming distance needs same length");
        a.bytes().zip(b.bytes()).filter(|(x, y)| x != y).count()
    }

    struct TwoFileCaseFixture {
        _workdir: TempDir,
        archive_dir: PathBuf,
        full_out_dir: PathBuf,
        targeted_out_dir: PathBuf,
        case_file: PathBuf,
        case_dump_full: PathBuf,
        case_dump_targeted: PathBuf,
        control_dump_full: PathBuf,
        control_dump_targeted: PathBuf,
    }

    fn setup_two_file_case_fixture(
        _mode: &str,
        case_filename: &str,
        case_seqs: &[String],
        control_filename: &str,
        control_seqs: &[String],
    ) -> TwoFileCaseFixture {
        let workdir = tempfile::tempdir().expect("create temp workdir");
        let archive_dir = workdir.path().join("archive");
        let full_out_dir = workdir.path().join("out_full");
        let targeted_out_dir = workdir.path().join("out_targeted");
        fs::create_dir_all(&archive_dir).expect("create archive directory");
        fs::create_dir_all(&full_out_dir).expect("create full output directory");
        fs::create_dir_all(&targeted_out_dir).expect("create targeted output directory");

        let case_file = workdir.path().join(case_filename);
        let control_file = workdir.path().join(control_filename);
        let case_refs: Vec<&str> = case_seqs.iter().map(String::as_str).collect();
        let control_refs: Vec<&str> = control_seqs.iter().map(String::as_str).collect();
        write_fasta(&case_file, &case_refs);
        write_fasta(&control_file, &control_refs);

        // Run streaming compression instead of write_records_archive
        run_compression(&[case_file.clone(), control_file.clone()], &archive_dir, K)
            .expect("compression failed");

        TwoFileCaseFixture {
            _workdir: workdir,
            archive_dir,
            full_out_dir: full_out_dir.clone(),
            targeted_out_dir: targeted_out_dir.clone(),
            case_file: case_file.clone(),
            case_dump_full: make_dump_path(&full_out_dir, &case_file),
            case_dump_targeted: make_dump_path(&targeted_out_dir, &case_file),
            control_dump_full: make_dump_path(&full_out_dir, &control_file),
            control_dump_targeted: make_dump_path(&targeted_out_dir, &control_file),
        }
    }

    fn run_full_decompression(archive_dir: &Path, out_dir: &Path) {
        let archive_dir_s = format!("{}/", archive_dir.display());
        let out_dir_s = format!("{}/", out_dir.display());
        decompress::decompress(
            &String::from("bucket_sizes.txt"),
            &String::from("id_to_color_id.txt.zst"),
            &String::from("tigs_kloe.fa"),
            &String::from("positions_kloe.bin"),
            &String::from("filenames_id.txt"),
            &out_dir_s,
            &String::from(""),
            archive_dir_s,
        )
        .expect("decompress full archive");
    }

    fn run_targeted_decompression(archive_dir: &Path, out_dir: &Path, wanted_files: &[PathBuf]) {
        let wanted = archive_dir.join("wanted.txt");
        let mut wanted_content = String::new();
        for file in wanted_files {
            wanted_content.push_str(&format!("{}\n", file.display()));
        }
        fs::write(&wanted, wanted_content).expect("write wanted files list");

        let archive_dir_s = format!("{}/", archive_dir.display());
        let out_dir_s = format!("{}/", out_dir.display());
        decompress::decompress(
            &String::from("bucket_sizes.txt"),
            &String::from("id_to_color_id.txt.zst"),
            &String::from("tigs_kloe.fa"),
            &String::from("positions_kloe.bin"),
            &String::from("filenames_id.txt"),
            &out_dir_s,
            &wanted.display().to_string(),
            archive_dir_s,
        )
        .expect("decompress targeted archive");
    }

    #[test]
    fn multi_file_full_and_targeted_preserve_kmer_content_all_modes() {
        for mode in MODES {
            let fixture = setup_multi_fixture(mode, multi_case_records());

            run_full_decompression(&fixture.archive_dir, &fixture.full_out_dir);
            assert_kmer_equivalent(&fixture.expected1, &fixture.dump1_full, K);
            assert_kmer_equivalent(&fixture.expected2, &fixture.dump2_full, K);
            assert_kmer_equivalent(&fixture.expected3, &fixture.dump3_full, K);

            run_targeted_decompression(
                &fixture.archive_dir,
                &fixture.targeted_out_dir,
                &[fixture.file1.clone(), fixture.file3.clone()],
            );
            assert_kmer_equivalent(&fixture.expected1, &fixture.dump1_targeted, K);
            assert_kmer_equivalent(&fixture.expected3, &fixture.dump3_targeted, K);
            assert!(
                !fixture.dump2_targeted.exists(),
                "mode={mode}: targeted decompression should not create non-target files"
            );
        }
    }

    #[test]
    fn reverse_complement_case_preserves_content_in_full_and_targeted_decompression() {
        for mode in MODES {
            let seq = "ACGTTAGCCATGATCGTACCGTTAGGCTAACCGTTAACGA".to_string();
            let rc = reverse_complement(&seq);
            assert_ne!(
                seq, rc,
                "test sequence should not be self reverse-complement"
            );
            let case_expected = vec![seq.clone(), rc.clone()];
            let control_expected = vec!["GGTACCGATCGTACCGATCGTACCGATCGTACCGAT".to_string()];

            let fixture = setup_two_file_case_fixture(
                mode,
                "reverse_complement.fa",
                &case_expected,
                "reverse_complement_control.fa",
                &control_expected,
            );

            run_full_decompression(&fixture.archive_dir, &fixture.full_out_dir);
            assert_kmer_equivalent(&case_expected, &fixture.case_dump_full, K);
            assert_kmer_equivalent(&control_expected, &fixture.control_dump_full, K);

            run_targeted_decompression(
                &fixture.archive_dir,
                &fixture.targeted_out_dir,
                &[fixture.case_file.clone()],
            );
            assert_kmer_equivalent(&case_expected, &fixture.case_dump_targeted, K);
            assert!(
                !fixture.control_dump_targeted.exists(),
                "mode={mode}: targeted decompression should not create non-target files"
            );
        }
    }

    #[test]
    fn single_nucleotide_difference_case_preserves_content_in_full_and_targeted_decompression() {
        for mode in MODES {
            let seq1 = "TTGACCGTTAACCGGTTAACCGGTTAACCGGTTAACCGGT".to_string();
            let mut seq2_bytes = seq1.as_bytes().to_vec();
            let idx = 20;
            seq2_bytes[idx] = if seq2_bytes[idx] == b'A' { b'C' } else { b'A' };
            let seq2 = String::from_utf8(seq2_bytes).expect("build mutated sequence");
            assert_eq!(
                hamming_distance(&seq1, &seq2),
                1,
                "mutated sequence should differ by one nucleotide"
            );
            let case_expected = vec![seq1.clone(), seq2.clone()];
            let control_expected = vec!["AACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTT".to_string()];

            let fixture = setup_two_file_case_fixture(
                mode,
                "single_nt_diff.fa",
                &case_expected,
                "single_nt_diff_control.fa",
                &control_expected,
            );

            run_full_decompression(&fixture.archive_dir, &fixture.full_out_dir);
            assert_kmer_equivalent(&case_expected, &fixture.case_dump_full, K);
            assert_kmer_equivalent(&control_expected, &fixture.control_dump_full, K);

            run_targeted_decompression(
                &fixture.archive_dir,
                &fixture.targeted_out_dir,
                &[fixture.case_file.clone()],
            );
            assert_kmer_equivalent(&case_expected, &fixture.case_dump_targeted, K);
            assert!(
                !fixture.control_dump_targeted.exists(),
                "mode={mode}: targeted decompression should not create non-target files"
            );
        }
    }

    #[test]
    fn mutation_plus_reverse_complement_case_preserves_content_in_full_and_targeted_decompression()
    {
        for mode in MODES {
            let seq = "GCTAACCGTTAACCGGTTACCGATGCTAACCGTTAACCGG".to_string();
            let perfect_rc = reverse_complement(&seq);
            let mut mutated_rc_bytes = perfect_rc.as_bytes().to_vec();
            let idx = 15;
            mutated_rc_bytes[idx] = if mutated_rc_bytes[idx] == b'A' {
                b'C'
            } else {
                b'A'
            };
            let mutated_rc =
                String::from_utf8(mutated_rc_bytes).expect("build mutation+reverse-complement");

            assert_eq!(
                hamming_distance(&perfect_rc, &mutated_rc),
                1,
                "mutated reverse complement should differ by one nucleotide"
            );

            let case_expected = vec![seq.clone(), mutated_rc.clone()];
            let control_expected = vec!["TTGCCATTGCCATTGCCATTGCCATTGCCATTGCCAT".to_string()];

            let fixture = setup_two_file_case_fixture(
                mode,
                "mut_plus_rc.fa",
                &case_expected,
                "mut_plus_rc_control.fa",
                &control_expected,
            );

            run_full_decompression(&fixture.archive_dir, &fixture.full_out_dir);
            assert_kmer_equivalent(&case_expected, &fixture.case_dump_full, K);
            assert_kmer_equivalent(&control_expected, &fixture.control_dump_full, K);

            run_targeted_decompression(
                &fixture.archive_dir,
                &fixture.targeted_out_dir,
                &[fixture.case_file.clone()],
            );
            assert_kmer_equivalent(&case_expected, &fixture.case_dump_targeted, K);
            assert!(
                !fixture.control_dump_targeted.exists(),
                "mode={mode}: targeted decompression should not create non-target files"
            );
        }
    }

    #[test]
    fn order_invariance_preserves_full_and_targeted_decompression_all_modes() {
        for mode in MODES {
            let baseline = setup_multi_fixture(mode, multi_case_records());
            let permuted = setup_multi_fixture(mode, multi_case_records_permuted());

            run_full_decompression(&baseline.archive_dir, &baseline.full_out_dir);
            run_full_decompression(&permuted.archive_dir, &permuted.full_out_dir);

            assert_kmer_maps_equal(
                &baseline.dump1_full,
                &permuted.dump1_full,
                K,
                &format!("mode={mode}: full decompression differs for sample1"),
            );
            assert_kmer_maps_equal(
                &baseline.dump2_full,
                &permuted.dump2_full,
                K,
                &format!("mode={mode}: full decompression differs for sample2"),
            );
            assert_kmer_maps_equal(
                &baseline.dump3_full,
                &permuted.dump3_full,
                K,
                &format!("mode={mode}: full decompression differs for sample3"),
            );

            run_targeted_decompression(
                &baseline.archive_dir,
                &baseline.targeted_out_dir,
                &[baseline.file1.clone(), baseline.file3.clone()],
            );
            run_targeted_decompression(
                &permuted.archive_dir,
                &permuted.targeted_out_dir,
                &[permuted.file1.clone(), permuted.file3.clone()],
            );

            assert_kmer_maps_equal(
                &baseline.dump1_targeted,
                &permuted.dump1_targeted,
                K,
                &format!("mode={mode}: targeted decompression differs for sample1"),
            );
            assert_kmer_maps_equal(
                &baseline.dump3_targeted,
                &permuted.dump3_targeted,
                K,
                &format!("mode={mode}: targeted decompression differs for sample3"),
            );

            assert!(
                !baseline.dump2_targeted.exists(),
                "mode={mode}: targeted baseline should not create non-target files"
            );
            assert!(
                !permuted.dump2_targeted.exists(),
                "mode={mode}: targeted permuted should not create non-target files"
            );
        }
    }

    #[test]
    fn large_100_files_len100_preserve_full_and_targeted_kmers_all_modes() {
        for mode in MODES {
            let fixture = setup_large_fixture(mode, 100);

            run_full_decompression(&fixture.archive_dir, &fixture.full_out_dir);
            assert_eq!(
                count_dump_files(&fixture.full_out_dir),
                100,
                "mode={mode}: full decompression should emit one file per input"
            );
            for i in 0..fixture.expected_seqs.len() {
                let expected = vec![fixture.expected_seqs[i].clone()];
                assert_kmer_equivalent(&expected, &fixture.dump_full[i], K);
            }

            let wanted_files: Vec<PathBuf> = fixture
                .targeted_indices
                .iter()
                .map(|&i| fixture.source_files[i].clone())
                .collect();
            run_targeted_decompression(
                &fixture.archive_dir,
                &fixture.targeted_out_dir,
                &wanted_files,
            );
            assert_eq!(
                count_dump_files(&fixture.targeted_out_dir),
                fixture.targeted_indices.len(),
                "mode={mode}: targeted decompression should only emit requested files"
            );
            for &i in &fixture.targeted_indices {
                let expected = vec![fixture.expected_seqs[i].clone()];
                assert_kmer_equivalent(&expected, &fixture.dump_targeted[i], K);
            }
            assert!(
                !fixture.dump_targeted[2].exists(),
                "mode={mode}: non-target file should not be emitted"
            );
        }
    }

    #[test]
    fn various_k_below_31_preserved_in_full_and_targeted_all_modes() {
        let ks = [3_usize, 5, 7, 11, 17, 23, 29, 30];

        for mode in MODES {
            let fixture = setup_multi_fixture(mode, multi_case_records());
            run_full_decompression(&fixture.archive_dir, &fixture.full_out_dir);
            run_targeted_decompression(
                &fixture.archive_dir,
                &fixture.targeted_out_dir,
                &[fixture.file1.clone(), fixture.file3.clone()],
            );

            for k in ks {
                assert_kmer_equivalent(&fixture.expected1, &fixture.dump1_full, k);
                assert_kmer_equivalent(&fixture.expected2, &fixture.dump2_full, k);
                assert_kmer_equivalent(&fixture.expected3, &fixture.dump3_full, k);
                assert_kmer_equivalent(&fixture.expected1, &fixture.dump1_targeted, k);
                assert_kmer_equivalent(&fixture.expected3, &fixture.dump3_targeted, k);
            }

            assert!(
                !fixture.dump2_targeted.exists(),
                "mode={mode}: targeted decompression should not create non-target files"
            );
        }
    }

    #[test]
    fn merge_two_archives_preserves_per_file_kmer_content() {
        let workdir = tempfile::tempdir().expect("create temp workdir");
        let archive_a = workdir.path().join("archive_a");
        let archive_b = workdir.path().join("archive_b");
        let merged_archive = workdir.path().join("archive_merged");
        let merged_out = workdir.path().join("merged_out");
        fs::create_dir_all(&archive_a).expect("create archive_a directory");
        fs::create_dir_all(&archive_b).expect("create archive_b directory");
        fs::create_dir_all(&merged_archive).expect("create merged archive directory");
        fs::create_dir_all(&merged_out).expect("create merged output directory");

        let shared = "ACGTACGTTGCAACGTACGTTGCAACGTACGTACGTA".to_string();
        let a1_only = "TTGCAACGTACGTTTGCAACGTACGTTTGCAACGTAC".to_string();
        let a2_only = "GGATCCGGATCCGGATCCGGATCCGGATCCGGATCCA".to_string();
        let b2_only = "CCGTAACCGTAACCGTAACCGTAACCGTAACCGTAAC".to_string();

        let a_file1 = workdir.path().join("a_sample1.fa");
        let a_file2 = workdir.path().join("a_sample2.fa");
        let b_file1 = workdir.path().join("b_sample1.fa");
        let b_file2 = workdir.path().join("b_sample2.fa");

        write_fasta(&a_file1, &[shared.as_str(), a1_only.as_str()]);
        write_fasta(&a_file2, &[a2_only.as_str()]);
        write_fasta(&b_file1, &[shared.as_str()]);
        write_fasta(&b_file2, &[b2_only.as_str()]);

        run_compression(&[a_file1.clone(), a_file2.clone()], &archive_a, K)
            .expect("compress archive A");
        run_compression(&[b_file1.clone(), b_file2.clone()], &archive_b, K)
            .expect("compress archive B");

        merge::merge_archives(
            &archive_a.display().to_string(),
            &archive_b.display().to_string(),
            &merged_archive.display().to_string(),
            K,
            1,
        )
        .expect("merge archives");

        run_full_decompression(&merged_archive, &merged_out);

        assert_kmer_equivalent(
            &[shared.clone(), a1_only.clone()],
            &make_dump_path(&merged_out, &a_file1),
            K,
        );
        assert_kmer_equivalent(
            &[a2_only.clone()],
            &make_dump_path(&merged_out, &a_file2),
            K,
        );
        assert_kmer_equivalent(&[shared.clone()], &make_dump_path(&merged_out, &b_file1), K);
        assert_kmer_equivalent(
            &[b2_only.clone()],
            &make_dump_path(&merged_out, &b_file2),
            K,
        );
    }
}
