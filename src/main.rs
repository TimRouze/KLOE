#![allow(dead_code)]

mod abundance;
mod compress;
mod decompress;
mod merge;
mod packed_tigs;
mod records;
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
    ///Minimizer size (< k), default = 12
    #[arg(short, long, default_value_t = 12)]
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
    /// Produce unitigs (default for compression, decompression and graph merge)
    #[arg(long = "unitig", default_value_t = false)]
    unitig: bool,
    /// Produce matchtigs (experimental; not yet benchmarked)
    #[arg(long = "matchtig", default_value_t = false)]
    matchtig: bool,
    /// Produce eulertigs (experimental; not yet benchmarked)
    #[arg(long = "eulertig", default_value_t = false)]
    eulertig: bool,
    /// Produce simplitigs using embedded GGCAT (experimental; not yet benchmarked)
    #[arg(long = "simplitig", default_value_t = false)]
    simplitig: bool,
    /// Output stored tigs without rebuilding (required to restore abundance)
    #[arg(long = "stored-tigs", default_value_t = false)]
    stored_tigs: bool,
    /// Compatibility option: decompression already rebuilds each result by default
    #[arg(long = "ggcat-rebuild", default_value_t = false)]
    ggcat_rebuild: bool,
    /// Compression used for decompressed FASTA output
    #[arg(
        long = "output-compression",
        alias = "compression",
        value_enum,
        default_value_t = decompress::OutputCompression::Zstd
    )]
    output_compression: decompress::OutputCompression,
    /// Emit one FASTA for a set operation on colors listed with -Q (difference: first minus rest)
    #[arg(long = "color-set-operation", alias = "colors-operation", value_enum)]
    color_set_operation: Option<decompress::ColorSetOperation>,
    /// Join archive containers without rebuilding the graph or removing shared redundancy
    #[arg(long = "structural-merge", default_value_t = false)]
    structural_merge: bool,
    /// Compress: preserve abundance fields; decompress: restore them in FASTA headers
    #[arg(long = "abundance", default_value_t = false)]
    abundance: bool,
    /// Decompress: keep only k-mers with decoded abundance strictly above this integer
    #[arg(long = "abundance-threshold")]
    abundance_threshold: Option<u64>,
    /// Compress: preserve integer abundance exactly, without quantization or saturation
    #[arg(long = "abundance-exact", conflicts_with_all = ["abundance_bits", "abundance_max"])]
    abundance_exact: bool,
    /// Store rounded abundance directly as u16; reject values rounding above 65535
    #[arg(long = "abundance-u16", conflicts_with_all = ["abundance_exact", "abundance_bits", "abundance_max"])]
    abundance_u16: bool,
    /// Abundance bits: exact integers if the maximum fits, otherwise automatic logarithmic encoding
    #[arg(long = "abundance-bits", default_value_t = 8, requires = "abundance", value_parser = clap::value_parser!(u8).range(2..=32))]
    abundance_bits: u8,
    /// Maximum abundance; larger values are capped before storage
    #[arg(long = "abundance-max", alias = "max-abundance", default_value_t = 50000, requires = "abundance", value_parser = clap::value_parser!(u64).range(1..))]
    abundance_max: u64,
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
    let tig_flags_set = [use_unitigs, use_matchtigs, use_eulertigs, args.simplitig]
        .iter()
        .filter(|&&f| f)
        .count();
    if tig_flags_set > 1 {
        eprintln!("Error: only one of --unitig, --simplitig, --matchtig, --eulertig can be set at a time.");
        std::process::exit(1);
    }
    if let Some(do_decompress) = args.decompress {
        if args.abundance_threshold.is_some() && do_decompress != "decompress" {
            eprintln!("Error: --abundance-threshold is a decompression option.");
            std::process::exit(2);
        }
        if (args.abundance_exact || args.abundance_u16) && do_decompress != "compress" {
            eprintln!("Error: --abundance-exact and --abundance-u16 are compression options; decompress with --stored-tigs --abundance.");
            std::process::exit(2);
        }
        let Some(memory) = args.memory else {
            eprintln!("Error: please provide a memory budget with -r/--memory <GB>.");
            std::process::exit(2);
        };
        if do_decompress == "decompress" {
            if args.stored_tigs && (tig_flags_set > 0 || ggcat_rebuild) {
                eprintln!("Error: --stored-tigs cannot be combined with a rebuild mode.");
                std::process::exit(2);
            }
            if args.color_set_operation.is_some() && wanted_path.is_empty() {
                eprintln!(
                    "Error: --color-set-operation requires -Q/--wanted-files with the selected colors."
                );
                std::process::exit(2);
            }
            if args.color_set_operation.is_some() && args.abundance {
                eprintln!(
                    "Error: --abundance cannot be combined with --color-set-operation because one combined output has no single per-dataset abundance."
                );
                std::process::exit(2);
            }
            println!("Checking archive integrity...");
            is_compressed_dir_complete(input_dir.clone());
            let rebuild_cfg = decompress::GgcatRebuildConfig {
                enabled: !args.stored_tigs,
                threads,
                memory_gb: memory,
                k,
                temp_dir: temp_dir.clone(),
                use_unitigs,
                use_matchtigs,
                use_eulertigs,
                use_simplitigs: args.simplitig,
                restore_abundance: args.abundance,
                abundance_threshold: args.abundance_threshold,
                output_compression: args.output_compression,
                color_set_operation: args.color_set_operation,
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
            if args.stored_tigs {
                eprintln!("Error: --stored-tigs is a decompression option.");
                std::process::exit(2);
            }
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
                abundance_log_base: None,
                abundance_bits: args.abundance_bits,
                abundance_max: (args.abundance && !args.abundance_exact && !args.abundance_u16)
                    .then_some(args.abundance_max),
                abundance_exact: args.abundance_exact,
                abundance_u16: args.abundance_u16,
                use_simplitigs: args.simplitig,
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
            if args.stored_tigs {
                eprintln!("Error: --stored-tigs is a decompression option.");
                std::process::exit(2);
            }
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
                structural: args.structural_merge,
                use_simplitigs: args.simplitig,
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
    } else if !Path::new(&format!("{input_dir}/id_to_color_id.txt.zst")).exists()
        && !Path::new(&format!("{input_dir}/{}", compress::CID_TO_DATASET_FILE)).exists()
        && !Path::new(&format!("{input_dir}/{}", compress::DATASET_TO_CID_FILE)).exists()
    {
        panic!("archive membership index not found");
    } else if !Path::new(&format!("{input_dir}/tigs_kloe.fa")).exists() {
        panic!("Tigs file not found");
    } else if Path::new(&format!("{input_dir}/{}", compress::DATASET_TO_CID_FILE)).exists() {
        for index in [
            compress::ARCHIVE_MANIFEST_FILE,
            compress::FILENAME_INDEX_FILE,
            compress::POSITIONS_INDEX_FILE,
            compress::BUCKET_SIZES_INDEX_FILE,
            compress::TIGS_INDEX_FILE,
        ] {
            if !Path::new(&format!("{input_dir}/{index}")).exists() {
                panic!("KLOE v4 index not found: {input_dir}/{index}");
            }
        }
    }
    println!("Archive complete, starting decompression...");
}

#[cfg(test)]
mod tests {
    use super::{compress, decompress, merge};
    use std::collections::BTreeSet;
    use std::fs::{self, File};
    use std::io::{Read, Write};
    use std::path::{Path, PathBuf};
    use std::sync::{Mutex, MutexGuard, OnceLock};
    use tempfile::TempDir;

    const K: usize = 31;
    const MODES: [&str; 4] = ["simplitigs", "unitigs", "matchtigs", "eulertigs"];
    const SEQ_SHARED_12: &str = "ACGTACGTTGCAACGTACGTTGCAACGTACGTACGTA";
    const SEQ_ONLY_1: &str = "TTGCAACGTACGTTTGCAACGTACGTTTGCAACGTAC";
    const SEQ_ONLY_2: &str = "GGATCCGGATCCGGATCCGGATCCGGATCCGGATCCA";
    const SEQ_SHARED_23: &str = "CCGTAACCGTAACCGTAACCGTAACCGTAACCGTAAC";
    const SEQ_ONLY_3: &str = "ATGCTTATGCTTATGCTTATGCTTATGCTTATGCTTA";

    fn ggcat_test_lock() -> MutexGuard<'static, ()> {
        static LOCK: OnceLock<Mutex<()>> = OnceLock::new();
        LOCK.get_or_init(|| Mutex::new(())).lock().unwrap()
    }

    fn write_fasta(path: &Path, seqs: &[&str]) {
        let mut fasta = String::new();
        for (i, seq) in seqs.iter().enumerate() {
            fasta.push_str(&format!(">seq{}\n{}\n", i + 1, seq));
        }
        fs::write(path, fasta).expect("write source fasta");
    }

    fn parse_fasta_sequences(path: &Path) -> Vec<String> {
        if fs::metadata(path).expect("stat FASTA").len() == 0 {
            return Vec::new();
        }
        let file = File::open(path).expect("open fasta");
        let (mut reader, _) =
            niffler::get_reader(Box::new(file)).expect("detect fasta compression");
        let mut content = String::new();
        reader.read_to_string(&mut content).expect("read fasta");
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
        run_compression_with_simplitigs(input_files, output_dir, k, false)
    }

    fn run_compression_with_simplitigs(
        input_files: &[PathBuf],
        output_dir: &Path,
        k: usize,
        use_simplitigs: bool,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let _ggcat_guard = ggcat_test_lock();
        // Create a file-of-files list
        let fof_path = output_dir.join("input.fof");
        let mut fof = File::create(&fof_path)?;
        for file in input_files {
            writeln!(fof, "{}", file.display())?;
        }
        drop(fof);

        // Run the streaming compression
        let output_dir_str = output_dir.display().to_string() + "/";
        compress::compress_with_ggcat(
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
            compress::GgcatCompressionConfig {
                use_simplitigs,
                ..compress::GgcatCompressionConfig::default()
            },
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
        file2: PathBuf,
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
            file2: file2.clone(),
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
        decompress::decompress_with_options(
            &String::from("bucket_sizes.txt"),
            &String::from("id_to_color_id.txt.zst"),
            &String::from("tigs_kloe.fa"),
            &String::from("positions_kloe.bin"),
            &String::from("filenames_id.txt"),
            &out_dir_s,
            &String::from(""),
            archive_dir_s,
            decompress::GgcatRebuildConfig {
                enabled: false,
                output_compression: decompress::OutputCompression::Fasta,
                ..decompress::GgcatRebuildConfig::default()
            },
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
        decompress::decompress_with_options(
            &String::from("bucket_sizes.txt"),
            &String::from("id_to_color_id.txt.zst"),
            &String::from("tigs_kloe.fa"),
            &String::from("positions_kloe.bin"),
            &String::from("filenames_id.txt"),
            &out_dir_s,
            &wanted.display().to_string(),
            archive_dir_s,
            decompress::GgcatRebuildConfig {
                enabled: false,
                output_compression: decompress::OutputCompression::Fasta,
                ..decompress::GgcatRebuildConfig::default()
            },
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
    fn simplitig_archive_preserves_full_and_targeted_dataset_kmers() {
        let work = tempfile::tempdir().unwrap();
        let archive = work.path().join("archive");
        let full = work.path().join("full");
        let targeted = work.path().join("targeted");
        for path in [&archive, &full, &targeted] {
            fs::create_dir_all(path).unwrap();
        }
        let (first, second, third) = multi_case_expected();
        let expected = [first, second, third];
        let inputs: Vec<_> = expected
            .iter()
            .enumerate()
            .map(|(dataset, sequences)| {
                let file = work.path().join(format!("dataset{dataset}.fa"));
                write_fasta(
                    &file,
                    &sequences.iter().map(String::as_str).collect::<Vec<_>>(),
                );
                file
            })
            .collect();
        run_compression_with_simplitigs(&inputs, &archive, K, true).unwrap();
        run_full_decompression(&archive, &full);
        run_targeted_decompression(&archive, &targeted, &[inputs[1].clone()]);
        for (dataset, input) in inputs.iter().enumerate() {
            assert_kmer_equivalent(&expected[dataset], &make_dump_path(&full, input), K);
        }
        assert_kmer_equivalent(&expected[1], &make_dump_path(&targeted, &inputs[1]), K);
        assert_eq!(count_dump_files(&targeted), 1);
    }

    #[test]
    fn selected_color_union_and_intersection_emit_one_deduplicated_fasta() {
        let fixture = setup_multi_fixture("simplitigs", multi_case_records());
        let wanted = fixture.archive_dir.join("wanted_color_set.txt");
        fs::write(
            &wanted,
            format!("{}\n{}\n", fixture.file1.display(), fixture.file2.display()),
        )
        .unwrap();

        let run = |operation: decompress::ColorSetOperation, output: &Path| {
            fs::create_dir_all(output).unwrap();
            decompress::decompress_with_options(
                &"bucket_sizes.txt".to_string(),
                &"id_to_color_id.txt.zst".to_string(),
                &"tigs_kloe.fa".to_string(),
                &"positions_kloe.bin".to_string(),
                &"filenames_id.txt".to_string(),
                &(output.display().to_string() + "/"),
                &wanted.display().to_string(),
                fixture.archive_dir.display().to_string() + "/",
                decompress::GgcatRebuildConfig {
                    enabled: false,
                    output_compression: decompress::OutputCompression::Fasta,
                    color_set_operation: Some(operation),
                    ..decompress::GgcatRebuildConfig::default()
                },
            )
            .unwrap();
        };

        let union_out = fixture._workdir.path().join("out_union");
        run(decompress::ColorSetOperation::Union, &union_out);
        let mut union_expected = fixture.expected1.clone();
        union_expected.extend(fixture.expected2.clone());
        assert_kmer_equivalent(&union_expected, &union_out.join("Dump_union.fa"), K);

        let intersection_out = fixture._workdir.path().join("out_intersection");
        run(
            decompress::ColorSetOperation::Intersection,
            &intersection_out,
        );
        assert_kmer_equivalent(
            &[SEQ_SHARED_12.to_string()],
            &intersection_out.join("Dump_intersection.fa"),
            K,
        );
    }

    #[test]
    fn abundance_threshold_requires_metadata_and_nonnegative_integer() {
        use super::Args;
        use clap::Parser;
        assert!(Args::try_parse_from(["kloe", "decompress", "--abundance-threshold", "-1"]).is_err());
        assert!(Args::try_parse_from(["kloe", "decompress", "--abundance-threshold", "1.5"]).is_err());
        assert_eq!(Args::try_parse_from(["kloe", "decompress", "--abundance-threshold", "0"])
            .unwrap().abundance_threshold, Some(0));
        let fixture = setup_multi_fixture("simplitigs", multi_case_records());
        let wanted = fixture.archive_dir.join("wanted_threshold.txt");
        fs::write(&wanted, format!("{}\n", fixture.file1.display())).unwrap();
        for operation in [None, Some(decompress::ColorSetOperation::Union)] {
            let output = fixture._workdir.path().join(format!("no_metadata_{operation:?}"));
            let result = decompress::decompress_with_options(
                &"bucket_sizes.txt".to_string(), &"id_to_color_id.txt.zst".to_string(),
                &"tigs_kloe.fa".to_string(), &"positions_kloe.bin".to_string(),
                &"filenames_id.txt".to_string(), &(output.display().to_string() + "/"),
                &wanted.display().to_string(), fixture.archive_dir.display().to_string() + "/",
                decompress::GgcatRebuildConfig {
                    enabled: false, abundance_threshold: Some(0),
                    color_set_operation: operation,
                    ..decompress::GgcatRebuildConfig::default()
                },
            );
            assert!(result.unwrap_err().to_string().contains("abundance metadata"));
        }
    }

    #[test]
    fn selected_color_differences_match_kmer_set_algebra() {
        use decompress::{ColorSetOperation as Op, OutputCompression};
        let fixture = setup_multi_fixture("simplitigs", multi_case_records());
        let _ggcat_guard = ggcat_test_lock();
        let files = [&fixture.file1, &fixture.file2, &fixture.file3];
        let sets = [&fixture.expected1, &fixture.expected2, &fixture.expected3]
            .map(|seqs| kmer_set_from_seqs(seqs, K));
        // Check operand order, duplicate operands, and more than two operands,
        // both as stored tigs and after rebuilding the selected graph.
        for (case, operands) in [vec![0, 1], vec![1, 0], vec![0, 0], vec![0, 1, 2]]
            .into_iter().enumerate()
        {
            let wanted = fixture.archive_dir.join(format!("wanted_difference_{case}.txt"));
            fs::write(&wanted, operands.iter()
                .map(|&i| format!("{}\n\n", files[i].display())).collect::<String>()).unwrap();
            for op in [Op::Difference, Op::SymmetricDifference] {
                let mut expected = sets[operands[0]].clone();
                for &i in &operands[1..] {
                    expected = match op {
                        Op::Difference => expected.difference(&sets[i]).cloned().collect(),
                        Op::SymmetricDifference => expected.symmetric_difference(&sets[i]).cloned().collect(),
                        _ => unreachable!(),
                    };
                }
                for rebuilt in [false, true] {
                    let output = fixture._workdir.path().join(format!("out_{op}_{case}_{rebuilt}"));
                    fs::create_dir_all(&output).unwrap();
                    decompress::decompress_with_options(
                        &"bucket_sizes.txt".to_string(),
                        &"id_to_color_id.txt.zst".to_string(),
                        &"tigs_kloe.fa".to_string(),
                        &"positions_kloe.bin".to_string(),
                        &"filenames_id.txt".to_string(),
                        &(output.display().to_string() + "/"),
                        &wanted.display().to_string(),
                        fixture.archive_dir.display().to_string() + "/",
                        decompress::GgcatRebuildConfig {
                            enabled: rebuilt, threads: 1, memory_gb: 1, k: K,
                            temp_dir: output.join("tmp").display().to_string(),
                            use_unitigs: true,
                            output_compression: OutputCompression::Zstd,
                            color_set_operation: Some(op),
                            ..decompress::GgcatRebuildConfig::default()
                        },
                    ).unwrap();
                    let suffix = if rebuilt { ".unitigs.fa.zst" } else { ".fa.zst" };
                    let observed = parse_fasta_sequences(&output.join(format!("Dump_{op}{suffix}")));
                    assert_eq!(kmer_set_from_seqs(&observed, K), expected,
                        "operation={op}, operands={operands:?}, rebuilt={rebuilt}");
                }
            }
        }
    }

    #[test]
    fn decompression_unitig_mode_rebuilds_actual_unitigs() {
        let fixture = setup_multi_fixture("simplitigs", multi_case_records());
        let _ggcat_guard = ggcat_test_lock();
        for (compression, label, suffix) in [
            (decompress::OutputCompression::Zstd, "zstd", ".fa.zst"),
            (decompress::OutputCompression::Xz, "xz", ".fa.xz"),
        ] {
            let output = fixture
                ._workdir
                .path()
                .join(format!("out_rebuilt_unitigs_{label}"));
            let temp = fixture
                ._workdir
                .path()
                .join(format!("unitig_rebuild_tmp_{label}"));
            fs::create_dir_all(&output).unwrap();
            decompress::decompress_with_options(
                &"bucket_sizes.txt".to_string(),
                &"id_to_color_id.txt.zst".to_string(),
                &"tigs_kloe.fa".to_string(),
                &"positions_kloe.bin".to_string(),
                &"filenames_id.txt".to_string(),
                &(output.display().to_string() + "/"),
                &String::new(),
                fixture.archive_dir.display().to_string() + "/",
                decompress::GgcatRebuildConfig {
                    enabled: true,
                    threads: 1,
                    memory_gb: 1,
                    k: K,
                    temp_dir: temp.display().to_string(),
                    use_unitigs: true,
                    output_compression: compression,
                    ..decompress::GgcatRebuildConfig::default()
                },
            )
            .unwrap();

            for (id, expected) in [&fixture.expected1, &fixture.expected2, &fixture.expected3]
                .into_iter()
                .enumerate()
            {
                let rebuilt =
                    output.join(format!("Dump_{}_sample{}.unitigs{}", id, id + 1, suffix));
                assert!(
                    rebuilt.is_file(),
                    "{label} dataset {id} unitig output was not created"
                );
                assert_kmer_equivalent(expected, &rebuilt, K);
            }
            assert_eq!(fs::read_dir(&output).unwrap().count(), 3);
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
    fn indexed_abundance_archive_restores_targeted_headers_without_cid_sidecar() {
        check_abundance_archive_roundtrip(8, 1.05);
    }

    #[test]
    fn abundance_precision_modes_roundtrip() {
        for bits in [2, 4, 5, 6, 7, 8, 16, 32] {
            check_abundance_archive_roundtrip(bits, compress::default_abundance_log_base(bits));
        }
    }

    #[test]
    fn exact_abundance_archive_roundtrip() {
        check_abundance_archive_roundtrip(8, crate::abundance::EXACT_BASE);
    }

    #[test]
    fn u16_abundance_archive_roundtrip() {
        for simplitig in [false, true] {
            check_abundance_archive_roundtrip_mode(16, crate::abundance::EXACT_BASE, simplitig);
        }
    }

    #[test]
    fn u16_abundance_cli_rejects_conflicting_encodings() {
        use clap::Parser;
        assert!(super::Args::try_parse_from(["kloe", "compress", "--abundance-u16"]).unwrap().abundance_u16);
        for extra in [vec!["--abundance-exact"], vec!["--abundance-bits", "16"], vec!["--abundance-max", "500"]] {
            let mut args = vec!["kloe", "compress", "--abundance-u16"];
            args.extend(extra);
            assert!(super::Args::try_parse_from(args).is_err());
        }
    }

    #[test]
    fn exact_abundance_cli_rejects_quantization_options() {
        use clap::Parser;
        assert!(
            super::Args::try_parse_from(["kloe", "compress", "--abundance-exact"])
                .unwrap()
                .abundance_exact
        );
        for option in [["--abundance-bits", "8"], ["--abundance-max", "500"]] {
            assert!(super::Args::try_parse_from([
                "kloe",
                "compress",
                "--abundance-exact",
                option[0],
                option[1]
            ])
            .is_err());
        }
    }

    #[test]
    fn exact_abundance_many_values_and_merge_preserve_targeted_kmers() {
        check_many_abundances_and_merge(false, false, None);
    }

    #[test]
    fn simplitig_exact_abundance_and_merge_preserve_targeted_kmers() {
        check_many_abundances_and_merge(true, false, None);
    }

    #[test]
    fn u16_abundance_merge_preserves_targeted_kmers() {
        check_many_abundances_and_merge(false, true, None);
    }

    #[test]
    fn dense_quantized_abundance_merge_preserves_decoded_values() {
        for bits in [8, 14] {
            check_many_abundances_and_merge(false, false, Some(bits));
        }
    }

    fn check_many_abundances_and_merge(use_simplitigs: bool, u16_mode: bool, log_bits: Option<u8>) {
        use std::collections::BTreeMap;
        let _ggcat_guard = ggcat_test_lock();
        let work = tempfile::tempdir().unwrap();
        let first = work.path().join("first.fa");
        let second = work.path().join("second.fa");
        let mut expected = [BTreeMap::new(), BTreeMap::new()];
        let mut inputs = [String::new(), String::new()];
        for index in 0..300usize {
            let middle: String = (0..5)
                .map(|shift| b"ACGT"[(index >> (2 * shift)) & 3] as char)
                .collect();
            let sequence = format!("AAA{middle}CCCCCCC");
            for dataset in 0..2 {
                if dataset == 1 && index >= 100 {
                    continue;
                }
                let value = if dataset == 0 {
                    index as u64 + 1
                } else if log_bits.is_some() {
                    50_000 - index as u64
                } else if u16_mode {
                    u16::MAX as u64 - index as u64
                } else {
                    u64::MAX - index as u64
                };
                inputs[dataset].push_str(&format!(">r{index} ka:f:{value}\n{sequence}\n"));
                let decoded = log_bits.map_or(value, |bits| {
                    let base = compress::default_abundance_log_base(bits);
                    compress::decode_abundance(compress::encode_abundance_with_bits(value, base, bits), base)
                });
                expected[dataset].insert(canonical_kmer(&sequence), decoded);
            }
        }
        fs::write(&first, &inputs[0]).unwrap();
        fs::write(&second, &inputs[1]).unwrap();
        let build = |name: &str, files: &[&Path]| {
            let archive = work.path().join(name);
            fs::create_dir_all(&archive).unwrap();
            let fof = work.path().join(format!("{name}.fof"));
            fs::write(
                &fof,
                files
                    .iter()
                    .map(|file| format!("{}\n", file.display()))
                    .collect::<String>(),
            )
            .unwrap();
            compress::compress_with_ggcat(
                &(archive.display().to_string() + "/"),
                &fof.display().to_string(),
                1,
                15,
                5,
                10,
                false,
                false,
                false,
                false,
                false,
                compress::GgcatCompressionConfig {
                    memory_gb: 1,
                    temp_dir: work.path().join("tmp").display().to_string(),
                    abundance_exact: !u16_mode && log_bits.is_none(),
                    abundance_u16: u16_mode,
                    abundance_bits: log_bits.unwrap_or(8),
                    abundance_max: log_bits.map(|_| 50_000),
                    use_simplitigs,
                    ..compress::GgcatCompressionConfig::default()
                },
            )
            .unwrap();
            archive
        };
        let direct = build("direct", &[&first, &second]);
        let a = build("a", &[&first]);
        let b = build("b", &[&second]);
        let merged = work.path().join("merged");
        merge::merge_archives_with_config(
            &a.display().to_string(),
            &b.display().to_string(),
            &merged.display().to_string(),
            15,
            1,
            merge::MergeConfig {
                minimizer_size: 5,
                temp_dir: work.path().join("tmp").display().to_string(),
                use_simplitigs,
                ..merge::MergeConfig::default()
            },
        )
        .unwrap();
        for (archive_index, archive) in [direct, merged].iter().enumerate() {
            for (dataset, source) in [&first, &second].iter().enumerate() {
                let output = work.path().join(format!("out_{archive_index}_{dataset}"));
                fs::create_dir_all(&output).unwrap();
                let wanted = work.path().join("wanted.txt");
                fs::write(&wanted, format!("{}\n", source.display())).unwrap();
                decompress::decompress_with_options(
                    &"bucket_sizes.txt".to_string(),
                    &"id_to_color_id.txt.zst".to_string(),
                    &"tigs_kloe.fa".to_string(),
                    &"positions_kloe.bin".to_string(),
                    &"filenames_id.txt".to_string(),
                    &(output.display().to_string() + "/"),
                    &wanted.display().to_string(),
                    archive.display().to_string() + "/",
                    decompress::GgcatRebuildConfig {
                        enabled: false,
                        restore_abundance: true,
                        output_compression: decompress::OutputCompression::Fasta,
                        ..decompress::GgcatRebuildConfig::default()
                    },
                )
                .unwrap();
                let text = fs::read_to_string(make_dump_path(&output, source)).unwrap();
                let mut actual = BTreeMap::new();
                let mut value = None;
                for line in text.lines() {
                    if line.starts_with('>') {
                        value = Some(crate::abundance::parse_header(line.as_bytes()).unwrap());
                    } else {
                        for kmer in line.as_bytes().windows(15) {
                            actual.insert(
                                canonical_kmer(std::str::from_utf8(kmer).unwrap()),
                                value.unwrap(),
                            );
                        }
                    }
                }
                assert_eq!(actual, expected[dataset]);
                assert_eq!(fs::read_dir(&output).unwrap().count(), 1);
            }
        }
    }

    #[test]
    fn automatic_abundance_roundtrips_exact_logarithmic_and_capped_values() {
        let _ggcat_guard = ggcat_test_lock();
        let work = tempfile::tempdir().unwrap();
        for (bits, maximum) in [(8, 50_000u64), (10, 50_000), (12, 50_000), (14, 50_000), (4, 500), (8, 255), (8, 256),
                                (9, 500), (16, 50_000), (32, 1_000_000), (2, 3)] {
            let dir = work.path().join(format!("{bits}-{maximum}"));
            let archive = dir.join("archive");
            fs::create_dir_all(&archive).unwrap();
            let first = dir.join("first.fa");
            let second = dir.join("second.fa");
            let value = (maximum / 3).max(1);
            fs::write(&first, format!(">a ka:f:{value}.25\nACGTACGTACGTACGT\n")).unwrap();
            fs::write(&second, format!(">b ka:f:{}\nACGTACGTACGTACGT\n", maximum + 100)).unwrap();
            let fof = dir.join("input.fof");
            fs::write(&fof, format!("{}\n{}\n", first.display(), second.display())).unwrap();
            compress::compress_with_ggcat(
                &(archive.display().to_string() + "/"), &fof.display().to_string(),
                1, 7, 5, 10, false, false, false, false, false,
                compress::GgcatCompressionConfig {
                    memory_gb: 1, temp_dir: work.path().join("tmp").display().to_string(),
                    abundance_bits: bits, abundance_max: Some(maximum),
                    ..compress::GgcatCompressionConfig::default()
                },
            ).unwrap();
            let output = dir.join("out");
            decompress::decompress_with_options(
                &"bucket_sizes.txt".to_string(), &"id_to_color_id.txt.zst".to_string(),
                &"tigs_kloe.fa".to_string(), &"positions_kloe.bin".to_string(),
                &"filenames_id.txt".to_string(), &(output.display().to_string() + "/"),
                &String::new(), archive.display().to_string() + "/",
                decompress::GgcatRebuildConfig {
                    enabled: false, restore_abundance: true,
                    output_compression: decompress::OutputCompression::Fasta,
                    ..decompress::GgcatRebuildConfig::default()
                },
            ).unwrap();
            let encoding = crate::abundance::Encoding::automatic(bits, maximum).unwrap();
            let expected_first = if encoding.base == crate::abundance::EXACT_BASE { value }
                else { compress::decode_abundance(compress::encode_abundance_with_bits(value, encoding.base, bits), encoding.base) };
            for (source, expected) in [(&first, expected_first), (&second, maximum)] {
                let dump = make_dump_path(&output, source);
                assert_kmer_equivalent(&["ACGTACGTACGTACGT".to_string()], &dump, 7);
                for header in fs::read_to_string(dump).unwrap().lines().filter(|line| line.starts_with('>')) {
                    assert_eq!(crate::abundance::parse_header(header.as_bytes()).unwrap(), expected,
                        "bits={bits}, maximum={maximum}, source={source:?}");
                }
            }
        }
    }

    #[test]
    fn abundance_precision_cli_accepts_only_supported_bits() {
        use clap::Parser;
        let defaults = super::Args::try_parse_from(["kloe", "compress", "--abundance"]).unwrap();
        assert_eq!(defaults.abundance_bits, 8);
        assert_eq!(defaults.abundance_max, 50_000);
        assert!(super::Args::try_parse_from(["kloe", "compress", "--abundance", "--abundance-max", "0"]).is_err());
        assert!(super::Args::try_parse_from(["kloe", "compress", "--abundance-bits", "16"]).is_err());
        for bits in [2, 4, 5, 6, 7, 8, 16, 32] {
            let args = super::Args::try_parse_from([
                "kloe",
                "compress",
                "--abundance",
                "--abundance-bits",
                &bits.to_string(),
                "--abundance-max",
                "500",
            ])
            .unwrap();
            assert_eq!(args.abundance_bits, bits);
            assert_eq!(args.abundance_max, 500);
        }
        for bits in ["0", "1", "33", "255"] {
            assert!(super::Args::try_parse_from([
                "kloe",
                "compress",
                "--abundance",
                "--abundance-bits",
                bits,
            ])
            .is_err());
        }
    }

    fn check_abundance_archive_roundtrip(bits: u8, base: f64) {
        check_abundance_archive_roundtrip_mode(bits, base, false);
    }

    #[test]
    fn simplitig_quantized_abundance_roundtrip() {
        for bits in [2, 4, 5, 6, 7, 8, 16, 32] {
            check_abundance_archive_roundtrip_mode(
                bits,
                compress::default_abundance_log_base(bits),
                true,
            );
        }
    }

    fn check_abundance_archive_roundtrip_mode(bits: u8, base: f64, use_simplitigs: bool) {
        let exact = base == crate::abundance::EXACT_BASE;
        let u16_mode = exact && bits == 16;
        let _ggcat_guard = ggcat_test_lock();
        let workdir = tempfile::tempdir().expect("create temp workdir");
        let archive = workdir.path().join("archive");
        let output = workdir.path().join("output");
        fs::create_dir_all(&archive).unwrap();
        fs::create_dir_all(&output).unwrap();
        let first = workdir.path().join("abundance_a.fa");
        let second = workdir.path().join("abundance_b.fa");
        let first_value = if u16_mode { 6143 } else if exact { 9_007_199_254_740_993u64 } else { 17 };
        let second_value = if u16_mode { u16::MAX as u64 } else if exact { u64::MAX } else { 50_000 };
        let first_header = if u16_mode { "6142.75".to_string() } else { first_value.to_string() };
        fs::write(&first, format!(">a ka:f:{first_header}\nACGTACGTACGTACGT\n")).unwrap();
        fs::write(
            &second,
            format!(">b ka:f:{second_value}\nACGTACGTACGTACGT\n"),
        )
        .unwrap();
        let fof = workdir.path().join("input.fof");
        fs::write(&fof, format!("{}\n{}\n", first.display(), second.display())).unwrap();
        compress::compress_with_ggcat(
            &(archive.display().to_string() + "/"),
            &fof.display().to_string(),
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
                memory_gb: 1,
                temp_dir: workdir.path().join("tmp").display().to_string(),
                abundance_log_base: (!exact).then_some(base),
                abundance_bits: bits,
                abundance_max: None,
                abundance_exact: exact && !u16_mode,
                abundance_u16: u16_mode,
                use_simplitigs,
            },
        )
        .unwrap();
        assert!(!archive.join(compress::CID_TO_DATASET_FILE).exists());

        let wanted = workdir.path().join("wanted.txt");
        fs::write(
            &wanted,
            format!("{}\n{}\n", first.display(), second.display()),
        )
        .unwrap();
        decompress::decompress_with_options(
            &"bucket_sizes.txt".to_string(),
            &"id_to_color_id.txt.zst".to_string(),
            &"tigs_kloe.fa".to_string(),
            &"positions_kloe.bin".to_string(),
            &"filenames_id.txt".to_string(),
            &(output.display().to_string() + "/"),
            &wanted.display().to_string(),
            archive.display().to_string() + "/",
            decompress::GgcatRebuildConfig {
                enabled: false,
                restore_abundance: true,
                output_compression: decompress::OutputCompression::Fasta,
                ..decompress::GgcatRebuildConfig::default()
            },
        )
        .unwrap();
        let first_dump = fs::read_to_string(make_dump_path(&output, &first)).unwrap();
        let second_dump = fs::read_to_string(make_dump_path(&output, &second)).unwrap();
        let restored_17 = if exact {
            first_value
        } else {
            compress::decode_abundance(
                compress::encode_abundance_with_bits(first_value, base, bits),
                base,
            )
        };
        let restored_max = if exact {
            second_value
        } else {
            compress::decode_abundance(
                compress::encode_abundance_with_bits(second_value, base, bits),
                base,
            )
        };
        assert!(
            first_dump.contains(&format!("ka:f:{restored_17}")),
            "{first_dump}"
        );
        assert!(
            second_dump.contains(&format!("ka:f:{restored_max}")),
            "{second_dump}"
        );
        // Both colors share their k-mers but have different abundances. Filtering
        // must restore each residual with its CID base before comparing values.
        for (case, (targeted, headers, rebuilt, threshold)) in [
            (true, true, false, restored_17),
            (false, false, false, restored_17),
            (true, false, true, restored_17),
            (true, false, false, restored_max),
        ].into_iter().enumerate() {
            let filtered = workdir.path().join(format!("filtered_{case}"));
            let wanted_path = if targeted { wanted.display().to_string() } else { String::new() };
            decompress::decompress_with_options(
                &"bucket_sizes.txt".to_string(), &"id_to_color_id.txt.zst".to_string(),
                &"tigs_kloe.fa".to_string(), &"positions_kloe.bin".to_string(),
                &"filenames_id.txt".to_string(), &(filtered.display().to_string() + "/"),
                &wanted_path, archive.display().to_string() + "/",
                decompress::GgcatRebuildConfig {
                    enabled: rebuilt, threads: 1, memory_gb: 1, k: 7,
                    temp_dir: workdir.path().join("filter_tmp").display().to_string(),
                    restore_abundance: headers, abundance_threshold: Some(threshold),
                    output_compression: decompress::OutputCompression::Fasta,
                    ..decompress::GgcatRebuildConfig::default()
                },
            ).unwrap();
            for (id, source) in [&first, &second].into_iter().enumerate() {
                let path = if rebuilt {
                    filtered.join(format!("Dump_{id}_{}.unitigs.fa", source.file_stem().unwrap().to_str().unwrap()))
                } else { make_dump_path(&filtered, source) };
                let sequences = parse_fasta_sequences(&path);
                let expected = if id == 1 && restored_max > threshold {
                    vec!["ACGTACGTACGTACGT".to_string()]
                } else { vec![] };
                assert_eq!(kmer_set_from_seqs(&sequences, 7), kmer_set_from_seqs(&expected, 7));
                if !headers { assert!(!fs::read_to_string(path).unwrap().contains("ka:f:")); }
            }
        }
        for op in [decompress::ColorSetOperation::Union, decompress::ColorSetOperation::Intersection,
                   decompress::ColorSetOperation::Difference, decompress::ColorSetOperation::SymmetricDifference] {
            let filtered = workdir.path().join(format!("filtered_{op}"));
            decompress::decompress_with_options(
                &"bucket_sizes.txt".to_string(), &"id_to_color_id.txt.zst".to_string(),
                &"tigs_kloe.fa".to_string(), &"positions_kloe.bin".to_string(),
                &"filenames_id.txt".to_string(), &(filtered.display().to_string() + "/"),
                &wanted.display().to_string(), archive.display().to_string() + "/",
                decompress::GgcatRebuildConfig {
                    enabled: false, abundance_threshold: Some(restored_17),
                    color_set_operation: Some(op),
                    output_compression: decompress::OutputCompression::Fasta,
                    ..decompress::GgcatRebuildConfig::default()
                },
            ).unwrap();
            let observed = parse_fasta_sequences(&filtered.join(format!("Dump_{op}.fa")));
            let expected = match op {
                decompress::ColorSetOperation::Union | decompress::ColorSetOperation::SymmetricDifference =>
                    vec!["ACGTACGTACGTACGT".to_string()],
                _ => vec![],
            };
            assert_eq!(kmer_set_from_seqs(&observed, 7), kmer_set_from_seqs(&expected, 7), "filtered {op}");
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
        for archive in [&archive_a, &archive_b] {
            assert!(archive.join(compress::DATASET_TO_CID_FILE).is_file());
            assert!(
                !archive.join(compress::CID_TO_DATASET_FILE).exists()
                    && !archive.join("id_to_color_id.txt.zst").exists(),
                "new archives must not duplicate the membership relation"
            );
        }

        {
            let _ggcat_guard = ggcat_test_lock();
            merge::merge_archives(
                &archive_a.display().to_string(),
                &archive_b.display().to_string(),
                &merged_archive.display().to_string(),
                K,
                1,
            )
            .expect("merge archives");
        }
        assert!(merged_archive.join(compress::DATASET_TO_CID_FILE).is_file());
        assert!(!merged_archive.join("id_to_color_id.txt.zst").exists());

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
