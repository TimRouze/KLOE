use bio::io::fasta;
use parking_lot::Mutex;
use simd_minimizers::packed_seq::{PackedSeqVec, SeqVec};
use std::collections::HashMap;
use std::env;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::thread;
use zstd::stream::write::Encoder as ZstdEncoder;

fn process_file(
    file_path: &Path,
    file_id: usize,
    partition_files: &Arc<Mutex<HashMap<u64, ZstdEncoder<File>>>>,
    k: usize,
    m: usize,
) -> io::Result<()> {
    println!("[Thread {:?}] Processing file {} (Path: \"{}\")", thread::current().id(), file_id, file_path.display());
    
    let k_super = k;
    let k_mini = m as u8;
    let w = (k_super - k_mini as usize + 1) as u8;
    let k_minus_m = k_super - k_mini as usize;

    let file = File::open(file_path)?;
    let reader = fasta::Reader::new(file);
    let mut local_buffers: HashMap<u64, Vec<u8>> = HashMap::new();
    let num_partitions = 4u64.pow(m as u32);

    for result in reader.records() {
        let record = result.expect("Failed to parse FASTA record");
        let seq = record.seq();

        let seq_preview = if seq.len() > 100 { &seq[..100] } else { seq };
        // println!("\n[DEBUG] Input Sequence ID: {} (len={}), Preview: {}...", record.id(), seq.len(), String::from_utf8_lossy(seq_preview));

        if seq.len() < k_super {
            continue;
        }

        let packed_seq = PackedSeqVec::from_ascii(seq);
        let mut minimizer_pos = Vec::new();
        minimizer_pos = simd_minimizers::canonical_minimizer_positions(packed_seq.as_slice(), k_mini as usize, w as usize);//, &mut minimizer_pos);
        
        if minimizer_pos.is_empty() {
            continue;
        }

        //simd_minimizers::iter_canonical_minimizer_values(packed_seq.as_slice(), k_mini as usize, &minimizer_pos).collect();
        let minimizer_values: Vec<_> = simd_minimizers::canonical_minimizers(k_mini as usize, w as usize).run_once(packed_seq.as_slice());
        //let minimizer_values: Vec<_> = simd_minimizers::canonical_minimizers(k_mini as usize, w as usize).;//, &minimizer_pos).collect();
        let minimizers: Vec<_> = minimizer_pos.into_iter().zip(minimizer_values.into_iter()).collect();

        if minimizers.is_empty() {
            continue;
        }
        
        // println!("[DEBUG] Found {} minimizers: {:?}", minimizers.len(), minimizers);
        
        for (current_pos, current_hash) in minimizers.iter() {
            let start_pos = (*current_pos as usize).saturating_sub(k_minus_m);
            let end_pos = std::cmp::min(seq.len(), (*current_pos as usize) + k_super);
            
            if end_pos > start_pos {
                let superkmer_slice = &seq[start_pos..end_pos];
                let partition_id = (*current_hash as u64 % num_partitions) as u64;

                let sk_preview = if superkmer_slice.len() > 60 { &superkmer_slice[..60] } else { superkmer_slice };
                // println!("[DEBUG]   > Minimizer at pos {}: Spawning Super-k-mer for partition {}, len={}, seq[{}..{}], Preview={}...", current_pos, partition_id, superkmer_slice.len(), start_pos, end_pos, String::from_utf8_lossy(sk_preview));

                let buffer = local_buffers.entry(partition_id).or_insert_with(Vec::new);
                write!(buffer, ">{}\n", file_id).unwrap();
                buffer.extend_from_slice(superkmer_slice);
                buffer.extend_from_slice(b"\n");
            }
        }
    }
    println!("[Thread {:?}] Flushing final buffers for file {}", thread::current().id(), file_id);
    let mut locked_files = partition_files.lock();
    for (partition_id, buffer) in local_buffers.iter_mut().filter(|(_, b)| !b.is_empty()) {
        if let Some(encoder) = locked_files.get_mut(partition_id) {
            if encoder.write_all(buffer).is_err() {
                 eprintln!("Error during final flush to partition {}", partition_id);
            }
        }
        buffer.clear();
    }
    Ok(())
}


pub fn test_parser(input_fof: &String, output_dir: &String, k: usize, m: usize) -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    let threads = num_cpus::get();
    
    let num_partitions = 4u64.pow(m as u32);
    let required_limit = num_partitions + 50;

    let (soft_limit, hard_limit) = rlimit::getrlimit(rlimit::Resource::NOFILE)?;
    if soft_limit < required_limit {
        println!("Current open file limit (soft={}, hard={}) is less than required ({})", soft_limit, hard_limit, required_limit);
        let new_soft_limit = std::cmp::min(required_limit, hard_limit);
        println!("Attempting to raise soft limit to {}...", new_soft_limit);
        if let Err(e) = rlimit::setrlimit(rlimit::Resource::NOFILE, new_soft_limit, hard_limit) {
            eprintln!("Failed to raise file limit: {}. Please run `ulimit -n {}` manually before running this program.", e, required_limit);
            return Err(e.into());
        }
        let (new_soft, _) = rlimit::getrlimit(rlimit::Resource::NOFILE)?;
        if new_soft < required_limit {
             eprintln!("FATAL: Could not raise open file limit high enough. Current limit is {}, but {} is required.", new_soft, required_limit);
             eprintln!("Try running `sudo sysctl -w fs.file-max=100000` or raising the hard limit in /etc/security/limits.conf");
             return Err(io::Error::new(io::ErrorKind::Other, "Failed to set rlimit"));
        }
        println!("Successfully raised open file limit to {}", new_soft);
    }

    println!("Starting super-k-mer partitioning with k={}, m={}, threads={}", k, m, threads);
    fs::create_dir_all(output_dir)?;

    println!("Pre-creating {} partition files...", num_partitions);
    let mut partition_files = HashMap::new();
    for i in 0..num_partitions {
        let file_path = PathBuf::from(output_dir).join(format!("{}.fa.zst", i));
        let file = match File::create(&file_path) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("FATAL: Could not create partition file {:?}: {}", file_path, e);
                return Err(e);
            }
        };
        let encoder = ZstdEncoder::new(file, 3)?;
        partition_files.insert(i, encoder);
    }
    println!("All {} partition files pre-created successfully.", num_partitions);
    let partition_files_arc = Arc::new(Mutex::new(partition_files));
    let file = File::open(input_fof)?;
    let reader = BufReader::new(file);
    let file_paths: Vec<PathBuf> = reader
        .lines()
        .filter_map(Result::ok)
        .filter(|line| !line.trim().is_empty())
        .map(PathBuf::from)
        .collect();
    println!("Found {} files to process.", file_paths.len());
    let pool = threadpool::ThreadPool::new(threads);
    for (i, file_path) in file_paths.into_iter().enumerate() {
        let partition_files_clone = Arc::clone(&partition_files_arc);
        let file_id = i + 1;
        pool.execute(move || {
            if let Err(e) = process_file(&file_path, file_id, &partition_files_clone, k, m) {
                eprintln!("Failed to process file {}: {}", file_path.display(), e);
            }
        });
    }
    pool.join();
    println!("Finalizing all partition files...");
    let mut final_files = partition_files_arc.lock();
    for (id, encoder) in final_files.drain() {
        if let Err(e) = encoder.finish() {
            eprintln!("Error finalizing partition file {}: {}", id, e);
        }
    }
    println!("Partitioning complete.");
    Ok(())
}