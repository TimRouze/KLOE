use std::path::PathBuf;

use anyhow::Result;
use clap::{Parser, Subcommand};

const EXPLORATORY_WARNING: &str =
    "WARNING: KLOER is an exploratory project. Its archive format and command-line interface may change.";

#[derive(Parser, Debug)]
#[command(author, version, about, before_help = EXPLORATORY_WARNING)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Map unitigs to a reference and encode coordinates plus divergence.
    #[command(before_help = EXPLORATORY_WARNING)]
    Compress {
        #[arg(short, long)]
        reference: PathBuf,
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(short, long, default_value_t = std::thread::available_parallelism().map_or(1, usize::from))]
        threads: usize,
        #[arg(long, default_value_t = 0.20)]
        max_divergence: f64,
        #[arg(long, default_value_t = 0.80)]
        min_coverage: f64,
        #[arg(long, default_value_t = 512)]
        bandwidth: i32,
        #[arg(long, default_value_t = 3)]
        compression_level: i32,
    },
    /// Compress a file-of-files into one multi-dataset reference archive.
    #[command(before_help = EXPLORATORY_WARNING)]
    CompressList {
        #[arg(short, long)]
        reference: PathBuf,
        #[arg(short = 'f', long)]
        input_list: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(short, long, default_value_t = std::thread::available_parallelism().map_or(1, usize::from))]
        threads: usize,
        #[arg(long, default_value_t = 0.20)]
        max_divergence: f64,
        #[arg(long, default_value_t = 0.80)]
        min_coverage: f64,
        #[arg(long, default_value_t = 512)]
        bandwidth: i32,
        #[arg(long, default_value_t = 9)]
        compression_level: i32,
    },
    /// Recover the exact unitig sequences using the same reference.
    #[command(before_help = EXPLORATORY_WARNING)]
    Decompress {
        #[arg(short, long)]
        reference: PathBuf,
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
    },
    /// Recover every dataset from a multi-dataset archive.
    #[command(before_help = EXPLORATORY_WARNING)]
    DecompressList {
        #[arg(short, long)]
        reference: PathBuf,
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        output_dir: PathBuf,
    },
}

fn main() -> Result<()> {
    let cli = Cli::try_parse().unwrap_or_else(|error| {
        if error.use_stderr() && !error.to_string().contains(EXPLORATORY_WARNING) {
            eprintln!("{EXPLORATORY_WARNING}\n");
        }
        error.exit()
    });
    match cli.command {
        Command::Compress {
            reference,
            input,
            output,
            threads,
            max_divergence,
            min_coverage,
            bandwidth,
            compression_level,
        } => {
            let options = kloer::CompressOptions {
                threads,
                max_divergence,
                min_coverage,
                bandwidth,
                compression_level,
            };
            let stats = kloer::compress(&reference, &input, &output, &options)?;
            println!(
                "KLOER compression complete: records={}, mapped={}, literal={}, input_bases={}, mapped_bases={}, literal_bases={}, divergence_bases={}, unique_records={}, unique_mapped={}, unique_literal={}, unique_bases={}, unique_mapped_bases={}, unique_literal_bases={}",
                stats.records,
                stats.mapped,
                stats.literal,
                stats.input_bases,
                stats.mapped_bases,
                stats.literal_bases,
                stats.divergence_bases,
                stats.unique_records,
                stats.unique_mapped,
                stats.unique_literal,
                stats.unique_bases,
                stats.unique_mapped_bases,
                stats.unique_literal_bases
            );
        }
        Command::CompressList {
            reference,
            input_list,
            output,
            threads,
            max_divergence,
            min_coverage,
            bandwidth,
            compression_level,
        } => {
            let options = kloer::CompressOptions {
                threads,
                max_divergence,
                min_coverage,
                bandwidth,
                compression_level,
            };
            let stats = kloer::compress_list(&reference, &input_list, &output, &options)?;
            println!(
                "KLOER multi-dataset compression complete: records={}, mapped={}, literal={}, input_bases={}, mapped_bases={}, literal_bases={}, divergence_bases={}, unique_records={}, unique_mapped={}, unique_literal={}, unique_bases={}, unique_mapped_bases={}, unique_literal_bases={}",
                stats.records,
                stats.mapped,
                stats.literal,
                stats.input_bases,
                stats.mapped_bases,
                stats.literal_bases,
                stats.divergence_bases,
                stats.unique_records,
                stats.unique_mapped,
                stats.unique_literal,
                stats.unique_bases,
                stats.unique_mapped_bases,
                stats.unique_literal_bases
            );
        }
        Command::Decompress {
            reference,
            input,
            output,
        } => {
            let stats = kloer::decompress(&reference, &input, &output)?;
            println!(
                "KLOER decompression complete: records={}, bases={}",
                stats.records, stats.bases
            );
        }
        Command::DecompressList {
            reference,
            input,
            output_dir,
        } => {
            let stats = kloer::decompress_list(&reference, &input, &output_dir)?;
            println!(
                "KLOER multi-dataset decompression complete: records={}, bases={}",
                stats.records, stats.bases
            );
        }
    }
    Ok(())
}
