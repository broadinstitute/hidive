//! Hidive is a targeted genome co-assembler for biobank-scale long-read and short-read data.
//! It is designed to assemble high-diversity loci from long-read data and to co-assemble these
//! loci with short-read data. Hidive is built on top of the Skydive library, which provides
//! a series-parallel graph representation of long-read data and a linked de Bruijn graph
//! representation of short-read data.
//!
//! # Quick Start
//!
//! The following commands will compile the Rust and Python codebase.
//!
//! ## Install Rust
//!
//! ```sh
//! curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
//! ```
//!
//! ## Download and Build Hidive
//!
//! ```sh
//! git clone https://github.com/broadinstitute/hidive.git
//! cd hidive
//! cargo build --release
//! ```
//!
//! ## Configure Python Environment
//!
//! ```sh
//! python -m venv venv
//! . venv/bin/activate
//! pip install -r dev-requirements.txt
//! ```
//!
//! ## Build Hidive's Python Codebase, Pydive
//!
//! ```sh
//! cd src/pydive/
//! maturin develop --release
//! ```
//! # Prerequisites
//! Hidive is designed to access local files or data in Google Cloud Storage (GCS).
//! Within certain cloud-computing environments (i.e. Terra, All of Us Researcher Workbench),
//! access to GCS is already configured. For accessing files in GCS on your local machine,
//! you will also need to install the [Google Cloud CLI](https://cloud.google.com/sdk/docs/install-sdk).
//! Then, configure your [Application Default Credentials (ADC)](https://cloud.google.com/docs/authentication/provide-credentials-adc#local-dev).
//! If accessing [requester pays buckets](https://cloud.google.com/storage/docs/requester-pays)
//! set the following environment variable before running hidive commands:
//! ```sh
//!  export GCS_REQUESTER_PAYS_PROJECT=<Google Project ID>
//! ```

use std::path::PathBuf;

use clap::{Parser, Subcommand};

mod assemble;
mod build;
mod cluster;
mod coassemble;
mod fetch;
mod filter;
mod impute;
mod rescue;
mod train;
mod trim;

#[derive(Debug, Parser)]
#[clap(name = "hidive")]
#[clap(about = "Analysis of high-diversity loci through genome co-assembly of long/short reads.", long_about = None)]
// #[clap(author = "Kiran V Garimella (kiran@broadinstitute.org)")]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

const DEFAULT_KMER_SIZE: usize = 17;

#[derive(Debug, Subcommand)]
enum Commands {
    /// Train a graph-cleaning model using short- and long-read data with ground truth assemblies.
    #[clap(arg_required_else_help = true)]
    Train {
        /// Output path for trained model.
        #[clap(short, long, value_parser, default_value = "/dev/stdout")]
        output: PathBuf,

        /// Kmer-size.
        #[clap(short, long, value_parser, default_value_t = DEFAULT_KMER_SIZE)]
        kmer_size: usize,

        /// Number of training iterations.
        #[clap(short, long, value_parser, default_value_t = 50)]
        iterations: usize,

        /// Indexed WGS BAM, CRAM, or FASTA files from which to extract relevant sequences.
        #[clap(short, long, value_parser, required = true)]
        long_read_seq_paths: Vec<PathBuf>,

        /// Indexed WGS BAM, CRAM, or FASTA files from which to extract relevant sequences.
        #[clap(short, long, value_parser, required = true)]
        short_read_seq_paths: Vec<PathBuf>,

        /// Indexed BAM files to use as ground truth (usually from ultra-high-quality assemblies).
        #[clap(value_parser, required = true)]
        truth_seq_paths: Vec<PathBuf>,

        /// Turn on debug mode.
        #[clap(short, long, value_parser)]
        debug: bool,
    },

    /// Stream selected loci from FASTA and long-read WGS BAM files stored locally or in Google Cloud Storage.
    #[clap(arg_required_else_help = true)]
    Fetch {
        /// Output path for FASTA file with reads spanning locus of interest.
        #[clap(short, long, value_parser, default_value = "/dev/stdout")]
        output: PathBuf,

        /// Zero or more genomic loci ("contig:start-stop") to extract from WGS BAM files.
        #[clap(short, long, value_parser, required = false)]
        loci: Vec<String>,

        /// Include unmapped reads.
        #[clap(short, long, value_parser)]
        unmapped: bool,

        /// Indexed WGS BAM, CRAM, or FASTA files from which to extract relevant sequences.
        #[clap(required = true, value_parser)]
        seq_paths: Vec<PathBuf>,
    },

    /// Find more sequences (aligned or unaligned) overlapping previously fetched reads.
    #[clap(arg_required_else_help = true)]
    Rescue {
        /// Output path for FASTA file with reads spanning locus of interest.
        #[clap(short, long, value_parser, default_value = "/dev/stdout")]
        output: PathBuf,

        /// Kmer-size.
        #[clap(short, long, value_parser, default_value_t = DEFAULT_KMER_SIZE)]
        kmer_size: usize,

        /// Minimum number of k-mers to require before examining a read more carefully.
        #[clap(short, long, value_parser, default_value_t = 10)]
        min_kmers: usize,

        /// FASTA files with reads to use as a filter for finding more reads.
        #[clap(short, long, value_parser, required = true)]
        fasta_paths: Vec<PathBuf>,

        /// Indexed WGS BAM, CRAM, or FASTA files from which to extract relevant sequences.
        #[clap(required = true, value_parser)]
        seq_paths: Vec<PathBuf>,
    },

    /// Filter rescued reads to those most closely matching the long-read data.
    #[clap(arg_required_else_help = true)]
    Filter {
        /// Output path for filtered short-read sequences.
        #[clap(short, long, value_parser, default_value = "/dev/stdout")]
        output: PathBuf,

        /// Kmer-size
        #[clap(short, long, value_parser, default_value_t = DEFAULT_KMER_SIZE)]
        kmer_size: usize,

        /// Minimum percentage of bases in short-read sequences covered by the matched kmers.
        #[clap(short, long, value_parser, default_value_t = 90)]
        min_score_pct: usize,

        /// FASTA files with short-read sequences (may contain one or more samples).
        #[clap(required = true, value_parser)]
        short_read_fasta_paths: Vec<PathBuf>,

        /// FASTA files with long-read sequences (may contain one or more samples).
        #[clap(short, long, required = true, value_parser)]
        long_read_fasta_paths: Vec<PathBuf>,
    },

    /// Cluster sequences based on k-mer presence/absence.
    #[clap(arg_required_else_help = true)]
    Cluster {
        /// Output path for clustered sequences.
        #[clap(short, long, value_parser, default_value = "/dev/stdout")]
        output: PathBuf,

        /// Kmer-size
        #[clap(short, long, value_parser, default_value_t = DEFAULT_KMER_SIZE)]
        kmer_size: usize,

        /// Multi-sample FASTA file with reads spanning locus of interest.
        #[clap(required = true, value_parser)]
        fasta_path: PathBuf,
    },

    /// Trim reads to a specific window around locus.
    #[clap(arg_required_else_help = true)]
    Trim {
        /// Output path for trimmed BAM.
        #[clap(short, long, value_parser, default_value = "/dev/stdout")]
        output: PathBuf,

        /// One or more genomic loci ("contig:start-stop") to extract from WGS BAM files.
        #[clap(short, long, value_parser)]
        loci: Vec<String>,

        /// Multi-sample BAM file with reads spanning locus of interest.
        #[clap(required = true, value_parser)]
        bam_path: PathBuf,
    },

    /// Build series-parallel graph from long-read data in multi-sample FASTA file.
    #[clap(arg_required_else_help = true)]
    Build {
        /// Output path for series-parallel graph.
        #[clap(short, long, value_parser, default_value = "/dev/stdout")]
        output: PathBuf,

        /// Kmer-size
        #[clap(short, long, value_parser, default_value_t = DEFAULT_KMER_SIZE)]
        kmer_size: usize,

        /// Name of sequence to use as reference.
        #[clap(short, long, value_parser, required = true)]
        reference_name: String,

        /// Multi-sample FASTA file with reads spanning locus of interest.
        #[clap(required = true, value_parser)]
        fasta_path: PathBuf,
    },

    /// Cluster edge matrix and impute missing edges.
    #[clap(arg_required_else_help = true)]
    Impute {
        /// Output path for series-parallel graph with imputed missing edges.
        #[clap(short, long, value_parser, default_value = "/dev/stdout")]
        output: PathBuf,

        /// Series-parallel graph.
        #[clap(required = true, value_parser)]
        graph: PathBuf,
    },

    /// Assemble target locus from long-read data in anchor-based series-parallel graph.
    #[clap(arg_required_else_help = true)]
    Assemble {
        /// Output path for assembled long-read sequences.
        #[clap(short, long, value_parser, default_value = "/dev/stdout")]
        output: PathBuf,

        /// Series-parallel graph.
        #[clap(required = true, value_parser)]
        graph: PathBuf,
    },

    /// Co-assemble target locus from long-read and short-read data using a linked de Bruijn graph.
    #[clap(arg_required_else_help = true)]
    Coassemble {
        /// Output path for assembled short-read sequences.
        #[clap(short, long, value_parser, default_value = "/dev/stdout")]
        output: PathBuf,

        /// Kmer-size
        #[clap(short, long, value_parser, default_value_t = DEFAULT_KMER_SIZE)]
        kmer_size: usize,

        /// Trained error-cleaning model.
        #[clap(short, long, required = true, value_parser)]
        model_path: PathBuf,

        /// FASTA files with short-read sequences (may contain one or more samples).
        #[clap(short, long, required = false, value_parser)]
        short_read_fasta_paths: Vec<PathBuf>,

        /// FASTA files with long-read sequences (may contain one or more samples).
        #[clap(required = true, value_parser)]
        long_read_fasta_paths: Vec<PathBuf>,
    },
}

fn main() {
    let args = Cli::parse();

    skydive::elog!("Hidive version {}", env!("CARGO_PKG_VERSION"));
    skydive::elog!("{:?}", args);

    let start_time = std::time::Instant::now();

    match args.command {
        Commands::Train {
            output,
            kmer_size,
            iterations,
            long_read_seq_paths,
            short_read_seq_paths,
            truth_seq_paths,
            debug,
        } => {
            train::start(
                &output,
                kmer_size,
                iterations,
                &long_read_seq_paths,
                &short_read_seq_paths,
                &truth_seq_paths,
                debug,
            );
        }
        Commands::Fetch {
            output,
            loci,
            unmapped,
            seq_paths,
        } => {
            fetch::start(&output, &loci, unmapped, &seq_paths);
        }
        Commands::Rescue {
            output,
            kmer_size,
            min_kmers,
            fasta_paths,
            seq_paths,
        } => {
            rescue::start(&output, kmer_size, min_kmers, &fasta_paths, &seq_paths);
        }
        Commands::Filter {
            output,
            kmer_size,
            min_score_pct,
            long_read_fasta_paths,
            short_read_fasta_paths,
        } => {
            filter::start(
                &output,
                kmer_size,
                min_score_pct,
                &long_read_fasta_paths,
                &short_read_fasta_paths,
            );
        }
        Commands::Cluster {
            output,
            kmer_size,
            fasta_path,
        } => {
            cluster::start(&output, kmer_size, &fasta_path);
        }
        Commands::Trim {
            output,
            loci,
            bam_path,
        } => {
            trim::start(&output, &loci, &bam_path);
        }
        Commands::Build {
            output,
            kmer_size,
            fasta_path,
            reference_name,
        } => {
            build::start(&output, kmer_size, &fasta_path, reference_name);
        }
        Commands::Impute { output, graph } => {
            impute::start(&output, &graph);
        }
        Commands::Assemble { output, graph } => {
            assemble::start(&output, &graph);
        }
        Commands::Coassemble {
            output,
            model_path,
            kmer_size,
            long_read_fasta_paths,
            short_read_fasta_paths,
        } => {
            coassemble::start(
                &output,
                kmer_size,
                &model_path,
                &long_read_fasta_paths,
                &short_read_fasta_paths,
            );
        }
    }

    skydive::elog!("Complete. Elapsed time: {}.", elapsed_time(start_time));
}

fn elapsed_time(start_time: std::time::Instant) -> String {
    let end_time = std::time::Instant::now();
    let elapsed_time = end_time.duration_since(start_time);

    let elapsed_secs = elapsed_time.as_secs_f64();
    let elapsed_str = if elapsed_secs < 60.0 {
        format!("{:.2} seconds", elapsed_secs)
    } else if elapsed_secs < 3600.0 {
        format!("{:.2} minutes", elapsed_secs / 60.0)
    } else if elapsed_secs < 86400.0 {
        format!("{:.2} hours", elapsed_secs / 3600.0)
    } else {
        format!("{:.2} days", elapsed_secs / 86400.0)
    };

    elapsed_str
}
