use std::path::PathBuf;

use clap::{Parser, Subcommand};

mod assemble;
mod build;
mod cluster;
mod coassemble;
mod fetch;
mod impute;
mod trim;

#[derive(Debug, Parser)] // requires `derive` feature
#[clap(name = "hidive")]
#[clap(about = "Analysis of high-diversity loci through genome co-assembly of long/short reads.", long_about = None)]
// #[clap(author = "Kiran V Garimella (kiran@broadinstitute.org)")]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Stream selected loci from FASTA and long-read WGS BAM files stored locally or in Google Cloud Storage.
    #[clap(arg_required_else_help = true)]
    Fetch {
        /// Output path for multi-sample BAM file with reads spanning locus of interest.
        #[clap(short, long, value_parser, default_value = "/dev/stdout")]
        output: PathBuf,

        /// One or more genomic loci ("contig:start-stop") to extract from WGS BAM files.
        #[clap(short, long, value_parser)]
        loci: Vec<String>,

        /// Indexed WGS BAM, CRAM, or FASTA files from which to extract relevant sequences.
        #[clap(required = true, value_parser)]
        seq_paths: Vec<PathBuf>,
    },

    /// Cluster sequences based on k-mer presence/absence.
    #[clap(arg_required_else_help = true)]
    Cluster {
        /// Output path for clustered sequences.
        #[clap(short, long, value_parser, default_value = "/dev/stdout")]
        output: PathBuf,

        /// Kmer-size
        #[clap(short, long, value_parser, default_value = "11")]
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
        #[clap(short, long, value_parser, default_value = "11")]
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

    match args.command {
        Commands::Fetch {
            output,
            loci,
            seq_paths,
        } => {
            fetch::start(&output, &loci, &seq_paths);
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
            long_read_fasta_paths,
            short_read_fasta_paths,
        } => {
            coassemble::start(&output, &long_read_fasta_paths, &short_read_fasta_paths);
        }
    }

    skydive::elog!("Complete.");
}
