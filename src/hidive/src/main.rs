use std::path::PathBuf;

use clap::{Parser, Subcommand};

mod prepare;
mod assemble;
mod join;

#[derive(Debug, Parser)] // requires `derive` feature
#[clap(name = "hidive")]
#[clap(about = "Analysis of high-diversity loci through genome co-assembly of long/short reads.", long_about = None)]
#[clap(author = "Kiran V Garimella (kiran@broadinstitute.org)")]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Extract target haplotypes from FASTA file(s).
    #[clap(arg_required_else_help = true)]
    Prepare {
        /// Output graph path.
        #[clap(short, long, value_parser)]
        output: PathBuf,
        
        /// Loci to extract from FASTA files.
        #[clap(short, long, value_parser)]
        locus: Option<Vec<String>>,

        /// Gene features in GFF2/GFF3/GTF format.
        #[clap(short, long, value_parser)]
        gff: Option<PathBuf>,

        /// Indexed FASTA files from which to extract haplotypes.
        #[clap(required = true, value_parser)]
        fasta: PathBuf,
    },
    /// Join two or more hidive graphs into a single graph.
    #[clap(arg_required_else_help = true)]
    Join {
        /// Output graph path.
        #[clap(short, long, value_parser)]
        output: PathBuf,

        /// Hidive graphs to combine.
        #[clap(required = true, value_parser)]
        graph: Vec<PathBuf>,
    },
    /// Assemble target locus from long/short-read data in a BAM/CRAM file.
    #[clap(arg_required_else_help = true)]
    Assemble {
        /// Output graph path.
        #[clap(short, long, value_parser)]
        output: PathBuf,
        
        /// Loci to extract from BAM/CRAM files.
        #[clap(short, long, value_parser)]
        locus: Option<Vec<String>>,

        /// Aligned reads in BAM/CRAM format.
        #[clap(required = true, value_parser)]
        reads: PathBuf,
    },
}

fn main() {
    let args = Cli::parse();

    match args.command {
        Commands::Prepare { output, locus, gff, fasta } => {
            prepare::start(output, &locus, gff, fasta);
        }
        Commands::Join { output, graph } => {
            join::start(output, graph);
        }
        Commands::Assemble { output, locus, reads } => {
            assemble::start(output, locus, reads);
        }
    }
}