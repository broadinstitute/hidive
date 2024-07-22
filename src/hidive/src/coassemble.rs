use std::path::PathBuf;
use std::{fs::File, io::{BufWriter, Write}};

use rayon::iter::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;

use skydive;
use skydive::ldbg::LdBG;
use skydive::mldbg::MLdBG;

pub fn start(
    output: &PathBuf,
    kmer_size: usize,
    long_read_fasta_paths: &Vec<PathBuf>,
    short_read_fasta_paths: &Vec<PathBuf>,
) {
    let long_read_seq_urls = skydive::parse::parse_file_names(long_read_fasta_paths);
    let short_read_seq_urls = skydive::parse::parse_file_names(short_read_fasta_paths);

    // Create a multi-color linked de Bruijn graph
    let mut g = MLdBG::new(kmer_size, true);

    // Construct a graph from long reads, and add it to g.
    for long_read_seq_url in long_read_seq_urls {
        let basename = skydive::utils::basename_without_extension(&long_read_seq_url, &[".fasta.gz", ".fa.gz", ".fasta", ".fa"]);
        let fasta_path = long_read_seq_url.to_file_path().unwrap();

        skydive::elog!("Processing long read sample {}...", basename);
        g.append_from_file(basename, &fasta_path);
    }

    // Construct a graph from the short read data, filtered for k-mer overlap with the existing graph, and add it to g as a separate color.
    for short_read_seq_url in short_read_seq_urls {
        let basename = skydive::utils::basename_without_extension(&short_read_seq_url, &[".fasta.gz", ".fa.gz", ".fasta", ".fa"]);
        let fasta_path = short_read_seq_url.to_file_path().unwrap();

        skydive::elog!("Processing short read sample {}...", basename);
        g.append_from_filtered_file(basename, &fasta_path, |r, kmer_union| {
            let num_contained = r.seq().windows(kmer_size)
                .filter(|kmer| kmer_union.contains(&skydive::ldbg::LdBG::canonicalize_kmer(kmer)))
                .count();
            num_contained > 3
        });
    }
}
