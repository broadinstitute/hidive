use std::collections::{HashMap, HashSet};
use std::{fs::File, path::PathBuf, io::Write};

use needletail::Sequence;
use petgraph::graph::NodeIndex;
use rayon::prelude::*;
use indicatif::ParallelProgressIterator;

use skydive::ldbg::LdBG;
use skydive::mldbg::MLdBG;
use skydive::utils::*;

pub fn start(
    output: &PathBuf,
    gfa_output: Option<PathBuf>,
    kmer_size: usize,
    model_path: &PathBuf,
    long_read_fasta_path: &PathBuf,
    short_read_fasta_path: &PathBuf,
) {
    let long_read_seq_urls = skydive::parse::parse_file_names(&[long_read_fasta_path.clone()]);
    let short_read_seq_urls = skydive::parse::parse_file_names(&[short_read_fasta_path.clone()]);

    // Read all long reads.
    skydive::elog!("Processing long-read samples {:?}...", long_read_seq_urls.iter().map(|url| url.as_str()).collect::<Vec<&str>>());
    let all_lr_seqs = skydive::utils::read_fasta(&vec![long_read_fasta_path.clone()]);

    // Read all short reads.
    skydive::elog!("Processing short-read samples {:?}...", short_read_seq_urls.iter().map(|url| url.as_str()).collect::<Vec<&str>>());
    let all_sr_seqs = skydive::utils::read_fasta(&vec![short_read_fasta_path.clone()]);

    let l1 = LdBG::from_sequences("lr".to_string(), kmer_size, &all_lr_seqs);
    let s1 = LdBG::from_sequences("sr".to_string(), kmer_size, &all_sr_seqs);

    let m = MLdBG::from_ldbgs(vec![l1, s1])
        .score_kmers(model_path)
        .collapse()
        .clean(0.1)
        .build_links(&all_lr_seqs, false);

    skydive::elog!("Built MLdBG with {} k-mers.", m.kmers.len());

    skydive::elog!("Correcting reads...");
    let corrected_seqs = m.correct_seqs(&all_lr_seqs);

    let mut fa_file = File::create(output).unwrap();
    for (i, corrected_seq) in corrected_seqs.iter().enumerate() {
        let _ = writeln!(fa_file, ">corrected_{}\n{}", i, String::from_utf8(corrected_seq.clone()).unwrap());
    }

    if let Some(gfa_output) = gfa_output {
        skydive::elog!("Writing GFA to {}", gfa_output.display());

        let g = m.traverse_all_kmers();

        let _ = write_gfa(&mut File::create(gfa_output.clone()).unwrap(), &g);

        let csv_output = gfa_output.with_extension("csv");
        let mut csv_file = File::create(&csv_output).unwrap();

        writeln!(csv_file, "node,label,kmer,cov,entropy").unwrap();

        for (node_index, node_label) in g.node_indices().zip(g.node_weights()) {
            let kmer = node_label.as_bytes();
            let cn_kmer = skydive::utils::canonicalize_kmer(kmer);
            let score = (100.0 * *m.scores.get(&cn_kmer).unwrap_or(&0.0)) as u32;
            let cov = if m.kmers.contains_key(&cn_kmer) { m.kmers.get(&cn_kmer).unwrap().coverage() } else { 0 };
            let entropy = skydive::utils::shannon_entropy(kmer);
            let sources = m.sources.get(&cn_kmer).unwrap_or(&vec![]).clone();

            let source = if sources.len() == 1 { sources[0] } else { 2 };

            writeln!(
                csv_file,
                "{},{},{},{},{}",
                node_index.index(),
                format!("{} ({})", source, score),
                node_label,
                cov,
                entropy,
            )
            .unwrap();
        }
    }
}