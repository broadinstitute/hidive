use std::{fs::File, path::PathBuf, io::Write};

use rayon::prelude::*;
use indicatif::ParallelProgressIterator;

use skydive::ldbg::LdBG;
use skydive::mldbg::MLdBG;
use skydive::utils::*;

pub fn start(
    output: &PathBuf,
    kmer_size: usize,
    model_path: &PathBuf,
    long_read_fasta_paths: &Vec<PathBuf>,
    short_read_fasta_paths: &Vec<PathBuf>,
) {
    let long_read_seq_urls = skydive::parse::parse_file_names(long_read_fasta_paths);
    let short_read_seq_urls = skydive::parse::parse_file_names(short_read_fasta_paths);

    // Read all long reads.
    skydive::elog!("Processing long-read samples {:?}...", long_read_seq_urls.iter().map(|url| url.as_str()).collect::<Vec<&str>>());
    let all_lr_seqs = skydive::utils::read_fasta(long_read_fasta_paths);

    // Read all short reads.
    skydive::elog!("Processing short-read samples {:?}...", short_read_seq_urls.iter().map(|url| url.as_str()).collect::<Vec<&str>>());
    let all_sr_seqs = skydive::utils::read_fasta(short_read_fasta_paths);

    let l1 = LdBG::from_sequences("lr".to_string(), kmer_size, &all_lr_seqs);
    let s1 = LdBG::from_sequences("sr".to_string(), kmer_size, &all_sr_seqs);

    let m = MLdBG::from_ldbgs(vec![l1, s1])
        .score_kmers(model_path)
        .collapse()
        .clean_color_specific_paths(1, 0.2)
        .clean_branches(0.01)
        .clean_tips(3*kmer_size, 0.01)
        .clean_contigs(100)
        .build_links(&all_lr_seqs);

    // let superbubbles = m.identify_superbubbles();
    // for superbubble in superbubbles {
    //     skydive::elog!("Superbubble: {:?}", superbubble);
    // }

    let progress_bar = skydive::utils::default_bounded_progress_bar("Correcting reads", all_lr_seqs.len() as u64);

    // let read = b"GAACAGCTCAGGTGGAGAAGGGGTGAAGGGTGGGGTCTGAGATTTCTTGTCTCACTGAGGGTTCCAAGGCCCCAGCTAGAAATGTGCCCTGTCTCATTACTGGGAAGCACCATCCACAATCATGGGCCGACCCAGCCTGGGCCCTGTGTGCCAGCACTTACTCTTTTGTAAAGCACCTGAGCAATGAAGGACAGATTTATCACCTTGATTATGGCGGTGATGGGACCTGATCCCAGCAGTCACAAGTCACAGGGGAAGGTCCCTGACGACAGATCTCAGGAGGGCGATTGGTCCAGGGCCCACATCTGCTTTCTTCATGTTTCCTGATCCTGCCCTGGGTCTGCAGTCACACATTTCTGGAAACTTCTCTGGGGTCCAAGACTAGGAGGTTCCTCTAGGACCTTAAGGCCCTGGCTCCTTTCTGTATCTCACAGGACATTTTCTTCCCACAGATAGAAAAGGAGGGAGCTACTCTCAGGCTGCAAGTAAGTATGAAGGAGGCTGATGCCTGAGGTCCTTGGGATATTGTGTTTGGGAGCCCATGGGGGAGCT";

    let corrected_seqs =
        // vec![read.to_vec()]
        // .iter()
        all_lr_seqs
        .par_iter()
        .progress_with(progress_bar)
        .map(|seq| m.correct_seq(seq))
        .flatten()
        .collect::<Vec<Vec<u8>>>();

    for (i, corrected_seq) in corrected_seqs.iter().enumerate() {
        println!(">corrected_{}\n{}", i, String::from_utf8(corrected_seq.clone()).unwrap());
    }

    let g = m.traverse_all_kmers();
    // let g = m.traverse_all_contigs();

    let _ = write_gfa(&mut File::create(output).unwrap(), &g);

    let csv_output = output.with_extension("csv");
    let mut csv_file = File::create(&csv_output).unwrap();

    for (node_index, node_label) in g.node_indices().zip(g.node_weights()) {
        let kmer = node_label.as_bytes();
        let cn_kmer = skydive::utils::canonicalize_kmer(kmer);
        let score = (100.0 * *m.scores.get(&cn_kmer).unwrap()) as u32;
        let sources = m.sources.get(&cn_kmer).unwrap();

        let source = if sources.len() == 1 { sources[0] } else { 2 };

        writeln!(
            csv_file,
            "{},{}",
            node_index.index(),
            format!("{} ({})", source, score)
        )
        .unwrap();
    }
}