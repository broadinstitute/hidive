use std::{fs::File, path::PathBuf, io::Write};

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
        .clean(0.2, 0.01)
        .build_links(&all_lr_seqs, true)
        ;

    skydive::elog!("Built MLdBG with {} k-mers.", m.kmers.len());

    // let contigs = m.assemble_at_bubbles();

    // let mut fa_file = File::create(&output).unwrap();
    // for (i, contig) in contigs.iter().enumerate() {
    //     let _ = writeln!(fa_file, ">contig_{}\n{}", i, String::from_utf8(contig.clone()).unwrap());
    // }

    let progress_bar = skydive::utils::default_bounded_progress_bar("Correcting reads", all_lr_seqs.len() as u64);

    let corrected_seqs = all_lr_seqs
        .par_iter()
        .progress_with(progress_bar)
        .map(|seq| m.correct_seq(seq))
        .flatten()
        .collect::<Vec<Vec<u8>>>();

    skydive::elog!("Writing reads to {}", output.display());
    
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

        writeln!(csv_file, "node,label").unwrap();

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
}