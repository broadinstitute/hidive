use std::fs::File;
use std::path::PathBuf;

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
        .clean_paths(0.5)
        .clean_tips(3*kmer_size)
        .clean_contigs(200);

    let _ = write_gfa(&mut File::create(output).unwrap(), &m.traverse_all_kmers());
}