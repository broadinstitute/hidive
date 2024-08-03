use std::collections::HashSet;
use std::path::PathBuf;
use std::{
    fs::File,
    io::{BufWriter, Write},
};

use bio::alignment::sparse::*;

use skydive::ldbg::LdBG;

pub fn start(
    output: &PathBuf,
    kmer_size: usize,
    min_score_pct: usize,
    long_read_fasta_paths: &Vec<PathBuf>,
    short_read_fasta_paths: &Vec<PathBuf>,
) {
    let long_read_seq_urls = skydive::parse::parse_file_names(long_read_fasta_paths);
    let short_read_seq_urls = skydive::parse::parse_file_names(short_read_fasta_paths);

    // Read all long reads.
    let mut all_lr_seqs: Vec<Vec<u8>> = Vec::new();
    for long_read_seq_url in &long_read_seq_urls {
        let basename = skydive::utils::basename_without_extension(
            &long_read_seq_url,
            &[".fasta.gz", ".fa.gz", ".fasta", ".fa"],
        );
        let fasta_path = long_read_seq_url.to_file_path().unwrap();

        skydive::elog!("Processing long-read sample {}...", basename);

        let reader = bio::io::fasta::Reader::from_file(&fasta_path).unwrap();
        let all_reads: Vec<bio::io::fasta::Record> = reader.records().map(|r| r.unwrap()).collect();

        all_lr_seqs.extend(
            all_reads
                .iter()
                .map(|r| r.seq().to_vec())
                .collect::<Vec<Vec<u8>>>(),
        );
    }

    // Read all short reads.
    let mut all_sr_seqs: Vec<Vec<u8>> = Vec::new();
    for short_read_seq_url in &short_read_seq_urls {
        let basename = skydive::utils::basename_without_extension(
            &short_read_seq_url,
            &[".fasta.gz", ".fa.gz", ".fasta", ".fa"],
        );
        let fasta_path = short_read_seq_url.to_file_path().unwrap();

        skydive::elog!("Processing short-read sample {}...", basename);

        let reader = bio::io::fasta::Reader::from_file(&fasta_path).unwrap();
        let all_reads: Vec<bio::io::fasta::Record> = reader.records().map(|r| r.unwrap()).collect();

        all_sr_seqs.extend(
            all_reads
                .iter()
                .map(|r| r.seq().to_vec())
                .collect::<Vec<Vec<u8>>>(),
        );
    }

    skydive::elog!("{} long-read sequences.", all_lr_seqs.len());
    skydive::elog!("{} short-read sequences.", all_sr_seqs.len());

    let filtered_sr_seqs = all_sr_seqs.iter().filter(|sr_seq| {
        let best_score = all_lr_seqs.iter().map(|lr_seq| {
            let matches = find_kmer_matches(*sr_seq, lr_seq, kmer_size);
            let sparse_al = lcskpp(&matches, kmer_size);

            (100.0 * sparse_al.score as f64 / sr_seq.len() as f64) as usize
        }).max();

        best_score.unwrap() >= min_score_pct
    }).collect::<Vec<&Vec<u8>>>();

    skydive::elog!("{} filtered short-read sequences sparse alignment score >= {}%.", filtered_sr_seqs.len(), min_score_pct);

    let mut writer = BufWriter::new(File::create(&output).unwrap());
    for (i, sr_seq) in filtered_sr_seqs.iter().enumerate() {
        writer.write_all(format!(">sr_{}\n{}\n", i, String::from_utf8(sr_seq.to_vec()).unwrap()).as_bytes()).unwrap();
    }

    skydive::elog!("Wrote {} filtered short-read sequences to {}.", filtered_sr_seqs.len(), output.to_str().unwrap());
}
