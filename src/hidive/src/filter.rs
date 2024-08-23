use std::borrow::Borrow;
use std::collections::HashSet;
use std::io::Read;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::{
    fs::File,
    io::{BufWriter, Write},
};

use bio::alignment::sparse::*;

use flate2::read::GzDecoder;
use gbdt::decision_tree::{Data, DataVec};
use gbdt::gradient_boost::GBDT;
use num_format::{Locale, ToFormattedString};
use rayon::iter::{IntoParallelRefIterator, ParallelBridge, ParallelIterator};
use skydive::ldbg::LdBG;

pub fn start(
    output: &PathBuf,
    kmer_size: usize,
    // min_score_pct: usize,
    model_path: &PathBuf,
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

    // Assemble contigs.
    let l1 = LdBG::from_sequences("l1".to_string(), kmer_size, &all_lr_seqs)
        .score_kmers(model_path)
        .clean_paths(0.5)
        .clean_tips(2*kmer_size);

    // skydive::elog!(
    //     "K-mers with p < 0.5: {} / {} ({:.2}%)",
    //     num_below_threshold,
    //     num_total,
    //     num_below_threshold as f32 / num_total as f32 * 100.0
    // );

    skydive::elog!("Removed {} k-mers in {} paths", l1.cleaned_path_kmers, l1.cleaned_paths);
    skydive::elog!("Removed {} k-mers in {} tips", l1.cleaned_tip_kmers, l1.cleaned_tips);
    skydive::elog!("{} k-mers remaining", l1.kmers.len());

    // let all_corrected_seqs = l1.correct_seqs(&vec![all_lr_seqs[67].clone()]);
    let all_corrected_seqs = l1.correct_seqs(&all_lr_seqs);
    let l2 = l1.build_links(&all_corrected_seqs);
    let contigs = l2.assemble_all();

    skydive::elog!("Corrected sequences: {}", all_corrected_seqs.len());

    // Assemble contigs.
    // l3.links = LdBG::build_links(kmer_size, &all_seqs, &l3.kmers);
    // l3.links = LdBG::build_links(kmer_size, &all_lr_seqs2, &l3.kmers);

    for contig in contigs {
        skydive::elog!("Contig: {}", contig.len());
    }

    /*
    // Read and filter short reads.
    let mut filtered_sr_seqs: Vec<Vec<u8>> = Vec::new();
    for short_read_seq_url in &short_read_seq_urls {
        let basename = skydive::utils::basename_without_extension(
            &short_read_seq_url,
            &[".fasta.gz", ".fa.gz", ".fasta", ".fa"],
        );
        let fasta_path = short_read_seq_url.to_file_path().unwrap();

        skydive::elog!("Processing short-read sample {}...", basename);

        let file = File::open(&fasta_path).unwrap();
        let reader: Box<dyn Read> = if output.to_str().unwrap().ends_with(".gz") {
            Box::new(GzDecoder::new(file))
        } else {
            Box::new(file)
        };

        let fasta_reader = bio::io::fasta::Reader::new(reader);
        let all_reads = fasta_reader.records().flatten().collect::<Vec<_>>();

        let progress_bar = skydive::utils::default_unbounded_progress_bar(format!(
            "Filtering short reads with score < {}% (0 retained)",
            min_score_pct
        ));

        // Create some thread-safe counters.
        let found_items = Arc::new(AtomicUsize::new(0));
        let processed_items = Arc::new(AtomicUsize::new(0));
        let progress_bar = Arc::new(progress_bar);

        const UPDATE_FREQUENCY: usize = 1_000_000;

        filtered_sr_seqs.extend(
            all_reads
                .par_iter()
                .filter_map(|r| {
                    let sr_seq = r.seq().to_vec();

                    // Increment progress counters and update the progress bar every once in a while.
                    let current_processed = processed_items.fetch_add(1, Ordering::Relaxed);
                    if current_processed % UPDATE_FREQUENCY == 0 {
                        progress_bar.set_message(format!(
                            "Filtering short reads with score < {}% ({} retained)",
                            min_score_pct,
                            found_items
                                .load(Ordering::Relaxed)
                                .to_formatted_string(&Locale::en)
                        ));
                        progress_bar.inc(UPDATE_FREQUENCY as u64);
                    }

                    let best_score = all_lr_seqs
                        .iter()
                        .map(|lr_seq| {
                            let matches = find_kmer_matches(&sr_seq, lr_seq, kmer_size);
                            let sparse_al = lcskpp(&matches, kmer_size);

                            (100.0 * sparse_al.score as f64 / sr_seq.len() as f64) as usize
                        })
                        .max();

                    if best_score.unwrap() >= min_score_pct {
                        found_items.fetch_add(1, Ordering::Relaxed);
                        Some(sr_seq)
                    } else {
                        None
                    }
                })
                .collect::<Vec<Vec<u8>>>(),
        );

        // Final update to ensure the last counts are displayed
        let final_found = found_items.load(Ordering::Relaxed);
        let final_processed = processed_items.load(Ordering::Relaxed);
        progress_bar.set_message(format!(
            "Filtered {} short reads, retained {} with score >= {}%.",
            final_processed.to_formatted_string(&Locale::en),
            final_found.to_formatted_string(&Locale::en),
            min_score_pct
        ));
        progress_bar.finish();
    }

    skydive::elog!("{} long-read sequences.", all_lr_seqs.len());
    skydive::elog!("{} short-read sequences.", filtered_sr_seqs.len());

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
    */
}
