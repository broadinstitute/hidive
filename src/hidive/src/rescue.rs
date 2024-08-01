use std::fs::File;
use std::io::Write;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::{collections::HashSet, path::PathBuf};

use bio::io::fasta::Reader;
use num_format::{Locale, ToFormattedString};

use rayon::iter::ParallelBridge;
use rayon::iter::ParallelIterator;
use rayon::prelude::*;

use rust_htslib::bam::{FetchDefinition, Read};

// Import the skydive module, which contains the necessary functions for staging data
use skydive;

pub fn start(output: &PathBuf, fasta_paths: &Vec<PathBuf>, seq_paths: &Vec<PathBuf>) {
    let fasta_urls = skydive::parse::parse_file_names(fasta_paths);
    let seq_urls = skydive::parse::parse_file_names(seq_paths);

    let k = 15;
    let num_threshold = 10;

    // Get the system's temporary directory path
    let cache_path = std::env::temp_dir();
    skydive::elog!("Intermediate data will be stored at {:?}.", cache_path);

    // Read the FASTA files and prepare a hashmap of k-mers to search for in the other files.
    let mut kmer_set = HashSet::new();
    for fasta_url in fasta_urls {
        let reader = Reader::from_file(fasta_url.to_file_path().unwrap()).unwrap();

        for record in reader.records().flatten() {
            let fw_seq = record.seq();

            let kmers = fw_seq
                .windows(k)
                .map(skydive::ldbg::LdBG::canonicalize_kmer)
                .collect::<HashSet<Vec<u8>>>();

            kmer_set.extend(kmers);
        }
    }

    // Read the CRAM files and search for the k-mers in each read.
    let mut all_records = Vec::new();
    for seq_url in seq_urls {
        let mut reader = skydive::stage::open_bam(&seq_url).unwrap();

        let progress_bar =
            skydive::utils::default_unbounded_progress_bar("Searching for similar reads (0 found)");

        // Create some thread-safe counters.
        let found_items = Arc::new(AtomicUsize::new(0));
        let processed_items = Arc::new(AtomicUsize::new(0));
        let progress_bar = Arc::new(progress_bar);

        const UPDATE_FREQUENCY: usize = 1_000_000;

        // Iterate over the records in parallel.
        reader
            .fetch(FetchDefinition::All)
            .expect("Failed to fetch reads");
        let records: Vec<_> = reader
            .records()
            .par_bridge()
            .flat_map(|record| record.ok())
            .filter_map(|read| {
                // Increment progress counters and update the progress bar every once in a while.
                let current_processed = processed_items.fetch_add(1, Ordering::Relaxed);
                if current_processed % UPDATE_FREQUENCY == 0 {
                    progress_bar.set_message(format!(
                        "Searching for similar reads ({} found)",
                        found_items.load(Ordering::Relaxed)
                    ));
                    progress_bar.inc(UPDATE_FREQUENCY as u64);
                }

                // Count the number of k-mers found in our hashset from the long reads.
                let seq = read.seq().as_bytes();
                let num_kmers = seq
                    .par_windows(k)
                    .map(|kmer| kmer_set.contains(&skydive::ldbg::LdBG::canonicalize_kmer(kmer)))
                    .filter(|&contains| contains)
                    .count();

                // If the number of k-mers found is greater than the threshold, add the read to the list of similar reads.
                if num_kmers > num_threshold {
                    let current_found = found_items.fetch_add(1, Ordering::Relaxed);
                    if current_found % UPDATE_FREQUENCY == 0 {
                        progress_bar.set_message(format!(
                            "Searching for similar reads ({} found)",
                            (current_found + 1).to_formatted_string(&Locale::en)
                        ));
                    }
                    Some(read)
                } else {
                    None
                }
            })
            .collect();

        // Final update to ensure the last counts are displayed
        let final_found = found_items.load(Ordering::Relaxed);
        let final_processed = processed_items.load(Ordering::Relaxed);
        progress_bar.set_message(format!(
            "Found {} reads in {} reads total.",
            final_found.to_formatted_string(&Locale::en),
            final_processed.to_formatted_string(&Locale::en)
        ));
        progress_bar.finish();

        // Add the reads to the list of all records.
        all_records.extend(records);
    }

    let mut file = File::create(output).expect("Could not open file for writing.");
    for (i, record) in all_records.iter().enumerate() {
        writeln!(file, ">read_{}", i).expect("Could not write to file");
        writeln!(
            file,
            "{}",
            String::from_utf8_lossy(&record.seq().as_bytes())
        )
        .expect("Could not write to file");
    }

    skydive::elog!(
        "Wrote {} reads to {}.",
        all_records.len(),
        output.to_str().unwrap()
    );
}
