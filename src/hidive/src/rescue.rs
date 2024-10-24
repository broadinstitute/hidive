use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::{collections::HashSet, path::PathBuf};

use bio::io::fasta::Reader;
use flate2::write::GzEncoder;
use flate2::Compression;
use minimap2::Aligner;
use num_format::{Locale, ToFormattedString};

use rayon::iter::ParallelBridge;
use rayon::iter::ParallelIterator;
use rayon::prelude::*;

use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{FetchDefinition, Read};

// Import the skydive module, which contains the necessary functions for staging data

pub fn start(
    output: &PathBuf,
    kmer_size: usize,
    min_kmers_pct: usize,
    contigs: &Vec<String>,
    ref_path: Option<PathBuf>,
    fasta_paths: &Vec<PathBuf>,
    seq_paths: &Vec<PathBuf>,
) {
    let fasta_urls = skydive::parse::parse_file_names(fasta_paths);
    let seq_urls = skydive::parse::parse_file_names(seq_paths);

    // let contigs = contigs.iter().map(|c| c.to_string()).collect::<HashSet<String>>();

    // Get the system's temporary directory path
    let cache_path = std::env::temp_dir();
    skydive::elog!("Intermediate data will be stored at {:?}.", cache_path);

    // Define the fetch definitions for the aligner.
    let mut fetches = HashSet::new();

    // Read the FASTA files and prepare a hashmap of k-mers to search for in the other files.
    let mut kmer_set = HashSet::new();
    for fasta_url in fasta_urls {
        let reader = Reader::from_file(fasta_url.to_file_path().unwrap()).unwrap();

        let mut reads = Vec::new();
        for record in reader.records().flatten() {
            let fw_seq = record.seq();
            let rl_seq = skydive::utils::homopolymer_compressed(fw_seq);

            let kmers = rl_seq
                .windows(kmer_size)
                .map(skydive::utils::canonicalize_kmer)
                .collect::<HashSet<Vec<u8>>>();

            kmer_set.extend(kmers);

            reads.push(fw_seq.to_vec());
        }

        if let Some(ref ref_path) = ref_path {
            skydive::elog!("Searching for relevant loci in reference genome...");

            let aligner = Aligner::builder()
                .sr()
                .with_sam_hit_only()
                .with_index(ref_path, None)
                .expect("Unable to build index");

            let mappings = reads
                .par_iter()
                .map(|seq| {
                    seq
                        .windows(150)
                        .map(|window| {
                            aligner.map(window, false, false, None, None).unwrap()
                        })
                        .flatten()
                        .collect::<Vec<_>>()
                })
                .flatten()
                .collect::<Vec<_>>();

            let sub_fetches = mappings.iter().filter_map(|m| {
                if let Some(target_name) = &m.target_name {
                    Some(target_name.to_string().as_bytes().to_vec())
                } else {
                    None
                }
            }).collect::<HashSet<_>>();

            fetches.extend(sub_fetches);
        } else {
            skydive::elog!("Searching for relevant loci in {:?}...", contigs);

            fetches.extend(contigs.iter().map(|c| c.as_bytes().to_vec()));
        }

        fetches.insert("*".as_bytes().to_vec());
    }

    // Read the CRAM files and search for the k-mers in each read.
    let mut tid_to_chrom: HashMap<i32, String> = HashMap::new();

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

        // Create an immutable map of tid to chromosome name
        tid_to_chrom.extend(reader
            .header()
            .target_names()
            .iter()
            .enumerate()
            .map(|(tid, &name)| (tid as i32, String::from_utf8_lossy(name).into_owned()))
            .collect::<HashMap<_, _>>());

        for fetch in &fetches {
            let fetch_definition = if String::from_utf8_lossy(fetch) != "*" {
                FetchDefinition::String(fetch)
            } else {
                FetchDefinition::Unmapped
            };

            // Iterate over the records in parallel.
            reader
                .fetch(fetch_definition)
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
                            "Searching for similar reads ({} found, most recent at {}:{})",
                            found_items
                                .load(Ordering::Relaxed)
                                .to_formatted_string(&Locale::en),
                            tid_to_chrom.get(&read.tid()).unwrap(),
                            read.reference_start().to_formatted_string(&Locale::en)
                        ));
                        progress_bar.inc(UPDATE_FREQUENCY as u64);
                    }

                    // Count the number of k-mers found in our hashset from the long reads.
                    let fw_seq = read.seq().as_bytes();
                    let rl_seq = skydive::utils::homopolymer_compressed(&fw_seq);

                    let num_kmers = rl_seq
                        .par_windows(kmer_size)
                        .map(|kmer| kmer_set.contains(&skydive::utils::canonicalize_kmer(kmer)))
                        .filter(|&contains| contains)
                        .count();

                    // If the number of k-mers found is greater than the threshold, add the read to the list of similar reads.
                    if num_kmers as f32 / rl_seq.len() as f32 > min_kmers_pct as f32 / 100.0 {
                        let current_found = found_items.fetch_add(1, Ordering::Relaxed);
                        if current_found % UPDATE_FREQUENCY == 0 {
                            progress_bar.set_message(format!(
                                "Searching for similar reads ({} found, most recent at {}:{})",
                                (current_found + 1).to_formatted_string(&Locale::en),
                                tid_to_chrom.get(&read.tid()).unwrap(),
                                read.reference_start().to_formatted_string(&Locale::en)
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
    }

    let file = File::create(output).expect("Could not open file for writing.");
    let mut writer: Box<dyn Write> = if output.to_str().unwrap().ends_with(".gz") {
        Box::new(GzEncoder::new(file, Compression::default()))
    } else {
        Box::new(file)
    };

    for record in all_records.iter() {
        writeln!(writer, ">read_{}_{}_{}_{}", String::from_utf8_lossy(&record.qname()), tid_to_chrom.get(&record.tid()).unwrap_or(&String::from("*")), record.reference_start(), record.reference_end()).expect("Could not write to file");
        writeln!(
            writer,
            "{}",
            String::from_utf8_lossy(&record.seq().as_bytes())
        )
        .expect("Could not write to file");
    }

    if !all_records.is_empty() {
        skydive::elog!(
            "Wrote {} reads to {}.",
            all_records.len().to_formatted_string(&Locale::en),
            output.to_str().unwrap()
        );
    } else {
        skydive::elog!("No reads were found in the CRAM files. Aborting.");
        std::process::exit(1);
    }
}
