use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::File;
use std::io::Write;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::{collections::HashSet, path::PathBuf};

use bio::data_structures::interval_tree::IntervalTree;
use bio::io::fasta::Reader;
use bio::utils::Interval;
use flate2::write::GzEncoder;
use flate2::Compression;
use minimap2::Aligner;
use needletail::Sequence;
use num_format::{Locale, ToFormattedString};

use rayon::iter::ParallelBridge;
use rayon::iter::ParallelIterator;
use rayon::prelude::*;

use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{FetchDefinition, Read};

// Import the skydive module, which contains the necessary functions for staging data
use skydive;

pub fn start(
    output: &PathBuf,
    kmer_size: usize,
    min_kmers_pct: usize,
    ref_path: Option<PathBuf>,
    fasta_paths: &Vec<PathBuf>,
    seq_paths: &Vec<PathBuf>,
) {
    let fasta_urls = skydive::parse::parse_file_names(fasta_paths);
    let seq_urls = skydive::parse::parse_file_names(seq_paths);

    // Get the system's temporary directory path
    let cache_path = std::env::temp_dir();
    skydive::elog!("Intermediate data will be stored at {:?}.", cache_path);

    // Define the fetch definitions for the aligner.
    let mut fetches = Vec::new();

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
            detect_relevant_loci(ref_path, reads, &mut fetches);
        }

        fetches.push(("*".to_string(), Interval::new(0..1).unwrap()));

        let total_length = fetches.iter().map(|(_, interval)| interval.end - interval.start).sum::<i32>();
        skydive::elog!(
            " -- will search unaligned reads and {} bases in {} contigs.",
            total_length.to_formatted_string(&Locale::en),
            fetches.len().to_formatted_string(&Locale::en)
        );

        for (contig, interval) in &fetches {
            skydive::elog!(" -- {}:{:?}", contig, interval);
        }
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

        for (contig, interval) in &fetches {
            let fetch_definition = if contig != "*" {
                FetchDefinition::RegionString(contig.as_bytes(), interval.start as i64, interval.end as i64)
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
                    // let read_seq = if read.is_reverse() { read.seq().as_bytes().reverse_complement() } else { read.seq().as_bytes() };
                    let read_seq = read.seq().as_bytes();
                    let rl_seq = skydive::utils::homopolymer_compressed(&read_seq);

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
        let fw_seq = record.seq().as_bytes();
        let rc_seq = fw_seq.reverse_complement();

        writeln!(writer, ">read_{}_{}_{}_{}", String::from_utf8_lossy(&record.qname()), tid_to_chrom.get(&record.tid()).unwrap_or(&"*".to_string()), record.reference_start(), record.reference_end()).expect("Could not write to file");
        writeln!(
            writer,
            "{}",
            if record.is_reverse() { String::from_utf8_lossy(&rc_seq) } else { String::from_utf8_lossy(&fw_seq) }
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

fn detect_relevant_loci(ref_path: &PathBuf, reads: Vec<Vec<u8>>, fetches: &mut Vec<(String, Interval<i32>)>) {
    let lr_aligner = Aligner::builder()
        .map_hifi()
        .with_sam_hit_only()
        .with_index(ref_path, None)
        .expect("Unable to build index");

    let lr_mappings = reads
        .par_iter()
        .map(|seq| {
            lr_aligner.map(seq, false, false, None, None).unwrap()
        })
        .flatten()
        .collect::<Vec<_>>();

    for lr_mapping in &lr_mappings {
        skydive::elog!("lr mapping {:?}:{}-{}", lr_mapping.target_name, lr_mapping.target_start, lr_mapping.target_end);
    }

    let sr_aligner = Aligner::builder()
        .sr()
        .with_sam_hit_only()
        .with_index(ref_path, None)
        .expect("Unable to build index");

    let sr_mappings = reads
        .par_iter()
        .map(|seq| {
            seq
                .windows(150)
                .map(|window| sr_aligner.map(window, false, false, None, None).unwrap())
                .flatten()
                .collect::<Vec<_>>()
        })
        .flatten()
        .collect::<Vec<_>>();

    let mut contig_lengths = HashMap::new();
    let loci = lr_mappings.iter().chain(sr_mappings.iter()).filter_map(|m| {
        if let Some(target_name) = &m.target_name {
            contig_lengths.insert(target_name.to_string(), m.target_len);

            Some((target_name.to_string().as_bytes().to_vec(), m.target_start, m.target_end))
        } else {
            None
        }
    }).collect::<HashSet<_>>();

    // First pass - put all the intervals into a tree.
    let mut trees = BTreeMap::new();
    for (seq, start, end) in &loci {
        let contig_name = String::from_utf8_lossy(&seq).to_string();
        let contig_length = contig_lengths.get(&contig_name).unwrap();

        if !trees.contains_key(&contig_name) {
            trees.insert(contig_name.clone(), IntervalTree::new());
        }

        let interval = Interval::new(start.saturating_sub(50000).max(0)..end.saturating_add(50000).min(*contig_length)).unwrap();
        trees.get_mut(&contig_name).unwrap().insert(interval, ());
    }

    // Second pass - figure out which intervals overlap, and merge them.
    for (seq, start, end) in &loci {
        let contig_name = String::from_utf8_lossy(&seq).to_string();
        let contig_length = contig_lengths.get(&contig_name).unwrap();

        let tree = trees.get_mut(&contig_name).unwrap();

        let current_interval = Interval::new(start.saturating_sub(50000).max(0)..end.saturating_add(50000).min(*contig_length)).unwrap();
        let overlaps = tree.find(current_interval.clone()).collect::<Vec<_>>();

        let mut merged_intervals = IntervalTree::new();

        if !overlaps.is_empty() {
            // Calculate the new merged interval
            let min_start = overlaps.iter()
                .map(|o| o.interval().start)
                .min()
                .unwrap()
                .min(*start);
            let max_end = overlaps.iter()
                .map(|o| o.interval().end)
                .max()
                .unwrap()
                .max(*end);
            
            // Insert the new merged interval
            let merged_interval = Interval::new(min_start..max_end).unwrap();
            merged_intervals.insert(merged_interval, ());
        } else {
            merged_intervals.insert(current_interval, ());
        }

        trees.insert(contig_name.clone(), merged_intervals);
    }

    let sub_fetches: HashSet<_> = trees.iter()
        .flat_map(|(contig_name, tree)| {
            tree.find(Interval::new(0..i32::MAX).unwrap())
                .map(move |entry| (contig_name.clone(), entry.interval().clone()))
        })
        .collect();

    fetches.extend(sub_fetches.into_iter());
}
