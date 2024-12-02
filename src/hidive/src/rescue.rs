use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::{collections::HashSet, path::PathBuf};
use std::cmp::PartialEq;
use clap::{ValueEnum};
use bio::data_structures::interval_tree::IntervalTree;
use bio::io::fasta::Reader;
use bio::utils::Interval;
use flate2::write::GzEncoder;
use flate2::Compression;
use minimap2::{Aligner, Mapping};
use needletail::Sequence;
use num_format::{Locale, ToFormattedString};


use indicatif::ProgressBar;
use rayon::iter::ParallelBridge;
use rayon::iter::ParallelIterator;
use rayon::prelude::*;

use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{FetchDefinition, IndexedReader, Read, Record as BamRecord};
use std::option::Option;

// Import the skydive module, which contains the necessary functions for staging data
use skydive;

#[derive(Clone, PartialEq, ValueEnum, Debug)]
pub(crate) enum SearchOption {
    All,
    Contig,
    ContigAndInterval,
    Unmapped,
}

pub fn start(
    output: &PathBuf,
    kmer_size: usize,
    min_kmers_pct: usize,
    search_option: SearchOption,
    ref_path: Option<PathBuf>,
    loci: Option<Vec<String>>,
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

        // Set fetches with the input loci or the detected relevant loci.
        if search_option == SearchOption::Contig || search_option == SearchOption::ContigAndInterval {
            match populate_fetches(&search_option, &loci, &ref_path, &reads, &mut fetches) {
                Ok(_) => {},
                Err(e) => {
                    skydive::elog!("Error: {}", e);
                    std::process::exit(1);
                }
            }
            print_fetches_info(&fetches);

            } else {
                assert!(
                    (search_option == SearchOption::All || search_option == SearchOption::Unmapped) && loci.is_none(),
                    "Assertion failed: search_option is 'All' or 'Unmapped' and loci is NOT None"
                );
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

        if search_option == SearchOption::Contig || search_option == SearchOption::ContigAndInterval {
            // Use the fetches to retrieve the records.
            retrieve_records(
                &mut reader,
                Some(&fetches),
                &search_option,
                kmer_size,
                min_kmers_pct,
                &kmer_set,
                &tid_to_chrom,
                &progress_bar,
                &found_items,
                &processed_items,
                &mut all_records,
            );
        } else {
            // Retrieve all or unmapped reads.
            retrieve_records(
                &mut reader,
                None,
                &search_option,
                kmer_size,
                min_kmers_pct,
                &kmer_set,
                &tid_to_chrom,
                &progress_bar,
                &found_items,
                &processed_items,
                &mut all_records,
            );
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
}


////////////////////////////////////////////////////////////////
//// print_fetches_info and its helper functions
////
////////////////////////////////////////////////////////////////

/// This function will populate the fetches with the input loci or the detected relevant loci.
/// Function should only be called if the search option is set to 'Contig' or 'ContigAndInterval'.
///
/// # Arguments
/// - `search_option` - The search option.
/// - `loci` - The loci to search for.
/// - `ref_path` - The path to the reference genome.
/// - `reads` - The reads to map.
/// - `fetches` - The fetches to populate.
///
/// # Panics
/// If the search option is set to 'Contig' or 'ContigAndInterval' and no reference genome is provided.
fn populate_fetches(
    search_option: &SearchOption,
    loci: &Option<Vec<String>>,
    ref_path: &Option<PathBuf>,
    reads: &Vec<Vec<u8>>,
    fetches: &mut Vec<(String, Interval<i32>)>,
) -> Result<(), String>{
    if !matches!(search_option, SearchOption::Contig | SearchOption::ContigAndInterval) {
        return Err("Function should only be called if the search option is 'Contig' or 'ContigAndInterval'".to_string());
    }

    if let Some(loci) = loci {
        skydive::elog!("Setting loci from input...");
        let loci = skydive::parse::parse_loci(loci, 0);
        for (contig, start, stop, _) in loci {
            match Interval::new(start as i32..stop as i32) {
                Ok(interval) => fetches.push((contig, interval)),
                Err(e) => {
                    skydive::elog!("Invalid interval for contig {}: {:?}, start: {}, stop: {}. Error: {}", contig, start..stop, start, stop, e);
                    return Err(format!("Error creating interval for contig {}: {:?}", contig, e));
                }
            }
        }
    } else if let Some(ref ref_path) = ref_path {
        skydive::elog!("Searching for relevant loci in reference genome...");
        detect_relevant_loci(ref_path, reads, fetches);
        fetches.push(("*".to_string(), Interval::new(0..1).unwrap()));
    } else {
        skydive::elog!("No reference genome provided but search option is set to 'Contig' or 'ContigAndInterval'. Aborting.");
        return Err("No reference genome provided.".to_string());
    }

    Ok(())
}

/// This function will print the fetches' info.
/// # Arguments
/// - `fetches` - The fetches to print.
///
/// # Panics
///
/// If the total length of the fetches cannot be calculated.
fn print_fetches_info(fetches: &[(String, Interval<i32>)]) {
    let total_length = fetches.iter().map(|(_, interval)| interval.end - interval.start).sum::<i32>();
    skydive::elog!(
        " -- will search unaligned reads and {} bases in {} contigs.",
        total_length.to_formatted_string(&Locale::en),
        fetches.len().to_formatted_string(&Locale::en)
    );

    for (contig, interval) in fetches {
        skydive::elog!(" -- {}:{:?}", contig, interval);
    }
}

//////////////////////////////////////////////////////////////////
//// detect_relevant_loci and its helper functions
////
//////////////////////////////////////////////////////////////////


/// This function detects the relevant loci in the reference genome.
///
/// # Arguments
///
/// * `ref_path` - The path to the reference genome.
/// * `reads` - The reads to map.
/// * `fetches` - The fetches to populate.
///
/// # Panics
///
/// If no mappings are found for the provided sequences.
fn detect_relevant_loci(ref_path: &PathBuf, reads: &Vec<Vec<u8>>, fetches: &mut Vec<(String, Interval<i32>)>) {
    let lr_aligner = initialize_aligner(ref_path, false);

    let lr_mappings = map_sequences(&lr_aligner, &reads, false);

    for lr_mapping in &lr_mappings {
        skydive::elog!("lr mapping {:?}:{}-{}", lr_mapping.target_name, lr_mapping.target_start, lr_mapping.target_end);
    }

    let sr_aligner = initialize_aligner(ref_path, true);

    let sr_mappings = map_sequences(&sr_aligner, &reads, true);

    let mut contig_lengths = HashMap::new();
    let loci = collect_loci(&lr_mappings, &sr_mappings, &mut contig_lengths);

    // First pass - put all the intervals into a tree.
    let mut trees = BTreeMap::new();
    populate_interval_trees(&mut trees, &loci, &contig_lengths);

    // Second pass
    merge_overlapping_intervals(&mut trees, &loci, &contig_lengths);

    let sub_fetches: HashSet<_> = trees.iter()
        .flat_map(|(contig_name, tree)| {
            tree.find(Interval::new(0..i32::MAX).unwrap())
                .map(move |entry| (contig_name.clone(), entry.interval().clone()))
        })
        .collect();

    fetches.extend(sub_fetches.into_iter());
}

/// This function initializes the aligner.
///
/// # Arguments
///
/// * `ref_path` - The path to the reference genome.
/// * `is_sr` - Whether the reads are short reads or not.
///
/// # Returns
///
/// An aligner object.
///
/// # Panics
///
/// If the index cannot be built.
fn initialize_aligner(ref_path: &PathBuf, is_sr: bool) -> Aligner {
    if is_sr {
        Aligner::builder()
            .sr()
            .with_sam_hit_only()
            .with_index(ref_path, None)
            .expect("Unable to build index")
    } else {
        Aligner::builder()
            .map_hifi()
            .with_sam_hit_only()
            .with_index(ref_path, None)
            .expect("Unable to build index")
    }
}

/// This function maps the sequences to the reference genome using the aligner.
///
/// # Arguments
///
/// * `aligner` - The aligner object to use for mapping.
/// * `reads` - The reads to map.
/// * `is_sr` - Whether the reads are short reads or not.
///
/// # Returns
///
/// A vector of mappings.
///
/// # Panics
///
/// If no mappings are found for the provided sequences.
fn map_sequences(aligner: &Aligner, reads: &[Vec<u8>], is_sr: bool) -> Vec<Mapping> {
    let mappings: Vec<Mapping> = reads.par_iter().flat_map(|seq| {
        if is_sr {
            seq.windows(150).flat_map(|window| aligner.map(window, false, false, None, None).unwrap_or_else(|_| vec![])).collect::<Vec<_>>()
        } else {
            aligner.map(seq, false, false, None, None).unwrap_or_else(|_| vec![])
        }
    }).collect();

    if mappings.is_empty() {
        panic!("No mappings found for the provided sequences.");
    }

    mappings
}

/// This function collects the loci from the mappings.
///
/// # Arguments
///
/// * `lr_mappings` - The long read mappings.
/// * `sr_mappings` - The short read mappings.
///
/// # Returns
///
/// A hashset of loci.
///
/// # Panics
///
/// If the target name is not found in the mappings.
fn collect_loci(
    lr_mappings: &[Mapping],
    sr_mappings: &[Mapping],
    contig_lengths: &mut HashMap<String, i32>,
) -> HashSet<(Vec<u8>, i32, i32)> {
    lr_mappings.iter().chain(sr_mappings.iter()).filter_map(|m| {
        if let Some(target_name) = &m.target_name {
            contig_lengths.insert(target_name.to_string(), m.target_len);
            Some((target_name.to_string().as_bytes().to_vec(), m.target_start, m.target_end))
        } else {
            None
        }
    }).collect::<HashSet<_>>()
}

/// This function populates the interval trees with the loci.
///
/// # Arguments
///
/// * `trees` - The interval trees to populate.
/// * `loci` - The loci to insert into the trees.
/// * `contig_lengths` - The lengths of the contigs.
///
/// # Panics
///
/// If the contig name is not found in the contig lengths.
fn populate_interval_trees(
    trees: &mut BTreeMap<String, IntervalTree<i32, ()>>,
    loci: &HashSet<(Vec<u8>, i32, i32)>,
    contig_lengths: &HashMap<String, i32>
) {
    for (seq, start, end) in loci {
        let contig_name = String::from_utf8_lossy(&seq).to_string();
        let contig_length = contig_lengths.get(&contig_name).unwrap();

        if !trees.contains_key(&contig_name) {
            trees.insert(contig_name.clone(), IntervalTree::new());
        }

        let interval = Interval::new(start.saturating_sub(50000).max(0)..end.saturating_add(50000).min(*contig_length)).unwrap();
        trees.get_mut(&contig_name).unwrap().insert(interval, ());
    }
}

/// This function merges overlapping intervals.
///
/// # Arguments
///
/// * `trees` - The interval trees to merge.
/// * `loci` - The loci to merge.
/// * `contig_lengths` - The lengths of the contigs.
///
/// # Panics
///
/// If the contig name is not found in the trees.
fn merge_overlapping_intervals(
    trees: &mut BTreeMap<String, IntervalTree<i32, ()>>,
    loci: &HashSet<(Vec<u8>, i32, i32)>,
    contig_lengths: &HashMap<String, i32>
) {
    for (seq, start, end) in loci {
        let contig_name = String::from_utf8_lossy(&seq).to_string();
        let contig_length = contig_lengths.get(&contig_name).unwrap();
        let tree = trees.get_mut(&contig_name).unwrap();

        let current_interval = Interval::new(start.saturating_sub(50000).max(0)..end.saturating_add(50000).min(*contig_length)).unwrap();
        let overlaps = tree.find(current_interval.clone()).collect::<Vec<_>>();

        let mut merged_intervals = IntervalTree::new();

        if !overlaps.is_empty() {
            let min_start = overlaps.iter().map(|o| o.interval().start).min().unwrap().min(*start);
            let max_end = overlaps.iter().map(|o| o.interval().end).max().unwrap().max(*end);
            merged_intervals.insert(Interval::new(min_start..max_end).unwrap(), ());
        } else {
            merged_intervals.insert(current_interval, ());
        }

        trees.insert(contig_name.clone(), merged_intervals);
    }
}

////////////////////////////////////////////////////////////////
//// process_fetches and its helper functions
////
////////////////////////////////////////////////////////////////

/// This function retrieves the records from the reader.
///
/// # Arguments
/// - `fetches` - The fetches to process.
/// - `reader` - The reader to use for fetching reads.
/// - `kmer_size` - The size of the k-mers to search for.
/// - `min_kmers_pct` - The minimum percentage of k-mers to search for.
/// - `search_option` - Whether to search all reads or not.
/// - `kmer_set` - The set of k-mers to search for.
/// - `tid_to_chrom` - The mapping of TID to chromosome name.
/// - `progress_bar` - The progress bar to update.
/// - `found_items` - The number of found items.
/// - `processed_items` - The number of processed items.
/// - `all_records` - The vector of all records.
///
fn retrieve_records(
    reader: &mut IndexedReader,
    fetches: Option<&[(String, Interval<i32>)]>,
    search_option: &SearchOption,
    kmer_size: usize,
    min_kmers_pct: usize,
    kmer_set: &HashSet<Vec<u8>>,
    tid_to_chrom: &HashMap<i32, String>,
    progress_bar: &Arc<ProgressBar>,
    found_items: &Arc<AtomicUsize>,
    processed_items: &Arc<AtomicUsize>,
    all_records: &mut Vec<BamRecord>,
) {
    const UPDATE_FREQUENCY: usize = 1_000_000;

    if let Some(fetches) = fetches {
        let fetches_updated = prepare_fetches(fetches, search_option);

        for (contig, interval) in fetches_updated {
            let fetch_definition = create_fetch_definition(search_option, Some(&contig), Some(interval));
            reader.fetch(fetch_definition).expect("Failed to fetch reads");

            let records = filter_and_collect_records(
                reader,
                kmer_size,
                min_kmers_pct,
                kmer_set,
                tid_to_chrom,
                progress_bar,
                found_items,
                processed_items,
            );

            all_records.extend(records);
        }
    } else {
        let fetch_definition = create_fetch_definition(search_option, None, None);
        reader.fetch(fetch_definition).expect("Failed to fetch reads");

        let records = filter_and_collect_records(
            reader,
            kmer_size,
            min_kmers_pct,
            kmer_set,
            tid_to_chrom,
            progress_bar,
            found_items,
            processed_items,
        );

        all_records.extend(records);
    }

    finalize_progress(progress_bar, found_items, processed_items);
}

/// This function processes all or unmapped reads.
/// # Arguments
/// - `reader` - The reader to use for fetching reads.
/// - `kmer_size` - The size of the k-mers to search for.
/// - `min_kmers_pct` - The minimum percentage of k-mers to search for.
/// - `kmer_set` - The set of k-mers to search for.
/// - `tid_to_chrom` - The mapping of TID to chromosome name.
/// - `progress_bar` - The progress bar to update.
/// - `found_items` - The number of found items.
/// - `processed_items` - The number of processed items.
///
/// # Panics
/// If the search option is not 'All' or 'Unmapped'.
///
/// # Returns
/// A vector of all records.
fn filter_and_collect_records(
    reader: &mut IndexedReader,
    kmer_size: usize,
    min_kmers_pct: usize,
    kmer_set: &HashSet<Vec<u8>>,
    tid_to_chrom: &HashMap<i32, String>,
    progress_bar: &Arc<ProgressBar>,
    found_items: &Arc<AtomicUsize>,
    processed_items: &Arc<AtomicUsize>,
) -> Vec<BamRecord> {
    const UPDATE_FREQUENCY: usize = 1_000_000;

    reader
        .records()
        .par_bridge()
        .flat_map(|record| record.ok())
        .filter_map(|read| {
            update_processed_progress(processed_items, UPDATE_FREQUENCY, tid_to_chrom, &read, &progress_bar);

            if is_valid_read(&read, kmer_size, min_kmers_pct, kmer_set) {
                update_found_progress(found_items, UPDATE_FREQUENCY, tid_to_chrom, &read, &progress_bar);
                Some(read)
            } else {
                None
            }
        })
        .collect()
}

/// This function prepares the fetches for searching.
/// If the search option is 'Contig', it will extract the unique contigs from the fetches and add a dummy interval.
/// If the search option is 'ContigAndInterval', it will return the fetches as is.
/// In both cases, the fetches will be returned as a vector of contigs and intervals.
///
/// # Arguments
/// - `fetches` - The fetches to prepare.
/// - `search_option` - Whether to search all reads or not.
///
/// # Returns
/// A vector of fetches to search.
///
/// # Panics
/// If the contig is not valid.
fn prepare_fetches(
    fetches: &[(String, Interval<i32>)], search_option: &SearchOption
) -> Vec<(String, Interval<i32>)> {
    if *search_option == SearchOption::Contig {
        let uniq_contigs: HashSet<_> = fetches.iter().map(|(contig, _)| contig.clone()).collect();
        uniq_contigs.into_iter().map(|contig| (contig, Interval::new(0..1).unwrap())).collect()
    } else if *search_option == SearchOption::ContigAndInterval {
        fetches.iter()
            .map(|(contig, interval)| (contig.clone(), interval.clone()))
            .collect()
    } else {
        panic!("Invalid search option");
    }
}

/// This function creates a fetch definition.
///
/// # Arguments
/// - `search_option` - Whether to search all reads or not.
/// - `contig` - The contig to fetch.
/// - `interval` - The interval to fetch.
///
/// # Returns
/// A fetch definition.
fn create_fetch_definition<'a>(
    search_option: &'a SearchOption,
    contig: Option<&'a String>,
    interval: Option<Interval<i32>>,
) -> FetchDefinition<'a> {
    match (contig, search_option) {
        (None, SearchOption::All) => FetchDefinition::All,
        (None, SearchOption::Unmapped) => FetchDefinition::Unmapped,
        (Some(contig), SearchOption::Contig) if contig != "*" => {
            FetchDefinition::String(contig.as_bytes())
        },
        (Some(contig), SearchOption::ContigAndInterval) if contig != "*" => {
            let interval = interval.unwrap_or_else(|| Interval::new(0..1).unwrap());
            FetchDefinition::RegionString(contig.as_bytes(), interval.start as i64, interval.end as i64)
        },
        (Some(contig), _) if contig == "*" => FetchDefinition::Unmapped,
        _ => panic!(
            "Invalid fetch definition settings. \
                If you are using the 'Contig' or 'ContigAndInterval' search options, \
                you must provide a contig name. \
                If you are using the 'All' or 'Unmapped' search options, \
                you must NOT provide a contig name. \
             search option: '{:?}' , contig: '{:?}', interval: '{:?}'",
            search_option,
            contig.as_deref().unwrap_or(&"None".to_string()),
            interval.map(|i| format!("{:?}", i)).unwrap_or_else(|| "None".to_string())
        ),
    }
}

/// This function updates the processed progress.
///
/// # Arguments
/// - `processed_items` - The number of processed items.
/// - `update_freq` - The update frequency.
/// - `tid_to_chrom` - The mapping of TID to chromosome name.
/// - `read` - The read to process.
/// - `progress_bar` - The progress bar to update.
///
/// # Panics
/// If the chromosome name is not found in the mapping.
fn update_processed_progress(processed_items: &Arc<AtomicUsize>, update_freq: usize, tid_to_chrom: &HashMap<i32, String>, read: &BamRecord, progress_bar: &Arc<ProgressBar>) {
    let current_processed = processed_items.fetch_add(1, Ordering::Relaxed);
    if current_processed % update_freq == 0 {
        let unknown_chrom = "Unknown".to_string();
        let chrom_name = tid_to_chrom.get(&read.tid()).unwrap_or(&unknown_chrom);
        progress_bar.set_message(format!(
            "Searching for similar reads ({} processed, most recent at {}:{})",
            current_processed.to_formatted_string(&Locale::en),
            chrom_name,
            read.reference_start().to_formatted_string(&Locale::en)
        ));
        progress_bar.inc(update_freq as u64);
    }
}

/// This function updates the found progress.
///
/// # Arguments
/// - `found_items` - The number of found items.
/// - `update_freq` - The update frequency.
/// - `tid_to_chrom` - The mapping of TID to chromosome name.
/// - `read` - The read to process.
/// - `progress_bar` - The progress bar to update.
///
/// # Panics
/// If the chromosome name is not found in the mapping.
fn update_found_progress(found_items: &Arc<AtomicUsize>, update_freq: usize, tid_to_chrom: &HashMap<i32, String>, read: &BamRecord, progress_bar: &Arc<ProgressBar>) {
    let current_found = found_items.fetch_add(1, Ordering::Relaxed);
    if current_found % update_freq == 0 {
        let unknown_chrom = "Unknown".to_string();
        let chrom_name = tid_to_chrom.get(&read.tid()).unwrap_or(&unknown_chrom);
        progress_bar.set_message(format!(
            "Found {} similar reads (most recent at {}:{})",
            current_found.to_formatted_string(&Locale::en),
            chrom_name,
            read.reference_start().to_formatted_string(&Locale::en)
        ));
    }
}

/// This function checks if a read is valid.
///
/// # Arguments
/// - `read` - The read to check.
/// - `kmer_size` - The size of the k-mers to search for.
/// - `min_kmers_pct` - The minimum percentage of k-mers to search for.
/// - `kmer_set` - The set of k-mers to search for.
///
/// # Returns
/// Whether the read is valid or not.
fn is_valid_read(read: &BamRecord, kmer_size: usize, min_kmers_pct: usize, kmer_set: &HashSet<Vec<u8>>) -> bool {
    let read_seq = read.seq().as_bytes();
    let rl_seq = skydive::utils::homopolymer_compressed(&read_seq);
    let num_kmers = rl_seq.par_windows(kmer_size)
        .filter(|kmer| kmer_set.contains(&skydive::utils::canonicalize_kmer(kmer)))
        .count();
    (num_kmers as f32 / rl_seq.len() as f32) > (min_kmers_pct as f32 / 100.0)
}

/// This function finalizes the progress.
///
/// # Arguments
/// - `progress_bar` - The progress bar to finalize.
/// - `found_items` - The number of found items.
/// - `processed_items` - The number of processed items.
fn finalize_progress(progress_bar: &Arc<ProgressBar>, found_items: &Arc<AtomicUsize>, processed_items: &Arc<AtomicUsize>) {
    let final_found = found_items.load(Ordering::Relaxed);
    let final_processed = processed_items.load(Ordering::Relaxed);
    progress_bar.set_message(format!(
        "Found {} reads in {} reads total.",
        final_found.to_formatted_string(&Locale::en),
        final_processed.to_formatted_string(&Locale::en)
    ));
    progress_bar.finish();
}
