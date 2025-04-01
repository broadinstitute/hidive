use std::collections::{HashMap, HashSet};
use std::{fs::File, io::Write, path::PathBuf};

use indicatif::{ParallelProgressIterator, ProgressIterator};
use linked_hash_set::LinkedHashSet;
use rayon::prelude::*;
// use rayon::iter::IntoParallelRefIterator;
// use rayon::iter::ParallelIterator;

use skydive::ldbg::LdBG;
use skydive::mldbg::MLdBG;

pub fn start(
    output: &PathBuf,
    loci_list: &Vec<String>,
    kmer_size: usize,
    window: usize,
    model_path: &PathBuf,
    long_read_fasta_path: &PathBuf,
    short_read_fasta_path: &PathBuf,
) {
    // Get the system's temporary directory path
    let cache_path = std::env::temp_dir();
    skydive::elog!("Intermediate data will be stored at {:?}.", cache_path);

    let model_absolute_path = if model_path.is_absolute() {
        model_path.clone()
    } else {
        std::env::current_dir().unwrap().join(model_path)
    };

    // Load datasets
    let long_read_seq_urls = skydive::parse::parse_file_names(&[long_read_fasta_path.clone()]);

    skydive::elog!(
        "Processing short-read sample {}...",
        short_read_fasta_path.display()
    );
    let all_sr_seqs = skydive::utils::read_fasta(&vec![short_read_fasta_path.clone()]);
    skydive::elog!(" - {} short reads loaded", all_sr_seqs.len());

    // Iterate over loci
    let loci = skydive::parse::parse_loci(loci_list, window as u64)
        .into_iter()
        .collect::<Vec<_>>();

    // let progress_bar = skydive::utils::default_bounded_progress_bar("Processing loci", loci.len() as u64);

    let fa_file = std::sync::Mutex::new(File::create(output).unwrap());
    loci
        // .par_iter()
        // .progress_with(progress_bar)
        .iter()
        .inspect(|(chrom, start, end, name)| {
            skydive::elog!("Processing locus {} ({}:{}-{})", name, chrom, start, end);
        })
        .for_each(|(chrom, start, end, name)| {
            let (padded_start, padded_end) = pad_interval(start, end, window);

            // let mut corrected_reads_old: HashMap<String, Vec<(Vec<u8>, HashMap<Vec<u8>, f32>)>> = HashMap::new();
            // for window_start in (padded_start..padded_end).step_by(window) {

            let corrected_reads: HashMap<String, Vec<(Vec<u8>, HashMap<Vec<u8>, f32>)>> =
                (padded_start..padded_end)
                    .step_by(window)
                    .collect::<Vec<_>>()
                    .into_par_iter()
                    .filter_map(|window_start| {
                        let window_end = window_start + window as u64;

                        let mut locus = LinkedHashSet::new();
                        locus.insert((chrom.clone(), window_start, window_end, name.clone()));

                        let r = skydive::stage::stage_data_in_memory(
                            &locus,
                            &long_read_seq_urls,
                            false,
                            &cache_path,
                        );

                        if let Ok(reads) = r {
                            Some(reads)
                        } else {
                            None
                        }
                    })
                    .map(|reads| {
                        let long_reads: HashMap<String, Vec<u8>> = reads
                            .into_iter()
                            .map(|read| (read.id().to_string(), read.seq().to_vec()))
                            .collect();

                        let cn_kmers = long_reads
                            .values()
                            .map(|seq| seq.windows(kmer_size).collect::<Vec<_>>())
                            .flatten()
                            .map(|kmer| skydive::utils::canonicalize_kmer(kmer))
                            .collect::<HashSet<_>>();

                        let sr_seqs = all_sr_seqs
                            .iter()
                            .filter_map(|seq| {
                                let kmers = seq.windows(kmer_size).collect::<Vec<_>>();
                                let read_kmers = kmers
                                    .iter()
                                    .map(|kmer| skydive::utils::canonicalize_kmer(kmer))
                                    .collect::<HashSet<_>>();

                                let count = cn_kmers.intersection(&read_kmers).count();

                                if count > (0.5 * read_kmers.len() as f64) as usize {
                                    Some(seq)
                                } else {
                                    None
                                }
                            })
                            .cloned()
                            .collect::<Vec<_>>();

                        let lr_seqs = long_reads.values().cloned().collect::<Vec<_>>();
                        let l1 = LdBG::from_sequences("lr".to_string(), kmer_size, &lr_seqs);
                        let s1 = LdBG::from_sequences("sr".to_string(), kmer_size, &sr_seqs);

                        let m = MLdBG::from_ldbgs(vec![l1, s1])
                            .score_kmers(&model_absolute_path)
                            .collapse()
                            .clean(0.1)
                            .build_links(&lr_seqs, false);

                        let g = m.traverse_all_kmers();

                        let mut window_corrections = HashMap::new();
                        for (id, seq) in long_reads {
                            let corrected_seq = m.correct_seq(&g, &seq);
                            window_corrections
                                .entry(id)
                                .or_insert_with(Vec::new)
                                .push((corrected_seq, m.scores.clone()));
                        }

                        window_corrections
                    })
                    .reduce(
                        || HashMap::new(),
                        |mut acc, window_map| {
                            for (id, corrections) in window_map {
                                acc.entry(id).or_insert_with(Vec::new).extend(corrections);
                            }
                            acc
                        },
                    );

            let mut file = fa_file.lock().unwrap();
            for id in corrected_reads.keys() {
                let pieces = corrected_reads.get(id).unwrap();

                let mut joined_seq = Vec::<u8>::new();
                let mut joined_scores = HashMap::new();
                for (piece, scores) in pieces {
                    joined_seq.extend(piece);
                    joined_scores.extend(scores);
                }

                let mut prob = vec![1.0; joined_seq.len()];
                let mut count = vec![1; joined_seq.len()];
                for (i, kmer) in joined_seq.windows(kmer_size).enumerate() {
                    let cn_kmer = skydive::utils::canonicalize_kmer(kmer);
                    let score = **joined_scores.get(&cn_kmer).unwrap_or(&&1.0);

                    for j in i..std::cmp::min(i + kmer_size, prob.len()) {
                        prob[j] *= score;
                        count[j] += 1;
                    }
                }

                let mut qual = vec![0; joined_seq.len()];
                for i in 0..prob.len() {
                    prob[i] = prob[i].powf(1.0 / count[i] as f32);
                    qual[i] = ((-10.0 * (1.0 - (prob[i] - 0.0001).max(0.0)).log10()) as u8 + 33)
                        .clamp(33, 73);
                }

                let _ = writeln!(
                    file,
                    "@{}\n{}\n+\n{}",
                    id,
                    String::from_utf8(joined_seq).unwrap(),
                    String::from_utf8_lossy(&qual)
                );
            }
        });
}

fn pad_interval(start: &u64, end: &u64, window: usize) -> (u64, u64) {
    // Calculate how much padding is needed to make the interval cleanly divisible by 1000
    let interval_length = end - start;
    let remainder = interval_length % window as u64;

    let (padded_start, padded_end) = if remainder > 0 {
        let padding_needed = window as u64 - remainder;
        let half_padding = padding_needed / 2;

        // Add padding evenly to start and end
        let new_start = start.saturating_sub(half_padding);
        let new_end = end + (padding_needed - half_padding);

        (new_start, new_end)
    } else {
        (*start, *end)
    };

    (padded_start, padded_end)
}
