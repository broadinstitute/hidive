use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::io::Read;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::{
    fs::File,
    io::{BufWriter, Write},
};

use flate2::read::GzDecoder;
use indicatif::ParallelProgressIterator;
use num_format::{Locale, ToFormattedString};

use rayon::prelude::*;

use spoa::AlignmentType;

pub fn start(
    output: &PathBuf,
    long_read_fasta_paths: &Vec<PathBuf>,
    short_read_fasta_paths: &Vec<PathBuf>,
) {
    let long_read_seq_urls = skydive::parse::parse_file_names(long_read_fasta_paths);
    let short_read_seq_urls = skydive::parse::parse_file_names(short_read_fasta_paths);

    // Read all long reads.
    skydive::elog!("Processing long-read samples {:?}...", long_read_seq_urls.iter().map(|url| url.as_str()).collect::<Vec<&str>>());
    let all_lr_seqs = long_read_fasta_paths
        .iter()
        .map(|p| {
            let reader = bio::io::fasta::Reader::from_file(p).expect("Failed to open file");
            reader.records().filter_map(|r| r.ok()).map(|r| r.seq().to_vec()).collect::<Vec<Vec<u8>>>()
        })
        .flatten()
        .collect::<Vec<Vec<u8>>>();

    let mut la = spoa::AlignmentEngine::new(AlignmentType::kOV, 5, -4, -8, -6, -10, -4);

    let mut sg = spoa::Graph::new();
    for lr_seq in &all_lr_seqs {
        let seq_cstr = std::ffi::CString::new(lr_seq.clone()).unwrap();
        let seq_qual = std::ffi::CString::new(vec![b'I'; lr_seq.len()]).unwrap();
        let a = la.align(seq_cstr.as_ref(), &sg);

        sg.add_alignment(&a, seq_cstr.as_ref(), seq_qual.as_ref());
    }

    skydive::elog!("Processing short-read samples {:?}...", short_read_seq_urls.iter().map(|url| url.as_str()).collect::<Vec<&str>>());
    let all_sr_seqs = short_read_fasta_paths
        .iter()
        .map(|p| {
            let reader = bio::io::fasta::Reader::from_file(p).expect("Failed to open file");
            reader.records().filter_map(|r| r.ok()).map(|r| r.seq().to_vec()).collect::<Vec<Vec<u8>>>()
        })
        .flatten()
        .collect::<Vec<Vec<u8>>>();

    let mut sa = spoa::AlignmentEngine::new(AlignmentType::kSW, 5, -10, -16, -12, -20, -8);

    let progress_bar = skydive::utils::default_bounded_progress_bar("Processing short reads", all_sr_seqs.len() as u64);
    for sr_seq in all_sr_seqs {
        let seq_cstr = std::ffi::CString::new(sr_seq.clone()).unwrap();
        let seq_qual = std::ffi::CString::new(vec![b'I'; sr_seq.len()]).unwrap();
        let a = sa.align(seq_cstr.as_ref(), &sg);

        sg.add_alignment(&a, seq_cstr.as_ref(), seq_qual.as_ref());
        progress_bar.inc(1);
    }

    let msa_cstrs = sg.multiple_sequence_alignment(false);
    let msa_strings = msa_cstrs.iter().map(|cstr| cstr.to_str().unwrap().to_string()).collect::<Vec<String>>();

    let mut lr_msas = msa_strings.iter().take(all_lr_seqs.len()).map(|msa| trim_unaligned_bases(msa)).collect::<Vec<String>>();
    let sr_msas = msa_strings.iter().skip(all_lr_seqs.len()).cloned().collect::<Vec<String>>();

    let rightmost_base_position = lr_msas.iter().map(|msa| {
        msa.rfind(|c| c == 'A' || c == 'C' || c == 'G' || c == 'T').unwrap_or(0)
    }).max().unwrap_or(0);

    progress_bar.reset();
    progress_bar.set_message("Filtering short reads");
    progress_bar.set_length(sr_msas.len() as u64);

    let mut filtered_sr_msas = sr_msas.par_iter().progress_with(progress_bar).filter_map(|sr_msa| {
        let num_good_bases = sr_msa.chars().enumerate().filter_map(|(i, sr_char)| {
            let lr_chars = &lr_msas
                .iter()
                .map(|lr_msa| lr_msa.chars().nth(i).unwrap())
                .filter(|c| c.is_ascii_alphabetic())
                .collect::<HashSet<char>>();

            if lr_chars.contains(&sr_char) {
                Some(1)
            } else {
                Some(0)
            }
        }).sum::<usize>();

        let len = sr_msa.replace("-", "").len();
        let score = 100.0 * num_good_bases as f64 / len as f64;

        if score > 90.0 {
            Some(trim_n_bases(&trim_unaligned_bases(&sr_msa[..=rightmost_base_position].to_string()), 10))
        } else {
            None
        }
    }).collect::<Vec<String>>();

    lr_msas.iter().for_each(|msa| {
        println!("{}", msa);
    });

    let length = filtered_sr_msas.first().unwrap().len();

    for column in 0..length {
        let mut lr_chars = BTreeMap::new();
        let mut sr_chars = BTreeMap::new();

        for row in 0..lr_msas.len() {
            let char = lr_msas[row].chars().nth(column).unwrap();
            *lr_chars.entry(char).or_insert(0) += 1;
        }

        lr_chars.remove(&' ');

        for row in 0..filtered_sr_msas.len() {
            let char = filtered_sr_msas[row].chars().nth(column).unwrap();
            *sr_chars.entry(char).or_insert(0) += 1;
        }

        sr_chars.remove(&' ');

        if lr_chars.len() > 0 && sr_chars.len() > 0 {
            let lr_mod_char = lr_chars.iter()
                .filter(|(_, &count)| count == 1)
                .map(|(&ch, _)| ch)
                .collect::<Vec<char>>();

            if lr_mod_char.len() > 0 {
                let lr_mod_char_not_in_sr = lr_mod_char.iter()
                    .filter(|&ch| !sr_chars.contains_key(ch))
                    .cloned()
                    .collect::<Vec<char>>();

                if lr_mod_char_not_in_sr.len() > 0 {
                    println!("{} {:?} {:?} {:?} {:?}", column, lr_chars, lr_mod_char, sr_chars, lr_mod_char_not_in_sr);

                    for base in lr_mod_char_not_in_sr {
                        for lr_msa in lr_msas.iter_mut() {
                            if let Some(ch) = lr_msa.chars().nth(column) {
                                if ch == base {
                                    let mut chars: Vec<char> = lr_msa.chars().collect();

                                    let most_common_sr_base = sr_chars.iter()
                                        .max_by_key(|&(_, count)| count)
                                        .map(|(base, _)| *base);
                                    if let Some(sr_base) = most_common_sr_base {
                                        chars[column] = sr_base;
                                    } else {
                                        chars[column] = '.';
                                    }

                                    // println!(" --> {} {} {}", chars[column], ch, base);

                                    // chars[column] = '.';
                                    *lr_msa = chars.into_iter().collect();

                                }
                            }
                        }
                    }
                }
            }
        }
    }

    println!();

    lr_msas.iter().for_each(|msa| {
        println!("{}", msa);
    });

    filtered_sr_msas.sort_by(|a, b| b.cmp(a));

    filtered_sr_msas.iter().for_each(|msa| {
        println!("{}", msa);
    });

    let output_file = File::create(output).unwrap();
    let mut writer = BufWriter::new(output_file);

    lr_msas.iter().enumerate().for_each(|(i, msa)| {
        let filtered_msa: String = msa.chars()
            .filter(|&c| c != '-' && c != ' ')
            .collect();

        let _ = writeln!(writer, ">{}\n{}", i, filtered_msa);
    });
}

fn trim_unaligned_bases(aligned_seq: &String) -> String {
    let trimmed_sr_msa_left: String = aligned_seq.chars().enumerate()
        .map(|(i, c)| {
            if i < aligned_seq.find(|ch: char| ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T').unwrap_or(aligned_seq.len()) {
                ' '
            } else {
                c
            }
        })
        .collect();

    let trimmed_sr_msa: String = trimmed_sr_msa_left.chars().rev()
        .enumerate()
        .map(|(i, c)| {
            if i < trimmed_sr_msa_left.chars().rev()
                .position(|ch| ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T')
                .unwrap_or(trimmed_sr_msa_left.len())
            {
                ' '
            } else {
                c
            }
        })
        .collect::<String>()
        .chars().rev().collect();

    trimmed_sr_msa
}

fn trim_n_bases(aligned_seq: &String, num_bases: usize) -> String {
    let is_base = |c: char| matches!(c, 'A' | 'C' | 'G' | 'T');
    let replace_bases = |s: &str| {
        s.chars()
            .scan(0, |count, c| {
                if *count < num_bases {
                    if is_base(c) {
                        *count += 1;
                    }
                    Some(c)
                } else {
                    Some(c)
                }
            })
            .collect::<String>()
    };

    let left_trimmed = replace_bases(aligned_seq);
    let right_trimmed = replace_bases(&left_trimmed.chars().rev().collect::<String>());
    right_trimmed.chars().rev().collect()
}
