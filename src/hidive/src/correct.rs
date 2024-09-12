use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::path::PathBuf;
use std::{
    fs::File,
    io::{BufWriter, Write},
};

use bio::alignment::pairwise::Scoring;
use bio::alignment::poa::Aligner;
use indicatif::ParallelProgressIterator;
use parquet::data_type::AsBytes;
use petgraph::dot::Dot;
use petgraph::visit::EdgeRef;
use petgraph::Graph;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;

use itertools::Itertools;

use spoa::AlignmentType;

pub fn start(
    output: &PathBuf,
    long_read_fasta_paths: &Vec<PathBuf>,
    short_read_fasta_paths: &Vec<PathBuf>,
) {
    let long_read_seq_urls = skydive::parse::parse_file_names(long_read_fasta_paths);
    let short_read_seq_urls = skydive::parse::parse_file_names(short_read_fasta_paths);

    // Read all long reads.
    skydive::elog!(
        "Processing long-read samples {:?}...",
        long_read_seq_urls
            .iter()
            .map(|url| url.as_str())
            .collect::<Vec<&str>>()
    );
    let all_lr_seqs = long_read_fasta_paths
        .iter()
        .map(|p| {
            let reader = bio::io::fasta::Reader::from_file(p).expect("Failed to open file");
            reader
                .records()
                .filter_map(|r| r.ok())
                .map(|r| r.seq().to_vec())
                .collect::<Vec<Vec<u8>>>()
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

    skydive::elog!(
        "Processing short-read samples {:?}...",
        short_read_seq_urls
            .iter()
            .map(|url| url.as_str())
            .collect::<Vec<&str>>()
    );
    let all_sr_seqs = short_read_fasta_paths
        .iter()
        .map(|p| {
            let reader = bio::io::fasta::Reader::from_file(p).expect("Failed to open file");
            reader
                .records()
                .filter_map(|r| r.ok())
                .map(|r| r.seq().to_vec())
                .collect::<Vec<Vec<u8>>>()
        })
        .flatten()
        .collect::<Vec<Vec<u8>>>();

    let mut sa = spoa::AlignmentEngine::new(AlignmentType::kSW, 5, -10, -16, -12, -20, -8);

    let progress_bar = skydive::utils::default_bounded_progress_bar(
        "Processing short reads",
        all_sr_seqs.len() as u64,
    );
    for sr_seq in all_sr_seqs {
        let seq_cstr = std::ffi::CString::new(sr_seq.clone()).unwrap();
        let seq_qual = std::ffi::CString::new(vec![b'I'; sr_seq.len()]).unwrap();
        let a = sa.align(seq_cstr.as_ref(), &sg);

        sg.add_alignment(&a, seq_cstr.as_ref(), seq_qual.as_ref());
        progress_bar.inc(1);
    }

    let msa_cstrs = sg.multiple_sequence_alignment(false);
    let msa_strings = msa_cstrs
        .iter()
        .map(|cstr| cstr.to_str().unwrap().to_string())
        .collect::<Vec<String>>();

    let mut lr_msas = msa_strings
        .iter()
        .take(all_lr_seqs.len())
        .map(|msa| trim_unaligned_bases(msa))
        .collect::<Vec<String>>();
    let sr_msas = msa_strings
        .iter()
        .skip(all_lr_seqs.len())
        .cloned()
        .collect::<Vec<String>>();

    // let leftmost_base_position = lr_msas
    //     .iter()
    //     .map(|msa| {
    //         msa.find(|c| c == 'A' || c == 'C' || c == 'G' || c == 'T')
    //             .unwrap_or(lr_msas.first().unwrap().len())
    //     })
    //     .min()
    //     .unwrap_or(0);

    let rightmost_base_position = lr_msas
        .iter()
        .map(|msa| {
            msa.rfind(|c| c == 'A' || c == 'C' || c == 'G' || c == 'T')
                .unwrap_or(0)
        })
        .max()
        .unwrap_or(0);

    progress_bar.reset();
    progress_bar.set_message("Filtering short reads");
    progress_bar.set_length(sr_msas.len() as u64);

    let mut filtered_sr_msas = sr_msas
        .par_iter()
        .progress_with(progress_bar)
        .filter_map(|sr_msa| {
            let num_good_bases = sr_msa
                .chars()
                .enumerate()
                .filter_map(|(i, sr_char)| {
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
                })
                .sum::<usize>();

            let len = sr_msa.replace("-", "").len();
            let score = 100.0 * num_good_bases as f64 / len as f64;

            if score > 90.0 {
                Some(trim_n_bases(
                    &trim_unaligned_bases(&sr_msa[..=rightmost_base_position].to_string()),
                    10,
                ))
            } else {
                None
            }
        })
        .collect::<Vec<String>>();

    // lr_msas.iter().for_each(|msa| {
    //     println!("{}", msa);
    // });

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
            let lr_mod_char = lr_chars
                .iter()
                .filter(|(_, &count)| count == 1)
                .map(|(&ch, _)| ch)
                .collect::<Vec<char>>();

            if lr_mod_char.len() > 0 {
                let lr_mod_char_not_in_sr = lr_mod_char
                    .iter()
                    .filter(|&ch| !sr_chars.contains_key(ch))
                    .cloned()
                    .collect::<Vec<char>>();

                if lr_mod_char_not_in_sr.len() > 0 {
                    // println!("{} {:?} {:?} {:?} {:?}", column, lr_chars, lr_mod_char, sr_chars, lr_mod_char_not_in_sr);

                    for base in lr_mod_char_not_in_sr {
                        for lr_msa in lr_msas.iter_mut() {
                            if let Some(ch) = lr_msa.chars().nth(column) {
                                if ch == base {
                                    let mut chars: Vec<char> = lr_msa.chars().collect();

                                    let most_common_sr_base = sr_chars
                                        .iter()
                                        .max_by_key(|&(_, count)| count)
                                        .map(|(base, _)| *base);
                                    if let Some(sr_base) = most_common_sr_base {
                                        chars[column] = sr_base;
                                    } else {
                                        chars[column] = '.';
                                    }

                                    *lr_msa = chars.into_iter().collect();
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    filtered_sr_msas.sort_by(|a, b| b.cmp(a));

    lr_msas.iter().for_each(|msa| {
        println!("{:5} {}", "", msa);
    });

    filtered_sr_msas.iter().for_each(|msa| {
        println!("{:5} {}", "", msa);
    });

    let output_file = File::create(output).unwrap();
    let mut writer = BufWriter::new(output_file);

    // lr_msas.iter().enumerate().for_each(|(i, msa)| {
    //     let (new_msa1, msas1) = recruit(&msa.to_lowercase(), &filtered_sr_msas, 149);
    //     let (new_msa2, _) = recruit(&msa.to_lowercase(), &msas1, 148);

    //     // let _ = writeln!(writer, ">{}_{}\n{}", i, 0, msa.chars().filter(|&c| c != '-' && c != ' ').collect::<String>());
    //     // let _ = writeln!(writer, ">{}_{}\n{}", i, 1, new_msa1.chars().filter(|&c| c != '-' && c != ' ').collect::<String>());
    //     // let _ = writeln!(writer, ">{}_{}\n{}", i, 2, new_msa2.chars().filter(|&c| c != '-' && c != ' ').collect::<String>());

    //     let _ = writeln!(writer, ">{}_{}\n{}", i, 0, msa);
    //     let _ = writeln!(writer, ">{}_{}\n{}", i, 1, new_msa1);
    //     let _ = writeln!(writer, ">{}_{}\n{}", i, 2, new_msa2);
    // });

    let lr_clusters = cluster_sequences(&lr_msas, &lr_msas, 200, 0.90);
    let lr_groups = group_clusters(&lr_clusters);

    for (lr_group_index, lr_group) in lr_groups.iter().enumerate() {
        println!("Cluster: {:?}", lr_group);

        for &index in lr_group {
            println!("  MSA {:5}: {}", index, lr_msas[index]);
        }
        println!();

        let group_lr_msas = lr_group.iter().map(|i| lr_msas[*i].clone()).collect::<Vec<String>>();
        let sr_clusters = cluster_sequences(&group_lr_msas, &filtered_sr_msas, 100, 0.90);
        let sr_groups = group_clusters(&sr_clusters);

        for sr_group in sr_groups {
            println!("  Group: {:?}", sr_group);

            let group_sr_msas = sr_group.iter().map(|i| filtered_sr_msas[*i].clone()).collect::<Vec<String>>();

            for (i, group_sr_msa) in group_sr_msas.iter().enumerate() {
                println!("  msa {:5}: {}", i, group_sr_msa);
            }
            println!();

            let g = create_graph(&group_lr_msas, &group_sr_msas, 0);
            let hap1 = haplotype(&g, true);
            let hap2 = haplotype(&g, false);

            let scoring = Scoring::from_scores(-8, -6, 5, -4);
            let mut aligner1 = Aligner::new(scoring, hap1.as_bytes());
            let mut aligner2 = Aligner::new(scoring, hap2.as_bytes());

            let total_score = all_lr_seqs
                .iter()
                .map(|lr_seq| {
                    let score1 = aligner1.local(&lr_seq.as_bytes()).alignment().score;
                    let score2 = aligner2.local(&lr_seq.as_bytes()).alignment().score;

                    println!("scores: {} {}", score1, score2);

                    score1 + score2
                })
                .sum::<i32>();

            println!("Total score: {}", total_score);

            println!("Consensus 1: {}", hap1);
            println!("Consensus 2: {}", hap2);

            let dot = Dot::new(&g);
            let mut dot_file = File::create(format!("graph_{}.dot", lr_group_index)).unwrap();
            write!(dot_file, "{:?}", dot).unwrap();

            // Write graph as GFA to a string
            let mut gfa_output = Vec::new();
            skydive::utils::write_gfa(&mut gfa_output, &g).unwrap();

            // Print GFA string (commented out for test)
            let gfa_string = String::from_utf8(gfa_output).unwrap();
            let _ = writeln!(writer, "{}", gfa_string);
        }
    }
}

fn haplotype(g: &Graph<String, f32>, hap1: bool) -> String {
    // Find the start node (a node with no incoming edges)
    let start_node = g.node_indices().find(|&n| g.edges_directed(n, petgraph::Direction::Incoming).count() == 0)
        .expect("No start node found");

    // Initialize the path with the start node
    let mut path = vec![start_node];
    let mut current_node = start_node;

    // Traverse the graph
    while let Some(next_edge) = g.edges_directed(current_node, petgraph::Direction::Outgoing)
        .filter(|e| if hap1 { *g.edge_weight(e.id()).unwrap() > -0.1 } else { *g.edge_weight(e.id()).unwrap() < 0.1 })
        .next()
    {
        let next_node = next_edge.target();
        path.push(next_node);
        current_node = next_node;

        // Break if we've reached an end node (no outgoing edges)
        if g.edges_directed(current_node, petgraph::Direction::Outgoing).count() == 0 {
            break;
        }
    }

    // Construct the consensus sequence from the path
    let consensus = path.iter()
        .map(|&node| g.node_weight(node).unwrap().clone())
        .collect::<String>();
    consensus
}

fn create_graph(lr_msas: &Vec<String>, sr_msas: &Vec<String>, seed: u64) -> Graph<String, f32> {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

    let mut graph = Graph::new();

    let mut prev_nodes = Vec::new();

    let mut contig = String::new();
    let mut index1 = 0;
    while index1 < lr_msas[0].len() {
        let combined_base_counts = pileup(lr_msas, sr_msas, index1);

        if combined_base_counts.len() == 0 {
            // contig.push('.');

            index1 += 1;
        } else if combined_base_counts.len() == 1 {
            let base = combined_base_counts.keys().next().unwrap();
            contig.push(*base);

            index1 += 1;
        } else {
            let node = graph.add_node(contig.clone());

            prev_nodes.iter().for_each(|prev_node| {
                graph.add_edge(*prev_node, node, 0.0);
            });
            prev_nodes.clear();

            let mut index2 = index1;
            let mut allele_base_counts = pileup(lr_msas, sr_msas, index2);

            while index2 < lr_msas[0].len() && allele_base_counts.len() > 1 {
                index2 += 1;
                allele_base_counts = pileup(lr_msas, sr_msas, index2);
            }

            let allele_counts = alleles(lr_msas, sr_msas, index1, index2)
                .into_iter()
                .filter(|&(_, count)| count >= 5)
                .sorted_by(|a, b| b.1.cmp(&a.1))
                .take(2)
                .collect::<BTreeMap<_, _>>();

            let filtered_allele_counts = allele_counts.into_iter().filter(|(_, count)| *count > 5).collect::<BTreeMap<_, _>>();

            if filtered_allele_counts.len() == 1 {
                for (allele, _) in filtered_allele_counts {
                    let new_node = graph.add_node(allele.clone());

                    graph.add_edge(node, new_node,  0.0);
                    prev_nodes.push(new_node);
                }
            } else {
                let mut phase = false;
                for (allele, _) in filtered_allele_counts {
                    let new_node = graph.add_node(allele.clone());

                    // let phase = rng.gen_bool(0.5);
                    graph.add_edge(node, new_node, if phase { 1.0 } else { -1.0 });
                    phase = !phase;

                    prev_nodes.push(new_node);
                }
            }

            index1 = index2;

            contig = String::new();
        }
    }

    graph
}

fn pileup(lr_msas: &Vec<String>, sr_msas: &Vec<String>, index: usize) -> BTreeMap<char, i32> {
    let combined_base_counts = lr_msas.iter().chain(sr_msas.iter())
        .map(|msa| msa.chars().nth(index).unwrap_or(' '))
        .filter(|&c| c != ' ')
        .fold(BTreeMap::new(), |mut counts, base| {
            *counts.entry(base).or_insert(0) += 1;
            counts
        });
    combined_base_counts
}

fn alleles(lr_msas: &Vec<String>, sr_msas: &Vec<String>, index1: usize, index2: usize) -> BTreeMap<String, i32> {
    let combined_allele_counts = lr_msas.iter().chain(sr_msas.iter())
        .map(|msa| msa[index1..index2].to_string().replace(" ", ""))
        .filter(|allele| allele.len() > 0)
        .fold(BTreeMap::new(), |mut counts, base| {
            *counts.entry(base).or_insert(0) += 1;
            counts
        });
    combined_allele_counts
}

fn group_clusters(clusters: &Vec<Vec<bool>>) -> Vec<BTreeSet<usize>> {
    let mut groups = Vec::new();
    let mut used = HashSet::new();
    for (i, cluster) in clusters.iter().enumerate() {
        if !used.contains(&i) {
            let mut group = BTreeSet::new();
            group.insert(i);

            for (j, similar) in cluster.iter().enumerate() {
                if i < j && *similar {
                    group.insert(j);
                }
            }

            used.extend(group.clone());
            groups.push(group);
        }
    }

    groups
}

fn cluster_sequences(msas1: &Vec<String>, msas2: &Vec<String>, min_overlap: usize, min_similarity: f64) -> Vec<Vec<bool>> {
    // let n = msas1.len();
    let mut similarity_matrix = vec![vec![false; msas2.len()]; msas1.len()];

    for (i, lr_msa1) in msas1.iter().enumerate() {
        for (j, lr_msa2) in msas2.iter().enumerate() {
            let mut matching_positions = 0;
            let mut total_positions = 0;

            for (c1, c2) in lr_msa1.chars().zip(lr_msa2.chars()) {
                if c1 != ' ' && c2 != ' ' {
                    total_positions += 1;
                    if c1 == c2 {
                        matching_positions += 1;
                    }
                }
            }

            let similarity = if total_positions > 0 {
                matching_positions as f64 / total_positions as f64
            } else {
                0.0
            };

            // println!("{} {} {} {}", i, j, matching_positions, similarity);
            similarity_matrix[i][j] = (matching_positions >= min_overlap) && (similarity >= min_similarity);

            // similarity_matrix[j][i] = (matching_positions >= min_overlap) && (similarity >= min_similarity);
        }
    }

    similarity_matrix
}

fn allele_counts(msas: &Vec<String>, column: usize) -> BTreeMap<char, usize> {
    let mut counts = BTreeMap::new();

    for msa in msas {
        let char = msa.chars().nth(column).unwrap();
        if char != ' ' {
            *counts.entry(char).or_insert(0) += 1;
        }
    }

    counts
}

fn recruit(lr_msa: &String, sr_msas: &Vec<String>, min_score: usize) -> (String, Vec<String>) {
    let scores = sr_msas.iter().map(|sr_msa| similarity(lr_msa, sr_msa)).collect::<Vec<usize>>();
    let mut score_msa_pairs: Vec<_> = scores.into_iter().zip(sr_msas.iter()).collect();
    score_msa_pairs.sort_by(|a, b| b.0.cmp(&a.0));

    let sorted_scores: Vec<usize> = score_msa_pairs.iter().map(|(score, _)| *score).collect();
    let sorted_sr_msas: Vec<&String> = score_msa_pairs.iter().map(|(_, msa)| *msa).collect();

    let mut msas = vec![];
    let mut msa_scores = vec![];
    for (score, sr_msa) in sorted_scores.iter().zip(sorted_sr_msas) {
        // println!("{:5} {}", *score, sr_msa);

        if *score >= min_score {
            msas.push(sr_msa.clone());
            msa_scores.push(*score);
        }
    }

    // Sort msa_scores and msas reverse-lexicographically based on msas
    let mut score_msa_pairs: Vec<_> = msa_scores.into_iter().zip(msas.iter_mut()).collect();
    score_msa_pairs.sort_by(|(_, a), (_, b)| b.cmp(a));
    
    // for (score, msa) in &score_msa_pairs {
    //     println!("{:5} {}", *score, msa);
    // }

    let mut used_positions: BTreeMap<usize, Vec<(usize, &String)>> = BTreeMap::new();
    for (score, msa) in &score_msa_pairs {
        if let Some(first_valid_pos) = msa.find(|c| c == 'A' || c == 'C' || c == 'G' || c == 'T') {
            if !used_positions.contains_key(&first_valid_pos) {
                used_positions.insert(first_valid_pos, vec![]);
            }

            used_positions.get_mut(&first_valid_pos).unwrap().push((*score, msa));
        }
    }

    let mut used = HashSet::new();

    let mut new_msa = lr_msa.clone();
    for (_, pos_vec) in used_positions {
        let mut pos_vec = pos_vec.into_iter().collect::<Vec<_>>();
        pos_vec.sort_by(|a, b| b.0.cmp(&a.0));

        let best_sr_msa = pos_vec.first().unwrap().1;

        for (i, c) in best_sr_msa.chars().enumerate() {
            if c != '-' && c != ' ' && matches!(new_msa.chars().nth(i), Some('a' | 'c' | 'g' | 't')) {
                new_msa.replace_range(i..=i, c.to_string().as_str());
            }
        }

        used.insert(best_sr_msa);
    }

    // println!("{:5} {}", "", new_msa);

    // Remove from sr_msas all entries in used
    let sr_msas: Vec<String> = sr_msas.into_iter().filter(|msa| !used.contains(msa)).cloned().collect();

    (new_msa, sr_msas)
}

fn similarity(lr_msa: &String, sr_msa: &String) -> usize {
    lr_msa.chars().zip(sr_msa.chars())
        .filter(|(a, b)| {
            matches!(a.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T') &&
            matches!(b.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T') &&
            a.to_ascii_uppercase() == b.to_ascii_uppercase()
        })
        .count()
}

fn trim_unaligned_bases(aligned_seq: &String) -> String {
    let trimmed_sr_msa_left: String = aligned_seq
        .chars()
        .enumerate()
        .map(|(i, c)| {
            if i < aligned_seq
                .find(|ch: char| ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T')
                .unwrap_or(aligned_seq.len())
            {
                ' '
            } else {
                c
            }
        })
        .collect();

    let trimmed_sr_msa: String = trimmed_sr_msa_left
        .chars()
        .rev()
        .enumerate()
        .map(|(i, c)| {
            if i < trimmed_sr_msa_left
                .chars()
                .rev()
                .position(|ch| ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T')
                .unwrap_or(trimmed_sr_msa_left.len())
            {
                ' '
            } else {
                c
            }
        })
        .collect::<String>()
        .chars()
        .rev()
        .collect();

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
