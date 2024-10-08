use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::path::PathBuf;

use indicatif::ParallelProgressIterator;
use petgraph::graph::NodeIndex;
use petgraph::Graph;
use rayon::prelude::*;

use itertools::Itertools;

use skydive::wmec::*;
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

    // let output_file = File::create(output).unwrap();
    // let mut writer = BufWriter::new(output_file);

    let lr_clusters = cluster_sequences(&lr_msas, &lr_msas, 200, 0.90);
    let lr_groups = group_clusters(&lr_clusters);

    for (lr_group_index, lr_group) in lr_groups.iter().enumerate() {
        // println!("Cluster: {:?}", lr_group);

        // for &index in lr_group {
        //     println!("  MSA {:5}: {}", index, lr_msas[index]);
        // }
        // println!();

        let group_lr_msas = lr_group.iter().map(|i| lr_msas[*i].clone()).collect::<Vec<String>>();
        let sr_clusters = cluster_sequences(&group_lr_msas, &filtered_sr_msas, 100, 0.90);
        let sr_groups = group_clusters(&sr_clusters);

        for sr_group in sr_groups {
            // println!("  Group: {:?}", sr_group);

            let group_sr_msas = sr_group.iter().map(|i| filtered_sr_msas[*i].clone()).collect::<Vec<String>>();

            // for (i, group_sr_msa) in group_sr_msas.iter().enumerate() {
            //     println!("  msa {:5}: {}", i, group_sr_msa);
            // }
            // println!();

            let g = create_graph(&group_lr_msas, &group_sr_msas);

            let m = create_read_allele_matrix(&group_lr_msas, &group_sr_msas);

            if m.len() > 0 {
                let wmec = create_wmec_matrix(&m);

                let (h1, h2) = phase(&wmec);
                let (hap1, hap2) = create_fully_phased_haplotypes(&group_lr_msas, &group_sr_msas, &h1);

                println!(">cluster_{}_hap1", lr_group_index);
                println!("{}", hap1.replace("-", ""));
                println!(">cluster_{}_hap2", lr_group_index);
                println!("{}", hap2.replace("-", ""));
            }

            /*
            // Print matrix m as a TSV file
            let tsv_file = File::create(format!("matrix_cluster_{}.tsv", lr_group_index)).unwrap();
            let mut tsv_writer = BufWriter::new(tsv_file);

            // Write header row with variant indices
            writeln!(tsv_writer, "read_index\t{}", (0..m.len()).join("\t")).unwrap();

            // Collect all unique read indices
            let all_read_indices: BTreeSet<_> = m.iter().flat_map(|col| col.keys()).collect();

            // Write each row
            for &read_index in &all_read_indices {
                write!(tsv_writer, "{}", read_index).unwrap();
                for col in &m {
                    write!(tsv_writer, "\t{}", col.get(read_index).unwrap_or(&String::from("."))).unwrap();
                }
                writeln!(tsv_writer).unwrap();
            }

            skydive::elog!("Creating het graph");
            let (h, source_node, sink_node) = create_het_graph(&group_lr_msas, &group_sr_msas);

            if lr_group_index > 0 {
                skydive::elog!("Writing het graph");

                // Write graph h to a dot file
                let dot_file = File::create(format!("graph_h_cluster_{}.dot", lr_group_index)).unwrap();
                let mut dot_writer = BufWriter::new(dot_file);
                write!(dot_writer, "{:?}", Dot::new(&h)).unwrap();
            }

            let dot = Dot::new(&g);
            let mut dot_file = File::create(format!("graph_{}.dot", lr_group_index)).unwrap();
            write!(dot_file, "{:?}", dot).unwrap();

            // Write graph as GFA to a string
            // let mut gfa_output = Vec::new();
            // skydive::utils::write_gfa(&mut gfa_output, &g).unwrap();

            // Print GFA string (commented out for test)
            // let gfa_string = String::from_utf8(gfa_output).unwrap();
            // let _ = writeln!(writer, "{}", gfa_string);
            */
        }
    }
}

fn create_read_allele_matrix(lr_msas: &Vec<String>, sr_msas: &Vec<String>) -> Vec<BTreeMap<usize, String>> {
    let mut matrix = Vec::new();

    let mut index1 = 0;
    while index1 < lr_msas[0].len() {
        let combined_base_counts = pileup_counts(lr_msas, sr_msas, index1);

        if combined_base_counts.len() > 1 {
            let mut index2 = index1;
            let mut allele_base_counts = pileup_counts(lr_msas, sr_msas, index2);

            while index2 < lr_msas[0].len() && allele_base_counts.len() > 1 {
                index2 += 1;
                allele_base_counts = pileup_counts(lr_msas, sr_msas, index2);
            }

            let allele_counts = allele_counts(lr_msas, sr_msas, index1, index2)
                .into_iter()
                .filter(|&(_, count)| count >= 5)
                .sorted_by(|a, b| b.1.cmp(&a.1))
                .take(2)
                .collect::<BTreeMap<_, _>>();

            let filtered_allele_counts = allele_counts.into_iter().filter(|(_, count)| *count > 5).collect::<BTreeMap<_, _>>();

            if filtered_allele_counts.len() == 2 {
                let mut column = BTreeMap::new();

                let alleles = get_allele_indices(lr_msas, sr_msas, index1, index2);
                let lr_alleles = get_allele_indices(lr_msas, &Vec::new(), index1, index2);

                let mut allele_index = 0;
                for (allele, _) in &filtered_allele_counts {
                    lr_alleles.iter().enumerate().for_each(|(i, a)| {
                        if *a == *allele {
                            // column.insert(i, allele.clone());
                            column.insert(i, String::from(allele_index.to_string()));
                        }
                    });

                    allele_index += 1;
                }

                matrix.push(column);
            }

            index1 = index2;
        } else {
            index1 += 1;
        }
    }

    matrix
}

fn create_wmec_matrix(matrix: &Vec<BTreeMap<usize, String>>) -> WMECData {
    let num_snps = matrix.len();
    let num_reads = matrix.iter().map(|m| m.keys().max().unwrap_or(&0) + 1).max().unwrap_or(0);

    let mut reads = vec![vec![None; num_snps]; num_reads];
    let mut confidences = vec![vec![None; num_snps]; num_reads];

    for (snp_idx, column) in matrix.iter().enumerate() {
        for (&read_idx, allele) in column {
            reads[read_idx][snp_idx] = Some(allele.parse::<u8>().unwrap());
            confidences[read_idx][snp_idx] = Some(32);
        }
    }

    WMECData::new(reads, confidences)
}

fn create_het_graph(lr_msas: &Vec<String>, sr_msas: &Vec<String>) -> (Graph<String, usize>, NodeIndex, NodeIndex) {
    let mut graph = Graph::new();

    let mut prev_nodes = Vec::new();
    let mut prev_read_sets= Vec::new();

    // Add a source dummy node to the graph
    let source_node = graph.add_node(String::from("^"));
    prev_nodes.push(source_node);
    // prev_read_sets.push((0..lr_msas.len()).collect::<HashSet<_>>());
    prev_read_sets.push((0..(lr_msas.len() + sr_msas.len())).collect::<HashSet<_>>());

    let mut index1 = 0;
    while index1 < lr_msas[0].len() {
        let combined_base_counts = pileup_counts(lr_msas, sr_msas, index1);

        if combined_base_counts.len() > 1 {
            let mut index2 = index1;
            let mut allele_base_counts = pileup_counts(lr_msas, sr_msas, index2);

            while index2 < lr_msas[0].len() && allele_base_counts.len() > 1 {
                index2 += 1;
                allele_base_counts = pileup_counts(lr_msas, sr_msas, index2);
            }

            let allele_counts = allele_counts(lr_msas, sr_msas, index1, index2)
                .into_iter()
                .filter(|&(_, count)| count >= 5)
                .sorted_by(|a, b| b.1.cmp(&a.1))
                .take(2)
                .collect::<BTreeMap<_, _>>();

            let filtered_allele_counts = allele_counts.into_iter().filter(|(_, count)| *count > 5).collect::<BTreeMap<_, _>>();

            if filtered_allele_counts.len() == 2 {
                let alleles = get_allele_indices(lr_msas, sr_msas, index1, index2);

                let mut new_nodes = Vec::new();
                let mut new_read_sets = Vec::new();

                for (allele, _) in &filtered_allele_counts {
                    let node = graph.add_node(allele.clone());

                    let allele_indices = &alleles.iter().enumerate().filter_map(|(i, a)| {
                        if *a == *allele {
                            Some(i)
                        } else {
                            None
                        }
                    })
                    .collect::<HashSet<_>>();

                    let lr_alleles = get_allele_indices(lr_msas, &Vec::new(), index1, index2);
                    let lr_allele_indices = lr_alleles.iter().enumerate().filter_map(|(i, a)| {
                        if *a == *allele {
                            Some(i)
                        } else {
                            None
                        }
                    })
                    .collect::<HashSet<_>>();
                    
                    for (prev_node, prev_read_set) in prev_nodes.iter().zip(&prev_read_sets) {
                        // let count = prev_read_set.intersection(&lr_allele_indices).count();
                        let count = prev_read_set.intersection(&allele_indices).count();
                        if count > 0 {
                            graph.add_edge(*prev_node, node, count);
                        }
                    }

                    new_nodes.push(node);
                    // new_read_sets.push(lr_allele_indices.clone());
                    new_read_sets.push(allele_indices.clone());
                }

                prev_nodes = new_nodes;
                prev_read_sets = new_read_sets;
            }

            index1 = index2;
        } else {
            index1 += 1;
        }
    }

    let sink_node = graph.add_node(String::from("$"));
    for prev_node in prev_nodes {
        graph.add_edge(prev_node, sink_node, lr_msas.len() + sr_msas.len());
    }

    (graph, source_node, sink_node)
}

fn create_graph(lr_msas: &Vec<String>, sr_msas: &Vec<String>) -> Graph<String, (HashSet<usize>, f32)> {
    let mut graph = Graph::new();

    let mut prev_nodes = Vec::new();

    let mut phase_matrix = Vec::new();

    let mut contig = String::new();
    let mut index1 = 0;
    while index1 < lr_msas[0].len() {
        let combined_base_counts = pileup_counts(lr_msas, sr_msas, index1);
        let bases = pileup_indices(lr_msas, sr_msas, index1);

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
                graph.add_edge(*prev_node, node, (HashSet::new(), 0.0));
            });
            prev_nodes.clear();

            let mut index2 = index1;
            let mut allele_base_counts = pileup_counts(lr_msas, sr_msas, index2);

            while index2 < lr_msas[0].len() && allele_base_counts.len() > 1 {
                index2 += 1;
                allele_base_counts = pileup_counts(lr_msas, sr_msas, index2);
            }

            let allele_counts = allele_counts(lr_msas, sr_msas, index1, index2)
                .into_iter()
                .filter(|&(_, count)| count >= 5)
                .sorted_by(|a, b| b.1.cmp(&a.1))
                .take(2)
                .collect::<BTreeMap<_, _>>();

            let filtered_allele_counts = allele_counts.into_iter().filter(|(_, count)| *count > 5).collect::<BTreeMap<_, _>>();

            if filtered_allele_counts.len() == 1 {
                for (allele, _) in filtered_allele_counts {
                    let new_node = graph.add_node(allele.clone());

                    graph.add_edge(node, new_node, (HashSet::new(), 0.0));
                    prev_nodes.push(new_node);
                }
            } else {
                let alleles = get_allele_indices(lr_msas, sr_msas, index1, index2);

                let mut dual_edges = BTreeMap::new();

                let mut phase = false;
                for (allele, _) in filtered_allele_counts {
                    let new_node = graph.add_node(allele.clone());

                    let allele_indices = alleles.iter().enumerate().filter_map(|(i, a)| {
                        if *a == allele {
                            Some(i)
                        } else {
                            None
                        }
                    })
                    .collect::<HashSet<_>>();

                    // let phase = rng.gen_bool(0.5);
                    let edge = graph.add_edge(node, new_node, (allele_indices, if phase { 1.0 } else { -1.0 }));
                    phase = !phase;

                    dual_edges.insert(edge, phase);

                    prev_nodes.push(new_node);
                }

                if dual_edges.len() == 2 {
                    phase_matrix.push(dual_edges);
                }
            }

            index1 = index2;

            contig = String::new();
        }
    }

    for i in 1..phase_matrix.len() {
        let entry = &phase_matrix[i];

        let (prev_edge0, prev_phase0) = phase_matrix[i - 1].first_key_value().unwrap();
        let (prev_edge1, prev_phase1) = phase_matrix[i - 1].last_key_value().unwrap();
        let prev_weight0 = graph.edge_weight(*prev_edge0).unwrap().0.clone();
        let prev_weight1 = graph.edge_weight(*prev_edge1).unwrap().0.clone();

        let (this_edge0, this_phase0) = entry.first_key_value().unwrap();
        let (this_edge1, this_phase1) = entry.last_key_value().unwrap();
        let this_weight0 = graph.edge_weight(*this_edge0).unwrap().0.clone();

        let intersection0 = prev_weight0.intersection(&this_weight0).count();
        let intersection1 = prev_weight1.intersection(&this_weight0).count();

        let (new_phase0, new_phase1) = if intersection0 > intersection1 {
            (*prev_phase0, *prev_phase1)
        } else {
            (*prev_phase1, *prev_phase0)
        };

        let mut graph = graph.clone();
        graph.edge_weight_mut(*this_edge0).unwrap().1 = if new_phase0 { 1.0 } else { -1.0 };
        graph.edge_weight_mut(*this_edge1).unwrap().1 = if new_phase1 { 1.0 } else { -1.0 };

        let mut new_entry = BTreeMap::new();
        new_entry.insert(*this_edge0, new_phase0);
        new_entry.insert(*this_edge1, new_phase1);
        phase_matrix[i] = new_entry;
    }

    graph
}

fn create_fully_phased_haplotypes(lr_msas: &Vec<String>, sr_msas: &Vec<String>, h1: &Vec<u8>) -> (String, String) {
    let mut index1 = 0;
    let mut hap_index = 0;

    let mut hap1 = String::new();
    let mut hap2 = String::new();

    while index1 < lr_msas[0].len() {
        let combined_base_counts = pileup_counts(lr_msas, sr_msas, index1);
        let bases = pileup_indices(lr_msas, sr_msas, index1);

        if combined_base_counts.len() == 0 {
            index1 += 1;
        } else if combined_base_counts.len() == 1 {
            let base = combined_base_counts.keys().next().unwrap();
            hap1.push(*base);
            hap2.push(*base);

            index1 += 1;
        } else {
            let mut index2 = index1;
            let mut allele_base_counts = pileup_counts(lr_msas, sr_msas, index2);

            while index2 < lr_msas[0].len() && allele_base_counts.len() > 1 {
                index2 += 1;
                allele_base_counts = pileup_counts(lr_msas, sr_msas, index2);
            }

            let allele_counts = allele_counts(lr_msas, sr_msas, index1, index2)
                .into_iter()
                .filter(|&(_, count)| count >= 5)
                .sorted_by(|a, b| b.1.cmp(&a.1))
                .take(2)
                .collect::<BTreeMap<_, _>>();

            let filtered_allele_counts = allele_counts.into_iter().filter(|(_, count)| *count > 5).collect::<BTreeMap<_, _>>();

            if filtered_allele_counts.len() == 1 {
                for (allele, _) in filtered_allele_counts {
                    hap1.extend(allele.chars());
                    hap2.extend(allele.chars());
                }
            } else {
                // let alleles = get_allele_indices(lr_msas, sr_msas, index1, index2);

                if hap_index < h1.len() {
                    let mut phase = h1[hap_index] == 1;
                    for (allele, _) in filtered_allele_counts {
                        if !phase {
                            hap1.extend(allele.chars());
                        } else {
                            hap2.extend(allele.chars());
                        }

                        phase = !phase;
                    }

                    hap_index += 1;
                }
            }

            index1 = index2;
        }
    }

    (hap1, hap2)
}

fn pileup_indices(lr_msas: &Vec<String>, sr_msas: &Vec<String>, index: usize) -> Vec<char> {
    let bases = lr_msas.iter().chain(sr_msas.iter())
        .map(|msa| msa.chars().nth(index).unwrap_or(' '))
        .collect::<Vec<char>>();

    bases
}

fn pileup_counts(lr_msas: &Vec<String>, sr_msas: &Vec<String>, index: usize) -> BTreeMap<char, i32> {
    let combined_base_counts = lr_msas.iter().chain(sr_msas.iter())
        .map(|msa| msa.chars().nth(index).unwrap_or(' '))
        .filter(|&c| c != ' ')
        .fold(BTreeMap::new(), |mut counts, base| {
            *counts.entry(base).or_insert(0) += 1;
            counts
        });
    combined_base_counts
}

fn get_allele_indices(lr_msas: &Vec<String>, sr_msas: &Vec<String>, index1: usize, index2: usize) -> Vec<String> {
    let alleles = lr_msas.iter().chain(sr_msas.iter())
        .map(|msa| msa[index1..index2].to_string().replace(" ", ""))
        .collect::<Vec<String>>();
    alleles
}

fn allele_counts(lr_msas: &Vec<String>, sr_msas: &Vec<String>, index1: usize, index2: usize) -> BTreeMap<String, i32> {
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

// fn allele_counts(msas: &Vec<String>, column: usize) -> BTreeMap<char, usize> {
//     let mut counts = BTreeMap::new();

//     for msa in msas {
//         let char = msa.chars().nth(column).unwrap();
//         if char != ' ' {
//             *counts.entry(char).or_insert(0) += 1;
//         }
//     }

//     counts
// }

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
