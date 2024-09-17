use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet, VecDeque};
use std::path::PathBuf;
use std::{
    fs::File,
    io::{BufWriter, Write},
};

use bio::alignment::pairwise::Scoring;
use bio::alignment::poa::Aligner;
use grb::expr::LinExpr;
use indicatif::ParallelProgressIterator;
use parquet::data_type::AsBytes;
use petgraph::algo::ford_fulkerson;
use petgraph::dot::Dot;
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::visit::{EdgeRef, IntoEdgesDirected};
use petgraph::Graph;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;

use itertools::Itertools;

use spoa::AlignmentType;

use grb::prelude::*;

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

            let g = create_graph(&group_lr_msas, &group_sr_msas);

            let m = create_read_allele_matrix(&group_lr_msas, &group_sr_msas);

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

            let (haplotype1, haplotype2) = wmec_ilp(&m, &group_lr_msas, &group_sr_msas).unwrap();

            println!("Haplotype 1: {:?}", haplotype1);
            println!("Haplotype 2: {:?}", haplotype2);

            skydive::elog!("Creating het graph");
            let (h, source_node, sink_node) = create_het_graph(&group_lr_msas, &group_sr_msas);

            // h.edges_directed(source_node, petgraph::Direction::Outgoing).for_each(|e| {
            //     let flow = ford_fulkerson(&h, e.target(), sink_node);

            //     skydive::elog!("Flow: {:?}", flow);
            // });

            if lr_group_index > 0 {
                skydive::elog!("Writing het graph");

                // Write graph h to a dot file
                let dot_file = File::create(format!("graph_h_cluster_{}.dot", lr_group_index)).unwrap();
                let mut dot_writer = BufWriter::new(dot_file);
                write!(dot_writer, "{:?}", Dot::new(&h)).unwrap();
            }

            // let mut intersection_matrix: BTreeMap<EdgeIndex, BTreeMap<EdgeIndex, usize>> = BTreeMap::new();

            // for edge_index1 in g.edge_indices() {
            //     let weight1 = g.edge_weight(edge_index1).unwrap();

            //     for edge_index2 in g.edge_indices() {
            //         let weight2 = g.edge_weight(edge_index2).unwrap();

            //         if (weight1.len() > 0) && (weight2.len() > 0) {
            //             let intersection = weight1.intersection(weight2).count();

            //             intersection_matrix
            //                 .entry(edge_index1)
            //                 .or_insert_with(BTreeMap::new)
            //                 .insert(edge_index2, intersection);
            //         }
            //     }
            // }

            // let p = phase_graph(g, &intersection_matrix);

            // println!("Intersection Matrix:");
            // println!("{:<10} {:<10} {:<10}", "Edge 1", "Edge 2", "Intersection");
            // println!("{:-<32}", "");
            // for (edge1, intersections) in &intersection_matrix {
            //     for (edge2, count) in intersections {
            //         if *count > 0 {
            //             println!("{:<10} {:<10} {:<10}", edge1.index(), edge2.index(), count);
            //         }
            //     }
            // }
            // println!();

            // let (hap1, hap2) = haplotypes(&g);
            let hap1 = haplotype(&g, true);
            let hap2 = haplotype(&g, false);

            println!("Consensus 1: {}", hap1);
            println!("Consensus 2: {}", hap2);

            /*

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

fn max_flow(g: &Graph<String, usize>, source_node: NodeIndex, sink_node: NodeIndex) -> Vec<NodeIndex> {
    // use std::collections::{VecDeque, HashMap};
    // use petgraph::graph::NodeIndex;
    // use petgraph::visit::EdgeRef;

    // Initialize residual graph
    let mut residual_graph = g.clone();
    let mut parent_map = HashMap::new();
    let mut max_flow = 0;

    loop {
        // BFS to find augmenting path
        parent_map.clear();
        let mut queue = VecDeque::new();
        queue.push_back(source_node);
        parent_map.insert(source_node, None);

        while let Some(current) = queue.pop_front() {
            if current == sink_node {
                break;
            }

            for edge in residual_graph.edges(current) {
                let next = edge.target();
                let capacity = *edge.weight();
                if capacity > 0 && !parent_map.contains_key(&next) {
                    parent_map.insert(next, Some((current, edge.id())));
                    queue.push_back(next);
                }
            }
        }

        // If sink is not reached, we're done
        if !parent_map.contains_key(&sink_node) {
            break;
        }

        // Find minimum residual capacity along the path
        let mut path_flow = usize::MAX;
        let mut current = sink_node;
        while let Some((parent, edge_id)) = parent_map[&current] {
            path_flow = path_flow.min(residual_graph[edge_id]);
            current = parent;
        }

        // Update residual capacities
        let mut current = sink_node;
        while let Some((parent, edge_id)) = parent_map[&current] {
            residual_graph[edge_id] -= path_flow;
            if let Some(reverse_edge) = residual_graph.find_edge(current, parent) {
                residual_graph[reverse_edge] += path_flow;
            } else {
                residual_graph.add_edge(current, parent, path_flow);
            }
            current = parent;
        }

        max_flow += path_flow;
    }

    // Reconstruct the path
    let mut path = Vec::new();
    let mut current = sink_node;
    while current != source_node {
        path.push(current);
        for edge in g.edges_directed(current, petgraph::Direction::Incoming) {
            if residual_graph[edge.id()] < g[edge.id()] {
                current = edge.source();
                break;
            }
        }
    }

    path.push(source_node);
    path.reverse();

    path
}

fn haplotype(g: &Graph<String, (HashSet<usize>, f32)>, hap1: bool) -> String {
    // Find the start node (a node with no incoming edges)
    let start_node = g.node_indices().find(|&n| g.edges_directed(n, petgraph::Direction::Incoming).count() == 0)
        .expect("No start node found");

    // Initialize the path with the start node
    let mut path = vec![start_node];
    let mut current_node = start_node;

    // Traverse the graph
    while let Some(next_edge) = g.edges_directed(current_node, petgraph::Direction::Outgoing)
        .filter(|e| if hap1 { g.edge_weight(e.id()).unwrap().1 > -0.1 } else { g.edge_weight(e.id()).unwrap().1 < 0.1 })
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

fn wmec_ilp(m: &Vec<BTreeMap<usize, String>>, lr_msas: &Vec<String>, sr_msas: &Vec<String>) -> Result<(Vec<usize>, Vec<usize>), grb::Error> {
    let n_variants = m.len();
    let n_reads = lr_msas.len();

    // Create a new Gurobi model
    let mut model = Model::new("wmec")?;

    // Add binary variables for each read
    let vars: Vec<Var> = (0..n_reads)
        // .map(|i| model.add_var(&format!("x_{}", i), Binary, 0.0, 0.0, 1.0, &[]))
        .map(|i| add_binvar!(model, name: &format!("x_{}", i), bounds: ..).unwrap())
        .collect();

    // Add constraints
    for (i, variant) in m.iter().enumerate() {
        let mut expr = LinExpr::new();
        for (&read, allele) in variant {
            let coeff = if allele == "0" { 1.0 } else { -1.0 };
            expr.add_term(coeff, vars[read]); // Use add_term to build the expression
        }
        // model.add_constr(&format!("c1_{}", i), expr.clone() <= 0.0)?;
        // model.add_constr(&format!("c2_{}", i), expr >= -1.0)?;

        model.add_constr(&format!("c1_{}", i), c!(expr.clone() <=  0.0)).unwrap();
        model.add_constr(&format!("c2_{}", i), c!(expr.clone() >= -1.0)).unwrap();
    }

    // Set objective function (minimize conflicts)
    let mut obj = LinExpr::new();
    for variant in m {
        for (&read1, allele1) in variant {
            for (&read2, allele2) in variant {
                if read1 < read2 && allele1 != allele2 {
                    obj.add_term(1.0, vars[read1]);
                    obj.add_term(1.0, vars[read2]);
                }
            }
        }
    }
    model.set_objective(obj, Minimize).unwrap();

    // Optimize the model
    model.optimize().unwrap();

    // Get solution
    let solution: Vec<f64> = vars
        .iter()
        .map(|var| model.get_obj_attr(attr::X, var).unwrap())
        .collect();

    // Partition reads based on solution
    let mut haplotype1 = Vec::new();
    let mut haplotype2 = Vec::new();
    for i in 0..n_reads {
        if solution[i] > 0.5 {
            haplotype1.push(i);
        } else {
            haplotype2.push(i);
        }
    }

    Ok((haplotype1, haplotype2))
}

fn wmec(m: Vec<BTreeMap<usize, String>>, lr_msas: &Vec<String>, sr_msas: &Vec<String>) -> (Vec<usize>, Vec<usize>) {
    let n_variants = m.len();
    let n_reads = lr_msas.len();

    skydive::elog!("n_variants: {}", n_variants);
    skydive::elog!("n_reads: {}", n_reads);

    // Initialize dynamic programming table
    let mut dp = vec![vec![f64::MAX; 1 << n_reads]; n_variants + 1];
    dp[0][0] = 0.0;

    skydive::elog!("allocated dp");

    // Initialize backtracking table
    let mut backtrack = vec![vec![0; 1 << n_reads]; n_variants + 1];

    skydive::elog!("allocated backtrack");

    // Fill the dynamic programming table
    for i in 1..=n_variants {
        for mask in 0..(1 << n_reads) {
            for read in 0..n_reads {
                if mask & (1 << read) != 0 {
                    continue;
                }
                
                let new_mask = mask | (1 << read);
                let mut cost = 0.0;

                if let Some(allele) = m[i-1].get(&read) {
                    for (other_read, other_allele) in &m[i-1] {
                        if mask & (1 << other_read) != 0 && allele != other_allele {
                            cost += 1.0;
                        }
                    }
                }

                let new_cost = dp[i-1][mask] + cost;
                if new_cost < dp[i][new_mask] {
                    dp[i][new_mask] = new_cost;
                    backtrack[i][new_mask] = mask;
                }
            }
        }
    }

    // Find the optimal partition
    let mut optimal_mask = 0;
    let mut min_cost = f64::MAX;
    for mask in 0..(1 << n_reads) {
        if dp[n_variants][mask] < min_cost {
            min_cost = dp[n_variants][mask];
            optimal_mask = mask;
        }
    }

    // Backtrack to reconstruct the optimal partition
    let mut haplotype1 = Vec::new();
    let mut haplotype2 = Vec::new();
    let mut current_mask = optimal_mask;

    for i in (1..=n_variants).rev() {
        let prev_mask = backtrack[i][current_mask];
        let changed_read = (current_mask ^ prev_mask).trailing_zeros() as usize;
        
        if current_mask & (1 << changed_read) != 0 {
            haplotype1.push(changed_read);
        } else {
            haplotype2.push(changed_read);
        }

        current_mask = prev_mask;
    }

    haplotype1.reverse();
    haplotype2.reverse();

    (haplotype1, haplotype2)
}

#[derive(Debug)]
struct WMECData {
    reads: Vec<Vec<Option<u8>>>, // Reads matrix where None represents missing data
    confidences: Vec<Vec<Option<u32>>>, // Confidence degrees matrix
    num_snps: usize,             // Number of SNP positions
}

impl WMECData {
    // Initialize the data structure with given reads and confidences
    fn new(reads: Vec<Vec<Option<u8>>>, confidences: Vec<Vec<Option<u32>>>) -> Self {
        let num_snps = reads[0].len();
        WMECData { reads, confidences, num_snps }
    }
    
    // Function to compute W^0(j, R) and W^1(j, R)
    // Cost to set all fragments in set R to 0 or 1 at SNP j
    fn compute_costs(&self, snp: usize, set_r: &BTreeSet<usize>) -> (u32, u32) {
        let mut w0 = 0; // Cost for setting to 0
        let mut w1 = 0; // Cost for setting to 1

        for &read_index in set_r {
            if let Some(allele) = self.reads[read_index][snp] {
                if let Some(confidence) = self.confidences[read_index][snp] {
                    if allele == 0 {
                        w1 += confidence; // Cost of flipping 0 -> 1
                    } else {
                        w0 += confidence; // Cost of flipping 1 -> 0
                    }
                }
            }
        }

        (w0, w1)
    }
    
    // Calculate minimum correction cost Delta C(j, (R, S))
    fn delta_c(&self, snp: usize, r: &BTreeSet<usize>, s: &BTreeSet<usize>) -> u32 {
        let (w0_r, w1_r) = self.compute_costs(snp, r);
        let (w0_s, w1_s) = self.compute_costs(snp, s);

        std::cmp::min(w0_r, w1_r) + std::cmp::min(w0_s, w1_s)
    }
}

// Function to generate all bipartitions of a set
fn generate_bipartitions(set: &BTreeSet<usize>) -> Vec<(BTreeSet<usize>, BTreeSet<usize>)> {
    let mut partitions = vec![];
    let set_vec: Vec<_> = set.iter().collect();

    let num_partitions = 1 << set_vec.len(); // 2^|set|
    for i in 0..num_partitions {
        let mut r = BTreeSet::new();
        let mut s = BTreeSet::new();

        for (j, &elem) in set_vec.iter().enumerate() {
            if i & (1 << j) == 0 {
                r.insert(*elem);
            } else {
                s.insert(*elem);
            }
        }

        partitions.push((r, s));
    }

    partitions
}

// Function to initialize the DP table for SNP 0
fn initialize_dp(data: &WMECData) -> (HashMap<(usize, BTreeSet<usize>, BTreeSet<usize>), u32>, HashMap<(usize, BTreeSet<usize>, BTreeSet<usize>), Option<(BTreeSet<usize>, BTreeSet<usize>)>>) {
    let mut dp: HashMap<(usize, BTreeSet<usize>, BTreeSet<usize>), u32> = HashMap::new();
    let mut backtrack: HashMap<(usize, BTreeSet<usize>, BTreeSet<usize>), Option<(BTreeSet<usize>, BTreeSet<usize>)>> = HashMap::new();

    let active_fragments: BTreeSet<usize> = data.reads.iter().enumerate()
        .filter(|(_, read)| read[0].is_some()) // Only consider fragments covering SNP 0
        .map(|(index, _)| index)
        .collect();

    let partitions = generate_bipartitions(&active_fragments);
    for (r, s) in partitions {
        let cost = data.delta_c(0, &r, &s);
        dp.insert((0, r.clone(), s.clone()), cost);
        backtrack.insert((0, r.clone(), s.clone()), None);
    }

    (dp, backtrack)
}

// Function to update the DP table for each SNP position
fn update_dp(data: &WMECData, dp: &mut HashMap<(usize, BTreeSet<usize>, BTreeSet<usize>), u32>, backtrack: &mut HashMap<(usize, BTreeSet<usize>, BTreeSet<usize>), Option<(BTreeSet<usize>, BTreeSet<usize>)>>, snp: usize) {
    let active_fragments: BTreeSet<usize> = data.reads.iter().enumerate()
        .filter(|(_, read)| read[snp].is_some()) // Only consider fragments covering SNP
        .map(|(index, _)| index)
        .collect();
    let partitions = generate_bipartitions(&active_fragments);

    for (r, s) in &partitions {
        let delta_cost = data.delta_c(snp, r, s);
        let mut min_cost = u32::MAX;
        let mut best_bipartition = None;

        let prev_active_fragments: BTreeSet<usize> = data.reads.iter().enumerate()
            .filter(|(_, read)| read[snp - 1].is_some()) // Fragments covering the previous SNP
            .map(|(index, _)| index)
            .collect();
        
        for (prev_r, prev_s) in generate_bipartitions(&prev_active_fragments) {
            // let r_intersect = prev_r.intersection(&active_fragments).cloned().collect::<BTreeSet<_>>();
            // let s_intersect = prev_s.intersection(&active_fragments).cloned().collect::<BTreeSet<_>>();

            // if r_intersect.is_subset(&prev_r) && s_intersect.is_subset(&prev_s) {
            let r_compatible = r.intersection(&prev_active_fragments).all(|&x| prev_r.contains(&x));
            let s_compatible = s.intersection(&prev_active_fragments).all(|&x| prev_s.contains(&x));

            if r_compatible && s_compatible {
                if let Some(&prev_cost) = dp.get(&(snp - 1, prev_r.clone(), prev_s.clone())) {
                    let current_cost = delta_cost + prev_cost;

                    // println!("{:?} {:?} {} {}", r_intersect, s_intersect, delta_cost, prev_cost);

                    if current_cost < min_cost {
                        min_cost = current_cost;
                        best_bipartition = Some((prev_r.clone(), prev_s.clone()));
                    }
                }
            }
        }

        dp.insert((snp, r.clone(), s.clone()), min_cost);
        backtrack.insert((snp, r.clone(), s.clone()), best_bipartition);
    }
}

/*
// Function to backtrack and find the optimal haplotypes
fn backtrack_haplotypes(data: &WMECData, dp: &HashMap<(usize, BTreeSet<usize>, BTreeSet<usize>), u32>, backtrack: &HashMap<(usize, BTreeSet<usize>, BTreeSet<usize>), Option<(BTreeSet<usize>, BTreeSet<usize>)>>) -> (Vec<u8>, Vec<u8>) {
    let mut best_cost = u32::MAX;
    let mut best_bipartition = None;
    let final_active_fragments: BTreeSet<usize> = data.reads.iter().enumerate()
        .filter(|(_, read)| read[data.num_snps - 1].is_some()) // Only consider fragments covering final SNP
        .map(|(index, _)| index)
        .collect();

    for (r, s) in generate_bipartitions(&final_active_fragments) {
        if let Some(&cost) = dp.get(&(data.num_snps - 1, r.clone(), s.clone())) {
            if cost < best_cost {
                best_cost = cost;
                best_bipartition = Some((r, s));
            }
        }
    }

    let mut haplotype1 = vec![0; data.num_snps];
    let mut haplotype2 = vec![1; data.num_snps];
    let mut current_bipartition = best_bipartition.unwrap();

    for snp in (0..data.num_snps).rev() {
        let (r, s) = current_bipartition;
        let (w0_r, w1_r) = data.compute_costs(snp, &r);
        if w0_r < w1_r {
            haplotype1[snp] = 0;
        } else {
            haplotype1[snp] = 1;
        }

        let (w0_s, w1_s) = data.compute_costs(snp, &s);
        if w0_s < w1_s {
            haplotype2[snp] = 0;
        } else {
            haplotype2[snp] = 1;
        }

        if let Some(next_bipartition) = backtrack.get(&(snp, r, s)).and_then(|x| x.as_ref()) {
            current_bipartition = next_bipartition.clone();
        } else {
            break;
        }
    }

    (haplotype1, haplotype2)
}
*/

fn backtrack_haplotypes(data: &WMECData, dp: &HashMap<(usize, BTreeSet<usize>, BTreeSet<usize>), u32>, backtrack: &HashMap<(usize, BTreeSet<usize>, BTreeSet<usize>), Option<(BTreeSet<usize>, BTreeSet<usize>)>>) -> (Vec<u8>, Vec<u8>) {
    let mut best_cost = u32::MAX;
    let mut best_bipartition = None;
    let final_active_fragments: BTreeSet<usize> = data.reads.iter().enumerate()
        .filter(|(_, read)| read[data.num_snps - 1].is_some())
        .map(|(index, _)| index)
        .collect();

    for (r, s) in generate_bipartitions(&final_active_fragments) {
        if let Some(&cost) = dp.get(&(data.num_snps - 1, r.clone(), s.clone())) {
            if cost < best_cost {
                best_cost = cost;
                best_bipartition = Some((r, s));
            }
        }
    }

    let mut haplotype1 = vec![0; data.num_snps];
    let mut haplotype2 = vec![0; data.num_snps];
    let mut current_bipartition = best_bipartition.unwrap();

    for snp in (0..data.num_snps).rev() {
        let (r, s) = &current_bipartition;
        let (w0_r, w1_r) = data.compute_costs(snp, r);
        let (w0_s, w1_s) = data.compute_costs(snp, s);

        if w0_r + w1_s <= w1_r + w0_s {
            haplotype1[snp] = 0;
            haplotype2[snp] = 1;
        } else {
            haplotype1[snp] = 1;
            haplotype2[snp] = 0;
        }

        if snp > 0 {
            if let Some(prev_bipartition) = backtrack.get(&(snp, r.clone(), s.clone())).and_then(|x| x.as_ref()) {
                current_bipartition = prev_bipartition.clone();
            } else {
                // This should not happen if the DP table is correctly filled
                panic!("No valid previous bipartition found for SNP {}", snp);
            }
        }
    }

    (haplotype1, haplotype2)
}

// Main function to perform WMEC using dynamic programming
fn perform_wmec(data: &WMECData) -> (Vec<u8>, Vec<u8>) {
    let (mut dp, mut backtrack) = initialize_dp(data);

    for snp in 1..data.num_snps {
        update_dp(data, &mut dp, &mut backtrack, snp);
    }

    backtrack_haplotypes(data, &dp, &backtrack)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_initialization_example() {
        // Define the dataset as a matrix of reads
        let reads = vec![
            vec![Some(0), None],    // f0
            vec![Some(1), Some(0)], // f1
            vec![Some(1), Some(1)], // f2
            vec![None,    Some(0)], // f3
        ];

        // Define the confidence degrees
        let confidences = vec![
            vec![Some(5), None],    // f0
            vec![Some(3), Some(2)], // f1
            vec![Some(6), Some(1)], // f2
            vec![None,    Some(2)], // f3
        ];

        // Initialize the dataset
        let data = WMECData::new(reads, confidences);

        // Define the expected values for C(1, ·)
        let expected_values = vec![
            (vec![0, 1, 2], vec![], 5), // C(1, ({f0, f1, f2}, ∅)) = 5
            (vec![0, 1], vec![2], 3),   // C(1, ({f0, f1}, {f2})) = 3
            (vec![0, 2], vec![1], 5),   // C(1, ({f0, f2}, {f1})) = 5
            (vec![1, 2], vec![0], 0),   // C(1, ({f1, f2}, {f0})) = 0
        ];

        // Check that the values in the dynamic programming table match the expected values
        for (r, s, expected_cost) in expected_values {
            let r_set: BTreeSet<usize> = r.into_iter().collect();
            let s_set: BTreeSet<usize> = s.into_iter().collect();
            let cost = data.delta_c(0, &r_set, &s_set);
            assert_eq!(cost, expected_cost, "Cost for partition ({:?}, {:?}) does not match expected", r_set, s_set);
        }
    }

    #[test]
    fn test_example_2() {
        // Define the dataset as a matrix of reads
        let reads = vec![
            vec![Some(0), None],    // f0
            vec![Some(1), Some(0)], // f1
            vec![Some(1), Some(1)], // f2
            vec![None,    Some(0)], // f3
        ];

        // Define the confidence degrees
        let confidences = vec![
            vec![Some(5), None],    // f0
            vec![Some(3), Some(2)], // f1
            vec![Some(6), Some(1)], // f2
            vec![None,    Some(2)], // f3
        ];

        // Initialize the dataset
        let data = WMECData::new(reads, confidences);

        // Define the expected values for C(2, ·)
        let expected_values = vec![
            (vec![1, 2, 3], vec![], 1), // C(2, ({f1, f2, f3}, ∅)) = 1
            (vec![1, 2], vec![3], 1),   // C(2, ({f1, f2}, {f3})) = 1
            (vec![1, 3], vec![2], 3),   // C(2, ({f1, f3}, {f2})) = 3
            (vec![2, 3], vec![1], 4),   // C(2, ({f2, f3}, {f1})) = 4
        ];

        let (mut dp, mut backtrack) = initialize_dp(&data);

        for snp in 1..data.num_snps {
            update_dp(&data, &mut dp, &mut backtrack, snp);
        }

        // Verify that the results comport with expected_values
        for (r, s, expected_cost) in expected_values {
            let r_set: BTreeSet<usize> = r.into_iter().collect();
            let s_set: BTreeSet<usize> = s.into_iter().collect();
            let actual_cost = dp.get(&(1, r_set.clone(), s_set.clone())).unwrap_or(&u32::MAX);
            assert_eq!(*actual_cost, expected_cost, "Cost for partition ({:?}, {:?}) does not match expected", r_set, s_set);
        }
    }

    #[test]
    fn test_whatshap_example() {
        // Define the dataset as a matrix of reads
        let reads = vec![
            vec![Some(0), None,    None,    Some(0), Some(1)], // f0
            vec![Some(1), Some(0), Some(0), None,    None   ], // f1
            vec![Some(1), Some(1), Some(0), None,    None   ], // f2
            vec![None,    Some(0), Some(0), Some(1), None   ], // f3
            vec![None,    None,    Some(1), Some(0), Some(1)], // f4
            vec![None,    None,    None,    Some(0), Some(1)], // f5
        ];

        // Define the confidence degrees
        let confidences = vec![
            vec![Some(32), None,     None,     Some(34), Some(17)], // f0
            vec![Some(15), Some(25), Some(13), None,     None    ], // f1
            vec![Some(7),  Some(3),  Some(15), None,     None    ], // f2
            vec![None,     Some(12), Some(23), Some(29), Some(31)], // f3
            vec![None,     None,     Some(25), Some(17), Some(19)], // f4
            vec![None,     None,     None,     Some(20), Some(10)], // f5
        ];

        // Initialize the dataset
        let data = WMECData::new(reads, confidences);

        // Perform the WMEC algorithm using dynamic programming
        let (haplotype1, haplotype2) = perform_wmec(&data);

        let expected_haplotype1 = vec![0, 1, 1, 0, 1];
        let expected_haplotype2 = vec![1, 0, 0, 1, 0];

        assert_eq!(haplotype1, expected_haplotype1, "Haplotype 1 does not match expected");
        assert_eq!(haplotype2, expected_haplotype2, "Haplotype 2 does not match expected");
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
                skydive::elog!("filt {:?} {}", filtered_allele_counts, prev_nodes.len());
                
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

                    skydive::elog!("Adding node {} {:?} {} {:?}", index1, filtered_allele_counts, allele, allele_indices);

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

    for (i, row) in phase_matrix.iter().enumerate() {
        println!("{} {:?}", i, row);

        if i > 5 { break; }
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

    // for (i, row) in phase_matrix.iter().enumerate() {
    //     println!("{} {:?}", i, row);

    //     if i > 5 { break; }
    // }

    graph
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
