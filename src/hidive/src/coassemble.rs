use std::collections::{BTreeMap, HashMap, HashSet};
use std::{fs::File, path::PathBuf, io::Write};

use itertools::Itertools;
use linked_hash_map::LinkedHashMap;
use minimap2::Aligner;
use needletail::Sequence;
use petgraph::graph::NodeIndex;
use rayon::prelude::*;
use indicatif::ParallelProgressIterator;

use skydive::ldbg::LdBG;
use skydive::mldbg::MLdBG;

use skydive::wmec::*;
use spoa::AlignmentType;

pub fn start(
    output: &PathBuf,
    kmer_size: usize,
    model_path: &PathBuf,
    reference_fasta_paths: &Vec<PathBuf>,
    long_read_fasta_path: PathBuf,
    short_read_fasta_path: PathBuf,
) {
    let long_read_seq_urls = skydive::parse::parse_file_names(&[long_read_fasta_path.clone()]);
    let short_read_seq_urls = skydive::parse::parse_file_names(&[short_read_fasta_path.clone()]);
    let reference_seq_urls = skydive::parse::parse_file_names(reference_fasta_paths);

    // Read all long reads.
    skydive::elog!("Processing long-read samples {:?}...", long_read_seq_urls.iter().map(|url| url.as_str()).collect::<Vec<&str>>());
    let all_lr_seqs = skydive::utils::read_fasta(&vec![long_read_fasta_path.clone()]);

    // Read all short reads.
    skydive::elog!("Processing short-read samples {:?}...", short_read_seq_urls.iter().map(|url| url.as_str()).collect::<Vec<&str>>());
    let all_sr_seqs = skydive::utils::read_fasta(&vec![short_read_fasta_path]);

    // Read all reference subsequences.
    skydive::elog!("Processing reference {:?}...", reference_seq_urls.iter().map(|url| url.as_str()).collect::<Vec<&str>>());
    let all_ref_seqs = skydive::utils::read_fasta(reference_fasta_paths);

    let l1 = LdBG::from_sequences("lr".to_string(), kmer_size, &all_lr_seqs);
    let s1 = LdBG::from_sequences("sr".to_string(), kmer_size, &all_sr_seqs);

    let m = MLdBG::from_ldbgs(vec![l1, s1])
        .score_kmers(model_path)
        .collapse()
        .clean(0.1)
        .build_links(&all_lr_seqs, false);

    skydive::elog!("Built MLdBG with {} k-mers.", m.kmers.len());

    skydive::elog!("Correcting reads...");
    let corrected_lr_seqs = m.correct_seqs(&all_lr_seqs);

    skydive::elog!("Clustering reads...");
    let (reads_hap1, reads_hap2) = cluster_reads(&m, &corrected_lr_seqs);

    skydive::elog!("Assembling haplotype 1...");
    let asm1 = assemble_haplotype(&all_ref_seqs, &reads_hap1);
    // let hap1 = refine_haplotype(asm1, reference_fasta_paths[0].clone(), all_ref_seqs.get(0).unwrap());

    skydive::elog!("Assembling haplotype 2...");
    let asm2 = assemble_haplotype(&all_ref_seqs, &reads_hap2);
    // let hap2 = refine_haplotype(asm2, reference_fasta_paths[0].clone(), all_ref_seqs.get(0).unwrap());

    let mut output = File::create(output).expect("Failed to create output file");
    writeln!(output, ">hap1").expect("Failed to write to output file");
    writeln!(output, "{}", asm1).expect("Failed to write to output file");
    writeln!(output, ">hap2").expect("Failed to write to output file");
    writeln!(output, "{}", asm2).expect("Failed to write to output file");

}

fn refine_haplotype(asm: String, ref_path: PathBuf, ref_seq: &Vec<u8>) -> String {
    let aligner = Aligner::builder()
        .map_hifi()
        .with_cigar()
        .with_index(ref_path.clone(), None)
        .expect("Failed to build aligner");

    let results = vec![asm]
        .iter()
        .map(|hap| aligner.map(hap.as_bytes(), true, false, None, None, None).unwrap())
        .collect::<Vec<_>>();

    let mut alt_seq: Vec<u8> = Vec::new();
    for result in results {
        for mapping in result {
            if let Some(alignment) = &mapping.alignment {
                if mapping.is_primary && !mapping.is_supplementary {
                    if let Some(cs) = &alignment.cs {
                        let re = regex::Regex::new(r"([:\*\+\-])(\w+)").unwrap();
                        let mut cursor = 0;
                        for cap in re.captures_iter(cs) {
                            let op = cap.get(1).unwrap().as_str();
                            let seq = cap.get(2).unwrap().as_str();

                            match op {
                                ":" => {
                                    // Match (copy from reference)
                                    let length = seq.parse::<usize>().unwrap();
                                    alt_seq.extend(&ref_seq[cursor..cursor + length]);
                                    cursor += length;
                                }
                                "*" => {
                                    // Substitution
                                    // let ref_base = seq.as_bytes()[0];
                                    let alt_base = seq.to_uppercase().as_bytes()[1];
                                    alt_seq.push(alt_base);
                                    cursor += 1;
                                }
                                "+" => {
                                    // Insertion
                                    alt_seq.extend(seq.to_uppercase().as_bytes());
                                }
                                "-" => {
                                    // Deletion
                                    // alt_seq.extend(vec![b'-'; seq.len()]);
                                    cursor += seq.len();
                                }
                                _ => unreachable!("Invalid CIGAR operation"),
                            }
                        }
                    }
                }

                // skydive::elog!("ref {}", String::from_utf8_lossy(ref_seq));
                // skydive::elog!("alt {}", String::from_utf8_lossy(&alt_seq));
            }
        }
    }

    String::from_utf8_lossy(&alt_seq).to_string()
}

fn assemble_haplotype(ref_seqs: &Vec<Vec<u8>>, reads: &Vec<Vec<u8>>) -> String {
    let mut sg = spoa::Graph::new();

    let oriented_reads = orient_reads(ref_seqs, reads);

    let mut la = spoa::AlignmentEngine::new(AlignmentType::kOV, 5, -4, -8, -6, -8, -4);
    oriented_reads.iter().filter(|read| read.len() > 500).for_each(|lr_seq| {
        let seq = lr_seq.clone();
        let seq_cstr = std::ffi::CString::new(seq.clone()).unwrap();
        let seq_qual = std::ffi::CString::new(vec![b'I'; seq.len()]).unwrap();

        let a = la.align(seq_cstr.as_ref(), &sg);
        sg.add_alignment(&a, seq_cstr.as_ref(), seq_qual.as_ref());
    });

    let consensus_cstr = sg.consensus();

    consensus_cstr.to_str().unwrap().to_string()
}

fn orient_reads(ref_seqs: &Vec<Vec<u8>>, reads: &Vec<Vec<u8>>) -> Vec<Vec<u8>> {
    let mut ref_kmers = HashSet::new();
    for ref_seq in ref_seqs {
        for kmer in ref_seq.windows(17) {
            ref_kmers.insert(kmer);
        }
    }

    let mut sorted_reads = reads.clone();
    sorted_reads.sort_by(|a, b| b.len().cmp(&a.len()));

    let mut oriented_reads = Vec::new();

    for i in 0..sorted_reads.len() {
        if oriented_reads.is_empty() {
            oriented_reads.push(sorted_reads[i].clone());
        } else {
            let last_read_fw = oriented_reads[oriented_reads.len() - 1].clone();
            let last_read_kmers = last_read_fw.windows(17).collect::<Vec<_>>();
            // let last_read_kmers = &ref_kmers;

            let this_read_fw = &sorted_reads[i];
            let this_read_rc = this_read_fw.reverse_complement();

            let fw_overlaps = this_read_fw.windows(17).filter(|kmer| last_read_kmers.contains(kmer)).count();
            let rc_overlaps = this_read_rc.windows(17).filter(|kmer| last_read_kmers.contains(kmer)).count();

            if fw_overlaps > rc_overlaps {
                oriented_reads.push(this_read_fw.clone());
            } else {
                oriented_reads.push(this_read_rc);
            }
        }
    }

    oriented_reads
}

fn cluster_reads(l: &LdBG, all_seqs: &[Vec<u8>]) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let g = l.traverse_all_kmers();
    let bubbles = skydive::ldbg::find_all_superbubbles(&g);

    skydive::elog!("Found {} bubbles.", bubbles.len());

    let lr_mat = assign_reads_to_bubbles(&bubbles, &all_seqs, l, &g);
    let (h1, h2) = skydive::wmec::phase(&lr_mat);

    let all_mat = assign_reads_to_bubbles(&bubbles, &all_seqs, l, &g);

    let mut reads_hap1 = Vec::new();
    let mut reads_hap2 = Vec::new();
    for (read_index, read) in all_mat.reads.iter().enumerate() {
        let read_u8 = read.iter().map(|x| if x.is_some() { x.unwrap() } else { 2 }).collect::<Vec<u8>>();

        let mut h1_matches = 0;
        let mut h2_matches = 0;
        let mut total_valid = 0;

        for (i, &allele) in read_u8.iter().enumerate() {
            if allele == 0 || allele == 1 {
                total_valid += 1;
                if allele == h1[i] {
                    h1_matches += 1;
                }
                if allele == h2[i] {
                    h2_matches += 1;
                }
            }
        }

        let h1_similarity = if total_valid > 0 { h1_matches as f64 / total_valid as f64 } else { 0.0 };
        let h2_similarity = if total_valid > 0 { h2_matches as f64 / total_valid as f64 } else { 0.0 };

        if h1_similarity > 0.8 || h2_similarity > 0.8 {
            if h1_similarity > h2_similarity {
                reads_hap1.push(all_seqs[read_index].clone());
            } else {
                reads_hap2.push(all_seqs[read_index].clone());
            }
        }
    }

    (reads_hap1, reads_hap2)
}

fn assign_reads_to_bubbles(bubbles: &LinkedHashMap<(NodeIndex, NodeIndex), Vec<NodeIndex>>, lr_seqs: &[Vec<u8>], l: &LdBG, g: &petgraph::Graph<String, f32>) -> WMECData {
    let mut reads = vec![vec![None; bubbles.len()]; lr_seqs.len()];
    let mut confidences = vec![vec![None; bubbles.len()]; lr_seqs.len()];

    let cn_kmers: HashMap<Vec<u8>, Vec<usize>> = lr_seqs
        .iter()
        .enumerate()
        .flat_map(|(read_index, seq)| {
            seq.windows(l.kmer_size)
                .map(move |kmer| (skydive::utils::canonicalize_kmer(kmer), read_index))
        })
        .fold(HashMap::new(), |mut acc, (kmer, read_idx)| {
            acc.entry(kmer).or_default().push(read_idx);
            acc
        });

    for (bubble_index, ((in_node, out_node), interior)) in bubbles.iter().enumerate() {
        let paths_fwd = petgraph::algo::all_simple_paths::<Vec<_>, _>(&g, *in_node, *out_node, 0, Some(interior.len()));
        let paths_rev = petgraph::algo::all_simple_paths::<Vec<_>, _>(&g, *out_node, *in_node, 0, Some(interior.len()));

        let mut path_info = paths_fwd.chain(paths_rev).map(|path| {
            let path_kmers = path
                .iter()
                .map(|node| g.node_weight(*node).unwrap().as_bytes().to_vec())
                .map(|kmer| skydive::utils::canonicalize_kmer(&kmer))
                .collect::<Vec<Vec<u8>>>();

            let mut read_index_counts = HashMap::new();

            let mut path_score = 1.0;
            for kmer in &path_kmers {
                if let Some(read_indices) = cn_kmers.get(kmer) {
                    for &read_index in read_indices {
                        *read_index_counts.entry(read_index).or_insert(0) += 1;
                    }
                }

                let cn_kmer = skydive::utils::canonicalize_kmer(kmer);
                path_score *= l.scores.get(&cn_kmer).unwrap_or(&1.0);
            }

            path_score = path_score.powf(1.0 / path_kmers.len() as f32);
            if path_score.is_nan() {
                path_score = 0.0;
            }

            let path_length = path_kmers.len();
            let read_index_counts_filtered: BTreeMap<usize, i32> = read_index_counts.iter()
                .filter(|(_, &count)| count as f32 > 0.8 * path_length as f32)
                .map(|(&index, &count)| (index, count))
                .collect();

            let read_indices = read_index_counts_filtered.keys().cloned().collect::<Vec<_>>();

            (path_kmers, read_indices, path_score)
        })
        .collect::<Vec<_>>();

        path_info.sort_by(|(_, _, score_a), (_, _, score_b)| score_b.partial_cmp(score_a).unwrap());

        for (path_index, (path_kmers, read_indices, path_score)) in path_info.iter().take(2).enumerate() {
            let mut allele = String::new();
            for path_kmer in path_kmers {
                let kmer = std::str::from_utf8(path_kmer).unwrap();
                if allele.is_empty() {
                    allele.push_str(kmer);
                } else {
                    allele.push_str(&kmer.chars().last().unwrap().to_string());
                }
            }

            let path_qual = (-10.0 * (1.0 - path_score).log10()) as u32;

            // crate::elog!("[{} {}]: {:.3} {} {} {:?}", bubble_index, path_index, path_score, path_qual, allele, read_indices);

            for read_index in read_indices {
                reads[*read_index][bubble_index] = Some(path_index as u8);
                confidences[*read_index][bubble_index] = Some(path_qual);
            }
        }
    }

    let mat = WMECData::new(reads, confidences);

    mat
}

fn create_fully_phased_haplotypes(lr_msas: &Vec<String>, h1: &Vec<u8>) -> (String, String) {
    let mut index1 = 0;
    let mut hap_index = 0;

    let mut hap1 = String::new();
    let mut hap2 = String::new();

    while index1 < lr_msas[0].len() {
        let combined_base_counts = allele_counts(lr_msas, index1, index1+1);
        // let bases = allele_indices(lr_msas, index1, index1+1);

        if combined_base_counts.is_empty() {
            index1 += 1;
        } else if combined_base_counts.len() == 1 {
            let base = combined_base_counts.keys().next().unwrap();
            hap1.push(base.chars().next().unwrap());
            hap2.push(base.chars().next().unwrap());

            index1 += 1;
        } else {
            let mut index2 = index1;
            let mut allele_base_counts = allele_counts(lr_msas, index2, index2+1);

            while index2 < lr_msas[0].len() && allele_base_counts.len() > 1 {
                index2 += 1;
                allele_base_counts = allele_counts(lr_msas, index2, index2+1);
            }

            let allele_counts = allele_counts(lr_msas, index1, index2)
                .into_iter()
                .sorted_by(|a, b| b.1.cmp(&a.1))
                .take(2)
                .collect::<BTreeMap<_, _>>();

            if allele_counts.len() == 1 {
                for (allele, _) in allele_counts {
                    hap1.push_str(&allele);
                    hap2.push_str(&allele);
                }
            } else {
                // let alleles = get_allele_indices(lr_msas, sr_msas, index1, index2);

                if hap_index < h1.len() {
                    let mut phase = h1[hap_index] == 1;
                    for (allele, _) in allele_counts {
                        if !phase {
                            hap1.push_str(&allele);
                        } else {
                            hap2.push_str(&allele);
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

fn create_read_allele_matrix(lr_msas: &Vec<String>) -> Vec<BTreeMap<usize, String>> {
    let mut matrix = Vec::new();

    let mut index1 = 0;
    while index1 < lr_msas[0].len() {
        let combined_base_counts = allele_counts(lr_msas, index1, index1+1);

        if combined_base_counts.len() > 1 {
            let mut index2 = index1;
            let mut allele_base_counts = allele_counts(lr_msas, index2, index2+1);

            while index2 < lr_msas[0].len() && allele_base_counts.len() > 1 {
                index2 += 1;
                allele_base_counts = allele_counts(lr_msas, index2, index2+1);
            }

            let allele_counts = allele_counts(lr_msas, index1, index2)
                .into_iter()
                .sorted_by(|a, b| b.1.cmp(&a.1))
                .take(2)
                .collect::<BTreeMap<_, _>>();

            // println!("{} {} {:?}", index1, index2, allele_counts);

            if allele_counts.len() == 2 {
                let mut column = BTreeMap::new();

                let alleles = allele_indices(lr_msas, index1, index2);

                let mut allele_index = 0;
                for (allele, _) in &allele_counts {
                    alleles.iter().enumerate().for_each(|(i, a)| {
                        if *a == *allele {
                            // column.insert(i, allele.clone());
                            column.insert(i, allele_index.to_string());
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

fn allele_indices(lr_msas: &Vec<String>, index1: usize, index2: usize) -> Vec<String> {
    let alleles = lr_msas.iter()
        .map(|msa| msa[index1..index2].to_string().replace(" ", ""))
        .collect::<Vec<String>>();
    alleles
}

fn allele_counts(lr_msas: &Vec<String>, index1: usize, index2: usize) -> BTreeMap<String, i32> {
    let combined_allele_counts = lr_msas.iter()
        .map(|msa| msa[index1..index2].to_string().replace(" ", ""))
        .filter(|allele| !allele.is_empty())
        .fold(BTreeMap::new(), |mut counts, base| {
            *counts.entry(base).or_insert(0) += 1;
            counts
        });
    combined_allele_counts
}
