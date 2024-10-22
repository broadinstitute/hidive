use std::collections::BTreeMap;
use std::{fs::File, path::PathBuf, io::Write};

use itertools::Itertools;
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
        .clean(0.2, 0.01)
        .build_links(&all_lr_seqs, true);

    skydive::elog!("Built MLdBG with {} k-mers.", m.kmers.len());

    let progress_bar = skydive::utils::default_bounded_progress_bar("Correcting reads", all_lr_seqs.len() as u64);

    let corrected_seqs =
        all_lr_seqs
        .par_iter()
        .progress_with(progress_bar)
        .map(|seq| m.correct_seq(seq))
        .flatten()
        .collect::<Vec<Vec<u8>>>();

    // let contigs = m.assemble_at_bubbles();

    // let mut fa_file = File::create(&output).unwrap();
    // for (i, contig) in corrected_seqs.iter().chain(contigs.iter()).enumerate() {
    //     let _ = writeln!(fa_file, ">contig_{}\n{}", i, String::from_utf8(contig.clone()).unwrap());
    // }

    skydive::elog!("Assembling diplotypes...");

    let mut la = spoa::AlignmentEngine::new(AlignmentType::kSW, 5, -10, -16, -12, -20, -8);

    let mut sg = spoa::Graph::new();
    for lr_seq in all_ref_seqs.iter().chain(corrected_seqs.iter()) {
        let seq_cstr = std::ffi::CString::new(lr_seq.clone()).unwrap();
        let seq_qual = std::ffi::CString::new(vec![b'I'; lr_seq.len()]).unwrap();
        let a = la.align(seq_cstr.as_ref(), &sg);

        sg.add_alignment(&a, seq_cstr.as_ref(), seq_qual.as_ref());
    }

    let msa_cstrs = sg.multiple_sequence_alignment(false);
    let msa_strings = msa_cstrs
        .iter()
        .map(|cstr| cstr.to_str().unwrap().to_string())
        .map(|msa_string| {
            let leading_dashes = msa_string.chars().take_while(|&c| c == '-').count();
            let trailing_dashes = msa_string.chars().rev().take_while(|&c| c == '-').count();
            let middle = &msa_string[leading_dashes..msa_string.len() - trailing_dashes];
            format!("{}{}{}", 
                " ".repeat(leading_dashes), 
                middle, 
                " ".repeat(trailing_dashes)
            )
        })
        .collect::<Vec<String>>();

    // for msa_string in &msa_strings {
    //     println!("{}", msa_string);
    // }

    let matrix = create_read_allele_matrix(&msa_strings);
    let wmec = create_wmec_matrix(&matrix);

    let (h1, _) = skydive::wmec::phase(&wmec);
    let (hap1, hap2) = create_fully_phased_haplotypes(&msa_strings, &h1);

    let mut output = File::create(output).expect("Failed to create output file");
    writeln!(output, ">hap1").expect("Failed to write to output file");
    writeln!(output, "{}", hap1.replace("-", "")).expect("Failed to write to output file");
    writeln!(output, ">hap2").expect("Failed to write to output file");
    writeln!(output, "{}", hap2.replace("-", "")).expect("Failed to write to output file");
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
