use std::collections::{BTreeMap, HashMap, HashSet};
use std::{fs::File, path::PathBuf, io::Write};

use itertools::Itertools;
use linked_hash_map::LinkedHashMap;
use minimap2::Aligner;
use needletail::Sequence;
use petgraph::graph::NodeIndex;
use rayon::prelude::*;
use indicatif::ParallelProgressIterator;

use rust_htslib::bam::{FetchDefinition, Read};
use skydive::ldbg::LdBG;
use skydive::mldbg::MLdBG;

use skydive::wmec::*;
use spoa::AlignmentType;

pub fn start(
    output: &PathBuf,
    loci_list: &Vec<String>,
    reference_fasta_path: &PathBuf,
    bam_path: &PathBuf,
) {
    // Get the system's temporary directory path
    let cache_path = std::env::temp_dir();
    skydive::elog!("Intermediate data will be stored at {:?}.", cache_path);

    // Parse reference sequence file path.
    let reference_seq_urls = skydive::parse::parse_file_names(&[reference_fasta_path.clone()]);
    let reference_seq_url = reference_seq_urls.iter().next().unwrap();

    // Parse BAM file path.
    let bam_urls = skydive::parse::parse_file_names(&[bam_path.clone()]);
    let bam_url = bam_urls.iter().next().unwrap();

    // Iterate over loci
    let loci = skydive::parse::parse_loci(loci_list, 0).into_iter().collect::<Vec<_>>();

    for (chr, start, stop, name) in loci {
        let fasta = skydive::stage::open_fasta(&reference_seq_url).unwrap();
        let mut bam = skydive::stage::open_bam(&bam_url).unwrap();

        skydive::elog!("Processing locus {} ({}:{}-{})...", name, chr, start, stop);
        // skydive::elog!("Reference sequence: {}", seq);

        let mut read_ids = HashMap::new();
        let mut matrix = Vec::new();
        let mut metadata: BTreeMap<u32, String> = BTreeMap::new();

        let _ = bam.fetch(FetchDefinition::RegionString(chr.as_bytes(), start as i64, stop as i64));
        for (cursor, p) in bam.pileup().enumerate() {
            let pileup = p.unwrap();

            let ref_str = fasta.fetch_seq_string(&chr, usize::try_from(pileup.pos()).unwrap(), usize::try_from(pileup.pos()).unwrap()).unwrap();
            let ref_base = ref_str.as_bytes()[0];

            let mut is_variant = false;

            let mut allele_map: BTreeMap<usize, (String, u8)> = BTreeMap::new();

            for alignment in pileup.alignments() {
                let record = alignment.record();
                let qname = String::from_utf8(record.qname().to_vec()).unwrap();

                if !alignment.is_del() && !alignment.is_refskip() {
                    let base = record.seq()[alignment.qpos().unwrap()];
                    let seq = vec![base];
                    let qual = vec![record.qual()[alignment.qpos().unwrap()]];
                    let q = *qual.iter().min().unwrap();

                    let len = read_ids.len();
                    read_ids.entry(qname.clone()).or_insert(len);
                    allele_map.insert(*read_ids.get(&qname).unwrap(), (String::from_utf8_lossy(&seq).to_string(), q));

                    if base != ref_base {
                        is_variant = true;
                    }
                }

                match alignment.indel() {
                    rust_htslib::bam::pileup::Indel::Ins(len) => {
                        let qpos_start = alignment.qpos().unwrap() + 1;
                        let qpos_stop = alignment.qpos().unwrap() + 1 + (len as usize);
                        let seq = record.seq().as_bytes()[qpos_start..qpos_stop].to_vec();
                        let qual = record.qual()[qpos_start..qpos_stop].to_vec();
                        let q = *qual.iter().min().unwrap();

                        let len = read_ids.len();
                        read_ids.entry(qname.clone()).or_insert(len);
                        allele_map.insert(*read_ids.get(&qname).unwrap(), (String::from_utf8_lossy(&seq).to_string(), q));

                        is_variant = true;
                    },
                    rust_htslib::bam::pileup::Indel::Del(len) => {
                        let seq = vec![b'-'; len as usize];
                        let qual = vec![record.qual()[alignment.qpos().unwrap()]];
                        let q = *qual.iter().min().unwrap();

                        let len = read_ids.len();
                        read_ids.entry(qname.clone()).or_insert(len);
                        allele_map.insert(*read_ids.get(&qname).unwrap(), (String::from_utf8_lossy(&seq).to_string(), q));

                        is_variant = true;
                    },
                    rust_htslib::bam::pileup::Indel::None => ()
                }
            }

            if is_variant {
                // skydive::elog!("Variant found at {}:{} {} ({:?})", chr, pileup.pos() + 1, ref_base, alleles);
                matrix.push(allele_map);
                // metadata.push((pileup.pos(), ref_base));
                // metadata.insert(cursor, ref_base);
                metadata.insert(pileup.pos() + 1, String::from_utf8(vec![ref_base]).unwrap());
            }
        }

        let (h1, h2) = phase_variants(&matrix);

        skydive::elog!("Haplotype 1: {:?} ({})", h1, h1.len());
        skydive::elog!("Haplotype 2: {:?} ({})", h2, h2.len());
        skydive::elog!("Metadata: {:?}", metadata);

        let fasta = skydive::stage::open_fasta(&reference_seq_url).unwrap();

        /*
        let seq = fasta
            .fetch_seq_string(&chr, usize::try_from(start).unwrap(), usize::try_from(stop - 1).unwrap())
            .unwrap()
            .as_bytes()
            .to_vec();
        */

        let mut hap1 = Vec::new();
        let mut hap2 = Vec::new();

        let mut cursor = 0;
        for pos in start as u32..stop as u32 {
            let ref_str = fasta.fetch_seq_string(&chr, usize::try_from(pos).unwrap(), usize::try_from(pos).unwrap()).unwrap();

            if metadata.contains_key(&pos) {
                let a1 = h1[cursor].clone();
                let a2 = h2[cursor].clone();

                hap1.push(a1.unwrap_or(ref_str.clone()));
                hap2.push(a2.unwrap_or(ref_str.clone()));

                cursor += 1;
            } else {
                hap1.push(ref_str.clone());
                hap2.push(ref_str.clone());
            }
        }


        /*
        let mut cursor = 0;
        for (pos, &base) in seq.iter().enumerate() {
            if metadata.contains_key(&pos) {
                let a1 = h1[cursor].clone();
                let a2 = h2[cursor].clone();

                hap1.push(a1.unwrap_or(String::from_utf8(vec![base]).unwrap()));
                hap2.push(a2.unwrap_or(String::from_utf8(vec![base]).unwrap()));

                cursor += 1;
            } else {
                hap1.push(String::from_utf8(vec![base]).unwrap());
                hap2.push(String::from_utf8(vec![base]).unwrap());
            }
        }
        */

        let h1 = hap1.join("").replace("-", "");
        let h2 = hap2.join("").replace("-", "");

        let mut output = File::create(output).expect("Failed to create output file");
        writeln!(output, ">h1").expect("Failed to write to output file");
        writeln!(output, "{}", h1).expect("Failed to write to output file");
        writeln!(output, ">h2").expect("Failed to write to output file");
        writeln!(output, "{}", h2).expect("Failed to write to output file");
    }
}

fn phase_variants(matrix: &Vec<BTreeMap<usize, (String, u8)>>) -> (Vec<Option<String>>, Vec<Option<String>>) {
    let num_snps = matrix.len();
    let num_reads = matrix.iter().map(|m| m.keys().max().unwrap_or(&0) + 1).max().unwrap_or(0);

    let mut reads = vec![vec![None; num_snps]; num_reads];
    let mut confidences = vec![vec![None; num_snps]; num_reads];
    let mut all_alleles = vec![HashMap::new(); num_snps];

    for (snp_idx, column) in matrix.iter().enumerate() {
        // Count frequency of each allele in this column
        let allele_counts = column
            .values()
            .map(|(allele, _)| allele)
            .fold(HashMap::new(), |mut counts, allele| {
                *counts.entry(allele.clone()).or_insert(0) += 1;
                counts
            })
            .into_iter()
            .sorted_by(|a, b| b.1.cmp(&a.1))
            .take(2)
            .collect::<Vec<_>>();

        let alleles = allele_counts 
            .iter()
            .map(|(a, _)| a)
            .cloned()
            .collect::<HashSet<_>>();

        let allele_map = alleles.iter().enumerate().map(|(i, a)| (a, i as u8)).collect::<HashMap<_, _>>();

        let mut index_map = HashMap::new();

        for (&read_idx, (allele, qual)) in column {
            if let Some(allele_idx) = allele_map.get(allele) {
                reads[read_idx][snp_idx] = Some(*allele_idx);
                confidences[read_idx][snp_idx] = Some(*qual as u32);

                index_map.insert(*allele_idx, allele.clone());
            }
        }

        all_alleles[snp_idx] = index_map;
    }

    let wmec_matrix = WMECData::new(reads, confidences);

    let (p1, p2) = skydive::wmec::phase(&wmec_matrix);

    let mut h1 = Vec::new();
    let mut h2 = Vec::new();
    for i in 0..p1.len() {
        // h1.push(all_alleles[i].get(&p1[i]).unwrap().clone());
        // h2.push(all_alleles[i].get(&p2[i]).unwrap().clone());

        let a1 = all_alleles[i].get(&p1[i]).cloned();
        let a2 = all_alleles[i].get(&p2[i]).cloned();

        h1.push(a1);
        h2.push(a2);
    }

    (h1, h2)
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
