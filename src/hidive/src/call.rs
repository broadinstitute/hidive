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

        let mut read_ids = HashMap::new();
        let mut matrix = Vec::new();
        let mut metadata: BTreeMap<u32, String> = BTreeMap::new();
        let mut mloci = Vec::new();

        let _ = bam.fetch(FetchDefinition::RegionString(chr.as_bytes(), start as i64, stop as i64));
        for p in bam.pileup() {
            let pileup = p.unwrap();

            if pileup.pos() + 1 < start as u32 || pileup.pos() + 1 > stop as u32 {
                continue;
            }

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
                        let qpos_start = alignment.qpos().unwrap();
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
                        let full_seq = [&[ref_base], seq.as_slice()].concat();
                        allele_map.insert(*read_ids.get(&qname).unwrap(), (String::from_utf8_lossy(&full_seq).to_string(), q));

                        is_variant = true;
                    },
                    rust_htslib::bam::pileup::Indel::None => ()
                }
            }

            if is_variant {
                matrix.push(allele_map);
                metadata.insert(pileup.pos() + 1, String::from_utf8(vec![ref_base]).unwrap());
                mloci.push(pileup.pos());
            }
        }

        skydive::elog!("Phasing {} variants...", mloci.len());

        let (h1, h2) = phase_variants(&matrix);

        let mut hap_alleles = HashMap::new();
        for ((pos, a1), a2) in mloci.iter().zip(h1.iter()).zip(h2.iter()) {
            if a1.is_some() && a2.is_some() {
                hap_alleles.insert(pos, (a1.clone().unwrap(), a2.clone().unwrap()));
            } else if a1.is_some() && a2.is_none() {
                hap_alleles.insert(pos, (a1.clone().unwrap(), a1.clone().unwrap()));
            } else if a1.is_none() && a2.is_some() {
                hap_alleles.insert(pos, (a2.clone().unwrap(), a2.clone().unwrap()));
            }
        }

        let fasta = skydive::stage::open_fasta(&reference_seq_url).unwrap();

        let mut hap1 = Vec::new();
        let mut hap2 = Vec::new();

        for pos in start as u32..stop as u32 {
            let ref_str = fasta.fetch_seq_string(&chr, usize::try_from(pos).unwrap(), usize::try_from(pos).unwrap()).unwrap();

            if hap_alleles.contains_key(&pos) {
                let a1 = hap_alleles.get(&pos).unwrap().0.clone();
                let a2 = hap_alleles.get(&pos).unwrap().1.clone();

                hap1.push(a1);
                hap2.push(a2);
            } else {
                hap1.push(ref_str.clone());
                hap2.push(ref_str.clone());
            }
        }

        for i in 0..hap1.len() {
            if hap1[i].contains("-") {
                let len = hap1[i].len();

                for j in 1..len {
                    hap1[i+j] = "".to_string();
                }
            }

            if hap2[i].contains("-") {
                let len = hap2[i].len();

                for j in 1..len {
                    hap2[i+j] = "".to_string();
                }
            }
        }

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

    let (p1, p2) = skydive::wmec::phase_all(&wmec_matrix, 30, 20);

    let mut h1 = Vec::new();
    let mut h2 = Vec::new();
    for i in 0..p1.len() {
        let a1 = all_alleles[i].get(&p1[i]).cloned();
        let a2 = all_alleles[i].get(&p2[i]).cloned();

        h1.push(a1);
        h2.push(a2);
    }

    (h1, h2)
}