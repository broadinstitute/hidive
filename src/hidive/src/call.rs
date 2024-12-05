use std::collections::{BTreeMap, HashMap, HashSet};
use std::{fs::File, path::PathBuf, io::Write};

use itertools::Itertools;

use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::{Format, Header, Writer};
use rust_htslib::{bam, bam::{FetchDefinition, Read}};
use minimap2::Aligner;

use skydive::wmec::*;

pub fn start(
    output: &PathBuf,
    sample_name: &String,
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

    // Initialize VCF header
    let mut vcf_header = Header::new();

    let mut fasta = skydive::stage::open_fasta(&reference_seq_url).unwrap();
    for seq_name in fasta.seq_names().unwrap() {
        let header_contig_line = format!("##contig=<ID={},length={}>", seq_name, fasta.fetch_seq_len(&seq_name));
        vcf_header.push_record(header_contig_line.as_bytes());
    }

    let header_gt_line = r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#;
    vcf_header.push_record(header_gt_line.as_bytes());
    vcf_header.push_sample(sample_name.as_bytes());

    let mut vcf = Writer::from_path(output, &vcf_header, true, Format::Vcf).unwrap();

    for (chr, start, stop, name) in loci {
        skydive::elog!("Processing locus {} ({}:{}-{})...", name, chr, start, stop);

        // The BAM reader gets renewed for each locus, but it's fast to open.
        let bam = skydive::stage::open_bam(&bam_url).unwrap();

        let (matrix, mloci) = prepare_matrix(bam, &chr, start, stop, &fasta);

        skydive::elog!(" - phasing {} variants...", mloci.len());

        let (h1, h2) = phase_variants(&matrix);

        // let (hap1, hap2) = build_haplotypes(mloci, h1, h2, reference_seq_url, start, stop, chr);
        // write_haplotypes(&mut output_file, name, hap1, hap2);

        add_variant_records(mloci, h1, h2, &mut fasta, chr, start, stop, &mut vcf);
    }
}

fn prepare_matrix(mut bam: rust_htslib::bam::IndexedReader, chr: &String, start: u64, stop: u64, fasta: &rust_htslib::faidx::Reader) -> (Vec<BTreeMap<usize, (String, u8)>>, Vec<u32>) {
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

        let ref_str = fasta.fetch_seq_string(chr, usize::try_from(pileup.pos()).unwrap(), usize::try_from(pileup.pos()).unwrap()).unwrap();
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
    (matrix, mloci)
}

fn write_haplotypes(output_file: &mut File, name: String, hap1: String, hap2: String) {
    writeln!(output_file, ">hap1_{}", name).expect("Failed to write to output file");
    writeln!(output_file, "{}_{}", hap1, name).expect("Failed to write to output file");
    writeln!(output_file, ">hap2_{}", name).expect("Failed to write to output file");
    writeln!(output_file, "{}_{}", hap2, name).expect("Failed to write to output file");
}

fn build_haplotypes(mloci: Vec<u32>, h1: Vec<Option<String>>, h2: Vec<Option<String>>, reference_seq_url: &url::Url, start: u64, stop: u64, chr: String) -> (String, String) {
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
    (h1, h2)
}

fn add_variant_records(mloci: Vec<u32>, h1: Vec<Option<String>>, h2: Vec<Option<String>>, fasta: &mut rust_htslib::faidx::Reader, chr: String, start: u64, stop: u64, vcf: &mut rust_htslib::bcf::Writer) {
    // let fasta = skydive::stage::open_fasta(&reference_seq_url).unwrap();

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

    for pos in start as u32..stop as u32 {
        let ref_base = fasta.fetch_seq_string(&chr, usize::try_from(pos).unwrap(), usize::try_from(pos).unwrap()).unwrap();

        if hap_alleles.contains_key(&pos) {
            let mut record = vcf.empty_record();

            // Set contig name and record position
            let rid = vcf.header().name2rid(chr.as_bytes()).unwrap();
            record.set_rid(Some(rid));
            record.set_pos(pos as i64);

            // Get haplotype alleles
            let mut a1 = hap_alleles.get(&pos).unwrap().0.clone();
            let mut a2 = hap_alleles.get(&pos).unwrap().1.clone();

            // First determine the full length of the reference allele
            let mut ref_allele = ref_base.clone();
            if a1.contains("-") {
                let ref_full_str = fasta.fetch_seq_string(&chr, usize::try_from(pos).unwrap(), usize::try_from(pos + (a1.len() as u32) - 1).unwrap()).unwrap();
                if ref_full_str.len() > ref_allele.len() {
                    ref_allele = ref_full_str;
                }
            }
            if a2.contains("-") {
                let ref_full_str = fasta.fetch_seq_string(&chr, usize::try_from(pos).unwrap(), usize::try_from(pos + (a2.len() as u32) - 1).unwrap()).unwrap();
                if ref_full_str.len() > ref_allele.len() {
                    ref_allele = ref_full_str;
                }
            }

            // Rewrite insertions relative to the reference allele
            if ref_allele.len() > 1 && !a1.contains("-") && a1.len() > 1 {
                a1 = ref_base.clone() + &a1[1..] + &ref_allele[1..];
            }

            if ref_allele.len() > 1 && !a2.contains("-") && a2.len() > 1 {
                a2 = ref_base.clone() + &a2[1..] + &ref_allele[1..];
            }

            // If the alternate alleles was the same as the reference allele, replace them with the full length reference allele
            if ref_allele.len() > 1 && a1 == ref_base {
                a1 = ref_allele.clone();
            }

            if ref_allele.len() > 1 && a2 == ref_base {
                a2 = ref_allele.clone();
            }

            let mut alt_alleles = HashSet::new();
            alt_alleles.insert(a1.clone());
            alt_alleles.insert(a2.clone());

            let mut allele_set = alt_alleles.iter().map(|s| s.replace("-", "").into_bytes()).collect::<HashSet<_>>();
            allele_set.insert(ref_allele.to_string().into_bytes());

            // Make sure the reference allele is first
            let mut alleles_vec = allele_set.into_iter().collect::<Vec<_>>();
            alleles_vec.sort_by(|a, b| {
                if a == ref_allele.as_bytes() {
                    std::cmp::Ordering::Less
                } else if b == ref_allele.as_bytes() {
                    std::cmp::Ordering::Greater
                } else {
                    a.cmp(b)
                }
            });

            let alleles: Vec<&[u8]> = alleles_vec.iter().map(|v| v.as_slice()).collect();
            record.set_alleles(&alleles).unwrap();

            // Get allele indices
            let i1 = record.alleles().iter().position(|&a| a == a1.replace("-", "").as_bytes()).unwrap() as i32;
            let i2 = record.alleles().iter().position(|&a| a == a2.replace("-", "").as_bytes()).unwrap() as i32;

            // Set record genotype - note first allele is always unphased
            let genotype = &[GenotypeAllele::Unphased(i1), GenotypeAllele::Phased(i2)];
            record.push_genotypes(genotype).unwrap();

            // Trim unused alleles
            record.trim_alleles();

            // Write record
            vcf.write(&record).unwrap()
        } else {
            // hap1.push(ref_str.clone());
            // hap2.push(ref_str.clone());
        }
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

    let (p1, p2) = skydive::wmec::phase_all(&wmec_matrix, 40, 20);

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