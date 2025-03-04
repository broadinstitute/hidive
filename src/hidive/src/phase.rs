use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::{fs::File, path::PathBuf, io::Write};

use itertools::Itertools;

use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::{Format, Header, Writer};
use rust_htslib::{bam, bam::{FetchDefinition, Read}};
use minimap2::Aligner;

use skydive::ldbg::LdBG;
use skydive::wmec::*;
use spoa::AlignmentType;

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
    let mut fasta = skydive::stage::open_fasta(&reference_seq_url).unwrap();

    // Parse BAM file path.
    let bam_urls = skydive::parse::parse_file_names(&[bam_path.clone()]);
    let bam_url = bam_urls.iter().next().unwrap();

    // Iterate over loci
    let loci = skydive::parse::parse_loci(loci_list, 0).into_iter().collect::<Vec<_>>();

    // Initialize VCF header
    // let vcf_header1 = initialize_vcf_header(&fasta, sample_name);
    // let mut vcf1 = Writer::from_path(PathBuf::from("phased1.vcf"), &vcf_header1, true, Format::Vcf).unwrap();

    // let vcf_header2 = initialize_vcf_header(&fasta, sample_name);
    // let mut vcf2 = Writer::from_path(PathBuf::from("phased2.vcf"), &vcf_header2, true, Format::Vcf).unwrap();

    // Initialize BAM header and writer for output
    let input_bam = skydive::stage::open_bam(&bam_url).unwrap();
    let mut bam1_writer = bam::Writer::from_path(format!("{}.phased1.bam", output.display()), &bam::Header::from_template(&input_bam.header()), bam::Format::Bam).unwrap();
    let mut bam2_writer = bam::Writer::from_path(format!("{}.phased2.bam", output.display()), &bam::Header::from_template(&input_bam.header()), bam::Format::Bam).unwrap();

    for (chr, start, stop, name) in loci {
        skydive::elog!("Processing locus {} ({}:{}-{})...", name, chr, start, stop);

        // The BAM reader gets renewed for each locus, but it's fast to open.
        let bam = skydive::stage::open_bam(&bam_url).unwrap();

        let (matrix, mloci, read_map) = prepare_matrix(&chr, start, stop, bam, &fasta);

        // for (read_idx, record) in &read_map {
        //     skydive::elog!("read: {} {}", read_idx, String::from_utf8_lossy(record.qname()));
        // }

        skydive::elog!(" - phasing {} variants...", mloci.len());

        let (h1_reads, h2_reads, h1, h2) = phase_variants(&matrix);

        for read in h1_reads.iter().filter(|idx| read_map.contains_key(idx)).map(|idx| read_map.get(idx).unwrap()) {
            bam1_writer.write(read).unwrap();
        }

        for read in h2_reads.iter().filter(|idx| read_map.contains_key(idx)).map(|idx| read_map.get(idx).unwrap()) {
            bam2_writer.write(read).unwrap();
        }

        // let h1_filtered = filter_variants(&mloci, &chr, &h1, &fasta);
        // let h2_filtered = filter_variants(&mloci, &chr, &h2, &fasta);

        // let mut h1_repaired = Vec::new();
        // let mut h2_repaired = Vec::new();

        // for ((l, a), b) in mloci.iter().zip(h1_filtered.iter()).zip(h2_filtered.iter()) {
        //     if a.is_some() && b.is_none() {
        //         h1_repaired.push(a.clone());
        //         h2_repaired.push(a.clone());
        //     } else if a.is_none() && b.is_some() {
        //         h1_repaired.push(b.clone());
        //         h2_repaired.push(b.clone());
        //     } else {
        //         h1_repaired.push(a.clone());
        //         h2_repaired.push(b.clone());
        //     }
        // }

        // let none_vec = vec![None; h1_filtered.len()];

        // add_variant_records(mloci.clone(), h1_repaired, none_vec.clone(), &mut fasta, chr.clone(), start, stop, &mut vcf1);
        // add_variant_records(mloci.clone(), none_vec.clone(), h2_repaired, &mut fasta, chr.clone(), start, stop, &mut vcf2);
    }
}

fn clean_variants(mloci: &Vec<u32>, h1: &Vec<Option<String>>, h2: &Vec<Option<String>>) -> (Vec<u32>, Vec<Option<String>>, Vec<Option<String>>) {
    let mut nloci: Vec<u32> = Vec::new();
    let mut h1_cleaned: Vec<Option<String>> = Vec::new();
    let mut h2_cleaned: Vec<Option<String>> = Vec::new();

    for i in 0..mloci.len() {
        if h1[i].is_some() || h2[i].is_some() {
            nloci.push(mloci[i]);
            h1_cleaned.push(h1[i].clone());
            h2_cleaned.push(h2[i].clone());
        }
    }

    (nloci, h1_cleaned, h2_cleaned)
}

/// Keeps only the longest insertion when multiple insertions are near each other
fn filter_variants(mloci: &Vec<u32>, chr: &String, variants: &[Option<String>], fasta: &rust_htslib::faidx::Reader) -> Vec<Option<String>> {
    let mut filtered = variants.to_vec();
    let mut include = vec![true; variants.len()];

    for i in 0..variants.len() {
        if include[i] && filtered[i].is_some() && filtered[i].as_ref().unwrap().len() > 100 {
            for j in (i + 1)..variants.len() {
                if include[j] && filtered[j].is_some() {
                    let i_len = filtered[i].as_ref().unwrap().len();
                    let j_len = filtered[j].as_ref().unwrap().len();
                    if (j_len as f64) >= (i_len as f64) * 0.5 && (j_len as f64) <= (i_len as f64) * 1.5 {
                        if filtered[j].as_ref().unwrap().len() > filtered[i].as_ref().unwrap().len() {
                            include[i] = false;
                        } else {
                            include[j] = false;
                        }
                    }
                }
            }
        }

        if include[i] && filtered[i].is_some() && filtered[i].as_ref().unwrap().len() == 2 {
            let ref_pos = mloci[i];
            let prev_ref_str = fasta.fetch_seq_string(chr, usize::try_from(ref_pos - 2).unwrap(), usize::try_from(ref_pos - 2).unwrap()).unwrap();
            let next_ref_str = fasta.fetch_seq_string(chr, usize::try_from(ref_pos).unwrap(), usize::try_from(ref_pos).unwrap()).unwrap();

            let first_allele_str = filtered[i].as_ref().unwrap().chars().next().unwrap().to_string();
            let last_allele_str = filtered[i].as_ref().unwrap().chars().last().unwrap().to_string();

            // skydive::elog!(" - examine {} {} {} {} {}", ref_pos, prev_ref_str, first_allele_str, last_allele_str, next_ref_str);

            if prev_ref_str == first_allele_str || next_ref_str == last_allele_str {
                include[i] = false;

                // skydive::elog!(" - filtered");
            }
        }
    }

    for i in 0..variants.len() {
        if !include[i] {
            filtered[i] = None;
        }
    }

    filtered
}

fn assemble_haplotype(reads: &Vec<Vec<u8>>) -> String {
    let mut sg = spoa::Graph::new();

    // let oriented_reads = orient_reads(ref_seqs, reads);

    let mut la = spoa::AlignmentEngine::new(AlignmentType::kOV, 5, -4, -8, -6, -8, -4);
    reads.iter().filter(|read| read.len() > 500).for_each(|lr_seq| {
        let seq = lr_seq.clone();
        let seq_cstr = std::ffi::CString::new(seq.clone()).unwrap();
        let seq_qual = std::ffi::CString::new(vec![b'I'; seq.len()]).unwrap();

        let a = la.align(seq_cstr.as_ref(), &sg);
        sg.add_alignment(&a, seq_cstr.as_ref(), seq_qual.as_ref());
    });

    let consensus_cstr = sg.consensus();

    consensus_cstr.to_str().unwrap().to_string()
}

fn initialize_vcf_header(fasta: &rust_htslib::faidx::Reader, sample_name: &String) -> Header {
    let mut vcf_header = Header::new();

    for seq_name in fasta.seq_names().unwrap() {
        let header_contig_line = format!("##contig=<ID={},length={}>", seq_name, fasta.fetch_seq_len(&seq_name));
        vcf_header.push_record(header_contig_line.as_bytes());
    }

    let header_end_line = r#"##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">"#;
    vcf_header.push_record(header_end_line.as_bytes());

    let header_gq_line = r#"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">"#;
    vcf_header.push_record(header_gq_line.as_bytes());

    let header_gt_line = r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#;
    vcf_header.push_record(header_gt_line.as_bytes());

    vcf_header.push_sample(sample_name.as_bytes());
    vcf_header
}

fn prepare_matrix_old(chr: &String, start: u64, stop: u64, mut bam: rust_htslib::bam::IndexedReader, fasta: &rust_htslib::faidx::Reader) -> (Vec<BTreeMap<usize, (String, u8)>>, Vec<u32>) {
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

                    // skydive::elog!("{} {}", pileup.pos(), len);

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

fn prepare_matrix(chr: &String, start: u64, stop: u64, mut bam: rust_htslib::bam::IndexedReader, fasta: &rust_htslib::faidx::Reader) -> (Vec<BTreeMap<usize, (String, u8)>>, Vec<u32>, BTreeMap<usize, Record>) {
    //                            pos           read    allele  qual
    let mut alleles_map: BTreeMap<u32, BTreeMap<usize, (String, u8)>> = BTreeMap::new();
    let mut read_ids = HashMap::new();
    let mut read_map = BTreeMap::new();

    let _ = bam.fetch(FetchDefinition::RegionString(chr.as_bytes(), start as i64, stop as i64));
    for (read_index, read) in bam.records().filter(|r| r.is_ok()).flatten().enumerate() {
        read_ids.insert(String::from_utf8_lossy(&read.qname()).to_string(), read_index);
        read_map.insert(read_index, read.clone());

        let mut read_pos = 0;
        let mut ref_pos = read.pos() as u32;

        for ce in read.cigar().iter() {
            match ce {
                Cigar::Ins(len) => {
                    if start as u32 <= ref_pos + 1 && ref_pos + 1 < stop as u32 {
                        let ref_str = fasta.fetch_seq_string(chr, usize::try_from(ref_pos - 1).unwrap(), usize::try_from(ref_pos - 1).unwrap()).unwrap();
                        let ref_base = ref_str.as_bytes()[0];

                        let qpos_start = read_pos as usize;
                        let qpos_stop = read_pos as usize + *len as usize;
                        let qual = read.qual()[qpos_start..qpos_stop].to_vec();
                        // let q = *qual.iter().min().unwrap();
                        let q = 10;
                        let seq = [&[ref_base], &read.seq().as_bytes()[qpos_start..qpos_stop]].concat();

                        let allele_map = alleles_map.entry(ref_pos).or_insert_with(BTreeMap::new);

                        if let Some((existing_seq, _)) = allele_map.get(&read_index) {
                            if seq.len() > existing_seq.len() {
                                allele_map.insert(read_index, (String::from_utf8_lossy(&seq).to_string(), q));
                            }
                        } else {
                            allele_map.insert(read_index, (String::from_utf8_lossy(&seq).to_string(), q));
                        }

                        // read_ids.insert(String::from_utf8_lossy(&read.qname()).to_string(), read_index);
                        // read_map.insert(read_index, read.clone());
                    }

                    read_pos += len;
                }
                Cigar::Del(len) => {
                    if start as u32 <= ref_pos + 1 && ref_pos + 1 < stop as u32 {
                        let ref_str = fasta.fetch_seq_string(chr, usize::try_from(ref_pos - 1).unwrap(), usize::try_from(ref_pos - 1).unwrap()).unwrap();
                        let ref_base = ref_str.as_bytes()[0];

                        let qpos_start = read_pos as usize;
                        // let q = read.qual()[qpos_start];
                        let q = 10;
                        let seq = vec![ref_base; 1].into_iter().chain(std::iter::repeat(b'-').take(*len as usize)).collect::<Vec<_>>();

                        let allele_map = alleles_map.entry(ref_pos - 1).or_insert_with(BTreeMap::new);

                        if let Some((existing_seq, _)) = allele_map.get(&read_index) {
                            if seq.len() > existing_seq.len() {
                                allele_map.insert(read_index, (String::from_utf8_lossy(&seq).to_string(), q));
                            }
                        } else {
                            allele_map.insert(read_index, (String::from_utf8_lossy(&seq).to_string(), q));
                        }

                        // read_ids.insert(String::from_utf8_lossy(&read.qname()).to_string(), read_index);
                        // read_map.insert(read_index, read.clone());
                    }

                    ref_pos += len;
                }
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                    for i in 0..*len as usize {
                        if start as u32 <= ref_pos + i as u32 + 1 && ref_pos + i as u32 + 1 < stop as u32 {
                            let ref_str = fasta.fetch_seq_string(chr, usize::try_from(ref_pos + i as u32).unwrap(), usize::try_from(ref_pos + i as u32).unwrap()).unwrap();
                            let ref_base = ref_str.as_bytes()[0];

                            let read_base = read.seq().as_bytes()[read_pos as usize + i];
                            // let q = read.qual()[read_pos as usize + i];
                            let q = 10;

                            if ref_base != read_base {
                                let seq = vec![read_base];
                                let allele_map = alleles_map.entry(ref_pos + i as u32).or_insert_with(BTreeMap::new);

                                if let Some((existing_seq, _)) = allele_map.get(&read_index) {
                                    if seq.len() > existing_seq.len() {
                                        allele_map.insert(read_index, (String::from_utf8_lossy(&seq).to_string(), q));
                                    }
                                } else {
                                    allele_map.insert(read_index, (String::from_utf8_lossy(&seq).to_string(), q));
                                }

                                // read_ids.insert(String::from_utf8_lossy(&read.qname()).to_string(), read_index);
                                // read_map.insert(read_index, read.clone());
                            }
                        }
                    }

                    read_pos += len;
                    ref_pos += len;
                }
                Cigar::SoftClip(len) => {
                    if start as u32 <= ref_pos + 1 && ref_pos + 1 < stop as u32 {
                        let ref_str = fasta.fetch_seq_string(chr, usize::try_from(ref_pos - 1).unwrap(), usize::try_from(ref_pos - 1).unwrap()).unwrap();
                        let ref_base = ref_str.as_bytes()[0];

                        let qpos_start = read_pos as usize;
                        let qpos_stop = read_pos as usize + *len as usize;
                        let qual = read.qual()[qpos_start..qpos_stop].to_vec();
                        // let q = *qual.iter().min().unwrap();
                        let q = 10;
                        let seq = [&[ref_base], &read.seq().as_bytes()[qpos_start..qpos_stop]].concat();

                        let allele_map = alleles_map.entry(ref_pos).or_insert_with(BTreeMap::new);

                        if let Some((existing_seq, _)) = allele_map.get(&read_index) {
                            if seq.len() > existing_seq.len() {
                                allele_map.insert(read_index, (String::from_utf8_lossy(&seq).to_string(), q));
                            }
                        } else {
                            allele_map.insert(read_index, (String::from_utf8_lossy(&seq).to_string(), q));
                        }

                        // read_ids.insert(String::from_utf8_lossy(&read.qname()).to_string(), read_index);
                        // read_map.insert(read_index, read.clone());
                    }

                    read_pos += len;
                }
                Cigar::HardClip(_) => { }
                Cigar::RefSkip(len) => {
                    ref_pos += len;
                }
                Cigar::Pad(_) => { }
            }
        }
    }

    let _ = bam.fetch(FetchDefinition::RegionString(chr.as_bytes(), start as i64, stop as i64));
    for p in bam.pileup() {
        let pileup = p.unwrap();

        for alignment in pileup.alignments() {
            let record = alignment.record();
            let qname = String::from_utf8(record.qname().to_vec()).unwrap();

            if read_ids.contains_key(&qname) {
                let read_index = *read_ids.get(&qname).unwrap();

                if alleles_map.contains_key(&(pileup.pos())) && !alleles_map.get(&(pileup.pos())).unwrap().contains_key(&read_index) {
                    if alignment.qpos().is_none() || alignment.qpos().unwrap() >= record.seq().len() {
                        continue;
                    }

                    let base = record.seq()[alignment.qpos().unwrap()];
                    let seq = vec![base];
                    // let q = record.qual()[alignment.qpos().unwrap()];
                    let q = 10;

                    let allele_map = alleles_map.get_mut(&pileup.pos()).expect("Position not found in alleles_map");
                    allele_map.insert(read_index, (String::from_utf8_lossy(&seq).to_string(), q));
                }
            }
        }
    }

    let mut matrix = Vec::new();
    let mut mloci = Vec::new();

    let mut last_pos = 0;
    for (pos, mut allele_map) in alleles_map {
        // skydive::elog!("before {} {:?}", pos, allele_map);

        for id in read_map.keys() {
            let mut allele_set = HashMap::new();
            for (_, (allele, _)) in &allele_map {
                if !allele_set.contains_key(allele) {
                    allele_set.insert(allele.clone(), 0);
                }

                if let Some(x) = allele_set.get_mut(allele) {
                    *x += 1;
                }
            }

            let mut most_frequent_allele = allele_set.keys().next().unwrap().clone();
            let mut most_frequent_count = 0;
            for (allele, count) in allele_set {
                if count > most_frequent_count {
                    most_frequent_count = count;
                    most_frequent_allele = allele;
                }
            }

            if !allele_map.contains_key(id) {
                allele_map.insert(*id, (most_frequent_allele, 0));
            }
        }

        // skydive::elog!("after  {} {:?}", pos, allele_map);

        // let score = allele_map.values().map(|(a, s)| *s as f32).sum::<f32>() / allele_map.len() as f32;

        let allele_set = allele_map.values().map(|(a, _)| a.clone()).collect::<HashSet<String>>();

        // skydive::elog!("{:?}", allele_set);

        // if allele_map.len() == read_map.len() {
        // if score > 5.0 {
        if allele_set.len() == 2 && allele_set.iter().map(|a| a.len() - 1).sum::<usize>() == 0 && pos - last_pos > 20 {
            matrix.push(allele_map);
            mloci.push(pos);
        }

        last_pos = pos;
    }

    (matrix, mloci, read_map)
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
    let mut hap_alleles = HashMap::new();
    for ((pos, a1), a2) in mloci.iter().zip(h1.iter()).zip(h2.iter()) {
        let ref_base = fasta.fetch_seq_string(&chr, usize::try_from(*pos).unwrap(), usize::try_from(*pos).unwrap()).unwrap();

        if a1.is_some() && a2.is_some() {
            hap_alleles.insert(pos, (a1.clone().unwrap(), a2.clone().unwrap()));
        } else if a1.is_some() && a2.is_none() {
            hap_alleles.insert(pos, (a1.clone().unwrap(), a1.clone().unwrap()));
        } else if a1.is_none() && a2.is_some() {
            hap_alleles.insert(pos, (a2.clone().unwrap(), a2.clone().unwrap()));
        } else {
            // skydive::elog!(" - no variant for {}", pos);
        }
    }

    let mut ref_start: Option<(String, u32)> = None;

    for pos in start as u32..stop as u32 {
        let ref_base = fasta.fetch_seq_string(&chr, usize::try_from(pos).unwrap(), usize::try_from(pos).unwrap()).unwrap();

        if hap_alleles.contains_key(&pos) {
            /*
            if ref_start.is_some() {
                let (ref_str, ref_pos) = ref_start.unwrap();

                let mut ref_record = vcf.empty_record();
                let rid = vcf.header().name2rid(chr.as_bytes()).unwrap();
                ref_record.set_rid(Some(rid));
                ref_record.set_pos(ref_pos as i64);

                ref_record.push_info_integer("END".as_bytes(), &[pos as i32]).unwrap();
                ref_record.push_format_integer("GQ".as_bytes(), &[40]).unwrap();

                ref_record.set_alleles(&[ref_str.to_string().as_bytes(), "<*>".as_bytes()]).unwrap();

                let genotype = &[GenotypeAllele::Unphased(0), GenotypeAllele::Phased(0)];
                ref_record.push_genotypes(genotype).unwrap();

                vcf.write(&ref_record).unwrap();

                ref_start = None;
            }
            */

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

            if alleles.len() > 1 {
                // alleles.push("<*>".as_bytes());
                record.set_alleles(&alleles).unwrap();

                record.push_format_integer("GQ".as_bytes(), &[40]).unwrap();

                // Get allele indices
                let i1 = record.alleles().iter().position(|&a| a == a1.replace("-", "").as_bytes()).unwrap() as i32;
                let i2 = record.alleles().iter().position(|&a| a == a2.replace("-", "").as_bytes()).unwrap() as i32;

                // Set record genotype - note first allele is always unphased
                let genotype = &[GenotypeAllele::Unphased(i1), GenotypeAllele::Phased(i2)];
                record.push_genotypes(genotype).unwrap();

                // Trim unused alleles
                // record.trim_alleles();

                // Write record
                vcf.write(&record).unwrap()
            }
        } else {
            if ref_start.is_none() {
                ref_start = Some((ref_base.clone(), pos));
            }
        }
    }
}

fn phase_variants(matrix: &Vec<BTreeMap<usize, (String, u8)>>) -> (BTreeSet<usize>, BTreeSet<usize>, Vec<Option<String>>, Vec<Option<String>>) {
    let num_snps = matrix.len();
    let num_reads = matrix.iter().map(|m| m.keys().max().unwrap_or(&0) + 1).max().unwrap_or(0);

    skydive::elog!("num_reads: {}", num_reads);
    skydive::elog!("num_snps: {}", num_snps);

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

    let wmec_matrix = WMECData::new(reads.clone(), confidences);

    let (p1, p2, r1, r2) = skydive::wmec::phase_all(&wmec_matrix, 60, 15);

    // skydive::elog!("aa: {:?}", all_alleles);
    skydive::elog!("p1: {:?}", p1);
    skydive::elog!("p2: {:?}", p2);
    skydive::elog!("r1: {:?}", r1);
    skydive::elog!("r2: {:?}", r2);

    // Convert phased variants to allele strings as before
    let mut h1 = Vec::new();
    let mut h2 = Vec::new();
    for i in 0..p1.len() {
        let a1 = all_alleles[i].get(&p1[i]).cloned();
        let a2 = all_alleles[i].get(&p2[i]).cloned();

        h1.push(a1);
        h2.push(a2);
    }

    // Return both the read assignments and the phased variants
    (r1, r2, h1, h2)
}