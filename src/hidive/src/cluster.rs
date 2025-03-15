use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

use bio::io::fasta;
use itertools::Itertools;

use minimap2::Aligner;
use needletail::Sequence;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::Cigar;
use rust_htslib::{bam, bam::{FetchDefinition, Read}};

use hdbscan::{DistanceMetric, Hdbscan, HdbscanHyperParams};
use spoa::AlignmentType;

pub fn start(
    output: &PathBuf,
    sample_name: &String,
    from_loci_list: &Vec<String>,
    to_loci_list: &Vec<String>,
    reference_fasta_path: &PathBuf,
    bam_path: &PathBuf,
) {
    // Get the system's temporary directory path
    let cache_path = std::env::temp_dir();
    skydive::elog!("Intermediate data will be stored at {:?}.", cache_path);

    // Parse reference sequence file path.
    let reference_seq_urls = skydive::parse::parse_file_names(&[reference_fasta_path.clone()]);
    let reference_seq_url = reference_seq_urls.iter().next().unwrap();
    let fasta = skydive::stage::open_fasta(&reference_seq_url).unwrap();

    // Parse BAM file path.
    let bam_urls = skydive::parse::parse_file_names(&[bam_path.clone()]);
    let bam_url = bam_urls.iter().next().unwrap();

    // Iterate over loc
    let from_loci = skydive::parse::parse_loci(from_loci_list, 0).into_iter().collect::<Vec<_>>();
    let to_loci = skydive::parse::parse_loci(to_loci_list, 0).into_iter().collect::<Vec<_>>();

    // Open FASTA file for writing.
    let mut buf_writer = BufWriter::new(File::create(output).unwrap());
    let mut fasta_writer = fasta::Writer::new(&mut buf_writer);

    for ((from_chr, from_start, from_stop, from_name), (to_chr, to_start, to_stop, to_name)) in from_loci.iter().zip(to_loci.iter()) {
        skydive::elog!("Found locus {} ({}:{}-{}) -> {} ({}:{}-{})...", from_name, from_chr, from_start, from_stop, to_name, to_chr, to_start, to_stop);
    }

    for ((from_chr, from_start, from_stop, from_name), (to_chr, to_start, to_stop, to_name)) in from_loci.iter().zip(to_loci.iter()) {
        skydive::elog!("Processing locus {} ({}:{}-{}) -> {} ({}:{}-{})...", from_name, from_chr, from_start, from_stop, to_name, to_chr, to_start, to_stop);

        // Extract the reference sequence for the current locus.
        let reference_seq = fasta.fetch_seq_string(to_chr.clone(), usize::try_from(*to_start).unwrap(), usize::try_from(*to_stop).unwrap()).unwrap();

        // Initialize aligner
        let mut aligner = Aligner::builder()
            .map_hifi()
            .with_cigar()
            .with_seq(reference_seq.as_bytes())
            .unwrap();

        aligner.mapopt.flag |= minimap2::ffi::MM_F_EQX as i64;

        // The BAM reader gets renewed for each locus, but it's fast to open.
        let mut bam = skydive::stage::open_bam(&bam_url).unwrap();

        let mut read_vectors = Vec::new();
        let mut reads = Vec::new();

        let _ = bam.fetch(FetchDefinition::RegionString(from_chr.as_bytes(), *from_start as i64, *from_stop as i64));
        for read in bam.records().filter(|r| r.is_ok()).flatten().filter(|r| !r.is_secondary() && !r.is_supplementary()) {
            let mappings = aligner
                .map(&read.seq().as_bytes(), false, false, None, None, Some(read.qname()))
                .unwrap();

            for mapping in mappings {
                // skydive::elog!("{:?}", mapping);

                let mut variants = BTreeMap::new();
                let mut ref_pos = mapping.target_start as i64 + 1;

                let cigar = mapping.alignment.unwrap().cigar_str.unwrap();
                let mut trim_left = 0;
                let mut trim_right = 0;

                let re = regex::Regex::new(r"(\d+)([MIDSNX=])").unwrap();
                for (cap_idx, cap) in re.captures_iter(&cigar).enumerate() {
                    let cigar_len = cap[1].parse::<u32>().unwrap();
                    match &cap[2] {
                        "M" | "=" => { ref_pos += cigar_len as i64; },
                        "I" => {
                            if cigar_len > 5 {
                                variants.insert(ref_pos, 1);
                            }
                        },
                        "D" => {
                            if cigar_len > 5 {
                                variants.insert(ref_pos, 1);
                            }
                            ref_pos += cigar_len as i64;
                        },
                        "X" => {
                            variants.insert(ref_pos, 1);
                            ref_pos += cigar_len as i64;
                        },
                        "S" => {
                            if cap_idx == 0 {
                                trim_left = cigar_len as usize;
                            } else {
                                trim_right = cigar_len as usize;
                            }
                        }
                        _ => {}
                    }
                }

                let seq = if mapping.strand == minimap2::Strand::Forward { read.seq().as_bytes() } else { read.seq().as_bytes().reverse_complement() };
                let subseq = seq[trim_left..seq.len() - trim_right].to_vec();

                read_vectors.push(variants);
                reads.push(subseq);
            }
        }

        let all_positions = read_vectors.iter()
            .flat_map(|r| r.keys())
            .collect::<HashSet<_>>()
            .into_iter()
            .sorted()
            .cloned()
            .collect::<Vec<_>>();

        let pos_to_idx = all_positions.iter()
            .enumerate()
            .map(|(i, p)| (p, i))
            .collect::<BTreeMap<_, _>>();

        let mut data = Vec::with_capacity(read_vectors.len());
        for read in read_vectors.iter() {
            let mut row = vec![0.0; all_positions.len()];
            for (&pos, &value) in read {
                if let Some(&idx) = pos_to_idx.get(&pos) {
                    row[idx] = value as f32;
                }
            }
            data.push(row);
        }

        // Configure HDBSCAN parameters
        let params = HdbscanHyperParams::builder()
            .min_cluster_size(2)
            .allow_single_cluster(true)
            .dist_metric(DistanceMetric::Euclidean)
            .build();

        // Create and run clusterer
        let clusterer = Hdbscan::new(&data, params);
        let cluster = clusterer.cluster();

        if let Ok(labels) = cluster {
            let unique_labels = labels.iter().unique().collect::<Vec<_>>();

            for unique_label in unique_labels {
                let mut sg = spoa::Graph::new();
                let mut la = spoa::AlignmentEngine::new(AlignmentType::kOV, 5, -4, -8, -6, -8, -4);

                for (label, read) in labels.iter().zip(reads.iter()) {
                    if label == unique_label {
                        let subseq = read.clone();

                        let seq_cstr = std::ffi::CString::new(subseq.clone()).unwrap();
                        let seq_qual = std::ffi::CString::new(vec![b'I'; subseq.len()]).unwrap();
                        let a = la.align(seq_cstr.as_ref(), &sg);
                        sg.add_alignment(&a, seq_cstr.as_ref(), seq_qual.as_ref());
                    }
                }

                if *unique_label >= 0 {
                    let consensus_cstr = sg.consensus();
                    let name = format!("{}.{}.{}", from_name, sample_name, unique_label);

                    let record = fasta::Record::with_attrs(name.as_str(), None, consensus_cstr.to_str().unwrap().as_bytes());
                    fasta_writer.write_record(&record).unwrap();
                }
            }
        }
    }
}