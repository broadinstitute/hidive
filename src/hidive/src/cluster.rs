use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::PathBuf;

use itertools::Itertools;

use minimap2::Aligner;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::Cigar;
use rust_htslib::{bam, bam::{FetchDefinition, Read}};

use hdbscan::{DistanceMetric, Hdbscan, HdbscanHyperParams};
use spoa::AlignmentType;

pub fn start(
    output: &PathBuf,
    // sample_name: &String,
    from_loci_list: &Vec<String>,
    to_loci_list: &Vec<String>,
    reference_fasta_path: &PathBuf,
    bam_path: &PathBuf,
) {
    // Get the system's temporary directory path
    let cache_path = std::env::temp_dir();
    skydive::elog!("Intermediate data will be stored at {:?}.", cache_path);

    // Parse reference sequence file path.
    // let reference_seq_urls = skydive::parse::parse_file_names(&[reference_fasta_path.clone()]);
    // let reference_seq_url = reference_seq_urls.iter().next().unwrap();
    // let mut fasta = skydive::stage::open_fasta(&reference_seq_url).unwrap();
    let mut aligner = Aligner::builder()
        .map_hifi()
        .with_cigar()
        .with_index(reference_fasta_path, None)
        .unwrap();

    // Parse BAM file path.
    let bam_urls = skydive::parse::parse_file_names(&[bam_path.clone()]);
    let bam_url = bam_urls.iter().next().unwrap();

    // Iterate over loc
    let from_loci = skydive::parse::parse_loci(from_loci_list, 0).into_iter().collect::<Vec<_>>();
    let to_loci = skydive::parse::parse_loci(to_loci_list, 0).into_iter().collect::<Vec<_>>();

    for ((from_chr, from_start, from_stop, from_name), (to_chr, to_start, to_stop, to_name)) in from_loci.iter().zip(to_loci.iter()) {
        skydive::elog!("Processing locus {} ({}:{}-{}) -> {} ({}:{}-{})...", from_name, from_chr, from_start, from_stop, to_name, to_chr, to_start, to_stop);

        // The BAM reader gets renewed for each locus, but it's fast to open.
        let mut bam = skydive::stage::open_bam(&bam_url).unwrap();

        let mut read_vectors = Vec::new();
        let mut reads = Vec::new();

        let _ = bam.fetch(FetchDefinition::RegionString(from_chr.as_bytes(), *from_start as i64, *from_stop as i64));
        for read in bam.records().filter(|r| r.is_ok()).flatten().filter(|r| !r.is_secondary() && !r.is_supplementary()) {
            let mut variants = HashMap::new();
            let mut ref_pos = read.pos();

            let alignments = aligner
                .map(&read.seq().as_bytes(), false, false, None, None, Some(read.qname()))
                .unwrap();

            skydive::elog!("{:?}", alignments);

            for cigar_element in read.cigar().iter() {
                match cigar_element {
                    Cigar::Match(len) => {
                        ref_pos += *len as i64;
                    }
                    Cigar::Equal(len) => {
                        ref_pos += *len as i64;
                    }
                    Cigar::Diff(len) => {
                        variants.insert(ref_pos, 1);
                        ref_pos += *len as i64;
                    }
                    Cigar::Ins(len) => {
                        if *len > 5 {
                            variants.insert(ref_pos, 1);
                        }
                    }
                    Cigar::Del(len) => {
                        if *len > 5 {
                            variants.insert(ref_pos, 1);
                        }
                        ref_pos += *len as i64;
                    }
                    _ => {}
                }
            }

            read_vectors.push(variants);
            reads.push(read);
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
            .dist_metric(DistanceMetric::Euclidean) // Since we're passing a distance matrix
            .build();

        // Create and run clusterer
        let clusterer = Hdbscan::new(&data, params);
        let labels = clusterer.cluster().unwrap();

        // Initialize BAM writer for output
        let mut bam_writer = bam::Writer::from_path("clustered.bam", &bam::Header::from_template(&bam.header()), bam::Format::Bam).unwrap();

        for (label, read) in labels.iter().zip(reads.iter()) {
            let mut tagged_read = read.clone();
            let _ = tagged_read.push_aux(b"CL", bam::record::Aux::I32(*label));

            skydive::elog!("{} {} {:?}", label, String::from_utf8_lossy(read.qname()), tagged_read.aux(b"CL"));

            bam_writer.write(&tagged_read).unwrap();
        }

        let unique_labels = labels.iter().unique().collect::<Vec<_>>();

        for unique_label in unique_labels {
            let mut sg = spoa::Graph::new();
            let mut la = spoa::AlignmentEngine::new(AlignmentType::kOV, 5, -4, -8, -6, -8, -4);

            for (label, read) in labels.iter().zip(reads.iter()) {
                if label == unique_label {
                    let seq = read.seq().as_bytes();

                    let mut trim_left = 0;
                    // let first_cigar_element = read.cigar().first().unwrap();
                    match read.cigar().first().unwrap() {
                        Cigar::SoftClip(len) => {
                            trim_left = *len as usize;
                        }
                        _ => {}
                    }

                    let mut trim_right = 0;
                    match read.cigar().last().unwrap() {
                        Cigar::SoftClip(len) => {
                            trim_right = *len as usize;
                        }
                        _ => {}
                    }

                    let subseq = seq[trim_left..seq.len() - trim_right].to_vec();

                    let seq_cstr = std::ffi::CString::new(subseq.clone()).unwrap();
                    let seq_qual = std::ffi::CString::new(vec![b'I'; subseq.len()]).unwrap();
                    let a = la.align(seq_cstr.as_ref(), &sg);
                    sg.add_alignment(&a, seq_cstr.as_ref(), seq_qual.as_ref());
                }
            }

            if *unique_label >= 0 {
                let consensus_cstr = sg.consensus();
                println!(">{}\n{}", unique_label, consensus_cstr.to_str().unwrap());
            }
        }
    }
}