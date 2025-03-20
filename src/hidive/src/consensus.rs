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
    loci_list: &Vec<String>,
    reference_fasta_path: &PathBuf,
    short_read_fasta_path: &PathBuf,
    bam_path: &PathBuf,
) {
    // Get the system's temporary directory path
    let cache_path = std::env::temp_dir();
    skydive::elog!("Intermediate data will be stored at {:?}.", cache_path);

    // Parse reference sequence file path.
    // let reference_seq_urls = skydive::parse::parse_file_names(&[reference_fasta_path.clone()]);
    // let reference_seq_url = reference_seq_urls.iter().next().unwrap();
    // let fasta = skydive::stage::open_fasta(&reference_seq_url).unwrap();

    skydive::elog!("Processing short-read sample {}...", short_read_fasta_path.display());
    let all_sr_seqs = skydive::utils::read_fasta(&vec![short_read_fasta_path.clone()]);
    skydive::elog!(" - {} short reads loaded", all_sr_seqs.len());

    // Parse BAM file path.
    let bam_urls = skydive::parse::parse_file_names(&[bam_path.clone()]);
    let bam_url = bam_urls.iter().next().unwrap();

    // Iterate over loc
    let loci = skydive::parse::parse_loci(loci_list, 0).into_iter().collect::<Vec<_>>();

    // Open FASTA file for writing.
    let mut buf_writer = BufWriter::new(File::create(output).unwrap());
    let mut fasta_writer = fasta::Writer::new(&mut buf_writer);

    for (chr, start_pos, stop_pos, name) in loci {
        skydive::elog!("Processing locus {} ({}:{}-{})...", name, chr, start_pos, stop_pos);

        let mut sg = spoa::Graph::new();

        // The BAM reader gets renewed for each locus, but it's fast to open.
        let mut bam = skydive::stage::open_bam(&bam_url).unwrap();
        let _ = bam.fetch(FetchDefinition::RegionString(chr.as_bytes(), start_pos as i64, stop_pos as i64));

        for read in bam.records().filter(|r| r.is_ok()).flatten() {
            let mut la = spoa::AlignmentEngine::new(AlignmentType::kOV, 5, -4, -8, -6, -8, -4);

            let subseq = read.seq().as_bytes();

            let seq_cstr = std::ffi::CString::new(subseq.clone()).unwrap();
            let seq_qual = std::ffi::CString::new(vec![b'I'; subseq.len()]).unwrap();
            let a = la.align(seq_cstr.as_ref(), &sg);
            sg.add_alignment(&a, seq_cstr.as_ref(), seq_qual.as_ref());
        }

        let consensus_cstr = sg.consensus();

        let record = fasta::Record::with_attrs(name.as_str(), None, consensus_cstr.to_str().unwrap().as_bytes());
        fasta_writer.write_record(&record).unwrap();
    }

    // let _ = fasta_writer.flush();
}