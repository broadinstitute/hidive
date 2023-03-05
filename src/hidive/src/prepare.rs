use std::fs::File;
use std::path::PathBuf;
use std::str;
use regex::Regex;

use bio::io::fasta::IndexedReader;
use bio::utils::Interval;
// use bio::io::gff;

use debruijn::dna_string::DnaString;
use debruijn::*;
use debruijn::kmer::*;

use skydive;
use skydive::dbg::DeBruijnGraph;

const TOTAL_INTERVAL_LENGTH_LIMIT: u64 = 100_000_000;

pub fn start(_output: PathBuf, locus: &Option<Vec<String>>, _gff: Option<PathBuf>, fasta: PathBuf) {
    // Load sequences and intervals (if any).
    let faidx = IndexedReader::from_file(&fasta).unwrap();
    let (chrs, intervals) = get_intervals(locus, &faidx);

    // Convert sequences to a form that we can use, and break sequences at 'N's.
    let fwd_seqs = get_sequences(&chrs, &intervals, faidx);

    // Add sequences to graph.
    let mut g = DeBruijnGraph::new();
    g.add_all(&fwd_seqs);

    // Assemble contiguous sequences.
    println!("{:?}", g);

    // let mut buffer = File::create(output).unwrap();
    // dbg.write_gfa(&mut buffer).unwrap();

    /*
    [x] vec loci (chrs and intervals (equal length))
    [x] - if locus vector is empty:
    [x]  - add all loci to loci vec
    [x] - else
    [x]  - add all contig coordinates to loci vec

    [x] abort if we try to add more than 100 megabases of data to the vector

    [x] vec sequences
    [x] for locus in loci:
    [ ] - if gff provided:
    [ ]   - get gff records
    [x] - get sequence
    [x] - construct graph
    [ ] - add links to graph
    */
}

// fn get_gff_records(chrs: &[String], intervals: &[Interval<u64>], gff: Option<PathBuf>) {
//     if !gff.is_none() {
//         let mut reader = gff::Reader::from_file(gff.unwrap(), gff::GffType::GFF3).unwrap();
// 
//         for record in reader.records() {
//             println!("{}", record.unwrap().feature_type());
//         }
//     }
// }

fn get_sequences(chrs: &[String], intervals: &[Interval<u64>], mut faidx: IndexedReader<std::fs::File>) -> Vec<DnaString> {
    let mut seqs: Vec<DnaString> = Vec::new();

    for (_, (chr, interval)) in chrs.iter().zip(intervals.iter()).enumerate() {
        let mut seq: Vec<u8> = Vec::new();

        faidx
            .fetch(chr, interval.start, interval.end)
            .expect("Couldn't fetch interval");

        faidx
            .read(&mut seq)
            .expect("Couldn't read the interval");

        for split_seq in seq.split(|ch| *ch == b'N') {
            if split_seq.len() > 0 {
                let split_seq_two_bit: Vec<u8> = split_seq.to_vec().into_iter().map(debruijn::base_to_bits).collect();
                seqs.push(DnaString::from_bytes(&split_seq_two_bit));
            }
        }
    }

    seqs
}

fn get_intervals(locus: &Option<Vec<String>>, faidx: &IndexedReader<std::fs::File>) -> (Vec<String>, Vec<Interval<u64>>) {
    let mut chrs = Vec::new();
    let mut intervals = Vec::new();

    if locus.as_ref().unwrap_or(&Vec::new()).is_empty() {
        for seq in faidx.index.sequences() {
            let interval = Interval::new(0..seq.len).unwrap();

            chrs.push(seq.name);
            intervals.push(interval);
        }
    } else {
        for l in locus.as_ref().unwrap() {
            let sep = Regex::new(r"[:-]").unwrap();
            let m = l.replace(',', "");
            let pieces: Vec<&str> = sep.split(&m).collect();

            let interval = Interval::new(pieces[1].parse::<u64>().unwrap()..pieces[2].parse::<u64>().unwrap()).unwrap();

            chrs.push(String::from(pieces[0]));
            intervals.push(interval);
        };
    }

    check_interval_length_limit(&intervals);

    (chrs, intervals)
}

fn check_interval_length_limit(intervals: &[Interval<u64>]) {
    let total_interval_length = sum_interval_lengths(intervals);
    if total_interval_length > TOTAL_INTERVAL_LENGTH_LIMIT {
        panic!("Total interval length is too large ({} bp); must be less than {} bp.", total_interval_length, TOTAL_INTERVAL_LENGTH_LIMIT);
    }
}

fn sum_interval_lengths(intervals: &[Interval<u64>]) -> u64 {
    let mut sum: u64 = 0;

    for interval in intervals {
        sum += interval.end - interval.start + 1;
    }

    sum
}