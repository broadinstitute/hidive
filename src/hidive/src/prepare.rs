use std::path::PathBuf;
use std::str;
use regex::Regex;

use bio::io::fasta::IndexedReader;
use bio::io::gff;
use bio::utils::Interval;

use debruijn::*;
use debruijn::kmer::*;

use skydive;

const TOTAL_INTERVAL_LENGTH_LIMIT: u64 = 100_000_000;

pub fn start(output: PathBuf, locus: &Option<Vec<String>>, gff: Option<PathBuf>, fasta: PathBuf) {
    let faidx = IndexedReader::from_file(&fasta).unwrap();

    let (chrs, intervals) = get_intervals(locus, &faidx);

    let total_interval_length = sum_interval_lengths(&intervals);
    if total_interval_length > TOTAL_INTERVAL_LENGTH_LIMIT {
        panic!("Total interval length is too large ({} bp); must be less than {} bp.", total_interval_length, TOTAL_INTERVAL_LENGTH_LIMIT);
    }

    let seqs = get_sequences(&chrs, &intervals, faidx);
    for seq in &seqs {
        // let s = match str::from_utf8(seq) {
        //     Ok(v) => v,
        //     Err(e) => panic!("Invalid UTF-8 sequence: {}", e),
        // };
        // println!("  -- len: {}", s.len());
        // println!("  -- seq: {}", s);

        construct_ldbg(seq);
    }

    // get_gff_records(&chrs, &intervals, gff);

    /*
    [x] vec loci (chrs and intervals (equal length))
    [x] - If locus vector is empty:
    [x]  - Add all loci to loci vec
    [x] - Else
    [x]  - Add all contig coordinates to loci vec

    [x] Abort if we try to add more than 100 megabases of data to the vector

    [x] vec sequences
    [x] for locus in loci:
    [ ] - if gff provided:
    [ ]   - get gff records
    [x] - get sequence
    [ ] - construct graph


    */
}

// fn get_gff_records(chrs: &[String], intervals: &[Interval<u64>], gff: Option<PathBuf>) {
//     if !gff.is_none() {
//         let mut reader = gff::Reader::from_file(gff.unwrap(), gff::GffType::GFF3).unwrap();

//         for record in reader.records() {
//             println!("{}", record.unwrap().feature_type());
//         }
//     }
// }

fn construct_ldbg(seq: &Vec<u8>) {
    let all_kmers = Kmer15::kmers_from_ascii(seq);

    for kmer in all_kmers {
        println!("{}", kmer.to_string());
    }

    // let dbg = DebruijnGraph<Kmer15, u16>::
}

fn get_sequences(chrs: &[String], intervals: &[Interval<u64>], mut faidx: IndexedReader<std::fs::File>) -> Vec<Vec<u8>> {
    let mut seqs: Vec<Vec<u8>> = Vec::new();

    for (i, (chr, interval)) in chrs.iter().zip(intervals.iter()).enumerate() {
        println!("{}: ({}:{}-{})", i, chr, interval.start, interval.end);

        let mut seq: Vec<u8> = Vec::new();

        faidx
            .fetch(chr, interval.start, interval.end)
            .expect("Couldn't fetch interval");

        faidx
            .read(&mut seq)
            .expect("Couldn't read the interval");

        seqs.push(seq);
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
            let m = l.replace(",", "");
            let pieces: Vec<&str> = sep.split(&m).collect();

            let interval = Interval::new(pieces[1].parse::<u64>().unwrap()..pieces[2].parse::<u64>().unwrap()).unwrap();

            chrs.push(String::from(pieces[0]));
            intervals.push(interval);
        };
    }

    (chrs, intervals)
}

fn sum_interval_lengths(intervals: &[Interval<u64>]) -> u64 {
    let mut sum: u64 = 0;

    for interval in intervals {
        sum += interval.end - interval.start + 1;
    }

    sum
}