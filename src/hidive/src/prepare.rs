use std::path::PathBuf;
use std::str;
use regex::Regex;

use bio::io::fasta::IndexedReader;
use bio::utils::Interval;

use needletail::Sequence;

use skydive;
use skydive::ldbg::LdBG;

const TOTAL_INTERVAL_LENGTH_LIMIT: u64 = 100_000_000;

pub fn start(_output: PathBuf, locus: &Option<Vec<String>>, _gff: Option<PathBuf>, fasta: PathBuf) {
    // Load sequences and intervals (if any).
    let faidx = IndexedReader::from_file(&fasta).expect("Unable to open specified FASTA file.");
    let (chrs, intervals) = get_intervals(locus, &faidx);

    // Convert sequences to a form that we can use, and break sequences at 'N's.
    let fwd_seqs = get_sequences(&chrs, &intervals, faidx);

    // Construct a linked de Bruijn graph from the sequences.
    let g: LdBG<15> = LdBG::from_sequences(fwd_seqs);

    for (kmer, record) in g.kmers {
        if record.coverage() > 1 {
            println!("{:?} {:?}", String::from_utf8(kmer), record.coverage());
        }
    }

    // println!("{:?}", g.kmers.len());

    // for fwd_seq in fwd_seqs {
    //     let s = String::from_utf8(fwd_seq).unwrap();
    //     println!("{}", s);
    // }

    // Assemble contiguous sequences.
    // println!("{:?}", g);

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
    [ ] - construct graph
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

/// Load sequence(s) from a FASTA file (optionally only those from specified loci), splitting on 'N's, and return them.
fn get_sequences(chrs: &[String], intervals: &[Interval<u64>], mut faidx: IndexedReader<std::fs::File>) -> Vec<Vec<u8>> {
    let mut seqs: Vec<Vec<u8>> = Vec::new();

    for (_, (chr, interval)) in chrs.iter().zip(intervals.iter()).enumerate() {
        let mut seq: Vec<u8> = Vec::new();

        faidx
            .fetch(chr, interval.start, interval.end)
            .expect("Couldn't fetch interval");

        faidx
            .read(&mut seq)
            .expect("Couldn't read the interval");

        // Normalize sequences (replace IUPAC bases other than {A,C,G,T} with 'N's) and split on 'N's.
        for split_seq in seq.normalize(false).split(|ch| *ch == b'N') {
            if split_seq.len() > 0 {
                seqs.push(split_seq.to_vec());
            }
        }
    }

    seqs
}

/// Parse optional locus specifications provided on the command line.
fn get_intervals(locus: &Option<Vec<String>>, faidx: &IndexedReader<std::fs::File>) -> (Vec<String>, Vec<Interval<u64>>) {
    let mut chrs = Vec::new();
    let mut intervals = Vec::new();

    if locus.as_ref().unwrap_or(&Vec::new()).is_empty() {
        // If no loci are specified, use the entire FASTA file.
        for seq in faidx.index.sequences() {
            let interval = Interval::new(0..seq.len).unwrap();

            chrs.push(seq.name);
            intervals.push(interval);
        }
    } else {
        // If loci are specified, parse the provided strings (removing any formatting) and store them.
        for l in locus.as_ref().unwrap() {
            let sep = Regex::new(r"[:-]").unwrap();
            let m = l.replace(',', "");
            let pieces: Vec<&str> = sep.split(&m).collect();

            let interval = Interval::new(
                pieces[1].parse::<u64>().expect("Could not parse start range of interval.")
                ..
                pieces[2].parse::<u64>().expect("Could not parse second range of interval.")
            ).expect("Could not construct Interval object for user-specified interval.");

            chrs.push(String::from(pieces[0]));
            intervals.push(interval);
        };
    }

    // Do not allow processing of very long (> 100 Mbp) regions.
    // This is a sanity check to prevent someone from accidentally specifying the entire
    // human genome to process, as this tool is intended only for assembling targeted regions.
    check_interval_length_limit(&intervals);

    (chrs, intervals)
}

/// Check that specified intervals do not exceed TOTAL_INTERVAL_LENGTH_LIMIT in total length.
fn check_interval_length_limit(intervals: &[Interval<u64>]) {
    let total_interval_length = sum_interval_lengths(intervals);
    if total_interval_length > TOTAL_INTERVAL_LENGTH_LIMIT {
        panic!("Total interval length is too large ({} bp); must be less than {} bp.", total_interval_length, TOTAL_INTERVAL_LENGTH_LIMIT);
    }
}

/// Sum the length of user-specified intervals.
fn sum_interval_lengths(intervals: &[Interval<u64>]) -> u64 {
    let mut sum: u64 = 0;

    for interval in intervals {
        sum += interval.end - interval.start + 1;
    }

    sum
}