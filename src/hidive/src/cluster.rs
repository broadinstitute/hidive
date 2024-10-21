// Import necessary standard library modules
use std::collections::HashSet;
use std::io::Write;
use std::path::PathBuf;

use bio::io::fasta::{Reader, Record};
use needletail::Sequence;

use gaoya::minhash::{MinHashIndex, MinHasher, MinHasher32};

// Import the skydive module, which contains the necessary functions for building graphs
// use skydive;

pub fn start(output: &PathBuf, k: usize, jaccard_threshold: f64, fasta_path: &PathBuf) {
    // Read the reads records (name and sequence) into a vector.
    let reader = Reader::from_file(fasta_path).unwrap();
    let all_reads: Vec<Record> = reader.records().map(|r| r.unwrap()).collect();

    // Get and store sample name.
    let mut sample_names = Vec::new();
    for record in &all_reads {
        let id = record.id().to_string();
        let v: Vec<String> = id.split('|').map(|s| s.to_string()).collect();

        // sample_names.push(v[1].clone());
        sample_names.push(id);
    }

    // Iterate over all reads.
    let sequences: Vec<String> = all_reads
        .iter()
        .map(|record| String::from_utf8_lossy(record.seq()).to_string().to_uppercase())
        .collect();

    let (num_bands, band_width) = (42, 3);
    let minhasher = MinHasher32::new(num_bands * band_width);
    let mut index = MinHashIndex::new(num_bands, band_width, jaccard_threshold);
    for (i, doc) in sequences.iter().enumerate() {
        index.insert(i, minhasher.create_signature(doc.as_bytes().kmers(k as u8)));
    }

    // Open a new file to write the k-mer matrix to.
    let file = std::fs::File::create(output).unwrap();
    let mut writer = std::io::BufWriter::new(&file);

    let mut used = HashSet::new();
    for (i, sequence) in sequences.iter().enumerate() {
        if used.contains(&i) {
            continue;
        }

        let indices = index.query_owned(&minhasher.create_signature(sequence.as_bytes().kmers(k as u8)));

        // println!("{} {:?}", sample_names[i], a);
        for index in indices {
            // println!(">cluster_{}_sequence_{}\n{}", i, index, sequences[index]);
            let _ = writeln!(writer, ">cluster{}_sequence{}\n{}", i, index, sequences[index]);

            used.insert(index);
        }
    }
}
