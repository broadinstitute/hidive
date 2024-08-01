// Import necessary standard library modules
use std::collections::HashMap;
use std::io::Write;
use std::path::PathBuf;

use bio::io::fasta::{Reader, Record};
// Import the skydive module, which contains the necessary functions for building graphs
// use skydive;

pub fn start(output: &PathBuf, k: usize, fasta_path: &PathBuf) {
    // Read the reads records (name and sequence) into a vector.
    let reader = Reader::from_file(fasta_path).unwrap();
    let all_reads: Vec<Record> = reader.records().map(|r| r.unwrap()).collect();

    // Get and store sample name.
    let mut sample_names = HashMap::new();
    let mut sample_index = 0;
    for record in &all_reads {
        let id = record.id().to_string();
        let v: Vec<String> = id.split('|').map(|s| s.to_string()).collect();

        if !sample_names.contains_key(&v[1]) {
            sample_names.insert(v[1].clone(), sample_index);
            sample_index += 1;
        }
    }

    // Create a HashMap to store the k-mer presence/absence matrix.
    let mut kmer_matrix: HashMap<String, Vec<bool>> = HashMap::new();

    // Iterate over all reads.
    for record in &all_reads {
        // Get the sample index of the current read.
        let id = record.id().to_string();
        let v: Vec<String> = id.split('|').map(|s| s.to_string()).collect();
        let sample_index = *sample_names.get(&v[1]).unwrap();

        // Get the sequence of the current read.
        let sequence = String::from_utf8_lossy(record.seq()).to_string();

        // Iterate over the sequence with a sliding window of size k.
        for i in 0..sequence.len() - k + 1 {
            // Get the current k-mer.
            let kmer = &sequence[i..i + k];

            // Check if the k-mer is already in the matrix.
            if let Some(presence_vector) = kmer_matrix.get_mut(kmer) {
                // If the k-mer is already in the matrix, add a presence (true) for the current read.
                presence_vector[sample_index] = true;
            } else {
                // If the k-mer is not in the matrix, create a new presence vector for it.
                let mut presence_vector = vec![false; sample_names.len()];
                presence_vector[sample_index] = true;
                kmer_matrix.insert(kmer.to_string(), presence_vector);
            }
        }
    }

    // Open a new file to write the k-mer matrix to.
    let file = std::fs::File::create(output).unwrap();
    let mut writer = std::io::BufWriter::new(&file);

    // Write the header row (sample names).
    let header = sample_names.keys().cloned().collect::<Vec<_>>().join(",");
    writer.write_all(b"kmer,").unwrap();
    writer.write_all(header.as_bytes()).unwrap();
    writer.write_all(b"\n").unwrap();

    // Write the presence/absence matrix.
    for (kmer, presence_vector) in &kmer_matrix {
        let row_str = presence_vector
            .iter()
            .map(|b| if *b { "1" } else { "0" })
            .collect::<Vec<&str>>()
            .join(",");
        writer.write_all(kmer.as_bytes()).unwrap();
        writer.write_all(b",").unwrap();
        writer.write_all(row_str.as_bytes()).unwrap();
        writer.write_all(b"\n").unwrap();
    }
}
