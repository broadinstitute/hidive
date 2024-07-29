use std::collections::HashSet;
use std::path::PathBuf;
use std::{fs::File, io::{BufWriter, Write}};

use skydive::ldbg::LdBG;
use skydive::mldbg::MLdBG;

pub fn start(
    output: &PathBuf,
    kmer_size: usize,
    long_read_fasta_paths: &Vec<PathBuf>,
    short_read_fasta_paths: &Vec<PathBuf>,
) {
    let long_read_seq_urls = skydive::parse::parse_file_names(long_read_fasta_paths);
    let short_read_seq_urls = skydive::parse::parse_file_names(short_read_fasta_paths);

    // Read all long reads.
    let mut all_lr_seqs: Vec<Vec<u8>> = Vec::new();
    for long_read_seq_url in &long_read_seq_urls {
        let basename = skydive::utils::basename_without_extension(&long_read_seq_url, &[".fasta.gz", ".fa.gz", ".fasta", ".fa"]);
        let fasta_path = long_read_seq_url.to_file_path().unwrap();

        skydive::elog!("Processing long-read sample {}...", basename);

        let reader = bio::io::fasta::Reader::from_file(&fasta_path).unwrap();
        let all_reads: Vec<bio::io::fasta::Record> = reader.records().map(|r| r.unwrap()).collect();

        all_lr_seqs.extend(all_reads
            .iter()
            .map(|r| r.seq().to_vec())
            .collect::<Vec<Vec<u8>>>());
    }

    let mut l1 = LdBG::from_sequences(String::from("l1"), kmer_size, &all_lr_seqs, false, false);

    // Read all short reads.
    let mut all_sr_seqs: Vec<Vec<u8>> = Vec::new();
    for short_read_seq_url in &short_read_seq_urls {
        let basename = skydive::utils::basename_without_extension(&short_read_seq_url, &[".fasta.gz", ".fa.gz", ".fasta", ".fa"]);
        let fasta_path = short_read_seq_url.to_file_path().unwrap();

        skydive::elog!("Processing short-read sample {}...", basename);

        let reader = bio::io::fasta::Reader::from_file(&fasta_path).unwrap();
        let all_reads: Vec<bio::io::fasta::Record> = reader.records().map(|r| r.unwrap()).collect();

        all_sr_seqs.extend(all_reads
            .iter()
            .map(|r| r.seq().to_vec())
            .collect::<Vec<Vec<u8>>>());
    }

    let mut s1 = LdBG::from_sequences(String::from("s1"), kmer_size, &all_sr_seqs, false, false);

    println!("l1 len before {}", &l1.kmers.len());
    println!("s1 len before {}", &s1.kmers.len());

    // Union of kmers from l1 and s1.
    let kmers: HashSet<_> = l1.kmers.keys().chain(s1.kmers.keys()).cloned().collect();

    for kmer in kmers {
        let lr = l1.kmers.get(&kmer);
        let sr = s1.kmers.get(&kmer);

        if lr.is_none() && sr.is_some() {
            s1.remove(&kmer);
        }
    }

    s1.infer_edges();

    println!("s1 len after {}", &s1.kmers.len());

    let filtered_sr_seqs = &all_sr_seqs
        .iter()
        .cloned()
        .filter(|read| {
            let num_contained = read.windows(kmer_size)
                .filter(|kmer| s1.kmers.contains_key(&skydive::ldbg::LdBG::canonicalize_kmer(kmer)))
                .count();

            // println!("{} {}", num_contained, read.len() - kmer_size + 1);

            // num_contained > read.len() / 4
            num_contained == read.len() - kmer_size + 1
        })
        .collect::<Vec<Vec<u8>>>();

    println!("{} {}", &all_sr_seqs.len(), &filtered_sr_seqs.len());

    // Print contigs.
    let contigs = s1.assemble_all();
    let output_file = File::create(output).unwrap();
    let mut writer = BufWriter::new(output_file);

    for (i, read) in filtered_sr_seqs.iter().enumerate() {
        writeln!(writer, ">sr_read_{}\n{}", i, String::from_utf8(read.clone()).unwrap());
    }

    // for (i, contig) in contigs.iter().enumerate() {
    //     writeln!(writer, ">sr_unitig_{}\n{}", i, String::from_utf8(contig.clone()).unwrap()).unwrap();
    // }

    /*
    // Create a linked de Bruijn graph
    // let lr_path = &long_read_seq_urls.iter().next().unwrap().to_file_path().unwrap();
    // let l1 = LdBG::from_file(String::from("l1"), kmer_size, &lr_path, false, true);
    let l1 = LdBG::from_sequences(String::from("l1"), kmer_size, &all_fwd_seqs, false, true);

    println!("number of reads {}", all_fwd_seqs.len());
    println!("l1 len {}", l1.kmers.len());

    // Create a multi-color linked de Bruijn graph
    let mut g = MLdBG::new(kmer_size, false, false);
    g.ldbgs.push(l1);

    // Construct a graph from the short read data, filtered for k-mer overlap with the existing graph, and add it to g as a separate color.
    for short_read_seq_url in short_read_seq_urls {
        let basename = skydive::utils::basename_without_extension(&short_read_seq_url, &[".fasta.gz", ".fa.gz", ".fasta", ".fa"]);
        let fasta_path = short_read_seq_url.to_file_path().unwrap();

        let reader = bio::io::fasta::Reader::from_file(&fasta_path).unwrap();
        let all_reads: Vec<bio::io::fasta::Record> = reader.records().map(|r| r.unwrap()).collect();
        println!("Number of short reads before filtering {}", all_reads.len());

        skydive::elog!("Filtering short-read sample {}...", basename);
        all_fwd_seqs.extend(
            g.filter_reads(&fasta_path, |r, kmer_union| {
            let num_contained = r.seq().windows(kmer_size)
                .filter(|kmer| kmer_union.contains(&skydive::ldbg::LdBG::canonicalize_kmer(kmer)))
                .count();

            num_contained > r.seq().len() / 4
        }));

        println!("Number of total reads after filtering {}", all_fwd_seqs.len());
    }

    for read in &all_fwd_seqs {
        println!("{}", read.len());
    }

    // Print more stats.
    // println!("Num reads: {}", all_fwd_seqs.len());

    // Create a new linked de Bruijn graph
    let l2 = LdBG::from_sequences(String::from("l2"), kmer_size, &all_fwd_seqs, true, true);

    // Print more stats.
    // println!("Num colors: {}", g.len());
    // for l in g.iter() {
    //     println!(" -- {}: {}", l.name(), l.kmers.len());
    // }

    // Print contigs.
    let contigs = l2.assemble_all();
    // let contigs = g.get(0).unwrap().assemble_all();
    let output_file = File::create(output).unwrap();
    let mut writer = BufWriter::new(output_file);

    for (i, contig) in contigs.iter().enumerate() {
        writeln!(writer, ">unitig_{}\n{}", i, String::from_utf8(contig.clone()).unwrap()).unwrap();
    }
    */
}
