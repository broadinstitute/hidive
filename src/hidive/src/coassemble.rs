use std::collections::HashSet;
use std::path::PathBuf;
use std::{
    fs::File,
    io::{BufWriter, Write},
};

use gbdt::config::{loss2string, Config, Loss};
use gbdt::decision_tree::{Data, DataVec};
use gbdt::gradient_boost::GBDT;

use skydive::ldbg::LdBG;

pub fn start(
    output: &PathBuf,
    kmer_size: usize,
    model_path: &PathBuf,
    long_read_fasta_paths: &Vec<PathBuf>,
    short_read_fasta_paths: &Vec<PathBuf>,
) {
    let long_read_seq_urls = skydive::parse::parse_file_names(long_read_fasta_paths);
    let short_read_seq_urls = skydive::parse::parse_file_names(short_read_fasta_paths);

    // Read all long reads.
    let mut all_lr_seqs: Vec<Vec<u8>> = Vec::new();
    for long_read_seq_url in &long_read_seq_urls {
        let basename = skydive::utils::basename_without_extension(
            &long_read_seq_url,
            &[".fasta.gz", ".fa.gz", ".fasta", ".fa"],
        );
        let fasta_path = long_read_seq_url.to_file_path().unwrap();

        skydive::elog!("Processing long-read sample {}...", basename);

        let reader = bio::io::fasta::Reader::from_file(&fasta_path).unwrap();
        let all_reads: Vec<bio::io::fasta::Record> = reader.records().map(|r| r.unwrap()).collect();

        all_lr_seqs.extend(
            all_reads
                .iter()
                .map(|r| r.seq().to_vec())
                .collect::<Vec<Vec<u8>>>(),
        );
    }

    let l1 = LdBG::from_sequences(String::from("l1"), kmer_size, &all_lr_seqs, false, false);

    // Read all short reads.
    let mut all_sr_seqs: Vec<Vec<u8>> = Vec::new();
    for short_read_seq_url in &short_read_seq_urls {
        let basename = skydive::utils::basename_without_extension(
            &short_read_seq_url,
            &[".fasta.gz", ".fa.gz", ".fasta", ".fa"],
        );
        let fasta_path = short_read_seq_url.to_file_path().unwrap();

        skydive::elog!("Processing short-read sample {}...", basename);

        let reader = bio::io::fasta::Reader::from_file(&fasta_path).unwrap();
        let all_reads: Vec<bio::io::fasta::Record> = reader.records().map(|r| r.unwrap()).collect();

        all_sr_seqs.extend(
            all_reads
                .iter()
                .map(|r| r.seq().to_vec())
                .collect::<Vec<Vec<u8>>>(),
        );
    }

    let mut s1 = LdBG::from_sequences(String::from("s1"), kmer_size, &all_sr_seqs, false, false);

    // Union of kmers from l1 and s1.
    let kmers: HashSet<_> = l1.kmers.keys().chain(s1.kmers.keys()).cloned().collect();
    for kmer in kmers {
        let lr = l1.kmers.get(&kmer);
        let sr = s1.kmers.get(&kmer);

        if lr.is_none() && sr.is_some() {
            s1.remove(&kmer);
        }
    }

    // s1.infer_edges();

    // Filter out the irrelevant short reads.
    let filtered_sr_seqs = &all_sr_seqs
        .iter()
        .cloned()
        .filter(|read| {
            let num_contained = read
                .windows(kmer_size)
                .filter(|kmer| {
                    s1.kmers
                        .contains_key(&skydive::ldbg::LdBG::canonicalize_kmer(kmer))
                })
                .count();

            num_contained == read.len() - kmer_size + 1
        })
        .collect::<Vec<Vec<u8>>>();

    // Report some stats.
    skydive::elog!("Long reads: {}", &all_lr_seqs.len());
    skydive::elog!("Short reads (before filtering): {}", &all_sr_seqs.len());
    skydive::elog!("Short reads (after filtering): {}", &filtered_sr_seqs.len());

    let all_seqs = all_lr_seqs
        .into_iter()
        .chain(filtered_sr_seqs.iter().cloned())
        .collect::<Vec<Vec<u8>>>();
    skydive::elog!("Combined reads: {}", &all_seqs.len());

    // Assemble contigs.
    let mut l3 = LdBG::from_sequences(String::from("l3"), kmer_size, &all_seqs, true, true);

    // Filter k-mers.
    let mut eval_data: DataVec = Vec::new();
    let graph_kmers = &l3.kmers.keys().cloned().collect::<Vec<_>>();
    for kmer in graph_kmers {
        let lcov = l1.kmers.get(kmer).map_or(0, |record| record.coverage());
        let scov = s1.kmers.get(kmer).map_or(0, |record| record.coverage());
        let compressed_len = skydive::utils::run_length_encoded(kmer).len();

        let data = Data::new_training_data(
            vec![
                lcov as f32,
                scov as f32,
                (kmer.len() - compressed_len) as f32,
            ],
            1.0,
            0.0,
            None,
        );

        eval_data.push(data);

        if lcov + scov < 10 {
            // l3.remove(kmer);
        }
    }

    let gbdt = GBDT::load_model(model_path.to_str().unwrap()).unwrap();
    let predictions = gbdt.predict(&eval_data);
    let mut num_below_threshold = 0;
    let mut num_total = 0;

    for (p, data) in predictions.iter().zip(eval_data.iter()) {
        // println!("pred={:.2} data={:?}", p, data);

        if *p < 0.5 {
            num_below_threshold += 1;
            // l3.remove(&graph_kmers[i]);
        }
        num_total += 1;
    }

    skydive::elog!(
        "Number of data points with p < 0.5: {} / {}",
        num_below_threshold,
        num_total
    );

    // l3.infer_edges();
    // l3.links = LdBG::build_links(kmer_size, &all_seqs, &l3.kmers);

    // Write contigs to disk.
    let contigs = l3.assemble_all();
    let output_file = File::create(output).unwrap();
    let mut writer = BufWriter::new(output_file);

    skydive::elog!("Writing {} contigs to disk...", &contigs.len());
    for (i, contig) in contigs.iter().enumerate() {
        writeln!(
            writer,
            ">unitig_{}\n{}",
            i,
            String::from_utf8(contig.clone()).unwrap()
        )
        .unwrap();
    }
}
