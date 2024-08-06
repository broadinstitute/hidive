use std::collections::HashSet;
use std::path::PathBuf;
use std::{
    fs::File,
    io::{BufWriter, Write},
};

use gbdt::config::{loss2string, Config, Loss};
use gbdt::decision_tree::{Data, DataVec};
use gbdt::gradient_boost::GBDT;

use needletail::Sequence;
use parquet::data_type::AsBytes;
use petgraph::dot::Dot;

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
    let mut all_lr_seqs2: Vec<Vec<u8>> = Vec::new();

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

        all_lr_seqs2.extend(
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

    let all_seqs = &all_lr_seqs
        .into_iter()
        .chain(filtered_sr_seqs.iter().cloned())
        .collect::<Vec<Vec<u8>>>();
    skydive::elog!("Combined reads: {}", &all_seqs.len());

    // Assemble contigs.
    let mut l3 = LdBG::from_sequences(String::from("l3"), kmer_size, &all_seqs, false, false);

    // l3.mark_tips(10*l3.kmer_size);

    // let mut num_below_threshold = 0;
    // for cn_kmer in l3.scores.keys() {
    //     if l3.scores.get(cn_kmer).unwrap() < &0.5 {
    //         num_below_threshold += 1;
    //     }
    // }

    // skydive::elog!("{} {}", num_below_threshold, l3.scores.len());

    // Filter k-mers.
    let mut eval_data: DataVec = Vec::new();
    let graph_kmers = &l3.kmers.keys().cloned().collect::<Vec<_>>();
    for kmer in graph_kmers {
        let lcov = l1.kmers.get(kmer).map_or(0, |record| record.coverage());
        let scov = s1.kmers.get(kmer).map_or(0, |record| record.coverage());
        let compressed_len = skydive::utils::homopolymer_compressed(kmer).len();

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
    }

    let gbdt = GBDT::load_model(model_path.to_str().unwrap()).unwrap();
    let predictions = gbdt.predict(&eval_data);
    let mut num_below_threshold = 0;
    let mut num_total = 0;

    for (p, cn_kmer) in predictions.iter().zip(graph_kmers.iter()) {
        l3.scores.insert(cn_kmer.clone(), p.clamp(0.0, 1.0));

        if *p < 0.5 {
            num_below_threshold += 1;
        }
        num_total += 1;
    }

    let (cleaned_kmers, cleaned_paths) = l3.clean_paths(0.4);
    let (cleaned_tips_kmers, cleaned_tips_paths) = l3.clean_tips(2*kmer_size);

    skydive::elog!("K-mers with p < 0.5: {} / {} ({:.2}%)", num_below_threshold, num_total, num_below_threshold as f32 / num_total as f32 * 100.0);
    skydive::elog!("Removed {} k-mers in {} paths", cleaned_kmers, cleaned_paths);
    skydive::elog!("Removed {} k-mers in {} tips", cleaned_tips_kmers, cleaned_tips_paths);

    let mut all_corrected_seqs = Vec::new();
    for lr_seq in all_seqs {
        let corrected_seqs = l3.correct_seq(lr_seq);

        all_corrected_seqs.extend(corrected_seqs);
    }

    skydive::elog!("Corrected sequences: {}", all_corrected_seqs.len());

    // Assemble contigs.
    // l3.links = LdBG::build_links(kmer_size, &all_seqs, &l3.kmers);
    // l3.links = LdBG::build_links(kmer_size, &all_lr_seqs2, &l3.kmers);
    l3.links = LdBG::build_links(kmer_size, &all_corrected_seqs, &l3.kmers);
    let contigs = l3.assemble_all();

    // Write contigs to disk.
    skydive::elog!("Writing {} contigs to disk...", &contigs.len());

    let output_file = File::create(output).unwrap();
    let mut writer = BufWriter::new(output_file);

    for (i, contig) in contigs.iter().enumerate() {
        writeln!(writer, ">contig_{}\n{}", i, String::from_utf8(contig.clone()).unwrap()).unwrap();
    }

    // let graph = l3.traverse_all_kmers();
    // let graph = l3.traverse_all_contigs();

    // let mut gfa_output = Vec::new();
    // skydive::utils::write_graph_as_gfa(&mut gfa_output, &graph).unwrap();
    // let gfa_string = String::from_utf8(gfa_output).unwrap();
    // writeln!(writer, "{}", gfa_string).unwrap();
}
