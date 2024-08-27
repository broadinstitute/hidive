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
use skydive::mldbg::MLdBG;

pub fn start(
    output: &PathBuf,
    kmer_size: usize,
    model_path: &PathBuf,
    long_read_fasta_paths: &Vec<PathBuf>,
    short_read_fasta_paths: &Vec<PathBuf>,
) {
    let long_read_seq_urls = skydive::parse::parse_file_names(long_read_fasta_paths);
    let short_read_seq_urls = skydive::parse::parse_file_names(short_read_fasta_paths);

    let m = MLdBG::new(kmer_size);

    // Read all long reads.
    skydive::elog!("Processing long-read samples {:?}...", long_read_seq_urls);
    // let l1 = LdBG::from_files(String::from("l1"), kmer_size, &long_read_fasta_paths);
    

    // Read all short reads.
    skydive::elog!("Processing short-read samples {:?}...", short_read_seq_urls);
    let s1 = LdBG::from_files(String::from("s1"), kmer_size, &short_read_fasta_paths);

    /*
    // Union of kmers from l1 and s1.
    let kmers: HashSet<_> = l1.kmers.keys().chain(s1.kmers.keys()).cloned().collect();
    for kmer in kmers {
        let lr = l1.kmers.get(&kmer);
        let sr = s1.kmers.get(&kmer);

        if lr.is_none() && sr.is_some() {
            s1.remove(&kmer);
        }
    }

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
    let mut l3 = LdBG::from_sequences(String::from("l3"), kmer_size, &all_seqs);

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
    let (cleaned_tips_kmers, cleaned_tips_paths) = l3.clean_tips(2 * kmer_size);

    skydive::elog!(
        "K-mers with p < 0.5: {} / {} ({:.2}%)",
        num_below_threshold,
        num_total,
        num_below_threshold as f32 / num_total as f32 * 100.0
    );
    skydive::elog!(
        "Removed {} k-mers in {} paths",
        cleaned_kmers,
        cleaned_paths
    );
    skydive::elog!(
        "Removed {} k-mers in {} tips",
        cleaned_tips_kmers,
        cleaned_tips_paths
    );

    let mut all_corrected_seqs = Vec::new();
    for lr_seq in all_seqs {
        let corrected_seqs = l3.correct_seq(lr_seq);

        all_corrected_seqs.extend(corrected_seqs);
    }

    skydive::elog!("Corrected sequences: {}", all_corrected_seqs.len());

    // Assemble contigs.
    l3.links = LdBG::build_links(kmer_size, &all_corrected_seqs, &l3.kmers);
    let contigs = l3.assemble_all();

    // Write contigs to disk.
    skydive::elog!("Writing {} contigs to disk...", &contigs.len());

    let output_file = File::create(output).unwrap();
    let mut writer = BufWriter::new(output_file);

    for (i, contig) in contigs.iter().enumerate() {
        writeln!(
            writer,
            ">contig_{}\n{}",
            i,
            String::from_utf8(contig.clone()).unwrap()
        )
        .unwrap();
    }
    */
}
