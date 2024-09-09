use num_format::{Locale, ToFormattedString};
use rand::SeedableRng;
use std::path::PathBuf;

use gbdt::config::{loss2string, Config, Loss};
use gbdt::decision_tree::{Data, DataVec};
use gbdt::gradient_boost::GBDT;

// Import the skydive module, which contains the necessary functions for staging data
use skydive::ldbg::LdBG;

pub fn start(
    output: &PathBuf,
    kmer_size: usize,
    iterations: usize,
    test_split: f32,
    long_read_seq_paths: &Vec<PathBuf>,
    short_read_seq_paths: &Vec<PathBuf>,
    truth_seq_paths: &Vec<PathBuf>,
    debug: bool,
) {
    let long_read_seq_urls = skydive::parse::parse_file_names(&long_read_seq_paths);
    let short_read_seq_urls = skydive::parse::parse_file_names(&short_read_seq_paths);
    let truth_seq_urls = skydive::parse::parse_file_names(&truth_seq_paths);

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

    let l1 = LdBG::from_sequences(String::from("l1"), kmer_size, &all_lr_seqs);

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

    let s1 = LdBG::from_sequences(String::from("s1"), kmer_size, &all_sr_seqs);

    // Read all truth sequences.
    let mut all_truth_seqs: Vec<Vec<u8>> = Vec::new();
    for truth_seq_url in &truth_seq_urls {
        let basename = skydive::utils::basename_without_extension(
            &truth_seq_url,
            &[".fasta.gz", ".fa.gz", ".fasta", ".fa"],
        );
        let fasta_path = truth_seq_url.to_file_path().unwrap();

        skydive::elog!("Processing truth sample {}...", basename);

        let reader = bio::io::fasta::Reader::from_file(&fasta_path).unwrap();
        let all_reads: Vec<bio::io::fasta::Record> = reader.records().map(|r| r.unwrap()).collect();

        all_truth_seqs.extend(
            all_reads
                .iter()
                .map(|r| r.seq().to_vec())
                .collect::<Vec<Vec<u8>>>(),
        );
    }

    let t1 = LdBG::from_sequences(String::from("s1"), kmer_size, &all_truth_seqs);

    // Configure GBDT.
    let mut cfg = Config::new();
    cfg.set_feature_size(3);
    cfg.set_max_depth(4);
    cfg.set_min_leaf_size(1);
    cfg.set_loss(&loss2string(&Loss::BinaryLogitraw));
    cfg.set_shrinkage(0.1);
    cfg.set_iterations(iterations);
    cfg.set_debug(debug);

    skydive::elog!("Training GBDT model with:");
    skydive::elog!(" - feature_size={}", cfg.feature_size);
    skydive::elog!(" - max_depth={}", cfg.max_depth);
    skydive::elog!(" - min_leaf_size={}", cfg.min_leaf_size);
    skydive::elog!(" - loss={}", loss2string(&cfg.loss));
    skydive::elog!(" - shrinkage={}", cfg.shrinkage);
    skydive::elog!(" - iterations={}", cfg.iterations);
    skydive::elog!(" - train/test={}/{}", 1.0 - test_split, test_split);

    // Initialize GBDT algorithm.
    let mut gbdt = GBDT::new(&cfg);

    let mut training_data: DataVec = Vec::new();
    let mut test_data: DataVec = Vec::new();

    let mut rng = rand::rngs::StdRng::seed_from_u64(0);

    let kmers = l1
        .kmers
        .keys()
        .chain(s1.kmers.keys())
        .chain(t1.kmers.keys());
    for kmer in kmers {
        let compressed_len = skydive::utils::homopolymer_compressed(kmer).len();

        let lcov = l1.kmers.get(kmer).map_or(0, |lr| lr.coverage());
        let scov = s1.kmers.get(kmer).map_or(0, |sr| sr.coverage());
        let tcov = t1.kmers.get(kmer).map_or(0, |tr| tr.coverage());

        let data = Data::new_training_data(
            vec![
                lcov as f32,
                scov as f32,
                (kmer.len() - compressed_len) as f32,
            ],
            1.0,
            if tcov > 0 { 1.0 } else { 0.0 },
            None,
        );

        // Save 20% of the data for testing, remaining for training.
        if rand::Rng::gen::<f32>(&mut rng) < test_split {
            test_data.push(data);
        } else {
            training_data.push(data);
        }
    }

    // Train the decision trees.
    skydive::elog!(
        "Training GBDT model with {} training points...",
        training_data.len().to_formatted_string(&Locale::en)
    );
    gbdt.fit(&mut training_data);

    // Predict the test data.
    skydive::elog!("Computing accuracy on test data...");
    let p = gbdt.predict(&test_data);

    // Evaluate accuracy of the model on the test data.
    let mut num_correct = 0u32;
    let mut num_total = 0u32;
    for (data, pred) in test_data.iter().zip(p.iter()) {
        let truth = data.label;
        let call = if *pred > 0.5 { 1.0 } else { 0.0 };
        if (call - truth).abs() < f32::EPSILON {
            num_correct += 1;
        }
        num_total += 1;
    }

    skydive::elog!(
        "Predicted accuracy: {}/{} ({:.2}%)",
        num_correct.to_formatted_string(&Locale::en),
        num_total.to_formatted_string(&Locale::en),
        100.0 * num_correct as f32 / num_total as f32
    );

    gbdt.save_model(output.to_str().unwrap())
        .expect("Unable to save model");

    skydive::elog!("Model saved to {}", output.to_str().unwrap());
}
