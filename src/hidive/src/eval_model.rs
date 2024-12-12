use std::collections::hash_map::Keys;
use std::collections::HashMap;
use num_format::{Locale, ToFormattedString};
use std::path::PathBuf;


use crate::train::{distance_to_a_contig_end, plot_roc_curve, process_reads};
use gbdt::gradient_boost::GBDT;
use gbdt::decision_tree::{Data, DataVec};
use skydive::ldbg::LdBG;
use std::io::Write;
use std::iter::Chain;

use gbdt::config::{loss2string, Config, Loss};


use plotters::prelude::*;
use skydive::record::Record;
use url::Url;
use skydive::nn_model::{KmerData, KmerDataVec};

//////////////////////////////
// Todo: Use NN model instead of GBDT
//
/////////////////////////////

pub fn start(
    output: &PathBuf,
    kmer_size: usize,
    long_read_seq_paths: &Vec<PathBuf>,
    short_read_seq_paths: &Vec<PathBuf>,
    truth_seq_paths: &Vec<PathBuf>,
    model_path: &PathBuf,
) {
    let long_read_seq_urls = skydive::parse::parse_file_names(long_read_seq_paths);
    let short_read_seq_urls = skydive::parse::parse_file_names(short_read_seq_paths);
    let truth_seq_urls = skydive::parse::parse_file_names(truth_seq_paths);


    // Read all long reads.
    let all_lr_seqs: Vec<Vec<u8>> = process_reads(&long_read_seq_urls, "long");
    let l1 = LdBG::from_sequences(String::from("l1"), kmer_size, &all_lr_seqs);

    // Read all short reads.
    let all_sr_seqs: Vec<Vec<u8>> = process_reads(&short_read_seq_urls, "short");
    let s1 = LdBG::from_sequences(String::from("s1"), kmer_size, &all_sr_seqs);

    // Read all truth sequences.
    let all_truth_seqs: Vec<Vec<u8>> = process_reads(&truth_seq_urls, "truth");
    let t1 = LdBG::from_sequences(String::from("s1"), kmer_size, &all_truth_seqs);

    let lr_contigs = l1.assemble_all();
    let lr_distances = distance_to_a_contig_end(&lr_contigs, kmer_size);

    let sr_contigs = s1.assemble_all();
    let sr_distances = distance_to_a_contig_end(&sr_contigs, kmer_size);


    // load model
    skydive::elog!("Loading GBDT model from {}...", model_path.to_str().unwrap());
    let gbdt = GBDT::load_model(model_path.to_str().unwrap()).expect("Unable to load model");

    // Prepare test data.
    let test_kmers = l1
        .kmers
        .keys()
        .chain(s1.kmers.keys())
        .chain(t1.kmers.keys());
    let test_data: DataVec = create_dataset_for_model(
        test_kmers,
        &lr_distances,
        &sr_distances,
        &l1,
        &s1,
        &t1,
    );

    // Predict the test data.
    skydive::elog!("Computing accuracy on test data...");
    let prediction = gbdt.predict(&test_data);
    let pred_threshold = 0.5;

    // Evaluate accuracy of the model on the test data.
    let mut num_correct = 0u32;
    let mut num_total = 0u32;
    for (data, pred) in test_data.iter().zip(prediction.iter()) {
        let truth = data.label;
        let call = if *pred > pred_threshold { 1.0 } else { 0.0 };
        if (call - truth).abs() < f32::EPSILON {
            num_correct += 1;
        }
        num_total += 1;
    }

    // Precision, Recall, and F1 score calculations.
    let (precision, recall, f1_score) = compute_precision_recall_f1(&test_data, &prediction, pred_threshold);
    skydive::elog!("Prediction threshold: {:.2}", pred_threshold);
    skydive::elog!("Precision: {:.2}%", 100.0 * precision);
    skydive::elog!("Recall: {:.2}%", 100.0 * recall);
    skydive::elog!("F1 score: {:.2}%", 100.0 * f1_score);

    // TPR and FPR calculations at various thresholds.
    let fpr_tpr = compute_fpr_tpr(&test_data, &prediction);

    // Save TPR and FPR at various thresholds to a file.
    let csv_output = output.with_extension("csv");
    let mut writer = std::fs::File::create(&csv_output).expect("Unable to create file");
    for (tpr, fpr) in &fpr_tpr {
        writeln!(writer, "{},{}", tpr, fpr).expect("Unable to write data");
    }
    skydive::elog!("TPR and FPR at various thresholds saved to {}", csv_output.to_str().unwrap());

    // Create a ROC curve.
    let png_output = output.with_extension("png");
    plot_roc_curve(&png_output, &fpr_tpr).expect("Unable to plot ROC curve");
    skydive::elog!("ROC curve saved to {}", png_output.to_str().unwrap());

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

/// Creates a dataset for the ML model.
pub fn create_dataset_for_model(
    kmers: Chain<Chain<Keys<Vec<u8>, Record>, Keys<Vec<u8>, Record>>, Keys<Vec<u8>, Record>>,
    lr_distances: &HashMap<Vec<u8>, usize>,
    sr_distances: &HashMap<Vec<u8>, usize>,
    l1: &LdBG,
    s1: &LdBG,
    t1: &LdBG,
) -> Vec<Data> {
    let mut dataset = Vec::new();

    for kmer in kmers {
        let compressed_len = skydive::utils::homopolymer_compressed(kmer).len();
        let compressed_len_diff = (kmer.len() - compressed_len) as f32;
        let entropy = skydive::utils::shannon_entropy(kmer);
        let gc_content = skydive::utils::gc_content(kmer);
        // let lr_distance = *lr_distances.get(kmer).unwrap_or(&0) as f32;
        // let sr_distance = *sr_distances.get(kmer).unwrap_or(&0) as f32;

        let lcov = l1.kmers.get(kmer).map_or(0, |lr| lr.coverage());

        let scov_fw = s1.kmers.get(kmer).map_or(0, |sr| sr.fw_coverage());
        let scov_rc = s1.kmers.get(kmer).map_or(0, |sr| sr.rc_coverage());
        let scov_total = scov_fw + scov_rc;
        let strand_ratio = if scov_total > 0 {
            (scov_fw as f32).max(scov_rc as f32) / scov_total as f32
        } else {
            0.5
        };

        // let expected = (scov_fw + scov_rc) as f32 / 2.0;
        // let chi_square = if expected > 0.0 {
        //     ((scov_fw as f32 - expected).powi(2) +
        //     (scov_rc as f32 - expected).powi(2)) / expected
        // } else {
        //     0.0
        // };

        let tcov = t1.kmers.get(kmer).map_or(0, |tr| tr.coverage());

        let data = Data::new_training_data(
            vec![
                if lcov > 0 { 1.0 } else { 0.0 },    // present in long reads
                scov_total as f32,                   // coverage in short reads
                strand_ratio as f32,                 // measure of strand bias (0.5 = balanced, 1.0 = all on one strand)
                compressed_len_diff,                 // homopolymer compression length difference
                entropy,                             // shannon entropy
                gc_content,                          // gc content
                // lr_distance,                         // distance to nearest long read contig end
                // sr_distance,                         // distance to nearest short read contig end
            ],
            1.0,
            if tcov > 0 { 1.0 } else { 0.0 }, // present in truth
            None,
        );

        dataset.push(data);
    }

    dataset
}

pub fn compute_precision_recall_f1(test_data: &DataVec, p: &Vec<f32>, pred_threshold: f32) -> (f32, f32, f32) {
    skydive::elog!("Computing precision, recall, and F1 score on test data...");
    let (mut num_true_positives, mut num_false_positives, mut num_false_negatives) = (0u32, 0u32, 0u32);

    for (data, pred) in test_data.iter().zip(p.iter()) {
        let truth = data.label;
        let call = if *pred > pred_threshold { 1.0 } else { 0.0 };
        if (call - truth).abs() < f32::EPSILON {
            if truth > 0.5 { num_true_positives += 1; }
        } else {
            if truth < 0.5 { num_false_positives += 1; }
            if truth > 0.5 { num_false_negatives += 1; }
        }
    }

    let precision = num_true_positives as f32 / (num_true_positives + num_false_positives) as f32;
    let recall = num_true_positives as f32 / (num_true_positives + num_false_negatives) as f32;
    let f1_score = 2.0 * precision * recall / (precision + recall);

    (precision, recall, f1_score)
}


/// Computes FPR and TPR at various thresholds.
pub fn compute_fpr_tpr(test_data: &DataVec, p: &Vec<f32>) -> Vec<(f32, f32)> {
    skydive::elog!("Computing FPR and TPR at various thresholds...");
    let mut fpr_tpr: Vec<(f32, f32)> = Vec::new();

    // Iterate over different thresholds
    for threshold in 0..100 {
        let pred_threshold = threshold as f32 / 100.0;
        let (mut num_true_positives, mut num_false_positives, mut num_false_negatives, mut num_true_negatives) = (0u32, 0u32, 0u32, 0u32);

        // Iterate over test data and predictions
        for (data, pred) in test_data.into_iter().zip(p.iter()) {
            let truth = data.label;
            let call = if *pred > pred_threshold { 1.0 } else { 0.0 };

            if (call - truth).abs() < f32::EPSILON {
                // Correct classification
                if truth > 0.5 { num_true_positives += 1; }  // True Positive
                else { num_true_negatives += 1; }  // True Negative
            } else {
                // Incorrect classification
                if truth < 0.5 { num_false_positives += 1; }  // False Positive
                if truth > 0.5 { num_false_negatives += 1; }  // False Negative
            }
        }

        // Calculate FPR (False Positive Rate)
        let fpr = if num_false_positives + num_true_negatives > 0 {
            num_false_positives as f32 / (num_false_positives + num_true_negatives) as f32
        } else {
            0.0
        };

        // Calculate TPR (True Positive Rate)
        let tpr = if num_true_positives + num_false_negatives > 0 {
            num_true_positives as f32 / (num_true_positives + num_false_negatives) as f32
        } else {
            0.0
        };

        fpr_tpr.push((fpr, tpr));
    }

    fpr_tpr
}

