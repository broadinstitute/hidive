use num_format::{Locale, ToFormattedString};
use std::path::PathBuf;


use crate::train::{create_dataset_for_model, distance_to_a_contig_end, plot_roc_curve, process_reads};
use gbdt::decision_tree::DataVec;
use gbdt::gradient_boost::GBDT;
use skydive::ldbg::LdBG;
use std::io::Write;

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
    let (precision, recall, f1_score) = crate::train::compute_precision_recall_f1(&test_data, &prediction, pred_threshold);
    skydive::elog!("Prediction threshold: {:.2}", pred_threshold);
    skydive::elog!("Precision: {:.2}%", 100.0 * precision);
    skydive::elog!("Recall: {:.2}%", 100.0 * recall);
    skydive::elog!("F1 score: {:.2}%", 100.0 * f1_score);

    // TPR and FPR calculations at various thresholds.
    let fpr_tpr = crate::train::compute_fpr_tpr(&test_data, &prediction);

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


