use num_format::{Locale, ToFormattedString};
// use rand::SeedableRng;
// use skydive::mldbg::MLdBG;
use skydive::utils::canonicalize_kmer;
use std::collections::HashMap;
use std::collections::HashSet;
use std::path::PathBuf;
use std::iter::Chain;
use std::collections::hash_map::Keys;

use gbdt::config::{loss2string, Config, Loss};
use gbdt::decision_tree::{Data, DataVec};
use gbdt::gradient_boost::GBDT;

use plotters::prelude::*;
use skydive::ldbg::LdBG;
use skydive::record::Record;
use std::io::Write;
use url::Url;

pub fn start(
    output: &PathBuf,
    kmer_size: usize,
    iterations: usize,
    // test_split: f32,
    long_read_seq_paths: &Vec<PathBuf>,
    short_read_seq_paths: &Vec<PathBuf>,
    truth_seq_paths: &Vec<PathBuf>,
    test_long_read_seq_paths: &Vec<PathBuf>,
    test_short_read_seq_paths: &Vec<PathBuf>,
    test_truth_seq_paths: &Vec<PathBuf>,
    debug: bool,
) {
    let long_read_seq_urls = skydive::parse::parse_file_names(&long_read_seq_paths);
    let short_read_seq_urls = skydive::parse::parse_file_names(&short_read_seq_paths);
    let truth_seq_urls = skydive::parse::parse_file_names(&truth_seq_paths);
    let test_long_read_seq_urls = skydive::parse::parse_file_names(&test_long_read_seq_paths);
    let test_short_read_seq_urls = skydive::parse::parse_file_names(&test_short_read_seq_paths);
    let test_truth_seq_urls = skydive::parse::parse_file_names(&test_truth_seq_paths);

    // Read all long reads.
    let all_lr_seqs: Vec<Vec<u8>> = process_reads(&long_read_seq_urls, "long");
    let l1 = LdBG::from_sequences(String::from("l1"), kmer_size, &all_lr_seqs);

    // Read all short reads.
    let all_sr_seqs: Vec<Vec<u8>> = process_reads(&short_read_seq_urls, "short");
    let s1 = LdBG::from_sequences(String::from("s1"), kmer_size, &all_sr_seqs);

    // Read all truth sequences.
    let all_truth_seqs: Vec<Vec<u8>> = process_reads(&truth_seq_urls, "truth");
    let t1 = LdBG::from_sequences(String::from("s1"), kmer_size, &all_truth_seqs);

    // let m = MLdBG::from_ldbgs(vec![l1.clone(), s1.clone()])
    //     .collapse()
    //     .assemble_all();

    let lr_contigs = l1.assemble_all();
    let lr_distances = distance_to_a_contig_end(&lr_contigs, kmer_size);

    let sr_contigs = s1.assemble_all();
    let sr_distances = distance_to_a_contig_end(&sr_contigs, kmer_size);

    // Read all test long reads.
    let all_test_lr_seqs: Vec<Vec<u8>> = process_reads(&test_long_read_seq_urls, "test long");
    let test_l1 = LdBG::from_sequences(String::from("test_l1"), kmer_size, &all_test_lr_seqs);

    // Read all test short reads.
    let all_test_sr_seqs: Vec<Vec<u8>> = process_reads(&test_short_read_seq_urls, "test short");
    let test_s1 = LdBG::from_sequences(String::from("test_s1"), kmer_size, &all_test_sr_seqs);

    // Read all test truth sequences.
    let all_test_truth_seqs: Vec<Vec<u8>> = process_reads(&test_truth_seq_urls, "test truth");
    let test_t1 = LdBG::from_sequences(String::from("test_t1"), kmer_size, &all_test_truth_seqs);

    let test_lr_contigs = test_l1.assemble_all();
    let test_lr_distances = distance_to_a_contig_end(&test_lr_contigs, kmer_size);

    let test_sr_contigs = test_s1.assemble_all();
    let test_sr_distances = distance_to_a_contig_end(&test_sr_contigs, kmer_size);


    // Configure GBDT.
    let mut cfg = Config::new();
    cfg.set_feature_size(7);
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
    // skydive::elog!(" - train/test={}/{}", 1.0 - test_split, test_split); //**** Commented out

    // Initialize GBDT algorithm.
    let mut gbdt = GBDT::new(&cfg);

    // let mut training_data: DataVec = Vec::new();
    // let mut test_data: DataVec = Vec::new();

    // let mut rng = rand::rngs::StdRng::seed_from_u64(0);

    let kmers = l1
        .kmers
        .keys()
        .chain(s1.kmers.keys())
        .chain(t1.kmers.keys());

    let mut training_data: DataVec = create_dataset_for_model(
        kmers,
        &lr_distances,
        &sr_distances,
        &l1,
        &s1,
        &t1,
    );

    // Train the decision trees.
    skydive::elog!(
        "Training GBDT model with {} training points...",
        training_data.len().to_formatted_string(&Locale::en)
    );
    gbdt.fit(&mut training_data);

    // Prepare test data.
    let test_kmers = test_l1
        .kmers
        .keys()
        .chain(test_s1.kmers.keys())
        .chain(test_t1.kmers.keys());

    let mut test_data: DataVec = create_dataset_for_model(
        test_kmers,
        &test_lr_distances,
        &test_sr_distances,
        &test_l1,
        &test_s1,
        &test_t1,
    );

    // Predict the test data.
    skydive::elog!("Computing accuracy on test data...");
    let prediction = gbdt.predict(&test_data);
    let pred_threshold: f32 = 0.5;

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
    skydive::elog!("Computing TPR and FPR at various thresholds...");
    let fpr_tpr: Vec<(f32, f32)> = compute_fpr_tpr(&test_data, &prediction);

    // Save TPR and FPR at various thresholds to a file.
    let csv_output = output.with_extension("csv");
    let mut writer = std::fs::File::create(csv_output).expect("Unable to create file");
    for (tpr, fpr) in &fpr_tpr {
        writeln!(writer, "{},{}", tpr, fpr).expect("Unable to write data");
    }
    skydive::elog!("TPR and FPR at various thresholds saved to {}", output.to_str().unwrap());

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

pub fn distance_to_a_contig_end(contigs: &Vec<Vec<u8>>, kmer_size: usize) -> HashMap<Vec<u8>, usize> {
    let mut distances = HashMap::new();

    for contig in contigs {
        for (distance_from_start, cn_kmer) in contig.windows(kmer_size).map(canonicalize_kmer).enumerate() {
            let distance_from_end = contig.len() - distance_from_start - kmer_size;

            distances.insert(cn_kmer, if distance_from_start < distance_from_end { distance_from_start } else { distance_from_end });
        }
    }

    distances
}
pub fn process_reads(read_seq_urls: &HashSet<Url>, read_type: &str)  -> Vec<Vec<u8>>{
    let mut all_seqs: Vec<Vec<u8>> = Vec::new();
    for read_seq_url in read_seq_urls {
        let basename = skydive::utils::basename_without_extension(
            &read_seq_url,
            &[".fasta.gz", ".fa.gz", ".fasta", ".fa"],
        );
        let fasta_path = read_seq_url.to_file_path().unwrap();

        skydive::elog!("Processing {}-read sample {}...", read_type, basename);

        let reader = bio::io::fasta::Reader::from_file(&fasta_path).unwrap();
        let all_reads: Vec<bio::io::fasta::Record> = reader.records().map(|r| r.unwrap()).collect();

        all_seqs.extend(
            all_reads
                .iter()
                .map(|r| r.seq().to_vec())
                .collect::<Vec<Vec<u8>>>(),
        );
    }

    all_seqs
}

pub fn plot_roc_curve(output: &PathBuf, roc_points: &[(f32, f32)]) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(output, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("ROC Curve", ("sans-serif", 50).into_font())
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(50)
        .build_cartesian_2d(0.0..1.0, 0.0..1.0)?;

    chart
        .configure_mesh()
        .x_desc("False Positive Rate")
        .y_desc("True Positive Rate")
        .draw()?;

    chart.draw_series(LineSeries::new(
        roc_points.iter().map(|(fpr, tpr)| (*fpr as f64, *tpr as f64)),
        &RED,
    ))?
        .label("ROC Curve")
        .legend(|(x, y)| PathElement::new([(x, y), (x + 20, y)], &RED));

    // Draw a diagonal dotted line for reference.
    chart.draw_series(LineSeries::new(
        [(0.0, 0.0), (1.0, 1.0)].iter().cloned(),
        &BLACK,
    ))?
        .label("Random")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLACK));

    chart.configure_series_labels().background_style(&WHITE).draw()?;

    Ok(())
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
        for (data, pred) in test_data.iter().zip(p.iter()) {
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

/// Computes precision, recall, and F1 score on test data.
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
        let lr_distance = *lr_distances.get(kmer).unwrap_or(&0) as f32;
        let sr_distance = *sr_distances.get(kmer).unwrap_or(&0) as f32;

        let lcov = l1.kmers.get(kmer).map_or(0, |lr| lr.coverage());
        let scov = s1.kmers.get(kmer).map_or(0, |sr| sr.coverage());
        let tcov = t1.kmers.get(kmer).map_or(0, |tr| tr.coverage());

        let data = Data::new_training_data(
            vec![
                if lcov > 0 { 1.0 } else { 0.0 },
                scov as f32,
                compressed_len_diff,
                entropy,
                gc_content,
                lr_distance,
                sr_distance,
            ],
            1.0,
            if tcov > 0 { 1.0 } else { 0.0 },
            None,
        );

        dataset.push(data);
    }

    dataset
}
