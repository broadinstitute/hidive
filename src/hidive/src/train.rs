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
use gbdt::decision_tree::{Data, DataVec, PredVec};
use gbdt::gradient_boost::GBDT;

use plotters::prelude::*;
use skydive::ldbg::LdBG;
use skydive::record::Record;
use std::io::Write;
use rand::seq::SliceRandom;
use url::Url;

use skydive::elog;

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
    let long_read_seq_urls = skydive::parse::parse_file_names(long_read_seq_paths);
    let short_read_seq_urls = skydive::parse::parse_file_names(short_read_seq_paths);
    let truth_seq_urls = skydive::parse::parse_file_names(truth_seq_paths);
    let test_long_read_seq_urls = skydive::parse::parse_file_names(test_long_read_seq_paths);
    let test_short_read_seq_urls = skydive::parse::parse_file_names(test_short_read_seq_paths);
    let test_truth_seq_urls = skydive::parse::parse_file_names(test_truth_seq_paths);

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
    cfg.set_feature_size(6);
    cfg.set_max_depth(4);
    cfg.set_min_leaf_size(5);
    cfg.set_loss(&loss2string(&Loss::BinaryLogitraw));
    cfg.set_shrinkage(0.05);
    cfg.set_iterations(iterations);
    // cfg.set_debug(debug);

    elog!("Training GBDT model with:");
    elog!(" - feature_size={}", cfg.feature_size);
    elog!(" - max_depth={}", cfg.max_depth);
    elog!(" - min_leaf_size={}", cfg.min_leaf_size);
    elog!(" - loss={}", loss2string(&cfg.loss));
    elog!(" - shrinkage={}", cfg.shrinkage);
    elog!(" - iterations={}", cfg.iterations);
    // elog!(" - train/test={}/{}", 1.0 - test_split, test_split); //**** Commented out

    // Initialize GBDT algorithm.
    // let mut gbdt = GBDT::new(&cfg);

    // There are three types of data: training, validation, test.
    // The training data is used to train the model.
    // The validation data is used to tune the model.
    // The test data is used to evaluate the model.

    // Create a dataset for the ML model.
    let kmers = l1
        .kmers
        .keys()
        .chain(s1.kmers.keys())
        .chain(t1.kmers.keys());

    let pre_kmer_with_feature_data: DataVec = create_dataset_for_model(
        kmers,
        &lr_distances,
        &sr_distances,
        &l1,
        &s1,
        &t1,
    );
    let kmer_with_feature_data: DataVec = undersample_classes(&pre_kmer_with_feature_data);


    // Save the data to a TSV file
    elog!("Saving data to file...");
    let mut tsv_data: Vec<String> = Vec::new();
    let headers = vec![
        "lr-present",
        "sr-coverage",
        "strand-bias",
        "hc-length-difference",
        "shannon-entropy",
        "gc-content",
        "truth-present",
    ];
    tsv_data.insert(0, headers.join("\t"));
    for data in &kmer_with_feature_data {
        tsv_data.push(
            data.feature
                .iter()
                .map(|f| f.to_string())
                .chain(vec![data.label.to_string()])
                .collect::<Vec<String>>()
                .join("\t"),
        );
    }
    let data_output = output.with_file_name("kmer_feature_data.tsv");
    std::fs::write(&data_output, tsv_data.join("\n")).expect("Unable to write data to file");


    // split data 80 : 20 for training and testing
    let mut training_data: DataVec;
    let mut validation_data: DataVec;
    (training_data, validation_data) = split_data(kmer_with_feature_data, 0.8);

    elog!("Training data size: {}", &training_data.len().to_formatted_string(&Locale::en));
    elog!("Validation data size: {}", &validation_data.len().to_formatted_string(&Locale::en));


    // Train the decision trees.
    elog!("Training GBDT model with {} training points...",
        &training_data.len().to_formatted_string(&Locale::en)
    );

    ////////////////
    //// This section is for training the model with different depths.
    //// The model is trained with different depths and the accuracy and roc are calculated and
    //// written to a file.
    //// This can be removed in production.
    ///////////////

    // varialbe to hold validation evaluation
    let mut training_validation_accuracy: Vec<Vec<f32>> = Vec::new();
    let mut training_validation_roc: Vec<Vec<f32>> = Vec::new();
    // let mut training_validation_loss: Vec<PredVec> = Vec::new();
    let validation_truth = &validation_data.iter().map(|d| d.label).collect::<Vec<f32>>();
    let training_data_truth = &training_data.iter().map(|d| d.label).collect::<Vec<f32>>();
    for depth in 1..10 {
        let mut cfg = Config::new();
        cfg.set_feature_size(6);
        cfg.set_min_leaf_size(5);
        cfg.set_loss(&loss2string(&Loss::BinaryLogitraw));
        cfg.set_shrinkage(0.05);
        cfg.set_iterations(iterations);
        cfg.set_max_depth(depth);

        let mut gbdt = GBDT::new(&cfg);
        gbdt.fit(&mut training_data);
        let training_validation_prediction = gbdt.predict(&validation_data);
        let training_data_prediction = gbdt.predict(&training_data);

        // Evaluate the model on the validation data
        let (thresholds, accuracies, precisions, recalls, f1_scores) = evaluate_thresholds(&training_validation_prediction, &validation_truth);
        let (precision, recall, f1_score) = compute_precision_recall_f1(&validation_truth, &training_validation_prediction, 0.5);
        let val_roc_auc = roc_auc_score(&training_validation_prediction, &validation_truth);
        let train_roc_auc = roc_auc_score(&training_data_prediction, &training_data_truth);
        elog!("Validation data evaluation at depth {}:", depth);

        training_validation_accuracy.push(accuracies);
        training_validation_roc.push(vec![val_roc_auc, train_roc_auc]);

        //Todo: retrieve loss

    }

    // Save accuracy to a TSV file
    let mut accuracy_tsv_data: Vec<String> = Vec::new();
    let accuracy_headers = vec!["accuracy_0.1", "accuracy_0.2", "accuracy_0.3", "accuracy_0.4", "accuracy_0.5", "accuracy_0.6", "accuracy_0.7", "accuracy_0.8", "accuracy_0.9"];
    accuracy_tsv_data.push(accuracy_headers.join("\t"));
    for accuracy in training_validation_accuracy {
        accuracy_tsv_data.push(accuracy.iter().map(|a| a.to_string()).collect::<Vec<String>>().join("\t"));
    }
    let accuracy_data_output = output.with_file_name("training_validation_accuracy.tsv");
    std::fs::write(&accuracy_data_output, accuracy_tsv_data.join("\n")).expect("Unable to write data to file");

    // Save roc to a TSV file
    let roc_tsv_output = output.with_file_name("training_validation_roc.tsv");
    let roc_headers = vec!["validation_roc_auc", "training_roc_auc"];
    let mut roc_tsv_data: Vec<String> = Vec::new();
    roc_tsv_data.push(roc_headers.join("\t"));
    for roc in training_validation_roc {
        roc_tsv_data.push(roc.iter().map(|a| a.to_string()).collect::<Vec<String>>().join("\t"));
    }
    std::fs::write(&roc_tsv_output, roc_tsv_data.join("\n")).expect("Unable to write data to file");
    ////////////////////////
    //// End of training with different depths
    ////////////////////////


    // Initialize GBDT algorithm.
    let mut gbdt = GBDT::new(&cfg);
    gbdt.fit(&mut training_data);

    // Run model on validation data
    elog!("Computing accuracy on validation data...");
    let validation_prediction = gbdt.predict(&validation_data);
    let validation_truth = &validation_data.iter().map(|d| d.label).collect::<Vec<f32>>();

    // Evaluate the model on the test data
    elog!("Computing precision, recall, and F1 score on test data...");
    let (thresholds, accuracies, precisions, recalls, f1_scores) = evaluate_thresholds(&validation_prediction, &validation_truth);

    log_evaluation_metrics(&thresholds, &precisions, &recalls, &f1_scores, &accuracies, &validation_prediction, &validation_truth);

    compute_and_save_tpr_fpr_roc(&output, &validation_data, &validation_prediction, "validation"
    ).expect("Unable to compute and save TPR, FPR, and ROC curve");


    // Prepare test data.
    let test_kmers = test_l1
        .kmers
        .keys()
        .chain(test_s1.kmers.keys())
        .chain(test_t1.kmers.keys());

    let test_data: DataVec = create_dataset_for_model(
        test_kmers,
        &test_lr_distances,
        &test_sr_distances,
        &test_l1,
        &test_s1,
        &test_t1,
    );

    // Predict the test data.
    elog!("Computing accuracy on test data...");
    let test_prediction = gbdt.predict(&test_data);
    let test_truth = test_data.iter().map(|d| d.label).collect::<Vec<f32>>();

    // Evaluate the model on the test data
    elog!("Computing precision, recall, and F1 score on test data...");
    let (thresholds, accuracies, precisions, recalls, f1_scores) = evaluate_thresholds(&test_prediction, &test_truth);

    log_evaluation_metrics(&thresholds, &precisions, &recalls, &f1_scores, &accuracies, &test_prediction, &test_truth);

    compute_and_save_tpr_fpr_roc(&output, &test_data, &test_prediction, "test").expect(
        "Unable to compute and save TPR, FPR, and ROC curve"
    );

    // Save the model to a file.
    gbdt.save_model(output.to_str().unwrap()).expect("Unable to save model");
    elog!("Model saved to {}", output.to_str().unwrap());
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
            read_seq_url,
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
        .legend(|(x, y)| PathElement::new([(x, y), (x + 20, y)], RED));

    // Draw a diagonal dotted line for reference.
    chart.draw_series(LineSeries::new(
        [(0.0, 0.0), (1.0, 1.0)].iter().cloned(),
        &BLACK,
    ))?
        .label("Random")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLACK));

    chart.configure_series_labels().background_style(WHITE).draw()?;

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

/// Computes and saves TPR and FPR at various thresholds and plots the ROC curve.
///
/// # Arguments
/// * `output` - The output path for the ROC curve image.
/// * `test_data` - The test data used for evaluation.
/// * `prediction` - The predicted values from the model.
///
/// # Returns
/// None
fn compute_and_save_tpr_fpr_roc(output: &PathBuf, test_data: &DataVec, prediction: &Vec<f32>, basename: &str )-> Result<(), Box<dyn std::error::Error>> {
    // TPR and FPR calculations at various thresholds.
    elog!("Computing TPR and FPR at various thresholds...");
    let fpr_tpr: Vec<(f32, f32)> = compute_fpr_tpr(test_data, prediction);

    // Save TPR and FPR at various thresholds to a file.
    let csv_output = output.with_file_name(format!("{}.data_tpr_fpr.csv", basename));
    let mut writer = std::fs::File::create(&csv_output).expect("Unable to create file");
    for (tpr, fpr) in &fpr_tpr {
        writeln!(writer, "{},{}", tpr, fpr).expect("Unable to write data");
    }
    elog!("TPR and FPR at various thresholds saved to {}", &csv_output.to_str().unwrap());

    // Create a ROC curve.
    let png_output = output.with_file_name(format!("{}.data_roc_curve.png", basename));
    plot_roc_curve(&png_output, &fpr_tpr).expect("Unable to plot ROC curve");
    elog!("ROC curve saved to {}", &png_output.to_str().unwrap());

    Ok(())
}



/// Creates a dataset for the ML model.
///
/// The dataset is created by iterating over all kmers and computing the features for each kmer.
/// The features are:
/// - present in long reads
/// - coverage in short reads
/// - measure of strand bias (0.5 = balanced, 1.0 = all on one strand)
/// - homopolymer compression length difference
/// - shannon entropy
/// - gc content
/// - distance to nearest long read contig end
/// - distance to nearest short read contig end
///
/// returns a vector of Data objects.
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

/// Splits the data into training and testing sets.
fn split_data(data: DataVec, ratio: f32) -> (DataVec, DataVec) {
    let split_index = (data.len() as f32 * ratio) as usize;
    let (train_data, test_data) = data.split_at(split_index);
    (train_data.to_vec(), test_data.to_vec())
}


/// Calculate the ROC AUC score. This is the area under the ROC curve.
/// The ROC curve is a plot of the true positive rate (TPR) against the false positive rate (FPR)
/// for the different possible thresholds of a binary classifier.
pub fn roc_auc_score(preds: &Vec<f32>, labels: &Vec<f32>) -> f32 {
    let mut roc_points: Vec<(f32, f32)> = Vec::new();
    let mut num_true_positives = 0u32;
    let mut num_false_positives = 0u32;
    let mut num_true_negatives = 0u32;
    let mut num_false_negatives = 0u32;

    // Iterate over different thresholds
    for threshold in 0..100 {
        let pred_threshold = threshold as f32 / 100.0;
        for (pred, truth) in preds.iter().zip(labels.iter()) {
            let call = if *pred > pred_threshold { 1.0 } else { 0.0 };
            if (call - *truth).abs() < f32::EPSILON {
                if *truth > 0.5 { num_true_positives += 1; }  // True Positive
                else { num_true_negatives += 1; }  // True Negative
            } else {
                if *truth < 0.5 { num_false_positives += 1; }  // False Positive
                if *truth > 0.5 { num_false_negatives += 1; }  // False Negative
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

        roc_points.push((fpr, tpr));
    }

    // Sort the ROC points by FPR
    roc_points.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    // Calculate the area under the ROC curve
    let mut auc = 0.0;
    for i in 1..roc_points.len() {
        let (fpr1, tpr1) = roc_points[i - 1];
        let (fpr2, tpr2) = roc_points[i];
        auc += (tpr1 + tpr2) * (fpr2 - fpr1) / 2.0;
    }

    auc
}

/// Evaluates the model on the test data.
fn evaluate_thresholds(preds: &[f32], truth_lables: &Vec<f32>) -> (Vec<f32>, Vec<f32>, Vec<f32>, Vec<f32>, Vec<f32>) {
    let mut thresholds = vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
    let mut accuracies = Vec::new();
    let mut precisions = Vec::new();
    let mut recalls = Vec::new();
    let mut f1_scores = Vec::new();

    for threshold in thresholds.iter() {
        let mut num_correct = 0u32;
        let mut num_total = 0u32;
        let mut p: Vec<f32> = Vec::new();

        for (pred, truth) in preds.iter().zip(truth_lables.iter()) {
            p.push(*pred);
            let call: f32 = if *pred > *threshold { 1.0 } else { 0.0 };
            if (call as f32 - truth).abs() < f32::EPSILON {
                num_correct += 1;
            }
            num_total += 1;
        }

        // elog!(
        // "Predicted accuracy: {}/{} ({:.2}%)",
        // num_correct.to_formatted_string(&Locale::en),
        // num_total.to_formatted_string(&Locale::en),
        // 100.0 * num_correct as f32 / num_total as f32
        // );

        let accuracy = num_correct as f32 / num_total as f32;
        let (precision, recall, f1_score) = compute_precision_recall_f1(&truth_lables, &p, *threshold);

        accuracies.push(accuracy);
        precisions.push(precision);
        recalls.push(recall);
        f1_scores.push(f1_score);
    }

    (thresholds, accuracies, precisions, recalls, f1_scores)
}


fn log_evaluation_metrics(
    thresholds: &Vec<f32>,
    precisions: &Vec<f32>,
    recalls: &Vec<f32>,
    f1_scores: &Vec<f32>,
    accuracies: &Vec<f32>,
    prediction: &Vec<f32>,
    test_truth: &Vec<f32>,
) {
    elog!("Prediction thresholds: {:?}", thresholds);
    elog!("Precision: {:?}", precisions.iter().map(|&x| format!("{:.2}%", 100.0 * x)).collect::<Vec<String>>());
    elog!("Recall:    {:?}", recalls.iter().map(|&x| format!("{:.2}%", 100.0 * x)).collect::<Vec<String>>());
    elog!("F1 score:  {:?}", f1_scores.iter().map(|&x| format!("{:.2}%", 100.0 * x)).collect::<Vec<String>>());
    elog!("Accuracy:  {:?}", accuracies.iter().map(|&x| format!("{:.2}%", 100.0 * x)).collect::<Vec<String>>());
    elog!("ROC AUC:   {:.2}", roc_auc_score(&prediction, &test_truth));
}


/// Computes precision, recall, and F1 score on test data.
pub fn compute_precision_recall_f1(truth_labels: &Vec<f32>, p: &Vec<f32>, pred_threshold: f32) -> (f32, f32, f32) {
    let (mut num_true_positives, mut num_false_positives, mut num_false_negatives) = (0u32, 0u32, 0u32);

    for (truth, pred) in truth_labels.iter().zip(p.iter()) {
        let call = if *pred > pred_threshold { 1.0 } else { 0.0 };
        if (call - truth).abs() < f32::EPSILON {
            if truth > &0.5 { num_true_positives += 1; }
        } else {
            if truth < &0.5 { num_false_positives += 1; }
            if truth > &0.5 { num_false_negatives += 1; }
        }
    }

    let precision = num_true_positives as f32 / (num_true_positives + num_false_positives) as f32;
    let recall = num_true_positives as f32 / (num_true_positives + num_false_negatives) as f32;
    let f1_score = 2.0 * precision * recall / (precision + recall);

    (precision, recall, f1_score)
}

pub fn undersample_classes(data: &DataVec) -> DataVec {
    let mut rng = rand::thread_rng();
    let mut class_0_samples: Vec<&Data> = data.iter().filter(|d| d.label == 0.0).collect();
    let mut class_1_samples: Vec<&Data> = data.iter().filter(|d| d.label == 1.0).collect();

    // Shuffle the samples
    class_0_samples.shuffle(&mut rng);
    class_1_samples.shuffle(&mut rng);

    let num_class_0 = class_0_samples.len();
    let num_class_1 = class_1_samples.len();

    // Take the minimum number of samples from each class
    let num_samples = num_class_0.min(num_class_1);

    // combined dataset of the two classes using the minimum number of samples from each class
    let mut undersampled_data: DataVec = Vec::new();
    undersampled_data.extend(class_0_samples.into_iter().take(num_samples).cloned());
    undersampled_data.extend(class_1_samples.into_iter().take(num_samples).cloned());
    undersampled_data.shuffle(&mut rng);

    undersampled_data
}


