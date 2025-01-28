use skydive::utils::canonicalize_kmer;
use std::collections::HashMap;
use std::collections::HashSet;
use std::path::PathBuf;
use std::iter::Chain;
use std::collections::hash_map::Keys;

use skydive::nn_model::{KmerData, KmerDataVec, KmerNN, train_model, evaluate_model, prepare_tensors, split_data, one_hot_encode_from_dna_ascii, calculate_class_weights};
use candle_nn::{Module, Optimizer, VarBuilder, VarMap};
use candle_core::{DType, Device, Tensor};
use num_format::Locale::el;

const DEVICE: Device = Device::Cpu;


use plotters::prelude::*;
use skydive::ldbg::LdBG;
use skydive::record::Record;
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

    // Load the training data
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

    // Load the test data
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


    // Train model

    let kmers = l1
        .kmers
        .keys()
        .chain(s1.kmers.keys())
        .chain(t1.kmers.keys());

    let data: KmerDataVec = create_dataset_for_model(
        kmers,
        &lr_distances,
        &sr_distances,
        &l1,
        &s1,
        &t1,
    );


    // Save the data to a TSV file
    elog!("Saving data to file...");
    let mut tsv_data: Vec<String> = KmerData::to_vec_string(&data);
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
    // Write the data to a TSV file
    let data_output = output.with_extension("tsv");
    std::fs::write(&data_output, tsv_data.join("\n")).expect("Unable to write data to file");


    elog!("Splitting data into training and validation sets...");
    let (training_data, validation_data) = split_data(&data);
    // get dimensions of the training data
    let  n_features = training_data[0].feature.len();

    // f_train  is the features and l_train is the labels
    elog!("Preparing tensors...");
    let (f_train, l_train) = prepare_tensors(&training_data, &DEVICE, true).unwrap();
    let (f_validation, l_validation) = prepare_tensors(&*KmerData::undersample_classes(&validation_data), &DEVICE, true).unwrap();



    // Create the model
    // varmap is used to store the variables, which are the weights and biases of the model
    let varmap = VarMap::new();
    // vb is used to create the variables, used to create and manage the weights and biases of the model
    let vb = VarBuilder::from_varmap(&varmap, DType::F32, &DEVICE);

    elog!("Creating model...");
    let model = match KmerNN::new(n_features, vb) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("Error creating model: {}", e);
            return;
        }
    };

    // Create the Optimizer
    let optim_config = candle_nn::ParamsAdamW{
        lr: 1e-1,
        ..Default::default()
    };

    elog!("Creating optimizer...");
    let mut optimizer = match candle_nn::AdamW::new(varmap.all_vars(), optim_config) {
        Ok(opt) => opt,
        Err(e) => {
            eprintln!("Error creating optimizer: {}", e);
            return;
        }
    };

    // let batch_size = 32;
    let batch_size = 512;
    let epochs = 100;

    // Calculate weights for the loss function
    let class_weights = calculate_class_weights(&training_data);
    eprintln!("Training Class weights (0, 1): {:?}", class_weights);


    elog!("Training model...");
    if let Err(e) = train_model(
        &model, &f_train, &l_train, &mut optimizer, epochs, batch_size, class_weights
    ) {
        eprintln!("Error training model: {}", e);
        return;
    }


    // Calculate weights for the loss function
    let eval_class_weights = calculate_class_weights(&*KmerData::undersample_classes(&validation_data));
    eprintln!("Validation Class weights (0, 1): {:?}", eval_class_weights);
    elog!("Evaluating model...");
    if let Err(e) = evaluate_model(&model, &f_validation, &l_validation, eval_class_weights) {
        eprintln!("Error evaluating model: {}", e);
        return;
    }

    // Evaluate accuracy of the model on the test data.
    let test_kmers = l1
        .kmers
        .keys()
        .chain(s1.kmers.keys())
        .chain(t1.kmers.keys());

    let test_data: KmerDataVec = create_dataset_for_model(
        test_kmers,
        &lr_distances,
        &sr_distances,
        &l1,
        &s1,
        &t1,
    );

    // f_test  is the features and l_test is the labels
    elog!("Preparing tensors...");
    // let mut num_correct = 0u32;
    // let mut num_total = 0u32;
    let mut p: Vec<f32> = Vec::new();
    // let pred_threshold: f32 = 0.5;

    let (f_test, l_test) = prepare_tensors(&test_data, &DEVICE, true).unwrap();

    let model_test_output = model.forward(&f_test).unwrap();
    let preds: Vec<f32> = model_test_output.squeeze(1).unwrap().to_vec1().unwrap();

    // for (pred, truth) in preds.iter().zip(l_test.to_vec1::<f32>().unwrap().iter()) {
    //     p.push(*pred);
    //     let call: f32 = if *pred > 0.5 { 1.0 } else { 0.0 };
    //     if (call as f32 - *truth as f32).abs() < f32::EPSILON {
    //         num_correct += 1;
    //     }
    //     num_total += 1;
    // }
    //
    // Precision, Recall, and F1 score calculations.
    // let (precision, recall, f1_score) = compute_precision_recall_f1(&test_data, &p, pred_threshold);
    // elog!("Prediction threshold: {:.2}", pred_threshold);
    // elog!("Precision: {:.2}%", 100.0 * precision);
    // elog!("Recall: {:.2}%", 100.0 * recall);
    // elog!("F1 score: {:.2}%", 100.0 * f1_score);
    // elog!("Accuracy: {:.2}%", 100.0 * num_correct as f32 / num_total as f32);

    // Evaluate the model on the test data
    elog!("Computing precision, recall, and F1 score on test data...");
    let (thresholds, accuracies, precisions, recalls, f1_scores) = evaluate_thresholds(&preds, &l_test, &test_data);

    elog!("Prediction thresholds: {:?}", thresholds);
    elog!("Precision: {:?}", precisions.iter().map(|&x| format!("{:.2}%", 100.0 * x)).collect::<Vec<String>>());
    elog!("Recall:    {:?}", recalls.iter().map(|&x| format!("{:.2}%", 100.0 * x)).collect::<Vec<String>>());
    elog!("F1 score:  {:?}", f1_scores.iter().map(|&x| format!("{:.2}%", 100.0 * x)).collect::<Vec<String>>());
    elog!("Accuracy:  {:?}", accuracies.iter().map(|&x| format!("{:.2}%", 100.0 * x)).collect::<Vec<String>>());
    elog!("ROC AUC:   {:.2}", roc_auc_score(&preds, &l_test.to_vec1::<f32>().unwrap()) );

    // Save the model to file
    let tensor_output = output.with_extension("safetensor");
    varmap.save(&tensor_output).expect("Unable to save tensor");
    elog!("Model saved to {}", tensor_output.display());

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

        elog!("Processing {}-read sample {}...", read_type, basename);

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
pub fn compute_fpr_tpr(test_data: &KmerDataVec, p: &Vec<f32>) -> Vec<(f32, f32)> {
    elog!("Computing FPR and TPR at various thresholds...");
    let mut fpr_tpr: Vec<(f32, f32)> = Vec::new();

    // Iterate over different thresholds
    for threshold in 0..100 {
        let pred_threshold = threshold as f32 / 100.0;
        let (mut num_true_positives, mut num_false_positives, mut num_false_negatives, mut num_true_negatives) = (0u32, 0u32, 0u32, 0u32);

        // Iterate over test data and predictions
        for (data, pred) in test_data.iter().zip(p.iter()) {
            let truth = data.label;
            let call = if *pred > pred_threshold { 1.0 } else { 0.0 };

            if (call - truth).abs() < f64::EPSILON {
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
pub fn compute_precision_recall_f1(test_data: &Vec<KmerData>, p: &Vec<f32>, pred_threshold: f32) -> (f32, f32, f32) {
    let (mut num_true_positives, mut num_false_positives, mut num_false_negatives) = (0u32, 0u32, 0u32);

    for (data, pred) in test_data.iter().zip(p.iter()) {
        let truth = data.label;
        let call = if *pred > pred_threshold { 1.0 } else { 0.0 };
        if (call - truth).abs() < f64::EPSILON {
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
) -> Vec<KmerData> {
    elog!("Creating dataset for model...");
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

        let mut features = Vec::new();
        features.extend(vec![
            if lcov > 0 { 1.0 } else { 0.0 },    // present in long reads
            scov_total as f32,                   // coverage in short reads
            strand_ratio as f32,                 // measure of strand bias (0.5 = balanced, 1.0 = all on one strand)
            compressed_len_diff,                 // homopolymer compression length difference
            entropy,                             // shannon entropy
            gc_content,                          // gc content
            // lr_distance,                         // distance to nearest long read contig end
            // sr_distance,                         // distance to nearest short read contig end
        ]);

        let data = KmerData::load_data(
            features,
            1.0,
            if tcov > 0 { 1.0 } else { 0.0 }, // present in truth
        );

        dataset.push(data);
    }

    dataset
}


/// Evaluates the model on the test data.
fn evaluate_thresholds(preds: &[f32], l_test: &Tensor, test_data: &KmerDataVec) -> (Vec<f32>, Vec<f32>, Vec<f32>, Vec<f32>, Vec<f32>) {
    let mut thresholds = vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
    let mut accuracies = Vec::new();
    let mut precisions = Vec::new();
    let mut recalls = Vec::new();
    let mut f1_scores = Vec::new();

    for threshold in thresholds.iter() {
        let mut num_correct = 0u32;
        let mut num_total = 0u32;
        let mut p: Vec<f32> = Vec::new();

        for (pred, truth) in preds.iter().zip(l_test.to_vec1::<f32>().unwrap().iter()) {
            p.push(*pred);
            let call: f32 = if *pred > *threshold { 1.0 } else { 0.0 };
            if (call as f32 - *truth as f32).abs() < f32::EPSILON {
                num_correct += 1;
            }
            num_total += 1;
        }

        let accuracy = num_correct as f32 / num_total as f32;
        let (precision, recall, f1_score) = compute_precision_recall_f1(&test_data, &p, *threshold);

        accuracies.push(accuracy);
        precisions.push(precision);
        recalls.push(recall);
        f1_scores.push(f1_score);
    }

    (thresholds, accuracies, precisions, recalls, f1_scores)
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


