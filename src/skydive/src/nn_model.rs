use candle_core::{DType, Device, Tensor};
use candle_nn::{linear, Linear, Module, Optimizer, VarBuilder, LayerNorm, ops::sigmoid};
use candle_nn::loss::cross_entropy;
use rand::seq::SliceRandom;
use rand::thread_rng;

use crate::elog;



/// First we make struct of our model and mention our model architecture layers.
/// Think of it as all data you would put in a __init__() function of python with self keyword.
/// So we could access it later as a attribute of the model which we will do in forward method.
#[derive(Debug)]
pub struct KmerNN{
    fc1: Linear,
    ln1: LayerNorm,
    fc2: Linear,
    ln2: LayerNorm,
    fc3: Linear,
}


/// First argument is training data dimension which in our case will be 6 features and the second
/// one is varbuilder which is necessary to name of layers so we can load or save model
/// weights to/from python from/to rust.
// This function instantiates a new Model
impl KmerNN {
    pub fn new(in_dim: usize, vb: VarBuilder) -> candle_core::Result<Self> {
        let fc1 = linear(in_dim, 512, vb.pp("fc1"))?;

        // LayerNorm is a normalization layer which is used to normalize the input data. the parameters are weight: Tensor, bias: Tensor, eps: f64
        let ln1 = LayerNorm::new(
            Tensor::ones(&[512], DType::F32, vb.device())?,
            Tensor::zeros(&[512],  DType::F32,  vb.device())?, 1e-5
        );
        let fc2 = linear(512, 128, vb.pp("fc2"))?;
        let ln2 = LayerNorm::new(
            Tensor::ones(&[128], DType::F32, vb.device())?,
            Tensor::zeros(&[128],  DType::F32,  vb.device())?, 1e-5
        );
        let fc3 = linear(128, 1, vb.pp("fc3"))?;
        Ok(Self { fc1, ln1, fc2, ln2, fc3})
    }
}


/// Then we will implement Module trait on our model. Which is simple forward method for our model.
// forward pass of our model using Module trait
impl Module for KmerNN {
    fn forward(&self, xs: &Tensor) -> candle_core::Result<Tensor> {
        let x = self.fc1.forward(xs)?;
        let x = self.ln1.forward(&x)?;
        let x = x.relu()?;
        let x = self.fc2.forward(&x)?;
        let x = self.ln2.forward(&x)?;
        let x = x.relu()?;
        let x = self.fc3.forward(&x)?;
        let x = sigmoid(&x)?;
        Ok(x)
    }
}

// This function trains the model
pub fn train_model(
    model: &KmerNN,
    x_train: &Tensor,
    y_train: &Tensor,
    optimizer: &mut candle_nn::AdamW,
    epochs: usize,
    batch_size: usize,
    class_weights: (f32, f32),
) -> anyhow::Result<()> {
    let f_batches = create_batches(&x_train, batch_size);
    let l_batches = create_batches(&y_train, batch_size);

    for epoch in 0..epochs {
        // Todo:  encapsulate epoch-loss computation in a helper function to keep the loop cleaner
        let mut epoch_loss = Tensor::zeros(&[], DType::F32, x_train.device())?; // Initialize epoch loss, used to print the loss

        for (f_batch, l_batch) in f_batches.iter().zip(l_batches.iter()) {
            let output = model.forward(f_batch)?;
            // let loss = candle_nn::loss::mse(&output.squeeze(1)?, l_batch)?;
            let loss = weighted_mse_loss(&output.squeeze(1)?, l_batch, class_weights)?;

            optimizer.backward_step(&loss)?;
            epoch_loss = loss.clone();
        }

        if epoch % 10 == 0 {
            elog!("Epoch: {}, Loss: {}", epoch, epoch_loss.to_scalar::<f32>()?);
        }

    }
    Ok(())
}

// This function evaluates the model
pub fn evaluate_model(model: &KmerNN, x_test: &Tensor, y_test: &Tensor, x_test_table: KmerDataVec, class_weights: (f32, f32) ) -> anyhow::Result<()> {
    let output = model.forward(x_test)?;
    // let loss = candle_nn::loss::mse(&output.squeeze(1)?, y_test)?;
    let loss = weighted_mse_loss(&output.squeeze(1)?, y_test, class_weights)?;
    elog!("Test Loss: {}", loss.to_scalar::<f32>()?);

    // Convert the output tensor to a vector of predictions
    let preds: Vec<f32> = output.squeeze(1)?.to_vec1::<f32>()?.to_vec();

    // Evaluate the model on the test data
    elog!("Computing precision, recall, and F1 score on test data...");
    let (thresholds, accuracies, precisions, recalls, f1_scores) = evaluate_thresholds(&preds, &y_test, &x_test_table);

    elog!("Prediction thresholds: {:?}", thresholds);
    elog!("Precision: {:?}", precisions.iter().map(|&x| format!("{:.2}%", 100.0 * x)).collect::<Vec<String>>());
    elog!("Recall:    {:?}", recalls.iter().map(|&x| format!("{:.2}%", 100.0 * x)).collect::<Vec<String>>());
    elog!("F1 score:  {:?}", f1_scores.iter().map(|&x| format!("{:.2}%", 100.0 * x)).collect::<Vec<String>>());
    elog!("Accuracy:  {:?}", accuracies.iter().map(|&x| format!("{:.2}%", 100.0 * x)).collect::<Vec<String>>());
    elog!("ROC AUC:   {:.2}", roc_auc_score(&preds, &y_test.to_vec1::<f32>().unwrap()) );

    Ok(())
}

/// Calculate class weights based on the frequency of each class in the dataset.
pub fn calculate_class_weights(data: &[KmerData]) -> (f32, f32) {
    let total_samples = data.len() as f32;
    let class_0_count = data.iter().filter(|d| d.label == 0.0).count() as f32;
    let class_1_count = data.iter().filter(|d| d.label == 1.0).count() as f32;

    let weight_0 = total_samples / (2.0 * class_0_count);
    let weight_1 = total_samples / (2.0 * class_1_count);

    (weight_0, weight_1)
}


/// Modified loss function to include class weights.
///
/// output: The tensor with the predicted values from the model
/// target: The tensor with the actual values
/// weights: A tuple of the class weights
fn weighted_mse_loss(output: &Tensor, target: &Tensor, weights: (f32, f32)) -> candle_core::Result<Tensor> {
    // Compute the difference and squared difference
    let diff = (output - target)?;
    let squared_diff = diff.sqr()?;

    // Generate a tensor with weights based on the target values
    let zero_mask = target.eq(0.0)?.to_dtype(DType::F32)?; // Mask for where target == 0.0, cast to F32
    let one_mask = target.eq(1.0)?.to_dtype(DType::F32)?;  // Mask for where target == 1.0, cast to F32

    // Convert weights to tensors and broadcast them to the shape of `target`
    let weight_0 = Tensor::new(weights.0, target.device())?.broadcast_as(target.shape())?;
    let weight_1 = Tensor::new(weights.1, target.device())?.broadcast_as(target.shape())?;

    // Create the weighted tensor using the masks
    let weights_tensor = (zero_mask * weight_0)? + (one_mask * weight_1)?;

    // Compute the weighted loss
    let weighted_loss = (squared_diff * weights_tensor)?;
    weighted_loss.mean_all()
}


use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmerData {
    /// Features of the data
    pub feature: Vec<f32>,
    /// sample's weight. Used in training.
    pub weight: f64,
    /// The label value of the sample.
    pub label: f64,

}
/// The vector of the samples
pub type KmerDataVec = Vec<KmerData>;

impl KmerData {
    pub fn load_data(feature: Vec<f32>, weight: f64, label: f64) -> Self {
        KmerData {
            feature,
            weight,
            label,
        }
    }

    /// This function converts the data to a string. Can be used to save data to TSV file.
    pub fn to_vec_string(data: &KmerDataVec) -> Vec<String> {
        data.iter().map(|d| {
            format!(
                "{}\t{}",
                d.feature.iter().map(|f| f.to_string()).collect::<Vec<String>>().join("\t"),
                d.label
            )
        }).collect()
    }

    /// This function undersamples the majority class to balance the dataset.
    pub fn undersample_classes(data: &KmerDataVec) -> KmerDataVec {
        let mut rng = rand::thread_rng();
        let mut class_0_samples: Vec<&KmerData> = data.iter().filter(|d| d.label == 0.0).collect();
        let mut class_1_samples: Vec<&KmerData> = data.iter().filter(|d| d.label == 1.0).collect();

        class_0_samples.shuffle(&mut rng);
        class_1_samples.shuffle(&mut rng);

        let num_class_0 = class_0_samples.len();
        let num_class_1 = class_1_samples.len();

        let num_samples = num_class_0.min(num_class_1);
        class_0_samples.iter().take(num_samples)
            .chain(class_1_samples.iter().take(num_samples))
            .map(|&d| d.clone())
            .collect()
    }
}

/// This function prepares the data for training by converting it to tensors
/// and returning the feature and label tensors as a tuple. It also optionally
/// normalizes the feature tensor.
pub fn prepare_tensors(data: &[KmerData], device: &Device, normalize: bool) -> candle_core::Result<(Tensor, Tensor)> {
    // Convert the data to vectors while flattening the features and labels
    let features: Vec<f32> = data.iter().flat_map(|d| d.feature.clone()).collect();
    let labels: Vec<f32> = data.iter().map(|d| d.label as f32).collect();

    // Convert the data to tensors
    let feature_tensor = Tensor::from_slice(&features, (data.len(), data[0].feature.len()), device)?;
    let label_tensor = Tensor::from_slice(&labels, data.len(), device)?;

    // Normalize the feature tensor if the normalize parameter is true
    let feature_tensor = if normalize {
        let mean = feature_tensor.mean(0)?.broadcast_as(feature_tensor.shape())?;
        let var = feature_tensor.var(0)?.broadcast_as(feature_tensor.shape())?;
        let std = var.sqrt()?;
        (feature_tensor - mean)? / std
    } else {
        Ok(feature_tensor)
    };

    Ok((feature_tensor?, label_tensor))
}



/// Split the data into training and testing data randomly.
/// The split ratio is used to determine the percentage of data to use for training.
/// The function returns a tuple of training and testing data.
pub fn split_data(
    data: &[KmerData],
    split_ratio: f64 ,
) -> (Vec<KmerData>, Vec<KmerData>) {
    let mut rng = thread_rng();
    let mut data = data.to_vec();
    let split = (data.len() as f64 * split_ratio) as usize;

    data.shuffle(&mut rng);
    let split = (data.len() as f64 * split_ratio) as usize;
    let (train, test) = data.split_at(split);

    (train.to_vec(), test.to_vec())
}


/// Create batches of tensor data for training
pub fn create_batches(data: &Tensor, batch_size: usize) -> Vec<Tensor> {
    let mut batches = Vec::new();
    let num_samples = data.shape().dims()[0];
    let num_batches = (num_samples + batch_size - 1) / batch_size; // Round up for partial batch

    for i in 0..num_batches {
        let start = i * batch_size;
        let length = batch_size.min(num_samples - start); // Ensure we don't exceed remaining samples
        let batch = data.narrow(0, start, length).unwrap(); // Slice with start and length
        batches.push(batch);
    }
    batches
}



/// One hot encode the seqnence data
/// The function returns a vector of one hot encoded dna sequence
pub fn one_hot_encode_from_dna_string(seq: &str) -> Vec<f32> {
    let mut one_hot = vec![0.0; 5 * seq.len()];
    for (i, c) in seq.chars().enumerate() {
        match c.to_uppercase().next().unwrap() {
            'A' => one_hot[i * 5] = 1.0,
            'C' => one_hot[i * 5 + 1] = 1.0,
            'G' => one_hot[i * 5 + 2] = 1.0,
            'T' => one_hot[i * 5 + 3] = 1.0,
            'N' => one_hot[i * 5 + 4] = 1.0,
            _ => {}
        }
    }
    one_hot
}

/// One hot encode the seqnence data from vec<u8>
pub fn one_hot_encode_from_dna_ascii(seq: &[u8]) -> Vec<f32> {
    let mut one_hot = vec![0.0; 5 * seq.len()];
    for (i, c) in seq.iter().enumerate() {
        match c.to_ascii_uppercase() {
            b'A' => one_hot[i * 5] = 1.0,
            b'C' => one_hot[i * 5 + 1] = 1.0,
            b'G' => one_hot[i * 5 + 2] = 1.0,
            b'T' => one_hot[i * 5 + 3] = 1.0,
            b'N' => one_hot[i * 5 + 4] = 1.0,
            _ => {}
        }
    }
    one_hot
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
