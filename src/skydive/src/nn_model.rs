use candle_core::{DType, Device, Tensor};
use candle_nn::{linear, Linear, Module, Optimizer, VarBuilder, LayerNorm};
use rand::seq::SliceRandom;
use rand::thread_rng;

const DEVICE: Device = Device::Cpu;


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
        let fc1 = linear(in_dim, 128, vb.pp("fc1"))?;

        // LayerNorm is a normalization layer which is used to normalize the input data. the parameters are weight: Tensor, bias: Tensor, eps: f64
        let ln1 = LayerNorm::new(
            Tensor::ones(&[128], DType::F32, vb.device())?,
            Tensor::zeros(&[128],  DType::F32,  vb.device())?, 1e-5
        );
        let fc2 = linear(128, 64, vb.pp("fc2"))?;
        let ln2 = LayerNorm::new(
            Tensor::ones(&[64], DType::F32, vb.device())?,
            Tensor::zeros(&[64],  DType::F32,  vb.device())?, 1e-5
        );
        let fc3 = linear(64, 1, vb.pp("fc3"))?;
        Ok(Self { fc1, ln1, fc2, ln2, fc3 })
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
        Ok(x)
    }
}

// This function trains the model
pub fn train_model(
    model: &KmerNN,
    x_train: &Tensor,
    y_train: &Tensor,
    optimizer: &mut candle_nn::AdamW,
    epochs: usize
) -> anyhow::Result<()> {
    for epoch in 0..epochs {

        let output = model.forward(x_train)?;
        let loss = candle_nn::loss::mse(&output.squeeze(1)?, y_train)?;

        optimizer.backward_step(&loss)?;
        if (epoch + 1) % 10 == 0 {
            println!("Epoch: {} | Train Loss: {} ", epoch + 1, loss.to_scalar::<f32>()?);
        }
    }
    Ok(())
}

// This function evaluates the model
pub fn evaluate_model(model: &KmerNN, x_test: &Tensor, y_test: &Tensor) -> anyhow::Result<()> {
    let output = model.forward(x_test)?;
    let loss = candle_nn::loss::mse(&output.squeeze(1)?, y_test)?;
    println!("Test Loss: {}", loss.to_scalar::<f32>()?);
    Ok(())
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
}


/// This function prepares the data for training by converting it to tensors
/// and returning the feature and label tensors as a tuple.
pub fn prepare_tensors(data: &[KmerData], device: &Device) -> (Tensor, Tensor) {
    // Convert the data to vectors while flattening the features and labels
    let features: Vec<f32> = data.iter().flat_map(|d| d.feature.clone()).collect();
    let labels: Vec<f32> = data.iter().map(|d| d.label as f32).collect();

    // Convert the data to tensors
    let feature_tensor = Tensor::from_slice(
        &features, (data.len(), data[0].feature.len()), device
    ).unwrap();
    let label_tensor = Tensor::from_slice(
        &labels, data.len(), device
    ).unwrap();


    (feature_tensor, label_tensor)
}




/// Split the data into training and testing data randomly.
/// The data is split into 80% training and 20% testing data.
/// The function returns a tuple of training and testing data.
pub fn split_data(data: &[KmerData]) -> (Vec<KmerData>, Vec<KmerData>) {
    let mut rng = thread_rng();
    let mut data = data.to_vec();

    data.shuffle(&mut rng);
    let split = (data.len() as f64 * 0.8) as usize;
    let (train, test) = data.split_at(split);

    (train.to_vec(), test.to_vec())
}


/// One hot encode the seqnence data
/// The function returns a vector of one hot encoded dna sequence
pub fn one_hot_encode_from_dna_string(seq: &str) -> Vec<f32> {
    let mut one_hot = vec![0.0; 4 * seq.len()];
    for (i, c) in seq.chars().enumerate() {
        match c {
            'A' => one_hot[i * 4] = 1.0,
            'C' => one_hot[i * 4 + 1] = 1.0,
            'G' => one_hot[i * 4 + 2] = 1.0,
            'T' => one_hot[i * 4 + 3] = 1.0,
            _ => {}
        }
    }
    one_hot
}

/// One hot encode the seqnence data from vec<u8>
pub fn one_hot_encode_from_dna_ascii(seq: &[u8]) -> Vec<f32> {
    let mut one_hot = vec![0.0; 4 * seq.len()];
    for (i, c) in seq.iter().enumerate() {
        match c {
            65 => one_hot[i * 4] = 1.0,
            67 => one_hot[i * 4 + 1] = 1.0,
            84 => one_hot[i * 4 + 2] = 1.0,
            71 => one_hot[i * 4 + 3] = 1.0,
            _ => {}
        }
    }
    one_hot
}

