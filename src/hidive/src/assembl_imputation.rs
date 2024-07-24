// Import necessary standard library modules
use std::collections::HashSet;
use std::collections::HashMap;
use std::path::PathBuf;
use serde_json::Value;
// Import the Absolutize trait to convert relative paths to absolute paths
use path_absolutize::Absolutize;

extern crate ndarray;

use ndarray::{Array,Array1, Array2, array};
use ndarray::Axis;
use linfa::traits::Transformer;
use std::f64::NAN;

// Import the Url type to work with URLs
use url::Url;

// Import the skydive module, which contains the necessary functions for staging data
use skydive;

pub fn find_all_read (graph: &skydive::agg::GraphicalGenome) -> HashSet<String> {
    let mut read_sets = HashSet::new();
    let edgelist = graph.edges.keys();
    for item in edgelist{
        if let Some(edge_data) = graph.edges.get(item) {
            if let Some(readlist) =  edge_data.get("reads"){
                if let Value::Array(array) = readlist {
                    for read in array.iter() {
                        if let Value::String(read_str) = read {
                            read_sets.insert(read_str.clone());
                        }
                    }
                }
            }               

        }
    }
    read_sets
}

fn construct_anchor_table(graph: &skydive::agg::GraphicalGenome) -> (Vec<Vec<Option<f64>>>,Vec<String>) {
    let read_sets = find_all_read(graph);
    let mut sorted_read_sets: Vec<String> = read_sets.into_iter().collect();
    sorted_read_sets.sort();

    let mut r: HashMap<&String, usize> = HashMap::new();
    for (index, read) in sorted_read_sets.iter().enumerate() {
        r.insert(read, index);
    }

    let mut anchorlist: Vec<&String> = graph.anchor.keys().collect();
    anchorlist.sort();

    let mut vector_matrix: Vec<Vec<Option<f64>>> = vec![vec![Some(NAN); sorted_read_sets.len()]; anchorlist.len()];
    let empty_vec: Vec<String> = Vec::new();
    for (anchor_index, anchor) in anchorlist.iter().enumerate() {
        let outgoinglist = graph.outgoing.get(*anchor).unwrap_or(&empty_vec);
        if outgoinglist.is_empty() {
            assert_eq!(*anchor, *anchorlist.last().unwrap());
        }

        let mut sorted_outgoinglist: Vec<&String> = outgoinglist.iter().collect();
        sorted_outgoinglist.sort();
        let d: HashMap<_, _> = sorted_outgoinglist.iter().enumerate().map(|(i, v)| (*v, i+1)).collect();

        for edge in outgoinglist {
            if let Some(reads) = graph.edges.get(edge).and_then(|e| e.get("reads")) {
                if let serde_json::Value::Array(reads_array) = reads {
                    for read in reads_array.iter() {
                        if let serde_json::Value::String(read_str) = read {
                            let read_index = *r.get(read_str).unwrap();
                            vector_matrix[anchor_index][read_index] = Some(*d.get(edge).unwrap() as f64);
                        }
                    }
                }
            }
        }
    }

    (vector_matrix, sorted_read_sets)
}

fn knn_impute(data: &Array2<f64>, k: usize) -> Array2<f64> {
    let mut imputed_data = data.clone();

    for ((row, col), value) in data.indexed_iter() {
        if value.is_nan() {
            let mut distances = Vec::new();
            for (other_row, other_value) in data.outer_iter().enumerate() {
                if other_row != row && !other_value[col].is_nan() {
                    let distance = count_non_equal_elements(&data.row(row).to_owned(), &data.row(other_row).to_owned());
                    distances.push((distance, other_value[col]));
                }
            }
            distances.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

            let sum: f64 = distances.iter().take(k).map(|&(_, val)| val).sum();
            imputed_data[(row, col)] = sum / k as f64;
        }
    }

    imputed_data
}

fn euclidean_distance(a: &Array1<f64>, b: &Array1<f64>) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| if x.is_nan() || y.is_nan() { 0.0 } else { (x - y).powi(2) })
        .sum::<f64>()
        .sqrt()
}

fn manhattan_distance(a: &Array1<f64>, b: &Array1<f64>) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| (x - y).abs()).sum::<f64>()
}

fn count_non_equal_elements(a: &Array1<f64>, b: &Array1<f64>) -> usize {
    // Count the number of elements that are not equal, ignoring NaNs
    a.iter()
     .zip(b.iter())
     .filter(|(&x, &y)| !x.is_nan() && !y.is_nan() && x != y)
     .count()
}



pub fn start(output: &PathBuf, graph_path: &PathBuf, k_nearest_neighbor: usize) {
    let graph = skydive::agg::GraphicalGenome::load_graph(graph_path.to_str().unwrap()).unwrap();
    let (vector_matrix, sorted_read_names) = construct_anchor_table(&graph);
    // println!("The answer is {:?} {:?}!", output, sorted_read_names);
    println!("columns:\n{:?}", sorted_read_names.len());
    let d = vector_matrix.len();
    println!("The dimensino of vector_matrix is {}", d);
    let vector_matrix_f64: Vec<Vec<f64>> = vector_matrix.iter()
        .map(|row| row.iter().map(|&x| x.unwrap_or(f64::NAN)).collect())
        .collect();

    let flat_data: Vec<f64> = vector_matrix_f64.iter().flatten().map(|&x| x as f64).collect();
    let rows = vector_matrix_f64.len();
    let cols = if !vector_matrix_f64.is_empty() { vector_matrix_f64[0].len() } else { 0 };
    let data = Array2::from_shape_vec((rows, cols), flat_data).unwrap();

    // Create the KNN imputer
    let data_t = data.t().to_owned();
    let imputed_data = knn_impute(&data_t, k_nearest_neighbor);
    println!("Input Data:\n{:?}", data);
    println!("Imputed Data:\n{:?}", imputed_data);

    //  find single sample matrix from the imputed data
    let sample = "HG002";
    let row_index: Vec<usize> = sorted_read_names.iter()
        .enumerate()
        .filter_map(|(index, readname)| {
            if readname.split('|').last().unwrap_or("") == sample{
                Some(index)
            }else{
                None
            }
        })
        .collect();
    let read_sets: Vec<&String> = row_index.iter()
        .map(|&i| &sorted_read_names[i])
        .collect();
    println!("Read_sets {:?}", read_sets);
    let vector_single_sample = imputed_data.select(Axis(0), &row_index);
    println!("Single sample Matrix:\n{:?}", vector_single_sample);

    let mut anchorlist: Vec<&String> = graph.anchor.keys().collect();
    anchorlist.sort();

    
}
   