// Import necessary standard library modules
use std::collections::HashSet;
use std::collections::HashMap;
use std::hash::Hash;
use std::path::PathBuf;
use serde_json::Value;
// Import the Absolutize trait to convert relative paths to absolute paths
use bio::io::fasta::{Reader, Record};

use minimap2::{Aligner, Preset};

extern crate ndarray;
extern crate gurobi;
use gurobi::*;

use rust_wfa2;

use ndarray::{Array,Array1, Array2, array};
use ndarray::Axis;
use skydive::agg::GraphicalGenome;
use std::f64::NAN;

// Import the Url type to work with URLs
use url::Url;

// Import the skydive module, which contains the necessary functions for staging data
use skydive;

<<<<<<< HEAD
<<<<<<< HEAD
pub fn start(output: &PathBuf, graph: &PathBuf) {
    let graph = skydive::agg::GraphicalGenome::load_graph(output.to_str().unwrap()).unwrap();

    println!("The answer is {:?} {:?}!", output, graph);
=======
=======

>>>>>>> 5a0bd5b (add read files)
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
>>>>>>> 7d233e6 (add single sample graph extraction based on knn imputation)
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
        let d: HashMap<_, _> = sorted_outgoinglist.iter().enumerate().map(|(i, v)| (*v, i+1)).collect(); // index + 1 in this function

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

pub fn find_targetseq_in_reads(read_seq: &str, source_kmer:&str, sink_kmer: &str) -> String{
    let source_rev = skydive::agg::reverse_complement(source_kmer);
    let sink_rev = skydive::agg::reverse_complement(sink_kmer);

    let spos = read_seq.find(&source_kmer);
    let epos = read_seq.find(&sink_kmer);
    let rspos = read_seq.find(&source_rev);
    let repos = read_seq.find(&sink_rev);

    match (spos, epos, rspos, repos) {
        (Some(spos), Some(epos), _, _) => {
            let end_index = epos + sink_kmer.chars().count();
            read_seq[spos..end_index].to_string()
        },
        (_, _, Some(rspos), Some(repos)) => {
            let end_index = rspos + sink_kmer.chars().count();
            skydive::agg::reverse_complement(&read_seq[repos..end_index])
        },
        (Some(spos), None, _, _) => {
            read_seq[spos..].to_string()
        },
        (None, Some(epos), _, _) => {
            read_seq[..epos].to_string()
        },
        (None, None, Some(rspos), None) => {
            skydive::agg::reverse_complement(&read_seq[..rspos])
        },
        (None, None, None, Some(repos)) => {
            skydive::agg::reverse_complement(&read_seq[repos..])
        },
        _ => "".to_string(),
    }

}

// pub fn alignment(pattern:&str, text:&str) -> i32 {
//     let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
//     let pat_len = pattern.as_bytes().len();
//     let text_len = text.as_bytes().len();
//     let mut aligner = Aligner::with_capacity(pattern.len(), text.len(), -5, -1, &score);
//     let alignment = aligner.global(pattern.as_bytes(), text.as_bytes());

//     let score = alignment.score;
//     score
// }

pub fn calculate_edit_distance(read_path:&PathBuf, sample:&str, single_sample_graph: GraphicalGenome, source_node_name:&str, source: &String, sink_node_name: &str, sink:&String) -> HashMap<usize, HashMap<String, String>>{
        // import read files
    let reader = Reader::from_file(read_path).unwrap();
    let all_reads: Vec<Record> = reader.records().map(|r| r.unwrap()).collect();
    

    // pairwise alignment
    let single_sample_readsets = find_all_read(&single_sample_graph);
    let path_list = skydive::agg::FindAllPathBetweenAnchors::new(&single_sample_graph, source_node_name, sink_node_name, single_sample_readsets.clone());
    println!( " candidate path number : {:?}", path_list.subpath.len());

    let mut data_info: HashMap<usize, HashMap<String, String>> = HashMap::new();
    let mut index:usize = 0;

    for (p, r) in path_list.subpath.iter() {
        let pseq = skydive::agg::reconstruct_path_seq(&single_sample_graph, p);
        // println!("path lenth:{:?}", pseq.len());
        // minimap2 let aligner = Aligner::builder().asm5().with_seq(pseq.as_bytes()).expect("Unable to build index");
        
        
        for record in all_reads.clone() {
            let h = String::from_utf8_lossy(record.id().as_bytes()).to_string();
            if single_sample_readsets.contains(&h) {
                let read_seq_upper = String::from_utf8(record.seq().to_ascii_uppercase()).expect("Invalid UTF-8 sequence");

                let target = find_targetseq_in_reads(&read_seq_upper, source, sink);
                // println!("target lenth:{:?}", target.len());

                
                if !target.is_empty() {
                    // minimap2 let hits =aligner.map(target.as_bytes(), false, false, None, None);
                    // minimap2 let score = hits.unwrap(); 

                    // rust_wfa2
                    let alignment_scope = rust_wfa2::aligner::AlignmentScope::Alignment;
                    let memory_model = rust_wfa2::aligner::MemoryModel::MemoryUltraLow;
                    let mut aligner = rust_wfa2::aligner::WFAlignerGapAffine::new(1, 5, 2, alignment_scope, memory_model);
                    let status = aligner.align_end_to_end(target.as_bytes(), pseq.as_bytes());
                    let score = aligner.score();

                    // bio::alignment too slow
                    // let score = alignment(&target, &pseq);


                    data_info.insert(index, {
                        let mut sub_map = HashMap::new();
                        sub_map.insert("sample".to_string(), sample.to_string());
                        sub_map.insert("read".to_string(), h.to_string());
                        sub_map.insert("path".to_string(), p.join(">"));
                        sub_map.insert("cost".to_string(), (-score).to_string());
                        sub_map
                    });

                    index += 1;

                    // println!("{:?}, {:?}", h, score);
                }else {
                    continue;
                }
            }

        }

    };

    data_info
}

pub fn start(output: &PathBuf, graph_path: &PathBuf, read_path:&PathBuf, k_nearest_neighbor: usize) {
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
    let read_sets: Vec<String> = row_index.iter()
        .map(|&i| sorted_read_names[i].clone())
        .collect();
    println!("Read_sets {:?}", read_sets);
    let vector_single_sample = imputed_data.select(Axis(0), &row_index);
    println!("Single sample Matrix:\n{:?}", vector_single_sample);
    
    // get anchorlist and start end anchor sequences
    let mut anchorlist: Vec<String> = graph.anchor.keys().cloned().collect();
    anchorlist.sort();
    println!("{:?}", anchorlist.last());
    let source_node_name = anchorlist.first().unwrap();
    let source = if let Some(anchor_data) = graph.anchor.get(source_node_name) {
        if let Some(serde_json::Value::String(seq)) = anchor_data.get("seq") {
            seq 
        } else {
            panic!("there is no sequences")
        }
    } else {
        panic!("source node name not found")
    };
    let source_rev = skydive::agg::reverse_complement(&source);
    let sink_node_name = anchorlist.last().unwrap();
    let sink = if let Some(anchor_data) = graph.anchor.get(sink_node_name) {
        if let Some(serde_json::Value::String(seq)) = anchor_data.get("seq") {
            seq
        } else {
            panic!("there is no anchor sequences")
        }
    }else{
        panic!("Sink node name not found")
    };
    let sink_rev = skydive::agg::reverse_complement(&sink);

    // Extract single sample graph
    let single_sample_graph = skydive::agg::GraphicalGenome::extract_single_sample_graph(&graph, &vector_single_sample, anchorlist.clone(), read_sets, sample).unwrap();
    println!("Single sample Graph:\n{:?}", graph.incoming);
 
    // pairwise alignment between reads and candidate paths from single_sample graph, 
    // construct HashMap for Ryan's optimizer
    let data_info = calculate_edit_distance(read_path, sample, single_sample_graph, source_node_name, source, sink_node_name, sink);
    
    //  Ryan's optimizer translated using gurobi
    let mut unique_samples = HashSet::new();
    let mut unique_paths = HashSet::new();
    let mut unique_reads = HashSet::new();

    for sub_map in data_info.values() {
        if let Some(sample_name) = sub_map.get("sample") {
            unique_samples.insert(sample_name.clone());
        }
        if let Some(path_name) = sub_map.get("path") {
            unique_paths.insert(path_name.clone());
        }
        if let Some(read_name) = sub_map.get("read") {
            unique_reads.insert(read_name.clone());
        }

    }

    let env = Env::new("logfile.log").unwrap();
    let mut model = env.new_model("hap").unwrap();
    

    
    





    

}
   