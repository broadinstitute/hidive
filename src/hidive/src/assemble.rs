// Import necessary standard library modules
use std::collections::HashSet;
use std::collections::HashMap;
use std::hash::Hash;

use std::iter::SkipWhile;
use std::path::PathBuf;
use rayon::string;


use serde_json::Value;
// Import the Absolutize trait to convert relative paths to absolute paths
use bio::io::fasta::{Reader, Record};

// use minimap2::{Aligner, Preset};


extern crate ndarray;

use rust_wfa2;
// use russcip::prelude::*;

// use good_lp::{variables, variable, scip, constraint, SolverModel, Solution};


use grb::prelude::*;

use ndarray::{Array,Array1, Array2, array};
use ndarray::Axis;
use skydive::agg::GraphicalGenome;
use std::f64::NAN;

// Import the Url type to work with URLs
use url::Url;
// import spoa
extern crate rust_spoa;
use rust_spoa::poa_consensus;

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

pub fn calculate_edit_distance(read_path:&PathBuf, sample:&str, single_sample_graph: &GraphicalGenome, source_node_name:&str, source: &String, sink_node_name: &str, sink:&String) -> HashMap<usize, HashMap<String, String>>{
        // import read files
    let reader = Reader::from_file(read_path).unwrap();
    let all_reads: Vec<Record> = reader.records().map(|r| r.unwrap()).collect();
    

    // pairwise alignment
    let single_sample_readsets = find_all_read(single_sample_graph);
    let path_list = skydive::agg::FindAllPathBetweenAnchors::new(single_sample_graph, source_node_name, sink_node_name, single_sample_readsets.clone());
    println!( " candidate path number : {:?}", path_list.subpath.len());

    let mut data_info: HashMap<usize, HashMap<String, String>> = HashMap::new();
    let mut index:usize = 0;

    for (p, r) in path_list.subpath.iter() {
        let pseq = skydive::agg::reconstruct_path_seq(single_sample_graph, p);
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

pub fn mip_optimization (data_info: HashMap<usize, HashMap<String, String>>) -> (HashMap<String, f64>, HashMap<String, f64>, HashMap<String, i32>) {
    let mut unique_samples = HashSet::new();
    let mut unique_paths = HashMap::new();
    let mut unique_reads = HashSet::new();
    let mut ind = 0;
    for sub_map in data_info.values() {
        if let Some(sample_name) = sub_map.get("sample") {
            unique_samples.insert(sample_name.clone());
        }
        if let Some(path_name) = sub_map.get("path") {
            unique_paths.insert(path_name.clone(), ind);
            ind += 1;
        }
        if let Some(read_name) = sub_map.get("read") {
            unique_reads.insert(read_name.clone());
        }

    }

    // create model
    let mut model = Model::new("haplotype").unwrap();
    
    let mut var_haps = HashMap::new();
    let mut var_sample_hap = HashMap::new();
    let mut var_flow = HashMap::new();

    for (hap_name, hap) in unique_paths.iter(){
        var_haps.insert(hap.clone(), add_binvar!(model, name:&hap.to_string()).unwrap());
    }

    let mut total_cost_expr = grb::expr::LinExpr::new();
    let mut flow_out = HashMap::new();
    let mut connected_haps_per_sample = HashMap::new();

    for (index, submap) in data_info.iter(){
        let s = submap.get("sample").unwrap();
        let r = submap.get("read").unwrap();
        let p = submap.get("path").unwrap();
        let c:f64 = submap.get("cost").unwrap().parse().expect("Not a valid integer for cost");

        let r_p_identifier = format!("{}_{}", r, unique_paths.get(p).unwrap());
        let s_p_identifier = format!("{}_{}", s, unique_paths.get(p).unwrap());
        var_flow.insert(r_p_identifier.clone(), add_var!(model, Continuous, name: &r_p_identifier, obj: 0.0, bounds: 0..1).unwrap());


        if !var_sample_hap.contains_key(&s_p_identifier){
            let inserted_value = add_binvar!(model, name:&s_p_identifier).unwrap();
            var_sample_hap.insert(s_p_identifier.clone(), inserted_value);
            connected_haps_per_sample.entry(s.clone()).or_insert_with(Vec::new).push(inserted_value);
        }



        // total flows:
        let flow_var = var_flow.get(&r_p_identifier).unwrap();
        let sample_hap_var = var_sample_hap.get(&s_p_identifier).unwrap();
        let hap_index = unique_paths.get(p).unwrap();
        let hap_var = var_haps.get(hap_index).unwrap();
        total_cost_expr.add_term( c, *flow_var);
        flow_out.entry(r).or_insert(grb::expr::LinExpr::new()).add_term(1.0, *flow_var) ;  

        // add constraints
        let name1 = (2*index).to_string();
        model.add_constr(&name1, c!(*flow_var - *sample_hap_var <= 0)).unwrap();
        let name2 = (2*index + 1).to_string();
        model.add_constr(&name2, c!( *sample_hap_var - *hap_var <= 0)).unwrap();      
    }

    for read in unique_reads.iter(){
        let mut flow_out_var = flow_out.get(read).unwrap();
        model.add_constr(read, c!(flow_out_var.clone() == 1));
    }

    for sample in unique_samples.iter() {
        let sample_haps = connected_haps_per_sample.get(sample).unwrap();
        let mut sum_expr: grb::expr::LinExpr = grb::expr::LinExpr::new();
        for sh in sample_haps.iter(){
            sum_expr.add_term(1.0, *sh);
        }
        model.add_constr(&format!("constraint_for_{}", sample), c!(sum_expr <= 2)).unwrap();
    }

    let mut var_num_haps = add_var!(model, Continuous, name: "var_num_haps", obj: 0.0).unwrap();
    let mut var_num_haps_expr = grb::expr::LinExpr::new();
    for (path_index, var) in var_haps.iter(){
        var_num_haps_expr.add_term(1.0, *var);
    }
    var_num_haps_expr.add_term(-1.0, var_num_haps);
    model.add_constr("constraints_of_var_num_haps", c!(var_num_haps_expr == 0));

    let mut var_total_cost = add_var!(model, Continuous, name: "var_total_cost", obj: 0.0).unwrap();
    total_cost_expr.add_term(-1.0, var_total_cost);
    model.add_constr("total_cost", c!(total_cost_expr == 0));


    // solving models
    // minimize total cost
    model.set_objective(var_total_cost, Minimize);
    model.optimize();
    assert_eq!(model.status().unwrap(), Status::Optimal);

    let least_cost = model.get_obj_attr(grb::attr::X, &var_total_cost).unwrap();
    println!("least_cost is  {:?}", least_cost);

    // minimize most haplotypes that minimizes total cost
    let aux = model.add_constr("total_cost_equals_leastcost", c!(var_total_cost == least_cost)).unwrap();
    model.set_objective(var_num_haps, Minimize);
    model.optimize();   
    let most_haps = model.get_obj_attr(grb::attr::X, &var_num_haps).unwrap();
    println!("most haplotype is {:?}", most_haps);

    // minimize haplotypes
    model.remove(aux);
    model.optimize();
    let least_haps = model.get_obj_attr(grb::attr::X, &var_num_haps).unwrap();
    println!("least haplotype is {:?}", least_haps);

    // find most_cost that minimizes haplotypes
    let aux1 = model.add_constr("least_haplotype", c!(var_num_haps == least_haps)).unwrap();
    model.set_objective(var_total_cost, Minimize);
    model.optimize();
    let most_cost = model.get_obj_attr(grb::attr::X, &var_total_cost).unwrap();
    model.remove(aux1);
    println!("most_cost = {:?}", most_cost);

    // final solution
    let range_cost = most_cost - least_cost;
    let range_haps = most_haps - least_haps;

    if range_cost == 0.0{
        assert_eq!( range_haps, 0.0);
        let optimal_value = 0;
        let optimal_cost = least_cost;
        let optimal_haps = least_haps;
    }else{
        let mut final_expr = grb::expr::QuadExpr::new();
        let coeff1 = 1.0/range_haps;
        let coeff2 = 1.0 / range_cost;
        let mut var1 = add_var!(model, Continuous, name: "var_num_hap_minus_leasthap", obj: 0.0).unwrap();
        model.add_constr("var1", c!(var_num_haps - least_haps - var1 == 0));
        // let var1 = var_num_haps - least_haps;
        let mut var2 = add_var!(model, Continuous, name: "var_total_cost_minus_leastcost", obj: 0.0).unwrap();
        model.add_constr("var2", c!(var_total_cost - least_cost - var2 == 0));
        // let var2 = var_total_cost - least_cost;
        final_expr.add_qterm(coeff1, var1.clone(), var1.clone());
        final_expr.add_qterm(coeff2, var2.clone(), var2.clone());

        model.set_objective(final_expr, Minimize);
        model.optimize();
        let optimal_value = model.get_attr(grb::attr::ObjVal).unwrap();
        let optimal_cost = model.get_obj_attr(grb::attr::X, &var_total_cost).unwrap();
        let optimal_haps = model.get_obj_attr(grb::attr::X, &var_num_haps).unwrap();
    }

    let mut var_sample_hap_values = HashMap::new();
    for (identifier, var) in var_sample_hap.iter() {
        let mut var_value = model.get_obj_attr(grb::attr::X, var).unwrap();
        var_sample_hap_values.insert(identifier.clone(), var_value); 
    }
    let mut var_flow_values = HashMap::new();
    for (ridentifier, var) in var_flow.iter() {
        let mut v = model.get_obj_attr(grb::attr::X, var).unwrap();
        var_flow_values.insert(ridentifier.clone(), v);
    }

    (var_sample_hap_values, var_flow_values, unique_paths)
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
    // println!("Single sample Matrix:\n{:?}", vector_single_sample);
    
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
    // println!("Single sample Graph:\n{:?}", graph.incoming);
 
    // pairwise alignment between reads and candidate paths from single_sample graph, 
    // construct HashMap for Ryan's optimizer
    let data_info = calculate_edit_distance(read_path, sample, &single_sample_graph, source_node_name, source, sink_node_name, sink);
    
    //  Ryan's optimizer translated using gurobi

    let (var_sample_hap_values, var_flow_values, unique_paths) = mip_optimization(data_info);
    

    // get path dictionary from index to path entities  
    let mut path_dict = HashMap::new();
    for (p, ind) in unique_paths.iter(){
        path_dict.insert(ind.to_string(), p);
    }
    // get path read pair
    let mut path_read_pair = HashMap::new();
    for (read_path_id, var_id) in var_flow_values.iter(){
        let path = read_path_id.split("_").last().unwrap();
        if *var_id == 1.0{
            path_read_pair.entry(path).or_insert_with(Vec::new).push(read_path_id.clone());
        }
    }
    println!("{:?}", path_read_pair);

    // get read sequences
    let reader = Reader::from_file(read_path).unwrap();
    let all_reads: Vec<Record> = reader.records().map(|r| r.unwrap()).collect();

    // calculate consensus
    for (sample_hap_id, var) in var_sample_hap_values.iter(){
        
        let hap_id = sample_hap_id.split("_").last().unwrap();
        let haplotype = path_dict.get(hap_id).expect("not in the haplotype dictionary");
        let path_absolute:Vec<String> = haplotype.split(">").map(|s| s.to_string()).collect();
        let mut input_seq = Vec::new();

        if *var == 1.0{
            let path_seq = skydive::agg::reconstruct_path_seq(&single_sample_graph, &path_absolute);
            println!("{:?}, {:?}, {:?}", sample_hap_id, hap_id, path_seq.len());
            input_seq.push(path_seq);
            let readlist = path_read_pair.get(hap_id).expect("haplotype not found");
            for read in readlist.iter(){
                for record in all_reads.clone() {
                    let h = String::from_utf8_lossy(record.id().as_bytes()).to_string();
                        if read.contains(&h) {
                            let read_seq = String::from_utf8(record.seq().to_ascii_uppercase()).expect("Invalid UTF-8 sequence");
                            let target = find_targetseq_in_reads(&read_seq, source, sink);
                            input_seq.push(target);
                        }
                }                
            }

            let mut sequences = vec![];
            for seq in input_seq.iter_mut(){
                seq.push('\0');
                sequences.push((*seq).bytes().map(|x| {x as u8}).collect::<Vec<u8>>());
            }
            let mut consensus = poa_consensus(&sequences, 10000000, 1, 5, -4, -3, -1);
            let consensus = String::from_utf8(consensus).expect("bytes to string");
            println!("{:?}, {:?}", hap_id, consensus);
        }
        

    }

}
