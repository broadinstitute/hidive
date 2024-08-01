use std::path::PathBuf;

use skydive;
use gfa;
use handlegraph::{
    handle::{Direction, Handle, NodeId},
    handlegraph::HandleGraph,
    hashgraph::HashGraph,
};
use recgraph::api::*;
use bio::io::fasta::{Reader, Record};

pub fn start(output: &PathBuf, graph: &PathBuf, seq_paths: &PathBuf) {
    println!("The answer is {:?} {:?}!", output, graph);
    let parser = gfa::parser::GFAParser::new();
    let gfa: gfa::gfa::GFA<usize, ()> = parser.parse_file(graph).unwrap();
    let g: HashGraph = HashGraph::from_gfa(&gfa);

    let reader = Reader::from_file(seq_paths).unwrap();
    let all_reads: Vec<Record> = reader.records().map(|r| r.unwrap()).collect();

    for record in all_reads.clone() {
        let h = String::from_utf8_lossy(record.id().as_bytes()).to_string();
        let read_seq_upper = String::from_utf8(record.seq().to_ascii_uppercase()).expect("Invalid UTF-8 sequence");
        let score_matrix = create_score_matrix_i32(Some(0), Some(1), None);
        // let alignment = align_global_no_gap(&read_seq_upper, &g, h, score_matrix, Some(0.1));
          
    }
    
}