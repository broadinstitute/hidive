// Import necessary standard library modules
use std::collections::HashSet;
use std::path::PathBuf;

use bio::io::fasta::{IndexedReader, Reader, Record};

use flate2::read::GzDecoder;
use serde_json::Value;

// Import the skydive module, which contains the necessary functions for building graphs
use skydive;

use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::sync::Mutex;

pub fn get_reference_kmer_profile(ref_hla: &str, k: usize) -> Vec<String> {
    let mut ref_kmer_profile: HashMap<String, i32> = HashMap::new();
    let r = ref_hla.len() as isize - k as isize + 1;

    for i in 0..r as usize {
        let kmer: String = ref_hla.chars().skip(i).take(k).collect();
        *ref_kmer_profile.entry(kmer).or_insert(0) += 1;
    }

    ref_kmer_profile
        .into_iter()
        .filter(|&(_, count)| count == 1)
        .map(|(kmer, _)| kmer)
        .collect()
}

pub fn reverse_complement(kmer: &str) -> String {
    kmer.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'N' => 'N',
            _ => panic!("Unexpected character: {}", c),
        })
        .collect()
}

pub fn find_sequences_between_sanchor_eanchor(
    reads: Vec<Record>,
    ref_hla: String,
    stem: &String,
) -> (HashSet<String>, HashMap<String, String>) {
    let mut hla_samples = HashSet::new();
    let mut hla_seq = HashMap::new();

    for read in reads.iter() {
        let h = String::from_utf8_lossy(read.id().as_bytes()).to_string();
        let seq_upper =
            String::from_utf8(read.seq().to_ascii_uppercase()).expect("Invalid UTF-8 sequence");

        let sample_identifier = h.split("|").last().unwrap_or_default().to_string();
        hla_samples.insert(sample_identifier.clone());
        hla_seq.insert(h.clone(), seq_upper);
    }

    hla_samples.insert(stem.to_owned());
    hla_seq.insert(stem.to_owned(), ref_hla);

    (hla_samples, hla_seq)
}

pub fn map_reference_unique_kmers_to_seq(
    unique_kmerlist: Vec<String>,
    hla_seq: &HashMap<String, String>,
    k: usize,
) -> (
    Mutex<HashMap<String, HashSet<String>>>,
    Mutex<HashMap<String, HashMap<String, Vec<usize>>>>,
) {
    let sample_dict: Mutex<HashMap<String, HashSet<String>>> = Mutex::new(HashMap::new());
    let position_dict: Mutex<HashMap<String, HashMap<String, Vec<usize>>>> =
        Mutex::new(HashMap::new());

    unique_kmerlist.iter().for_each(|kmer| {
        let kmer_upper = kmer.to_uppercase();
        let rev = reverse_complement(&kmer_upper);
        {
            let mut sample_dict = sample_dict.lock().unwrap();
            sample_dict.entry(kmer_upper.clone()).or_default();
            sample_dict.entry(rev.clone()).or_default();
        }
        {
            let mut position_dict = position_dict.lock().unwrap();
            position_dict.entry(kmer_upper.clone()).or_default();
            position_dict.entry(rev.clone()).or_default();
        }
    });

    hla_seq.par_iter().for_each(|(readname, seq)| {
        let seq_len = seq.len();

        if seq_len >= k {
            for i in 0..=seq_len - k {
                let kmer: String = seq[i..i + k].to_uppercase();
                // Direct slice instead of skip and take
                {
                    let mut sample_dict = sample_dict.lock().unwrap();
                    if let Some(samples) = sample_dict.get_mut(&kmer) {
                        samples.insert(readname.split('|').last().unwrap_or_default().to_string());
                    }
                }
                {
                    let mut position_dict = position_dict.lock().unwrap();
                    if let Some(positions) = position_dict.get_mut(&kmer) {
                        positions.entry(readname.clone()).or_default().push(i);
                    }
                }
            }
        } else {
            println!("Panic! Read sequences are too short!");
        }
    });

    (sample_dict, position_dict)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnchorInfo {
    seq: String,
    pos: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EdgeInfo {
    seq: String,
    reads: Vec<String>,
    samples: HashSet<String>,
    src: String,
    dst: String,
}

impl AnchorInfo {
    // Getter method for `seq`
    pub fn get_seq(&self) -> &String {
        &self.seq
    }
}

pub fn get_anchor_information(
    sample_dict: &HashMap<String, HashSet<String>>,
    hla_samples: &HashSet<String>,
    position_dict: &HashMap<String, HashMap<String, Vec<usize>>>,
    stem: &String,
) -> Vec<String> {
    let mut anchorlist = Vec::new();
    for (kmer, samplelist) in sample_dict.iter() {
        if samplelist.len() == hla_samples.len() {
            if let Some(positions) = position_dict.get(kmer) {
                if positions.len() > 1 && positions.get(stem).map_or(false, |v| v.len() == 1) {
                    anchorlist.push(kmer.clone());
                }
            }
        }
    }
    anchorlist
}

pub fn get_anchors(
    anchorlist: &Vec<String>,
    position_dict: &HashMap<String, HashMap<String, Vec<usize>>>,
    k: usize,
    stem: &String,
) -> HashMap<String, AnchorInfo> {
    let mut anchor_info = HashMap::new();
    for kmer in anchorlist.iter() {
        if let Some(positions) = position_dict.get(kmer) {
            if let Some(pos) = positions.get(stem).and_then(|v| v.first()) {
                let anchor_name = format!("A{:06}", pos / k + 1);
                anchor_info.insert(
                    anchor_name,
                    AnchorInfo {
                        seq: kmer.clone(),
                        pos: *pos,
                    },
                );
            }
        }
    }

    anchor_info
}

pub fn get_final_anchor(
    anchor_info: &HashMap<String, AnchorInfo>,
    k: usize,
) -> HashMap<String, &AnchorInfo> {
    let mut final_anchor = HashMap::new();
    let mut anchornames: Vec<&String> = anchor_info.keys().collect();
    anchornames.sort();

    let mut anchor_unadjacent_list = Vec::new();
    let mut index = 0;

    while index < anchornames.len() - 1 {
        let sanchor = anchornames[index];
        let mut danchor = sanchor;
        let mut next_index = index;
        for (i, next_anchor) in anchornames[index + 1..].iter().enumerate() {
            if anchor_info.get(*next_anchor).unwrap().pos
                > anchor_info.get(sanchor).unwrap().pos + k + 1
            {
                break;
            }
            next_index = index + i + 1;
            danchor = next_anchor;
        }
        index = next_index + 1;

        // index = anchornames.iter().enumerate().find(|&(_, ref e)| *e == &danchor).map(|(index, _)| index).unwrap_or(0);
        // println!("{}, {}",index, danchor);
        anchor_unadjacent_list.push(sanchor.clone());
        anchor_unadjacent_list.push(danchor.clone());
    }

    anchor_unadjacent_list.sort();
    anchor_unadjacent_list.dedup();

    for anchor_name in &anchor_unadjacent_list {
        if let Some(anchor) = anchor_info.get(anchor_name) {
            final_anchor.insert(anchor_name.clone(), anchor);
        }
    }

    final_anchor
}

pub fn mapping_info(
    final_anchor: &HashMap<String, &AnchorInfo>,
    contig: String,
    k: usize,
) -> (HashMap<String, Option<usize>>, HashMap<String, Vec<usize>>) {
    let anchor_list: Vec<String> = final_anchor.keys().cloned().collect();
    let seqlist: Vec<String> = anchor_list
        .iter()
        .filter_map(|anchor| final_anchor.get(anchor).map(|info| info.get_seq().clone()))
        .collect();

    let mut position_dict: HashMap<String, Vec<usize>> = HashMap::new();
    for anchor_seq in seqlist.into_iter() {
        let anchor_rev = reverse_complement(&anchor_seq);
        position_dict
            .entry(anchor_seq.clone())
            .or_insert_with(|| Vec::new());
        position_dict
            .entry(anchor_rev.clone())
            .or_insert_with(|| Vec::new());
    }

    for i in 0..contig.len() - k + 1 {
        let kmer: String = contig.chars().skip(i).take(k).collect();
        if let Some(positions) = position_dict.get_mut(&kmer) {
            positions.push(i);
        }
    }

    let mut a = HashMap::new();
    let mut svs = HashMap::new();

    for (anchor, info) in final_anchor.iter() {
        let anchor_seq = info.seq.clone();
        let positionlist = position_dict.get(&anchor_seq).unwrap_or(&vec![]).clone();
        if positionlist.len() == 1 {
            a.insert(anchor.clone(), positionlist.first().copied());
        } else {
            svs.insert(anchor.clone(), positionlist.clone());
        }
    }

    (a, svs)
}
pub fn construct_edges(
    src_pos: usize,
    dst_pos: usize,
    k: usize,
    contig: String,
    contigname: String,
    sample: String,
    Anchorseq: &HashMap<String, String>,
) -> EdgeInfo {
    let src_seq: String;
    let mut pr = false;
    let mut src = "".to_string();

    if src_pos == 0 {
        src = "SOURCE".to_string();
        src_seq = "".to_string();
        pr = true;
    } else {
        src_seq = contig.chars().skip(src_pos).take(k).collect();
        src = match Anchorseq.get(&src_seq) {
            Some(value) => value.clone(),
            None => {
                let reversed_src_seq = reverse_complement(&src_seq);
                pr = true;
                Anchorseq
                    .get(&reversed_src_seq)
                    .unwrap_or(&"".to_string())
                    .clone()
            }
        };
    }

    let mut sr = false;
    let mut dst = "".to_string();
    let dst_seq;
    if dst_pos == contig.len() {
        dst = "SINK".to_string();
        dst_seq = "".to_string();
        sr = true;
    } else {
        dst_seq = contig.chars().skip(dst_pos).take(k).collect::<String>();
        dst = match Anchorseq.get(&dst_seq) {
            Some(value) => value.clone(),
            None => {
                let reversed_dst_seq = reverse_complement(&dst_seq);
                sr = true;
                Anchorseq
                    .get(&reversed_dst_seq)
                    .unwrap_or(&"".to_string())
                    .clone()
            }
        };
    }

    let mut edge_seq = if src_pos == 0 {
        contig.get(0..dst_pos).unwrap_or_default().to_string()
    } else {
        contig
            .get(src_pos + k..dst_pos)
            .unwrap_or_default()
            .to_string()
    };

    if pr && sr {
        edge_seq = reverse_complement(&edge_seq);
        let node = src.clone();
        src = dst.clone();
        dst = node.clone();
    }

    EdgeInfo {
        seq: edge_seq,
        src: src,
        dst: dst,
        reads: vec![contigname],
        samples: vec![sample].into_iter().collect(),
    }
}

pub fn create_edge_file(
    hla_seq: &HashMap<String, String>,
    final_anchor: &HashMap<String, &AnchorInfo>,
    k: usize,
) -> (HashMap<String, EdgeInfo>, HashMap<String, Vec<String>>) {
    let mut edge_info: HashMap<String, EdgeInfo> = HashMap::new();
    let mut outgoing: HashMap<String, Vec<String>> = HashMap::new();
    // let mut edgenum_perread = Vec::new();

    let anchorseq: HashMap<_, _> = final_anchor
        .iter()
        .map(|(anchor, info)| (info.seq.clone(), anchor.clone()))
        .collect();
    let mut contig_index = 0;

    for (contig_name, contig) in hla_seq.iter() {
        contig_index += 1;
        let sample_name = contig_name
            .split('|')
            .last()
            .unwrap_or_default()
            .to_string();
        let (a, svs) = mapping_info(&final_anchor, contig.to_string(), k);
        let mut splitposlist: Vec<_> = a.values().filter_map(|&x| x).collect();
        splitposlist.sort();
        let mut edgeindex = 0;
        let mut src_pos = 0;
        for &dst_pos in &splitposlist {
            let e = construct_edges(
                src_pos,
                dst_pos,
                k,
                contig.to_string(),
                contig_name.to_string(),
                sample_name.clone(),
                &anchorseq,
            );
            let src = &e.src;
            let edgelist = outgoing.entry(src.clone()).or_insert_with(|| Vec::new());
            if let Some(pos) = edgelist
                .iter()
                .position(|edge| edge_info[edge].dst == e.dst && edge_info[edge].seq == e.seq)
            {
                let edge = edge_info.get_mut(&edgelist[pos]).unwrap();
                edge.reads.extend(e.reads);
                edge.samples = edge.samples.union(&e.samples).cloned().collect();
            } else {
                let edgename = format!("E{:05}.{:04}", contig_index, edgeindex);
                edge_info.insert(edgename.clone(), e);
                edgelist.push(edgename);
                edgeindex += 1;
            }
            src_pos = dst_pos;
        }

        let dst_pos = contig.len();
        let e = construct_edges(
            src_pos,
            dst_pos,
            k,
            contig.to_string(),
            contig_name.to_string(),
            sample_name.clone(),
            &anchorseq,
        );
        let src = &e.src;
        let edgelist = outgoing.entry(src.clone()).or_insert_with(|| Vec::new());
        if let Some(pos) = edgelist
            .iter()
            .position(|edge| edge_info[edge].dst == e.dst && edge_info[edge].seq == e.seq)
        {
            let edge = edge_info.get_mut(&edgelist[pos]).unwrap();
            edge.reads.extend(e.reads);
            edge.samples = edge.samples.union(&e.samples).cloned().collect();
        } else {
            let edgename = format!("E{:05}.{:04}", contig_index, edgeindex);
            edge_info.insert(edgename.clone(), e);
            edgelist.push(edgename);
            edgeindex += 1;
        }

        // edgenum_perread.push(edgeindex);
    }

    (edge_info, outgoing)
}

pub fn write_gfa(
    final_anchor: &HashMap<String, AnchorInfo>,
    edge_info: &HashMap<String, EdgeInfo>,
    output_filename: &PathBuf,
) -> Result<(), Box<dyn Error>> {
    let mut file = File::create(output_filename)?;

    writeln!(file, "H\tVN:Z:1.0")?;

    let mut anchor_output = Vec::new();
    for (anchor, info) in final_anchor.iter() {
        let seq = &info.seq;
        let mut anchor_info_clone = HashMap::new();
        anchor_info_clone.insert("pos".to_string(), info.pos);
        // anchor_info_clone.seq = String::new();
        let json_string =
            serde_json::to_string(&anchor_info_clone).unwrap_or_else(|_| "{}".to_string());
        let formatted_string = format!("S\t{}\t{}\tPG:J:{}", anchor, seq, json_string);
        anchor_output.push(formatted_string);
        // writeln!(file,"S\t{}\t{}\tPG:J:{}\n", anchor, seq, json_string)?;
    }
    let mut edge_output = Vec::new();
    let mut link_output = Vec::new();

    for (edge, edge_data) in edge_info {
        let seq = &edge_data.seq;
        let src = &edge_data.src;
        let dst = &edge_data.dst;
        let mut edge_data_clone = HashMap::new();
        edge_data_clone.insert("reads".to_string(), &edge_data.reads);
        let sample_vec: Vec<String> = edge_data.samples.iter().cloned().collect();
        let src_vec = vec![edge_data.src.clone()];
        let dst_vec = vec![edge_data.dst.clone()];
        edge_data_clone.insert("samples".to_string(), &sample_vec);
        edge_data_clone.insert("src".to_string(), &src_vec);
        edge_data_clone.insert("dst".to_string(), &dst_vec);

        let json_string =
            serde_json::to_string(&edge_data_clone).unwrap_or_else(|_| "{}".to_string());
        let formatted_string = if !edge_data.reads.is_empty() {
            format!(
                "S\t{}\t{}\tPG:J:{}\tRC:i:{}",
                edge,
                seq,
                json_string,
                edge_data.reads.len()
            )
            // writeln!(file, "S\t{}\t{}\tPG:J:{}\tRC:i:{}\n", edge, seq, json_string, edge_data.reads.len())?;
        } else {
            format!("S\t{}\t{}", edge, seq)
            // writeln!(file, "S\t{}\t{}\n", edge, seq)?;
        };
        edge_output.push(formatted_string);

        link_output.push(format!("L\t{}\t+\t{}\t+\t0M", src, edge));
        link_output.push(format!("L\t{}\t+\t{}\t+\t0M", edge, dst));
        // writeln!(file, "L\t{}\t+\t{}\t+\t0M", src, edge)?;
        // writeln!(file, "L\t{}\t+\t{}\t+\t0M", edge, dst)?;
    }

    for s in anchor_output {
        writeln!(file, "{}", s)?;
    }

    for s in edge_output {
        writeln!(file, "{}", s)?;
    }

    for l in link_output {
        writeln!(file, "{}", l)?;
    }

    Ok(())
}

pub fn filter_undersupported_edges(
    edge_info: &HashMap<String, EdgeInfo>,
    reference: &String,
    threshold: i32,
) -> HashMap<String, EdgeInfo> {
    let mut filtered_edges = edge_info.clone(); // Clone the input HashMap to create an owned version

    for (edge, data) in edge_info.iter() {
        if data.reads.contains(&reference.to_string()) || data.reads.len() >= threshold as usize {
            continue;
        }
        filtered_edges.remove(edge); // Remove the edge from the cloned HashMap
    }

    filtered_edges // Return the owned HashMap
}

pub struct GraphicalGenome {
    pub anchor: HashMap<String, Value>,
    pub edges: HashMap<String, Value>,
    pub outgoing: HashMap<String, Vec<String>>,
    pub incoming: HashMap<String, Vec<String>>,
}

impl GraphicalGenome {
    pub fn load_graph(filename: &str) -> io::Result<GraphicalGenome> {
        let file = File::open(filename)?;
        let reader: Box<dyn BufRead> = if filename.ends_with(".gz") {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };
        let mut anchor_dict = HashMap::new();
        let mut edge_dict = HashMap::new();
        let mut outgoing: HashMap<String, Vec<String>> = HashMap::new();
        let mut incoming: HashMap<String, Vec<String>> = HashMap::new();
        for line in reader.lines() {
            let line = line?.trim_end().to_string();
            if line.starts_with('S') {
                let itemlist: Vec<&str> = line.split('\t').collect();
                let name = itemlist[1];
                let seq = itemlist[2];
                let annotation = &itemlist[3][5..];

                let json_value: Value = serde_json::from_str(annotation).unwrap();

                let mut value = json_value;

                value["seq"] = Value::String(seq.to_string());

                if name.starts_with('A') {
                    anchor_dict.insert(name.to_string(), value);
                } else if name.starts_with('E') {
                    edge_dict.insert(name.to_string(), value);
                }
            } else if line.starts_with('L') {
                let itemlist: Vec<&str> = line.split('\t').collect();
                let src = itemlist[1];
                let dst = itemlist[3];
                outgoing
                    .entry(src.to_string())
                    .or_default()
                    .push(dst.to_string());
                incoming
                    .entry(dst.to_string())
                    .or_default()
                    .push(src.to_string());
            }
        }
        Ok(GraphicalGenome {
            anchor: anchor_dict,
            edges: edge_dict,
            outgoing: outgoing,
            incoming: incoming,
        })
    }
}

pub struct FindAllPathBetweenAnchors {
    pub subpath: Vec<(Vec<String>, HashSet<String>)>,
}

impl FindAllPathBetweenAnchors {
    pub fn new(graph: &GraphicalGenome, start: &str, end: &str, read_sets: HashSet<String>) -> Self {
        let mut finder = FindAllPathBetweenAnchors {
            subpath: Vec::new(),
        };
        finder.find_path(graph, start, end, Vec::new(), 0, read_sets);
        finder
    }

    pub fn find_path(&mut self, g: &GraphicalGenome, start: &str, end: &str, mut sofar: Vec<String>, depth: usize, readset: HashSet<String>) {
        if start == end {
            let mut sofar1 = sofar.clone();
            sofar1.push(end.to_string());
            if !readset.is_empty() {
                self.subpath.push((sofar1, readset));
            }
            return;
        }

        if readset.is_empty() || start == "SINK" {
            return;
        }

        if !g.outgoing.contains_key(start) {
            return;
        }

        let depth1 = depth + 1;

        if let Some(outgoing) = g.outgoing.get(start) {
            
            for dst in outgoing {
                let mut readset1 = readset.clone();
                if dst.starts_with("E") {
                    // let edge_reads: HashSet<String> = g.edges.get(dst)
                    //                                         .and_then(|edge| edge.get("reads").and_then(|r| r.as_array()))
                    //                                         .map(|reads| reads.iter().filter_map(|read| read.as_str().map(String::from)).collect())
                    //                                         .unwrap_or_default();

                    // let readset1: HashSet<String> = readset.intersection(&edge_reads).cloned().collect();
                    if let Some(edge_reads) = g.edges.get(dst).and_then(|e| e.get("reads").and_then(|r| r.as_array())) {
                        readset1.retain(|read| edge_reads.iter().any(|r| r.as_str() == Some(read)));
                    }
                    
                }
                let mut sofar1 = sofar.clone();
                sofar1.push(start.to_string());
                self.find_path(g, dst, end, sofar1, depth1, readset1);
            }
        }
    }
}

// Series parallele graph
pub fn reconstruct_path_seq(graph: &GraphicalGenome, path: &[String]) -> String {
    let mut seq = String::new();
    for item in path {
        if item.starts_with('A') {
            if let Some(anchor) = graph.anchor.get(item) {
                seq += &anchor["seq"].as_str().unwrap_or_default(); // Assuming `anchor` is a HashMap and "seq" is a key
                // println!("{:?}", anchor["seq"].as_str().unwrap_or_default());
            }
        } else if item.starts_with("E") {
            if let Some(edge) = graph.edges.get(item) {
                seq += &edge["seq"].as_str().unwrap_or_default(); // Assuming `edges` is a HashMap and "seq" is a key
            }
        }
    }
    seq
}

pub fn find_all_reads(graph: &GraphicalGenome) -> HashSet<String> {
    let mut read_sets = HashSet::new();
    for (item, edge) in &graph.edges {
        if let Some(readlist) = edge.get("reads").and_then(|r| r.as_array()) {
            for read in readlist {
                if let Some(read_str) = read.as_str() {
                    read_sets.insert(read_str.to_string());
                }
            }
        }
    }
    read_sets
}
pub struct GetSeriesParallelGraph {
    pub nodelist: Vec<String>,
    pub anchor: HashMap<String, Value>, // Assuming Value is a struct or type that represents the anchor data
    pub edges: HashMap<String, Value>,  // Assuming Value is a struct or type that represents the edge data
    pub outgoing: HashMap<String, Vec<String>>,
    pub incoming: HashMap<String, Vec<String>>,
}

impl GetSeriesParallelGraph {
    pub fn new(graph: &GraphicalGenome) -> GraphicalGenome {
        let nodelist = Self::series_parallel_graph_nodelist(graph);
        // println!("{:?}", nodelist);
        // let (anchor, edges, outgoing, incoming) = Self::series_parallel_graph(&nodelist, graph);
        let (anchor, edges, outgoing, incoming) = Self::series_parallel_graph(&nodelist, graph);
        GraphicalGenome {
            anchor,
            edges,
            outgoing,
            incoming,
        }
    }

    fn find_furthest_node(node_candidate: &[String], subgraph: &GraphicalGenome, start_node: &str) -> String {
        let mut max_distance = -1;
        let mut node = "".to_string();
        for n in node_candidate {
            if let (Some(pos_n), Some(pos_start)) = (subgraph.anchor[n]["pos"].as_i64(), subgraph.anchor[start_node]["pos"].as_i64()) {
                let d = (pos_n - pos_start).abs(); 
                if d > max_distance {
                    node = n.clone();
                    max_distance = d;
                }
            }
        }
        node
    }

    fn series_parallel_graph_nodelist( subgraph: &GraphicalGenome) -> Vec<String> {
        let mut keys: Vec<_> = subgraph.anchor.keys().cloned().collect();
        keys.sort();
        let start_node = keys[0].clone();
        let end_node = keys[keys.len() - 1].clone();
        println!("{}, {}", start_node, end_node);
        let mut nodelist = vec![start_node.clone()];

        let mut node_candidate: Vec<String> = subgraph.outgoing[&start_node]
            .iter()
            .flat_map(|edge| &subgraph.outgoing[edge])
            .filter(|node| subgraph.anchor.contains_key(*node))
            .cloned()
            .collect();
        node_candidate.sort();
        let mut node = Self::find_furthest_node(&node_candidate, subgraph, &start_node);
        nodelist.push(node.clone());

        while node != end_node && node != "" {
            let edgelist = &subgraph.outgoing[&node];
            let mut node_candidate: Vec<String> = Vec::new();
            for edge in edgelist {
                let nodelist = &subgraph.outgoing[edge];
                if nodelist.contains(&"SINK".to_string()) {
                    continue;
                }
                if nodelist[0] != "" && subgraph.anchor.contains_key(&nodelist[0]) && subgraph.outgoing.contains_key(&nodelist[0]) {
                    node_candidate.extend(nodelist.iter().cloned());
                }
            }
            node_candidate.sort();
            node = Self::find_furthest_node(&node_candidate, subgraph, &node);
            if nodelist.contains(&node) {
                nodelist.push(node.clone());
                break;
            }
            nodelist.push(node.clone());
        }

        nodelist
    }


    fn series_parallel_graph(nodelist: &[String], subgraph: &GraphicalGenome) -> (HashMap<String, Value>, HashMap<String, Value>, HashMap<String, Vec<String>>, HashMap<String, Vec<String>>) {
        let mut node_dict = HashMap::new();
        let mut edge_dict: HashMap<String, serde_json::Value>  = HashMap::new();
        let mut outgoing_dict:HashMap<String, Vec<String>> = HashMap::new();
        let mut incoming_dict:HashMap<String, Vec<String>> = HashMap::new();

        for (i, node) in nodelist.iter().enumerate().take(nodelist.len() - 1) {
            let start_node = node;
            let end_node = &nodelist[i + 1];
            node_dict.insert(start_node.clone(), subgraph.anchor[start_node].clone());
            let initial_set = find_all_reads(subgraph);
            let path = FindAllPathBetweenAnchors::new(subgraph, start_node, end_node, initial_set.clone());
            let mut index = 0;
            for (p, rs) in path.subpath.iter() {
                let edgename = format!("E{:05}.{:04}", start_node[1..].parse::<usize>().unwrap(), index);
                let seq = reconstruct_path_seq(subgraph, &p[1..p.len() - 1]);
                let edgelist = outgoing_dict.get(start_node).cloned().unwrap_or_else(Vec::new);
                let mut found = false;

                for edge in edgelist {
                    if edge_dict[&edge]["dst"] != end_node.as_str() {
                        continue;
                    }
                    if edge_dict[&edge]["seq"] == seq {
                        let reads: HashSet<_> = rs.iter().cloned().collect();
                        let existing_reads: HashSet<_> = edge_dict.get_mut(&edge).unwrap()["reads"].as_array().unwrap().iter().map(|v| v.as_str().unwrap().to_string()).collect();
                        let updated_reads: HashSet<_> = existing_reads.union(&reads).cloned().collect();
                        edge_dict.get_mut(&edge).unwrap()["reads"] = serde_json::to_value(updated_reads.iter().cloned().collect::<Vec<_>>()).unwrap();
                        edge_dict.get_mut(&edge).unwrap()["samples"] = serde_json::to_value(rs.iter().map(|item| item.split('|').last().unwrap().to_string()).collect::<HashSet<_>>()).unwrap();
                        found = true;
                        break;
                    }
                }
                if !found {
                    let edgename = format!("E{:05}.{:04}", start_node[1..].parse::<usize>().unwrap(), index);
                    let mut edge_info = serde_json::Map::new();
                    edge_info.insert("seq".to_string(), serde_json::to_value(&seq).unwrap());
                    edge_info.insert("src".to_string(), serde_json::to_value(start_node).unwrap());
                    edge_info.insert("dst".to_string(), serde_json::to_value(end_node).unwrap());
                    edge_info.insert("reads".to_string(), serde_json::to_value(rs).unwrap());
                    edge_info.insert("samples".to_string(), serde_json::to_value(rs.iter().map(|item| item.split('|').last().unwrap().to_string()).collect::<HashSet<_>>()).unwrap());

                    edge_dict.insert(edgename.clone(), serde_json::Value::Object(edge_info));
                    outgoing_dict.entry(start_node.clone()).or_default().push(edgename.clone());
                    outgoing_dict.insert(edgename.clone(), vec![end_node.clone()]);
                    incoming_dict.entry(end_node.clone()).or_default().push(edgename.clone());
                    incoming_dict.insert(edgename.clone(), vec![start_node.clone()]);
                    index += 1;
                }
            }
        }

        (node_dict, edge_dict, outgoing_dict, incoming_dict)
    }
}  

pub fn write_graph_from_graph(filename: &str, graph: &GraphicalGenome) -> std::io::Result<()> {
    let mut file = File::create(filename)?;

    writeln!(file, "H\tVN:Z:1.0")?;

    for (anchor, data) in &graph.anchor {
        let seq = data["seq"].as_str().unwrap_or_default();
        let mut data_clone = data.clone();
        data_clone.as_object_mut().unwrap().remove("seq");
        let json_string = serde_json::to_string(&data_clone).unwrap_or_else(|_| "{}".to_string());
        writeln!(file, "S\t{}\t{}\tPG:J:{}", anchor, seq, json_string)?;
    }
    let mut edge_output = Vec::new();
    let mut link_output = Vec::new();

    for (edge, edge_data) in &graph.edges {
        let seq = edge_data["seq"].as_str().unwrap_or_default();
        let src = graph.incoming[edge][0].clone();
        let dst = graph.outgoing[edge][0].clone();
        let mut edge_data_clone = edge_data.clone();
        edge_data_clone.as_object_mut().unwrap().remove("seq");
        let json_string = serde_json::to_string(&edge_data_clone).unwrap_or_else(|_| "{}".to_string());
        let formatted_string = if edge_data.get("reads").and_then(|r| r.as_array()).map_or(false, |arr| !arr.is_empty()) {
            format!("S\t{}\t{}\tPG:J:{}\tRC:i:{}", edge, seq, json_string, edge_data.get("reads").and_then(|r| r.as_array()).map_or(0, |arr| arr.len()))
        } else {
            format!("S\t{}\t{}", edge, seq)
        };
        edge_output.push(formatted_string);

        link_output.push(format!("L\t{}\t+\t{}\t+\t0M", src, edge));
        link_output.push(format!("L\t{}\t+\t{}\t+\t0M", edge, dst));
    }
    for s in edge_output {
        writeln!(file, "{}", s)?;
    }
    for l in link_output {
        writeln!(file, "{}", l)?;
    }
    Ok(())
}


pub fn start(
    output: &PathBuf,
    k: usize,
    fasta_path: &PathBuf,
    reference_name: String
) {
    // Use the basename of the reads fasta to be the reference sequence name in the graph.
    let stem = fasta_path
        .file_stem()
        .unwrap()
        .to_str()
        .unwrap()
        .to_string();

    // Read the reads records (name and sequence) into a vector.
    let reader = Reader::from_file(fasta_path).unwrap();
    let all_reads: Vec<Record> = reader.records().map(|r| r.unwrap()).collect();

    // Find the reference sequence in the read sequences and create a new vector excluding the reference.
    let mut reference = String::new();
    let mut reads: Vec<Record> = Vec::new();
    for record in all_reads {
        if record.id().ends_with(&reference_name) {
            skydive::elog!("Using '{}' as reference for anchor discovery.", record.id());

            reference = String::from_utf8_lossy(record.seq()).to_string();
        } else {
            reads.push(record);
        }
    }

    // Create the k-mer profile.
    let unique_kmer_list = get_reference_kmer_profile(&reference, k);
    let (hla_samples, hla_seq) =
        find_sequences_between_sanchor_eanchor(reads, reference, &stem);
    let (sample_dict, position_dict) =
        map_reference_unique_kmers_to_seq(unique_kmer_list, &hla_seq, k);

    // Compute anchors.
    let anchorlist = get_anchor_information(
        &sample_dict.lock().unwrap(),
        &hla_samples,
        &position_dict.lock().unwrap(),
        &stem,
    );
    let anchors = get_anchors(&anchorlist, &position_dict.lock().unwrap(), k, &stem);

    let final_anchor = get_final_anchor(&anchors, k);
    let (edge_info, outgoing) = create_edge_file(&hla_seq, &final_anchor, k);

    let dereferenced_final_anchor = final_anchor
        .iter()
        .map(|(k, v)| (k.clone(), (*v).clone()))
        .collect();

    // Filter edges with little read support.
    let filtered_edges = filter_undersupported_edges(&edge_info, &stem, 4);

    // Write final graph to disk.
    write_gfa(&dereferenced_final_anchor, &filtered_edges, output);

    let graph = GraphicalGenome::load_graph(output.to_str().unwrap()).unwrap();
    // println!("{:?}", graph.anchor)

    let mut sp_graph = GetSeriesParallelGraph::new(&graph);
    let outputfilename_str = output.with_extension("sp.gfa");
    println!("{:?}", outputfilename_str);
    write_graph_from_graph(outputfilename_str.to_str().unwrap(), &mut sp_graph);    println!("{:?}", outputfilename_str);
    // write_graph_from_graph("HLA-A.sp.gfa", &mut sp_graph);

}
