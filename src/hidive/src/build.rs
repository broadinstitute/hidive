// Import necessary standard library modules
use std::path::PathBuf;
use std::collections::HashSet;

// Import the Absolutize trait to convert relative paths to absolute paths
use path_absolutize::Absolutize;

// Import the Url type to work with URLs
use url::Url;

// Import the skydive module, which contains the necessary functions for building graphs
use skydive;

// From Hang Su

use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::ffi::FromVecWithNulError;
use std::fmt::DebugStruct;
use std::fs::File;
use std::io::Write; 
use std::hash::Hash;
use std::io::{self, BufRead, BufReader, Read};
use rayon::prelude::*;
use std::sync::Mutex;
use serde::{Serialize, Deserialize};
use serde_json::json;
use std::error::Error;

// use std::fs::File;

pub fn load_fasta(filename: &str) -> io::Result<(Vec<String>, Vec<String>)> {
    let file = File::open(filename)?;
    let reader: Box<dyn BufRead> = if filename.ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut headers = Vec::new();
    let mut sequences = Vec::new();
    let mut current_sequence = String::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_sequence.is_empty() {
                // Save the current sequence before starting a new one
                sequences.push(current_sequence.clone());
                current_sequence.clear();
            }
            // Save the header, trimming the leading '>'
            headers.push(line[1..].to_string());
        } else {
            // Append the line to the current sequence
            current_sequence.push_str(&line);
        }
    }

    // Don't forget to add the last sequence if the file doesn't end with a header
    if !current_sequence.is_empty() {
        sequences.push(current_sequence);
    }

    Ok((headers, sequences))
}

pub fn find_reference_and_source_sink_anchors(
    reference: &str,
    genestart: i32,
    geneend: i32,
    chromosome: i32,
) -> String {
    let (_headers, _sequences) = load_fasta(&reference).expect("fail to load reference file");
    println!("{}", "test 1");
    let chromo = (chromosome - 1) as usize;
    let contig = &_sequences[chromo];
    let end = geneend.min(contig.chars().count() as i32) as usize;
    let start = genestart as usize;
    let refseq: String = contig.chars().skip(start).take(end - start).collect();
    println!("{}", "test 2");
    return refseq.to_uppercase();
}
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

pub fn find_sequences_between_sanchor_eanchor(headers: Vec<String>, sequences: Vec<String>, ref_hla: String) -> (HashSet<String>, HashMap<String, String>){
    let mut hla_samples = HashSet::new();
    let mut hla_seq = HashMap::new();
    for (i, seq) in sequences.iter().enumerate() {
        let h = &headers[i];
        let seq_upper = seq.to_uppercase();
        let sample_identifier = h.split("|").last().unwrap_or_default().to_string();
        hla_samples.insert(sample_identifier.clone());
        hla_seq.insert(h.clone(), seq_upper);
    }
    hla_samples.insert("MHC-CHM13".to_string());
    hla_seq.insert("MHC-CHM13".to_string(), ref_hla);
    (hla_samples, hla_seq)
}

pub fn map_reference_unique_kmers_to_seq(unique_kmerlist: Vec<String>, hla_seq: &HashMap<String, String>, k: usize) -> (Mutex<HashMap<String, HashSet<String>>>, Mutex<HashMap<String, HashMap<String, Vec<usize>>>>){
    let sample_dict: Mutex<HashMap<String, HashSet<String>>> = Mutex::new(HashMap::new());
    let position_dict: Mutex<HashMap<String, HashMap<String, Vec<usize>>>> = Mutex::new(HashMap::new());    
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
pub struct EdgeInfo{
    seq: String,
    reads: Vec<String>,
    samples: HashSet<String>,
    src: String,
    dst:String,
}

impl AnchorInfo {
    // Getter method for `seq`
    pub fn get_seq(&self) -> &String {
        &self.seq
    }
}

pub fn get_anchor_information(
    sample_dict:&HashMap<String, HashSet<String>>, 
    hla_samples:&HashSet<String>, 
    position_dict:&HashMap<String, HashMap<String, Vec<usize>>>, 
    k: usize) 
    -> Vec<String> {
    let mut anchorlist = Vec::new();
    for (kmer, samplelist) in sample_dict.iter(){
        if samplelist.len() == hla_samples.len(){
            if let Some(positions) = position_dict.get(kmer){
                if positions.len()>1 && positions.get("MHC-CHM13").map_or(false, |v| v.len() == 1){
                    anchorlist.push(kmer.clone());
                }
            }

        }
    }
    anchorlist
}


pub fn get_anchors(anchorlist:&Vec<String>, position_dict: &HashMap<String, HashMap<String, Vec<usize>>>, k: usize) -> HashMap<String, AnchorInfo>{
    let mut anchor_info = HashMap::new();
    for kmer in anchorlist.iter(){
        if let Some(positions) = position_dict.get(kmer){
            if let Some(pos) = positions.get("MHC-CHM13").and_then(|v| v.first()){
                let anchor_name = format!("A{:06}", pos / k + 1);
                anchor_info.insert(anchor_name, AnchorInfo {
                    seq: kmer.clone(),
                    pos: *pos,
                });
            }
        }
    }

    anchor_info
}

pub fn get_final_anchor(anchor_info: & HashMap<String, AnchorInfo>, k: usize) -> HashMap<String, &AnchorInfo>{
    let mut final_anchor = HashMap::new();
    let mut anchornames: Vec<&String> = anchor_info.keys().collect();
    anchornames.sort();

    let mut anchor_unadjacent_list = Vec::new();
    let mut index = 0;

    while index < anchornames.len() - 1{
        let sanchor = anchornames[index];
        let mut danchor = sanchor;
        let mut next_index = index;
        for (i, next_anchor) in anchornames[index+1..].iter().enumerate(){
            if anchor_info.get(*next_anchor).unwrap().pos > anchor_info.get(sanchor).unwrap().pos + k + 1 {
                break
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

    for anchor_name in &anchor_unadjacent_list{
        if let Some(anchor) = anchor_info.get(anchor_name) {
            final_anchor.insert(anchor_name.clone(), anchor);
        }
    }

    final_anchor

}

pub fn mapping_info(final_anchor: &HashMap<String, &AnchorInfo>, contig: String, k: usize) -> (HashMap<String,Option<usize>>, HashMap<String, Vec<usize>>){
    let anchor_list: Vec<String> = final_anchor.keys().cloned().collect();
    let seqlist: Vec<String> = anchor_list.iter()
        .filter_map(|anchor| final_anchor.get(anchor).map(|info| info.get_seq().clone()))
        .collect();
    let mut position_dict:HashMap<String, Vec<usize>> = HashMap::new();
    for anchor_seq in seqlist.into_iter(){
        let anchor_rev = reverse_complement(&anchor_seq);
        position_dict.entry(anchor_seq.clone()).or_insert_with(|| Vec::new());
        position_dict.entry(anchor_rev.clone()).or_insert_with(|| Vec::new());
    }
    for i in 0..contig.len()-k+1 {
        let kmer: String = contig.chars().skip(i).take(k).collect();
        if let Some(positions) = position_dict.get_mut(&kmer) {
            positions.push(i);
        }
    }
    let mut a = HashMap::new();
    let mut svs = HashMap::new();
    for (anchor, info) in final_anchor.iter(){
        let anchor_seq = info.seq.clone();
        let positionlist = position_dict.get(&anchor_seq).unwrap_or(&vec![]).clone();
        if positionlist.len() == 1{
            a.insert(anchor.clone(), positionlist.first().copied());
        }else{
            svs.insert(anchor.clone(), positionlist.clone());
        }

    }
    (a, svs)

}
pub fn construct_edges(src_pos:usize, dst_pos: usize, k:usize, contig:String, contigname: String, sample:String, Anchorseq:&HashMap<String, String>) -> EdgeInfo {
    let src_seq: String;
    let mut pr = false;
    let mut src = "".to_string();

    if src_pos == 0 {
        src = "SOURCE".to_string();
        src_seq = "".to_string();
        pr = true;
    }else{
        src_seq = contig.chars().skip(src_pos).take(k).collect();
        src = match Anchorseq.get(&src_seq) {
            Some(value) => value.clone(),
            None => {
                let reversed_src_seq = reverse_complement(&src_seq);
                pr = true;
                Anchorseq.get(&reversed_src_seq).unwrap_or(&"".to_string()).clone()
            }
        };
    }
    
    let mut sr = false;
    let mut dst = "".to_string();
    let dst_seq;
    if dst_pos == contig.len(){
        dst = "SINK".to_string();
        dst_seq = "".to_string();
        sr = true;
    }else{
        dst_seq = contig.chars().skip(dst_pos).take(k).collect::<String>();
        dst = match Anchorseq.get(&dst_seq) {
            Some(value) => value.clone(),
            None => {
                let reversed_dst_seq = reverse_complement(&dst_seq);
                sr = true;
                Anchorseq.get(&reversed_dst_seq).unwrap_or(&"".to_string()).clone()
            }
        };

    }

    let mut edge_seq = if src_pos == 0 {
        contig.get(0..dst_pos).unwrap_or_default().to_string()
    } else {
        contig.get(src_pos + k..dst_pos).unwrap_or_default().to_string()
    };

    if pr && sr{
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

pub fn create_edge_file(hla_seq: &HashMap<String, String>, 
    final_anchor: &HashMap<String, &AnchorInfo>, 
    k: usize) -> (HashMap<String, EdgeInfo>, HashMap<String, Vec<String>>) {
        let mut edge_info:HashMap<String, EdgeInfo> = HashMap::new();
        let mut outgoing:HashMap<String, Vec<String>> = HashMap::new();
        // let mut edgenum_perread = Vec::new();
        let anchorseq: HashMap<_, _> = final_anchor.iter()
            .map(|(anchor, info)| (info.seq.clone(), anchor.clone()))
            .collect();
        let mut contig_index = 0;

        for (contig_name, contig) in hla_seq.iter(){
            contig_index += 1;
            let sample_name = contig_name.split('|').last().unwrap_or_default().to_string();
            let (a, svs) = mapping_info(&final_anchor, contig.to_string(), k);
            let mut splitposlist: Vec<_> = a.values().filter_map(|&x| x).collect();
            splitposlist.sort();
            let mut edgeindex = 0;
            let mut src_pos = 0;
            for &dst_pos in &splitposlist{
                let e = construct_edges(src_pos, dst_pos, k, contig.to_string(), contig_name.to_string(), sample_name.clone(), &anchorseq);
                let src = &e.src;
                let edgelist = outgoing.entry(src.clone()).or_insert_with(|| Vec::new());
                if let Some(pos) = edgelist.iter().position(|edge| {
                    edge_info[edge].dst == e.dst && edge_info[edge].seq == e.seq
                }) {
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
            let e = construct_edges(src_pos, dst_pos, k, contig.to_string(), contig_name.to_string(), sample_name.clone(), &anchorseq);
            let src = &e.src;
            let edgelist = outgoing.entry(src.clone()).or_insert_with(|| Vec::new());
            if let Some(pos) = edgelist.iter().position(|edge| {
                edge_info[edge].dst == e.dst && edge_info[edge].seq == e.seq
            }) {
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

pub fn write_gfa(final_anchor: &HashMap<String, AnchorInfo>, edge_info: &HashMap<String, EdgeInfo>, output_filename: &str) -> Result<(), Box<dyn Error>> {
    let mut file = File::create(output_filename)?;

    writeln!(file, "H\tVN:Z:1.0")?;

    let mut anchor_output = Vec::new();
    for (anchor, info) in final_anchor.iter(){
        let seq = &info.seq;
        let mut anchor_info_clone = HashMap::new();
        anchor_info_clone.insert("pos".to_string(), info.pos);
        // anchor_info_clone.seq = String::new();
        let json_string = serde_json::to_string(&anchor_info_clone).unwrap_or_else(|_| "{}".to_string());
        let formatted_string = format!("S\t{}\t{}\tPG:J:{}", anchor, seq, json_string);
        anchor_output.push(formatted_string);
        // writeln!(file,"S\t{}\t{}\tPG:J:{}\n", anchor, seq, json_string)?;
    }
    let mut edge_output = Vec::new();
    let mut link_output = Vec::new();

    for (edge, edge_data) in edge_info{
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

        let json_string = serde_json::to_string(&edge_data_clone).unwrap_or_else(|_| "{}".to_string());
        let formatted_string = if !edge_data.reads.is_empty() {
            format!("S\t{}\t{}\tPG:J:{}\tRC:i:{}", edge, seq, json_string, edge_data.reads.len())
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

pub fn filter_undersupported_edges(edge_info: &HashMap<String, EdgeInfo>, reference: &str, threshold: i32) -> HashMap<String, EdgeInfo> {
    let mut filtered_edges = edge_info.clone(); // Clone the input HashMap to create an owned version

    for (edge, data) in edge_info.iter() {
        if data.reads.contains(&reference.to_string()) || data.reads.len() >= threshold as usize {
            continue;
        }
        filtered_edges.remove(edge); // Remove the edge from the cloned HashMap
    }

    filtered_edges // Return the owned HashMap
}

pub fn start(output: &PathBuf, fasta_path: &PathBuf, reference_path: &PathBuf) {
    println!("Hello, Cargo!");
    // let filename = "/Users/suhang/Google Drive/My Drive/Analysis/HLAC/HLAE/HLA-E.aln.fa";
    let filename = fasta_path.to_str().unwrap();
    let (headers, sequences) = load_fasta(filename).expect("fail to load reads");
    println!("Hello, Cargo! 1");
    // let reference = "/Users/suhang/Analysis/chm13v2.0.ebv.fa.gz";
    let reference = reference_path.to_str().unwrap();
    // let genestart: i32 = 32408981;
    // let geneend: i32 = 32422402;
    let genestart = 29940532;
    let geneend = 29947870;

    // chr6:29,940,532-29,947,870

    let chromosome: i32 = 6;
    let reference_hla =
        find_reference_and_source_sink_anchors(reference, genestart, geneend, chromosome);
    println!("Hello, Cargo! 2");
    let k: usize = 11;
    let unique_kmer_list = get_reference_kmer_profile(&reference_hla, k);
    println!("Hello, Cargo! 3");
    let (hla_samples, hla_seq) = find_sequences_between_sanchor_eanchor(headers, sequences, reference_hla);
    println!("Hello, Cargo! 4");
    let (sample_dict, position_dict) = map_reference_unique_kmers_to_seq(unique_kmer_list,  &hla_seq, k);
    // println!("{:?}, {:?}", sample_dict, position_dict);
    println!("Hello, Cargo! 5");
    let anchorlist = get_anchor_information(&sample_dict.lock().unwrap(), &hla_samples, &position_dict.lock().unwrap(), k);
    println!("{}", anchorlist.len());
    println!("Hello, Cargo! 6");
    let anchors = get_anchors(&anchorlist, &position_dict.lock().unwrap(), k);
    println!("{}", anchors.len());
    let final_anchor = get_final_anchor(&anchors, k);
    println!("Hello, Cargo! 7");
    let (edge_info, outgoing) = create_edge_file(&hla_seq, &final_anchor, k);
    println!("Hello, Cargo! 8");


    let outputfilename = "test";
    let dereferenced_final_anchor: HashMap<_, AnchorInfo> = final_anchor.iter()
    .map(|(k, v)| (k.clone(), (*v).clone()))
    .collect();

    let reference = "MHC-CHM13";
    let filtered_edges = filter_undersupported_edges(&edge_info,reference, 4);

    write_gfa(&dereferenced_final_anchor, &filtered_edges, &format!("{}.trimmed.gfa", outputfilename));
      println!("Hello, Cargo! 9");

    // let anchor_list: Vec<String> = final_anchor.keys().cloned().collect();
    // let seqlist: Vec<String> = anchor_list.iter()
    //     .filter_map(|anchor| final_anchor.get(anchor).map(|info| info.get_seq().clone()))
    //     .collect();

    // println!("{:?}", seqlist);
}