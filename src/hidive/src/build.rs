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

pub fn start(
    output: &PathBuf,
    loci_list: &Vec<String>,
    k: usize,
    fasta_path: &PathBuf,
    reference_path: &PathBuf,
) {
    // Use the basename of the reads fasta to be the reference sequence name in the graph.
    let stem = fasta_path
        .file_stem()
        .unwrap()
        .to_str()
        .unwrap()
        .to_string();

    let mut faidx = IndexedReader::from_file(&reference_path).unwrap();
    for (chr, start, stop) in skydive::utils::parse_loci(loci_list).iter() {
        // Fetch the requested sequence from the reference.
        faidx
            .fetch(chr, *start, *stop)
            .expect("Couldn't fetch interval");

        // Read the reference bases
        let mut seq = Vec::new();
        faidx.read(&mut seq).expect("Couldn't read the interval");

        let reference =
            String::from_utf8(seq.to_ascii_uppercase()).expect("Invalid UTF-8 sequence");

        println!("{}", reference);

        // Read the reads records (name and sequence) into a vector.
        let reader = Reader::from_file(fasta_path).unwrap();
        let reads: Vec<Record> = reader.records().map(|r| r.unwrap()).collect();

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

        // let anchor_list: Vec<String> = final_anchor.keys().cloned().collect();
        // let seqlist: Vec<String> = anchor_list.iter()
        //     .filter_map(|anchor| final_anchor.get(anchor).map(|info| info.get_seq().clone()))
        //     .collect();

        // println!("{:?}", seqlist);

        // let graph_filename =
        let graph = GraphicalGenome::load_graph(output.to_str().unwrap()).unwrap();
        // println!("{:?}", graph.anchor)
    }
}
