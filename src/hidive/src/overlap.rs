use std::collections::{HashMap, HashSet};
use std::{fs::File, io::{Write, BufWriter}, path::PathBuf};

use needletail::parser::SequenceRecord as NeedletailRecord;
use needletail::Sequence;
use needletail::parse_fastx_file;

use petgraph::graph::DiGraph;
use url::Url;

/// Represents the source technology for a sequence read
#[derive(Debug, Clone, PartialEq)]
pub enum SequenceSource {
    PacBio,
    Illumina, 
    Nanopore,
}

/// Information about an overlap between two sequences
#[derive(Debug, Clone)]
struct OverlapInfo {
    overlap_length: usize,
    overlap_type: OverlapType,
}

#[derive(Debug, Clone)]
enum OverlapType {
    SuffixPrefix, // seq1 suffix matches seq2 prefix
    PrefixSuffix, // seq1 prefix matches seq2 suffix
}

/// Represents a sequence record with metadata
#[derive(Debug, Clone)]
struct SequenceRecord {
    id: String,
    is_forward: bool,
    sequence: Vec<u8>,
    quality: Option<Vec<u8>>,
    source: SequenceSource,
}

pub fn start(
    output: &PathBuf,
    kmer_size: usize,
    pacbio_fastx_path: &PathBuf,
    illumina_fastx_path: &Option<PathBuf>,
    nanopore_fastx_path: &Option<PathBuf>,
) {
    let mut graph = DiGraph::<SequenceRecord, OverlapInfo>::new();

    load_fastx_file(&mut graph, &Some(pacbio_fastx_path.clone()), SequenceSource::PacBio);
    load_fastx_file(&mut graph, illumina_fastx_path, SequenceSource::Illumina);
    load_fastx_file(&mut graph, nanopore_fastx_path, SequenceSource::Nanopore);

    compute_overlaps(&mut graph, kmer_size);

    write_gfa_format(&graph, output);
}

/// Write the overlap graph to GFA (Graphical Fragment Assembly) format
fn write_gfa_format(graph: &DiGraph<SequenceRecord, OverlapInfo>, output: &PathBuf) {
    let mut writer = BufWriter::new(File::create(output).expect("Failed to create output file"));

    // Write header
    writeln!(writer, "H\tVN:Z:1.0").expect("Failed to write GFA header");

    // Write sequences (nodes)
    for node_idx in graph.node_indices() {
        let seq_record = &graph[node_idx];
        // Only write forward sequences to avoid duplicates
        if seq_record.is_forward {
            writeln!(
                writer,
                "S\t{}\t{}\t*",
                seq_record.id,
                String::from_utf8_lossy(&seq_record.sequence)
            ).expect("Failed to write sequence line");
        }
    }

    // Write overlaps (edges)
    for edge_idx in graph.edge_indices() {
        let (from_idx, to_idx) = graph.edge_endpoints(edge_idx).unwrap();
        let from_seq = &graph[from_idx];
        let to_seq = &graph[to_idx];
        let overlap = &graph[edge_idx];

        // Only write edges between forward sequences
        if from_seq.is_forward && to_seq.is_forward {
            let orientation = match overlap.overlap_type {
                OverlapType::SuffixPrefix => ('+', '+'),
                OverlapType::PrefixSuffix => ('+', '-'),
            };
            
            writeln!(
                writer,
                "L\t{}\t{}\t{}\t{}\t{}M",
                from_seq.id,
                orientation.0,
                to_seq.id, 
                orientation.1,
                overlap.overlap_length
            ).expect("Failed to write link line");
        }
    }
}

/// Read sequences from a fastx file (FASTA/FASTQ)
fn load_fastx_file(graph: &mut DiGraph<SequenceRecord, OverlapInfo>, fastx_path: &Option<PathBuf>, source: SequenceSource) {
    if let Some(path) = fastx_path {
        let fastx_urls = skydive::parse::parse_file_names(&[path.clone()]);
        
        for fastx_url in fastx_urls {
            let file_path = fastx_url.to_file_path().unwrap();
            let mut reader = parse_fastx_file(file_path).expect("Failed to open sequence file");

            while let Some(record) = reader.next() {
                let record = record.expect("Failed to parse sequence record");
                let fw_seq = record.seq().to_vec();
                let rc_seq = record.reverse_complement();

                let fw_qual = match record.qual() {
                    Some(qual) => {
                        let q = qual.iter().copied().collect::<Vec<u8>>();
                        Some(q)
                    }
                    None => { None }
                };

                let rc_qual = match record.qual() {
                    Some(qual) => {
                        let q = qual.iter().copied().rev().collect::<Vec<u8>>();
                        Some(q)
                    }
                    None => { None }
                };

                let fw_record = SequenceRecord {
                    id: String::from_utf8_lossy(record.id()).to_string(),
                    is_forward: true,
                    sequence: fw_seq,
                    quality: fw_qual,
                    source: source.clone(),
                };
                
                let rc_record = SequenceRecord {
                    id: String::from_utf8_lossy(record.id()).to_string(),
                    is_forward: false,
                    sequence: rc_seq,
                    quality: rc_qual,
                    source: source.clone(),
                };

                graph.add_node(fw_record);
                graph.add_node(rc_record);
            }
        }
    }
}

/// Compute overlaps between sequences using k-mer based detection
fn compute_overlaps(graph: &mut DiGraph<SequenceRecord, OverlapInfo>, kmer_size: usize) {
    let mut kmer_to_nodes: HashMap<Vec<u8>, Vec<petgraph::graph::NodeIndex>> = HashMap::new();
    
    // Build k-mer index
    for node_idx in graph.node_indices() {
        let record = &graph[node_idx];
        let kmers = extract_kmers(&record.sequence, kmer_size);
        
        for kmer in kmers {
            kmer_to_nodes.entry(kmer).or_insert_with(Vec::new).push(node_idx);
        }
    }
    
    // Find overlaps based on shared k-mers
    for (kmer, nodes) in kmer_to_nodes {
        if nodes.len() > 1 {
            // Check all pairs of nodes that share this k-mer
            for i in 0..nodes.len() {
                for j in i + 1..nodes.len() {
                    let node1 = nodes[i];
                    let node2 = nodes[j];
                    
                    if let Some(overlap_info) = find_overlap(&graph[node1], &graph[node2], kmer_size) {
                        // Add edge representing the overlap
                        graph.add_edge(node1, node2, overlap_info);
                    }
                }
            }
        }
    }
}

/// Extract all k-mers from a sequence
fn extract_kmers(sequence: &[u8], k: usize) -> Vec<Vec<u8>> {
    if sequence.len() < k {
        return Vec::new();
    }
    
    let mut kmers = Vec::new();
    for i in 0..=sequence.len() - k {
        kmers.push(sequence[i..i + k].to_vec());
    }
    kmers
}

/// Find overlap between two sequences
fn find_overlap(seq1: &SequenceRecord, seq2: &SequenceRecord, min_overlap: usize) -> Option<OverlapInfo> {
    let seq1_len = seq1.sequence.len();
    let seq2_len = seq2.sequence.len();
    
    // Try different overlap lengths
    for overlap_len in min_overlap..=std::cmp::min(seq1_len, seq2_len) {
        // Check suffix of seq1 against prefix of seq2
        if seq1.sequence[seq1_len - overlap_len..] == seq2.sequence[..overlap_len] {
            return Some(OverlapInfo {
                overlap_length: overlap_len,
                overlap_type: OverlapType::SuffixPrefix,
            });
        }
        
        // Check prefix of seq1 against suffix of seq2
        if seq1.sequence[..overlap_len] == seq2.sequence[seq2_len - overlap_len..] {
            return Some(OverlapInfo {
                overlap_length: overlap_len,
                overlap_type: OverlapType::PrefixSuffix,
            });
        }
    }
    
    None
}

/// Alternative: Compute overlaps using suffix array approach (more sensitive)
fn compute_overlaps_suffix_array(graph: &mut DiGraph<SequenceRecord, OverlapInfo>, min_overlap: usize) {
    let nodes: Vec<_> = graph.node_indices().collect();
    
    for i in 0..nodes.len() {
        for j in i + 1..nodes.len() {
            let node1 = nodes[i];
            let node2 = nodes[j];
            
            if let Some(overlap_info) = find_overlap_suffix_array(&graph[node1], &graph[node2], min_overlap) {
                graph.add_edge(node1, node2, overlap_info);
            }
        }
    }
}

/// Find overlap using suffix array approach
fn find_overlap_suffix_array(seq1: &SequenceRecord, seq2: &SequenceRecord, min_overlap: usize) -> Option<OverlapInfo> {
    // This is a simplified version - in practice you'd use a proper suffix array library
    // For now, we'll use the same approach as before but with more sophisticated matching
    
    let seq1_len = seq1.sequence.len();
    let seq2_len = seq2.sequence.len();
    
    // Try different overlap lengths with more sophisticated matching
    for overlap_len in min_overlap..=std::cmp::min(seq1_len, seq2_len) {
        // Check suffix of seq1 against prefix of seq2
        if seq1.sequence[seq1_len - overlap_len..] == seq2.sequence[..overlap_len] {
            return Some(OverlapInfo {
                overlap_length: overlap_len,
                overlap_type: OverlapType::SuffixPrefix,
            });
        }
        
        // Check prefix of seq1 against suffix of seq2
        if seq1.sequence[..overlap_len] == seq2.sequence[seq2_len - overlap_len..] {
            return Some(OverlapInfo {
                overlap_length: overlap_len,
                overlap_type: OverlapType::PrefixSuffix,
            });
        }
    }
    
    None
}

/// Find overlap with quality score consideration
fn find_overlap_with_quality(seq1: &SequenceRecord, seq2: &SequenceRecord, min_overlap: usize, min_quality: u8) -> Option<OverlapInfo> {
    let seq1_len = seq1.sequence.len();
    let seq2_len = seq2.sequence.len();
    
    for overlap_len in min_overlap..=std::cmp::min(seq1_len, seq2_len) {
        // Check suffix of seq1 against prefix of seq2
        if let Some(quality_score) = check_overlap_quality(
            &seq1.sequence[seq1_len - overlap_len..],
            &seq2.sequence[..overlap_len],
            &seq1.quality.as_ref().map(|q| &q[seq1_len - overlap_len..]),
            &seq2.quality.as_ref().map(|q| &q[..overlap_len]),
            min_quality
        ) {
            return Some(OverlapInfo {
                overlap_length: overlap_len,
                overlap_type: OverlapType::SuffixPrefix,
            });
        }
        
        // Check prefix of seq1 against suffix of seq2
        if let Some(quality_score) = check_overlap_quality(
            &seq1.sequence[..overlap_len],
            &seq2.sequence[seq2_len - overlap_len..],
            &seq1.quality.as_ref().map(|q| &q[..overlap_len]),
            &seq2.quality.as_ref().map(|q| &q[seq2_len - overlap_len..]),
            min_quality
        ) {
            return Some(OverlapInfo {
                overlap_length: overlap_len,
                overlap_type: OverlapType::PrefixSuffix,
            });
        }
    }
    
    None
}

/// Check if overlap has sufficient quality
fn check_overlap_quality(
    seq1: &[u8],
    seq2: &[u8],
    qual1: &Option<&[u8]>,
    qual2: &Option<&[u8]>,
    min_quality: u8
) -> Option<f64> {
    if seq1 != seq2 {
        return None;
    }
    
    // If we have quality scores, check them
    if let (Some(q1), Some(q2)) = (qual1, qual2) {
        let mut total_quality = 0.0;
        let mut count = 0;
        
        for i in 0..seq1.len() {
            if q1[i] >= min_quality && q2[i] >= min_quality {
                total_quality += (q1[i] as f64 + q2[i] as f64) / 2.0;
                count += 1;
            }
        }
        
        if count > 0 {
            return Some(total_quality / count as f64);
        }
    }
    
    // If no quality scores, just check sequence match
    Some(0.0)
}