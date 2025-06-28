use std::collections::{HashMap, HashSet};
use std::{fs::File, io::{Write, BufWriter}, path::PathBuf};

use needletail::parser::SequenceRecord as NeedletailRecord;
use needletail::Sequence;
use needletail::parse_fastx_file;

use petgraph::graph::DiGraph;
use url::Url;

use bio::alignment::sparse::{find_kmer_matches, lcskpp, SparseAlignmentResult};

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
    Contained,    // seq1 is contained within seq2
    Container,    // seq1 contains seq2
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

    // Print overlapping sequence IDs and sequences for each PacBio read
    for node_idx in graph.node_indices() {
        let seq_record = &graph[node_idx];
        
        // Only process PacBio reads (and only forward orientation to avoid duplicates)
        if seq_record.source == SequenceSource::PacBio {
            // Collect unique neighbors to avoid duplicates
            let mut neighbors: Vec<_> = graph.neighbors(node_idx).collect();
            neighbors.sort(); // Sort to ensure consistent ordering
            neighbors.dedup(); // Remove any duplicates
            
            if !neighbors.is_empty() {
                // Create a POA graph of overlapping reads
                let mut sg = spoa::Graph::new();
                let mut la_ov = spoa::AlignmentEngine::new(spoa::AlignmentType::kOV, 5, -4, -8, -6, -8, -4);
                let mut la_sw = spoa::AlignmentEngine::new(spoa::AlignmentType::kSW, 5, -4, -8, -6, -8, -4);

                // Track sequence names in order they're added to the graph
                let mut sequence_names = Vec::new();

                // Add PacBio read first
                let pb_seq_cstr = std::ffi::CString::new(&seq_record.sequence[..]).unwrap();
                let pb_qual_cstr = std::ffi::CString::new(vec![b'I'; seq_record.sequence.len()]).unwrap();
                let a = la_ov.align(pb_seq_cstr.as_ref(), &sg);
                sg.add_alignment(&a, pb_seq_cstr.as_ref(), pb_qual_cstr.as_ref());
                sequence_names.push(seq_record.id.clone());

                // Process all neighbors in a single pass, using appropriate alignment method
                for neighbor in neighbors {
                    let neighbor_seq = &graph[neighbor];
                    
                    // Skip if this is the same sequence (shouldn't happen, but safety check)
                    if neighbor == node_idx {
                        continue;
                    }
                    
                    let seq_cstr = std::ffi::CString::new(&neighbor_seq.sequence[..]).unwrap();
                    let seq_qual = std::ffi::CString::new(vec![b'I'; neighbor_seq.sequence.len()]).unwrap();
                    
                    // Use overlap alignment for long reads, Smith-Waterman for short reads
                    let a = match neighbor_seq.source {
                        SequenceSource::PacBio | SequenceSource::Nanopore => {
                            la_ov.align(seq_cstr.as_ref(), &sg)
                        }
                        SequenceSource::Illumina => {
                            la_sw.align(seq_cstr.as_ref(), &sg)
                        }
                    };
                    
                    sg.add_alignment(&a, seq_cstr.as_ref(), seq_qual.as_ref());
                    sequence_names.push(neighbor_seq.id.clone());
                }

                // Get MSA
                let msa = sg.multiple_sequence_alignment(true);
                for (i, seq) in msa.iter().enumerate() {
                    let read_name = if i < sequence_names.len() {
                        &sequence_names[i]
                    } else {
                        "consensus"
                    };
                    skydive::elog!(
                        "\tMSA sequence: {} ({})", 
                        seq.to_str().unwrap(),
                        read_name,
                    );
                }
                skydive::elog!("");
            }
        }
    }
}

/// Write the overlap graph to GFA (Graphical Fragment Assembly) format
fn write_gfa_format(graph: &DiGraph<SequenceRecord, OverlapInfo>, output: &PathBuf) {
    let mut writer = BufWriter::new(File::create(output).expect("Failed to create output file"));

    // Write header
    writeln!(writer, "H\tVN:Z:1.0").expect("Failed to write GFA header");

    // Write sequences (nodes)
    for node_idx in graph.node_indices() {
        let seq_record = &graph[node_idx];
        writeln!(
            writer,
            "S\t{}\t{}\t*",
            seq_record.id,
            String::from_utf8_lossy(&seq_record.sequence)
        ).expect("Failed to write sequence line");
    }

    // Write overlaps (edges)
    for edge_idx in graph.edge_indices() {
        let (from_idx, to_idx) = graph.edge_endpoints(edge_idx).unwrap();
        let from_seq = &graph[from_idx];
        let to_seq = &graph[to_idx];
        let overlap = &graph[edge_idx];

        // Only write edges between forward sequences
        if from_seq.is_forward == to_seq.is_forward {
            let orientation = match overlap.overlap_type {
                OverlapType::SuffixPrefix => ('+', '+'),
                OverlapType::PrefixSuffix => ('+', '-'),
                OverlapType::Contained => ('+', '+'), // Contained reads have same orientation
                OverlapType::Container => ('+', '+'), // Container reads have same orientation
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
                    id: format!("{}-fw", String::from_utf8_lossy(record.id())),
                    is_forward: true,
                    sequence: fw_seq,
                    quality: fw_qual,
                    source: source.clone(),
                };
                
                let rc_record = SequenceRecord {
                    id: format!("{}-rc", String::from_utf8_lossy(record.id())),
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

/// Find overlap between two sequences using sparse alignment
fn find_overlap(seq1: &SequenceRecord, seq2: &SequenceRecord, min_overlap: usize) -> Option<OverlapInfo> {
    let matches = find_kmer_matches(&seq1.sequence, &seq2.sequence, min_overlap);
    let sparse_al = lcskpp(&matches, min_overlap);

    if sparse_al.score >= min_overlap as u32 {
        return Some(OverlapInfo {
            overlap_length: 1,
            overlap_type: OverlapType::SuffixPrefix, // Default assumption
        });
    }

    None
}