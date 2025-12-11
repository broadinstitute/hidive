//! CRF training module for pangenome graph-based haplotype inference.
//!
//! This module trains a Conditional Random Field (CRF) model on a pangenome graph
//! to learn how to score nodes and edges for haplotype inference.

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use std::num::NonZeroU32;

use crate::pangenome::{get_tier_weight, AssemblyTier};
use skydive::elog;
use rucrf::{FeatureProvider, FeatureSet, Lattice, Edge, Trainer, Model, RawModel};
use bincode::{Encode, Decode};

/// Features extracted for CRF training
#[derive(Debug, Clone)]
pub struct CRFFeatures {
    /// Read coverage (kmer counts from PacBio reads)
    pub read_coverage: f64,
    /// Tier confidence weight
    pub tier_confidence: f64,
    /// Graph structure metrics (connectivity, path consistency)
    pub graph_structure: GraphStructureFeatures,
    /// Sequence similarity to reference
    pub sequence_similarity: f64,
}

/// Graph structure features for CRF
#[derive(Debug, Clone)]
pub struct GraphStructureFeatures {
    /// Number of incoming edges
    pub in_degree: usize,
    /// Number of outgoing edges
    pub out_degree: usize,
    /// Path consistency score
    pub path_consistency: f64,
    /// Connectivity score
    pub connectivity: f64,
}

/// Configuration for CRF training
pub struct CRFTrainingConfig {
    /// Path to the pangenome graph (GFA format)
    pub graph_path: PathBuf,
    /// Paths to PacBio read files (FASTA/FASTQ)
    pub read_paths: Vec<PathBuf>,
    /// Paths to ground truth haplotype files (FASTA)
    pub truth_haplotype_paths: Vec<PathBuf>,
    /// Output path for trained CRF model
    pub output_path: PathBuf,
    /// K-mer size for feature extraction
    pub kmer_size: usize,
    /// Number of training iterations
    pub iterations: usize,
}

/// Train a CRF model on the pangenome graph
pub fn train_crf(config: &CRFTrainingConfig) -> Result<(), Box<dyn std::error::Error>> {
    elog!("Starting CRF training...");
    elog!("Graph: {}", config.graph_path.display());
    elog!("Reads: {} files", config.read_paths.len());
    elog!("Ground truth haplotypes: {} files", config.truth_haplotype_paths.len());

    // Step 1: Load the pangenome graph
    elog!("Loading pangenome graph...");
    let graph = load_pangenome_graph(&config.graph_path)?;

    // Step 2: Align reads to the graph
    elog!("Aligning reads to pangenome graph...");
    let read_alignments = align_reads_to_graph(&config.read_paths, &graph, config.kmer_size)?;

    // Step 3: Extract features from nodes and edges
    elog!("Extracting features from graph nodes and edges...");
    let node_features = extract_node_features(&graph, &read_alignments, config.kmer_size)?;
    let edge_features = extract_edge_features(&graph, &read_alignments)?;

    // Step 4: Load ground truth haplotypes and create training labels
    elog!("Loading ground truth haplotypes...");
    let truth_haplotypes = load_truth_haplotypes(&config.truth_haplotype_paths)?;

    // Step 5: Determine truth paths by aligning truth haplotypes to graph
    // For now, we'll use a simple approach: find nodes that match truth sequences
    elog!("Determining truth paths from haplotypes...");
    let truth_paths_with_paths = determine_truth_paths(&graph, &truth_haplotypes)?;

    // Step 6: Train the CRF model
    elog!("Training CRF model...");
    let model = train_crf_model(
        &graph,
        &node_features,
        &edge_features,
        &truth_paths_with_paths,
        config.iterations,
    )?;

    // Step 7: Save the model
    elog!("Saving CRF model to {}...", config.output_path.display());
    save_crf_model(&model, &config.output_path)?;

    elog!("CRF training complete!");
    Ok(())
}

use crate::pangenome::PangenomeGraph;

/// Load a pangenome graph from GFA format
fn load_pangenome_graph(path: &PathBuf) -> Result<PangenomeGraph, Box<dyn std::error::Error>> {
    // Use the public function from pangenome module
    crate::pangenome::load_pangenome_graph(path)
}

/// Read alignment information
pub struct ReadAlignment {
    /// Read sequence
    pub sequence: Vec<u8>,
    /// Aligned node IDs
    pub aligned_nodes: Vec<u64>,
    /// Alignment scores
    pub scores: Vec<f64>,
}

/// Align reads to the pangenome graph
fn align_reads_to_graph(
    read_paths: &[PathBuf],
    graph: &PangenomeGraph,
    kmer_size: usize,
) -> Result<Vec<ReadAlignment>, Box<dyn std::error::Error>> {
    elog!("Aligning reads to graph using k-mer matching (k={})", kmer_size);
    
    // Build k-mer index: kmer -> list of (node_id, position_in_node)
    let mut kmer_index: HashMap<Vec<u8>, Vec<(u64, usize)>> = HashMap::new();
    
    for (&node_id, node_seq) in &graph.nodes {
        if node_seq.len() < kmer_size {
            continue;
        }
        // Extract all k-mers from this node
        for i in 0..=node_seq.len().saturating_sub(kmer_size) {
            let kmer = node_seq[i..i + kmer_size].to_vec();
            kmer_index.entry(kmer).or_insert_with(Vec::new).push((node_id, i));
        }
    }
    
    elog!("Built k-mer index: {} unique k-mers across {} nodes", kmer_index.len(), graph.nodes.len());
    
    let mut alignments = Vec::new();
    
    for read_path in read_paths {
        let reader = bio::io::fasta::Reader::from_file(read_path)?;
        for record in reader.records() {
            let record = record?;
            let seq = record.seq().to_vec();
            
            // Find nodes matching this read using k-mer matching
            let (aligned_nodes, scores) = find_matching_nodes(&seq, &kmer_index, kmer_size);
            
            alignments.push(ReadAlignment {
                sequence: seq,
                aligned_nodes,
                scores,
            });
        }
    }
    
    elog!("Aligned {} reads to graph", alignments.len());
    Ok(alignments)
}

/// Find nodes matching a sequence using k-mer matching
fn find_matching_nodes(
    seq: &[u8],
    kmer_index: &HashMap<Vec<u8>, Vec<(u64, usize)>>,
    kmer_size: usize,
) -> (Vec<u64>, Vec<f64>) {
    if seq.len() < kmer_size {
        return (Vec::new(), Vec::new());
    }
    
    // Count how many times each node appears in k-mer matches
    let mut node_match_counts: HashMap<u64, usize> = HashMap::new();
    let mut node_positions: HashMap<u64, Vec<usize>> = HashMap::new();
    
    // Extract k-mers from the read and find matching nodes
    for i in 0..=seq.len().saturating_sub(kmer_size) {
        let kmer = seq[i..i + kmer_size].to_vec();
        
        if let Some(node_matches) = kmer_index.get(&kmer) {
            for &(node_id, node_pos) in node_matches {
                *node_match_counts.entry(node_id).or_insert(0) += 1;
                node_positions.entry(node_id).or_insert_with(Vec::new).push(i);
            }
        }
    }
    
    // Score nodes based on k-mer match count and coverage
    // Nodes with more k-mer matches are more likely to be part of the true alignment
    let mut node_scores: Vec<(u64, f64)> = node_match_counts
        .iter()
        .map(|(&node_id, &count)| {
            // Score = normalized match count (higher is better)
            // Also consider how well the matches span the read
            let positions = node_positions.get(&node_id).unwrap();
            let coverage_ratio = count as f64 / seq.len().saturating_sub(kmer_size).max(1) as f64;
            let score = count as f64 * (1.0 + coverage_ratio);
            (node_id, score)
        })
        .collect();
    
    // Sort by score (descending) and take top matches
    node_scores.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    
    // Take nodes that have at least 2 k-mer matches (to filter noise)
    let threshold = 2.0;
    let (nodes, scores): (Vec<u64>, Vec<f64>) = node_scores
        .into_iter()
        .filter(|(_, score)| *score >= threshold)
        .unzip();
    
    (nodes, scores)
}

/// Extract features for each node in the graph
fn extract_node_features(
    graph: &PangenomeGraph,
    alignments: &[ReadAlignment],
    kmer_size: usize,
) -> Result<HashMap<u64, CRFFeatures>, Box<dyn std::error::Error>> {
    // Calculate read coverage for each node
    // Count how many reads align to each node (not k-mer matches)
    let mut node_read_count: HashMap<u64, usize> = HashMap::new();
    let mut node_total_score: HashMap<u64, f64> = HashMap::new();
    
    for alignment in alignments {
        // Count each read only once per node (even if multiple k-mers match)
        let mut node_scores_this_read: HashMap<u64, f64> = HashMap::new();
        
        for (idx, &node_id) in alignment.aligned_nodes.iter().enumerate() {
            // Sum up alignment scores for this node from this read
            let score = alignment.scores.get(idx).copied().unwrap_or(1.0);
            *node_scores_this_read.entry(node_id).or_insert(0.0) += score;
        }
        
        // Count each node once per read, but weight by alignment score
        for (node_id, score) in node_scores_this_read {
            *node_read_count.entry(node_id).or_insert(0) += 1;
            *node_total_score.entry(node_id).or_insert(0.0) += score;
        }
    }
    
    // Extract features for each node
    let mut node_features = HashMap::new();
    
    // Include all nodes in the graph, even if they have no coverage
    for &node_id in graph.nodes.keys() {
        let read_count = node_read_count.get(&node_id).copied().unwrap_or(0);
        let total_score = node_total_score.get(&node_id).copied().unwrap_or(0.0);
        
        // Coverage = number of reads + weighted by alignment scores
        let read_coverage = read_count as f64 + total_score * 0.1;
        
        // Get tier confidence from path name (if available)
        let tier_confidence = get_tier_confidence_for_node(graph, node_id);
        
        // Calculate graph structure features
        let graph_structure = calculate_graph_structure_features(graph, node_id);
        
        // Calculate sequence similarity (placeholder)
        let sequence_similarity = 1.0; // TODO: Calculate actual similarity
        
        node_features.insert(
            node_id,
            CRFFeatures {
                read_coverage,
                tier_confidence,
                graph_structure,
                sequence_similarity,
            },
        );
    }
    
    Ok(node_features)
}

/// Extract features for each edge in the graph
fn extract_edge_features(
    graph: &PangenomeGraph,
    alignments: &[ReadAlignment],
) -> Result<HashMap<(u64, u64), CRFFeatures>, Box<dyn std::error::Error>> {
    // Count edge traversals from read alignments
    let mut edge_coverage: HashMap<(u64, u64), usize> = HashMap::new();
    for alignment in alignments {
        for window in alignment.aligned_nodes.windows(2) {
            let edge = (window[0], window[1]);
            *edge_coverage.entry(edge).or_insert(0) += 1;
        }
    }
    
    let mut edge_features = HashMap::new();
    for &(from, to) in &graph.edges {
        let coverage = edge_coverage.get(&(from, to)).copied().unwrap_or(0);
        let tier_confidence = (get_tier_confidence_for_node(graph, from)
            + get_tier_confidence_for_node(graph, to))
            / 2.0;
        
        edge_features.insert(
            (from, to),
            CRFFeatures {
                read_coverage: coverage as f64,
                tier_confidence,
                graph_structure: GraphStructureFeatures {
                    in_degree: 0,
                    out_degree: 0,
                    path_consistency: 1.0,
                    connectivity: 1.0,
                },
                sequence_similarity: 1.0,
            },
        );
    }
    
    Ok(edge_features)
}

/// Get tier confidence for a node
fn get_tier_confidence_for_node(graph: &PangenomeGraph, node_id: u64) -> f64 {
    // Try to find tier information from paths containing this node
    for (path_name, nodes) in &graph.paths {
        if nodes.contains(&node_id) {
            return crate::pangenome::get_tier_weight(path_name);
        }
    }
    // Default to lowest tier confidence
    crate::pangenome::TIER_3_WEIGHT
}

/// Calculate graph structure features for a node
fn calculate_graph_structure_features(
    graph: &PangenomeGraph,
    node_id: u64,
) -> GraphStructureFeatures {
    let in_degree = graph.edges.iter().filter(|(_, to)| *to == node_id).count();
    let out_degree = graph.edges.iter().filter(|(from, _)| *from == node_id).count();
    
    // Calculate path consistency (simplified)
    let path_consistency = if in_degree > 0 && out_degree > 0 {
        1.0
    } else {
        0.5
    };
    
    // Calculate connectivity (simplified)
    let connectivity = (in_degree + out_degree) as f64 / 10.0; // Normalize
    
    GraphStructureFeatures {
        in_degree,
        out_degree,
        path_consistency,
        connectivity,
    }
}

/// Ground truth haplotype paths
pub struct TruthHaplotype {
    /// Haplotype name
    pub name: String,
    /// Sequence
    pub sequence: Vec<u8>,
    /// Path through graph (node IDs)
    pub path: Vec<u64>,
}

/// Load ground truth haplotypes
fn load_truth_haplotypes(
    paths: &[PathBuf],
) -> Result<Vec<TruthHaplotype>, Box<dyn std::error::Error>> {
    let mut haplotypes = Vec::new();
    
    for path in paths {
        let reader = bio::io::fasta::Reader::from_file(path)?;
        for record in reader.records() {
            let record = record?;
            haplotypes.push(TruthHaplotype {
                name: record.id().to_string(),
                sequence: record.seq().to_vec(),
                path: Vec::new(), // Will be determined by aligning to graph
            });
        }
    }
    
    Ok(haplotypes)
}

/// Determine truth paths by aligning truth haplotypes to graph nodes
fn determine_truth_paths(
    graph: &PangenomeGraph,
    truth_haplotypes: &[TruthHaplotype],
) -> Result<Vec<TruthHaplotype>, Box<dyn std::error::Error>> {
    let mut path_sequences: Vec<(String, Vec<u64>, Vec<u8>)> = Vec::new();
    for (path_name, node_ids) in &graph.paths {
        if node_ids.is_empty() {
            continue;
        }
        let mut seq = Vec::new();
        for &node_id in node_ids {
            if let Some(node_seq) = graph.nodes.get(&node_id) {
                seq.extend_from_slice(node_seq);
            }
        }
        path_sequences.push((path_name.clone(), node_ids.clone(), seq));
    }
    
    let mut truth_paths_with_paths = Vec::new();
    
    for truth_hap in truth_haplotypes {
        let mut best_path = Vec::new();
        
        // Step 0: Attempt to match entire path sequence via k-mer overlap
        let revcomp_truth = revcomp(&truth_hap.sequence);
        let mut best_identity = 0.0;
        let mut best_seq_nodes: Option<Vec<u64>> = None;
        let mut matched_name: Option<String> = None;
        for (path_name, node_ids, seq) in &path_sequences {
            let identity = kmer_identity(&truth_hap.sequence, seq, 13)
                .max(kmer_identity(&revcomp_truth, seq, 13));
            if identity > best_identity {
                best_identity = identity;
                best_seq_nodes = Some(node_ids.clone());
                matched_name = Some(path_name.clone());
            }
        }
        if let Some(nodes) = best_seq_nodes.clone() {
            if best_identity >= 0.4 {
                elog!(
                    "Truth haplotype '{}' matched graph path '{}' via k-mer identity {:.3} ({} nodes)",
                    truth_hap.name,
                    matched_name.unwrap_or_else(|| "<unknown>".into()),
                    best_identity,
                    nodes.len()
                );
                best_path = nodes;
            }
        }
        
        // Fallback: find single best-matching node if we couldn't map an entire path
        if best_path.is_empty() {
            let mut best_node: Option<u64> = None;
            let mut best_score = 0.0;
            for (&node_id, node_seq) in &graph.nodes {
                let k = 11;
                if truth_hap.sequence.len() < k || node_seq.len() < k {
                    continue;
                }
                let score = kmer_identity(&truth_hap.sequence, node_seq, k);
                if score > best_score {
                    best_score = score;
                    best_node = Some(node_id);
                }
            }
            if let Some(node_id) = best_node {
                elog!(
                    "Truth haplotype '{}' fallback-matched node {} (k-mer identity {:.2})",
                    truth_hap.name,
                    node_id,
                    best_score
                );
                best_path.push(node_id);
            } else {
                elog!(
                    "Warning: Could not find matching path or node for truth haplotype '{}'",
                    truth_hap.name
                );
            }
        }
        
        truth_paths_with_paths.push(TruthHaplotype {
            name: truth_hap.name.clone(),
            sequence: truth_hap.sequence.clone(),
            path: best_path,
        });
    }
    
    Ok(truth_paths_with_paths)
}

/// Strip tier annotation (`|tier=TierX`) to get the base path name
fn base_name_without_tier(name: &str) -> &str {
    name.split('|').next().unwrap_or(name)
}

fn encode_kmer(seq: &[u8]) -> u64 {
    let mut val = 0u64;
    for &b in seq {
        val <<= 2;
        val |= match b {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 0,
        };
    }
    val
}

fn kmer_identity(a: &[u8], b: &[u8], k: usize) -> f64 {
    if a.len() < k || b.len() < k {
        return 0.0;
    }
    let mut set: HashSet<u64> = HashSet::with_capacity(a.len().saturating_sub(k) + 1);
    for i in 0..=a.len() - k {
        set.insert(encode_kmer(&a[i..i + k]));
    }
    let mut matches = 0usize;
    for i in 0..=b.len() - k {
        if set.contains(&encode_kmer(&b[i..i + k])) {
            matches += 1;
        }
    }
    let denom = set.len().max(1);
    (matches.min(denom)) as f64 / denom as f64
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => *b,
        })
        .collect()
}

/// Training example for CRF
pub struct TrainingExample {
    /// Node/edge ID
    pub id: String,
    /// Features
    pub features: CRFFeatures,
    /// Label (1 if in truth path, 0 otherwise)
    pub label: f64,
}

// Note: create_training_examples is no longer needed as we build lattices directly
// Keeping the struct for potential future use
#[allow(dead_code)]
fn create_training_examples(
    _graph: &PangenomeGraph,
    _node_features: &HashMap<u64, CRFFeatures>,
    _edge_features: &HashMap<(u64, u64), CRFFeatures>,
    _truth_haplotypes: &[TruthHaplotype],
) -> Result<Vec<TrainingExample>, Box<dyn std::error::Error>> {
    // This function is kept for compatibility but not used
    // Training now uses lattices directly
    Ok(Vec::new())
}

/// Trained CRF model
/// Note: RawModel already implements Encode/Decode from bincode 2.0
#[derive(Encode, Decode)]
pub struct CRFModel {
    /// CRF model for path inference (includes feature provider)
    pub model: RawModel,
    /// Feature mapping: (node_id or edge) -> feature_set_id
    /// Note: feature_provider is part of RawModel, so we don't need to store it separately
    pub node_feature_map: HashMap<u64, u32>,
    pub edge_feature_map: HashMap<(u64, u64), u32>,
}

/// Convert features to feature IDs for rucrf
/// rucrf uses integer feature IDs, so we quantize continuous features
fn allocate_feature_id(counter: &mut u32) -> Result<NonZeroU32, Box<dyn std::error::Error>> {
    if *counter == 0 {
        *counter = 1;
    }
    if *counter == u32::MAX {
        return Err("Feature ID overflow".into());
    }
    let id = NonZeroU32::new(*counter).ok_or("Invalid feature ID")?;
    *counter += 1;
    Ok(id)
}

/// Train a CRF model using rucrf
fn train_crf_model(
    graph: &PangenomeGraph,
    node_features: &HashMap<u64, CRFFeatures>,
    edge_features: &HashMap<(u64, u64), CRFFeatures>,
    truth_haplotypes: &[TruthHaplotype],
    iterations: usize,
) -> Result<CRFModel, Box<dyn std::error::Error>> {
    elog!("Training CRF model with {} iterations...", iterations);
    
    // Create feature provider
    let mut provider = FeatureProvider::new();
    let mut node_feature_map: HashMap<u64, u32> = HashMap::new();
    let mut edge_feature_map: HashMap<(u64, u64), u32> = HashMap::new();
    
    let mut feature_id_counter: u32 = 1;
    
    // Register features for all nodes
    for (&node_id, _features) in node_features {
        let unique_feature = allocate_feature_id(&mut feature_id_counter)?;
        let empty: [Option<NonZeroU32>; 0] = [];
        let feature_set = FeatureSet::new(
            &[unique_feature],
            &empty,
            &empty,
        );
        let feature_set_id = provider.add_feature_set(feature_set)?;
        node_feature_map.insert(node_id, feature_set_id.get());
    }
    
    // Register features for all edges
    for (&edge, _features) in edge_features {
        let unique_feature = allocate_feature_id(&mut feature_id_counter)?;
        let empty: [Option<NonZeroU32>; 0] = [];
        let feature_set = FeatureSet::new(
            &[unique_feature],
            &empty,
            &empty,
        );
        let feature_set_id = provider.add_feature_set(feature_set)?;
        edge_feature_map.insert(edge, feature_set_id.get());
    }
    
    // Build adjacency for negative examples
    let mut adjacency: HashMap<u64, Vec<u64>> = HashMap::new();
    for &(from, to) in &graph.edges {
        adjacency.entry(from).or_default().push(to);
    }
    
    // Create training lattices from truth haplotypes
    let mut lattices = Vec::new();
    
    for truth_hap in truth_haplotypes {
        if truth_hap.path.is_empty() {
            continue;
        }
        
        let node_seq = &truth_hap.path;
        let lattice_len = node_seq.len().checked_mul(2)
            .and_then(|v| v.checked_sub(1))
            .ok_or("Invalid truth path length")?;
        let mut lattice = Lattice::new(lattice_len)?;
        
        for (idx, &node_id) in node_seq.iter().enumerate() {
            let node_feature_id = node_feature_map
                .get(&node_id)
                .ok_or_else(|| format!("Missing node feature for node {}", node_id))?;
            let node_feature_nz = NonZeroU32::new(*node_feature_id)
                .ok_or_else(|| format!("Invalid node feature ID {} for node {}", node_feature_id, node_id))?;
            lattice.add_edge(2 * idx, Edge::new(2 * idx + 1, node_feature_nz))?;
            
            if idx + 1 < node_seq.len() {
                let current = node_id;
                let next_true = node_seq[idx + 1];
                let neighbors = adjacency
                    .get(&current)
                    .ok_or_else(|| format!("No outgoing edges registered for node {}", current))?;
                // create ordered list with true neighbor first
                let mut neighbor_list = neighbors.clone();
                if let Some(pos) = neighbor_list.iter().position(|&n| n == next_true) {
                    neighbor_list.swap(0, pos);
                } else {
                    // ensure true edge exists even if not in adjacency (shouldn't happen)
                    neighbor_list.insert(0, next_true);
                }
                for &neighbor in &neighbor_list {
                    let edge = (current, neighbor);
                    let edge_feature_id = edge_feature_map
                        .get(&edge)
                        .ok_or_else(|| format!("Missing edge feature for transition {} -> {}", edge.0, edge.1))?;
                    let edge_feature_nz = NonZeroU32::new(*edge_feature_id)
                        .ok_or_else(|| format!("Invalid edge feature ID {} for edge {:?}", edge_feature_id, edge))?;
                    lattice.add_edge(2 * idx + 1, Edge::new(2 * (idx + 1), edge_feature_nz))?;
                }
            }
        }
        
        lattices.push(lattice);
    }
    
    if lattices.is_empty() {
        elog!("Warning: No training lattices created, using default model");
        // Create a minimal model
        let trainer = Trainer::new();
        let model = trainer.max_iter(1)?.train(&[], provider);
        
        return Ok(CRFModel {
            model,
            node_feature_map,
            edge_feature_map,
        });
    }
    
    elog!("Created {} training lattices", lattices.len());
    
    // Train the CRF model
    // Note: provider is consumed by train(), so we can't reuse it
    // We'll need to reconstruct it during inference
    let trainer = Trainer::new()
        .max_iter(iterations as u64)?;
    let model = trainer.train(&lattices, provider);
    
    elog!("CRF training complete!");
    if let Ok(debug_model) = model.merge() {
        let total_weight: f64 = debug_model
            .feature_sets
            .iter()
            .map(|fs| fs.weight.abs())
            .sum();
        elog!("Merged model L1 weight sum: {:.6}", total_weight);
    }
    
    // Note: The feature provider is now part of RawModel (stored inside it),
    // so we don't need to reconstruct it separately. The node_feature_map and
    // edge_feature_map are already populated during feature registration.
    
    Ok(CRFModel {
        model,
        node_feature_map,
        edge_feature_map,
    })
}

/// Save a CRF model to disk using bincode
fn save_crf_model(model: &CRFModel, path: &PathBuf) -> Result<(), Box<dyn std::error::Error>> {
    use std::fs::File;
    use std::io::Write;
    
    // Serialize the entire model using bincode 2.0
    // RawModel implements Encode from bincode, so we can serialize it directly
    let encoded = bincode::encode_to_vec(model, bincode::config::standard())?;
    
    let mut file = File::create(path)?;
    file.write_all(&encoded)?;
    
    elog!("CRF model saved successfully using bincode ({} bytes)", encoded.len());
    
    Ok(())
}

/// Load a CRF model from disk using bincode
/// Supports both new bincode format and legacy JSON format for backward compatibility
pub fn load_crf_model(path: &PathBuf) -> Result<CRFModel, Box<dyn std::error::Error>> {
    use std::io::Read;
    
    let mut file = std::fs::File::open(path)?;
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)?;
    
    // Try to deserialize as bincode 2.0 first (new format)
    match bincode::decode_from_slice(&buffer, bincode::config::standard()) {
        Ok((model, _)) => {
            elog!("Loaded CRF model from bincode format");
            return Ok(model);
        }
        Err(_) => {
            // If bincode fails, try legacy JSON format for backward compatibility
            elog!("Bincode deserialization failed, trying legacy JSON format...");
        }
    }
    
    // Legacy JSON format support (for models saved before bincode implementation)
    let json: serde_json::Value = serde_json::from_slice(&buffer)?;
    
    let node_feature_map: HashMap<u64, u32> = json["node_feature_map"]
        .as_object()
        .ok_or("node_feature_map not found or not an object")?
        .iter()
        .map(|(k, v)| -> Result<(u64, u32), Box<dyn std::error::Error>> {
            let node_id = k.parse().map_err(|e| format!("Failed to parse node_id {}: {}", k, e))?;
            let feature_id = v.as_u64()
                .ok_or(format!("Feature ID for node {} is not a number", k))?
                as u32;
            Ok((node_id, feature_id))
        })
        .collect::<Result<HashMap<_, _>, _>>()?;
    
    let edge_feature_map: HashMap<(u64, u64), u32> = json["edge_feature_map"]
        .as_object()
        .ok_or("edge_feature_map not found or not an object")?
        .iter()
        .map(|(k, v)| -> Result<((u64, u64), u32), Box<dyn std::error::Error>> {
            let parts: Vec<&str> = k.split('_').collect();
            if parts.len() != 2 {
                return Err(format!("Invalid edge key format: {}", k).into());
            }
            let n1 = parts[0].parse().map_err(|e| format!("Failed to parse node1 in {}: {}", k, e))?;
            let n2 = parts[1].parse().map_err(|e| format!("Failed to parse node2 in {}: {}", k, e))?;
            let feature_id = v.as_u64()
                .ok_or(format!("Feature ID for edge {} is not a number", k))?
                as u32;
            Ok(((n1, n2), feature_id))
        })
        .collect::<Result<HashMap<_, _>, _>>()?;
    
    // Create a minimal model for legacy format (will need retraining)
    elog!("Warning: Loaded legacy JSON format - model will need retraining");
    let trainer = Trainer::new();
    let model = trainer.max_iter(1)?.train(&[], FeatureProvider::new());
    
    Ok(CRFModel {
        model,
        node_feature_map,
        edge_feature_map,
    })
}

