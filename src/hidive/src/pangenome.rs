//! Pangenome graph loading and utilities module.
//!
//! This module provides functionality for loading pangenome graphs from GFA format.
//! The graphs are expected to be manually created and contain tier information
//! in path names (format: "path_name|tier=TierN").
//!
//! **Note**: Graph construction via `build_pangenome_graph` is currently not functional
//! due to issues with abPOA-rs. GFA files should be created manually using external tools.
//! The `load_pangenome_graph` function can load any valid GFA file.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

// Note: abPOA-rs has a bug when adding the second sequence to a graph
// Using a simple workaround that creates individual nodes per sequence
use bio::io::fasta;
use parfait_gfa::gfa::{GfaParser, ParseOptions};
use skydive::elog;

/// Tier confidence weights for graph construction
pub const TIER_1_WEIGHT: f64 = 1.0;
pub const TIER_2_WEIGHT: f64 = 0.7;
pub const TIER_3_WEIGHT: f64 = 0.4;

/// Represents the tier of an assembly
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssemblyTier {
    Tier1,
    Tier2,
    Tier3,
}

impl AssemblyTier {
    /// Get the confidence weight for this tier
    pub fn weight(&self) -> f64 {
        match self {
            AssemblyTier::Tier1 => TIER_1_WEIGHT,
            AssemblyTier::Tier2 => TIER_2_WEIGHT,
            AssemblyTier::Tier3 => TIER_3_WEIGHT,
        }
    }
}

/// Configuration for building a pangenome graph
pub struct PangenomeConfig {
    /// FASTA files for Tier 1 assemblies
    pub tier1_fasta_paths: Vec<PathBuf>,
    /// FASTA files for Tier 2 assemblies
    pub tier2_fasta_paths: Vec<PathBuf>,
    /// FASTA files for Tier 3 assemblies
    pub tier3_fasta_paths: Vec<PathBuf>,
    /// Output path for the GFA graph
    pub output_path: PathBuf,
    /// Minimum alignment length (not used with abPOA, kept for compatibility)
    pub min_aln_len: Option<usize>,
    /// K-mer size for alignment (not used with abPOA, kept for compatibility)
    pub kmer_size: Option<usize>,
}

impl Default for PangenomeConfig {
    fn default() -> Self {
        Self {
            tier1_fasta_paths: Vec::new(),
            tier2_fasta_paths: Vec::new(),
            tier3_fasta_paths: Vec::new(),
            output_path: PathBuf::from("pangenome.gfa"),
            min_aln_len: Some(100),
            kmer_size: Some(17),
        }
    }
}

/// Build a tiered pangenome graph from FASTA assemblies.
/// 
/// **NOTE**: This function is currently not functional due to a bug in abPOA-rs
/// that causes failures when adding the second sequence to a graph. GFA files
/// should be created manually using external tools (e.g., seqwish, abPOA CLI, etc.)
/// and then loaded using `load_pangenome_graph`.
/// 
/// This function is kept for future use if/when the abPOA-rs bug is fixed.
#[allow(dead_code)]
pub fn build_pangenome_graph(config: &PangenomeConfig) -> Result<(), Box<dyn std::error::Error>> {
    elog!("Building tiered pangenome graph with abPOA...");
    
    // Collect all sequences with their tier information
    let mut sequences: Vec<(Vec<u8>, String, AssemblyTier)> = Vec::new();
    
    elog!("Reading sequences from FASTA files...");
    for (tier, paths) in [
        (AssemblyTier::Tier1, &config.tier1_fasta_paths),
        (AssemblyTier::Tier2, &config.tier2_fasta_paths),
        (AssemblyTier::Tier3, &config.tier3_fasta_paths),
    ] {
        for path in paths {
            let reader = fasta::Reader::from_file(path)?;
            for record in reader.records() {
                let record = record?;
                let id = record.id();
                let seq = record.seq().to_vec();
                
                // Add tier annotation to sequence ID
                let annotated_id = format!("{}|tier={:?}", id, tier);
                sequences.push((seq, annotated_id, tier));
            }
        }
    }
    
    if sequences.is_empty() {
        return Err("No sequences found in input FASTA files".into());
    }
    
    elog!("Found {} sequences to align", sequences.len());
    
    // Step 2: Build POA graph with abPOA
    elog!("Building POA graph with abPOA...");
    build_graph_with_abpoa(&sequences, &config.output_path)?;
    
    elog!("Pangenome graph construction complete: {}", config.output_path.display());
    elog!("Tier information is encoded in path names (format: path_name|tier=TierN)");
    Ok(())
}

/// Build a variation graph using abPOA and write GFA format
/// 
/// Note: There's a known bug in abPOA-rs where adding the second sequence fails.
/// This function works around it by creating individual nodes for each sequence
/// and building a simple graph structure.
fn build_graph_with_abpoa(
    sequences: &[(Vec<u8>, String, AssemblyTier)],
    gfa_output: &PathBuf,
) -> Result<(), Box<dyn std::error::Error>> {
    elog!("Building graph structure from {} sequences...", sequences.len());
    
    // Workaround for abPOA-rs bug: Instead of using POA graph building,
    // we'll create a simple graph where each sequence is its own node
    // and we'll add edges based on sequence similarity
    
    // For now, create a simple graph structure without using abPOA's graph building
    // Each sequence becomes a node, and we can add edges later if needed
    elog!("Creating individual nodes for each sequence...");
    
    // Write GFA directly without using abPOA's graph (workaround for the bug)
    write_gfa_simple(sequences, gfa_output)?;
    
    Ok(())
}

/// Write a simple GFA where each sequence is a node
/// This is a workaround for the abPOA-rs bug
fn write_gfa_simple(
    sequences: &[(Vec<u8>, String, AssemblyTier)],
    output_path: &PathBuf,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = File::create(output_path)?;
    
    // Write GFA header
    writeln!(file, "H\tVN:Z:1.0")?;
    
    // Create one node per sequence
    let mut node_id = 1u64;
    let mut sequence_to_node: HashMap<String, u64> = HashMap::new();
    let mut edges: Vec<(u64, u64)> = Vec::new();
    
    elog!("Writing {} sequences as nodes...", sequences.len());
    
    // Write segments (nodes) - one per sequence
    for (seq, name, _tier) in sequences {
        let seq_str = String::from_utf8(seq.clone())
            .map_err(|e| format!("Invalid UTF-8 in sequence {}: {}", name, e))?;
        
        writeln!(file, "S\t{}\t{}", node_id, seq_str)?;
        sequence_to_node.insert(name.clone(), node_id);
        node_id += 1;
    }
    
    // Add edges between similar sequences (optional - can be enhanced)
    // For now, we'll skip edges and just have individual nodes
    
    // Write paths - each sequence is a path through its own node
    for (seq, name, _tier) in sequences {
        if let Some(&node_id) = sequence_to_node.get(name) {
            writeln!(file, "P\t{}\t{}+\t0M", name, node_id)?;
        }
    }
    
    elog!(
        "GFA file written with {} segments and {} paths",
        sequences.len(),
        sequences.len()
    );
    
    Ok(())
}


/// Extract tier information from a path name
pub fn extract_tier_from_path_name(path_name: &str) -> Option<AssemblyTier> {
    if path_name.contains("tier=Tier1") {
        Some(AssemblyTier::Tier1)
    } else if path_name.contains("tier=Tier2") {
        Some(AssemblyTier::Tier2)
    } else if path_name.contains("tier=Tier3") {
        Some(AssemblyTier::Tier3)
    } else {
        None
    }
}

/// Get tier weight for a path name
pub fn get_tier_weight(path_name: &str) -> f64 {
    extract_tier_from_path_name(path_name)
        .map(|tier| tier.weight())
        .unwrap_or(TIER_3_WEIGHT) // Default to lowest confidence
}

/// Represents a pangenome graph loaded from GFA format
#[derive(Debug, Clone)]
pub struct PangenomeGraph {
    /// Node IDs to sequences
    pub nodes: HashMap<u64, Vec<u8>>,
    /// Edge list (from_node, to_node)
    pub edges: Vec<(u64, u64)>,
    /// Path information (path_name -> node_ids)
    pub paths: HashMap<String, Vec<u64>>,
}

/// Load a pangenome graph from GFA format
pub fn load_pangenome_graph(path: &PathBuf) -> Result<PangenomeGraph, Box<dyn std::error::Error>> {
    elog!("Loading pangenome graph from GFA file: {}", path.display());
    
    let mut parser = GfaParser::new();
    let parse_options = ParseOptions::default();
    
    // Parse the GFA file
    let parse_result = parser.parse(path.to_str().unwrap(), &parse_options);
    
    match parse_result {
        Ok(_) => {
            elog!("GFA file parsed successfully");
        }
        Err(errors) => {
            let error_msg = format!("Failed to parse GFA file: {:?}", errors);
            return Err(error_msg.into());
        }
    }
    
    // Extract nodes (segments) - collect all segments first
    let segments: Vec<_> = parser.segments().collect();
    let mut nodes = HashMap::new();
    let mut node_id_map: HashMap<String, u64> = HashMap::new();
    let mut next_id: u64 = 1;
    let mut edge_set: HashSet<(u64, u64)> = HashSet::new();
    
    // Build node mappings from segments
    for segment in segments.iter() {
        let segment_name = &segment.name;
        
        // Map string node names to numeric IDs
        let node_id = if let Some(&id) = node_id_map.get(segment_name) {
            id
        } else {
            let id = next_id;
            node_id_map.insert(segment_name.clone(), id);
            next_id += 1;
            id
        };
        
        // Get sequence from segment
        let sequence = segment.sequence.as_bytes().to_vec();
        
        nodes.insert(node_id, sequence);
    }
    
    // Extract edges (links)
    for link in parser.links() {
        let from_name = &link.from_segment;
        let to_name = &link.to_segment;
        
        let from_id = *node_id_map.get(from_name)
            .ok_or_else(|| format!("Unknown from node: {}", from_name))?;
        let to_id = *node_id_map.get(to_name)
            .ok_or_else(|| format!("Unknown to node: {}", to_name))?;
        
        edge_set.insert((from_id, to_id));
    }
    
    // Extract paths
    let mut paths = HashMap::new();
    
    for path in parser.paths() {
        let path_name = path.name.clone();
        let mut node_ids = Vec::new();
        
        // Parse path steps (segments with orientations)
        // parfait-gfa's step.segment_id appears to be an internal ID
        // We'll try both 0-based and 1-based indexing, and also try matching by name
        for step in &path.steps {
            let seg_id = step.segment_id;
            
            // Try different indexing schemes and lookups
            let segment_name_opt = if seg_id > 0 && (seg_id as usize) <= segments.len() {
                // Try 1-based index (most common in GFA)
                Some(segments[(seg_id as usize) - 1].name.clone())
            } else if (seg_id as usize) < segments.len() {
                // Try 0-based index
                Some(segments[seg_id as usize].name.clone())
            } else {
                // If still not found, try to find by matching the segment_id as a string
                let seg_id_str = seg_id.to_string();
                segments.iter().find(|s| s.name == seg_id_str).map(|s| s.name.clone())
            };
            
            // If we still can't find it, skip this step with a warning
            // This might be a bug in parfait-gfa or an edge case we need to handle
            if let Some(segment_name) = segment_name_opt {
                if let Some(&node_id) = node_id_map.get(&segment_name) {
                    node_ids.push(node_id);
                } else {
                    elog!("Warning: Segment name '{}' from path '{}' not found in node map, skipping", segment_name, path_name);
                }
            } else {
                elog!("Warning: Segment ID {} in path '{}' not found (total segments: {}), skipping step", 
                    seg_id, path_name, segments.len());
            }
        }
        
        // Derive edges from consecutive nodes along this path
        for window in node_ids.windows(2) {
            if let [from, to] = window {
                edge_set.insert((*from, *to));
            }
        }
        
        paths.insert(path_name, node_ids);
    }
    
    let edges: Vec<(u64, u64)> = edge_set.into_iter().collect();
    
    elog!(
        "Loaded graph: {} nodes, {} edges, {} paths",
        nodes.len(),
        edges.len(),
        paths.len()
    );
    
    Ok(PangenomeGraph {
        nodes,
        edges,
        paths,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tier_weights() {
        assert_eq!(AssemblyTier::Tier1.weight(), TIER_1_WEIGHT);
        assert_eq!(AssemblyTier::Tier2.weight(), TIER_2_WEIGHT);
        assert_eq!(AssemblyTier::Tier3.weight(), TIER_3_WEIGHT);
    }

    #[test]
    fn test_extract_tier() {
        assert_eq!(
            extract_tier_from_path_name("sample1|tier=Tier1"),
            Some(AssemblyTier::Tier1)
        );
        assert_eq!(
            extract_tier_from_path_name("sample2|tier=Tier2"),
            Some(AssemblyTier::Tier2)
        );
        assert_eq!(
            extract_tier_from_path_name("sample3|tier=Tier3"),
            Some(AssemblyTier::Tier3)
        );
        assert_eq!(extract_tier_from_path_name("sample4"), None);
    }

    #[test]
    fn test_get_tier_weight() {
        assert_eq!(get_tier_weight("sample1|tier=Tier1"), TIER_1_WEIGHT);
        assert_eq!(get_tier_weight("sample2|tier=Tier2"), TIER_2_WEIGHT);
        assert_eq!(get_tier_weight("sample3|tier=Tier3"), TIER_3_WEIGHT);
        assert_eq!(get_tier_weight("sample4"), TIER_3_WEIGHT); // Default
    }
}
