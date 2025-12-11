use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};

use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::algo::toposort;

pub fn start(gfa_path: &PathBuf) {
    skydive::elog!("Training CRF model from GFA file: {}", gfa_path.display());
    
    // Read the GFA file
    let g = read_gfa(gfa_path).unwrap();

    skydive::elog!("Loaded GFA with {} nodes", g.graph.node_count());
    
    // Topologically sort the graph
    let sorted_nodes = match toposort(&g.graph, None) {
        Ok(nodes) => nodes,
        Err(e) => {
            skydive::elog!("Graph contains a cycle: {:?}", e);
            return;
        }
    };

    // --- Toy longest-path DP ---
    let (best_path, best_score) = longest_path_dag(&g.graph, &sorted_nodes);
    skydive::elog!(
        "Best path node count: {}, score: {}",
        best_path.len(),
        best_score
    );

    let assembled = emit_sequence(&g.graph, &best_path);
    skydive::elog!("Assembled contig length: {}", assembled.len());

    for path in &g.paths {
        let s = seq_for_path(&g, path);
        skydive::elog!("Allele {}: length {}", path.name, s.len());
    }
}

/// Per-node data: sequence and metadata.
#[derive(Debug, Clone)]
pub struct NodeData {
    pub id: String,   // GFA segment ID (e.g. "S1", "seg123")
    pub seq: String,  // nucleotide sequence
    pub len: usize,   // cached length = seq.len()
}

/// Per-edge data: transition weight / score.
#[derive(Debug, Clone, Copy)]
pub struct EdgeData {
    pub weight: f32,  // you can treat this as log-pot or cost later
}

/// One named path (e.g. an HLA allele) through the graph.
#[derive(Debug, Clone)]
pub struct AllelePath {
    pub name: String,          // path name from P line (e.g. "HLA_A_01:01")
    pub nodes: Vec<NodeIndex>, // resolved sequence of nodes
}

/// Pangenome wrapper: graph + paths.
#[derive(Debug)]
pub struct Pangenome {
    pub graph: DiGraph<NodeData, EdgeData>,
    pub paths: Vec<AllelePath>,
}

/// Internal helper: raw path with segment IDs before we resolve them.
#[derive(Debug)]
struct RawPath {
    name: String,
    seg_ids: Vec<String>, // e.g. ["seg1", "seg2", "seg3"]
}

/// Reads a GFA file and constructs a directed graph + allele paths.
///
/// - S lines -> nodes
/// - L lines -> directed edges
/// - P lines -> named paths through the graph
pub fn read_gfa<Pth: AsRef<Path>>(path: Pth) -> io::Result<Pangenome> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut graph: DiGraph<NodeData, EdgeData> = DiGraph::new();
    let mut node_map: HashMap<String, NodeIndex> = HashMap::new();
    let mut raw_paths: Vec<RawPath> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.is_empty() {
            continue;
        }

        match fields[0] {
            // Segment line: S <id> <sequence> ...
            "S" => {
                if fields.len() < 3 {
                    continue; // malformed S line
                }
                let id = fields[1].to_string();
                let seq = fields[2].to_string();
                let len = seq.len();

                let node_data = NodeData { id: id.clone(), seq, len };
                let node_index = graph.add_node(node_data);
                node_map.insert(id, node_index);
            }

            // Link line: L <from> <from_orient> <to> <to_orient> <overlap> [opt...]
            "L" => {
                if fields.len() < 6 {
                    continue; // malformed L line
                }
                let from_id = fields[1];
                let _from_orient = fields[2];
                let to_id = fields[3];
                let _to_orient = fields[4];

                // Overlap is fields[5], but we don't really use it right now.
                // Optional: parse a weight from an optional field (e.g. "cg:f:1.0") if present.
                // For now, default to 1.0 unless you have a convention in the 6th/7th field.
                let weight = 1.0_f32;

                if let (Some(&from), Some(&to)) = (node_map.get(from_id), node_map.get(to_id)) {
                    graph.add_edge(from, to, EdgeData { weight });
                }
            }

            // Path line: P <name> <segment_ids> <overlaps> [opt...]
            "P" => {
                if fields.len() < 4 {
                    continue; // malformed P line
                }
                let name = fields[1].to_string();
                let segment_list = fields[2]; // e.g., "seg1+,seg2-,seg3+"

                let seg_ids: Vec<String> = segment_list
                    .split(',')
                    .filter_map(|token| {
                        if token.is_empty() {
                            return None;
                        }
                        // token is like "seg1+" or "seg2-"
                        let (seg_id, _orient) = token.split_at(token.len().saturating_sub(1));
                        if seg_id.is_empty() {
                            None
                        } else {
                            Some(seg_id.to_string())
                        }
                    })
                    .collect();

                raw_paths.push(RawPath { name, seg_ids });
            }

            _ => {
                // Ignore other record types (H, W, etc.)
            }
        }
    }

    // Resolve RawPath seg_ids -> NodeIndex using node_map
    let mut paths: Vec<AllelePath> = Vec::new();
    for raw in raw_paths {
        let mut node_indices = Vec::new();
        let mut missing = false;

        for seg_id in &raw.seg_ids {
            if let Some(&idx) = node_map.get(seg_id) {
                node_indices.push(idx);
            } else {
                // This path references a segment we never saw.
                // You can choose to skip this path entirely or keep partial.
                eprintln!("Warning: path '{}' references unknown segment '{}'; skipping path.",
                          raw.name, seg_id);
                missing = true;
                break;
            }
        }

        if !missing && !node_indices.is_empty() {
            paths.push(AllelePath {
                name: raw.name,
                nodes: node_indices,
            });
        }
    }

    Ok(Pangenome { graph, paths })
}

use petgraph::Direction;
use petgraph::visit::EdgeRef;
use std::f32;

fn longest_path_dag(
    graph: &DiGraph<NodeData, EdgeData>,
    topo: &[NodeIndex],
) -> (Vec<NodeIndex>, f32) {
    let n_nodes = graph.node_count();
    let mut score = vec![f32::NEG_INFINITY; n_nodes];
    let mut prev: Vec<Option<NodeIndex>> = vec![None; n_nodes];

    // initialize sources: nodes with no incoming edges get base score 0.0
    for &node in topo {
        let in_deg = graph
            .neighbors_directed(node, Direction::Incoming)
            .count();
        if in_deg == 0 {
            score[node.index()] = 0.0;
        }
    }

    // DP in topological order
    for &u in topo {
        let u_score = score[u.index()];
        if u_score == f32::NEG_INFINITY {
            continue; // unreachable
        }

        // simple node "emission": length of this node's sequence
        let emit = graph[u].len as f32;

        for edge in graph.edges(u) {
            let v = edge.target();
            let w = edge.weight().weight;

            let cand = u_score + emit + w;
            if cand > score[v.index()] {
                score[v.index()] = cand;
                prev[v.index()] = Some(u);
            }
        }
    }

    // find best sink: node with no outgoing edges and max score
    let mut best_node: Option<NodeIndex> = None;
    let mut best_score = f32::NEG_INFINITY;

    for &node in topo {
        let out_deg = graph
            .neighbors_directed(node, Direction::Outgoing)
            .count();
        if out_deg == 0 {
            let s = score[node.index()];
            if s > best_score {
                best_score = s;
                best_node = Some(node);
            }
        }
    }

    let end = best_node.expect("No sink node found in DAG");

    // backtrack path
    let mut path = Vec::new();
    let mut cur = Some(end);
    while let Some(nidx) = cur {
        path.push(nidx);
        cur = prev[nidx.index()];
    }
    path.reverse();

    (path, best_score)
}

fn emit_sequence(graph: &DiGraph<NodeData, EdgeData>, path: &[NodeIndex]) -> String {
    let mut seq = String::new();
    for &n in path {
        seq.push_str(&graph[n].seq);
    }
    seq
}

fn seq_for_path(g: &Pangenome, path: &AllelePath) -> String {
    let mut seq = String::new();
    for &n in &path.nodes {
        seq.push_str(&g.graph[n].seq);
    }
    seq
}