use petgraph::{graph::Graph, stable_graph::NodeIndex, stable_graph::EdgeIndex};

use debruijn::dna_string::DnaString;
use debruijn::{kmer::*, Vmer};

use std::collections::HashMap;

#[derive(Debug)]
pub struct DeBruijnGraph {
    g: Graph<Kmer15, u32>,
    idx: HashMap<Kmer15, NodeIndex>
}

impl DeBruijnGraph {
    pub fn new() -> Self {
        DeBruijnGraph {
            g: Graph::<Kmer15, u32>::new(),
            idx: HashMap::<Kmer15, NodeIndex>::new()
        }
    }

    pub fn add(&mut self, s: &DnaString) {
        let mut last_kmer = s.first_kmer::<Kmer15>();

        for cur_kmer in s.iter_kmers::<Kmer15>().skip(1) {
            let n1 = self.add_vertex(last_kmer);
            let n2 = self.add_vertex(cur_kmer);
            
            self.add_edge(n1, n2, 0);

            last_kmer = cur_kmer;
        }
    }

    pub fn has_vertex(&mut self, v: &Kmer15) -> bool {
        self.idx.contains_key(v)
    }

    pub fn get_index(&mut self, v: &Kmer15) -> Result<&NodeIndex, String> {
        if !self.has_vertex(v) {
            return Err("Vertex not found".to_string());
        }

        return Ok(self.idx.get(v).unwrap())
    }

    pub fn add_vertex(&mut self, v: Kmer15) -> NodeIndex {
        if self.has_vertex(&v) {
            return *self.get_index(&v).unwrap();
        }

        let ni = self.g.add_node(v);

        self.idx.insert(v, ni);

        ni
    }

    pub fn add_edge(&mut self, n1: NodeIndex, n2: NodeIndex, e: u32) -> EdgeIndex {
        self.g.add_edge(n1, n2, e)
    }
}