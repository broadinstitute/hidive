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

#[cfg(test)]
mod tests {
    use super::*;

    fn get_dna_string() -> DnaString {
        DnaString::from_dna_string("GCTATAGCATATGCTGTCATGCAA")
    }

    #[test]
    fn test_add() {
        let mut graph = DeBruijnGraph::new();

        assert_eq!(0, graph.g.node_count());
        assert_eq!(0, graph.idx.len());

        for (i, fw_kmer) in get_dna_string().iter_kmers::<Kmer15>().enumerate() {
            graph.add_vertex(fw_kmer);

            assert_eq!(i+1, graph.g.node_count());
            assert_eq!(i+1, graph.idx.len());
        }
    }

    #[test]
    fn test_has_vertex() {
        let mut graph = DeBruijnGraph::new();

        let good_vertex = get_dna_string().first_kmer::<Kmer15>();
        let bad_vertex = get_dna_string().last_kmer::<Kmer15>();

        graph.add_vertex(good_vertex);

        assert_eq!(graph.has_vertex(&good_vertex), true);
        assert_eq!(graph.has_vertex(&bad_vertex), false);
    }

    #[test]
    fn test_get_index() {
        let mut graph = DeBruijnGraph::new();

        for fw_kmer in get_dna_string().iter_kmers::<Kmer15>() {
            let index_in = graph.add_vertex(fw_kmer);
            let index_out = graph.g.node_indices().find(|i| graph.g[*i] == fw_kmer).unwrap();

            assert_eq!(index_in, index_out);
        }
    }

    #[test]
    fn test_add_vertex() {
        let mut graph = DeBruijnGraph::new();

        let first_vertex = get_dna_string().first_kmer::<Kmer15>();
        let last_vertex = get_dna_string().last_kmer::<Kmer15>();

        assert_eq!(0, graph.g.node_count());

        graph.add_vertex(first_vertex);

        assert_eq!(1, graph.g.node_count());

        graph.add_vertex(last_vertex);

        assert_eq!(2, graph.g.node_count());
    }

    #[test]
    fn test_add_edge() {
        let mut graph = DeBruijnGraph::new();

        let s = get_dna_string();
        let mut iter = s.iter_kmers::<Kmer15>();

        let first_kmer = iter.next().unwrap();
        let second_kmer = iter.next().unwrap();

        let n1 = graph.add_vertex(first_kmer);
        let n2 = graph.add_vertex(second_kmer);

        assert_eq!(0, graph.g.edge_count());

        graph.add_edge(n1, n2, 0);

        assert_eq!(1, graph.g.edge_count());
    }
}