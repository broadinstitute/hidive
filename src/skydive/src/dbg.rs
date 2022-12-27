use petgraph::{graph::Graph, stable_graph::NodeIndex, stable_graph::EdgeIndex};

use debruijn::dna_string::DnaString;
use debruijn::{kmer::*, Vmer};

use std::collections::HashMap;

/// Represents a DeBruijn graph with a built-in and fixed k-mer size of 15.
#[derive(Debug)]
pub struct DeBruijnGraph {
    graph: Graph<Kmer15, u32>,
    cov: HashMap<Kmer15, u32>,
    idx: HashMap<Kmer15, NodeIndex>
}

impl DeBruijnGraph {
    /// Create an empty de Bruijn graph.
    pub fn new() -> Self {
        DeBruijnGraph {
            graph: Graph::<Kmer15, u32>::new(),
            cov: HashMap::<Kmer15, u32>::new(),
            idx: HashMap::<Kmer15, NodeIndex>::new()
        }
    }

    /// Add a vector of DnaString sequences to the graph as k-mers and edges.
    pub fn add_all(&mut self, seqs: &Vec<DnaString>) {
        seqs.iter().for_each(|s| {
            self.add(s);
        });
    }

    /// Add a single DnaString sequence to the graph as k-mers and edges.
    pub fn add(&mut self, s: &DnaString) {
        let first_kmer = s.first_kmer::<Kmer15>();
        let mut n1 = self.add_vertex(first_kmer);

        for cur_kmer in s.iter_kmers::<Kmer15>().skip(1) {
            let n2 = self.add_vertex(cur_kmer);
            
            self.add_edge(n1, n2, 0);

            n1 = n2;
        }
    }

    /// Check if a given k-mer exists in the graph.
    pub fn has_vertex(&mut self, v: &Kmer15) -> bool {
        self.idx.contains_key(v)
    }

    /// Get the internal index associated with a particular k-mer.
    pub fn get_index(&mut self, v: &Kmer15) -> Result<&NodeIndex, String> {
        if !self.has_vertex(v) {
            return Err("Vertex not found".to_string());
        }

        return Ok(self.idx.get(v).unwrap())
    }

    /// Get coverage of k-mer in graph
    pub fn get_coverage(&mut self, v: &Kmer15) -> u32 {
        if self.has_vertex(v) {
            return *self.cov.get(v).unwrap();
        }

        0
    }

    /// Add a single k-mer to the graph.
    pub fn add_vertex(&mut self, v: Kmer15) -> NodeIndex {
        if self.has_vertex(&v) {
            self.cov.insert(v, self.cov.get(&v).unwrap() + 1);
            return *self.get_index(&v).unwrap();
        }

        let ni = self.graph.add_node(v);

        self.cov.insert(v, 1);
        self.idx.insert(v, ni);

        ni
    }

    /// Add a single edge between two nodes to the graph.
    pub fn add_edge(&mut self, n1: NodeIndex, n2: NodeIndex, e: u32) -> EdgeIndex {
        self.graph.add_edge(n1, n2, e)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn get_dna_string() -> DnaString {
        DnaString::from_dna_string("GCTATAGCATATGCTGTCATGCAA")
    }

    #[test]
    fn test_add_all() {
        let seqs = vec![
            DnaString::from_dna_string("ATCGTAGCTCCCACT"),
            DnaString::from_dna_string("GGCATAGCTCTAGCA"),
            DnaString::from_dna_string("CACCAACTGATCAGA"),
        ];

        let mut dbg = DeBruijnGraph::new();

        assert_eq!(0, dbg.graph.node_count());

        dbg.add_all(&seqs);

        assert_eq!(3, dbg.graph.node_count());
    }

    #[test]
    fn test_add() {
        let mut dbg = DeBruijnGraph::new();

        assert_eq!(0, dbg.graph.node_count());
        assert_eq!(0, dbg.idx.len());

        for (i, fw_kmer) in get_dna_string().iter_kmers::<Kmer15>().enumerate() {
            dbg.add_vertex(fw_kmer);

            assert_eq!(i+1, dbg.graph.node_count());
            assert_eq!(i+1, dbg.idx.len());
        }
    }

    #[test]
    fn test_count_added_kmers() {
        let kmer_count_1 = DnaString::from_dna_string("CACCAACTGATCAGA").first_kmer();
        let kmer_count_2 = DnaString::from_dna_string("ATCGTAGCTCCCACT").first_kmer();

        let mut dbg = DeBruijnGraph::new();

        assert_eq!(0, dbg.graph.node_count());

        dbg.add_vertex(kmer_count_1);
        dbg.add_vertex(kmer_count_2);
        dbg.add_vertex(kmer_count_2);

        assert_eq!(2, dbg.graph.node_count());

        assert_eq!(1, dbg.get_coverage(&kmer_count_1));
        assert_eq!(2, dbg.get_coverage(&kmer_count_2));
    }

    #[test]
    fn test_has_vertex() {
        let mut dbg = DeBruijnGraph::new();

        let good_vertex = get_dna_string().first_kmer::<Kmer15>();
        let bad_vertex = get_dna_string().last_kmer::<Kmer15>();

        dbg.add_vertex(good_vertex);

        assert_eq!(dbg.has_vertex(&good_vertex), true);
        assert_eq!(dbg.has_vertex(&bad_vertex), false);
    }

    #[test]
    fn test_get_index() {
        let mut dbg = DeBruijnGraph::new();

        for fw_kmer in get_dna_string().iter_kmers::<Kmer15>() {
            let index_in = dbg.add_vertex(fw_kmer);
            let index_out = dbg.graph.node_indices().find(|i| dbg.graph[*i] == fw_kmer).unwrap();

            assert_eq!(index_in, index_out);
        }
    }

    #[test]
    fn test_add_vertex() {
        let mut dbg = DeBruijnGraph::new();

        let first_vertex = get_dna_string().first_kmer::<Kmer15>();
        let last_vertex = get_dna_string().last_kmer::<Kmer15>();

        assert_eq!(0, dbg.graph.node_count());

        dbg.add_vertex(first_vertex);

        assert_eq!(1, dbg.graph.node_count());

        dbg.add_vertex(last_vertex);

        assert_eq!(2, dbg.graph.node_count());
    }

    #[test]
    fn test_add_edge() {
        let mut dbg = DeBruijnGraph::new();

        let s = get_dna_string();
        let mut iter = s.iter_kmers::<Kmer15>();

        let first_kmer = iter.next().unwrap();
        let second_kmer = iter.next().unwrap();

        let n1 = dbg.add_vertex(first_kmer);
        let n2 = dbg.add_vertex(second_kmer);

        assert_eq!(0, dbg.graph.edge_count());

        dbg.add_edge(n1, n2, 0);

        assert_eq!(1, dbg.graph.edge_count());
    }
}