use std::collections::BTreeMap;

use indextree::{Arena, NodeId};

use crate::record::Record;

fn print_junction_tree(root: NodeId, arena: &Arena<u8>, depth: usize) {
    for child in root.children(arena) {
        let node = arena.get(child).unwrap();

        println!("{} {}", depth, *node.get() as char);

        print_junction_tree(child, arena, depth + 1);
    }
}

// #[derive(Debug)]
// struct KmerRecordMap(BTreeMap<Vec<u8>, Record>);

/// Represents a linked de Bruijn graph with a variable k-mer size specified at construction time.
#[derive(Debug)]
pub struct LdBG<const K: usize> {
    pub kmers: BTreeMap<Vec<u8>, Record>
}

impl<const K: usize> LdBG<K> {
    /// Create a de Bruijn graph from a list of sequences.
    pub fn from_sequences(fwd_seqs: Vec<Vec<u8>>) -> Self {
        let mut kmers: BTreeMap<Vec<u8>, Record> = BTreeMap::new();

        // Iterate over sequences
        for fwd_seq in &fwd_seqs {
            // Iterate over k-mers
            for (i, fwd_kmer) in fwd_seq.windows(K).enumerate() {
                // If it's not already there, insert k-mer and empty record into k-mer map.
                if !kmers.contains_key(fwd_kmer) {
                    kmers.insert(fwd_kmer.to_owned(), Record::new(0));
                }

                // Increment the k-mer's coverage.
                kmers.get_mut(fwd_kmer).unwrap().increment_coverage();

                // Set incoming edge.
                if i > 0 {
                    kmers.get_mut(fwd_kmer).unwrap().set_incoming_edge(fwd_seq[i - 1]);
                }

                // Set outgoing edge.
                if i + K < fwd_seq.len() {
                    kmers.get_mut(fwd_kmer).unwrap().set_outgoing_edge(fwd_seq[i + K]);
                }
            }
        }

        let mut junction_tree: BTreeMap<Vec<u8>, Arena<u8>> = BTreeMap::new();

        // Iterate over sequences again to add in the link information.
        for fwd_seq in &fwd_seqs {
            println!("{:?}", fwd_seq.len());

            let mut anchor_kmer: Option<&[u8]> = None;

            // Iterate over k-mers.
            for (i, fwd_kmer) in fwd_seq.windows(K).enumerate() {
                // If this k-mer has out-degree > 1, and we have a place to anchor the links, start the link annotation process.
                if kmers.get(fwd_kmer).unwrap().out_degree() > 1 && anchor_kmer.is_some() {
                    // Create the junction tree.
                    if !junction_tree.contains_key(anchor_kmer.unwrap()) {
                        let arena = &mut Arena::new();
                        junction_tree.insert(anchor_kmer.unwrap().to_owned(), arena.to_owned());
                    }

                    let arena = junction_tree.get_mut(anchor_kmer.unwrap()).unwrap();
                    let root = arena.new_node('-' as u8);
                    let mut last_junction = root;

                    for (j, _) in fwd_seq
                        .windows(K)
                        .enumerate()
                        .skip(i)
                        .filter(|x| {
                            kmers.get((*x).1).unwrap().out_degree() > 1
                        }) {

                        // println!("{:?} {:?} {:?} {:?} {} {}", i, String::from_utf8(fwd_kmer.to_vec()), j, String::from_utf8(junction_kmer.to_vec()), fwd_seq[j + K] as char, kmers.get(junction_kmer).unwrap());

                        let this_junction = arena.new_node(fwd_seq[j + K]);
                        last_junction.append(this_junction, arena);
                        last_junction = this_junction;
                    }

                    // println!("{:?}", arena);
                    // print_junction_tree(root, arena, 0);
                }

                anchor_kmer = Some(fwd_kmer);
            }

            // print_junction_tree(root, arena, 0);
        }

        LdBG {
            kmers
        }
    }

}

#[cfg(test)]
mod tests {
    /*
    use parquet::file::writer::SerializedColumnWriter;

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

    // Trying some stuff

    #[test]
    fn sample_test() {
        use color_eyre::{Result};
        use polars::prelude::*;

        let mut df = df!(
            "foo" => [1, 2, 3],
            "bar" => [None, Some("bak"), Some("baz")],
        )
        .unwrap();
        
        let mut file = std::fs::File::create("path.parquet").unwrap();
        ParquetWriter::new(&mut file).finish(&mut df).unwrap();

        println!("{:?}", df);
    }
    */
}