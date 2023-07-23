use std::collections::BTreeMap;

use indextree::{Arena, NodeId};
use parquet::data_type::AsBytes;

use needletail::Sequence;
use needletail::sequence::complement;

use crate::record::Record;

type KmerGraph = BTreeMap<Vec<u8>, Record>;

/// Represents a linked de Bruijn graph with a k-mer size specified at construction time.
#[derive(Debug)]
pub struct LdBG {
    pub kmers: KmerGraph
}

impl LdBG {
    /// Create a de Bruijn graph from a list of sequences.
    pub fn from_sequences(k: usize, fwd_seqs: Vec<Vec<u8>>) -> Self {
        let graph = Self::build_graph(k, &fwd_seqs);
        // let links = build_links(kmer_size, &fwd_seqs, &graph);

        LdBG {
            kmers: graph
        }
    }

    fn print_kmer(kmer: &[u8], prev_base: u8, next_base: u8, record: Option<&Record>) {
        println!("{} {} {} {}",
            std::char::from_u32(prev_base as u32).unwrap_or('N'),
            std::str::from_utf8(kmer).unwrap(),
            std::char::from_u32(next_base as u32).unwrap_or('N'),
            record.map_or(String::from(""), |r| format!("{}", r))
        );
    }

    fn add_to_graph(graph: &mut KmerGraph, fw_kmer: &[u8], fw_prev_base: u8, fw_next_base: u8) {
        let rc_kmer_vec = fw_kmer.reverse_complement();
        let rc_kmer = rc_kmer_vec.as_bytes();

        let (can_kmer, can_prev_base, can_next_base) = 
            if fw_kmer < rc_kmer {
                (fw_kmer, fw_prev_base, fw_next_base)
            } else {
                (rc_kmer, complement(fw_next_base), complement(fw_prev_base))
            };

        // If it's not already there, insert k-mer and empty record into k-mer map.
        if !graph.contains_key(can_kmer) {
            graph.insert(can_kmer.to_owned(), Record::new(0));
        }

        // Increment the canonical k-mer's coverage.
        graph.get_mut(can_kmer).unwrap().increment_coverage();

        // Set incoming edge for canonical k-mer.
        graph.get_mut(can_kmer).unwrap().set_incoming_edge(can_prev_base);

        // Set outgoing edge for canonical k-mer.
        graph.get_mut(can_kmer).unwrap().set_outgoing_edge(can_next_base);

        // Print some information.
        Self::print_kmer(fw_kmer, fw_prev_base, fw_next_base, None);
        Self::print_kmer(can_kmer, can_prev_base, can_next_base, graph.get(can_kmer));
        println!("");
    }

    fn build_graph(k: usize, fwd_seqs: &Vec<Vec<u8>>) -> KmerGraph {
        let mut graph: KmerGraph = BTreeMap::new();

        // Iterate over sequences
        for fwd_seq in fwd_seqs {
            // Iterate over k-mers
            for i in 0..fwd_seq.len()-k+1 {
                let fw_kmer = &fwd_seq[i..i+k];
        
                let prev = fwd_seq.get(i.wrapping_sub(1));
                let fw_prev_base = *prev.unwrap_or(&b'N');

                let next = fwd_seq.get(i+k);
                let fw_next_base = *next.unwrap_or(&b'N');

                Self::add_to_graph(&mut graph, fw_kmer, fw_prev_base, fw_next_base);
            }
        }

        graph
    }

}

pub fn build_links(k: usize, fwd_seqs: &Vec<Vec<u8>>, graph: &KmerGraph) -> BTreeMap<Vec<u8>, Arena<u8>> {
    let mut junction_tree: BTreeMap<Vec<u8>, Arena<u8>> = BTreeMap::new();

    // Iterate over sequences again to add in the link information.
    for fwd_seq in fwd_seqs {
        println!("{:?}", fwd_seq.len());

        let mut anchor_kmer: Option<&[u8]> = None;

        // Iterate over k-mers.
        for (i, fwd_kmer) in fwd_seq.windows(k).enumerate() {
            // If this k-mer has out-degree > 1, and we have a place to anchor the links, start the link annotation process.
            if graph.get(fwd_kmer).unwrap().out_degree() > 1 && anchor_kmer.is_some() {
                // Create the junction tree.
                if !junction_tree.contains_key(anchor_kmer.unwrap()) {
                    let arena = &mut Arena::new();
                    junction_tree.insert(anchor_kmer.unwrap().to_owned(), arena.to_owned());
                }

                let arena = junction_tree.get_mut(anchor_kmer.unwrap()).unwrap();
                let root = arena.new_node('-' as u8);
                let mut last_junction = root;

                for (j, _) in fwd_seq
                    .windows(k)
                    .enumerate()
                    .skip(i)
                    .filter(|x| {
                        graph.get((*x).1).unwrap().out_degree() > 1
                    }) {

                    let this_junction = arena.new_node(fwd_seq[j + k]);
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

    junction_tree
}

fn print_junction_tree(root: NodeId, arena: &Arena<u8>, depth: usize) {
    for child in root.children(arena) {
        let node = arena.get(child).unwrap();

        println!("{} {}", depth, *node.get() as char);

        print_junction_tree(child, arena, depth + 1);
    }
}

#[cfg(test)]
mod tests {
    use std::io::Read;

    use super::*;

    fn get_test_genome() -> Vec<u8> {
        "ACTGATTTCGATGCGATGCGATGCCACGGTGG".as_bytes().to_vec()
        /*
         CCACCGTGGCATCGCATCGCATCGAAATCAGT
         */
    }

    fn get_test_read() -> Vec<u8> {
        "TTTCGATGCGATGCGATGCCACG".as_bytes().to_vec()
    }

    #[test]
    fn test_add_all() {
        let genome = get_test_genome();
        let rc = genome.reverse_complement();

        println!("{}", String::from_utf8(genome.to_owned()).unwrap());
        println!("{}", String::from_utf8(rc).unwrap());
        // let read = get_test_read();

        let g = LdBG::from_sequences(5, vec!(genome));

        for (kmer, record) in g.kmers {
            println!("{:?} {}", String::from_utf8(kmer), record);
        }
    }

    /*
AAATC 1 ..g.A...
AATCA 1 a.....G.
ACCGT 1 .c....G.
ACTGA 1 .......T
ATCAG 1 a......T
ATCGA 1 .c..A...
ATCGC 2 .c..A...
ATGCC 1 ..g.A...
ATGCG 2 ..g.A...
ATTTC 1 ..g...G.
CACCG 1 .c.....T
CACGG 1 .c.....T
CATCG 3 ..g.AC..
CCACC 1 ......G.
CCACG 1 ..g...G.
CGAAA 1 ...t...T
GATGC 3 .c...CG.
GCCAC 1 ...t..G.
TCGAA 1 a...A...
TCGCA 2 a......T
TGCCA 1 a....C..
     */
}