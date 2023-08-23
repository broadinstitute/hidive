use std::collections::BTreeMap;

use parquet::data_type::AsBytes;

use needletail::Sequence;
use needletail::sequence::complement;

use crate::record::Record;

type KmerGraph = BTreeMap<Vec<u8>, Record>;
type Links = BTreeMap<Vec<u8>, Vec<Vec<u8>>>;

/// Represents a linked de Bruijn graph with a k-mer size specified at construction time.
#[derive(Debug)]
pub struct LdBG {
    pub kmers: KmerGraph,
    pub fw_links: Links,
    pub rc_links: Links
}

impl LdBG {
    /// Create a de Bruijn graph from a list of sequences.
    pub fn from_sequences(k: usize, fwd_seqs: &Vec<Vec<u8>>) -> Self {
        let kmers = Self::build_graph(k, fwd_seqs);
        let fw_links = Self::build_links(k, fwd_seqs, &kmers, false);
        let rc_links = Self::build_links(k, fwd_seqs, &kmers, true);

        LdBG {
            kmers,
            fw_links,
            rc_links
        }
    }

    /// Add a k-mer, the preceeding base, and following base to the graph
    fn add_record_to_graph(graph: &mut KmerGraph, fw_kmer: &[u8], fw_prev_base: u8, fw_next_base: u8) {
        let rc_kmer_vec = fw_kmer.reverse_complement();
        let rc_kmer = rc_kmer_vec.as_bytes();

        let (cn_kmer, can_prev_base, can_next_base) = 
            if fw_kmer < rc_kmer {
                (fw_kmer, fw_prev_base, fw_next_base)
            } else {
                (rc_kmer, complement(fw_next_base), complement(fw_prev_base))
            };

        // If it's not already there, insert k-mer and empty record into k-mer map.
        if !graph.contains_key(cn_kmer) {
            graph.insert(cn_kmer.to_owned(), Record::new(0, None));
        }

        // Increment the canonical k-mer's coverage.
        graph.get_mut(cn_kmer).unwrap().increment_coverage();

        // Set incoming edge for canonical k-mer.
        graph.get_mut(cn_kmer).unwrap().set_incoming_edge(can_prev_base);

        // Set outgoing edge for canonical k-mer.
        graph.get_mut(cn_kmer).unwrap().set_outgoing_edge(can_next_base);
    }

    /// Build a de Bruijn graph from a vector of sequences.
    fn build_graph(k: usize, fwd_seqs: &Vec<Vec<u8>>) -> KmerGraph {
        let mut graph: KmerGraph = KmerGraph::new();

        // Iterate over sequences
        for fwd_seq in fwd_seqs {
            // Iterate over k-mers
            for i in 0..fwd_seq.len()-k+1 {
                let fw_kmer = &fwd_seq[i..i+k];
        
                let prev = fwd_seq.get(i.wrapping_sub(1));
                let fw_prev_base = *prev.unwrap_or(&b'.');

                let next = fwd_seq.get(i+k);
                let fw_next_base = *next.unwrap_or(&b'.');

                Self::add_record_to_graph(&mut graph, fw_kmer, fw_prev_base, fw_next_base);
            }
        }

        graph
    }

    /// Add all junction choices from a given sequence.
    fn add_record_to_links(links: &mut Links, fwd_seq: &Vec<u8>, k: usize, graph: &KmerGraph) {
        // Iterate over k-mers to find junctions
        let mut anchor_kmer = &fwd_seq[0..k];

        for i in 0..fwd_seq.len()-k+1 {
            let fw_kmer = &fwd_seq[i..i+k];

            if !LdBG::has_junction(&graph, fw_kmer, false) && !LdBG::has_junction(&graph, fw_kmer, true) {
                anchor_kmer = &fwd_seq[i..i+k];
            }

            if LdBG::has_junction(&graph, fw_kmer, true) {
                let mut junctions = Vec::new();

                // Iterate over k-mers after the current junction and record further junction choices
                for j in i..fwd_seq.len()-k+1 {
                    let next_kmer = &fwd_seq[j..j+k];

                    let has_junction = LdBG::has_junction(&graph, next_kmer, true);
                    if has_junction {
                        junctions.push(fwd_seq[j+k]);
                    }
                }

                let cn_kmer_vec = LdBG::canonicalize_kmer(anchor_kmer);
                let cn_kmer = cn_kmer_vec.as_bytes();

                if !links.contains_key(cn_kmer) {
                    links.insert(cn_kmer.to_owned(), Vec::new());
                }

                links.get_mut(cn_kmer).unwrap().push(junctions);
            }
        }
    }

    /// Build the links for a de Bruijn graph from a vector of sequences.
    fn build_links(k: usize, fwd_seqs: &Vec<Vec<u8>>, graph: &KmerGraph, is_reverse: bool) -> Links {
        let mut links = Links::new();

        // Iterate over sequences again to add links
        for fwd_seq in fwd_seqs {
            let seq = if is_reverse { fwd_seq.clone().reverse_complement() } else { fwd_seq.clone() };
            LdBG::add_record_to_links(&mut links, &seq, k, graph);
        }

        links
    }

    /// Get the canonical (lexicographically-lowest) version of a k-mer.
    fn canonicalize_kmer(kmer: &[u8]) -> Vec<u8> {
        std::cmp::min(kmer, &kmer.reverse_complement()).to_vec()
    }

    /// Check if the given k-mer represents a junction (in the orientation of the given k-mer).
    fn has_junction(graph: &KmerGraph, kmer: &[u8], outgoing: bool) -> bool {
        let cn_kmer_vec = LdBG::canonicalize_kmer(kmer).to_owned();
        let cn_kmer = cn_kmer_vec.as_bytes();

        if let Some(r) = graph.get(cn_kmer) {
            if kmer == cn_kmer {
                return r.out_degree() > 1 && outgoing || r.in_degree() > 1 && !outgoing;
            } else {
                return r.in_degree() > 1 && outgoing || r.out_degree() > 1 && !outgoing;
            }
        }

        false
    }

    /// Starting at a given k-mer, get the next k-mer (or return None if there isn't a single outgoing edge).
    fn next_kmer(&self, kmer: &[u8], links_in_scope: &mut Vec<Vec<u8>>) -> Option<Vec<u8>> {
        let cn_kmer_vec = LdBG::canonicalize_kmer(kmer).to_owned();
        let cn_kmer = cn_kmer_vec.as_bytes();

        let r = self.kmers.get(cn_kmer)?;

        let next_base = if cn_kmer == kmer {
            match r.out_degree() {
                0 => return None,
                1 => r.outgoing_edges()[0],
                _ => {
                    let consensus_junction_choice = links_in_scope.get(0)?[0];
                    if r.outgoing_edges().contains(&consensus_junction_choice) {
                        links_in_scope.iter_mut().for_each(|link| { link.remove(0); });
                        links_in_scope.retain(|link| !link.is_empty());
                        consensus_junction_choice
                    } else {
                        return None;
                    }
                }
            }
        } else {
            match r.in_degree() {
                0 => return None,
                1 => complement(r.incoming_edges()[0]),
                _ => return None,
            }
        };

        let next_kmer = [&kmer[1..kmer.len()], &[next_base]].concat();

        Some(next_kmer)
    }

    /// Starting at a given k-mer, get the previous k-mer (or return None if there isn't a single incoming edge).
    fn prev_kmer(&self, kmer: &[u8], links_in_scope: &mut Vec<Vec<u8>>) -> Option<Vec<u8>> {
        let cn_kmer_vec = LdBG::canonicalize_kmer(kmer).to_owned();
        let cn_kmer = cn_kmer_vec.as_bytes();

        let r = self.kmers.get(cn_kmer)?;

        let prev_base = if cn_kmer == kmer {
            match r.in_degree() {
                0 => return None,
                1 => r.incoming_edges()[0],
                _ => return None,
            }
        } else {
            match r.out_degree() {
                0 => return None,
                1 => complement(r.outgoing_edges()[0]),
                _ => {
                    let consensus_junction_choice = links_in_scope.get(0)?[0];
                    if r.outgoing_edges().contains(&consensus_junction_choice) {
                        links_in_scope.iter_mut().for_each(|link| { link.remove(0); });
                        links_in_scope.retain(|link| !link.is_empty());
                        complement(consensus_junction_choice)
                    } else {
                        return None;
                    }
                }
            }
        };

        let prev_kmer = [&[prev_base], &kmer[0..kmer.len()-1]].concat();

        Some(prev_kmer)
    }

    /// Starting at a given k-mer, assemble a contig.
    pub fn assemble(&self, kmer: &[u8]) -> Vec<u8> {
        let mut contig: Vec<u8> = kmer.to_vec();

        self.assemble_forward(&mut contig, kmer.to_vec());
        self.assemble_backward(&mut contig, kmer.to_vec());

        contig
    }

    /// Assemble a contig in the forward direction.
    fn assemble_forward(&self, contig: &mut Vec<u8>, start_kmer: Vec<u8>) {
        let mut links_in_scope: Vec<Vec<u8>> = Vec::new();
        let mut last_kmer = start_kmer.clone();
        loop {
            self.update_links(&mut links_in_scope, &last_kmer);

            match self.next_kmer(&last_kmer, &mut links_in_scope) {
                Some(this_kmer) => {
                    contig.push(this_kmer[this_kmer.len() - 1]);
                    last_kmer = this_kmer;
                },
                None => {
                    break;
                }
            }
        }
    }

    /// Assemble a contig in the backward direction.
    fn assemble_backward(&self, contig: &mut Vec<u8>, start_kmer: Vec<u8>) {
        let mut links_in_scope: Vec<Vec<u8>> = Vec::new();
        let mut last_kmer = start_kmer.clone();
        loop {
            self.update_links(&mut links_in_scope, &last_kmer);

            match self.prev_kmer(&last_kmer, &mut links_in_scope) {
                Some(this_kmer) => {
                    contig.insert(0, this_kmer[0]);
                    last_kmer = this_kmer;
                },
                None => {
                    break;
                }
            }
        }
    }

    /// Update available links.
    fn update_links(&self, links_in_scope: &mut Vec<Vec<u8>>, last_kmer: &Vec<u8>) {
        let cn_kmer_vec = LdBG::canonicalize_kmer(last_kmer.as_bytes()).to_owned();
        let cn_kmer = cn_kmer_vec.as_bytes();

        if self.fw_links.contains_key(cn_kmer.as_bytes()) {
            let mut new_links_in_scope = self.fw_links.get(cn_kmer.as_bytes()).unwrap().clone();
            links_in_scope.append(&mut new_links_in_scope);
        }

        if self.rc_links.contains_key(cn_kmer.as_bytes()) {
            let mut new_links_in_scope = self.rc_links.get(cn_kmer.as_bytes()).unwrap().clone();
            links_in_scope.append(&mut new_links_in_scope);
        }
    }

    fn print_kmer(kmer: &[u8], prev_base: u8, next_base: u8, record: Option<&Record>) {
        println!("{} {} {} {}",
            std::char::from_u32(prev_base as u32).unwrap_or('.'),
            std::str::from_utf8(kmer).unwrap(),
            std::char::from_u32(next_base as u32).unwrap_or('.'),
            record.map_or(String::from(""), |r| format!("{}", r))
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::edges::Edges;

    fn get_test_genome() -> Vec<u8> {
        "ACTGATTTCGATGCGATGCGATGCCACGGTGG".as_bytes().to_vec()
    }

    // fn get_test_read() -> Vec<u8> {
    //     "TTTCGATGCGATGCGATGCCACG".as_bytes().to_vec()
    // }

    #[test]
    fn test_from_sequences() {
        let genome = get_test_genome();
        let fwd_seqs = vec!(genome);

        let g = LdBG::from_sequences(5, &fwd_seqs);

        let mut exp_graph = KmerGraph::new();
        exp_graph.insert(b"AAATC".to_vec(), Record::new(1, Some(Edges::from_string("..g.A...".to_string()))));
        exp_graph.insert(b"AATCA".to_vec(), Record::new(1, Some(Edges::from_string("a.....G.".to_string()))));
        exp_graph.insert(b"ACCGT".to_vec(), Record::new(1, Some(Edges::from_string(".c....G.".to_string()))));
        exp_graph.insert(b"ACTGA".to_vec(), Record::new(1, Some(Edges::from_string(".......T".to_string()))));
        exp_graph.insert(b"ATCAG".to_vec(), Record::new(1, Some(Edges::from_string("a......T".to_string()))));
        exp_graph.insert(b"ATCGA".to_vec(), Record::new(1, Some(Edges::from_string(".c..A...".to_string()))));
        exp_graph.insert(b"ATCGC".to_vec(), Record::new(2, Some(Edges::from_string(".c..A...".to_string()))));
        exp_graph.insert(b"ATGCC".to_vec(), Record::new(1, Some(Edges::from_string("..g.A...".to_string()))));
        exp_graph.insert(b"ATGCG".to_vec(), Record::new(2, Some(Edges::from_string("..g.A...".to_string()))));
        exp_graph.insert(b"ATTTC".to_vec(), Record::new(1, Some(Edges::from_string("..g...G.".to_string()))));
        exp_graph.insert(b"CACCG".to_vec(), Record::new(1, Some(Edges::from_string(".c.....T".to_string()))));
        exp_graph.insert(b"CACGG".to_vec(), Record::new(1, Some(Edges::from_string(".c.....T".to_string()))));
        exp_graph.insert(b"CATCG".to_vec(), Record::new(3, Some(Edges::from_string("..g.AC..".to_string()))));
        exp_graph.insert(b"CCACC".to_vec(), Record::new(1, Some(Edges::from_string("......G.".to_string()))));
        exp_graph.insert(b"CCACG".to_vec(), Record::new(1, Some(Edges::from_string("..g...G.".to_string()))));
        exp_graph.insert(b"CGAAA".to_vec(), Record::new(1, Some(Edges::from_string("...t...T".to_string()))));
        exp_graph.insert(b"GATGC".to_vec(), Record::new(3, Some(Edges::from_string(".c...CG.".to_string()))));
        exp_graph.insert(b"GCCAC".to_vec(), Record::new(1, Some(Edges::from_string("...t..G.".to_string()))));
        exp_graph.insert(b"TCGAA".to_vec(), Record::new(1, Some(Edges::from_string("a...A...".to_string()))));
        exp_graph.insert(b"TCGCA".to_vec(), Record::new(2, Some(Edges::from_string("a......T".to_string()))));
        exp_graph.insert(b"TGCCA".to_vec(), Record::new(1, Some(Edges::from_string("a....C..".to_string()))));

        for (kmer, record) in g.kmers {
            assert!(exp_graph.contains_key(kmer.as_bytes()));

            let exp_record = exp_graph.get(kmer.as_bytes()).unwrap();
            assert!(record.coverage() == exp_record.coverage());
            assert!(record.incoming_edges() == exp_record.incoming_edges());
            assert!(record.outgoing_edges() == exp_record.outgoing_edges());
        }
    }

    #[test]
    fn test_assemble() {
        let fw_genome = get_test_genome();
        let rc_genome = fw_genome.reverse_complement();

        let fwd_seqs = vec!(fw_genome.clone());
        let g = LdBG::from_sequences(5, &fwd_seqs);

        // assembly outside cycle should recapitulate entire genome
        assert!(fw_genome == g.assemble(b"ACTGA"));
        assert!(fw_genome == g.assemble(b"TTCGA"));
        assert!(fw_genome == g.assemble(b"TGCCA"));
        assert!(fw_genome == g.assemble(b"GGTGG"));

        // assembly outside cycle should recapitulate entire genome
        assert!(rc_genome == g.assemble(b"ACTGA".to_vec().reverse_complement().as_bytes()));
        assert!(rc_genome == g.assemble(b"TTCGA".to_vec().reverse_complement().as_bytes()));
        assert!(rc_genome == g.assemble(b"TGCCA".to_vec().reverse_complement().as_bytes()));
        assert!(rc_genome == g.assemble(b"GGTGG".to_vec().reverse_complement().as_bytes()));

        // assembly inside cycle should not recapitulate entire genome
        assert!(fw_genome != g.assemble(b"TGCGA"));
        assert!(rc_genome != g.assemble(b"TGCGA".to_vec().reverse_complement().as_bytes()));
        assert!(fw_genome != g.assemble(b"CGATG"));
        assert!(rc_genome != g.assemble(b"CGATG".to_vec().reverse_complement().as_bytes()));
    }
}