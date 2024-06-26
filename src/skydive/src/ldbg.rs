use std::collections::{HashMap, HashSet};
// use std::sync::Mutex;

use indicatif::ProgressIterator;
use parquet::data_type::AsBytes;

use std::path::PathBuf;
use needletail::Sequence;
use needletail::sequence::complement;

use rayon::iter::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;
use indicatif::ParallelProgressIterator;

use crate::record::Record;
use crate::link::Link;

type KmerGraph = HashMap<Vec<u8>, Record>;
type Links = HashMap<Vec<u8>, HashMap<Link, u16>>;

/// Represents a linked de Bruijn graph with a k-mer size specified at construction time.
#[derive(Debug)]
pub struct LdBG {
    pub name: String,
    pub kmers: KmerGraph,
    pub links: Links,
    pub junctions: HashMap<Vec<u8>, bool>,
}

impl LdBG {
    /// Create a de Bruijn graph (and optional links) from a file path.
    pub fn from_file(name: String, k: usize, seq_path: &PathBuf, build_links: bool) -> Self {
        let reader = bio::io::fasta::Reader::from_file(seq_path).unwrap();
        let all_reads: Vec<bio::io::fasta::Record> = reader.records().map(|r| r.unwrap()).collect();

        let fwd_seqs: Vec<Vec<u8>> = all_reads
            .iter()
            .map(|r| r.seq().as_bytes().to_vec())
            .collect();

        let kmers = Self::build_graph(k, &fwd_seqs);

        let junctions = Self::find_junctions(&kmers);

        let links = match build_links {
            true => Self::build_links(k, &fwd_seqs, &kmers, &junctions),
            false => Links::new()
        };

        LdBG {
            name,
            kmers,
            links,
            junctions
        }
    }

    /// Create a de Bruijn graph (and optional links) from a list of sequences.
    pub fn from_sequences(name: String, k: usize, fwd_seqs: &Vec<Vec<u8>>, build_links: bool) -> Self {
        let kmers = Self::build_graph(k, fwd_seqs);

        let junctions = Self::find_junctions(&kmers);

        let links = match build_links {
            true => Self::build_links(k, fwd_seqs, &kmers, &junctions),
            false => Links::new()
        };

        LdBG {
            name,
            kmers,
            links,
            junctions
        }
    }

    /// Get name of graph.
    pub fn name(&self) -> &String {
        &self.name
    }

    /// Add a k-mer, the preceeding base, and following base to the graph.
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
    fn add_record_to_links(links: &mut Links, fwd_seq: &Vec<u8>, k: usize, graph: &KmerGraph, junctions: &HashMap<Vec<u8>, bool>, outgoing: bool) {
        let mut anchor_kmer = &fwd_seq[0..1];

        // Iterate over k-mers to find junctions.
        for i in 0..fwd_seq.len()-k+1 {
            let fw_kmer = &fwd_seq[i..i+k];

            if !LdBG::has_junction(&graph, junctions, fw_kmer, false) && !LdBG::has_junction(&graph, junctions, fw_kmer, true) {
                anchor_kmer = &fwd_seq[i..i+k];
            }

            if anchor_kmer.len() == k && LdBG::has_junction(&graph, junctions, fw_kmer, outgoing) {
                let cn_anchor_kmer_vec = LdBG::canonicalize_kmer(anchor_kmer);
                let cn_anchor_kmer = cn_anchor_kmer_vec.as_bytes();

                // Populate link.
                let mut link = Link::new(anchor_kmer == cn_anchor_kmer);

                for j in i..fwd_seq.len()-k {
                    let next_kmer = &fwd_seq[j..j+k];

                    let has_junction = LdBG::has_junction(&graph, junctions, next_kmer, true);
                    if has_junction {
                        link.push_back(fwd_seq[j+k]);
                    }
                }

                // Add link to links map.
                if !links.contains_key(cn_anchor_kmer) {
                    links.insert(cn_anchor_kmer.to_owned(), HashMap::new());
                }

                if !links.get(cn_anchor_kmer).unwrap().contains_key(&link) {
                    links.get_mut(cn_anchor_kmer)
                         .unwrap()
                         .insert(link, 1);
                } else {
                    let linkcov = *links.get_mut(cn_anchor_kmer).unwrap().get(&link).unwrap();

                    links.get_mut(cn_anchor_kmer)
                         .unwrap()
                         .insert(link, linkcov.saturating_add(1));
                }
            }
        }
    }

    /// Build the links for a de Bruijn graph from a vector of sequences.
    fn build_links(k: usize, fwd_seqs: &Vec<Vec<u8>>, graph: &KmerGraph, junctions: &HashMap<Vec<u8>, bool>) -> Links {
        let progress_bar_style = indicatif::ProgressStyle::default_bar()
            .template("Building links... [{elapsed_precise}] [{bar:40.white/white}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("#>-");

        let links: Links = fwd_seqs.par_iter().progress_with_style(progress_bar_style).map(|fwd_seq| {
            let mut local_links = Links::new();
            let fw_seq = fwd_seq.clone();
            let rc_seq = fw_seq.reverse_complement();

            LdBG::add_record_to_links(&mut local_links, &fw_seq, k, graph, junctions, true);
            LdBG::add_record_to_links(&mut local_links, &rc_seq, k, graph, junctions, false);

            local_links
        }).reduce(Links::new, |mut acc, local_links| {
            for (k, v) in local_links {
                acc.entry(k).or_insert_with(HashMap::new).extend(v);
            }
            acc
        });

        links
    }

    /// Get the canonical (lexicographically-lowest) version of a k-mer.
    // #[inline(always)]
    // fn canonicalize_kmer(kmer: &[u8]) -> Vec<u8> {
    //     std::cmp::min(kmer, &kmer.reverse_complement()).to_vec()
    // }

    /// Get the canonical (lexicographically-lowest) version of a k-mer.
    #[inline(always)]
    fn canonicalize_kmer(kmer: &[u8]) -> Vec<u8> {
        let rc_kmer = kmer.reverse_complement();
        if kmer < rc_kmer.as_bytes() {
            kmer.to_vec()
        } else {
            rc_kmer.as_bytes().to_vec()
        }
    }

    /// Find junctions in the graph.
    fn find_junctions(kmers: &KmerGraph) -> HashMap<Vec<u8>, bool> {
        kmers
            .par_iter()
            .flat_map(|(cn_kmer, r)| {
                let mut junction_kmers = HashMap::new();

                let fw_kmer = cn_kmer.clone();
                let rc_kmer = cn_kmer.reverse_complement();

                let (in_degree, out_degree) = (r.in_degree(), r.out_degree());

                if out_degree > 1 {
                    junction_kmers.insert(fw_kmer.clone(), true);
                    junction_kmers.insert(rc_kmer.clone(), false);
                }

                if in_degree > 1 {
                    junction_kmers.insert(fw_kmer.clone(), false);
                    junction_kmers.insert(rc_kmer.clone(), true);
                }

                // for outgoing in [ true, false ] {
                //     if (out_degree > 1 && outgoing) || (in_degree > 1 && !outgoing) {
                //         junction_kmers.insert(fw_kmer.clone(), outgoing);
                //     }

                //     if (in_degree > 1 && outgoing) || (out_degree > 1 && !outgoing) {
                //         junction_kmers.insert(rc_kmer.clone(), outgoing);
                //     }
                // }

                junction_kmers
            })
            .collect()
    }


    /// Check if the given k-mer represents a junction (in the orientation of the given k-mer).
    fn has_junction(graph: &KmerGraph, junctions: &HashMap<Vec<u8>, bool>, kmer: &[u8], outgoing: bool) -> bool {
        // let cn_kmer = LdBG::canonicalize_kmer(kmer);

        // if let Some(r) = graph.get(&cn_kmer) {
        //     let is_canonical = kmer == cn_kmer.as_slice();
        //     let (in_degree, out_degree) = (r.in_degree(), r.out_degree());

        //     return if is_canonical {
        //         (out_degree > 1 && outgoing) || (in_degree > 1 && !outgoing)
        //     } else {
        //         (in_degree > 1 && outgoing) || (out_degree > 1 && !outgoing)
        //     };
        // }

        // false

        junctions.contains_key(kmer) && *junctions.get(kmer).unwrap() == outgoing
    }

    /// Starting at a given k-mer, get the next k-mer (or return None if there isn't a single outgoing edge).
    fn next_kmer(&self, kmer: &[u8], links_in_scope: &mut Vec<Link>) -> Option<Vec<u8>> {
        let cn_kmer_vec = LdBG::canonicalize_kmer(kmer).to_owned();
        let cn_kmer = cn_kmer_vec.as_bytes();

        let r = self.kmers.get(cn_kmer)?;

        let next_base = if cn_kmer == kmer {
            match r.out_degree() {
                0 => return None,
                1 => r.outgoing_edges()[0],
                _ => {
                    let consensus_junction_choice = *links_in_scope.get(0)?.front().unwrap();

                    if r.outgoing_edges().contains(&consensus_junction_choice) {
                        links_in_scope.iter_mut().for_each(|link| { link.pop_front(); });
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
                _ => {
                    let consensus_junction_choice = *links_in_scope.get(0)?.front().unwrap();

                    if r.incoming_edges().contains(&complement(consensus_junction_choice)) {
                        links_in_scope.iter_mut().for_each(|link| { link.pop_front(); });
                        links_in_scope.retain(|link| !link.is_empty());
                        consensus_junction_choice
                    } else {
                        return None;
                    }
                }
            }
        };

        let next_kmer = [&kmer[1..kmer.len()], &[next_base]].concat();

        Some(next_kmer)
    }

    /// Starting at a given k-mer, get the previous k-mer (or return None if there isn't a single incoming edge).
    fn prev_kmer(&self, kmer: &[u8], links_in_scope: &mut Vec<Link>) -> Option<Vec<u8>> {
        let cn_kmer_vec = LdBG::canonicalize_kmer(kmer).to_owned();
        let cn_kmer = cn_kmer_vec.as_bytes();
        // let dbg_cn_kmer_str = String::from_utf8(cn_kmer_vec.clone()).unwrap();

        let r = self.kmers.get(cn_kmer)?;
        // let dbg_r_str = r.to_string();

        let prev_base = if cn_kmer == kmer {
            match r.in_degree() {
                0 => return None,
                1 => r.incoming_edges()[0],
                _ => {
                    // let consensus_junction_choice = *links_in_scope.get(0)?.front().unwrap();
                    let consensus_junction_choice = match links_in_scope.get(0)?.front() {
                        Some(choice) => *choice,
                        None => return None,
                    };

                    if r.incoming_edges().contains(&consensus_junction_choice) {
                        links_in_scope.iter_mut().for_each(|link| { link.pop_front(); });
                        links_in_scope.retain(|link| !link.is_empty());
                        consensus_junction_choice
                    } else {
                        return None;
                    }
                }
            }
        } else {
            match r.out_degree() {
                0 => return None,
                1 => complement(r.outgoing_edges()[0]),
                _ => {
                    let consensus_junction_choice = match links_in_scope.get(0)?.front() {
                        Some(choice) => *choice,
                        None => return None,
                    };

                    if r.outgoing_edges().contains(&complement(consensus_junction_choice)) {
                        links_in_scope.iter_mut().for_each(|link| { link.pop_front(); });
                        links_in_scope.retain(|link| !link.is_empty());
                        consensus_junction_choice
                    } else {
                        return None;
                    }
                }
            }
        };

        let prev_kmer = [&[prev_base], &kmer[0..kmer.len()-1]].concat();

        Some(prev_kmer)
    }

    /// Assemble all contigs from the linked de Bruijn graph.
    pub fn assemble_all(&self) -> Vec<Vec<u8>> {
        let mut contigs = Vec::new();

        let mut used_kmers = HashSet::new();
        let k = self.kmers.keys().next().unwrap().len();

        for cn_kmer in self.kmers.keys() {
            if !used_kmers.contains(cn_kmer) {
                let r = self.kmers.get(cn_kmer).unwrap();

                if r.in_degree() == 1 && r.out_degree() == 1 {
                    let contig = self.assemble(cn_kmer);
                    for kmer_in_contig in contig.windows(k) {
                        used_kmers.insert(Self::canonicalize_kmer(kmer_in_contig));
                    }
                    contigs.push(contig);
                }

                used_kmers.insert(cn_kmer.clone());
            }
        }

        contigs
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
        let mut links_in_scope: Vec<Link> = Vec::new();
        let mut used_links = HashSet::new();
        let mut last_kmer = start_kmer.clone();
        loop {
            self.update_links(&mut links_in_scope, &last_kmer, &mut used_links);

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
        let mut links_in_scope: Vec<Link> = Vec::new();
        let mut used_links = HashSet::new();
        let mut last_kmer = start_kmer.clone();
        loop {
            self.update_links(&mut links_in_scope, &last_kmer, &mut used_links);

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

    /// Update links available to inform navigation during graph traversal.
    fn update_links(&self, links_in_scope: &mut Vec<Link>, last_kmer: &Vec<u8>, used_links: &mut HashSet<Link>) {
        let cn_kmer_vec = LdBG::canonicalize_kmer(last_kmer.as_bytes()).to_owned();
        let cn_kmer = cn_kmer_vec.as_bytes();

        if self.links.contains_key(cn_kmer.as_bytes()) {
            let available_links = self.links.get(cn_kmer.as_bytes()).unwrap().clone();

            let record_orientation_matches_kmer = last_kmer == cn_kmer;

            for jv in available_links {
                let link_goes_forward = record_orientation_matches_kmer == jv.0.is_forward();
                let new_link = if link_goes_forward { jv.0.clone() } else { jv.0.complement() };

                if !used_links.contains(&new_link) {
                    used_links.insert(new_link.clone());
                    links_in_scope.push(new_link);
                }
            }
        }
    }

    /// Pretty-print k-mer, edge, and record information.
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

    use rand::{Rng, SeedableRng};

    /// Canonical example genome from https://academic.oup.com/bioinformatics/article/34/15/2556/4938484
    fn get_test_genome() -> Vec<u8> {
        "ACTGATTTCGATGCGATGCGATGCCACGGTGG".as_bytes().to_vec()
    }

    /// Canonical example read from https://academic.oup.com/bioinformatics/article/34/15/2556/4938484
    fn get_test_read() -> Vec<u8> {
        "TTTCGATGCGATGCGATGCCACG".as_bytes().to_vec()
    }

    fn generate_random_genome(length: usize, seed: u64) -> Vec<u8> {
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(seed);

        let alphabet = b"ACGT";
        let mut genome = Vec::with_capacity(length);

        for _ in 0..length {
            let idx = rng.gen_range(0..4);
            genome.push(alphabet[idx]);
        }

        genome
    }

    #[test]
    fn test_from_sequences() {
        let genome = get_test_genome();
        let fwd_seqs = vec!(genome);

        let g = LdBG::from_sequences(String::from("test"), 5, &fwd_seqs, true);

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
        let g = LdBG::from_sequences(String::from("test"), 5, &fwd_seqs, true);

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
    }

    #[test]
    fn test_assemble_random_genomes() {
        for length in (50..3000).step_by(50) {
            let random_genome = generate_random_genome(length, 0);
            for k in (11..21).step_by(2) {
                let fwd_seqs = vec!(random_genome.clone());

                println!("length={} k={}", length, k);

                let g = LdBG::from_sequences(String::from("test"), k, &fwd_seqs, true);
                let contigs = g.assemble_all();

                // println!("genome -- {:?}", std::str::from_utf8(random_genome.as_bytes()));

                for fw_contig in contigs {
                    assert!(!fw_contig.is_empty(), "Contig should not be empty (length={} kmer={})", length, k);
                    assert!(fw_contig.len() >= k, "Contig should be at least k-mer size (length={} kmer={})", length, k);

                    let rc_contig = fw_contig.reverse_complement();

                    // println!("contig -- {:?}", std::str::from_utf8(fw_contig.as_bytes()));
                    // println!("contig -- {:?}", std::str::from_utf8(rc_contig.as_bytes()));

                    assert!(random_genome == fw_contig || random_genome == rc_contig, "Contig should match the genome or its reverse complement (length={} kmer={})", length, k);
                }
            }
        }
    }

    #[test]
    fn test_assemble_with_and_without_links() {
        let length = 2700;
        let k = 11;

        let random_genome = generate_random_genome(length, 0);
        let fwd_seqs = vec!(random_genome.clone());

        let g1 = LdBG::from_sequences(String::from("without_links"), k, &fwd_seqs, false);
        let g2 = LdBG::from_sequences(String::from("with_links"), k, &fwd_seqs, true);

        let contigs1 = g1.assemble_all();
        let contigs2 = g2.assemble_all();

        assert!(contigs1.len() == 3);
        assert!(contigs2.len() == 1);
    }
}