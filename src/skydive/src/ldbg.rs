use std::collections::{HashMap, HashSet, BTreeMap};
use std::io::{stdout, Write, BufRead, BufReader};
use std::fs::File;
use flate2::read::GzDecoder;

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
    pub kmer_size: usize,
    pub kmers: KmerGraph,
    pub links: Links,
    pub junctions: HashMap<Vec<u8>, bool>,
}

impl LdBG {
    /// Create a de Bruijn graph (and optional links) from a file path.
    ///
    /// # Arguments
    ///
    /// * `name` - A string representing the name of the graph.
    /// * `k` - The k-mer size.
    /// * `seq_path` - A path to the sequence file.
    /// * `build_links` - A boolean indicating whether to build links.
    ///
    /// # Returns
    ///
    /// A new instance of `LdBG`.
    pub fn from_file(name: String, kmer_size: usize, seq_path: &PathBuf, build_links: bool) -> Self {
        let reader = bio::io::fasta::Reader::from_file(seq_path).unwrap();
        let all_reads: Vec<bio::io::fasta::Record> = reader.records().map(|r| r.unwrap()).collect();

        let fwd_seqs: Vec<Vec<u8>> = all_reads
            .iter()
            .map(|r| r.seq().as_bytes().to_vec())
            .collect();

        let kmers = Self::build_graph(kmer_size, &fwd_seqs);

        let junctions = Self::find_junctions(&kmers);

        let links = match build_links {
            true => Self::build_links(kmer_size, &fwd_seqs, &kmers, &junctions),
            false => Links::new()
        };

        LdBG {
            name,
            kmer_size,
            kmers,
            links,
            junctions
        }
    }

    /// Create a de Bruijn graph (and optional links) from a list of sequences.
    ///
    /// # Arguments
    ///
    /// * `name` - A string representing the name of the graph.
    /// * `k` - The k-mer size.
    /// * `fwd_seqs` - A vector of forward sequences.
    /// * `build_links` - A boolean indicating whether to build links.
    ///
    /// # Returns
    ///
    /// A new instance of `LdBG`.
    pub fn from_sequences(name: String, k: usize, fwd_seqs: &Vec<Vec<u8>>, build_links: bool) -> Self {
        let kmers = Self::build_graph(k, fwd_seqs);

        let junctions = Self::find_junctions(&kmers);

        let links = match build_links {
            true => Self::build_links(k, fwd_seqs, &kmers, &junctions),
            false => Links::new()
        };

        LdBG {
            name,
            kmer_size: k,
            kmers,
            links,
            junctions
        }
    }

    /// Get the name of the graph.
    ///
    /// # Returns
    ///
    /// A reference to the name of the graph.
    pub fn name(&self) -> &String {
        &self.name
    }

    /// Add a k-mer, the preceding base, and following base to the graph.
    ///
    /// # Arguments
    ///
    /// * `graph` - A mutable reference to the k-mer graph.
    /// * `fw_kmer` - A slice representing the forward k-mer.
    /// * `fw_prev_base` - The preceding base of the forward k-mer.
    /// * `fw_next_base` - The following base of the forward k-mer.
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
    ///
    /// # Arguments
    ///
    /// * `k` - The k-mer size.
    /// * `fwd_seqs` - A vector of forward sequences.
    ///
    /// # Returns
    ///
    /// A k-mer graph.
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

    /// Find an anchor k-mer in a sequence.
    ///
    /// # Arguments
    ///
    /// * `index` - The starting index in the sequence.
    /// * `seq` - A reference to the sequence.
    /// * `k` - The k-mer size.
    /// * `graph` - A reference to the k-mer graph.
    /// * `junctions` - A reference to the junctions map.
    /// * `reverse` - A boolean indicating whether to search in reverse.
    ///
    /// # Returns
    ///
    /// An optional tuple containing the anchor k-mer and its index.
    fn find_anchor_kmer(index: usize, seq: &Vec<u8>, k: usize, graph: &KmerGraph, junctions: &HashMap<Vec<u8>, bool>, reverse: bool) -> Option<(Vec<u8>, usize)> {
        if index > 0 {
            let mut index = index - 1;
            loop {
                if index == 0 { break; }

                let anchor_kmer = &seq[index..index+k];

                if !LdBG::has_junction(&graph, junctions, anchor_kmer, !reverse) || LdBG::has_junction(&graph, junctions, anchor_kmer, reverse) {
                    return Some((anchor_kmer.to_vec(), index));
                }

                index -= 1;
            }
        }

        None
    }

    /// Add all junction choices from a given sequence.
    ///
    /// # Arguments
    ///
    /// * `links` - A mutable reference to the links map.
    /// * `seq` - A reference to the sequence.
    /// * `k` - The k-mer size.
    /// * `graph` - A reference to the k-mer graph.
    /// * `junctions` - A reference to the junctions map.
    /// * `reverse` - A boolean indicating whether to search in reverse.
    /// * `fw` - A boolean indicating whether the sequence is forward.
    fn add_record_to_links(links: &mut Links, seq: &Vec<u8>, k: usize, graph: &KmerGraph, junctions: &HashMap<Vec<u8>, bool>) {
        let range = (0..seq.len()-k+1).collect::<Vec<_>>();

        // Iterate over k-mers to find junctions.
        for i in range {
            let fw_kmer = &seq[i..i+k];

            if LdBG::has_junction(&graph, junctions, fw_kmer, true) {
                if let Some((anchor_kmer_vec, index)) = Self::find_anchor_kmer(i, seq, k, graph, junctions, true) {
                    let anchor_kmer = anchor_kmer_vec.as_bytes();

                    let cn_anchor_kmer_vec = LdBG::canonicalize_kmer(anchor_kmer);
                    let cn_anchor_kmer = cn_anchor_kmer_vec.as_bytes();

                    // Populate link.
                    let mut link = Link::new(anchor_kmer == cn_anchor_kmer);

                    let sub_range = (index..seq.len()-k).collect::<Vec<_>>();

                    for j in sub_range {
                        let next_kmer = &seq[j..j+k];

                        let has_junction = LdBG::has_junction(&graph, junctions, next_kmer, false);
                        if has_junction {
                            let choice = seq[j+k];
                            link.push_back(choice);
                        }
                    }

                    if link.junctions.len() > 0 {
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
        }
    }

    /// Build the links for a de Bruijn graph from a vector of sequences.
    ///
    /// # Arguments
    ///
    /// * `k` - The k-mer size.
    /// * `fwd_seqs` - A vector of forward sequences.
    /// * `graph` - A reference to the k-mer graph.
    /// * `junctions` - A reference to the junctions map.
    ///
    /// # Returns
    ///
    /// A map of links.
    fn build_links(k: usize, fwd_seqs: &Vec<Vec<u8>>, graph: &KmerGraph, junctions: &HashMap<Vec<u8>, bool>) -> Links {
        let progress_bar_style = indicatif::ProgressStyle::default_bar()
            .template("Building links... [{elapsed_precise}] [{bar:40.white/white}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("#>-");

        let links: Links = fwd_seqs.par_iter().progress_with_style(progress_bar_style).map(|fwd_seq| {
        // let links: Links = fwd_seqs.par_iter().map(|fwd_seq| {
            let mut local_links = Links::new();
            let fw_seq = fwd_seq.clone();
            let rc_seq = fw_seq.reverse_complement();

            LdBG::add_record_to_links(&mut local_links, &fw_seq, k, graph, junctions);
            LdBG::add_record_to_links(&mut local_links, &rc_seq, k, graph, junctions);

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
    ///
    /// # Arguments
    ///
    /// * `kmer` - A slice representing the k-mer.
    ///
    /// # Returns
    ///
    /// A vector containing the canonical k-mer.
    #[inline(always)]
    pub fn canonicalize_kmer(kmer: &[u8]) -> Vec<u8> {
        let rc_kmer = kmer.reverse_complement();
        if kmer < rc_kmer.as_bytes() {
            kmer.to_vec()
        } else {
            rc_kmer.as_bytes().to_vec()
        }
    }

    /// Find junctions in the graph.
    ///
    /// # Arguments
    ///
    /// * `kmers` - A reference to the k-mer graph.
    ///
    /// # Returns
    ///
    /// A map of junctions.
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

                junction_kmers
            })
            .collect()
    }

    /// Check if the given k-mer represents a junction (in the orientation of the given k-mer).
    ///
    /// # Arguments
    ///
    /// * `graph` - A reference to the k-mer graph.
    /// * `junctions` - A reference to the junctions map.
    /// * `kmer` - A slice representing the k-mer.
    /// * `reverse` - A boolean indicating whether to check in reverse.
    ///
    /// # Returns
    ///
    /// A boolean indicating whether the k-mer is a junction.
    fn has_junction(graph: &KmerGraph, junctions: &HashMap<Vec<u8>, bool>, kmer: &[u8], reverse: bool) -> bool {
        let cn_kmer = LdBG::canonicalize_kmer(kmer);

        if let Some(r) = graph.get(&cn_kmer) {
            let is_canonical = String::from_utf8_lossy(kmer).to_string() == String::from_utf8_lossy(&cn_kmer).to_string();
            let (in_degree, out_degree) = (r.in_degree(), r.out_degree());

            if is_canonical {
                if !reverse {
                    return out_degree > 1;
                } else {
                    return in_degree > 1;
                }
            } else {
                if !reverse {
                    return in_degree > 1;
                } else {
                    return out_degree > 1;
                }
            }
        }

        false
    }

    /// Starting at a given k-mer, get the next k-mer (or return None if there isn't a single outgoing edge).
    ///
    /// # Arguments
    ///
    /// * `kmer` - A slice representing the current k-mer.
    /// * `links_in_scope` - A mutable reference to the links in scope.
    ///
    /// # Returns
    ///
    /// An optional vector containing the next k-mer.
    fn next_kmer(&self, kmer: &[u8], links_in_scope: &mut Vec<Link>) -> Option<Vec<u8>> {
        let cn_kmer_vec = LdBG::canonicalize_kmer(kmer).to_owned();
        let cn_kmer = cn_kmer_vec.as_bytes();

        let ru = self.kmers.get(cn_kmer);
        let r = ru?;

        let next_base = if cn_kmer == kmer {
            match r.out_degree() {
                0 => return None,
                1 => r.outgoing_edges()[0],
                _ => {
                    links_in_scope.retain(|link| {
                        if let Some(&next_char) = link.front() {
                            r.outgoing_edges().contains(&next_char)
                        } else {
                            false
                        }
                    });

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
                    links_in_scope.retain(|link| {
                        if let Some(&next_char) = link.front() {
                            r.incoming_edges().contains(&complement(next_char))
                        } else {
                            false
                        }
                    });

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
    ///
    /// # Arguments
    ///
    /// * `kmer` - A slice representing the current k-mer.
    /// * `links_in_scope` - A mutable reference to the links in scope.
    ///
    /// # Returns
    ///
    /// An optional vector containing the previous k-mer.
    fn prev_kmer(&self, kmer: &[u8], links_in_scope: &mut Vec<Link>) -> Option<Vec<u8>> {
        let cn_kmer_vec = LdBG::canonicalize_kmer(kmer).to_owned();
        let cn_kmer = cn_kmer_vec.as_bytes();

        let ru = self.kmers.get(cn_kmer);
        let r = ru?;

        let prev_base = if cn_kmer == kmer {
            match r.in_degree() {
                0 => return None,
                1 => r.incoming_edges()[0],
                _ => {
                    links_in_scope.retain(|link| {
                        if let Some(&next_char) = link.front() {
                            r.incoming_edges().contains(&next_char)
                        } else {
                            false
                        }
                    });

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
                    links_in_scope.retain(|link| {
                        if let Some(&next_char) = link.front() {
                            r.outgoing_edges().contains(&complement(next_char))
                        } else {
                            false
                        }
                    });

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
    ///
    /// # Returns
    ///
    /// A vector of contigs.
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
    ///
    /// # Arguments
    ///
    /// * `kmer` - A slice representing the starting k-mer.
    ///
    /// # Returns
    ///
    /// A vector containing the assembled contig.
    pub fn assemble(&self, kmer: &[u8]) -> Vec<u8> {
        let mut contig: Vec<u8> = kmer.to_vec();

        assert!(
            kmer.len() == self.kmer_size,
            "kmer length {} does not match expected length {}",
            kmer.len(),
            self.kmer_size
        );

        self.assemble_forward(&mut contig, kmer.to_vec());
        self.assemble_backward(&mut contig, kmer.to_vec());

        contig
    }

    /// Assemble a contig in the forward direction.
    ///
    /// # Arguments
    ///
    /// * `contig` - A mutable reference to the contig being assembled.
    /// * `start_kmer` - A vector representing the starting k-mer.
    fn assemble_forward(&self, contig: &mut Vec<u8>, start_kmer: Vec<u8>) {
        let mut links_in_scope: Vec<Link> = Vec::new();
        let mut used_links = HashSet::new();
        let mut last_kmer = start_kmer.clone();
        loop {
            self.update_links(&mut links_in_scope, &last_kmer, &mut used_links, true);

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
    ///
    /// # Arguments
    ///
    /// * `contig` - A mutable reference to the contig being assembled.
    /// * `start_kmer` - A vector representing the starting k-mer.
    fn assemble_backward(&self, contig: &mut Vec<u8>, start_kmer: Vec<u8>) {
        let mut links_in_scope: Vec<Link> = Vec::new();
        let mut used_links = HashSet::new();
        let mut last_kmer = start_kmer.clone();
        loop {
            self.update_links(&mut links_in_scope, &last_kmer, &mut used_links, false);

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
    ///
    /// # Arguments
    ///
    /// * `links_in_scope` - A mutable reference to the links in scope.
    /// * `last_kmer` - A reference to the last k-mer.
    /// * `used_links` - A mutable reference to the set of used links.
    /// * `forward` - A boolean indicating whether to update forward links.
    fn update_links(&self, links_in_scope: &mut Vec<Link>, last_kmer: &Vec<u8>, used_links: &mut HashSet<Link>, forward: bool) {
        let cn_kmer_vec = LdBG::canonicalize_kmer(last_kmer.as_bytes()).to_owned();
        let cn_kmer = cn_kmer_vec.as_bytes();

        if self.links.contains_key(cn_kmer.as_bytes()) {
            let available_links = self.links.get(cn_kmer.as_bytes()).unwrap().clone();

            let record_orientation_matches_kmer = last_kmer == cn_kmer;

            for jv in available_links {
                let link_goes_forward = record_orientation_matches_kmer == jv.0.is_forward();
                let new_link = if link_goes_forward { jv.0.clone() } else { jv.0.complement() };

                if forward == link_goes_forward && !used_links.contains(&new_link) {
                    used_links.insert(new_link.clone());
                    links_in_scope.push(new_link);
                }
            }
        }
    }

    /// Pretty-print k-mer, edge, and record information.
    ///
    /// # Arguments
    ///
    /// * `kmer` - A slice representing the k-mer.
    /// * `prev_base` - The preceding base.
    /// * `next_base` - The following base.
    /// * `record` - An optional reference to the record.
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

    use std::process::{Command, Stdio};
    use std::env;
    use std::fs::File;
    use std::io::Write;

    use proptest::prelude::*;

    /// Canonical example genome from https://academic.oup.com/bioinformatics/article/34/15/2556/4938484
    fn get_test_genome() -> Vec<u8> {
        "ACTGATTTCGATGCGATGCGATGCCACGGTGG".as_bytes().to_vec()
    }

    /// Canonical example read from https://academic.oup.com/bioinformatics/article/34/15/2556/4938484
    fn get_test_read() -> Vec<u8> {
        "TTTCGATGCGATGCGATGCCACG".as_bytes().to_vec()
    }

    /// Generate a genome sequence with tandem repeats.
    ///
    /// # Arguments
    ///
    /// * `flank_length` - The length of the flanking regions.
    /// * `repeat_length` - The length of the repeat region.
    /// * `num_repeats` - The number of repeats.
    /// * `seed` - The seed for the random number generator.
    ///
    /// # Returns
    ///
    /// A vector containing the genome sequence with tandem repeats.
    fn generate_genome_with_tandem_repeats(flank_length: usize, repeat_length: usize, num_repeats: usize, seed: u64) -> Vec<u8> {
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(seed);
        let alphabet = b"ACGT";

        let mut left_flank = Vec::with_capacity(flank_length);
        for _ in 0..flank_length {
            let idx = rng.gen_range(0..4);
            left_flank.push(alphabet[idx]);
        }

        let mut repeat = Vec::with_capacity(repeat_length);
        for _ in 0..repeat_length {
            let idx = rng.gen_range(0..4);
            repeat.push(alphabet[idx]);
        }

        let mut right_flank = Vec::with_capacity(flank_length);
        for _ in 0..flank_length {
            let idx = rng.gen_range(0..4);
            right_flank.push(alphabet[idx]);
        }

        let mut genome = Vec::new();
        genome.extend_from_slice(&left_flank);
        for _ in 0..num_repeats {
            genome.extend_from_slice(&repeat);
        }
        genome.extend_from_slice(&right_flank);

        genome
    }

    /// Assemble using the reference implementation, McCortex.
    fn assemble_with_mccortex(k: usize, input_genome: &Vec<u8>, build_links: bool, use_cache: bool) -> (BTreeMap<String, String>, Links) {
        let current_dir = env::current_dir().unwrap();
        let test_dir = current_dir.join("tests/test_data");
        let test_dir_str = test_dir.to_str().unwrap();

        let copied_genome = input_genome.clone();
        let md5_hash = format!("{:x}", md5::compute(&copied_genome));

        let genome_path = test_dir.join(format!("test.{}.k_{}.l_{}.fa", md5_hash, k, build_links));
        let links_path = test_dir.join(format!("links.{}.k_{}.l_{}.ctp.gz", md5_hash, k, build_links));
        let contigs_path = test_dir.join(format!("contigs.{}.k_{}.l_{}.fa", md5_hash, k, build_links));

        if !use_cache || !contigs_path.exists() || !links_path.exists() {
            let mut genome_file = File::create(genome_path).expect("Unable to create file");

            writeln!(genome_file, ">genome").expect("Unable to write to file");
            writeln!(genome_file, "{}", String::from_utf8(copied_genome).unwrap()).expect("Unable to write to file");
        
            Command::new("docker")
                .args(&[
                    "run", "-v", &format!("{}:/data", test_dir_str), "-it", "us.gcr.io/broad-dsp-lrma/lr-mccortex:1.0.0",
                    "mccortex31", "build", "-f", "-m", "1g", "-k", &format!("{}", k), "-s", "test",
                    "-1", &format!("/data/test.{}.k_{}.l_{}.fa", md5_hash, k, build_links),
                    &format!("/data/test.{}.k_{}.l_{}.ctx", md5_hash, k, build_links),
                ])
                .stdout(Stdio::null())
                .stderr(Stdio::null())
                .status()
                .expect("failed to execute process");

            if build_links {
                Command::new("docker")
                    .args(&[
                        "run", "-v", &format!("{}:/data", test_dir_str), "-it", "us.gcr.io/broad-dsp-lrma/lr-mccortex:1.0.0",
                        "mccortex31", "thread", "-f", "-m 1g", "-w",
                        "-o", &format!("/data/links.{}.k_{}.l_{}.ctp.gz", md5_hash, k, build_links),
                        "-1", &format!("/data/test.{}.k_{}.l_{}.fa", md5_hash, k, build_links),
                        &format!("/data/test.{}.k_{}.l_{}.ctx", md5_hash, k, build_links),
                    ])
                    .stdout(Stdio::null())
                    .stderr(Stdio::null())
                    .status()
                    .expect("failed to execute process");
        
                Command::new("docker")
                    .args(&[
                        "run", "-v", &format!("{}:/data", test_dir_str), "-it", "us.gcr.io/broad-dsp-lrma/lr-mccortex:1.0.0",
                        "mccortex31", "contigs", "-f", "-m", "1g",
                        "-p", &format!("/data/links.{}.k_{}.l_{}.ctp.gz", md5_hash, k, build_links),
                        "-o", &format!("/data/contigs.{}.k_{}.l_{}.fa", md5_hash, k, build_links),
                        &format!("/data/test.{}.k_{}.l_{}.ctx", md5_hash, k, build_links),
                    ])
                    .stdout(Stdio::null())
                    .stderr(Stdio::null())
                    .status()
                    .expect("failed to execute process");
            } else {
                Command::new("docker")
                    .args(&[
                        "run", "-v", &format!("{}:/data", test_dir_str), "-it", "us.gcr.io/broad-dsp-lrma/lr-mccortex:1.0.0",
                        "mccortex31", "contigs", "-f", "-m", "1g",
                        "-o", &format!("/data/contigs.{}.k_{}.l_{}.fa", md5_hash, k, build_links),
                        &format!("/data/test.{}.k_{}.l_{}.ctx", md5_hash, k, build_links),
                    ])
                    .stdout(Stdio::null())
                    .stderr(Stdio::null())
                    .status()
                    .expect("failed to execute process");
            }
        }

        let contig_map = parse_mccortex_fasta(contigs_path);
        let links = parse_mccortex_ctp_file(links_path);

        (contig_map, links)
    }

    /// Parse the McCortex contigs file, grabbing both the assembly seed and assembled sequence.
    fn parse_mccortex_fasta(contigs_path: PathBuf) -> BTreeMap<String, String> {
        let contigs = bio::io::fasta::Reader::from_file(contigs_path).unwrap();
        let contig_map: BTreeMap<_, _> = contigs.records().into_iter()
            .filter_map(Result::ok)
            .map(|r| {
                let pieces = r.desc().unwrap().split_whitespace()
                    .flat_map(|s| s.split("=").collect::<Vec<&str>>())
                    .collect::<Vec<&str>>();
    
                let seed = pieces
                    .chunks(2)
                    .find(|pair| pair[0] == "seed")
                    .unwrap()[1]
                    .to_string();
    
                let seq = String::from_utf8_lossy(r.seq()).to_string();
        
                (seed, seq)
            })
            .collect();

        contig_map
    }
    
    /// Parse the McCortex links file.
    fn parse_mccortex_ctp_file(links_path: PathBuf) -> Links {
        let mut links = Links::new();
    
        let file = File::open(links_path).expect("Unable to open file");
        let decoder = GzDecoder::new(file);
        let reader = BufReader::new(decoder);
    
        let lines: Vec<Vec<String>> = reader.lines().map(|line| {
            line.unwrap().split_whitespace().map(|s| s.to_string()).collect()
        }).collect();
    
        for i in 0..lines.len()-1 {
            let pieces1 = &lines[i];

            if pieces1.len() < 1 { continue; }
    
            if let Some(first_char_1) = pieces1.get(0).unwrap().chars().next() {
                if "ACGT".contains(first_char_1) {
                    let anchor_kmer = pieces1.get(0).unwrap().as_bytes();
                    let cn_anchor_vec = LdBG::canonicalize_kmer(anchor_kmer);
                    let cn_anchor_kmer = cn_anchor_vec.as_bytes();
    
                    for j in i+1..lines.len() {
                        let pieces2 = &lines[j];

                        if pieces2.len() < 1 { continue; }
    
                        if let Some(first_char_2) = pieces2.get(0).unwrap().chars().next() {
                            if "ACGT".contains(first_char_2) {
                                break;
                            } else if "FR".contains(first_char_2) {
                                let mut link = Link::new(pieces2.get(0).unwrap() == "F");
    
                                let junctions = pieces2.get(3).unwrap();
                                for junction in junctions.chars() {
                                    link.push_back(junction as u8);
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
                }
            }
        }

        links
    }

    #[test]
    fn test_canonicalize_kmer() {
        let kmer1 = b"CGTA";
        let kmer2 = b"TACG";
        let kmer3 = b"AAAA";
        let kmer4 = b"TTTT";

        // Test canonical k-mer for kmer1 and kmer2
        assert_eq!(LdBG::canonicalize_kmer(kmer1), b"CGTA".to_vec());
        assert_eq!(LdBG::canonicalize_kmer(kmer2), b"CGTA".to_vec());

        // Test canonical k-mer for kmer3 and kmer4
        assert_eq!(LdBG::canonicalize_kmer(kmer3), b"AAAA".to_vec());
        assert_eq!(LdBG::canonicalize_kmer(kmer4), b"AAAA".to_vec());
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

        let mut exp_links = Links::new();
        exp_links.insert(b"ATGCC".to_vec(), HashMap::from([(Link::from_junctions(false, b"CCA"), 1)]));
        exp_links.insert(b"ATCGC".to_vec(), HashMap::from([
            (Link::from_junctions(false, b"GC"), 1),
            (Link::from_junctions(false, b"C"), 1),
        ]));
        exp_links.insert(b"ATGCG".to_vec(), HashMap::from([
            (Link::from_junctions(false, b"CA"), 1),
            (Link::from_junctions(false, b"A"), 1),
        ]));
        exp_links.insert(b"ATCGA".to_vec(), HashMap::from([(Link::from_junctions(false, b"GGC"), 1)]));

        for (kmer, record) in g.kmers {
            assert!(exp_graph.contains_key(kmer.as_bytes()));

            let exp_record = exp_graph.get(kmer.as_bytes()).unwrap();
            assert!(record.coverage() == exp_record.coverage());
            assert!(record.incoming_edges() == exp_record.incoming_edges());
            assert!(record.outgoing_edges() == exp_record.outgoing_edges());
        }

        for (kmer, link_map) in &g.links {
            assert!(exp_links.contains_key(kmer));
            assert!(link_map.len() == exp_links.get(kmer).unwrap().len());
            for (link, count) in link_map {
                assert!(exp_links.get(kmer).unwrap().contains_key(link));
                assert!(count == exp_links.get(kmer).unwrap().get(link).unwrap());
            }
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

    proptest! {
        #[test]
        fn test_assemble_random_genomes(
            repeat_length in (3..18usize).prop_filter("Increment by 3", |v| v % 3 == 0),
            num_repeats in (2..=3usize),
            k in (11..=17usize).prop_filter("Must be odd", |v| v % 2 == 1)
        ) {
            let random_genome = generate_genome_with_tandem_repeats(500, repeat_length, num_repeats, 0);

            let g = LdBG::from_sequences(String::from("test"), k, &vec!(random_genome.clone()), true);
            let hd_links = g.links.clone();

            let (mc_contigs, mc_links) = assemble_with_mccortex(k, &random_genome, true, true);
            let mc_links: BTreeMap<Vec<u8>, BTreeMap<Link, u16>> = mc_links.into_iter()
                .map(|(k, v)| (k, v.into_iter().collect()))
                .collect();

            let all_kmers: HashSet<_> = mc_links.keys().chain(hd_links.keys()).collect();

            for kmer in all_kmers {
                let mc_links_map = mc_links.get(kmer);
                let hd_links_map = hd_links.get(kmer);

                match (mc_links_map, hd_links_map) {
                    (Some(mc_map), Some(hd_map)) => {
                        assert_eq!(mc_map.len(), hd_map.len(), "Maps have different lengths for kmer: {}", String::from_utf8_lossy(kmer));
                        for mc_link in mc_map.keys() {
                            assert!(hd_map.contains_key(mc_link), "Link {} not in hd map for kmer: {}", mc_link, String::from_utf8_lossy(kmer));
                        }
                    }
                    (Some(mc_map), None) => {
                        for (link, count) in mc_map {
                            println!("mc only: {} {} {}", String::from_utf8_lossy(kmer), link, count);
                        }
                        panic!("There are links unique to reference implementation.");
                    }
                    (None, Some(hd_map)) => {
                        for (link, count) in hd_map {
                            println!("hd only: {} {} {}", String::from_utf8_lossy(kmer), link, count);
                        }
                        panic!("There are links unique to skydive implementation.");
                    }
                    (None, None) => {}
                }
            }

            for (seed, mc_contig) in mc_contigs {
                let hd_contig = String::from_utf8(g.assemble(seed.as_bytes())).expect("Invalid UTF-8 sequence");

                assert_eq!(mc_contig, hd_contig);
            }
        }
    }
}