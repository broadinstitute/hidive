use std::collections::{HashMap, HashSet, VecDeque};
use std::io::{self, Write};

use indicatif::ProgressIterator;
use itertools::Itertools;
use parquet::data_type::AsBytes;

use needletail::sequence::complement;
use needletail::Sequence;
use std::path::PathBuf;

use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::visit::{EdgeRef, NodeIndexable, NodeRef};

use indicatif::ParallelProgressIterator;
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;

use gbdt::config::{loss2string, Config, Loss};
use gbdt::decision_tree::{Data, DataVec};
use gbdt::gradient_boost::GBDT;

use crate::edges::Edges;
use crate::link::Link;
use crate::record::Record;

type KmerGraph = HashMap<Vec<u8>, Record>;
type KmerScores = HashMap<Vec<u8>, f32>;
type Links = HashMap<Vec<u8>, HashMap<Link, u16>>;

/// Represents a linked de Bruijn graph with a k-mer size specified at construction time.
#[derive(Debug)]
pub struct LdBG {
    pub name: String,
    pub kmer_size: usize,
    pub kmers: KmerGraph,
    pub scores: KmerScores,
    pub links: Links,
}

impl LdBG {
    /// Create a de Bruijn graph (and optional links) from a file path.
    ///
    /// # Arguments
    ///
    /// * `name` - A string representing the name of the graph.
    /// * `kmer_size` - The k-mer size.
    /// * `seq_path` - A path to the sequence file.
    /// * `build_links` - A boolean indicating whether to build links.
    ///
    /// # Returns
    ///
    /// A new instance of `LdBG`.
    pub fn from_file(
        name: String,
        kmer_size: usize,
        seq_path: &PathBuf,
        clean: bool,
        build_links: bool,
    ) -> Self {
        let reader = bio::io::fasta::Reader::from_file(seq_path).unwrap();
        let all_reads: Vec<bio::io::fasta::Record> = reader.records().map(|r| r.unwrap()).collect();

        let fwd_seqs: Vec<Vec<u8>> = all_reads
            .iter()
            .map(|r| r.seq().as_bytes().to_vec())
            .collect();

        let mut kmers = Self::build_graph(kmer_size, &fwd_seqs);

        let scores: KmerScores = kmers.keys().map(|k| (k.clone(), 1.0)).collect();

        if clean {
            kmers = Self::clean_graph(&kmers, &scores);
        }

        let links = match build_links {
            true => Self::build_links(kmer_size, &fwd_seqs, &kmers),
            false => Links::new(),
        };

        LdBG {
            name,
            kmer_size,
            kmers,
            scores,
            links,
        }
    }

    /// Create a de Bruijn graph (and optional links) from a list of sequences.
    ///
    /// # Arguments
    ///
    /// * `name` - A string representing the name of the graph.
    /// * `kmer_size` - The k-mer size.
    /// * `fwd_seqs` - A vector of forward sequences.
    /// * `build_links` - A boolean indicating whether to build links.
    ///
    /// # Returns
    ///
    /// A new instance of `LdBG`.
    pub fn from_sequences(
        name: String,
        kmer_size: usize,
        fwd_seqs: &Vec<Vec<u8>>,
        clean: bool,
        build_links: bool,
    ) -> Self {
        let mut kmers = Self::build_graph(kmer_size, fwd_seqs);

        let scores: KmerScores = kmers.keys().map(|k| (k.clone(), 1.0)).collect();

        if clean {
            kmers = Self::clean_graph(&kmers, &scores);
        }

        let links = match build_links {
            true => Self::build_links(kmer_size, fwd_seqs, &kmers),
            false => Links::new(),
        };

        LdBG {
            name,
            kmer_size,
            kmers,
            scores,
            links,
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
    fn add_record_to_graph(
        graph: &mut KmerGraph,
        fw_kmer: &[u8],
        fw_prev_base: u8,
        fw_next_base: u8,
    ) {
        let rc_kmer_vec = fw_kmer.reverse_complement();
        let rc_kmer = rc_kmer_vec.as_bytes();

        let (cn_kmer, can_prev_base, can_next_base) = if fw_kmer < rc_kmer {
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
        graph
            .get_mut(cn_kmer)
            .unwrap()
            .set_incoming_edge(can_prev_base);

        // Set outgoing edge for canonical k-mer.
        graph
            .get_mut(cn_kmer)
            .unwrap()
            .set_outgoing_edge(can_next_base);
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
            for i in 0..fwd_seq.len() - k + 1 {
                let fw_kmer = &fwd_seq[i..i + k];

                let prev = fwd_seq.get(i.wrapping_sub(1));
                let fw_prev_base = *prev.unwrap_or(&b'.');

                let next = fwd_seq.get(i + k);
                let fw_next_base = *next.unwrap_or(&b'.');

                Self::add_record_to_graph(&mut graph, fw_kmer, fw_prev_base, fw_next_base);
            }
        }

        graph
    }

    fn clean_graph(graph: &KmerGraph, scores: &KmerScores) -> KmerGraph {
        let threshold = 0.2;
        let mut cleaned_graph = KmerGraph::new();

        for (cn_kmer, record) in graph {
            if scores.get(cn_kmer).unwrap() > &threshold {
                cleaned_graph.insert(cn_kmer.clone(), record.clone());
            }
        }

        for (cn_kmer, record) in graph {
            if cleaned_graph.contains_key(cn_kmer) {
                let mut new_record = Record::new(record.coverage(), Some(Edges::empty()));

                for e in record.incoming_edges() {
                    let prev_kmer = std::iter::once(e)
                        .chain(cn_kmer[0..cn_kmer.len() - 1].iter().cloned())
                        .collect::<Vec<u8>>();

                    if cleaned_graph.contains_key(&LdBG::canonicalize_kmer(&prev_kmer)) {
                        new_record.set_incoming_edge(e);
                    }
                }

                for e in record.outgoing_edges() {
                    let next_kmer = cn_kmer[1..]
                        .iter()
                        .cloned()
                        .chain(std::iter::once(e))
                        .collect::<Vec<u8>>();

                    if cleaned_graph.contains_key(&LdBG::canonicalize_kmer(&next_kmer)) {
                        new_record.set_outgoing_edge(e);
                    }
                }

                cleaned_graph.insert(cn_kmer.clone(), new_record);
            }
        }

        cleaned_graph
    }

    pub fn mark_tips(&mut self, max_tip_length: usize) {
        for cn_kmer in self.kmers.keys() {
            if self.kmers.get(cn_kmer).unwrap().out_degree() > 1 {
                let next_kmers = self.next_kmers(cn_kmer);

                for cur_kmer in next_kmers {
                    let mut fw_contig = Vec::new();

                    self.assemble_forward(&mut fw_contig, cur_kmer);

                    if fw_contig.len() < max_tip_length {
                        for kmer in fw_contig.kmers(self.kmer_size as u8) {
                            self.scores.insert(LdBG::canonicalize_kmer(&kmer), 0.0);
                        };
                    }
                }
            }

            if self.kmers.get(cn_kmer).unwrap().in_degree() > 1 {
                let prev_kmers = self.prev_kmers(cn_kmer);

                for cur_kmer in prev_kmers {
                    let mut rv_contig = Vec::new();

                    self.assemble_backward(&mut rv_contig, cur_kmer);

                    if rv_contig.len() < max_tip_length {
                        for kmer in rv_contig.kmers(self.kmer_size as u8) {
                            self.scores.insert(LdBG::canonicalize_kmer(&kmer), 0.0);
                        };
                    }
                }
            }
        }
    }

    pub fn remove(&mut self, kmer: &[u8]) -> Option<Record> {
        let cn_kmer = LdBG::canonicalize_kmer(kmer);
        self.kmers.remove(&cn_kmer)
    }

    pub fn infer_edges(&mut self) {
        /*
        let mut kmers = self.kmers.clone();

        kmers.iter_mut().for_each(|(cn_kmer, record)| {
            let mut new_record = Record::new(record.coverage(), Some(Edges::empty()));

            for e in record.incoming_edges() {
                let prev_kmer = std::iter::once(e)
                    .chain(cn_kmer[0..cn_kmer.len() - 1].iter().cloned())
                    .collect::<Vec<u8>>();

                if self
                    .kmers
                    .contains_key(&LdBG::canonicalize_kmer(&prev_kmer))
                {
                    new_record.set_incoming_edge(e);
                }
            }

            for e in record.outgoing_edges() {
                let next_kmer = cn_kmer[1..]
                    .iter()
                    .cloned()
                    .chain(std::iter::once(e))
                    .collect::<Vec<u8>>();

                if self
                    .kmers
                    .contains_key(&LdBG::canonicalize_kmer(&next_kmer))
                {
                    new_record.set_outgoing_edge(e);
                }
            }

            self.kmers.insert(cn_kmer.clone(), new_record);
        });
        */

        let mut kmers = KmerGraph::new();

        self.kmers.iter().for_each(|(cn_kmer, record)| {
            let mut new_record = Record::new(record.coverage(), Some(Edges::empty()));

            for next_kmer in self.next_kmers(cn_kmer) {
                let next_cn_kmer = LdBG::canonicalize_kmer(&next_kmer);

                if self.kmers.contains_key(&next_cn_kmer) {
                    new_record.set_outgoing_edge(next_kmer[next_kmer.len() - 1]);
                }
            }

            for prev_kmer in self.prev_kmers(cn_kmer) {
                let prev_cn_kmer = LdBG::canonicalize_kmer(&prev_kmer);

                if self.kmers.contains_key(&prev_cn_kmer) {
                    new_record.set_incoming_edge(prev_kmer[0]);
                }
            }

            // if record.in_degree() + record.out_degree() != new_record.in_degree() + new_record.out_degree() {
            //     println!("{} {}", String::from_utf8_lossy(cn_kmer), record);
            //     println!("{} {}", String::from_utf8_lossy(cn_kmer), new_record);
            // }

            kmers.insert(cn_kmer.clone(), new_record);
        });

        self.kmers = kmers;
    }

    /// Find an anchor k-mer in a sequence.
    ///
    /// # Arguments
    ///
    /// * `index` - The starting index in the sequence.
    /// * `seq` - A reference to the sequence.
    /// * `k` - The k-mer size.
    /// * `graph` - A reference to the k-mer graph.
    /// * `reverse` - A boolean indicating whether to search in reverse.
    ///
    /// # Returns
    ///
    /// An optional tuple containing the anchor k-mer and its index.
    fn find_anchor_kmer(
        index: usize,
        seq: &[u8],
        k: usize,
        graph: &KmerGraph,
        reverse: bool,
    ) -> Option<(Vec<u8>, usize)> {
        if index > 0 {
            let mut index = index - 1;
            loop {
                if index == 0 {
                    break;
                }

                let anchor_kmer = &seq[index..index + k];

                if !LdBG::has_junction(graph, anchor_kmer, !reverse)
                    || LdBG::has_junction(graph, anchor_kmer, reverse)
                {
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
    /// * `reverse` - A boolean indicating whether to search in reverse.
    /// * `fw` - A boolean indicating whether the sequence is forward.
    fn add_record_to_links(links: &mut Links, seq: &[u8], k: usize, graph: &KmerGraph) {
        let range = (0..seq.len() - k + 1).collect::<Vec<_>>();

        // Iterate over k-mers to find junctions.
        for i in range {
            let fw_kmer = &seq[i..i + k];

            if LdBG::has_junction(graph, fw_kmer, true) {
                if let Some((anchor_kmer_vec, index)) =
                    Self::find_anchor_kmer(i, seq, k, graph, true)
                {
                    let anchor_kmer = anchor_kmer_vec.as_bytes();

                    let cn_anchor_kmer_vec = LdBG::canonicalize_kmer(anchor_kmer);
                    let cn_anchor_kmer = cn_anchor_kmer_vec.as_bytes();

                    // Populate link.
                    let mut link = Link::new(anchor_kmer == cn_anchor_kmer);

                    let sub_range = (index..seq.len() - k).collect::<Vec<_>>();

                    for j in sub_range {
                        let next_kmer = &seq[j..j + k];

                        let has_junction = LdBG::has_junction(graph, next_kmer, false);
                        if has_junction {
                            let choice = seq[j + k];
                            link.push_back(choice);
                        }
                    }

                    if !link.junctions.is_empty() {
                        // Add link to links map.
                        if !links.contains_key(cn_anchor_kmer) {
                            links.insert(cn_anchor_kmer.to_owned(), HashMap::new());
                        }

                        if !links.get(cn_anchor_kmer).unwrap().contains_key(&link) {
                            links.get_mut(cn_anchor_kmer).unwrap().insert(link, 1);
                        } else {
                            let linkcov =
                                *links.get_mut(cn_anchor_kmer).unwrap().get(&link).unwrap();

                            links
                                .get_mut(cn_anchor_kmer)
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
    ///
    /// # Returns
    ///
    /// A map of links.
    pub fn build_links(k: usize, fwd_seqs: &Vec<Vec<u8>>, graph: &KmerGraph) -> Links {
        let progress_bar =
            crate::utils::default_bounded_progress_bar("Building links", fwd_seqs.len() as u64);

        let links: Links = fwd_seqs
            .par_iter()
            .progress_with(progress_bar)
            .map(|fwd_seq| {
                let mut local_links = Links::new();

                let fw_seq = fwd_seq.clone();
                let rc_seq = fw_seq.reverse_complement();

                LdBG::add_record_to_links(&mut local_links, &fw_seq, k, graph);
                LdBG::add_record_to_links(&mut local_links, &rc_seq, k, graph);

                local_links
            })
            .reduce(Links::new, |mut acc, local_links| {
                for (k, v) in local_links {
                    acc.entry(k).or_default().extend(v);
                }
                acc
            });

        links
    }

    pub fn correct_seq(&self, seq: &[u8]) -> Vec<u8> {
        let mut corrected_seq = seq.to_owned();

        let mut uncorrected_regions = Vec::new();
        let mut start = 0;
        let mut in_uncorrected = false;

        for (i, kmer) in seq.windows(self.kmer_size).enumerate() {
            let cn_kmer = Self::canonicalize_kmer(kmer);

            if !self.kmers.contains_key(&cn_kmer) {
                if !in_uncorrected {
                    start = i;
                    in_uncorrected = true;
                }
            } else if in_uncorrected {
                uncorrected_regions.push((start, i + self.kmer_size - 1));
                in_uncorrected = false;
            }
        }

        if in_uncorrected {
            uncorrected_regions.push((start, seq.len()));
        }

        for (start, end) in uncorrected_regions.iter().rev() {
            println!("Read {}, uncorrected region: {} to {}", seq.len(), start, end);

            // Print the uncorrected region of the read, padded by the k-mer size
            let pad_start = start.saturating_sub(self.kmer_size);
            let pad_end = (end + self.kmer_size).min(seq.len());

            let padded_region = &seq[pad_start..pad_end];
            println!("Uncorrected region (padded): {}", String::from_utf8_lossy(padded_region));

            let start_kmer = padded_region.kmers(self.kmer_size as u8).next().unwrap();
            let end_kmer = padded_region.kmers(self.kmer_size as u8).last().unwrap();

            println!("Start k-mer: {}", String::from_utf8_lossy(start_kmer));
            println!("End k-mer: {}", String::from_utf8_lossy(end_kmer));

            let mut fw_contig = start_kmer.to_vec();
            self.assemble_forward(&mut fw_contig, start_kmer.to_vec());

            println!("Fw contig: {}", String::from_utf8_lossy(&fw_contig));

            let mut rv_contig = start_kmer.to_vec();
            self.assemble_backward(&mut rv_contig, end_kmer.to_vec());

            println!("Rc contig: {}", String::from_utf8_lossy(&rv_contig));
        }

        corrected_seq
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

    /// Check if the given k-mer represents a junction (in the orientation of the given k-mer).
    ///
    /// # Arguments
    ///
    /// * `graph` - A reference to the k-mer graph.
    /// * `kmer` - A slice representing the k-mer.
    /// * `reverse` - A boolean indicating whether to check in reverse.
    ///
    /// # Returns
    ///
    /// A boolean indicating whether the k-mer is a junction.
    fn has_junction(graph: &KmerGraph, kmer: &[u8], reverse: bool) -> bool {
        let cn_kmer = LdBG::canonicalize_kmer(kmer);

        if let Some(r) = graph.get(&cn_kmer) {
            let is_canonical = String::from_utf8_lossy(kmer) == String::from_utf8_lossy(&cn_kmer);
            let (in_degree, out_degree) = (r.in_degree(), r.out_degree());

            return if is_canonical {
                if !reverse {
                    out_degree > 1
                } else {
                    in_degree > 1
                }
            } else if !reverse {
                in_degree > 1
            } else {
                out_degree > 1
            };
        }

        false
    }

    fn next_kmers(&self, kmer: &[u8]) -> Vec<Vec<u8>> {
        let cn_kmer_vec = LdBG::canonicalize_kmer(kmer).to_owned();
        let cn_kmer = cn_kmer_vec.as_bytes();

        let ru = self.kmers.get(cn_kmer);
        if ru.is_none() {
            return vec![];
        }

        let r = ru.unwrap();

        let next_kmers = if cn_kmer == kmer {
            r.outgoing_edges().iter().map(|&e| {
                let mut next_kmer = kmer[1..].to_owned();
                next_kmer.push(e);
                next_kmer
            }).collect::<Vec<Vec<u8>>>()
        } else {
            r.incoming_edges().iter().map(|&e| {
                let mut next_kmer = kmer[1..].to_owned();
                next_kmer.push(complement(e));
                next_kmer
            }).collect::<Vec<Vec<u8>>>()
        };

        next_kmers
    }

    fn prev_kmers(&self, kmer: &[u8]) -> Vec<Vec<u8>> {
        let cn_kmer_vec = LdBG::canonicalize_kmer(kmer).to_owned();
        let cn_kmer = cn_kmer_vec.as_bytes();

        let ru = self.kmers.get(cn_kmer);
        if ru.is_none() {
            return vec![];
        }

        let r = ru.unwrap();

        let prev_kmers = if cn_kmer == kmer {
            r.incoming_edges().iter().map(|&e| {
                let mut prev_kmer = kmer[0..kmer.len() - 1].to_owned();
                prev_kmer.insert(0, e);
                prev_kmer
            }).collect::<Vec<Vec<u8>>>()
        } else {
            r.outgoing_edges().iter().map(|&e| {
                let mut prev_kmer = kmer[0..kmer.len() - 1].to_owned();
                prev_kmer.insert(0, complement(e));
                prev_kmer
            }).collect::<Vec<Vec<u8>>>()
        };

        prev_kmers
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

                    let consensus_junction_choice = *links_in_scope.first()?.front().unwrap();

                    if r.outgoing_edges().contains(&consensus_junction_choice) {
                        links_in_scope.iter_mut().for_each(|link| {
                            link.pop_front();
                        });
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

                    let consensus_junction_choice = *links_in_scope.first()?.front().unwrap();

                    if r.incoming_edges()
                        .contains(&complement(consensus_junction_choice))
                    {
                        links_in_scope.iter_mut().for_each(|link| {
                            link.pop_front();
                        });
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

                    let consensus_junction_choice = match links_in_scope.first()?.front() {
                        Some(choice) => *choice,
                        None => return None,
                    };

                    if r.incoming_edges().contains(&consensus_junction_choice) {
                        links_in_scope.iter_mut().for_each(|link| {
                            link.pop_front();
                        });
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

                    let consensus_junction_choice = match links_in_scope.first()?.front() {
                        Some(choice) => *choice,
                        None => return None,
                    };

                    if r.outgoing_edges()
                        .contains(&complement(consensus_junction_choice))
                    {
                        links_in_scope.iter_mut().for_each(|link| {
                            link.pop_front();
                        });
                        links_in_scope.retain(|link| !link.is_empty());
                        consensus_junction_choice
                    } else {
                        return None;
                    }
                }
            }
        };

        let prev_kmer = [&[prev_base], &kmer[0..kmer.len() - 1]].concat();

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
        let k = self.kmer_size;

        let progress_bar = crate::utils::default_bounded_progress_bar(
            "Assembling contigs",
            self.kmers.len() as u64,
        );

        for cn_kmer in self.kmers.keys().progress_with(progress_bar) {
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

        // println!("fw: {}", String::from_utf8_lossy(kmer));
        self.assemble_forward(&mut contig, kmer.to_vec());

        // println!("rc: {}", String::from_utf8_lossy(kmer));
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
                }
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
                }
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
    fn update_links(
        &self,
        links_in_scope: &mut Vec<Link>,
        last_kmer: &Vec<u8>,
        used_links: &mut HashSet<Link>,
        forward: bool,
    ) {
        let cn_kmer_vec = LdBG::canonicalize_kmer(last_kmer.as_bytes()).to_owned();
        let cn_kmer = cn_kmer_vec.as_bytes();

        if self.links.contains_key(cn_kmer.as_bytes()) {
            let available_links = self.links.get(cn_kmer.as_bytes()).unwrap().clone();

            let record_orientation_matches_kmer = last_kmer == cn_kmer;

            for jv in available_links {
                let link_goes_forward = record_orientation_matches_kmer == jv.0.is_forward();
                let new_link = if link_goes_forward {
                    jv.0.clone()
                } else {
                    jv.0.complement()
                };

                if forward == link_goes_forward && !used_links.contains(&new_link) {
                    used_links.insert(new_link.clone());
                    links_in_scope.push(new_link);
                }
            }
        }
    }

    // fn print_kmer(kmer: &[u8], prev_base: u8, next_base: u8, record: Option<&Record>) {
    //     println!("{} {} {} {}",
    //         std::char::from_u32(prev_base as u32).unwrap_or('.'),
    //         std::str::from_utf8(kmer).unwrap(),
    //         std::char::from_u32(next_base as u32).unwrap_or('.'),
    //         record.map_or(String::from(""), |r| format!("{}", r))
    //     );
    // }

    fn first_kmer(&self, contig: &[u8]) -> Option<Vec<u8>> {
        if contig.len() < self.kmer_size as usize {
            return None;
        }

        Some(contig[0..self.kmer_size].to_vec())
    }

    fn last_kmer(&self, contig: &[u8]) -> Option<Vec<u8>> {
        if contig.len() < self.kmer_size as usize {
            return None;
        }

        Some(contig[contig.len() - self.kmer_size as usize..].to_vec())
    }

    pub fn clean_paths(&mut self) -> (usize, usize) {
        let bad_cn_kmers = self.kmers
            .keys()
            .cloned()
            .filter(|cn_kmer| self.scores.get(cn_kmer).unwrap_or(&1.0) < &0.5)
            .filter(|cn_kmer| self.kmers.get(cn_kmer).unwrap().in_degree() == 1 && self.kmers.get(cn_kmer).unwrap().out_degree() == 1)
            .collect::<Vec<Vec<u8>>>();

        // let bad_cn_kmers = vec![b"CTGAGCCGGCCATGTCC".to_vec()];

        let mut to_remove = HashSet::new();
        let mut bad_paths: usize = 0;

        for bad_cn_kmer in bad_cn_kmers {
            if !to_remove.contains(&bad_cn_kmer) {
                let mut seen_kmers = HashSet::new();
                let mut score_sum = 0.0;

                let mut fw_contig = bad_cn_kmer.clone();
                self.assemble_forward(&mut fw_contig, bad_cn_kmer.clone());

                for kmer in fw_contig.windows(self.kmer_size) {
                    let cn_kmer = LdBG::canonicalize_kmer(kmer);

                    if let Some(r) = self.kmers.get(&cn_kmer) {
                        // crate::elog!("fw {} {} {}", String::from_utf8_lossy(&kmer), String::from_utf8_lossy(&cn_kmer), r);

                        if r.in_degree() == 1 && r.out_degree() == 1 {
                            seen_kmers.insert(cn_kmer.clone());
                            score_sum += self.scores.get(&cn_kmer).unwrap_or(&1.0);
                        } else {
                            break;
                        }
                    }
                }

                let mut rc_contig = bad_cn_kmer.clone();
                self.assemble_backward(&mut rc_contig, bad_cn_kmer.clone());

                for kmer in rc_contig.windows(self.kmer_size).rev().skip(1) {
                    let cn_kmer = LdBG::canonicalize_kmer(kmer);

                    if let Some(r) = self.kmers.get(&cn_kmer) {
                        // crate::elog!("rv {} {} {}", String::from_utf8_lossy(&kmer), String::from_utf8_lossy(&cn_kmer), r);

                        if r.in_degree() == 1 && r.out_degree() == 1 {
                            seen_kmers.insert(cn_kmer.clone());
                            score_sum += self.scores.get(&cn_kmer).unwrap_or(&1.0);
                        } else {
                            break;
                        }
                    }
                }

                let weight = score_sum / seen_kmers.len() as f32;

                if weight < 0.5 {
                    // crate::elog!("score {} {} {} {} {}", String::from_utf8_lossy(&bad_cn_kmer), String::from_utf8_lossy(&bad_cn_kmer.reverse_complement()), seen_kmers.len(), score_sum, weight);

                    to_remove.extend(seen_kmers);
                    bad_paths += 1;
                }
            }
        }

        for cn_kmer in &to_remove {
            self.kmers.remove(cn_kmer);
            self.scores.remove(cn_kmer);
        }

        self.infer_edges();

        (to_remove.len(), bad_paths)
    }

    pub fn clean_tips(&mut self, max_tip_length: usize) -> (usize, usize) {
        let mut to_remove = HashSet::new();
        let mut bad_paths: usize = 0;

        // let my_kmer = b"AAATTCGTATCTGTCAA".to_vec();
        // let my_cn_kmer = LdBG::canonicalize_kmer(&my_kmer);

        for cn_kmer in self.kmers.keys() {
            if self.kmers.get(cn_kmer).unwrap().out_degree() > 1 {
                let next_kmers = self.next_kmers(cn_kmer);

                for cur_kmer in next_kmers {
                    let mut fw_contig = cur_kmer.to_vec();
                    self.assemble_forward(&mut fw_contig, cur_kmer.clone());

                    let last_kmer = fw_contig.kmers(self.kmer_size as u8).last().unwrap();
                    let num_next = self.next_kmers(last_kmer).len();

                    if fw_contig.len() <= max_tip_length && num_next == 0 {
                        // crate::elog!("fw {} {} {}", String::from_utf8_lossy(&cur_kmer), String::from_utf8_lossy(&fw_contig), fw_contig.len());

                        for kmer in fw_contig.kmers(self.kmer_size as u8) {
                            to_remove.insert(LdBG::canonicalize_kmer(&kmer));
                        }

                        bad_paths += 1;
                    }
                }
            }

            if self.kmers.get(cn_kmer).unwrap().in_degree() > 1 {
                let prev_kmers = self.prev_kmers(cn_kmer);

                for cur_kmer in prev_kmers {
                    let mut rv_contig = cur_kmer.to_vec();
                    self.assemble_backward(&mut rv_contig, cur_kmer.clone());

                    let first_kmer = rv_contig.kmers(self.kmer_size as u8).next().unwrap();
                    let num_prev = self.prev_kmers(first_kmer).len();

                    if rv_contig.len() <= max_tip_length && num_prev == 0 {
                        // crate::elog!("rv {} {} {}", String::from_utf8_lossy(&cur_kmer), String::from_utf8_lossy(&rv_contig), rv_contig.len());

                        for kmer in rv_contig.kmers(self.kmer_size as u8) {
                            to_remove.insert(LdBG::canonicalize_kmer(&kmer));
                        }

                        bad_paths += 1;
                    }
                }
            }
        }

        for cn_kmer in &to_remove {
            self.kmers.remove(cn_kmer);
            self.scores.remove(cn_kmer);
        }

        self.infer_edges();

        (to_remove.len(), bad_paths)
    }

    pub fn traverse_kmers(&self, start_kmer: Vec<u8>) -> DiGraph<String, f32> {
        let mut graph = DiGraph::new();
        let mut visited = HashMap::new();

        let start_node = graph.add_node(String::from_utf8_lossy(&start_kmer).to_string());

        // Traverse forward
        let mut fwd_stack = vec![start_node];
        while let Some(node) = fwd_stack.pop() {
            let this_kmer = graph.node_weight(node).unwrap().as_bytes();

            for next_kmer in self.next_kmers(this_kmer) {
                let next_node = graph.add_node(String::from_utf8_lossy(&next_kmer).to_string());
                // graph.add_edge(node, next_node, 1.0);
                graph.add_edge(node, next_node, *self.scores.get(&LdBG::canonicalize_kmer(&next_kmer)).unwrap_or(&1.0));

                if !visited.contains_key(&next_kmer) {
                    visited.insert(next_kmer.clone(), next_node);

                    fwd_stack.push(next_node);
                } else {
                    for end_kmer in self.next_kmers(&next_kmer) {
                        let end_node = visited.get(&end_kmer).unwrap();
                        // graph.add_edge(next_node, *end_node, 1.0);
                        graph.add_edge(next_node, *end_node, *self.scores.get(&LdBG::canonicalize_kmer(&end_kmer)).unwrap_or(&1.0));
                    }
                }
            }
        }

        // Traverse reverse
        let mut rev_stack = vec![start_node];
        while let Some(node) = rev_stack.pop() {
            let this_kmer = graph.node_weight(node).unwrap().as_bytes();

            for prev_kmer in self.prev_kmers(this_kmer) {
                let prev_node = graph.add_node(String::from_utf8_lossy(&prev_kmer).to_string());
                // graph.add_edge(prev_node, node, 1.0);
                graph.add_edge(prev_node, node, *self.scores.get(&LdBG::canonicalize_kmer(&prev_kmer)).unwrap_or(&1.0));

                if !visited.contains_key(&prev_kmer) {
                    visited.insert(prev_kmer.clone(), prev_node);

                    rev_stack.push(prev_node);
                } else {
                    for end_kmer in self.prev_kmers(&prev_kmer) {
                        let end_node = visited.get(&end_kmer).unwrap();
                        // graph.add_edge(*end_node, prev_node, 1.0);
                        graph.add_edge(*end_node, prev_node, *self.scores.get(&LdBG::canonicalize_kmer(&end_kmer)).unwrap_or(&1.0));
                    }
                }
            }
        }

        graph
    }

    pub fn traverse_contigs(&self, start_kmer: Vec<u8>) -> DiGraph<String, f32> {
        let mut graph = DiGraph::new();
        let mut visited = HashMap::new();

        let start_contig = self.assemble(&start_kmer);
        let start_node = graph.add_node(String::from_utf8_lossy(&start_contig).to_string());

        // Traverse forward
        let mut fwd_stack = vec![start_node];
        while let Some(node) = fwd_stack.pop() {
            let this_contig = graph.node_weight(node).unwrap().as_bytes();

            if let Some(last_kmer) = self.last_kmer(this_contig) {
                for next_kmer in self.next_kmers(&last_kmer) {
                    if !visited.contains_key(&LdBG::canonicalize_kmer(&next_kmer)) {
                        let mut next_contig = next_kmer.clone();
                        self.assemble_forward(&mut next_contig, next_kmer.clone());

                        let weight = next_contig.as_bytes().kmers(self.kmer_size as u8).into_iter().map(|x| self.scores.get(x).unwrap_or(&1.0)).fold(f32::INFINITY, |acc, x| acc.min(*x));

                        let next_node = graph.add_node(String::from_utf8_lossy(&next_contig).to_string());
                        graph.add_edge(node, next_node, weight);

                        visited.insert(LdBG::canonicalize_kmer(&next_kmer), next_node);

                        fwd_stack.push(next_node);
                    } else {
                        let next_node = visited.get(&LdBG::canonicalize_kmer(&next_kmer)).unwrap();
                        graph.add_edge(node, *next_node, 1.0);
                    }
                }
            }
        }

        // Traverse backward
        let mut rev_stack = vec![start_node];
        while let Some(node) = rev_stack.pop() {
            let this_contig = graph.node_weight(node).unwrap().as_bytes();

            if let Some(first_kmer) = self.first_kmer(this_contig) {
                for prev_kmer in self.prev_kmers(&first_kmer) {
                    if !visited.contains_key(&LdBG::canonicalize_kmer(&prev_kmer)) {
                        let mut prev_contig = prev_kmer.clone();
                        self.assemble_backward(&mut prev_contig, prev_kmer.clone());

                        let weight = prev_contig.as_bytes().kmers(self.kmer_size as u8).into_iter().map(|x| self.scores.get(x).unwrap_or(&1.0)).fold(f32::INFINITY, |acc, x| acc.min(*x));

                        let prev_node = graph.add_node(String::from_utf8_lossy(&prev_contig).to_string());
                        graph.add_edge(node, prev_node, weight);

                        visited.insert(LdBG::canonicalize_kmer(&prev_kmer), prev_node);

                        rev_stack.push(prev_node);
                    } else {
                        let prev_node = visited.get(&LdBG::canonicalize_kmer(&prev_kmer)).unwrap();
                        graph.add_edge(node, *prev_node, 1.0);
                    }
                }
            }
        }

        graph
    }

    pub fn traverse_all_kmers(&self) -> DiGraph<String, f32> {
        let cn_kmers = self.kmers
            .iter()
            .filter(|(_, record)| record.in_degree() <= 1 && record.out_degree() <= 1)
            .map(|(cn_kmer, _)| cn_kmer.clone())
            .collect::<Vec<Vec<u8>>>();

        let mut visited = HashSet::new();
        let mut graph = DiGraph::new();

        for cn_kmer in cn_kmers {
            if !visited.contains(&cn_kmer) {
                let g = self.traverse_kmers(cn_kmer);

                g.node_weights().for_each(|node| {
                    let cn_kmer = LdBG::canonicalize_kmer(node.as_bytes());
                    visited.insert(cn_kmer);
                });

                // Add all nodes from g to graph
                for node_index in g.node_indices() {
                    let node_weight = g.node_weight(node_index).unwrap();
                    graph.add_node(node_weight.clone());
                }

                // Add all edges from g to graph
                for edge in g.edge_references() {
                    let (source, target) = (edge.source(), edge.target());
                    let source_weight = g.node_weight(source).unwrap();
                    let target_weight = g.node_weight(target).unwrap();
                    
                    let new_source = graph.node_indices()
                        .find(|&i| graph[i] == *source_weight)
                        .unwrap();
                    let new_target = graph.node_indices()
                        .find(|&i| graph[i] == *target_weight)
                        .unwrap();
                    
                    graph.add_edge(new_source, new_target, *edge.weight());
                }
            }
        }

        graph
    }

    pub fn traverse_all_contigs(&self) -> DiGraph<String, f32> {
        let cn_kmers = self.kmers
            .iter()
            .filter(|(_, record)| record.in_degree() <= 1 && record.out_degree() <= 1)
            .map(|(cn_kmer, _)| cn_kmer.clone())
            .collect::<Vec<Vec<u8>>>();

        let mut visited = HashSet::new();
        let mut graph = DiGraph::new();

        for cn_kmer in cn_kmers {
            if !visited.contains(&cn_kmer) {
                let g = self.traverse_contigs(cn_kmer);

                g.node_weights().for_each(|node| {
                    for kmer in node.as_bytes().kmers(self.kmer_size as u8) {
                        visited.insert(LdBG::canonicalize_kmer(&kmer));
                    }
                });

                // Add all nodes from g to graph
                for node_index in g.node_indices() {
                    let node_weight = g.node_weight(node_index).unwrap();
                    graph.add_node(node_weight.clone());
                }

                // Add all edges from g to graph
                for edge in g.edge_references() {
                    let (source, target) = (edge.source(), edge.target());
                    let source_weight = g.node_weight(source).unwrap();
                    let target_weight = g.node_weight(target).unwrap();
                    
                    let new_source = graph.node_indices()
                        .find(|&i| graph[i] == *source_weight)
                        .unwrap();
                    let new_target = graph.node_indices()
                        .find(|&i| graph[i] == *target_weight)
                        .unwrap();
                    
                    graph.add_edge(new_source, new_target, *edge.weight());
                }
            }
        }

        graph
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::edges::Edges;

    use rand::{Rng, SeedableRng};

    use std::collections::BTreeMap;
    use std::env;
    use std::fs::File;
    use std::io::{BufRead, BufReader, Write};
    use std::process::{Command, Stdio};

    use flate2::read::GzDecoder;

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
    fn generate_genome_with_tandem_repeats(
        flank_length: usize,
        repeat_length: usize,
        num_repeats: usize,
        seed: u64,
    ) -> Vec<u8> {
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
    fn assemble_with_mccortex(
        k: usize,
        input_genome: &Vec<u8>,
        build_links: bool,
        use_cache: bool,
    ) -> (BTreeMap<String, String>, Links) {
        let current_dir = env::current_dir().unwrap();
        let test_dir = current_dir.join("tests/test_data");
        let test_dir_str = test_dir.to_str().unwrap();

        let copied_genome = input_genome.clone();
        let md5_hash = format!("{:x}", md5::compute(&copied_genome));

        let genome_path = test_dir.join(format!("test.{}.k_{}.l_{}.fa", md5_hash, k, build_links));
        let links_path = test_dir.join(format!(
            "links.{}.k_{}.l_{}.ctp.gz",
            md5_hash, k, build_links
        ));
        let contigs_path =
            test_dir.join(format!("contigs.{}.k_{}.l_{}.fa", md5_hash, k, build_links));

        if !use_cache || !contigs_path.exists() || !links_path.exists() {
            let mut genome_file = File::create(genome_path).expect("Unable to create file");

            writeln!(genome_file, ">genome").expect("Unable to write to file");
            writeln!(genome_file, "{}", String::from_utf8(copied_genome).unwrap())
                .expect("Unable to write to file");

            Command::new("docker")
                .args(&[ "run", "-v", &format!("{}:/data", test_dir_str), "-it", "us.gcr.io/broad-dsp-lrma/lr-mccortex:1.0.0", "mccortex31",
                    "build", "-f", "-m", "1g", "-k", &format!("{}", k), "-s", "test",
                    "-1", &format!("/data/test.{}.k_{}.l_{}.fa", md5_hash, k, build_links),
                    &format!("/data/test.{}.k_{}.l_{}.ctx", md5_hash, k, build_links),
                ])
                .stdout(Stdio::null())
                .stderr(Stdio::null())
                .status()
                .expect("failed to execute process");

            if build_links {
                Command::new("docker")
                    .args(&[ "run", "-v", &format!("{}:/data", test_dir_str), "-it", "us.gcr.io/broad-dsp-lrma/lr-mccortex:1.0.0", "mccortex31",
                        "thread", "-f", "-m", "1g", "-w", "-o", &format!("/data/links.{}.k_{}.l_{}.ctp.gz", md5_hash, k, build_links),
                        "-1", &format!("/data/test.{}.k_{}.l_{}.fa", md5_hash, k, build_links),
                        &format!("/data/test.{}.k_{}.l_{}.ctx", md5_hash, k, build_links),
                    ])
                    .stdout(Stdio::null())
                    .stderr(Stdio::null())
                    .status()
                    .expect("failed to execute process");

                Command::new("docker")
                    .args(&[ "run", "-v", &format!("{}:/data", test_dir_str), "-it", "us.gcr.io/broad-dsp-lrma/lr-mccortex:1.0.0", "mccortex31",
                        "contigs", "-f", "-m", "1g",
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
                    .args(&[ "run", "-v", &format!("{}:/data", test_dir_str), "-it", "us.gcr.io/broad-dsp-lrma/lr-mccortex:1.0.0", "mccortex31",
                        "contigs", "-f", "-m", "1g",
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
        let contig_map: BTreeMap<_, _> = contigs
            .records()
            .into_iter()
            .filter_map(Result::ok)
            .map(|r| {
                let pieces = r
                    .desc()
                    .unwrap()
                    .split_whitespace()
                    .flat_map(|s| s.split("=").collect::<Vec<&str>>())
                    .collect::<Vec<&str>>();

                let seed = pieces.chunks(2).find(|pair| pair[0] == "seed").unwrap()[1].to_string();

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

        let lines: Vec<Vec<String>> = reader
            .lines()
            .map(|line| {
                line.unwrap()
                    .split_whitespace()
                    .map(|s| s.to_string())
                    .collect()
            })
            .collect();

        for i in 0..lines.len() - 1 {
            let pieces1 = &lines[i];

            if pieces1.len() < 1 {
                continue;
            }

            if let Some(first_char_1) = pieces1.get(0).unwrap().chars().next() {
                if "ACGT".contains(first_char_1) {
                    let anchor_kmer = pieces1.get(0).unwrap().as_bytes();
                    let cn_anchor_vec = LdBG::canonicalize_kmer(anchor_kmer);
                    let cn_anchor_kmer = cn_anchor_vec.as_bytes();

                    for j in i + 1..lines.len() {
                        let pieces2 = &lines[j];

                        if pieces2.len() < 1 {
                            continue;
                        }

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
                                    links.get_mut(cn_anchor_kmer).unwrap().insert(link, 1);
                                } else {
                                    let linkcov =
                                        *links.get_mut(cn_anchor_kmer).unwrap().get(&link).unwrap();

                                    links
                                        .get_mut(cn_anchor_kmer)
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
    fn test_traverse_kmers() {
        let genome = get_test_genome();
        let fwd_seqs = vec![genome];

        let g = LdBG::from_sequences(String::from("test"), 5, &fwd_seqs, false, false);
        let graph = g.traverse_kmers(b"ATTTC".to_vec());

        // println!("{}", Dot::with_config(&graph, &[petgraph::dot::Config::EdgeNoLabel]));

        // Check if the graph is consistent with the expected structure
        assert_eq!(graph.node_count(), 22);
        assert_eq!(graph.edge_count(), 22);

        // Check specific nodes and their labels
        let node_labels: Vec<_> = graph.node_weights().cloned().collect();
        assert!(node_labels.contains(&"ATTTC".to_string()));
        assert!(node_labels.contains(&"TTTCG".to_string()));
        assert!(node_labels.contains(&"TTCGA".to_string()));
        assert!(node_labels.contains(&"TCGAT".to_string()));
        assert!(node_labels.contains(&"CGATG".to_string()));
        assert!(node_labels.contains(&"GATGC".to_string()));
        assert!(node_labels.contains(&"ATGCC".to_string()));
        assert!(node_labels.contains(&"ATGCG".to_string()));
        assert!(node_labels.contains(&"TGCGA".to_string()));
        assert!(node_labels.contains(&"GCGAT".to_string()));
        assert!(node_labels.contains(&"TGCCA".to_string()));
        assert!(node_labels.contains(&"GCCAC".to_string()));
        assert!(node_labels.contains(&"CCACG".to_string()));
        assert!(node_labels.contains(&"CACGG".to_string()));
        assert!(node_labels.contains(&"ACGGT".to_string()));
        assert!(node_labels.contains(&"CGGTG".to_string()));
        assert!(node_labels.contains(&"GGTGG".to_string()));
        assert!(node_labels.contains(&"GATTT".to_string()));
        assert!(node_labels.contains(&"TGATT".to_string()));
        assert!(node_labels.contains(&"CTGAT".to_string()));
        assert!(node_labels.contains(&"ACTGA".to_string()));

        // Check specific edges
        let edges = graph.edge_indices().collect::<Vec<_>>();
        assert_eq!(edges.len(), 22);

        // Helper function to find node index by label
        let find_node = |label: &str| graph.node_indices().find(|&i| graph[i] == label).unwrap();

        // Check some specific edges
        assert!(graph.contains_edge(find_node("ATTTC"), find_node("TTTCG")));
        assert!(graph.contains_edge(find_node("TTTCG"), find_node("TTCGA")));
        assert!(graph.contains_edge(find_node("GATGC"), find_node("ATGCC")));
        assert!(graph.contains_edge(find_node("GATGC"), find_node("ATGCG")));
        assert!(graph.contains_edge(find_node("CGATG"), find_node("GATGC")));
        assert!(graph.contains_edge(find_node("GCCAC"), find_node("CCACG")));
        assert!(graph.contains_edge(find_node("ACTGA"), find_node("CTGAT")));
    }

    #[test]
    fn test_traverse_all_kmers() {
        let genome = get_test_genome();
        let fwd_seqs = vec![genome];

        let g = LdBG::from_sequences(String::from("test"), 5, &fwd_seqs, false, false);
        let graph = g.traverse_all_kmers();

        // Write graph as GFA to a string
        let mut gfa_output = Vec::new();
        crate::utils::write_graph_as_gfa(&mut gfa_output, &graph).unwrap();

        // Print GFA string (commented out for test)
        let gfa_string = String::from_utf8(gfa_output).unwrap();
        println!("{}", gfa_string);
    }

    #[test]
    fn test_traverse_all_contigs() {
        let genome = get_test_genome();
        let fwd_seqs = vec![genome];

        let g = LdBG::from_sequences(String::from("test"), 5, &fwd_seqs, false, false);
        let graph = g.traverse_all_contigs();

        // Write graph as GFA to a string
        let mut gfa_output = Vec::new();
        crate::utils::write_graph_as_gfa(&mut gfa_output, &graph).unwrap();

        // Print GFA string (commented out for test)
        let gfa_string = String::from_utf8(gfa_output).unwrap();
        println!("{}", gfa_string);
    }

    #[test]
    fn test_traverse_contigs() {
        let genome = get_test_genome();
        let fwd_seqs = vec![genome];

        let g = LdBG::from_sequences(String::from("test"), 5, &fwd_seqs, false, false);
        let graph = g.traverse_contigs(b"CCACG".to_vec());

        // println!("{}", Dot::with_config(&graph, &[petgraph::dot::Config::EdgeNoLabel]));

        // Check if the graph is consistent with the expected structure
        assert_eq!(graph.node_count(), 3);
        assert_eq!(graph.edge_count(), 4);

        // Check specific nodes and their labels
        let node_labels: Vec<_> = graph.node_weights().cloned().collect();
        assert!(node_labels.contains(&"CGATGCCACGGTGG".to_string()));
        assert!(node_labels.contains(&"ACTGATTTCGAT".to_string()));
        assert!(node_labels.contains(&"CGATGCGAT".to_string()));

        // Check specific edges
        let edges = graph.edge_indices().collect::<Vec<_>>();
        assert_eq!(edges.len(), 4);

        // Helper function to find node index by label
        let find_node = |label: &str| graph.node_indices().find(|&i| graph[i] == label).unwrap();

        // Check some specific edges
        assert!(graph.contains_edge(find_node("CGATGCCACGGTGG"), find_node("ACTGATTTCGAT")));
        assert!(graph.contains_edge(find_node("CGATGCCACGGTGG"), find_node("CGATGCGAT")));
        assert!(graph.contains_edge(find_node("CGATGCGAT"), find_node("ACTGATTTCGAT")));
        assert!(graph.contains_edge(find_node("CGATGCGAT"), find_node("CGATGCGAT")));
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
    fn test_prev_kmers() {
        let genome = get_test_genome();
        let fwd_seqs = vec![genome];

        let g = LdBG::from_sequences(String::from("test"), 5, &fwd_seqs, false, false);

        let prev_kmers_1 = g.prev_kmers(b"TTTCG");
        assert_eq!(prev_kmers_1.len(), 1);
        assert_eq!(prev_kmers_1[0], b"ATTTC");

        let prev_kmers_2: HashSet<Vec<u8>> = g.prev_kmers(b"CGATG").into_iter().collect();
        assert_eq!(prev_kmers_2.len(), 2);
        assert!(prev_kmers_2.contains(&b"TCGAT".to_vec()));
        assert!(prev_kmers_2.contains(&b"GCGAT".to_vec()));

        let prev_kmers_3 = g.prev_kmers(b"ACTGA");
        assert_eq!(prev_kmers_3.len(), 0);
    }

    #[test]
    fn test_next_kmers() {
        let genome = get_test_genome();
        let fwd_seqs = vec![genome];

        let g = LdBG::from_sequences(String::from("test"), 5, &fwd_seqs, false, false);

        let next_kmers_1 = g.next_kmers(b"TTTCG");
        assert_eq!(next_kmers_1.len(), 1);
        assert_eq!(next_kmers_1[0], b"TTCGA");

        let next_kmers_2: HashSet<Vec<u8>> = g.next_kmers(b"GATGC").into_iter().collect();
        assert_eq!(next_kmers_2.len(), 2);
        assert!(next_kmers_2.contains(&b"ATGCG".to_vec()));
        assert!(next_kmers_2.contains(&b"ATGCC".to_vec()));

        let next_kmers_3 = g.next_kmers(b"GGTGG");
        assert_eq!(next_kmers_3.len(), 0);
    }


    #[test]
    fn test_from_sequences() {
        let genome = get_test_genome();
        let fwd_seqs = vec![genome];

        let g = LdBG::from_sequences(String::from("test"), 5, &fwd_seqs, false, true);

        let mut exp_graph = KmerGraph::new();
        exp_graph.insert(
            b"AAATC".to_vec(),
            Record::new(1, Some(Edges::from_string("..g.A...".to_string()))),
        );
        exp_graph.insert(
            b"AATCA".to_vec(),
            Record::new(1, Some(Edges::from_string("a.....G.".to_string()))),
        );
        exp_graph.insert(
            b"ACCGT".to_vec(),
            Record::new(1, Some(Edges::from_string(".c....G.".to_string()))),
        );
        exp_graph.insert(
            b"ACTGA".to_vec(),
            Record::new(1, Some(Edges::from_string(".......T".to_string()))),
        );
        exp_graph.insert(
            b"ATCAG".to_vec(),
            Record::new(1, Some(Edges::from_string("a......T".to_string()))),
        );
        exp_graph.insert(
            b"ATCGA".to_vec(),
            Record::new(1, Some(Edges::from_string(".c..A...".to_string()))),
        );
        exp_graph.insert(
            b"ATCGC".to_vec(),
            Record::new(2, Some(Edges::from_string(".c..A...".to_string()))),
        );
        exp_graph.insert(
            b"ATGCC".to_vec(),
            Record::new(1, Some(Edges::from_string("..g.A...".to_string()))),
        );
        exp_graph.insert(
            b"ATGCG".to_vec(),
            Record::new(2, Some(Edges::from_string("..g.A...".to_string()))),
        );
        exp_graph.insert(
            b"ATTTC".to_vec(),
            Record::new(1, Some(Edges::from_string("..g...G.".to_string()))),
        );
        exp_graph.insert(
            b"CACCG".to_vec(),
            Record::new(1, Some(Edges::from_string(".c.....T".to_string()))),
        );
        exp_graph.insert(
            b"CACGG".to_vec(),
            Record::new(1, Some(Edges::from_string(".c.....T".to_string()))),
        );
        exp_graph.insert(
            b"CATCG".to_vec(),
            Record::new(3, Some(Edges::from_string("..g.AC..".to_string()))),
        );
        exp_graph.insert(
            b"CCACC".to_vec(),
            Record::new(1, Some(Edges::from_string("......G.".to_string()))),
        );
        exp_graph.insert(
            b"CCACG".to_vec(),
            Record::new(1, Some(Edges::from_string("..g...G.".to_string()))),
        );
        exp_graph.insert(
            b"CGAAA".to_vec(),
            Record::new(1, Some(Edges::from_string("...t...T".to_string()))),
        );
        exp_graph.insert(
            b"GATGC".to_vec(),
            Record::new(3, Some(Edges::from_string(".c...CG.".to_string()))),
        );
        exp_graph.insert(
            b"GCCAC".to_vec(),
            Record::new(1, Some(Edges::from_string("...t..G.".to_string()))),
        );
        exp_graph.insert(
            b"TCGAA".to_vec(),
            Record::new(1, Some(Edges::from_string("a...A...".to_string()))),
        );
        exp_graph.insert(
            b"TCGCA".to_vec(),
            Record::new(2, Some(Edges::from_string("a......T".to_string()))),
        );
        exp_graph.insert(
            b"TGCCA".to_vec(),
            Record::new(1, Some(Edges::from_string("a....C..".to_string()))),
        );

        let mut exp_links = Links::new();
        exp_links.insert(
            b"ATGCC".to_vec(),
            HashMap::from([(Link::from_junctions(false, b"CCA"), 1)]),
        );
        exp_links.insert(
            b"ATCGC".to_vec(),
            HashMap::from([
                (Link::from_junctions(false, b"GC"), 1),
                (Link::from_junctions(false, b"C"), 1),
            ]),
        );
        exp_links.insert(
            b"ATGCG".to_vec(),
            HashMap::from([
                (Link::from_junctions(false, b"CA"), 1),
                (Link::from_junctions(false, b"A"), 1),
            ]),
        );
        exp_links.insert(
            b"ATCGA".to_vec(),
            HashMap::from([(Link::from_junctions(false, b"GGC"), 1)]),
        );

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

        let fwd_seqs = vec![fw_genome.clone()];
        let g = LdBG::from_sequences(String::from("test"), 5, &fwd_seqs, false, true);

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

            let g = LdBG::from_sequences(String::from("test"), k, &vec!(random_genome.clone()), false, true);
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
