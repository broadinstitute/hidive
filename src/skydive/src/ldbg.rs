use anyhow::Result;

use std::collections::{HashMap, HashSet};

use parquet::data_type::AsBytes;

use needletail::sequence::complement;
use needletail::Sequence;
use std::path::PathBuf;

use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::visit::EdgeRef;

use indicatif::{ProgressIterator, ParallelProgressIterator};
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;

use gbdt::decision_tree::Data;
use gbdt::gradient_boost::GBDT;

use bio::alignment::distance::*;

use crate::edges::Edges;
use crate::link::Link;
use crate::record::Record;

type KmerGraph = HashMap<Vec<u8>, Record>;
type KmerScores = HashMap<Vec<u8>, f32>;
type Links = HashMap<Vec<u8>, HashMap<Link, u16>>;
type Sources = HashMap<Vec<u8>, Vec<usize>>;

/// Represents a linked de Bruijn graph with a k-mer size specified at construction time.
#[derive(Debug, Clone)]
pub struct LdBG {
    pub name: String,
    pub kmer_size: usize,
    pub kmers: KmerGraph,
    pub scores: KmerScores,
    pub links: Links,
    pub sources: Sources,
    pub verbose: bool,
}

impl LdBG {
    pub fn new(name: String, kmer_size: usize) -> Self {
        LdBG {
            name,
            kmer_size,
            kmers: KmerGraph::new(),
            scores: KmerScores::new(),
            links: Links::new(),
            sources: Sources::new(),
            verbose: false,
        }
    }

    /// Create a de Bruijn graph (and optional links) from a file path.
    ///
    /// # Arguments
    ///
    /// * `name` - A string representing the name of the graph.
    /// * `kmer_size` - The k-mer size.
    /// * `seq_path` - A path to the sequence file.
    ///
    /// # Returns
    ///
    /// A new instance of `LdBG`.
    pub fn from_file(
        name: String,
        kmer_size: usize,
        seq_path: &PathBuf,
    ) -> Self {
        let reader = bio::io::fasta::Reader::from_file(seq_path).unwrap();
        let all_reads: Vec<bio::io::fasta::Record> = reader.records().map(|r| r.unwrap()).collect();

        let fwd_seqs: Vec<Vec<u8>> = all_reads
            .iter()
            .map(|r| r.seq().as_bytes().to_vec())
            .collect();

        let kmers = Self::build_graph(kmer_size, &fwd_seqs);
        let scores: KmerScores = kmers.keys().map(|k| (k.clone(), 1.0)).collect();
        let links: Links = Links::new();
        let sources: Sources = Sources::new();

        LdBG {
            name,
            kmer_size,
            kmers,
            scores,
            links,
            sources,
            verbose: false,
        }
    }

    /// Create a de Bruijn graph (and optional links) from many file paths.
    ///
    /// # Arguments
    ///
    /// * `name` - A string representing the name of the graph.
    /// * `kmer_size` - The k-mer size.
    /// * `seq_paths` - Paths to sequence files.
    ///
    /// # Returns
    ///
    /// A new instance of `LdBG`.
    pub fn from_files(
        name: String,
        kmer_size: usize,
        seq_paths: &Vec<PathBuf>,
    ) -> Self {
        let fwd_seqs = seq_paths
            .iter()
            .map(|p| {
                let reader = bio::io::fasta::Reader::from_file(p).expect("Failed to open file");
                reader.records().filter_map(|r| r.ok()).map(|r| r.seq().to_vec()).collect::<Vec<Vec<u8>>>()
            })
            .flatten()
            .collect();

        let kmers = Self::build_graph(kmer_size, &fwd_seqs);
        let scores: KmerScores = kmers.keys().map(|k| (k.clone(), 1.0)).collect();
        let links: Links = Links::new();
        let sources: Sources = Sources::new();

        LdBG {
            name,
            kmer_size,
            kmers,
            scores,
            links,
            sources,
            verbose: false,
        }
    }

    /// Create a de Bruijn graph (and optional links) from a sequence.
    ///
    /// # Arguments
    ///
    /// * `name` - A string representing the name of the graph.
    /// * `kmer_size` - The k-mer size.
    /// * `fwd_seq` - A forward sequence.
    ///
    /// # Returns
    ///
    /// A new instance of `LdBG`.
    pub fn from_sequence(
        name: String,
        kmer_size: usize,
        fwd_seq: &Vec<u8>,
    ) -> Self {
        let fwd_seqs = vec![fwd_seq.clone()];

        let kmers = Self::build_graph(kmer_size, &fwd_seqs);
        let scores: KmerScores = kmers.keys().map(|k| (k.clone(), 1.0)).collect();
        let links: Links = Links::new();
        let sources: Sources = Sources::new();

        LdBG {
            name,
            kmer_size,
            kmers,
            scores,
            links,
            sources,
            verbose: false,
        }
    }

    /// Create a de Bruijn graph (and optional links) from a list of sequences.
    ///
    /// # Arguments
    ///
    /// * `name` - A string representing the name of the graph.
    /// * `kmer_size` - The k-mer size.
    /// * `fwd_seqs` - A vector of forward sequences.
    ///
    /// # Returns
    ///
    /// A new instance of `LdBG`.
    pub fn from_sequences(
        name: String,
        kmer_size: usize,
        fwd_seqs: &Vec<Vec<u8>>,
    ) -> Self {
        let kmers = Self::build_graph(kmer_size, fwd_seqs);
        let scores: KmerScores = kmers.keys().map(|k| (k.clone(), 1.0)).collect();
        let links: Links = Links::new();
        let sources: Sources = Sources::new();

        LdBG {
            name,
            kmer_size,
            kmers,
            scores,
            links,
            sources,
            verbose: false,
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

    pub fn verbose(mut self, verbose: bool) -> Self {
        self.verbose = verbose;
        self
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
            if fwd_seq.len() < k + 1 {
                continue;
            }

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

    /// Build the links for a de Bruijn graph from a vector of sequences.
    ///
    /// # Arguments
    ///
    /// * `k` - The k-mer size.
    /// * `fwd_seqs` - A vector of forward sequences.
    ///
    /// # Returns
    ///
    /// A map of links.
    pub fn build_links(mut self, fwd_seqs: &Vec<Vec<u8>>, correct: bool) -> Self {
        let progress_bar = if self.verbose {
            crate::utils::default_bounded_progress_bar("Building links", fwd_seqs.len() as u64)
        } else {
            crate::utils::default_hidden_progress_bar()
        };

        let links: Links = fwd_seqs
            .par_iter()
            .progress_with(progress_bar)
            .map(|fwd_seq| if correct { self.correct_seq(fwd_seq) } else { vec![fwd_seq.clone()] })
            .flatten()
            .map(|fwd_seq| {
                let mut local_links = Links::new();

                let fw_seq = fwd_seq.clone();
                let rc_seq = fw_seq.reverse_complement();

                LdBG::add_record_to_links(&mut local_links, &fw_seq, self.kmer_size, &self.kmers);
                LdBG::add_record_to_links(&mut local_links, &rc_seq, self.kmer_size, &self.kmers);

                local_links
            })
            .reduce(Links::new, |mut acc, local_links| {
                for (k, v) in local_links {
                    acc.entry(k).or_default().extend(v);
                }
                acc
            });

        self.links = links;

        self
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
        if seq.len() < k + 1 { return; }

        let range = (0..seq.len() - k + 1).collect::<Vec<_>>();

        // Iterate over k-mers to find junctions.
        for i in range {
            let fw_kmer = &seq[i..i + k];

            if LdBG::has_junction(graph, fw_kmer, true) {
                if let Some((anchor_kmer_vec, index)) =
                    Self::find_anchor_kmer(i, seq, k, graph, true)
                {
                    let anchor_kmer = anchor_kmer_vec.as_bytes();

                    let cn_anchor_kmer_vec = crate::utils::canonicalize_kmer(anchor_kmer);
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
        let cn_kmer = crate::utils::canonicalize_kmer(kmer);

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
        let cn_kmer_vec = crate::utils::canonicalize_kmer(kmer).to_owned();
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
        let cn_kmer_vec = crate::utils::canonicalize_kmer(kmer).to_owned();
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
        let cn_kmer_vec = crate::utils::canonicalize_kmer(kmer).to_owned();
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
        let cn_kmer_vec = crate::utils::canonicalize_kmer(kmer).to_owned();
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

    pub fn remove(&mut self, kmer: &[u8]) -> Option<Record> {
        let cn_kmer = crate::utils::canonicalize_kmer(kmer);
        self.kmers.remove(&cn_kmer)
    }

    pub fn score_kmers(mut self, model_path: &PathBuf) -> Self {
        let gbdt = GBDT::load_model(model_path.to_str().unwrap()).unwrap();

        self.scores = self.kmers.keys().map(|cn_kmer| {
            let lcov = self.kmers.get(cn_kmer).unwrap().coverage();
            let scov: u16 = 0;
            let compressed_len = crate::utils::homopolymer_compressed(cn_kmer).len();

            let data = Data::new_test_data(vec![lcov as f32, scov as f32, (cn_kmer.len() - compressed_len) as f32], Some(0.0));
            let prediction = *gbdt.predict(&vec![data]).first().unwrap_or(&0.0);

            (cn_kmer.clone(), prediction)
        }).collect();

        self
    }

    pub fn infer_edges(&mut self) {
        let mut kmers = KmerGraph::new();

        self.kmers.iter().for_each(|(cn_kmer, record)| {
            let mut new_record = Record::new(record.coverage(), Some(Edges::empty()));

            for edge_base in vec![b'A', b'C', b'G', b'T'] {
                let next_kmer = [&cn_kmer[1..], &[edge_base]].concat();
                let next_cn_kmer = crate::utils::canonicalize_kmer(&next_kmer);

                if self.kmers.contains_key(&next_cn_kmer) {
                    new_record.set_outgoing_edge(edge_base);
                }

                let prev_kmer = [&[edge_base], &cn_kmer[0..cn_kmer.len() - 1]].concat();
                let prev_cn_kmer = crate::utils::canonicalize_kmer(&prev_kmer);

                if self.kmers.contains_key(&prev_cn_kmer) {
                    new_record.set_incoming_edge(edge_base);
                }
            }

            kmers.insert(cn_kmer.clone(), new_record);
        });

        self.kmers = kmers;
    }

    pub fn correct_seq(&self, seq: &[u8]) -> Vec<Vec<u8>> {
        let (mut graph, node_indices) = self.initialize_read_graph(seq);

        for window in node_indices.windows(2) {
            let (from, to) = (window[0], window[1]);

            let from_kmer = graph.node_weight(from).unwrap().as_bytes();
            let to_kmer = graph.node_weight(to).unwrap().as_bytes();

            let from_cn_kmer = crate::utils::canonicalize_kmer(from_kmer);
            let to_cn_kmer = crate::utils::canonicalize_kmer(to_kmer);

            let from_next_kmers = self.next_kmers(&from_kmer);
            let to_prev_kmers = self.prev_kmers(&to_kmer);

            if self.kmers.contains_key(&from_cn_kmer) && self.kmers.contains_key(&to_cn_kmer) && from_next_kmers.contains(&to_kmer.to_vec()) && to_prev_kmers.contains(&from_kmer.to_vec()) {
                graph.add_edge(from, to, 1.0);
            } else {
                if let Ok((contig, forward)) = self.try_simple_paths(from_kmer, to_kmer) {
                    self.add_to_read_graph(contig, &mut graph, forward);
                } else if let Ok((contig, forward)) = self.try_alternate_paths(&graph, from, to, 5) {
                    self.add_to_read_graph(contig, &mut graph, forward);
                } else if let Ok((contig, forward)) = self.try_overlapping_paths(&graph, from, to, 5) {
                    self.add_to_read_graph(contig, &mut graph, forward);
                } else if let Ok((contig, forward)) = self.try_exploratory_paths(&graph, from, to, seq, 5, 100) {
                    self.add_to_read_graph(contig, &mut graph, forward);
                }
            }
        }

        let initial_segments = traverse_read_graph(&graph, 0);

        if let Ok(connecting_segments) = self.try_connecting_all_segments(initial_segments) {
            for (contig, forward) in connecting_segments.into_iter() {
                self.add_to_read_graph(contig, &mut graph, forward);
            }
        }

        // if graph.node_count() > 0 {
        //     self.remove_short_branches(&mut graph, 2 * self.kmer_size);
        
        //     let dot = Dot::with_config(&graph, &[Config::EdgeNoLabel]);
        //     let mut file = File::create("read_graph.dot").expect("Unable to create file");
        //     write!(file, "{:?}", dot).expect("Unable to write data");
        // }

        let corrected_segments = traverse_read_graph(&graph, self.kmer_size);

        if corrected_segments.len() > 1 {
            vec![self.combine_read_segments(seq, &corrected_segments)]
        } else {
            corrected_segments
        }
    }

    fn combine_read_segments(&self, seq: &[u8], corrected_segments: &Vec<Vec<u8>>) -> Vec<u8> {
        use spoa::AlignmentType;
        let mut la = spoa::AlignmentEngine::new(AlignmentType::kSW, 5, -10, -16, -12, -20, -8);

        let mut sg = spoa::Graph::new();
        for lr_seq in vec![seq.to_vec()].iter().chain(corrected_segments.iter()) {
            let seq_cstr = std::ffi::CString::new(lr_seq.clone()).unwrap();
            let seq_qual = std::ffi::CString::new(vec![b'I'; lr_seq.len()]).unwrap();
            let a = la.align(seq_cstr.as_ref(), &sg);

            sg.add_alignment(&a, seq_cstr.as_ref(), seq_qual.as_ref());
        }

        let msa_cstrs = sg.multiple_sequence_alignment(false);
        let msa_strings = msa_cstrs
            .iter()
            .map(|cstr| cstr.to_str().unwrap().to_string())
            .map(|msa_string| {
                let leading_dashes = msa_string.chars().take_while(|&c| c == '-').count();
                let trailing_dashes = msa_string.chars().rev().take_while(|&c| c == '-').count();
                let middle = &msa_string[leading_dashes..msa_string.len() - trailing_dashes];
                format!("{}{}{}", 
                    " ".repeat(leading_dashes), 
                    middle, 
                    " ".repeat(trailing_dashes)
                )
            })
            .collect::<Vec<String>>();

        let mut combined_seq = vec![];
        for column in 0..msa_strings[0].len() {
            let mut pileup = vec![];
            for msa_string in msa_strings.iter().skip(1).chain(msa_strings.iter().take(1)) {
                let char_at_column = msa_string.chars().nth(column).unwrap();

                if char_at_column != ' ' {
                    pileup.push(char_at_column);
                }
            }

            combined_seq.push(pileup[0] as u8);
        }

        combined_seq.retain(|&x| x != b'-');

        combined_seq
    }

    fn remove_short_branches(&self, graph: &mut petgraph::Graph<String, f64>, min_size: usize) {
        // Find nodes with more than one outgoing edge and prune short branches
        let mut nodes_to_remove = HashSet::new();
        for node in graph.node_indices() {
            let out_degree = graph.edges_directed(node, petgraph::Outgoing).count();
            if out_degree > 1 {
                for edge in graph.edges_directed(node, petgraph::Outgoing) {
                    let mut branch = vec![edge.target()];
                    let mut current = edge.target();
                    let mut branch_length = 0;
    
                    while branch_length < 2 * self.kmer_size {
                        let out_edges: Vec<_> = graph.edges_directed(current, petgraph::Outgoing).collect();
                        if out_edges.len() != 1 {
                            break;
                        }
                        current = out_edges[0].target();
                        branch.push(current);
                        branch_length += 1;
                    }
    
                    if branch_length < min_size {
                        nodes_to_remove.extend(branch);
                    }
                }
            }
        }
    
        // Remove the identified nodes from the graph
        for node in nodes_to_remove {
            graph.remove_node(node);
        }
    }
    
    fn try_simple_paths(&self, from_kmer: &[u8], to_kmer: &[u8]) -> Result<(Vec<u8>, bool)> {
        let mut forward_contig = from_kmer.to_vec();
        let mut backward_contig = to_kmer.to_vec();

        if let Ok(_) = self.assemble_forward_until(&mut forward_contig, from_kmer.to_vec(), HashSet::from([to_kmer.to_vec()]), 100) {
            return Ok((forward_contig, true));
        } else if let Ok(_) = self.assemble_backward_until(&mut backward_contig, to_kmer.to_vec(), HashSet::from([from_kmer.to_vec()]), 100) {
            return Ok((backward_contig, false));
        }

        Err(anyhow::anyhow!("No simple path found"))
    }

    fn try_alternate_paths(&self, graph: &petgraph::Graph<String, f64>, from: NodeIndex, to: NodeIndex, window: u8) -> Result<(Vec<u8>, bool)> {
        for scan_left in 0..=window {
            for scan_right in 0..=window {
                if let (Some(new_from), Some(new_to)) = (navigate_backward(from, scan_left, graph), navigate_forward(to, scan_right, graph)) {
                    let new_from_kmer = graph.node_weight(new_from).unwrap().as_bytes();
                    let new_to_kmer = graph.node_weight(new_to).unwrap().as_bytes();

                    let mut forward_contig = new_from_kmer.to_vec();
                    let mut backward_contig = new_to_kmer.to_vec();

                    if let Ok(_) = self.assemble_forward_until(&mut forward_contig, new_from_kmer.to_vec(), HashSet::from([new_to_kmer.to_vec()]), 100) {
                        return Ok((forward_contig, true));
                    } else if let Ok(_) = self.assemble_backward_until(&mut backward_contig, new_to_kmer.to_vec(), HashSet::from([new_from_kmer.to_vec()]), 100) {
                        return Ok((backward_contig, false));
                    }
                }
            }
        }

        return Err(anyhow::anyhow!("No alternate path found"));
    }

    fn try_overlapping_paths(&self, graph: &petgraph::Graph<String, f64>, from: NodeIndex, to: NodeIndex, window: u8) -> Result<(Vec<u8>, bool)> {
        for scan_left in 0..=window {
            for scan_right in 0..=window {
                if let (Some(new_from), Some(new_to)) = (navigate_backward(from, scan_left, graph), navigate_forward(to, scan_right, graph)) {
                    let new_from_kmer = graph.node_weight(new_from).unwrap().as_bytes();
                    let new_to_kmer = graph.node_weight(new_to).unwrap().as_bytes();

                    let mut forward_contig = new_from_kmer.to_vec();
                    let mut backward_contig = new_to_kmer.to_vec();

                    let seqs = vec![forward_contig.clone(), backward_contig.clone()];
                    let l = LdBG::from_sequences("contig".to_string(), self.kmer_size, &seqs).build_links(&seqs, false);

                    if let Ok(_) = l.assemble_forward_until(&mut forward_contig, new_from_kmer.to_vec(), HashSet::from([new_to_kmer.to_vec()]), 100) {
                        return Ok((forward_contig, true));
                    } else if let Ok(_) = l.assemble_backward_until(&mut backward_contig, new_to_kmer.to_vec(), HashSet::from([new_from_kmer.to_vec()]), 100) {
                        return Ok((backward_contig, false));
                    }
                }
            }
        }

        return Err(anyhow::anyhow!("No overlapping path found"));
    }

    fn try_exploratory_paths(&self, graph: &petgraph::Graph<String, f64>, from: NodeIndex, to: NodeIndex, seq: &[u8], window: u8, max_intermediate_nodes: usize) -> Result<(Vec<u8>, bool)> {
        for scan_left in 0..=window {
            for scan_right in 0..=window {
                if let (Some(new_from), Some(new_to)) = (navigate_backward(from, scan_left, graph), navigate_forward(to, scan_right, graph)) {
                    let new_from_kmer = graph.node_weight(new_from).unwrap().as_bytes();
                    let new_to_kmer = graph.node_weight(new_to).unwrap().as_bytes();

                    if let Some(substring) = self.find_read_substring(seq, new_from_kmer, new_to_kmer) {
                        let mut subgraph = DiGraph::new();
                        let mut visited = HashMap::<String, NodeIndex>::new();

                        let start_label = String::from_utf8_lossy(&new_from_kmer).to_string();
                        let start_node = subgraph.add_node(start_label.clone());
                        visited.insert(start_label.clone(), start_node);

                        if let Ok(_) = self.traverse_forward_until(&mut subgraph, &mut visited, start_node, new_to_kmer) {
                            let start_label = String::from_utf8_lossy(&new_from_kmer).to_string();
                            let end_label = String::from_utf8_lossy(&new_to_kmer).to_string();

                            // Find the node indices for the start and end kmers
                            let start_node = subgraph.node_indices().find(|&n| subgraph[n] == start_label);
                            let end_node = subgraph.node_indices().find(|&n| subgraph[n] == end_label);

                            if let (Some(start), Some(end)) = (start_node, end_node) {
                                let mut best_contig = String::new();
                                let mut min_ldist = u32::MAX;

                                // List all paths from start_node to end_node
                                for path in petgraph::algo::all_simple_paths::<Vec<_>, _>(&subgraph, start, end, 0, Some(max_intermediate_nodes)) {
                                    // Convert the path to labels
                                    let path_labels: Vec<&String> = path.iter().map(|&n| &subgraph[n]).collect();

                                    // Create the full contig from the path_labels
                                    let mut new_contig = String::new();
                                    for kmer in path_labels.iter() {
                                        if new_contig.is_empty() {
                                            new_contig.push_str(kmer);
                                        } else {
                                            new_contig.push_str(&kmer.chars().last().unwrap().to_string());
                                        }
                                    }

                                    let ldist = levenshtein(&new_contig.as_bytes().to_vec(), &substring);
                                    if ldist < min_ldist {
                                        min_ldist = ldist;
                                        best_contig = new_contig;
                                    }
                                }

                                if min_ldist < (substring.len() / 10) as u32 {
                                    return Ok((best_contig.as_bytes().to_vec(), true));
                                }
                            }
                        }
                    }
                }
            }
        }

        return Err(anyhow::anyhow!("No alternate path found"));
    }

    fn try_connecting_all_segments(&self, segments: Vec<Vec<u8>>) -> Result<Vec<(Vec<u8>, bool)>> {
        let mut connecting_segments = Vec::new();

        for segment1 in segments.iter() {
            for segment2 in segments.iter() {
                if segment1 == segment2 {
                    continue;
                }

                if let Ok(result) = self.try_connecting_segments(segment1, segment2) {
                    connecting_segments.push(result);
                }
            }
        }

        if connecting_segments.len() == 0 {
            return Err(anyhow::anyhow!("No connecting path found"));
        }

        Ok(connecting_segments)
    }

    fn try_connecting_segments(&self, segment1: &[u8], segment2: &[u8]) -> Result<(Vec<u8>, bool)> {
        // let start_kmers = segment1.windows(self.kmer_size).skip(segment1.len() - 2*self.kmer_size).map(|kmer| kmer.to_vec()).collect::<HashSet<Vec<u8>>>();
        let start_kmers = segment1.windows(self.kmer_size).skip(segment1.len().saturating_sub(2*self.kmer_size)).map(|kmer| kmer.to_vec()).collect::<HashSet<Vec<u8>>>();
        let end_kmers = segment2.windows(self.kmer_size).take(10).map(|kmer| kmer.to_vec()).collect::<HashSet<Vec<u8>>>();

        for start_kmer in start_kmers.iter() {
            let mut forward_contig = start_kmer.to_vec();
            if let Ok(_) = self.assemble_forward_until(&mut forward_contig, start_kmer.to_vec(), end_kmers.clone(), 100) {
                return Ok((forward_contig, true));
            }
        }

        for end_kmer in end_kmers.iter() {
            let mut backward_contig = end_kmer.to_vec();
            if let Ok(_) = self.assemble_backward_until(&mut backward_contig, end_kmer.to_vec(), start_kmers.clone(), 100) {
                return Ok((backward_contig, false));
            }
        }

        return Err(anyhow::anyhow!("No connecting path found"));
    }
    
    fn find_read_substring(&self, seq: &[u8], from_kmer: &[u8], to_kmer: &[u8]) -> Option<Vec<u8>> {
        seq.windows(from_kmer.len())
            .position(|window| window == from_kmer)
            .and_then(|start_pos| {
                seq[start_pos + from_kmer.len()..]
                    .windows(to_kmer.len())
                    .position(|window| window == to_kmer)
                    .map(|end_pos| {
                        seq[start_pos..start_pos + from_kmer.len() + end_pos + to_kmer.len()].to_vec()
                    })
            })
    }

    fn add_to_read_graph(&self, contig: Vec<u8>, graph: &mut petgraph::Graph<String, f64>, forward: bool) {
        if forward {
            self.add_to_read_graph_forward(contig, graph);
        } else {
            self.add_to_read_graph_backward(contig, graph);
        }
    }

    fn add_to_read_graph_backward(&self, backward_contig: Vec<u8>, graph: &mut petgraph::Graph<String, f64>) {
        // If backward assembly succeeded, add path to graph
        for i in (0..backward_contig.len() - self.kmer_size).rev() {
            let current_kmer = &backward_contig[i..i + self.kmer_size];
            let prev_kmer = &backward_contig[i + 1..i + 1 + self.kmer_size];
    
            let prev_node = graph.node_indices().find(|&n| *graph[n] == String::from_utf8_lossy(current_kmer).to_string())
                .unwrap_or_else(|| graph.add_node(String::from_utf8_lossy(current_kmer).to_string()));
            let current_node = graph.node_indices().find(|&n| *graph[n] == String::from_utf8_lossy(prev_kmer).to_string())
                .unwrap_or_else(|| graph.add_node(String::from_utf8_lossy(prev_kmer).to_string()));
    
            if !graph.contains_edge(prev_node, current_node) {
                graph.add_edge(prev_node, current_node, 1.0);
            }
        }
    }
    
    fn add_to_read_graph_forward(&self, forward_contig: Vec<u8>, graph: &mut petgraph::Graph<String, f64>) {
        // If forward assembly succeeded, add path to graph
        for i in 0..forward_contig.len() - self.kmer_size {
            let current_kmer = &forward_contig[i..i + self.kmer_size];
            let next_kmer = &forward_contig[i + 1..i + 1 + self.kmer_size];
    
            let current_node = graph.node_indices().find(|&n| *graph[n] == String::from_utf8_lossy(current_kmer).to_string())
                .unwrap_or_else(|| graph.add_node(String::from_utf8_lossy(current_kmer).to_string()));
            let next_node = graph.node_indices().find(|&n| *graph[n] == String::from_utf8_lossy(next_kmer).to_string())
                .unwrap_or_else(|| graph.add_node(String::from_utf8_lossy(next_kmer).to_string()));
    
            if !graph.contains_edge(current_node, next_node) {
                graph.add_edge(current_node, next_node, 1.0);
            }
        }
    }
    
    fn initialize_read_graph(&self, seq: &[u8]) -> (petgraph::Graph<String, f64>, Vec<NodeIndex>) {
        let mut graph = DiGraph::<String, f64>::new();
    
        let mut node_indices = Vec::new();
        for kmer in seq.windows(self.kmer_size) {
            if self.kmers.contains_key(&crate::utils::canonicalize_kmer(kmer)) {
                let node_index = graph.add_node(String::from_utf8_lossy(kmer).to_string());
                node_indices.push(node_index);
            }
        }
        (graph, node_indices)
    }
    
    pub fn correct_seqs(&self, seqs: &Vec<Vec<u8>>) -> Vec<Vec<u8>> {
        let progress_bar = if self.verbose {
            crate::utils::default_bounded_progress_bar("Correcting reads", seqs.len() as u64)
        } else {
            crate::utils::default_hidden_progress_bar()
        };

        seqs
            .par_iter()
            .progress_with(progress_bar)
            .map(|seq| self.correct_seq(seq))
            .flatten()
            .collect::<Vec<Vec<u8>>>()
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
        let mut visited = HashSet::new();

        loop {
            self.update_links(&mut links_in_scope, &last_kmer, &mut used_links, true);

            match self.next_kmer(&last_kmer, &mut links_in_scope) {
                Some(this_kmer) => {
                    if links_in_scope.len() == 0 {
                        if visited.contains(&this_kmer) {
                            break;
                        }

                        visited.insert(this_kmer.clone());
                    }

                    contig.push(this_kmer[this_kmer.len() - 1]);
                    last_kmer = this_kmer;
                }
                None => {
                    break;
                }
            }
        }
    }

    /// Assemble a contig in the forward direction until a stopping k-mer is reached or the limit is hit.
    ///
    /// # Arguments
    ///
    /// * `contig` - A mutable reference to the contig being assembled.
    /// * `start_kmer` - A vector representing the starting k-mer.
    /// * `stop_kmers` - A HashSet of k-mers to stop at.
    /// * `limit` - The maximum number of nodes to traverse before giving up.
    fn assemble_forward_until(&self, contig: &mut Vec<u8>, start_kmer: Vec<u8>, stop_kmers: HashSet<Vec<u8>>, limit: usize) -> Result<Vec<u8>> {
        let mut links_in_scope: Vec<Link> = Vec::new();
        let mut used_links = HashSet::new();
        let mut last_kmer = start_kmer.clone();
        let mut nodes_traversed = 0;

        loop {
            if stop_kmers.contains(&last_kmer) {
                return Ok(last_kmer);
            }

            if nodes_traversed >= limit {
                return Err(anyhow::anyhow!("Reached limit without finding a stopping k-mer"));
            }

            self.update_links(&mut links_in_scope, &last_kmer, &mut used_links, true);

            match self.next_kmer(&last_kmer, &mut links_in_scope) {
                Some(this_kmer) => {
                    contig.push(this_kmer[this_kmer.len() - 1]);
                    last_kmer = this_kmer;
                    nodes_traversed += 1;
                }
                None => {
                    break;
                }
            }
        }

        Err(anyhow::anyhow!("No stopping k-mer found"))
    }

    fn assemble_forward_until_condition<F>(&self, contig: &mut Vec<u8>, start_kmer: Vec<u8>, limit: usize, stopping_condition: F) -> Result<Vec<u8>>
    where
        F: Fn(&[u8], usize, &Self) -> bool,
    {
        let initial_contig_length = contig.len();

        let mut links_in_scope: Vec<Link> = Vec::new();
        let mut used_links = HashSet::new();
        let mut last_kmer = start_kmer.clone();
    
        loop {
            if stopping_condition(&last_kmer, contig.len() - initial_contig_length, self) {
                return Ok(last_kmer);
            }

            if contig.len() - initial_contig_length >= limit {
                return Err(anyhow::anyhow!("Reached limit without finding a stopping k-mer"));
            }
    
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
    
        Err(anyhow::anyhow!("Stopping condition not satisfied"))
    }

    fn assemble_forward_limit(&self, contig: &mut Vec<u8>, start_kmer: Vec<u8>, limit: usize) -> Result<Vec<u8>> {
        let mut links_in_scope: Vec<Link> = Vec::new();
        let mut used_links = HashSet::new();
        let mut last_kmer = start_kmer.clone();
        loop {
            if contig.len() >= limit {
                return Ok(last_kmer);
            }

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

        Err(anyhow::anyhow!("No stopping k-mer found"))
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
        let mut visited = HashSet::new();

        loop {
            self.update_links(&mut links_in_scope, &last_kmer, &mut used_links, false);

            match self.prev_kmer(&last_kmer, &mut links_in_scope) {
                Some(this_kmer) => {
                    if links_in_scope.len() == 0 {
                        if visited.contains(&this_kmer) {
                            break;
                        }

                        visited.insert(this_kmer.clone());
                    }

                    contig.insert(0, this_kmer[0]);
                    last_kmer = this_kmer;
                }
                None => {
                    break;
                }
            }
        }
    }

    /// Assemble a contig in the backward direction until a stopping k-mer is reached or the limit is hit.
    ///
    /// # Arguments
    ///
    /// * `contig` - A mutable reference to the contig being assembled.
    /// * `start_kmer` - A vector representing the starting k-mer.
    /// * `stop_kmers` - A HashSet of vectors representing the stopping k-mers.
    /// * `limit` - The maximum number of nodes to traverse before giving up.
    fn assemble_backward_until(&self, contig: &mut Vec<u8>, start_kmer: Vec<u8>, stop_kmers: HashSet<Vec<u8>>, limit: usize) -> Result<Vec<u8>> {
        let mut links_in_scope: Vec<Link> = Vec::new();
        let mut used_links = HashSet::new();
        let mut last_kmer = start_kmer.clone();
        let mut nodes_traversed = 0;

        loop {
            if stop_kmers.contains(&last_kmer) {
                return Ok(last_kmer);
            }

            if nodes_traversed >= limit {
                return Err(anyhow::anyhow!("Limit reached before finding a stopping k-mer"));
            }

            self.update_links(&mut links_in_scope, &last_kmer, &mut used_links, false);

            match self.prev_kmer(&last_kmer, &mut links_in_scope) {
                Some(this_kmer) => {
                    contig.insert(0, this_kmer[0]);
                    last_kmer = this_kmer;
                    nodes_traversed += 1;
                }
                None => {
                    break;
                }
            }
        }

        Err(anyhow::anyhow!("No stopping k-mer found"))
    }

    fn assemble_backward_until_condition<F>(&self, contig: &mut Vec<u8>, start_kmer: Vec<u8>, limit: usize, stopping_condition: F) -> Result<Vec<u8>>
    where
        F: Fn(&[u8], usize, &Self) -> bool,
    {
        let initial_contig_length = contig.len();

        let mut links_in_scope: Vec<Link> = Vec::new();
        let mut used_links = HashSet::new();
        let mut last_kmer = start_kmer.clone();
    
        loop {
            if stopping_condition(&last_kmer, contig.len() - initial_contig_length, self) {
                return Ok(last_kmer);
            }

            if contig.len() - initial_contig_length >= limit {
                return Err(anyhow::anyhow!("Reached limit without finding a stopping k-mer"));
            }
    
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
    
        Err(anyhow::anyhow!("Stopping condition not satisfied"))
    }

    fn assemble_backward_limit(&self, contig: &mut Vec<u8>, start_kmer: Vec<u8>, limit: usize) -> Result<Vec<u8>>{
        let mut links_in_scope: Vec<Link> = Vec::new();
        let mut used_links = HashSet::new();
        let mut last_kmer = start_kmer.clone();
        loop {
            if contig.len() >= limit {
                return Ok(last_kmer);
            }

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

        Err(anyhow::anyhow!("No stopping k-mer found"))
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

    /// Assemble all contigs from the linked de Bruijn graph.
    ///
    /// # Returns
    ///
    /// A vector of contigs.
    pub fn assemble_all(&self) -> Vec<Vec<u8>> {
        let mut contigs = Vec::new();

        let mut used_kmers = HashSet::new();
        let k = self.kmer_size;

        let progress_bar = if self.verbose {
            crate::utils::default_bounded_progress_bar("Assembling contigs", self.kmers.len() as u64)
        } else {
            crate::utils::default_hidden_progress_bar()
        };

        for cn_kmer in self.kmers.keys().progress_with(progress_bar) {
            if !used_kmers.contains(cn_kmer) {
                let r = self.kmers.get(cn_kmer).unwrap();

                if r.in_degree() == 1 && r.out_degree() == 1 {
                    let contig = self.assemble(cn_kmer);
                    for kmer_in_contig in contig.windows(k) {
                        used_kmers.insert(crate::utils::canonicalize_kmer(kmer_in_contig));
                    }
                    contigs.push(contig);
                } else if r.in_degree() == 0 && r.out_degree() == 0 {
                    contigs.push(cn_kmer.to_vec());
                }

                used_kmers.insert(cn_kmer.clone());
            }
        }

        contigs
    }

    pub fn assemble_at_bubbles(&self) -> Vec<Vec<u8>> {
        let g = self.traverse_all_kmers();
        let bubbles = find_all_superbubbles(&g);

        let mut all_contigs = Vec::new();
        let mut visited = HashSet::new();

        for ((in_node, out_node), interior) in bubbles {
            let paths_fwd = petgraph::algo::all_simple_paths::<Vec<_>, _>(&g, in_node, out_node, 0, Some(interior.len()));
            let paths_rev = petgraph::algo::all_simple_paths::<Vec<_>, _>(&g, out_node, in_node, 0, Some(interior.len()));

            let paths: Vec<Vec<NodeIndex>> = paths_fwd.chain(paths_rev).collect();

            let mut unique_nodes = Vec::new();

            for (i, path) in paths.iter().enumerate() {
                let other_paths: HashSet<_> = paths.iter().enumerate()
                    .filter(|&(j, _)| j != i)
                    .flat_map(|(_, p)| p)
                    .collect();

                if let Some(unique_node) = path.iter().find(|&node| !other_paths.contains(node)) {
                    let kmer = g.node_weight(*unique_node).unwrap().as_bytes().to_vec();
                    let cn_kmer = crate::utils::canonicalize_kmer(&kmer);

                    if !visited.contains(&cn_kmer) {
                        unique_nodes.push(*unique_node);
                    }
                }
            }

            let contigs = unique_nodes.iter().map(|unique_node| {
                let kmer = g.node_weight(*unique_node).unwrap().as_bytes().to_vec();
                self.assemble(&kmer)
            }).collect::<Vec<_>>();

            all_contigs.extend(contigs.clone());

            for contig in contigs {
                for kmer in contig.windows(self.kmer_size) {
                    visited.insert(crate::utils::canonicalize_kmer(kmer));
                }
            }
        }

        all_contigs
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
        let cn_kmer_vec = crate::utils::canonicalize_kmer(last_kmer.as_bytes()).to_owned();
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

    pub fn clean(self, lenient_threshold: f32, strict_threshold: f32) -> Self {
        let kmer_size = self.kmer_size;

        self
            .clean_color_specific_paths(1, lenient_threshold)
            .clean_tangles(1, 100, lenient_threshold)
            .clean_branches(strict_threshold)
            .clean_tips(3*kmer_size, strict_threshold)
            // .clean_bubbles(2.0*lenient_threshold)
            .clean_contigs(100)
    }

    pub fn clean_color_specific_paths(mut self, color: usize, min_score: f32) -> Self {
        let bad_cn_kmers = self.kmers
            .keys()
            .cloned()
            .filter(|cn_kmer| self.sources.get(cn_kmer).unwrap_or(&vec![]) == &vec![color])
            .filter(|cn_kmer| self.scores.get(cn_kmer).unwrap_or(&1.0) < &min_score)
            .filter(|cn_kmer| self.kmers.get(cn_kmer).unwrap().in_degree() == 1 && self.kmers.get(cn_kmer).unwrap().out_degree() == 1)
            .collect::<Vec<Vec<u8>>>();

        let mut to_remove = HashSet::new();
        let mut bad_paths: usize = 0;

        for bad_cn_kmer in bad_cn_kmers {
            if !to_remove.contains(&bad_cn_kmer) {
                let mut seen_kmers = HashSet::new();
                let mut score_sum = 0.0;

                let mut fw_contig = bad_cn_kmer.clone();
                self.assemble_forward(&mut fw_contig, bad_cn_kmer.clone());

                for kmer in fw_contig.windows(self.kmer_size) {
                    let cn_kmer = crate::utils::canonicalize_kmer(kmer);

                    if let Some(r) = self.kmers.get(&cn_kmer) {
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
                    let cn_kmer = crate::utils::canonicalize_kmer(kmer);

                    if let Some(r) = self.kmers.get(&cn_kmer) {
                        if r.in_degree() == 1 && r.out_degree() == 1 {
                            seen_kmers.insert(cn_kmer.clone());
                            score_sum += self.scores.get(&cn_kmer).unwrap_or(&1.0);
                        } else {
                            break;
                        }
                    }
                }

                let weight = score_sum / seen_kmers.len() as f32;

                if weight < min_score {
                    to_remove.extend(seen_kmers);
                    bad_paths += 1;
                }
            }
        }

        for cn_kmer in &to_remove {
            self.kmers.remove(cn_kmer);
            self.scores.remove(cn_kmer);
        }

        crate::elog!(" -- Removed {} bad paths specific to color {} ({} kmers)", bad_paths, color, to_remove.len());

        self.infer_edges();

        self
    }

    pub fn clean_branches(mut self, min_score: f32) -> Self {
        let bad_cn_kmers = self.kmers
            .keys()
            .cloned()
            .filter(|cn_kmer| self.scores.get(cn_kmer).unwrap_or(&1.0) < &min_score)
            .filter(|cn_kmer| self.kmers.get(cn_kmer).unwrap().in_degree() == 1 && self.kmers.get(cn_kmer).unwrap().out_degree() == 1)
            .collect::<Vec<Vec<u8>>>();

        let mut to_remove = HashSet::new();
        let mut bad_paths: usize = 0;

        for bad_cn_kmer in bad_cn_kmers {
            if !to_remove.contains(&bad_cn_kmer) {
                let mut seen_kmers = HashSet::new();
                let mut score_sum = 0.0;

                let mut fw_contig = bad_cn_kmer.clone();
                self.assemble_forward(&mut fw_contig, bad_cn_kmer.clone());

                for kmer in fw_contig.windows(self.kmer_size) {
                    let cn_kmer = crate::utils::canonicalize_kmer(kmer);

                    if let Some(r) = self.kmers.get(&cn_kmer) {
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
                    let cn_kmer = crate::utils::canonicalize_kmer(kmer);

                    if let Some(r) = self.kmers.get(&cn_kmer) {
                        if r.in_degree() == 1 && r.out_degree() == 1 {
                            seen_kmers.insert(cn_kmer.clone());
                            score_sum += self.scores.get(&cn_kmer).unwrap_or(&1.0);
                        } else {
                            break;
                        }
                    }
                }

                let weight = score_sum / seen_kmers.len() as f32;

                if weight < min_score {
                    to_remove.extend(seen_kmers);
                    bad_paths += 1;
                }
            }
        }

        for cn_kmer in &to_remove {
            self.kmers.remove(cn_kmer);
            self.scores.remove(cn_kmer);
        }

        crate::elog!(" -- Removed {} bad paths ({} kmers)", bad_paths, to_remove.len());

        self.infer_edges();

        self
    }

    pub fn clean_tangles(mut self, color: usize, limit: usize, min_score: f32) -> Self {
        let mut to_remove = HashSet::new();
        let mut bad_tangles: usize = 0;

        let mut visited = HashSet::new();

        let stopping_condition = |kmer: &[u8], _: usize, g: &LdBG| {
            g.scores.get(&crate::utils::canonicalize_kmer(kmer)).unwrap_or(&1.0) < &min_score
        };

        for cn_kmer in self.kmers.keys() {
            if !visited.contains(cn_kmer) && self.sources.get(cn_kmer).unwrap_or(&vec![]) == &vec![color] && self.scores.get(cn_kmer).unwrap_or(&1.0) < &min_score && self.kmers.get(cn_kmer).unwrap().in_degree() + self.kmers.get(cn_kmer).unwrap().out_degree() >= 4 {
                let g = self.traverse_kmers_until_condition(cn_kmer.to_vec(), color, limit, stopping_condition);

                visited.insert(cn_kmer.clone());
                for node in g.node_indices() {
                    let current_kmer = g.node_weight(node).unwrap().as_bytes();
                    let current_cn_kmer = crate::utils::canonicalize_kmer(current_kmer);

                    to_remove.insert(current_cn_kmer.clone());

                    // Mark k-mers as visited, but avoid marking k-mers that are new potential starting points for tangles.
                    if !(self.scores.get(&current_cn_kmer).unwrap_or(&1.0) < &min_score && self.kmers.get(&current_cn_kmer).unwrap().in_degree() + self.kmers.get(&current_cn_kmer).unwrap().out_degree() >= 4) {
                        visited.insert(current_cn_kmer.clone());
                    }

                }

                bad_tangles += 1;
            }
        }

        for cn_kmer in &to_remove {
            self.kmers.remove(cn_kmer);
            self.scores.remove(cn_kmer);
        }

        crate::elog!(" -- Removed {} tangles ({} kmers)", bad_tangles, to_remove.len());

        self.infer_edges();

        self
    }

    pub fn clean_tips(mut self, limit: usize, min_score: f32) -> Self {
        let mut to_remove = HashSet::new();
        let mut bad_paths: usize = 0;

        let stopping_condition = |kmer: &[u8], _: usize, g: &LdBG| { g.prev_kmers(kmer).len() > 1 || g.next_kmers(kmer).len() > 1 };

        for cn_kmer in self.kmers.keys() {
            if self.kmers.get(cn_kmer).unwrap().in_degree() == 0 && self.kmers.get(cn_kmer).unwrap().out_degree() > 1 {
                let mut fw_contig = cn_kmer.to_vec();
                if let Ok(_) = self.assemble_forward_until_condition(&mut fw_contig, cn_kmer.clone(), limit, stopping_condition) {
                    let score_sum = fw_contig.kmers(self.kmer_size as u8)
                        .map(|kmer| self.scores.get(&crate::utils::canonicalize_kmer(&kmer)).unwrap_or(&1.0))
                        .sum::<f32>();
                    let weight = score_sum / (fw_contig.len() - self.kmer_size + 1) as f32;

                    if weight < min_score {
                        for kmer in fw_contig.kmers(self.kmer_size as u8) {
                            to_remove.insert(crate::utils::canonicalize_kmer(&kmer));
                        }

                        bad_paths += 1;
                    }
                }
            }

            if self.kmers.get(cn_kmer).unwrap().out_degree() == 0 && self.kmers.get(cn_kmer).unwrap().in_degree() > 1 {
                let mut rv_contig = cn_kmer.to_vec();
                if let Ok(_) = self.assemble_backward_until_condition(&mut rv_contig, cn_kmer.clone(), limit, stopping_condition) {
                    let score_sum = rv_contig.kmers(self.kmer_size as u8)
                        .map(|kmer| self.scores.get(&crate::utils::canonicalize_kmer(&kmer)).unwrap_or(&1.0))
                        .sum::<f32>();
                    let weight = score_sum / (rv_contig.len() - self.kmer_size + 1) as f32;

                    if weight < min_score {
                        for kmer in rv_contig.kmers(self.kmer_size as u8) {
                            to_remove.insert(crate::utils::canonicalize_kmer(&kmer));
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

        crate::elog!(" -- Removed {} bad tips ({} kmers)", bad_paths, to_remove.len());

        self.infer_edges();

        self
    }

    pub fn clean_bubbles(mut self, min_score: f32) -> Self {
        let mut to_remove = HashSet::new();
        let mut bad_bubbles: usize = 0;

        let g = self.traverse_all_kmers();
        let bubbles = find_all_superbubbles(&g);

        for (i, ((in_node, out_node), interior)) in bubbles.iter().enumerate() {
            let paths_fwd = petgraph::algo::all_simple_paths::<Vec<_>, _>(&g, *in_node, *out_node, 0, Some(interior.len()));
            let paths_rev = petgraph::algo::all_simple_paths::<Vec<_>, _>(&g, *out_node, *in_node, 0, Some(interior.len()));

            let mut paths = paths_fwd
                .chain(paths_rev)
                .map(|path| {
                    let score = path.iter()
                        .map(|node| {
                            let kmer = g.node_weight(*node).unwrap().as_bytes();
                            let cn_kmer = crate::utils::canonicalize_kmer(kmer);
                            *self.scores.get(&cn_kmer).unwrap_or(&1.0)
                        })
                        .sum::<f32>() / path.len() as f32;

                    (path, score)
                })
                .collect::<Vec<(Vec<NodeIndex>, f32)>>();

            paths.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

            if paths.len() > 2 {
                let mut to_keep = HashSet::new();
                let mut cleaned_bubble = false;

                for (j, (path, total_score)) in paths.iter().enumerate() {
                    // let mut contig = String::new();

                    for node in path {
                        let kmer = g.node_weight(*node).unwrap().as_bytes();
                        let cn_kmer = crate::utils::canonicalize_kmer(kmer);

                        if j < 2 {
                            to_keep.insert(cn_kmer.clone());
                        } else if !to_keep.contains(&cn_kmer) && *total_score < min_score {
                            to_remove.insert(cn_kmer.clone());
                            cleaned_bubble = true;
                        }

                        // if contig.is_empty() {
                        //     contig = String::from_utf8_lossy(kmer).to_string();
                        // } else {
                        //     contig.push_str(&String::from_utf8_lossy(&kmer[self.kmer_size - 1..]));
                        // }
                    }

                    // crate::elog!(" -- Bubble {}: {} (score: {}) (remove: {})", i, contig, total_score, num_remove);
                }

                // crate::elog!("");

                if cleaned_bubble {
                    bad_bubbles += paths.len() - 2;
                }
            }
        }

        for cn_kmer in &to_remove {
            self.kmers.remove(cn_kmer);
            self.scores.remove(cn_kmer);
        }

        crate::elog!(" -- Removed {} bad bubbles ({} kmers)", bad_bubbles, to_remove.len());

        self.infer_edges();

        self
    }

    pub fn clean_contigs(mut self, min_contig_length: usize) -> Self {
        let mut to_remove = HashSet::new();
        let mut bad_contigs: usize = 0;

        self.assemble_all().iter().for_each(|contig| {
            if contig.len() < min_contig_length {
                let mut process = true;

                // Check if the first k-mer has no incoming edges
                let first_kmer = contig.windows(self.kmer_size).next().unwrap();
                let prev_kmers = self.prev_kmers(first_kmer);
                if prev_kmers.len() > 0 {
                    process = false;
                }

                // Check if the last k-mer has no outgoing edges
                let last_kmer = contig.windows(self.kmer_size).last().unwrap();
                let next_kmers = self.next_kmers(last_kmer);
                if next_kmers.len() > 0 {
                    process = false;
                }

                if process {
                    for kmer in contig.windows(self.kmer_size) {
                        to_remove.insert(crate::utils::canonicalize_kmer(kmer));
                    }

                    bad_contigs += 1;
                }
            }
        });

        for cn_kmer in &to_remove {
            self.kmers.remove(cn_kmer);
            self.scores.remove(cn_kmer);
        }

        crate::elog!(" -- Removed {} bad contigs ({} kmers)", bad_contigs, to_remove.len());

        self.infer_edges();

        self
    }

    fn traverse_forward(&self, graph: &mut petgraph::Graph<String, f32>, visited: &mut HashMap<String, NodeIndex>, start_node: NodeIndex) {
        let mut fwd_stack = vec![start_node];
        while let Some(node) = fwd_stack.pop() {
            let this_kmer = graph.node_weight(node).unwrap().as_bytes();
    
            for next_kmer in self.next_kmers(this_kmer) {
                let next_label = String::from_utf8_lossy(&next_kmer).to_string();
                let next_node = if let Some(&existing_node) = visited.get(&next_label) {
                    existing_node
                } else {
                    let new_node = graph.add_node(next_label.clone());
                    visited.insert(next_label.clone(), new_node);
                    fwd_stack.push(new_node);
                    new_node
                };
    
                if !graph.contains_edge(node, next_node) {
                    graph.add_edge(
                        node,
                        next_node,
                        *self.scores.get(&crate::utils::canonicalize_kmer(&next_kmer))
                            .unwrap_or(&1.0),
                    );
                }
            }
        }
    }

    fn traverse_forward_until(&self, graph: &mut petgraph::Graph<String, f32>, visited: &mut HashMap<String, NodeIndex>, start_node: NodeIndex, stop_node: &[u8]) -> Result<()> {
        let mut found = false;

        let mut fwd_stack = vec![start_node];
        while let Some(node) = fwd_stack.pop() {
            let this_kmer = graph.node_weight(node).unwrap().as_bytes();
    
            for next_kmer in self.next_kmers(this_kmer) {
                let next_label = String::from_utf8_lossy(&next_kmer).to_string();
                let next_node = if let Some(&existing_node) = visited.get(&next_label) {
                    existing_node
                } else {
                    let new_node = graph.add_node(next_label.clone());
                    visited.insert(next_label.clone(), new_node);

                    if graph.node_weight(new_node).unwrap().as_bytes() != stop_node {
                        fwd_stack.push(new_node);
                    } else {
                        found = true;
                    }

                    new_node
                };
    
                if !graph.contains_edge(node, next_node) {
                    graph.add_edge(
                        node,
                        next_node,
                        *self.scores.get(&crate::utils::canonicalize_kmer(&next_kmer))
                            .unwrap_or(&1.0),
                    );
                }
            }
        }

        if found {
            Ok(())
        } else {
            Err(anyhow::anyhow!("Stop node not found"))
        }
    }

    fn traverse_forward_until_condition<F>(&self, graph: &mut petgraph::Graph<String, f32>, visited: &mut HashMap<String, NodeIndex>, color: usize, start_node: NodeIndex, limit: usize, stopping_condition: F) -> Result<()>
    where
        F: Fn(&[u8], usize, &Self) -> bool,
    {
        let mut found = false;
    
        // Use a vector of tuples (node, depth) instead of just nodes
        let mut fwd_stack = vec![(start_node, 0)];
        while let Some((node, depth)) = fwd_stack.pop() {
            if depth >= limit {
                continue; // Skip this node if we've reached the depth limit
            }
    
            let this_kmer = graph.node_weight(node).unwrap().as_bytes();
    
            for next_kmer in self.next_kmers(this_kmer) {
                let same_color = self.sources.get(&crate::utils::canonicalize_kmer(&next_kmer)).unwrap_or(&vec![]) == &vec![color];
                if !same_color {
                    continue;
                }

                let next_label = String::from_utf8_lossy(&next_kmer).to_string();
                let next_node = if let Some(&existing_node) = visited.get(&next_label) {
                    existing_node
                } else {
                    let new_node = graph.add_node(next_label.clone());
                    visited.insert(next_label.clone(), new_node);
    
                    if !stopping_condition(&next_kmer, depth + 1, self) {
                        // Push the new node with an incremented depth
                        fwd_stack.push((new_node, depth + 1));
                    } else {
                        found = true;
                    }
    
                    new_node
                };
    
                if !graph.contains_edge(node, next_node) {
                    graph.add_edge(
                        node,
                        next_node,
                        *self.scores.get(&crate::utils::canonicalize_kmer(&next_kmer))
                            .unwrap_or(&1.0),
                    );
                }
            }
        }
    
        if found {
            Ok(())
        } else {
            Err(anyhow::anyhow!("Stop condition not met within depth limit"))
        }
    }

    fn traverse_backward(&self, graph: &mut petgraph::Graph<String, f32>, visited: &mut HashMap<String, NodeIndex>, start_node: NodeIndex) {
        let mut rev_stack = vec![start_node];
        while let Some(node) = rev_stack.pop() {
            let this_kmer = graph.node_weight(node).unwrap().as_bytes();
    
            for prev_kmer in self.prev_kmers(this_kmer) {
                let prev_label = String::from_utf8_lossy(&prev_kmer).to_string();
                let prev_node = if let Some(&existing_node) = visited.get(&prev_label) {
                    existing_node
                } else {
                    let new_node = graph.add_node(prev_label.clone());
                    visited.insert(prev_label.clone(), new_node);
                    rev_stack.push(new_node);
                    new_node
                };
    
                if !graph.contains_edge(prev_node, node) {
                    graph.add_edge(
                        prev_node,
                        node,
                        *self.scores.get(&crate::utils::canonicalize_kmer(&prev_kmer))
                            .unwrap_or(&1.0),
                    );
                }
            }
        }
    }

    fn traverse_backward_until(&self, graph: &mut petgraph::Graph<String, f32>, visited: &mut HashMap<String, NodeIndex>, start_node: NodeIndex, stop_node: &[u8]) -> Result<()> {
        let mut found = false;

        let mut rev_stack = vec![start_node];
        while let Some(node) = rev_stack.pop() {
            let this_kmer = graph.node_weight(node).unwrap().as_bytes();
    
            for prev_kmer in self.prev_kmers(this_kmer) {
                let prev_label = String::from_utf8_lossy(&prev_kmer).to_string();
                let prev_node = if let Some(&existing_node) = visited.get(&prev_label) {
                    existing_node
                } else {
                    let new_node = graph.add_node(prev_label.clone());
                    visited.insert(prev_label.clone(), new_node);

                    if graph.node_weight(new_node).unwrap().as_bytes() != stop_node {
                        rev_stack.push(new_node);
                    } else {
                        found = true;
                    }

                    new_node
                };
    
                if !graph.contains_edge(prev_node, node) {
                    graph.add_edge(
                        prev_node,
                        node,
                        *self.scores.get(&crate::utils::canonicalize_kmer(&prev_kmer))
                            .unwrap_or(&1.0),
                    );
                }
            }
        }

        if found {
            Ok(())
        } else {
            Err(anyhow::anyhow!("Stop node not found"))
        }
    }

    fn traverse_backward_until_condition<F>(&self, graph: &mut petgraph::Graph<String, f32>, visited: &mut HashMap<String, NodeIndex>, color: usize, start_node: NodeIndex, limit: usize, stopping_condition: F) -> Result<()>
    where
        F: Fn(&[u8], usize, &Self) -> bool,
    {
        let mut found = false;
    
        // Use a vector of tuples (node, depth) instead of just nodes
        let mut rev_stack = vec![(start_node, 0)];
        while let Some((node, depth)) = rev_stack.pop() {
            if depth >= limit {
                continue; // Skip this node if we've reached the depth limit
            }
    
            let this_kmer = graph.node_weight(node).unwrap().as_bytes();
    
            for prev_kmer in self.prev_kmers(this_kmer) {
                let same_color = self.sources.get(&crate::utils::canonicalize_kmer(&prev_kmer)).unwrap_or(&vec![]) == &vec![color];
                if !same_color {
                    continue;
                }

                let prev_label = String::from_utf8_lossy(&prev_kmer).to_string();
                let prev_node = if let Some(&existing_node) = visited.get(&prev_label) {
                    existing_node
                } else {
                    let new_node = graph.add_node(prev_label.clone());
                    visited.insert(prev_label.clone(), new_node);
    
                    if !stopping_condition(&prev_kmer, depth + 1, self) {
                        // Push the new node with an incremented depth
                        rev_stack.push((new_node, depth + 1));
                    } else {
                        found = true;
                    }
    
                    new_node
                };
    
                if !graph.contains_edge(prev_node, node) {
                    graph.add_edge(
                        prev_node,
                        node,
                        *self.scores.get(&crate::utils::canonicalize_kmer(&prev_kmer))
                            .unwrap_or(&1.0),
                    );
                }
            }
        }
    
        if found {
            Ok(())
        } else {
            Err(anyhow::anyhow!("Stop condition not met within depth limit"))
        }
    }

    pub fn traverse_kmers_until_condition<F>(&self, start_kmer: Vec<u8>, color: usize, limit: usize, stopping_condition: F) -> DiGraph<String, f32>
    where
        F: Fn(&[u8], usize, &Self) -> bool,
    {
        let mut graph = DiGraph::new();
        let mut visited = HashMap::<String, NodeIndex>::new(); // Map node labels to indices

        let start_label = String::from_utf8_lossy(&start_kmer).to_string();
        let start_node = graph.add_node(start_label.clone());
        visited.insert(start_label.clone(), start_node);

        let _ = self.traverse_forward_until_condition(&mut graph, &mut visited, color, start_node, limit, &stopping_condition);
        let _ = self.traverse_backward_until_condition(&mut graph, &mut visited, color, start_node, limit, &stopping_condition);

        graph
    }

    /// Traverse kmers starting from a given kmer and build a graph.
    pub fn traverse_kmers(&self, start_kmer: Vec<u8>) -> DiGraph<String, f32> {
        let mut graph = DiGraph::new();
        let mut visited = HashMap::<String, NodeIndex>::new(); // Map node labels to indices

        let start_label = String::from_utf8_lossy(&start_kmer).to_string();
        let start_node = graph.add_node(start_label.clone());
        visited.insert(start_label.clone(), start_node);

        self.traverse_forward(&mut graph, &mut visited, start_node);
        self.traverse_backward(&mut graph, &mut visited, start_node);

        graph
    }

    pub fn traverse_contigs(&self, start_kmer: Vec<u8>) -> DiGraph<String, f32> {
        let mut graph = DiGraph::new();
        let mut visited = HashMap::new();
    
        let start_contig = self.assemble(&start_kmer);
        let start_node = graph.add_node(String::from_utf8_lossy(&start_contig).to_string());
        
        // Mark all k-mers in the start contig as visited
        for kmer in start_contig.windows(self.kmer_size) {
            visited.insert(crate::utils::canonicalize_kmer(kmer), start_node);
        }
    
        // Traverse forward
        let mut fwd_stack = vec![start_node];
        while let Some(node) = fwd_stack.pop() {
            let this_contig = graph.node_weight(node).unwrap().as_bytes();
    
            if let Some(last_kmer) = self.last_kmer(this_contig) {
                for next_kmer in self.next_kmers(&last_kmer) {
                    if !visited.contains_key(&crate::utils::canonicalize_kmer(&next_kmer)) {
                        let mut next_contig = next_kmer.clone();
                        self.assemble_forward(&mut next_contig, next_kmer.clone());
    
                        let weight = next_contig.windows(self.kmer_size)
                            .map(|x| self.scores.get(&crate::utils::canonicalize_kmer(x)).unwrap_or(&1.0))
                            .fold(f32::INFINITY, |acc, x| acc.min(*x));
    
                        let next_node = graph.add_node(String::from_utf8_lossy(&next_contig).to_string());
                        graph.add_edge(node, next_node, weight);
    
                        // Mark all k-mers in the new contig as visited
                        for kmer in next_contig.windows(self.kmer_size) {
                            visited.insert(crate::utils::canonicalize_kmer(kmer), next_node);
                        }
    
                        fwd_stack.push(next_node);
                    } else {
                        let next_node = *visited.get(&crate::utils::canonicalize_kmer(&next_kmer)).unwrap();
                        if !graph.contains_edge(node, next_node) {
                            let weight = *self.scores.get(&crate::utils::canonicalize_kmer(&next_kmer)).unwrap_or(&1.0);
                            graph.add_edge(node, next_node, weight);
                        }
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
                    if !visited.contains_key(&crate::utils::canonicalize_kmer(&prev_kmer)) {
                        let mut prev_contig = prev_kmer.clone();
                        self.assemble_backward(&mut prev_contig, prev_kmer.clone());
    
                        let weight = prev_contig.windows(self.kmer_size)
                            .map(|x| self.scores.get(&crate::utils::canonicalize_kmer(x)).unwrap_or(&1.0))
                            .fold(f32::INFINITY, |acc, x| acc.min(*x));
    
                        let prev_node = graph.add_node(String::from_utf8_lossy(&prev_contig).to_string());
                        graph.add_edge(prev_node, node, weight);  // Note the reversed edge direction
    
                        // Mark all k-mers in the new contig as visited
                        for kmer in prev_contig.windows(self.kmer_size) {
                            visited.insert(crate::utils::canonicalize_kmer(kmer), prev_node);
                        }
    
                        rev_stack.push(prev_node);
                    } else {
                        let prev_node = *visited.get(&crate::utils::canonicalize_kmer(&prev_kmer)).unwrap();
                        if !graph.contains_edge(prev_node, node) {
                            let weight = *self.scores.get(&crate::utils::canonicalize_kmer(&prev_kmer)).unwrap_or(&1.0);
                            graph.add_edge(prev_node, node, weight);  // Note the reversed edge direction
                        }
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
        let mut node_indices = HashMap::<String, NodeIndex>::new(); // Map node labels to indices

        for cn_kmer in cn_kmers {
            if !visited.contains(&cn_kmer) {
                let g = self.traverse_kmers(cn_kmer.clone());

                // Mark kmers as visited
                for node_label in g.node_weights() {
                    let cn_kmer = crate::utils::canonicalize_kmer(node_label.as_bytes());
                    visited.insert(cn_kmer);
                }

                // Add nodes to the main graph, avoiding duplicates
                for node_index in g.node_indices() {
                    let node_label = g.node_weight(node_index).unwrap().clone();
                    if !node_indices.contains_key(&node_label) {
                        let new_node_idx = graph.add_node(node_label.clone());
                        node_indices.insert(node_label.clone(), new_node_idx);
                    }
                }

                // Add edges, mapping nodes correctly
                for edge in g.edge_references() {
                    let (source, target) = (edge.source(), edge.target());
                    let source_label = g.node_weight(source).unwrap();
                    let target_label = g.node_weight(target).unwrap();

                    let new_source = *node_indices.get(source_label).expect("Source node should exist");
                    let new_target = *node_indices.get(target_label).expect("Target node should exist");

                    if !graph.contains_edge(new_source, new_target) {
                        graph.add_edge(new_source, new_target, *edge.weight());
                    }
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
                        visited.insert(crate::utils::canonicalize_kmer(&kmer));
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

fn traverse_read_graph(graph: &petgraph::Graph<String, f64>, min_length: usize) -> Vec<Vec<u8>> {
    // Populate corrected_seqs with all contiguous stretches of the graph more than 1 node long
    let mut corrected_seqs = Vec::new();
    let mut visited = HashSet::new();
    for node in graph.node_indices() {
        if !visited.contains(&node) {
            let mut current_seq = Vec::new();
            let mut current_node = node;
        
            loop {
                visited.insert(current_node);
                if current_seq.is_empty() {
                    current_seq.push(graph[current_node].clone());
                } else {
                    current_seq.push(graph[current_node].chars().last().unwrap().to_string());
                }
            
                let outgoing_edges: Vec<_> = graph.edges_directed(current_node, petgraph::Direction::Outgoing).collect();
                if outgoing_edges.len() == 1 {
                    let next_node = outgoing_edges[0].target();
                    if !visited.contains(&next_node) {
                        current_node = next_node;
                    } else {
                        break;
                    }
                } else {
                    break;
                }
            }
        
            if current_seq.len() > min_length {
                let contig = current_seq.join("").as_bytes().to_vec();
                corrected_seqs.push(contig);
            }
        }
    }
    corrected_seqs
}

fn navigate_backward(node: NodeIndex, steps: u8, graph: &petgraph::Graph<String, f64>) -> Option<NodeIndex> {
    let mut current_node = node;
    let mut steps_backward = 0;

    while steps_backward < steps {
        let incoming_edges: Vec<_> = graph.edges_directed(current_node, petgraph::Direction::Incoming).collect();
        if incoming_edges.is_empty() {
            return None;
        }
        current_node = incoming_edges[0].source();
        steps_backward += 1;
    }

    if steps_backward == steps {
        Some(current_node)
    } else {
        None
    }
}

fn navigate_forward(node: NodeIndex, steps: u8, graph: &petgraph::Graph<String, f64>) -> Option<NodeIndex> {
    let mut current_node = node;
    let mut steps_forward = 0;

    while steps_forward < steps {
        let outgoing_edges: Vec<_> = graph.edges_directed(current_node, petgraph::Direction::Outgoing).collect();
        if outgoing_edges.is_empty() {
            return None;
        }
        current_node = outgoing_edges[0].source();
        steps_forward += 1;
    }

    if steps_forward == steps {
        Some(current_node)
    } else {
        None
    }
}

// See reference implementation at
// https://github.com/fawaz-dabbaghieh/bubble_gun/blob/53b8d68d0c9d0c35252da4dc1b4bbb6af57b4631/BubbleGun/find_bubbles.py#L29
pub fn find_superbubble(graph: &DiGraph<String, f32>, s: NodeIndex, direction: petgraph::Direction) -> Option<(NodeIndex, NodeIndex, Vec<NodeIndex>)> {
    let mut seen = HashSet::new();
    let mut visited = HashSet::new();
    let mut nodes_inside = Vec::new();

    // seen.insert((s.index(), direction));
    seen.insert(s.index());

    let mut stack = vec![(s, direction)];
    while !stack.is_empty() {

        let (v, v_direction) = stack.pop().unwrap();
        visited.insert(v.index());

        nodes_inside.push(v);

        // seen.remove(&(v.index(), v_direction));
        seen.remove(&v.index());

        let children = match v_direction {
            petgraph::Direction::Incoming => graph.neighbors_directed(v, petgraph::Direction::Incoming).collect::<Vec<_>>(),
            petgraph::Direction::Outgoing => graph.neighbors_directed(v, petgraph::Direction::Outgoing).collect::<Vec<_>>(),
        };

        if children.is_empty() {
            break;
        }

        for u in children {
            let (u_child_direction, u_parents) = if v_direction == petgraph::Direction::Outgoing {
                (petgraph::Direction::Outgoing, graph.neighbors_directed(u, petgraph::Direction::Incoming).collect::<Vec<_>>())
            } else {
                (petgraph::Direction::Incoming, graph.neighbors_directed(u, petgraph::Direction::Outgoing).collect::<Vec<_>>())
            };

            if u.index() == s.index() {
                stack.clear();
                break;
            }

            seen.insert(u.index());

            if u_parents.iter().all(|&n| visited.contains(&n.index())) {
                stack.push((u, u_child_direction));
            }
        }

        if stack.len() == 1 && seen.len() == 1 {
            let (t, _) = stack.pop().unwrap();
            nodes_inside.push(t);

            if nodes_inside.len() == 2 {
                break;
            }

            nodes_inside.retain(|&x| x != s && x != t);
            return Some((s, t, nodes_inside));
        }
    }

    None
}

pub fn find_all_superbubbles(graph: &petgraph::Graph<String, f32>) -> HashMap<(NodeIndex, NodeIndex), Vec<NodeIndex>> {
    let mut bubbles = HashMap::new();

    let mut visited: HashSet<NodeIndex> = HashSet::new();
    for n in graph.node_indices() {
        if !visited.contains(&n) {
            for d in [petgraph::Direction::Outgoing, petgraph::Direction::Incoming] {
                let bubble = find_superbubble(&graph, n, d);

                if let Some(bubble) = bubble {
                    visited.extend(bubble.2.iter().cloned());

                    let key_fwd = (bubble.0, bubble.1);
                    let key_rev = (bubble.1, bubble.0);

                    if !bubbles.contains_key(&key_fwd) && !bubbles.contains_key(&key_rev) {
                        bubbles.insert((bubble.0, bubble.1), bubble.2);
                    }
                }
            }
        }
    }

    bubbles
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::edges::Edges;

    use petgraph::data::DataMap;
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
                    let cn_anchor_vec = crate::utils::canonicalize_kmer(anchor_kmer);
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

        let g = LdBG::from_sequences(String::from("test"), 5, &fwd_seqs);
        let graph = g.traverse_kmers(b"ATTTC".to_vec());

        // println!("{}", Dot::with_config(&graph, &[petgraph::dot::Config::EdgeNoLabel]));

        // Check if the graph is consistent with the expected structure
        assert_eq!(graph.node_count(), 21);
        assert_eq!(graph.edge_count(), 21);

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
        assert_eq!(edges.len(), 21);

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

    // TODO: This test is not working yet
    fn test_traverse_all_kmers() {
        let genome = get_test_genome();
        let fwd_seqs = vec![genome];

        let g = LdBG::from_sequences(String::from("test"), 5, &fwd_seqs);
        let graph = g.traverse_all_kmers();

        // Write graph as GFA to a string
        let mut gfa_output = Vec::new();
        crate::utils::write_gfa(&mut gfa_output, &graph).unwrap();

        // Print GFA string (commented out for test)
        let gfa_string = String::from_utf8(gfa_output).unwrap();
        println!("{}", gfa_string);
    }

    // TODO: This test is not working yet
    fn test_traverse_all_contigs() {
        let genome = get_test_genome();
        let fwd_seqs = vec![genome];

        let g = LdBG::from_sequences(String::from("test"), 5, &fwd_seqs);
        let graph = g.traverse_all_contigs();

        // Write graph as GFA to a string
        let mut gfa_output = Vec::new();
        crate::utils::write_gfa(&mut gfa_output, &graph).unwrap();

        // Print GFA string (commented out for test)
        let gfa_string = String::from_utf8(gfa_output).unwrap();
        println!("{}", gfa_string);
    }

    #[test]
    fn test_traverse_contigs() {
        let genome = get_test_genome();
        let fwd_seqs = vec![genome];

        let g = LdBG::from_sequences(String::from("test"), 5, &fwd_seqs);
        let graph = g.traverse_contigs(b"CCACG".to_vec());

        // Uncomment to manually verify structure of the graph
        // use petgraph::dot::Dot;
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
        assert!(graph.contains_edge(find_node("ACTGATTTCGAT"), find_node("CGATGCCACGGTGG")));
        assert!(graph.contains_edge(find_node("ACTGATTTCGAT"), find_node("CGATGCGAT")));
        assert!(graph.contains_edge(find_node("CGATGCGAT"), find_node("CGATGCCACGGTGG")));
        assert!(graph.contains_edge(find_node("CGATGCGAT"), find_node("CGATGCGAT")));
    }

    #[test]
    fn test_prev_kmers() {
        let genome = get_test_genome();
        let fwd_seqs = vec![genome];

        let g = LdBG::from_sequences(String::from("test"), 5, &fwd_seqs);

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

        let g = LdBG::from_sequences(String::from("test"), 5, &fwd_seqs);

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

        let g = LdBG::from_sequences(String::from("test"), 5, &fwd_seqs);

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

        let g = LdBG::from_sequences(String::from("test"), 5, &fwd_seqs)
            .build_links(&fwd_seqs, false);

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
    fn test_assemble_until() {
        let fw_genome = get_test_genome();
        let fwd_seqs = vec![fw_genome.clone()];

        let g = LdBG::from_sequences(String::from("test"), 5, &fwd_seqs)
            .build_links(&fwd_seqs, false);

        let mut contig1 = b"ACTGA".to_vec();
        let _ = g.assemble_forward_until(&mut contig1, b"ACTGA".to_vec(), HashSet::from([b"TTCGA".to_vec()]), 100);
        assert_eq!(contig1, b"ACTGATTTCGA".to_vec());

        let mut contig2 = b"GGTGG".to_vec();
        let _ = g.assemble_backward_until(&mut contig2, b"GGTGG".to_vec(), HashSet::from([b"TGCCA".to_vec()]), 100);
        assert_eq!(contig2, b"TGCCACGGTGG".to_vec());
    }

    #[test]
    fn test_correct_seq() {
        let genome = get_test_genome();
        let g = LdBG::from_sequence("test".to_string(), 5, &genome);

        let mut uncorrected_seq_map = BTreeMap::new();

        // Things that should get corrected successfully.
        uncorrected_seq_map.insert(b"ACTGATTTCGATGCGATGCGATGCCACGGTGG".to_vec(), vec![b"ACTGATTTCGATGCGATGCGATGCCACGGTGG".to_vec()]);
        uncorrected_seq_map.insert(b"ACTGAATTCGATGCGATGCGATGCCACGGTGG".to_vec(), vec![b"ACTGATTTCGATGCGATGCGATGCCACGGTGG".to_vec()]);
        uncorrected_seq_map.insert(b"ACTGATATCGATGCGATGCGATGCCACGGTGG".to_vec(), vec![b"ACTGATTTCGATGCGATGCGATGCCACGGTGG".to_vec()]);
        uncorrected_seq_map.insert(b"ACTGATTTCGATGCGATGCGATGCCATGGTGG".to_vec(), vec![b"ACTGATTTCGATGCGATGCGATGCCACGGTGG".to_vec()]);
        uncorrected_seq_map.insert(b"ACTGATTTCGATGCGATGCGATGCCATCGGTGG".to_vec(), vec![b"ACTGATTTCGATGCGATGCGATGCCACGGTGG".to_vec()]);
        uncorrected_seq_map.insert(b"ACTGATTTCGATGCGATGCGATGCCGGTGG".to_vec(), vec![b"ACTGATTTCGATGCGATGCGATGCCACGGTGG".to_vec()]);

        for (uncorrected_seq, expected_seqs) in uncorrected_seq_map.iter() {
            let corrected_seqs = g.correct_seq(&uncorrected_seq);

            assert_eq!(corrected_seqs.len(), expected_seqs.len());
            for (corrected_seq, expected_seq) in corrected_seqs.iter().zip(expected_seqs) {
                assert_eq!(*corrected_seq, *expected_seq);
            }
        }
    }

    proptest! {
        #[test]
        fn test_assemble_random_genomes(
            repeat_length in (3..18usize).prop_filter("Increment by 3", |v| v % 3 == 0),
            num_repeats in (2..=3usize),
            k in (11..=17usize).prop_filter("Must be odd", |v| v % 2 == 1)
        ) {
            let random_genome = generate_genome_with_tandem_repeats(500, repeat_length, num_repeats, 0);

            let g = LdBG::from_sequences(String::from("test"), k, &vec!(random_genome.clone()))
                .build_links(&vec!(random_genome.clone()), false);

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

    // From "BubbleGun: enumerating bubbles and superbubbles in genome graphs", Dabbaghie et al. 2022
    // https://academic.oup.com/bioinformatics/article/38/17/4217/6633304?login=false
    // https://github.com/fawaz-dabbaghieh/bubble_gun/blob/master/example/paper_example2.gfa
    fn create_bubblegun_graph() -> petgraph::Graph<String, f32> {
        // Create a new directed graph
        let mut graph = petgraph::Graph::<String, f32>::new();
        
        // Create a HashMap to store node indices
        let mut node_indices = std::collections::HashMap::new();
    
        // Add nodes to the graph
        let gfa_nodes = vec![
            ("1", "CCCAACAAGTG"),
            ("7", "AACAAGTGTACTCATTG"),
            ("25", "ACTCATTGACG"),
            ("31", "CATTGACGTAAATTTGGTGCGGGCCTCAAGGTGTCCA"),
            ("89", "GGTGTCCATTGGGGTC"),
            ("105", "TTGGGGTCAGCACAAAT"),
            ("123", "GCACAAATTGCCA"),
            ("163", "CATTGACGAGGCACCC"),
            ("179", "AGGCACCCGTGCATTAG"),
            ("197", "TGCATTAGGCAGGGTGT"),
            ("215", "CAGGGTGTCCA"),
            ("237", "TTGGGGTCTGCACAAAT"),
            ("271", "AACAAGTGAACTCATTG"),
            ("311", "AGGCACCCCTGCATTAG"),
            ("427", "CATTGACGCAACCGGCATTGAATACACAGGGTGT"),
        ];
    
        for (id, seq) in gfa_nodes {
            let node_index = graph.add_node(seq.to_string());
            node_indices.insert(id.to_string(), node_index);
        }
    
        // Add edges to the graph
        let gfa_edges = vec![
            ("1", "7"), ("1", "271"),
            ("7", "25"),
            ("25", "31"), ("25", "427"), ("25", "163"),
            ("31", "89"),
            ("89", "237"), ("89", "105"),
            ("105", "123"),
            ("163", "179"), ("163", "311"),
            ("179", "197"),
            ("197", "215"),
            ("215", "89"),
            ("237", "123"),
            ("271", "25"),
            ("311", "197"),
            ("427", "215"),
        ];
    
        for (from, to) in gfa_edges {
            if let (Some(&from_index), Some(&to_index)) = (node_indices.get(from), node_indices.get(to)) {
                graph.add_edge(from_index, to_index, 1.0);
            }
        }
        graph
    }

    #[test]
    fn test_find_all_superbubbles() {
        let graph = create_bubblegun_graph();

        let expected_bubbles = HashMap::from([
            ((b"CCCAACAAGTG".to_vec(),      b"ACTCATTGACG".to_vec()),      2usize),
            // ((b"ACTCATTGACG".to_vec(),      b"CCCAACAAGTG".to_vec()),      2usize),

            ((b"ACTCATTGACG".to_vec(),      b"GGTGTCCATTGGGGTC".to_vec()), 7usize),
            // ((b"GGTGTCCATTGGGGTC".to_vec(), b"ACTCATTGACG".to_vec()),      7usize),

            ((b"GGTGTCCATTGGGGTC".to_vec(), b"GCACAAATTGCCA".to_vec()),    2usize),
            // ((b"GCACAAATTGCCA".to_vec(),    b"GGTGTCCATTGGGGTC".to_vec()), 2usize),
        ]);

        let bubbles = find_all_superbubbles(&graph);
        assert_eq!(bubbles.len(), expected_bubbles.len());

        for ((in_node, out_node), interior) in bubbles {
            let in_seq = graph.node_weight(in_node).unwrap().as_bytes().to_vec();
            let out_seq = graph.node_weight(out_node).unwrap().as_bytes().to_vec();

            let key_fwd = (in_seq.clone(), out_seq.clone());
            let key_rev = (out_seq.clone(), in_seq.clone());

            assert!(expected_bubbles.contains_key(&key_fwd) || expected_bubbles.contains_key(&key_rev));

            if expected_bubbles.contains_key(&key_fwd) {
                assert_eq!(expected_bubbles.get(&key_fwd).unwrap(), &interior.len());
            } else {
                assert_eq!(expected_bubbles.get(&key_rev).unwrap(), &interior.len());
            }
        }
    }

    #[test]
    fn test_find_all_superbubbles_2() {
        let graph = crate::utils::read_gfa("/Users/kiran/repositories/hidive/test.gfa").unwrap();

        let bubbles = find_all_superbubbles(&graph);

        let start = b"GAGGGAACAGCGACTTC".to_vec();
        let end = b"TGGACACAGGCACCTGG".to_vec();
        for ((in_node, out_node), interior) in bubbles {
            let in_seq = graph.node_weight(in_node).unwrap().as_bytes().to_vec();
            let out_seq = graph.node_weight(out_node).unwrap().as_bytes().to_vec();

            if in_seq == start && out_seq == end {
                println!("Found bubble: {}", interior.len());

                for path in petgraph::algo::all_simple_paths::<Vec<_>, _>(&graph, in_node, out_node, 0, Some(interior.len() + 10)) {
                    println!("out path: {}", path.len());
                }

                for path in petgraph::algo::all_simple_paths::<Vec<_>, _>(&graph, out_node, in_node, 0, Some(interior.len() + 10)) {
                    println!("in path:  {}", path.len());

                    for node in path {
                        println!(" - {}", graph.node_weight(node).unwrap());
                    }
                }
            }
        }
    }
}
