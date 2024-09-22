use std::{collections::{HashMap, HashSet}, path::PathBuf};

use gbdt::{decision_tree::Data, gradient_boost::GBDT};

use crate::{ldbg::LdBG, record::Record};

/// Represents a multi-color linked de Bruijn graph, all built with same k-mer size.
#[derive(Debug)]
pub struct MLdBG {
    pub kmer_size: usize,
    pub ldbgs: Vec<LdBG>,
    pub scores: HashMap<Vec<u8>, f32>,
}

impl MLdBG {
    /// Create an empty multi-color LdBG.
    pub fn new(kmer_size: usize) -> Self {
        MLdBG {
            kmer_size,
            ldbgs: Vec::new(),
            scores: HashMap::new(),
        }
    }

    /// Create an MLdBG from a vector of LdBGs.
    pub fn from_ldbgs(ldbgs: Vec<LdBG>) -> Self {
        let kmer_size = ldbgs[0].kmer_size;

        for ldbg in &ldbgs {
            assert!(
                ldbg.kmer_size == kmer_size,
                "The k-mer size of the LdBG does not match the k-mer size of the MLdBG."
            );
        }

        MLdBG {
            kmer_size,
            ldbgs,
            scores: HashMap::new(),
        }
    }

    /// Add a LdBG to the MLdBG.
    pub fn push(&mut self, ldbg: LdBG) {
        assert!(
            ldbg.kmer_size == self.kmer_size,
            "The k-mer size of the LdBG does not match the k-mer size of the MLdBG."
        );

        self.ldbgs.push(ldbg);
    }

    /// Insert a LdBG at a specific position in the MLdBG.
    pub fn insert(&mut self, index: usize, ldbg: LdBG) {
        if index <= self.ldbgs.len() {
            self.ldbgs.insert(index, ldbg);
        }
    }

    /// Append a LdBG to the end of the MLdBG.
    pub fn append(&mut self, ldbg: LdBG) {
        self.ldbgs.push(ldbg);
    }

    /// Append an LdBG to the end of the MLdBG, created anew from a fasta file.
    pub fn append_from_file(&mut self, name: String, seq_path: &PathBuf) {
        let l = LdBG::from_file(name, self.kmer_size, seq_path);
        self.ldbgs.push(l);
    }

    pub fn append_from_filtered_file<F>(&mut self, name: String, seq_path: &PathBuf, filter: F)
    where
        F: Fn(&bio::io::fasta::Record, &HashSet<Vec<u8>>) -> bool,
    {
        let reader = bio::io::fasta::Reader::from_file(seq_path).unwrap();
        let all_reads: Vec<bio::io::fasta::Record> = reader.records().map(|r| r.unwrap()).collect();

        let kmer_union = self.union_of_kmers();

        let filtered_reads: Vec<Vec<u8>> = all_reads
            .into_iter()
            .filter(|r| filter(r, &kmer_union))
            .map(|r| r.seq().to_vec())
            .collect();

        let l = LdBG::from_sequences(
            name,
            self.kmer_size,
            &filtered_reads,
        );
        self.ldbgs.push(l);
    }

    pub fn score_kmers(mut self, model_path: &PathBuf) -> Self {
        let gbdt = GBDT::load_model(model_path.to_str().unwrap()).unwrap();

        self.scores = self.union_of_kmers().iter().map(|cn_kmer| {
            let mut features = vec![];

            for ldbg in &self.ldbgs {
                let cov = ldbg.kmers.get(cn_kmer).map_or(0, |record| record.coverage());
                features.push(cov as f32);
            }

            let compressed_len = crate::utils::homopolymer_compressed(cn_kmer).len();
            features.push(compressed_len as f32);

            let data = Data::new_test_data(features, Some(0.0));
            let prediction = *gbdt.predict(&vec![data]).first().unwrap_or(&0.0);

            (cn_kmer.clone(), prediction)
        }).collect();

        self
    }

    pub fn collapse(&mut self) -> LdBG {
        let mut ldbg = LdBG::new(self.ldbgs[0].name.clone(), self.kmer_size);

        for cn_kmer in self.union_of_kmers() {
            let coverage = self.ldbgs
                .iter()
                .map(|ldbg|
                    ldbg.kmers
                        .get(&cn_kmer)
                        .map_or(0, |record| record.coverage())
                )
                .sum::<u16>();

            ldbg.kmers.insert(cn_kmer, Record::new(coverage, None));
        }

        ldbg.scores = self.scores.clone();

        ldbg.infer_edges();

        ldbg
    }

    pub fn filter_reads<F>(&mut self, seq_path: &PathBuf, filter: F) -> Vec<Vec<u8>>
    where
        F: Fn(&bio::io::fasta::Record, &HashSet<Vec<u8>>) -> bool,
    {
        let reader = bio::io::fasta::Reader::from_file(seq_path).unwrap();
        let all_reads: Vec<bio::io::fasta::Record> = reader.records().map(|r| r.unwrap()).collect();

        let kmer_union = self.union_of_kmers();

        let filtered_reads: Vec<Vec<u8>> = all_reads
            .into_iter()
            .filter(|r| filter(r, &kmer_union))
            .map(|r| r.seq().to_vec())
            .collect();

        filtered_reads
    }

    /// Get the union of kmers from all LdBGs in the MLdBG.
    pub fn union_of_kmers(&self) -> HashSet<Vec<u8>> {
        let mut kmer_union = HashSet::new();

        for ldbg in &self.ldbgs {
            for kmer in ldbg.kmers.keys() {
                kmer_union.insert(kmer.clone());
            }
        }

        kmer_union
    }

    /// Get a reference to the LdBG at a specific index.
    pub fn get(&self, index: usize) -> Option<&LdBG> {
        self.ldbgs.get(index)
    }

    /// Returns an iterator over the LdBGs in the MLdBG.
    pub fn iter(&self) -> std::slice::Iter<LdBG> {
        self.ldbgs.iter()
    }

    /// Returns a mutable iterator over the LdBGs in the MLdBG.
    pub fn iter_mut(&mut self) -> std::slice::IterMut<LdBG> {
        self.ldbgs.iter_mut()
    }

    /// Clear all LdBGs from the MLdBG.
    pub fn clear(&mut self) {
        self.ldbgs.clear();
    }

    /// Remove a LdBG from the MLdBG by index.
    pub fn remove(&mut self, index: usize) -> Option<LdBG> {
        if index < self.ldbgs.len() {
            Some(self.ldbgs.remove(index))
        } else {
            None
        }
    }

    /// Returns the number of LdBGs in the MLdBG.
    pub fn len(&self) -> usize {
        self.ldbgs.len()
    }

    /// Check if the MLdBG is empty.
    pub fn is_empty(&self) -> bool {
        self.ldbgs.is_empty()
    }

    /// Remove and return the last LdBG from the MLdBG.
    pub fn pop(&mut self) -> Option<LdBG> {
        self.ldbgs.pop()
    }

    /// Remove and return the LdBG from the MLdBG if it matches a certain condition.
    pub fn pop_if<F>(&mut self, condition: F) -> Option<LdBG>
    where
        F: Fn(&LdBG) -> bool,
    {
        let index = self.ldbgs.iter().position(condition);
        if let Some(index) = index {
            Some(self.ldbgs.remove(index))
        } else {
            None
        }
    }
}
