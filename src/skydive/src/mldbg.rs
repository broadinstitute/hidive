use std::{collections::HashSet, path::PathBuf};

use crate::ldbg::LdBG;

/// Represents a multi-color linked de Bruijn graph, all built with same k-mer size.
#[derive(Debug)]
pub struct MLdBG {
    pub ldbgs: Vec<LdBG>,
    pub kmer_size: usize,
    pub build_links: bool
}

impl MLdBG {
    /// Create an empty multi-color LdBG.
    pub fn new(kmer_size: usize, build_links: bool) -> Self {
        MLdBG {
            ldbgs: Vec::new(),
            kmer_size,
            build_links
        }
    }

    /// Add a LdBG to the MLdBG.
    pub fn push(&mut self, ldbg: LdBG) {
        assert!(ldbg.kmer_size == self.kmer_size, "The k-mer size of the LdBG does not match the k-mer size of the MLdBG.");

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
        let l = LdBG::from_file(name, self.kmer_size, seq_path, self.build_links);
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

        let l = LdBG::from_sequences(name, self.kmer_size, &filtered_reads, self.build_links);
        self.ldbgs.push(l);
    }

    /// Get the union of kmers from all LdBGs in the MLdBG.
    fn union_of_kmers(&self) -> HashSet<Vec<u8>> {
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
        let index = self.ldbgs.iter().position(|x| condition(x));
        if let Some(index) = index {
            Some(self.ldbgs.remove(index))
        } else {
            None
        }
    }
}