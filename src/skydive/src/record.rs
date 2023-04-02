use std::fmt;

use crate::edges::Edges;

/// Represents a de Bruijn graph record.
#[derive(Debug)]
pub struct Record {
    coverage: u16,
    edges: Edges
}

impl Record {
    /// Create an empty de Bruijn graph record.
    pub fn new(coverage: u16) -> Self {
        Record {
            coverage: coverage,
            edges: Edges::empty()
        }
    }

    /// Return the Record's coverage.
    pub fn coverage(&self) -> u16 {
        self.coverage
    }

    /// Increment the coverage value by 1.
    pub fn increment_coverage(&mut self) {
        self.coverage = self.coverage.saturating_add(1);
    }

    /// Set the coverage value.
    pub fn set_coverage(&mut self, coverage: u16) {
        self.coverage = coverage;
    }

    /// Set incoming edge.
    pub fn set_incoming_edge(&mut self, nucleotide: u8) {
        match nucleotide {
            b'A' => self.edges.insert(Edges::FLAG_EDGE_IN_A),
            b'C' => self.edges.insert(Edges::FLAG_EDGE_IN_C),
            b'G' => self.edges.insert(Edges::FLAG_EDGE_IN_G),
            b'T' => self.edges.insert(Edges::FLAG_EDGE_IN_T),
            _ => println!("Invalid nucleotide: {}", nucleotide as char)
        }
    }

    /// Set outgoing edge.
    pub fn set_outgoing_edge(&mut self, nucleotide: u8) {
        match nucleotide {
            b'A' => self.edges.insert(Edges::FLAG_EDGE_OUT_A),
            b'C' => self.edges.insert(Edges::FLAG_EDGE_OUT_C),
            b'G' => self.edges.insert(Edges::FLAG_EDGE_OUT_G),
            b'T' => self.edges.insert(Edges::FLAG_EDGE_OUT_T),
            _ => println!("Invalid nucleotide: {}", nucleotide as char)
        }
    }
}

impl fmt::Display for Record {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}{}{}{} {}{}{}{}] {}",
            if self.edges.contains(Edges::FLAG_EDGE_IN_A)  { 'a' } else { '.' },
            if self.edges.contains(Edges::FLAG_EDGE_IN_C)  { 'c' } else { '.' },
            if self.edges.contains(Edges::FLAG_EDGE_IN_G)  { 'g' } else { '.' },
            if self.edges.contains(Edges::FLAG_EDGE_IN_T)  { 't' } else { '.' },

            if self.edges.contains(Edges::FLAG_EDGE_OUT_A) { 'A' } else { '.' },
            if self.edges.contains(Edges::FLAG_EDGE_OUT_C) { 'C' } else { '.' },
            if self.edges.contains(Edges::FLAG_EDGE_OUT_G) { 'G' } else { '.' },
            if self.edges.contains(Edges::FLAG_EDGE_OUT_T) { 'T' } else { '.' },

            self.coverage,
        )
    }
}