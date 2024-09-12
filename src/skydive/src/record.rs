use std::fmt;

use crate::edges::Edges;

/// Represents a de Bruijn graph record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    coverage: u16,
    edges: Edges,
}

impl Record {
    /// Create an empty de Bruijn graph record.
    #[must_use]
    pub fn new(coverage: u16, edges: Option<Edges>) -> Self {
        Record {
            coverage,
            edges: edges.unwrap_or(Edges::empty()),
        }
    }

    // /// Create a default de Bruijn graph record with 0 coverage and no edges.
    // pub fn default() -> Self {
    //     Record {
    //         coverage: 0,
    //         edges: Edges::empty(),
    //     }
    // }

    /// Return the Record's coverage.
    #[must_use]
    pub fn coverage(&self) -> u16 {
        self.coverage
    }

    // Return all edges.
    #[must_use]
    pub fn edges(&self) -> Edges {
        self.edges
    }

    /// Return incoming edges
    #[must_use]
    pub fn incoming_edges(&self) -> Vec<u8> {
        let mut edges = Vec::new();

        if self.edges.contains(Edges::FLAG_EDGE_IN_A) {
            edges.push(b'A');
        }
        if self.edges.contains(Edges::FLAG_EDGE_IN_C) {
            edges.push(b'C');
        }
        if self.edges.contains(Edges::FLAG_EDGE_IN_G) {
            edges.push(b'G');
        }
        if self.edges.contains(Edges::FLAG_EDGE_IN_T) {
            edges.push(b'T');
        }

        edges
    }

    /// Return outgoing edges
    #[must_use]
    pub fn outgoing_edges(&self) -> Vec<u8> {
        let mut edges = Vec::new();

        if self.edges.contains(Edges::FLAG_EDGE_OUT_A) {
            edges.push(b'A');
        }
        if self.edges.contains(Edges::FLAG_EDGE_OUT_C) {
            edges.push(b'C');
        }
        if self.edges.contains(Edges::FLAG_EDGE_OUT_G) {
            edges.push(b'G');
        }
        if self.edges.contains(Edges::FLAG_EDGE_OUT_T) {
            edges.push(b'T');
        }

        edges
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
            _ => (),
        }
    }

    /// Set outgoing edge.
    pub fn set_outgoing_edge(&mut self, nucleotide: u8) {
        match nucleotide {
            b'A' => self.edges.insert(Edges::FLAG_EDGE_OUT_A),
            b'C' => self.edges.insert(Edges::FLAG_EDGE_OUT_C),
            b'G' => self.edges.insert(Edges::FLAG_EDGE_OUT_G),
            b'T' => self.edges.insert(Edges::FLAG_EDGE_OUT_T),
            _ => (),
        }
    }

    /// The in-degree of a particular k-mer.
    #[must_use]
    pub fn in_degree(&self) -> u8 {
        let mut degree: u8 = 0;
        degree += if self.edges.contains(Edges::FLAG_EDGE_IN_A) {
            1
        } else {
            0
        };
        degree += if self.edges.contains(Edges::FLAG_EDGE_IN_C) {
            1
        } else {
            0
        };
        degree += if self.edges.contains(Edges::FLAG_EDGE_IN_G) {
            1
        } else {
            0
        };
        degree += if self.edges.contains(Edges::FLAG_EDGE_IN_T) {
            1
        } else {
            0
        };

        degree
    }

    /// The out-degree of a particular k-mer.
    #[must_use]
    pub fn out_degree(&self) -> u8 {
        let mut degree: u8 = 0;
        degree += if self.edges.contains(Edges::FLAG_EDGE_OUT_A) {
            1
        } else {
            0
        };
        degree += if self.edges.contains(Edges::FLAG_EDGE_OUT_C) {
            1
        } else {
            0
        };
        degree += if self.edges.contains(Edges::FLAG_EDGE_OUT_G) {
            1
        } else {
            0
        };
        degree += if self.edges.contains(Edges::FLAG_EDGE_OUT_T) {
            1
        } else {
            0
        };

        degree
    }

    /// Identifies junctions in the graph
    #[must_use]
    pub fn is_junction(&self) -> bool {
        self.in_degree() > 1 || self.out_degree() > 1
    }
}

impl fmt::Display for Record {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} {}{}{}{}{}{}{}{}",
            self.coverage,
            if self.edges.contains(Edges::FLAG_EDGE_IN_A) {
                'a'
            } else {
                '.'
            },
            if self.edges.contains(Edges::FLAG_EDGE_IN_C) {
                'c'
            } else {
                '.'
            },
            if self.edges.contains(Edges::FLAG_EDGE_IN_G) {
                'g'
            } else {
                '.'
            },
            if self.edges.contains(Edges::FLAG_EDGE_IN_T) {
                't'
            } else {
                '.'
            },
            if self.edges.contains(Edges::FLAG_EDGE_OUT_A) {
                'A'
            } else {
                '.'
            },
            if self.edges.contains(Edges::FLAG_EDGE_OUT_C) {
                'C'
            } else {
                '.'
            },
            if self.edges.contains(Edges::FLAG_EDGE_OUT_G) {
                'G'
            } else {
                '.'
            },
            if self.edges.contains(Edges::FLAG_EDGE_OUT_T) {
                'T'
            } else {
                '.'
            },
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // use proptest::prelude::*;

    #[test]
    fn test_get_coverage() {
        let r1 = Record::new(1, None);
        let r10 = Record::new(10, None);

        assert!(r1.coverage() == 1);
        assert!(r10.coverage() == 10);
    }

    #[test]
    fn test_set_coverage() {
        let mut r100 = Record::new(0, None);
        r100.set_coverage(100);

        let mut r1000 = Record::new(0, None);
        r1000.set_coverage(1000);

        assert!(r100.coverage() == 100);
        assert!(r1000.coverage() == 1000);
    }

    #[test]
    fn test_increment_coverage() {
        let mut r100 = Record::new(99, None);
        r100.increment_coverage();

        let mut r1000 = Record::new(999, None);
        r1000.increment_coverage();

        assert!(r100.coverage() == 100);
        assert!(r1000.coverage() == 1000);
    }

    #[test]
    fn test_increment_coverage_saturates() {
        let mut rmax = Record::new(u16::MAX, None);
        rmax.increment_coverage();

        assert!(rmax.coverage() == u16::MAX);
    }

    #[test]
    fn test_set_incoming_edge() {
        let mut r = Record::new(1, None);

        assert!(r.edges.is_empty());

        r.set_incoming_edge('A' as u8);
        assert!(r.edges.contains(Edges::FLAG_EDGE_IN_A));
        assert!(!r.edges.contains(Edges::FLAG_EDGE_IN_C));
        assert!(!r.edges.contains(Edges::FLAG_EDGE_IN_G));
        assert!(!r.edges.contains(Edges::FLAG_EDGE_IN_T));

        r.set_incoming_edge('C' as u8);
        assert!(r.edges.contains(Edges::FLAG_EDGE_IN_A));
        assert!(r.edges.contains(Edges::FLAG_EDGE_IN_C));
        assert!(!r.edges.contains(Edges::FLAG_EDGE_IN_G));
        assert!(!r.edges.contains(Edges::FLAG_EDGE_IN_T));

        r.set_incoming_edge('G' as u8);
        assert!(r.edges.contains(Edges::FLAG_EDGE_IN_A));
        assert!(r.edges.contains(Edges::FLAG_EDGE_IN_C));
        assert!(r.edges.contains(Edges::FLAG_EDGE_IN_G));
        assert!(!r.edges.contains(Edges::FLAG_EDGE_IN_T));

        r.set_incoming_edge('T' as u8);
        assert!(r.edges.contains(Edges::FLAG_EDGE_IN_A));
        assert!(r.edges.contains(Edges::FLAG_EDGE_IN_C));
        assert!(r.edges.contains(Edges::FLAG_EDGE_IN_G));
        assert!(r.edges.contains(Edges::FLAG_EDGE_IN_T));

        assert!(!r.edges.contains(Edges::FLAG_EDGE_OUT_A));
        assert!(!r.edges.contains(Edges::FLAG_EDGE_OUT_C));
        assert!(!r.edges.contains(Edges::FLAG_EDGE_OUT_G));
        assert!(!r.edges.contains(Edges::FLAG_EDGE_OUT_T));
    }

    #[test]
    fn test_set_outgoing_edge() {
        let mut r = Record::new(1, None);

        assert!(r.edges.is_empty());

        r.set_outgoing_edge('A' as u8);
        assert!(r.edges.contains(Edges::FLAG_EDGE_OUT_A));
        assert!(!r.edges.contains(Edges::FLAG_EDGE_OUT_C));
        assert!(!r.edges.contains(Edges::FLAG_EDGE_OUT_G));
        assert!(!r.edges.contains(Edges::FLAG_EDGE_OUT_T));

        r.set_outgoing_edge('C' as u8);
        assert!(r.edges.contains(Edges::FLAG_EDGE_OUT_A));
        assert!(r.edges.contains(Edges::FLAG_EDGE_OUT_C));
        assert!(!r.edges.contains(Edges::FLAG_EDGE_OUT_G));
        assert!(!r.edges.contains(Edges::FLAG_EDGE_OUT_T));

        r.set_outgoing_edge('G' as u8);
        assert!(r.edges.contains(Edges::FLAG_EDGE_OUT_A));
        assert!(r.edges.contains(Edges::FLAG_EDGE_OUT_C));
        assert!(r.edges.contains(Edges::FLAG_EDGE_OUT_G));
        assert!(!r.edges.contains(Edges::FLAG_EDGE_OUT_T));

        r.set_outgoing_edge('T' as u8);
        assert!(r.edges.contains(Edges::FLAG_EDGE_OUT_A));
        assert!(r.edges.contains(Edges::FLAG_EDGE_OUT_C));
        assert!(r.edges.contains(Edges::FLAG_EDGE_OUT_G));
        assert!(r.edges.contains(Edges::FLAG_EDGE_OUT_T));

        assert!(!r.edges.contains(Edges::FLAG_EDGE_IN_A));
        assert!(!r.edges.contains(Edges::FLAG_EDGE_IN_C));
        assert!(!r.edges.contains(Edges::FLAG_EDGE_IN_G));
        assert!(!r.edges.contains(Edges::FLAG_EDGE_IN_T));
    }

    #[test]
    fn test_in_degree() {
        let mut r = Record::new(1, None);

        assert!(r.in_degree() == 0);

        r.set_incoming_edge('A' as u8);
        assert!(r.in_degree() == 1);

        r.set_incoming_edge('C' as u8);
        assert!(r.in_degree() == 2);

        r.set_incoming_edge('G' as u8);
        assert!(r.in_degree() == 3);

        r.set_incoming_edge('T' as u8);
        assert!(r.in_degree() == 4);

        r.set_outgoing_edge('A' as u8);
        assert!(r.in_degree() == 4);
    }

    #[test]
    fn test_out_degree() {
        let mut r = Record::new(1, None);

        assert!(r.out_degree() == 0);

        r.set_outgoing_edge('A' as u8);
        assert!(r.out_degree() == 1);

        r.set_outgoing_edge('C' as u8);
        assert!(r.out_degree() == 2);

        r.set_outgoing_edge('G' as u8);
        assert!(r.out_degree() == 3);

        r.set_outgoing_edge('T' as u8);
        assert!(r.out_degree() == 4);

        r.set_outgoing_edge('A' as u8);
        assert!(r.out_degree() == 4);
    }

    #[test]
    fn test_incoming_junction() {
        let mut r = Record::new(1, None);

        assert!(r.is_junction() == false);

        r.set_incoming_edge('A' as u8);
        assert!(r.is_junction() == false);

        r.set_incoming_edge('C' as u8);
        assert!(r.is_junction() == true);
    }

    #[test]
    fn test_outgoing_junction() {
        let mut r = Record::new(1, None);

        assert!(r.is_junction() == false);

        r.set_outgoing_edge('A' as u8);
        assert!(r.is_junction() == false);

        r.set_outgoing_edge('C' as u8);
        assert!(r.is_junction() == true);
    }
}
