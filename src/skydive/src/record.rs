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
}