use std::fmt;

/// Represents metadata on a link (a series of junction choices in de Bruijn graph).
#[derive(Debug)]
pub struct LinkData {
    coverage: u16,
    is_forward: bool
}

impl LinkData {
    /// Create an empty de Bruijn graph record.
    pub fn new(coverage: u16, is_forward: bool) -> Self {
        LinkData {
            coverage,
            is_forward
        }
    }

    /// Return orientation of the link.
    pub fn is_forward(&self) -> bool {
        self.is_forward
    }

    /// Return the LinkData's coverage.
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

impl fmt::Display for LinkData {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}",
            if self.is_forward { "F" } else { "R" },
            self.coverage,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_coverage() {
        let r1 = LinkData::new(1, true);
        let r10 = LinkData::new(10, true);

        assert!(r1.coverage() == 1);
        assert!(r10.coverage() == 10);
    }

    #[test]
    fn test_set_coverage() {
        let mut r100 = LinkData::new(0, true);
        r100.set_coverage(100);

        let mut r1000 = LinkData::new(0, true);
        r1000.set_coverage(1000);

        assert!(r100.coverage() == 100);
        assert!(r1000.coverage() == 1000);
    }

    #[test]
    fn test_increment_coverage() {
        let mut r100 = LinkData::new(99, true);
        r100.increment_coverage();

        let mut r1000 = LinkData::new(999, true);
        r1000.increment_coverage();

        assert!(r100.coverage() == 100);
        assert!(r1000.coverage() == 1000);
    }

    #[test]
    fn test_increment_coverage_saturates() {
        let mut rmax = LinkData::new(u16::MAX, true);
        rmax.increment_coverage();

        assert!(rmax.coverage() == u16::MAX);
    }
}