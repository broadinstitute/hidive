use std::fmt;
use std::collections::VecDeque;

use needletail::sequence::complement;

use parquet::data_type::AsBytes;

/// Represents metadata on a link (a series of junction choices in de Bruijn graph).
#[derive(Debug, Hash, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Link {
    pub is_forward: bool,
    pub junctions: VecDeque<u8>
}

impl Link {
    /// Create an empty link record.
    pub fn new(is_forward: bool) -> Self {
        Link {
            is_forward,
            junctions: VecDeque::new()
        }
    }

    /// Create a link from a sequence of junction choices.
    pub fn from_junctions(is_forward: bool, seq: &[u8]) -> Self {
        Link {
            is_forward,
            junctions: VecDeque::from(seq.to_vec())
        }
    }

    /// Return orientation of the link.
    pub fn is_forward(&self) -> bool {
        self.is_forward
    }

    /// Return the number of junction choices in the link.
    pub fn len(&self) -> usize {
        self.junctions.len()
    }

    /// Indicate whether the list of junctions is empty.
    pub fn is_empty(&self) -> bool {
        self.junctions.is_empty()
    }

    /// Add a junction to the queue.
    pub fn push_back(&mut self, junction: u8) {
        self.junctions.push_back(junction);
    }

    /// Peek at a junction from the front of the queue.
    pub fn front(&self) -> Option<&u8> {
        self.junctions.front()
    }

    /// Take a junction from the front of the queue.
    pub fn pop_front(&mut self) -> Option<u8> {
        self.junctions.pop_front()
    }

    /// Return a new link with the junction choices complemented.
    pub fn complement(&self) -> Link {
        let mut new_link = Link::new(self.is_forward);

        for junction in &self.junctions {
            new_link.push_back(complement(*junction));
        }

        new_link
    }
}

impl fmt::Display for Link {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut jvec = Vec::new();
        for junction in &self.junctions {
            jvec.push(*junction);
        }
        
        write!(f, "{} {:?}",
            if self.is_forward { "F" } else { "R" },
            std::str::from_utf8(jvec.as_bytes())
        )
    }
}