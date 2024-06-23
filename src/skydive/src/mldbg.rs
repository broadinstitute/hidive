use crate::ldbg::LdBG;

/// Represents a multi-color linked de Bruijn graph, all built with same k-mer size.
#[derive(Debug)]
pub struct MLdBG {
    pub ldbgs: Vec<LdBG>
}

impl MLdBG {
    /// Create an empty multi-color LdBG.
    pub fn new() -> Self {
        MLdBG {
            ldbgs: Vec::new()
        }
    }

    /// Add a LdBG to the MLdBG.
    pub fn push(&mut self, ldbg: LdBG) {
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