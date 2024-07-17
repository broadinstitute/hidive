use bitflags::bitflags;

bitflags! {
    // #[derive(Debug)]
    pub struct Edges: u8 {
        const FLAG_EDGE_IN_A  = 0b0000_0001;
        const FLAG_EDGE_IN_C  = 0b0000_0010;
        const FLAG_EDGE_IN_G  = 0b0000_0100;
        const FLAG_EDGE_IN_T  = 0b0000_1000;
        const FLAG_EDGE_OUT_A = 0b0001_0000;
        const FLAG_EDGE_OUT_C = 0b0010_0000;
        const FLAG_EDGE_OUT_G = 0b0100_0000;
        const FLAG_EDGE_OUT_T = 0b1000_0000;
    }
}

impl Edges {
    pub fn from_string(s: String) -> Self {
        let mut edges = Edges::empty();

        if s.chars().nth(0) == Some('a') { edges.insert(Self::FLAG_EDGE_IN_A); }
        if s.chars().nth(1) == Some('c') { edges.insert(Self::FLAG_EDGE_IN_C); }
        if s.chars().nth(2) == Some('g') { edges.insert(Self::FLAG_EDGE_IN_G); }
        if s.chars().nth(3) == Some('t') { edges.insert(Self::FLAG_EDGE_IN_T); }

        if s.chars().nth(4) == Some('A') { edges.insert(Self::FLAG_EDGE_OUT_A); }
        if s.chars().nth(5) == Some('C') { edges.insert(Self::FLAG_EDGE_OUT_C); }
        if s.chars().nth(6) == Some('G') { edges.insert(Self::FLAG_EDGE_OUT_G); }
        if s.chars().nth(7) == Some('T') { edges.insert(Self::FLAG_EDGE_OUT_T); }

        edges
    }
}