use bitflags::bitflags;

bitflags! {
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