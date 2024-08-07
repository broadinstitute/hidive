//! Crate for building and manipulating a series-parallel graph representation
//! of long-read data and a linked de Bruijn graph representation of short-read data.
//!

pub mod edges;
pub mod ldbg;
pub mod link;
pub mod record;
pub mod agg;

pub mod mldbg;

pub mod env;
pub mod parse;
pub mod utils;

pub mod stage;
pub mod storage_gcs;
pub mod storage_local;

#[macro_export]
macro_rules! elog {
    ($($arg:tt)*) => {
        eprintln!("[{}] {}", chrono::Local::now().format("%Y-%m-%d %H:%M:%S"), format_args!($($arg)*));
    };
}
