pub mod ldbg;
pub mod link;
pub mod edges;
pub mod record;

pub mod mldbg;

pub mod env;
pub mod parse;
pub mod utils;

pub mod stage;
pub mod storage_gcs;
pub mod storage_local;
pub mod agg;

#[macro_export]
macro_rules! elog {
    ($($arg:tt)*) => {
        eprintln!("[{}] {}", chrono::Local::now().format("%Y-%m-%d %H:%M:%S"), format_args!($($arg)*));
    };
}