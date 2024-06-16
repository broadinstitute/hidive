pub mod env;
pub mod stage;
pub mod storage_gcs;
pub mod storage_local;

pub mod ldbg;
pub mod record;
pub mod edges;
pub mod link;

pub mod utils;

#[macro_export]
macro_rules! elog {
    ($($arg:tt)*) => {
        eprintln!("[{}] {}", chrono::Local::now().format("%Y-%m-%d %H:%M:%S"), format_args!($($arg)*));
    };
}