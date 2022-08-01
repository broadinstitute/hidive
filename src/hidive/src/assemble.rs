use std::path::PathBuf;

use skydive;

pub fn start(output: PathBuf, locus: Option<Vec<String>>, reads: PathBuf) {
    println!("The answer is {:?} {} {}!", locus, output.display(), reads.display());
}