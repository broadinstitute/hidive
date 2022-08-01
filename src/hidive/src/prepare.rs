use std::path::PathBuf;

use skydive;

pub fn start(output: PathBuf, locus: Option<Vec<String>>, _gff: Option<PathBuf>, fasta: Vec<PathBuf>) {
    println!("The answer is {:?} {} {:?}!", locus, output.display(), fasta)
}