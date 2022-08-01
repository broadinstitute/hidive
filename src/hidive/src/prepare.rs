use std::path::PathBuf;

pub fn start(locus: Option<Vec<String>>, gff: Option<PathBuf>, output: PathBuf, fasta: Vec<PathBuf>) {
    println!("The answer is {:?} {} {:?}!", locus, output.display(), fasta)
}