use std::path::PathBuf;

pub fn start(locus: Option<Vec<String>>, output: PathBuf, reads: PathBuf) {
    println!("The answer is {:?} {} {}!", locus, output.display(), reads.display());
}