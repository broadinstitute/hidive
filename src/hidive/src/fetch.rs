use std::path::PathBuf;

use skydive;

pub fn start(output: &PathBuf, locus: &Vec<String>, bams: &Vec<PathBuf>) {
    println!("{:?} {:?} {:?}!", output, locus, bams)
}