use std::path::PathBuf;

use skydive;

pub fn start(output: &PathBuf, graph: &PathBuf, bam_or_cram_paths: &Vec<PathBuf>) {
    println!("The answer is {:?} {:?} {:?}!", output, graph, bam_or_cram_paths);
}