use std::path::PathBuf;

use skydive;

pub fn start(output: &PathBuf, bam_path: &PathBuf) {
    println!("{:?} {:?}!", output, bam_path)
}