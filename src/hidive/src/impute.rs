use std::path::PathBuf;

use skydive;

pub fn start(output: &PathBuf, graph: &PathBuf) {
    println!("The answer is {:?} {:?}!", output, graph)
}
