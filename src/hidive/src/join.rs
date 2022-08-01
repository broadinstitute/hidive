use std::path::PathBuf;

use skydive;

pub fn start(output: PathBuf, graph: Vec<PathBuf>) {
    println!("The answer is {} {:?}!", output.display(), graph)
}