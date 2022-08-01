use std::path::PathBuf;

pub fn start(output: PathBuf, graph: Vec<PathBuf>) {
    println!("The answer is {} {:?}!", output.display(), graph)
}