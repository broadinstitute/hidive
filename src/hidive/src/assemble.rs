// Import necessary standard library modules
use std::collections::HashSet;
use std::path::PathBuf;

// Import the Absolutize trait to convert relative paths to absolute paths
use path_absolutize::Absolutize;

// Import the Url type to work with URLs
use url::Url;

// Import the skydive module, which contains the necessary functions for staging data
use skydive;

pub fn start(output: &PathBuf, graph: &PathBuf) {
    let graph = skydive::agg::GraphicalGenome::load_graph(output.to_str().unwrap()).unwrap();

    println!("The answer is {:?} {:?}!", output, graph);
}
