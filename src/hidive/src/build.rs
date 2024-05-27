// Import necessary standard library modules
use std::path::PathBuf;
use std::collections::HashSet;

// Import the Absolutize trait to convert relative paths to absolute paths
use path_absolutize::Absolutize;

// Import the Url type to work with URLs
use url::Url;

// Import the skydive module, which contains the necessary functions for building graphs
use skydive;

pub fn start(output: &PathBuf, bam_path: &PathBuf) {
    println!("{:?} {:?}!", output, bam_path)
}