// Import necessary standard library modules
use std::path::PathBuf;
use std::collections::HashSet;

// Import the Absolutize trait to convert relative paths to absolute paths
use path_absolutize::Absolutize;

// Import the Url type to work with URLs
use url::Url;

// Import the skydive module, which contains the necessary functions for staging data
use skydive;

/// Starts the data fetching process.
///
/// # Arguments
///
/// * `output` - A reference to a PathBuf that specifies the output file path.
/// * `loci_list` - A reference to a vector of strings, each representing a genomic locus.
/// * `bam_paths` - A reference to a vector of PathBufs, each representing a path to a BAM file.
/// * `require_spanning_reads` - Flag indicating that reads should fully span the requested locus.
///
/// # Panics
///
/// Panics if any locus in `loci_list` cannot be parsed.
pub fn start(output: &PathBuf, loci_list: &Vec<String>, bam_paths: &Vec<PathBuf>, require_spanning_reads: bool) {
    let loci = skydive::utils::parse_loci(loci_list);

    // Convert the list of BAM file paths into a HashSet of URLs
    let reads_urls: HashSet<Url> = bam_paths
        .iter()
        // Use filter_map to attempt to parse each path as a URL, and collect the successful ones
        .filter_map(|path| {
            let path_str = path.to_string_lossy();
            if path_str.starts_with("gs://") {
                Url::parse(&path_str).ok()
            } else {
                Url::from_file_path(path.absolutize().unwrap()).ok()
            }
        })
        .collect();

    // Get the system's temporary directory path
    let cache_path = std::env::temp_dir();

    // Convert the output path to an absolute path and own it
    let output_path = output.absolutize().unwrap().into_owned();

    // Call the stage_data function from the skydive module to process and stage the data
    let r = skydive::stage::stage_data(&output_path, &loci, &reads_urls, &cache_path, require_spanning_reads);

    match r {
        Ok(_) => {},
        Err(_) => { panic!("Failed to write multi-sample locus BAM.") },
    }
}