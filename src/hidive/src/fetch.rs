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
///
/// # Panics
///
/// Panics if any locus in `loci_list` cannot be parsed.
pub fn start(output: &PathBuf, loci_list: &Vec<String>, bam_paths: &Vec<PathBuf>) {
    // Initialize a HashSet to store unique loci after parsing
    let mut loci = HashSet::new();

    // Iterate over each locus in the provided list
    for locus in loci_list {
        // Attempt to parse the locus using a function from the skydive module
        match skydive::utils::parse_locus(locus.to_owned()) {
            Ok(l_fmt) => {
                // If parsing is successful, insert the formatted locus into the HashSet
                loci.insert(l_fmt);
            }
            Err(_) => {
                // If parsing fails, panic and terminate the program, providing an error message
                panic!("Could not parse locus '{}'.", locus);
            }
        }
    }

    // Convert the list of BAM file paths into a HashSet of URLs
    let reads_urls: HashSet<Url> = bam_paths
        .iter()
        // Use filter_map to attempt to parse each path as a URL, and collect the successful ones
        .filter_map(
            |path| Url::parse(&path.to_string_lossy()).ok()
        )
        .collect();

    // Get the system's temporary directory path
    let cache_path = std::env::temp_dir();

    // Convert the output path to an absolute path and own it
    let output_path = output.absolutize().unwrap().into_owned();

    // Call the stage_data function from the skydive module to process and stage the data
    let r = skydive::stage::stage_data(&output_path, &loci, &reads_urls, &cache_path);

    match r {
        Ok(_) => { eprintln!("Done!") },
        Err(_) => { panic!("Failed to write multi-sample locus BAM.") },
    }
}