// Import necessary standard library modules
use std::path::PathBuf;

// Import the Absolutize trait to convert relative paths to absolute paths
use path_absolutize::Absolutize;

// Import the skydive module, which contains the necessary functions for staging data
use skydive;

/// Starts the data fetching process.
///
/// # Arguments
///
/// * `output` - A reference to a PathBuf that specifies the output file path.
/// * `loci_list` - A reference to a vector of strings, each representing a genomic locus.
/// * `seq_paths` - A reference to a vector of PathBufs, each representing a path to a BAM or FASTA file.
///
/// # Panics
///
/// Panics if any locus in `loci_list` cannot be parsed.
pub fn start(output: &PathBuf, loci_list: &Vec<String>, seq_paths: &Vec<PathBuf>) {
    let loci = skydive::utils::parse_loci(loci_list);
    let seq_urls = skydive::utils::parse_file_names(seq_paths);

    // Get the system's temporary directory path
    let cache_path = std::env::temp_dir();
    eprintln!("[{}] Data will be cached to {:?}.", chrono::Local::now().format("%Y-%m-%d %H:%M:%S"), cache_path);

    // Convert the output path to an absolute path and own it
    let output_path = output.absolutize().unwrap().into_owned();

    // Call the stage_data function from the skydive module to process and stage the data
    eprintln!("[{}] Fetching data...", chrono::Local::now().format("%Y-%m-%d %H:%M:%S"));
    let r = skydive::stage::stage_data(&output_path, &loci, &seq_urls, &cache_path);

    match r {
        Ok(_) => {}
        Err(_) => {
            panic!("Failed to write multi-sample locus BAM.")
        }
    }
}
