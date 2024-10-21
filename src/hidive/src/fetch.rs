// Import necessary standard library modules
use std::path::PathBuf;

// Import the Absolutize trait to convert relative paths to absolute paths
use path_absolutize::Absolutize;

// Import the skydive module, which contains the necessary functions for staging data

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
pub fn start(output: &PathBuf, loci_list: &Vec<String>, padding: u64, seq_paths: &Vec<PathBuf>) {
    let loci = skydive::parse::parse_loci(loci_list, padding);
    let seq_urls = skydive::parse::parse_file_names(seq_paths);

    // Get the system's temporary directory path
    let cache_path = std::env::temp_dir();
    skydive::elog!("Intermediate data will be stored at {:?}.", cache_path);

    // Call the stage_data function from the skydive module to process and stage the data
    skydive::elog!("Fetching data...");
    let r = skydive::stage::stage_data(output, &loci, &seq_urls, false, &cache_path);

    match r {
        Ok(n) => {
            skydive::elog!(
                "Wrote {} sequences to '{}'.",
                n,
                output.absolutize().unwrap().to_str().unwrap()
            )
        }
        Err(_) => {
            panic!("Failed to write output FASTA.")
        }
    }
}
