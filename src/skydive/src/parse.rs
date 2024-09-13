use std::{collections::HashSet, io::BufRead, path::PathBuf};

use anyhow::Result;
use url::Url;

use path_absolutize::Absolutize;

/// Parse a list of loci into a `HashSet` of formatted loci.
///
/// # Arguments
///
/// * `loci_list` - A reference to a vector of strings representing loci.
///
/// # Returns
///
/// A `HashSet` containing tuples of chromosome, start, and stop positions.
///
/// # Panics
///
/// If a locus cannot be parsed.
#[must_use]
pub fn parse_loci(loci_list: &Vec<String>) -> HashSet<(String, u64, u64)> {
    // Initialize a HashSet to store unique loci after parsing
    let mut loci = HashSet::new();

    // Iterate over each locus in the provided list
    for locus in loci_list {
        // Attempt to parse the locus using a function from the skydive module
        match parse_locus(&locus.to_owned()) {
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

    loci
}

/// Parse a locus string into a tuple of chromosome, start, and stop positions.
/// The locus should be in the format `chr:start[-stop]`.
///
/// # Arguments
///
/// * `locus` - A string representing the locus.
///
/// # Returns
///
/// A tuple containing the chromosome, start, and stop positions.
///
/// # Errors
///
/// This function returns an error if the locus format is incorrect.
///
/// # Panics
///
/// This function panics if the start or stop positions cannot be parsed.
pub fn parse_locus(locus: &str) -> Result<(String, u64, u64)> {
    let l_fmt = locus.replace(',', "");
    let parts1: Vec<&str> = l_fmt.split(|c| c == ':').collect();
    let parts2: Vec<&str> = parts1[1].split(|c| c == '-').collect();

    let chr = parts1[0].to_string();

    if parts2.len() == 1 {
        let start = match parts2[0].parse::<u64>() {
            Ok(val) => val,
            Err(e) => {
                return Err(anyhow::Error::new(e));
            }
        };

        Ok((chr, start - 1000, start + 1000))
    } else if parts2.len() == 2 {
        let start = match parts2[0].parse::<u64>() {
            Ok(val) => val,
            Err(e) => {
                return Err(anyhow::Error::new(e));
            }
        };

        let stop = match parts2[1].parse::<u64>() {
            Ok(val) => val,
            Err(e) => {
                return Err(anyhow::Error::new(e));
            }
        };

        Ok((chr, start, stop))
    } else {
        anyhow::bail!(
            "Locus format for '{}' is incorrect. It should be 'chr:start[-stop]'.",
            locus
        );
    }
}

/// Parse a list of file paths into a `HashSet` of URLs.
///
/// # Arguments
///
/// * `bam_paths` - A reference to a vector of `PathBufs` representing file paths.
///
/// # Returns
///
/// A `HashSet` containing URLs for each file path.
///
/// # Panics
///
/// If a file path cannot be parsed as a URL.
pub fn parse_file_names(bam_paths: &[PathBuf]) -> HashSet<Url> {
    // Convert the list of BAM file paths into a HashSet of URLs
    let mut reads_urls: HashSet<Url> = bam_paths
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

    // If any of the files are a local file ending in .txt, assume it's a file of filenames.
    let mut local_file_contents = HashSet::new();
    let mut to_remove = HashSet::new();
    for url in &reads_urls {
        if url.scheme() == "file" {
            let path = url.to_file_path().unwrap();
            if path.extension().and_then(std::ffi::OsStr::to_str) == Some("txt") {
                if let Ok(file) = std::fs::File::open(&path) {
                    let reader = std::io::BufReader::new(file);
                    for line in reader.lines().map_while(Result::ok) {
                        let abs_path = PathBuf::from(line);
                        local_file_contents.insert(abs_path);
                    }
                }

                to_remove.insert(url.clone());
            }
        }
    }

    // Remove FOFN files from the set of BAM/CRAM files.
    to_remove.iter().for_each(|url| {
        let _ = reads_urls.remove(url);
    });

    // Add the files from the file of filenames to the full list of files.
    reads_urls.extend(local_file_contents.into_iter().filter_map(|path| {
        let path_str = path.to_string_lossy();
        if path_str.starts_with("gs://") {
            Url::parse(&path_str).ok()
        } else {
            Url::from_file_path(path.absolutize().unwrap()).ok()
        }
    }));

    reads_urls
}
