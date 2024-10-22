use std::{collections::HashSet, io::BufRead, path::PathBuf};

use anyhow::Result;
use url::Url;
use regex::Regex;

use path_absolutize::Absolutize;

/// Parse a list of loci strings into a `HashSet` of formatted loci.
/// If a locus string is a file on disk, read its contents and parse each line as a locus.
///
/// # Arguments
///
/// * `loci_list` - A list of strings representing loci.
/// * `padding` - A u64 representing the amount to extend the start and stop positions by.
///
/// # Returns
///
/// A `HashSet` of formatted loci.
///
/// # Example
///
/// ```
/// let loci_list = vec![
///     "chr1:1000-2000".to_string(),
///     "chr2:3000-4000".to_string(),
///     "chr3:5000-6000".to_string(),
/// ];
/// let padding = 0;
/// let loci = parse_loci(&loci_list, padding);
/// ```
///
/// # Panics
///
/// This function will panic if it cannot parse a given locus.
/// It will also panic if it cannot read a given file path.
///
#[must_use]
pub fn parse_loci(loci_list: &Vec<String>, padding: u64) -> HashSet<(String, u64, u64, String)> {
    // Initialize a HashSet to store unique loci after parsing
    let mut loci = HashSet::new();

    // Iterate over each locus in the provided list
    for locus in loci_list {
        // Check if the locus represents a file on disk
        let path = PathBuf::from(locus);
        if path.is_file() {
            // If it's a file, read its contents and parse each line as a locus
            let file = std::fs::File::open(&path).expect("Failed to open file");
            let reader = std::io::BufReader::new(file);
            for line in reader.lines() {
                let line = line.expect("Failed to read line");

                // Skip lines starting with '#'
                if line.trim().starts_with('#') {
                    continue;
                }

                match parse_locus(line.to_owned(), padding) {
                    Ok(l_fmt) => {
                        loci.insert(l_fmt);
                    }
                    Err(_) => {
                        panic!("Could not parse locus '{}' from file '{}'.", line, locus);
                    }
                }
            }
            // Skip the rest of the loop iteration
            continue;
        } else {
            // Attempt to parse the locus
            match parse_locus(locus.to_owned(), padding) {
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
    }

    loci
}

/// Parse a locus string into a tuple of contig name, start position, stop position, and optional name.
/// The locus string can be in the following formats:
/// - chr:start-stop
/// - chr:start-stop|name
/// - chr start stop
/// - chr start stop name
///
/// The start and stop positions are 1-based and inclusive.
/// The optional padding parameter can be used to extend the start and stop positions by a specified amount.
///
/// # Arguments
///
/// * `locus` - A string representing a locus.
/// * `padding` - A u64 representing the amount to extend the start and stop positions by.
///
/// # Returns
///
/// A tuple containing the contig name, start position, stop position, and optional name.
///
/// # Example
/// ```
/// let locus = "chr1:1000-2000".to_string();
/// let padding = 0;
/// let result = parse_locus(locus, padding);
/// ```
///
/// # Panics
///
/// This function will panic if the locus format is incorrect.
pub fn parse_locus(locus: String, padding: u64) -> Result<(String, u64, u64, String)> {
    // Regex to capture the contig name, start position, stop position, and optional name.
    // Accepts:
    // - chr:start-stop
    // - chr:start-stop|name
    // - chr start stop
    // - chr start stop name
    let re = Regex::new(r"(.*)[:\s]+(\d+)[-\s]+(\d+)(?:[|\s+](.*))?")?;

    // Remove commas from the locus string
    let locus = locus.replace(",", "");

    if let Some(captures) = re.captures(&locus) {
        let chr = captures.get(1).unwrap().as_str().to_string();
        let start = captures.get(2).unwrap().as_str().parse::<u64>()? - padding;
        let stop = captures.get(3).unwrap().as_str().parse::<u64>()? + padding;
        let name = captures.get(4).map_or_else(
            || format!("{}:{}-{}", chr, start, stop),
            |m| m.as_str().to_string()
        );

        if start > stop {
            anyhow::bail!("Locus format for '{}' is incorrect. Start position ({}) is greater than stop position ({}).", locus, start, stop);
        }

        Ok((chr, start, stop, name))
    } else {
        anyhow::bail!(
            "Locus format for '{}' is incorrect. It should be 'chr:start-stop', 'chr:start-stop|name', 'chr start stop', or 'chr start stop name'.",
            locus
        );
    }
}

/// Parse a list of BAM file paths into a `HashSet` of URLs.
/// If any of the files are a local file ending in .txt, assume it's a file of filenames.
///
/// # Arguments
///
/// * `bam_paths` - A list of BAM file paths.
///
/// # Returns
///
/// A `HashSet` of URLs.
///
/// # Example
///
/// ```
/// let bam_paths = vec![
///    PathBuf::from("gs://bucket/file1.bam"),
///    PathBuf::from("gs://bucket/file2.bam"),
///    PathBuf::from("gs://bucket/file3.txt"),
/// ];
/// let urls = parse_file_names(&bam_paths);
/// ```
///
/// # Panics
///
/// This function will panic if it cannot parse a given file path.
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_locus() {
        // Valid locus without padding
        let result = parse_locus("chr1:1000-2000".to_string(), 0);
        assert_eq!(result.ok(), Some(("chr1".to_string(), 1000 as u64, 2000 as u64, "chr1:1000-2000".to_string())));

        // Valid locus with padding
        let result = parse_locus("chr2:5000-6000".to_string(), 100);
        assert!(result.is_ok());
        assert_eq!(result.ok(), Some(("chr2".to_string(), 4900 as u64, 6100 as u64, "chr2:4900-6100".to_string())));

        // Valid locus with name
        let result = parse_locus("chr3:10000-20000|gene1".to_string(), 0);
        assert_eq!(result.ok(), Some(("chr3".to_string(), 10000 as u64, 20000 as u64, "gene1".to_string())));

        // Valid locus with commas
        let result = parse_locus("chr3:10,000-20,000|gene1".to_string(), 0);
        assert_eq!(result.ok(), Some(("chr3".to_string(), 10000 as u64, 20000 as u64, "gene1".to_string())));

        // Combination of space and colon separators
        let result = parse_locus("chr4 30000-40000".to_string(), 0);
        assert_eq!(result.ok(), Some(("chr4".to_string(), 30000 as u64, 40000 as u64, "chr4:30000-40000".to_string())));

        // Invalid format (non-numeric start position)
        let result = parse_locus("chr5:start-50000".to_string(), 0);
        assert!(result.is_err());

        // Invalid format (start position greater than end position)
        let result = parse_locus("chr6:60000-50000".to_string(), 0);
        assert!(result.is_err());

        // Valid locus with tab-separated fields
        let result = parse_locus("chr7\t70000\t80000".to_string(), 0);
        assert_eq!(result.ok(), Some(("chr7".to_string(), 70000 as u64, 80000 as u64, "chr7:70000-80000".to_string())));

        // Valid locus with tab-separated fields and name
        let result = parse_locus("chr8\t90000\t100000\tgene2".to_string(), 0);
        assert_eq!(result.ok(), Some(("chr8".to_string(), 90000 as u64, 100000 as u64, "gene2".to_string())));

        // Valid locus with mixed tab and colon separators
        let result = parse_locus("chr9:110000\t120000".to_string(), 0);
        assert_eq!(result.ok(), Some(("chr9".to_string(), 110000 as u64, 120000 as u64, "chr9:110000-120000".to_string())));

        // Contig name with dash in it
        let result = parse_locus("chr10-A:130000-140000|chr10-A".to_string(), 0);
        assert_eq!(result.ok(), Some(("chr10-A".to_string(), 130000 as u64, 140000 as u64, "chr10-A".to_string())));
    }
}

