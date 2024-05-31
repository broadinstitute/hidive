use std::{collections::HashSet, io::BufRead, path::PathBuf};

use anyhow::Result;
use url::Url;

use path_absolutize::Absolutize;

pub fn parse_loci(loci_list: &Vec<String>) -> HashSet<(String, u64, u64)> {
    // Initialize a HashSet to store unique loci after parsing
    let mut loci = HashSet::new();
    
    // Iterate over each locus in the provided list
    for locus in loci_list {
        // Attempt to parse the locus using a function from the skydive module
        match parse_locus(locus.to_owned()) {
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

pub fn parse_locus(locus: String) -> Result<(String, u64, u64)> {
    let l_fmt = locus.replace(",", "");
    let parts: Vec<&str> = l_fmt.split(|c| (c == ':' || c == '-')).collect();

    let chr = parts[0].to_string();

    if parts.len() == 2 {
        let start = match parts[1].parse::<u64>() {
            Ok(val) => val,
            Err(e) => {
                return Err(anyhow::Error::new(e));
            }
        };

        Ok((chr, start - 1000, start + 1000))
    } else if parts.len() == 3 {
        let start = match parts[1].parse::<u64>() {
            Ok(val) => val,
            Err(e) => {
                return Err(anyhow::Error::new(e));
            }
        };

        let stop = match parts[2].parse::<u64>() {
            Ok(val) => val,
            Err(e) => {
                return Err(anyhow::Error::new(e));
            }
        };

        Ok((chr, start, stop))
    } else {
        anyhow::bail!("Locus format for '{}' is incorrect. It should be 'chr:start[-stop]'.", locus);
    }
}

pub fn parse_file_names(bam_paths: &Vec<PathBuf>) -> HashSet<Url> {
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
                    for line in reader.lines() {
                        if let Ok(line_path) = line {
                            let abs_path = PathBuf::from(line_path);
                            local_file_contents.insert(abs_path);
                        }
                    }
                }

                to_remove.insert(url.clone());
            }
        }
    }

    // Remove FOFN files from the set of BAM/CRAM files.
    to_remove.iter().for_each(|url| { let _ = reads_urls.remove(url); });

    // Add the files from the file of filenames to the full list of files.
    reads_urls.extend(local_file_contents
        .into_iter()
        .filter_map(|path| {
            let path_str = path.to_string_lossy();
            if path_str.starts_with("gs://") {
                Url::parse(&path_str).ok()
            } else {
                Url::from_file_path(path.absolutize().unwrap()).ok()
            }
        })
    );

    reads_urls
}
