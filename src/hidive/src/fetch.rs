use std::path::PathBuf;
use std::collections::HashSet;
use path_absolutize::Absolutize;
use url::Url;

use skydive;

pub fn start(output: &PathBuf, loci_list: &Vec<String>, bam_paths: &Vec<PathBuf>) {
    let mut loci = HashSet::new();

    for locus in loci_list {
        match skydive::utils::parse_locus(locus.to_owned()) {
            Ok(l_fmt) => {
                loci.insert(l_fmt);
            }
            Err(_) => {
                panic!("Could not parse locus '{}'.", locus);
            }
        }
    }

    let reads_urls: HashSet<Url> = bam_paths
        .iter()
        .filter_map(
            |path| Url::parse(&path.to_string_lossy()).ok()
        )
        .collect();

    let cache_path = std::env::temp_dir();
    let output_path = output.absolutize().unwrap().into_owned();

    skydive::stage::stage_data(&output_path, &loci, &reads_urls, &cache_path);
}