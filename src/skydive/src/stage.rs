// Import the Result type from the anyhow crate for error handling.
use anyhow::Result;
use rust_htslib::bam::ext::BamRecordExtensions;

// Import various standard library collections.
use std::collections::{HashMap, HashSet};
use std::env;
use std::path::PathBuf;
use std::io::{self, Write};
use linear_map::LinearMap;

// Import the Url type to work with URLs.
use url::Url;

// Import ExponentialBackoff for retrying operations.
use backoff::ExponentialBackoff;

// Import Gag to suppress output to stderr.
use gag::Gag;

// Import rayon's parallel iterator traits.
use rayon::prelude::*;
use rayon::iter::{ParallelIterator, IntoParallelRefIterator};

// Import types from rust_htslib for working with BAM files.
use rust_htslib::bam::{ self, Header, IndexedReader, Read };
use rust_htslib::bam::header::HeaderRecord;
use bio_types::genome::Interval;

// Import functions for authorizing access to Google Cloud Storage.
use crate::env::{ gcs_authorize_data_access, local_guess_curl_ca_bundle };

// Function to open a BAM file from a URL and cache its contents locally.
fn open_bam(reads_url: &Url, cache_path: &PathBuf) -> Result<IndexedReader> {
    env::set_current_dir(cache_path).unwrap();

    // Try to open the BAM file from the URL, with retries for authorization.
    let bam = match IndexedReader::from_url(reads_url) {
        Ok(bam) => bam,
        Err(_) => {
            // If opening fails, try authorizing access to Google Cloud Storage.
            gcs_authorize_data_access();

            // Try opening the BAM file again.
            match IndexedReader::from_url(reads_url) {
                Ok(bam) => bam,
                Err(_) => {
                    // If it still fails, guess the cURL CA bundle path.
                    local_guess_curl_ca_bundle();

                    // Try one last time to open the BAM file.
                    IndexedReader::from_url(reads_url)?
                }
            }
        }
    };

    Ok(bam)
}

// Function to extract reads from a BAM file within a specified genomic region.
fn extract_reads(bam: &mut IndexedReader, chr: &String, start: &u64, stop: &u64) -> Result<Vec<bam::Record>> {
    let mut records = Vec::new();

    // Fetch the genomic region from the BAM file.
    let _ = bam.fetch(((*chr).as_bytes(), *start, *stop));
    for (_, r) in bam.records().enumerate() {
        let record = r?;

        records.push(record);
    }

    Ok(records)
}

// Function to stage data from a single BAM file.
fn stage_data_from_one_file(
    reads_url: &Url,
    loci: &HashSet<(String, u64, u64)>,
    cache_path: &PathBuf,
) -> Result<(Vec<LinearMap<String, String>>, Vec<bam::Record>)> {
    let mut bam = open_bam(reads_url, cache_path)?;
    let header_view = bam.header().to_owned();

    let header = Header::from_template(&header_view);
    let header_map = header.to_hashmap();
    let read_groups = header_map.get("RG").unwrap().clone();

    let mut all_reads = Vec::new();
    for (chr, start, stop) in loci.iter() {
        // Extract reads for the current locus.
        let reads = extract_reads(&mut bam, chr, start, stop).unwrap();

        // Extend the all_reads vector with the reads from the current locus.
        all_reads.extend(reads);
    }

    Ok((read_groups, all_reads))
}

// Function to stage data from multiple BAM files.
fn stage_data_from_all_files(
    reads_urls: &HashSet<Url>,
    loci: &HashSet<(String, u64, u64)>,
    cache_path: &PathBuf,
) -> Result<Vec<(Vec<LinearMap<std::string::String, std::string::String>>, Vec<rust_htslib::bam::Record>)>> {

    // Use a parallel iterator to process multiple BAM files concurrently.
    let all_data: Vec<_> = reads_urls
        .par_iter()
        .map(|reads_url| {
            // Define an operation to stage data from one file.
            let op = || {
                let (read_groups, reads) = stage_data_from_one_file(reads_url, loci, cache_path)?;
                Ok((read_groups, reads))
            };

            // Retry the operation with exponential backoff in case of failure.
            match backoff::retry(ExponentialBackoff::default(), op) {
                Ok((read_groups, reads)) => { (read_groups, reads) }
                Err(e) => {
                    // If all retries fail, panic with an error message.
                    panic!("Error: {}", e);
                }
            }
        })
        .collect();

    // Return a vector of vectors, each containing read groups and reads from one BAM file.
    Ok(all_data)
}

// Function to retrieve the header of a BAM file.
fn get_bam_header(
    reads_url: &Url,
    cache_path: &PathBuf,
) -> Result<bam::HeaderView> {
    // Open the BAM file.
    let bam = open_bam(reads_url, cache_path)?;

    // Return the header of the BAM file.
    Ok(bam.header().to_owned())
}

pub fn read_spans_locus(
    start: i64,
    end: i64,
    loci: &HashSet<(String, u64, u64)>,
) -> bool {
    loci.iter().any(|e| start <= e.1 as i64 && end >= e.2 as i64)
}

// Public function to stage data from multiple BAM files and write to an output file.
pub fn stage_data(
    output_path: &PathBuf,
    loci: &HashSet<(String, u64, u64)>,
    reads_urls: &HashSet<Url>,
    cache_path: &PathBuf,
    require_spanning_reads: bool,
) -> Result<()> {
    // Disable stderr from trying to open an IndexedReader a few times, so
    // that the Jupyter notebook user doesn't get confused by intermediate
    // error messages that are nothing to worry about. The gag will end
    // automatically when it goes out of scope at the end of the function.
    let mut _stderr_gag = Gag::stderr().unwrap();

    // Get the header from the first BAM file.
    let first_url = reads_urls.iter().next().unwrap();
    let first_header = match get_bam_header(first_url, cache_path) {
        Ok(first_header) => first_header,
        Err(e) => {
            drop(_stderr_gag);

            panic!("Error: {}", e);
        }
    };

    // Create a new header based on the first header.
    let mut header = bam::Header::from_template(&first_header);
    let header_map = header.to_hashmap();
    let rgs = header_map.get("RG").unwrap();
    let mut rg_ids: HashSet<String> = rgs
        .iter()
        .filter_map(|a| a.get("ID").cloned())
        .collect();

    // Stage data from all BAM files.
    let all_data = match stage_data_from_all_files(reads_urls, loci, cache_path) {
        Ok(all_data) => all_data,
        Err(e) => {
            drop(_stderr_gag);

            panic!("Error: {}", e);
        }
    };

    // Populate the output header with read group information, avoiding the read group that's already in the header.
    all_data
        .iter()
        .for_each(|(h, _)| {
            let mut hr = HeaderRecord::new("RG".as_bytes());

            h
                .iter()
                .filter(|a| {
                    let status = !rg_ids.contains(a.get("ID").unwrap());
                    rg_ids.insert(a.get("ID").unwrap().clone());

                    status
                })
                .flatten()
                .for_each(|(k, v)| {
                    hr.push_tag(k.as_bytes(), v);
                });

            header.push_record(&hr);
        });

    // Create a BAM writer to write the reads to the output file.
    let mut bam_writer = bam::Writer::from_path(output_path, &header, bam::Format::Bam)?;
    let mut seen_read_ids = HashSet::new();
    all_data
        .iter()
        .for_each(|(_, reads)| {
            reads
                .iter()
                .for_each(|read| {
                    if !seen_read_ids.contains(read.qname()) && read_spans_locus(read.reference_start(), read.reference_end(), loci) {
                        let _ = bam_writer.write(read);
                    }

                    seen_read_ids.insert(read.qname());
                });
        });

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use url::Url;
    use std::collections::HashSet;

    // This test may pass, but still print a message to stderr regarding its failure to access data. This is because
    // open_bam() tries a couple of authorization methods before accessing data, and the initial failures print a
    // message to stderr. Elsewhere in the code, we suppress such messages (i.e. in stage_data()), but here we don't.
    #[test]
    fn test_open_bam() {
        let reads_url = Url::parse(
            "gs://fc-8c3900db-633f-477f-96b3-fb31ae265c44/results/PBFlowcell/m84060_230907_210011_s2/reads/ccs/aligned/m84060_230907_210011_s2.bam"
        ).unwrap();
        let cache_path = std::env::temp_dir();

        let bam = open_bam(&reads_url, &cache_path);

        assert!(bam.is_ok(), "Failed to open bam file");
    }

    #[test]
    fn test_stage_data_from_one_file() {
        let reads_url = Url::parse(
            "gs://fc-8c3900db-633f-477f-96b3-fb31ae265c44/results/PBFlowcell/m84060_230907_210011_s2/reads/ccs/aligned/m84060_230907_210011_s2.bam"
        ).unwrap();
        let loci = HashSet::from([("chr15".to_string(), 23960193, 23963918)]);
        let cache_path = std::env::temp_dir();

        let result = stage_data_from_one_file(&reads_url, &loci, &cache_path);

        assert!(result.is_ok(), "Failed to stage data from one file");
    }

    #[test]
    fn test_stage_data() {
        let cache_path = std::env::temp_dir();
        let output_path = cache_path.join("test.bam");

        let loci = HashSet::from([("chr15".to_string(), 23960193, 23963918)]);
        let reads_url = Url::parse(
            "gs://fc-8c3900db-633f-477f-96b3-fb31ae265c44/results/PBFlowcell/m84060_230907_210011_s2/reads/ccs/aligned/m84060_230907_210011_s2.bam"
        ).unwrap();
        let reads_urls = HashSet::from([reads_url]);

        let result = stage_data(&output_path, &loci, &reads_urls, &cache_path, false);

        assert!(result.is_ok(), "Failed to stage data from file");

        println!("{:?}", result.unwrap());
    }

    #[test]
    fn test_stage_multiple_data() {
        let cache_path = std::env::temp_dir();
        let output_path = cache_path.join("test.bam");

        let reads_url_1 = Url::parse(
            "gs://fc-8c3900db-633f-477f-96b3-fb31ae265c44/results/PBFlowcell/m84060_230907_210011_s2/reads/ccs/aligned/m84060_230907_210011_s2.bam"
        ).unwrap();
        let reads_url_2 = Url::parse(
            "gs://fc-8c3900db-633f-477f-96b3-fb31ae265c44/results/PBFlowcell/m84043_230901_211947_s1/reads/ccs/aligned/m84043_230901_211947_s1.bam"
        ).unwrap();
        let loci = HashSet::from([("chr15".to_string(), 23960193, 23963918)]);
        let reads_urls = HashSet::from([ reads_url_1, reads_url_2 ]);

        let result = stage_data(&output_path, &loci, &reads_urls, &cache_path, false);

        println!("{:?}", result);

        assert!(result.is_ok(), "Failed to stage data from all files");
    }
}
