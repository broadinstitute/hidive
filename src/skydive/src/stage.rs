use anyhow::Result;
use std::collections::{ HashSet, HashMap };
use std::env;
use std::path::PathBuf;
use url::Url;

use backoff::ExponentialBackoff;
use gag::Gag;
use rayon::prelude::*;

use rust_htslib::bam::record::{ Aux, Cigar };
use rust_htslib::bam::{ self, Read, IndexedReader, ext::BamRecordExtensions };

use crate::env::{ gcs_authorize_data_access, local_guess_curl_ca_bundle };

fn open_bam(reads_url: &Url, cache_path: &PathBuf) -> Result<IndexedReader> {
    env::set_current_dir(cache_path).unwrap();

    let bam = match IndexedReader::from_url(reads_url) {
        Ok(bam) => bam,
        Err(_) => {
            gcs_authorize_data_access();

            match IndexedReader::from_url(reads_url) {
                Ok(bam) => bam,
                Err(_) => {
                    local_guess_curl_ca_bundle();

                    IndexedReader::from_url(reads_url)?
                }
            }
        }
    };

    Ok(bam)
}

fn extract_reads(bam: &mut IndexedReader, chr: &String, start: &u64, stop: &u64) -> Result<Vec<bam::Record>> {
    let mut records = Vec::new();
    let _ = bam.fetch(((*chr).as_bytes(), *start, *stop));
    for (_, r) in bam.records().enumerate() {
        let record = r?;

        records.push(record);
    }

    Ok(records)
}

fn stage_data_from_one_file(
    reads_url: &Url,
    loci: &HashSet<(String, u64, u64)>,
    cache_path: &PathBuf,
) -> Result<Vec<bam::Record>> {
    let mut bam = open_bam(reads_url, cache_path)?;

    let mut all_reads = Vec::new();
    for (chr, start, stop) in loci.iter() {
        let reads = extract_reads(&mut bam, chr, start, stop).unwrap();

        all_reads.extend(reads);
    }

    Ok(all_reads)
}

fn stage_data_from_all_files(
    reads_urls: &HashSet<Url>,
    loci: &HashSet<(String, u64, u64)>,
    cache_path: &PathBuf,
) -> Result<Vec<Vec<bam::Record>>> {
    let all_reads: Vec<_> = reads_urls
        .par_iter()
        .map(|reads_url| {
            let op = || {
                let reads = stage_data_from_one_file(reads_url, loci, cache_path)?;
                Ok(reads)
            };

            match backoff::retry(ExponentialBackoff::default(), op) {
                Ok(reads) => { reads }
                Err(e) => {
                    panic!("Error: {}", e);
                }
            }
        })
        .collect();

    Ok(all_reads)
}

fn get_bam_header(
    reads_url: &Url,
    cache_path: &PathBuf,
) -> Result<bam::HeaderView> {
    let bam = open_bam(reads_url, cache_path)?;

    Ok(bam.header().to_owned())
}

pub fn stage_data(
    output_path: &PathBuf,
    loci: &HashSet<(String, u64, u64)>,
    reads_urls: &HashSet<Url>,
    cache_path: &PathBuf,
) -> Result<()> {
    // Disable stderr from trying to open an IndexedReader a few times, so
    // that the Jupyter notebook user doesn't get confused by intermediate
    // error messages that are nothing to worry about. The gag will end
    // automatically when it goes out of scope at the end of the function.
    let _stderr_gag = Gag::stderr().unwrap();

    let first_url = reads_urls.iter().next().unwrap();
    let first_header = get_bam_header(first_url, cache_path)?;
    let header = bam::Header::from_template(&first_header);

    let all_reads = stage_data_from_all_files(reads_urls, loci, cache_path)?;

    let mut bam_writer = bam::Writer::from_path(output_path, &header, bam::Format::Bam)?;
    all_reads
        .iter()
        .flatten()
        .try_for_each(|read| bam_writer.write(read))?;

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

        let result = stage_data(&output_path, &loci, &reads_urls, &cache_path);

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

        let result = stage_data(&output_path, &loci, &reads_urls, &cache_path);

        println!("{:?}", result);

        assert!(result.is_ok(), "Failed to stage data from all files");
    }
}
