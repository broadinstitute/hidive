// Import the Result type from the anyhow crate for error handling.
use anyhow::Result;
use rust_htslib::bam::record::Aux;

// Import various standard library collections.
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::path::PathBuf;
use std::io::BufWriter;

// Import the Url type to work with URLs.
use url::Url;

// Import ExponentialBackoff for retrying operations.
use backoff::ExponentialBackoff;

// Import rayon's parallel iterator traits.
use rayon::prelude::*;
use rayon::iter::{ParallelIterator, IntoParallelRefIterator};

// Import types from rust_htslib for working with BAM files.
use rust_htslib::bam::{ self, IndexedReader, Read };
use rust_htslib::faidx::Reader;
use bio::io::fasta;

// Import functions for authorizing access to Google Cloud Storage.
use crate::env::{ gcs_authorize_data_access, local_guess_curl_ca_bundle };
use crate::elog;

// Function to open a BAM file from a URL and cache its contents locally.
fn open_bam(seqs_url: &Url, cache_path: &PathBuf) -> Result<IndexedReader> {
    env::set_current_dir(cache_path).unwrap();

    if env::var("GCS_OAUTH_TOKEN").is_err() {
        gcs_authorize_data_access();
    }

    // Try to open the BAM file from the URL, with retries for authorization.
    let bam = match IndexedReader::from_url(seqs_url) {
        Ok(bam) => bam,
        Err(_) => {
<<<<<<< HEAD
            elog!("Read '{}', attempt 2 (reauthorizing to GCS)", seqs_url);
=======
            elog!("[{}] Read '{}', attempt 2 (reauthorizing to GCS)", chrono::Local::now().format("%Y-%m-%d %H:%M:%S"), seqs_url);
>>>>>>> bfdb82b (Simplifying build command)

            // If opening fails, try authorizing access to Google Cloud Storage.
            gcs_authorize_data_access();

            // Try opening the BAM file again.
            match IndexedReader::from_url(seqs_url) {
                Ok(bam) => bam,
                Err(_) => {
<<<<<<< HEAD
                    elog!("Read '{}', attempt 3 (overriding cURL CA bundle)", seqs_url);
=======
                    elog!("[{}] Read '{}', attempt 3 (overriding cURL CA bundle)", chrono::Local::now().format("%Y-%m-%d %H:%M:%S"), seqs_url);
>>>>>>> bfdb82b (Simplifying build command)

                    // If it still fails, guess the cURL CA bundle path.
                    local_guess_curl_ca_bundle();

                    // Try one last time to open the BAM file.
                    IndexedReader::from_url(seqs_url)?
                }
            }
        }
    };

    Ok(bam)
}

// Function to open a FASTA file from a URL and cache its contents locally.
fn open_fasta(seqs_url: &Url, cache_path: &PathBuf) -> Result<Reader> {
    env::set_current_dir(cache_path).unwrap();

    if env::var("GCS_OAUTH_TOKEN").is_err() {
        gcs_authorize_data_access();
    }

    // Try to open the FASTA file from the URL, with retries for authorization.
    let fasta = match Reader::from_url(seqs_url) {
        Ok(fasta) => fasta,
        Err(_) => {
<<<<<<< HEAD
            elog!("Read '{}', attempt 2 (reauthorizing to GCS)", seqs_url);
=======
            elog!("[{}] Read '{}', attempt 2 (reauthorizing to GCS)", chrono::Local::now().format("%Y-%m-%d %H:%M:%S"), seqs_url);
>>>>>>> bfdb82b (Simplifying build command)

            // If opening fails, try authorizing access to Google Cloud Storage.
            gcs_authorize_data_access();

            // Try opening the BAM file again.
            match Reader::from_url(seqs_url) {
                Ok(bam) => bam,
                Err(_) => {
<<<<<<< HEAD
                    elog!("Read '{}', attempt 3 (overriding cURL CA bundle)", seqs_url);
=======
                    elog!("[{}] Read '{}', attempt 3 (overriding cURL CA bundle)", chrono::Local::now().format("%Y-%m-%d %H:%M:%S"), seqs_url);
>>>>>>> bfdb82b (Simplifying build command)

                    // If it still fails, guess the cURL CA bundle path.
                    local_guess_curl_ca_bundle();

                    // Try one last time to open the FASTA file.
                    Reader::from_url(seqs_url)?
                }
            }
        }
    };

    Ok(fasta)
}

// Function to get a mapping between read group and sample name from a BAM header.
fn get_rg_to_sm_mapping(bam: &IndexedReader) -> HashMap<String, String> {
    let header = bam::Header::from_template(bam.header());

    let rg_sm_map: HashMap<String, String> = header
        .to_hashmap()
        .into_iter()
        .flat_map(|(_, records)| records)
        .filter(|record| record.contains_key("ID") && record.contains_key("SM"))
        .map(|record| (record["ID"].to_owned(), record["SM"].to_owned()))
        .collect();

    rg_sm_map
}

fn get_sm_name_from_rg(alignment: &bam::pileup::Alignment, rg_sm_map: &HashMap<String, String>) -> Result<String> {
    let r = alignment.record();
    let rg = r.aux(b"RG")?;

    if let Aux::String(v) = rg {
        if let Some(sm) = rg_sm_map.get(v) {
            Ok(sm.to_owned())
        } else {
            Err(anyhow::Error::msg(format!("Sample name not found for read group: {}", v)))
        }
    } else {
        Err(anyhow::Error::msg("Read group is not a string"))
    }
}

// Function to extract seqs from a BAM file within a specified genomic region.
fn extract_bam_reads(basename: &String, bam: &mut IndexedReader, chr: &String, start: &u64, stop: &u64) -> Result<Vec<fasta::Record>> {
    let rg_sm_map = get_rg_to_sm_mapping(&bam);

    let mut bmap = HashMap::new();

    let _ = bam.fetch(((*chr).as_bytes(), *start, *stop));
    for p in bam.pileup() {
        let pileup = p.unwrap();

        if *start <= (pileup.pos() as u64) && (pileup.pos() as u64) < *stop {
            for alignment in pileup.alignments() {
                let qname = String::from_utf8_lossy(alignment.record().qname()).into_owned();
                let sm = get_sm_name_from_rg(&alignment, &rg_sm_map)?;

                let seq_name = format!("{}|{}", qname, sm);

                if !bmap.contains_key(&seq_name) {
                    bmap.insert(seq_name.to_owned(), String::new());
                }

                if !alignment.is_del() && !alignment.is_refskip() {
                    let a = alignment.record().seq()[alignment.qpos().unwrap()];

                    bmap.get_mut(&seq_name).unwrap().push(a as char);
                }

                match alignment.indel() {
                    bam::pileup::Indel::Ins(len) => {
                        let pos1 = alignment.qpos().unwrap() as usize;
                        let pos2 = pos1 + (len as usize);
                        for pos in pos1..pos2 {
                            let a = alignment.record().seq()[pos];

                            bmap.get_mut(&seq_name).unwrap().push(a as char);
                        }
                    }
                    _ => {}
                }
            }
        }
    }

    let records = bmap
        .iter()
        .map(|kv| fasta::Record::with_attrs(kv.0.as_str(), None, kv.1.as_bytes()))
        .collect();

    Ok(records)
}

// Function to extract seqs from a FASTA file within a specified genomic region.
fn extract_fasta_seqs(basename: &String, fasta: &mut Reader, chr: &String, start: &u64, stop: &u64) -> Result<Vec<fasta::Record>> {
    let id = format!("{}:{}-{}|{}", chr, start, stop, basename);
<<<<<<< HEAD
    let seq = fasta.fetch_seq_string(chr, *start as usize, (*stop - 1) as usize).unwrap();
=======
    let seq = fasta.fetch_seq_string(chr, *start as usize, *stop as usize).unwrap();
>>>>>>> bfdb82b (Simplifying build command)

    let records = vec![fasta::Record::with_attrs(id.as_str(), None, seq.as_bytes())];

    Ok(records)
}

// Function to stage data from a single BAM file.
fn stage_data_from_one_file(
    seqs_url: &Url,
    loci: &HashSet<(String, u64, u64)>,
    cache_path: &PathBuf,
) -> Result<Vec<fasta::Record>> {
    let mut all_seqs = Vec::new();

    let basename = seqs_url.path_segments().map(|c| c.collect::<Vec<_>>()).unwrap().last().unwrap().to_string();

    let seqs_str = seqs_url.as_str();
    if seqs_str.ends_with(".bam") {
        // Handle BAM file processing
        let basename = basename
            .trim_end_matches(".bam")
            .to_string();
        let mut bam = open_bam(seqs_url, cache_path)?;

        for (chr, start, stop) in loci.iter() {
            // Extract seqs for the current locus.
            let seqs = extract_bam_reads(&basename, &mut bam, chr, start, stop).unwrap();

            // Extend the all_seqs vector with the seqs from the current locus.
            all_seqs.extend(seqs);
        }
    } else if seqs_str.ends_with(".fa") || seqs_str.ends_with(".fasta") || seqs_str.ends_with(".fa.gz") || seqs_str.ends_with(".fasta.gz") {
        // Handle FASTA file processing
        let basename = basename
            .trim_end_matches(".fa")
            .trim_end_matches(".fasta")
            .trim_end_matches(".fa.gz")
            .trim_end_matches(".fasta.gz")
            .to_string();
        let mut fasta = open_fasta(seqs_url, cache_path)?;

        for (chr, start, stop) in loci.iter() {
            // Extract seqs for the current locus.
            let seqs = extract_fasta_seqs(&basename, &mut fasta, chr, start, stop).unwrap();

            // Extend the all_seqs vector with the seqs from the current locus.
            all_seqs.extend(seqs);
        }
    } else if seqs_str.ends_with(".cram") {
        // Handle CRAM file processing

        !unimplemented!();
    } else {
        // Handle unknown file extension
        return Err(anyhow::anyhow!("Unsupported file type: {}", seqs_url));
    }

    Ok(all_seqs)
}

// Function to stage data from multiple BAM files.
fn stage_data_from_all_files(
    seq_urls: &HashSet<Url>,
    loci: &HashSet<(String, u64, u64)>,
    cache_path: &PathBuf,
) -> Result<Vec<fasta::Record>> {

    // Use a parallel iterator to process multiple BAM files concurrently.
    let all_data: Vec<_> = seq_urls
        .par_iter()
        .map(|seqs_url| {
            // Define an operation to stage data from one file.
            let op = || {
                let seqs = stage_data_from_one_file(seqs_url, loci, cache_path)?;
                Ok(seqs)
            };

            // Retry the operation with exponential backoff in case of failure.
            match backoff::retry(ExponentialBackoff::default(), op) {
                Ok(seqs) => { seqs }
                Err(e) => {
                    // If all retries fail, panic with an error message.
                    panic!("Error: {}", e);
                }
            }
        })
        .collect();

    let flattened_data = all_data.into_iter().flatten().collect::<Vec<_>>();

    // Return a flattened vector of sequences
    Ok(flattened_data)
}

// Function to retrieve the header of a BAM file.
fn get_bam_header(
    seqs_url: &Url,
    cache_path: &PathBuf,
) -> Result<bam::HeaderView> {
    // Open the BAM file.
    let bam = open_bam(seqs_url, cache_path)?;

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
    seq_urls: &HashSet<Url>,
    cache_path: &PathBuf
) -> Result<()> {
    // Stage data from all BAM files.
    let all_data = match stage_data_from_all_files(seq_urls, loci, cache_path) {
        Ok(all_data) => all_data,
        Err(e) => {
            panic!("Error: {}", e);
        }
    };

    // Write to a FASTA file.
    let mut buf_writer = BufWriter::new(File::create(output_path)?);
    let mut fasta_writer = fasta::Writer::new(&mut buf_writer);

    for record in all_data.iter() {
        fasta_writer.write_record(record)?;
    }

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
        let seqs_url = Url::parse(
            "gs://fc-8c3900db-633f-477f-96b3-fb31ae265c44/results/PBFlowcell/m84060_230907_210011_s2/seqs/ccs/aligned/m84060_230907_210011_s2.bam"
        ).unwrap();
        let cache_path = std::env::temp_dir();

        let bam = open_bam(&seqs_url, &cache_path);

        assert!(bam.is_ok(), "Failed to open bam file");
    }

    #[test]
    fn test_stage_data_from_one_file() {
        let seqs_url = Url::parse(
            "gs://fc-8c3900db-633f-477f-96b3-fb31ae265c44/results/PBFlowcell/m84060_230907_210011_s2/seqs/ccs/aligned/m84060_230907_210011_s2.bam"
        ).unwrap();
        let loci = HashSet::from([("chr15".to_string(), 23960193, 23963918)]);
        let cache_path = std::env::temp_dir();

        let result = stage_data_from_one_file(&seqs_url, &loci, &cache_path);

        assert!(result.is_ok(), "Failed to stage data from one file");
    }

    #[test]
    fn test_stage_data() {
        let cache_path = std::env::temp_dir();
        let output_path = cache_path.join("test.bam");

        let loci = HashSet::from([("chr15".to_string(), 23960193, 23963918)]);
        let seqs_url = Url::parse(
            "gs://fc-8c3900db-633f-477f-96b3-fb31ae265c44/results/PBFlowcell/m84060_230907_210011_s2/seqs/ccs/aligned/m84060_230907_210011_s2.bam"
        ).unwrap();
        let seq_urls = HashSet::from([seqs_url]);

        let result = stage_data(&output_path, &loci, &seq_urls, &cache_path);

        assert!(result.is_ok(), "Failed to stage data from file");

        println!("{:?}", result.unwrap());
    }

    #[test]
    fn test_stage_multiple_data() {
        let cache_path = std::env::temp_dir();
        let output_path = cache_path.join("test.bam");

        let seqs_url_1 = Url::parse(
            "gs://fc-8c3900db-633f-477f-96b3-fb31ae265c44/results/PBFlowcell/m84060_230907_210011_s2/seqs/ccs/aligned/m84060_230907_210011_s2.bam"
        ).unwrap();
        let seqs_url_2 = Url::parse(
            "gs://fc-8c3900db-633f-477f-96b3-fb31ae265c44/results/PBFlowcell/m84043_230901_211947_s1/seqs/ccs/aligned/m84043_230901_211947_s1.bam"
        ).unwrap();
        let loci = HashSet::from([("chr15".to_string(), 23960193, 23963918)]);
        let seq_urls = HashSet::from([ seqs_url_1, seqs_url_2 ]);

        let result = stage_data(&output_path, &loci, &seq_urls, &cache_path);

        println!("{:?}", result);

        assert!(result.is_ok(), "Failed to stage data from all files");
    }
}
