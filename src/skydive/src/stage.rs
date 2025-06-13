// Import the Result type from the anyhow crate for error handling.
use anyhow::Result;
use linked_hash_set::LinkedHashSet;
use parquet::data_type::AsBytes;
use rust_htslib::bam::record::Aux;

// Import various standard library collections.
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

// Import the Url type to work with URLs.
use url::Url;

// Import ExponentialBackoff for retrying operations.
use backoff::ExponentialBackoff;

// Import rayon's parallel iterator traits.
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

// Import types from rust_htslib for working with BAM files.
use bio::io::fastq;
use rust_htslib::bam::{self, FetchDefinition, IndexedReader, Read};
use rust_htslib::faidx::Reader;

// Import functions for authorizing access to Google Cloud Storage.
use crate::env::{gcs_authorize_data_access, local_guess_curl_ca_bundle};

/// Function to open a BAM/CRAM file from a URL and cache its contents locally.
///
/// # Arguments
///
/// * `seqs_url` - A reference to a URL object representing the sequence file URL.
///
/// # Returns
///
/// An `IndexedReader` object representing the opened BAM/CRAM file.
///
/// # Errors
///
/// This function returns an error if the BAM/CRAM file cannot be opened.
///
/// # Panics
///
/// This function panics if the URL scheme is not recognized.
pub fn open_bam(seqs_url: &Url) -> Result<IndexedReader> {
    if seqs_url.to_string().starts_with("gs://") && env::var("GCS_OAUTH_TOKEN").is_err() {
        gcs_authorize_data_access();
    }

    // Try to open the BAM file from the URL, with retries for authorization.
    let bam = match IndexedReader::from_url(seqs_url) {
        Ok(bam) => bam,
        Err(_) => {
            crate::elog!("Read '{}', attempt 2 (reauthorizing to GCS)", seqs_url);

            // If opening fails, try authorizing access to Google Cloud Storage.
            gcs_authorize_data_access();

            // Try opening the BAM file again.
            match IndexedReader::from_url(seqs_url) {
                Ok(bam) => bam,
                Err(_) => {
                    crate::elog!("Read '{}', attempt 3 (overriding cURL CA bundle)", seqs_url);

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

/// Function to open a FASTA file from a URL and cache its contents locally.
///
/// # Arguments
///
/// * `seqs_url` - A reference to a URL object representing the sequence file URL.
///
/// # Returns
///
/// A `Reader` object representing the opened FASTA file.
///
/// # Errors
///
/// This function returns an error if the FASTA file cannot be opened.
///
/// # Panics
///
/// This function panics if the URL scheme is not recognized.
pub fn open_fasta(seqs_url: &Url) -> Result<Reader> {
    if seqs_url.to_string().starts_with("gs://") && env::var("GCS_OAUTH_TOKEN").is_err() {
        gcs_authorize_data_access();
    }

    // Try to open the FASTA file from the URL, with retries for authorization.
    let fasta = match Reader::from_url(seqs_url) {
        Ok(fasta) => fasta,
        Err(_) => {
            crate::elog!("Read '{}', attempt 2 (reauthorizing to GCS)", seqs_url);

            // If opening fails, try authorizing access to Google Cloud Storage.
            gcs_authorize_data_access();

            // Try opening the BAM file again.
            match Reader::from_url(seqs_url) {
                Ok(bam) => bam,
                Err(_) => {
                    crate::elog!("Read '{}', attempt 3 (overriding cURL CA bundle)", seqs_url);

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
        .map(|record| (record["ID"].clone(), record["SM"].clone()))
        .collect();

    rg_sm_map
}

fn get_sm_name_from_rg(read: &bam::Record, rg_sm_map: &HashMap<String, String>) -> Result<String> {
    let rg = read.aux(b"RG")?;

    if let Aux::String(v) = rg {
        if let Some(sm) = rg_sm_map.get(v) {
            Ok(sm.to_owned())
        } else {
            Err(anyhow::anyhow!(
                "Sample name not found for read group: {}",
                v
            ))
        }
    } else {
        Err(anyhow::anyhow!("Read group is not a string"))
    }
}

// Function to extract seqs from a BAM file within a specified genomic region.
pub fn extract_aligned_bam_reads(
    _basename: &str,
    bam: &mut IndexedReader,
    chr: &str,
    start: &u64,
    stop: &u64,
    name: &str,
    haplotype: Option<u8>,
) -> Result<Vec<fastq::Record>> {
    let rg_sm_map = get_rg_to_sm_mapping(bam);

    let mut bmap = HashMap::new();

    let _ = bam.fetch(((*chr).as_bytes(), *start, *stop));
    for p in bam.pileup() {
        let pileup = p.unwrap();

        if *start <= (pileup.pos() as u64) && (pileup.pos() as u64) < *stop {
            for (i, alignment) in pileup.alignments().enumerate().filter(|(_, a)| {
                haplotype.is_none()
                    || a.record()
                        .aux(b"HP")
                        .ok()
                        .map(|aux| match aux {
                            Aux::U8(v) => v == haplotype.unwrap(),
                            _ => false,
                        })
                        .unwrap_or(false)
            }) {
                let qname = String::from_utf8_lossy(alignment.record().qname()).into_owned();
                let sm = match get_sm_name_from_rg(&alignment.record(), &rg_sm_map) {
                    Ok(a) => a,
                    Err(_) => String::from("unknown"),
                };

                let is_secondary = alignment.record().is_secondary();
                let is_supplementary = alignment.record().is_supplementary();
                let seq_name = format!("{qname}|{name}|{sm}|{i}|{is_secondary}|{is_supplementary}");

                if !bmap.contains_key(&seq_name) {
                    bmap.insert(seq_name.clone(), (String::new(), Vec::new()));
                }

                if !alignment.is_del() && !alignment.is_refskip() {
                    let a = alignment.record().seq()[alignment.qpos().unwrap()];
                    let q = alignment.record().qual()[alignment.qpos().unwrap()];

                    bmap.get_mut(&seq_name).unwrap().0.push(a as char);
                    bmap.get_mut(&seq_name).unwrap().1.push(q + 33 as u8);
                }

                if let bam::pileup::Indel::Ins(len) = alignment.indel() {
                    if let Some(pos1) = alignment.qpos() {
                        let pos2 = pos1 + (len as usize);
                        for pos in pos1..pos2 {
                            let a = alignment.record().seq()[pos];
                            let q = alignment.record().qual()[pos];

                            bmap.get_mut(&seq_name).unwrap().0.push(a as char);
                            bmap.get_mut(&seq_name).unwrap().1.push(q + 33 as u8);
                        }
                    }
                }
            }
        }
    }

    let records = bmap
        .iter()
        .map(|kv| fastq::Record::with_attrs(kv.0.as_str(), None, kv.1.0.as_bytes(), kv.1.1.as_bytes()))
        .collect();

    Ok(records)
}

// Function to extract unaligned seqs from a BAM file
fn extract_unaligned_bam_reads(
    _basename: &str,
    bam: &mut IndexedReader,
) -> Result<Vec<fastq::Record>> {
    let rg_sm_map = get_rg_to_sm_mapping(bam);

    let _ = bam.fetch(FetchDefinition::Unmapped);
    let records = bam
        .records()
        .map(|r| {
            let read = r.unwrap();
            let qname = String::from_utf8_lossy(read.qname()).into_owned();
            let sm = get_sm_name_from_rg(&read, &rg_sm_map).unwrap();

            let seq_name = format!("{qname}|{sm}");

            let vseq = read.seq().as_bytes();
            let bseq = vseq.as_bytes();
            let qseq = read.qual().iter().map(|&q| q + 33).collect::<Vec<u8>>();

            let seq = fastq::Record::with_attrs(seq_name.as_str(), Some(""), bseq, &qseq);

            seq
        })
        .collect();

    Ok(records)
}

// Function to extract seqs from a FASTA file within a specified genomic region.
fn extract_fasta_seqs(
    basename: &String,
    fasta: &mut Reader,
    chr: &String,
    start: &u64,
    stop: &u64,
    name: &String,
) -> Result<Vec<fastq::Record>> {
    let id = format!("{chr}:{start}-{stop}|{name}|{basename}");
    let seq = fasta.fetch_seq_string(chr, usize::try_from(*start)?, usize::try_from(*stop - 1)?)?;

    if seq.len() > 0 {
        let records = vec![fastq::Record::with_attrs(id.as_str(), None, seq.as_bytes(), vec![30; seq.len()].as_slice())];

        return Ok(records);
    }

    Err(anyhow::anyhow!("No sequence found for locus: {}", id))
}

// Function to stage data from a single file.
fn stage_data_from_one_file(
    seqs_url: &Url,
    loci: &LinkedHashSet<(String, u64, u64, String)>,
    unmapped: bool,
    haplotype: Option<u8>,
) -> Result<Vec<fastq::Record>> {
    let mut all_seqs = Vec::new();

    let basename = seqs_url
        .path_segments()
        .map(|c| c.collect::<Vec<_>>())
        .unwrap()
        .last()
        .unwrap()
        .to_string();

    let seqs_str = seqs_url.as_str();
    if seqs_str.ends_with(".bam") || seqs_str.ends_with(".cram") {
        // Handle BAM/CRAM file processing
        let basename = basename
            .trim_end_matches(".bam")
            .trim_end_matches(".cram")
            .to_string();
        let mut bam = open_bam(seqs_url)?;

        // Extract seqs for the current locus.
        for (chr, start, stop, name) in loci.iter() {
            let aligned_seqs =
                extract_aligned_bam_reads(&basename, &mut bam, chr, start, stop, name, haplotype)
                    .unwrap();
            all_seqs.extend(aligned_seqs);
        }

        // Optionally extract unaligned reads.
        if unmapped {
            let unaligned_seqs = extract_unaligned_bam_reads(&basename, &mut bam).unwrap();
            all_seqs.extend(unaligned_seqs);
        }
    } else if seqs_str.ends_with(".fa")
        || seqs_str.ends_with(".fasta")
        || seqs_str.ends_with(".fa.gz")
        || seqs_str.ends_with(".fasta.gz")
    {
        // Handle FASTA file processing
        let basename = basename
            .trim_end_matches(".fasta.gz")
            .trim_end_matches(".fasta")
            .trim_end_matches(".fa.gz")
            .trim_end_matches(".fa")
            .to_string();
        let mut fasta = open_fasta(seqs_url)?;

        for (chr, start, stop, name) in loci.iter() {
            // Extract seqs for the current locus.
            let seqs = extract_fasta_seqs(&basename, &mut fasta, chr, start, stop, name)
                .map_or_else(|_| Vec::new(), |s| s);

            // Extend the all_seqs vector with the seqs from the current locus.
            all_seqs.extend(seqs);
        }
    } else {
        // Handle unknown file extension
        return Err(anyhow::anyhow!("Unsupported file type: {}", seqs_url));
    }

    Ok(all_seqs)
}

// Function to stage data from multiple BAM files.
fn stage_data_from_all_files(
    seq_urls: &HashSet<Url>,
    loci: &LinkedHashSet<(String, u64, u64, String)>,
    unmapped: bool,
    haplotype: Option<u8>,
) -> Result<Vec<fastq::Record>> {
    // Use a parallel iterator to process multiple BAM files concurrently.
    let all_data: Vec<_> = seq_urls
        .par_iter()
        .map(|seqs_url| {
            // Define an operation to stage data from one file.
            let op = || {
                let seqs = stage_data_from_one_file(seqs_url, loci, unmapped, haplotype)?;
                Ok(seqs)
            };

            // Retry the operation with exponential backoff in case of failure.
            match backoff::retry(ExponentialBackoff::default(), op) {
                Ok(seqs) => seqs,
                Err(e) => {
                    // If all retries fail, panic with an error message.
                    panic!("Error: {e}");
                }
            }
        })
        .collect();

    let flattened_data = all_data.into_iter().flatten().collect::<Vec<_>>();

    // Return a flattened vector of sequences
    Ok(flattened_data)
}

/// Checks if a given genomic range spans any of the loci in the provided set.
///
/// # Arguments
///
/// * `start` - The start position of the genomic range.
/// * `end` - The end position of the genomic range.
/// * `loci` - A reference to a `HashSet` of tuples representing the loci.
///
/// # Returns
///
/// A `Result` containing `true` if the range spans any loci, `false` otherwise. Returns an error if the positions are negative.
///
/// # Errors
///
/// Returns an error if the start or end positions are negative.
///
/// # Panics
///
/// This function does not panic.
pub fn read_spans_locus(
    start: i64,
    end: i64,
    loci: &HashSet<(String, i64, i64)>,
) -> Result<bool, String> {
    if start < 0 || end < 0 {
        return Err("Error: Negative genomic positions are not allowed.".to_string());
    }

    Ok(loci.iter().any(|e| start <= e.1 && end >= e.2))
}

/// Public function to stage data from multiple BAM files and write to an output file.
///
/// # Arguments
///
/// * `output_path` - A reference to a `PathBuf` representing the output file path.
/// * `loci` - A reference to a `HashSet` of tuples representing the loci to extract.
/// * `seq_urls` - A reference to a `HashSet` of URLs representing the sequence files.
/// * `unmapped` - A boolean indicating whether to extract unmapped reads.
/// * `cache_path` - A reference to a `PathBuf` representing the cache directory path.
///
/// # Returns
///
/// The number of records written to the output file.
///
/// # Errors
///
/// This function returns an error if the output file cannot be created.
///
/// # Panics
///
/// If an error occurs while staging data from the files.
///
pub fn stage_data(
    output_path: &PathBuf,
    loci: &LinkedHashSet<(String, u64, u64, String)>,
    seq_urls: &HashSet<Url>,
    unmapped: bool,
    haplotype: Option<u8>,
    cache_path: &PathBuf,
) -> Result<usize> {
    let current_dir = env::current_dir()?;
    env::set_current_dir(cache_path).unwrap();

    // Stage data from all BAM files.
    let all_data = match stage_data_from_all_files(seq_urls, loci, unmapped, haplotype) {
        Ok(all_data) => all_data,
        Err(e) => {
            panic!("Error: {e}");
        }
    };

    env::set_current_dir(current_dir).unwrap();

    // Write to a FASTQ file.
    let mut buf_writer = BufWriter::new(File::create(output_path)?);
    let mut fastq_writer = fastq::Writer::new(&mut buf_writer);

    for record in all_data.iter() {
        if record.seq().len() > 0 {
            fastq_writer.write_record(record)?;
        }
    }

    let _ = fastq_writer.flush();

    Ok(all_data.len())
}

pub fn stage_data_in_memory(
    loci: &LinkedHashSet<(String, u64, u64, String)>,
    seq_urls: &HashSet<Url>,
    unmapped: bool,
    cache_path: &PathBuf,
) -> Result<Vec<fastq::Record>> {
    let current_dir = env::current_dir()?;
    env::set_current_dir(cache_path).unwrap();

    // Stage data from all BAM files.
    let all_data = match stage_data_from_all_files(seq_urls, loci, unmapped, None) {
        Ok(all_data) => all_data,
        Err(e) => {
            panic!("Error: {e}");
        }
    };

    env::set_current_dir(current_dir).unwrap();

    Ok(all_data)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;
    use url::Url;

    // This test may pass, but still print a message to stderr regarding its failure to access data. This is because
    // open_bam() tries a couple of authorization methods before accessing data, and the initial failures print a
    // message to stderr.
    #[test]
    fn test_open_bam() {
        let seqs_url = Url::parse(
            "gs://fc-8c3900db-633f-477f-96b3-fb31ae265c44/results/PBFlowcell/m84060_230907_210011_s2/reads/ccs/aligned/m84060_230907_210011_s2.bam"
        ).unwrap();
        let bam = open_bam(&seqs_url);

        assert!(bam.is_ok(), "Failed to open bam file");
    }

    #[test]
    fn test_stage_data_from_one_file() {
        let seqs_url = Url::parse(
            "gs://fc-8c3900db-633f-477f-96b3-fb31ae265c44/results/PBFlowcell/m84060_230907_210011_s2/reads/ccs/aligned/m84060_230907_210011_s2.bam"
        ).unwrap();
        let mut loci = LinkedHashSet::new();
        loci.insert(("chr15".to_string(), 23960193, 23963918, "test".to_string()));

        let result = stage_data_from_one_file(&seqs_url, &loci, false, None);

        assert!(result.is_ok(), "Failed to stage data from one file");
    }

    #[test]
    fn test_stage_data() {
        let cache_path = std::env::temp_dir();
        let output_path = cache_path.join("test.bam");

        let mut loci = LinkedHashSet::new();
        loci.insert(("chr15".to_string(), 23960193, 23963918, "test".to_string()));

        let seqs_url = Url::parse(
            "gs://fc-8c3900db-633f-477f-96b3-fb31ae265c44/results/PBFlowcell/m84060_230907_210011_s2/reads/ccs/aligned/m84060_230907_210011_s2.bam"
        ).unwrap();
        let seq_urls = HashSet::from([seqs_url]);

        let result = stage_data(&output_path, &loci, &seq_urls, false, None, &cache_path);

        assert!(result.is_ok(), "Failed to stage data from file");
    }

    #[test]
    fn test_stage_multiple_data() {
        let cache_path = std::env::temp_dir();
        let output_path = cache_path.join("test.bam");

        let seqs_url_1 = Url::parse(
            "gs://fc-8c3900db-633f-477f-96b3-fb31ae265c44/results/PBFlowcell/m84060_230907_210011_s2/reads/ccs/aligned/m84060_230907_210011_s2.bam"
        ).unwrap();
        let seqs_url_2 = Url::parse(
            "gs://fc-8c3900db-633f-477f-96b3-fb31ae265c44/results/PBFlowcell/m84043_230901_211947_s1/reads/ccs/aligned/m84043_230901_211947_s1.hifi_reads.bc2080.bam"
        ).unwrap();
        let mut loci = LinkedHashSet::new();
        loci.insert(("chr15".to_string(), 23960193, 23963918, "test".to_string()));

        let seq_urls = HashSet::from([seqs_url_1, seqs_url_2]);

        let result = stage_data(&output_path, &loci, &seq_urls, false, None, &cache_path);

        assert!(result.is_ok(), "Failed to stage data from all files");
    }
}
