use std::borrow::Cow;

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::collections::HashMap;

use needletail::Sequence;
use parquet::data_type::AsBytes;
use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::visit::{EdgeRef, NodeIndexable, NodeRef};

/// This function takes a sequence URL and a list of possible extensions, and returns the base name of the file
/// without any of the provided extensions. It does this by first extracting the last segment of the URL path,
/// and then iteratively removing any of the specified extensions from the end of the base name.
///
/// # Arguments
///
/// * `seq_url` - A reference to a URL object representing the sequence file URL.
/// * `extensions` - A slice of string slices representing the possible file extensions to be removed.
///
/// # Returns
///
/// * A `String` containing the base name of the file without any of the specified extensions.
///
/// # Example
///
/// ```
/// let url = url::Url::parse("http://example.com/path/to/file.fasta.gz").unwrap();
/// let extensions = [".fasta.gz", ".fa.gz", ".fasta", ".fa"];
/// let basename = basename_without_extension(&url, &extensions);
/// assert_eq!(basename, "file");
/// ```
/// # Panics
/// This function will panic:
/// 1. If `seq_url.path_segments()` returns `None`, indicating that the URL does not have a path.
/// 2. If `seq_url.path_segments().map(|c| c.collect::<Vec<_>>()).unwrap().last()` returns `None`,
/// indicating that the path does not have any segments.
#[must_use]
pub fn basename_without_extension(seq_url: &url::Url, extensions: &[&str]) -> String {
    let mut basename = seq_url
        .path_segments()
        .map(|c| c.collect::<Vec<_>>())
        .unwrap()
        .last()
        .unwrap()
        .to_string();

    let mut sorted_extensions = extensions.to_vec();
    sorted_extensions.sort_by_key(|b| std::cmp::Reverse(b.len()));

    for ext in sorted_extensions {
        basename = basename.trim_end_matches(ext).to_string();
    }

    basename
}

/// Given fasta files this function will read and return a list of lists containing the contents
/// of the fasta files
///
/// # Arguments
/// * `path`: paths to fasta files
///
/// # Returns
/// A list of lists containing the contents of the fasta files
///
/// # Panics
/// This function will panic if it cannot read a given file path
#[must_use]
pub fn read_fasta(paths: &Vec<PathBuf>) -> Vec<Vec<u8>> {
    paths
        .iter()
        .map(|p| {
            let reader = bio::io::fasta::Reader::from_file(p).expect("Failed to open file");
            reader
                .records()
                .filter_map(|r| r.ok())
                .map(|r| r.seq().to_vec())
                .collect::<Vec<Vec<u8>>>()
        })
        .flatten()
        .collect::<Vec<Vec<u8>>>()
}

#[must_use]
pub fn default_hidden_progress_bar() -> indicatif::ProgressBar {
    indicatif::ProgressBar::hidden()
}

/// Create a new bounded progress bar with the specified message and length.
/// The progress bar will be a bar with a spinner.
/// The progress bar will display the elapsed time, the progress bar, the current position,
/// the total length, and the estimated time remaining.
///
/// # Arguments
/// * `msg`: The message to display on the progress bar.
/// * `len`: The total length of the progress bar.
///
/// # Returns
/// A new bounded progress bar.
///
/// # Example
/// ```
/// let progress_bar = default_bounded_progress_bar("Processing sequences", 100);
/// ```
/// This will create a new bounded progress bar with the message "Processing sequences" and a total length of 100.
/// The progress bar will be a bar with a spinner.
///
/// # Panics
/// This function will panic if the progress bar style cannot be created.
pub fn default_bounded_progress_bar(
    msg: impl Into<Cow<'static, str>>,
    len: u64,
) -> indicatif::ProgressBar {
    let progress_bar_style = indicatif::ProgressStyle::default_bar()
        .template(
            "{msg} ... [{elapsed_precise}] [{bar:40.white/white}] {human_pos}/{human_len} ({eta})",
        )
        .unwrap()
        .progress_chars("#>-");

    let progress_bar = indicatif::ProgressBar::new(len);
    progress_bar.set_style(progress_bar_style);
    progress_bar.set_message(msg);

    progress_bar
}

/// Create a new unbounded progress bar with the specified message.
///
/// # Arguments
/// `msg`: The message to display on the progress bar.
///
/// # Returns
/// A new unbounded progress bar.
///
/// # Example
/// ```
/// let progress_bar = default_unbounded_progress_bar("Processing sequences");
/// ```
///
/// This will create a new unbounded progress bar with the message "Processing sequences".
/// The progress bar will be a spinner.
///
/// # Panics
/// This function will panic if the progress bar style cannot be created.
///
pub fn default_unbounded_progress_bar(msg: impl Into<Cow<'static, str>>) -> indicatif::ProgressBar {
    let progress_bar_style = indicatif::ProgressStyle::default_bar()
        .template("{msg} ... [{elapsed_precise}] {human_pos}")
        .unwrap()
        .progress_chars("#>-");

    let progress_bar = indicatif::ProgressBar::new_spinner();
    progress_bar.set_style(progress_bar_style);
    progress_bar.set_message(msg);

    progress_bar
}

/// Get the canonical (lexicographically-lowest) version of a k-mer.
///
/// # Arguments
///
/// * `kmer` - A slice representing the k-mer.
///
/// # Returns
///
/// A vector containing the canonical k-mer.
#[inline(always)]
#[must_use]
pub fn canonicalize_kmer(kmer: &[u8]) -> Vec<u8> {
    let rc_kmer = kmer.reverse_complement();
    if kmer < rc_kmer.as_bytes() {
        kmer.to_vec()
    } else {
        rc_kmer.as_bytes().to_vec()
    }
}

#[must_use]
pub fn homopolymer_compressed(seq: &[u8]) -> Vec<u8> {
    let mut compressed = Vec::new();
    let mut prev = None;

    for &base in seq {
        if Some(base) != prev {
            compressed.push(base);
        }
        prev = Some(base);
    }

    compressed
}

#[must_use]
pub fn shannon_entropy(seq: &[u8]) -> f32 {
    let mut freq = HashMap::new();
    let len = seq.len() as f32;

    for &base in seq {
        *freq.entry(base).or_insert(0) += 1;
    }

    -freq.values().map(|&count| {
        let p = count as f32 / len;
        p * p.log2()
    }).sum::<f32>()
}

#[must_use]
pub fn gc_content(seq: &[u8]) -> f32 {
    let gc_count = seq.iter().filter(|&&base| base == b'G' || base == b'C').count();
    gc_count as f32 / seq.len() as f32
}

pub fn write_gfa<W: std::io::Write>(writer: &mut W, graph: &DiGraph<String, f32>) -> std::io::Result<()> {
    // Write header
    writeln!(writer, "H\tVN:Z:1.0")?;

    // Write segments
    for (node_index, node_label) in graph.node_indices().zip(graph.node_weights()) {
        writeln!(writer, "S\t{}\t{}", node_index.index(), node_label)?;
    }

    // Write links
    for edge in graph.edge_references() {
        let (from, to) = (edge.source().index(), edge.target().index());
        let weight = edge.weight();
        writeln!(writer, "L\t{}\t+\t{}\t+\t0M\tRC:f:{}", from, to, (100.0*weight).round() as u8)?;
    }

    Ok(())
}

pub fn read_gfa<P: AsRef<Path>>(path: P) -> std::io::Result<DiGraph<String, f32>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut graph = DiGraph::new();
    let mut node_map = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        
        match fields[0] {
            "S" => {
                let id = fields[1];
                let sequence = fields[2].to_string();
                let node_index = graph.add_node(sequence);
                node_map.insert(id.to_string(), node_index);
            },
            "L" => {
                if fields.len() < 6 {
                    continue; // Skip malformed lines
                }
                let from_id = fields[1];
                // let from_orient = fields[2];
                let to_id = fields[3];
                // let to_orient = fields[4];

                let weight = fields.get(5)
                    .and_then(|s| s.split(':').last())
                    .and_then(|s| s.parse::<f32>().ok())
                    .unwrap_or(1.0);

                if let (Some(&from), Some(&to)) = (node_map.get(from_id), node_map.get(to_id)) {
                    graph.add_edge(from, to, weight);
                }
            },
            _ => {} // Ignore other lines
        }
    }

    Ok(graph)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_canonicalize_kmer() {
        let kmer1 = b"CGTA";
        let kmer2 = b"TACG";
        let kmer3 = b"AAAA";
        let kmer4 = b"TTTT";

        // Test canonical k-mer for kmer1 and kmer2
        assert_eq!(canonicalize_kmer(kmer1), b"CGTA".to_vec());
        assert_eq!(canonicalize_kmer(kmer2), b"CGTA".to_vec());

        // Test canonical k-mer for kmer3 and kmer4
        assert_eq!(canonicalize_kmer(kmer3), b"AAAA".to_vec());
        assert_eq!(canonicalize_kmer(kmer4), b"AAAA".to_vec());
    }
}
