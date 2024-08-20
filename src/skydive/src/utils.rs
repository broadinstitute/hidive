use std::borrow::Cow;

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
pub fn canonicalize_kmer(kmer: &[u8]) -> Vec<u8> {
    let rc_kmer = kmer.reverse_complement();
    if kmer < rc_kmer.as_bytes() {
        kmer.to_vec()
    } else {
        rc_kmer.as_bytes().to_vec()
    }
}

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

pub fn write_graph_as_gfa<W: std::io::Write>(writer: &mut W, graph: &DiGraph<String, f32>) -> std::io::Result<()> {
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
