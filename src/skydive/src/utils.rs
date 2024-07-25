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