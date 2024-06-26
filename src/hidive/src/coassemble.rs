use std::path::PathBuf;

use skydive;
use skydive::ldbg::LdBG;
use skydive::mldbg::MLdBG;

pub fn start(output: &PathBuf, long_read_fasta_paths: &Vec<PathBuf>, short_read_fasta_paths: &Vec<PathBuf>) {
    let long_read_seq_urls = skydive::parse::parse_file_names(long_read_fasta_paths);
    let short_read_seq_urls = skydive::parse::parse_file_names(short_read_fasta_paths);

    let mut g = MLdBG::new();

    for seq_url in long_read_seq_urls {
        let basename = seq_url.path_segments().map(|c| c.collect::<Vec<_>>()).unwrap().last().unwrap().to_string();
        let basename = basename
            .trim_end_matches(".fasta.gz")
            .trim_end_matches(".fa.gz")
            .trim_end_matches(".fasta")
            .trim_end_matches(".fa")
            .to_string();

        let fasta_path = seq_url.to_file_path().unwrap();

        skydive::elog!("Processing {}...", basename);
        let l = LdBG::from_file(basename, 11, &fasta_path, true);

        g.push(l);
    }

    // sleep(Duration::from_secs(1000));
}
