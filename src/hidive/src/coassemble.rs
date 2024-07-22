use std::path::PathBuf;

use skydive;
use skydive::ldbg::LdBG;
use skydive::mldbg::MLdBG;

pub fn start(
    output: &PathBuf,
    long_read_fasta_paths: &Vec<PathBuf>,
    short_read_fasta_paths: &Vec<PathBuf>,
) {
    let long_read_seq_urls = skydive::parse::parse_file_names(long_read_fasta_paths);
    let short_read_seq_urls = skydive::parse::parse_file_names(short_read_fasta_paths);

    let mut g = MLdBG::new();

    for seq_url in long_read_seq_urls {
        let basename = skydive::utils::basename_without_extension(&seq_url, &[".fasta.gz", ".fa.gz", ".fasta", ".fa"]);

        let fasta_path = seq_url.to_file_path().unwrap();

        skydive::elog!("Processing {}...", basename);
        let l = LdBG::from_file(basename, 11, &fasta_path, true);

        let unitigs = l.assemble_all();
        for unitig in unitigs {
            println!("{}", String::from_utf8_lossy(&unitig));
        }

        g.push(l);
    }
}
