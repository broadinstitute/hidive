use std::fs::File;
use std::path::PathBuf;

use skydive::ldbg::LdBG;
use skydive::mldbg::MLdBG;
use skydive::utils::*;

use std::io::Write;

pub fn start(
    output: &PathBuf,
    kmer_size: usize,
    model_path: &PathBuf,
    long_read_fasta_paths: &Vec<PathBuf>,
    short_read_fasta_paths: &Vec<PathBuf>,
) {
    let long_read_seq_urls = skydive::parse::parse_file_names(long_read_fasta_paths);
    let short_read_seq_urls = skydive::parse::parse_file_names(short_read_fasta_paths);

    // Read all long reads.
    skydive::elog!("Processing long-read samples {:?}...", long_read_seq_urls.iter().map(|url| url.as_str()).collect::<Vec<&str>>());
    let all_lr_seqs = skydive::utils::read_fasta(long_read_fasta_paths);

    // Read all short reads.
    skydive::elog!("Processing short-read samples {:?}...", short_read_seq_urls.iter().map(|url| url.as_str()).collect::<Vec<&str>>());
    let all_sr_seqs = skydive::utils::read_fasta(short_read_fasta_paths);

    let l1 = LdBG::from_sequences("lr".to_string(), kmer_size, &all_lr_seqs);
    let s1 = LdBG::from_sequences("sr".to_string(), kmer_size, &all_sr_seqs);

    let m = MLdBG::from_ldbgs(vec![l1, s1])
        .score_kmers(model_path)
        .collapse()
        .clean_paths(0.01)
        .clean_tips(3*kmer_size)
        .clean_contigs(100)
        ;

    let contigs = m.assemble_all();

    for (i, all_lr_seq) in all_lr_seqs.iter().enumerate() {
        // println!(">{}\n{}", i, String::from_utf8(contig.clone()).unwrap());
        // println!(">{}\n{}", i, String::from_utf8(all_lr_seq.clone()).unwrap());

        let kmers = all_lr_seq
            .windows(kmer_size)
            .map(|kmer| {
                let cn_kmer = skydive::utils::canonicalize_kmer(kmer);
                if m.kmers.contains_key(&cn_kmer) {
                    kmer.to_vec()
                } else {
                    "N".repeat(kmer_size).as_bytes().to_vec()
                }
            })
            .collect::<Vec<Vec<u8>>>();

        let first_kmer = kmers.first().unwrap();
        let mut seq = String::from_utf8(first_kmer[0..kmer_size-1].to_vec()).unwrap();

        for kmer in kmers {
            seq.push(kmer[kmer_size-1] as char);
        }

        // println!(">{}\n{}", i, seq);

        // let corrected_seqs = m.correct_seq(all_lr_seq);
        // for (j, corrected_seq) in corrected_seqs.iter().enumerate() {
        //     println!(">corrected_{}_{}\n{}", i, j, String::from_utf8(corrected_seq.clone()).unwrap());
        // }
    }

    let orig_seq = b"GATTCTCCCCAGACGCCGAGGATGGCCGTCATGGCGCCCCGAACCCTCGTCCTGCTACTCTCGGGGGCTCTGGCCCTGACCCAGACCTGGGCGGGTGAGTGCGGGGTCGGGGAGGGAAACGGCCTCTGTGGGGAGAAGCAACGGGCCCGCCTGGGCGGGGGCGCAGGACCCGGGAAGCCGCGCCGGGAGGAGGGTCGGGCGGGTCTCAAGCCACTCCTCTCCCCAGGCTCTCACTCCATGAGGTATTTCTACACCTCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCAGTGGG";
    let corr_seqs = m.correct_seq(orig_seq);

    skydive::elog!("Corrected sequence: {}", corr_seqs.len());

    for (i, corr_seq) in corr_seqs.iter().enumerate() {
        println!(">corrected_{}\n{}", i, String::from_utf8(corr_seq.clone()).unwrap());
    }

    let g = m.traverse_all_kmers();
    let _ = write_gfa(&mut File::create(output).unwrap(), &g);

    let csv_output = output.with_extension("csv");
    let mut csv_file = File::create(&csv_output).unwrap();

    for (node_index, node_label) in g.node_indices().zip(g.node_weights()) {
        let kmer = node_label.as_bytes();
        let cn_kmer = skydive::utils::canonicalize_kmer(kmer);
        let score = (100.0 * *m.scores.get(&cn_kmer).unwrap()) as u32;

        writeln!(
            csv_file,
            "{},{}",
            node_index.index(),
            score
        )
        .unwrap();
    }
}