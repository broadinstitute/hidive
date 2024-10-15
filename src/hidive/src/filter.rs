use std::io::Read;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::{
    fs::File,
    io::{BufWriter, Write},
};

use flate2::read::GzDecoder;
use num_format::{Locale, ToFormattedString};

use rayon::prelude::*;

use minimap2::Aligner;
use tempfile::NamedTempFile;

pub fn start(output: &PathBuf, gfa_path: &PathBuf, short_read_fasta_paths: &Vec<PathBuf>) {
    let short_read_seq_urls = skydive::parse::parse_file_names(short_read_fasta_paths);

    let g = skydive::utils::read_gfa(gfa_path).unwrap();

    skydive::elog!("{} contigs", g.node_count());

    let temp_file = NamedTempFile::new().expect("Failed to create temporary file");
    let temp_path = temp_file.path().to_owned();

    {
        let file = File::create(&temp_path).expect("Failed to open temporary file for writing");
        let mut writer = BufWriter::new(file);

        for (index, node) in g.node_indices().enumerate() {
            let node_str = format!(">contig_{}\n{}\n", index, g[node]);
            writer
                .write_all(node_str.as_bytes())
                .expect("Failed to write to temporary file");
        }
    }

    skydive::elog!(
        "Wrote {} nodes to temporary file: {:?}",
        g.node_count(),
        temp_path
    );

    let aligner = Aligner::builder()
        .short()
        .with_cigar()
        .with_index(temp_path, None)
        .expect("Unable to build index");

    let mut filtered_sr_seqs: Vec<Vec<u8>> = Vec::new();
    for short_read_seq_url in &short_read_seq_urls {
        let basename = skydive::utils::basename_without_extension(
            &short_read_seq_url,
            &[".fasta.gz", ".fa.gz", ".fasta", ".fa"],
        );
        let fasta_path = short_read_seq_url.to_file_path().unwrap();

        skydive::elog!("Processing short-read sample {}...", basename);

        let file = File::open(&fasta_path).unwrap();
        let reader: Box<dyn Read> = if output.to_str().unwrap().ends_with(".gz") {
            Box::new(GzDecoder::new(file))
        } else {
            Box::new(file)
        };

        let fasta_reader = bio::io::fasta::Reader::new(reader);
        let all_reads = fasta_reader.records().flatten().collect::<Vec<_>>();

        let progress_bar = skydive::utils::default_unbounded_progress_bar(format!(
            "Filtering short reads (0 retained)",
        ));

        // Create some thread-safe counters.
        let found_items = Arc::new(AtomicUsize::new(0));
        let processed_items = Arc::new(AtomicUsize::new(0));
        let progress_bar = Arc::new(progress_bar);

        const UPDATE_FREQUENCY: usize = 1_000_000;

        filtered_sr_seqs.extend(
            all_reads
                .par_iter()
                .filter_map(|r| {
                    let sr_seq = r.seq().to_vec();

                    // Increment progress counters and update the progress bar every once in a while.
                    let current_processed = processed_items.fetch_add(1, Ordering::Relaxed);
                    if current_processed % UPDATE_FREQUENCY == 0 {
                        progress_bar.set_message(format!(
                            "Filtering short reads ({} retained)",
                            found_items
                                .load(Ordering::Relaxed)
                                .to_formatted_string(&Locale::en)
                        ));
                        progress_bar.inc(UPDATE_FREQUENCY as u64);
                    }

                    let count = aligner
                        .map(&sr_seq, false, false, None, None)
                        .iter()
                        .flatten()
                        .filter_map(|a| {
                            if a.target_name.is_some() {
                                Some(1)
                            } else {
                                None
                            }
                        })
                        .sum::<usize>();

                    if count > 0 {
                        found_items.fetch_add(1, Ordering::Relaxed);
                        Some(sr_seq)
                    } else {
                        None
                    }
                })
                .collect::<Vec<Vec<u8>>>(),
        );

        // Write filtered short read sequences to output
        let mut writer = std::io::BufWriter::new(std::fs::File::create(output).unwrap());
        for (i, seq) in filtered_sr_seqs.iter().enumerate() {
            let _ = writeln!(
                writer,
                ">filtered_read_{}\n{}",
                i,
                String::from_utf8(seq.to_vec()).unwrap()
            );
        }

        // Final update to ensure the last counts are displayed
        let final_found = found_items.load(Ordering::Relaxed);
        let final_processed = processed_items.load(Ordering::Relaxed);
        progress_bar.set_message(format!(
            "Filtered {} short reads, retained {}.",
            final_processed.to_formatted_string(&Locale::en),
            final_found.to_formatted_string(&Locale::en),
        ));
        progress_bar.finish();
    }
}
