// Import necessary standard library modules
use std::collections::HashSet;
use std::io::Write;
use std::{collections::HashMap, path::PathBuf};

// Import the skydive module, which contains the necessary functions for staging data
use skydive;

// Import types from rust_htslib for working with BAM files.
use rust_htslib::bam::{
    self


    , IndexedReader, Read,
};

pub fn start(output: &PathBuf, loci_list: &Vec<String>, bam_path: &PathBuf) {
    // Initialize a HashSet to store unique loci after parsing
    let mut loci = HashSet::new();

    // Iterate over each locus in the provided list
    for locus in loci_list {
        // Attempt to parse the locus using a function from the skydive module
        match skydive::parse::parse_locus(locus.to_owned()) {
            Ok(l_fmt) => {
                // If parsing is successful, insert the formatted locus into the HashSet
                loci.insert(l_fmt);
            }
            Err(_) => {
                // If parsing fails, panic and terminate the program, providing an error message
                panic!("Could not parse locus '{}'.", locus);
            }
        }
    }

    // Open output file for writing.
    let mut output_file = std::fs::File::create(output).expect("Unable to create output file");

    let mut bam = IndexedReader::from_path(bam_path).unwrap();
    let rg_sm_map = get_rg_to_sm_mapping(&bam);

    loci.iter().for_each(|(chr, start, stop)| {
        let _ = bam.fetch(((*chr).as_bytes(), *start, *stop));

        let mut bmap = HashMap::new();

        for p in bam.pileup() {
            let pileup = p.unwrap();

            if *start <= (pileup.pos() as u64) && (pileup.pos() as u64) < *stop {
                for alignment in pileup.alignments() {
                    let qname = String::from_utf8_lossy(alignment.record().qname()).into_owned();
                    if !bmap.contains_key(&qname) {
                        bmap.insert(qname.to_owned(), String::new());
                    }

                    if !alignment.is_del() && !alignment.is_refskip() {
                        let a = alignment.record().seq()[alignment.qpos().unwrap()];

                        bmap.get_mut(&qname).unwrap().push(a as char);
                    }

                    match alignment.indel() {
                        bam::pileup::Indel::Ins(len) => {
                            let pos1 = alignment.qpos().unwrap() as usize;
                            let pos2 = pos1 + (len as usize);
                            for pos in pos1..pos2 {
                                let a = alignment.record().seq()[pos];

                                bmap.get_mut(&qname).unwrap().push(a as char);
                            }
                        }
                        _ => {}
                    }
                }
            }
        }

        for kv in bmap {
            let output_str = format!(">{}\n{}", kv.0, kv.1);
            writeln!(output_file, "{}", output_str).expect("Unable to write to file");
        }

        /*
        for (_, r) in bam.records().enumerate() {
            let record = r.unwrap();

            let mut ref_pos: u32 = (record.reference_start() as u32) + 1;
            let mut read_pos: u32 = 0;

            let mut sequence = Vec::new();

            for c in record.cigar().iter() {
                let seq_bytes = record.seq().as_bytes();

                match c {
                    Cigar::Match(len) | Cigar::Equal(len) => {
                        // Handle Match and Equal cases (both consume query, ref)
                        for _ in 0..*len {
                            let cigar_seq: &[u8] = &[record.seq()[read_pos as usize]];

                            if *start <= (ref_pos as u64) && (ref_pos as u64) < *stop {
                                sequence.push(String::from_utf8_lossy(cigar_seq).into_owned());
                            }

                            ref_pos += 1;
                            read_pos += 1;
                        }
                    }
                    Cigar::Diff(len) => {
                        // Handle Difference case (consumes query, ref)
                        let cigar_seq: &[u8] = &[record.seq()[read_pos as usize]];

                        if *start <= ref_pos as u64 && ref_pos as u64 <= *stop {
                            sequence.push(String::from_utf8_lossy(cigar_seq).into_owned());
                        }

                        ref_pos += len;
                        read_pos += len;
                    }
                    Cigar::Ins(len) => {
                        // Handle Insertion case (consumes query)
                        let cigar_start = read_pos as usize;
                        let cigar_end = (read_pos + *len) as usize;
                        let cigar_seq = &seq_bytes[cigar_start..cigar_end];

                        if *start <= ref_pos as u64 && ref_pos as u64 <= *stop {
                            sequence.push(String::from_utf8_lossy(cigar_seq).into_owned());
                        }

                        read_pos += len;
                    }
                    Cigar::Del(len) => {
                        // Handle Deletion case (consumes ref)
                        let cigar_seq = &[];

                        if *start <= ref_pos as u64 && ref_pos as u64 <= *stop {
                            sequence.push(String::from_utf8_lossy(cigar_seq).into_owned());
                        }

                        ref_pos += len;
                    }
                    Cigar::RefSkip(len) => {
                        // Handle Reference Skip case (consumes ref)
                        ref_pos += len;
                    }
                    Cigar::SoftClip(len) => {
                        // Handle Soft Clip case (consumes query)
                        read_pos += len;
                    }
                    Cigar::HardClip(_) => {
                        // Handle Hard Clip case (consumes nothing)
                    }
                    Cigar::Pad(_) => {
                        // Handle Padding case (consumes nothing)
                    }
                }
            }

            if let Ok(Aux::String(rg)) = record.aux(b"RG") {
                let sample_name = rg_sm_map.get(rg).unwrap().to_owned();
                let sequence_str = sequence.join("");
                let output_str = format!(">{}.{}\n{}", sample_name, String::from_utf8_lossy(record.qname()), sequence_str);
                writeln!(output_file, "{}", output_str).expect("Unable to write to file");
            }
        }
        */
    });
}

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
