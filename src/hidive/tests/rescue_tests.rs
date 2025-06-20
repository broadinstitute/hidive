use std::path::PathBuf;
use std::process::Command;
use tempfile::TempDir;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use bio::io::{fasta, fastq};

/// Validates that a FASTQ file exists and contains valid records
fn validate_fastq_output(path: &PathBuf) {
    let reader = bio::io::fastq::Reader::from_file(path).expect("Failed to open output file");
    let mut records = Vec::new();

    for result in reader.records() {
        // Each record should be valid and parseable
        let record = result.expect("Error reading FASTQ record");
        records.push(record);
    }

    // Basic validation that we got some records
    assert!(!records.is_empty(), "No FASTQ records found in output file");

    // Validate each record has the expected components
    for record in records {
        // Check sequence name exists
        assert!(!record.id().is_empty(), "FASTQ record missing sequence ID");
        
        // Check sequence exists and is not empty
        assert!(!record.seq().is_empty(), "FASTQ record has empty sequence");
        
        // Check quality scores exist and match sequence length
        assert_eq!(record.seq().len(), record.qual().len(), 
            "Quality score length does not match sequence length");
    }
}

/// Validates that a FASTA file exists and contains valid records
fn validate_fasta_output(path: &PathBuf) {
    let reader = bio::io::fasta::Reader::from_file(path).expect("Failed to open output file");
    let mut records = Vec::new();

    for result in reader.records() {
        // Each record should be valid and parseable
        let record = result.expect("Error reading FASTA record");
        records.push(record);
    }

    // Basic validation that we got some records
    assert!(!records.is_empty(), "No FASTA records found in output file");

    // Validate each record has the expected components
    for record in records {
        // Check sequence name exists
        assert!(!record.id().is_empty(), "FASTA record missing sequence ID");
        
        // Check sequence exists and is not empty
        assert!(!record.seq().is_empty(), "FASTA record has empty sequence");
    }
}

/// Converts a FASTQ file to FASTA format
fn convert_fastq_to_fasta(fastq_path: &PathBuf, fasta_path: &PathBuf) {
    let reader = bio::io::fastq::Reader::from_file(fastq_path).expect("Failed to open FASTQ file");
    let mut writer = bio::io::fasta::Writer::new(BufWriter::new(File::create(fasta_path).expect("Failed to create FASTA file")));

    for result in reader.records() {
        let fastq_record = result.expect("Error reading FASTQ record");
        let fasta_record = bio::io::fasta::Record::with_attrs(
            fastq_record.id(),
            fastq_record.desc(),
            fastq_record.seq(),
        );
        writer.write_record(&fasta_record).expect("Failed to write FASTA record");
    }
}

#[test]
fn test_rescue_basic() {
    // Create a temporary directory for test files
    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("output.fastq");
    
    // Test files
    let test_bam = "gs://fc-8c3900db-633f-477f-96b3-fb31ae265c44/results/PBFlowcell/m84060_230907_210011_s2/reads/ccs/aligned/m84060_230907_210011_s2.bam";
    let test_cram = "gs://fc-1ee08173-e353-4494-ad28-7a3d7bd99734/working/HPRC/HG00438/raw_data/Illumina/child/HG00438.final.cram";
    
    // First fetch some reads to use as input for rescue
    let fetch_output = temp_dir.path().join("fetch_output.fastq");
    let fetch_cmd = Command::new("cargo")
        .args([
            "run",
            "--",
            "fetch",
            "--loci",
            "chr15:23960193-23963918",
            "--output",
            fetch_output.to_str().unwrap(),
            test_bam,
        ])
        .output()
        .expect("Failed to execute fetch command");

    assert!(fetch_cmd.status.success(), "Fetch command failed: {:?}", fetch_cmd);
    assert!(fetch_output.exists(), "Fetch output file was not created");

    // Now run the rescue command
    let output = Command::new("cargo")
        .args([
            "run",
            "--",
            "rescue",
            "--ref-path",
            "/Users/kiran/repositories/hidive/scratch/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa",
            "--fastx-paths",
            fetch_output.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
            test_cram,
        ])
        .output()
        .expect("Failed to execute rescue command");

    // Check if the command was successful
    assert!(output.status.success(), "Rescue command failed: {:?}", output);
    
    // Verify the output file exists
    assert!(output_path.exists(), "Output file was not created");

    // Validate the format of the output FASTQ file
    validate_fastq_output(&output_path);
}

#[test]
fn test_rescue_basic_fasta() {
    // Create a temporary directory for test files
    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("output.fastq");
    
    // Test files
    let test_bam = "gs://fc-8c3900db-633f-477f-96b3-fb31ae265c44/results/PBFlowcell/m84060_230907_210011_s2/reads/ccs/aligned/m84060_230907_210011_s2.bam";
    let test_cram = "gs://fc-1ee08173-e353-4494-ad28-7a3d7bd99734/working/HPRC/HG00438/raw_data/Illumina/child/HG00438.final.cram";
    
    // First fetch some reads to use as input for rescue
    let fetch_output = temp_dir.path().join("fetch_output.fastq");
    let fetch_cmd = Command::new("cargo")
        .args([
            "run",
            "--",
            "fetch",
            "--loci",
            "chr15:23960193-23963918",
            "--output",
            fetch_output.to_str().unwrap(),
            test_bam,
        ])
        .output()
        .expect("Failed to execute fetch command");

    assert!(fetch_cmd.status.success(), "Fetch command failed: {:?}", fetch_cmd);
    assert!(fetch_output.exists(), "Fetch output file was not created");

    // Convert FASTQ to FASTA for testing
    let fasta_input = temp_dir.path().join("input.fasta");
    convert_fastq_to_fasta(&fetch_output, &fasta_input);
    assert!(fasta_input.exists(), "FASTA conversion failed");
    
    // Validate the FASTA file was created correctly
    validate_fasta_output(&fasta_input);

    // Now run the rescue command with FASTA input
    let output = Command::new("cargo")
        .args([
            "run",
            "--",
            "rescue",
            "--ref-path",
            "/Users/kiran/repositories/hidive/scratch/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa",
            "--fastx-paths",
            fasta_input.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
            test_cram,
        ])
        .output()
        .expect("Failed to execute rescue command");

    // Check if the command was successful
    assert!(output.status.success(), "Rescue command failed: {:?}", output);
    
    // Verify the output file exists
    assert!(output_path.exists(), "Output file was not created");

    // Validate the format of the output FASTQ file
    validate_fastq_output(&output_path);
}

/*
#[test]
fn test_rescue_with_search_option() {
    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("output.fastq");
    
    // Test files
    let test_bam = "gs://fc-8c3900db-633f-477f-96b3-fb31ae265c44/results/PBFlowcell/m84060_230907_210011_s2/reads/ccs/aligned/m84060_230907_210011_s2.bam";
    let test_cram = "gs://fc-1ee08173-e353-4494-ad28-7a3d7bd99734/working/HPRC/HG00438/raw_data/Illumina/child/HG00438.final.cram";
    let ref_fasta = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta";
    
    // First fetch some reads
    let fetch_output = temp_dir.path().join("fetch_output.fastq");
    let fetch_cmd = Command::new("cargo")
        .args([
            "run",
            "--",
            "fetch"
            "--loci",
            "chr15:23960193-23963918",
            "--output",
            fetch_output.to_str().unwrap(),
            test_bam,
        ])
        .output()
        .expect("Failed to execute fetch command");

    assert!(fetch_cmd.status.success(), "Fetch command failed: {:?}", fetch_cmd);
    assert!(fetch_output.exists(), "Fetch output file was not created");

    // Run rescue with contig-and-interval search option
    let output = Command::new("cargo")
        .args([
            "run",
            "--",
            "rescue",
            "--output",
            output_path.to_str().unwrap(),
            "--ref",
            ref_fasta,
            "--fastx",
            fetch_output.to_str().unwrap(),
            "--search-option",
            "contig-and-interval",
            test_cram,
        ])
        .output()
        .expect("Failed to execute rescue command");

    assert!(output.status.success(), "Rescue command failed: {:?}", output);
    assert!(output_path.exists(), "Output file was not created");
    validate_fastq_output(&output_path);
}
*/

#[test]
fn test_rescue_missing_required_args() {
    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("output.fastq");
    
    let output = Command::new("cargo")
        .args([
            "run",
            "--",
            "rescue",
            "--output",
            output_path.to_str().unwrap(),
            // Missing required --ref, --fasta arguments and input file
        ])
        .output()
        .expect("Failed to execute command");

    // This should fail because required arguments are missing
    assert!(!output.status.success(), "Command should have failed");
} 