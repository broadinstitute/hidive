use std::path::PathBuf;
use std::process::Command;
use tempfile::TempDir;

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

#[test]
fn test_fetch_basic() {
    // Create a temporary directory for test files
    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("output.fastq");
    
    // Test BAM file
    let test_bam = "gs://fc-8c3900db-633f-477f-96b3-fb31ae265c44/results/PBFlowcell/m84060_230907_210011_s2/reads/ccs/aligned/m84060_230907_210011_s2.bam";
    
    // Run the fetch command
    let output = Command::new("cargo")
        .args([
            "run",
            "--",
            "fetch",
            "--loci",
            "chr15:23960193-23963918",
            "--output",
            output_path.to_str().unwrap(),
            test_bam,
        ])
        .output()
        .expect("Failed to execute command");

    // Check if the command was successful
    assert!(output.status.success(), "Command failed: {:?}", output);
    
    // Verify the output file exists
    assert!(output_path.exists(), "Output file was not created");

    // Validate the format of the output FASTQ file
    validate_fastq_output(&output_path);
}

#[test]
fn test_fetch_with_padding() {
    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("output.fastq");
    let test_bam = "gs://fc-8c3900db-633f-477f-96b3-fb31ae265c44/results/PBFlowcell/m84060_230907_210011_s2/reads/ccs/aligned/m84060_230907_210011_s2.bam";
    
    let output = Command::new("cargo")
        .args([
            "run",
            "--",
            "fetch",
            "--loci",
            "chr15:23960193-23963918",
            "--padding",
            "100",
            "--output",
            output_path.to_str().unwrap(),
            test_bam
        ])
        .output()
        .expect("Failed to execute command");

    assert!(output.status.success(), "Command failed: {:?}", output);
    assert!(output_path.exists(), "Output file was not created");
    validate_fastq_output(&output_path);
}

// #[test]
// fn test_fetch_with_haplotype() {
//     let temp_dir = TempDir::new().unwrap();
//     let output_path = temp_dir.path().join("output.fastq");
//     let test_bam = "gs://fc-8c3900db-633f-477f-96b3-fb31ae265c44/results/PBFlowcell/m84060_230907_210011_s2/reads/ccs/aligned/m84060_230907_210011_s2.bam";
    
//     let output = Command::new("cargo")
//         .args([
//             "run",
//             "--",
//             "fetch",
//             "--loci",
//             "chr15:23960193-23963918",
//             "--haplotype",
//             "1",
//             "--output",
//             output_path.to_str().unwrap(),
//             test_bam,
//         ])
//         .output()
//         .expect("Failed to execute command");

//     assert!(output.status.success(), "Command failed: {:?}", output);
//     assert!(output_path.exists(), "Output file was not created");
//     validate_fastq_output(&output_path);
// }

#[test]
fn test_fetch_invalid_locus() {
    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("output.fastq");
    let test_bam = "gs://fc-8c3900db-633f-477f-96b3-fb31ae265c44/results/PBFlowcell/m84060_230907_210011_s2/reads/ccs/aligned/m84060_230907_210011_s2.bam";
    
    let output = Command::new("cargo")
        .args([
            "run",
            "--",
            "fetch",
            "--loci",
            "invalid_locus",
            "--output",
            output_path.to_str().unwrap(),
            test_bam,
        ])
        .output()
        .expect("Failed to execute command");

    // This should fail because the locus format is invalid
    assert!(!output.status.success(), "Command should have failed");
}

#[test]
fn test_fetch_missing_required_args() {
    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("output.fastq");
    
    let output = Command::new("cargo")
        .args([
            "run",
            "--",
            "fetch",
            "--output",
            output_path.to_str().unwrap(),
            // Missing required --loci argument and BAM file
        ])
        .output()
        .expect("Failed to execute command");

    // This should fail because required arguments are missing
    assert!(!output.status.success(), "Command should have failed");
} 