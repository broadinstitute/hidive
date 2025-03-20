version 1.0

workflow HidiveSummarize {
    input {
        Array[File] hap1_bams
        Array[File] hap2_bams
    }

    call CombineContigs as CombineHTT { input: hap1_bams = hap1_bams, hap2_bams = hap2_bams, gene_name = "HTT" }

    call SummarizeHTT { input: fasta = CombineHTT.fasta }
    
    output {
        File report_html = SummarizeHTT.report_html
        File legend_html = SummarizeHTT.legend_html
    }
}

task CombineContigs {
    input {
        Array[File] hap1_bams
        Array[File] hap2_bams

        String gene_name
    }

    Int disk_size_gb = 1 + 4*ceil(size(hap1_bams, "GB"))

    command <<<
        python3 <<CODE
import pysam

# Iterate over both hap1 and hap2 bams together
hap1_files = "~{sep=',' hap1_bams}".split(",")
hap2_files = "~{sep=',' hap2_bams}".split(",")

print("Hap1 files:", hap1_files)
print("Hap2 files:", hap2_files)

for hap1_bam, hap2_bam in zip(hap1_files, hap2_files): 
    # Extract sample name from hap1 bam filename
    sample_name = hap1_bam.split("/")[-1].split(".")[0]
    print(f"Processing sample: {sample_name}")

    # Open both BAM files and output fasta file
    with open("out.fasta", "a") as fasta_out:
        with pysam.AlignmentFile(hap1_bam, "rb") as h1, pysam.AlignmentFile(hap2_bam, "rb") as h2:
            # Process hap1 reads
            for read in h1:
                if read.query_name.startswith("~{gene_name}"):
                    fasta_out.write(f">{sample_name}_hap1\n{read.query_sequence}\n")
                    
            # Process hap2 reads  
            for read in h2:
                if read.query_name.startswith("~{gene_name}"):
                    fasta_out.write(f">{sample_name}_hap2\n{read.query_sequence}\n")
CODE
    >>>

    output {
        File fasta = "out.fasta"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:0.1.107"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}

task SummarizeHTT {
    input {
        File fasta
    }

    Int disk_size_gb = 1 + 2*ceil(size(fasta, "GB"))

    command <<<
        set -euxo pipefail

        mv ~{fasta} htt.fa

        python /scripts/make_viz.py /scripts/htt_params.json
    >>>

    output {
        File report_html = "HTT_viz.html"
        File legend_html = "HTT_viz_legend.html"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_call_star_alleles"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}