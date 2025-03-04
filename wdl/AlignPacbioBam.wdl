version 1.0

workflow AlignPacbioBam {
    input {
        File bam
        File reference
    }

    call Align {
        input:
            bam = bam,
            ref_fasta = reference
    }

    output {
        File pb_aligned_bam = Align.aligned_bam
        File pb_aligned_bai = Align.aligned_bai
    }
}

task Align {
    input {
        File bam
        File ref_fasta

        String sample_name
        String prefix = "out"

        Int num_cpus = 16 
    }

    Int disk_size_gb = 1 + 4*ceil(size([bam, ref_fasta], "GB"))
    Int memory_gb = 6*num_cpus

    command <<<
        set -euxo pipefail

        pbmm2 align ~{bam} ~{ref_fasta} ~{prefix}.bam \
            --sample ~{sample_name} \
            --sort \
            --unmapped
    >>>

    output {
        File aligned_bam = "~{prefix}.bam"
        File aligned_bai = "~{prefix}.bam.bai"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-smrttools:11.0.0.146107"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}