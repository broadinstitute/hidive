version 1.0

workflow SubsetAndAlignBam {
    input {
        File bam
        File reference

        String? to_locus
        File? to_loci
    }

    if (defined(to_locus)) { call PrepareLocus as PrepareToLocus { input: locus = select_first([to_locus]) } }
    File to_bed = select_first([PrepareToLocus.bed, to_loci])

    call SubsetAndAlign {
        input:
            bam = bam,
            ref_fasta = reference,
            to_bed = to_bed
    }

    output {
        File pb_aligned_bam = SubsetAndAlign.aligned_bam
        File pb_aligned_bai = SubsetAndAlign.aligned_bai
    }
}

task PrepareLocus {
    input {
        String locus
    }

    command <<<
        set -euxo pipefail

        echo "~{locus}" | \
            sed 's/,//g' | \
            sed 's/[:\|]/\t/g' | \
            sed 's/-/\t/' > locus.bed
    >>>

    output {
        File bed = "locus.bed"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:0.1.101"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 1 SSD"
    }
}

task SubsetAndAlign {
    input {
        File bam
        File ref_fasta

        File to_bed

        String prefix = "out"

        Int num_cpus = 16 
    }

    Int disk_size_gb = 1 + 4*ceil(size([bam, ref_fasta, to_bed], "GB"))
    Int memory_gb = 6*num_cpus

    command <<<
        set -euxo pipefail

        samtools view -b -L ~{to_bed} ~{bam} | \
            samtools fastq - | \
            minimap2 -ayYL --eqx -x map-hifi ~{ref_fasta} - | \
            samtools sort --write-index -O BAM -o ~{prefix}.bam
    >>>

    output {
        File aligned_bam = "~{prefix}.bam"
        File aligned_bai = "~{prefix}.bam.csi"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_eval"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}