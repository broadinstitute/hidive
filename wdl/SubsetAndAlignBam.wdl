version 1.0

workflow SubsetAndAlignBam {
    input {
        File bam
        File bai
        File reference

        String? from_locus
        File? from_loci

        String? to_locus
        File? to_loci
    }

    if (defined(from_locus)) { call PrepareLocus as PrepareFromLocus { input: locus = select_first([from_locus]) } }
    File from_bed = select_first([PrepareFromLocus.bed, from_loci])

    if (defined(to_locus)) { call PrepareLocus as PrepareToLocus { input: locus = select_first([to_locus]) } }
    File to_bed = select_first([PrepareToLocus.bed, to_loci])

    call SubsetAndAlign {
        input:
            bam = bam,
            bai = bai,
            ref_fasta = reference,
            from_bed = from_bed,
            to_bed = to_bed
    }

    output {
        File subset_fa = SubsetAndAlign.subset_fa
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
        File bai
        File ref_fasta

        File from_bed
        File to_bed

        String prefix = "subset"

        Int num_cpus = 16 
    }

    Int disk_size_gb = 1 + 4*ceil(size([bam, ref_fasta, to_bed], "GB"))
    Int memory_gb = 6*num_cpus

    command <<<
        set -euxo pipefail

        awk '{ print $1 ":" $2 "-" $3 }' ~{to_bed} | samtools faidx -r - ~{ref_fasta} > ~{prefix}.fa

        hidive fetch -l ~{from_bed} ~{bam} | \
            minimap2 -ayYL --eqx -x map-hifi ~{prefix}.fa - | \
            samtools sort --write-index -O BAM -o ~{prefix}.bam
    >>>

    output {
        File subset_fa = "~{prefix}.fa"
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