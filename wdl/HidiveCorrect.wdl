version 1.0

workflow Hidive {
    input {
        File long_reads_bam
        File long_reads_bai
        File short_reads_cram
        File short_reads_crai

        File reference
        String locus
        File model

        Int padding = 500
    }

    call Fetch {
        input:
            bam = long_reads_bam,
            locus = locus,
            padding = padding,
    }

    call Rescue {
        input:
            long_reads_fasta = Fetch.fasta,
            short_reads_cram = short_reads_cram,
            short_reads_crai = short_reads_crai,
    }

    call Correct {
        input:
            model = model,
            long_read_fasta = Fetch.fasta,
            short_read_fasta = Rescue.fasta,
    }

    call Align {
        input:
            reference = reference,
            sequences = Correct.fasta,
    }

    output {
        File corrected_fa = Correct.fasta
        File aligned_bam = Align.aligned_bam
        File aligned_bai = Align.aligned_bai
    }
}

task Fetch {
    input {
        String bam
        String locus
        Int padding

        String prefix = "out"

        Int disk_size_gb = 2
        Int num_cpus = 4
    }

    command <<<
        set -euxo pipefail

        hidive fetch -l "~{locus}" -p ~{padding} ~{bam} > ~{prefix}.fa
    >>>

    output {
        File fasta = "~{prefix}.fa"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:0.1.77"
        memory: "2 GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}

task Rescue {
    input {
        File long_reads_fasta
        File short_reads_cram
        File short_reads_crai

        File ref_fa_with_alt
        File ref_fai_with_alt

        String prefix = "out"

        Int num_cpus = 4
    }

    Int disk_size_gb = 1 + 2*ceil(size([long_reads_fasta, short_reads_cram, short_reads_crai], "GB"))

    command <<<
        set -euxo pipefail

        export REF_PATH="$(dirname ~{ref_fa_with_alt})"

        hidive rescue -f ~{long_reads_fasta} ~{short_reads_cram} > ~{prefix}.fa
    >>>

    output {
        File fasta = "~{prefix}.fa"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:0.1.77"
        memory: "2 GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}

task Correct {
    input {
        File model
        File long_read_fasta
        File short_read_fasta

        String prefix = "out"

        Int num_cpus = 4
    }

    Int disk_size_gb = 1 + 2*ceil(size([model, long_read_fasta, short_read_fasta], "GB"))

    command <<<
        set -euxo pipefail

        hidive correct -m ~{model} -s ~{short_read_fasta} ~{long_read_fasta} > ~{prefix}.fa
    >>>

    output {
        File fasta = "~{prefix}.fa"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:latest"
        memory: "4 GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}

task Align {
    input {
        File reference
        File sequences

        String preset = "map-hifi"
        String prefix = "out"

        Int num_cpus = 8
    }

    command <<<
        set -euxo pipefail

        minimap2 -t ~{num_cpus} -ayYL -x ~{preset} ~{reference} ~{sequences} | samtools sort --write-index -O BAM -o ~{prefix}.bam
    >>>

    output {
        File aligned_bam = "~{prefix}.bam"
        File aligned_bai = "~{prefix}.bam.bai"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:latest"
        memory: "32 GB"
        cpu: num_cpus
        disks: "local-disk 100 SSD"
    }
}