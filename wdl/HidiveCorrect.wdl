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
            long_reads_bam = long_reads_bam,
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
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_phase"
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
        File ref_cache_tar_gz

        String prefix = "out"
        String? contig

        Int num_cpus = 16
    }

    Int disk_size_gb = 1 + 2*ceil(size([long_reads_fasta, short_reads_cram, short_reads_crai, ref_fa_with_alt, ref_fai_with_alt, ref_cache_tar_gz], "GB"))

    command <<<
        set -euxo pipefail

        mv ~{ref_fa_with_alt} Homo_sapiens_assembly38.fasta
        mv ~{ref_fai_with_alt} Homo_sapiens_assembly38.fasta.fai 
        mv ~{ref_cache_tar_gz} Homo_sapiens_assembly38.ref_cache.tar.gz

        tar xzf Homo_sapiens_assembly38.ref_cache.tar.gz >/dev/null 2>&1

        export REF_PATH="$(pwd)/ref/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s"
        export REF_CACHE="$(pwd)/ref/cache/%2s/%2s/%s"

        hidive rescue \
            ~{true='-c' false='' defined(contig)} ~{select_first([contig, ""])} \
            -f ~{long_reads_fasta} \
            ~{short_reads_cram} \
            > ~{prefix}.fa
    >>>

    output {
        File fasta = "~{prefix}.fa"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_phase"
        memory: "~{num_cpus} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}

task Correct {
    input {
        File model
        String long_reads_bam
        File short_read_fasta

        String locus
        String prefix = "out"

        Int num_cpus = 4
    }

    Int disk_size_gb = 1 + 2*ceil(size([model, short_read_fasta], "GB"))

    command <<<
        set -euxo pipefail

        hidive correct -l "~{locus}" -m ~{model} ~{long_reads_bam} ~{short_read_fasta} > ~{prefix}.fa
    >>>

    output {
        File fasta = "~{prefix}.fa"
        File gfa = "~{prefix}.gfa"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_phase"
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
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_phase"
        memory: "32 GB"
        cpu: num_cpus
        disks: "local-disk 100 SSD"
    }
}