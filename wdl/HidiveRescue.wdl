version 1.0

workflow HidiveRescue {
    input {
        File long_reads_bam
        File long_reads_bai

        File short_reads_cram
        File short_reads_crai

        String? locus
        File? loci

        String sample_name

        File ref_fa_with_alt
        File ref_fai_with_alt
        File ref_cache_tar_gz

        Int padding = 500
    }

    if (defined(locus)) { call PrepareLocus { input: locus = select_first([locus]) } }
    File bed = select_first([PrepareLocus.bed, loci])

    call Fetch {
        input:
            bam = long_reads_bam,
            loci = bed,
            padding = padding,
            prefix = sample_name
    }

    call Rescue {
        input:
            long_reads_fastq = Fetch.fastq,
            short_reads_cram = short_reads_cram,
            short_reads_crai = short_reads_crai,
            ref_fa_with_alt = ref_fa_with_alt,
            ref_fai_with_alt = ref_fai_with_alt,
            ref_cache_tar_gz = ref_cache_tar_gz,
            prefix = sample_name
    }

    output {
        File rescued_fastq_gz = Rescue.fastq_gz
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
            sed 's/\|/\t/g' | \
            sed 's/\(.*\):/\1\t/' | \
            sed 's/\(.*\)-\(.*\)/\1\t\2/' > locus.bed
    >>>

    output {
        File bed = "locus.bed"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:0.1.113"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 1 SSD"
    }
}

task Fetch {
    input {
        String bam
        String? locus
        File? loci
        Int padding

        String prefix = "out"

        Int disk_size_gb = 2
        Int num_cpus = 4
    }

    command <<<
        set -euxo pipefail

        hidive fetch -l ~{select_first([locus, loci])} -p ~{padding} ~{bam} > ~{prefix}.fq
    >>>

    output {
        File fastq = "~{prefix}.fq"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:0.1.113"
        memory: "2 GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}

task Rescue {
    input {
        File long_reads_fastq
        File short_reads_cram
        File short_reads_crai

        File ref_fa_with_alt
        File ref_fai_with_alt
        File ref_cache_tar_gz

        String prefix = "out"

        Int num_cpus = 16
    }

    Int disk_size_gb = 1 + 2*ceil(size([long_reads_fastq, short_reads_cram, short_reads_crai, ref_fa_with_alt, ref_fai_with_alt, ref_cache_tar_gz], "GB"))
    Int memory_gb = 2*num_cpus

    command <<<
        set -euxo pipefail

        mv ~{ref_fa_with_alt} Homo_sapiens_assembly38.fasta
        mv ~{ref_fai_with_alt} Homo_sapiens_assembly38.fasta.fai 
        mv ~{ref_cache_tar_gz} Homo_sapiens_assembly38.ref_cache.tar.gz

        tar xzf Homo_sapiens_assembly38.ref_cache.tar.gz >/dev/null 2>&1

        export REF_PATH="$(pwd)/ref/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s"
        export REF_CACHE="$(pwd)/ref/cache/%2s/%2s/%s"

        hidive rescue -r Homo_sapiens_assembly38.fasta -f ~{long_reads_fastq} ~{short_reads_cram} | gzip > ~{prefix}.fq.gz
    >>>

    output {
        File fastq_gz = "~{prefix}.fq.gz"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:0.1.113"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
        maxRetries: 0
    }
}