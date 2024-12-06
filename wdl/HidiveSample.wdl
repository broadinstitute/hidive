version 1.0

workflow HidiveSample {
    input {
        File long_reads_bam
        File long_reads_bai
        File short_reads_cram
        File short_reads_crai

        File reference
        String locus
        File model
        String sample_name

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
            locus = locus,
            model = model,
            reference = reference,
            long_reads_bam = long_reads_bam,
            short_read_fasta = Rescue.fasta,
    }

    call Call {
        input:
            locus = locus,
            reference = reference,
            sample_name = sample_name,
            aligned_reads_bam = Correct.corrected_bam,
            aligned_reads_csi = Correct.corrected_csi,
    }

    call Consensus {
        input:
            locus = locus,
            reference = reference,
            calls_vcf = Call.calls_vcf,
            calls_tbi = Call.calls_tbi,
    }

    output {
        File corrected_bam = Correct.corrected_bam
        File corrected_csi = Correct.corrected_csi

        File calls_vcf = Call.calls_vcf
        File calls_tbi = Call.calls_tbi

        File h1_bam = Consensus.h1_bam
        File h1_csi = Consensus.h1_csi
        File h2_bam = Consensus.h2_bam
        File h2_csi = Consensus.h2_csi
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
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:0.1.98"
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

        Int num_cpus = 16
    }

    Int disk_size_gb = 1 + 2*ceil(size([long_reads_fasta, short_reads_cram, short_reads_crai, ref_fa_with_alt, ref_fai_with_alt, ref_cache_tar_gz], "GB"))
    Int memory_gb = 2*num_cpus

    command <<<
        set -euxo pipefail

        mv ~{ref_fa_with_alt} Homo_sapiens_assembly38.fasta
        mv ~{ref_fai_with_alt} Homo_sapiens_assembly38.fasta.fai 
        mv ~{ref_cache_tar_gz} Homo_sapiens_assembly38.ref_cache.tar.gz

        tar xzf Homo_sapiens_assembly38.ref_cache.tar.gz >/dev/null 2>&1

        export REF_PATH="$(pwd)/ref/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s"
        export REF_CACHE="$(pwd)/ref/cache/%2s/%2s/%s"

        hidive rescue -r Homo_sapiens_assembly38.fasta -f ~{long_reads_fasta} ~{short_reads_cram} > ~{prefix}.fa
    >>>

    output {
        File fasta = "~{prefix}.fa"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_call"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
        maxRetries: 2
    }
}

task Correct {
    input {
        File short_read_fasta
        String long_reads_bam

        String locus
        File model
        File reference

        String prefix = "out"

        Int num_cpus = 8
    }

    Int disk_size_gb = 1 + 2*ceil(size([model, reference, short_read_fasta], "GB"))
    Int memory_gb = 2*num_cpus

    command <<<
        set -x

        hidive correct -l "~{locus}" -m ~{model} ~{long_reads_bam} ~{short_read_fasta} | \
            minimap2 -ayYL -x map-hifi ~{reference} - | \
            samtools sort --write-index -O BAM -o ~{prefix}.bam
    >>>

    output {
        File corrected_bam = "~{prefix}.bam"
        File corrected_csi = "~{prefix}.bam.csi"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:0.1.98"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}

task Call {
    input {
        String locus
        File reference
        String sample_name

        File aligned_reads_bam
        File aligned_reads_csi

        String prefix = "out"

        Int num_cpus = 8
    }

    Int disk_size_gb = 1 + 2*ceil(size([reference, aligned_reads_bam, aligned_reads_csi], "GB"))
    Int memory_gb = 2*num_cpus

    command <<<
        set -x

        hidive call -s "~{sample_name}" -l "~{locus}" -r ~{reference} ~{aligned_reads_bam} | bgzip -c > ~{prefix}.vcf.bgz
        tabix -p vcf ~{prefix}.vcf.bgz
    >>>

    output {
        File calls_vcf = "~{prefix}.vcf.bgz"
        File calls_tbi = "~{prefix}.vcf.bgz.tbi"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_call"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}

task Consensus {
    input {
        File reference
        File calls_vcf
        File calls_tbi

        String locus

        Int num_cpus = 8
    }

    Int disk_size_gb = 1 + 2*ceil(size([reference, calls_vcf, calls_tbi], "GB"))
    Int memory_gb = 2*num_cpus

    command <<<
        set -x

        LOCUS=$(echo "~{locus}" | sed 's/,//g' | sed 's/|.*//')

        samtools faidx ~{reference} ${LOCUS} | bcftools consensus -H 1pIu ~{calls_vcf} | minimap2 -ayYL -x map-hifi ~{reference} - | samtools sort --write-index -O BAM -o h1.bam
        samtools faidx ~{reference} ${LOCUS} | bcftools consensus -H 2pIu ~{calls_vcf} | minimap2 -ayYL -x map-hifi ~{reference} - | samtools sort --write-index -O BAM -o h2.bam
    >>>

    output {
        File h1_bam = "h1.bam"
        File h1_csi = "h1.bam.csi"
        File h2_bam = "h2.bam"
        File h2_csi = "h2.bam.csi"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_call"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}