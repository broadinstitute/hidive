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
            aligned_reads_bam = Correct.corrected_bam,
            aligned_reads_csi = Correct.corrected_csi,
    }

    call Phase {
        input:
            reference = reference,
            called_bam = Call.called_bam,
            called_csi = Call.called_csi,
            sample_name = sample_name,
    }

    output {
        File corrected_bam = Correct.corrected_bam
        File corrected_csi = Correct.corrected_csi
        File called_bam = Call.called_bam
        File called_csi = Call.called_csi
        File phased_gvcf = Phase.phased_gvcf
        File phased_gvcf_tbi = Phase.phased_gvcf_tbi
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
        # maxRetries: 1
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
        File reference
        File aligned_reads_bam
        File aligned_reads_csi

        String locus
        String prefix = "out"

        Int num_cpus = 8
    }

    Int disk_size_gb = 1 + 2*ceil(size([reference, aligned_reads_bam, aligned_reads_csi], "GB"))
    Int memory_gb = 2*num_cpus

    command <<<
        set -x

        hidive call -l "~{locus}" -r ~{reference} ~{aligned_reads_bam} | \
            minimap2 -ayYL -x map-hifi ~{reference} - | \
            samtools sort --write-index -O BAM -o ~{prefix}.bam
    >>>

    output {
        File called_bam = "~{prefix}.bam"
        File called_csi = "~{prefix}.bam.csi"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:0.1.98"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}

task Phase {
    input {
        File called_bam
        File called_csi
        File reference

        String sample_name
        String prefix = "out"

        Int num_cpus = 2
    }

    Int disk_size_gb = 1 + 2*ceil(size([reference, called_bam, called_csi], "GB"))
    Int memory_gb = 2*num_cpus

    command <<<
        set -x

        samtools view -h ~{called_bam} | grep -e '^@' -e '^h1' | samtools view -bS --write-index -o h1.bam
        samtools view -h ~{called_bam} | grep -e '^@' -e '^h2' | samtools view -bS --write-index -o h2.bam

        bcftools mpileup -f ~{reference} -g 0 -Ou -a DP -m 1 --indel-size 30000 h1.bam | bcftools call -m --ploidy 1 -g 0 -Oz -W -o h1.variants.gvcf.bgz
        bcftools mpileup -f ~{reference} -g 0 -Ou -a DP -m 1 --indel-size 30000 h2.bam | bcftools call -m --ploidy 1 -g 0 -Oz -W -o h2.variants.gvcf.bgz

        bcftools merge -g ~{reference} -0 -Ov h1.variants.gvcf.bgz h2.variants.gvcf.bgz | \
            awk -F'\t' '
            $0 ~ /^##/ {print; next}
            $0 ~ /^#CHROM/ {
                for(i=1; i<=9; i++) printf "%s\t", $i
                print "~{sample_name}"
                next
            }
            {
                h1 = ($10 ~ /^\./) ? "0" : substr($10,1,1)
                h2 = ($11 ~ /^\./) ? "0" : substr($11,1,1)
                for(i=1; i<=8; i++) printf "%s\t", $i
                printf "GT\t%s|%s\n", h1, h2
            }' > ~{prefix}.phased.gvcf

        bgzip -c ~{prefix}.phased.gvcf > ~{prefix}.phased.gvcf.bgz
        tabix -p vcf ~{prefix}.phased.gvcf.bgz
    >>>

    output {
        File phased_gvcf = "~{prefix}.phased.gvcf.bgz"
        File phased_gvcf_tbi = "~{prefix}.phased.gvcf.bgz.tbi"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_call"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}