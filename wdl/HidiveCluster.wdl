version 1.0

workflow HidiveCluster {
    input {
        File long_reads_bam
        File long_reads_bai
        File short_reads_cram
        File short_reads_crai

        String? locus
        File? loci

        File model
        String sample_name

        File reference
        File ref_fa_with_alt
        File ref_fai_with_alt
        File ref_cache_tar_gz

        Int padding = 500
    }

    if (defined(locus)) {
        call PrepareLocus { input: locus = select_first([locus]) }
    }

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
            long_reads_fasta = Fetch.fasta,
            short_reads_cram = short_reads_cram,
            short_reads_crai = short_reads_crai,
            ref_fa_with_alt = ref_fa_with_alt,
            ref_fai_with_alt = ref_fai_with_alt,
            ref_cache_tar_gz = ref_cache_tar_gz,
            prefix = sample_name
    }

    call Correct {
        input:
            loci = bed,
            model = model,
            reference = reference,
            long_reads_bam = long_reads_bam,
            short_read_fasta = Rescue.fasta,
            prefix = sample_name
    }

    call Phase {
        input:
            loci = bed,
            reference = reference,
            sample_name = sample_name,
            aligned_reads_bam = Correct.corrected_bam,
            aligned_reads_csi = Correct.corrected_csi,
            prefix = sample_name
    }

    call Cluster as Cluster1 {
        input:
            loci = bed,
            reference = reference,
            hap_bam = Phase.hap1_bam,
            prefix = sample_name + ".hap1"
    }

    call Cluster as Cluster2 {
        input:
            loci = bed,
            reference = reference,
            hap_bam = Phase.hap2_bam,
            prefix = sample_name + ".hap2"
    }

    output {
        File corrected_bam = Correct.corrected_bam
        File corrected_csi = Correct.corrected_csi

        File hap1_bam = Phase.hap1_bam
        File hap2_bam = Phase.hap2_bam

        File hap1_fa = Cluster1.cluster_fa
        File hap2_fa = Cluster2.cluster_fa
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

task Fetch {
    input {
        String bam
        File loci
        Int padding

        String prefix = "out"

        Int disk_size_gb = 2
        Int num_cpus = 4
    }

    command <<<
        set -euxo pipefail

        hidive fetch -l ~{loci} -p ~{padding} ~{bam} > ~{prefix}.fa
    >>>

    output {
        File fasta = "~{prefix}.fa"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:0.1.101"
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
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:0.1.101"
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

        File loci
        File model
        File reference

        String prefix = "out"

        Int num_cpus = 8
    }

    Int disk_size_gb = 1 + 2*ceil(size([model, reference, short_read_fasta], "GB"))
    Int memory_gb = 2*num_cpus

    command <<<
        set -x

        hidive correct -l ~{loci} -m ~{model} ~{long_reads_bam} ~{short_read_fasta} | \
            minimap2 -ayYL -x map-hifi ~{reference} - | \
            samtools sort --write-index -O BAM -o ~{prefix}.bam
    >>>

    output {
        File corrected_bam = "~{prefix}.bam"
        File corrected_csi = "~{prefix}.bam.csi"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:0.1.101"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}

task Phase {
    input {
        File loci
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

        hidive phase -s "~{sample_name}" -l ~{loci} -r ~{reference} -o ~{prefix} ~{aligned_reads_bam}
    >>>

    output {
        File hap1_bam = "~{prefix}.phased1.bam"
        File hap2_bam = "~{prefix}.phased2.bam"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_eval"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}

task Cluster {
    input {
        File reference
        File hap_bam

        File loci
        String prefix

        Int num_cpus = 8
    }

    Int disk_size_gb = 1 + 2*ceil(size([reference, hap_bam], "GB"))
    Int memory_gb = 2*num_cpus

    command <<<
        set -x

        samtools index ~{hap_bam}

        hidive cluster -f ~{loci} -t ~{loci} -r ~{reference} -o ~{prefix} ~{hap_bam} > ~{prefix}.fa
    >>>

    output {
        File cluster_fa = "~{prefix}.fa"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_eval"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}