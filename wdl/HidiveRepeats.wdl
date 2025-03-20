version 1.0

workflow HidiveRepeats {
    input {
        File long_reads_bam
        File long_reads_bai
        File short_reads_cram
        File short_reads_crai

        String? locus
        File? loci

        Int low_k
        File low_k_model

        Int high_k
        File high_k_model

        String sample_name

        File reference
        File ref_fa_with_alt
        File ref_fai_with_alt
        File ref_cache_tar_gz

        Int padding = 500

        File? mat_aln_bam
        File? pat_aln_bam
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
            k = low_k,
            model = low_k_model,
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

    call Consensus as Consensus1 {
        input:
            loci = bed,
            reference = reference,
            aligned_reads_bam = Phase.hap1_bam,
            aligned_reads_csi = Phase.hap1_bai,
            prefix = sample_name + ".hap1"
    }

    call Correct as Correct1 {
        input:
            loci = bed,
            k = high_k,
            model = high_k_model,
            reference = reference,
            long_reads_bam = Consensus1.consensus_bam,
            short_read_fasta = Rescue.fasta,
            prefix = sample_name + ".hap1"
    }

    call Consensus as Consensus2 {
        input:
            loci = bed,
            reference = reference,
            aligned_reads_bam = Phase.hap2_bam,
            aligned_reads_csi = Phase.hap2_bai,
            prefix = sample_name + ".hap2"
    }

    call Correct as Correct2 {
        input:
            loci = bed,
            k = high_k,
            model = high_k_model,
            reference = reference,
            long_reads_bam = Consensus2.consensus_bam,
            short_read_fasta = Rescue.fasta,
            prefix = sample_name + ".hap2"
    }

    if (defined(mat_aln_bam)) {
        call Fetch as FetchMat {
            input:
                bam = select_first([mat_aln_bam]),
                loci = bed,
                padding = padding,
                prefix = sample_name
        }

        call Align as AlignMat {
            input:
                fasta = FetchMat.fasta,
                reference = reference,
                prefix = sample_name
        }
    }

    if (defined(pat_aln_bam)) {
        call Fetch as FetchPat {
            input:
                bam = select_first([pat_aln_bam]),
                loci = bed,
                padding = padding,
                prefix = sample_name
        }

        call Align as AlignPat {
            input:
                fasta = FetchPat.fasta,
                reference = reference,
                prefix = sample_name
        }
    }

    output {
        File repeats_hap1_bam = Correct1.corrected_bam
        File repeats_hap1_csi = Correct1.corrected_csi

        File repeats_hap2_bam = Correct2.corrected_bam
        File repeats_hap2_csi = Correct2.corrected_csi

        File? repeats_mat_bam = AlignMat.cluster_bam
        File? repeats_mat_csi = AlignMat.cluster_csi

        File? repeats_pat_bam = AlignPat.cluster_bam
        File? repeats_pat_csi = AlignPat.cluster_csi
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
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:0.1.101"
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

        hidive fetch -l ~{select_first([locus, loci])} -p ~{padding} ~{bam} > ~{prefix}.fa
    >>>

    output {
        File fasta = "~{prefix}.fa"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:0.1.107"
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
        Int k
        File model
        File reference

        String prefix = "out"

        Int num_cpus = 8
    }

    Int disk_size_gb = 1 + 2*ceil(size([model, reference, short_read_fasta], "GB"))
    Int memory_gb = 2*num_cpus

    command <<<
        set -euxo pipefail

        hidive correct -l ~{loci} -k ~{k} -m ~{model} ~{long_reads_bam} ~{short_read_fasta} | \
            minimap2 -ayYL -x map-hifi ~{reference} - | \
            samtools sort --write-index -O BAM -o ~{prefix}.bam
    >>>

    output {
        File corrected_bam = "~{prefix}.bam"
        File corrected_csi = "~{prefix}.bam.csi"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_dep_fix"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
        maxRetries: 2
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
        set -euxo pipefail

        hidive phase -s "~{sample_name}" -l ~{loci} -r ~{reference} -o ~{prefix} ~{aligned_reads_bam}

        samtools index ~{prefix}.phased1.bam
        samtools index ~{prefix}.phased2.bam
    >>>

    output {
        File hap1_bam = "~{prefix}.phased1.bam"
        File hap1_bai = "~{prefix}.phased1.bam.bai"

        File hap2_bam = "~{prefix}.phased2.bam"
        File hap2_bai = "~{prefix}.phased2.bam.bai"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_call_star_alleles"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}

task Consensus {
    input {
        File loci

        File reference

        File aligned_reads_bam
        File aligned_reads_csi

        String prefix = "out"

        Int num_cpus = 8
    }

    Int disk_size_gb = 1 + 2*ceil(size([reference, aligned_reads_bam, aligned_reads_csi], "GB"))
    Int memory_gb = 2*num_cpus

    command <<<
        set -x

        hidive consensus -l ~{loci} -r ~{reference} -o ~{prefix}.fa ~{aligned_reads_bam}

        minimap2 -ayYL -x map-hifi ~{reference} ~{prefix}.fa | \
            samtools sort --write-index -O BAM -o ~{prefix}.bam

        ls -lah *
    >>>

    output {
        File consensus_bam = "~{prefix}.bam"
        File consensus_csi = "~{prefix}.bam.csi"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_call_star_alleles"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}

task BamToFasta {
    input {
        File bam

        String prefix
    }

    Int disk_size_gb = 1 + 2*ceil(size([bam], "GB"))
    Int num_cpus = 1
    Int memory_gb = 2

    command <<<
        set -euxo pipefail

        samtools fasta ~{bam} > ~{prefix}.fa
    >>>

    output {
        File fasta = "~{prefix}.fa"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_eval"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}

task Align {
    input {
        File reference
        File fasta

        String prefix
    }

    Int disk_size_gb = 1 + 2*ceil(size([reference, fasta], "GB"))
    Int num_cpus = 1
    Int memory_gb = 2

    command <<<
        set -euxo pipefail

        minimap2 -ayYL --eqx -x asm20 -R '@RG\tID:~{prefix}\tSM:~{prefix}' ~{reference} ~{fasta} > out.sam
        samtools sort --write-index -O BAM -o ~{prefix}.bam out.sam
    >>>

    output {
        File cluster_bam = "~{prefix}.bam"
        File cluster_csi = "~{prefix}.bam.csi"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_eval"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}

task Plot {
    input {
        File cluster_fa
        String prefix
    }

    Int disk_size_gb = 1 + 2*ceil(size([cluster_fa], "GB"))
    Int num_cpus = 1
    Int memory_gb = 2

    command <<<
        set -euxo pipefail

        python3 /scripts/divide_fasta_into_kmers.py ~{cluster_fa} hidive_clusters_hprc_kmers.fa
        minimap2 -a -Y -x map-hifi -A 2 -B 4 /scripts/CYP2D6_and_7_T2T_w_grch38_rep_and_spacer.fa hidive_clusters_hprc_kmers.fa > kmer_align.sam
        python3 /scripts/make_kmer_plot.py kmer_align.sam /scripts/cyp2d6_plot_opts.json 10 CYP2D cyp2d6_plot.~{prefix}.html
    >>>

    output {
        File cluster_html = "cyp2d6_plot.~{prefix}.html"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_dep_fix"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}

task MergeAlignments {
    input {
        Array[File] bams
        String prefix
    }

    Int disk_size_gb = 1 + 2*ceil(size(bams, "GB"))
    Int num_cpus = 4
    Int memory_gb = 8 

    command <<<
        set -euxo pipefail

        samtools merge --write-index -O BAM -o ~{prefix}.merged.bam ~{sep=' ' bams}
    >>>

    output {
        File merged_bam = "~{prefix}.merged.bam"
        File merged_csi = "~{prefix}.merged.bam.csi"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_dep_fix"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}

task Type {
    input {
        File cluster_fa

        File cyp2d6_catalog
        File cyp2d6_haplotypes

        String prefix
    }

    Int disk_size_gb = 1 + 2*ceil(size([cluster_fa, cyp2d6_catalog, cyp2d6_haplotypes], "GB"))
    Int num_cpus = 4
    Int memory_gb = 8 

    command <<<
        set -euxo pipefail

        unzip ~{cyp2d6_catalog}
        CATALOG_DIR=$(basename ~{cyp2d6_catalog} .zip)

        minimap2 -x asm20 -a -Y -P --MD --eqx --sam-hit-only ~{cluster_fa} ~{cyp2d6_haplotypes} > star_alleles_aligned_haps.sam
        python /scripts/get_star_allele_site_matches.py star_alleles_aligned_haps.sam $CATALOG_DIR/GRCh38 site_matches.tsv
        python /scripts/get_final_star_alleles.py star_alleles_aligned_haps.sam site_matches.tsv $CATALOG_DIR/GRCh38 ~{prefix}.star_alleles.tsv
    >>>

    output {
        File star_alleles = "~{prefix}.star_alleles.tsv"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_call_star_alleles"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}