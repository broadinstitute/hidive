version 1.0

workflow HidiveCluster {
    input {
        File long_reads_bam
        File long_reads_bai
        File short_reads_cram
        File short_reads_crai

        String? from_locus
        File? from_loci

        String? to_locus
        File? to_loci

        File model
        String sample_name

        File reference
        File ref_fa_with_alt
        File ref_fai_with_alt
        File ref_cache_tar_gz

        File cyp2d6_catalog
        File cyp2d6_haplotypes

        Int padding = 500

        File? mat_aln_bam
        File? pat_aln_bam
    }

    if (defined(from_locus)) { call PrepareLocus as PrepareFromLocus { input: locus = select_first([from_locus]) } }
    File from_bed = select_first([PrepareFromLocus.bed, from_loci])

    if (defined(to_locus)) { call PrepareLocus as PrepareToLocus { input: locus = select_first([to_locus]) } }
    File to_bed = select_first([PrepareToLocus.bed, to_loci])

    call Fetch {
        input:
            bam = long_reads_bam,
            loci = from_bed,
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
            loci = from_bed,
            model = model,
            reference = reference,
            long_reads_bam = long_reads_bam,
            short_read_fasta = Rescue.fasta,
            prefix = sample_name
    }

    call Phase {
        input:
            loci = from_bed,
            reference = reference,
            sample_name = sample_name,
            aligned_reads_bam = Correct.corrected_bam,
            aligned_reads_csi = Correct.corrected_csi,
            prefix = sample_name
    }

    call Cluster as Cluster1 {
        input:
            sample_name = sample_name,
            from_loci = from_bed,
            to_loci = to_bed,
            reference = reference,
            hap_bam = Phase.hap1_bam,
            hap_bai = Phase.hap1_bai,
            prefix = sample_name + ".hap1"
    }

    call Cluster as Cluster2 {
        input:
            sample_name = sample_name,
            from_loci = from_bed,
            to_loci = to_bed,
            reference = reference,
            hap_bam = Phase.hap2_bam,
            hap_bai = Phase.hap2_bai,
            prefix = sample_name + ".hap2"
    }

    call SubsetReference {
        input:
            reference = reference,
            to_loci = to_bed
    }

    call BamToFasta as BamToFasta1 { input: bam = Phase.hap1_bam, prefix = sample_name + ".reads1" }
    call BamToFasta as BamToFasta2 { input: bam = Phase.hap2_bam, prefix = sample_name + ".reads2" }

    call Align as AlignReads1 {
        input:
            reference = SubsetReference.ref_subset_fa,
            fasta = BamToFasta1.fasta,
            prefix = sample_name + ".reads1"
    }

    call Align as AlignReads2 {
        input:
            reference = SubsetReference.ref_subset_fa,
            fasta = BamToFasta2.fasta,
            prefix = sample_name + ".reads2"
    }

    call Align as AlignCluster1 {
        input:
            reference = SubsetReference.ref_subset_fa,
            fasta = Cluster1.cluster_fa,
            prefix = sample_name + ".hap1"
    }

    call Align as AlignCluster2 {
        input:
            reference = SubsetReference.ref_subset_fa,
            fasta = Cluster2.cluster_fa,
            prefix = sample_name + ".hap2"
    }

    call Plot as PlotCluster1 {
        input:
            cluster_fa = Cluster1.cluster_fa,
            prefix = sample_name + ".hap1"
    }

    call Plot as PlotCluster2 {
        input:
            cluster_fa = Cluster2.cluster_fa,
            prefix = sample_name + ".hap2"
    }

    call Type as TypeCluster1 {
        input:
            cluster_fa = Cluster1.cluster_fa,
            cyp2d6_catalog = cyp2d6_catalog,
            cyp2d6_haplotypes = cyp2d6_haplotypes,
            prefix = sample_name + ".hap1",
    }

    call Type as TypeCluster2 {
        input:
            cluster_fa = Cluster2.cluster_fa,
            cyp2d6_catalog = cyp2d6_catalog,
            cyp2d6_haplotypes = cyp2d6_haplotypes,
            prefix = sample_name + ".hap2",
    }

    if (defined(mat_aln_bam) && defined(pat_aln_bam)) {
        call Fetch as FetchMat {
            input:
                bam = select_first([mat_aln_bam]),
                loci = from_bed,
                padding = padding,
                prefix = sample_name + ".mat"
        }

        call Align as AlignMatReads {
            input:
                reference = SubsetReference.ref_subset_fa,
                fasta = FetchMat.fasta,
                prefix = sample_name + ".mat"
        }

        call PrepareLocus as PrepareSubsetLocus { input: locus = "chr22:42121531-42135680:1-14150" }

        call Fetch as FetchMatSubset {
            input:
                bam = AlignMatReads.cluster_bam,
                loci = PrepareSubsetLocus.bed,
                padding = padding,
                prefix = sample_name + ".subset.mat"
        }

        call Type as TypeClusterMat {
            input:
                cluster_fa = FetchMatSubset.fasta,
                cyp2d6_catalog = cyp2d6_catalog,
                cyp2d6_haplotypes = cyp2d6_haplotypes,
                prefix = sample_name + ".subset.mat",
        }

        call Fetch as FetchPat {
            input:
                bam = select_first([pat_aln_bam]),
                loci = from_bed,
                padding = padding,
                prefix = sample_name + ".pat"
        }

        call Align as AlignPatReads {
            input:
                reference = SubsetReference.ref_subset_fa,
                fasta = FetchPat.fasta,
                prefix = sample_name + ".pat"
        }

        call Fetch as FetchPatSubset {
            input:
                bam = AlignPatReads.cluster_bam,
                loci = PrepareSubsetLocus.bed,
                padding = padding,
                prefix = sample_name + ".subset.pat"
        }

        call Type as TypeClusterPat {
            input:
                cluster_fa = FetchPatSubset.fasta,
                cyp2d6_catalog = cyp2d6_catalog,
                cyp2d6_haplotypes = cyp2d6_haplotypes,
                prefix = sample_name + ".subset.pat",
        }

        call MergeAlignments {
            input:
                bams = [
                    AlignMatReads.cluster_bam,
                    AlignPatReads.cluster_bam,
                    AlignCluster1.cluster_bam,
                    AlignCluster2.cluster_bam,
                    AlignReads1.cluster_bam,
                    AlignReads2.cluster_bam
                ],
                prefix = sample_name
        }
    }
    
    output {
        File corrected_bam = Correct.corrected_bam
        File corrected_csi = Correct.corrected_csi

        File hap1_bam = Phase.hap1_bam
        File hap1_bai = Phase.hap1_bai

        File hap2_bam = Phase.hap2_bam
        File hap2_bai = Phase.hap2_bai

        File ref_subset_fa = SubsetReference.ref_subset_fa
        File ref_subset_fai = SubsetReference.ref_subset_fai

        File cluster1_fa = Cluster1.cluster_fa
        File cluster1_bam = AlignCluster1.cluster_bam
        File cluster1_csi = AlignCluster1.cluster_csi

        File cluster2_fa = Cluster2.cluster_fa
        File cluster2_bam = AlignCluster2.cluster_bam
        File cluster2_csi = AlignCluster2.cluster_csi

        File cluster1_html = PlotCluster1.cluster_html
        File cluster2_html = PlotCluster2.cluster_html

        File cluster1_star_alleles = TypeCluster1.star_alleles
        File cluster2_star_alleles = TypeCluster2.star_alleles

        File? cluster_mat_star_alleles = TypeClusterMat.star_alleles
        File? cluster_pat_star_alleles = TypeClusterPat.star_alleles

        File? cluster_mat_bam = AlignMatReads.cluster_bam
        File? cluster_mat_csi = AlignMatReads.cluster_csi

        File? cluster_pat_bam = AlignPatReads.cluster_bam
        File? cluster_pat_csi = AlignPatReads.cluster_csi

        File? merged_bam = MergeAlignments.merged_bam
        File? merged_csi = MergeAlignments.merged_csi
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
        set -euxo pipefail

        hidive correct -l ~{loci} -m ~{model} ~{long_reads_bam} ~{short_read_fasta} | \
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
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_dep_fix"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}

task Cluster {
    input {
        File reference
        File hap_bam
        File hap_bai

        String sample_name

        File from_loci
        File to_loci
        String prefix

        Int num_cpus = 8
    }

    Int disk_size_gb = 1 + 2*ceil(size([reference, hap_bam, hap_bai, from_loci, to_loci], "GB"))
    Int memory_gb = 4*num_cpus

    command <<<
        set -euxo pipefail

        hidive cluster -s "~{sample_name}" -f ~{from_loci} -t ~{to_loci} -r ~{reference} -o ~{prefix} ~{hap_bam} > ~{prefix}.fa
    >>>

    output {
        File cluster_fa = "~{prefix}.fa"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_dep_fix"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}

task SubsetReference {
    input {
        File reference
        File to_loci
    }

    Int disk_size_gb = 1 + 2*ceil(size([reference, to_loci], "GB"))
    Int num_cpus = 1
    Int memory_gb = 2

    command <<<
        set -euxo pipefail

        awk '{ print $1 ":" $2 "-" $3 }' ~{to_loci} | samtools faidx -r - ~{reference} > ref.subset.fa
        samtools faidx ref.subset.fa
    >>>

    output {
        File ref_subset_fa = "ref.subset.fa"
        File ref_subset_fai = "ref.subset.fa.fai"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_eval"
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

        minimap2 -ayYL --eqx -x asm20 -R '@RG\tID:~{prefix}\tSM:~{prefix}' ~{reference} ~{fasta} | \
            samtools sort --write-index -O BAM -o ~{prefix}.bam
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