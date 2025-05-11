version 1.0

workflow HidiveMSA {
    input {
        File repeat_hap1_bam
        File repeat_hap2_bam

        File? hprc_hap1_bam
        File? hprc_hap2_bam

        File? ont_bam

        String sample_name
    }

    String locus_htt = "chr4:3,074,832-3,074,978|HTT"
    String locus_fmr1 = "chrX:147,911,969-147,912,124|FMR1"

    call Fetch as FetchHTT {
        input:
            bams = select_all([repeat_hap1_bam, repeat_hap2_bam, hprc_hap1_bam, hprc_hap2_bam]),
            ont_bam = ont_bam,
            locus = locus_htt,
            prefix = "HTT",
            sample_name = sample_name
    }

    call Abpoa as AbpoaHTT { input: fasta = FetchHTT.fasta }

    call Fetch as FetchFMR1 {
        input:
            bams = select_all([repeat_hap1_bam, repeat_hap2_bam, hprc_hap1_bam, hprc_hap2_bam]),
            ont_bam = ont_bam,
            locus = locus_fmr1,
            prefix = "FMR1",
            sample_name = sample_name
    }

    call Abpoa as AbpoaFMR1 { input: fasta = FetchFMR1.fasta }

    output {
        File msa_htt = AbpoaHTT.msa_fa
        File msa_fmr1 = AbpoaFMR1.msa_fa
    }
}

task Fetch {
    input {
        Array[String] bams
        String? ont_bam

        String locus
        String sample_name

        String prefix = "out"

        Int disk_size_gb = 2
        Int num_cpus = 4
    }

    command <<<
        set -euxo pipefail

        # Write bams to a file
        echo '~{sep="," bams}' | tr ',' '\n' > bams.txt

        hidive fetch -l '~{locus}' bams.txt | sed 's/>/>~{sample_name}#/' > ~{prefix}.fa 

        # If ONT BAM is provided, append its reads to the FASTA
        if [ -n "~{ont_bam}" ]; then
            hidive fetch -l '~{locus}' --haplotype 1 ~{ont_bam} | sed 's/>/>~{sample_name}_ont_1#/' >> ~{prefix}.fa
            hidive fetch -l '~{locus}' --haplotype 2 ~{ont_bam} | sed 's/>/>~{sample_name}_ont_2#/' >> ~{prefix}.fa
        fi
    >>>

    output {
        File fasta = "~{prefix}.fa"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_fix_consensus"
        memory: "2 GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}

task Abpoa {
    input {
        File fasta

        String prefix = "out"

        Int disk_size_gb = 2
        Int num_cpus = 4
    }

    command <<<
        set -euxo pipefail

        abpoa -m 0 -r 1 ~{fasta} > ~{prefix}.fa
    >>>

    output {
        File msa_fa = "~{prefix}.fa"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_fix_consensus"
        memory: "2 GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}