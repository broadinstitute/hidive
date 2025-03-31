version 1.0

workflow HidiveSummarizeMSA {
    input {
        Array[File] all_msa_htt
        Array[File] all_msa_fmr1
    }

    call AbpoaAll as AbpoaAllHTT { input: all_msa = all_msa_htt, prefix = "HTT" }
    call AbpoaAll as AbpoaAllFMR1 { input: all_msa = all_msa_fmr1, prefix = "FMR1" }

    output {
        File msa_htt = AbpoaAllHTT.msa_fa
        File msa_fmr1 = AbpoaAllFMR1.msa_fa
    }
}

task AbpoaAll {
    input {
        Array[File] all_msa

        String prefix = "out"

        Int disk_size_gb = 2
        Int num_cpus = 4
    }

    command <<<
        set -euxo pipefail

        sed 's/-//g' ~{sep="\n" all_msa} > all.fa

        abpoa -m 0 -r 1 all.fa > ~{prefix}.fa
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