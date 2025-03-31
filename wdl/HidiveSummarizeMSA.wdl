version 1.0

workflow HidiveSummarizeMSA {
    input {
        Array[File] all_msa_htt
        Array[File] all_msa_fmr1
    }

    call CombineFastas as CombineFastasHTT { input: all_msa = all_msa_htt, prefix = "HTT" }
    call AbpoaAll as AbpoaAllHTT { input: all_msa = CombineFastasHTT.msa_fa, prefix = "HTT" }

    call CombineFastas as CombineFastasFMR1 { input: all_msa = all_msa_fmr1, prefix = "FMR1" }
    call AbpoaAll as AbpoaAllFMR1 { input: all_msa = CombineFastasFMR1.msa_fa, prefix = "FMR1" }

    output {
        File msa_htt = AbpoaAllHTT.msa_fa
        File msa_fmr1 = AbpoaAllFMR1.msa_fa
    }
}

task CombineFastas {
    input {
        Array[String] all_msa

        String prefix = "out"

        Int disk_size_gb = 2
        Int num_cpus = 4
    }

    command <<<
        set -euxo pipefail

        echo '~{sep="\n" all_msa}' > file_list.txt

        while read -r file; do
            gsutil cat "$file" | sed 's/-//g' >> ~{prefix}.fa
        done < file_list.txt
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

task AbpoaAll {
    input {
        File all_msa

        String prefix = "out"

        Int disk_size_gb = 2
        Int num_cpus = 4
        Int memory_gb = 16
    }

    command <<<
        set -euxo pipefail

        abpoa -m 0 -r 1 ~{all_msa} > ~{prefix}.fa
    >>>

    output {
        File msa_fa = "~{prefix}.fa"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_fix_consensus"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}