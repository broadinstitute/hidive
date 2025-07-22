version 1.0

workflow BuildLocityperDB {
    input {
        File ref_fa_with_alt
        File ref_fai_with_alt
        File counts_jf

        File vcf
        File vcf_tbi
        File bed
    }

    call GenerateDBFromVCF {
        input:
            reference = ref_fa_with_alt,
            reference_index = ref_fai_with_alt,
            counts_jf = counts_jf,
            vcf = vcf,
            vcf_tbi = vcf_tbi,
            bed = bed
    }

    output {
        File locityper_db_tar_gz = GenerateDBFromVCF.db_tar
    }
}

task GenerateDBFromVCF {
    input {
        File reference
        File reference_index
        File counts_jf
        File vcf
        File vcf_tbi
        File bed
    }

    Int disk_size = 1 + 10*ceil(size([reference, vcf, counts_jf, bed], "GiB"))
    String output_tar = "vcf_db.tar.gz"

    command <<<
        set -euxo pipefail

        nthreads=$(nproc)
        echo "using ${nthreads} threads"

        gunzip -c ~{reference} > reference.fa
        samtools faidx reference.fa

        locityper add \
            -@ ${nthreads} \
            -d vcf_db \
            -v ~{vcf} \
            -r reference.fa \
            -j ~{counts_jf} \
            -L ~{bed}

        echo "compressing DB"
        tar -czf ~{output_tar} vcf_db

        echo "done compressing DB"
    >>>

    output {
        File db_tar = output_tar
    }

    runtime {
        memory: "8 GB"
        cpu: "4"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
        docker: "eichlerlab/locityper:0.19.1"
    }
}