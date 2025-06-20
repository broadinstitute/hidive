version 1.0

workflow RescueAndLocityper {
    input {
        String sample_id
        File cram
        File crai

        File ref_fa_with_alt
        File ref_fai_with_alt
        File ref_cache_tar_gz
        File counts_jf

        # File reference
        # File reference_index
        # File db_path

        String locus_name
        String locus_coordinates

        File alleles_fa
        # File weights_file

        String locityper_docker = "eichlerlab/locityper:0.19.1"
    }

    call GenerateDB {
        input:
            reference = ref_fa_with_alt,
            reference_index = ref_fai_with_alt,
            counts_jf = counts_jf,
            locus_name = locus_name,
            locus_coordinates = locus_coordinates,
            alleles_fa = alleles_fa,
            docker = locityper_docker
    }

    call Rescue {
        input:
            long_reads_fastx = alleles_fa,
            short_reads_cram = cram,
            short_reads_crai = crai,
            ref_fa_with_alt = ref_fa_with_alt,
            ref_fai_with_alt = ref_fai_with_alt,
            ref_cache_tar_gz = ref_cache_tar_gz,
            prefix = sample_id
    }

    call LocityperPreprocessAndGenotype {
        input:
            sample_id = sample_id,
            input_fq1 = Rescue.fastq_gz,
            db_targz = GenerateDB.db_tar,
            counts_file = counts_jf,
            reference = ref_fa_with_alt,
            locus_name = locus_name,
            reference_index = ref_fai_with_alt,
            docker = locityper_docker,
            locityper_n_cpu = 4,
            locityper_mem_gb = 32
    }

    output {
        File results = LocityperPreprocessAndGenotype.genotype_tar
    }
}

task SubsetVCF {
    input {
        File vcf
        File bed
    }

    Int disk_size = 1 + ceil(size([vcf, bed], "GiB"))

    command <<<
        set -euxo pipefail

        bcftools view -R ~{bed} ~{vcf} | bgzip > subset.vcf.gz
        tabix -p vcf subset.vcf.gz
    >>>

    output {
        File subset_vcf_gz = "subset.vcf.gz"
        File subset_vcf_gz_tbi = "subset.vcf.gz.tbi"
    }

    runtime {
        memory: "8 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
        docker: "staphb/bcftools:1.22"
    }
}

task GenerateDBFromVCF {
    input {
        File reference
        File reference_index
        File counts_jf
        File vcf
        File bed
        String docker
    }

    Int disk_size = 10
    String output_tar = "vcf_db.tar.gz"

    command <<<
        set -euxo pipefail

        gunzip -c ~{reference} > reference.fa
        samtools faidx reference.fa

        locityper add -d vcf_db \
            -v ~{vcf} \
            -r reference.fa \
            -j ~{counts_jf} \
            -L ~{bed}

        echo "compressing DB"
        tar -czf ~{output_tar} vcf_db
        echo "done compressing DB"
    >>>

    runtime {
        memory: "8 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
        docker: docker
    }

    output {
        File db_tar = output_tar
    }
}

task GenerateDB {
    input {
        File reference
        File reference_index
        File counts_jf
        String locus_name
        String locus_coordinates
        File alleles_fa
        String docker
    }

    Int disk_size = 10
    String output_tar = locus_name + "db.tar.gz"

    command <<<
        set -euxo pipefail

        gunzip -c ~{reference} > reference.fa
        samtools faidx reference.fa

        locityper add -d ~{locus_name}.db \
            -r reference.fa \
            -j ~{counts_jf} \
            -l ~{locus_name} ~{locus_coordinates} ~{alleles_fa}

        find ~{locus_name}.db -type f -exec ls -lah {} \;
        
        echo "compressing DB"
        tar -czf ~{output_tar} ~{locus_name}.db
        echo "done compressing DB"
    >>>

    runtime {
        memory: "8 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
        docker: docker
    }

    output {
        File db_tar = output_tar
    }
}

task Rescue {
    input {
        File long_reads_fastx
        File short_reads_cram
        File short_reads_crai

        File ref_fa_with_alt
        File ref_fai_with_alt
        File ref_cache_tar_gz

        String prefix = "out"

        Int num_cpus = 16
    }

    Int disk_size_gb = 1 + 2*ceil(size([long_reads_fastx, short_reads_cram, short_reads_crai, ref_fa_with_alt, ref_fai_with_alt, ref_cache_tar_gz], "GB"))
    Int memory_gb = 3*num_cpus

    command <<<
        set -euxo pipefail

        mv ~{ref_fa_with_alt} Homo_sapiens_assembly38.fasta
        mv ~{ref_fai_with_alt} Homo_sapiens_assembly38.fasta.fai 
        mv ~{ref_cache_tar_gz} Homo_sapiens_assembly38.ref_cache.tar.gz

        tar xzf Homo_sapiens_assembly38.ref_cache.tar.gz >/dev/null 2>&1

        export REF_PATH="$(pwd)/ref/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s"
        export REF_CACHE="$(pwd)/ref/cache/%2s/%2s/%s"

        hidive rescue -r Homo_sapiens_assembly38.fasta -f ~{long_reads_fastx} ~{short_reads_cram} | gzip > ~{prefix}.fq.gz

        hidive fetch -l "chr17:72062001-76562000" -p 10000 Homo_sapiens_assembly38.fasta > chr17.fa
        hidive rescue -r Homo_sapiens_assembly38.fasta -f chr17.fa ~{short_reads_cram} | gzip >> ~{prefix}.fq.gz
    >>>

    output {
        File fastq_gz = "~{prefix}.fq.gz"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_interlaced"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
        maxRetries: 0
    }
}

task LocityperPreprocessAndGenotype {
    input {
        File input_fq1
        File? input_fq2
        File counts_file
        File reference
        File reference_index
        # File weights_file
        String docker
        File db_targz
        String sample_id
        String locus_name

        Int locityper_n_cpu
        Int locityper_mem_gb
    }

    Int disk_size = 80 + ceil(size(select_all([input_fq1, input_fq2]), "GiB"))
    String output_tar = sample_id + "." + locus_name + ".tar.gz"

    command <<<
        set -euxo pipefail

        gunzip -c ~{reference} > reference.fa
        samtools faidx reference.fa

        nthreads=$(nproc)
        echo "using ${nthreads} threads"

        mkdir -p locityper_prepoc
        locityper preproc -i ~{sep=" " select_all([input_fq1, input_fq2])} \
            --interleaved \
            -j ~{counts_file} \
            -@ ${nthreads} \
            --technology illumina \
            -r reference.fa \
            -o locityper_prepoc

        mkdir -p db
        tar --strip-components 1 -C db -xvzf ~{db_targz}

        mkdir -p out_dir
        locityper genotype -i ~{sep=" " select_all([input_fq1, input_fq2])} \
            -d db \
            -p locityper_prepoc \
            -@ ${nthreads} \
            --debug 2 \
            -o out_dir

        tar -czf ~{output_tar} out_dir
    >>>

    runtime {
        memory: "~{locityper_mem_gb} GB"
        cpu: locityper_n_cpu
        disks: "local-disk ~{disk_size} HDD"
        preemptible: 3
        docker: docker
    }

    output {
        File genotype_tar = output_tar
    }
}