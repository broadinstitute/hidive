version 1.0

workflow HidiveTrain {
    input {
        File mat_assembly_bam
        File pat_assembly_bam

        File long_reads_bam

        File short_reads_cram
        File short_reads_crai

        File ref_fa_with_alt
        File ref_fai_with_alt
        File ref_cache_tar_gz

        String train_locus = "chr19:53,000,000-57,000,000" # KIR
        String test_locus  = "chr6:29,000,000-30,000,000"  # HLA-A
    }

    call PrepareLocus as PrepareTrainLocus { input: locus = train_locus }
    call PrepareLocus as PrepareTestLocus { input: locus = test_locus }

    call Fetch as FetchMatTrain { input: bam = mat_assembly_bam, loci = PrepareTrainLocus.bed }
    call Fetch as FetchPatTrain { input: bam = pat_assembly_bam, loci = PrepareTrainLocus.bed }

    call Fetch as FetchMatTest { input: bam = mat_assembly_bam, loci = PrepareTestLocus.bed }
    call Fetch as FetchPatTest { input: bam = pat_assembly_bam, loci = PrepareTestLocus.bed }

    call Fetch as FetchSampleTrain { input: bam = long_reads_bam, loci = PrepareTrainLocus.bed }
    call Fetch as FetchSampleTest { input: bam = long_reads_bam, loci = PrepareTestLocus.bed }

    call Rescue as RescueSampleTrain {
        input:
            long_reads_fasta = FetchSampleTrain.fasta,
            short_reads_cram = short_reads_cram,
            short_reads_crai = short_reads_crai,
            ref_fa_with_alt = ref_fa_with_alt,
            ref_fai_with_alt = ref_fai_with_alt,
            ref_cache_tar_gz = ref_cache_tar_gz,
    }

    call Rescue as RescueSampleTest {
        input:
            long_reads_fasta = FetchSampleTest.fasta,
            short_reads_cram = short_reads_cram,
            short_reads_crai = short_reads_crai,
            ref_fa_with_alt = ref_fa_with_alt,
            ref_fai_with_alt = ref_fai_with_alt,
            ref_cache_tar_gz = ref_cache_tar_gz,
    }

    scatter (kmer_size in [17, 37, 57, 77, 97, 117]) {
        call Train {
            input:
                long_read_train_fa = FetchSampleTrain.fasta,
                long_read_test_fa = FetchSampleTest.fasta,
                short_read_train_fa = RescueSampleTrain.fasta,
                short_read_test_fa = RescueSampleTest.fasta,
                mat_train_fa = FetchMatTrain.fasta,
                mat_test_fa = FetchMatTest.fasta,
                pat_train_fa = FetchPatTrain.fasta,
                pat_test_fa = FetchPatTest.fasta,
                kmer_size = kmer_size
        }
    }

    output {
        File model_k17  = Train.model[0]
        File model_k37  = Train.model[1]
        File model_k57  = Train.model[2]
        File model_k77  = Train.model[3]
        File model_k97  = Train.model[4]
        File model_k117 = Train.model[5]
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
        Int padding = 0

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

task Train {
    input {
        File long_read_train_fa
        File long_read_test_fa

        File short_read_train_fa
        File short_read_test_fa

        File mat_train_fa
        File mat_test_fa

        File pat_train_fa
        File pat_test_fa

        Int kmer_size
    }

    Int disk_size_gb = 1 + 2*ceil(size([long_read_train_fa, long_read_test_fa, short_read_train_fa, short_read_test_fa, mat_train_fa, mat_test_fa, pat_train_fa, pat_test_fa], "GB"))

    command <<<
        set -euxo pipefail

        hidive train \
            --debug \
            --long-read-seq-paths ~{long_read_train_fa} \
            --short-read-seq-paths ~{short_read_train_fa} \
            --truth-seq-paths ~{mat_train_fa} ~{pat_train_fa} \
            --test-long-read-seq-paths ~{long_read_test_fa} \
            --test-short-read-seq-paths ~{short_read_test_fa} \
            --test-truth-seq-paths ~{mat_test_fa} ~{pat_test_fa} \
            -o model.k~{kmer_size}.json
    >>>

    output {
        File model = "model.k~{kmer_size}.json"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:0.1.101"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}