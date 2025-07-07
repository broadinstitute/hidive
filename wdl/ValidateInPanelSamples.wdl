version 1.0

workflow ValidateInPanelSamples {
    input {
        String sample_id

        File long_reads_bam

        File cram
        File crai

        File ref_fa_with_alt
        File ref_fai_with_alt
        File ref_cache_tar_gz
        File counts_jf

        File vcf
        File vcf_tbi
        File bed
        File filtered_bed

        Int N = 20
    }

    call SubsetVCF {
        input:
            vcf = vcf,
            vcf_tbi = vcf_tbi,
            bed = bed,
            sample_id = sample_id
    }

    call GenerateDBFromVCF {
        input:
            reference = ref_fa_with_alt,
            reference_index = ref_fai_with_alt,
            counts_jf = counts_jf,
            vcf = SubsetVCF.subset_vcf,
            bed = bed,
    }

    call Fetch {
        input:
            bam = long_reads_bam,
            loci = bed,
            padding = 50000,
            prefix = sample_id
    }

    call Rescue {
        input:
            long_reads_fastx = Fetch.fastq,
            short_reads_cram = cram,
            ref_fa_with_alt = ref_fa_with_alt,
            ref_fai_with_alt = ref_fai_with_alt,
            ref_cache_tar_gz = ref_cache_tar_gz,
            prefix = sample_id
    }

    call DeduplicateFastq {
        input:
            fastq = Rescue.fastq_gz,
            prefix = sample_id
    }

    call SplitBedNames { input: bed = filtered_bed, N = N }

    scatter (names_file in SplitBedNames.name_parts) {
        call LocityperPreprocessAndGenotype {
            input:
                sample_id = sample_id,
                input_fq1 = DeduplicateFastq.fastq_gz,
                db_targz = GenerateDBFromVCF.db_tar,
                counts_file = counts_jf,
                reference = ref_fa_with_alt,
                locus_names = read_lines(names_file),
                reference_index = ref_fai_with_alt,
                locityper_n_cpu = 4,
                locityper_mem_gb = 32
        }
    }

    call CombineTarFiles {
        input:
            tar_files = LocityperPreprocessAndGenotype.genotype_tar,
            sample_id = sample_id
    }

    call Summarize { input: sample_id = sample_id, genotype_tar = CombineTarFiles.combined_tar }

    output {
        File results = Summarize.summary_csv
    }
}

task SubsetVCF {
    input {
        File vcf
        File vcf_tbi
        File bed
        String? sample_id
    }

    Int disk_size = 1 + ceil(size([vcf, bed], "GiB"))

    command <<<
        set -euxo pipefail

        bcftools view ~{if defined(sample_id) then "-s " + sample_id else ""} -R ~{bed} ~{vcf} > subset.vcf
    >>>

    output {
        File subset_vcf = "subset.vcf"
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
    }

    Int disk_size = 1 + 4*ceil(size([reference, vcf, counts_jf, bed], "GiB"))
    String output_tar = "vcf_db.tar.gz"

    command <<<
        set -euxo pipefail

        gunzip -c ~{reference} > reference.fa
        samtools faidx reference.fa

        mv ~{vcf} subset.vcf
        bgzip subset.vcf
        tabix -p vcf subset.vcf.gz

        locityper add -d vcf_db \
            -v subset.vcf.gz \
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
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
        docker: "eichlerlab/locityper:0.19.1"
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

        hidive fetch -l ~{select_first([locus, loci])} -p ~{padding} ~{bam} > ~{prefix}.fq
    >>>

    output {
        File fastq = "~{prefix}.fq"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:0.1.122"
        memory: "2 GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
    }
}

task Rescue {
    input {
        File long_reads_fastx
        String short_reads_cram
        # File short_reads_crai

        File ref_fa_with_alt
        File ref_fai_with_alt
        File ref_cache_tar_gz

        String prefix = "out"

        Int num_cpus = 16
    }

    Int disk_size_gb = 1 + 2*ceil(size([long_reads_fastx, short_reads_cram, ref_fa_with_alt, ref_fai_with_alt, ref_cache_tar_gz], "GB"))
    Int memory_gb = 8*num_cpus

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
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:kvg_locityper_modes"
        memory: "~{memory_gb} GB"
        cpu: num_cpus
        disks: "local-disk ~{disk_size_gb} SSD"
        maxRetries: 2
    }
}

task DeduplicateFastq {
    input {
        File fastq
        String prefix = "out"
    }

    Int disk_size_gb = 1 + 2*ceil(size([fastq], "GB"))
    Int memory_gb = 2

    command <<<
        set -x

        python3 <<EOF
import gzip
from collections import defaultdict

seen_reads = defaultdict(int)
written_reads = set()

with gzip.open("~{fastq}", "rt") as f_in, gzip.open("~{prefix}.fq.gz", "wt") as f_out:
    while True:
        # Try to read the 4 lines that make up a FASTQ record
        name = f_in.readline().strip()
        if not name:  # EOF
            break
            
        seq = f_in.readline()
        plus = f_in.readline()
        qual = f_in.readline()
        
        # Get just the read name without /1 or /2
        read_name = name.split()[0]
        
        # Count this occurrence
        seen_reads[read_name] += 1
        
        # Only write if we haven't written this read name yet
        # and we've seen it exactly twice (paired)
        if seen_reads[read_name] <= 2 and read_name not in written_reads:
            f_out.write(f"{name}\n{seq}{plus}{qual}")
            if seen_reads[read_name] == 2:
                written_reads.add(read_name)

EOF
    >>>

    output {
        File fastq_gz = "~{prefix}.fq.gz"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:0.1.122"
        memory: "~{memory_gb} GB"
        cpu: 2
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
        File db_targz
        String sample_id
        Array[String] locus_names

        Int locityper_n_cpu
        Int locityper_mem_gb

        String docker = "eichlerlab/locityper:0.19.1"
    }

    Int disk_size = 1 + 1*length(locus_names) + 2*ceil(size(select_all([input_fq1, input_fq2, counts_file, reference, reference_index, db_targz]), "GiB"))
    String output_tar = sample_id + ".locityper.tar.gz"

    command <<<
        set -euxo pipefail
        
        df -h

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
            --interleaved \
            -d db \
            -p locityper_prepoc \
            -@ ${nthreads} \
            --debug 2 \
            --subset-loci ~{sep=" " locus_names} \
            -o out_dir

        tar -czf ~{output_tar} out_dir
    >>>

    runtime {
        memory: "~{locityper_mem_gb} GB"
        cpu: locityper_n_cpu
        disks: "local-disk ~{disk_size} HDD"
        preemptible: 0
        docker: docker
    }

    output {
        File genotype_tar = output_tar
    }
}

task SplitBedNames {
    input {
        File bed
        Int N = 20
    }

    Int disk_size = 1 + ceil(size(bed, "GiB"))

    command <<<
        set -euxo pipefail

        cut -f4 ~{bed} | split -l ~{N} - split_part_ && wc -l split_part_*
    >>>

    output {
        Array[File] name_parts = glob("split_part_*")
    }

    runtime {
        memory: "4 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
        docker: "staphb/bcftools:1.22"
    }
}

task CombineTarFiles {
    input {
        Array[File] tar_files
        String sample_id
    }

    Int disk_size = 1 + 2*ceil(size(tar_files, "GiB"))
    String combined_tar = sample_id + ".combined.tar.gz"

    command <<<
        set -euxo pipefail

        # Get the first tar file to start with
        first_tar=$(echo ~{sep=" " tar_files} | cut -d' ' -f1)
        echo "Starting with first tar file: $first_tar"
        
        # Copy the first tar file as our base
        cp "$first_tar" ~{combined_tar}
        
        # Append the rest of the tar files (skip the first one)
        remaining_tars=$(echo ~{sep=" " tar_files} | cut -d' ' -f2-)
        for tar_file in $remaining_tars; do
            echo "Appending $tar_file"
            tar -Af ~{combined_tar} "$tar_file"
        done
        
        echo "Combined tar file created successfully"
        ls -lh ~{combined_tar}
    >>>

    output {
        File combined_tar = combined_tar
    }

    runtime {
        memory: "8 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:0.1.122"
    }
}

task Summarize {
    input {
        String sample_id
        File genotype_tar
    }

    Int disk_size = 1 + 10*ceil(size(genotype_tar, "GiB"))

    command <<<
        set -euxo pipefail

        tar -xzvf ~{genotype_tar}

        mv out_dir ~{sample_id}
        python3 /locityper/extra/into_csv.py -i ./~{sample_id} -o gts.csv
    >>>

    output {
        File summary_csv = "gts.csv"
    }

    runtime {
        memory: "4 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
        docker: "us.gcr.io/broad-dsp-lrma/lr-hidive:0.1.122"
    }
}