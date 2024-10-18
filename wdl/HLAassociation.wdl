version 1.0

workflow HLAAssociation {

    input {
        File genotype_vcf
        File reference_vcf
        File reference_vcf_tbi
        File phenotype
        File covariate
        String prefix
        String referenceprefix
        String imputationprefix
        String association_prefix
        String trait_name
        Int makereferencememory
    }

    call FormatTransition { input:
        vcf = genotype_vcf,
        outputprefix = prefix
    }

    call MakeReference { input:
        reference_vcf = reference_vcf,
        reference_vcf_tbi = reference_vcf_tbi,
        outputprefix = referenceprefix,
        memory = makereferencememory
    }

    call SNP2HLA { input:
        genotype = FormatTransition.genotype_data,
        genotype_prefix = prefix,
        reference = MakeReference.reference_data,
        reference_prefix = referenceprefix,
        outputprefix = imputationprefix

    }

    call Association { input:
        prefix = imputationprefix,
        imputed_data = SNP2HLA.outputfiles,
        phenotype = phenotype,
        covariate = covariate,
        pheno_name = trait_name,
        outputprefix = association_prefix,
        memory = 4

    }

    output {
        Array[File] imputed_data = SNP2HLA.outputfiles
        Array[File] association_results = Association.outputfiles
    }
}

task FormatTransition {
    input {
        File vcf
        String outputprefix
        Int memory
    }

    String docker_dir = "/"
    String work_dir = "/cromwell_root/"

    command <<<

        cp ~{docker_dir}/get_duprem_var.py ~{work_dir}/get_duprem_var.py
        cd ~{work_dir}

        plink --vcf ~{vcf} --make-bed --out ~{outputprefix} --vcf-half-call m 
        # remove duplicated SNPs
        plink --bfile ~{outputprefix} --missing --out ~{outputprefix}
        python get_duprem_var.py ~{outputprefix}
        plink --bfile ~{outputprefix}  --exclude ~{outputprefix}.remdup.snp --make-bed --out ~{outputprefix}

    >>>

    output {
        # File bed = "~{outputprefix}.bed"
        # File bim = "~{outputprefix}.bim"
        # File fam = "~{outputprefix}.fam"
        # File log = "~{outputprefix}.log"
        # File nosex = "~{outputprefix}.nosex"

        Array[File] genotype_data = glob("~{outputprefix}.*")
    }
    ###################
    runtime {
        cpu: 4
        memory: memory + " GiB"
        disks: "local-disk 375 LOCAL"
        bootDiskSizeGb: 10
        preemptible_tries:     1
        max_retries:           0
        docker:"hangsuunc/hla-tapas:v1"
    }
}

task MakeReference {
    input {
        File reference_vcf
        File reference_vcf_tbi
        String outputprefix
        Int memory
    }
    command <<<
        
        plink --vcf ~{reference_vcf} --keep-allele-order --make-bed --out ~{outputprefix}
        plink --bfile ~{outputprefix} --freq --out ~{outputprefix}.FRQ
        zcat ~{reference_vcf} | grep -v "#" | awk '{print $3,$2,$4,$5}' > ~{outputprefix}.markers
        mv ~{reference_vcf} ~{outputprefix}.bgl.phased.vcf.gz
        mv ~{reference_vcf_tbi} ~{outputprefix}.bgl.phased.vcf.gz.tbi

        # plink --bfile ~{outputprefix} --hardy --allow-no-sex --out ~{outputprefix}.miss.frq.diff
        # cat ~{outputprefix}.miss.frq.diff.hwe | awk '/ALL/{split($6,V,"/"); if(V[1]+V[2]+V[3]==0){Freq=0}else{Freq=(V[2]/2+V[1])/(V[1]+V[2]+V[3])}print $2, $4, $5, Freq}' > ~{outputprefix}.miss.frq.diff.hwe.freq 
        # cat ~{outputprefix}.miss.frq.diff.hwe.freq | sort -k1,1 | join - <(cat ~{outputprefix}.bim | awk '{print $2, $1 "_" $4}' | sort -k1,1) | awk '{print $5, $2, $3, $4}' | sort -k1,1 > ~{outputprefix}.Ref.Frq.chr_pos_allele

    >>>

    output {
        # File bed = "~{outputprefix}.bed"
        # File bim = "~{outputprefix}.bim"
        # File fam = "~{outputprefix}.fam"
        # File log = "~{outputprefix}.log"
        # File nosex = "~{outputprefix}.nosex"
        # File hwe = "~{outputprefix}.miss.frq.diff.hwe"
        # File freq = "~{outputprefix}.miss.frq.diff.hwe.freq"
        # File pos_allele = "~{outputprefix}.Ref.Frq.chr_pos_allele"

        Array[File] reference_data = glob("~{outputprefix}.*")
    }
    ###################
    runtime {
        cpu: 4
        memory: memory + " GiB"
        disks: "local-disk 375 LOCAL"
        bootDiskSizeGb: 10
        preemptible_tries:     1
        max_retries:           0
        docker:"hangsuunc/hla-tapas:v1"
    }
}


task SNP2HLA {

    input {
        Array[File] genotype
        String genotype_prefix
        Array[File] reference
        String reference_prefix
        String outputprefix
        Int memory
    }

    parameter_meta {
        genotype: {localization_optional: true}
        reference: {localization_optional: true}
    }
    
    String docker_dir = "/"
    String work_dir = "/cromwell_root/"
    
    command <<<
        set -euxo pipefail

        cp ~{docker_dir}/rename_bim.py ~{work_dir}/rename_bim.py
        cd ~{work_dir}

        # we do single-sample phased VCFs localization ourselves
        time \
        gcloud storage cp ~{sep=" " genotype} /cromwell_root/

        time \
        gcloud storage cp ~{sep=" " reference} /cromwell_root/


        mv ~{genotype_prefix}.bim ~{genotype_prefix}.bim.old
        python rename_bim.py ~{reference_prefix}.bim ~{genotype_prefix}.bim

        cd /HLA-TAPAS

        python /HLA-TAPAS/SNP2HLA.py \
        -i "/cromwell_root/~{genotype_prefix}" \
        -o ~{outputprefix} \
        -rf "/cromwell_root/~{reference_prefix}" \
        --nthreads 10 \
        --java-mem=80g --tolerated-diff=0.5 --dependency "/HLA-TAPAS/dependency"

        # python -m SNP2HLA \
        # --target "/cromwell_root/~{genotype_prefix}" \
        # --out ~{outputprefix} \
        # --reference "/cromwell_root/~{reference_prefix}" \
        # --nthreads 2 \
        # --mem 4g

        cp /HLA-TAPAS/~{outputprefix}* /cromwell_root/

    >>>

    output {
        Array[File] outputfiles = glob("~{outputprefix}*")
        # Array[File] vcf_by_sample_tbi = glob("output/*vcf.gz.tbi")

    }
    ###################
    runtime {
        cpu: 4
        memory: memory + " GiB"
        disks: "local-disk 375 LOCAL"
        bootDiskSizeGb: 10
        preemptible_tries:     1
        max_retries:           0
        docker:"hangsuunc/hla-tapas:v1"
    }
}

task Association {
    input {
        String prefix
        Array[File] imputed_data
        File phenotype
        File covariate
        String pheno_name
        String outputprefix
        Int memory
    }

    command <<<
        set -euxo pipefail

        # Created ${imputed}.{psam, pgen, pvar}
        plink2 \
            --vcf ~{prefix}.bgl.phased.vcf.gz dosage=DS \
            --make-pgen --out ~{prefix}

        plink2 --pfile ${imputed} \
            --pheno ~{phenotype} \
            --covar ~{covariate} \
            --pheno-name ~{pheno_name} \
            --glm omit-ref hide-covar cols=chrom,pos,ref,alt,test,nobs,beta,se,ci,tz,p,a1freqcc,a1freq \
            --ci 0.95 \
            --out ~{outputprefix} \
            --covar-variance-standardize
    >>>

    output {
        Array[File] outputfiles = glob("~{outputprefix}*")
    }

    runtime {
        cpu: 4
        memory: memory + " GiB"
        disks: "local-disk 375 LOCAL"
        bootDiskSizeGb: 10
        preemptible_tries:     1
        max_retries:           0
        docker:"hangsuunc/hla-tapas:v1"
    }
}