#!/bin/bash

set -euxo pipefail

LOCUS="chr6:25,726,063-33,410,226"

cargo build --release

./target/release/hidive fetch \
	-l $LOCUS \
	-o HG00438.LR.HLA-A.fasta \
	gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/HPRC_grch38/PBCCSWholeGenome/HG00438/alignments/HG00438.bam

./target/release/hidive fetch \
	-l $LOCUS \
	-o HG00438.SR.HLA-A.fasta \
	gs://fc-1ee08173-e353-4494-ad28-7a3d7bd99734/working/HPRC/HG00438/raw_data/Illumina/child/HG00438.final.cram

minimap2 -ayYL --MD --eqx -x map-hifi -R '@RG\tID:lr\tSM:lr' GCA_000001405.15_GRCh38_no_alt_analysis_set.fa HG00438.LR.HLA-A.fasta | \
	samtools sort -O BAM -o HG00438.LR.HLA-A.bam && \
	samtools index HG00438.LR.HLA-A.bam

minimap2 -ayYL --MD --eqx -x map-hifi -R '@RG\tID:sr\tSM:sr' GCA_000001405.15_GRCh38_no_alt_analysis_set.fa HG00438.SR.HLA-A.fasta | \
	samtools sort -O BAM -o HG00438.SR.HLA-A.bam && \
	samtools index HG00438.SR.HLA-A.bam

minimap2 -ayYL --MD --eqx -x map-hifi -R '@RG\tID:ar\tSM:ar' HG00438.LR.HLA-A.fasta HG00438.SR.HLA-A.fasta | \
	samtools sort -O BAM -o HG00438.AR.HLA-A.bam && \
	samtools index HG00438.AR.HLA-A.bam

#./target/release/hidive build -r GCA_000001405.15_GRCh38_no_alt_analysis_set -o HLA-A.gfa HLA-A.fasta
