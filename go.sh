#!/bin/bash

set -euxo pipefail

#LOCUS="chr6:25,726,063-33,410,226"
LOCUS="chr6:29,941,260-29,949,572"
NAME="HLA-A"

cargo build --release

./target/release/hidive fetch \
	-l $LOCUS \
	-o HG00438.LR.$NAME.fasta \
	gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/HPRC_grch38/PBCCSWholeGenome/HG00438/alignments/HG00438.bam

./target/release/hidive fetch \
	-l $LOCUS \
	-o HG00438.SR.$NAME.fasta \
	gs://fc-1ee08173-e353-4494-ad28-7a3d7bd99734/working/HPRC/HG00438/raw_data/Illumina/child/HG00438.final.cram

minimap2 -ayYL --MD --eqx -x map-hifi -R '@RG\tID:lr\tSM:lr' GCA_000001405.15_GRCh38_no_alt_analysis_set.fa HG00438.LR.$NAME.fasta | \
	samtools sort -O BAM -o HG00438.LR.$NAME.bam && \
	samtools index HG00438.LR.$NAME.bam

minimap2 -ayYL --MD --eqx -x map-hifi -R '@RG\tID:sr\tSM:sr' GCA_000001405.15_GRCh38_no_alt_analysis_set.fa HG00438.SR.$NAME.fasta | \
	samtools sort -O BAM -o HG00438.SR.$NAME.bam && \
	samtools index HG00438.SR.$NAME.bam

minimap2 -ayYL --MD --eqx -x map-hifi -R '@RG\tID:ar\tSM:ar' HG00438.LR.$NAME.fasta HG00438.SR.$NAME.fasta | \
	samtools sort -O BAM -o HG00438.AR.$NAME.bam && \
	samtools index HG00438.AR.$NAME.bam

#./target/release/hidive build -r GCA_000001405.15_GRCh38_no_alt_analysis_set -o $NAME.gfa $NAME.fasta

./target/release/hidive coassemble -o coassemble.fa -s HG00438.SR.$NAME.fasta HG00438.LR.$NAME.fasta

minimap2 -ayYL --MD --eqx -x map-hifi -R '@RG\tID:asm\tSM:asm' GCA_000001405.15_GRCh38_no_alt_analysis_set.fa coassemble.fa | \
	samtools sort -O BAM -o coassemble.bam && \
	samtools index coassemble.bam
