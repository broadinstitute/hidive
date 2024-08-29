#!/bin/bash

set -euxo pipefail

LOCUS="chr6:32,438,878-32,446,046"
NAME="HLA-DRA"

OUTPUT="scratch/$NAME"

cargo build --release

./target/release/hidive fetch \
	-l chr6:29,940,532-29,947,870 \
	-o HLA-A.fasta \
	gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/HPRC_grch38/PBCCSWholeGenome/HG002/alignments/HG002.bam \
	gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/HPRC_grch38/PBCCSWholeGenome/HG00438/alignments/HG00438.bam \
	gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/HPRC_grch38/PBCCSWholeGenome/HG02559/alignments/HG02559.bam \
	gs://broad-dsde-methods-long-reads/resources/references/grch38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

./target/release/hidive build -r GCA_000001405.15_GRCh38_no_alt_analysis_set -o HLA-A.gfa HLA-A.fasta


./target/release/hidive assemble --output haplotype.fasta ./HLA-A.sp.gfa ./HLA-A.fasta 3 "HG00438"

minimap2 -ayYL --MD --eqx -x map-hifi -R '@RG\tID:asm\tSM:asm' ./MHC-HG002.1.fa haplotype.fasta | \
	samtools sort -O BAM --write-index -o $OUTPUT/assemble.bam