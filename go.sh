#!/bin/bash

set -euxo pipefail

cargo build --release

./target/release/hidive fetch \
	-r \
	-l chr6:29,940,532-29,947,870 \
	-o HG002.HLA-A.bam \
	gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/HPRC_grch38/PBCCSWholeGenome/HG002/alignments/HG002.bam

./target/release/hidive fetch \
	-r \
	-l chr6:29,940,532-29,947,870 \
	-o HG00438.HLA-A.bam \
	gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/HPRC_grch38/PBCCSWholeGenome/HG00438/alignments/HG00438.bam

./target/release/hidive fetch \
	-r \
	-l chr6:29,940,532-29,947,870 \
	-o HG02559.HLA-A.bam \
	gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/HPRC_grch38/PBCCSWholeGenome/HG02559/alignments/HG02559.bam

samtools fasta HG002.HLA-A.bam   | sed '/^>/s/$/\|HG002/' > HG002.HLA-A.fasta
samtools fasta HG00438.HLA-A.bam | sed '/^>/s/$/\|HG00438/' > HG00438.HLA-A.fasta
samtools fasta HG02559.HLA-A.bam | sed '/^>/s/$/\|HG02559/' > HG02559.HLA-A.fasta

cat HG002.HLA-A.fasta HG00438.HLA-A.fasta HG02559.HLA-A.fasta > HLA-A.fasta

./target/release/hidive build -l chr6:29,940,532-29,947,870 -o HLA-A.gfa HLA-A.fasta GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
