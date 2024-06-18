#!/bin/bash

set -euxo pipefail

cargo build --release

./target/release/hidive fetch \
	-l chr6:29,940,532-29,947,870 \
	-o HLA-A.fasta \
	hprc_bams.txt \
	gs://broad-dsde-methods-long-reads/resources/references/grch38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

./target/release/hidive build -r GCA_000001405.15_GRCh38_no_alt_analysis_set -o HLA-A.gfa HLA-A.fasta
