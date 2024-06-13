#!/bin/bash

set -euxo pipefail

cargo build --release >/dev/null 2>&1

./target/debug/hidive fetch \
	-r \
	-o HG002.HLA-A.bam \
	-l chr6:29,940,532-29,947,870 \
	gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/HPRC_grch38/PBCCSWholeGenome/HG002/alignments/HG002.bam

samtools index HG002.HLA-A.bam 
