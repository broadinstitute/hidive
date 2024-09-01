#!/bin/bash

set -euxo pipefail

#LOCUS="chr6:25,726,063-33,410,226"
#NAME="MHC"

LOCUS="chr20:40,000,000-60,000,000"
NAME="chr20_v2"

#LOCUS="chr6:29,941,260-29,949,572"
#NAME="HLA-A"

OUTPUT="scratch/training/$NAME"

export GCS_REQUESTER_PAYS_PROJECT="broad-dsp-lrma"

cargo build --release

mkdir -p $OUTPUT

# maternal haplotype

gsutil -u $GCS_REQUESTER_PAYS_PROJECT cp gs://fc-4310e737-a388-4a10-8c9e-babe06aaf0cf/working/HPRC/HG00438/assemblies/year1_f1_assembly_v2_genbank/alignment/assembly-to-reference/HG00438.maternal.GRCh38_no_alt.bam $OUTPUT/HG00438.maternal.GRCh38_no_alt.bam
gsutil -u $GCS_REQUESTER_PAYS_PROJECT cp gs://fc-4310e737-a388-4a10-8c9e-babe06aaf0cf/working/HPRC/HG00438/assemblies/year1_f1_assembly_v2_genbank/alignment/assembly-to-reference/HG00438.maternal.GRCh38_no_alt.bam.bai $OUTPUT/HG00438.maternal.GRCh38_no_alt.bam.bai

./target/release/hidive fetch \
	-l $LOCUS \
	-o $OUTPUT/HG00438.maternal.GRCh38_no_alt.$NAME.fasta \
	$OUTPUT/HG00438.maternal.GRCh38_no_alt.bam

# paternal haplotype

#gsutil -u $GCS_REQUESTER_PAYS_PROJECT cp gs://fc-4310e737-a388-4a10-8c9e-babe06aaf0cf/working/HPRC/HG00438/assemblies/year1_f1_assembly_v2_genbank/alignment/assembly-to-reference/HG00438.paternal.GRCh38_no_alt.bam $OUTPUT/HG00438.paternal.GRCh38_no_alt.bam
#gsutil -u $GCS_REQUESTER_PAYS_PROJECT cp gs://fc-4310e737-a388-4a10-8c9e-babe06aaf0cf/working/HPRC/HG00438/assemblies/year1_f1_assembly_v2_genbank/alignment/assembly-to-reference/HG00438.paternal.GRCh38_no_alt.bam.bai $OUTPUT/HG00438.paternal.GRCh38_no_alt.bam.bai

./target/release/hidive fetch \
	-l $LOCUS \
	-o $OUTPUT/HG00438.paternal.GRCh38_no_alt.$NAME.fasta \
	$OUTPUT/HG00438.paternal.GRCh38_no_alt.bam

# fetch long reads

./target/release/hidive fetch \
	-l $LOCUS \
	-o $OUTPUT/HG00438.LR.$NAME.fasta \
	gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/HPRC_grch38/PBCCSWholeGenome/HG00438/alignments/HG00438.bam

# fetch short reads

./target/release/hidive fetch \
	-l $LOCUS \
	-o $OUTPUT/HG00438.SR.$NAME.fasta \
	HG00438.final.cram

# train

./target/release/hidive train \
	-o $OUTPUT/training.long_and_short_reads.$NAME.json \
	-l $OUTPUT/HG00438.LR.$NAME.fasta \
	-s $OUTPUT/HG00438.SR.$NAME.fasta \
	$OUTPUT/HG00438.maternal.GRCh38_no_alt.$NAME.fasta \
	$OUTPUT/HG00438.paternal.GRCh38_no_alt.$NAME.fasta

./target/release/hidive train \
	-o $OUTPUT/training.long_reads.$NAME.json \
	-l $OUTPUT/HG00438.LR.$NAME.fasta \
	$OUTPUT/HG00438.maternal.GRCh38_no_alt.$NAME.fasta \
	$OUTPUT/HG00438.paternal.GRCh38_no_alt.$NAME.fasta
