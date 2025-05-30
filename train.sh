#!/bin/bash

set -euxo pipefail

export GCS_REQUESTER_PAYS_PROJECT="broad-dsp-lrma"

#TRAIN_LOCUS="chr20:40,000,000-60,000,000"
#TRAIN_NAME="chr20"

TRAIN_LOCUS="chr19:53,000,000-57,000,000"
TRAIN_NAME="KIR"

#TRAIN_LOCUS="chr6:27,100,000-30,500,000"
#TRAIN_NAME="MHC" # https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1207-1

#TEST_LOCUS="chr6:29,941,260-29,949,572"
TEST_LOCUS="chr6:29,000,000-30,000,000"
TEST_NAME="HLA-A"

SAMPLE="HG00438"
REF="Homo_sapiens_assembly38.fasta"
OUTPUT="scratch/training/per_sample_model/$SAMPLE/$TRAIN_NAME"

export GCS_REQUESTER_PAYS_PROJECT="broad-dsp-lrma"
export RUST_BACKTRACE=full

cargo build
cargo build --release

#HIDIVE_DEBUG="./target/debug/hidive"
HIDIVE="./target/release/hidive"

mkdir -p $OUTPUT

# maternal haplotype training
if [ ! -f $OUTPUT/$SAMPLE.maternal.GRCh38_no_alt.bam ]; then
  gsutil -u $GCS_REQUESTER_PAYS_PROJECT cp gs://fc-4310e737-a388-4a10-8c9e-babe06aaf0cf/working/HPRC/$SAMPLE/assemblies/year1_f1_assembly_v2_genbank/alignment/assembly-to-reference/$SAMPLE.maternal.GRCh38_no_alt.bam $OUTPUT/$SAMPLE.maternal.GRCh38_no_alt.bam
  gsutil -u $GCS_REQUESTER_PAYS_PROJECT cp gs://fc-4310e737-a388-4a10-8c9e-babe06aaf0cf/working/HPRC/$SAMPLE/assemblies/year1_f1_assembly_v2_genbank/alignment/assembly-to-reference/$SAMPLE.maternal.GRCh38_no_alt.bam.bai $OUTPUT/$SAMPLE.maternal.GRCh38_no_alt.bam.bai
fi

# if file not in folder fetch
if [ ! -f $OUTPUT/$SAMPLE.maternal.train.GRCh38_no_alt.$TRAIN_NAME.fasta ]; then
  $HIDIVE fetch \
    -l $TRAIN_LOCUS \
    -o $OUTPUT/$SAMPLE.maternal.train.GRCh38_no_alt.$TRAIN_NAME.fasta \
    $OUTPUT/$SAMPLE.maternal.GRCh38_no_alt.bam
fi

# maternal haplotype testing
if [ ! -f $OUTPUT/$SAMPLE.maternal.test.GRCh38_no_alt.$TEST_NAME.fasta ]; then
  $HIDIVE fetch \
    -l $TEST_LOCUS \
    -o $OUTPUT/$SAMPLE.maternal.test.GRCh38_no_alt.$TEST_NAME.fasta \
    $OUTPUT/$SAMPLE.maternal.GRCh38_no_alt.bam
fi

# paternal haplotype train
if [ ! -f $OUTPUT/$SAMPLE.paternal.GRCh38_no_alt.bam ]; then
  gsutil -u $GCS_REQUESTER_PAYS_PROJECT cp gs://fc-4310e737-a388-4a10-8c9e-babe06aaf0cf/working/HPRC/$SAMPLE/assemblies/year1_f1_assembly_v2_genbank/alignment/assembly-to-reference/$SAMPLE.paternal.GRCh38_no_alt.bam $OUTPUT/$SAMPLE.paternal.GRCh38_no_alt.bam
  gsutil -u $GCS_REQUESTER_PAYS_PROJECT cp gs://fc-4310e737-a388-4a10-8c9e-babe06aaf0cf/working/HPRC/$SAMPLE/assemblies/year1_f1_assembly_v2_genbank/alignment/assembly-to-reference/$SAMPLE.paternal.GRCh38_no_alt.bam.bai $OUTPUT/$SAMPLE.paternal.GRCh38_no_alt.bam.bai
fi

if [ ! -f $OUTPUT/$SAMPLE.paternal.train.GRCh38_no_alt.$TRAIN_NAME.fasta ]; then
  $HIDIVE fetch \
    -l $TRAIN_LOCUS \
    -o $OUTPUT/$SAMPLE.paternal.train.GRCh38_no_alt.$TRAIN_NAME.fasta \
    $OUTPUT/$SAMPLE.paternal.GRCh38_no_alt.bam
fi

# paternal haplotype testing
if [ ! -f $OUTPUT/$SAMPLE.paternal.test.GRCh38_no_alt.$TEST_NAME.fasta ]; then
  $HIDIVE fetch \
    -l $TEST_LOCUS \
    -o $OUTPUT/$SAMPLE.paternal.test.GRCh38_no_alt.$TEST_NAME.fasta \
    $OUTPUT/$SAMPLE.paternal.GRCh38_no_alt.bam
fi

# fetch long reads training
if [ ! -f $OUTPUT/$SAMPLE.LR.$TRAIN_NAME.fasta ]; then
  $HIDIVE fetch \
    -l $TRAIN_LOCUS \
    -o $OUTPUT/$SAMPLE.LR.$TRAIN_NAME.fasta \
    gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/HPRC_grch38/$SAMPLE/$SAMPLE.bam
fi

# fetch long reads testing
if [ ! -f $OUTPUT/$SAMPLE.LR.$TEST_NAME.fasta ]; then
  $HIDIVE fetch \
    -l $TEST_LOCUS \
    -o $OUTPUT/$SAMPLE.LR.$TEST_NAME.fasta \
    gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/HPRC_grch38/$SAMPLE/$SAMPLE.bam
fi

# fetch short reads training
if [ ! -f $OUTPUT/$SAMPLE.SR.$TRAIN_NAME.fasta ]; then
  $HIDIVE rescue \
    -o $OUTPUT/$SAMPLE.SR.$TRAIN_NAME.fasta \
    -r $REF \
    -f $OUTPUT/$SAMPLE.LR.$TRAIN_NAME.fasta \
    $SAMPLE.final.cram
fi

# fetch short reads testing
if [ ! -f $OUTPUT/$SAMPLE.SR.$TEST_NAME.fasta ]; then
  $HIDIVE rescue \
    -o $OUTPUT/$SAMPLE.SR.$TEST_NAME.fasta \
    -r $REF \
    -f $OUTPUT/$SAMPLE.LR.$TEST_NAME.fasta \
    $SAMPLE.final.cram
fi

# train

if [ ! -f $OUTPUT/training.long_and_short_reads.$TRAIN_NAME.json ]; then
  $HIDIVE train -d \
    -o $OUTPUT/training.long_and_short_reads.$TRAIN_NAME.json \
    --long-read-seq-paths $OUTPUT/$SAMPLE.LR.$TRAIN_NAME.fasta \
    --short-read-seq-paths $OUTPUT/$SAMPLE.SR.$TRAIN_NAME.fasta \
    --truth-seq-paths	$OUTPUT/$SAMPLE.maternal.train.GRCh38_no_alt.$TRAIN_NAME.fasta $OUTPUT/$SAMPLE.paternal.train.GRCh38_no_alt.$TRAIN_NAME.fasta \
    --test-long-read-seq-paths $OUTPUT/$SAMPLE.LR.$TEST_NAME.fasta \
    --test-short-read-seq-paths $OUTPUT/$SAMPLE.SR.$TEST_NAME.fasta \
    --test-truth-seq-paths $OUTPUT/$SAMPLE.maternal.test.GRCh38_no_alt.$TEST_NAME.fasta $OUTPUT/$SAMPLE.paternal.test.GRCh38_no_alt.$TEST_NAME.fasta
fi

minimap2 -t 16 -ayYL -x sr Homo_sapiens_assembly38.fasta $OUTPUT/$SAMPLE.SR.$TRAIN_NAME.fasta | samtools sort --write-index -O BAM -o $OUTPUT/$SAMPLE.SR.$TRAIN_NAME.bam
