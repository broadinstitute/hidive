#!/bin/bash

set -euxo pipefail

#LOCUS="chr6:25,726,063-33,410,226"
#NAME="MHC"

#TRAIN_LOCUS="chr20:40,000,000-60,000,000"
#TRAIN_NAME="chr20"

#TEST_LOCUS="chr6:40,000,000-55,000,000"
#TEST_NAME="chr6"

TEST_LOCUS="chr6:29,941,260-29,949,572"
TEST_NAME="HLA-A"

#LOCUS="chr6:29,941,260-29,949,572"
#NAME="HLA-A"

SAMPLE="HG01071"

MODEL_PATH="/Users/bshifaw/IdeaProjects/hidive/scratch/training/testtrain/KIR/training.long_and_short_reads.KIR.json"

OUTPUT="scratch/training/eval/$TEST_NAME"

export GCS_REQUESTER_PAYS_PROJECT="broad-dsp-lrma"

cargo build --release

mkdir -p $OUTPUT

# maternal haplotype
if [ ! -f $OUTPUT/$SAMPLE.maternal.GRCh38_no_alt.bam ]; then
  gsutil -u $GCS_REQUESTER_PAYS_PROJECT cp gs://fc-4310e737-a388-4a10-8c9e-babe06aaf0cf/working/HPRC/$SAMPLE/assemblies/year1_f1_assembly_v2_genbank/alignment/assembly-to-reference/$SAMPLE.maternal.GRCh38_no_alt.bam $OUTPUT/$SAMPLE.maternal.GRCh38_no_alt.bam
  gsutil -u $GCS_REQUESTER_PAYS_PROJECT cp gs://fc-4310e737-a388-4a10-8c9e-babe06aaf0cf/working/HPRC/$SAMPLE/assemblies/year1_f1_assembly_v2_genbank/alignment/assembly-to-reference/$SAMPLE.maternal.GRCh38_no_alt.bam.bai $OUTPUT/$SAMPLE.maternal.GRCh38_no_alt.bam.bai
fi

if [ ! -f $OUTPUT/$SAMPLE.maternal.GRCh38_no_alt.$TEST_NAME.fasta ]; then
  ./target/release/hidive fetch \
    -l $TEST_LOCUS \
    -o $OUTPUT/$SAMPLE.maternal.GRCh38_no_alt.$TEST_NAME.fasta \
    $OUTPUT/$SAMPLE.maternal.GRCh38_no_alt.bam
fi

# paternal haplotype train
if [ ! -f $OUTPUT/$SAMPLE.paternal.GRCh38_no_alt.bam ]; then
  gsutil -u $GCS_REQUESTER_PAYS_PROJECT cp gs://fc-4310e737-a388-4a10-8c9e-babe06aaf0cf/working/HPRC/$SAMPLE/assemblies/year1_f1_assembly_v2_genbank/alignment/assembly-to-reference/$SAMPLE.paternal.GRCh38_no_alt.bam $OUTPUT/$SAMPLE.paternal.GRCh38_no_alt.bam
  gsutil -u $GCS_REQUESTER_PAYS_PROJECT cp gs://fc-4310e737-a388-4a10-8c9e-babe06aaf0cf/working/HPRC/$SAMPLE/assemblies/year1_f1_assembly_v2_genbank/alignment/assembly-to-reference/$SAMPLE.paternal.GRCh38_no_alt.bam.bai $OUTPUT/$SAMPLE.paternal.GRCh38_no_alt.bam.bai
fi

if [ ! -f $OUTPUT/$SAMPLE.paternal.GRCh38_no_alt.$TEST_NAME.fasta ]; then
  ./target/release/hidive fetch \
    -l $TEST_LOCUS \
    -o $OUTPUT/$SAMPLE.paternal.GRCh38_no_alt.$TEST_NAME.fasta \
    $OUTPUT/$SAMPLE.paternal.GRCh38_no_alt.bam
fi

# fetch long reads testing
if [ ! -f $OUTPUT/$SAMPLE.LR.$TEST_NAME.fasta ]; then
  ./target/release/hidive fetch \
    -l $TEST_LOCUS \
    -o $OUTPUT/$SAMPLE.LR.$TEST_NAME.fasta \
    gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/HPRC_grch38/PBCCSWholeGenome/$SAMPLE/alignments/$SAMPLE.bam
#    gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/HPRC_grch38/PBCCSWholeGenome/$SAMPLE/alignments/$SAMPLE.bam
fi

# fetch short reads testing
if [ ! -f $OUTPUT/$SAMPLE.SR.$TEST_NAME.fasta ]; then
  ./target/release/hidive fetch \
    -l $TEST_LOCUS \
    -o $OUTPUT/$SAMPLE.SR.$TEST_NAME.fasta \
    gs://fc-1ee08173-e353-4494-ad28-7a3d7bd99734/working/HPRC/$SAMPLE/raw_data/Illumina/child/$SAMPLE.final.cram
#    gs://fc-1ee08173-e353-4494-ad28-7a3d7bd99734/working/HPRC/$SAMPLE/raw_data/Illumina/child/$SAMPLE.final.cram

fi

# eval

./target/release/hidive eval-model \
	-o $OUTPUT/training.long_and_short_reads.$TEST_NAME.json \
	--long-read-seq-paths $OUTPUT/$SAMPLE.LR.$TEST_NAME.fasta \
	--short-read-seq-paths $OUTPUT/$SAMPLE.SR.$TEST_NAME.fasta \
  --truth-seq-paths	$OUTPUT/$SAMPLE.maternal.GRCh38_no_alt.$TEST_NAME.fasta $OUTPUT/$SAMPLE.paternal.GRCh38_no_alt.$TEST_NAME.fasta \
  --model-path $MODEL_PATH

#./target/release/hidive eval-model \
#-o scratch/training/combined_sample_model/copy.training.long_and_short_reads.KIR.json \
#--long-read-seq-paths scratch/training/testtrain/KIR/HG00438.LR.HLA-A.fasta \
#--short-read-seq-paths scratch/training/testtrain/KIR/HG00438.SR.HLA-A.fasta \
#--truth-seq-paths scratch/training/testtrain/KIR/HG00438.maternal.test.GRCh38_no_alt.HLA-A.fasta scratch/training/testtrain/KIR/HG00438.paternal.test.GRCh38_no_alt.HLA-A.fasta \
#--model-path scratch/training/combined_sample_model/training.long_and_short_reads.KIR.json


./target/release/hidive eval-model \
-o scratch/training/combined_sample_model/MHC/copy.training.long_and_short_reads.MHC.json \
--long-read-seq-paths scratch/training/per_sample_model/HG00438/KIR/HG00438.LR.HLA-A.fasta \
--short-read-seq-paths scratch/training/per_sample_model/HG00438/KIR/HG00438.SR.HLA-A.fasta \
--truth-seq-paths scratch/training/per_sample_model/HG00438/KIR/HG00438.maternal.test.GRCh38_no_alt.HLA-A.fasta scratch/training/per_sample_model/HG00438/KIR/HG00438.paternal.test.GRCh38_no_alt.HLA-A.fasta \
--model-path scratch/training/combined_sample_model/MHC/training.long_and_short_reads.MHC.json