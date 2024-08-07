#!/bin/bash

set -euxo pipefail

#LOCUS="chr6:25,726,063-33,410,226"
#NAME="MHC"

#LOCUS="chr6:29,941,260-29,949,572"
#NAME="HLA-A"

#LOCUS="chr6:31,352,872-31,368,067"
#NAME="HLA-B"

LOCUS="chr6:32,438,878-32,446,046"
NAME="HLA-DRA"

OUTPUT="scratch/$NAME"

cargo build --release

mkdir -p $OUTPUT

# maternal haplotype

./target/release/hidive fetch \
	-l $LOCUS \
	-o $OUTPUT/HG00438.maternal.GRCh38_no_alt.$NAME.fasta \
	gs://fc-4310e737-a388-4a10-8c9e-babe06aaf0cf/working/HPRC/HG00438/assemblies/year1_f1_assembly_v2_genbank/alignment/assembly-to-reference/HG00438.maternal.GRCh38_no_alt.bam

minimap2 -ayYL --MD --eqx -x map-hifi -R '@RG\tID:maternal\tSM:maternal' $OUTPUT/HG00438.maternal.GRCh38_no_alt.$NAME.fasta $OUTPUT/HG00438.maternal.GRCh38_no_alt.$NAME.fasta | \
	samtools sort -O BAM --write-index -o $OUTPUT/HG00438.maternal.GRCh38_no_alt.$NAME.bam

# paternal haplotype

./target/release/hidive fetch \
	-l $LOCUS \
	gs://fc-4310e737-a388-4a10-8c9e-babe06aaf0cf/working/HPRC/HG00438/assemblies/year1_f1_assembly_v2_genbank/alignment/assembly-to-reference/HG00438.paternal.GRCh38_no_alt.bam | \
	minimap2 -ayYL --MD --eqx -x map-hifi -R '@RG\tID:paternal\tSM:paternal' $OUTPUT/HG00438.maternal.GRCh38_no_alt.$NAME.fasta - | \
	samtools sort -O BAM --write-index -o $OUTPUT/HG00438.paternal.GRCh38_no_alt.$NAME.bam

# fetch

./target/release/hidive fetch \
	-l $LOCUS \
	-o $OUTPUT/HG00438.LR.$NAME.fasta \
	gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/HPRC_grch38/PBCCSWholeGenome/HG00438/alignments/HG00438.bam

minimap2 -ayYL --MD --eqx -x map-hifi -R '@RG\tID:lr\tSM:lr' $OUTPUT/HG00438.maternal.GRCh38_no_alt.$NAME.fasta $OUTPUT/HG00438.LR.$NAME.fasta | \
	samtools sort -O BAM --write-index -o $OUTPUT/HG00438.LR.$NAME.bam

# rescue

./target/release/hidive rescue \
	-o $OUTPUT/HG00438.SR.$NAME.fasta \
	-f $OUTPUT/HG00438.LR.$NAME.fasta \
	gs://fc-1ee08173-e353-4494-ad28-7a3d7bd99734/working/HPRC/HG00438/raw_data/Illumina/child/HG00438.final.cram

minimap2 -ayYL --MD --eqx -x sr -R '@RG\tID:sr\tSM:sr' $OUTPUT/HG00438.maternal.GRCh38_no_alt.$NAME.fasta $OUTPUT/HG00438.SR.$NAME.fasta | \
	samtools sort -O BAM --write-index -o $OUTPUT/HG00438.SR.$NAME.bam

# filter

./target/release/hidive filter \
	-o $OUTPUT/HG00438.SR.$NAME.filtered.fasta \
	-l $OUTPUT/HG00438.LR.$NAME.fasta \
	$OUTPUT/HG00438.SR.$NAME.fasta

minimap2 -ayYL --MD --eqx -x sr -R '@RG\tID:sr-filt\tSM:sr-filt' $OUTPUT/HG00438.maternal.GRCh38_no_alt.$NAME.fasta $OUTPUT/HG00438.SR.$NAME.filtered.fasta | \
	samtools sort -O BAM --write-index -o $OUTPUT/HG00438.SR.$NAME.filtered.bam

# co-assemble

./target/release/hidive coassemble -m training.json -s $OUTPUT/HG00438.SR.$NAME.filtered.fasta -o $OUTPUT/coassemble.fasta $OUTPUT/HG00438.LR.$NAME.fasta

# align
cat $OUTPUT/coassemble.fasta | \
	minimap2 -ayYL --MD --eqx -x map-hifi -R '@RG\tID:asm\tSM:asm' $OUTPUT/HG00438.maternal.GRCh38_no_alt.$NAME.fasta - | \
	samtools sort -O BAM --write-index -o $OUTPUT/coassemble.bam

