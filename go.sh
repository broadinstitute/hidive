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
	-o $OUTPUT/HG002.maternal.GRCh38_no_alt.$NAME.fasta \
	gs://fc-4310e737-a388-4a10-8c9e-babe06aaf0cf/working/HPRC_PLUS/HG002/assemblies/year1_f1_assembly_v2_genbank/alignment/assembly-to-reference/HG002.maternal.CHM13Y_EBV.bam

minimap2 -ayYL --MD --eqx -x map-hifi -R '@RG\tID:maternal\tSM:maternal' /Users/suhang/Analysis/hidive/MHC-HG002.1.fa $OUTPUT/HG00438.maternal.GRCh38_no_alt.$NAME.fasta | \
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
	-o $OUTPUT/HG002.LR.$NAME.fasta \
	gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/HPRC_grch38/HG002/HG002.bam

minimap2 -ayYL --MD --eqx -x map-hifi -R '@RG\tID:lr\tSM:lr' /Users/suhang/Analysis/hidive/MHC-HG002.1.fa $OUTPUT/HG002.LR.$NAME.fasta | \
	samtools sort -O BAM --write-index -o $OUTPUT/HG002.LR.$NAME.bam

# rescue

./target/release/hidive rescue \
	-r /Users/suhang/Analysis/hidive_backup/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
	-o $OUTPUT/HG002.SR.mrna.$NAME.fasta \
	-f $OUTPUT/HG002.LR.$NAME.fasta \
	/Users/suhang/Analysis/hg002_gm24385.mrna.grch38.cram

minimap2 -ayYL --MD --eqx -x sr -R '@RG\tID:sr\tSM:sr' /Users/suhang/Analysis/hidive/MHC-HG002.1.fa $OUTPUT/HG002.SR.mrna.$NAME.fasta | \
	samtools sort -O BAM --write-index -o $OUTPUT/HG002.SR.mrna.$NAME.bam

# filter

./target/release/hidive filter \
	-o $OUTPUT/HG00438.SR.$NAME.filtered.fasta \
	-l $OUTPUT/HG00438.LR.$NAME.fasta \
	$OUTPUT/HG00438.SR.$NAME.fasta

minimap2 -ayYL --MD --eqx -x sr -R '@RG\tID:sr-filt\tSM:sr-filt' $OUTPUT/HG002.maternal.GRCh38_no_alt.$NAME.fasta $OUTPUT/HG002.SR.$NAME.filtered.fasta | \
	samtools sort -O BAM --write-index -o $OUTPUT/HG00438.SR.$NAME.filtered.bam

# co-assemble

./target/release/hidive coassemble -m training.json -s $OUTPUT/HG00438.SR.$NAME.filtered.fasta -o $OUTPUT/coassemble.fasta $OUTPUT/HG00438.LR.$NAME.fasta

# align
cat $OUTPUT/coassemble.fasta | \
	minimap2 -ayYL --MD --eqx -x map-hifi -R '@RG\tID:asm\tSM:asm' $OUTPUT/HG00438.maternal.GRCh38_no_alt.$NAME.fasta - | \
	samtools sort -O BAM --write-index -o $OUTPUT/coassemble.bam

