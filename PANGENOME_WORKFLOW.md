# Pangenome Graph Workflow

This document describes the workflow for building and using pangenome graphs with hidive.

## Overview

hidive uses GFA (Graphical Fragment Assembly) format for pangenome graphs. The `build-pangenome` command outputs GFA format, which can be loaded directly by `train-crf` and `infer-haplotypes` commands. No conversion step is needed.

## Workflow

### Step 1: Build Pangenome Graph (GFA)

```bash
hidive build-pangenome \
    --tier1-fasta-path tier1_assembly1.fa \
    --tier1-fasta-path tier1_assembly2.fa \
    --tier2-fasta-path tier2_assembly1.fa \
    --tier3-fasta-path tier3_assembly1.fa \
    --output pangenome.gfa \
    --kmer-size 17 \
    --min-aln-len 100
```

This will:
- Combine all FASTA files with tier annotations in sequence names
- Build a variation graph using `seqwish`
- Output a GFA file with tier information encoded in path names

**Tier annotation format**: Path names will contain `|tier=TierN` (e.g., `sample1|tier=Tier1`)

### Step 2: Train CRF Model

```bash
hidive train-crf \
    --graph pangenome.gfa \
    --reads sample_reads.fq \
    --truth-haplotypes truth_hap1.fa \
    --truth-haplotypes truth_hap2.fa \
    --output crf_model.json \
    --kmer-size 17
```

### Step 3: Infer Haplotypes

```bash
hidive infer-haplotypes \
    --graph pangenome.gfa \
    --model crf_model.json \
    --reads sample_reads.fq \
    --output haplotypes.fa \
    --kmer-size 17
```

## Why GFA Format?

1. **Simplicity**: GFA is a text-based format that's easy to parse and debug
2. **No External Dependencies**: No need for ODGI conversion or command-line tools
3. **Direct Loading**: GFA files can be loaded directly using the `parfait-gfa` crate
4. **Standard Format**: GFA is a widely-used standard for variation graphs

## Tier Information

Tier information is preserved in path names during graph construction:
- **Tier 1**: `path_name|tier=Tier1` (weight: 1.0)
- **Tier 2**: `path_name|tier=Tier2` (weight: 0.7)
- **Tier 3**: `path_name|tier=Tier3` (weight: 0.4)

The CRF model will use these tier weights as features during training and inference.

## Dummy Test Dataset

A small end-to-end dataset lives in `examples/dummy_data/`. To run the entire pipeline on it:

```bash
chmod +x scripts/run_dummy_pipeline.sh
scripts/run_dummy_pipeline.sh  # optional argument: output directory
```

The script will:
- Build a GFA graph from the dummy tiered FASTAs
- Train the placeholder CRF with the bundled reads and truth haplotypes
- Infer haplotypes and write them to `examples/dummy_data/output/haplotypes.fa`

## Dependencies

- `seqwish`: Used by `build-pangenome` to construct variation graphs (included as Rust dependency)
- `parfait-gfa`: Used for parsing GFA files (included as Rust dependency)

