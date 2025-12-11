#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DATA_DIR="${REPO_ROOT}/examples/dummy_data"
OUT_DIR="${1:-${DATA_DIR}/output}"

GFA_PATH="${OUT_DIR}/pangenome.gfa"
MODEL_PATH="${OUT_DIR}/crf_model.json"
HAP_PATH="${OUT_DIR}/haplotypes.fa"

mkdir -p "${OUT_DIR}"

if ! command -v seqwish >/dev/null 2>&1; then
  echo "Error: seqwish is not available in PATH. Install it before running this script." >&2
  exit 1
fi

echo "=== Building dummy pangenome graph ==="
cargo run --bin hidive -- build-pangenome \
  --tier1-fasta-paths "${DATA_DIR}/tier1.fa" \
  --tier2-fasta-paths "${DATA_DIR}/tier2.fa" \
  --tier3-fasta-paths "${DATA_DIR}/tier3.fa" \
  --output "${GFA_PATH}" \
  --kmer-size 11 \
  --min-aln-len 30

echo "=== Training CRF on dummy data ==="
cargo run --bin hidive -- train-crf \
  --graph "${GFA_PATH}" \
  --reads "${DATA_DIR}/reads.fa" \
  --truth-haplotypes "${DATA_DIR}/truth_hap1.fa" \
  --truth-haplotypes "${DATA_DIR}/truth_hap2.fa" \
  --output "${MODEL_PATH}" \
  --kmer-size 11 \
  --iterations 5

echo "=== Inferring haplotypes on dummy data ==="
cargo run --bin hidive -- infer-haplotypes \
  --graph "${GFA_PATH}" \
  --model "${MODEL_PATH}" \
  --reads "${DATA_DIR}/reads.fa" \
  --output "${HAP_PATH}" \
  --kmer-size 11

echo "Dummy pipeline complete."
echo "GFA:      ${GFA_PATH}"
echo "Model:    ${MODEL_PATH}"
echo "Haplotypes: ${HAP_PATH}"

