#!/usr/bin/env bash
# Thin launcher for Phase 2 PPI & Pathway Network Context Analysis.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Support default flat output directory if no arguments provided
if [ "$#" -eq 0 ]; then
  exec uv run python "${ROOT}/scripts/run_network_analysis.py" \
    --gene-mechanism "${ROOT}/data/processed/track2/02_phase1_gene_mechanism_matrix.csv" \
    --variant-matrix "${ROOT}/data/processed/track2/01_phase1_variant_mechanism_matrix.csv" \
    --output-dir "${ROOT}/data/processed/track2/" \
    --flat-names \
    --string-score 700
fi

exec uv run python "${ROOT}/scripts/run_network_analysis.py" "$@"
