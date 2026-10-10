#!/usr/bin/env bash
# Thin launcher for Phase 1 Mechanism & Directionality Matrix Generation.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Support default flat output directory if no arguments provided
if [ "$#" -eq 0 ]; then
  exec uv run python "${ROOT}/scripts/build_mechanism_matrix.py" \
    --candidates "${ROOT}/data/processed/track2/00_phase0_consolidated_candidates.csv" \
    --output-variant-matrix "${ROOT}/data/processed/track2/01_phase1_variant_mechanism_matrix.csv" \
    --output-gene-matrix "${ROOT}/data/processed/track2/02_phase1_gene_mechanism_matrix.csv" \
    --output-provenance "${ROOT}/data/processed/track2/02_phase1_provenance.json"
fi

exec uv run python "${ROOT}/scripts/build_mechanism_matrix.py" "$@"
