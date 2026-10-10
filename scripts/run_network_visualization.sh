#!/usr/bin/env bash
# Thin launcher for Phase 2 PPI & Pathway Network Visualization Generation.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Support default flat output directory if no arguments provided
if [ "$#" -eq 0 ]; then
  exec uv run python "${ROOT}/scripts/run_network_visualization.py" \
    --gene-context "${ROOT}/data/processed/track2/07_phase2_network_gene_context.csv" \
    --edges "${ROOT}/data/processed/track2/04_phase2_string_edges.csv" \
    --enrichment "${ROOT}/data/processed/track2/06_phase2_pathway_enrichment.csv" \
    --output-dir "${ROOT}/data/processed/track2/" \
    --flat-names
fi

exec uv run python "${ROOT}/scripts/run_network_visualization.py" "$@"
