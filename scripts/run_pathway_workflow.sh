#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec uv run python "${ROOT}/scripts/run_pathway_workflow.py" \
  --gene-mechanism "${ROOT}/data/processed/track2/02_phase1_gene_mechanism_matrix.csv" \
  --network-context "${ROOT}/data/processed/track2/07_phase2_network_gene_context.csv" \
  --network-panel "${ROOT}/data/processed/track2/18_phase3_query_gene_panel.csv" \
  --network-drugs "${ROOT}/data/processed/track2/14_phase3_drug_evidence_matrix.csv" \
  --output-dir "${ROOT}/data/processed/track2/pathway" "$@"
