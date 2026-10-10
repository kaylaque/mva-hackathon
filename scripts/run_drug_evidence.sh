#!/usr/bin/env bash
# ==============================================================================
# Phase 3: Multi-Source Drug Repurposing Evidence and Safety Review
# ==============================================================================
# Usage:
#   bash scripts/run_drug_evidence.sh [OPTIONS]
#
# Default behavior uses Track 2 flat layout paths:
#   Input gene mechanism:  data/processed/track2/02_phase1_gene_mechanism_matrix.csv
#   Input network context: data/processed/track2/07_phase2_network_gene_context.csv
#   Output directory:      data/processed/track2/
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

GENE_MECHANISM="${1:-"${REPO_ROOT}/data/processed/track2/02_phase1_gene_mechanism_matrix.csv"}"
NETWORK_CONTEXT="${2:-"${REPO_ROOT}/data/processed/track2/07_phase2_network_gene_context.csv"}"
OUTPUT_DIR="${3:-"${REPO_ROOT}/data/processed/track2"}"

shift 3 2>/dev/null || true

echo "=== Running Phase 3 Drug Repurposing Evidence Workflow ==="
echo "Gene Mechanism Matrix: ${GENE_MECHANISM}"
echo "Network Gene Context:  ${NETWORK_CONTEXT}"
echo "Output Directory:      ${OUTPUT_DIR}"

uv run python "${REPO_ROOT}/scripts/run_drug_evidence.py" \
    --gene-mechanism "${GENE_MECHANISM}" \
    --network-context "${NETWORK_CONTEXT}" \
    --output-dir "${OUTPUT_DIR}" \
    --flat-names \
    "$@"
