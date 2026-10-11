#!/usr/bin/env bash
# ==============================================================================
# Phase 3: Query Gene Panel Generation Engine
# ==============================================================================
# Usage:
#   bash scripts/run_query_gene_panel.sh [TRACK2_DIR] [OUTPUT_DIR] [ADDITIONAL_ARGS]
#
# Default behavior uses Track 2 flat layout paths:
#   Track 2 input dir:  data/processed/track2/
#   Output directory:   data/processed/track2/
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

TRACK2_DIR="${1:-"${REPO_ROOT}/data/processed/track2"}"
OUTPUT_DIR="${2:-"${REPO_ROOT}/data/processed/track2"}"

shift 2 2>/dev/null || true

echo "=== Running Phase 3 Query Gene Panel Generation ==="
echo "Track 2 Input Directory: ${TRACK2_DIR}"
echo "Output Directory:        ${OUTPUT_DIR}"

uv run python "${REPO_ROOT}/scripts/run_query_gene_panel.py" \
    --track2-dir "${TRACK2_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --flat-names \
    "$@"
