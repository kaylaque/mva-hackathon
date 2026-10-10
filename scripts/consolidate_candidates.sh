#!/usr/bin/env bash
# Thin launcher for Phase 0 Candidate Consolidation.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Support default flat output directory if no arguments provided
if [ "$#" -eq 0 ]; then
  exec uv run python "${ROOT}/scripts/consolidate_candidates.py" \
    --evidence "${ROOT}/data/processed/variant_evidence_master.csv" \
    --exomiser-variants "${ROOT}/data/interim/exomiser/mva_patient_exomiser.variants.tsv" \
    --exomiser-genes "${ROOT}/data/interim/exomiser/mva_patient_exomiser.genes.tsv" \
    --output-candidates "${ROOT}/data/processed/track2/00_phase0_consolidated_candidates.csv" \
    --output-provenance "${ROOT}/data/processed/track2/00_phase0_provenance.json"
fi

exec uv run python "${ROOT}/scripts/consolidate_candidates.py" "$@"
