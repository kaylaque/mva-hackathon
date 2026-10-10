#!/usr/bin/env bash
# Thin launcher for Phase 0 Candidate Consolidation.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec uv run python "${ROOT}/scripts/consolidate_candidates.py" "$@"
