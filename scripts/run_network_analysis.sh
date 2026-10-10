#!/usr/bin/env bash
# Thin launcher for Phase 2 PPI & Pathway Network Context Analysis.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec uv run python "${ROOT}/scripts/run_network_analysis.py" "$@"
