#!/usr/bin/env bash
# Thin launcher for Phase 1 Mechanism & Directionality Matrix Generation.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec uv run python "${ROOT}/scripts/build_mechanism_matrix.py" "$@"
