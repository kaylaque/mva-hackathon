#!/usr/bin/env bash
# Thin launcher: runs only when explicitly invoked by the researcher.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec uv run python "${ROOT}/scripts/run_enrichment.py" --root "${ROOT}" "$@"
