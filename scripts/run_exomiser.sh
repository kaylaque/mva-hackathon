#!/usr/bin/env bash
# Thin launcher: prepare by default; Exomiser executes only with --execute.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec uv run python "${ROOT}/scripts/run_exomiser.py" --root "${ROOT}" "$@"
