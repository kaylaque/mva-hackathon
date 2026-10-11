#!/usr/bin/env bash
# Phase 4 Signed-Pathway/Reaction Adjudication & Literature Validation Runner
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

echo "==> Running Phase 4 Adjudication Engine..."
uv run python "${REPO_ROOT}/scripts/run_phase4_adjudication.py" "$@"
