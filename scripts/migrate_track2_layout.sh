#!/usr/bin/env bash
# Thin launcher for Track 2 Artifact Migration.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec uv run python "${ROOT}/scripts/migrate_track2_layout.py" "$@"
