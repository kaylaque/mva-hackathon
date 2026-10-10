"""Unit tests for Track 2 artifact migration script using synthetic non-patient fixtures."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile

REPO_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import pytest
from migrate_track2_layout import (
    TRACK2_ARTIFACT_MAPPINGS,
    check_sources_exist,
    migrate_track2_artifacts,
)


def create_synthetic_track2_tree(base_dir: Path) -> dict[str, int]:
    """Create synthetic empty/small files matching the expected Track 2 structure."""
    expected_sizes = {}
    for mapping in TRACK2_ARTIFACT_MAPPINGS:
        src = base_dir / mapping.source_rel
        src.parent.mkdir(parents=True, exist_ok=True)
        if src.suffix == ".json":
            content = json.dumps({"test_stage": mapping.source_rel, "synthetic": True})
        else:
            content = "GENE,SCORE,MOCK\nTEST,1.0,SYNTHETIC\n"
        src.write_text(content, encoding="utf-8")
        expected_sizes[mapping.source_rel] = src.stat().st_size
    return expected_sizes


def test_check_sources_exist_complete():
    with tempfile.TemporaryDirectory() as tmp_dir:
        base_dir = Path(tmp_dir)
        create_synthetic_track2_tree(base_dir)
        all_present, missing = check_sources_exist(base_dir)
        assert all_present is True
        assert len(missing) == 0


def test_check_sources_exist_missing():
    with tempfile.TemporaryDirectory() as tmp_dir:
        base_dir = Path(tmp_dir)
        # Empty directory
        all_present, missing = check_sources_exist(base_dir)
        assert all_present is False
        expected_missing = [m for m in TRACK2_ARTIFACT_MAPPINGS if not m.optional]
        assert len(missing) == len(expected_missing)


def test_migration_dry_run_leaves_files_untouched():
    with tempfile.TemporaryDirectory() as tmp_dir:
        base_dir = Path(tmp_dir)
        create_synthetic_track2_tree(base_dir)

        results = migrate_track2_artifacts(base_dir=base_dir, dry_run=True)
        assert len(results) == len(TRACK2_ARTIFACT_MAPPINGS)
        for mapping in TRACK2_ARTIFACT_MAPPINGS:
            dst = base_dir / mapping.target_rel
            assert not dst.exists(), f"Target {dst} should not exist after dry run"


def test_migration_copy_and_validation():
    with tempfile.TemporaryDirectory() as tmp_dir:
        base_dir = Path(tmp_dir)
        sizes = create_synthetic_track2_tree(base_dir)

        results = migrate_track2_artifacts(base_dir=base_dir, dry_run=False)
        assert len(results) == len(TRACK2_ARTIFACT_MAPPINGS)

        for mapping in TRACK2_ARTIFACT_MAPPINGS:
            src = base_dir / mapping.source_rel
            dst = base_dir / mapping.target_rel
            assert src.is_file(), f"Source {src} must remain intact (duplicates kept)"
            assert dst.is_file(), f"Target {dst} must exist after copy"
            assert src.stat().st_size == dst.stat().st_size
            assert results[mapping.target_rel]["action"] == "COPIED"

        # Now test validate_only mode
        val_results = migrate_track2_artifacts(base_dir=base_dir, validate_only=True)
        assert all(v["valid"] is True for v in val_results.values())
        assert all(v["content_match"] is True for v in val_results.values())


def test_validation_rejects_same_size_content_mismatch():
    with tempfile.TemporaryDirectory() as tmp_dir:
        base_dir = Path(tmp_dir)
        create_synthetic_track2_tree(base_dir)
        migrate_track2_artifacts(base_dir=base_dir)

        target = base_dir / "00_phase0_consolidated_candidates.csv"
        source = base_dir / "phase0/consolidated_candidates.csv"
        original = target.read_bytes()
        target.write_bytes(b"X" + original[1:])
        assert target.stat().st_size == source.stat().st_size

        val_results = migrate_track2_artifacts(base_dir=base_dir, validate_only=True)
        assert val_results[target.name]["size_match"] is True
        assert val_results[target.name]["content_match"] is False
        assert val_results[target.name]["valid"] is False


def test_migration_fails_loudly_on_missing_source():
    with tempfile.TemporaryDirectory() as tmp_dir:
        base_dir = Path(tmp_dir)
        create_synthetic_track2_tree(base_dir)

        # Remove one source file
        (base_dir / "phase2/string_edges.csv").unlink()

        with pytest.raises(FileNotFoundError, match="missing required source files"):
            migrate_track2_artifacts(base_dir=base_dir, dry_run=False)
