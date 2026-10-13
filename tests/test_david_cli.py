"""CLI isolation tests for the DAVID native-export handoff."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch

import pandas as pd

REPO_ROOT = Path(__file__).parents[1]
SCRIPT = REPO_ROOT / "scripts" / "run_pathway_workflow.py"
sys.path.insert(0, str(REPO_ROOT / "src"))

spec = importlib.util.spec_from_file_location("run_pathway_workflow", SCRIPT)
assert spec and spec.loader
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)



def _full_workflow(events: list[str]) -> MagicMock:
    workflow = MagicMock()
    normalized = pd.DataFrame(
        [{
            "SYMBOL": "GENE",
            "TOTAL_PASSING_CANDIDATE_VARIANTS": 1,
            "TOTAL_VARIANTS_EVALUATED": 1,
            "IN_FOREGROUND": True,
        }]
    )
    mapping = pd.DataFrame(
        [{"SYMBOL": "GENE", "HGNC_ID": "1", "ENSEMBL_ID": "ENSG000001", "STATUS": "MAPPED"}]
    )
    empty = pd.DataFrame(columns=["SOURCE", "TERM_ID"])
    workflow.run_p0.return_value = (normalized, {"status": "completed", "count": 1})
    workflow.mapping_from_gprofiler.return_value = mapping
    workflow.run_p1.return_value = empty
    workflow.run_p2.return_value = empty
    workflow.run_p3.return_value = empty
    workflow.run_p4.return_value = empty
    workflow.run_p5.return_value = empty
    workflow.run_p6.side_effect = lambda *_: events.append("p6") or {"status": "completed", "count": 0}
    workflow.run_p7.side_effect = lambda *_: events.append("p7")
    workflow.stage_summary.side_effect = lambda frame, _: {"status": "completed", "count": len(frame)}
    workflow.provenance.side_effect = lambda stages: {"stages": dict(stages)}
    return workflow
class DavidCliTests(unittest.TestCase):
    def test_david_requires_no_legacy_gene_input_or_provider_workflow(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            baseline = root / "baseline"
            baseline.mkdir()
            output = root / "external-output"
            with patch.object(runner, "run_david_addition", return_value=2) as david, patch.object(
                runner, "PathwayWorkflow", side_effect=AssertionError("legacy workflow must not start")
            ):
                with patch.object(sys, "argv", [str(SCRIPT), "--david", "--reuse-pathway-dir", str(baseline), "--output-dir", str(output)]):
                    self.assertEqual(runner.main(), 2)
            david.assert_called_once_with(baseline, output, None, force=False)
            self.assertFalse(output.exists())

    def test_david_rejects_non_david_stage(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with patch.object(sys, "argv", [str(SCRIPT), "--david", "--stage", "p1", "--reuse-pathway-dir", str(root), "--output-dir", str(root / "output")]):
                with self.assertRaisesRegex(ValueError, "stage"):
                    runner.main()


    def test_with_david_runs_after_p6_and_records_only_public_summary(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture_dir = root / "fixture"
            fixture_dir.mkdir()
            (fixture_dir / "gprofiler.json").write_text("{}", encoding="utf-8")
            output_dir = root / "pathway"
            david_output_dir = root / "david"
            events: list[str] = []
            workflow = _full_workflow(events)
            summary = {
                "status": "completed",
                "run_key": "run-key",
                "import_version": "receipt-hash",
                "external_output_dir": str(david_output_dir),
                "private": "must not persist",
            }
            with (
                patch.object(runner, "PathwayWorkflow", return_value=workflow),
                patch.object(runner, "run_david_addition", side_effect=lambda *_args, **_kwargs: events.append("david") or 0) as david,
                patch.object(runner, "david_run_summary", return_value=summary) as david_summary,
                patch.object(sys, "argv", [
                    str(SCRIPT),
                    "--gene-mechanism", str(root / "gene_mechanism.csv"),
                    "--offline-fixture", str(fixture_dir),
                    "--with-david",
                    "--david-output-dir", str(david_output_dir),
                    "--output-dir", str(output_dir),
                ]),
            ):
                self.assertEqual(runner.main(), 0)

            self.assertEqual(events, ["p6", "david", "p7"])
            david.assert_called_once_with(output_dir, david_output_dir, None, force=False)
            david_summary.assert_called_once_with(david_output_dir, 0)
            self.assertEqual(
                workflow.run_p7.call_args.args[1]["david"],
                {key: value for key, value in summary.items() if key != "private"},
            )
            self.assertEqual(
                workflow.provenance.call_args.args[0]["david"],
                {key: value for key, value in summary.items() if key != "private"},
            )

    def test_with_david_import_error_returns_combined_partial_exit(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture_dir = root / "fixture"
            fixture_dir.mkdir()
            (fixture_dir / "gprofiler.json").write_text("{}", encoding="utf-8")
            workflow = _full_workflow([])
            with (
                patch.object(runner, "PathwayWorkflow", return_value=workflow),
                patch.object(runner, "run_david_addition", side_effect=ValueError("INVALID_RECEIPT")),
                patch.object(runner, "david_run_summary", return_value={"status": "failed"}) as david_summary,
                patch.object(sys, "argv", [
                    str(SCRIPT),
                    "--gene-mechanism", str(root / "gene_mechanism.csv"),
                    "--offline-fixture", str(fixture_dir),
                    "--with-david",
                    "--david-output-dir", str(root / "david"),
                    "--output-dir", str(root / "pathway"),
                ]),
            ):
                self.assertEqual(runner.main(), 2)

            david_summary.assert_called_once_with(root / "david", 2, reason="INVALID_RECEIPT")

    def test_with_david_validates_flags_before_workflow_construction(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with (
                patch.object(runner, "PathwayWorkflow", side_effect=AssertionError("workflow constructed")),
                patch.object(sys, "argv", [
                    str(SCRIPT),
                    "--with-david",
                    "--stage", "p6",
                    "--david-output-dir", str(root / "david"),
                    "--output-dir", str(root / "pathway"),
                ]),
            ):
                with self.assertRaisesRegex(ValueError, "stage"):
                    runner.main()

    def test_non_david_rejects_reused_pathway_directory_before_workflow_construction(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with (
                patch.object(runner, "PathwayWorkflow", side_effect=AssertionError("workflow constructed")),
                patch.object(sys, "argv", [
                    str(SCRIPT),
                    "--gene-mechanism", str(root / "gene_mechanism.csv"),
                    "--reuse-pathway-dir", str(root / "baseline"),
                    "--output-dir", str(root / "pathway"),
                ]),
            ):
                with self.assertRaisesRegex(ValueError, "reuse-pathway-dir"):
                    runner.main()

    def test_p7_renders_compact_david_addendum(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            report_path = Path(temporary) / "p7.md"
            runner.PathwayWorkflow().run_p7(
                report_path,
                {
                    "p0": {"status": "completed", "foreground_count": 1, "background_count": 2},
                    "david": {
                        "status": "completed",
                        "run_key": "run-key",
                        "import_version": "receipt-hash",
                        "external_output_dir": "/external/david",
                        "private": "must not render",
                    },
                },
            )
            report = report_path.read_text(encoding="utf-8")

        self.assertIn("## DAVID addendum", report)
        self.assertIn("- Import version: receipt-hash", report)
        self.assertNotIn("must not render", report)

if __name__ == "__main__":
    unittest.main()
