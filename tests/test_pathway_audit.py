"""Synthetic tests for aggregate-only pathway validation reporting."""

from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd

REPO_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from mva_hackathon.pathway_audit import audit_pathway_outputs


class PathwayAuditTests(unittest.TestCase):
    def test_audit_explains_null_baseline_and_threshold_two_signal(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            pd.DataFrame(
                [
                    {"STATUS": "mapped", "IN_FOREGROUND": True},
                    {"STATUS": "unmapped", "IN_FOREGROUND": False},
                ]
            ).to_csv(output_dir / "p0_gene_qc.csv", index=False)
            pd.DataFrame(
                [
                    {"SOURCE": "REAC", "RAW_P": 0.01, "BH_Q": 0.2},
                    {"SOURCE": "HP", "RAW_P": 0.2, "BH_Q": 0.3},
                ]
            ).to_csv(output_dir / "p1_enrichment_results.csv", index=False)
            pd.DataFrame(columns=["TERM_ID", "BH_Q"]).to_csv(
                output_dir / "p2_prioritized_pathways.csv", index=False
            )
            pd.DataFrame(
                [{"SOURCE": "REAC", "BH_Q": 0.01, "TERM_ID": "R-HSA-1"}]
            ).to_csv(output_dir / "p2_sensitivity_min_passing_2.csv", index=False)
            pd.DataFrame(columns=["TERM_ID", "BH_Q"]).to_csv(
                output_dir / "p2_sensitivity_hp_excluded.csv", index=False
            )
            pd.DataFrame(
                [{"shared_targets": 0, "pathway_only_targets": 0, "string_only_targets": 3}]
            ).to_csv(output_dir / "p6_evaluation_summary.csv", index=False)

            summary, report_path = audit_pathway_outputs(output_dir)
            report = report_path.read_text()

        self.assertEqual(summary["baseline"]["provider_bh_le_0_05"], 0)
        self.assertEqual(summary["threshold_two"]["prioritized_pathways"], 1)
        self.assertIn("No baseline term passed provider BH-FDR", report)
        self.assertIn("Threshold-2 sensitivity produced 1 prioritized pathway", report)
        self.assertNotIn("VARIANT_KEY", report)

        self.assertIn("## Threshold-2 prioritized pathways", report)
        self.assertIn("NOT_APPLICABLE_NO_ELIGIBLE_INPUT", report)
        self.assertIn("Drug comparison: UNAVAILABLE", report)


    def test_audit_does_not_mark_nonempty_baseline_as_not_applicable(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            pd.DataFrame([{"STATUS": "mapped", "IN_FOREGROUND": True}]).to_csv(
                output_dir / "p0_gene_qc.csv", index=False
            )
            pd.DataFrame([{"SOURCE": "REAC", "RAW_P": 0.01, "BH_Q": 0.01}]).to_csv(
                output_dir / "p1_enrichment_results.csv", index=False
            )
            pd.DataFrame([{"TERM_ID": "R-HSA-1", "BH_Q": 0.01}]).to_csv(
                output_dir / "p2_prioritized_pathways.csv", index=False
            )
            pd.DataFrame([{"shared_targets": 1, "pathway_only_targets": 0, "string_only_targets": 0}]).to_csv(
                output_dir / "p6_evaluation_summary.csv", index=False
            )

            _, report_path = audit_pathway_outputs(output_dir)
            report = report_path.read_text()

        self.assertNotIn("P3 mechanism, P4 target expansion, and P5 drug validation: NOT_APPLICABLE", report)
if __name__ == "__main__":
    unittest.main()
