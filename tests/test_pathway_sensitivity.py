"""Synthetic tests for the isolated pathway sensitivity workflow."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch

import pandas as pd

REPO_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from mva_hackathon.pathway_sensitivity import run_sensitivities, verify_sensitivity_artifacts


class PathwaySensitivityTests(unittest.TestCase):
    def test_preparation_writes_exactly_two_gene_only_requests_without_http(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline = root / "baseline"
            baseline.mkdir()
            pd.DataFrame(
                [
                    {"SYMBOL": "GENE1", "ENSEMBL_ID": "ENSG1", "STATUS": "mapped", "IN_FOREGROUND": "true", "PASSING_VARIANTS": "private"},
                    {"SYMBOL": "GENE2", "ENSEMBL_ID": "ENSG2", "STATUS": "mapped", "IN_FOREGROUND": "true", "PASSING_VARIANTS": "private"},
                    {"SYMBOL": "GENE3", "ENSEMBL_ID": "ENSG3", "STATUS": "mapped", "IN_FOREGROUND": "false", "PASSING_VARIANTS": "private"},
                    {"SYMBOL": "GENE4", "ENSEMBL_ID": "ENSG4", "STATUS": "mapped", "IN_FOREGROUND": "false", "PASSING_VARIANTS": "private"},
                    {"SYMBOL": "UNMAPPED", "ENSEMBL_ID": "", "STATUS": "unmapped", "IN_FOREGROUND": "true", "PASSING_VARIANTS": "private"},
                ]
            ).to_csv(baseline / "p0_gene_qc.csv", index=False)
            pd.DataFrame(columns=["TERM_ID", "TERM_NAME", "SOURCE", "RAW_P", "BH_Q", "FOLD_ENRICHMENT"]).to_csv(
                baseline / "p1_enrichment_results.csv", index=False
            )
            multihit = root / "multihit.csv"
            pd.DataFrame({"SYMBOL": ["gene1", "GENE1"]}).to_csv(multihit, index=False)
            output = root / "sensitivity"

            self.assertEqual(run_sensitivities(baseline, multihit, output), 0)

            self.assertEqual(sorted(path.name for path in output.iterdir()), ["hp_excluded", "multihit", "sensitivity_report.md"])
            multihit_request = json.loads((output / "multihit" / "request.json").read_text())
            hp_excluded_request = json.loads((output / "hp_excluded" / "request.json").read_text())
            self.assertEqual(multihit_request["payload"]["query"], ["ENSG1"])
            self.assertEqual(multihit_request["payload"]["background"], ["ENSG1", "ENSG2", "ENSG3", "ENSG4"])
            self.assertEqual(hp_excluded_request["payload"]["query"], ["ENSG1", "ENSG2"])
            self.assertNotIn("HP", hp_excluded_request["payload"]["sources"])
            self.assertNotIn("private", (output / "multihit" / "request.json").read_text())


    def test_offline_results_preserve_provider_q_and_keep_hp_as_annotation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline = root / "baseline"
            baseline.mkdir()
            pd.DataFrame(
                [
                    {"SYMBOL": "GENE1", "ENSEMBL_ID": "ENSG1", "STATUS": "mapped", "IN_FOREGROUND": "true"},
                    {"SYMBOL": "GENE2", "ENSEMBL_ID": "ENSG2", "STATUS": "mapped", "IN_FOREGROUND": "true"},
                    {"SYMBOL": "GENE3", "ENSEMBL_ID": "ENSG3", "STATUS": "mapped", "IN_FOREGROUND": "false"},
                    {"SYMBOL": "GENE4", "ENSEMBL_ID": "ENSG4", "STATUS": "mapped", "IN_FOREGROUND": "false"},
                ]
            ).to_csv(baseline / "p0_gene_qc.csv", index=False)
            pd.DataFrame(columns=["TERM_ID", "TERM_NAME", "SOURCE", "RAW_P", "BH_Q", "FOLD_ENRICHMENT"]).to_csv(
                baseline / "p1_enrichment_results.csv", index=False
            )
            multihit = root / "multihit.csv"
            pd.DataFrame({"SYMBOL": ["GENE1"]}).to_csv(multihit, index=False)
            fixture = root / "fixture"
            fixture.mkdir()
            response = {
                "meta": {"organism": "hsapiens", "domain_scope": "custom", "significance_threshold_method": "fdr"},
                "result": [
                    {"source": "REAC", "native": "R-HSA-1", "name": "Synthetic pathway", "effective_domain_size": 4, "term_size": 2, "query_size": 1, "intersection_size": 1, "intersections": ["ENSG1"], "p_value": 0.04},
                    {"source": "HP", "native": "HP:1", "name": "Synthetic annotation", "effective_domain_size": 4, "term_size": 2, "query_size": 1, "intersection_size": 1, "intersections": ["ENSG1"], "p_value": 0.06},
                    {"source": "REAC", "native": "R-HSA-2", "name": "Not significant", "effective_domain_size": 4, "term_size": 2, "query_size": 1, "intersection_size": 1, "intersections": ["ENSG1"], "p_value": 0.2},
                ],
            }
            for name in ("multihit", "hp_excluded"):
                (fixture / f"{name}.json").write_text(json.dumps(response), encoding="utf-8")
            output = root / "sensitivity"

            self.assertEqual(run_sensitivities(baseline, multihit, output, offline_fixture=fixture), 0)

            prioritized = pd.read_csv(output / "multihit" / "prioritized.csv")
            self.assertEqual(prioritized["TERM_ID"].tolist(), ["R-HSA-1"])
            enrichment = pd.read_csv(output / "multihit" / "enrichment.csv")
            annotation = enrichment.loc[enrichment["TERM_ID"].eq("HP:1")].iloc[0]
            self.assertEqual(annotation["SOURCE"], "HP")
            report = (output / "sensitivity_report.md").read_text()
            self.assertIn("Provider-adjusted q-values are not re-corrected.", report)
            self.assertIn("multihit: COMPLETED", report)
            self.assertNotIn("mechanistic pathway", (output / "sensitivity_report.md").read_text().lower())
    def test_cli_runs_preparation_and_offline_fixture(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline = root / "baseline"
            baseline.mkdir()
            pd.DataFrame([
                {"SYMBOL": "GENE1", "ENSEMBL_ID": "ENSG1", "STATUS": "mapped", "IN_FOREGROUND": "true"},
                {"SYMBOL": "GENE2", "ENSEMBL_ID": "ENSG2", "STATUS": "mapped", "IN_FOREGROUND": "false"},
            ]).to_csv(baseline / "p0_gene_qc.csv", index=False)
            pd.DataFrame(columns=["TERM_ID", "TERM_NAME", "SOURCE", "RAW_P", "BH_Q", "FOLD_ENRICHMENT"]).to_csv(
                baseline / "p1_enrichment_results.csv", index=False
            )
            multihit = root / "multihit.csv"
            pd.DataFrame({"SYMBOL": ["GENE1"]}).to_csv(multihit, index=False)
            output = root / "sensitivity"
            import subprocess

            prepared = subprocess.run(
                [sys.executable, str(REPO_ROOT / "scripts" / "run_pathway_sensitivity.py"), "--baseline-dir", str(baseline), "--multihit-genes", str(multihit), "--output-dir", str(output)],
                capture_output=True, text=True, check=False,
            )
            self.assertEqual(prepared.returncode, 0, prepared.stderr)
            self.assertIn("PREPARED_NOT_EXECUTED", prepared.stdout)
    def test_live_provider_failure_preserves_other_analysis(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline = root / "baseline"
            baseline.mkdir()
            pd.DataFrame([
                {"SYMBOL": "GENE1", "ENSEMBL_ID": "ENSG1", "STATUS": "mapped", "IN_FOREGROUND": "true"},
                {"SYMBOL": "GENE2", "ENSEMBL_ID": "ENSG2", "STATUS": "mapped", "IN_FOREGROUND": "false"},
            ]).to_csv(baseline / "p0_gene_qc.csv", index=False)
            pd.DataFrame(columns=["TERM_ID", "TERM_NAME", "SOURCE", "RAW_P", "BH_Q", "FOLD_ENRICHMENT"]).to_csv(
                baseline / "p1_enrichment_results.csv", index=False
            )
            multihit = root / "multihit.csv"
            pd.DataFrame({"SYMBOL": ["GENE1"]}).to_csv(multihit, index=False)
            response = MagicMock(status_code=200)
            response.json.return_value = {"meta": {}, "result": []}
            client = MagicMock()
            client.post.side_effect = [response, None]
            output = root / "sensitivity"

            with patch("mva_hackathon.pathway_sensitivity.HttpClient", return_value=client):
                self.assertEqual(run_sensitivities(baseline, multihit, output, execute=True), 2)

            self.assertIn("UNVERIFIED_PROVIDER_METADATA", (output / "multihit" / "manifest.json").read_text())
            self.assertIn("UNAVAILABLE_PROVIDER_RESPONSE", (output / "sensitivity_report.md").read_text())
    def test_invalid_provider_counts_are_unavailable_not_negative_results(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline = root / "baseline"
            baseline.mkdir()
            pd.DataFrame([
                {"SYMBOL": "GENE1", "ENSEMBL_ID": "ENSG1", "STATUS": "mapped", "IN_FOREGROUND": "true"},
                {"SYMBOL": "GENE2", "ENSEMBL_ID": "ENSG2", "STATUS": "mapped", "IN_FOREGROUND": "false"},
            ]).to_csv(baseline / "p0_gene_qc.csv", index=False)
            pd.DataFrame(columns=["TERM_ID", "TERM_NAME", "SOURCE", "RAW_P", "BH_Q", "FOLD_ENRICHMENT"]).to_csv(
                baseline / "p1_enrichment_results.csv", index=False
            )
            multihit = root / "multihit.csv"
            pd.DataFrame({"SYMBOL": ["GENE1"]}).to_csv(multihit, index=False)
            fixture = root / "fixture"
            fixture.mkdir()
            invalid = {"meta": {}, "result": [{"source": "REAC", "native": "R-HSA-1", "name": "Invalid", "effective_domain_size": 2, "term_size": 3, "query_size": 1, "intersection_size": 1, "p_value": 0.01}]}
            for name in ("multihit", "hp_excluded"):
                (fixture / f"{name}.json").write_text(json.dumps(invalid), encoding="utf-8")
            output = root / "sensitivity"

            self.assertEqual(run_sensitivities(baseline, multihit, output, offline_fixture=fixture), 2)
            self.assertIn("UNAVAILABLE_INVALID_PROVIDER_RESPONSE", (output / "sensitivity_report.md").read_text())
            self.assertIn("UNAVAILABLE_INVALID_PROVIDER_RESPONSE", (output / "multihit" / "manifest.json").read_text())
    def test_identity_hashes_p1_rows_and_records_unmapped_multihit_exclusions(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            multihit = root / "multihit.csv"
            pd.DataFrame({"SYMBOL": ["GENE1", "UNMAPPED"]}).to_csv(multihit, index=False)
            manifests = []
            for index, raw_p in enumerate((0.1, 0.2)):
                baseline = root / f"baseline-{index}"
                baseline.mkdir()
                pd.DataFrame([
                    {"SYMBOL": "GENE1", "ENSEMBL_ID": "ENSG1", "STATUS": "mapped", "IN_FOREGROUND": "true"},
                    {"SYMBOL": "GENE2", "ENSEMBL_ID": "ENSG2", "STATUS": "mapped", "IN_FOREGROUND": "false"},
                    {"SYMBOL": "UNMAPPED", "ENSEMBL_ID": "", "STATUS": "unmapped", "IN_FOREGROUND": "true"},
                ]).to_csv(baseline / "p0_gene_qc.csv", index=False)
                pd.DataFrame([{"TERM_ID": "R-HSA-1", "TERM_NAME": "Term", "SOURCE": "REAC", "RAW_P": raw_p, "BH_Q": 0.2, "FOLD_ENRICHMENT": 1.0}]).to_csv(
                    baseline / "p1_enrichment_results.csv", index=False
                )
                output = root / f"out-{index}"
                self.assertEqual(run_sensitivities(baseline, multihit, output), 0)
                manifests.append(json.loads((output / "multihit" / "manifest.json").read_text()))

            self.assertNotEqual(manifests[0]["baseline_projection_sha256"], manifests[1]["baseline_projection_sha256"])
            self.assertEqual(manifests[0]["unmapped_multihit_symbols"], ["UNMAPPED"])
    def test_tampered_completed_output_fails_integrity_check(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline = root / "baseline"
            baseline.mkdir()
            pd.DataFrame([
                {"SYMBOL": "GENE1", "ENSEMBL_ID": "ENSG1", "STATUS": "mapped", "IN_FOREGROUND": "true"},
                {"SYMBOL": "GENE2", "ENSEMBL_ID": "ENSG2", "STATUS": "mapped", "IN_FOREGROUND": "false"},
            ]).to_csv(baseline / "p0_gene_qc.csv", index=False)
            pd.DataFrame(columns=["TERM_ID", "TERM_NAME", "SOURCE", "RAW_P", "BH_Q", "FOLD_ENRICHMENT"]).to_csv(
                baseline / "p1_enrichment_results.csv", index=False
            )
            multihit = root / "multihit.csv"
            pd.DataFrame({"SYMBOL": ["GENE1"]}).to_csv(multihit, index=False)
            fixture = root / "fixture"
            fixture.mkdir()
            for name in ("multihit", "hp_excluded"):
                (fixture / f"{name}.json").write_text(json.dumps({"meta": {}, "result": []}), encoding="utf-8")
            output = root / "sensitivity"
            self.assertEqual(run_sensitivities(baseline, multihit, output, offline_fixture=fixture), 0)
            self.assertTrue(verify_sensitivity_artifacts(output))
            (output / "multihit" / "enrichment.csv").write_text("tampered", encoding="utf-8")
            self.assertFalse(verify_sensitivity_artifacts(output))
if __name__ == "__main__":
    unittest.main()
