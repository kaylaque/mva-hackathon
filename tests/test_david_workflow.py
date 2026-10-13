"""Synthetic privacy-safe tests for the resumable native DAVID addendum."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd

REPO_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from mva_hackathon.david_workflow import _comparison, _stable_id, david_run_summary, run_david_addition


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class DavidWorkflowTests(unittest.TestCase):
    def _baseline(self, root: Path, *, foreground_count: int = 1) -> Path:
        baseline = root / "baseline"
        baseline.mkdir()
        genes = [
            {
                "SYMBOL": "GENE_A",
                "ENSEMBL_ID": "ENSG_A",
                "STATUS": "mapped",
                "IN_FOREGROUND": "True" if foreground_count else "False",
                "PRIVATE_VARIANT": "must-not-propagate",
            },
            {
                "SYMBOL": "GENE_B",
                "ENSEMBL_ID": "ENSG_B",
                "STATUS": "mapped",
                "IN_FOREGROUND": "False",
                "PRIVATE_VARIANT": "must-not-propagate",
            },
            {
                "SYMBOL": "UNMAPPED",
                "ENSEMBL_ID": None,
                "STATUS": "unmapped",
                "IN_FOREGROUND": "False",
                "PRIVATE_VARIANT": "must-not-propagate",
            },
        ]
        if foreground_count > 1:
            genes = [
                {
                    "SYMBOL": f"GENE_{index}",
                    "ENSEMBL_ID": f"ENSG_{index}",
                    "STATUS": "mapped",
                    "IN_FOREGROUND": "True",
                    "PRIVATE_VARIANT": "must-not-propagate",
                }
                for index in range(foreground_count)
            ] + genes[1:]
        pd.DataFrame(genes).to_csv(baseline / "p0_gene_qc.csv", index=False)
        pd.DataFrame([
            {
                "SOURCE": "GO:BP",
                "TERM_ID": "GO:0000001",
                "TERM_NAME": "Synthetic term",
                "BH_Q": 0.04,
                "PROVIDER_ADJUSTED_P": 0.04,
                "PROVIDER_CORRECTION": "BH",
                "RAW_P": 0.01,
                "FOLD_ENRICHMENT": 2.0,
                "CONTRIBUTING_GENES": "GENE_A",
                "BACKGROUND_LABEL": "custom",
                "PRIVATE_PATIENT": "must-not-propagate",
            },
        ]).to_csv(baseline / "p1_enrichment_results.csv", index=False)
        return baseline

    def _prepared(self, root: Path, *, foreground_count: int = 1) -> tuple[Path, Path, dict[str, object]]:
        baseline = self._baseline(root, foreground_count=foreground_count)
        output = root / "david-output"
        self.assertEqual(run_david_addition(baseline, output, None), 2)
        return baseline, output, json.loads((output / "david_manifest.json").read_text())

    def _exports(self, root: Path, manifest: dict[str, object], *, receipt_patch: dict[str, object] | None = None, chart: str | None = None, clusters: str | None = None) -> Path:
        exports = root / "exports"
        exports.mkdir(parents=True)
        chart_path = exports / "chart.csv"
        clusters_path = exports / "clusters.csv"
        chart_path.write_text(chart or "Category,Term,Count,P-Value,Benjamini,FDR,Fisher Exact,User Ids\nGOTERM_BP_DIRECT,GO:0000001~Synthetic,2,0.01,0.04,0.03,0.02,ENSG_A\n", encoding="utf-8")
        clusters_path.write_text(clusters or "Cluster,Cluster Enrichment Score,Category,Term,Count,P-Value,Benjamini\n1,2.1,GOTERM_BP_DIRECT,GO:0000001~Synthetic,2,0.01,0.04\n", encoding="utf-8")
        receipt: dict[str, object] = {
            "schema_version": 1,
            "run_key": manifest["run_key"],
            "list_sha256": manifest["list_sha256"],
            "namespace": "ENSEMBL_GENE_ID",
            "taxon": 9606,
            "population_confirmed": True,
            "timestamps": {"exported_at": "2026-10-13T00:00:00Z"},
            "selected_population": "custom uploaded background",
            "foreground_mapping": {"submitted": 1, "mapped": 1, "unique": 1},
            "background_mapping": {"mapped": 2, "reason": "DAVID mapping summary"},
            "interface": {
                "selected_categories": ["GOTERM_BP_DIRECT"],
                "settings": {
                    "chart": {"Count": 2, "EASE": 0.1},
                    "clustering": {
                        "displayed": "Medium",
                        "actual_categories_and_thresholds": {"overlap": 3, "similarity": 50},
                    },
                },
            },
            "reports": {
                "chart": {"complete": True, "complete_unfiltered": True, "empty": False},
                "clusters": {"complete": True, "complete_unfiltered": True, "empty": False},
            },
            "export_sha256": {"chart.csv": _sha256(chart_path), "clusters.csv": _sha256(clusters_path)},
        }
        if receipt_patch:
            receipt.update(receipt_patch)
        (exports / "receipt.json").write_text(json.dumps(receipt), encoding="utf-8")
        return exports

    def test_preparation_uses_boolean_string_mapping_and_excludes_private_columns(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline = self._baseline(root)
            before = {
                name: (baseline / name).read_bytes()
                for name in ("p0_gene_qc.csv", "p1_enrichment_results.csv")
            }
            output = root / "david-output"
            self.assertEqual(run_david_addition(baseline, output, None), 2)
            manifest = json.loads((output / "david_manifest.json").read_text())
            self.assertEqual((output / "foreground.txt").read_text(), "ENSG_A\n")
            self.assertEqual((output / "background.txt").read_text(), "ENSG_A\nENSG_B\n")
            self.assertEqual(manifest["status"], "awaiting_web_results")
            self.assertNotIn("PRIVATE_VARIANT", json.dumps(manifest))
            self.assertNotIn("PRIVATE_PATIENT", json.dumps(manifest))
            self.assertEqual(
                {name: (baseline / name).read_bytes() for name in before},
                before,
            )

    def test_foreground_over_3000_is_rejected_without_truncating_background(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline = self._baseline(root, foreground_count=3001)
            with self.assertRaisesRegex(ValueError, "3000"):
                run_david_addition(baseline, root / "out", None)

    def test_background_size_has_no_foreground_limit(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline = self._baseline(root)
            p0 = pd.read_csv(baseline / "p0_gene_qc.csv")
            extra = pd.DataFrame([
                {"SYMBOL": f"BG_{index}", "ENSEMBL_ID": f"ENSG_BG_{index}", "STATUS": "mapped", "IN_FOREGROUND": "False", "PRIVATE_VARIANT": "private"}
                for index in range(3001)
            ])
            pd.concat([p0, extra], ignore_index=True).to_csv(baseline / "p0_gene_qc.csv", index=False)
            output = root / "out"
            self.assertEqual(run_david_addition(baseline, output, None), 2)
            self.assertEqual(len((output / "background.txt").read_text().splitlines()), 3003)

    def test_receipt_requires_complete_mapping_interface_and_unfiltered_reports(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline, output, manifest = self._prepared(root)
            exports = self._exports(root, manifest, receipt_patch={"population_confirmed": False})
            with self.assertRaisesRegex(ValueError, "population_confirmed"):
                run_david_addition(baseline, output, exports)
            receipt = json.loads((exports / "receipt.json").read_text())
            receipt["population_confirmed"] = True
            receipt["reports"]["chart"]["complete_unfiltered"] = False
            (exports / "receipt.json").write_text(json.dumps(receipt), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "complete_unfiltered"):
                run_david_addition(baseline, output, exports)
            receipt["reports"]["chart"]["complete_unfiltered"] = True
            receipt["reports"]["clusters"]["complete"] = False
            (exports / "receipt.json").write_text(json.dumps(receipt), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "incomplete"):
                run_david_addition(baseline, output, exports)
            receipt["reports"]["clusters"]["complete"] = True
            del receipt["interface"]["selected_categories"]
            (exports / "receipt.json").write_text(json.dumps(receipt), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "selected_categories"):
                run_david_addition(baseline, output, exports)
            receipt["interface"]["selected_categories"] = ["GOTERM_BP_DIRECT"]
            del receipt["interface"]["settings"]
            (exports / "receipt.json").write_text(json.dumps(receipt), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "settings"):
                run_david_addition(baseline, output, exports)

    def test_receipt_requires_valid_timestamp_exact_chart_settings_and_recorded_categories(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline, output, manifest = self._prepared(root)
            exports = self._exports(root, manifest)
            receipt_path = exports / "receipt.json"
            receipt = json.loads(receipt_path.read_text())
            del receipt["timestamps"]
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "timestamps"):
                run_david_addition(baseline, output, exports)
            receipt["timestamps"] = {"exported_at": "not-a-timestamp"}
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "timestamp"):
                run_david_addition(baseline, output, exports)
            receipt["timestamps"] = {"exported_at": "2026-10-13T00:00:00Z"}
            receipt["interface"]["settings"]["chart"]["Count"] = 1
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Count"):
                run_david_addition(baseline, output, exports)
            receipt["interface"]["settings"]["chart"]["Count"] = 2
            receipt["interface"]["settings"]["chart"]["EASE"] = 0.2
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "EASE"):
                run_david_addition(baseline, output, exports)
            receipt["interface"]["settings"]["chart"]["EASE"] = 0.1
            receipt["interface"]["settings"]["clustering"]["actual_categories_and_thresholds"] = {}
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "clustering"):
                run_david_addition(baseline, output, exports)
            receipt["interface"]["settings"]["clustering"]["actual_categories_and_thresholds"] = {"overlap": 3}
            receipt["interface"]["selected_categories"] = ["OTHER_CATEGORY"]
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Category"):
                run_david_addition(baseline, output, exports)

    def test_receipt_rejects_each_required_root_field_and_invalid_mapping_counts(self) -> None:
        required = (
            "schema_version", "run_key", "list_sha256", "namespace", "taxon",
            "population_confirmed", "timestamps", "selected_population", "foreground_mapping",
            "background_mapping", "interface", "reports", "export_sha256",
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline, output, manifest = self._prepared(root)
            for field in required:
                with self.subTest(field=field):
                    exports = self._exports(root / field, manifest)
                    receipt_path = exports / "receipt.json"
                    receipt = json.loads(receipt_path.read_text())
                    del receipt[field]
                    receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
                    with self.assertRaises(ValueError):
                        run_david_addition(baseline, output, exports)
            exports = self._exports(root / "mapping", manifest)
            receipt_path = exports / "receipt.json"
            receipt = json.loads(receipt_path.read_text())
            receipt["foreground_mapping"]["submitted"] = 0
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "submitted count"):
                run_david_addition(baseline, output, exports)
            receipt["foreground_mapping"]["submitted"] = 1
            receipt["foreground_mapping"]["mapped"] = 2
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "exceeds"):
                run_david_addition(baseline, output, exports)
            receipt["foreground_mapping"]["mapped"] = 1
            receipt["background_mapping"] = {"mapped": None, "reason": ""}
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "background mapping reason"):
                run_david_addition(baseline, output, exports)
            hash_exports = self._exports(root / "hash", manifest)
            receipt_path = hash_exports / "receipt.json"
            receipt = json.loads(receipt_path.read_text())
            receipt["export_sha256"]["chart.csv"] = "not-a-sha256"
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "SHA-256"):
                run_david_addition(baseline, output, hash_exports)


    def test_receipt_rejects_unique_foreground_count_above_mapped_count(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline, output, manifest = self._prepared(root, foreground_count=2)
            exports = self._exports(root, manifest)
            receipt_path = exports / "receipt.json"
            receipt = json.loads(receipt_path.read_text())
            receipt["foreground_mapping"]["submitted"] = 2
            receipt["foreground_mapping"]["unique"] = 2
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unique count exceeds mapped"):
                run_david_addition(baseline, output, exports)

    def test_header_only_reports_require_explicit_completed_empty_receipt(self) -> None:

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline, output, manifest = self._prepared(root)
            exports = self._exports(
                root,
                manifest,
                chart="Category,Term,Count,P-Value,Benjamini\n",
                clusters="Cluster,Cluster Enrichment Score,Category,Term,Count,P-Value,Benjamini\n",
            )
            with self.assertRaisesRegex(ValueError, "header-only"):
                run_david_addition(baseline, output, exports)
            receipt = json.loads((exports / "receipt.json").read_text())
            receipt["reports"]["chart"]["empty"] = True
            receipt["reports"]["clusters"]["empty"] = True
            chart = exports / "chart.csv"
            clusters = exports / "clusters.csv"
            receipt["export_sha256"] = {"chart.csv": _sha256(chart), "clusters.csv": _sha256(clusters)}
            (exports / "receipt.json").write_text(json.dumps(receipt), encoding="utf-8")
            self.assertEqual(run_david_addition(baseline, output, exports), 0)

    def test_import_creates_immutable_version_and_identical_replay_does_not_write(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline, output, manifest = self._prepared(root)
            exports = self._exports(root, manifest)
            self.assertEqual(run_david_addition(baseline, output, exports), 0)
            summary = david_run_summary(output, 0)
            version = output / "imports" / str(summary["run_key"]) / str(summary["import_version"])
            self.assertTrue(version.is_dir())
            self.assertEqual((version / "receipt.json").read_bytes(), (exports / "receipt.json").read_bytes())
            before = {path.name: path.read_bytes() for path in version.iterdir() if path.is_file()}
            self.assertEqual(run_david_addition(baseline, output, exports), 0)
            after = {path.name: path.read_bytes() for path in version.iterdir() if path.is_file()}
            self.assertEqual(after, before)
            self.assertEqual(set(summary), {"status", "run_key", "import_version", "external_output_dir"})

    def test_changed_prepared_lists_or_baseline_projection_fail_without_force(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline, output, manifest = self._prepared(root)
            (output / "foreground.txt").write_text("ENSG_TAMPERED\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "STALE_OR_INVALID_PROVENANCE"):
                run_david_addition(baseline, output, None)
            self.assertEqual(run_david_addition(baseline, output, None, force=True), 2)
            self.assertEqual((output / "foreground.txt").read_text(), "ENSG_A\n")

            p1_path = baseline / "p1_enrichment_results.csv"
            added = pd.read_csv(p1_path).iloc[[0]].assign(TERM_ID="GO:0000002", TERM_NAME="Added term")
            pd.concat([pd.read_csv(p1_path), added], ignore_index=True).to_csv(p1_path, index=False)

            with self.assertRaisesRegex(ValueError, "STALE_OR_INVALID_PROVENANCE"):
                run_david_addition(baseline, output, None)

    def test_completed_import_resume_verifies_and_detects_tampered_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline, output, manifest = self._prepared(root)
            exports = self._exports(root, manifest)
            self.assertEqual(run_david_addition(baseline, output, exports), 0)

            completed = json.loads((output / "david_manifest.json").read_text())
            self.assertEqual(completed["status"], "completed")
            self.assertIsNotNone(completed["result_sha256"])
            self.assertEqual(run_david_addition(baseline, output, None), 0)

            version_dir = output / "imports" / completed["run_key"] / completed["import_version"]
            (version_dir / "david_chart.csv").write_text("Category,Term\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "STALE_OR_INVALID_PROVENANCE"):
                run_david_addition(baseline, output, None)

    def test_comparison_keeps_each_unrecognized_provider_term_unmapped(self) -> None:
        p1 = pd.DataFrame([{"SOURCE": "GO:BP", "TERM_ID": "not-a-stable-id", "CONTRIBUTING_GENES": "GENE_A"}])
        chart = pd.DataFrame([{"Category": "GOTERM_BP_DIRECT", "Term": "also-not-a-stable-id", "PROVIDER": "DAVID"}])
        comparison = _comparison(p1, chart, {"GENE_A": "ENSG_A"}, ["ENSG_A"])
        self.assertEqual(list(comparison["MATCH_STATUS"]), ["unmapped_term", "unmapped_term"])
        self.assertEqual(comparison.iloc[0]["GPROFILER_ROWS_JSON"] != "[]", True)
        self.assertEqual(comparison.iloc[1]["DAVID_ROWS_JSON"] != "[]", True)
        self.assertEqual(_stable_id("hsa04110~Cell cycle"), "hsa04110")

    def test_comparison_matches_native_prefixed_ids_and_full_terms(self) -> None:
        p1 = pd.DataFrame([
            {"SOURCE": "REAC", "TERM_ID": "REAC:R-HSA-141444", "CONTRIBUTING_GENES": "GENE_A", "BH_Q": 0.01},
            {"SOURCE": "WP", "TERM_ID": "WP:WP5344", "CONTRIBUTING_GENES": "GENE_A", "BH_Q": 0.02},
            {"SOURCE": "GO:BP", "TERM_ID": "GO:0000012", "CONTRIBUTING_GENES": "GENE_B", "BH_Q": 0.03},
        ])
        chart = pd.DataFrame([
            {"Category": "REACTOME_PATHWAY", "Term": "Mitotic Spindle Checkpoint", "Full Term": "R-HSA-141444~Mitotic Spindle Checkpoint", "User Ids": "ENSG_A", "PROVIDER": "DAVID"},
            {"Category": "KEGG_PATHWAY", "Term": "Cell cycle", "Full Term": "hsa04110:Cell cycle", "User Ids": "ENSG_A", "PROVIDER": "DAVID"},
            {"Category": "GOTERM_BP_DIRECT", "Term": "protein-containing complex assembly", "Full Term": "GO:0065003~protein-containing complex assembly", "User Ids": "ENSG_A", "PROVIDER": "DAVID"},
            {"Category": "SMART", "Term": "DSRM", "Full Term": "SM00358:DSRM", "User Ids": "ENSG_A", "PROVIDER": "DAVID"},
        ])

        comparison = _comparison(p1, chart, {"GENE_A": "ENSG_A", "GENE_B": "ENSG_B"}, ["ENSG_A", "ENSG_B"])
        status = {
            row["TERM_ID"]: (row["MATCH_STATUS"], row["GENE_JACCARD"])
            for _, row in comparison.iterrows()
            if row["TERM_ID"] is not None
        }

        self.assertEqual(status["R-HSA-141444"], ("shared", 1.0))
        self.assertEqual(status["WP5344"][0], "gprofiler_only")
        self.assertEqual(status["GO:0000012"][0], "gprofiler_only")
        self.assertEqual(status["GO:0065003"][0], "david_only")
        self.assertEqual(status["hsa04110"][0], "david_only")
        self.assertEqual(int(comparison["MATCH_STATUS"].eq("unmapped_term").sum()), 1)
        self.assertEqual(_stable_id("REAC:R-HSA-141444"), "R-HSA-141444")
        self.assertEqual(_stable_id("WP:WP5344"), "WP5344")
        self.assertEqual(_stable_id("KEGG:04110"), "hsa04110")
        self.assertIsNone(_stable_id("SM00358:DSRM"))
        self.assertIsNone(_stable_id("KW-1274~Primary mitochondrial disease"))
        self.assertIsNone(_stable_id("Primary mitochondrial disease"))


    def test_import_preserves_native_statistic_meanings(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline, output, manifest = self._prepared(root)
            self.assertEqual(run_david_addition(baseline, output, self._exports(root, manifest)), 0)
            summary = david_run_summary(output, 0)
            version = output / "imports" / str(summary["run_key"]) / str(summary["import_version"])
            chart = pd.read_csv(version / "david_chart.csv")
            clusters = pd.read_csv(version / "david_clusters.csv")
            self.assertEqual(chart.loc[0, "P-Value"], 0.01)
            self.assertEqual(chart.loc[0, "Benjamini"], 0.04)
            self.assertEqual(chart.loc[0, "FDR"], 0.03)
            self.assertEqual(clusters.loc[0, "Cluster Enrichment Score"], 2.1)


if __name__ == "__main__":
    unittest.main()
