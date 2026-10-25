"""Synthetic privacy-safe tests for the independent pathway workflow."""

from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import MagicMock

import pandas as pd

REPO_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from mva_hackathon.pathway_workflow import PathwayWorkflow, PathwayWorkflowConfig


class PathwayInputQcTests(unittest.TestCase):
    def test_p0_reads_only_allowlisted_gene_columns_and_marks_foreground(self) -> None:
        """P0 must not propagate variant-level fields into its normalized table."""
        source = pd.DataFrame(
            [
                {
                    "SYMBOL": "BUB1B",
                    "TOTAL_VARIANTS_EVALUATED": 2,
                    "TOTAL_PASSING_CANDIDATE_VARIANTS": 1,
                    "THERAPEUTIC_DIRECTIONALITY": "RESTORE_CHECKPOINT_FUNCTION",
                    "TOP_VARIANT_KEYS": "private-variant-key",
                },
                {
                    "SYMBOL": "GENE0",
                    "TOTAL_VARIANTS_EVALUATED": 1,
                    "TOTAL_PASSING_CANDIDATE_VARIANTS": 0,
                    "TOP_VARIANT_KEYS": "private-variant-key",
                },
            ]
        )

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "genes.csv"
            source.to_csv(path, index=False)
            workflow = PathwayWorkflow(PathwayWorkflowConfig(min_passing=1, offline=True))
            normalized, report = workflow.run_p0(path)

        self.assertEqual(list(normalized.columns), [
            "SYMBOL",
            "TOTAL_VARIANTS_EVALUATED",
            "TOTAL_PASSING_CANDIDATE_VARIANTS",
            "THERAPEUTIC_DIRECTIONALITY",
            "PASSING_CONSEQUENCE_SPECTRUM",
            "IN_FOREGROUND",
        ])
        self.assertNotIn("TOP_VARIANT_KEYS", normalized.columns)
        self.assertEqual(normalized.loc[0, "IN_FOREGROUND"], True)
        self.assertEqual(report["foreground_count"], 1)
        self.assertEqual(report["background_count"], 2)

    def test_p0_rejects_conflicting_duplicate_counts(self) -> None:
        """Conflicting duplicates must fail instead of silently summing candidates."""
        source = pd.DataFrame(
            [
                {"SYMBOL": "BUB1B", "TOTAL_VARIANTS_EVALUATED": 1, "TOTAL_PASSING_CANDIDATE_VARIANTS": 1},
                {"SYMBOL": "BUB1B", "TOTAL_VARIANTS_EVALUATED": 2, "TOTAL_PASSING_CANDIDATE_VARIANTS": 1},
            ]
        )

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "genes.csv"
            source.to_csv(path, index=False)
            workflow = PathwayWorkflow(PathwayWorkflowConfig(offline=True))
            with self.assertRaisesRegex(ValueError, "conflicting duplicate"):
                workflow.run_p0(path)


    def test_hgnc_mapping_only_uses_symbols_and_marks_unmapped(self) -> None:
        """HGNC requests must send normalized gene symbols, never patient-derived fields."""
        response = MagicMock()
        response.status_code = 200
        response.json.return_value = {
            "response": {"docs": [{"hgnc_id": "HGNC:1148", "ensembl_gene_id": "ENSG000001", "symbol": "BUB1B"}]}
        }
        unmapped_response = MagicMock()
        unmapped_response.status_code = 404
        client = MagicMock()
        client.get.side_effect = [response, unmapped_response]

        qc = PathwayWorkflow(PathwayWorkflowConfig(offline=True)).map_hgnc(["BUB1B", "UNKNOWN"], client)

        self.assertEqual(qc.loc[0, "ENSEMBL_ID"], "ENSG000001")
        self.assertEqual(qc.loc[1, "STATUS"], "unmapped")
        self.assertEqual(client.get.call_args_list[0].args[0], "https://rest.genenames.org/fetch/symbol/BUB1B")

    def test_gprofiler_payload_contains_only_resolved_gene_ids(self) -> None:
        """Provider payloads must exclude candidate counts and clinical/variant fields."""
        genes = pd.DataFrame(
            [
                {"ENSEMBL_ID": "ENSG1", "IN_FOREGROUND": True, "TOTAL_PASSING_CANDIDATE_VARIANTS": 2},
                {"ENSEMBL_ID": "ENSG2", "IN_FOREGROUND": False, "TOTAL_PASSING_CANDIDATE_VARIANTS": 0},
            ]
        )
        payload = PathwayWorkflow(PathwayWorkflowConfig()).gprofiler_payload(genes)

        self.assertEqual(payload["query"], ["ENSG1"])
        self.assertEqual(payload["background"], ["ENSG1", "ENSG2"])
        self.assertNotIn("TOTAL_PASSING_CANDIDATE_VARIANTS", str(payload))

    def test_request_cache_key_is_stable_for_reordered_payload_keys(self) -> None:
        """Cache identity is endpoint plus canonical payload, never secrets or source rows."""
        workflow = PathwayWorkflow(PathwayWorkflowConfig())
        first = workflow.request_cache_key("https://example.test/api", {"query": ["ENSG1"], "background": ["ENSG1", "ENSG2"]})
        second = workflow.request_cache_key("https://example.test/api", {"background": ["ENSG1", "ENSG2"], "query": ["ENSG1"]})
        self.assertEqual(first, second)


    def test_gprofiler_payload_can_exclude_hp_for_sensitivity(self) -> None:
        """HP-excluded sensitivity must alter provider source family, not post-filter rows."""
        genes = pd.DataFrame(
            [
                {"ENSEMBL_ID": "ENSG1", "IN_FOREGROUND": True},
                {"ENSEMBL_ID": "ENSG2", "IN_FOREGROUND": False},
            ]
        )

        payload = PathwayWorkflow(PathwayWorkflowConfig()).gprofiler_payload(
            genes, sources=["GO:BP", "REAC", "WP"]
        )

        self.assertEqual(payload["significance_threshold_method"], "fdr")
        self.assertNotIn("HP", payload["sources"])
class PathwayEnrichmentTests(unittest.TestCase):
    def test_p1_preserves_provider_bh_and_decodes_intersections(self) -> None:
        """Provider-adjusted p-values remain authoritative and identifiers map by query mapping."""
        normalized = pd.DataFrame(
            [
                {"SYMBOL": "BUB1B", "TOTAL_VARIANTS_EVALUATED": 2, "TOTAL_PASSING_CANDIDATE_VARIANTS": 1, "THERAPEUTIC_DIRECTIONALITY": "RESTORE", "PASSING_CONSEQUENCE_SPECTRUM": pd.NA, "IN_FOREGROUND": True},
                {"SYMBOL": "MAD1L1", "TOTAL_VARIANTS_EVALUATED": 1, "TOTAL_PASSING_CANDIDATE_VARIANTS": 1, "THERAPEUTIC_DIRECTIONALITY": "RESTORE", "PASSING_CONSEQUENCE_SPECTRUM": pd.NA, "IN_FOREGROUND": True},
                {"SYMBOL": "GENE0", "TOTAL_VARIANTS_EVALUATED": 1, "TOTAL_PASSING_CANDIDATE_VARIANTS": 0, "THERAPEUTIC_DIRECTIONALITY": pd.NA, "PASSING_CONSEQUENCE_SPECTRUM": pd.NA, "IN_FOREGROUND": False},
            ]
        )
        response = {
            "meta": {"query_metadata": {"queries": {"query_1": {"mapping": {"ENSG1": ["BUB1B"], "ENSG2": ["MAD1L1"]}}}}},
            "result": [{
                "source": "REAC",
                "native": "R-HSA-141444",
                "name": "Mitotic Spindle Checkpoint",
                "term_size": 10,
                "query_size": 2,
                "intersection_size": 2,
                "effective_domain_size": 100,
                "intersections": ["ENSG2", "ENSG1"],
                "p_value": 0.01,
            }],
        }

        results = PathwayWorkflow(PathwayWorkflowConfig(offline=True)).run_p1(normalized, response)

        row = results.iloc[0]
        self.assertEqual(row["BH_Q"], 0.01)
        self.assertEqual(row["PROVIDER_ADJUSTED_P"], 0.01)
        self.assertEqual(row["PROVIDER_CORRECTION"], "BH")
        self.assertEqual(row["CONTRIBUTING_GENES"], "BUB1B;MAD1L1")
        self.assertLess(row["RAW_P"], row["BH_Q"])

    def test_gprofiler_request_uses_gene_only_payload(self) -> None:
        """Live enrichment transport receives the exact privacy-filtered payload."""
        client = MagicMock()
        response = MagicMock()
        response.status_code = 200
        response.json.return_value = {"result": []}
        client.post.return_value = response
        genes = pd.DataFrame([
            {"ENSEMBL_ID": "ENSG1", "IN_FOREGROUND": True},
            {"ENSEMBL_ID": "ENSG2", "IN_FOREGROUND": False},
        ])

        PathwayWorkflow(PathwayWorkflowConfig()).fetch_gprofiler(genes, client)

        payload = client.post.call_args.kwargs["json_data"]
        self.assertEqual(payload["query"], ["ENSG1"])
        self.assertEqual(payload["background"], ["ENSG1", "ENSG2"])




class PathwayPrioritizationTests(unittest.TestCase):
    def test_p2_keeps_significance_order_and_clusters_redundant_terms(self) -> None:
        """P2 ranks by provider statistics alone and never adds an MVA relevance score."""
        enrichment = pd.DataFrame(
            [
                {"SOURCE": "REAC", "TERM_ID": "R1", "TERM_NAME": "Mitotic spindle checkpoint", "BH_Q": 0.01, "FOLD_ENRICHMENT": 5.0, "INTERSECTION_SIZE": 2, "CONTRIBUTING_GENES": "BUB1B;MAD1L1"},
                {"SOURCE": "GO:BP", "TERM_ID": "G1", "TERM_NAME": "cell cycle", "BH_Q": 0.02, "FOLD_ENRICHMENT": 4.0, "INTERSECTION_SIZE": 2, "CONTRIBUTING_GENES": "BUB1B;MAD1L1"},
            ]
        )

        results = PathwayWorkflow(PathwayWorkflowConfig(offline=True)).run_p2(enrichment)

        self.assertEqual(list(results["TERM_ID"]), ["R1", "G1"])
        self.assertNotIn("MVA_RELEVANT", results.columns)
        self.assertNotIn("CORE_MVA_GENE_OVERLAP", results.columns)
        self.assertEqual(results.loc[results["TERM_ID"] == "G1", "REDUNDANT_OF"].item(), "R1")


class PathwayTargetAndDrugTests(unittest.TestCase):
    def test_p4_retains_full_pathway_members_as_exploratory_targets(self) -> None:
        """A pathway member absent from Phase 1 remains an explicit exploratory target."""
        pathways = pd.DataFrame([{"TERM_ID": "R-HSA-141444", "BH_Q": 0.01}])
        memberships = {"R-HSA-141444": [{"SYMBOL": "BUB1B", "ENSEMBL_ID": "ENSG1"}, {"SYMBOL": "AURKB", "ENSEMBL_ID": "ENSG2"}]}
        phase1 = pd.DataFrame([{"SYMBOL": "BUB1B", "THERAPEUTIC_DIRECTIONALITY": "RESTORE_CHECKPOINT_FUNCTION"}])

        targets = PathwayWorkflow(PathwayWorkflowConfig(offline=True)).run_p4(pathways, memberships, phase1)

        aurkb = targets.loc[targets["SYMBOL"] == "AURKB"].iloc[0]
        self.assertEqual(aurkb["TARGET_SUPPORT"], "PATHWAY_EXPLORATORY")
        self.assertFalse(aurkb["VARIANT_SUPPORTED"])

    def test_p5_marks_downregulator_as_a_contradiction(self) -> None:
        """All inhibitor-like action labels use the pathway safety rule."""
        targets = pd.DataFrame([{"SYMBOL": "BUB1B", "ENSEMBL_ID": "ENSG1", "THERAPEUTIC_DIRECTIONALITY": "RESTORE_CHECKPOINT_FUNCTION"}])
        evidence = PathwayWorkflow(PathwayWorkflowConfig(offline=True)).run_p5(
            targets,
            [{"SOURCE": "ChEMBL", "TARGET_SYMBOL": "BUB1B", "ACTION_TYPE": "DOWNREGULATOR"}],
        )
        self.assertTrue(evidence.loc[0, "CONTRADICTION_FLAG"])

    def test_p5_marks_inhibitors_of_restoration_targets_as_hazards(self) -> None:
        """Inhibiting a target whose Phase 1 direction requires restoration is contradictory."""
        targets = pd.DataFrame([{"SYMBOL": "BUB1B", "ENSEMBL_ID": "ENSG1", "THERAPEUTIC_DIRECTIONALITY": "RESTORE_CHECKPOINT_FUNCTION"}])
        records = [{"SOURCE": "ChEMBL", "TARGET_SYMBOL": "BUB1B", "ACTION_TYPE": "INHIBITOR", "DRUG_NAME": "synthetic"}]

        evidence = PathwayWorkflow(PathwayWorkflowConfig(offline=True)).run_p5(targets, records)

        self.assertTrue(evidence.loc[0, "CONTRADICTION_FLAG"])
        self.assertEqual(evidence.loc[0, "DIRECTIONALITY_LABEL"], "CONTRADICTION")


class PathwayEvaluationTests(unittest.TestCase):
    def test_p6_compares_against_the_queried_string_panel(self) -> None:
        """P6 reports descriptive overlap against query-panel targets, not all network rows."""
        pathway_targets = pd.DataFrame([{"SYMBOL": "BUB1B"}, {"SYMBOL": "AURKB"}])
        queried_panel = pd.DataFrame([{"QUERY_GENE": "BUB1B"}, {"QUERY_GENE": "TRIP13"}])

        summary = PathwayWorkflow(PathwayWorkflowConfig(offline=True)).run_p6(pathway_targets, queried_panel)

        self.assertEqual(summary["shared_targets"], 1)
        self.assertEqual(summary["pathway_only_targets"], 1)
        self.assertEqual(summary["string_only_targets"], 1)
        self.assertEqual(summary["jaccard"], 1 / 3)


class PathwayMechanismAndReportTests(unittest.TestCase):
    def test_p3_marks_pathway_only_gene_direction_as_unresolved(self) -> None:
        """Mechanism rows distinguish Phase 1 support from pathway-only membership."""
        phase1 = pd.DataFrame([{"SYMBOL": "BUB1B", "THERAPEUTIC_DIRECTIONALITY": "RESTORE_CHECKPOINT_FUNCTION"}])
        evidence = [{"TERM_ID": "R-HSA-1", "SYMBOL": "AURKB", "ROLE": "participant", "EVIDENCE_TIER": "EXPERIMENTAL"}]

        result = PathwayWorkflow(PathwayWorkflowConfig(offline=True)).run_p3(phase1, evidence)

        self.assertEqual(result.loc[0, "THERAPEUTIC_DIRECTIONALITY"], "UNRESOLVED_DIRECTION")

    def test_p7_reports_partial_stages_without_fabricated_counts(self) -> None:
        """P7 renders explicit unavailable values and the background limitation."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "p7_pathway_report.md"
            PathwayWorkflow(PathwayWorkflowConfig(offline=True)).run_p7(
                path,
                {"p0": {"status": "completed", "foreground_count": 2}, "p1": {"status": "failed"}},
            )
            report = path.read_text()

        self.assertIn("p1 | failed | UNAVAILABLE", report)
        self.assertIn("Background is Phase 1 evaluated genes", report)

    def test_p7_includes_required_markdown_sections(self) -> None:
        """Report content remains complete even though its format is Markdown-only."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "p7_pathway_report.md"
            PathwayWorkflow(PathwayWorkflowConfig(offline=True)).run_p7(
                path,
                {
                    "p0": {"status": "completed", "foreground_count": 2, "background_count": 12},
                    "p6": {"status": "partial", "reason": "STRING_BASELINE_UNAVAILABLE"},
                },
            )
            report = path.read_text()

        for heading in (
            "## Run header",
            "## Top pathways",
            "## Top targets",
            "## Top drug hypotheses",
            "## STRING comparison",
            "## Next actions",
        ):
            self.assertIn(heading, report)



class PathwayEmptyResultTests(unittest.TestCase):
    def test_empty_p2_has_explicit_no_significant_enrichment_status(self) -> None:
        """No BH-significant terms is a scientific result, not a completed pathway panel."""
        summary = PathwayWorkflow(PathwayWorkflowConfig(offline=True)).stage_summary(
            pd.DataFrame(), "NO_SIGNIFICANT_ENRICHMENT"
        )

        self.assertEqual(
            summary,
            {
                "status": "completed",
                "outcome": "NO_SIGNIFICANT_ENRICHMENT",
                "count": 0,
            },
        )

    def test_p7_reports_completed_null_outcome(self) -> None:
        """Scientific nulls must not be rendered as unavailable execution failures."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "p7_pathway_report.md"
            PathwayWorkflow(PathwayWorkflowConfig(offline=True)).run_p7(
                path,
                {"p2": {"status": "completed", "outcome": "NO_SIGNIFICANT_ENRICHMENT", "count": 0}},
            )

            report = path.read_text()

        self.assertIn("NO_SIGNIFICANT_ENRICHMENT", report)

    def test_empty_downstream_evidence_has_schema(self) -> None:
        """Zero eligible pathways must produce readable header-only evidence tables."""
        workflow = PathwayWorkflow(PathwayWorkflowConfig(offline=True))

        self.assertIn("SYMBOL", workflow.run_p3(pd.DataFrame(columns=["SYMBOL"]), []).columns)
        targets = pd.DataFrame(columns=["SYMBOL", "THERAPEUTIC_DIRECTIONALITY"])
        self.assertIn("TARGET_SYMBOL", workflow.run_p5(targets, []).columns)

    def test_source_native_evidence_fields_survive_normalization(self) -> None:
        """Pathway evidence must retain source-native identifiers and measurements."""
        targets = pd.DataFrame(
            [{"SYMBOL": "BUB1B", "THERAPEUTIC_DIRECTIONALITY": "RESTORE_CHECKPOINT_FUNCTION"}]
        )
        evidence = PathwayWorkflow(PathwayWorkflowConfig(offline=True)).run_p5(
            targets,
            [{
                "SOURCE": "ChEMBL",
                "TARGET_SYMBOL": "BUB1B",
                "ACTION_TYPE": "INHIBITOR",
                "ACTIVITY_ID": "activity-1",
                "STANDARD_VALUE": 10,
            }],
        )

        self.assertEqual(evidence.loc[0, "ACTIVITY_ID"], "activity-1")
        self.assertEqual(evidence.loc[0, "STANDARD_VALUE"], 10)

    def test_p7_explains_empty_downstream_sections(self) -> None:
        """Null pathway enrichment must explain absent targets and drugs."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "p7_pathway_report.md"
            PathwayWorkflow(PathwayWorkflowConfig(offline=True)).run_p7(
                path,
                {"p2": {"status": "completed", "outcome": "NO_SIGNIFICANT_ENRICHMENT", "count": 0}},
                pd.DataFrame(),
                pd.DataFrame(),
                pd.DataFrame(),
            )
            report = path.read_text()

        self.assertIn("No terms passed provider BH-FDR", report)
        self.assertIn("No eligible pathway targets", report)
if __name__ == "__main__":
    unittest.main()
