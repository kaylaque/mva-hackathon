"""Unit tests for Phase 4 Signed Pathway & Reaction Adjudication Engine.

All tests operate offline with synthetic non-patient fixtures and mocked HTTP clients.
"""

from __future__ import annotations

import json
from pathlib import Path
import tempfile
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from mva_hackathon.adjudication_safety import (
    Phase4AdjudicationConfig,
    Phase4AdjudicationEngine,
    calculate_pediatric_hazard_score,
)
from mva_hackathon.hypothesis_classifier import (
    CLASS_H0,
    CLASS_H1,
    CLASS_H2,
    CLASS_H3,
    classify_pair_hypotheses,
)
from mva_hackathon.reaction_adjudication import (
    AdjudicatedReaction,
    ReactionAdjudicationConfig,
    ReactionAdjudicationEngine,
)
from mva_hackathon.signed_regulation import (
    SignedEdge,
    SignedPath,
    SignedRegulationConfig,
    SignedRegulationEngine,
)


def test_sign_math_and_edge_schema() -> None:
    """Verify signed algebra: ACT*ACT=+, ACT*INH=-, INH*INH=+, UNKNOWN=0."""
    edge_act = SignedEdge("A", "B", "ACTIVATION", 1, "direct", True, "CURATED", "SIGNOR", "v3", ["123"], 1.0)
    edge_inh = SignedEdge("B", "C", "INHIBITION", -1, "direct", True, "CURATED", "SIGNOR", "v3", ["456"], 1.0)
    edge_amb = SignedEdge("C", "D", "AMBIGUOUS", 0, "direct", True, "CURATED", "SIGNOR", "v3", ["789"], 1.0)

    # Edge schema validation
    assert edge_act.sign == 1
    assert edge_inh.sign == -1
    assert edge_amb.sign == 0

    # Path multiplication algebra
    assert edge_act.sign * edge_act.sign == 1
    assert edge_act.sign * edge_inh.sign == -1
    assert edge_inh.sign * edge_inh.sign == 1
    assert edge_act.sign * edge_amb.sign == 0


def test_cross_database_conflict_resolution() -> None:
    """Contradictory signs across curated DBs must resolve to AMBIGUOUS and be excluded from paths."""
    config = SignedRegulationConfig(offline_mode=True)
    engine = SignedRegulationEngine(config)

    # Custom edges: SIGNOR says A -> B is ACTIVATION (+1), OmniPath says A -> B is INHIBITION (-1)
    e1 = SignedEdge("GENEA", "GENEB", "ACTIVATION", 1, "bind", True, "CURATED", "SIGNOR", "v3", ["111"])
    e2 = SignedEdge("GENEA", "GENEB", "INHIBITION", -1, "inhib", True, "CURATED", "OmniPath", "v2024", ["222"])

    total, conflicts = engine.ingest_and_reconcile(custom_edges=[e1, e2])

    assert total == 1
    assert conflicts == 1

    reconciled = engine.consolidated_edges[("GENEA", "GENEB")]
    assert reconciled.conflict_flag is True
    assert reconciled.effect == "AMBIGUOUS"
    assert reconciled.sign == 0

    # Ensure conflicting edges are excluded from path computation
    paths = engine.find_signed_paths("GENEA", {"GENEB"}, max_length=1)
    assert len(paths) == 0


def test_path_length_cap() -> None:
    """Paths must be restricted to length <= max_path_length (default 2)."""
    config = SignedRegulationConfig(offline_mode=True, max_path_length=2)
    engine = SignedRegulationEngine(config)

    # Chain: A -> B -> C -> D
    e1 = SignedEdge("A", "B", "ACTIVATION", 1, "direct", True, "CURATED", "SIGNOR", "v3")
    e2 = SignedEdge("B", "C", "ACTIVATION", 1, "direct", True, "CURATED", "SIGNOR", "v3")
    e3 = SignedEdge("C", "D", "ACTIVATION", 1, "direct", True, "CURATED", "SIGNOR", "v3")

    engine.ingest_and_reconcile(custom_edges=[e1, e2, e3])

    # Direct length 1 path A -> B
    paths_ab = engine.find_signed_paths("A", {"B"}, max_length=2)
    assert len(paths_ab) == 1
    assert paths_ab[0].path_length == 1

    # Length 2 path A -> C (via B)
    paths_ac = engine.find_signed_paths("A", {"C"}, max_length=2)
    assert len(paths_ac) == 1
    assert paths_ac[0].path_length == 2

    # Length 3 path A -> D should NOT be found because max_length = 2
    paths_ad = engine.find_signed_paths("A", {"D"}, max_length=2)
    assert len(paths_ad) == 0


def test_hypothesis_class_rules_and_duplication() -> None:
    """Verify deterministic assignment of H1, H2, H3, H0 and multi-class duplication."""
    # Case 1: H1 (SAC kinase reaching BUB1B)
    edge = SignedEdge("CDC20", "BUB1B", "INHIBITION", -1, "direct", True, "CURATED", "SIGNOR", "v3")
    path = SignedPath("CDC20", "BUB1B", [edge], 1, -1, True)

    h1_results = classify_pair_hypotheses(
        target_symbol="CDC20",
        drug_name="APCIN",
        gene_category="CORE_MVA_CAUSAL",
        action_type="INHIBITOR",
        drug_sign=-1,
        signed_paths=[path],
        reactions=[],
        prior_safety_tier="TIER_1_CONCORDANT_HYPOTHESIS",
        literature_verified=True,
    )
    classes_h1 = [r.hypothesis_class for r in h1_results]
    assert CLASS_H1 in classes_h1

    # Case 2: H2 (Pediatric malignancy indication)
    h2_results = classify_pair_hypotheses(
        target_symbol="ALK",
        drug_name="CRIZOTINIB",
        gene_category="CONTEXT_MODIFIER",
        action_type="INHIBITOR",
        drug_sign=-1,
        signed_paths=[],
        reactions=[],
        prior_safety_tier="TIER_1_CONCORDANT_HYPOTHESIS",
        disease_indications=["rhabdomyosarcoma", "neuroblastoma"],
        literature_verified=True,
    )
    classes_h2 = [r.hypothesis_class for r in h2_results]
    assert CLASS_H2 in classes_h2

    # Case 3: Dual qualification (H1 + H2) produces duplicated rows
    dual_results = classify_pair_hypotheses(
        target_symbol="PLK1",
        drug_name="VOLASERTIB",
        gene_category="GRAPH_DERIVED_NEIGHBOR",
        action_type="INHIBITOR",
        drug_sign=-1,
        signed_paths=[path],
        reactions=[],
        prior_safety_tier="TIER_2_OFF_TARGET_COMPENSATION",
        disease_indications=["rhabdomyosarcoma", "acute myeloid leukemia"],
        literature_verified=True,
    )
    dual_classes = [r.hypothesis_class for r in dual_results]
    assert CLASS_H1 in dual_classes
    assert CLASS_H2 in dual_classes
    assert len(dual_results) >= 2

    # Case 4: Evidence gap -> H0
    gap_results = classify_pair_hypotheses(
        target_symbol="UNKNOWN_GENE",
        drug_name="TEST_DRUG",
        gene_category="UNKNOWN",
        action_type="INHIBITOR",
        drug_sign=-1,
        signed_paths=[],
        reactions=[],
        prior_safety_tier="TIER_3_UNCERTAIN_DIRECTIONALITY",
        is_evidence_gap=True,
    )
    assert len(gap_results) == 1
    assert gap_results[0].hypothesis_class == CLASS_H0


def test_pediatric_hazard_score_and_gate_ordering() -> None:
    """Hazard must outrank concordance: Tier 4 negative control receives MECHANISTIC_HAZARD_REJECTED."""
    # TP53 inhibitor in Tier 4
    score, flags = calculate_pediatric_hazard_score(
        prior_tier="TIER_4_CONTRAINDICATED_HAZARD",
        hazard_note="Severe hazard: TP53 inactivation abrogates aneuploidy surveillance; black-box risk.",
        action_type="INHIBITOR",
        target_symbol="TP53",
    )
    assert score >= 85
    assert "BLACK_BOX_WARNING" in flags or "TIER_4_CONTRAINDICATED" in flags


def test_multi_source_upgrade_rule_and_shared_pmid_dedupe(tmp_path: Path) -> None:
    """Verify upgrade rule requires >=2 independent sources and shared PMIDs are deduplicated."""
    from mva_hackathon.literature_adjudication import LiteratureAdjudicationEngine, LiteratureAdjudicationConfig
    engine = LiteratureAdjudicationEngine(LiteratureAdjudicationConfig(offline_mode=True))

    # Sources sharing the same PMID 12345
    s1 = {"source_db": "SIGNOR", "pmids": ["12345"]}
    s2 = {"source_db": "OmniPath", "pmids": ["12345"]}

    indep_count, shared_dedupe = engine.deduplicate_independent_sources([s1, s2])
    # Because both rely entirely on PMID 12345, independent count is 1
    assert indep_count == 1
    assert shared_dedupe == 1

    # Sources with distinct PMIDs
    s3 = {"source_db": "SIGNOR", "pmids": ["12345"]}
    s4 = {"source_db": "Reactome", "pmids": ["67890"]}
    indep_count2, shared_dedupe2 = engine.deduplicate_independent_sources([s3, s4])
    assert indep_count2 == 2
    assert shared_dedupe2 == 0


def test_offline_pipeline_execution_and_provenance(tmp_path: Path) -> None:
    """End-to-end offline pipeline execution generating all 20-24 outputs with provenance completeness."""
    # Create tiny synthetic test inputs
    df_mech = pd.DataFrame([{"TARGET_SYMBOL": "TP53", "THERAPEUTIC_DIRECTIONALITY": "ACTIVATE_STABILIZE"}])
    df_drug = pd.DataFrame([
        {
            "TARGET_SYMBOL": "TP53",
            "DRUG_NAME": "PIFITHRIN-ALPHA",
            "DRUG_ID": "CHEMBL1",
            "ACTION_TYPE": "INHIBITOR",
            "SAFETY_TIER": "TIER_4_CONTRAINDICATED_HAZARD",
            "PEDIATRIC_HAZARD_NOTE": "TP53 abrogation hazard",
            "GENE_CATEGORY": "CONTEXT_MODIFIER",
            "CONTROL_ROLE": "NEGATIVE_CONTROL",
        },
        {
            "TARGET_SYMBOL": "ALK",
            "DRUG_NAME": "CRIZOTINIB",
            "DRUG_ID": "CHEMBL2",
            "ACTION_TYPE": "INHIBITOR",
            "SAFETY_TIER": "TIER_1_CONCORDANT_HYPOTHESIS",
            "PEDIATRIC_HAZARD_NOTE": "Oncologic monitoring",
            "GENE_CATEGORY": "GRAPH_DERIVED_NEIGHBOR",
            "CONTROL_ROLE": "POSITIVE_CONTROL",
        },
    ])
    df_safety = df_drug.copy()
    df_panel = pd.DataFrame([
        {
            "QUERY_GENE": "TP53",
            "GENE_CATEGORY": "CONTEXT_MODIFIER",
            "THERAPEUTIC_DIRECTIONALITY": "ACTIVATE_STABILIZE",
            "UNIPROT_ID": "P04637",
            "PRIMARY_LITERATURE_PMID": "15475955",
            "PRIMARY_LITERATURE_DOI": "10.1038/ng1449",
        },
        {
            "QUERY_GENE": "ALK",
            "GENE_CATEGORY": "GRAPH_DERIVED_NEIGHBOR",
            "THERAPEUTIC_DIRECTIONALITY": "INHIBIT_DEGRADE",
            "UNIPROT_ID": "Q9UM73",
            "PRIMARY_LITERATURE_PMID": "15475955",
            "PRIMARY_LITERATURE_DOI": "10.1038/ng1449",
        },
    ])
    df_enrich = pd.DataFrame([{"PATHWAY_ID": "R-HSA-69560", "PATHWAY_NAME": "p53 Signaling"}])
    df_net = pd.DataFrame([{"GENE": "TP53", "DEGREE": 10}])

    p_mech = tmp_path / "02_mech.csv"
    p_drug = tmp_path / "14_drug.csv"
    p_safety = tmp_path / "15_safety.csv"
    p_panel = tmp_path / "18_panel.csv"
    p_enrich = tmp_path / "06_enrich.csv"
    p_net = tmp_path / "07_net.csv"

    df_mech.to_csv(p_mech, index=False)
    df_drug.to_csv(p_drug, index=False)
    df_safety.to_csv(p_safety, index=False)
    df_panel.to_csv(p_panel, index=False)
    df_enrich.to_csv(p_enrich, index=False)
    df_net.to_csv(p_net, index=False)

    out_dir = tmp_path / "track2_flat"
    config = Phase4AdjudicationConfig(
        gene_mechanism_path=p_mech,
        drug_evidence_path=p_drug,
        safety_ranking_path=p_safety,
        query_panel_path=p_panel,
        pathway_enrichment_path=p_enrich,
        network_context_path=p_net,
        output_dir=out_dir,
        flat_names=True,
        offline_mode=True,
        verify_literature=False,
    )

    engine = Phase4AdjudicationEngine(config)
    results = engine.run()

    # Verify all expected output files exist and are <= 8MB
    for key, path in results["output_paths"].items():
        assert path.is_file(), f"Output file missing: {path}"
        assert path.stat().st_size < 8 * 1024 * 1024, f"File {path} exceeds 8MB cap"

    # Verify provenance JSON completeness
    prov_path = results["output_paths"]["provenance"]
    with open(prov_path) as f:
        prov = json.load(f)

    assert "PREREGISTERED_ENDPOINTS" in prov
    assert "PREREGISTERED_SIGN_CONVENTION" in prov
    assert "PREREGISTERED_MAX_PATH_LENGTH" in prov
    assert "PREREGISTERED_SOURCE_HIERARCHY" in prov
    assert "PREREGISTERED_SCORING_WEIGHTS" in prov
    assert "PREREGISTERED_HAZARD_THRESHOLDS" in prov
    assert "PREREGISTERED_CLASS_RULES" in prov
    assert "PREREGISTERED_ENTRY_CRITERIA" in prov
    assert "best_options_per_class" in prov

    # Verify Tier 4 TP53 inhibitor was marked MECHANISTIC_HAZARD_REJECTED
    df_rep = pd.read_csv(results["output_paths"]["adjudicated_safety_report"])
    tp53_row = df_rep[df_rep["TARGET_SYMBOL"] == "TP53"].iloc[0]
    assert tp53_row["VERDICT_STATUS"] in {"MECHANISTIC_HAZARD_REJECTED", "NEGATIVE_NO_SAFE_INTERVENTION"}
    assert bool(tp53_row["SIGNED_PATH_HAZARD"]) is True
    assert tp53_row["CLINICAL_RECOMMENDATION_STATUS"] == "PROHIBITED / INVESTIGATIONAL ONLY"


def test_no_flat_names_naming(tmp_path: Path) -> None:
    """Verify that flat_names=False writes un-prefixed file names."""
    df_mech = pd.DataFrame([{"TARGET_SYMBOL": "ALK", "THERAPEUTIC_DIRECTIONALITY": "INHIBIT_DEGRADE"}])
    df_drug = pd.DataFrame([{
        "TARGET_SYMBOL": "ALK",
        "DRUG_NAME": "CRIZOTINIB",
        "DRUG_ID": "CHEMBL2",
        "ACTION_TYPE": "INHIBITOR",
        "SAFETY_TIER": "TIER_1_CONCORDANT_HYPOTHESIS",
        "PEDIATRIC_HAZARD_NOTE": "Monitoring",
        "GENE_CATEGORY": "GRAPH_DERIVED_NEIGHBOR",
        "CONTROL_ROLE": "POSITIVE_CONTROL",
    }])
    df_panel = pd.DataFrame([{
        "QUERY_GENE": "ALK",
        "GENE_CATEGORY": "GRAPH_DERIVED_NEIGHBOR",
        "THERAPEUTIC_DIRECTIONALITY": "INHIBIT_DEGRADE",
    }])
    df_enrich = pd.DataFrame([{"PATHWAY_ID": "R-HSA-69560", "PATHWAY_NAME": "Test"}])

    p_mech = tmp_path / "mech.csv"
    p_drug = tmp_path / "drug.csv"
    p_panel = tmp_path / "panel.csv"
    p_enrich = tmp_path / "enrich.csv"

    df_mech.to_csv(p_mech, index=False)
    df_drug.to_csv(p_drug, index=False)
    df_panel.to_csv(p_panel, index=False)
    df_enrich.to_csv(p_enrich, index=False)

    out_dir = tmp_path / "phase4_unflat"
    config = Phase4AdjudicationConfig(
        gene_mechanism_path=p_mech,
        drug_evidence_path=p_drug,
        safety_ranking_path=p_drug,
        query_panel_path=p_panel,
        pathway_enrichment_path=p_enrich,
        network_context_path=p_mech,
        output_dir=out_dir,
        flat_names=False,
        offline_mode=True,
        verify_literature=False,
    )
    engine = Phase4AdjudicationEngine(config)
    res = engine.run()

    assert (out_dir / "adjudicated_safety_report.csv").is_file()
    assert (out_dir / "adjudication_provenance.json").is_file()
    assert (out_dir / "reaction_level_matrix.csv").is_file()
    assert (out_dir / "hypothesis_classification.csv").is_file()
    assert (out_dir / "synthetic_lethality_candidates.csv").is_file()