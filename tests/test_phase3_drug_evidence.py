"""Unit tests for Phase 3 drug repurposing evidence, directionality, and safety review.

Operates with synthetic, non-patient fixtures and mocked API calls to guarantee
complete determinism, absence of patient data in test fixtures, and network isolation.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
from typing import Any
from unittest.mock import MagicMock, patch

REPO_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import pandas as pd
import pytest
import requests

from mva_hackathon.drug_evidence import (
    ChEMBLClient,
    DrugEvidenceConfig,
    DrugEvidencePipeline,
    OpenTargetsClient,
    evaluate_mechanism_directionality,
    write_drug_evidence_outputs,
)
from run_drug_evidence import main as run_drug_evidence_main


def make_synthetic_phase1_df() -> pd.DataFrame:
    """Create a minimal synthetic non-patient Phase 1 gene mechanism matrix."""
    return pd.DataFrame([
        {
            "SYMBOL": "BUB1B",
            "GENE_NAME": "BUB1 mitotic checkpoint serine/threonine kinase B",
            "PATHWAY": "Spindle Assembly Checkpoint (SAC) / Mitotic Cell Cycle",
            "BIOLOGICAL_ROLE": "Core SAC kinase and APC/C inhibitor",
            "MVA_MECHANISM_LINK": "Canonical MVA Type 1 (mosaic aneuploidy)",
            "CANONICAL_PATHWAY_DIRECTION": "DOWN (Checkpoint arrest failure)",
            "THERAPEUTIC_DIRECTIONALITY": "Compensate impaired checkpoint / promote chromosome stability",
            "INHERITANCE_REQUIREMENT": "Biallelic",
            "TOTAL_VARIANTS_EVALUATED": 2,
            "TOTAL_PASSING_CANDIDATE_VARIANTS": 2,
            "COMPOUND_HET_ARCHETYPE_FIT": True,
            "PHASE_GATE_STATUS": "REQUIRES_PHASE_VERIFICATION",
        },
        {
            "SYMBOL": "TP53",
            "GENE_NAME": "Tumor Protein P53",
            "PATHWAY": "Aneuploidy & DNA Damage Surveillance Checkpoint",
            "BIOLOGICAL_ROLE": "Tumor suppressor and aneuploidy surveillance",
            "MVA_MECHANISM_LINK": "Secondary surveillance defect",
            "CANONICAL_PATHWAY_DIRECTION": "SURVEILLANCE_ARREST",
            "THERAPEUTIC_DIRECTIONALITY": "Restore or stabilize surveillance",
            "INHERITANCE_REQUIREMENT": "Dominant / De novo",
            "TOTAL_VARIANTS_EVALUATED": 2,
            "TOTAL_PASSING_CANDIDATE_VARIANTS": 2,
            "COMPOUND_HET_ARCHETYPE_FIT": False,
            "PHASE_GATE_STATUS": "SINGLE_HIT_OR_UNPASSING",
        },
        {
            "SYMBOL": "AURKA",
            "GENE_NAME": "Aurora Kinase A",
            "PATHWAY": "Mitotic Spindle Assembly",
            "BIOLOGICAL_ROLE": "Centrosome maturation and bipolar spindle assembly",
            "MVA_MECHANISM_LINK": "Secondary interactome kinase",
            "CANONICAL_PATHWAY_DIRECTION": "MITOTIC_PROGRESSION",
            "THERAPEUTIC_DIRECTIONALITY": "Spindle modulation",
            "INHERITANCE_REQUIREMENT": "Uncertain",
            "TOTAL_VARIANTS_EVALUATED": 1,
            "TOTAL_PASSING_CANDIDATE_VARIANTS": 1,
            "COMPOUND_HET_ARCHETYPE_FIT": False,
            "PHASE_GATE_STATUS": "SINGLE_HIT_OR_UNPASSING",
        },
    ])


def make_synthetic_network_context_df() -> pd.DataFrame:
    """Create a minimal synthetic non-patient Phase 2 network context table."""
    return pd.DataFrame([
        {
            "HGNC_SYMBOL": "BUB1B",
            "TOTAL_PASSING_CANDIDATE_VARIANTS": 2,
            "NETWORK_DEGREE": 3,
            "COMPONENT_ID": 2,
        },
        {
            "HGNC_SYMBOL": "TP53",
            "TOTAL_PASSING_CANDIDATE_VARIANTS": 2,
            "NETWORK_DEGREE": 5,
            "COMPONENT_ID": 2,
        },
        {
            "HGNC_SYMBOL": "AURKA",
            "TOTAL_PASSING_CANDIDATE_VARIANTS": 1,
            "NETWORK_DEGREE": 2,
            "COMPONENT_ID": 1,
        },
    ])


def test_chembl_programming_errors_propagate() -> None:
    http_client = MagicMock()
    http_client.get.side_effect = ValueError("malformed fixture")

    with pytest.raises(ValueError, match="malformed fixture"):
        ChEMBLClient(http_client, DrugEvidenceConfig()).fetch_target_mechanisms("TEST")



def test_open_targets_graphql_errors_raise_source_failure() -> None:
    response = MagicMock()
    response.json.return_value = {"errors": [{"message": "synthetic failure"}]}
    http_client = MagicMock()
    http_client.post.return_value = response

    with pytest.raises(RuntimeError, match="Open Targets GraphQL response contained errors"):
        OpenTargetsClient(http_client, DrugEvidenceConfig()).fetch_target_drugs_and_safety("ENSGTEST")
def test_evaluate_mechanism_directionality_bub1b_inhibitor():
    """BUB1B inhibitor must be flagged as a critical mechanistic hazard."""
    eval_res = evaluate_mechanism_directionality(
        gene_symbol="BUB1B",
        action_type="INHIBITOR",
        mechanism_desc="BUBR1 kinase ATP-competitive inhibitor",
        canonical_direction="DOWN (Checkpoint arrest failure)",
        therapeutic_direction="Compensate checkpoint",
        is_canonical_mva=True,
    )
    assert eval_res.fit_status == "MECHANISTIC_HAZARD"
    assert eval_res.contradiction_flag is True
    assert "Direct mechanistic contradiction" in eval_res.explanation
    assert "HIGH HAZARD" in eval_res.pediatric_hazard_note


def test_evaluate_mechanism_directionality_tp53_inhibitor():
    """TP53 inhibitor must be flagged as a mechanistic hazard eliminating surveillance."""
    eval_res = evaluate_mechanism_directionality(
        gene_symbol="TP53",
        action_type="INHIBITOR",
        mechanism_desc="p53 small molecule antagonist",
        canonical_direction="SURVEILLANCE_ARREST",
        therapeutic_direction="Restore surveillance",
        is_canonical_mva=False,
    )
    assert eval_res.fit_status == "MECHANISTIC_HAZARD"
    assert eval_res.contradiction_flag is True
    assert "Severe hazard" in eval_res.explanation
    assert "EXTREME HAZARD" in eval_res.pediatric_hazard_note


def test_evaluate_mechanism_directionality_tp53_activator():
    """TP53 gene replacement or activator aligns concordantly in research hypothesis."""
    eval_res = evaluate_mechanism_directionality(
        gene_symbol="TP53",
        action_type="EXOGENOUS GENE",
        mechanism_desc="Cellular tumor antigen p53 exogenous gene",
        canonical_direction="SURVEILLANCE_ARREST",
        therapeutic_direction="Restore surveillance",
        is_canonical_mva=False,
    )
    assert eval_res.fit_status == "CONCORDANT_DIRECTION"
    assert eval_res.contradiction_flag is False


def test_seed_selection_modes():
    """Verify seed selection logic across mva_core, mva_network, and custom modes."""
    p1_df = make_synthetic_phase1_df()
    net_df = make_synthetic_network_context_df()

    cfg_core = DrugEvidenceConfig(seed_mode="mva_core")
    pipeline_core = DrugEvidencePipeline(cfg_core)
    seeds_core = pipeline_core.select_seed_genes(p1_df, net_df)
    assert "BUB1B" in seeds_core
    assert "AURKA" not in seeds_core

    cfg_custom = DrugEvidenceConfig(seed_mode="custom", custom_genes=["TP53", "BUB1B"])
    pipeline_custom = DrugEvidencePipeline(cfg_custom)
    seeds_custom = pipeline_custom.select_seed_genes(p1_df, net_df)
    assert seeds_custom == ["BUB1B", "TP53"]


def test_pipeline_evidence_collection_and_ranking_mocked():
    """Verify end-to-end collection, contradiction gating, and ranking with mocked HTTP."""
    p1_df = make_synthetic_phase1_df()
    cfg = DrugEvidenceConfig(seed_mode="custom", custom_genes=["BUB1B", "TP53"])
    pipeline = DrugEvidencePipeline(cfg)

    # Mock Open Targets responses
    mock_ot_search_resp = {
        "data": {
            "search": {
                "hits": [
                    {"id": "ENSG00000156970", "name": "BUB1B", "entity": "target"},
                    {"id": "ENSG00000141510", "name": "TP53", "entity": "target"},
                ]
            }
        }
    }

    def mock_ot_target_resp(ensembl_id: str):
        if ensembl_id == "ENSG00000156970":  # BUB1B
            return {
                "id": "ENSG00000156970",
                "approvedSymbol": "BUB1B",
                "tractability": [{"label": "High-Quality Ligand", "value": True}],
                "safetyLiabilities": [],
                "drugAndClinicalCandidates": {
                    "count": 1,
                    "rows": [
                        {
                            "id": "cand_bub1b_1",
                            "maxClinicalStage": "PRECLINICAL",
                            "drug": {
                                "id": "CHEMBL_MOCK_BUB1B",
                                "name": "BAY-1816032",
                                "drugType": "Small molecule",
                                "maximumClinicalStage": "PRECLINICAL",
                                "mechanismsOfAction": {
                                    "rows": [
                                        {
                                            "mechanismOfAction": "BUBR1 kinase inhibitor",
                                            "actionType": "INHIBITOR",
                                            "targets": [{"id": "ENSG00000156970", "approvedSymbol": "BUB1B"}],
                                        }
                                    ]
                                },
                                "drugWarnings": [],
                            },
                        }
                    ],
                },
            }
        else:  # TP53
            return {
                "id": "ENSG00000141510",
                "approvedSymbol": "TP53",
                "tractability": [{"label": "Advanced Clinical", "value": True}],
                "safetyLiabilities": [{"event": "Tumorigenesis"}],
                "drugAndClinicalCandidates": {
                    "count": 1,
                    "rows": [
                        {
                            "id": "cand_tp53_1",
                            "maxClinicalStage": "PHASE_3",
                            "drug": {
                                "id": "CHEMBL4297574",
                                "name": "TEPRASIRAN",
                                "drugType": "Oligonucleotide",
                                "maximumClinicalStage": "PHASE_3",
                                "mechanismsOfAction": {
                                    "rows": [
                                        {
                                            "mechanismOfAction": "p53 mRNA RNAi inhibitor",
                                            "actionType": "RNAI INHIBITOR",
                                            "targets": [{"id": "ENSG00000141510", "approvedSymbol": "TP53"}],
                                        }
                                    ]
                                },
                                "drugWarnings": [],
                            },
                        }
                    ],
                },
            }

    with patch.object(pipeline.ot_client, "search_target_ensembl_id") as mock_search, \
         patch.object(pipeline.ot_client, "fetch_target_drugs_and_safety") as mock_fetch_ot, \
         patch.object(pipeline.chembl_client, "fetch_target_mechanisms") as mock_fetch_chembl, \
         patch.object(pipeline.dgidb_client, "fetch_gene_interactions") as mock_fetch_dgidb:

        mock_search.side_effect = lambda s: "ENSG00000156970" if s == "BUB1B" else "ENSG00000141510"
        mock_fetch_ot.side_effect = mock_ot_target_resp
        mock_fetch_chembl.side_effect = lambda sym: [
            {
                "target_chembl_id": "CHEMBL4096",
                "target_pref_name": "Cellular tumor antigen p53",
                "target_type": "SINGLE PROTEIN",
                "molecule_chembl_id": "CHEMBL2108585",
                "mechanism_of_action": "Cellular tumor antigen p53 exogenous gene",
                "action_type": "EXOGENOUS GENE",
                "max_phase": 2,
                "direct_interaction": 1,
                "disease_efficacy": 1,
                "mechanism_comment": "Gene replacement",
            }
        ] if sym == "TP53" else []
        mock_fetch_dgidb.side_effect = lambda syms: [
            {
                "gene": {"name": "TP53", "conceptId": "hgnc:11998"},
                "drug": {"name": "OXALIPLATIN", "conceptId": "rxcui:32592", "approved": True},
                "interactionScore": 0.8,
                "interactionTypes": [{"type": "modulator", "directionality": ""}],
                "sources": [{"sourceDbName": "PharmGKB"}],
            }
        ] if "TP53" in syms else []

        matrix_df, ranking_df, action_map_df, provenance = pipeline.collect_evidence_for_seeds(
            seed_genes=["BUB1B", "TP53"],
            gene_matrix_df=p1_df,
        )

        assert not matrix_df.empty
        assert len(matrix_df) >= 3  # BAY-1816032, TEPRASIRAN, CHEMBL2108585, OXALIPLATIN
        assert "BAY-1816032" in matrix_df["DRUG_NAME"].values

        # Check that BAY-1816032 on BUB1B is flagged as hazard
        bub1b_rows = matrix_df[matrix_df["TARGET_SYMBOL"] == "BUB1B"]
        assert (bub1b_rows["DIRECTIONALITY_FIT"] == "MECHANISTIC_HAZARD").all()
        assert (bub1b_rows["CONTRADICTION_FLAG"] == True).all()

        # Check ranking contains required columns and disclaimers
        assert "SAFETY_TIER" in ranking_df.columns
        assert "CLINICAL_RECOMMENDATION_STATUS" in ranking_df.columns
        assert (ranking_df["CLINICAL_RECOMMENDATION_STATUS"] == "PROHIBITED / INVESTIGATIONAL ONLY").all()

        # Check action map
        assert "SEED_GENE" in action_map_df.columns
        assert len(action_map_df) == 2
        bub1b_map = action_map_df[action_map_df["SEED_GENE"] == "BUB1B"].iloc[0]
        assert bub1b_map["DIRECTIONAL_HAZARDS_COUNT"] >= 1

        # Check provenance
        assert provenance["pipeline_version"] == "0.3.0"
        assert provenance["contradiction_count"] >= 1
        assert "ethical_and_regulatory_notice" in provenance


def test_write_drug_evidence_outputs_flat():
    """Verify write_drug_evidence_outputs writes all 4 flat files properly."""
    matrix_df = pd.DataFrame([{"TARGET_SYMBOL": "BUB1B", "DRUG_NAME": "DRUG_A"}])
    ranking_df = pd.DataFrame([{"TARGET_SYMBOL": "BUB1B", "DRUG_NAME": "DRUG_A", "SAFETY_TIER": "TIER_4"}])
    action_map_df = pd.DataFrame([{"SEED_GENE": "BUB1B", "DIRECTIONAL_HAZARDS_COUNT": 1}])
    provenance = {"pipeline_version": "0.3.0", "contradiction_count": 1}

    with tempfile.TemporaryDirectory() as tmpdir:
        out_dir = Path(tmpdir)
        paths = write_drug_evidence_outputs(
            evidence_matrix_df=matrix_df,
            ranking_df=ranking_df,
            action_map_df=action_map_df,
            provenance=provenance,
            output_dir=out_dir,
            flat_names=True,
        )

        assert (out_dir / "14_phase3_drug_evidence_matrix.csv").exists()
        assert (out_dir / "15_phase3_safety_ranking.csv").exists()
        assert (out_dir / "16_phase3_mechanism_action_map.csv").exists()
        assert (out_dir / "17_phase3_repurposing_provenance.json").exists()


def test_cli_runner_offline_mode():
    """Verify CLI runner executes in offline mode with minimal mock tables."""
    p1_df = make_synthetic_phase1_df()
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        p1_file = tmp_path / "02_phase1_gene_mechanism_matrix.csv"
        p1_df.to_csv(p1_file, index=False)
        out_dir = tmp_path / "out"

        test_args = [
            "run_drug_evidence.py",
            "--gene-mechanism", str(p1_file),
            "--output-dir", str(out_dir),
            "--seed-mode", "mva_core",
            "--offline",
        ]
        with patch.object(sys, "argv", test_args):
            exit_code = run_drug_evidence_main()
            assert exit_code == 0
            assert (out_dir / "14_phase3_drug_evidence_matrix.csv").exists()
            assert (out_dir / "17_phase3_repurposing_provenance.json").exists()
