"""Unit tests for Phase 3 graph-derived query gene panel builder and drug evidence wiring.

Operates exclusively on synthetic non-patient fixtures, deterministic mock tables,
and mocked HTTP responses to guarantee zero patient data leakage and offline execution.
"""

from __future__ import annotations

import json
from pathlib import Path
import re
import sys
import tempfile
from typing import Any
from unittest.mock import MagicMock, patch

REPO_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import pandas as pd
import pytest

from mva_hackathon.drug_evidence import (
    DrugEvidenceConfig,
    DrugEvidencePipeline,
)
from mva_hackathon.mva_gene_sets import (
    CANONICAL_MVA_GENES,
    CORE_MVA_CAUSAL_GENES,
    EXTENDED_MVA_GENES,
)
from mva_hackathon.query_gene_panel import (
    DOI_PATTERN,
    PMID_PATTERN,
    QueryGenePanelBuilder,
    QueryGenePanelConfig,
    verify_literature_via_pubmed,
    write_query_panel_outputs,
)
from run_query_gene_panel import main as run_query_gene_panel_main


def make_synthetic_track2_fixtures(tmp_dir: Path) -> dict[str, Path]:
    """Create minimal synthetic Track 2 tables for testing panel derivation."""
    # 1. Phase 1 gene mechanism matrix (BUB1B present, CEP57 & TRIP13 intentionally absent)
    p1_df = pd.DataFrame([
        {
            "SYMBOL": "BUB1B",
            "GENE_NAME": "BUB1 mitotic checkpoint serine/threonine kinase B",
            "PATHWAY": "Spindle Assembly Checkpoint",
            "BIOLOGICAL_ROLE": "Core SAC kinase",
            "MVA_MECHANISM_LINK": "Canonical MVA Type 1",
            "CANONICAL_PATHWAY_DIRECTION": "DOWN (Checkpoint arrest failure)",
            "THERAPEUTIC_DIRECTIONALITY": "Compensate impaired checkpoint",
            "INHERITANCE_REQUIREMENT": "Biallelic",
            "TOTAL_VARIANTS_EVALUATED": 2,
            "TOTAL_PASSING_CANDIDATE_VARIANTS": 2,
            "COMPOUND_HET_ARCHETYPE_FIT": True,
            "ARCHETYPE_DESCRIPTION": "CLASSIC_MVA_ARCHETYPE",
            "PHASE_GATE_STATUS": "REQUIRES_PHASE_VERIFICATION",
        },
        {
            "SYMBOL": "MAD1L1",
            "GENE_NAME": "MAD1 mitotic arrest deficient 1 like 1",
            "PATHWAY": "Kinetochore Signaling",
            "BIOLOGICAL_ROLE": "Kinetochore scaffold",
            "MVA_MECHANISM_LINK": "Reported systemic mosaic aneuploidy",
            "CANONICAL_PATHWAY_DIRECTION": "DOWN",
            "THERAPEUTIC_DIRECTIONALITY": "Strengthen checkpoint",
            "INHERITANCE_REQUIREMENT": "Biallelic",
            "TOTAL_VARIANTS_EVALUATED": 2,
            "TOTAL_PASSING_CANDIDATE_VARIANTS": 2,
            "COMPOUND_HET_ARCHETYPE_FIT": True,
            "ARCHETYPE_DESCRIPTION": "CLASSIC_MVA_ARCHETYPE",
            "PHASE_GATE_STATUS": "REQUIRES_PHASE_VERIFICATION",
        },
        {
            "SYMBOL": "CAND1",
            "GENE_NAME": "Candidate Gene 1",
            "PATHWAY": "Mitotic Cell Cycle",
            "BIOLOGICAL_ROLE": "Spindle organization",
            "MVA_MECHANISM_LINK": "Secondary candidate",
            "CANONICAL_PATHWAY_DIRECTION": "DOWN",
            "THERAPEUTIC_DIRECTIONALITY": "Mitotic stabilization",
            "INHERITANCE_REQUIREMENT": "Recessive",
            "TOTAL_VARIANTS_EVALUATED": 3,
            "TOTAL_PASSING_CANDIDATE_VARIANTS": 2,
            "COMPOUND_HET_ARCHETYPE_FIT": True,
            "ARCHETYPE_DESCRIPTION": "CLASSIC_MVA_ARCHETYPE",
            "PHASE_GATE_STATUS": "REQUIRES_PHASE_VERIFICATION",
        },
    ])

    # 2. Phase 2 seed genes (testing symbol hygiene)
    p2_seeds_df = pd.DataFrame([
        {"INPUT_SYMBOL": "BUB1B", "HGNC_SYMBOL": "BUB1B", "HGNC_STATUS": "mapped", "DESCRIPTION": "mitotic kinase"},
        {"INPUT_SYMBOL": "MAD1L1", "HGNC_SYMBOL": "MAD1L1", "HGNC_STATUS": "mapped", "DESCRIPTION": "checkpoint scaffold"},
        {"INPUT_SYMBOL": "NEIGHBOR_A", "HGNC_SYMBOL": "NEIGHBOR_A", "HGNC_STATUS": "mapped", "DESCRIPTION": "centrosomal interactor"},
        {"INPUT_SYMBOL": "NEIGHBOR_B", "HGNC_SYMBOL": "NEIGHBOR_B", "HGNC_STATUS": "mapped", "DESCRIPTION": "spindle member"},
        {"INPUT_SYMBOL": "CAND1", "HGNC_SYMBOL": "CAND1", "HGNC_STATUS": "mapped", "DESCRIPTION": "mitotic regulator", "TOTAL_PASSING_CANDIDATE_VARIANTS": 2, "COMPOUND_HET_ARCHETYPE_FIT": True},
        # Hygiene drop cases
        {"INPUT_SYMBOL": "UNMAPPED_GENE", "HGNC_SYMBOL": "UNMAPPED_GENE", "HGNC_STATUS": "unmapped_retained", "DESCRIPTION": "unmapped locus"},
        {"INPUT_SYMBOL": "HLA-DRA", "HGNC_SYMBOL": "HLA-DRA", "HGNC_STATUS": "mapped", "DESCRIPTION": "major histocompatibility complex"},
        {"INPUT_SYMBOL": "GENE1-GENE2", "HGNC_SYMBOL": "GENE1-GENE2", "HGNC_STATUS": "mapped", "DESCRIPTION": "readthrough transcript"},
        {"INPUT_SYMBOL": "TESTP1", "HGNC_SYMBOL": "TESTP1", "HGNC_STATUS": "mapped", "DESCRIPTION": "ribosomal protein pseudogene"},
        # Sensory receptor negative control candidates
        {"INPUT_SYMBOL": "OR10A1", "HGNC_SYMBOL": "OR10A1", "HGNC_STATUS": "mapped", "DESCRIPTION": "olfactory receptor family 10"},
        {"INPUT_SYMBOL": "TAS2R4", "HGNC_SYMBOL": "TAS2R4", "HGNC_STATUS": "mapped", "DESCRIPTION": "taste 2 receptor member 4"},
    ])

    # 3. Phase 2 string edges (NEIGHBOR_A touches BUB1B with curated support; NEIGHBOR_B touches NEIGHBOR_A)
    edges_df = pd.DataFrame([
        {
            "PREFERRED_NAME_A": "BUB1B",
            "PREFERRED_NAME_B": "NEIGHBOR_A",
            "SCORE": 0.850,
            "CURATED_SUPPORT": "YES",
        },
        {
            "PREFERRED_NAME_A": "NEIGHBOR_A",
            "PREFERRED_NAME_B": "NEIGHBOR_B",
            "SCORE": 0.750,
            "CURATED_SUPPORT": "NO",
        },
    ])

    # 4. Phase 2 pathway enrichment
    pw_df = pd.DataFrame([
        {
            "PATHWAY_NAME": "Mitotic Spindle Assembly Checkpoint",
            "HIT_GENES": "BUB1B;NEIGHBOR_A;CAND1",
            "FDR": 0.001,
            "PATHWAY_BROWSER_URL": "https://reactome.org/test",
        }
    ])

    # 5. Phase 2 network gene context
    ctx_df = pd.DataFrame([
        {"HGNC_SYMBOL": "BUB1B", "NETWORK_DEGREE": 4, "COMPONENT_ID": 5, "IS_ISOLATED": False, "STRING_INTERACTORS": "NEIGHBOR_A"},
        {"HGNC_SYMBOL": "MAD1L1", "NETWORK_DEGREE": 2, "COMPONENT_ID": 5, "IS_ISOLATED": False, "STRING_INTERACTORS": ""},
        {"HGNC_SYMBOL": "NEIGHBOR_A", "NETWORK_DEGREE": 3, "COMPONENT_ID": 5, "IS_ISOLATED": False, "STRING_INTERACTORS": "BUB1B;NEIGHBOR_B"},
        {"HGNC_SYMBOL": "NEIGHBOR_B", "NETWORK_DEGREE": 1, "COMPONENT_ID": 5, "IS_ISOLATED": False, "STRING_INTERACTORS": "NEIGHBOR_A"},
        {"HGNC_SYMBOL": "CAND1", "NETWORK_DEGREE": 2, "COMPONENT_ID": 8, "IS_ISOLATED": False, "STRING_INTERACTORS": ""},
        {"HGNC_SYMBOL": "OR10A1", "NETWORK_DEGREE": 0, "COMPONENT_ID": 99, "IS_ISOLATED": True, "STRING_INTERACTORS": ""},
        {"HGNC_SYMBOL": "TAS2R4", "NETWORK_DEGREE": 0, "COMPONENT_ID": 100, "IS_ISOLATED": True, "STRING_INTERACTORS": ""},
    ])

    # 6. Phase 2 network nodes (Component 5 is rank 4 size 4 in this toy graph!)
    nodes_df = pd.DataFrame([
        {"HGNC_SYMBOL": "BUB1B", "COMPONENT_RANK": 4, "NETWORK_DEGREE": 4},
        {"HGNC_SYMBOL": "MAD1L1", "COMPONENT_RANK": 4, "NETWORK_DEGREE": 2},
        {"HGNC_SYMBOL": "NEIGHBOR_A", "COMPONENT_RANK": 4, "NETWORK_DEGREE": 3},
        {"HGNC_SYMBOL": "NEIGHBOR_B", "COMPONENT_RANK": 4, "NETWORK_DEGREE": 1},
        {"HGNC_SYMBOL": "CAND1", "COMPONENT_RANK": 7, "NETWORK_DEGREE": 2},
    ])

    # 7. Phase 2 component summary (Rank 4 contains BUB1B and MAD1L1)
    comp_df = pd.DataFrame([
        {
            "COMPONENT_RANK": 4,
            "COMPONENT_SIZE": 4,
            "EDGE_COUNT": 2,
            "HAS_MVA_KNOWN_GENE": True,
            "MVA_GENES": "BUB1B;MAD1L1",
            "GENE_LIST": "BUB1B;MAD1L1;NEIGHBOR_A;NEIGHBOR_B",
        },
        {
            "COMPONENT_RANK": 7,
            "COMPONENT_SIZE": 1,
            "EDGE_COUNT": 0,
            "HAS_MVA_KNOWN_GENE": False,
            "MVA_GENES": "",
            "GENE_LIST": "CAND1",
        },
    ])

    paths = {
        "02_phase1_gene_mechanism_matrix.csv": p1_df,
        "03_phase2_seed_genes.csv": p2_seeds_df,
        "04_phase2_string_edges.csv": edges_df,
        "06_phase2_pathway_enrichment.csv": pw_df,
        "07_phase2_network_gene_context.csv": ctx_df,
        "09_phase2_network_nodes.csv": nodes_df,
        "10_phase2_component_summary.csv": comp_df,
    }

    res = {}
    for fname, df in paths.items():
        fp = tmp_dir / fname
        df.to_csv(fp, index=False)
        res[fname] = fp

    return res


def test_component_detection_anchor_based():
    """Verify component detection finds the anchor component dynamically rather than hard-coding rank 2."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        make_synthetic_track2_fixtures(tmp_path)

        cfg = QueryGenePanelConfig(track2_dir=tmp_path)
        builder = QueryGenePanelBuilder(cfg)
        tables = builder.load_track2_tables()

        comp_info = builder.detect_mva_focus_component(tables["p2_components"], tables["p2_nodes"])
        # In our toy fixture, the anchor component was intentionally assigned rank 4!
        assert comp_info["component_rank"] == 4
        assert "BUB1B" in comp_info["anchors"]
        assert comp_info["size"] == 4
        assert comp_info["edge_count"] == 2


def test_core_injection_and_gap_flags():
    """Verify all 7 core MVA genes are included and missing ones (CEP57/TRIP13) have PHASE1_ABSENT."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        make_synthetic_track2_fixtures(tmp_path)

        cfg = QueryGenePanelConfig(track2_dir=tmp_path)
        builder = QueryGenePanelBuilder(cfg)
        panel_df, _, prov = builder.build_panel()

        for cg in CORE_MVA_CAUSAL_GENES:
            assert cg in panel_df["QUERY_GENE"].values

        # Check CEP57 & TRIP13 injection
        cep57 = panel_df[panel_df["QUERY_GENE"] == "CEP57"].iloc[0]
        assert bool(cep57["PHASE1_PRESENT"]) is False
        assert "PHASE1_ABSENT" in cep57["EVIDENCE_GAP_FLAGS"]
        assert cep57["GENE_CATEGORY"] == "CORE_MVA_CAUSAL"
        assert cep57["MVA_GENE_CLASS"] == "CANONICAL"

        # Check BUB1B
        bub1b = panel_df[panel_df["QUERY_GENE"] == "BUB1B"].iloc[0]
        assert bool(bub1b["PHASE1_PRESENT"]) is True
        assert "PHASE1_ABSENT" not in bub1b["EVIDENCE_GAP_FLAGS"]


def test_one_hop_vs_component_member_tiering():
    """Verify ONE_HOP_OF_CORE vs TWO_HOP_OR_COMPONENT_MEMBER tiering in graph neighbors."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        make_synthetic_track2_fixtures(tmp_path)

        cfg = QueryGenePanelConfig(track2_dir=tmp_path)
        builder = QueryGenePanelBuilder(cfg)
        panel_df, _, _ = builder.build_panel()

        neighbor_a = panel_df[panel_df["QUERY_GENE"] == "NEIGHBOR_A"].iloc[0]
        assert neighbor_a["GENE_CATEGORY"] == "GRAPH_DERIVED_NEIGHBOR"
        assert neighbor_a["NEIGHBOR_HOP"] == "ONE_HOP_OF_CORE"
        assert neighbor_a["CURATED_EDGE_SUPPORT"] == "YES"
        assert neighbor_a["BEST_EDGE_PARTNER"] == "BUB1B"

        neighbor_b = panel_df[panel_df["QUERY_GENE"] == "NEIGHBOR_B"].iloc[0]
        assert neighbor_b["GENE_CATEGORY"] == "GRAPH_DERIVED_NEIGHBOR"
        assert neighbor_b["NEIGHBOR_HOP"] == "TWO_HOP_OR_COMPONENT_MEMBER"


def test_hygiene_quarantine():
    """Verify Gate 0 quarantines unmapped, HLA, readthrough, and pseudogenes."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        make_synthetic_track2_fixtures(tmp_path)

        cfg = QueryGenePanelConfig(track2_dir=tmp_path)
        builder = QueryGenePanelBuilder(cfg)
        panel_df, _, prov = builder.build_panel()

        symbols_in_panel = set(panel_df["QUERY_GENE"].values)
        assert "UNMAPPED_GENE" not in symbols_in_panel
        assert "HLA-DRA" not in symbols_in_panel
        assert "GENE1-GENE2" not in symbols_in_panel
        assert "TESTP1" not in symbols_in_panel

        reasons = prov["gate0_dropped_ledger_counts"]
        assert "UNMAPPED_HGNC_STATUS" in reasons
        assert "HLA_COMPLEX_NON_QUERY" in reasons
        assert "READTHROUGH_OR_NON_CODING_TRANSCRIPT" in reasons
        assert "PSEUDOGENE" in reasons


def test_negative_control_derivation():
    """Verify programmatically derived negative controls from isolated sensory receptors."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        make_synthetic_track2_fixtures(tmp_path)

        cfg = QueryGenePanelConfig(track2_dir=tmp_path, control_negative_max=2)
        builder = QueryGenePanelBuilder(cfg)
        panel_df, _, _ = builder.build_panel()

        neg_controls = panel_df[panel_df["CONTROL_ROLE"] == "NEGATIVE"]
        assert len(neg_controls) >= 1
        for _, r in neg_controls.iterrows():
            assert r["GENE_CATEGORY"] == "DRUG_TARGET_CONTROL"
            assert r["CONTROL_EXPECTATION"] == "EXPECTED_NO_HITS"
            assert r["IS_ISOLATED"] is True


def test_doi_pmid_strict_validation():
    """Verify that all non-empty DOIs and PMIDs match strict regex and contain no placeholders."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        make_synthetic_track2_fixtures(tmp_path)

        cfg = QueryGenePanelConfig(track2_dir=tmp_path)
        builder = QueryGenePanelBuilder(cfg)
        panel_df, citations_df, _ = builder.build_panel()

        for doi in panel_df["PRIMARY_LITERATURE_DOI"]:
            if doi != "":
                assert DOI_PATTERN.match(str(doi)), f"Invalid DOI format: {doi}"
                assert "NA" not in doi and "TBD" not in doi

        for pmid in panel_df["PRIMARY_LITERATURE_PMID"]:
            if pmid != "":
                assert PMID_PATTERN.match(str(pmid)), f"Invalid PMID format: {pmid}"

        # Empty string not fabricated policy
        for col in ["PRIMARY_LITERATURE_PMID", "PRIMARY_LITERATURE_DOI", "MVA_DISEASE_OMIM"]:
            val = panel_df[col].iloc[0]
            assert isinstance(val, str)


def test_write_outputs_and_provenance():
    """Verify write_query_panel_outputs writes flat files within 8MB limit."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        make_synthetic_track2_fixtures(tmp_path)

        cfg = QueryGenePanelConfig(track2_dir=tmp_path, output_dir=tmp_path)
        builder = QueryGenePanelBuilder(cfg)
        panel_df, citations_df, prov = builder.build_panel()

        written = write_query_panel_outputs(
            panel_df=panel_df,
            provenance=prov,
            citations_df=citations_df,
            output_dir=tmp_path,
            flat_names=True,
        )

        assert written["query_panel"].exists()
        assert written["provenance"].exists()
        assert written["citations"].exists()

        assert written["query_panel"].name == "18_phase3_query_gene_panel.csv"
        assert written["provenance"].name == "18_phase3_query_panel_provenance.json"
        assert written["citations"].name == "19_phase3_query_panel_citations.csv"

        for p in written.values():
            assert p.stat().st_size < 8 * 1024 * 1024


def test_verify_literature_via_pubmed_mocked():
    """Verify NCBI E-utilities lookup with mocked matching and non-matching responses."""
    mock_esearch_resp = json.dumps({
        "esearchresult": {"idlist": ["12345678"]}
    }).encode("utf-8")

    mock_esummary_matching = json.dumps({
        "result": {
            "12345678": {
                "title": "Aneuploidy and chromosome instability in BUB1B syndrome.",
                "source": "Nat Genet",
                "pubdate": "2020 Nov",
                "articleids": [{"idtype": "doi", "value": "10.1038/ng.1234"}],
            }
        }
    }).encode("utf-8")

    mock_esummary_non_matching = json.dumps({
        "result": {
            "12345678": {
                "title": "Unrelated cardiovascular review in elderly patients.",
                "source": "Heart",
                "pubdate": "2020",
                "articleids": [{"idtype": "doi", "value": "10.1038/heart.1234"}],
            }
        }
    }).encode("utf-8")

    # 1. Matching case
    with patch("urllib.request.urlopen") as mock_url:
        mock_resp1 = MagicMock()
        mock_resp1.read.return_value = mock_esearch_resp
        mock_resp1.__enter__.return_value = mock_resp1

        mock_resp2 = MagicMock()
        mock_resp2.read.return_value = mock_esummary_matching
        mock_resp2.__enter__.return_value = mock_resp2

        mock_url.side_effect = [mock_resp1, mock_resp2]

        res = verify_literature_via_pubmed(
            term="test term",
            expected_title_words=["aneuploidy", "bub1b"],
            polite_delay=0.0,
        )
        assert res["verified"] is True
        assert res["pmid"] == "12345678"
        assert res["doi"] == "10.1038/ng.1234"

    # 2. Non-matching title case
    with patch("urllib.request.urlopen") as mock_url:
        mock_resp1 = MagicMock()
        mock_resp1.read.return_value = mock_esearch_resp
        mock_resp1.__enter__.return_value = mock_resp1

        mock_resp2 = MagicMock()
        mock_resp2.read.return_value = mock_esummary_non_matching
        mock_resp2.__enter__.return_value = mock_resp2

        mock_url.side_effect = [mock_resp1, mock_resp2]

        res = verify_literature_via_pubmed(
            term="test term",
            expected_title_words=["aneuploidy", "bub1b"],
            polite_delay=0.0,
        )
        assert res["verified"] is False
        assert res["pmid"] == ""


def test_drug_evidence_query_panel_wiring_mocked():
    """Verify that DrugEvidencePipeline correctly loads a query panel and passes metadata through."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        make_synthetic_track2_fixtures(tmp_path)

        # Build synthetic panel
        cfg_panel = QueryGenePanelConfig(track2_dir=tmp_path, output_dir=tmp_path)
        builder = QueryGenePanelBuilder(cfg_panel)
        panel_df, citations_df, prov = builder.build_panel()
        written = write_query_panel_outputs(panel_df, prov, citations_df, tmp_path, flat_names=True)

        p1_df = pd.read_csv(tmp_path / "02_phase1_gene_mechanism_matrix.csv")

        # Configure DrugEvidencePipeline with --query-panel
        d_cfg = DrugEvidenceConfig(
            seed_mode="query_panel",
            query_panel_path=written["query_panel"],
            offline_mode=True,
        )
        pipeline = DrugEvidencePipeline(d_cfg)
        seeds = pipeline.select_seed_genes(p1_df)

        assert len(seeds) == len(panel_df)
        assert "BUB1B" in seeds
        assert "CEP57" in seeds  # Core gene injected even though absent from Phase 1!

        # Mock collection with empty offline data
        mat_df, rank_df, act_df, d_prov = pipeline.collect_evidence_for_seeds(
            seed_genes=seeds,
            gene_matrix_df=p1_df,
        )

        # Action map must include all queried genes and report NO_DRUG_EVIDENCE for zero-hit genes
        assert len(act_df) == len(seeds)
        assert "GENE_CATEGORY" in act_df.columns
        assert "DRUG_ROWS_BY_SOURCE" in act_df.columns
        cep57_act = act_df[act_df["SEED_GENE"] == "CEP57"].iloc[0]
        assert cep57_act["RESEARCH_STATUS"] == "NO_DRUG_EVIDENCE"
        assert cep57_act["DRUG_ROWS_BY_SOURCE"] == "OPEN_TARGETS:0;CHEMBL:0;DGIDB:0"

        # Provenance sanity reporting
        assert d_prov["seed_mode"] == "query_panel"
        assert "genes_with_zero_hits" in d_prov
        assert len(d_prov["genes_with_zero_hits"]) == len(seeds)
        assert "negative_control_sanity" in d_prov


def test_cli_runner_offline():
    """Verify scripts/run_query_gene_panel.py executes cleanly offline."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        make_synthetic_track2_fixtures(tmp_path)
        out_dir = tmp_path / "out"

        test_args = [
            "run_query_gene_panel.py",
            "--track2-dir", str(tmp_path),
            "--output-dir", str(out_dir),
            "--flat-names",
        ]
        with patch.object(sys, "argv", test_args):
            exit_code = run_query_gene_panel_main()
            assert exit_code == 0
            assert (out_dir / "18_phase3_query_gene_panel.csv").exists()
            assert (out_dir / "18_phase3_query_panel_provenance.json").exists()
            assert (out_dir / "19_phase3_query_panel_citations.csv").exists()
