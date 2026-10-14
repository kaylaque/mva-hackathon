"""Unit tests for Phase 2 PPI and pathway network context engine.

All tests operate exclusively with synthetic fixtures and mocked external APIs
to guarantee determinism, data privacy, and absence of external network dependency.
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

from mva_hackathon.network_analysis import (
    HttpClient,
    NetworkAnalysisConfig,
    NetworkContextEngine,
    write_network_outputs,
)
from run_network_analysis import main as run_network_analysis_main


def make_synthetic_gene_mechanism_df() -> pd.DataFrame:
    """Create a minimal synthetic non-patient gene mechanism matrix."""
    return pd.DataFrame([
        {
            "SYMBOL": "BUB1B",
            "GENE_NAME": "BUB1 mitotic checkpoint serine/threonine kinase",
            "PATHWAY": "Spindle Assembly Checkpoint (SAC) / Mitotic Cell Cycle",
            "BIOLOGICAL_ROLE": "Core SAC kinase",
            "MVA_MECHANISM_LINK": "Biallelic LoF causes MVA type 1 (OMIM 257300)",
            "CANONICAL_PATHWAY_DIRECTION": "REDUCED_SAC_FUNCTION",
            "THERAPEUTIC_DIRECTIONALITY": "RESTORE_CHECKPOINT_FUNCTION",
            "INHERITANCE_REQUIREMENT": "Biallelic",
            "TOTAL_VARIANTS_EVALUATED": 2,
            "TOTAL_PASSING_CANDIDATE_VARIANTS": 2,
            "COMPOUND_HET_ARCHETYPE_FIT": True,
            "ARCHETYPE_DESCRIPTION": "CLASSIC_MVA_ARCHETYPE (Null/Trun + Missense/Splice)",
            "PHASE_GATE_STATUS": "REQUIRES_PHASE_VERIFICATION",
        },
        {
            "SYMBOL": "MAD1L1",
            "GENE_NAME": "MAD1 mitotic checkpoint component",
            "PATHWAY": "SAC Initiation & Kinetochore Signaling",
            "BIOLOGICAL_ROLE": "SAC complex component",
            "MVA_MECHANISM_LINK": "Reported in mosaic variegated aneuploidy",
            "CANONICAL_PATHWAY_DIRECTION": "REDUCED_SAC_FUNCTION",
            "THERAPEUTIC_DIRECTIONALITY": "RESTORE_CHECKPOINT_FUNCTION",
            "INHERITANCE_REQUIREMENT": "Biallelic or dominant-negative",
            "TOTAL_VARIANTS_EVALUATED": 1,
            "TOTAL_PASSING_CANDIDATE_VARIANTS": 1,
            "COMPOUND_HET_ARCHETYPE_FIT": False,
            "ARCHETYPE_DESCRIPTION": "NO_ARCHETYPE",
            "PHASE_GATE_STATUS": "SINGLE_HIT_OR_UNPASSING",
        },
        {
            "SYMBOL": "TRIP13",
            "GENE_NAME": "thyroid hormone receptor interactor 13",
            "PATHWAY": "SAC Silencing and Mitotic Checkpoint Complex Turnover",
            "BIOLOGICAL_ROLE": "Pachytene checkpoint / MCC turnover",
            "MVA_MECHANISM_LINK": "Canonical MVA gene",
            "CANONICAL_PATHWAY_DIRECTION": "REDUCED_SAC_FUNCTION",
            "THERAPEUTIC_DIRECTIONALITY": "RESTORE_CHECKPOINT_FUNCTION",
            "INHERITANCE_REQUIREMENT": "Autosomal recessive",
            "TOTAL_VARIANTS_EVALUATED": 1,
            "TOTAL_PASSING_CANDIDATE_VARIANTS": 1,
            "COMPOUND_HET_ARCHETYPE_FIT": False,
            "ARCHETYPE_DESCRIPTION": "NO_ARCHETYPE",
            "PHASE_GATE_STATUS": "SINGLE_HIT_OR_UNPASSING",
        },
        {
            "SYMBOL": "GENE_ISOLATED",
            "GENE_NAME": "synthetic isolated gene",
            "PATHWAY": "General Cellular Maintenance",
            "BIOLOGICAL_ROLE": "Unknown",
            "MVA_MECHANISM_LINK": "NOT_IN_MVA_KNOWLEDGE_BASE",
            "CANONICAL_PATHWAY_DIRECTION": "UNKNOWN",
            "THERAPEUTIC_DIRECTIONALITY": "UNKNOWN",
            "INHERITANCE_REQUIREMENT": "Uncertain",
            "TOTAL_VARIANTS_EVALUATED": 1,
            "TOTAL_PASSING_CANDIDATE_VARIANTS": 0,
            "COMPOUND_HET_ARCHETYPE_FIT": False,
            "ARCHETYPE_DESCRIPTION": "NO_ARCHETYPE",
            "PHASE_GATE_STATUS": "SINGLE_HIT_OR_UNPASSING",
        },
    ])


def mock_response(status_code: int = 200, json_data: Any = None, text: str = "") -> requests.Response:
    resp = requests.Response()
    resp.status_code = status_code
    if json_data is not None:
        resp._content = json.dumps(json_data).encode("utf-8")
    else:
        resp._content = text.encode("utf-8")
    return resp


class MockNetworkHttpClient(HttpClient):
    """Synthetic HTTP client with canned responses for STRING, Reactome, and g:Profiler."""

    def __init__(self, config: NetworkAnalysisConfig) -> None:
        super().__init__(config=config)
        self.call_log: list[tuple[str, str]] = []

    def get(self, url: str, params: dict[str, Any] | None = None) -> requests.Response | None:
        self.call_log.append(("GET", url))
        # Route on exact endpoint suffixes so the version-pinned network URL
        # (https://version-12-5.string-db.org/api/json/network) is not caught
        # by a naive "version" substring match.
        if url.endswith("/api/json/version"):
            return mock_response(200, json_data=[{"string_version": "12.5"}])
        if url.endswith("/api/json/network"):
            # Simulate real STRING behavior: HTTP 404 when NO identifier in the
            # request batch resolves (keeps the skip-failed-batch path tested).
            identifiers = str((params or {}).get("identifiers", ""))
            batch_symbols = {s.strip().upper() for s in identifiers.split("\n") if s.strip()}
            known_symbols = {"BUB1B", "MAD1L1", "TRIP13"}
            if batch_symbols and not (batch_symbols & known_symbols):
                missing = sorted(batch_symbols)[0]
                return mock_response(
                    404,
                    text=(
                        f'[{{ "Error" : "not found", "ErrorMessage" : '
                        f'"<p>Sorry, STRING did not find a protein called \'{missing}\' '
                        f"in the taxon '9606'.</p>\" }}"
                    ),
                )
            # Synthetic STRING interaction edges between BUB1B, MAD1L1, TRIP13
            edges = [
                {
                    "stringId_A": "9606.ENSP00000371077",
                    "stringId_B": "9606.ENSP00000262104",
                    "preferredName_A": "BUB1B",
                    "preferredName_B": "TRIP13",
                    "ncbiTaxonId": 9606,
                    "score": 0.962,
                    "nscore": 0.0,
                    "fscore": 0.0,
                    "pscore": 0.0,
                    "ascore": 0.35,
                    "escore": 0.82,
                    "dscore": 0.90,
                    "tscore": 0.85,
                },
                {
                    "stringId_A": "9606.ENSP00000371077",
                    "stringId_B": "9606.ENSP00000288820",
                    "preferredName_A": "BUB1B",
                    "preferredName_B": "MAD1L1",
                    "ncbiTaxonId": 9606,
                    "score": 0.985,
                    "nscore": 0.0,
                    "fscore": 0.0,
                    "pscore": 0.0,
                    "ascore": 0.40,
                    "escore": 0.95,
                    "dscore": 0.92,
                    "tscore": 0.88,
                },
            ]
            return mock_response(200, json_data=edges)
        return mock_response(404, text="Not Found")

    def post(
        self,
        url: str,
        json_data: Any | None = None,
        data: Any | None = None,
        headers: dict[str, str] | None = None,
        params: dict[str, Any] | None = None,
    ) -> requests.Response | None:
        self.call_log.append(("POST", url))
        if "convert" in url:
            # g:Profiler convert response
            res = [
                {"incoming": "BUB1B", "converted": "BUB1B", "name": "BUB1B", "description": "BUB1 mitotic checkpoint serine/threonine kinase B"},
                {"incoming": "MAD1L1", "converted": "MAD1L1", "name": "MAD1L1", "description": "MAD1 mitotic arrest deficient like 1"},
                {"incoming": "TRIP13", "converted": "TRIP13", "name": "TRIP13", "description": "thyroid hormone receptor interactor 13"},
                {"incoming": "GENE_ISOLATED", "converted": "GENE_ISOLATED", "name": "GENE_ISOLATED", "description": "synthetic gene"},
            ]
            return mock_response(200, json_data={"result": res})
        if "projection" in url:
            # Reactome response
            pathways = [
                {
                    "stId": "R-HSA-141424",
                    "name": "Amplification of signal from the kinetochores",
                    "entities": {"found": 2, "total": 15, "pValue": 0.0001, "fdr": 0.001},
                    "species": {"name": "Homo sapiens"},
                },
                {
                    "stId": "R-HSA-141444",
                    "name": "Mitotic Spindle Checkpoint",
                    "entities": {"found": 3, "total": 28, "pValue": 0.00005, "fdr": 0.0008},
                    "species": {"name": "Homo sapiens"},
                },
            ]
            return mock_response(200, json_data={"pathways": pathways})
        if "profile" in url:
            # g:Profiler g:GOSt response
            gost_res = [
                {
                    "native": "GO:0007094",
                    "name": "mitotic spindle assembly checkpoint signaling",
                    "source": "GO:BP",
                    "p_value": 5.7e-6,
                    "term_size": 45,
                    "intersection_size": 3,
                    "intersections": [True, True, True, False],
                }
            ]
            return mock_response(200, json_data={"result": gost_res})
        return mock_response(404, text="Not Found")


def test_network_analysis_end_to_end_mocked() -> None:
    """Verify Phase 2 network analysis pipeline execution and metric generation."""
    config = NetworkAnalysisConfig(string_required_score=700, string_additional_nodes=0)
    mock_client = MockNetworkHttpClient(config)
    engine = NetworkContextEngine(config=config, http_client=mock_client)

    df_genes = make_synthetic_gene_mechanism_df()
    results = engine.run_analysis(df_genes)

    # Check presence of all expected tables
    for key in [
        "seed_genes_df",
        "string_edges_df",
        "network_summary_df",
        "pathway_enrichment_df",
        "network_gene_context_df",
        "provenance",
    ]:
        assert key in results, f"Missing output key: {key}"

    # Verify seed genes mapping
    seed_df = results["seed_genes_df"]
    assert len(seed_df) == 4
    assert set(seed_df["HGNC_SYMBOL"]) == {"BUB1B", "MAD1L1", "TRIP13", "GENE_ISOLATED"}

    # Verify STRING edges
    edges_df = results["string_edges_df"]
    assert len(edges_df) == 2
    assert "BUB1B" in edges_df["PREFERRED_NAME_A"].values or "BUB1B" in edges_df["PREFERRED_NAME_B"].values
    assert (edges_df["SCORE"] >= 0.700).all()

    # Verify Network Summary topology metrics
    summary_df = results["network_summary_df"]
    assert len(summary_df) == 1
    summary = summary_df.iloc[0]
    assert summary["TOTAL_SEED_GENES"] == 4
    assert summary["TOTAL_PPI_EDGES"] == 2
    assert summary["ISOLATED_GENES_COUNT"] == 1  # GENE_ISOLATED has degree 0
    assert summary["CONNECTED_COMPONENTS_COUNT"] == 2  # {BUB1B, TRIP13, MAD1L1} and {GENE_ISOLATED}
    assert summary["LARGEST_COMPONENT_SIZE"] == 3

    # Verify Pathway Enrichment
    enrichment_df = results["pathway_enrichment_df"]
    assert len(enrichment_df) >= 3
    pathway_names = enrichment_df["PATHWAY_NAME"].tolist()
    assert any("Mitotic Spindle Checkpoint" in name for name in pathway_names)
    assert any("mitotic spindle assembly checkpoint signaling" in name for name in pathway_names)

    # Verify Network Gene Context
    context_df = results["network_gene_context_df"]
    assert len(context_df) == 4
    bub1b_row = context_df[context_df["HGNC_SYMBOL"] == "BUB1B"].iloc[0]
    assert bub1b_row["NETWORK_DEGREE"] == 2
    assert not bub1b_row["IS_ISOLATED"]
    assert "TRIP13" in bub1b_row["STRING_INTERACTORS"]
    assert "MAD1L1" in bub1b_row["STRING_INTERACTORS"]

    isolated_row = context_df[context_df["HGNC_SYMBOL"] == "GENE_ISOLATED"].iloc[0]
    assert isolated_row["NETWORK_DEGREE"] == 0
    assert isolated_row["IS_ISOLATED"]


def test_score_channel_filtering() -> None:
    """Verify that score channel flags (min_escore, min_dscore) correctly filter PPI edges."""
    config = NetworkAnalysisConfig(
        string_required_score=700,
        min_escore=0.90,  # BUB1B-TRIP13 has escore 0.82, BUB1B-MAD1L1 has escore 0.95
    )
    mock_client = MockNetworkHttpClient(config)
    engine = NetworkContextEngine(config=config, http_client=mock_client)

    df_genes = make_synthetic_gene_mechanism_df()
    results = engine.run_analysis(df_genes)

    edges_df = results["string_edges_df"]
    # Only BUB1B-MAD1L1 should pass min_escore=0.90
    assert len(edges_df) == 1
    assert "MAD1L1" in [edges_df.iloc[0]["PREFERRED_NAME_A"], edges_df.iloc[0]["PREFERRED_NAME_B"]]


def test_tiny_and_zero_edge_lists() -> None:
    """Verify graceful handling for 1-gene input and disconnected networks."""
    config = NetworkAnalysisConfig()

    class DisconnectedClient(MockNetworkHttpClient):
        def get(self, url: str, params: dict[str, Any] | None = None) -> requests.Response | None:
            if "network" in url:
                return mock_response(200, json_data=[])  # Zero edges found
            return super().get(url, params)

    # Case 1: 1-gene input
    single_gene_df = pd.DataFrame([{"SYMBOL": "BUB1B", "MVA_MECHANISM_LINK": "Biallelic LoF causes MVA type 1"}])
    engine_single = NetworkContextEngine(config=config, http_client=MockNetworkHttpClient(config))
    res_single = engine_single.run_analysis(single_gene_df)
    assert len(res_single["string_edges_df"]) == 0
    assert res_single["network_summary_df"].iloc[0]["TOTAL_SEED_GENES"] == 1
    assert res_single["network_summary_df"].iloc[0]["TOTAL_PPI_EDGES"] == 0

    # Case 2: Multi-gene with zero edges returned
    engine_disco = NetworkContextEngine(config=config, http_client=DisconnectedClient(config))
    res_disco = engine_disco.run_analysis(make_synthetic_gene_mechanism_df())
    assert len(res_disco["string_edges_df"]) == 0
    assert res_disco["network_summary_df"].iloc[0]["TOTAL_PPI_EDGES"] == 0
    assert res_disco["network_summary_df"].iloc[0]["ISOLATED_GENES_COUNT"] == 4


def test_http_retry_and_failure_fallback() -> None:
    """Verify that HTTP failures fall back gracefully to empty tables and identity mapping."""
    config = NetworkAnalysisConfig(max_retries=2, retry_backoff=0.01, polite_delay=0.0)

    class FailingClient(HttpClient):
        def get(self, url: str, params: dict[str, Any] | None = None) -> requests.Response | None:
            return mock_response(503, text="Service Unavailable")

        def post(self, url: str, **kwargs: Any) -> requests.Response | None:
            return mock_response(500, text="Internal Server Error")

    engine = NetworkContextEngine(config=config, http_client=FailingClient(config))
    df_genes = make_synthetic_gene_mechanism_df()
    results = engine.run_analysis(df_genes)

    # Ensure all tables are produced without throwing exceptions
    assert len(results["seed_genes_df"]) == 4
    # Identity mapping fallback used
    assert list(results["seed_genes_df"]["HGNC_SYMBOL"]) == ["BUB1B", "MAD1L1", "TRIP13", "GENE_ISOLATED"]
    assert len(results["string_edges_df"]) == 0
    assert len(results["pathway_enrichment_df"]) == 0
    assert len(results["network_gene_context_df"]) == 4


def test_write_network_outputs_to_dir() -> None:
    """Verify file persistence of all Phase 2 result artifacts."""
    config = NetworkAnalysisConfig()
    engine = NetworkContextEngine(config=config, http_client=MockNetworkHttpClient(config))
    results = engine.run_analysis(make_synthetic_gene_mechanism_df())

    with tempfile.TemporaryDirectory() as tmp_dir:
        out_dir = Path(tmp_dir) / "phase2_results"
        paths = write_network_outputs(results, out_dir)

        assert (out_dir / "seed_genes.csv").is_file()
        assert (out_dir / "string_edges.csv").is_file()
        assert (out_dir / "network_summary.csv").is_file()
        assert (out_dir / "pathway_enrichment.csv").is_file()
        assert (out_dir / "network_gene_context.csv").is_file()
        assert (out_dir / "phase2_network_provenance.json").is_file()

        # Check provenance JSON content
        with (out_dir / "phase2_network_provenance.json").open("r", encoding="utf-8") as f:
            prov = json.load(f)
            assert "workflow_phase" in prov
            assert prov["counts"]["mapped_hgnc_count"] == 4
            assert any("directionality" in c for c in prov["interpretation_caveats"])

        # Also test flat_prefix=True
        flat_dir = Path(tmp_dir) / "flat_results"
        flat_paths = write_network_outputs(results, flat_dir, flat_prefix=True)
        assert (flat_dir / "03_phase2_seed_genes.csv").is_file()
        assert (flat_dir / "04_phase2_string_edges.csv").is_file()
        assert (flat_dir / "05_phase2_network_summary.csv").is_file()
        assert (flat_dir / "06_phase2_pathway_enrichment.csv").is_file()
        assert (flat_dir / "07_phase2_network_gene_context.csv").is_file()
        assert (flat_dir / "08_phase2_network_provenance.json").is_file()


def test_cli_runner_with_mocks() -> None:
    """Verify that scripts/run_network_analysis.py CLI entrypoint executes cleanly."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        gene_mech_path = tmp_path / "gene_mechanism_matrix.csv"
        make_synthetic_gene_mechanism_df().to_csv(gene_mech_path, index=False)
        out_dir = tmp_path / "cli_network_output"

        cli_args = [
            "run_network_analysis.py",
            "--gene-mechanism",
            str(gene_mech_path),
            "--output-dir",
            str(out_dir),
            "--string-score",
            "700",
        ]

        with patch("sys.argv", cli_args):
            with patch("mva_hackathon.network_analysis.HttpClient", MockNetworkHttpClient):
                run_network_analysis_main()

        assert (out_dir / "network_gene_context.csv").is_file()
        assert (out_dir / "phase2_network_provenance.json").is_file()


def test_seed_filter_min_passing_variants() -> None:
    """Genes without passing Phase 0 candidates are excluded from the seed set."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        df = make_synthetic_gene_mechanism_df()  # GENE_ISOLATED has 0 passing variants
        with patch("mva_hackathon.network_analysis.HttpClient", MockNetworkHttpClient):
            # Engine must be constructed INSIDE the patch: __init__ binds the client.
            engine = NetworkContextEngine(
                config=NetworkAnalysisConfig(min_passing_variants=1, polite_delay=0.0)
            )
            results = engine.run_analysis(df, None)
        seed_symbols = set(results["seed_genes_df"]["HGNC_SYMBOL"])
        assert len(seed_symbols) == 3
        assert "GENE_ISOLATED" not in seed_symbols


def test_string_batching_preserves_edges() -> None:
    """STRING GET requests are chunked without losing edges across batch merges.

    The mock returns 404 for a batch in which no identifier resolves (real
    STRING behavior), so this also exercises the skip-failed-batch path.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        df = make_synthetic_gene_mechanism_df()
        with patch("mva_hackathon.network_analysis.HttpClient", MockNetworkHttpClient):
            engine = NetworkContextEngine(
                config=NetworkAnalysisConfig(string_batch_size=1, polite_delay=0.0)
            )
            results = engine.run_analysis(df, None)
        edges = results["string_edges_df"]
        assert len(edges) == 2
        pairs = {
            tuple(sorted([e.PREFERRED_NAME_A, e.PREFERRED_NAME_B]))
            for e in edges.itertuples()
        }
        assert ("BUB1B", "TRIP13") in pairs
        assert ("BUB1B", "MAD1L1") in pairs


def test_gprofiler_never_attributes_submitted_genes_without_an_intersection_vector() -> None:
    """An empty or unusable intersection vector must not become all-input membership."""
    engine = NetworkContextEngine()
    captured: dict[str, Any] = {}

    def respond(*args: Any, **kwargs: Any) -> Any:
        response = MagicMock()
        response.status_code = 200
        response.json.return_value = {"result": captured["terms"]}
        return response

    engine.client.post = respond  # type: ignore[method-assign]

    captured["terms"] = [
        {"native": "R-HSA-69620", "name": "Cell Cycle Checkpoints", "source": "REAC", "p_value": 0.001, "term_size": 100, "intersection_size": 0, "intersections": [False, False, False]},
        {"native": "GO:0007059", "name": "chromosome segregation", "source": "GO:BP", "p_value": 0.002, "term_size": 200, "intersection_size": 2, "intersections": [True, False, True]},
        {"native": "hsa04110", "name": "Cell cycle", "source": "KEGG", "p_value": 0.003, "term_size": 124, "intersection_size": 1, "intersections": "malformed"},
    ]
    frame, status = engine.query_gprofiler(["BUB1B", "TRIP13", "MAD1L1"])

    assert status == "success"
    by_id = {row["PATHWAY_ID"]: row for _, row in frame.iterrows()}
    assert by_id["R-HSA-69620"]["HIT_GENES"] == ""
    assert by_id["R-HSA-69620"]["MEMBERSHIP_STATUS"] == "UNAVAILABLE"
    assert by_id["GO:0007059"]["HIT_GENES"] == "BUB1B;MAD1L1"
    assert by_id["GO:0007059"]["MEMBERSHIP_STATUS"] == "OBSERVED_INTERSECTION"
    assert by_id["hsa04110"]["MEMBERSHIP_STATUS"] == "UNAVAILABLE"
    assert "TRIP13" not in by_id["hsa04110"]["HIT_GENES"]
