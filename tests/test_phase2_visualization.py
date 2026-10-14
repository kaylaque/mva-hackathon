"""Unit tests for Phase 2 network visualization engine and CLI runner.

All tests operate strictly with synthetic fixtures to guarantee determinism,
data privacy, and absence of external network dependency.
"""

from __future__ import annotations

from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

REPO_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import pandas as pd
import pytest

from mva_hackathon.network_visualization import (
    NetworkVisualizer,
    fruchterman_reingold_layout,
)
from run_network_visualization import main as run_network_visualization_main


def make_synthetic_network_fixtures() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Create synthetic non-patient gene context, edges, and enrichment frames."""
    # Gene context (nodes)
    nodes_df = pd.DataFrame([
        {
            "HGNC_SYMBOL": "BUB1B",
            "INPUT_SYMBOL": "BUB1B",
            "PATHWAY": "Spindle Assembly Checkpoint (SAC) / Mitotic Cell Cycle",
            "MVA_MECHANISM_LINK": "Biallelic LoF causes MVA type 1",
            "CANONICAL_PATHWAY_DIRECTION": "REDUCED_SAC_FUNCTION",
            "THERAPEUTIC_DIRECTIONALITY": "RESTORE_CHECKPOINT_FUNCTION",
            "INHERITANCE_REQUIREMENT": "Biallelic",
            "TOTAL_VARIANTS_EVALUATED": 2,
            "TOTAL_PASSING_CANDIDATE_VARIANTS": 2,
            "COMPOUND_HET_ARCHETYPE_FIT": True,
            "PHASE_GATE_STATUS": "REQUIRES_PHASE_VERIFICATION",
            "NETWORK_DEGREE": 2,
            "CLUSTERING_COEFFICIENT": 1.0,
            "COMPONENT_ID": 1,
            "IS_ISOLATED": False,
            "STRING_INTERACTORS": "MAD1L1;TRIP13",
            "TOP_ENRICHED_PATHWAYS": "Mitotic Checkpoint",
            "NETWORK_ROBUSTNESS_NOTE": "Direct PPI cluster",
        },
        {
            "HGNC_SYMBOL": "MAD1L1",
            "INPUT_SYMBOL": "MAD1L1",
            "PATHWAY": "SAC Initiation & Kinetochore Signaling",
            "MVA_MECHANISM_LINK": "Reported in mosaic variegated aneuploidy",
            "CANONICAL_PATHWAY_DIRECTION": "REDUCED_SAC_FUNCTION",
            "THERAPEUTIC_DIRECTIONALITY": "RESTORE_CHECKPOINT_FUNCTION",
            "INHERITANCE_REQUIREMENT": "Biallelic",
            "TOTAL_VARIANTS_EVALUATED": 1,
            "TOTAL_PASSING_CANDIDATE_VARIANTS": 1,
            "COMPOUND_HET_ARCHETYPE_FIT": False,
            "PHASE_GATE_STATUS": "SINGLE_HIT_OR_UNPASSING",
            "NETWORK_DEGREE": 2,
            "CLUSTERING_COEFFICIENT": 1.0,
            "COMPONENT_ID": 1,
            "IS_ISOLATED": False,
            "STRING_INTERACTORS": "BUB1B;TRIP13",
            "TOP_ENRICHED_PATHWAYS": "Mitotic Checkpoint",
            "NETWORK_ROBUSTNESS_NOTE": "Direct PPI cluster",
        },
        {
            "HGNC_SYMBOL": "TRIP13",
            "INPUT_SYMBOL": "TRIP13",
            "PATHWAY": "SAC Silencing",
            "MVA_MECHANISM_LINK": "Canonical MVA gene",
            "CANONICAL_PATHWAY_DIRECTION": "REDUCED_SAC_FUNCTION",
            "THERAPEUTIC_DIRECTIONALITY": "RESTORE_CHECKPOINT_FUNCTION",
            "INHERITANCE_REQUIREMENT": "Autosomal recessive",
            "TOTAL_VARIANTS_EVALUATED": 1,
            "TOTAL_PASSING_CANDIDATE_VARIANTS": 1,
            "COMPOUND_HET_ARCHETYPE_FIT": False,
            "PHASE_GATE_STATUS": "SINGLE_HIT_OR_UNPASSING",
            "NETWORK_DEGREE": 2,
            "CLUSTERING_COEFFICIENT": 1.0,
            "COMPONENT_ID": 1,
            "IS_ISOLATED": False,
            "STRING_INTERACTORS": "BUB1B;MAD1L1",
            "TOP_ENRICHED_PATHWAYS": "Mitotic Checkpoint",
            "NETWORK_ROBUSTNESS_NOTE": "Direct PPI cluster",
        },
        {
            "HGNC_SYMBOL": "RPL13",
            "INPUT_SYMBOL": "RPL13",
            "PATHWAY": "Translation / Ribosome",
            "MVA_MECHANISM_LINK": "Non-canonical",
            "CANONICAL_PATHWAY_DIRECTION": "UNKNOWN",
            "THERAPEUTIC_DIRECTIONALITY": "UNKNOWN",
            "INHERITANCE_REQUIREMENT": "Uncertain",
            "TOTAL_VARIANTS_EVALUATED": 1,
            "TOTAL_PASSING_CANDIDATE_VARIANTS": 1,
            "COMPOUND_HET_ARCHETYPE_FIT": False,
            "PHASE_GATE_STATUS": "SINGLE_HIT_OR_UNPASSING",
            "NETWORK_DEGREE": 1,
            "CLUSTERING_COEFFICIENT": 0.0,
            "COMPONENT_ID": 2,
            "IS_ISOLATED": False,
            "STRING_INTERACTORS": "RPS16",
            "TOP_ENRICHED_PATHWAYS": "Ribosome",
            "NETWORK_ROBUSTNESS_NOTE": "Direct PPI cluster",
        },
        {
            "HGNC_SYMBOL": "RPS16",
            "INPUT_SYMBOL": "RPS16",
            "PATHWAY": "Translation / Ribosome",
            "MVA_MECHANISM_LINK": "Non-canonical",
            "CANONICAL_PATHWAY_DIRECTION": "UNKNOWN",
            "THERAPEUTIC_DIRECTIONALITY": "UNKNOWN",
            "INHERITANCE_REQUIREMENT": "Uncertain",
            "TOTAL_VARIANTS_EVALUATED": 1,
            "TOTAL_PASSING_CANDIDATE_VARIANTS": 1,
            "COMPOUND_HET_ARCHETYPE_FIT": False,
            "PHASE_GATE_STATUS": "SINGLE_HIT_OR_UNPASSING",
            "NETWORK_DEGREE": 1,
            "CLUSTERING_COEFFICIENT": 0.0,
            "COMPONENT_ID": 2,
            "IS_ISOLATED": False,
            "STRING_INTERACTORS": "RPL13",
            "TOP_ENRICHED_PATHWAYS": "Ribosome",
            "NETWORK_ROBUSTNESS_NOTE": "Direct PPI cluster",
        },
        {
            "HGNC_SYMBOL": "GENE_ISOLATED",
            "INPUT_SYMBOL": "GENE_ISOLATED",
            "PATHWAY": "General Maintenance",
            "MVA_MECHANISM_LINK": "Non-canonical",
            "CANONICAL_PATHWAY_DIRECTION": "UNKNOWN",
            "THERAPEUTIC_DIRECTIONALITY": "UNKNOWN",
            "INHERITANCE_REQUIREMENT": "Uncertain",
            "TOTAL_VARIANTS_EVALUATED": 1,
            "TOTAL_PASSING_CANDIDATE_VARIANTS": 1,
            "COMPOUND_HET_ARCHETYPE_FIT": False,
            "PHASE_GATE_STATUS": "SINGLE_HIT_OR_UNPASSING",
            "NETWORK_DEGREE": 0,
            "CLUSTERING_COEFFICIENT": 0.0,
            "COMPONENT_ID": -1,
            "IS_ISOLATED": True,
            "STRING_INTERACTORS": "",
            "TOP_ENRICHED_PATHWAYS": "None",
            "NETWORK_ROBUSTNESS_NOTE": "Singleton",
        },
    ])

    # STRING edges
    edges_df = pd.DataFrame([
        {
            "PREFERRED_NAME_A": "BUB1B",
            "PREFERRED_NAME_B": "MAD1L1",
            "STRING_ID_A": "ENSP1",
            "STRING_ID_B": "ENSP2",
            "SCORE": 0.985,
            "EXPERIMENTAL_SCORE": 0.95,
            "DATABASE_SCORE": 0.92,
            "TEXTMINING_SCORE": 0.88,
            "CURATED_SUPPORT": "YES",
            "TEXTMINING_ONLY": "NO",
        },
        {
            "PREFERRED_NAME_A": "BUB1B",
            "PREFERRED_NAME_B": "TRIP13",
            "STRING_ID_A": "ENSP1",
            "STRING_ID_B": "ENSP3",
            "SCORE": 0.962,
            "EXPERIMENTAL_SCORE": 0.82,
            "DATABASE_SCORE": 0.90,
            "TEXTMINING_SCORE": 0.85,
            "CURATED_SUPPORT": "YES",
            "TEXTMINING_ONLY": "NO",
        },
        {
            "PREFERRED_NAME_A": "MAD1L1",
            "PREFERRED_NAME_B": "TRIP13",
            "STRING_ID_A": "ENSP2",
            "STRING_ID_B": "ENSP3",
            "SCORE": 0.910,
            "EXPERIMENTAL_SCORE": 0.75,
            "DATABASE_SCORE": 0.80,
            "TEXTMINING_SCORE": 0.80,
            "CURATED_SUPPORT": "YES",
            "TEXTMINING_ONLY": "NO",
        },
        {
            "PREFERRED_NAME_A": "RPL13",
            "PREFERRED_NAME_B": "RPS16",
            "STRING_ID_A": "ENSP4",
            "STRING_ID_B": "ENSP5",
            "SCORE": 0.999,
            "EXPERIMENTAL_SCORE": 0.99,
            "DATABASE_SCORE": 0.99,
            "TEXTMINING_SCORE": 0.99,
            "CURATED_SUPPORT": "YES",
            "TEXTMINING_ONLY": "NO",
        },
    ])

    # Pathway enrichment
    enrichment_df = pd.DataFrame([
        {
            "SOURCE_DATABASE": "Reactome",
            "PATHWAY_ID": "R-HSA-141444",
            "PATHWAY_NAME": "Mitotic Spindle Checkpoint",
            "ENTITIES_FOUND": 3,
            "ENTITIES_TOTAL": 28,
            "P_VALUE": 1.2e-5,
            "FDR": 1.5e-4,
            "PATHWAY_SPECIES": "Homo sapiens",
            "HIT_GENES": "BUB1B;MAD1L1;TRIP13",
        },
        {
            "SOURCE_DATABASE": "g:Profiler (GO:BP)",
            "PATHWAY_ID": "GO:0007094",
            "PATHWAY_NAME": "mitotic spindle assembly checkpoint signaling",
            "ENTITIES_FOUND": 3,
            "ENTITIES_TOTAL": 45,
            "P_VALUE": 5.7e-6,
            "FDR": 2.1e-5,
            "PATHWAY_SPECIES": "hsapiens",
            "HIT_GENES": "BUB1B;MAD1L1;TRIP13",
        },
        {
            "SOURCE_DATABASE": "g:Profiler (KEGG)",
            "PATHWAY_ID": "KEGG:03010",
            "PATHWAY_NAME": "Ribosome",
            "ENTITIES_FOUND": 2,
            "ENTITIES_TOTAL": 130,
            "P_VALUE": 0.002,
            "FDR": 0.015,
            "PATHWAY_SPECIES": "hsapiens",
            "HIT_GENES": "RPL13;RPS16",
        },
    ])

    return nodes_df, edges_df, enrichment_df


def test_deterministic_spring_layout() -> None:
    """Verify that layout coordinates are 100% reproducible with fixed random seed."""
    nodes = ["A", "B", "C", "D"]
    edges = [("A", "B"), ("B", "C"), ("C", "D")]

    pos1 = fruchterman_reingold_layout(nodes, edges, seed=42)
    pos2 = fruchterman_reingold_layout(nodes, edges, seed=42)

    for n in nodes:
        assert pos1[n] == pytest.approx(pos2[n], abs=1e-9)


def test_build_network_node_table_and_components() -> None:
    """Verify network node calculation, degree assignment, and component ranking."""
    nodes_df, edges_df, _ = make_synthetic_network_fixtures()
    visualizer = NetworkVisualizer(seed=42)

    network_nodes_df, comp_summary_df, components = visualizer.build_network_node_table(
        nodes_df, edges_df
    )

    # Component 1 (size 3: BUB1B, MAD1L1, TRIP13), component 2 (size 2: RPL13, RPS16),
    # and a singleton component for the node that carries context but no edge.
    assert len(components) == 3
    assert len(components[0]) == 3
    assert len(components[1]) == 2
    assert components[2] == {"GENE_ISOLATED"}

    # Verify component summary table
    assert len(comp_summary_df) == 3
    c1 = comp_summary_df.iloc[0]
    assert c1["COMPONENT_RANK"] == 1
    assert c1["COMPONENT_SIZE"] == 3
    assert c1["EDGE_COUNT"] == 3
    assert bool(c1["HAS_MVA_KNOWN_GENE"]) is True
    assert "BUB1B" in c1["MVA_GENES"]

    c2 = comp_summary_df.iloc[1]
    assert c2["COMPONENT_RANK"] == 2
    assert c2["COMPONENT_SIZE"] == 2
    assert c2["EDGE_COUNT"] == 1
    assert bool(c2["HAS_MVA_KNOWN_GENE"]) is False

    c3 = comp_summary_df.iloc[2]
    assert c3["COMPONENT_RANK"] == 3
    assert c3["COMPONENT_SIZE"] == 1
    assert c3["EDGE_COUNT"] == 0

    # Verify network nodes table contains all 6 input genes
    assert len(network_nodes_df) == 6
    bub1b_row = network_nodes_df[network_nodes_df["HGNC_SYMBOL"] == "BUB1B"].iloc[0]
    assert bub1b_row["NETWORK_DEGREE"] == 2
    assert bub1b_row["COMPONENT_RANK"] == 1
    assert bool(bub1b_row["MVA_CANONICAL"]) is True
    assert bool(bub1b_row["IS_ISOLATED"]) is False

    iso_row = network_nodes_df[network_nodes_df["HGNC_SYMBOL"] == "GENE_ISOLATED"].iloc[0]
    assert iso_row["NETWORK_DEGREE"] == 0
    # An isolated node now carries its own component rank; the rank is not a node rank.
    assert iso_row["COMPONENT_RANK"] == 3
    assert bool(iso_row["IS_ISOLATED"]) is True


def test_visualization_pipeline_file_generation() -> None:
    """Verify all 5 visualization output files are generated correctly with flat names."""
    nodes_df, edges_df, enrichment_df = make_synthetic_network_fixtures()

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        p_nodes_in = tmp_path / "07_phase2_network_gene_context.csv"
        p_edges_in = tmp_path / "04_phase2_string_edges.csv"
        p_enrich_in = tmp_path / "06_phase2_pathway_enrichment.csv"

        nodes_df.to_csv(p_nodes_in, index=False)
        edges_df.to_csv(p_edges_in, index=False)
        enrichment_df.to_csv(p_enrich_in, index=False)

        out_dir = tmp_path / "viz_out"
        visualizer = NetworkVisualizer(seed=42)
        paths = visualizer.run_visualization_pipeline(
            nodes_context_path=p_nodes_in,
            edges_path=p_edges_in,
            enrichment_path=p_enrich_in,
            output_dir=out_dir,
            flat_prefix=True,
        )

        expected_files = [
            "09_phase2_network_nodes.csv",
            "10_phase2_component_summary.csv",
            "11_phase2_largest_component.html",
            "12_phase2_mva_focus_network.html",
            "13_phase2_pathway_enrichment.html",
        ]

        for fname in expected_files:
            file_path = out_dir / fname
            assert file_path.is_file(), f"Expected file {fname} not created!"
            # Ensure HTML files are non-empty and well under 8 MB
            size = file_path.stat().st_size
            assert size > 50, f"File {fname} is suspiciously small ({size} bytes)"
            assert size < 8 * 1024 * 1024, f"File {fname} exceeds 8MB ({size} bytes)!"


def test_cli_runner_visualization() -> None:
    """Verify CLI entrypoint scripts/run_network_visualization.py executes cleanly."""
    nodes_df, edges_df, enrichment_df = make_synthetic_network_fixtures()

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        p_nodes_in = tmp_path / "nodes.csv"
        p_edges_in = tmp_path / "edges.csv"
        p_enrich_in = tmp_path / "enrichment.csv"
        out_dir = tmp_path / "cli_out"

        nodes_df.to_csv(p_nodes_in, index=False)
        edges_df.to_csv(p_edges_in, index=False)
        enrichment_df.to_csv(p_enrich_in, index=False)

        cli_args = [
            "run_network_visualization.py",
            "--gene-context",
            str(p_nodes_in),
            "--edges",
            str(p_edges_in),
            "--enrichment",
            str(p_enrich_in),
            "--output-dir",
            str(out_dir),
            "--flat-names",
        ]

        with patch("sys.argv", cli_args):
            run_network_visualization_main()

        assert (out_dir / "09_phase2_network_nodes.csv").is_file()
        assert (out_dir / "10_phase2_component_summary.csv").is_file()
        assert (out_dir / "11_phase2_largest_component.html").is_file()
        assert (out_dir / "12_phase2_mva_focus_network.html").is_file()
        assert (out_dir / "13_phase2_pathway_enrichment.html").is_file()
