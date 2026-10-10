"""Unit tests for Phase 0 and Phase 1 candidate consolidation and mechanism matrix."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile

REPO_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import pandas as pd
import pytest

from mva_hackathon.candidate_consolidation import (
    CandidateConsolidator,
    ConsolidationGates,
    classify_mva_gene_category,
    normalize_chromosome,
)
from mva_hackathon.mechanism_matrix import MechanismMatrixBuilder
from build_mechanism_matrix import main as build_mechanism_matrix_main


def make_synthetic_evidence_df() -> pd.DataFrame:
    """Create a minimal synthetic non-patient evidence table."""
    return pd.DataFrame([
        {
            "CHROM": "chr15",
            "POS": 40263435,
            "REF": "A",
            "ALT": "C",
            "SYMBOL": "BUB1B",
            "CONSEQUENCE": "stop_gained",
            "IMPACT": "HIGH",
            "HGVS_C": "c.2210T>G",
            "HGVS_P": "p.Leu737*",
            "AF_POPULATION": 0.000098,
            "GNOMAD_AF": 0.000098,
            "GNOMAD_STATUS": "PASS",
            "GNOMAD_NHOMALT": 0,
            "CADD_PHRED": 36.0,
            "CLIN_SIG": "Pathogenic",
            "CANONICAL": "YES",
            "PRIORITY_SCORE": 80,
        },
        {
            "CHROM": "15",
            "POS": 40280000,
            "REF": "T",
            "ALT": "G",
            "SYMBOL": "BUB1B",
            "CONSEQUENCE": "missense_variant",
            "IMPACT": "MODERATE",
            "HGVS_C": "c.3006T>G",
            "HGVS_P": "p.Asn1002Lys",
            "AF_POPULATION": 0.000001,
            "GNOMAD_AF": 0.000001,
            "GNOMAD_STATUS": "PASS",
            "GNOMAD_NHOMALT": 0,
            "CADD_PHRED": 24.5,
            "REVEL_SCORE": 0.472,
            "AM_SCORE": 0.92,
            "AM_CLASS": "likely_pathogenic",
            "CANONICAL": "YES",
            "PRIORITY_SCORE": 70,
        },
        {
            "CHROM": "chr7",
            "POS": 10000,
            "REF": "G",
            "ALT": "A",
            "SYMBOL": "MAD1L1",
            "CONSEQUENCE": "missense_variant",
            "IMPACT": "MODERATE",
            "HGVS_C": "c.175C>T",
            "HGVS_P": "p.Arg59Cys",
            "AF_POPULATION": 0.0006,
            "GNOMAD_AF": 0.0006,
            "GNOMAD_STATUS": "PASS",
            "GNOMAD_NHOMALT": 0,
            "CADD_PHRED": 22.0,
            "REVEL_SCORE": 0.55,
            "AM_SCORE": 0.65,
            "CLIN_SIG": "conflicting",
            "CANONICAL": "YES",
            "PRIORITY_SCORE": 55,
        },
        {
            "CHROM": "1",
            "POS": 20000,
            "REF": "C",
            "ALT": "T",
            "SYMBOL": "GENE_COMMON",
            "CONSEQUENCE": "missense_variant",
            "IMPACT": "MODERATE",
            "AF_POPULATION": 0.05,
            "GNOMAD_AF": 0.05,
            "GNOMAD_STATUS": "PASS",
            "GNOMAD_NHOMALT": 10,
            "CADD_PHRED": 12.0,
            "CANONICAL": "YES",
            "PRIORITY_SCORE": 10,
        },
    ])


def make_synthetic_exomiser_variants_df() -> pd.DataFrame:
    """Create a minimal synthetic Exomiser variants export."""
    return pd.DataFrame([
        {
            "chrom": "15",
            "pos": 40263435,
            "ref": "A",
            "alt": "C",
            "variantScore": 0.95,
            "geneScore": 0.98,
            "phenotypeScore": 0.92,
            "rank": 1,
        },
        {
            "chrom": "15",
            "pos": 40280000,
            "ref": "T",
            "alt": "G",
            "variantScore": 0.88,
            "geneScore": 0.98,
            "phenotypeScore": 0.92,
            "rank": 2,
        },
    ])


def make_synthetic_exomiser_genes_df() -> pd.DataFrame:
    """Create a minimal synthetic Exomiser genes export."""
    return pd.DataFrame([
        {
            "geneSymbol": "BUB1B",
            "geneScore": 0.98,
            "phenotypeScore": 0.92,
            "rank": 1,
        },
        {
            "geneSymbol": "MAD1L1",
            "geneScore": 0.75,
            "phenotypeScore": 0.70,
            "rank": 5,
        },
    ])


def test_chromosome_normalization():
    assert normalize_chromosome("chr1") == "1"
    assert normalize_chromosome("CHR22") == "22"
    assert normalize_chromosome("chrX") == "X"
    assert normalize_chromosome("23") == "X"
    assert normalize_chromosome("24") == "Y"
    assert normalize_chromosome("M") == "MT"


def test_classify_mva_gene_category():
    assert classify_mva_gene_category("BUB1B") == "CANONICAL_MVA"
    assert classify_mva_gene_category("CEP57") == "CANONICAL_MVA"
    assert classify_mva_gene_category("TRIP13") == "CANONICAL_MVA"
    assert classify_mva_gene_category("MAD1L1") == "REPORTED_MVA_CA"
    assert classify_mva_gene_category("SLF2") == "OVERLAPPING_ATELIS"
    assert classify_mva_gene_category("UNKNOWN_GENE") == "GENOME_WIDE"


def test_phase_0_consolidation_and_gating():
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        ev_file = tmp_path / "evidence.csv"
        exo_v_file = tmp_path / "exomiser_var.tsv"
        exo_g_file = tmp_path / "exomiser_gene.tsv"
        out_cand = tmp_path / "consolidated_candidates.csv"

        make_synthetic_evidence_df().to_csv(ev_file, index=False)
        make_synthetic_exomiser_variants_df().to_csv(exo_v_file, sep="\t", index=False)
        make_synthetic_exomiser_genes_df().to_csv(exo_g_file, sep="\t", index=False)

        consolidator = CandidateConsolidator(ConsolidationGates(max_af=0.001))
        df, prov = consolidator.consolidate(
            evidence_path=ev_file,
            exomiser_variants_path=exo_v_file,
            exomiser_genes_path=exo_g_file,
            output_path=out_cand,
        )

        assert out_cand.is_file()
        assert len(df) == 4
        # Common variant should fail rarity gate
        common_var = df[df["SYMBOL"] == "GENE_COMMON"].iloc[0]
        assert not common_var["GATE_RARITY"]
        assert not common_var["GATE_OVERALL_PASS"]

        # BUB1B variants should pass gates and form a multi-hit candidate
        bub1b_vars = df[df["SYMBOL"] == "BUB1B"]
        assert len(bub1b_vars) == 2
        assert bub1b_vars["GATE_OVERALL_PASS"].all()
        assert bub1b_vars["MULTI_HIT_CANDIDATE"].all()
        assert bub1b_vars["GENE_QUALIFYING_HITS"].iloc[0] == 2

        # Check provenance
        assert prov["total_variants_ingested"] == 4
        assert prov["total_variants_passing_gates"] == 3
        assert prov["total_multihit_genes"] == 1


def test_phase_1_mechanism_matrix():
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        ev_file = tmp_path / "evidence.csv"
        make_synthetic_evidence_df().to_csv(ev_file, index=False)

        consolidator = CandidateConsolidator()
        cand_df, _ = consolidator.consolidate(evidence_path=ev_file)

        matrix_builder = MechanismMatrixBuilder()
        v_matrix, g_matrix, prov = matrix_builder.build_matrix(cand_df)

        assert len(v_matrix) == 4
        assert len(g_matrix) == 3

        # Check BUB1B in gene matrix
        bub1b_row = g_matrix[g_matrix["SYMBOL"] == "BUB1B"].iloc[0]
        assert bub1b_row["TOTAL_PASSING_CANDIDATE_VARIANTS"] == 2
        assert bool(bub1b_row["COMPOUND_HET_ARCHETYPE_FIT"]) is True
        assert bub1b_row["PHASE_GATE_STATUS"] == "REQUIRES_PHASE_VERIFICATION"
        assert "Spindle Assembly Checkpoint" in bub1b_row["PATHWAY"]

        # Check variant consequence classes
        bub1b_v1 = v_matrix[(v_matrix["SYMBOL"] == "BUB1B") & (v_matrix["CONSEQUENCE"] == "stop_gained")].iloc[0]
        assert bub1b_v1["VARIANT_CONSEQUENCE_CLASS"] == "PREDICTED_LOSS_OF_FUNCTION_NULL"
        assert "LOSS_OF_FUNCTION" in bub1b_v1["PREDICTED_MOLECULAR_DIRECTION"]


def test_phase_1_cli_requires_phase_0_columns(monkeypatch):
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        raw_evidence = tmp_path / "raw_evidence.csv"
        output_variant = tmp_path / "variant_matrix.csv"
        output_gene = tmp_path / "gene_matrix.csv"
        make_synthetic_evidence_df().to_csv(raw_evidence, index=False)

        monkeypatch.setattr(
            sys,
            "argv",
            [
                "build_mechanism_matrix.py",
                "--candidates", str(raw_evidence),
                "--output-variant-matrix", str(output_variant),
                "--output-gene-matrix", str(output_gene),
            ],
        )

        with pytest.raises(ValueError, match="requires the output of consolidate_candidates.py"):
            build_mechanism_matrix_main()
