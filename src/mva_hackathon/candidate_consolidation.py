"""Phase 0: Candidate Intake, Quality & Evidence Gates, and Candidate Consolidation.

Consumes caller-supplied existing evidence tables and optional Exomiser outputs
via CLI paths without rerunning upstream QC, VEP, snpEff, or Exomiser.
Applies transparent research gates and outputs a consolidated candidate table
with full provenance.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


# Gene sets relevant to MVA and related chromosome-instability / clinical phenotypes
CANONICAL_MVA_GENES: frozenset[str] = frozenset({"BUB1B", "CEP57", "TRIP13"})
REPORTED_MVA_CA_GENES: frozenset[str] = frozenset({"CENATAC", "MAD1L1"})
ATELIS_SYNDROME_GENES: frozenset[str] = frozenset({"SLF2", "SMC5"})
SAC_CORE_GENES: frozenset[str] = frozenset({
    "BUB1", "BUB3", "MAD2L1", "MAD2L2", "CDC20", "AURKA", "AURKB",
    "PLK1", "CENPE", "TTK", "MPS1", "ZW10", "ZWILCH", "KNTC1", "KNL1", "NDC80"
})
CHROMOSOME_COHESIN_GENES: frozenset[str] = frozenset({
    "STAG2", "SMC1A", "SMC3", "RAD21", "NIPBL", "NCAPD2", "NCAPD3",
    "NCAPH", "NCAPG", "SMC2", "SMC4"
})
CILIA_CENTROSOME_GENES: frozenset[str] = frozenset({
    "CEP290", "OFD1", "PKD1", "PKD2", "TMEM216", "TMEM67", "BBS1", "BBS2",
    "BBS4", "BBS7", "BBS9", "CEP152", "PCNT", "PLK4", "B9D1", "B9D2",
    "NPHP1", "NPHP3", "NPHP4", "INPP5E", "CC2D2A", "RPGRIP1L"
})
NEPHROCALCINOSIS_ION_GENES: frozenset[str] = frozenset({
    "CLDN16", "CLDN19", "CASR", "SLC12A1", "KCNJ1", "BSND", "SLC4A1",
    "SLC34A1", "SLC34A3", "ATP6V0A4", "ATP6V1B1", "CYP24A1", "AGXT",
    "GRHPR", "HOGA1", "SLC22A12", "SLC2A9", "MAGED2", "CLCN5", "OCRL",
    "VDR", "AP2S1", "GNA11"
})
RMS_CANCER_GENES: frozenset[str] = frozenset({
    "TP53", "DICER1", "WT1", "NF1", "HRAS", "KRAS", "NRAS", "FGFR4",
    "PIK3CA", "FBXW7", "BCOR", "CTNNB1", "ALK", "MYOD1", "PAX3", "PAX7"
})

ALL_MVA_RELEVANT_GENES: frozenset[str] = (
    CANONICAL_MVA_GENES
    | REPORTED_MVA_CA_GENES
    | ATELIS_SYNDROME_GENES
    | SAC_CORE_GENES
    | CHROMOSOME_COHESIN_GENES
)


@dataclass(frozen=True)
class ConsolidationGates:
    """Configurable thresholds for candidate filtering and research triage."""

    max_af: float = 0.001
    min_cadd_phred: float = 20.0
    min_spliceai_ds: float = 0.2
    min_revel: float = 0.5
    min_alphamissense: float = 0.5
    reference_build: str = "GRCh38"
    pipeline_version: str = "0.2.0"


def normalize_chromosome(val: object) -> str:
    """Standardize chromosome representation (e.g. 'chr1' -> '1', '23' -> 'X')."""
    if pd.isna(val):
        return ""
    chrom = str(val).strip()
    chrom = chrom.removeprefix("chr").removeprefix("CHR")
    mapping = {"23": "X", "24": "Y", "M": "MT"}
    return mapping.get(chrom, chrom).upper()


def build_variant_key(df: pd.DataFrame) -> pd.Series:
    """Construct canonical VARIANT_KEY (CHROM:POS:REF:ALT)."""
    if "VARIANT_KEY" in df.columns and df["VARIANT_KEY"].notna().all():
        return df["VARIANT_KEY"].astype(str)

    chrom = df["CHROM"].map(normalize_chromosome)
    pos = pd.to_numeric(df["POS"], errors="coerce")
    ref = df["REF"].astype(str).str.strip().str.upper()
    alt = df["ALT"].astype(str).str.strip().str.upper()

    valid = chrom.ne("") & pos.notna() & ref.ne("") & alt.ne("") & ref.ne("NAN") & alt.ne("NAN")
    keys = pd.Series(pd.NA, index=df.index, dtype="string")
    keys.loc[valid] = (
        chrom.loc[valid]
        + ":" + pos.loc[valid].astype(int).astype(str)
        + ":" + ref.loc[valid]
        + ":" + alt.loc[valid]
    )
    return keys


def classify_mva_gene_category(symbol: object) -> str:
    """Assign structural biological category to candidate gene."""
    if pd.isna(symbol):
        return "UNKNOWN"
    sym = str(symbol).strip().upper()
    if sym in CANONICAL_MVA_GENES:
        return "CANONICAL_MVA"
    if sym in REPORTED_MVA_CA_GENES:
        return "REPORTED_MVA_CA"
    if sym in ATELIS_SYNDROME_GENES:
        return "OVERLAPPING_ATELIS"
    if sym in SAC_CORE_GENES:
        return "SAC_CORE"
    if sym in CHROMOSOME_COHESIN_GENES:
        return "CHROMOSOME_COHESIN"
    if sym in CILIA_CENTROSOME_GENES:
        return "CILIA_CENTROSOME"
    if sym in NEPHROCALCINOSIS_ION_GENES:
        return "NEPHROCALCINOSIS_ION"
    if sym in RMS_CANCER_GENES:
        return "RMS_CANCER_PREDISPOSITION"
    return "GENOME_WIDE"


class CandidateConsolidator:
    """Executes Phase 0 candidate intake, verification, gating, and consolidation."""

    REQUIRED_GENOMIC_COLUMNS = frozenset({"CHROM", "POS", "REF", "ALT"})

    def __init__(self, gates: ConsolidationGates | None = None) -> None:
        self.gates = gates or ConsolidationGates()

    def load_evidence_table(self, path: Path | str) -> pd.DataFrame:
        """Load and normalize caller-supplied variant evidence table."""
        source = Path(path)
        if not source.is_file():
            raise FileNotFoundError(f"Evidence table not found: {source}")

        sep = "\t" if source.suffix.lower() in {".tsv", ".txt"} else ","
        df = pd.read_csv(source, sep=sep, low_memory=False)
        df.columns = df.columns.str.strip()

        # Check genomic coordinates
        if not self.REQUIRED_GENOMIC_COLUMNS.issubset(df.columns):
            if "VARIANT_KEY" in df.columns:
                # Parse CHROM, POS, REF, ALT from VARIANT_KEY if missing
                split_keys = df["VARIANT_KEY"].astype(str).str.split(":", expand=True)
                if split_keys.shape[1] == 4:
                    df["CHROM"] = split_keys[0]
                    df["POS"] = pd.to_numeric(split_keys[1], errors="coerce")
                    df["REF"] = split_keys[2]
                    df["ALT"] = split_keys[3]
                else:
                    raise ValueError(f"{source} lacks required genomic columns (CHROM, POS, REF, ALT)")
            else:
                raise ValueError(f"{source} lacks required genomic columns: {sorted(self.REQUIRED_GENOMIC_COLUMNS - set(df.columns))}")

        df["CHROM_NORM"] = df["CHROM"].map(normalize_chromosome)
        df["VARIANT_KEY"] = build_variant_key(df)
        if df["VARIANT_KEY"].isna().any():
            raise ValueError(f"{source} contains records without valid genomic coordinates")

        # Resolve representative transcript if multiline annotations exist
        canonical = df.get("CANONICAL", pd.Series("", index=df.index)).astype(str).str.upper().isin({"YES", "Y", "TRUE", "1"})
        priority = pd.to_numeric(df.get("PRIORITY_SCORE", pd.Series(0, index=df.index)), errors="coerce").fillna(0)
        df = df.assign(_canonical=canonical, _priority=priority)
        df = (
            df.sort_values(["VARIANT_KEY", "_canonical", "_priority"], ascending=[True, False, False])
            .drop_duplicates("VARIANT_KEY", keep="first")
            .drop(columns=["_canonical", "_priority"])
            .reset_index(drop=True)
        )
        return df

    def load_exomiser_variants(self, path: Path | str) -> pd.DataFrame:
        """Load and normalize caller-supplied Exomiser variants export."""
        source = Path(path)
        if not source.is_file():
            raise FileNotFoundError(f"Exomiser variants file not found: {source}")

        sep = "\t" if source.suffix.lower() in {".tsv", ".txt"} else ","
        df = pd.read_csv(source, sep=sep, low_memory=False)
        df.columns = df.columns.str.strip()

        aliases = {
            "variantScore": "EXOMISER_VARIANT_SCORE",
            "geneScore": "EXOMISER_GENE_SCORE",
            "phenotypeScore": "EXOMISER_PHENOTYPE_SCORE",
            "rank": "EXOMISER_RANK",
            "chrom": "CHROM",
            "pos": "POS",
            "ref": "REF",
            "alt": "ALT",
        }
        df = df.rename(columns={k: v for k, v in aliases.items() if k in df.columns})

        if not self.REQUIRED_GENOMIC_COLUMNS.issubset(df.columns):
            if "VARIANT_KEY" in df.columns:
                pass
            else:
                raise ValueError(f"{source} lacks genomic columns required to join Exomiser variants")

        df["VARIANT_KEY"] = build_variant_key(df)
        keep = ["VARIANT_KEY"]
        for col in ["EXOMISER_VARIANT_SCORE", "EXOMISER_GENE_SCORE", "EXOMISER_PHENOTYPE_SCORE", "EXOMISER_RANK"]:
            if col in df.columns:
                keep.append(col)
        return df[keep].drop_duplicates("VARIANT_KEY", keep="first")

    def load_exomiser_genes(self, path: Path | str) -> pd.DataFrame:
        """Load and normalize caller-supplied Exomiser genes export."""
        source = Path(path)
        if not source.is_file():
            raise FileNotFoundError(f"Exomiser genes file not found: {source}")

        sep = "\t" if source.suffix.lower() in {".tsv", ".txt"} else ","
        df = pd.read_csv(source, sep=sep, low_memory=False)
        df.columns = df.columns.str.strip()

        aliases = {
            "geneSymbol": "SYMBOL",
            "GENE_SYMBOL": "SYMBOL",
            "gene": "SYMBOL",
            "score": "EXOMISER_GENE_SCORE",
            "geneScore": "EXOMISER_GENE_SCORE",
            "phenotypeScore": "EXOMISER_PHENOTYPE_SCORE",
            "rank": "EXOMISER_GENE_RANK",
        }
        df = df.rename(columns={k: v for k, v in aliases.items() if k in df.columns})
        if "SYMBOL" not in df.columns:
            raise ValueError(f"{source} lacks gene symbol column for gene-level merge")

        df["SYMBOL"] = df["SYMBOL"].astype(str).str.strip().str.upper()
        keep = ["SYMBOL"]
        for col in ["EXOMISER_GENE_SCORE", "EXOMISER_PHENOTYPE_SCORE", "EXOMISER_GENE_RANK"]:
            if col in df.columns:
                keep.append(col)
        return df[keep].drop_duplicates("SYMBOL", keep="first")

    def apply_gates(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply transparent research quality, population rarity, and impact gates."""
        res = df.copy()

        # 1. Coordinate Integrity Gate
        pos_numeric = pd.to_numeric(res["POS"], errors="coerce")
        coord_valid = (
            res["CHROM_NORM"].ne("")
            & pos_numeric.notna()
            & pos_numeric.gt(0)
            & res["REF"].astype(str).ne("")
            & res["ALT"].astype(str).ne("")
        )
        res["GATE_COORDINATE"] = coord_valid

        # 2. Population Rarity Gate
        gnomad_af = pd.to_numeric(res.get("GNOMAD_AF", pd.Series(np.nan, index=res.index)), errors="coerce")
        pop_af = pd.to_numeric(res.get("AF_POPULATION", pd.Series(np.nan, index=res.index)), errors="coerce")
        best_af = gnomad_af.combine_first(pop_af)
        res["CONSOLIDATED_AF"] = best_af

        gnomad_status = res.get("GNOMAD_STATUS", pd.Series(pd.NA, index=res.index)).astype(str).str.lower()
        is_absent = gnomad_status.eq("not_found")
        is_rare_af = best_af.le(self.gates.max_af)
        is_unknown_af = best_af.isna() & ~is_absent

        rarity_pass = is_absent | is_rare_af
        res["GATE_RARITY"] = rarity_pass
        res["RARITY_TIER"] = np.select(
            [is_absent, best_af.le(1e-5), best_af.le(1e-4), best_af.le(1e-3), best_af.le(1e-2), best_af.gt(1e-2)],
            ["ABSENT", "ULTRA_RARE (<=1e-5)", "VERY_RARE (<=1e-4)", "RARE (<=1e-3)", "LOW_FREQ (<=0.01)", "COMMON (>0.01)"],
            default="UNKNOWN"
        )

        # Flag homozygotes in population
        nhomalt = pd.to_numeric(res.get("GNOMAD_NHOMALT", pd.Series(0, index=res.index)), errors="coerce").fillna(0)
        res["HOMOZYGOTES_IN_POPULATION"] = nhomalt.gt(0)

        # 3. Functional Impact & Deleteriousness Gate
        impact = res.get("IMPACT", pd.Series("", index=res.index)).astype(str).str.upper()
        consequence = res.get("CONSEQUENCE", pd.Series("", index=res.index)).astype(str).str.lower()
        clin_sig = res.get("CLIN_SIG", pd.Series("", index=res.index)).astype(str).str.lower()

        cadd_phred = pd.to_numeric(res.get("CADD_PHRED", pd.Series(np.nan, index=res.index)), errors="coerce")
        spliceai_max = pd.to_numeric(res.get("SPLICEAI_MAX_DS", pd.Series(np.nan, index=res.index)), errors="coerce")
        revel = pd.to_numeric(res.get("REVEL_SCORE", pd.Series(np.nan, index=res.index)), errors="coerce")
        am_score = pd.to_numeric(res.get("AM_SCORE", pd.Series(np.nan, index=res.index)), errors="coerce")
        am_class = res.get("AM_CLASS", pd.Series("", index=res.index)).astype(str).str.lower()

        is_high_impact = impact.eq("HIGH")
        is_pvs1_type = consequence.str.contains("stop_gained|frameshift|splice_donor|splice_acceptor", regex=True)
        is_clinvar_path = clin_sig.str.contains("pathogenic|likely_pathogenic", regex=True)
        is_cadd_high = cadd_phred.ge(self.gates.min_cadd_phred)
        is_splice_high = spliceai_max.ge(self.gates.min_spliceai_ds)
        is_revel_high = revel.ge(self.gates.min_revel)
        is_am_high = am_score.ge(self.gates.min_alphamissense) | am_class.str.contains("pathogenic", regex=True)

        impact_pass = (
            is_high_impact
            | is_pvs1_type
            | is_clinvar_path
            | is_cadd_high
            | is_splice_high
            | is_revel_high
            | is_am_high
        )
        res["GATE_IMPACT"] = impact_pass

        # 4. MVA Biological Relevance Category
        symbol = res.get("SYMBOL", pd.Series("", index=res.index)).astype(str).str.strip().str.upper()
        res["MVA_BIOLOGICAL_CATEGORY"] = symbol.map(classify_mva_gene_category)
        res["GATE_MVA_CONTEXT"] = res["MVA_BIOLOGICAL_CATEGORY"].ne("GENOME_WIDE")

        # 5. Multi-hit / Candidate Recessive Pattern (Count rare qualifying hits per gene)
        qualifying_hit = coord_valid & rarity_pass & impact_pass
        gene_hit_counts = res.loc[qualifying_hit].groupby("SYMBOL")["VARIANT_KEY"].nunique()
        res["GENE_QUALIFYING_HITS"] = res["SYMBOL"].map(gene_hit_counts).fillna(0).astype(int)
        res["MULTI_HIT_CANDIDATE"] = res["GENE_QUALIFYING_HITS"].ge(2)
        res["PHASE_STATUS"] = np.where(
            res["MULTI_HIT_CANDIDATE"],
            "UNCONFIRMED_PHASE (Trans/Cis verification required)",
            "SINGLE_HIT"
        )

        # 6. Overall Candidate Pass Gate
        # A candidate passes if coordinate is valid, rarity passes (or unknown with high consequence in MVA gene),
        # and impact passes.
        overall_pass = coord_valid & (rarity_pass | (is_unknown_af & res["GATE_MVA_CONTEXT"])) & impact_pass
        res["GATE_OVERALL_PASS"] = overall_pass

        # Reason string detailing passed/failed criteria
        reasons = []
        for _, row in res.iterrows():
            parts = []
            if not row["GATE_COORDINATE"]:
                parts.append("FAIL_COORD")
            if not row["GATE_RARITY"]:
                parts.append("FAIL_RARITY")
            if not row["GATE_IMPACT"]:
                parts.append("FAIL_IMPACT")
            if row["GATE_OVERALL_PASS"]:
                cat = row["MVA_BIOLOGICAL_CATEGORY"]
                if cat != "GENOME_WIDE":
                    parts.append(f"PASS_{cat}")
                else:
                    parts.append("PASS_GENOME_WIDE")
                if row["MULTI_HIT_CANDIDATE"]:
                    parts.append("MULTI_HIT")
            res_str = "; ".join(parts) if parts else "REVIEW"
            reasons.append(res_str)
        res["GATE_SUMMARY"] = reasons

        # 7. Composite Prioritization Score (Transparent integer rank for research triage)
        score = pd.Series(0, index=res.index, dtype=int)
        score += np.where(res["GATE_OVERALL_PASS"], 20, 0)
        score += np.where(res["MVA_BIOLOGICAL_CATEGORY"].eq("CANONICAL_MVA"), 50, 0)
        score += np.where(res["MVA_BIOLOGICAL_CATEGORY"].eq("REPORTED_MVA_CA"), 35, 0)
        score += np.where(res["MVA_BIOLOGICAL_CATEGORY"].eq("SAC_CORE"), 25, 0)
        score += np.where(res["MULTI_HIT_CANDIDATE"], 30, 0)
        score += np.where(is_high_impact | is_pvs1_type, 20, 0)
        score += np.where(is_clinvar_path, 25, 0)
        score += np.where(is_cadd_high, 10, 0)
        score += np.where(is_splice_high, 10, 0)
        score += np.where(is_revel_high | is_am_high, 10, 0)
        score += np.where(is_absent, 15, np.where(best_af.le(1e-4), 10, 0))

        # Add Exomiser score contribution if present
        if "EXOMISER_VARIANT_SCORE" in res.columns:
            exo_var = pd.to_numeric(res["EXOMISER_VARIANT_SCORE"], errors="coerce").fillna(0)
            score += (exo_var * 20).astype(int)
        if "EXOMISER_GENE_SCORE" in res.columns:
            exo_gene = pd.to_numeric(res["EXOMISER_GENE_SCORE"], errors="coerce").fillna(0)
            score += (exo_gene * 20).astype(int)

        res["CONSOLIDATION_PRIORITY_SCORE"] = score
        return res.sort_values(
            ["CONSOLIDATION_PRIORITY_SCORE", "GATE_OVERALL_PASS", "GENE_QUALIFYING_HITS"],
            ascending=[False, False, False]
        ).reset_index(drop=True)

    def consolidate(
        self,
        evidence_path: Path | str,
        exomiser_variants_path: Path | str | None = None,
        exomiser_genes_path: Path | str | None = None,
        output_path: Path | str | None = None,
    ) -> tuple[pd.DataFrame, dict[str, Any]]:
        """Consolidate evidence, merge optional Exomiser tables, and apply gates."""
        evidence_df = self.load_evidence_table(evidence_path)

        # Merge Exomiser variant table if provided
        if exomiser_variants_path:
            exo_var_df = self.load_exomiser_variants(exomiser_variants_path)
            # Avoid duplicate column names if already in evidence_df
            drop_cols = [c for c in exo_var_df.columns if c != "VARIANT_KEY" and c in evidence_df.columns]
            if drop_cols:
                evidence_df = evidence_df.drop(columns=drop_cols)
            evidence_df = evidence_df.merge(exo_var_df, on="VARIANT_KEY", how="left")

        # Merge Exomiser gene table if provided
        if exomiser_genes_path:
            exo_gene_df = self.load_exomiser_genes(exomiser_genes_path)
            # Prefix to avoid collision if variant-level gene scores exist
            rename_map = {
                c: f"GENE_LEVEL_{c}" for c in exo_gene_df.columns
                if c != "SYMBOL" and c in evidence_df.columns
            }
            if rename_map:
                exo_gene_df = exo_gene_df.rename(columns=rename_map)
            evidence_df = evidence_df.merge(exo_gene_df, on="SYMBOL", how="left")

        # Apply gates & ranking
        consolidated = self.apply_gates(evidence_df)

        # Append provenance
        now_utc = datetime.now(timezone.utc).isoformat()
        consolidated["PROVENANCE_EVIDENCE_SOURCE"] = str(Path(evidence_path).name)
        consolidated["PROVENANCE_EXOMISER_VAR_SOURCE"] = str(Path(exomiser_variants_path).name) if exomiser_variants_path else "NONE"
        consolidated["PROVENANCE_EXOMISER_GENE_SOURCE"] = str(Path(exomiser_genes_path).name) if exomiser_genes_path else "NONE"
        consolidated["PROVENANCE_GENOME_BUILD"] = self.gates.reference_build
        consolidated["PROVENANCE_PIPELINE_VERSION"] = self.gates.pipeline_version
        consolidated["PROVENANCE_TIMESTAMP"] = now_utc

        # Summary provenance dictionary
        provenance = {
            "pipeline_stage": "Phase 0 (Candidate Intake, Gating, and Consolidation)",
            "timestamp_utc": now_utc,
            "genome_build": self.gates.reference_build,
            "pipeline_version": self.gates.pipeline_version,
            "gates_config": asdict(self.gates),
            "input_evidence": str(Path(evidence_path).resolve()),
            "input_exomiser_variants": str(Path(exomiser_variants_path).resolve()) if exomiser_variants_path else None,
            "input_exomiser_genes": str(Path(exomiser_genes_path).resolve()) if exomiser_genes_path else None,
            "total_variants_ingested": len(consolidated),
            "total_variants_passing_gates": int(consolidated["GATE_OVERALL_PASS"].sum()),
            "total_mva_core_variants": int(consolidated["MVA_BIOLOGICAL_CATEGORY"].isin(["CANONICAL_MVA", "REPORTED_MVA_CA"]).sum()),
            "total_multihit_genes": int(consolidated.loc[consolidated["MULTI_HIT_CANDIDATE"], "SYMBOL"].nunique()),
            "disclaimer": "Research candidate exploration only; non-diagnostic; not medical advice.",
        }

        if output_path:
            out = Path(output_path)
            out.parent.mkdir(parents=True, exist_ok=True)
            sep = "\t" if out.suffix.lower() in {".tsv", ".txt"} else ","
            consolidated.to_csv(out, sep=sep, index=False)

        return consolidated, provenance
