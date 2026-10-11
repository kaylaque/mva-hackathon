"""Phase 1: Mechanism & Directionality Matrix Engine.

Constructs gene- and variant-level biological mechanism profiles:
- Functional consequence classification (Loss-of-Function, Missense/Hypomorphic, Splicing disruption, Regulatory)
- Cellular process and spindle assembly checkpoint (SAC) mapping
- Pathway directionality (impaired SAC arrest, precocious mitotic exit, chromosome missegregation, aneuploidy)
- Actionable target directionality (e.g. SAC compensation vs downstream aneuploidy stress)
- Rigorous gating and evidence provenance tracking
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


# Biological and pathway knowledge bases for MVA and mitotic / SAC pathways
PATHWAY_KNOWLEDGE_BASE: dict[str, dict[str, str]] = {
    "BUB1B": {
        "gene_name": "BUB1 Mitotic Checkpoint Serine/Threonine Kinase B (BUBR1)",
        "biological_role": "Core component of the Mitotic Checkpoint Complex (MCC); monitors kinetochore-microtubule attachment and delays anaphase onset via CDC20 inhibition.",
        "pathway": "Spindle Assembly Checkpoint (SAC) / Mitotic Cell Cycle",
        "mva_phenotype_link": "Canonical MVA Type 1 (mosaic aneuploidy, microcephaly, growth deficiency, embryonal rhabdomyosarcoma/Wilms tumor risk).",
        "canonical_mechanism": "Loss-of-function / hypomorphic allele combination causes premature APC/C activation and chromosome missegregation.",
        "pathway_directionality": "DOWN (Checkpoint arrest failure -> accelerated mitotic exit -> random chromosome gains/losses).",
        "therapeutic_direction": "Compensate impaired checkpoint / promote chromosome stability or selectively exploit aneuploidy-induced proteotoxic/metabolic stress.",
        "inheritance_requirement": "Biallelic (typically compound heterozygous LoF + missense/splice; constitutive homozygous null is embryonic lethal).",
    },
    "CEP57": {
        "gene_name": "Centrosomal Protein 57",
        "biological_role": "Centrosome architecture and microtubule attachment; anchors kinetochore fibers and regulates spindle pole organization.",
        "pathway": "Centrosome & Spindle Pole Integrity",
        "mva_phenotype_link": "Canonical MVA Type 2 (mosaic aneuploidy, rhizomelic micromelia, short stature, microcephaly, low tumor rate compared to BUB1B).",
        "canonical_mechanism": "Loss of centrosome microtubule anchoring leads to multipolar spindles and aneuploidy.",
        "pathway_directionality": "DOWN (Centrosome structural destabilization -> multipolar spindle formation).",
        "therapeutic_direction": "Stabilize spindle architecture / microtubule dynamics.",
        "inheritance_requirement": "Autosomal recessive (biallelic loss of function).",
    },
    "TRIP13": {
        "gene_name": "Thyroid Hormone Receptor Interactor 13 (AAA+ ATPase)",
        "biological_role": "Converts closed active MAD2 (C-MAD2) to open inactive MAD2 (O-MAD2), dismantling the MCC to allow mitotic exit when kinetochores attach properly.",
        "pathway": "SAC Silencing and Mitotic Checkpoint Complex Turnover",
        "mva_phenotype_link": "Canonical MVA Type 3 (mosaic aneuploidy, Wilms tumor predisposition, microcephaly).",
        "canonical_mechanism": "Loss of TRIP13 impairs SAC silencing and MCC remodeling, disrupting mitotic timing and chromosome fidelity.",
        "pathway_directionality": "DYSREGULATED (Loss of MCC remodeling kinetics -> aberrant mitotic arrest and missegregation).",
        "therapeutic_direction": "Restore SAC turnover kinetics or target checkpoint-dependent mitotic vulnerability.",
        "inheritance_requirement": "Autosomal recessive (biallelic loss of function).",
    },
    "MAD1L1": {
        "gene_name": "MAD1 Mitotic Arrest Deficient 1 Like 1",
        "biological_role": "Kinetochore scaffold that recruits MAD2 to unattached kinetochores and catalyzes C-MAD2 formation.",
        "pathway": "SAC Initiation & Kinetochore Signaling",
        "mva_phenotype_link": "Reported systemic mosaic variegated aneuploidy and childhood cancer predisposition.",
        "canonical_mechanism": "Hypomorphic/LoF alleles reduce MAD2 kinetochore recruitment, attenuating the 'wait anaphase' signal.",
        "pathway_directionality": "DOWN (Impaired SAC signal amplification -> premature anaphase).",
        "therapeutic_direction": "Strengthen checkpoint signaling or target aneuploid survival pathways.",
        "inheritance_requirement": "Biallelic or dominant-negative reports; requires segregation confirmation.",
    },
    "CENATAC": {
        "gene_name": "Centromere-Associated Protein CENATAC (CCDC84 Minor Spliceosome Component)",
        "biological_role": "Core component of the minor spliceosome (U12-dependent) required for accurate splicing of cell-cycle and mitotic regulatory transcripts.",
        "pathway": "U12 Minor Spliceosome & Mitotic Transcript Processing",
        "mva_phenotype_link": "Reported mosaic variegated aneuploidy and microcephaly due to mis-splicing of essential mitotic transcripts.",
        "canonical_mechanism": "Defective minor intron excision causes mis-splicing of critical mitotic regulators (e.g. CENPC, AURKB), impairing chromosome segregation and inducing aneuploidy.",
        "pathway_directionality": "DOWN (Minor spliceosome loss-of-function -> mitotic transcript mis-splicing -> aneuploidy).",
        "therapeutic_direction": "Mitigate spliceosomal deficiency, restore minor intron splicing fidelity, or target resulting mitotic vulnerabilities.",
        "inheritance_requirement": "Autosomal recessive (biallelic loss of function).",
    },
    "SLF2": {
        "gene_name": "SMC5-SMC6 Complex Localization Factor 2",
        "biological_role": "Recruits SMC5-SMC6 cohesin-like complex to DNA lesions and replication forks during cell division.",
        "pathway": "Replication Fork Protection & Cohesion",
        "mva_phenotype_link": "Atelis syndrome (mosaic variegated hyperploidy, microcephaly, growth impairment).",
        "canonical_mechanism": "Loss of SMC5-SMC6 recruitment leads to unresolved replication stress and hyperploid mosaicism.",
        "pathway_directionality": "DOWN (Replication checkpoint failure -> incomplete sister chromatid separation).",
        "therapeutic_direction": "Mitigate replication stress / ATR pathway modulation.",
        "inheritance_requirement": "Autosomal recessive.",
    },
    "SMC5": {
        "gene_name": "Structural Maintenance of Chromosomes 5",
        "biological_role": "Core component of the SMC5-SMC6 complex mediating sister chromatid resolution and chromosome stability.",
        "pathway": "Replication Fork Protection & Cohesion",
        "mva_phenotype_link": "Atelis syndrome with overlapping mosaic variegated hyperploidy.",
        "canonical_mechanism": "Failure to resolve topological entanglements during mitosis leads to chromosome bridges and aneuploidy.",
        "pathway_directionality": "DOWN (Defective chromatid resolution -> chromosome fragmentation and aneuploidy).",
        "therapeutic_direction": "Modulate cohesion turnover and DNA repair machinery.",
        "inheritance_requirement": "Autosomal recessive.",
    },
}

DEFAULT_MECHANISM = {
    "gene_name": "Genome-wide Candidate",
    "biological_role": "Secondary/co-occurring variant in genome-wide screening.",
    "pathway": "General Cellular / Genomic Maintenance",
    "mva_phenotype_link": "Non-canonical or potential dual-diagnosis / secondary clinical finding.",
    "canonical_mechanism": "Variant effect requires manual correlation with patient phenotype.",
    "pathway_directionality": "UNCERTAIN / VARIANT_SPECIFIC",
    "therapeutic_direction": "Assess per variant and symptom profile.",
    "inheritance_requirement": "Assess based on reported OMIM / ClinVar gene inheritance mode.",
}


def infer_variant_consequence_type(row: pd.Series) -> str:
    """Infer functional alteration class from variant annotation fields."""
    consequence = str(row.get("CONSEQUENCE", "")).lower()
    impact = str(row.get("IMPACT", "")).upper()
    splice_ds = float(row.get("SPLICEAI_MAX_DS", 0) or 0)

    # 1. Splice disruption
    if splice_ds >= 0.5 or "splice_donor" in consequence or "splice_acceptor" in consequence:
        return "HIGH_CONFIDENCE_SPLICE_ALTERATION"
    if splice_ds >= 0.2:
        return "POTENTIAL_SPLICE_ALTERATION"

    # 2. Loss of Function (Protein Truncation)
    if any(k in consequence for k in ["stop_gained", "frameshift_variant", "start_lost", "stop_lost"]):
        return "PREDICTED_LOSS_OF_FUNCTION_NULL"

    # 3. Missense / In-frame InDel
    if "missense_variant" in consequence or "inframe" in consequence:
        revel = float(row.get("REVEL_SCORE", 0) or 0)
        cadd = float(row.get("CADD_PHRED", 0) or 0)
        am_score = float(row.get("AM_SCORE", 0) or 0)
        am_class = str(row.get("AM_CLASS", "")).lower()
        if revel >= 0.7 or am_score >= 0.75 or "pathogenic" in am_class:
            return "DELETERIOUS_MISSENSE_HYPOMORPHIC"
        if revel >= 0.5 or am_score >= 0.5 or cadd >= 20.0:
            return "PREDICTED_DAMAGING_MISSENSE"
        return "MISSENSE_UNCERTAIN_SIGNIFICANCE"

    # 4. Regulatory / Non-coding
    if impact == "MODIFIER" or any(k in consequence for k in ["intron_variant", "regulatory_region", "3_prime_utr", "5_prime_utr"]):
        return "REGULATORY_OR_NON_CODING"

    return "OTHER_CONSEQUENCE"


def infer_direction_of_effect(variant_type: str, gene_info: dict[str, str]) -> str:
    """Infer cellular directionality (e.g. checkpoint suppression vs activation)."""
    if "LOSS_OF_FUNCTION" in variant_type or "SPLICE" in variant_type:
        return "LOSS_OF_FUNCTION (Checkpoint/Pathway Down-Regulation)"
    if "DELETERIOUS_MISSENSE" in variant_type or "HYPOMORPHIC" in variant_type:
        return "PARTIAL_LOSS_OF_FUNCTION (Hypomorphic / Sub-threshold Activity)"
    if "DAMAGING_MISSENSE" in variant_type:
        return "SUSPECTED_LOSS_OF_FUNCTION (Reduced Activity or Protein Destabilization)"
    return "UNKNOWN_OR_UNRESOLVED_DIRECTIONALITY"


class MechanismMatrixBuilder:
    """Builds gene-level and variant-level Mechanism and Directionality Matrices."""

    def __init__(self, knowledge_base: dict[str, dict[str, str]] | None = None) -> None:
        self.kb = knowledge_base or PATHWAY_KNOWLEDGE_BASE

    def build_matrix(self, consolidated_candidates: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
        """Construct both variant-level and gene-level mechanism and directionality matrices.

        Returns:
            variant_matrix: Detailed variant mechanism, directionality, and actionability matrix
            gene_matrix: Aggregated gene-level biological mechanism and multi-hit summary
            provenance: Execution metadata and gate verification
        """
        df = consolidated_candidates.copy()

        # Build variant-level matrix
        v_rows = []
        for _, row in df.iterrows():
            sym = str(row.get("SYMBOL", "")).strip().upper()
            gene_meta = self.kb.get(sym, DEFAULT_MECHANISM)

            var_type = infer_variant_consequence_type(row)
            dir_effect = infer_direction_of_effect(var_type, gene_meta)

            entry = {
                "VARIANT_KEY": row.get("VARIANT_KEY", ""),
                "SYMBOL": sym,
                "CHROM": row.get("CHROM", ""),
                "POS": row.get("POS", ""),
                "REF": row.get("REF", ""),
                "ALT": row.get("ALT", ""),
                "HGVS_C": row.get("HGVS_C", "-"),
                "HGVS_P": row.get("HGVS_P", "-"),
                "CONSEQUENCE": row.get("CONSEQUENCE", "-"),
                "IMPACT": row.get("IMPACT", "-"),
                "VARIANT_CONSEQUENCE_CLASS": var_type,
                "PREDICTED_MOLECULAR_DIRECTION": dir_effect,
                "PATHWAY": gene_meta["pathway"],
                "BIOLOGICAL_ROLE": gene_meta["biological_role"],
                "MVA_MECHANISM_LINK": gene_meta["mva_phenotype_link"],
                "CANONICAL_PATHWAY_DIRECTION": gene_meta["pathway_directionality"],
                "THERAPEUTIC_DIRECTIONALITY": gene_meta["therapeutic_direction"],
                "INHERITANCE_REQUIREMENT": gene_meta["inheritance_requirement"],
                "GATE_OVERALL_PASS": bool(row.get("GATE_OVERALL_PASS", False)),
                "GATE_RARITY": bool(row.get("GATE_RARITY", False)),
                "GATE_IMPACT": bool(row.get("GATE_IMPACT", False)),
                "CONSOLIDATED_AF": row.get("CONSOLIDATED_AF", np.nan),
                "RARITY_TIER": row.get("RARITY_TIER", "UNKNOWN"),
                "HOMOZYGOTES_IN_POPULATION": bool(row.get("HOMOZYGOTES_IN_POPULATION", False)),
                "CADD_PHRED": row.get("CADD_PHRED", np.nan),
                "SPLICEAI_MAX_DS": row.get("SPLICEAI_MAX_DS", np.nan),
                "REVEL_SCORE": row.get("REVEL_SCORE", np.nan),
                "AM_SCORE": row.get("AM_SCORE", np.nan),
                "AM_CLASS": row.get("AM_CLASS", "-"),
                "CLIN_SIG": row.get("CLIN_SIG", "-"),
                "ACMG_CLASSIFICATION": row.get("ACMG_CLASSIFICATION", "-"),
                "GENE_QUALIFYING_HITS": int(row.get("GENE_QUALIFYING_HITS", 0)),
                "MULTI_HIT_CANDIDATE": bool(row.get("MULTI_HIT_CANDIDATE", False)),
                "PHASE_STATUS": row.get("PHASE_STATUS", "SINGLE_HIT"),
                "CONSOLIDATION_PRIORITY_SCORE": int(row.get("CONSOLIDATION_PRIORITY_SCORE", 0)),
            }
            v_rows.append(entry)

        variant_matrix = pd.DataFrame(v_rows)
        variant_matrix = variant_matrix.sort_values(
            ["CONSOLIDATION_PRIORITY_SCORE", "GATE_OVERALL_PASS", "SYMBOL"],
            ascending=[False, False, True]
        ).reset_index(drop=True)

        # Build gene-level matrix
        gene_groups = variant_matrix.groupby("SYMBOL")
        g_rows = []
        for symbol, g_df in gene_groups:
            gene_meta = self.kb.get(symbol, DEFAULT_MECHANISM)
            passing_vars = g_df[g_df["GATE_OVERALL_PASS"]]
            qualifying_hits = len(passing_vars)

            all_consequences = "; ".join(sorted(g_df["VARIANT_CONSEQUENCE_CLASS"].unique()))
            passing_consequences = "; ".join(sorted(passing_vars["VARIANT_CONSEQUENCE_CLASS"].unique())) if qualifying_hits > 0 else "NONE_PASSING"

            has_null = passing_vars["VARIANT_CONSEQUENCE_CLASS"].eq("PREDICTED_LOSS_OF_FUNCTION_NULL").any()
            has_missense = passing_vars["VARIANT_CONSEQUENCE_CLASS"].str.contains("MISSENSE").any()
            has_splice = passing_vars["VARIANT_CONSEQUENCE_CLASS"].str.contains("SPLICE").any()

            # Compound heterozygous plausibility check
            compound_het_archetype = False
            archetype_description = "N/A"
            if qualifying_hits >= 2:
                if has_null and (has_missense or has_splice):
                    compound_het_archetype = True
                    archetype_description = "CLASSIC_MVA_ARCHETYPE (Null/Trun + Missense/Splice)"
                elif has_null and qualifying_hits >= 2:
                    compound_het_archetype = True
                    archetype_description = "BIALLELIC_NULL_WARNING (High lethality risk if homozygous)"
                elif has_missense:
                    compound_het_archetype = True
                    archetype_description = "BIALLELIC_HYPOMORPHIC (Compound missense candidates)"

            entry = {
                "SYMBOL": symbol,
                "GENE_NAME": gene_meta["gene_name"],
                "PATHWAY": gene_meta["pathway"],
                "BIOLOGICAL_ROLE": gene_meta["biological_role"],
                "MVA_MECHANISM_LINK": gene_meta["mva_phenotype_link"],
                "CANONICAL_PATHWAY_DIRECTION": gene_meta["pathway_directionality"],
                "THERAPEUTIC_DIRECTIONALITY": gene_meta["therapeutic_direction"],
                "INHERITANCE_REQUIREMENT": gene_meta["inheritance_requirement"],
                "TOTAL_VARIANTS_EVALUATED": len(g_df),
                "TOTAL_PASSING_CANDIDATE_VARIANTS": qualifying_hits,
                "PASSING_CONSEQUENCE_SPECTRUM": passing_consequences,
                "ALL_CONSEQUENCE_SPECTRUM": all_consequences,
                "COMPOUND_HET_ARCHETYPE_FIT": compound_het_archetype,
                "ARCHETYPE_DESCRIPTION": archetype_description,
                "PHASE_GATE_STATUS": "REQUIRES_PHASE_VERIFICATION" if qualifying_hits >= 2 else "SINGLE_HIT_OR_UNPASSING",
                "MAX_CONSOLIDATION_PRIORITY_SCORE": int(g_df["CONSOLIDATION_PRIORITY_SCORE"].max()),
                "TOP_VARIANT_KEYS": "; ".join(g_df.loc[g_df["GATE_OVERALL_PASS"], "VARIANT_KEY"].tolist()),
            }
            g_rows.append(entry)

        gene_matrix = pd.DataFrame(g_rows)
        gene_matrix = gene_matrix.sort_values(
            ["TOTAL_PASSING_CANDIDATE_VARIANTS", "COMPOUND_HET_ARCHETYPE_FIT", "MAX_CONSOLIDATION_PRIORITY_SCORE"],
            ascending=[False, False, False]
        ).reset_index(drop=True)

        now_utc = datetime.now(timezone.utc).isoformat()
        provenance = {
            "pipeline_stage": "Phase 1 (Mechanism & Directionality Matrix)",
            "timestamp_utc": now_utc,
            "total_variants_evaluated": len(variant_matrix),
            "total_genes_evaluated": len(gene_matrix),
            "canonical_mva_genes_evaluated": int(gene_matrix["SYMBOL"].isin(list(PATHWAY_KNOWLEDGE_BASE.keys())).sum()),
            "multihit_genes_passing_gates": int(gene_matrix[gene_matrix["TOTAL_PASSING_CANDIDATE_VARIANTS"] >= 2]["SYMBOL"].count()),
            "disclaimer": "Research candidate exploration only; non-diagnostic; not medical advice.",
        }

        return variant_matrix, gene_matrix, provenance
