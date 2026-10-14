"""Phase 4 Hypothesis Classification Engine.

Deterministically classifies drug-target candidate pairs into four hypothesis classes:
- H1: DISEASE_MECHANISM_CORRECTION (signed path/reaction moves patient-deficit endpoint
      toward physiological direction: SAC restoration/stabilization, KT-MT attachment
      stabilization, minor-intron splicing rescue for CENATAC axis)
- H2: CANCER_COMPLICATION_MANAGEMENT (mitotic/aneuploidy-selective action OR curated
      indication in MVA-associated pediatric malignancies: rhabdomyosarcoma, Wilms tumor,
      leukemia; carries mandatory pediatric-toxicity caveat)
- H3: SUPPORTIVE_CARE_SYMPTOM_MANAGEMENT (renal/nephrocalcinosis/ion-transport or
      neurodevelopmental supportive axis with NO signed path to mitosis/SAC)
- H0: NO_CLASS (no resolvable mechanism)

Pure, deterministic, unit-testable functions. A single pair may receive multiple
classifications, emitting duplicated records with per-class justification.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

from mva_hackathon.mva_gene_sets import (
    CORE_MVA_CAUSAL_GENES,
    NEPHROCALCINOSIS_ION_GENES,
    RMS_CANCER_GENES,
    SAC_CORE_GENES,
)

RESEARCH_ONLY_DISCLAIMER = (
    "RESEARCH EXPLORATION ONLY — STRICTLY NOT CLINICAL OR TREATMENT ADVICE. "
    "Hypotheses generated computationally from public databases for mechanistic inquiry."
)

CLASS_H1 = "H1_DISEASE_MECHANISM_CORRECTION"
CLASS_H2 = "H2_CANCER_COMPLICATION_MANAGEMENT"
CLASS_H3 = "H3_SUPPORTIVE_CARE_SYMPTOM_MANAGEMENT"
CLASS_H0 = "H0_NO_CLASS"

# MVA-associated pediatric malignancies (EFO / text matching)
MVA_MALIGNANCY_KEYWORDS = {
    "rhabdomyosarcoma",
    "wilms",
    "nephroblastoma",
    "leukemia",
    "acute lymphoblastic leukemia",
    "acute myeloid leukemia",
    "embryonal rhabdomyosarcoma",
    "sarcoma",
}

# Endpoints for H1
SAC_ENDPOINTS = {
    "BUB1B", "MAD1L1", "CDC20", "BUB1", "BUB3", "MAD2L1", "MAD2L2",
    "AURKA", "PLK1", "CENPE", "KNL1", "NDC80", "TTK",
}
SPLICING_ENDPOINTS = {"CENATAC", "CCDC84", "RNU11", "RNU12"}
CHROMOSOME_ENDPOINTS = {"CEP57", "SLF2", "SMC5", "TP53"}


@dataclass
class HypothesisClassificationResult:
    """Classified hypothesis row with evidence metrics and justifications."""

    target_symbol: str
    drug_name: str
    gene_category: str
    hypothesis_class: str
    class_justification: str
    endpoint_hit: str
    signed_path_count_curated: int
    signed_path_count_inferred: int
    net_flux_agreement: bool
    indication_evidence: str
    indication_efo_ids: str
    prior_safety_tier: str
    mechanism_confidence: float
    literature_support_status: str
    evidence_gap_flags: str
    pediatric_toxicity_caveat: str = ""
    research_only_disclaimer: str = RESEARCH_ONLY_DISCLAIMER

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def classify_pair_hypotheses(
    target_symbol: str,
    drug_name: str,
    gene_category: str,
    action_type: str,
    drug_sign: int,
    signed_paths: list[Any],  # list of SignedPath
    reactions: list[Any],     # list of AdjudicatedReaction
    prior_safety_tier: str,
    disease_indications: list[str] | None = None,
    efo_ids: list[str] | None = None,
    literature_verified: bool = False,
    is_evidence_gap: bool = False,
) -> list[HypothesisClassificationResult]:
    """Deterministically classify a drug-target pair into one or more hypothesis classes.

    Returns a list of HypothesisClassificationResult (duplicated rows if multiple classes qualify).
    """
    target = target_symbol.upper()
    indications = [ind.lower() for ind in (disease_indications or [])]
    efo_list = efo_ids or []

    # Count curated vs inferred paths
    curated_paths = [p for p in signed_paths if getattr(p, "is_curated", True)]
    inferred_paths = [p for p in signed_paths if not getattr(p, "is_curated", True)]
    all_endpoints = {getattr(p, "end_node", "").upper() for p in signed_paths}

    results: list[HypothesisClassificationResult] = []

    if is_evidence_gap:
        results.append(
            HypothesisClassificationResult(
                target_symbol=target,
                drug_name=drug_name,
                gene_category=gene_category,
                hypothesis_class=CLASS_H0,
                class_justification="EVIDENCE_GAP: Target lacks resolvable signed regulatory edges or biochemical reactions in curated repositories.",
                endpoint_hit="NONE",
                signed_path_count_curated=0,
                signed_path_count_inferred=0,
                net_flux_agreement=False,
                indication_evidence="NONE",
                indication_efo_ids="",
                prior_safety_tier=prior_safety_tier or "TIER_3_UNCERTAIN_DIRECTIONALITY",
                mechanism_confidence=0.10,
                literature_support_status="PENDING/UNVERIFIED",
                evidence_gap_flags="NO_SIGNED_PATHS;NO_REACTIONS",
            )
        )
        return results

    # -------------------------------------------------------------
    # Check H1: DISEASE_MECHANISM_CORRECTION
    # -------------------------------------------------------------
    # Condition: Target connects via signed path or reaction to SAC, KT-MT, or splicing
    h1_endpoints = (SAC_ENDPOINTS | SPLICING_ENDPOINTS | CHROMOSOME_ENDPOINTS)
    h1_hits = (all_endpoints & h1_endpoints) or ({target} if target in h1_endpoints else set())

    # Check for direct or path-mediated correction
    h1_qualifies = False
    h1_endpoint_name = "NONE"
    h1_net_agreement = False
    h1_confidence = 0.50

    if h1_hits:
        h1_endpoint_name = ";".join(sorted(h1_hits))
        # Directional logic:
        # If target is SAC inhibitor (e.g. CDC20) and drug is inhibitor (drug_sign == -1) ->
        # net effect = (-1)*(-1) = +1 (relieves inhibition -> checkpoint prolonging / stabilization).
        # If target is SAC kinase (e.g. PLK1/AURKA) and drug is inhibitor:
        # while commonly cytotoxic, partial modulation can alter tension, but complete inhibition causes mitotic exit failure.
        # Check net signs of paths reaching SAC endpoints
        concordant_paths = [p for p in signed_paths if getattr(p, "net_sign", 0) * drug_sign > 0 and getattr(p, "end_node", "") in h1_endpoints]
        if concordant_paths or target in CORE_MVA_CAUSAL_GENES:
            h1_qualifies = True
            h1_net_agreement = bool(concordant_paths)
            h1_confidence = 0.85 if (curated_paths and literature_verified) else 0.65
        elif target in SAC_CORE_GENES and drug_sign != 0:
            # Kinase / checkpoint modulator acting on mitotic spindle
            h1_qualifies = True
            h1_net_agreement = False
            h1_confidence = 0.55

    if h1_qualifies:
        justification = (
            f"Signed path/reaction links target {target} to core MVA deficit endpoints ({h1_endpoint_name}). "
            f"Perturbation modifies spindle assembly checkpoint (SAC) kinetics or kinetochore attachment fidelity."
        )
        results.append(
            HypothesisClassificationResult(
                target_symbol=target,
                drug_name=drug_name,
                gene_category=gene_category,
                hypothesis_class=CLASS_H1,
                class_justification=justification,
                endpoint_hit=h1_endpoint_name,
                signed_path_count_curated=len(curated_paths),
                signed_path_count_inferred=len(inferred_paths),
                net_flux_agreement=h1_net_agreement,
                indication_evidence=";".join(disease_indications) if disease_indications else "NONE",
                indication_efo_ids=";".join(efo_list),
                prior_safety_tier=prior_safety_tier,
                mechanism_confidence=h1_confidence,
                literature_support_status="VERIFIED" if literature_verified else "PENDING/UNVERIFIED",
                evidence_gap_flags="NONE",
            )
        )

    # -------------------------------------------------------------
    # Check H2: CANCER_COMPLICATION_MANAGEMENT
    # -------------------------------------------------------------
    # Condition: Curated indication in MVA-associated pediatric malignancy OR target in RMS/cancer set
    has_mva_cancer_indication = any(
        any(kw in ind for kw in MVA_MALIGNANCY_KEYWORDS) for ind in indications
    )
    is_cancer_gene = target in RMS_CANCER_GENES or target in {"ALK", "KIT", "PARP1", "TP53"}

    if has_mva_cancer_indication or is_cancer_gene:
        matched_malignancies = [
            ind for ind in indications if any(kw in ind for kw in MVA_MALIGNANCY_KEYWORDS)
        ]
        mal_str = ";".join(matched_malignancies) if matched_malignancies else "Pediatric oncologic vulnerability / aneuploidy surveillance"
        justification = (
            f"Target {target} is implicated in pediatric malignancy management ({mal_str}) "
            f"relevant to MVA cancer predisposition (rhabdomyosarcoma, Wilms tumor, leukemia)."
        )
        caveat = (
            "MANDATORY PEDIATRIC CAVEAT: Constitutional chromosomal instability in MVA creates heightened "
            "vulnerability to myelosuppression, growth retardation, and secondary treatment-induced malignancies."
        )
        h2_confidence = 0.80 if (has_mva_cancer_indication and literature_verified) else 0.60
        results.append(
            HypothesisClassificationResult(
                target_symbol=target,
                drug_name=drug_name,
                gene_category=gene_category,
                hypothesis_class=CLASS_H2,
                class_justification=justification,
                endpoint_hit="PEDIATRIC_MALIGNANCY_RMS_WILMS_LEUKEMIA",
                signed_path_count_curated=len(curated_paths),
                signed_path_count_inferred=len(inferred_paths),
                net_flux_agreement=True,
                indication_evidence=";".join(disease_indications) if disease_indications else "ONCOLOGIC_TARGET",
                indication_efo_ids=";".join(efo_list),
                prior_safety_tier=prior_safety_tier,
                mechanism_confidence=h2_confidence,
                literature_support_status="VERIFIED" if literature_verified else "PENDING/UNVERIFIED",
                evidence_gap_flags="NONE",
                pediatric_toxicity_caveat=caveat,
            )
        )

    # -------------------------------------------------------------
    # Check H3: SUPPORTIVE_CARE_SYMPTOM_MANAGEMENT
    # -------------------------------------------------------------
    # Condition: Target in renal/nephrocalcinosis/ion transport OR non-mitotic supportive care
    # with NO signed path to mitosis/SAC
    is_nephro_gene = target in NEPHROCALCINOSIS_ION_GENES
    has_renal_indication = any("kidney" in ind or "renal" in ind or "nephro" in ind for ind in indications)
    no_sac_paths = not bool(all_endpoints & SAC_ENDPOINTS)

    if (is_nephro_gene or has_renal_indication) and no_sac_paths:
        justification = (
            f"Target {target} acts in supportive care / renal tubular transport / nephrocalcinosis axis "
            f"without directly perturbing the mitotic spindle assembly checkpoint."
        )
        results.append(
            HypothesisClassificationResult(
                target_symbol=target,
                drug_name=drug_name,
                gene_category=gene_category,
                hypothesis_class=CLASS_H3,
                class_justification=justification,
                endpoint_hit="RENAL_SUPPORTIVE_AXIS",
                signed_path_count_curated=len(curated_paths),
                signed_path_count_inferred=len(inferred_paths),
                net_flux_agreement=True,
                indication_evidence=";".join(disease_indications) if disease_indications else "SUPPORTIVE_CARE",
                indication_efo_ids=";".join(efo_list),
                prior_safety_tier=prior_safety_tier,
                mechanism_confidence=0.70 if literature_verified else 0.50,
                literature_support_status="VERIFIED" if literature_verified else "PENDING/UNVERIFIED",
                evidence_gap_flags="NONE",
            )
        )

    # -------------------------------------------------------------
    # Fallback to H0 if no class met
    # -------------------------------------------------------------
    if not results:
        results.append(
            HypothesisClassificationResult(
                target_symbol=target,
                drug_name=drug_name,
                gene_category=gene_category,
                hypothesis_class=CLASS_H0,
                class_justification="Target-drug pair does not satisfy specific criteria for H1 (mechanism correction), H2 (cancer management), or H3 (supportive care).",
                endpoint_hit="NONE",
                signed_path_count_curated=len(curated_paths),
                signed_path_count_inferred=len(inferred_paths),
                net_flux_agreement=False,
                indication_evidence=";".join(disease_indications) if disease_indications else "NONE",
                indication_efo_ids=";".join(efo_list),
                prior_safety_tier=prior_safety_tier or "TIER_3_UNCERTAIN_DIRECTIONALITY",
                mechanism_confidence=0.30,
                literature_support_status="PENDING/UNVERIFIED",
                evidence_gap_flags="UNCERTAIN_MECHANISM",
            )
        )

    return results
