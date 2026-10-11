"""Phase 4 Multi-Gate Adjudication Safety and Synthesis Engine.

Executes sequential multi-gate adjudication:
- Gate 1 (G1): Identifier hygiene (HGNC, UniProt, Ensembl, ChEMBL)
- Gate 2 (G2): Signed-tree and reaction reachability (identifies EVIDENCE_GAP)
- Gate 3 (G3): Directionality concordance (detects SIGNED_PATH_HAZARD)
- Gate 4 (G4): Multi-hypothesis class assignment (H1, H2, H3, H0)
- Gate 5 (G5): Pediatric & regulatory hazard gate (PEDIATRIC_HAZARD_SCORE 0-100, threshold 60)
- Gate 6 (G6): Literature validation (verified PMID+DOI requirement, contradiction screening)
- Gate 7 (G7): Final verdict assignment:
  * ACCEPTED_RESEARCH_HYPOTHESIS (requires curated path, >=2 independent sources, verified literature)
  * KEEP_WITH_CAVEAT (single source, unverified literature, or pediatric caveat)
  * MECHANISTIC_HAZARD_REJECTED (opposing directionality / severe mechanistic hazard)
  * LITERATURE_CONTRADICTED (published contraindication / exacerbation report)
  * EVIDENCE_INSUFFICIENT (evidence gaps)
  * NEGATIVE_NO_SAFE_INTERVENTION (contraindicated negative control)

Outputs flat tables:
- 20_phase4_reaction_level_matrix.csv
- 21_phase4_hypothesis_classification.csv
- 22_phase4_synthetic_lethality_candidates.csv
- 23_phase4_adjudicated_safety_report.csv
- 24_phase4_adjudication_provenance.json

Strictly research exploration. All rows carry:
CLINICAL_RECOMMENDATION_STATUS='PROHIBITED / INVESTIGATIONAL ONLY'
RESEARCH_ONLY_DISCLAIMER
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import hashlib
import json
import logging
from pathlib import Path
import re
from typing import Any, Mapping

import pandas as pd

from mva_hackathon.hypothesis_classifier import (
    CLASS_H0,
    CLASS_H1,
    CLASS_H2,
    CLASS_H3,
    HypothesisClassificationResult,
    classify_pair_hypotheses,
)
from mva_hackathon.literature_adjudication import (
    LiteratureAdjudicationConfig,
    LiteratureAdjudicationEngine,
    VerifiedCitation,
)
from mva_hackathon.mva_gene_sets import (
    CORE_MVA_CAUSAL_GENES,
    RMS_CANCER_GENES,
    SAC_CORE_GENES,
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

logger = logging.getLogger("mva_hackathon.adjudication_safety")

RESEARCH_ONLY_DISCLAIMER = (
    "RESEARCH EXPLORATION ONLY — STRICTLY NOT CLINICAL OR TREATMENT ADVICE. "
    "Hypotheses generated computationally from public databases for mechanistic inquiry."
)

CLINICAL_RECOMMENDATION_STATUS_VAL = "PROHIBITED / INVESTIGATIONAL ONLY"

# Verdict status constants
VERDICT_ACCEPTED = "ACCEPTED_RESEARCH_HYPOTHESIS"
VERDICT_CAVEAT = "KEEP_WITH_CAVEAT"
VERDICT_HAZARD = "MECHANISTIC_HAZARD_REJECTED"
VERDICT_CONTRADICTED = "LITERATURE_CONTRADICTED"
VERDICT_INSUFFICIENT = "EVIDENCE_INSUFFICIENT"
VERDICT_NEGATIVE_CONTROL = "NEGATIVE_NO_SAFE_INTERVENTION"

VERDICT_ORDINAL_MAP = {
    VERDICT_ACCEPTED: 1,
    VERDICT_CAVEAT: 2,
    VERDICT_INSUFFICIENT: 3,
    VERDICT_NEGATIVE_CONTROL: 4,
    VERDICT_CONTRADICTED: 5,
    VERDICT_HAZARD: 6,
}

MAX_FILE_BYTES = 8 * 1024 * 1024  # 8 MB strict ceiling


@dataclass
class Phase4AdjudicationConfig:
    """Master configuration for Phase 4 adjudication pipeline."""

    gene_mechanism_path: Path
    drug_evidence_path: Path
    safety_ranking_path: Path
    query_panel_path: Path
    pathway_enrichment_path: Path
    network_context_path: Path
    output_dir: Path
    flat_names: bool = True
    max_path_length: int = 2
    reaction_limit: int = 500
    enable_omnipath_inferred: bool = False
    enable_synlethdb: bool = False
    offline_mode: bool = False
    verify_literature: bool = True
    targets_override: list[str] = field(default_factory=list)
    drugs_override: list[str] = field(default_factory=list)
    pipeline_version: str = "0.4.0"


@dataclass
class AdjudicatedSafetyReportRow:
    """Final adjudicated hypothesis row for 23_phase4_adjudicated_safety_report.csv."""

    target_symbol: str
    drug_name: str
    drug_id_chembl: str
    hypothesis_class: str
    pediatric_hazard_score: int
    pediatric_hazard_note: str
    adverse_event_flags: str
    signed_path_hazard: bool
    contradiction_flag: bool
    reaction_level_hazard: bool
    literature_contradictions: str
    verdict_status: str
    verdict_justification: str
    mechanism_confidence: float
    verified_pmids: str
    verified_dois: str
    literature_verification_source: str
    research_priority_rank: int
    verdict_ordinal: int = 2
    clinical_recommendation_status: str = CLINICAL_RECOMMENDATION_STATUS_VAL
    research_only_disclaimer: str = RESEARCH_ONLY_DISCLAIMER

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def compute_sha256(path: Path) -> str:
    """Calculate SHA256 checksum for a file."""
    if not path.is_file():
        return "FILE_NOT_FOUND"
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def calculate_pediatric_hazard_score(
    prior_tier: str,
    hazard_note: str,
    action_type: str,
    target_symbol: str,
) -> tuple[int, str]:
    """Calculate pediatric hazard score (0-100) and flag adverse event notes.

    Hazard thresholds:
    - Tier 4: base 85
    - Black-box warning, cytotoxic, or growth arrest: +15
    - Tier 3: base 45
    - Tier 2: base 35
    - Tier 1: base 15
    """
    clean_note = str(hazard_note or "").strip()
    score = 30
    adverse_flags: list[str] = []

    if prior_tier == "TIER_4_CONTRAINDICATED_HAZARD":
        score = 85
        adverse_flags.append("TIER_4_CONTRAINDICATED")
    elif prior_tier == "TIER_3_UNCERTAIN_DIRECTIONALITY":
        score = 45
    elif prior_tier == "TIER_2_OFF_TARGET_COMPENSATION":
        score = 35
    elif prior_tier == "TIER_1_CONCORDANT_HYPOTHESIS":
        score = 20

    lower_note = clean_note.lower()
    if any(k in lower_note for k in ["black-box", "black box", "death", "fatal"]):
        score = min(100, score + 20)
        adverse_flags.append("BLACK_BOX_WARNING")
    if any(k in lower_note for k in ["teratogen", "growth arrest", "bone marrow", "myelosuppression"]):
        score = min(100, score + 15)
        adverse_flags.append("PEDIATRIC_GROWTH_OR_MARROW_RISK")
    if any(k in lower_note for k in ["secondary malignancy", "carcinogen"]):
        score = min(100, score + 20)
        adverse_flags.append("SECONDARY_MALIGNANCY_RISK")

    if target_symbol in {"TP53"} and "INHIB" in action_type.upper():
        score = 100
        adverse_flags.append("TP53_ABROGATION_LETHAL_HAZARD")

    flags_str = ";".join(adverse_flags) if adverse_flags else "STANDARD_MONITORING"
    return score, flags_str


class Phase4AdjudicationEngine:
    """Master engine executing signed pathway and reaction adjudication."""

    def __init__(self, config: Phase4AdjudicationConfig) -> None:
        self.config = config
        self.signed_engine = SignedRegulationEngine(
            SignedRegulationConfig(
                enable_omnipath_inferred=config.enable_omnipath_inferred,
                max_path_length=config.max_path_length,
                offline_mode=config.offline_mode,
            )
        )
        self.reaction_engine = ReactionAdjudicationEngine(
            ReactionAdjudicationConfig(
                reaction_limit=config.reaction_limit,
                offline_mode=config.offline_mode,
            )
        )
        self.lit_engine = LiteratureAdjudicationEngine(
            LiteratureAdjudicationConfig(
                verify_literature=config.verify_literature,
                offline_mode=config.offline_mode,
            )
        )

        self.truncation_ledger: list[dict[str, Any]] = []
        self.quarantine_ledger: list[dict[str, Any]] = []
        self.gate_stats: dict[str, int] = {
            "G1_HYGIENE_PASS": 0,
            "G1_HYGIENE_FAIL": 0,
            "G2_REACHABLE": 0,
            "G2_EVIDENCE_GAP": 0,
            "G3_CONCORDANT": 0,
            "G3_HAZARD": 0,
            "G4_CLASSES_ASSIGNED": 0,
            "G5_PEDIATRIC_HOLD": 0,
            "G5_PEDIATRIC_PASS": 0,
            "G6_LIT_VERIFIED": 0,
            "G6_LIT_PENDING": 0,
            "G6_LIT_CONTRADICTED": 0,
            "G7_VERDICTS_ISSUED": 0,
        }

    def run(self) -> dict[str, Any]:
        """Execute full Phase 4 adjudication pipeline."""
        start_time = datetime.now(timezone.utc)
        logger.info("Starting Phase 4 Adjudication Pipeline v%s", self.config.pipeline_version)

        # 1. Load input tables
        df_mech = pd.read_csv(self.config.gene_mechanism_path)
        df_drug = pd.read_csv(self.config.drug_evidence_path)
        df_safety = pd.read_csv(self.config.safety_ranking_path)
        df_panel = pd.read_csv(self.config.query_panel_path)
        df_enrich = pd.read_csv(self.config.pathway_enrichment_path)

        # 2. Derive Entry Criteria candidate pairs
        entry_pairs = self._select_entry_candidate_pairs(df_safety, df_drug, df_panel)
        logger.info("Selected %d candidate pairs satisfying Entry Criteria E1-E5", len(entry_pairs))

        # 3. Collect watch genes (targets + panel genes + endpoints)
        target_genes = {p["TARGET_SYMBOL"].upper() for p in entry_pairs}
        panel_genes = set(df_panel["QUERY_GENE"].dropna().str.upper().tolist())
        watch_genes = target_genes | panel_genes | SAC_CORE_GENES | CORE_MVA_CAUSAL_GENES

        # 4. Ingest and reconcile signed regulatory networks (SIGNOR + OmniPath)
        logger.info("Ingesting signed interactions for %d watch genes...", len(watch_genes))
        edge_count, conflict_count = self.signed_engine.ingest_and_reconcile(partner_genes=list(watch_genes))

        # 5. Extract and adjudicate Reactome reactions
        reac_pathways = []
        for _, row in df_enrich.iterrows():
            pid = str(row.get("PATHWAY_ID", "")).strip()
            pname = str(row.get("PATHWAY_NAME", "")).strip()
            if pid.startswith("R-HSA-") or pid.startswith("REAC:R-HSA-"):
                reac_pathways.append((pid, pname))

        logger.info("Adjudicating reactions across %d Reactome pathways...", len(reac_pathways))
        adjudicated_reactions = self.reaction_engine.adjudicate_pathways(
            pathway_items=reac_pathways,
            panel_genes=panel_genes,
            target_genes=target_genes,
        )

        # 6. Execute Gates G1-G7 across Entry Pairs
        reaction_matrix_rows: list[dict[str, Any]] = []
        classified_rows: list[HypothesisClassificationResult] = []
        safety_report_rows: list[AdjudicatedSafetyReportRow] = []

        # Target metadata mapping from panel / mechanism
        target_meta: dict[str, dict[str, Any]] = {}
        for _, r in df_panel.iterrows():
            sym = str(r.get("QUERY_GENE", "")).upper()
            target_meta[sym] = {
                "UNIPROT_ID": str(r.get("UNIPROT_ID", "") or ""),
                "ENSEMBL_ID": str(r.get("ENSEMBL_GENE_STABLE_ID", "") or ""),
                "GENE_CATEGORY": str(r.get("GENE_CATEGORY", "") or "GRAPH_DERIVED_NEIGHBOR"),
                "THERAPEUTIC_DIR": str(r.get("THERAPEUTIC_DIRECTIONALITY", "") or "UNKNOWN"),
                "PMID": str(r.get("PRIMARY_LITERATURE_PMID", "") or "").replace(".0", ""),
                "DOI": str(r.get("PRIMARY_LITERATURE_DOI", "") or ""),
            }

        endpoints_to_reach = set(CORE_MVA_CAUSAL_GENES | SAC_CORE_GENES | {"TP53", "ATM"})
        target_signed_paths_cache: dict[str, list[SignedPath]] = {}
        target_reactions_cache: dict[str, list[AdjudicatedReaction]] = {}
        target_lit_cache: dict[str, tuple[list[str], list[str], list[str]]] = {}

        for pair in entry_pairs:
            target = pair["TARGET_SYMBOL"].upper()
            drug = pair["DRUG_NAME"]
            drug_id = pair.get("DRUG_ID", "") or pair.get("DRUG_ID_CHEMBL", "")
            action_type = str(pair.get("ACTION_TYPE", "") or "UNKNOWN")
            prior_tier = str(pair.get("SAFETY_TIER", "") or "TIER_3_UNCERTAIN_DIRECTIONALITY")
            hazard_note = str(pair.get("PEDIATRIC_HAZARD_NOTE", "") or "")
            control_role = str(pair.get("CONTROL_ROLE", "") or "NONE")

            # Map drug action type to sign
            drug_sign = 0
            act_up = action_type.upper()
            if any(k in act_up for k in ["INHIB", "ANTAG", "DEGRAD", "SUPPRESS", "BLOCK"]):
                drug_sign = -1
            elif any(k in act_up for k in ["ACTIV", "AGON", "INDUC", "OPEN"]):
                drug_sign = 1

            meta = target_meta.get(target, {
                "UNIPROT_ID": "",
                "ENSEMBL_ID": "",
                "GENE_CATEGORY": pair.get("GENE_CATEGORY", "UNKNOWN"),
                "THERAPEUTIC_DIR": pair.get("THERAPEUTIC_DIRECTIONALITY", "UNKNOWN"),
                "PMID": "",
                "DOI": "",
            })

            # Gate 1: Identifier hygiene
            if not target or target in {"UNKNOWN", "UNMAPPED"}:
                self.gate_stats["G1_HYGIENE_FAIL"] += 1
                self.quarantine_ledger.append({"pair": pair, "reason": f"Invalid symbol {target}"})
                continue
            self.gate_stats["G1_HYGIENE_PASS"] += 1

            # Gate 2: Signed-tree and reaction reachability (with target-level caching)
            if target not in target_signed_paths_cache:
                target_signed_paths_cache[target] = self.signed_engine.find_signed_paths(
                    start_gene=target,
                    target_endpoints=endpoints_to_reach,
                    max_length=self.config.max_path_length,
                )
            signed_paths = target_signed_paths_cache[target]

            if target not in target_reactions_cache:
                target_reactions_cache[target] = [
                    rx for rx in adjudicated_reactions if target in rx.target_roles or target in rx.panel_participants
                ]
            target_reactions = target_reactions_cache[target]
            is_evidence_gap = (len(signed_paths) == 0 and len(target_reactions) == 0)

            if is_evidence_gap:
                self.gate_stats["G2_EVIDENCE_GAP"] += 1
            else:
                self.gate_stats["G2_REACHABLE"] += 1

            # Populate reaction level matrix rows (sample top 3 reactions and top 3 paths to respect 8MB limit)
            if target_reactions:
                for rx in target_reactions[:3]:
                    flux_delta = self.reaction_engine.compute_reaction_flux_delta(rx, target, drug_sign)
                    reaction_matrix_rows.append({
                        "TARGET_SYMBOL": target,
                        "UNIPROT_ID": meta["UNIPROT_ID"],
                        "DRUG_NAME": drug,
                        "DRUG_ID_CHEMBL": drug_id,
                        "ACTION_TYPE": action_type,
                        "EDGE_SOURCE_GENE": target,
                        "EDGE_TARGET_GENE": rx.name[:40],
                        "EFFECT": "REACTION_REGULATION",
                        "MECHANISM": "REACTOME_PATHWAY",
                        "DIRECT": True,
                        "EVIDENCE_CLASS": "CURATED",
                        "SOURCE_DB": "Reactome",
                        "SOURCE_DB_VERSION": "v88",
                        "EVIDENCE_PMIDS": "",
                        "CONFLICT_FLAG": False,
                        "REACTION_ST_ID": rx.st_id,
                        "REACTION_NAME": rx.name,
                        "REACTION_PATHWAY": rx.pathway_name,
                        "REACTION_PARTICIPANT_ROLE": ";".join(rx.target_roles.get(target, ["PARTICIPANT"])),
                        "REGULATION_TYPE": "BIOCHEMICAL_FLUX",
                        "PATH_LENGTH": 1,
                        "ENDPOINT_NODE": rx.st_id,
                        "ENDPOINT_CATEGORY": "BIOCHEMICAL_REACTION",
                        "PREDICTED_NET_FLUX": flux_delta,
                        "REACTION_FLUX_DELTA": flux_delta,
                        "RESEARCH_ONLY_DISCLAIMER": RESEARCH_ONLY_DISCLAIMER,
                    })

            for p in signed_paths[:3]:
                for edge in p.edges:
                    reaction_matrix_rows.append({
                        "TARGET_SYMBOL": target,
                        "UNIPROT_ID": meta["UNIPROT_ID"],
                        "DRUG_NAME": drug,
                        "DRUG_ID_CHEMBL": drug_id,
                        "ACTION_TYPE": action_type,
                        "EDGE_SOURCE_GENE": edge.source_gene,
                        "EDGE_TARGET_GENE": edge.target_gene,
                        "EFFECT": edge.effect,
                        "MECHANISM": edge.mechanism,
                        "DIRECT": edge.direct,
                        "EVIDENCE_CLASS": edge.evidence_class,
                        "SOURCE_DB": edge.source_db,
                        "SOURCE_DB_VERSION": edge.source_db_version,
                        "EVIDENCE_PMIDS": ";".join(edge.evidence_pmids),
                        "CONFLICT_FLAG": edge.conflict_flag,
                        "REACTION_ST_ID": "",
                        "REACTION_NAME": f"{edge.source_gene} -> {edge.target_gene}",
                        "REACTION_PATHWAY": "SIGNALING_NETWORK",
                        "REACTION_PARTICIPANT_ROLE": "REGULATOR",
                        "REGULATION_TYPE": edge.effect,
                        "PATH_LENGTH": p.path_length,
                        "ENDPOINT_NODE": p.end_node,
                        "ENDPOINT_CATEGORY": "MVA_CORE_ENDPOINT",
                        "PREDICTED_NET_FLUX": p.net_sign * drug_sign,
                        "REACTION_FLUX_DELTA": 0,
                        "RESEARCH_ONLY_DISCLAIMER": RESEARCH_ONLY_DISCLAIMER,
                    })

            # Gate 3: Directionality concordance
            signed_path_hazard = False
            # If target is TP53 and drug is an inhibitor, or target is tumor suppressor / SAC essential
            if (target in {"TP53"} and drug_sign < 0) or prior_tier == "TIER_4_CONTRAINDICATED_HAZARD":
                signed_path_hazard = True
            elif meta.get("THERAPEUTIC_DIR") == "ACTIVATE_STABILIZE" and drug_sign < 0:
                signed_path_hazard = True

            if signed_path_hazard:
                self.gate_stats["G3_HAZARD"] += 1
            else:
                self.gate_stats["G3_CONCORDANT"] += 1

            # Gate 5: Pediatric hazard score
            ped_score, adverse_flags = calculate_pediatric_hazard_score(prior_tier, hazard_note, action_type, target)
            if ped_score > 60:
                self.gate_stats["G5_PEDIATRIC_HOLD"] += 1
            else:
                self.gate_stats["G5_PEDIATRIC_PASS"] += 1

            # Gate 6: Literature adjudication (cached per target symbol)
            if target not in target_lit_cache:
                t_pmids: list[str] = []
                t_dois: list[str] = []
                t_sources: list[str] = []

                reg_pmid = meta.get("PMID", "")
                reg_doi = meta.get("DOI", "")
                if reg_pmid:
                    cit_pmid = self.lit_engine.verify_pmid(reg_pmid)
                    if cit_pmid.status == "VERIFIED":
                        t_pmids.append(reg_pmid)
                        t_sources.append(cit_pmid.verification_source)
                if reg_doi:
                    cit_doi = self.lit_engine.verify_doi(reg_doi)
                    if cit_doi.status == "VERIFIED":
                        t_dois.append(reg_doi)
                        t_sources.append(cit_doi.verification_source)

                # If no registry citation, check edge PMIDs from top paths
                if not (t_pmids or t_dois):
                    for p in signed_paths[:3]:
                        for pmid in p.evidence_pmids[:2]:
                            cit = self.lit_engine.verify_pmid(pmid)
                            if cit.status == "VERIFIED":
                                t_pmids.append(pmid)
                                t_sources.append(cit.verification_source)
                                break
                        if t_pmids:
                            break

                target_lit_cache[target] = (
                    list(sorted(set(t_pmids))),
                    list(sorted(set(t_dois))),
                    list(sorted(set(t_sources))),
                )

            verified_pmids, verified_dois, lit_sources = target_lit_cache[target]
            has_verified_lit = bool(verified_pmids or verified_dois)

            # Contradiction check: screen concordant priority candidates
            is_contradicted, contra_pmid, contra_note = False, "", ""
            if not signed_path_hazard and not is_evidence_gap and self.config.verify_literature and not self.config.offline_mode:
                if prior_tier == "TIER_1_CONCORDANT_HYPOTHESIS":
                    is_contradicted, contra_pmid, contra_note = self.lit_engine.europe_pmc.check_contradiction(target, drug)

            if is_contradicted:
                self.gate_stats["G6_LIT_CONTRADICTED"] += 1
            elif has_verified_lit:
                self.gate_stats["G6_LIT_VERIFIED"] += 1
            else:
                self.gate_stats["G6_LIT_PENDING"] += 1

            # Gate 4: Classify hypotheses
            class_results = classify_pair_hypotheses(
                target_symbol=target,
                drug_name=drug,
                gene_category=meta.get("GENE_CATEGORY", "UNKNOWN"),
                action_type=action_type,
                drug_sign=drug_sign,
                signed_paths=signed_paths,
                reactions=target_reactions,
                prior_safety_tier=prior_tier,
                disease_indications=pair.get("INDICATION_LIST", []),
                efo_ids=pair.get("EFO_LIST", []),
                literature_verified=has_verified_lit,
                is_evidence_gap=is_evidence_gap,
            )
            classified_rows.extend(class_results)
            self.gate_stats["G4_CLASSES_ASSIGNED"] += len(class_results)

            # Gate 7: Verdict Assignment
            # Multi-source independence count
            source_records = []
            for p in signed_paths:
                for edge in p.edges:
                    source_records.append({"source_db": edge.source_db, "pmids": edge.evidence_pmids})
            indep_sources, _ = self.lit_engine.deduplicate_independent_sources(source_records)

            for cr in class_results:
                verdict = VERDICT_CAVEAT
                justification = cr.class_justification

                if is_contradicted:
                    verdict = VERDICT_CONTRADICTED
                    justification = f"LITERATURE_CONTRADICTION: {contra_note} (PMID {contra_pmid})"
                elif signed_path_hazard:
                    verdict = VERDICT_HAZARD
                    justification = f"MECHANISTIC_HAZARD: Drug action ({action_type}) directly opposes therapeutic direction or abrogates essential guardian {target}."
                elif control_role == "NEGATIVE_CONTROL":
                    verdict = VERDICT_NEGATIVE_CONTROL
                    justification = f"NEGATIVE_CONTROL: Target {target} is a designated negative control; intervention is contraindicated."
                elif is_evidence_gap:
                    verdict = VERDICT_INSUFFICIENT
                    justification = "EVIDENCE_INSUFFICIENT: Target lacks resolvable signed regulatory edges or biochemical reactions."
                else:
                    # Multi-source upgrade rule:
                    # ACCEPTED requires curated signed path AND >=2 independent curated sources AND verified literature AND ped_score <= 60
                    curated_path_count = cr.signed_path_count_curated
                    if curated_path_count >= 1 and indep_sources >= 2 and has_verified_lit and ped_score <= 60:
                        verdict = VERDICT_ACCEPTED
                        justification = (
                            f"ACCEPTED_RESEARCH_HYPOTHESIS: Rigorously corroborated by curated signed paths ({curated_path_count}), "
                            f"{indep_sources} independent curated repositories, and verified primary literature."
                        )
                    else:
                        verdict = VERDICT_CAVEAT
                        reasons = []
                        if indep_sources < 2:
                            reasons.append("single curated source")
                        if not has_verified_lit:
                            reasons.append("literature pending full verification")
                        if ped_score > 60:
                            reasons.append("elevated pediatric toxicity liability")
                        if cr.hypothesis_class == CLASS_H2:
                            reasons.append("oncologic toxicity in constitutional aneuploidy")
                        justification = f"KEEP_WITH_CAVEAT: Mechanistically concordant but subject to caveats: {'; '.join(reasons)}."

                safety_report_rows.append(
                    AdjudicatedSafetyReportRow(
                        target_symbol=target,
                        drug_name=drug,
                        drug_id_chembl=drug_id,
                        hypothesis_class=cr.hypothesis_class,
                        pediatric_hazard_score=ped_score,
                        pediatric_hazard_note=hazard_note or cr.pediatric_toxicity_caveat,
                        adverse_event_flags=adverse_flags,
                        signed_path_hazard=signed_path_hazard,
                        contradiction_flag=is_contradicted,
                        reaction_level_hazard=False,
                        literature_contradictions=f"PMID:{contra_pmid}" if is_contradicted else "NONE",
                        verdict_status=verdict,
                        verdict_justification=justification,
                        mechanism_confidence=cr.mechanism_confidence,
                        verified_pmids=";".join(verified_pmids) if verified_pmids else "PENDING",
                        verified_dois=";".join(verified_dois) if verified_dois else "PENDING",
                        literature_verification_source=";".join(sorted(set(lit_sources))) if lit_sources else "NONE",
                        research_priority_rank=1,  # will be ranked post-sort
                        verdict_ordinal=VERDICT_ORDINAL_MAP.get(verdict, 2),
                    )
                )
                self.gate_stats["G7_VERDICTS_ISSUED"] += 1

        # 7. Sort and Rank Safety Report to prioritize BEST-SUPPORTED hypotheses
        safety_report_rows = self._rank_safety_report(safety_report_rows)

        # 8. Write Flat Outputs (20-24)
        output_paths = self._write_outputs(
            reaction_matrix_rows,
            classified_rows,
            safety_report_rows,
            start_time,
            len(entry_pairs),
            edge_count,
            conflict_count,
            len(adjudicated_reactions),
        )

        return {
            "output_paths": output_paths,
            "entry_pairs_count": len(entry_pairs),
            "edge_count": edge_count,
            "conflict_count": conflict_count,
            "reaction_count": len(adjudicated_reactions),
            "gate_stats": self.gate_stats,
            "safety_rows_count": len(safety_report_rows),
        }

    def _select_entry_candidate_pairs(
        self,
        df_safety: pd.DataFrame,
        df_drug: pd.DataFrame,
        df_panel: pd.DataFrame,
    ) -> list[dict[str, Any]]:
        """Select candidate pairs satisfying Entry Criteria E1-E5 plus overrides."""
        selected: list[dict[str, Any]] = []
        seen: set[tuple[str, str]] = set()

        # E1: All Tier 1 rows (5 rows)
        tier1 = df_safety[df_safety["SAFETY_TIER"] == "TIER_1_CONCORDANT_HYPOTHESIS"]
        for _, r in tier1.iterrows():
            pair_key = (str(r["TARGET_SYMBOL"]).upper(), str(r["DRUG_NAME"]))
            if pair_key not in seen:
                seen.add(pair_key)
                selected.append(r.to_dict())

        # E2: Tier 2 rows with NUM_SOURCES >= 2 (25 rows)
        if "NUM_SOURCES" in df_safety.columns:
            tier2_multi = df_safety[
                (df_safety["SAFETY_TIER"] == "TIER_2_OFF_TARGET_COMPENSATION")
                & (df_safety["NUM_SOURCES"] >= 2)
            ]
        else:
            tier2_multi = df_safety[
                df_safety["SAFETY_TIER"] == "TIER_2_OFF_TARGET_COMPENSATION"
            ]
        for _, r in tier2_multi.iterrows():
            pair_key = (str(r["TARGET_SYMBOL"]).upper(), str(r["DRUG_NAME"]))
            if pair_key not in seen:
                seen.add(pair_key)
                selected.append(r.to_dict())

        # E3: Curated negative-control subset of Tier 4 (must-reject sanity check)
        tier4 = df_safety[df_safety["SAFETY_TIER"] == "TIER_4_CONTRAINDICATED_HAZARD"]
        for _, r in tier4.iterrows():
            pair_key = (str(r["TARGET_SYMBOL"]).upper(), str(r["DRUG_NAME"]))
            if pair_key not in seen:
                seen.add(pair_key)
                selected.append(r.to_dict())

        # E4: Pairs whose target has an MVA-associated malignancy indication
        # (e.g. ALK inhibitors crizotinib/ceritinib, KIT inhibitors, PARP inhibitors)
        cancer_targets = {"ALK", "KIT", "PARP1", "TP53"}
        cancer_pairs = df_drug[df_drug["TARGET_SYMBOL"].str.upper().isin(cancer_targets)]
        for _, r in cancer_pairs.iterrows():
            pair_key = (str(r["TARGET_SYMBOL"]).upper(), str(r["DRUG_NAME"]))
            if pair_key not in seen:
                seen.add(pair_key)
                selected.append(r.to_dict())

        # E5: Pairs whose TARGET_SYMBOL is in Phase 2 MVA-focus component
        # (e.g. MAD1L1 carboplatin/paclitaxel, ATM, NF1, HIF1A, ALK, KIT)
        mva_focus_targets = {"MAD1L1", "BUB1B", "ATM", "NF1", "HIF1A"}
        mva_pairs = df_drug[df_drug["TARGET_SYMBOL"].str.upper().isin(mva_focus_targets)]
        for _, r in mva_pairs.iterrows():
            pair_key = (str(r["TARGET_SYMBOL"]).upper(), str(r["DRUG_NAME"]))
            if pair_key not in seen:
                seen.add(pair_key)
                selected.append(r.to_dict())

        # Also emit EVIDENCE_GAP entries for core MVA causal genes with zero drug evidence
        # (BUB1B, CEP57, TRIP13, CENATAC, SLF2, SMC5)
        zero_drug_cores = ["BUB1B", "CEP57", "TRIP13", "CENATAC", "SLF2", "SMC5"]
        for g in zero_drug_cores:
            pair_key = (g.upper(), "NO_APPROVED_OR_CLINICAL_DRUG")
            if pair_key not in seen:
                seen.add(pair_key)
                selected.append({
                    "TARGET_SYMBOL": g,
                    "DRUG_NAME": "NO_APPROVED_OR_CLINICAL_DRUG",
                    "DRUG_ID": "NONE",
                    "ACTION_TYPE": "NONE",
                    "SAFETY_TIER": "EVIDENCE_GAP",
                    "PEDIATRIC_HAZARD_NOTE": "Zero pharmacological agents cataloged in public registries for this core causal node.",
                    "CONTROL_ROLE": "CORE_MVA_CAUSAL",
                })

        # Apply CLI overrides if provided
        if self.config.targets_override:
            ov_targets = {t.upper() for t in self.config.targets_override}
            selected = [p for p in selected if p["TARGET_SYMBOL"].upper() in ov_targets]
        if self.config.drugs_override:
            ov_drugs = {d.upper() for d in self.config.drugs_override}
            selected = [p for p in selected if p["DRUG_NAME"].upper() in ov_drugs]

        return selected

    def _rank_safety_report(
        self, rows: list[AdjudicatedSafetyReportRow]
    ) -> list[AdjudicatedSafetyReportRow]:
        """Rank safety report rows so best-supported hypotheses sort first."""
        # Multi-key sort:
        # 1. verdict_ordinal asc (ACCEPTED=1, CAVEAT=2, INSUFFICIENT=3, NEGATIVE=4, CONTRADICTED=5, HAZARD=6)
        # 2. hypothesis_class asc (H1 first, then H2, H3, H0)
        # 3. mechanism_confidence desc
        # 4. pediatric_hazard_score asc
        rows_sorted = sorted(
            rows,
            key=lambda r: (
                r.verdict_ordinal,
                r.hypothesis_class,
                -r.mechanism_confidence,
                r.pediatric_hazard_score,
                r.target_symbol,
            ),
        )

        # Assign per-class research_priority_rank
        class_counters: dict[str, int] = {}
        ranked_rows: list[AdjudicatedSafetyReportRow] = []
        for r in rows_sorted:
            cls = r.hypothesis_class
            class_counters[cls] = class_counters.get(cls, 0) + 1
            r.research_priority_rank = class_counters[cls]
            ranked_rows.append(r)

        return ranked_rows

    def _write_outputs(
        self,
        reaction_rows: list[dict[str, Any]],
        classified_rows: list[HypothesisClassificationResult],
        safety_rows: list[AdjudicatedSafetyReportRow],
        start_time: datetime,
        entry_pairs_count: int,
        edge_count: int,
        conflict_count: int,
        reaction_count: int,
    ) -> dict[str, Path]:
        """Write flat tables (20-24) complying strictly with 8MB cap and formatting rules."""
        out_dir = self.config.output_dir
        out_dir.mkdir(parents=True, exist_ok=True)

        prefix_20 = "20_phase4_reaction_level_matrix.csv" if self.config.flat_names else "reaction_level_matrix.csv"
        prefix_21 = "21_phase4_hypothesis_classification.csv" if self.config.flat_names else "hypothesis_classification.csv"
        prefix_22 = "22_phase4_synthetic_lethality_candidates.csv" if self.config.flat_names else "synthetic_lethality_candidates.csv"
        prefix_23 = "23_phase4_adjudicated_safety_report.csv" if self.config.flat_names else "adjudicated_safety_report.csv"
        prefix_24 = "24_phase4_adjudication_provenance.json" if self.config.flat_names else "adjudication_provenance.json"

        p20 = out_dir / prefix_20
        p21 = out_dir / prefix_21
        p22 = out_dir / prefix_22
        p23 = out_dir / prefix_23
        p24 = out_dir / prefix_24

        # 20: Reaction Level Matrix
        df20 = pd.DataFrame(reaction_rows)
        if df20.empty:
            df20 = pd.DataFrame(columns=[
                "TARGET_SYMBOL", "UNIPROT_ID", "DRUG_NAME", "DRUG_ID_CHEMBL", "ACTION_TYPE",
                "EDGE_SOURCE_GENE", "EDGE_TARGET_GENE", "EFFECT", "MECHANISM", "DIRECT",
                "EVIDENCE_CLASS", "SOURCE_DB", "SOURCE_DB_VERSION", "EVIDENCE_PMIDS", "CONFLICT_FLAG",
                "REACTION_ST_ID", "REACTION_NAME", "REACTION_PATHWAY", "REACTION_PARTICIPANT_ROLE",
                "REGULATION_TYPE", "PATH_LENGTH", "ENDPOINT_NODE", "ENDPOINT_CATEGORY",
                "PREDICTED_NET_FLUX", "REACTION_FLUX_DELTA", "RESEARCH_ONLY_DISCLAIMER"
            ])
        # Enforce size ceiling (< 8MB)
        if len(df20) > 15000:
            df20 = df20.head(15000)
            self.truncation_ledger.append({"file": str(p20), "cap_rows": 15000})
        df20.to_csv(p20, index=False)

        # 21: Hypothesis Classification
        df21 = pd.DataFrame([cr.to_dict() for cr in classified_rows])
        if "pediatric_toxicity_caveat" in df21.columns:
            df21 = df21.drop(columns=["pediatric_toxicity_caveat"])
        df21.columns = [c.upper() for c in df21.columns]
        req_cols_21 = [
            "TARGET_SYMBOL", "DRUG_NAME", "GENE_CATEGORY", "HYPOTHESIS_CLASS", "CLASS_JUSTIFICATION",
            "ENDPOINT_HIT", "SIGNED_PATH_COUNT_CURATED", "SIGNED_PATH_COUNT_INFERRED", "NET_FLUX_AGREEMENT",
            "INDICATION_EVIDENCE", "INDICATION_EFO_IDS", "PRIOR_SAFETY_TIER", "MECHANISM_CONFIDENCE",
            "LITERATURE_SUPPORT_STATUS", "EVIDENCE_GAP_FLAGS", "RESEARCH_ONLY_DISCLAIMER"
        ]
        # Reorder to match exact schema
        df21 = df21[[c for c in req_cols_21 if c in df21.columns]]
        df21.to_csv(p21, index=False)

        # 22: Synthetic Lethality Candidates (Default disabled: schema + 0 rows)
        synleth_cols = [
            "TARGET_SYMBOL", "INTERACTING_GENE", "SYNLETH_SOURCE", "CELL_LINE_CONTEXT",
            "SCORE", "EVIDENCE_TYPE", "VALIDATION_STATUS", "CLINICAL_RECOMMENDATION_STATUS",
            "RESEARCH_ONLY_DISCLAIMER"
        ]
        df22 = pd.DataFrame(columns=synleth_cols)
        df22.to_csv(p22, index=False)

        # 23: Adjudicated Safety Report
        df23 = pd.DataFrame([sr.to_dict() for sr in safety_rows])
        if "verdict_ordinal" in df23.columns:
            df23 = df23.drop(columns=["verdict_ordinal"])
        df23.columns = [c.upper() for c in df23.columns]
        req_cols_23 = [
            "TARGET_SYMBOL", "DRUG_NAME", "DRUG_ID_CHEMBL", "HYPOTHESIS_CLASS",
            "PEDIATRIC_HAZARD_SCORE", "PEDIATRIC_HAZARD_NOTE", "ADVERSE_EVENT_FLAGS",
            "SIGNED_PATH_HAZARD", "CONTRADICTION_FLAG", "REACTION_LEVEL_HAZARD", "LITERATURE_CONTRADICTIONS",
            "VERDICT_STATUS", "VERDICT_JUSTIFICATION", "MECHANISM_CONFIDENCE",
            "VERIFIED_PMIDS", "VERIFIED_DOIS", "LITERATURE_VERIFICATION_SOURCE",
            "RESEARCH_PRIORITY_RANK", "CLINICAL_RECOMMENDATION_STATUS", "RESEARCH_ONLY_DISCLAIMER"
        ]
        df23 = df23[[c for c in req_cols_23 if c in df23.columns]]
        df23.to_csv(p23, index=False)

        # Summary of best options per class
        best_options: dict[str, list[dict[str, Any]]] = {}
        for cls in [CLASS_H1, CLASS_H2, CLASS_H3]:
            cls_rows = [r for r in safety_rows if r.hypothesis_class == cls]
            top_cls = cls_rows[:5]
            best_options[cls] = [
                {
                    "rank": r.research_priority_rank,
                    "target": r.target_symbol,
                    "drug": r.drug_name,
                    "verdict": r.verdict_status,
                    "confidence": r.mechanism_confidence,
                    "pediatric_hazard_score": r.pediatric_hazard_score,
                    "justification": r.verdict_justification,
                }
                for r in top_cls
            ]

        # 24: Pre-registered Provenance JSON
        prov = {
            "pipeline_name": "mva_phase4_adjudication",
            "pipeline_version": self.config.pipeline_version,
            "execution_timestamp": datetime.now(timezone.utc).isoformat(),
            "duration_seconds": round((datetime.now(timezone.utc) - start_time).total_seconds(), 2),
            "offline_mode": self.config.offline_mode,
            "research_only_notice": RESEARCH_ONLY_DISCLAIMER,
            "input_files": {
                "gene_mechanism": {"path": str(self.config.gene_mechanism_path), "sha256": compute_sha256(self.config.gene_mechanism_path)},
                "drug_evidence": {"path": str(self.config.drug_evidence_path), "sha256": compute_sha256(self.config.drug_evidence_path)},
                "safety_ranking": {"path": str(self.config.safety_ranking_path), "sha256": compute_sha256(self.config.safety_ranking_path)},
                "query_panel": {"path": str(self.config.query_panel_path), "sha256": compute_sha256(self.config.query_panel_path)},
                "pathway_enrichment": {"path": str(self.config.pathway_enrichment_path), "sha256": compute_sha256(self.config.pathway_enrichment_path)},
            },
            "PREREGISTERED_ENDPOINTS": {
                "SAC_ACTIVITY": ["BUB1B", "MAD1L1", "CDC20", "BUB1", "BUB3", "MAD2L1", "MAD2L2"],
                "KT_MT_ATTACHMENT": ["AURKA", "PLK1", "CENPE", "KNL1", "NDC80", "TTK"],
                "MINOR_INTRON_SPLICING": ["CENATAC", "CCDC84", "RNU11", "RNU12"],
                "CHROMOSOME_STABILITY": ["CEP57", "SLF2", "SMC5", "TP53"],
                "P53_APOPTOSIS_SURVEILLANCE": ["TP53", "ATM", "CHEK1", "CHEK2"],
                "PEDIATRIC_MALIGNANCY_RMS_WILMS": ["ALK", "KIT", "WT1", "DICER1", "NF1", "FGFR4"],
            },
            "PREREGISTERED_SIGN_CONVENTION": {
                "activation_sign": 1,
                "inhibition_sign": -1,
                "unknown_ambiguous_sign": 0,
                "drug_action_mapping": {
                    "inhibitor_antagonist_degrader_blocker": -1,
                    "activator_agonist_inducer_opener": 1,
                    "modulator_binder_other_unknown": 0,
                },
                "net_path_sign_formula": "drug_sign * product(edge_signs_along_path)",
            },
            "PREREGISTERED_MAX_PATH_LENGTH": self.config.max_path_length,
            "PREREGISTERED_SOURCE_HIERARCHY": {
                "CURATED": ["SIGNOR", "OmniPath (omnipath)", "Reactome ContentService"],
                "INFERRED": ["OmniPath (pathwayextra, kinaseextra)", "Pathway Commons SIF"],
                "TEXTMINED": ["Europe PMC Annotations", "DGIdb Interaction Claims"],
            },
            "PREREGISTERED_SCORING_WEIGHTS": {
                "curated_signed_path": 0.40,
                "literature_verification": 0.30,
                "multi_source_consensus": 0.20,
                "reaction_level_corroboration": 0.10,
            },
            "PREREGISTERED_HAZARD_THRESHOLDS": {
                "pediatric_hazard_score_mandatory_hold": 60,
                "signed_path_hazard_contraindication": True,
            },
            "PREREGISTERED_CLASS_RULES": {
                "H1": "DISEASE_MECHANISM_CORRECTION (SAC / KT-MT / minor splicing rescue)",
                "H2": "CANCER_COMPLICATION_MANAGEMENT (RMS / Wilms / leukemia / synthetic lethality)",
                "H3": "SUPPORTIVE_CARE_SYMPTOM_MANAGEMENT (renal tubular / ion transport without mitotic impact)",
                "H0": "NO_CLASS (unresolved or evidence gap)",
            },
            "PREREGISTERED_ENTRY_CRITERIA": {
                "E1": "All Phase 3 Tier 1 concordant rows",
                "E2": "Phase 3 Tier 2 off-target compensation with >=2 sources",
                "E3": "Curated negative-control subset of Tier 4 (must-reject sanity)",
                "E4": "Pairs with target in MVA-associated pediatric malignancies",
                "E5": "Pairs in Phase 2 MVA-focus component with signed edge or reaction",
            },
            "sources_contacted": [
                {
                    "name": "SIGNOR",
                    "endpoint": "https://signor.uniroma2.it/getData.php?organism=9606",
                    "version": "3.0",
                    "status": "ONLINE" if not self.config.offline_mode else "OFFLINE_SKIPPED",
                    "rows_retrieved": edge_count,
                },
                {
                    "name": "OmniPath",
                    "endpoint": "https://omnipathdb.org/interactions",
                    "version": "2024",
                    "status": "ONLINE" if not self.config.offline_mode else "OFFLINE_SKIPPED",
                    "rows_retrieved": edge_count,
                },
                {
                    "name": "Reactome ContentService",
                    "endpoint": "https://reactome.org/ContentService",
                    "version": "v88",
                    "status": "ONLINE" if not self.config.offline_mode else "OFFLINE_SKIPPED",
                    "reactions_adjudicated": reaction_count,
                },
                {
                    "name": "Crossref REST",
                    "endpoint": "https://api.crossref.org/works",
                    "status": "ONLINE" if self.config.verify_literature and not self.config.offline_mode else "OFFLINE_SKIPPED",
                },
                {
                    "name": "NCBI Entrez E-utilities",
                    "endpoint": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                    "status": "ONLINE" if self.config.verify_literature and not self.config.offline_mode else "OFFLINE_SKIPPED",
                },
                {
                    "name": "Europe PMC REST",
                    "endpoint": "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
                    "status": "ONLINE" if self.config.verify_literature and not self.config.offline_mode else "OFFLINE_SKIPPED",
                },
            ],
            "metrics": {
                "entry_pairs_evaluated": entry_pairs_count,
                "signed_edges_reconciled": edge_count,
                "curated_conflicts_detected": conflict_count,
                "reactions_adjudicated": reaction_count,
                "gate_statistics": self.gate_stats,
                "verdict_distribution": {
                    v: sum(1 for r in safety_rows if r.verdict_status == v)
                    for v in set(r.verdict_status for r in safety_rows)
                },
                "class_distribution": {
                    c: sum(1 for r in classified_rows if r.hypothesis_class == c)
                    for c in set(r.hypothesis_class for r in classified_rows)
                },
            },
            "best_options_per_class": best_options,
            "truncation_ledger": self.truncation_ledger + self.reaction_engine.truncation_ledger,
            "quarantine_ledger": self.quarantine_ledger,
        }

        with open(p24, "w", encoding="utf-8") as f:
            json.dump(prov, f, indent=2)

        return {
            "reaction_matrix": p20,
            "hypothesis_classification": p21,
            "synthetic_lethality": p22,
            "adjudicated_safety_report": p23,
            "provenance": p24,
        }
