"""Phase 4 Reaction-Level Adjudication Engine.

Queries Reactome ContentService (https://reactome.org/ContentService) for
enriched pathways from Phase 2 (06_phase2_pathway_enrichment.csv).
Identifies ReactionLikeEvents containing candidate drug targets or MVA panel genes,
captures reaction inputs, outputs, catalysts, and regulators (positive/negative),
and calculates reaction flux deltas based on drug perturbation.

Applies a hard cap via --reaction-limit with an explicit truncation ledger.
Research-only: strictly investigational.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import logging
from pathlib import Path
import time
from typing import Any

import requests

logger = logging.getLogger("mva_hackathon.reaction_adjudication")

RESEARCH_ONLY_DISCLAIMER = (
    "RESEARCH EXPLORATION ONLY — STRICTLY NOT CLINICAL OR TREATMENT ADVICE. "
    "Hypotheses generated computationally from public databases for mechanistic inquiry."
)


@dataclass
class ReactionAdjudicationConfig:
    """Configuration for Reactome ContentService reaction queries."""

    content_service_url: str = "https://reactome.org/ContentService"
    reaction_limit: int = 500
    user_agent: str = "mva_hackathon_phase4/0.4.0 (research; mailto:mva_research@example.org)"
    request_timeout: int = 25
    max_retries: int = 3
    retry_backoff: float = 1.0
    polite_delay: float = 0.2
    offline_mode: bool = False


@dataclass
class ReactionParticipant:
    """A physical participant or regulator in a biochemical reaction."""

    display_name: str
    role: str  # INPUT, OUTPUT, CATALYST, POSITIVE_REGULATOR, NEGATIVE_REGULATOR, REGULATOR
    gene_symbol: str = ""
    st_id: str = ""


@dataclass
class AdjudicatedReaction:
    """A Reactome reaction with participants, regulators, and flux delta evaluation."""

    st_id: str
    name: str
    pathway_id: str
    pathway_name: str
    species: str
    inputs: list[str] = field(default_factory=list)
    outputs: list[str] = field(default_factory=list)
    catalysts: list[str] = field(default_factory=list)
    positive_regulators: list[str] = field(default_factory=list)
    negative_regulators: list[str] = field(default_factory=list)
    panel_participants: list[str] = field(default_factory=list)
    target_roles: dict[str, list[str]] = field(default_factory=dict)
    reaction_flux_delta: int = 0
    research_only_disclaimer: str = RESEARCH_ONLY_DISCLAIMER

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def extract_gene_symbol(name: str) -> str:
    """Heuristic extraction of standard HGNC symbol from Reactome entity name."""
    if not name:
        return ""
    # E.g. "p-S15,S20-TP53 Tetramer [nucleoplasm]" -> "TP53"
    # E.g. "BUB1B [cytosol]" -> "BUB1B"
    # Remove bracketed compartments
    clean = name.split("[")[0].strip()
    # Remove phosphorylation/modification prefixes
    if "-" in clean and len(clean.split("-")[0]) <= 10 and any(c in clean.split("-")[0] for c in ["p-", "S", "T", "Y", "Ac-"]):
        clean = clean.split("-", 1)[1].strip()
    # Match standard alphanumeric gene symbol token
    tokens = clean.split()
    for tok in tokens:
        tok_clean = tok.strip(":,()[]{}")
        if tok_clean.isalnum() and tok_clean.isupper() and 2 <= len(tok_clean) <= 12:
            return tok_clean
    return tokens[0].strip(":,()[]{}") if tokens else ""


class ReactomeReactionClient:
    """Client for Reactome ContentService."""

    def __init__(self, config: ReactionAdjudicationConfig) -> None:
        self.config = config

    def fetch_contained_events(
        self, pathway_id: str, session: requests.Session | None = None
    ) -> list[dict[str, Any]]:
        """Fetch contained events for a given Reactome pathway."""
        if self.config.offline_mode:
            return []

        clean_id = pathway_id.replace("REAC:", "").strip()
        url = f"{self.config.content_service_url}/data/pathway/{clean_id}/containedEvents"
        sess = session or requests.Session()
        sess.headers.update({"User-Agent": self.config.user_agent})

        for attempt in range(1, self.config.max_retries + 1):
            try:
                if self.config.polite_delay > 0:
                    time.sleep(self.config.polite_delay)
                resp = sess.get(url, timeout=self.config.request_timeout)
                if resp.status_code == 404:
                    return []
                resp.raise_for_status()
                data = resp.json()
                if isinstance(data, list):
                    return [ev for ev in data if isinstance(ev, dict)]
                return []
            except Exception as exc:
                if attempt == self.config.max_retries:
                    logger.warning("Failed fetching events for pathway %s: %s", clean_id, exc)
                    return []
                time.sleep(self.config.retry_backoff * attempt)
        return []

    def fetch_event_query(
        self, st_id: str, session: requests.Session | None = None
    ) -> dict[str, Any] | None:
        """Fetch full record for a reaction event."""
        if self.config.offline_mode:
            return None

        url = f"{self.config.content_service_url}/data/query/{st_id}"
        sess = session or requests.Session()
        sess.headers.update({"User-Agent": self.config.user_agent})

        for attempt in range(1, self.config.max_retries + 1):
            try:
                if self.config.polite_delay > 0:
                    time.sleep(self.config.polite_delay)
                resp = sess.get(url, timeout=self.config.request_timeout)
                if resp.status_code == 404:
                    return None
                resp.raise_for_status()
                data = resp.json()
                if isinstance(data, dict):
                    return data
                return None
            except Exception as exc:
                if attempt == self.config.max_retries:
                    logger.warning("Failed fetching query for reaction %s: %s", st_id, exc)
                    return None
                time.sleep(self.config.retry_backoff * attempt)
        return None


class ReactionAdjudicationEngine:
    """Orchestrates reaction extraction, participant filtering, and flux delta evaluation."""

    def __init__(self, config: ReactionAdjudicationConfig) -> None:
        self.config = config
        self.client = ReactomeReactionClient(config)
        self.truncation_ledger: list[dict[str, Any]] = []

    def adjudicate_pathways(
        self,
        pathway_items: list[tuple[str, str]],  # [(pathway_id, pathway_name)]
        panel_genes: set[str],
        target_genes: set[str],
        session: requests.Session | None = None,
    ) -> list[AdjudicatedReaction]:
        """Extract and adjudicate reactions involving panel or target genes across pathways."""
        all_reactions: list[AdjudicatedReaction] = []
        seen_reaction_ids: set[str] = set()
        panel_upper = {g.upper() for g in panel_genes}
        target_upper = {g.upper() for g in target_genes}
        all_watch_genes = panel_upper | target_upper

        if self.config.offline_mode:
            logger.info("Reaction adjudication in offline mode; skipping live Reactome calls")
            return all_reactions

        sess = session or requests.Session()
        sess.headers.update({"User-Agent": self.config.user_agent})

        for pw_id, pw_name in pathway_items:
            if len(all_reactions) >= self.config.reaction_limit:
                self.truncation_ledger.append({
                    "reason": "reaction_limit_reached",
                    "pathway_id": pw_id,
                    "cap": self.config.reaction_limit,
                })
                break

            events = self.client.fetch_contained_events(pw_id, session=sess)
            for ev in events:
                if len(all_reactions) >= self.config.reaction_limit:
                    self.truncation_ledger.append({
                        "reason": "reaction_limit_reached",
                        "last_reaction_id": ev.get("stId"),
                        "cap": self.config.reaction_limit,
                    })
                    break

                st_id = ev.get("stId", "")
                schema_class = ev.get("schemaClass", "")
                if not st_id or st_id in seen_reaction_ids:
                    continue

                # Filter for reaction-like events
                if schema_class not in {"Reaction", "ReactionLikeEvent", "BlackBoxEvent", "Polymerisation", "Depolymerisation"}:
                    continue

                disp_name = ev.get("displayName", "")
                disp_upper = disp_name.upper()
                # Quick check if reaction name mentions any of the watch genes
                name_words = {w.strip(":,()[]{}").upper() for w in disp_name.split()}
                candidate_overlap = (name_words & all_watch_genes) or {g for g in all_watch_genes if f" {g} " in f" {disp_upper} " or f"-{g}" in disp_upper or f"{g}-" in disp_upper}
                if not candidate_overlap:
                    continue

                detail = self.client.fetch_event_query(st_id, session=sess)
                if not detail:
                    continue

                # Parse inputs, outputs, catalysts, regulators
                inputs = [item.get("displayName", "") for item in detail.get("input", []) if isinstance(item, dict)]
                outputs = [item.get("displayName", "") for item in detail.get("output", []) if isinstance(item, dict)]
                catalysts = [item.get("displayName", "") for item in detail.get("catalystActivity", []) if isinstance(item, dict)]
                regulators = [item for item in detail.get("regulatedBy", []) if isinstance(item, dict)]

                pos_regs: list[str] = []
                neg_regs: list[str] = []
                for reg in regulators:
                    reg_name = reg.get("displayName", "")
                    reg_class = reg.get("className", "") or reg.get("schemaClass", "")
                    if "Positive" in reg_class or "positive" in reg_name.lower():
                        pos_regs.append(reg_name)
                    else:
                        neg_regs.append(reg_name)

                # Collect all participant gene symbols
                all_participant_texts = inputs + outputs + catalysts + pos_regs + neg_regs
                participant_genes: set[str] = set()
                for p_text in all_participant_texts:
                    sym = extract_gene_symbol(p_text)
                    if sym:
                        participant_genes.add(sym.upper())

                relevant_genes = (participant_genes | candidate_overlap) & all_watch_genes
                if not relevant_genes:
                    continue

                # Identify participant roles for targets
                target_roles: dict[str, list[str]] = {}
                for t in target_upper:
                    roles: list[str] = []
                    if any(t in p.upper() for p in inputs):
                        roles.append("INPUT")
                    if any(t in p.upper() for p in outputs):
                        roles.append("OUTPUT")
                    if any(t in p.upper() for p in catalysts):
                        roles.append("CATALYST")
                    if any(t in p.upper() for p in pos_regs):
                        roles.append("POSITIVE_REGULATOR")
                    if any(t in p.upper() for p in neg_regs):
                        roles.append("NEGATIVE_REGULATOR")
                    if roles:
                        target_roles[t] = roles

                reaction = AdjudicatedReaction(
                    st_id=st_id,
                    name=disp_name,
                    pathway_id=pw_id,
                    pathway_name=pw_name,
                    species=detail.get("speciesName", "Homo sapiens"),
                    inputs=inputs,
                    outputs=outputs,
                    catalysts=catalysts,
                    positive_regulators=pos_regs,
                    negative_regulators=neg_regs,
                    panel_participants=list(sorted(relevant_genes)),
                    target_roles=target_roles,
                    reaction_flux_delta=0,
                )
                seen_reaction_ids.add(st_id)
                all_reactions.append(reaction)

        logger.info("Adjudicated %d relevant reactions across %d pathways", len(all_reactions), len(pathway_items))
        return all_reactions

    def compute_reaction_flux_delta(
        self,
        reaction: AdjudicatedReaction,
        target_symbol: str,
        drug_action_sign: int,  # +1 (activator) or -1 (inhibitor)
    ) -> int:
        """Compute the net predicted reaction flux delta under drug action.

        - If drug inhibits a negative regulator -> +1 (disinhibition / activation)
        - If drug inhibits a catalyst or positive regulator -> -1 (suppression)
        - If drug activates a positive regulator or catalyst -> +1 (stimulation)
        - If drug activates a negative regulator -> -1 (suppression)
        """
        roles = reaction.target_roles.get(target_symbol.upper(), [])
        if not roles:
            return 0

        flux_delta = 0
        if "NEGATIVE_REGULATOR" in roles:
            if drug_action_sign < 0:
                flux_delta += 1  # Inhibiting an inhibitor restores reaction flux
            elif drug_action_sign > 0:
                flux_delta -= 1  # Activating an inhibitor shuts down reaction flux
        if "POSITIVE_REGULATOR" in roles or "CATALYST" in roles:
            if drug_action_sign < 0:
                flux_delta -= 1  # Inhibiting an essential positive element dampens flux
            elif drug_action_sign > 0:
                flux_delta += 1  # Activating an essential positive element boosts flux

        return flux_delta
