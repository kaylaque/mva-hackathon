"""Phase 4 Signed Pathway Regulation Engine.

Integrates curated and inferred directional signed interaction resources:
- SIGNOR (CausalTAB TSV human 9606)
- OmniPath (curated 'omnipath' dataset + optional 'pathwayextra'/'kinaseextra' inferred)
- Conflict resolution across curated databases (contradictory signs marked AMBIGUOUS and excluded)
- Graph traversal for signed path search (length <= 2) with net flux / sign math:
  ACT * ACT = +1, ACT * INH = -1, INH * INH = +1, AMBIGUOUS/UNKNOWN = 0.

Strict research exploration: no clinical claims.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import json
import logging
from pathlib import Path
import re
import time
from typing import Any, Mapping

import requests

logger = logging.getLogger("mva_hackathon.signed_regulation")

RESEARCH_ONLY_DISCLAIMER = (
    "RESEARCH EXPLORATION ONLY — STRICTLY NOT CLINICAL OR TREATMENT ADVICE. "
    "Hypotheses generated computationally from public databases for mechanistic inquiry."
)

SIGN_ACTIVATION = 1
SIGN_INHIBITION = -1
SIGN_UNKNOWN = 0


@dataclass
class SignedRegulationConfig:
    """Configuration for Signed Pathway Regulation queries and parsing."""

    signor_url: str = "https://signor.uniroma2.it/getData.php?organism=9606"
    omnipath_url: str = "https://omnipathdb.org/interactions"
    enable_omnipath_inferred: bool = False
    max_path_length: int = 2
    user_agent: str = "mva_hackathon_phase4/0.4.0 (research; mailto:mva_research@example.org)"
    request_timeout: int = 30
    max_retries: int = 3
    retry_backoff: float = 1.0
    polite_delay: float = 0.35
    offline_mode: bool = False
    cache_dir: Path | None = None


@dataclass
class SignedEdge:
    """A directed signed interaction edge between two biological entities."""

    source_gene: str
    target_gene: str
    effect: str  # ACTIVATION, INHIBITION, AMBIGUOUS, UNKNOWN
    sign: int  # +1, -1, 0
    mechanism: str
    direct: bool
    evidence_class: str  # CURATED, INFERRED, TEXTMINED
    source_db: str  # SIGNOR, OmniPath, Reactome
    source_db_version: str
    evidence_pmids: list[str] = field(default_factory=list)
    score: float = 1.0
    conflict_flag: bool = False
    source_uniprot: str = ""
    target_uniprot: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class SignedPath:
    """A directed signed path from a target gene to an endpoint gene."""

    start_node: str
    end_node: str
    edges: list[SignedEdge]
    path_length: int
    net_sign: int
    is_curated: bool
    evidence_pmids: list[str] = field(default_factory=list)
    sources: list[str] = field(default_factory=list)


def parse_signor_effect(effect_str: str) -> tuple[str, int]:
    """Map SIGNOR causal statements to standardized effect and sign (+1, -1, 0)."""
    clean = effect_str.strip().lower()
    if clean in {
        "up-regulates",
        "up-regulates activity",
        "up-regulates quantity",
        "up-regulates quantity by expression",
        "up-regulates quantity by stabilization",
        "stabilizes",
        "activates",
    }:
        return "ACTIVATION", SIGN_ACTIVATION
    if clean in {
        "down-regulates",
        "down-regulates activity",
        "down-regulates quantity",
        "down-regulates quantity by destabilization",
        "down-regulates quantity by repression",
        "destabilizes",
        "inhibits",
    }:
        return "INHIBITION", SIGN_INHIBITION
    return "UNKNOWN", SIGN_UNKNOWN


def parse_omnipath_sign(
    is_stim: bool, is_inhib: bool, cons_stim: bool | None = None, cons_inhib: bool | None = None
) -> tuple[str, int]:
    """Determine effect and sign from OmniPath stimulation/inhibition flags."""
    stim = bool(is_stim or cons_stim)
    inhib = bool(is_inhib or cons_inhib)
    if stim and not inhib:
        return "ACTIVATION", SIGN_ACTIVATION
    if inhib and not stim:
        return "INHIBITION", SIGN_INHIBITION
    if stim and inhib:
        return "AMBIGUOUS", SIGN_UNKNOWN
    return "UNKNOWN", SIGN_UNKNOWN


class SignorClient:
    """Client for SIGNOR Signaling Network Open Resource CausalTAB dataset."""

    def __init__(self, config: SignedRegulationConfig) -> None:
        self.config = config
        self.version = "SIGNOR 3.0"

    def fetch_edges(self, session: requests.Session | None = None) -> list[SignedEdge]:
        """Fetch and parse human causal interactions from SIGNOR."""
        if self.config.offline_mode:
            logger.info("Offline mode active; skipping live SIGNOR download")
            return []

        sess = session or requests.Session()
        sess.headers.update({"User-Agent": self.config.user_agent})

        edges: list[SignedEdge] = []
        try:
            time.sleep(self.config.polite_delay)
            resp = sess.get(self.config.signor_url, timeout=self.config.request_timeout)
            resp.raise_for_status()
            text = resp.text
            edges = self.parse_tsv(text)
            logger.info("Parsed %d signed edges from SIGNOR", len(edges))
        except Exception as exc:
            logger.warning("SIGNOR download/parse failed: %s", exc)
            if self.config.offline_mode:
                return []
            raise
        return edges

    def parse_tsv(self, tsv_text: str) -> list[SignedEdge]:
        """Parse SIGNOR CausalTAB formatted TSV text into SignedEdge records."""
        edges: list[SignedEdge] = []
        for line in tsv_text.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 10:
                continue

            entity_a = parts[0].strip()
            id_a = parts[2].strip() if len(parts) > 2 else ""
            entity_b = parts[4].strip() if len(parts) > 4 else ""
            id_b = parts[6].strip() if len(parts) > 6 else ""
            effect_raw = parts[8].strip() if len(parts) > 8 else ""
            mechanism = parts[9].strip() if len(parts) > 9 else ""
            pmid_raw = parts[21].strip() if len(parts) > 21 else ""
            direct_raw = parts[22].strip() if len(parts) > 22 else "t"
            score_raw = parts[27].strip() if len(parts) > 27 else "1.0"

            if not entity_a or not entity_b:
                continue

            effect, sign = parse_signor_effect(effect_raw)
            if sign == SIGN_UNKNOWN and effect == "UNKNOWN":
                continue

            direct = direct_raw.lower() == "t"
            pmids = [p.strip() for p in pmid_raw.split(";") if p.strip() and p.strip().isdigit()]
            try:
                score = float(score_raw) if score_raw else 1.0
            except ValueError:
                score = 1.0

            edge = SignedEdge(
                source_gene=entity_a,
                target_gene=entity_b,
                effect=effect,
                sign=sign,
                mechanism=mechanism or "causal_regulation",
                direct=direct,
                evidence_class="CURATED",
                source_db="SIGNOR",
                source_db_version=self.version,
                evidence_pmids=pmids,
                score=score,
                conflict_flag=False,
                source_uniprot=id_a,
                target_uniprot=id_b,
            )
            edges.append(edge)
        return edges


class OmniPathClient:
    """Client for OmniPath causal interactions service."""

    def __init__(self, config: SignedRegulationConfig) -> None:
        self.config = config
        self.version = "OmniPath 2024"

    def fetch_edges(
        self,
        partner_genes: list[str] | None = None,
        session: requests.Session | None = None,
    ) -> list[SignedEdge]:
        """Fetch signed interactions from OmniPath."""
        if self.config.offline_mode:
            logger.info("Offline mode active; skipping live OmniPath download")
            return []

        sess = session or requests.Session()
        sess.headers.update({"User-Agent": self.config.user_agent})

        datasets = "omnipath"
        if self.config.enable_omnipath_inferred:
            datasets = "omnipath,pathwayextra,kinaseextra"

        params: dict[str, Any] = {
            "genesymbols": "yes",
            "datasets": datasets,
            "organisms": "9606",
            "fields": "sources,references,curation_effort",
        }
        if partner_genes:
            params["partners"] = ",".join(sorted(set(partner_genes)))

        edges: list[SignedEdge] = []
        try:
            time.sleep(self.config.polite_delay)
            resp = sess.get(
                self.config.omnipath_url,
                params=params,
                timeout=self.config.request_timeout,
            )
            resp.raise_for_status()
            edges = self.parse_tsv(resp.text)
            logger.info("Parsed %d signed edges from OmniPath", len(edges))
        except Exception as exc:
            logger.warning("OmniPath download/parse failed: %s", exc)
            if self.config.offline_mode:
                return []
            raise
        return edges

    def parse_tsv(self, tsv_text: str) -> list[SignedEdge]:
        """Parse OmniPath interactions TSV text into SignedEdge records."""
        edges: list[SignedEdge] = []
        lines = tsv_text.splitlines()
        if not lines:
            return edges

        header = [h.strip() for h in lines[0].split("\t")]
        col_map = {name: idx for idx, name in enumerate(header)}

        src_sym_idx = col_map.get("source_genesymbol")
        tgt_sym_idx = col_map.get("target_genesymbol")
        src_up_idx = col_map.get("source")
        tgt_up_idx = col_map.get("target")
        is_stim_idx = col_map.get("is_stimulation")
        is_inhib_idx = col_map.get("is_inhibition")
        cons_stim_idx = col_map.get("consensus_stimulation")
        cons_inhib_idx = col_map.get("consensus_inhibition")
        sources_idx = col_map.get("sources")
        refs_idx = col_map.get("references")
        curation_idx = col_map.get("curation_effort")

        if src_sym_idx is None or tgt_sym_idx is None:
            return edges

        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) <= max(src_sym_idx, tgt_sym_idx):
                continue

            src_gene = parts[src_sym_idx].strip()
            tgt_gene = parts[tgt_sym_idx].strip()
            if not src_gene or not tgt_gene:
                continue

            src_up = parts[src_up_idx].strip() if src_up_idx is not None and len(parts) > src_up_idx else ""
            tgt_up = parts[tgt_up_idx].strip() if tgt_up_idx is not None and len(parts) > tgt_up_idx else ""

            is_stim = parts[is_stim_idx].strip() == "True" if is_stim_idx is not None and len(parts) > is_stim_idx else False
            is_inhib = parts[is_inhib_idx].strip() == "True" if is_inhib_idx is not None and len(parts) > is_inhib_idx else False
            cons_stim = parts[cons_stim_idx].strip() == "True" if cons_stim_idx is not None and len(parts) > cons_stim_idx else None
            cons_inhib = parts[cons_inhib_idx].strip() == "True" if cons_inhib_idx is not None and len(parts) > cons_inhib_idx else None

            effect, sign = parse_omnipath_sign(is_stim, is_inhib, cons_stim, cons_inhib)
            if sign == SIGN_UNKNOWN and effect != "AMBIGUOUS":
                continue

            # Extract PMIDs from references field (e.g. SPIKE:15531928;Lit-BM-17:16169070)
            pmids: list[str] = []
            if refs_idx is not None and len(parts) > refs_idx:
                raw_refs = parts[refs_idx].strip()
                for ref in raw_refs.split(";"):
                    if ":" in ref:
                        ref = ref.split(":")[-1]
                    ref = ref.strip()
                    if ref.isdigit():
                        pmids.append(ref)

            curation_effort = 1.0
            if curation_idx is not None and len(parts) > curation_idx:
                try:
                    curation_effort = float(parts[curation_idx].strip() or 1.0)
                except ValueError:
                    curation_effort = 1.0

            sources_val = parts[sources_idx].strip() if sources_idx is not None and len(parts) > sources_idx else "OmniPath"
            evidence_class = "CURATED"

            edge = SignedEdge(
                source_gene=src_gene,
                target_gene=tgt_gene,
                effect=effect,
                sign=sign,
                mechanism="activity_flow",
                direct=True,
                evidence_class=evidence_class,
                source_db="OmniPath",
                source_db_version=self.version,
                evidence_pmids=list(sorted(set(pmids))),
                score=curation_effort,
                conflict_flag=False,
                source_uniprot=src_up,
                target_uniprot=tgt_up,
            )
            edges.append(edge)

        return edges


class SignedRegulationEngine:
    """Consolidates signed regulatory edges, resolves conflicts, and searches signed paths."""

    def __init__(self, config: SignedRegulationConfig) -> None:
        self.config = config
        self.signor_client = SignorClient(config)
        self.omnipath_client = OmniPathClient(config)
        self.curated_edges: list[SignedEdge] = []
        self.conflict_records: list[dict[str, Any]] = []
        self.edge_map: dict[tuple[str, str], list[SignedEdge]] = {}
        self.consolidated_edges: dict[tuple[str, str], SignedEdge] = {}

    def ingest_and_reconcile(
        self,
        partner_genes: list[str] | None = None,
        custom_edges: list[SignedEdge] | None = None,
        session: requests.Session | None = None,
    ) -> tuple[int, int]:
        """Ingest edges from SIGNOR and OmniPath, resolve cross-database conflicts.

        Returns (total_reconciled_edges, conflict_count).
        """
        raw_edges: list[SignedEdge] = []
        if custom_edges:
            raw_edges.extend(custom_edges)
        else:
            raw_signor = self.signor_client.fetch_edges(session=session)
            raw_omni = self.omnipath_client.fetch_edges(partner_genes=partner_genes, session=session)
            raw_edges.extend(raw_signor)
            raw_edges.extend(raw_omni)

        self.curated_edges = raw_edges

        # Group by ordered pair (source_gene, target_gene)
        pair_groups: dict[tuple[str, str], list[SignedEdge]] = {}
        for edge in raw_edges:
            pair = (edge.source_gene.upper(), edge.target_gene.upper())
            pair_groups.setdefault(pair, []).append(edge)

        self.edge_map = pair_groups
        self.consolidated_edges = {}
        self.conflict_records = []

        for pair, edges in pair_groups.items():
            signs = {e.sign for e in edges if e.evidence_class == "CURATED" and e.sign != SIGN_UNKNOWN}
            if SIGN_ACTIVATION in signs and SIGN_INHIBITION in signs:
                # Direct contradiction across curated databases
                conflict_dbs = list({e.source_db for e in edges})
                self.conflict_records.append({
                    "source_gene": pair[0],
                    "target_gene": pair[1],
                    "databases": conflict_dbs,
                    "signs": list(signs),
                    "pmids": list(sorted({p for e in edges for p in e.evidence_pmids})),
                })
                # Reconcile as AMBIGUOUS with sign 0
                all_pmids = list(sorted({p for e in edges for p in e.evidence_pmids}))
                primary = edges[0]
                reconciled = SignedEdge(
                    source_gene=pair[0],
                    target_gene=pair[1],
                    effect="AMBIGUOUS",
                    sign=SIGN_UNKNOWN,
                    mechanism="conflicting_curation",
                    direct=any(e.direct for e in edges),
                    evidence_class="CURATED",
                    source_db=";".join(sorted(set(conflict_dbs))),
                    source_db_version="reconciled",
                    evidence_pmids=all_pmids,
                    score=0.0,
                    conflict_flag=True,
                    source_uniprot=primary.source_uniprot,
                    target_uniprot=primary.target_uniprot,
                )
                self.consolidated_edges[pair] = reconciled
            else:
                # Concordant signs
                chosen_sign = SIGN_ACTIVATION if SIGN_ACTIVATION in signs else (
                    SIGN_INHIBITION if SIGN_INHIBITION in signs else SIGN_UNKNOWN
                )
                effect = "ACTIVATION" if chosen_sign == SIGN_ACTIVATION else (
                    "INHIBITION" if chosen_sign == SIGN_INHIBITION else "UNKNOWN"
                )
                all_dbs = ";".join(sorted({e.source_db for e in edges}))
                all_pmids = list(sorted({p for e in edges for p in e.evidence_pmids}))
                primary = edges[0]
                reconciled = SignedEdge(
                    source_gene=pair[0],
                    target_gene=pair[1],
                    effect=effect,
                    sign=chosen_sign,
                    mechanism=primary.mechanism,
                    direct=any(e.direct for e in edges),
                    evidence_class="CURATED" if any(e.evidence_class == "CURATED" for e in edges) else "INFERRED",
                    source_db=all_dbs,
                    source_db_version=primary.source_db_version,
                    evidence_pmids=all_pmids,
                    score=max(e.score for e in edges),
                    conflict_flag=False,
                    source_uniprot=primary.source_uniprot,
                    target_uniprot=primary.target_uniprot,
                )
                self.consolidated_edges[pair] = reconciled

        logger.info(
            "Consolidated %d distinct directed pairs with %d conflicts",
            len(self.consolidated_edges),
            len(self.conflict_records),
        )
        return len(self.consolidated_edges), len(self.conflict_records)

    def find_signed_paths(
        self,
        start_gene: str,
        target_endpoints: set[str],
        max_length: int = 2,
    ) -> list[SignedPath]:
        """Find directed signed paths of length <= max_length from start_gene to endpoints.

        Excludes edges with sign == 0 or conflict_flag == True.
        Net sign math: product of edge signs along the path.
        """
        start = start_gene.upper()
        endpoints = {ep.upper() for ep in target_endpoints}
        paths: list[SignedPath] = []

        # Length 1 paths: start -> endpoint
        for ep in endpoints:
            if ep == start:
                continue
            pair = (start, ep)
            edge = self.consolidated_edges.get(pair)
            if edge and edge.sign != SIGN_UNKNOWN and not edge.conflict_flag:
                paths.append(
                    SignedPath(
                        start_node=start,
                        end_node=ep,
                        edges=[edge],
                        path_length=1,
                        net_sign=edge.sign,
                        is_curated=(edge.evidence_class == "CURATED"),
                        evidence_pmids=list(edge.evidence_pmids),
                        sources=[edge.source_db],
                    )
                )

        if max_length >= 2:
            # Length 2 paths: start -> intermediate -> endpoint
            # First find all outgoing edges from start
            for (s, mid), edge1 in self.consolidated_edges.items():
                if s != start or edge1.sign == SIGN_UNKNOWN or edge1.conflict_flag:
                    continue
                if mid == start or mid in endpoints:
                    # direct edges to endpoint already handled
                    continue

                for ep in endpoints:
                    if ep == mid or ep == start:
                        continue
                    pair2 = (mid, ep)
                    edge2 = self.consolidated_edges.get(pair2)
                    if edge2 and edge2.sign != SIGN_UNKNOWN and not edge2.conflict_flag:
                        net_sign = edge1.sign * edge2.sign
                        combined_pmids = list(sorted(set(edge1.evidence_pmids + edge2.evidence_pmids)))
                        combined_sources = list(sorted(set([edge1.source_db, edge2.source_db])))
                        is_curated = (edge1.evidence_class == "CURATED" and edge2.evidence_class == "CURATED")
                        paths.append(
                            SignedPath(
                                start_node=start,
                                end_node=ep,
                                edges=[edge1, edge2],
                                path_length=2,
                                net_sign=net_sign,
                                is_curated=is_curated,
                                evidence_pmids=combined_pmids,
                                sources=combined_sources,
                            )
                        )

        return paths
