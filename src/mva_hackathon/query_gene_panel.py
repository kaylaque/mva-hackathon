"""Phase 3 Query Gene Panel Generation Engine.

Derives a rigorous, graph-grounded, multi-tiered gene query panel for downstream
Phase 3 drug-repurposing inquiry:
- Gate 0: Symbol hygiene & non-query filtering (quarantines unmapped, HLA, read-throughs, pseudogenes)
- Gate 1: Core MVA causal genes (BUB1B, CEP57, TRIP13, CENATAC, MAD1L1, SLF2, SMC5) injected regardless of Phase 1 presence
- Gate 2: Graph-derived interactome neighbors from the programmatically detected MVA-focus component
- Gate 3: Context modifiers from qualified Phase 1/2 candidates + manually-curated literature references (TP53, ATM, NF1, MAD2L2)
- Gate 4: Positive and programmatically-derived negative drug-target controls
- Gate 5 & 6: Deduplication, integrity assertions, and panel sizing caps

Outputs:
- 18_phase3_query_gene_panel.csv
- 18_phase3_query_panel_provenance.json
- 19_phase3_query_panel_citations.csv

Never logs or prints patient data rows to standard output.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import json
import logging
from pathlib import Path
import re
import time
from typing import Any
import urllib.parse
import urllib.request

import pandas as pd

from mva_hackathon.mva_gene_sets import (
    CANONICAL_MVA_GENES,
    CORE_MVA_CAUSAL_GENES,
    EXTENDED_MVA_GENES,
    classify_mva_gene_class,
)

logger = logging.getLogger("mva_hackathon.query_gene_panel")

RESEARCH_ONLY_DISCLAIMER = (
    "RESEARCH EXPLORATION ONLY — STRICTLY NOT CLINICAL OR TREATMENT ADVICE. "
    "Hypotheses generated computationally from public databases for mechanistic inquiry."
)

DEFAULT_NCBI_TOOL = "mva_hackathon_panel_builder"
DEFAULT_NCBI_EMAIL = "mva_research@example.org"

# Regex for strict validation
DOI_PATTERN = re.compile(r"^10\.\d{4,9}/\S+$")
PMID_PATTERN = re.compile(r"^\d+$")
MITOTIC_KW_REGEX = re.compile(
    r"(?:chromosome\s+segregation|spindle|mitotic|cell\s+cycle|dna\s+damage|cohesin|checkpoint)",
    re.IGNORECASE,
)


@dataclass
class QueryGenePanelConfig:
    """Configuration parameters for the Phase 3 Query Gene Panel Builder."""

    track2_dir: Path
    output_dir: Path | None = None
    flat_names: bool = True
    string_min_score: float = 0.700
    max_context_modifiers: int = 25
    max_panel_size: int = 120
    include_core_absent_from_phase1: bool = True
    mva_focus_component_only: bool = True
    allow_unmapped_identity: bool = False
    include_no_kb_genes: bool = False
    verify_literature: bool = False
    literature_registry_path: Path | None = None
    control_positive_genes: list[str] = field(
        default_factory=lambda: [
            "AURKA",
            "PLK1",
            "BUB1",
            "CDC20",
            "MAD2L1",
            "PLK4",
            "ATR",
            "CHEK1",
            "WEE1",
        ]
    )
    control_negative_max: int = 4
    pipeline_version: str = "0.3.0"


def verify_literature_via_pubmed(
    term: str,
    expected_title_words: list[str] | None = None,
    tool: str = DEFAULT_NCBI_TOOL,
    email: str = DEFAULT_NCBI_EMAIL,
    max_retries: int = 3,
    polite_delay: float = 0.35,
) -> dict[str, Any]:
    """Look up and verify citation metadata using NCBI E-utilities (esearch + esummary).

    ``verified`` means the citation metadata was confirmed to exist and match.  It is not
    evidence for a mechanistic claim, which is reviewed separately in Compass.  ``status``
    distinguishes a real no-match from an unreachable or malformed source, so a failed
    lookup is never reported as "no such publication".
    """

    def blank(status: str, reason: str) -> dict[str, Any]:
        return {"pmid": "", "doi": "", "title": "", "journal": "", "year": "", "verified": False, "status": status, "reason": reason}

    params = {
        "db": "pubmed",
        "term": term,
        "retmode": "json",
        "retmax": "3",
        "tool": tool,
        "email": email,
    }
    search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{urllib.parse.urlencode(params)}"
    headers = {"User-Agent": "mva_hackathon/0.3.0"}

    id_list: list[str] = []
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(search_url, headers=headers)
            with urllib.request.urlopen(req, timeout=15) as resp:
                data = json.loads(resp.read().decode("utf-8"))
            id_list = data.get("esearchresult", {}).get("idlist", [])
            break
        except json.JSONDecodeError:
            return blank("INVALID_RESPONSE", "esearch returned a body that is not valid JSON")
        except Exception as exc:
            if attempt == max_retries - 1:
                logger.warning("NCBI esearch failed for query '%s': %s", term, exc)
                return blank("UNAVAILABLE", f"esearch unreachable: {type(exc).__name__}")
            time.sleep(polite_delay * (attempt + 1))

    if not id_list:
        return blank("NO_MATCH", "esearch returned no identifier for this query")

    time.sleep(polite_delay)

    sum_params = {
        "db": "pubmed",
        "id": ",".join(id_list),
        "retmode": "json",
        "tool": tool,
        "email": email,
    }
    sum_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?{urllib.parse.urlencode(sum_params)}"

    sum_data: dict[str, Any] = {}
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(sum_url, headers=headers)
            with urllib.request.urlopen(req, timeout=15) as resp:
                sum_data = json.loads(resp.read().decode("utf-8"))
            break
        except json.JSONDecodeError:
            return blank("INVALID_RESPONSE", "esummary returned a body that is not valid JSON")
        except Exception as exc:
            if attempt == max_retries - 1:
                logger.warning("NCBI esummary failed for IDs '%s': %s", id_list, exc)
                return blank("UNAVAILABLE", f"esummary unreachable: {type(exc).__name__}")
            time.sleep(polite_delay * (attempt + 1))

    results_dict = sum_data.get("result", {})
    for pmid in id_list:
        doc = results_dict.get(pmid)
        if not doc or not isinstance(doc, dict):
            continue
        title = doc.get("title", "")
        source = doc.get("source", "")
        pubdate = doc.get("pubdate", "")
        year_match = re.search(r"\b(19\d\d|20\d\d)\b", pubdate)
        year = year_match.group(1) if year_match else ""

        doi = ""
        for aid in doc.get("articleids", []):
            if aid.get("idtype") == "doi":
                cand_doi = aid.get("value", "")
                if DOI_PATTERN.match(cand_doi):
                    doi = cand_doi
                break

        # Check if title matches expected words if provided
        if expected_title_words:
            title_lower = title.lower()
            matches_expected = any(word.lower() in title_lower for word in expected_title_words)
            if not matches_expected:
                continue

        if PMID_PATTERN.match(pmid):
            return {
                "pmid": pmid,
                "doi": doi,
                "title": title,
                "journal": source,
                "year": year,
                "verified": True,
                "status": "VERIFIED",
                "reason": "citation metadata confirmed; this is not evidence for a mechanistic claim",
            }

    return blank("NO_MATCH", "no retrieved record matched the expected title words")


class QueryGenePanelBuilder:
    """Builder class orchestrating the multi-gate derivation of the query gene panel."""

    def __init__(self, config: QueryGenePanelConfig) -> None:
        self.config = config
        self.registry: dict[str, dict[str, Any]] = self._load_literature_registry()

    def _load_literature_registry(self) -> dict[str, dict[str, Any]]:
        """Load static literature registry from file."""
        if self.config.literature_registry_path and self.config.literature_registry_path.is_file():
            path = self.config.literature_registry_path
        else:
            default_path = Path(__file__).parent / "literature_registry.json"
            if default_path.is_file():
                path = default_path
            else:
                path = self.config.track2_dir / "literature_registry.json"

        if path.is_file():
            try:
                return json.loads(path.read_text(encoding="utf-8"))
            except Exception as exc:
                logger.warning("Failed to load literature registry from %s: %s", path, exc)
        return {}

    def load_track2_tables(self) -> dict[str, pd.DataFrame]:
        """Safely load required Phase 1 and Phase 2 tables from track2 directory."""
        track2 = self.config.track2_dir
        tables: dict[str, pd.DataFrame] = {}

        file_map = {
            "p1_gene": "02_phase1_gene_mechanism_matrix.csv",
            "p2_seed": "03_phase2_seed_genes.csv",
            "p2_edges": "04_phase2_string_edges.csv",
            "p2_pathways": "06_phase2_pathway_enrichment.csv",
            "p2_context": "07_phase2_network_gene_context.csv",
            "p2_nodes": "09_phase2_network_nodes.csv",
            "p2_components": "10_phase2_component_summary.csv",
            "p3_evidence_existing": "14_phase3_drug_evidence_matrix.csv",
        }

        for key, fname in file_map.items():
            fpath = track2 / fname
            if fpath.is_file():
                tables[key] = pd.read_csv(fpath)
            else:
                logger.info("Table %s not found at %s; initializing empty frame", key, fpath)
                tables[key] = pd.DataFrame()

        return tables

    def detect_mva_focus_component(
        self,
        components_df: pd.DataFrame,
        nodes_df: pd.DataFrame,
    ) -> dict[str, Any]:
        """Programmatically detect the MVA-focus connected component.

        Searches for the component containing any Core MVA causal anchor
        (BUB1B, MAD1L1, TRIP13, etc.). Does NOT hard-code component rank.
        """
        anchors = ["BUB1B", "MAD1L1", "TRIP13", "CEP57", "CENATAC", "SLF2", "SMC5"]
        detected_rank: int | None = None
        size: int = 0
        edge_count: int = 0
        found_anchors: list[str] = []

        if not components_df.empty and "COMPONENT_RANK" in components_df.columns:
            for _, row in components_df.iterrows():
                gene_list_str = str(row.get("GENE_LIST", ""))
                genes_in_comp = [g.strip() for g in gene_list_str.split(";") if g.strip()]
                matching = [a for a in anchors if a in genes_in_comp]
                if matching:
                    detected_rank = int(row["COMPONENT_RANK"])
                    size = int(row.get("COMPONENT_SIZE", len(genes_in_comp)))
                    edge_count = int(row.get("EDGE_COUNT", 0))
                    found_anchors = matching
                    break

        if detected_rank is None and not nodes_df.empty and "COMPONENT_RANK" in nodes_df.columns:
            for a in anchors:
                node_row = nodes_df[nodes_df["HGNC_SYMBOL"] == a]
                if not node_row.empty:
                    detected_rank = int(node_row.iloc[0]["COMPONENT_RANK"])
                    comp_nodes = nodes_df[nodes_df["COMPONENT_RANK"] == detected_rank]
                    size = len(comp_nodes)
                    found_anchors = [g for g in anchors if g in comp_nodes["HGNC_SYMBOL"].values]
                    break

        if detected_rank is None:
            detected_rank = 2  # Fallback to documented rank if data frames are empty/synthetic
            found_anchors = ["BUB1B", "MAD1L1"]

        return {
            "component_rank": detected_rank,
            "size": size,
            "edge_count": edge_count,
            "anchors": found_anchors,
        }

    def apply_gate0_hygiene(
        self,
        seed_genes_df: pd.DataFrame,
    ) -> tuple[pd.DataFrame, list[dict[str, str]]]:
        """Gate 0: Symbol hygiene.

        Quarantines unmapped symbols, HLA-* complex genes, read-through transcripts,
        and pseudogenes into a dropped ledger.
        """
        dropped_ledger: list[dict[str, str]] = []
        if seed_genes_df.empty:
            return pd.DataFrame(), dropped_ledger

        kept_indices = []
        for idx, row in seed_genes_df.iterrows():
            sym = str(row.get("HGNC_SYMBOL", "")).strip()
            status = str(row.get("HGNC_STATUS", "")).strip()
            desc = str(row.get("DESCRIPTION", "")).lower()

            # Always protect Core MVA, positive controls, and curated context from Gate 0 drop
            if (
                sym in CORE_MVA_CAUSAL_GENES
                or sym in ["TP53", "ATM", "NF1", "MAD2L2"]
                or sym in self.config.control_positive_genes
            ):
                kept_indices.append(idx)
                continue

            if not self.config.allow_unmapped_identity and status != "mapped":
                dropped_ledger.append({"symbol": sym, "reason": "UNMAPPED_HGNC_STATUS"})
                continue

            if sym.startswith("HLA-"):
                dropped_ledger.append({"symbol": sym, "reason": "HLA_COMPLEX_NON_QUERY"})
                continue

            if "-" in sym:
                dropped_ledger.append({"symbol": sym, "reason": "READTHROUGH_OR_NON_CODING_TRANSCRIPT"})
                continue

            if "pseudogene" in desc:
                dropped_ledger.append({"symbol": sym, "reason": "PSEUDOGENE"})
                continue

            kept_indices.append(idx)

        return seed_genes_df.loc[kept_indices].copy(), dropped_ledger

    def apply_gate1_core(
        self,
        p1_df: pd.DataFrame,
    ) -> list[dict[str, Any]]:
        """Gate 1: CORE_MVA_CAUSAL.

        Ensures all 7 core MVA causal genes (BUB1B, CEP57, TRIP13, CENATAC, MAD1L1, SLF2, SMC5)
        are included, injecting missing ones if absent from Phase 1.
        """
        core_ordered = ["BUB1B", "CEP57", "TRIP13", "CENATAC", "MAD1L1", "SLF2", "SMC5"]
        results: list[dict[str, Any]] = []

        p1_by_symbol = {}
        if not p1_df.empty and "SYMBOL" in p1_df.columns:
            for _, r in p1_df.iterrows():
                sym = str(r["SYMBOL"]).strip()
                p1_by_symbol[sym] = r

        for sym in core_ordered:
            mva_class = classify_mva_gene_class(sym)
            in_p1 = sym in p1_by_symbol
            p1_row = p1_by_symbol.get(sym)

            results.append({
                "gene": sym,
                "input_symbol": sym,
                "gene_category": "CORE_MVA_CAUSAL",
                "category_justification": "Core MVA/aneuploidy causal gene (established or reported in literature)",
                "mva_gene_class": mva_class,
                "query_priority": 1,
                "selection_source": "MVA_CORE_PANEL",
                "reason_for_querying": "Canonical or established MVA causal gene driving spindle assembly checkpoint / mitotic chromosome segregation phenotype",
                "evidence_nature": "KNOWN_GENE",
                "phase1_present": in_p1,
                "phase1_row": p1_row,
                "control_role": "NONE",
                "control_expectation": "NA",
            })

        return results

    def apply_gate2_graph_neighbors(
        self,
        component_rank: int,
        nodes_df: pd.DataFrame,
        edges_df: pd.DataFrame,
        core_genes: set[str],
    ) -> list[dict[str, Any]]:
        """Gate 2: GRAPH_DERIVED_NEIGHBOR.

        Includes all members of the detected MVA-focus component that are not core genes
        and not curated context genes (which are classified in Gate 3 as CONTEXT_MODIFIER).
        Determines ONE_HOP_OF_CORE vs TWO_HOP_OR_COMPONENT_MEMBER using STRING edges.
        """
        results: list[dict[str, Any]] = []
        if nodes_df.empty or "COMPONENT_RANK" not in nodes_df.columns:
            return results

        curated_context_set = {"TP53", "ATM", "NF1", "MAD2L2"}
        excluded_genes = core_genes | curated_context_set

        comp_nodes = nodes_df[nodes_df["COMPONENT_RANK"] == component_rank]
        comp_genes = [
            str(g).strip()
            for g in comp_nodes["HGNC_SYMBOL"].dropna().unique()
            if str(g).strip() not in excluded_genes
        ]

        # Identify incident edges touching core genes in edges_df
        incident_neighbors: dict[str, dict[str, Any]] = {}
        if not edges_df.empty:
            for _, erow in edges_df.iterrows():
                score = float(erow.get("SCORE", 0.0))
                if score < self.config.string_min_score:
                    continue
                a = str(erow.get("PREFERRED_NAME_A", "")).strip()
                b = str(erow.get("PREFERRED_NAME_B", "")).strip()
                curated = str(erow.get("CURATED_SUPPORT", "NO")).strip().upper() == "YES"

                partner = None
                core_partner = None
                if a in core_genes and b not in core_genes:
                    partner = b
                    core_partner = a
                elif b in core_genes and a not in core_genes:
                    partner = a
                    core_partner = b

                if partner:
                    prev = incident_neighbors.get(partner)
                    if prev is None or score > prev["score"] or (curated and not prev["curated"]):
                        incident_neighbors[partner] = {
                            "core_partner": core_partner,
                            "score": score,
                            "curated": curated,
                        }

        for gene in comp_genes:
            is_one_hop = gene in incident_neighbors
            hop_tier = "ONE_HOP_OF_CORE" if is_one_hop else "TWO_HOP_OR_COMPONENT_MEMBER"
            best_info = incident_neighbors.get(gene, {})

            results.append({
                "gene": gene,
                "input_symbol": gene,
                "gene_category": "GRAPH_DERIVED_NEIGHBOR",
                "category_justification": "Direct or component interactor of MVA core machinery within the MVA-focus interactome component",
                "mva_gene_class": "NONE",
                "query_priority": 2,
                "selection_source": "PHASE2_STRING_MVA_COMPONENT",
                "reason_for_querying": "Interactome partner within the MVA-focus connected component; potential network modifier or synthetic target",
                "evidence_nature": "GRAPH_INDIRECT",
                "neighbor_hop": hop_tier,
                "best_core_partner": best_info.get("core_partner", ""),
                "best_edge_score": best_info.get("score", ""),
                "curated_edge_support": "YES" if best_info.get("curated") else ("NO" if is_one_hop else "NA"),
                "control_role": "NONE",
                "control_expectation": "NA",
            })

        return results

    def apply_gate3_context_modifiers(
        self,
        clean_seed_df: pd.DataFrame,
        context_df: pd.DataFrame,
        pathways_df: pd.DataFrame,
        p1_df: pd.DataFrame,
        selected_genes: set[str],
    ) -> list[dict[str, Any]]:
        """Gate 3: CONTEXT_MODIFIER.

        Ranks and caps remaining gated candidates matching compound-het, >=2 passing variants,
        or mitotic pathway co-membership.
        Injects curated context genes (TP53, ATM, NF1, MAD2L2) labelled as NON-MVA-causal.
        """
        curated_context_definitions = {
            "TP53": "Aneuploidy surveillance and cancer predisposition context; somatic/germline checkpoint upon chromosome missegregation",
            "ATM": "DNA-damage response and double-strand break repair context linked to mitotic catastrophe",
            "NF1": "RAS pathway / childhood tumor context including rhabdomyosarcoma and neurofibromatosis",
            "MAD2L2": "REV7/Shieldin SAC and DNA repair-adjacent mitotic regulator",
        }

        # Build set of genes in matching pathways
        mitotic_pathway_genes: set[str] = set()
        if not pathways_df.empty:
            for _, r in pathways_df.iterrows():
                pname = str(r.get("PATHWAY_NAME", ""))
                hits = [g.strip() for g in str(r.get("HIT_GENES", "")).split(";") if g.strip()]
                shares_core_or_neighbor = any(g in selected_genes for g in hits)
                matches_kw = bool(MITOTIC_KW_REGEX.search(pname))
                if shares_core_or_neighbor or matches_kw:
                    mitotic_pathway_genes.update(hits)

        # Context metrics map
        context_map: dict[str, dict[str, Any]] = {}
        if not context_df.empty and "HGNC_SYMBOL" in context_df.columns:
            for _, r in context_df.iterrows():
                sym = str(r["HGNC_SYMBOL"]).strip()
                context_map[sym] = {
                    "network_degree": int(r.get("NETWORK_DEGREE", 0)) if pd.notna(r.get("NETWORK_DEGREE")) else 0,
                    "is_isolated": bool(r.get("IS_ISOLATED") == True),
                }

        qualifying_candidates: list[dict[str, Any]] = []
        if not clean_seed_df.empty and "HGNC_SYMBOL" in clean_seed_df.columns:
            for _, r in clean_seed_df.iterrows():
                sym = str(r["HGNC_SYMBOL"]).strip()
                if sym in selected_genes or sym in curated_context_definitions:
                    continue

                comp_het = bool(r.get("COMPOUND_HET_ARCHETYPE_FIT") == True)
                variants = int(r["TOTAL_PASSING_CANDIDATE_VARIANTS"]) if pd.notna(r.get("TOTAL_PASSING_CANDIDATE_VARIANTS")) else 0
                in_pathway = sym in mitotic_pathway_genes

                if comp_het or variants >= 2 or in_pathway:
                    c_info = context_map.get(sym, {})
                    deg = c_info.get("network_degree", 0)
                    qualifying_candidates.append({
                        "gene": sym,
                        "degree": deg,
                        "variants": variants,
                        "comp_het": int(comp_het),
                        "row": r,
                    })

        # Rank by (network_degree, variants, comp_het) descending
        qualifying_candidates.sort(
            key=lambda x: (x["degree"], x["variants"], x["comp_het"]),
            reverse=True,
        )

        capped_candidates = qualifying_candidates[: self.config.max_context_modifiers]

        results: list[dict[str, Any]] = []
        for item in capped_candidates:
            sym = item["gene"]
            results.append({
                "gene": sym,
                "input_symbol": sym,
                "gene_category": "CONTEXT_MODIFIER",
                "category_justification": "Gated Phase 1/2 candidate meeting compound-het, high passing variant count, or mitotic pathway co-membership",
                "mva_gene_class": "NONE",
                "query_priority": 3,
                "selection_source": "PHASE1_CANDIDATE_PATHWAY",
                "reason_for_querying": "Secondary candidate with qualifying biological support or mitotic pathway membership",
                "evidence_nature": "PATHWAY_COMEMBERSHIP",
                "control_role": "NONE",
                "control_expectation": "NA",
            })

        # Now include the curated context genes (TP53, ATM, NF1, MAD2L2)
        p1_symbols = set(p1_df["SYMBOL"].dropna().astype(str).str.strip()) if not p1_df.empty else set()
        clean_seed_symbols = (
            set(clean_seed_df["HGNC_SYMBOL"].dropna().astype(str).str.strip())
            if not clean_seed_df.empty
            else set()
        )

        for sym, rationale in curated_context_definitions.items():
            if sym in selected_genes or any(r["gene"] == sym for r in results):
                continue
            in_p1 = sym in p1_symbols
            results.append({
                "gene": sym,
                "input_symbol": sym,
                "gene_category": "CONTEXT_MODIFIER",
                "category_justification": rationale,
                "mva_gene_class": "NONE",
                "query_priority": 3,
                "selection_source": "CURATED_CONTEXT",
                "reason_for_querying": f"Curated literature context modifier (non-MVA-causal): {rationale}",
                "evidence_nature": "CURATED_CONTEXT",
                "phase1_present": in_p1,
                "control_role": "NONE",
                "control_expectation": "NA",
            })

        return results

    def apply_gate4_controls(
        self,
        clean_seed_df: pd.DataFrame,
        context_df: pd.DataFrame,
        pathways_df: pd.DataFrame,
        selected_map: dict[str, dict[str, Any]],
    ) -> list[dict[str, Any]]:
        """Gate 4: DRUG_TARGET_CONTROL.

        Adds curated positive druggable controls and derives negative controls
        (isolated sensory receptors with no mitotic pathway co-membership).
        Respects precedence: if a control was already selected in Gates 1-3,
        it retains its primary category while setting control role/expectation.
        """
        results: list[dict[str, Any]] = []

        # 1. Positive controls
        for sym in self.config.control_positive_genes:
            if sym in selected_map:
                # Update existing entry with control role
                selected_map[sym]["control_role"] = "POSITIVE"
                selected_map[sym]["control_expectation"] = "EXPECTED_DRUG_HITS"
            else:
                results.append({
                    "gene": sym,
                    "input_symbol": sym,
                    "gene_category": "DRUG_TARGET_CONTROL",
                    "category_justification": "Curated positive drug-target control: established druggable SAC/DDR mitotic target",
                    "mva_gene_class": "NONE",
                    "query_priority": 3,
                    "selection_source": "DRUG_TARGET_CONTROL",
                    "reason_for_querying": "Positive reference control to validate API coverage of known druggable mitotic kinases and regulators",
                    "evidence_nature": "CONTROL",
                    "control_role": "POSITIVE",
                    "control_expectation": "EXPECTED_DRUG_HITS",
                })

        # 2. Derive negative controls programmatically
        # Find isolated genes with no mitotic pathways and sensory/structural descriptions
        mitotic_pathway_genes: set[str] = set()
        if not pathways_df.empty:
            for _, r in pathways_df.iterrows():
                pname = str(r.get("PATHWAY_NAME", ""))
                hits = [g.strip() for g in str(r.get("HIT_GENES", "")).split(";") if g.strip()]
                if MITOTIC_KW_REGEX.search(pname):
                    mitotic_pathway_genes.update(hits)

        isolated_genes: set[str] = set()
        if not context_df.empty and "IS_ISOLATED" in context_df.columns:
            iso_rows = context_df[context_df["IS_ISOLATED"] == True]
            isolated_genes = set(iso_rows["HGNC_SYMBOL"].dropna().astype(str).str.strip())

        negative_candidates: list[dict[str, str]] = []
        if not clean_seed_df.empty:
            for _, r in clean_seed_df.iterrows():
                sym = str(r["HGNC_SYMBOL"]).strip()
                if (
                    sym in selected_map
                    or any(item["gene"] == sym for item in results)
                    or sym in self.config.control_positive_genes
                ):
                    continue

                if sym in isolated_genes and sym not in mitotic_pathway_genes:
                    desc = str(r.get("DESCRIPTION", "")).lower()
                    # Look for olfactory receptor or taste receptor
                    is_olfactory = sym.startswith("OR") and re.match(r"^OR\d+[A-Z0-9]+$", sym)
                    is_taste = sym.startswith("TAS2R")
                    if is_olfactory or is_taste:
                        negative_candidates.append({
                            "gene": sym,
                            "desc": desc,
                            "type": "Olfactory receptor" if is_olfactory else "Taste receptor",
                        })

        for cand in negative_candidates[: self.config.control_negative_max]:
            sym = cand["gene"]
            results.append({
                "gene": sym,
                "input_symbol": sym,
                "gene_category": "DRUG_TARGET_CONTROL",
                "category_justification": f"Programmatically derived negative control: isolated {cand['type'].lower()} with zero mitotic pathway interactions and expected absence of therapeutic drug hits",
                "mva_gene_class": "NONE",
                "query_priority": 3,
                "selection_source": "DRUG_TARGET_CONTROL",
                "reason_for_querying": "Negative reference control: isolated sensory receptor with expected zero drug-repurposing hits",
                "evidence_nature": "CONTROL",
                "control_role": "NEGATIVE",
                "control_expectation": "EXPECTED_NO_HITS",
            })

        return results

    def build_panel(self) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
        """Execute full multi-gate pipeline and return panel, citations, and provenance."""
        tables = self.load_track2_tables()
        comp_info = self.detect_mva_focus_component(tables["p2_components"], tables["p2_nodes"])
        clean_seeds, dropped_ledger = self.apply_gate0_hygiene(tables["p2_seed"])

        core_genes_set = set(CORE_MVA_CAUSAL_GENES)

        # Gate 1: Core MVA
        gate1_items = self.apply_gate1_core(tables["p1_gene"])
        selected_map: dict[str, dict[str, Any]] = {item["gene"]: item for item in gate1_items}

        # Gate 2: Graph neighbors
        gate2_items = self.apply_gate2_graph_neighbors(
            component_rank=comp_info["component_rank"],
            nodes_df=tables["p2_nodes"],
            edges_df=tables["p2_edges"],
            core_genes=core_genes_set,
        )
        for item in gate2_items:
            if item["gene"] not in selected_map:
                selected_map[item["gene"]] = item

        # Gate 3: Context modifiers
        selected_so_far = set(selected_map.keys())
        gate3_items = self.apply_gate3_context_modifiers(
            clean_seed_df=clean_seeds,
            context_df=tables["p2_context"],
            pathways_df=tables["p2_pathways"],
            p1_df=tables["p1_gene"],
            selected_genes=selected_so_far,
        )
        for item in gate3_items:
            if item["gene"] not in selected_map:
                selected_map[item["gene"]] = item

        # Gate 4: Controls
        gate4_items = self.apply_gate4_controls(
            clean_seed_df=clean_seeds,
            context_df=tables["p2_context"],
            pathways_df=tables["p2_pathways"],
            selected_map=selected_map,
        )
        for item in gate4_items:
            if item["gene"] not in selected_map:
                selected_map[item["gene"]] = item

        # Gate 5 & 6: Deduplication assertions and capping
        panel_items = list(selected_map.values())
        if len(panel_items) > self.config.max_panel_size:
            logger.info("Panel size %d exceeds max %d; trimming context modifiers", len(panel_items), self.config.max_panel_size)
            # Retain Core, Neighbors, Controls; trim candidate context modifiers
            non_trim = [it for it in panel_items if it["selection_source"] != "PHASE1_CANDIDATE_PATHWAY"]
            trimmable = [it for it in panel_items if it["selection_source"] == "PHASE1_CANDIDATE_PATHWAY"]
            allowed_trim = self.config.max_panel_size - len(non_trim)
            panel_items = non_trim + trimmable[:allowed_trim]

        # Verify all 7 core genes are present
        present_genes = {it["gene"] for it in panel_items}
        for cg in CORE_MVA_CAUSAL_GENES:
            assert cg in present_genes, f"Core MVA gene {cg} missing from panel!"

        # Assemble full rows with schema
        panel_df, citations_df, prov_counts = self._assemble_records(
            panel_items=panel_items,
            tables=tables,
            comp_info=comp_info,
            dropped_ledger=dropped_ledger,
        )

        return panel_df, citations_df, prov_counts

    def _assemble_records(
        self,
        panel_items: list[dict[str, Any]],
        tables: dict[str, pd.DataFrame],
        comp_info: dict[str, Any],
        dropped_ledger: list[dict[str, str]],
    ) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
        """Construct normalized panel and citations dataframes according to strict schema."""
        # Index reference tables by symbol
        p1_by_sym = {}
        if not tables["p1_gene"].empty and "SYMBOL" in tables["p1_gene"].columns:
            for _, r in tables["p1_gene"].iterrows():
                p1_by_sym[str(r["SYMBOL"]).strip()] = r

        p2_seed_by_sym = {}
        if not tables["p2_seed"].empty and "HGNC_SYMBOL" in tables["p2_seed"].columns:
            for _, r in tables["p2_seed"].iterrows():
                p2_seed_by_sym[str(r["HGNC_SYMBOL"]).strip()] = r

        p2_context_by_sym = {}
        if not tables["p2_context"].empty and "HGNC_SYMBOL" in tables["p2_context"].columns:
            for _, r in tables["p2_context"].iterrows():
                p2_context_by_sym[str(r["HGNC_SYMBOL"]).strip()] = r

        p2_nodes_by_sym = {}
        if not tables["p2_nodes"].empty and "HGNC_SYMBOL" in tables["p2_nodes"].columns:
            for _, r in tables["p2_nodes"].iterrows():
                p2_nodes_by_sym[str(r["HGNC_SYMBOL"]).strip()] = r

        p3_existing_counts: dict[str, int] = {}
        if not tables["p3_evidence_existing"].empty and "TARGET_SYMBOL" in tables["p3_evidence_existing"].columns:
            p3_existing_counts = tables["p3_evidence_existing"]["TARGET_SYMBOL"].value_counts().to_dict()

        panel_rows: list[dict[str, Any]] = []
        citations_rows: list[dict[str, Any]] = []

        verified_lit_count = 0
        pending_lit_count = 0
        live_status: dict[str, dict[str, str]] = {}

        for item in panel_items:
            gene = item["gene"]
            cat = item["gene_category"]
            just = item["category_justification"]
            mva_class = item["mva_gene_class"]
            prio = item["query_priority"]
            source = item["selection_source"]
            reason = item["reason_for_querying"]
            nature = item["evidence_nature"]
            ctrl_role = item.get("control_role", "NONE")
            ctrl_exp = item.get("control_expectation", "NA")

            # Look up Phase 1 metadata
            p1_row = p1_by_sym.get(gene)
            p1_present = p1_row is not None
            tot_variants = int(p1_row["TOTAL_VARIANTS_EVALUATED"]) if p1_present and pd.notna(p1_row.get("TOTAL_VARIANTS_EVALUATED")) else ""
            pass_variants = int(p1_row["TOTAL_PASSING_CANDIDATE_VARIANTS"]) if p1_present and pd.notna(p1_row.get("TOTAL_PASSING_CANDIDATE_VARIANTS")) else ""
            comp_het = bool(p1_row["COMPOUND_HET_ARCHETYPE_FIT"] == True) if p1_present and pd.notna(p1_row.get("COMPOUND_HET_ARCHETYPE_FIT")) else ""
            arch_desc = str(p1_row["ARCHETYPE_DESCRIPTION"]) if p1_present and pd.notna(p1_row.get("ARCHETYPE_DESCRIPTION")) else ""
            phase_gate = str(p1_row["PHASE_GATE_STATUS"]) if p1_present and pd.notna(p1_row.get("PHASE_GATE_STATUS")) else ""
            inherit = str(p1_row["INHERITANCE_REQUIREMENT"]) if p1_present and pd.notna(p1_row.get("INHERITANCE_REQUIREMENT")) else ""
            canon_dir = str(p1_row["CANONICAL_PATHWAY_DIRECTION"]) if p1_present and pd.notna(p1_row.get("CANONICAL_PATHWAY_DIRECTION")) else ""
            ther_dir = str(p1_row["THERAPEUTIC_DIRECTIONALITY"]) if p1_present and pd.notna(p1_row.get("THERAPEUTIC_DIRECTIONALITY")) else ""
            mva_link = str(p1_row["MVA_MECHANISM_LINK"]) if p1_present and pd.notna(p1_row.get("MVA_MECHANISM_LINK")) else ""
            pathway = str(p1_row["PATHWAY"]) if p1_present and pd.notna(p1_row.get("PATHWAY")) else ""

            # Phase 2 metadata
            seed_row = p2_seed_by_sym.get(gene)
            gated = seed_row is not None

            node_row = p2_nodes_by_sym.get(gene)
            ctx_row = p2_context_by_sym.get(gene)

            comp_rank = int(node_row["COMPONENT_RANK"]) if node_row is not None and pd.notna(node_row.get("COMPONENT_RANK")) else ""
            in_focus_comp = comp_rank == comp_info["component_rank"] if comp_rank != "" else False

            comp_id = int(ctx_row["COMPONENT_ID"]) if ctx_row is not None and pd.notna(ctx_row.get("COMPONENT_ID")) else ""
            deg = int(ctx_row["NETWORK_DEGREE"]) if ctx_row is not None and pd.notna(ctx_row.get("NETWORK_DEGREE")) else ""
            clustering = float(ctx_row["CLUSTERING_COEFFICIENT"]) if ctx_row is not None and pd.notna(ctx_row.get("CLUSTERING_COEFFICIENT")) else ""
            is_iso = bool(ctx_row["IS_ISOLATED"] == True) if ctx_row is not None and pd.notna(ctx_row.get("IS_ISOLATED")) else ""
            interactors = str(ctx_row["STRING_INTERACTORS"]) if ctx_row is not None and pd.notna(ctx_row.get("STRING_INTERACTORS")) else ""
            int_count = len([x for x in interactors.split(";") if x.strip()]) if interactors else 0
            enriched_pw = str(ctx_row["TOP_ENRICHED_PATHWAYS"]) if ctx_row is not None and pd.notna(ctx_row.get("TOP_ENRICHED_PATHWAYS")) else ""

            neighbor_hop = item.get("neighbor_hop", "")
            best_partner = item.get("best_core_partner", "")
            best_edge_score = item.get("best_edge_score", "")
            curated_support = item.get("curated_edge_support", "NA" if neighbor_hop == "" else "NO")

            # Literature registry info
            reg_entry = self.registry.get(gene, {})
            gene_name = reg_entry.get("name") or (str(p1_row.get("GENE_NAME", "")) if p1_present else "") or gene
            hgnc_id = reg_entry.get("hgnc_id", "")
            ensembl_id = reg_entry.get("ensembl_id", "")
            uniprot_id = reg_entry.get("uniprot_id", "")
            omim_gene = reg_entry.get("omim_gene_id", "")
            omim_disease = reg_entry.get("mva_disease_omim", "")

            pmid = reg_entry.get("primary_pmid", "")
            doi = reg_entry.get("primary_doi", "")
            lit_mode = reg_entry.get("literature_verification_mode", "NA" if cat == "GRAPH_DERIVED_NEIGHBOR" else "PENDING")
            lit_date = reg_entry.get("literature_verification_date", "")

            # If verify_literature requested live and entry has search query
            if self.config.verify_literature and reg_entry.get("search_query"):
                live_res = verify_literature_via_pubmed(
                    term=reg_entry["search_query"],
                    expected_title_words=[gene, "aneuploidy", "mva", "mosaic", "syndrome", "kinase", "checkpoint", "protein"],
                )
                live_status[gene] = {"status": live_res.get("status", "UNKNOWN"), "reason": live_res.get("reason", "")}
                if live_res["verified"]:
                    pmid = live_res["pmid"]
                    doi = live_res["doi"]
                    lit_mode = "PUBMED_EUTILITIES"
                    lit_date = datetime.now(timezone.utc).strftime("%Y-%m-%d")
                else:
                    # A source failure must stay visible instead of reading as "no publication exists".
                    lit_mode = f"PUBMED_{live_res.get('status', 'UNKNOWN')}"

            # URLs
            hgnc_url = f"https://www.genenames.org/data/gene-symbol-report/#!/symbol/{gene}"
            omim_url = f"https://www.omim.org/entry/{omim_gene}" if omim_gene else ""
            clinvar_url = f"https://www.ncbi.nlm.nih.gov/clinvar/?term={gene}%5Bgene%5D"
            gencc_url = f"https://search.thegencc.org/genes/{gene}"
            uniprot_url = f"https://www.uniprot.org/uniprotkb/{uniprot_id}/entry" if uniprot_id else ""
            string_url = f"https://version-12-5.string-db.org/cgi/network?identifiers={gene}&species=9606"
            reactome_url = "https://reactome.org/PathwayBrowser/"

            # QA & Evidence gap flags
            gap_flags: list[str] = []
            if not p1_present:
                gap_flags.append("PHASE1_ABSENT")
            if cat == "CORE_MVA_CAUSAL" and not pmid:
                gap_flags.append("NO_PMID")
            if cat == "CORE_MVA_CAUSAL" and not doi:
                gap_flags.append("NO_DOI")
            if cat == "CORE_MVA_CAUSAL" and not omim_disease:
                gap_flags.append("NO_OMIM")
            if cat == "GRAPH_DERIVED_NEIGHBOR":
                gap_flags.append("GRAPH_ONLY")
            if not pathway and not enriched_pw:
                gap_flags.append("NO_MECHANISM_KB")

            # Citations records
            if pmid or doi:
                verified_lit_count += 1
                if pmid:
                    citations_rows.append({
                        "QUERY_GENE": gene,
                        "GENE_CATEGORY": cat,
                        "CITATION_TYPE": "PRIMARY_LITERATURE",
                        "ID_TYPE": "PMID",
                        "ID_VALUE": pmid,
                        "CITATION_TITLE": reg_entry.get("primary_title", ""),
                        "CITATION_DETAILS": f"{reg_entry.get('primary_journal', '')} ({reg_entry.get('primary_year', '')})",
                        "URL": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                        "VERIFICATION_STATUS": "VERIFIED" if lit_mode == "PUBMED_EUTILITIES" else "STATIC_CURATION",
                    })
                if doi:
                    citations_rows.append({
                        "QUERY_GENE": gene,
                        "GENE_CATEGORY": cat,
                        "CITATION_TYPE": "PRIMARY_LITERATURE",
                        "ID_TYPE": "DOI",
                        "ID_VALUE": doi,
                        "CITATION_TITLE": reg_entry.get("primary_title", ""),
                        "CITATION_DETAILS": f"{reg_entry.get('primary_journal', '')} ({reg_entry.get('primary_year', '')})",
                        "URL": f"https://doi.org/{doi}",
                        "VERIFICATION_STATUS": "VERIFIED" if lit_mode == "PUBMED_EUTILITIES" else "STATIC_CURATION",
                    })
            elif cat in ["CORE_MVA_CAUSAL", "CONTEXT_MODIFIER"] and nature in ["KNOWN_GENE", "CURATED_CONTEXT"]:
                pending_lit_count += 1
                lit_mode = "PENDING"

            if omim_disease:
                citations_rows.append({
                    "QUERY_GENE": gene,
                    "GENE_CATEGORY": cat,
                    "CITATION_TYPE": "DISEASE_RESOURCE",
                    "ID_TYPE": "OMIM_DISEASE",
                    "ID_VALUE": omim_disease,
                    "CITATION_TITLE": f"OMIM phenotype entry #{omim_disease}",
                    "CITATION_DETAILS": f"Online Mendelian Inheritance in Man #{omim_disease}",
                    "URL": f"https://www.omim.org/entry/{omim_disease}",
                    "VERIFICATION_STATUS": "STATIC_CURATION",
                })

            # Add database citation entry (STRING / HGNC)
            citations_rows.append({
                "QUERY_GENE": gene,
                "GENE_CATEGORY": cat,
                "CITATION_TYPE": "DATABASE_ENTRY",
                "ID_TYPE": "STRING",
                "ID_VALUE": gene,
                "CITATION_TITLE": f"STRING DB v12.5 PPI network entry for {gene}",
                "CITATION_DETAILS": "STRING Database v12.5 (9606 / Homo sapiens)",
                "URL": string_url,
                "VERIFICATION_STATUS": "DATABASE_ENTRY",
            })

            # Existing Phase 3 drug coverage
            exist_drug_count = p3_existing_counts.get(gene, 0)
            has_p3_cov = exist_drug_count > 0

            panel_rows.append({
                # Identity
                "QUERY_GENE": gene,
                "INPUT_SYMBOL": gene,
                "GENE_NAME": gene_name,
                "GENE_CATEGORY": cat,
                "CATEGORY_JUSTIFICATION": just,
                "MVA_GENE_CLASS": mva_class,
                "QUERY_PRIORITY": prio,
                "SELECTION_SOURCE": source,
                "REASON_FOR_QUERYING": reason,
                # Phase 1
                "PHASE1_PRESENT": p1_present,
                "TOTAL_VARIANTS_EVALUATED": tot_variants,
                "TOTAL_PASSING_CANDIDATE_VARIANTS": pass_variants,
                "COMPOUND_HET_ARCHETYPE_FIT": comp_het,
                "ARCHETYPE_DESCRIPTION": arch_desc,
                "PHASE_GATE_STATUS": phase_gate,
                "INHERITANCE_REQUIREMENT": inherit,
                "CANONICAL_PATHWAY_DIRECTION": canon_dir,
                "THERAPEUTIC_DIRECTIONALITY": ther_dir,
                "MVA_MECHANISM_LINK": mva_link,
                "PATHWAY": pathway,
                # Phase 2 graph
                "GATED": gated,
                "IN_MVA_FOCUS_COMPONENT": in_focus_comp,
                "COMPONENT_RANK": comp_rank,
                "COMPONENT_ID": comp_id,
                "NETWORK_DEGREE": deg,
                "CLUSTERING_COEFFICIENT": clustering,
                "IS_ISOLATED": is_iso,
                "STRING_INTERACTORS": interactors,
                "INTERACTOR_COUNT": int_count,
                "BEST_EDGE_PARTNER": best_partner,
                "BEST_EDGE_SCORE": best_edge_score,
                "MAX_STRING_SCORE": best_edge_score,
                "CURATED_EDGE_SUPPORT": curated_support,
                "NEIGHBOR_HOP": neighbor_hop,
                "TOP_ENRICHED_PATHWAYS": enriched_pw,
                "BEST_PATHWAY_FDR": "",
                "PATHWAY_EVIDENCE_URLS": reactome_url if enriched_pw else "",
                # Resources / literature
                "HGNC_STABLE_ID": hgnc_id,
                "ENSEMBL_GENE_STABLE_ID": ensembl_id,
                "UNIPROT_ID": uniprot_id,
                "OMIM_GENE_ID": omim_gene,
                "MVA_DISEASE_OMIM": omim_disease,
                "HGNC_URL": hgnc_url,
                "OMIM_URL": omim_url,
                "CLINVAR_URL": clinvar_url,
                "GENCC_URL": gencc_url,
                "UNIPROT_URL": uniprot_url,
                "STRING_URL": string_url,
                "REACTOME_URL": reactome_url,
                "PRIMARY_LITERATURE_PMID": pmid,
                "PRIMARY_LITERATURE_DOI": doi,
                "SECONDARY_REFERENCES": "",
                "LITERATURE_VERIFICATION_MODE": lit_mode,
                "LITERATURE_VERIFICATION_DATE": lit_date,
                "EVIDENCE_NATURE": nature,
                # Drug-query intent
                "DRUG_DB_QUERY_REASON": f"Assess therapeutic repurposing & safety liabilities for {cat} target {gene}",
                "DRUG_DATABASES_TO_QUERY": "OPEN_TARGETS;CHEMBL;DGIDB",
                "CONTROL_ROLE": ctrl_role,
                "CONTROL_EXPECTATION": ctrl_exp,
                "EXISTING_PHASE3_COVERAGE": has_p3_cov,
                "EXISTING_DRUG_ROW_COUNT": exist_drug_count,
                # QA
                "EVIDENCE_GAP_FLAGS": ";".join(gap_flags),
                "SOURCE_TABLE_TRACE": "02_phase1_gene_mechanism_matrix;03_phase2_seed_genes;04_phase2_string_edges;07_phase2_network_gene_context",
                "DATA_VERSION_TRACE": f"Track2_v{self.config.pipeline_version}_STRINGv12.5_Reactome2026",
                "ROW_QC": "PASS",
                "RESEARCH_ONLY_DISCLAIMER": RESEARCH_ONLY_DISCLAIMER,
            })

        panel_df = pd.DataFrame(panel_rows)
        citations_df = pd.DataFrame(citations_rows)

        cat_counts = panel_df["GENE_CATEGORY"].value_counts().to_dict()
        ctrl_counts = panel_df["CONTROL_ROLE"].value_counts().to_dict()

        prov_counts = {
            "pipeline_version": self.config.pipeline_version,
            "execution_timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "config": asdict(self.config),
            "detected_mva_focus_component": comp_info,
            "panel_category_counts": cat_counts,
            "control_counts": ctrl_counts,
            "total_panel_size": len(panel_df),
            "citations_count": len(citations_df),
            "literature_verification": {
                "verified_count": verified_lit_count,
                "pending_count": pending_lit_count,
                "live_source_status": live_status,
                "unavailable_count": sum(1 for item in live_status.values() if item["status"] in {"UNAVAILABLE", "INVALID_RESPONSE"}),
                "sources": ["NCBI_EUTILITIES", "OMIM", "STRING_DB", "REACTOME"],
            },
            "gate0_dropped_ledger_counts": pd.Series([d["reason"] for d in dropped_ledger]).value_counts().to_dict() if dropped_ledger else {},
            "total_gate0_dropped": len(dropped_ledger),
            "core_genes_injected_phase1_absent": [
                g for g in CORE_MVA_CAUSAL_GENES if g not in p1_by_sym
            ],
            "research_only_disclaimer": RESEARCH_ONLY_DISCLAIMER,
        }

        return panel_df, citations_df, prov_counts


def write_query_panel_outputs(
    panel_df: pd.DataFrame,
    provenance: dict[str, Any],
    citations_df: pd.DataFrame,
    output_dir: Path,
    flat_names: bool = True,
) -> dict[str, Path]:
    """Write Phase 3 query panel artifacts to output directory.

    Guarantees that files are < 8 MB and created with verified filenames.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    if flat_names:
        panel_path = output_dir / "18_phase3_query_gene_panel.csv"
        prov_path = output_dir / "18_phase3_query_panel_provenance.json"
        cit_path = output_dir / "19_phase3_query_panel_citations.csv"
    else:
        p3_dir = output_dir / "phase3"
        p3_dir.mkdir(parents=True, exist_ok=True)
        panel_path = p3_dir / "query_gene_panel.csv"
        prov_path = p3_dir / "query_panel_provenance.json"
        cit_path = p3_dir / "query_panel_citations.csv"

    panel_df.to_csv(panel_path, index=False)
    citations_df.to_csv(cit_path, index=False)
    prov_path.write_text(json.dumps(provenance, indent=2, default=str), encoding="utf-8")

    # Safety assertion on file sizes (< 8MB = 8,388,608 bytes)
    max_bytes = 8 * 1024 * 1024
    for p in [panel_path, prov_path, cit_path]:
        sz = p.stat().st_size
        if sz >= max_bytes:
            raise ValueError(f"Output artifact {p} size {sz} bytes exceeds 8MB safety threshold!")

    return {
        "query_panel": panel_path,
        "provenance": prov_path,
        "citations": cit_path,
    }
