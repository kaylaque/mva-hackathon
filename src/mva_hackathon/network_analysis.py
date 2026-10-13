"""Phase 2: Protein-Protein Interaction (PPI) and Pathway Network Context.

Consumes caller-supplied Phase 1 gene mechanism matrices (and optional variant matrices),
querying external biomedical APIs at runtime:
- g:Profiler (g:Convert): Official HGNC symbol mapping and normalization
- STRING Database v12.5: High-confidence physical and functional PPI networks (default score: 700)
- Reactome Analysis Service: Pathway overrepresentation and projection
- g:Profiler (g:GOSt): Multi-source functional enrichment (GO:BP, REAC, KEGG, HP, WP)

Provides:
- Graceful error handling (retries, timeouts, fallback handling for offline/mocked environments)
- Polite rate-limiting and User-Agent identification
- Zero network expansion by default (preserves candidate focus without diluting seed list)
- Configurable STRING interaction score channels and thresholds
- Pure-Python graph topology analysis (degree, clustering coefficient, connected components)
- Output tables: seed_genes, string_edges, network_summary, pathway_enrichment, network_gene_context, provenance
- Strict non-diagnostic disclaimer: graph topology does not infer directionality or drug suitability.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import json
import logging
from pathlib import Path
import time
from typing import Any, Mapping, Sequence

import pandas as pd
import requests

logger = logging.getLogger("mva_hackathon.network_analysis")

# Default API endpoints
DEFAULT_GPROFILER_CONVERT_URL = "https://biit.cs.ut.ee/gprofiler/api/convert/convert/"
DEFAULT_GPROFILER_GOST_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
# Version-pinned STRING endpoint so results are reproducible across releases.
DEFAULT_STRING_API_URL = "https://version-12-5.string-db.org/api/json/network"
DEFAULT_STRING_VERSION_URL = "https://version-12-5.string-db.org/api/json/version"
DEFAULT_REACTOME_ANALYSIS_URL = "https://reactome.org/AnalysisService/identifiers/projection"

DEFAULT_USER_AGENT = "mva-hackathon/0.2.0 (SageBio MVA Hackathon 2026; research exploration; non-diagnostic)"


@dataclass
class NetworkAnalysisConfig:
    """Configuration for Phase 2 PPI & Pathway Network Analysis."""

    # STRING Database parameters
    string_species: int = 9606  # Homo sapiens
    string_required_score: int = 700  # High confidence primary threshold (0-1000)
    string_additional_nodes: int = 0  # No network expansion by default (keep seeds only)
    string_network_flavor: str = "confidence"  # "confidence" or "evidence"

    # Score channel filtering (None or minimum threshold 0.0-1.0)
    # nscore: gene neighborhood, fscore: gene fusion, pscore: phylogenetic co-occurrence
    # ascore: coexpression, escore: experimental, dscore: database, tscore: textmining
    min_escore: float | None = None  # experimental
    min_dscore: float | None = None  # database curated
    min_tscore: float | None = None  # textmining
    min_ascore: float | None = None  # co-expression

    # Seed filtering: keep genes carrying at least N passing Phase 0 candidate
    # variants (0 disables filtering and uses the full input gene matrix).
    min_passing_variants: int = 0

    # STRING GET requests are batched so large seed lists stay within URL limits.
    string_batch_size: int = 300

    # Reactome parameters
    reactome_page_size: int = 50

    # g:Profiler parameters
    gprofiler_organism: str = "hsapiens"
    gprofiler_target: str = "HGNC"
    gprofiler_sources: list[str] = field(
        default_factory=lambda: ["GO:BP", "REAC", "KEGG", "HP", "WP"]
    )
    gprofiler_user_threshold: float = 0.05

    # HTTP client behavior
    request_timeout: float = 15.0
    max_retries: int = 3
    retry_backoff: float = 1.0
    polite_delay: float = 0.25  # seconds between HTTP calls
    user_agent: str = DEFAULT_USER_AGENT

    # Endpoint URLs
    string_api_url: str = DEFAULT_STRING_API_URL
    string_version_url: str = DEFAULT_STRING_VERSION_URL
    reactome_api_url: str = DEFAULT_REACTOME_ANALYSIS_URL
    gprofiler_convert_url: str = DEFAULT_GPROFILER_CONVERT_URL
    gprofiler_gost_url: str = DEFAULT_GPROFILER_GOST_URL


class HttpClient:
    """Polite, retry-aware HTTP client for external biomedical APIs."""

    def __init__(self, config: NetworkAnalysisConfig, session: requests.Session | None = None) -> None:
        self.config = config
        self.session = session or requests.Session()
        self.session.headers.update({"User-Agent": config.user_agent})

    def get(self, url: str, params: dict[str, Any] | None = None) -> requests.Response | None:
        """Execute GET request with exponential backoff retry."""
        return self._request_with_retry("GET", url, params=params)

    def post(
        self,
        url: str,
        json_data: Any | None = None,
        data: Any | None = None,
        headers: dict[str, str] | None = None,
        params: dict[str, Any] | None = None,
    ) -> requests.Response | None:
        """Execute POST request with exponential backoff retry."""
        return self._request_with_retry(
            "POST", url, json_data=json_data, data=data, headers=headers, params=params
        )

    def _request_with_retry(
        self,
        method: str,
        url: str,
        params: dict[str, Any] | None = None,
        json_data: Any | None = None,
        data: Any | None = None,
        headers: dict[str, str] | None = None,
    ) -> requests.Response | None:
        for attempt in range(1, self.config.max_retries + 1):
            if self.config.polite_delay > 0:
                time.sleep(self.config.polite_delay)
            try:
                resp = self.session.request(
                    method=method,
                    url=url,
                    params=params,
                    json=json_data,
                    data=data,
                    headers=headers,
                    timeout=self.config.request_timeout,
                )
                if resp.status_code in {200, 204}:
                    return resp
                if resp.status_code in {429, 500, 502, 503, 504} and attempt < self.config.max_retries:
                    sleep_time = self.config.retry_backoff * (2 ** (attempt - 1))
                    logger.warning(
                        "HTTP %d from %s (attempt %d/%d). Retrying in %.1fs...",
                        resp.status_code,
                        url,
                        attempt,
                        self.config.max_retries,
                        sleep_time,
                    )
                    time.sleep(sleep_time)
                    continue
                logger.error("HTTP error %d from %s: %s", resp.status_code, url, resp.text[:200])
                return resp
            except (requests.RequestException, TimeoutError) as exc:
                if attempt < self.config.max_retries:
                    sleep_time = self.config.retry_backoff * (2 ** (attempt - 1))
                    logger.warning(
                        "Network error %s contacting %s (attempt %d/%d). Retrying in %.1fs...",
                        exc,
                        url,
                        attempt,
                        self.config.max_retries,
                        sleep_time,
                    )
                    time.sleep(sleep_time)
                else:
                    logger.error("Failed to reach %s after %d attempts: %s", url, self.config.max_retries, exc)
                    return None
        return None


class NetworkContextEngine:
    """Downstream PPI & Pathway Network Analysis Engine for Phase 2."""

    def __init__(
        self,
        config: NetworkAnalysisConfig | None = None,
        http_client: HttpClient | None = None,
    ) -> None:
        self.config = config or NetworkAnalysisConfig()
        self.client = http_client or HttpClient(self.config)

    def run_analysis(
        self,
        gene_mechanism_df: pd.DataFrame,
        variant_matrix_df: pd.DataFrame | None = None,
    ) -> dict[str, Any]:
        """Execute end-to-end network context analysis on candidate genes.

        Returns a dictionary containing:
        - seed_genes_df: Table of seed genes with HGNC mapping and input status
        - string_edges_df: Table of PPI network edges with scores and channel breakdowns
        - network_summary_df: Single-row or multi-metric topological network statistics
        - pathway_enrichment_df: Combined Reactome and g:Profiler enrichment terms
        - network_gene_context_df: Gene-level synthesis combining Phase 1 mechanism and Phase 2 network metrics
        - provenance: Detailed metadata dictionary of versions, parameters, endpoints, and timestamps
        """
        start_time = datetime.now(timezone.utc).isoformat()

        # 1. Restrict seeds to gated candidates, then extract and sanitize symbols.
        #    A gene qualifies if it carries at least min_passing_variants passing
        #    Phase 0 candidates (0 disables filtering = full input gene matrix).
        if (
            self.config.min_passing_variants > 0
            and "TOTAL_PASSING_CANDIDATE_VARIANTS" in gene_mechanism_df.columns
        ):
            passing_numeric = pd.to_numeric(
                gene_mechanism_df["TOTAL_PASSING_CANDIDATE_VARIANTS"], errors="coerce"
            ).fillna(0)
            keep_mask = passing_numeric >= self.config.min_passing_variants
            filtered_df = gene_mechanism_df.loc[keep_mask].reset_index(drop=True)
            logger.info(
                "Seed filter: kept %d of %d gene rows with >= %d passing candidate variants",
                len(filtered_df),
                len(gene_mechanism_df),
                self.config.min_passing_variants,
            )
            gene_mechanism_df = filtered_df
        input_symbols = self._extract_symbols(gene_mechanism_df)
        logger.info("Extracted %d unique candidate symbols from input table", len(input_symbols))

        # 2. HGNC Symbol normalization via g:Profiler Convert API
        hgnc_mapping, gprofiler_convert_status = self.map_hgnc_symbols(input_symbols)

        # 3. Create Seed Genes DataFrame
        seed_genes_df = self._build_seed_genes_df(
            input_symbols, hgnc_mapping, gene_mechanism_df, variant_matrix_df
        )

        valid_symbols = [
            row["HGNC_SYMBOL"]
            for _, row in seed_genes_df.iterrows()
            if pd.notna(row["HGNC_SYMBOL"]) and row["HGNC_SYMBOL"] != ""
        ]

        # 4. Query STRING DB v12.5 PPI network
        string_version = self.get_string_version()
        string_edges_df, string_status = self.query_string_ppi(valid_symbols)

        # 5. Compute network topology
        network_summary_df, node_topologies = self.compute_network_topology(
            valid_symbols, string_edges_df
        )

        # 6. Query Reactome pathway enrichment
        reactome_df, reactome_status = self.query_reactome(valid_symbols)

        # 7. Query g:Profiler (g:GOSt) functional enrichment
        gprofiler_df, gprofiler_status = self.query_gprofiler(valid_symbols)

        # 8. Combine pathway enrichment
        pathway_enrichment_df = self._combine_enrichment(reactome_df, gprofiler_df)

        # 9. Synthesize Network Gene Context
        network_gene_context_df = self._build_gene_context_df(
            seed_genes_df, node_topologies, pathway_enrichment_df
        )

        # 10. Compile Provenance
        end_time = datetime.now(timezone.utc).isoformat()
        provenance = {
            "workflow_phase": "Phase 2: Protein-Protein Interaction (PPI) and Pathway Network Context",
            "timestamp_start": start_time,
            "timestamp_end": end_time,
            "config": asdict(self.config),
            "external_resources": {
                "string_db": {
                    "version": string_version,
                    "endpoint": self.config.string_api_url,
                    "species": self.config.string_species,
                    "required_score": self.config.string_required_score,
                    "additional_nodes": self.config.string_additional_nodes,
                    "network_flavor": self.config.string_network_flavor,
                    "status": string_status,
                },
                "reactome": {
                    "endpoint": self.config.reactome_api_url,
                    "status": reactome_status,
                },
                "gprofiler": {
                    "convert_endpoint": self.config.gprofiler_convert_url,
                    "gost_endpoint": self.config.gprofiler_gost_url,
                    "organism": self.config.gprofiler_organism,
                    "convert_status": gprofiler_convert_status,
                    "gost_status": gprofiler_status,
                },
            },
            "counts": {
                "input_genes_count": len(input_symbols),
                "mapped_hgnc_count": len(valid_symbols),
                "string_edges_count": len(string_edges_df),
                "pathways_count": len(pathway_enrichment_df),
            },
            "interpretation_caveats": [
                "Research exploration only. Not for clinical diagnosis or treatment guidance.",
                "Network topology and co-membership DO NOT infer biological directionality.",
                "Direct edges in STRING reflect physical or functional association, not causality or drug suitability.",
                "Enrichment p-values reflect term overrepresentation among candidates, not disease-causing mechanisms.",
            ],
        }

        return {
            "seed_genes_df": seed_genes_df,
            "string_edges_df": string_edges_df,
            "network_summary_df": network_summary_df,
            "pathway_enrichment_df": pathway_enrichment_df,
            "network_gene_context_df": network_gene_context_df,
            "provenance": provenance,
        }

    def _extract_symbols(self, df: pd.DataFrame) -> list[str]:
        """Extract unique gene symbols from input mechanism or candidate table."""
        col = None
        for candidate in ["SYMBOL", "GENE_SYMBOL", "geneSymbol", "gene", "Gene"]:
            if candidate in df.columns:
                col = candidate
                break
        if col is None:
            raise ValueError(f"Input dataframe must have a gene symbol column (looked for SYMBOL, GENE_SYMBOL). Columns: {list(df.columns)}")

        symbols = [
            str(s).strip()
            for s in df[col].dropna()
            if str(s).strip() and str(s).strip().upper() not in {"NAN", "NONE", "."}
        ]
        # Preserve input order while eliminating duplicates
        seen = set()
        ordered = []
        for s in symbols:
            s_clean = s.upper()
            if s_clean not in seen:
                seen.add(s_clean)
                ordered.append(s)
        return ordered

    def map_hgnc_symbols(self, symbols: list[str]) -> tuple[dict[str, dict[str, Any]], str]:
        """Map and validate incoming symbols to official HGNC symbols using g:Profiler."""
        if not symbols:
            return {}, "empty_input"

        mapping: dict[str, dict[str, Any]] = {}
        payload = {
            "organism": self.config.gprofiler_organism,
            "target": self.config.gprofiler_target,
            "query": symbols,
        }

        resp = self.client.post(self.config.gprofiler_convert_url, json_data=payload)
        if resp is None or resp.status_code != 200:
            status = f"http_failure_{resp.status_code if resp else 'no_response'}"
            # Fallback: identity mapping
            for s in symbols:
                mapping[s] = {
                    "incoming": s,
                    "converted": s,
                    "description": "",
                    "status": "fallback_identity",
                }
            return mapping, status

        try:
            data = resp.json().get("result", [])
            for item in data:
                inc = str(item.get("incoming", "")).strip()
                conv = item.get("converted")
                name = item.get("name")
                desc = item.get("description", "")
                if inc:
                    mapping[inc] = {
                        "incoming": inc,
                        "converted": conv if conv and conv != "None" else inc,
                        "name": name if name and name != "None" else inc,
                        "description": desc or "",
                        "status": "mapped" if conv and conv != "None" else "unmapped_retained",
                    }
            # Catch any missing from result
            for s in symbols:
                if s not in mapping:
                    mapping[s] = {
                        "incoming": s,
                        "converted": s,
                        "name": s,
                        "description": "",
                        "status": "fallback_retained",
                    }
            return mapping, "success"
        except Exception as exc:
            logger.error("Failed to parse g:Profiler convert response: %s", exc)
            for s in symbols:
                mapping[s] = {"incoming": s, "converted": s, "status": "error_fallback"}
            return mapping, f"parse_error_{type(exc).__name__}"

    def get_string_version(self) -> str:
        """Query STRING API version endpoint to capture current version."""
        resp = self.client.get(self.config.string_version_url)
        if resp is not None and resp.status_code == 200:
            try:
                data = resp.json()
                if isinstance(data, list) and data:
                    return str(data[0].get("string_version", "12.5"))
            except Exception:
                pass
        return "12.5 (version-pinned endpoint; live version lookup unavailable)"

    def query_string_ppi(self, symbols: list[str]) -> tuple[pd.DataFrame, str]:
        """Query STRING v12.5 PPI network for input symbols."""
        columns = [
            "PREFERRED_NAME_A",
            "PREFERRED_NAME_B",
            "STRING_ID_A",
            "STRING_ID_B",
            "SCORE",
            "EXPERIMENTAL_SCORE",
            "DATABASE_SCORE",
            "TEXTMINING_SCORE",
            "COEXPRESSION_SCORE",
            "NEIGHBORHOOD_SCORE",
            "FUSION_SCORE",
            "COOCCURRENCE_SCORE",
            "TAXON_ID",
            "CURATED_SUPPORT",
            "TEXTMINING_ONLY",
        ]
        if len(symbols) < 2:
            return pd.DataFrame(columns=columns), "tiny_list_under_2_genes"

        # Batch the GET requests so large seed lists stay within URL limits.
        batch_size = max(1, int(self.config.string_batch_size))
        batches = [symbols[i : i + batch_size] for i in range(0, len(symbols), batch_size)]

        all_rows: list[dict[str, Any]] = []
        batch_statuses: list[str] = []
        for batch in batches:
            params = {
                "identifiers": "\n".join(batch),
                "species": self.config.string_species,
                "required_score": self.config.string_required_score,
                "add_nodes": self.config.string_additional_nodes,
                "network_flavor": self.config.string_network_flavor,
            }

            resp = self.client.get(self.config.string_api_url, params=params)
            if resp is None or resp.status_code != 200:
                # STRING returns 404 when NO identifier in a batch resolves.
                # A failed batch must not discard results from other batches.
                batch_status = f"http_failure_{resp.status_code if resp else 'no_response'}"
                batch_statuses.append(batch_status)
                logger.warning(
                    "STRING batch %d/%d failed (%s); skipping batch",
                    batches.index(batch) + 1,
                    len(batches),
                    batch_status,
                )
                continue

            try:
                data = resp.json()
            except Exception as exc:
                logger.error("Failed to parse STRING response: %s", exc)
                batch_statuses.append(f"parse_error_{type(exc).__name__}")
                continue
            if not isinstance(data, list):
                batch_statuses.append("unexpected_format")
                continue
            if not data:
                batch_statuses.append("zero_edges_found")
                continue

            for item in data:
                node_a = str(item.get("preferredName_A", "")).strip().upper()
                node_b = str(item.get("preferredName_B", "")).strip().upper()
                score = float(item.get("score", 0.0) or 0.0)

                # Channel scores (0.0 to 1.0)
                escore = float(item.get("escore", 0.0) or 0.0)
                dscore = float(item.get("dscore", 0.0) or 0.0)
                tscore = float(item.get("tscore", 0.0) or 0.0)
                ascore = float(item.get("ascore", 0.0) or 0.0)
                nscore = float(item.get("nscore", 0.0) or 0.0)
                fscore = float(item.get("fscore", 0.0) or 0.0)
                pscore = float(item.get("pscore", 0.0) or 0.0)

                # Apply score channel filters if configured
                if self.config.min_escore is not None and escore < self.config.min_escore:
                    continue
                if self.config.min_dscore is not None and dscore < self.config.min_dscore:
                    continue
                if self.config.min_tscore is not None and tscore < self.config.min_tscore:
                    continue
                if self.config.min_ascore is not None and ascore < self.config.min_ascore:
                    continue

                all_rows.append({
                    "PREFERRED_NAME_A": node_a,
                    "PREFERRED_NAME_B": node_b,
                    "STRING_ID_A": item.get("stringId_A", ""),
                    "STRING_ID_B": item.get("stringId_B", ""),
                    "SCORE": score,
                    "EXPERIMENTAL_SCORE": escore,
                    "DATABASE_SCORE": dscore,
                    "TEXTMINING_SCORE": tscore,
                    "COEXPRESSION_SCORE": ascore,
                    "NEIGHBORHOOD_SCORE": nscore,
                    "FUSION_SCORE": fscore,
                    "COOCCURRENCE_SCORE": pscore,
                    "TAXON_ID": item.get("ncbiTaxonId", self.config.string_species),
                    # Curated support = experimental or curated-database evidence.
                    # Text-mining-only edges are associations, not physical interactions.
                    "CURATED_SUPPORT": "YES" if max(escore, dscore) > 0.0 else "NO",
                    "TEXTMINING_ONLY": "YES" if (tscore > 0.0 and escore == 0.0 and dscore == 0.0) else "NO",
                })

        rows = all_rows
        df = pd.DataFrame(rows, columns=columns)
        if df.empty:
            failed_statuses = [
                s for s in batch_statuses
                if s.startswith("http_failure") or s.startswith("parse_error")
            ]
            if failed_statuses and "success" not in batch_statuses:
                return df, "all_batches_failed"
            return df, "zero_edges_found"

        # Deduplicate undirected edges (e.g. A-B vs B-A)
        df["EDGE_PAIR"] = df.apply(
            lambda r: tuple(sorted([r["PREFERRED_NAME_A"], r["PREFERRED_NAME_B"]])), axis=1
        )
        df = df.drop_duplicates(subset=["EDGE_PAIR"]).drop(columns=["EDGE_PAIR"])
        df = df.sort_values(by="SCORE", ascending=False).reset_index(drop=True)

        failed_count = sum(
            1 for s in batch_statuses
            if s.startswith("http_failure") or s.startswith("parse_error")
        )
        status = (
            "success"
            if failed_count == 0
            else f"partial_batch_failures_{failed_count}_of_{len(batch_statuses)}"
        )
        return df, status

    def compute_network_topology(
        self,
        nodes: list[str],
        edges_df: pd.DataFrame,
    ) -> tuple[pd.DataFrame, dict[str, dict[str, Any]]]:
        """Compute pure-Python graph metrics without requiring heavy C extensions."""
        upper_nodes = sorted(list({n.upper() for n in nodes if n}))
        adj: dict[str, set[str]] = {n: set() for n in upper_nodes}

        if not edges_df.empty:
            for _, r in edges_df.iterrows():
                u = str(r["PREFERRED_NAME_A"]).strip().upper()
                v = str(r["PREFERRED_NAME_B"]).strip().upper()
                if u not in adj:
                    adj[u] = set()
                if v not in adj:
                    adj[v] = set()
                if u != v:
                    adj[u].add(v)
                    adj[v].add(u)

        all_nodes = sorted(list(adj.keys()))
        n_nodes = len(all_nodes)
        n_edges = sum(len(neighbors) for neighbors in adj.values()) // 2

        # 1. Connected components (BFS)
        visited = set()
        components = []
        for n in all_nodes:
            if n not in visited:
                comp = set()
                queue = [n]
                visited.add(n)
                while queue:
                    curr = queue.pop(0)
                    comp.add(curr)
                    for nbr in adj[curr]:
                        if nbr not in visited:
                            visited.add(nbr)
                            queue.append(nbr)
                components.append(comp)

        n_components = len(components)
        largest_comp_size = max(len(c) for c in components) if components else 0

        # 2. Local clustering coefficients
        node_metrics: dict[str, dict[str, Any]] = {}
        clustering_coeffs = []
        for n in all_nodes:
            neighbors = adj[n]
            k = len(neighbors)
            if k < 2:
                cc = 0.0
            else:
                links_between_neighbors = 0
                nbr_list = list(neighbors)
                for i in range(k):
                    for j in range(i + 1, k):
                        if nbr_list[j] in adj[nbr_list[i]]:
                            links_between_neighbors += 1
                cc = (2.0 * links_between_neighbors) / (k * (k - 1))
            clustering_coeffs.append(cc)

            # Find which component this node belongs to
            comp_id = -1
            for idx, c in enumerate(components):
                if n in c:
                    comp_id = idx + 1
                    break

            node_metrics[n] = {
                "DEGREE": k,
                "CLUSTERING_COEFFICIENT": round(cc, 4),
                "COMPONENT_ID": comp_id,
                "IS_ISOLATED": k == 0,
                "CONNECTED_PARTNERS": sorted(list(neighbors)),
            }

        avg_degree = (2.0 * n_edges) / n_nodes if n_nodes > 0 else 0.0
        avg_clustering = sum(clustering_coeffs) / n_nodes if n_nodes > 0 else 0.0
        density = (2.0 * n_edges) / (n_nodes * (n_nodes - 1)) if n_nodes > 1 else 0.0

        summary_df = pd.DataFrame([
            {
                "TOTAL_SEED_GENES": len(upper_nodes),
                "TOTAL_NETWORK_NODES": n_nodes,
                "TOTAL_PPI_EDGES": n_edges,
                "AVERAGE_DEGREE": round(avg_degree, 3),
                "AVERAGE_CLUSTERING_COEFFICIENT": round(avg_clustering, 4),
                "NETWORK_DENSITY": round(density, 4),
                "CONNECTED_COMPONENTS_COUNT": n_components,
                "LARGEST_COMPONENT_SIZE": largest_comp_size,
                "ISOLATED_GENES_COUNT": sum(1 for m in node_metrics.values() if m["IS_ISOLATED"]),
                "EXPANSION_APPLIED": self.config.string_additional_nodes > 0,
                "STRING_SCORE_THRESHOLD": self.config.string_required_score,
            }
        ])

        return summary_df, node_metrics

    def query_reactome(self, symbols: list[str]) -> tuple[pd.DataFrame, str]:
        """Query Reactome Analysis Service for pathway overrepresentation.

        Uses the recommended '#Genes' header line for the text/plain body.
        The projection summary reports per-pathway entity counts only; it does
        NOT identify which submitted symbols hit each pathway, so HIT_GENES is
        left empty here (per-pathway hit genes require token-based follow-up).
        The analysis token is retained per row for later retrieval (valid 7 days).
        """
        columns = [
            "SOURCE_DATABASE",
            "PATHWAY_ID",
            "PATHWAY_NAME",
            "ENTITIES_FOUND",
            "ENTITIES_TOTAL",
            "P_VALUE",
            "FDR",
            "PATHWAY_SPECIES",
            "HIT_GENES",
            "MEMBERSHIP_STATUS",
            "REACTOME_TOKEN",
            "PATHWAY_BROWSER_URL",
        ]
        if not symbols:
            return pd.DataFrame(columns=columns), "empty_input"

        headers = {"Content-Type": "text/plain"}
        body = "#Genes\n" + "\n".join(symbols)
        params = {
            "pageSize": self.config.reactome_page_size,
            "page": 1,
            "species": "Homo sapiens",
        }

        resp = self.client.post(
            self.config.reactome_api_url,
            data=body,
            headers=headers,
            params=params,
        )
        if resp is None or resp.status_code != 200:
            status = f"http_failure_{resp.status_code if resp else 'no_response'}"
            return pd.DataFrame(columns=columns), status

        try:
            data = resp.json()
            pathways = data.get("pathways", [])
            token = str(data.get("summary", {}).get("token", "") or "")
            pathways_found = data.get("pathwaysFound", len(pathways))
            identifiers_not_found = data.get("identifiersNotFound", 0)
            if identifiers_not_found:
                logger.info(
                    "Reactome: %d submitted identifier(s) not mapped: %s",
                    identifiers_not_found,
                    data.get("identifiersNotFoundList", []),
                )
            if not pathways:
                return pd.DataFrame(columns=columns), "no_pathways_found"
            logger.info("Reactome analysis token: %s (%d pathways found)", token, pathways_found)

            rows = []
            for p in pathways:
                entities = p.get("entities", {})
                pathway_id = p.get("stId", "")
                browser_url = (
                    f"https://reactome.org/PathwayBrowser/#{pathway_id}&DTAB=AN&ANALYSIS={token}"
                    if token and pathway_id
                    else ""
                )
                rows.append({
                    "SOURCE_DATABASE": "Reactome",
                    "PATHWAY_ID": pathway_id,
                    "PATHWAY_NAME": p.get("name", ""),
                    "ENTITIES_FOUND": entities.get("found", 0),
                    "ENTITIES_TOTAL": entities.get("total", 0),
                    "P_VALUE": entities.get("pValue", 1.0),
                    "FDR": entities.get("fdr", 1.0),
                    "PATHWAY_SPECIES": p.get("species", {}).get("name", "Homo sapiens"),
                    "HIT_GENES": "",
                    # The projection summary reports counts only; it cannot name the contributors.
                    "MEMBERSHIP_STATUS": "UNAVAILABLE",
                    "REACTOME_TOKEN": token,
                    "PATHWAY_BROWSER_URL": browser_url,
                })

            df = pd.DataFrame(rows, columns=columns)
            df = df.sort_values(by=["FDR", "P_VALUE"]).reset_index(drop=True)
            return df, "success"
        except Exception as exc:
            logger.error("Failed to parse Reactome response: %s", exc)
            return pd.DataFrame(columns=columns), f"parse_error_{type(exc).__name__}"

    def query_gprofiler(self, symbols: list[str]) -> tuple[pd.DataFrame, str]:
        """Query g:Profiler g:GOSt API for multi-source enrichment."""
        columns = [
            "SOURCE_DATABASE",
            "PATHWAY_ID",
            "PATHWAY_NAME",
            "ENTITIES_FOUND",
            "ENTITIES_TOTAL",
            "P_VALUE",
            "FDR",
            "PATHWAY_SPECIES",
            "HIT_GENES",
            "MEMBERSHIP_STATUS",
        ]
        if not symbols:
            return pd.DataFrame(columns=columns), "empty_input"

        payload = {
            "organism": self.config.gprofiler_organism,
            "query": symbols,
            "sources": self.config.gprofiler_sources,
            "user_threshold": self.config.gprofiler_user_threshold,
        }

        resp = self.client.post(self.config.gprofiler_gost_url, json_data=payload)
        if resp is None or resp.status_code != 200:
            status = f"http_failure_{resp.status_code if resp else 'no_response'}"
            return pd.DataFrame(columns=columns), status

        try:
            data = resp.json().get("result", [])
            if not data:
                return pd.DataFrame(columns=columns), "no_terms_found"

            rows = []
            for term in data:
                # The intersection vector is only meaningful when it aligns with the submitted
                # query order.  A missing or mismatched vector stays UNAVAILABLE: submitted
                # genes are never attributed to a term as a fallback.
                hit_genes, membership_status = self._decode_intersections(term, symbols)
                rows.append({
                    "SOURCE_DATABASE": f"g:Profiler ({term.get('source', '')})",
                    "PATHWAY_ID": term.get("native", ""),
                    "PATHWAY_NAME": term.get("name", ""),
                    "ENTITIES_FOUND": term.get("intersection_size", 0),
                    "ENTITIES_TOTAL": term.get("term_size", 0),
                    "P_VALUE": term.get("p_value", 1.0),
                    "FDR": term.get("p_value", 1.0),  # g:Profiler returns multiple-testing corrected p-value
                    "PATHWAY_SPECIES": self.config.gprofiler_organism,
                    "HIT_GENES": ";".join(hit_genes),
                    "MEMBERSHIP_STATUS": membership_status,
                })

            df = pd.DataFrame(rows, columns=columns)
            df = df.sort_values(by="P_VALUE").reset_index(drop=True)
            return df, "success"
        except Exception as exc:
            logger.error("Failed to parse g:Profiler g:GOSt response: %s", exc)
            return pd.DataFrame(columns=columns), f"parse_error_{type(exc).__name__}"

    @staticmethod
    def _decode_intersections(term: Mapping[str, Any], symbols: Sequence[str]) -> tuple[list[str], str]:
        """Decode an explicit intersection vector, or report membership as unavailable."""
        intersections = term.get("intersections")
        if not isinstance(intersections, list) or len(intersections) != len(symbols):
            return [], "UNAVAILABLE"
        if any(not isinstance(flag, bool) for flag in intersections):
            return [], "UNAVAILABLE"
        decoded = [symbol for symbol, is_hit in zip(symbols, intersections) if is_hit]
        return decoded, "OBSERVED_INTERSECTION" if decoded else "UNAVAILABLE"

    def _combine_enrichment(
        self,
        reactome_df: pd.DataFrame,
        gprofiler_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """Combine and normalize enrichment results across engines."""
        combined = pd.concat([reactome_df, gprofiler_df], ignore_index=True)
        if combined.empty:
            return pd.DataFrame(
                columns=[
                    "SOURCE_DATABASE",
                    "PATHWAY_ID",
                    "PATHWAY_NAME",
                    "ENTITIES_FOUND",
                    "ENTITIES_TOTAL",
                    "P_VALUE",
                    "FDR",
                    "PATHWAY_SPECIES",
                    "HIT_GENES",
                    "MEMBERSHIP_STATUS",
                ]
            )

        combined["P_VALUE"] = pd.to_numeric(combined["P_VALUE"], errors="coerce").fillna(1.0)
        combined["FDR"] = pd.to_numeric(combined["FDR"], errors="coerce").fillna(1.0)
        return combined.sort_values(by=["FDR", "P_VALUE"]).reset_index(drop=True)

    def _build_seed_genes_df(
        self,
        symbols: list[str],
        hgnc_mapping: dict[str, dict[str, Any]],
        gene_mech_df: pd.DataFrame,
        var_matrix_df: pd.DataFrame | None = None,
    ) -> pd.DataFrame:
        """Build seed gene tracking table with HGNC status and Phase 1 context."""
        # Index gene mechanism table by symbol
        mech_by_sym = {}
        for _, row in gene_mech_df.iterrows():
            sym = str(row.get("SYMBOL", "")).strip().upper()
            if sym:
                mech_by_sym[sym] = row

        rows = []
        for s in symbols:
            s_upper = s.upper()
            map_info = hgnc_mapping.get(s, {})
            hgnc_sym = map_info.get("converted", s_upper)
            if not hgnc_sym or hgnc_sym == "None":
                hgnc_sym = s_upper

            # Align with the real Phase 1 gene-matrix schema (mechanism_matrix.py).
            mech_row = mech_by_sym.get(s_upper, {})
            pathway = mech_row.get("PATHWAY", "")
            mva_link = mech_row.get("MVA_MECHANISM_LINK", "")
            pathway_dir = mech_row.get("CANONICAL_PATHWAY_DIRECTION", "")
            therapeutic_dir = mech_row.get("THERAPEUTIC_DIRECTIONALITY", "")
            inheritance = mech_row.get("INHERITANCE_REQUIREMENT", "")
            n_vars = mech_row.get("TOTAL_VARIANTS_EVALUATED", 0)
            n_passing = mech_row.get("TOTAL_PASSING_CANDIDATE_VARIANTS", 0)
            archetype_fit = mech_row.get("COMPOUND_HET_ARCHETYPE_FIT", False)
            archetype_desc = mech_row.get("ARCHETYPE_DESCRIPTION", "")
            phase_gate = mech_row.get("PHASE_GATE_STATUS", "")

            rows.append({
                "INPUT_SYMBOL": s,
                "HGNC_SYMBOL": hgnc_sym,
                "HGNC_STATUS": map_info.get("status", "unmapped"),
                "DESCRIPTION": map_info.get("description", ""),
                "PATHWAY": pathway,
                "MVA_MECHANISM_LINK": mva_link,
                "CANONICAL_PATHWAY_DIRECTION": pathway_dir,
                "THERAPEUTIC_DIRECTIONALITY": therapeutic_dir,
                "INHERITANCE_REQUIREMENT": inheritance,
                "TOTAL_VARIANTS_EVALUATED": n_vars,
                "TOTAL_PASSING_CANDIDATE_VARIANTS": n_passing,
                "COMPOUND_HET_ARCHETYPE_FIT": archetype_fit,
                "ARCHETYPE_DESCRIPTION": archetype_desc,
                "PHASE_GATE_STATUS": phase_gate,
            })

        return pd.DataFrame(rows)

    def _build_gene_context_df(
        self,
        seed_genes_df: pd.DataFrame,
        node_topologies: dict[str, dict[str, Any]],
        enrichment_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """Synthesize Phase 1 mechanism with Phase 2 topological and pathway context."""
        rows = []
        for _, s_row in seed_genes_df.iterrows():
            hgnc = s_row["HGNC_SYMBOL"]
            topo = node_topologies.get(hgnc, {})

            degree = topo.get("DEGREE", 0)
            cc = topo.get("CLUSTERING_COEFFICIENT", 0.0)
            comp_id = topo.get("COMPONENT_ID", -1)
            is_iso = topo.get("IS_ISOLATED", True)
            partners = ";".join(topo.get("CONNECTED_PARTNERS", []))

            # Find top enriched pathways containing this gene.
            # Only rows with concrete per-pathway hit genes count (g:Profiler
            # intersections). Reactome projection rows carry empty HIT_GENES by
            # design (the summary API does not report per-pathway hit symbols).
            gene_pathways = []
            if not enrichment_df.empty:
                for _, p_row in enrichment_df.iterrows():
                    hits_raw = str(p_row.get("HIT_GENES", "") or "").strip()
                    if not hits_raw:
                        continue
                    hits = hits_raw.split(";")
                    if hgnc in [h.strip().upper() for h in hits]:
                        p_name = p_row.get("PATHWAY_NAME", "")
                        fdr = p_row.get("FDR", 1.0)
                        gene_pathways.append(f"{p_name} (FDR: {fdr:.2e})")
                        if len(gene_pathways) >= 3:
                            break

            rows.append({
                "HGNC_SYMBOL": hgnc,
                "INPUT_SYMBOL": s_row["INPUT_SYMBOL"],
                "PATHWAY": s_row["PATHWAY"],
                "MVA_MECHANISM_LINK": s_row["MVA_MECHANISM_LINK"],
                "CANONICAL_PATHWAY_DIRECTION": s_row["CANONICAL_PATHWAY_DIRECTION"],
                "THERAPEUTIC_DIRECTIONALITY": s_row["THERAPEUTIC_DIRECTIONALITY"],
                "INHERITANCE_REQUIREMENT": s_row["INHERITANCE_REQUIREMENT"],
                "TOTAL_VARIANTS_EVALUATED": s_row["TOTAL_VARIANTS_EVALUATED"],
                "TOTAL_PASSING_CANDIDATE_VARIANTS": s_row["TOTAL_PASSING_CANDIDATE_VARIANTS"],
                "COMPOUND_HET_ARCHETYPE_FIT": s_row["COMPOUND_HET_ARCHETYPE_FIT"],
                "PHASE_GATE_STATUS": s_row["PHASE_GATE_STATUS"],
                "NETWORK_DEGREE": degree,
                "CLUSTERING_COEFFICIENT": cc,
                "COMPONENT_ID": comp_id,
                "IS_ISOLATED": is_iso,
                "STRING_INTERACTORS": partners,
                "TOP_ENRICHED_PATHWAYS": "; ".join(gene_pathways) if gene_pathways else "None",
                "NETWORK_ROBUSTNESS_NOTE": "Direct PPI cluster" if degree > 0 else "Singleton / unlinked in primary threshold",
            })

        df = pd.DataFrame(rows)
        if not df.empty:
            df = df.sort_values(by=["NETWORK_DEGREE", "TOTAL_PASSING_CANDIDATE_VARIANTS"], ascending=[False, False]).reset_index(drop=True)
        return df


def write_network_outputs(
    results: Mapping[str, Any],
    output_dir: Path,
    flat_prefix: bool = False,
) -> dict[str, Path]:
    """Write all Phase 2 result tables and provenance record to output directory.

    If flat_prefix is True, writes ordered Track 2 filenames:
    - 03_phase2_seed_genes.csv
    - 04_phase2_string_edges.csv
    - 05_phase2_network_summary.csv
    - 06_phase2_pathway_enrichment.csv
    - 07_phase2_network_gene_context.csv
    - 08_phase2_network_provenance.json
    Otherwise writes standard un-prefixed phase filenames.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}

    table_mapping = [
        ("seed_genes_df", "03_phase2_seed_genes.csv" if flat_prefix else "seed_genes.csv"),
        ("string_edges_df", "04_phase2_string_edges.csv" if flat_prefix else "string_edges.csv"),
        ("network_summary_df", "05_phase2_network_summary.csv" if flat_prefix else "network_summary.csv"),
        ("pathway_enrichment_df", "06_phase2_pathway_enrichment.csv" if flat_prefix else "pathway_enrichment.csv"),
        ("network_gene_context_df", "07_phase2_network_gene_context.csv" if flat_prefix else "network_gene_context.csv"),
    ]

    for key, filename in table_mapping:
        df: pd.DataFrame = results[key]
        dest = output_dir / filename
        df.to_csv(dest, index=False)
        paths[key] = dest

    prov_filename = "08_phase2_network_provenance.json" if flat_prefix else "phase2_network_provenance.json"
    prov_dest = output_dir / prov_filename
    with prov_dest.open("w", encoding="utf-8") as f:
        json.dump(results["provenance"], f, indent=2)
    paths["provenance"] = prov_dest

    return paths
