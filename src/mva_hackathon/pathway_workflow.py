"""Independent, privacy-preserving pathway workflow.

Only normalized gene identifiers cross an external-service boundary.  This is
research exploration, not diagnostic or therapeutic guidance.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
import hashlib
import json
from pathlib import Path
import re
from typing import Any

import pandas as pd

from mva_hackathon.mva_gene_sets import CORE_MVA_CAUSAL_GENES

_REQUIRED_COLUMNS = (
    "SYMBOL",
    "TOTAL_VARIANTS_EVALUATED",
    "TOTAL_PASSING_CANDIDATE_VARIANTS",
)
_OPTIONAL_COLUMNS = (
    "THERAPEUTIC_DIRECTIONALITY",
    "PASSING_CONSEQUENCE_SPECTRUM",
)
_NORMALIZED_COLUMNS = (*_REQUIRED_COLUMNS, *_OPTIONAL_COLUMNS, "IN_FOREGROUND")
_GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
_HGNC_URL = "https://rest.genenames.org/fetch/symbol/"
_REACTOME_ANALYSIS_URL = "https://reactome.org/AnalysisService/identifiers/projection"
_REACTOME_CONTENT_URL = "https://reactome.org/ContentService/data/query/"
_UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search"
_QUICKGO_URL = "https://www.ebi.ac.uk/QuickGO/services/annotation/search"
_INHIBITOR_TOKENS = (
    "INHIBITOR",
    "ANTAGONIST",
    "BLOCKER",
    "DOWNREGULATOR",
    "DEGRADER",
    "RNAI",
    "SUPPRESSOR",
)


@dataclass(frozen=True)
class PathwayWorkflowConfig:
    """Configuration that affects reproducible pathway computations."""

    min_passing: int = 1
    q_threshold: float = 0.05
    top_pathways: int = 15
    offline: bool = False

    def __post_init__(self) -> None:
        if self.min_passing < 0:
            raise ValueError("min_passing must be non-negative")
        if not 0 < self.q_threshold <= 1:
            raise ValueError("q_threshold must be in (0, 1]")
        if self.top_pathways < 1:
            raise ValueError("top_pathways must be positive")


class PathwayWorkflow:
    """Run P0–P7 while retaining source-native evidence and partial failures."""

    def __init__(self, config: PathwayWorkflowConfig | None = None) -> None:
        self.config = config or PathwayWorkflowConfig()

    @staticmethod
    def request_cache_key(endpoint: str, payload: dict[str, Any]) -> str:
        """Hash an endpoint and canonical non-secret payload."""
        encoded = json.dumps(
            {"endpoint": endpoint, "payload": payload},
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
        return hashlib.sha256(encoded).hexdigest()

    @staticmethod
    def atomic_write_csv(frame: pd.DataFrame, path: Path) -> None:
        """Write a completed table atomically, never exposing partial output."""
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(f"{path.suffix}.tmp")
        frame.to_csv(temporary, index=False)
        temporary.replace(path)

    @staticmethod
    def atomic_write_json(payload: dict[str, Any], path: Path) -> None:
        """Write provenance/checkpoint JSON atomically."""
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(f"{path.suffix}.tmp")
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")
        temporary.replace(path)

    @staticmethod
    def file_hash(path: Path) -> str:
        """Hash a local artifact for checkpoint invalidation without parsing it."""
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()

    def map_hgnc(self, symbols: Iterable[str], client: Any) -> pd.DataFrame:
        """Resolve exact HGNC symbols; ambiguous results are deliberately unmapped."""
        records: list[dict[str, Any]] = []
        for symbol in sorted({str(item).strip().upper() for item in symbols if str(item).strip()}):
            response = client.get(f"{_HGNC_URL}{symbol}")
            try:
                docs = (
                    response.json().get("response", {}).get("docs", [])
                    if response is not None and response.status_code == 200
                    else []
                )
            except (ValueError, TypeError):
                docs = []
            exact = [doc for doc in docs if doc.get("symbol") == symbol]
            ensembl_ids = {doc.get("ensembl_gene_id") for doc in exact if doc.get("ensembl_gene_id")}
            if len(exact) == 1 and len(ensembl_ids) == 1:
                records.append(
                    {
                        "SYMBOL": symbol,
                        "HGNC_ID": exact[0].get("hgnc_id"),
                        "ENSEMBL_ID": next(iter(ensembl_ids)),
                        "STATUS": "mapped",
                    }
                )
            else:
                records.append(
                    {
                        "SYMBOL": symbol,
                        "HGNC_ID": pd.NA,
                        "ENSEMBL_ID": pd.NA,
                        "STATUS": "unmapped",
                    }
                )
        return pd.DataFrame(records, columns=["SYMBOL", "HGNC_ID", "ENSEMBL_ID", "STATUS"])

    def mapping_from_gprofiler(self, normalized: pd.DataFrame, response: dict[str, Any]) -> pd.DataFrame:
        """Build an offline-safe mapping from GOSt's explicit query mapping only."""
        mapping = (
            response.get("meta", {})
            .get("query_metadata", {})
            .get("queries", {})
            .get("query_1", {})
            .get("mapping", {})
        )
        symbol_to_ids: dict[str, set[str]] = {}
        for ensembl_id, symbols in mapping.items():
            for symbol in symbols or []:
                symbol_to_ids.setdefault(str(symbol).upper(), set()).add(str(ensembl_id))
        records: list[dict[str, Any]] = []
        for symbol in normalized["SYMBOL"]:
            ids = symbol_to_ids.get(symbol, set())
            records.append(
                {
                    "SYMBOL": symbol,
                    "HGNC_ID": pd.NA,
                    "ENSEMBL_ID": next(iter(ids)) if len(ids) == 1 else pd.NA,
                    "STATUS": "mapped" if len(ids) == 1 else "unmapped",
                }
            )
        return pd.DataFrame(records)

    def gprofiler_payload(self, mapped_genes: pd.DataFrame, domain_scope: str = "custom") -> dict[str, Any]:
        """Build GOSt request with gene IDs only and no patient-level data."""
        status = mapped_genes.get("STATUS", pd.Series("mapped", index=mapped_genes.index))
        resolved = mapped_genes.loc[
            status.eq("mapped") & mapped_genes["ENSEMBL_ID"].notna()
        ].copy()
        background = sorted(resolved["ENSEMBL_ID"].astype(str).unique())
        query = sorted(resolved.loc[resolved["IN_FOREGROUND"], "ENSEMBL_ID"].astype(str).unique())
        if not background or not query:
            raise ValueError("resolved foreground and background are required for enrichment")
        return {
            "organism": "hsapiens",
            "query": query,
            "sources": ["GO:BP", "REAC", "WP", "HP"],
            "domain_scope": domain_scope,
            "background": background,
            "significance_threshold_method": "fdr",
            "all_results": True,
            "no_evidences": False,
            "user_threshold": self.config.q_threshold,
        }

    def fetch_gprofiler(self, mapped_genes: pd.DataFrame, client: Any, domain_scope: str = "custom") -> dict[str, Any]:
        """Call GOSt and reject transport failures rather than treating them as zero terms."""
        response = client.post(_GPROFILER_URL, json_data=self.gprofiler_payload(mapped_genes, domain_scope))
        if response is None or response.status_code != 200:
            raise RuntimeError("gProfiler request failed")
        return response.json()

    def fetch_reactome_analysis(self, symbols: Iterable[str], client: Any) -> dict[str, Any]:
        """Obtain secondary Reactome method-comparison results."""
        body = "#Genes\n" + "\n".join(sorted({str(symbol).upper() for symbol in symbols})) + "\n"
        response = client.post(_REACTOME_ANALYSIS_URL, data=body, headers={"Content-Type": "text/plain"})
        if response is None or response.status_code != 200:
            raise RuntimeError("Reactome analysis request failed")
        return response.json()

    def run_p0(self, gene_mechanism_path: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
        """Read only the Phase 1 gene-level allowlist and normalize it."""
        separator = "\t" if gene_mechanism_path.suffix.lower() in {".tsv", ".txt"} else ","
        headers = pd.read_csv(gene_mechanism_path, sep=separator, nrows=0).columns
        missing = [column for column in _REQUIRED_COLUMNS if column not in headers]
        if missing:
            raise ValueError(f"missing required columns: {', '.join(missing)}")
        columns = [*_REQUIRED_COLUMNS, *(column for column in _OPTIONAL_COLUMNS if column in headers)]
        frame = pd.read_csv(gene_mechanism_path, sep=separator, usecols=columns)
        for column in _OPTIONAL_COLUMNS:
            if column not in frame:
                frame[column] = pd.NA
        frame["SYMBOL"] = frame["SYMBOL"].fillna("").astype(str).str.strip().str.upper()
        quarantined = int(frame["SYMBOL"].eq("").sum())
        frame = frame.loc[frame["SYMBOL"].ne("")].copy()
        for column in _REQUIRED_COLUMNS[1:]:
            values = pd.to_numeric(frame[column], errors="coerce")
            if values.isna().any() or (values < 0).any():
                raise ValueError(f"invalid non-negative count in {column}")
            frame[column] = values.astype(int)
        duplicate_count = int(frame.duplicated().sum())
        duplicates = frame.loc[frame.duplicated("SYMBOL", keep=False)]
        if not duplicates.empty:
            conflicts = duplicates.groupby("SYMBOL").filter(lambda rows: len(rows.drop_duplicates()) > 1)
            if not conflicts.empty:
                raise ValueError("conflicting duplicate rows for SYMBOL")
            frame = frame.drop_duplicates("SYMBOL", keep="first")
        frame["IN_FOREGROUND"] = frame["TOTAL_PASSING_CANDIDATE_VARIANTS"] >= self.config.min_passing
        frame = frame.loc[:, _NORMALIZED_COLUMNS].sort_values("SYMBOL").reset_index(drop=True)
        foreground_count = int(frame["IN_FOREGROUND"].sum())
        report = {
            "stage": "p0",
            "status": "completed" if len(frame) and foreground_count else "failed",
            "background_count": len(frame),
            "foreground_count": foreground_count,
            "duplicate_count": duplicate_count,
            "quarantined_blank_symbols": quarantined,
            "background_limitation": "Background is Phase 1 evaluated genes, not whole-exome callable genes.",
        }
        return frame, report

    @staticmethod
    def _decoded_intersection(response: dict[str, Any], identifiers: Iterable[Any]) -> list[str]:
        mapping = (
            response.get("meta", {})
            .get("query_metadata", {})
            .get("queries", {})
            .get("query_1", {})
            .get("mapping", {})
        )
        return sorted(
            {
                str(symbol)
                for identifier in identifiers
                for symbol in mapping.get(str(identifier), [])
            }
        )

    def run_p1(self, normalized: pd.DataFrame, response: dict[str, Any]) -> pd.DataFrame:
        """Normalize GOSt output, preserving provider BH values without re-correction."""
        from scipy.stats import hypergeom

        records: list[dict[str, Any]] = []
        for term in response.get("result", []):
            m, k, n, x = (
                term.get("effective_domain_size"),
                term.get("term_size"),
                term.get("query_size"),
                term.get("intersection_size"),
            )
            valid = all(isinstance(value, int) and value >= 0 for value in (m, k, n, x))
            valid = valid and bool(m) and k <= m and n <= m and x <= min(k, n)
            raw_p = hypergeom.sf(x - 1, m, k, n) if valid else None
            fold = ((x / n) / (k / m)) if valid and k and n else None
            records.append(
                {
                    "SOURCE": term.get("source"),
                    "TERM_ID": term.get("native"),
                    "TERM_NAME": term.get("name"),
                    "PATHWAY_DB": term.get("source"),
                    "INTERSECTION_SIZE": x,
                    "TERM_SIZE": k,
                    "QUERY_SIZE": n,
                    "BACKGROUND_SIZE": m,
                    "RAW_P": raw_p,
                    "BH_Q": term.get("p_value"),
                    "PROVIDER_ADJUSTED_P": term.get("p_value"),
                    "PROVIDER_CORRECTION": "BH",
                    "FOLD_ENRICHMENT": fold,
                    "CONTRIBUTING_GENES": ";".join(self._decoded_intersection(response, term.get("intersections") or [])),
                    "BACKGROUND_LABEL": "Phase 1 evaluated genes",
                    "LOCAL_STATS_REASON": pd.NA if valid else "PROVIDER_EFFECTIVE_COUNTS_UNAVAILABLE",
                }
            )
        return pd.DataFrame(records)

    def run_p2(self, enrichment: pd.DataFrame) -> pd.DataFrame:
        """Prioritize significant terms using stated components, never a hidden score."""
        if enrichment.empty:
            return enrichment.copy()
        results = enrichment.loc[enrichment["BH_Q"].le(self.config.q_threshold)].copy()
        if results.empty:
            return results
        results["CONTRIBUTING_GENE_SET"] = results["CONTRIBUTING_GENES"].fillna("").map(
            lambda value: frozenset(gene for gene in value.split(";") if gene)
        )
        pattern = re.compile(r"MITOT|SPINDLE|CHROMOSOME SEGREGAT|ANEUPLOID|CHECKPOINT|COHESIN|CENTROSOME", re.I)
        results["CORE_MVA_GENE_OVERLAP"] = results["CONTRIBUTING_GENE_SET"].map(
            lambda genes: ";".join(sorted(genes & CORE_MVA_CAUSAL_GENES))
        )
        results["TERM_NAME_MVA_MATCH"] = results["TERM_NAME"].fillna("").map(lambda value: bool(pattern.search(value)))
        results["MVA_RELEVANT"] = results["CORE_MVA_GENE_OVERLAP"].ne("") | results["TERM_NAME_MVA_MATCH"]
        results["OBSERVATIONAL_PHENOTYPE_ANNOTATION"] = results["TERM_ID"].fillna("").str.startswith("HP:")
        results["REDUNDANT_OF"] = pd.NA
        representatives: list[tuple[str, frozenset[str]]] = []
        for index, row in results.sort_values(["BH_Q", "SOURCE", "TERM_ID"], kind="stable").iterrows():
            genes = row["CONTRIBUTING_GENE_SET"]
            redundant = next(
                (
                    term_id
                    for term_id, representative in representatives
                    if genes and representative and len(genes & representative) / len(genes | representative) > 0.6
                ),
                None,
            )
            if redundant is None:
                representatives.append((row["TERM_ID"], genes))
            else:
                results.at[index, "REDUNDANT_OF"] = redundant
        return results.sort_values(
            ["MVA_RELEVANT", "BH_Q", "FOLD_ENRICHMENT", "INTERSECTION_SIZE", "SOURCE", "TERM_ID"],
            ascending=[False, True, False, False, True, True],
            kind="stable",
        ).drop(columns="CONTRIBUTING_GENE_SET").reset_index(drop=True)

    def reactome_members(self, stable_id: str, client: Any) -> list[dict[str, Any]]:
        """Retrieve complete Reactome event members, preserving only gene identifiers."""
        response = client.get(f"{_REACTOME_CONTENT_URL}{stable_id}")
        if response is None or response.status_code != 200:
            raise RuntimeError("Reactome content request failed")
        payload = response.json()
        members: list[dict[str, Any]] = []
        for entity in payload.get("participants", []) + payload.get("hasEvent", []):
            reference = entity.get("referenceEntity") or entity
            symbol = reference.get("geneName") or reference.get("displayName")
            if symbol:
                members.append({"SYMBOL": str(symbol).upper(), "ENSEMBL_ID": reference.get("identifier")})
        return members


    def mechanism_evidence(self, prioritized: pd.DataFrame, client: Any) -> list[dict[str, Any]]:
        """Fetch lightweight function/GO evidence for top Reactome intersections."""
        records: list[dict[str, Any]] = []
        reactome_terms = prioritized.loc[prioritized["SOURCE"].eq("REAC")].head(self.config.top_pathways)
        for term in reactome_terms.itertuples(index=False):
            contributed = term.CONTRIBUTING_GENES if pd.notna(term.CONTRIBUTING_GENES) else ""
            for symbol in str(contributed).split(";"):
                if not symbol:
                    continue
                protein = client.get(
                    _UNIPROT_URL,
                    params={
                        "query": f"gene_exact:{symbol} AND organism_id:9606",
                        "format": "tsv",
                        "fields": "accession,function,cc_subunit",
                        "size": 1,
                    },
                )
                if protein is None or protein.status_code != 200:
                    raise RuntimeError("UniProt request failed")
                rows = protein.text.splitlines()
                accession = rows[1].split("\t", 1)[0] if len(rows) > 1 else None
                evidence_tier = "UNAVAILABLE"
                if accession:
                    annotations = client.get(
                        _QUICKGO_URL,
                        params={"geneProductId": f"UniProtKB:{accession}", "limit": 100},
                    )
                    if annotations is None or annotations.status_code != 200:
                        raise RuntimeError("QuickGO request failed")
                    codes = {
                        str(item.get("goEvidence") or item.get("evidenceCode") or "").upper()
                        for item in annotations.json().get("results", [])
                    }
                    evidence_tier = "EXPERIMENTAL" if "IDA" in codes else "COMPUTATIONAL" if codes else "UNAVAILABLE"
                records.append(
                    {
                        "TERM_ID": term.TERM_ID,
                        "SYMBOL": symbol,
                        "ROLE": "enrichment_intersection",
                        "UNIPROT_ACCESSION": accession,
                        "EVIDENCE_TIER": evidence_tier,
                    }
                )
        return records

    def run_p3(self, phase1: pd.DataFrame, evidence_rows: list[dict[str, Any]]) -> pd.DataFrame:
        """Attach Phase 1 direction only where a pathway gene is variant-supported."""
        directions = phase1.drop_duplicates("SYMBOL").set_index("SYMBOL").get("THERAPEUTIC_DIRECTIONALITY", pd.Series(dtype=object))
        records: list[dict[str, Any]] = []
        for evidence in evidence_rows:
            record = dict(evidence)
            symbol = str(record.get("SYMBOL", "")).upper()
            record.update(
                {
                    "SYMBOL": symbol,
                    "THERAPEUTIC_DIRECTIONALITY": directions.get(symbol, "UNRESOLVED_DIRECTION"),
                    "VARIANT_SUPPORTED": symbol in directions.index,
                }
            )
            records.append(record)
        return pd.DataFrame(records)

    def run_p4(self, prioritized: pd.DataFrame, memberships: dict[str, list[dict[str, Any]]], phase1: pd.DataFrame) -> pd.DataFrame:
        """Expand selected pathways to full supplied memberships, not intersections alone."""
        directions = phase1.drop_duplicates("SYMBOL").set_index("SYMBOL").get("THERAPEUTIC_DIRECTIONALITY", pd.Series(dtype=object))
        records: list[dict[str, Any]] = []
        for term_id in prioritized["TERM_ID"].dropna().head(self.config.top_pathways):
            for member in memberships.get(term_id, []):
                symbol = str(member["SYMBOL"]).upper()
                supported = symbol in directions.index
                records.append(
                    {
                        "TERM_ID": term_id,
                        "SYMBOL": symbol,
                        "ENSEMBL_ID": member.get("ENSEMBL_ID"),
                        "VARIANT_SUPPORTED": supported,
                        "TARGET_SUPPORT": "VARIANT_SUPPORTED" if supported else "PATHWAY_EXPLORATORY",
                        "THERAPEUTIC_DIRECTIONALITY": directions.get(symbol, "UNRESOLVED_DIRECTION"),
                    }
                )
        columns = ["TERM_ID", "SYMBOL", "ENSEMBL_ID", "VARIANT_SUPPORTED", "TARGET_SUPPORT", "THERAPEUTIC_DIRECTIONALITY"]
        return pd.DataFrame(records, columns=columns).drop_duplicates(["TERM_ID", "SYMBOL"]).sort_values(["TERM_ID", "SYMBOL"]).reset_index(drop=True)

    def run_p5(self, targets: pd.DataFrame, source_records: list[dict[str, Any]]) -> pd.DataFrame:
        """Preserve source records and label contradiction hazards independently of efficacy."""
        directions = targets.drop_duplicates("SYMBOL").set_index("SYMBOL")["THERAPEUTIC_DIRECTIONALITY"]
        records: list[dict[str, Any]] = []
        for source_record in source_records:
            record = dict(source_record)
            symbol = str(record.get("TARGET_SYMBOL", "")).upper()
            action = str(record.get("ACTION_TYPE", "")).upper()
            direction = str(directions.get(symbol, "UNRESOLVED_DIRECTION")).upper()
            contradiction = "RESTOR" in direction and any(token in action for token in _INHIBITOR_TOKENS)
            record.update(
                {
                    "TARGET_SYMBOL": symbol,
                    "THERAPEUTIC_DIRECTIONALITY": direction,
                    "CONTRADICTION_FLAG": contradiction,
                    "DIRECTIONALITY_LABEL": "CONTRADICTION" if contradiction else "UNKNOWN",
                    "RESEARCH_ONLY": True,
                }
            )
            records.append(record)
        return pd.DataFrame(records)

    def run_p6(self, pathway_targets: pd.DataFrame, network_panel: pd.DataFrame | None) -> dict[str, Any]:
        """Provide descriptive target overlap; no unjustified universe-based test."""
        pathway_symbols = {str(symbol).upper() for symbol in pathway_targets.get("SYMBOL", pd.Series(dtype=str)).dropna()}
        if network_panel is None:
            return {"status": "partial", "reason": "STRING_BASELINE_UNAVAILABLE", "count": len(pathway_symbols)}
        column = "QUERY_GENE" if "QUERY_GENE" in network_panel else "SYMBOL"
        string_symbols = {str(symbol).upper() for symbol in network_panel[column].dropna()}
        union = pathway_symbols | string_symbols
        return {
            "status": "completed",
            "pathway_only_targets": len(pathway_symbols - string_symbols),
            "string_only_targets": len(string_symbols - pathway_symbols),
            "shared_targets": len(pathway_symbols & string_symbols),
            "combined_targets": len(union),
            "jaccard": len(pathway_symbols & string_symbols) / len(union) if union else None,
            "overlap_test": "UNAVAILABLE: no prespecified eligible target universe",
        }

    def provenance(self, stages: dict[str, dict[str, Any]]) -> dict[str, Any]:
        """Create a compact provenance payload with no source table values."""
        return {
            "research_only": True,
            "privacy_boundary": "Only normalized gene identifiers were disclosed to external services.",
            "generated_at": datetime.now(UTC).isoformat(),
            "configuration": asdict(self.config),
            "stages": stages,
        }

    def run_p7(
        self,
        report_path: Path,
        stage_reports: dict[str, dict[str, Any]],
        pathways: pd.DataFrame | None = None,
        targets: pd.DataFrame | None = None,
        drugs: pd.DataFrame | None = None,
    ) -> None:
        """Write Markdown-only cumulative reporting with explicit partial states."""
        p0 = stage_reports.get("p0", {})
        lines = [
            "# Pathway workflow report",
            "",
            "Research exploration only; enrichment is not efficacy and pharmacology is not biological rescue.",
            "",
            "## Run header",
            "",
            f"- Status: {p0.get('status', 'UNAVAILABLE')}",
            f"- Foreground genes: {p0.get('foreground_count', 'UNAVAILABLE')}",
            f"- Background genes: {p0.get('background_count', 'UNAVAILABLE')}",
            "",
            "## Stage QC",
            "",
            "| Stage | Status | Count |",
            "| --- | --- | --- |",
        ]
        for stage in ("p0", "p1", "p2", "p3", "p4", "p5", "p6"):
            report = stage_reports.get(stage, {})
            lines.append(f"| {stage} | {report.get('status', 'UNAVAILABLE')} | {report.get('count', report.get('foreground_count', 'UNAVAILABLE'))} |")
        lines.extend(["", "## Top pathways", ""])
        lines.extend(self._markdown_table(pathways, ["TERM_ID", "TERM_NAME", "BH_Q", "FOLD_ENRICHMENT"]) if pathways is not None and not pathways.empty else ["UNAVAILABLE"])
        lines.extend(["", "## Top targets", ""])
        lines.extend(self._markdown_table(targets, ["TERM_ID", "SYMBOL", "TARGET_SUPPORT", "THERAPEUTIC_DIRECTIONALITY"]) if targets is not None and not targets.empty else ["UNAVAILABLE"])
        lines.extend(["", "## Top drug hypotheses", ""])
        lines.extend(self._markdown_table(drugs, ["TARGET_SYMBOL", "DRUG_NAME", "SOURCE", "DIRECTIONALITY_LABEL"]) if drugs is not None and not drugs.empty else ["UNAVAILABLE"])
        lines.extend(
            [
                "",
                "## STRING comparison",
                "",
                str(stage_reports.get("p6", {}).get("reason", stage_reports.get("p6", {}).get("jaccard", "UNAVAILABLE"))),
                "",
                "## Limitations",
                "",
                "- Background is Phase 1 evaluated genes, not whole-exome callable genes.",
                "- Pathway coverage, direction uncertainty, missing assays, and enrichment limitations remain.",
                "- HP terms are observational gene-annotation enrichment, not observed phenotype.",
                "",
                "## Next actions",
                "",
                "- Validate prioritized pathways with independent experimental evidence.",
                "- Verify literature and assay provenance before biological interpretation.",
            ]
        )
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    @staticmethod
    def _markdown_table(frame: pd.DataFrame, columns: list[str]) -> list[str]:
        available = [column for column in columns if column in frame]
        if not available:
            return ["UNAVAILABLE"]
        rows = ["| " + " | ".join(available) + " |", "| " + " | ".join(["---"] * len(available)) + " |"]
        for _, row in frame.loc[:, available].head(10).iterrows():
            rows.append("| " + " | ".join(str(row[column]).replace("|", "\\|") for column in available) + " |")
        return rows
