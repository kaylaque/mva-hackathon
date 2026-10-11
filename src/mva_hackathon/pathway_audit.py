"""Aggregate-only validation reporting for pathway workflow artifacts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
from typing import Any

import pandas as pd


_SENSITIVITY_PROVENANCE_FIELDS = {
    "baseline_p0_sha256",
    "background_count",
    "query_count",
    "sources",
    "fdr_method",
    "min_passing",
    "generated_at",
    "cache_key",
    "outcome",
}


_FDR_METHOD = "fdr"


def _file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()



def _require_nonignored_paths(paths: list[Path]) -> None:
    repo_root = Path(__file__).parents[2]
    for path in paths:
        check = subprocess.run(
            ["git", "check-ignore", "--quiet", str(path.resolve())],
            cwd=repo_root,
            check=False,
            capture_output=True,
        )
        if check.returncode == 0:
            raise ValueError(f"ignored artifact path is forbidden: {path}")

def _read_csv(path: Path) -> pd.DataFrame:
    """Read required audit table with a path-specific error."""
    if not path.is_file():
        raise FileNotFoundError(f"missing pathway artifact: {path.name}")
    return pd.read_csv(path)


def _baseline_summary(output_dir: Path) -> dict[str, int]:
    p0 = _read_csv(output_dir / "p0_gene_qc.csv")
    p1 = _read_csv(output_dir / "p1_enrichment_results.csv")
    p2 = _read_csv(output_dir / "p2_prioritized_pathways.csv")
    q_values = pd.to_numeric(p1["BH_Q"], errors="coerce")
    raw_p = pd.to_numeric(p1["RAW_P"], errors="coerce")
    return {
        "background_genes": len(p0),
        "foreground_genes": int(p0["IN_FOREGROUND"].sum()),
        "mapped_genes": int(p0["STATUS"].eq("mapped").sum()),
        "unmapped_genes": int((~p0["STATUS"].eq("mapped")).sum()),
        "enrichment_terms": len(p1),
        "nominal_raw_p_lt_0_05": int(raw_p.lt(0.05).sum()),
        "provider_bh_le_0_05": int(q_values.le(0.05).sum()),
        "provider_bh_le_0_10": int(q_values.le(0.10).sum()),
        "prioritized_pathways": len(p2),
    }


def _sensitivity_summary(output_dir: Path, name: str, baseline_hash: str) -> dict[str, int | str]:
    path = output_dir / name
    manifest_path = output_dir / f"{path.stem}.provenance.json"
    if not path.is_file() or not manifest_path.is_file():
        return {"status": "UNAVAILABLE", "reason": "MISSING_PROVENANCE"}
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (
        not _SENSITIVITY_PROVENANCE_FIELDS.issubset(manifest)
        or manifest["baseline_p0_sha256"] != baseline_hash
        or manifest["fdr_method"] != _FDR_METHOD
    ):
        return {"status": "UNAVAILABLE", "reason": "STALE_OR_INVALID_PROVENANCE"}
    frame = _read_csv(path)
    return {"status": "AVAILABLE", "prioritized_pathways": len(frame), "outcome": str(manifest["outcome"])}

def audit_pathway_outputs(output_dir: Path) -> tuple[dict[str, Any], Path]:
    """Create local aggregate report from provenance-bound, non-ignored artifacts."""
    report_path = output_dir / "pathway_validation_report.md"
    _require_nonignored_paths(
        [
            output_dir / name
            for name in (
                "p0_gene_qc.csv",
                "p1_enrichment_results.csv",
                "p2_prioritized_pathways.csv",
                "p6_evaluation_summary.csv",
                "p2_sensitivity_min_passing_2.csv",
                "p2_sensitivity_min_passing_2.provenance.json",
                "p2_sensitivity_hp_excluded.csv",
                "p2_sensitivity_hp_excluded.provenance.json",
                report_path.name,
            )
        ]
    )
    baseline_hash = _file_sha256(output_dir / "p0_gene_qc.csv")
    baseline = _baseline_summary(output_dir)
    threshold_two = _sensitivity_summary(output_dir, "p2_sensitivity_min_passing_2.csv", baseline_hash)
    hp_excluded = _sensitivity_summary(output_dir, "p2_sensitivity_hp_excluded.csv", baseline_hash)
    p6 = _read_csv(output_dir / "p6_evaluation_summary.csv")
    comparison = {key: int(p6.iloc[0][key]) for key in ("pathway_only_targets", "string_only_targets", "shared_targets") if key in p6}
    summary = {
        "baseline": baseline,
        "threshold_two": threshold_two,
        "hp_excluded": hp_excluded,
        "comparison": comparison,
        "privacy": {
            "transmitted_identifier_types": ["HGNC symbol", "ENSG", "UniProt accession", "source target ID"],
            "excluded_patient_fields": ["variants", "genotypes", "coordinates", "phenotype narratives", "patient identifiers"],
        },
    }
    report_path = output_dir / "pathway_validation_report.md"
    lines = [
        "# Pathway validation report",
        "",
        "Research-only. Enrichment and database co-occurrence do not establish therapeutic efficacy.",
        "",
        "## Baseline enrichment QC",
        "",
        f"- Background genes: {baseline['background_genes']}",
        f"- Foreground genes: {baseline['foreground_genes']}",
        f"- ENSG mapped/unmapped: {baseline['mapped_genes']}/{baseline['unmapped_genes']}",
        f"- GOSt terms: {baseline['enrichment_terms']}",
        f"- Nominal raw-p < 0.05: {baseline['nominal_raw_p_lt_0_05']}",
        f"- Provider BH-FDR ≤ 0.05: {baseline['provider_bh_le_0_05']}",
        "",
        "No baseline term passed provider BH-FDR ≤ 0.05." if baseline["provider_bh_le_0_05"] == 0 else "Baseline provider BH-significant terms exist.",
        "",
        "## Sensitivity",
        "",
        _sensitivity_line("Threshold-2", threshold_two),
        _sensitivity_line("HP-excluded", hp_excluded),
        "- Gene-length bias: UNMEASURED; no documented length covariate was supplied.",
        *_sensitivity_table(output_dir, "p2_sensitivity_min_passing_2.csv", threshold_two),
        "",
        *_sensitivity_table(output_dir, "p2_sensitivity_hp_excluded.csv", hp_excluded, "HP-excluded prioritized pathways"),
        "",
        "## Downstream evidence applicability",
        "",
        *(
            [
                "P3 mechanism, P4 target expansion, and P5 drug validation: NOT_APPLICABLE_NO_ELIGIBLE_INPUT.",
                "No baseline pathway met provider BH-FDR ≤ 0.05, so no pathway-expanded targets were eligible for Reactome/UniProt/QuickGO or ChEMBL/DGIdb/Open Targets validation.",
            ]
            if baseline["prioritized_pathways"] == 0
            else [
                "Baseline pathways were eligible for downstream evidence assessment.",
                "Review P3/P4/P5 source status and source-native records before interpreting any target or drug hypothesis.",
            ]
        ),
        "Source failure versus source-empty semantics remain a required pipeline correction before any non-empty P5 interpretation.",
        "## STRING comparison",
        "",
        f"- Pathway-only targets: {comparison.get('pathway_only_targets', 'UNAVAILABLE')}",
        f"- STRING-only targets: {comparison.get('string_only_targets', 'UNAVAILABLE')}",
        f"- Shared targets: {comparison.get('shared_targets', 'UNAVAILABLE')}",
        "",
        "## Interpretation limits",
        "",
        "- Drug comparison: UNAVAILABLE. Current P6 consumes target-panel overlap only; it does not consume network drug evidence or retrieve union-panel evidence.",
        *(
            ["- No baseline pathway target or drug hypothesis is eligible when no provider BH-significant pathway exists."]
            if baseline["prioritized_pathways"] == 0
            else []
        ),
        "- Threshold-2 findings are sensitivity results, not confirmation of baseline enrichment.",
        "- ChEMBL activity/assay evidence remains unavailable unless retrieved with source-native records.",
        "- Directional pharmacology does not prove biological rescue.",
        "",
        "## Privacy inventory",
        "",
        "- External identifiers: HGNC symbols, ENSG, UniProt accessions, source target IDs.",
        "- Never external: variants, genotypes, coordinates, phenotype narratives, patient identifiers.",
    ]
    temporary = report_path.with_suffix(".md.tmp")
    temporary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    temporary.replace(report_path)
    return summary, report_path




def _sensitivity_table(output_dir: Path, name: str, summary: dict[str, int | str], heading: str = "Threshold-2 prioritized pathways") -> list[str]:
    """Render only provenance-bound local pathway evidence."""
    lines = [f"## {heading}", ""]
    if summary["status"] != "AVAILABLE":
        return [*lines, f"UNAVAILABLE: {summary['reason']}"]
    frame = _read_csv(output_dir / name)
    columns = [column for column in ("TERM_ID", "TERM_NAME", "SOURCE", "BH_Q", "FOLD_ENRICHMENT", "CONTRIBUTING_GENES") if column in frame]
    if frame.empty or not columns:
        return [*lines, "No prioritized pathways."]
    lines.extend([
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ])
    for _, row in frame.loc[:, columns].head(15).iterrows():
        lines.append("| " + " | ".join(str(row[column]).replace("|", "\\|") for column in columns) + " |")
    return lines
def _sensitivity_line(label: str, summary: dict[str, int | str]) -> str:
    """Format aggregate sensitivity availability without pathway identities."""
    if summary["status"] == "UNAVAILABLE":
        return f"- {label}: UNAVAILABLE ({summary['reason']})"
    count = summary["prioritized_pathways"]
    noun = "pathway" if count == 1 else "pathways"
    return f"- {label} sensitivity produced {count} prioritized {noun}."
