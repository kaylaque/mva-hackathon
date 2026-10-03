"""Variant-level evidence consolidation for existing enrichment stage outputs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd

from .config import PipelinePaths


@dataclass(frozen=True)
class EvidenceStage:
    """One optional annotation table and the fields retained from it."""

    name: str
    path: Path
    columns: tuple[str, ...]


class EvidencePipeline:
    """Build one variant-level evidence table from VEP and completed stages."""

    REQUIRED_VEP_COLUMNS = frozenset({"CHROM", "POS", "REF", "ALT"})
    STAGE_COLUMNS = {
        "cadd": ("CADD_RAW", "CADD_PHRED", "CADD_STATUS", "CADD_SOURCE"),
        "gnomad": (
            "GNOMAD_STATUS", "GNOMAD_AF", "GNOMAD_AC", "GNOMAD_AN",
            "GNOMAD_NHOMALT", "GNOMAD_GRPMAX_AF", "GNOMAD_RARITY",
            "GNOMAD_RARE_FLAG",
        ),
        "revel": ("REVEL_SCORE",),
        "alphamissense": (
            "AM_SCORE", "AM_CLASS", "AM_UNIPROT", "AM_TRANSCRIPT",
            "AM_PROTEIN_VARIANT",
        ),
        "spliceai": (
            "SPLICEAI_GENE", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG",
            "DP_AL", "DP_DG", "DP_DL", "SPLICEAI_MAX_DS",
        ),
        "alphagenome": (
            "ALPHAGENOME_AVI_RAW", "ALPHAGENOME_AVI_QUANTILE",
            "ALPHAGENOME_AVI_PHRED", "ALPHAGENOME_TOP_FEATURE",
            "ALPHAGENOME_TOP_FEATURE_VALUE", "ALPHAGENOME_STATUS",
        ),
        "exomiser": (
            "EXOMISER_VARIANT_SCORE", "EXOMISER_GENE_SCORE",
            "EXOMISER_PHENOTYPE_SCORE", "EXOMISER_RANK",
        ),
    }

    def __init__(self, paths: PipelinePaths) -> None:
        self.paths = paths

    @staticmethod
    def _chromosome(value: object) -> str:
        chrom = str(value).strip()
        chrom = chrom.removeprefix("chr").removeprefix("CHR")
        return {"23": "X", "24": "Y"}.get(chrom, chrom)

    @classmethod
    def _variant_key(cls, frame: pd.DataFrame) -> pd.Series:
        positions = pd.to_numeric(frame["POS"], errors="coerce")
        valid = positions.notna() & frame["REF"].notna() & frame["ALT"].notna()
        result = pd.Series(pd.NA, index=frame.index, dtype="string")
        result.loc[valid] = (
            frame.loc[valid, "CHROM"].map(cls._chromosome)
            + ":" + positions.loc[valid].astype(int).astype(str)
            + ":" + frame.loc[valid, "REF"].astype(str).str.upper()
            + ":" + frame.loc[valid, "ALT"].astype(str).str.upper()
        )
        return result
    @staticmethod
    def _numeric_column(frame: pd.DataFrame, name: str) -> pd.Series:
        return pd.to_numeric(
            frame.get(name, pd.Series(float("nan"), index=frame.index)),
            errors="coerce",
        )


    def load_vep(self, path: Path | None = None) -> pd.DataFrame:
        source = path or self.paths.vep_input
        frame = pd.read_csv(source, low_memory=False)
        frame.columns = frame.columns.str.strip()
        missing = self.REQUIRED_VEP_COLUMNS.difference(frame.columns)
        if missing:
            raise ValueError(f"{source} is missing required VEP columns: {sorted(missing)}")
        frame["VARIANT_KEY"] = self._variant_key(frame)
        if frame["VARIANT_KEY"].isna().any():
            raise ValueError(f"{source} contains rows without a usable variant key")
        canonical = frame.get("CANONICAL", pd.Series("", index=frame.index))
        priority = self._numeric_column(frame, "PRIORITY_SCORE")
        frame = frame.assign(
            _canonical=canonical.astype(str).str.upper().isin({"YES", "Y", "TRUE", "1"}),
            _priority=priority,
        )
        return (
            frame.sort_values(["VARIANT_KEY", "_canonical", "_priority"], ascending=[True, False, False])
            .drop_duplicates("VARIANT_KEY", keep="first")
            .drop(columns=["_canonical", "_priority"])
            .reset_index(drop=True)
        )

    def available_stages(self) -> Iterable[EvidenceStage]:
        for name, columns in self.STAGE_COLUMNS.items():
            path = self.paths.stages[name]
            if path.is_file():
                yield EvidenceStage(name, path, columns)

    def merge_completed_stages(self, variants: pd.DataFrame) -> pd.DataFrame:
        evidence = variants.copy()
        for stage in self.available_stages():
            annotation = pd.read_csv(stage.path, low_memory=False)
            annotation.columns = annotation.columns.str.strip()
            if "VARIANT_KEY" not in annotation:
                if self.REQUIRED_VEP_COLUMNS.issubset(annotation.columns):
                    annotation["VARIANT_KEY"] = self._variant_key(annotation)
                else:
                    raise ValueError(f"{stage.path} has no VARIANT_KEY or genomic-coordinate columns")
            keep = ["VARIANT_KEY", *(field for field in stage.columns if field in annotation)]
            annotation = annotation[keep].drop_duplicates("VARIANT_KEY", keep="first")
            evidence = evidence.merge(annotation, on="VARIANT_KEY", how="left", validate="one_to_one")
        return evidence

    @staticmethod
    def _rarity(af: object, gnomad_status: object) -> str:
        if not pd.isna(gnomad_status) and gnomad_status == "not_found":
            return "ABSENT"
        if pd.isna(af):
            return "UNKNOWN"
        if af <= 1e-5:
            return "ULTRA_RARE"
        if af <= 1e-4:
            return "VERY_RARE"
        if af <= 1e-3:
            return "RARE"
        if af <= 1e-2:
            return "LOW_FREQUENCY"
        return "COMMON"

    def prioritize(self, evidence: pd.DataFrame) -> pd.DataFrame:
        """Add transparent research triage flags; this is not ACMG classification."""
        result = evidence.copy()
        population_af = self._numeric_column(result, "AF_POPULATION")
        gnomad_af = self._numeric_column(result, "GNOMAD_AF")
        result["BEST_AF"] = gnomad_af.combine_first(population_af)
        gnomad_status = result.get("GNOMAD_STATUS", pd.Series(pd.NA, index=result.index))
        result["BEST_RARITY"] = [self._rarity(af, status) for af, status in zip(result["BEST_AF"], gnomad_status)]
        rare = result["BEST_AF"].le(0.001) | gnomad_status.eq("not_found")
        cadd_high = self._numeric_column(result, "CADD_PHRED").ge(20)
        splice_high = self._numeric_column(result, "SPLICEAI_MAX_DS").ge(0.2)
        revel_high = self._numeric_column(result, "REVEL_SCORE").ge(0.5)
        am_high = self._numeric_column(result, "AM_SCORE").ge(0.5)
        result["RARE_AND_CADD_HIGH"] = rare & cadd_high
        result["PRIORITIZATION_SIGNALS"] = rare.astype(int) + cadd_high.astype(int) + splice_high.astype(int) + revel_high.astype(int) + am_high.astype(int)
        return result.sort_values("PRIORITIZATION_SIGNALS", ascending=False, kind="stable")

    def consolidate(self, vep_path: Path | None = None) -> pd.DataFrame:
        self.paths.make_output_directories()
        evidence = self.prioritize(self.merge_completed_stages(self.load_vep(vep_path)))
        evidence.to_csv(self.paths.master_evidence, index=False)
        evidence[evidence["PRIORITIZATION_SIGNALS"] > 0].to_csv(self.paths.prioritized_candidates, index=False)
        return evidence
