"""Filesystem contract for private MVA analysis artifacts."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class PipelinePaths:
    """Resolve repository-local paths without embedding a workstation path."""

    root: Path

    @classmethod
    def from_root(cls, root: str | Path) -> "PipelinePaths":
        return cls(Path(root).resolve())

    @property
    def data(self) -> Path:
        return self.root / "data"

    @property
    def interim(self) -> Path:
        return self.data / "interim"

    @property
    def processed(self) -> Path:
        return self.data / "processed"

    @property
    def references(self) -> Path:
        return self.data / "input" / "refs"

    @property
    def vep_input(self) -> Path:
        return self.interim / "vep" / "vep_high_moderate.csv"

    @property
    def stages(self) -> dict[str, Path]:
        return {
            "cadd": self.interim / "cadd" / "cadd_scores.csv",
            "gnomad": self.interim / "gnomad" / "gnomad_candidates.csv",
            "revel": self.interim / "revel" / "revel.csv",
            "alphamissense": self.interim / "alphamissense" / "alphamissense.csv",
            "spliceai": self.interim / "spliceai" / "spliceai.csv",
            "alphagenome": self.interim / "alphagenome" / "alphagenome_avi.csv",
            "exomiser": self.interim / "exomiser" / "exomiser_variants.csv",
        }

    @property
    def master_evidence(self) -> Path:
        return self.processed / "variant_evidence_master.csv"

    @property
    def prioritized_candidates(self) -> Path:
        return self.processed / "prioritized_candidates.csv"

    def make_output_directories(self) -> None:
        """Create only output directories, never reference or input directories."""
        self.processed.mkdir(parents=True, exist_ok=True)
        for path in self.stages.values():
            path.parent.mkdir(parents=True, exist_ok=True)
