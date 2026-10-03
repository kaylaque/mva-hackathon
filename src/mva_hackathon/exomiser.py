"""Create Exomiser inputs and normalize its tabular variant output."""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .config import PipelinePaths
from .enrichment import EvidencePipeline


@dataclass(frozen=True)
class ExomiserRun:
    jar: Path | None
    phenotype_terms: tuple[str, ...]
    assembly: str = "hg38"


class ExomiserRunner:
    """Own all Exomiser-specific file formats; do not run it from notebooks."""

    def __init__(self, paths: PipelinePaths, run: ExomiserRun) -> None:
        self.paths = paths
        self.run = run
        self.output_dir = paths.stages["exomiser"].parent

    def write_input(self, variants: pd.DataFrame) -> tuple[Path, Path]:
        """Write minimal VCF and YAML config from the one-row-per-variant evidence table."""
        self.output_dir.mkdir(parents=True, exist_ok=True)
        vcf_path = self.output_dir / "candidates.vcf"
        yaml_path = self.output_dir / "analysis.yml"
        with vcf_path.open("w", encoding="utf-8") as handle:
            handle.write("##fileformat=VCFv4.2\n##reference=GRCh38\n")
            handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for row in variants.itertuples(index=False):
                handle.write(f"{row.CHROM}\t{int(row.POS)}\t.\t{row.REF}\t{row.ALT}\t.\tPASS\t.\n")
        terms = "\n".join(f"  - {term}" for term in self.run.phenotype_terms)
        yaml_path.write_text(
            f"analysis:\n"
            f"  genomeAssembly: {self.run.assembly}\n"
            f"  vcf: {vcf_path}\n"
            f"  hpoIds:\n{terms}\n"
            f"  outputPrefix: {self.output_dir / 'exomiser'}\n",
            encoding="utf-8",
        )
        return vcf_path, yaml_path

    def command(self, yaml_path: Path) -> list[str]:
        if self.run.jar is None:
            raise ValueError("An Exomiser JAR is required to execute an analysis")
        return ["java", "-jar", str(self.run.jar), "--analysis", str(yaml_path)]

    def execute(self, yaml_path: Path) -> None:
        """Run only on explicit CLI request; all long-running work is opt-in."""
        subprocess.run(self.command(yaml_path), check=True)

    def normalize_results(self, source: Path) -> pd.DataFrame:
        """Retain stable score fields from a user-supplied Exomiser TSV/CSV export."""
        separator = "\t" if source.suffix.lower() in {".tsv", ".txt"} else ","
        table = pd.read_csv(source, sep=separator, low_memory=False)
        aliases = {
            "variantScore": "EXOMISER_VARIANT_SCORE",
            "geneScore": "EXOMISER_GENE_SCORE",
            "phenotypeScore": "EXOMISER_PHENOTYPE_SCORE",
            "rank": "EXOMISER_RANK",
            "chrom": "CHROM",
            "pos": "POS",
            "ref": "REF",
            "alt": "ALT",
        }
        table = table.rename(columns={key: value for key, value in aliases.items() if key in table.columns})
        required = {"CHROM", "POS", "REF", "ALT"}
        if not required.issubset(table.columns):
            raise ValueError(f"{source} lacks genomic columns required for a variant-level Exomiser merge")
        table["VARIANT_KEY"] = EvidencePipeline._variant_key(table)
        output = self.paths.stages["exomiser"]
        output.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(output, index=False)
        return table
