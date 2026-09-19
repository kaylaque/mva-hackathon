# mva-hackathon

Early technical exploration of whole-genome sequencing (WGS) and VCF data for the SageBio **Rare Disease, Real Kid: MVA Hackathon 2026**. The project profiles VCF, BAM/CRAM, and FASTQ inputs, summarizes quality and genotype signals, and produces exploratory plots and compact tables. Results support research and candidate discovery only; they are not diagnostic or medical advice.

## Installation

Requires Python 3.13 and [`uv`](https://docs.astral.sh/uv/).

```bash
uv sync
brew install bcftools samtools
uv run python -m ipykernel install --user \
  --name mva-hackathon \
  --display-name "Python (mva-hackathon)"
```

The Homebrew command installs the local HTS utilities used for VCF and alignment inspection. The Python environment includes `cyvcf2`, `pysam`, `scikit-allel`, `numpy`, `pandas`, `matplotlib`, `seaborn`, `plotly`, `upsetplot`, JupyterLab, and related project dependencies.

## Run the exploration notebook

```bash
uv run jupyter lab MVA_WGS_VCF_Early_Exploration.ipynb
```

Use the `Python (mva-hackathon)` kernel. On Kaggle, the notebook reads from `/kaggle/input` and writes to `/kaggle/working`. Locally, place source files under `input/`; outputs are written to `results/`.

Supported inputs include `.vcf`, `.vcf.gz`, `.bam`, `.cram`, `.fastq`, and `.fastq.gz`. Large patient/genomic files are intentionally excluded from Git by `.gitignore`; do not force-add them.

## Project guidance

See [`AGENTS.md`](AGENTS.md) for repository architecture, data-handling rules, VCF analysis conventions, and QA guidance.
