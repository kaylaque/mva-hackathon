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

## Project structure

```text
mva-hackathon/
├── data/
│   ├── raw/                          # source data (ignored): WGS VCF + tabix index, clinical phenotype DOCX
│   └── interim/                      # pipeline intermediates (ignored by default)
│       ├── qc/                       # bcftools QC summaries            (committed - curated)
│       ├── filtered/                 # PASS-quality filtered VCF        (ignored)
│       ├── snpeff/                   # snpEff annotation outputs        (ignored)
│       │   └── filtered_csv/         # snpEff candidate tables          (committed - curated)
│       └── vep/                      # VEP annotation outputs           (ignored)
│           └── filtered_csv/         # VEP candidate tables             (committed - curated)
├── docs/
│   └── README.md                     # project overview, methods brief, case background
├── envs/                             # conda environments: mva_qc, mva_snpeff, mva_vep
├── scripts/
│   ├── 01_run_qc.sh                  # bcftools QC + PASS-quality filtering
│   ├── 02_run_snpeff.sh              # snpEff functional annotation
│   ├── 03_run_ensembl_vep.sh         # offline Ensembl VEP annotation
│   └── py/
│       ├── filter_variants_snpeff.py # reduce snpEff output to candidate CSVs
│       └── filter_variants_vep.py    # reduce VEP output to candidate CSVs
├── mutation_analysis.ipynb           # variant exploration notebook
├── quality_check_analysis.ipynb      # QC review notebook
├── AGENTS.md                         # repo conventions and data-handling rules
└── pyproject.toml / uv.lock          # Python 3.13 environment for notebooks
```

Raw and large derived data (anything over 8 MB) never enter Git history — enforced by `.gitignore` and a pre-commit hook.

## WGS variant annotation pipeline

`scripts/` contains a three-step bcftools -> snpEff -> Ensembl VEP workflow for the WGS VCF, driven by conda environments in `envs/` (bioconda tools that the uv environment does not provide):

```bash
conda env create -f envs/mva_qc.yml      # bcftools, htslib
conda env create -f envs/mva_snpeff.yml  # snpEff, SnpSift, OpenJDK 17
conda env create -f envs/mva_vep.yml     # ensembl-vep, htslib

bash scripts/01_run_qc.sh             # QC + PASS-quality filtering -> data/interim/filtered/
bash scripts/02_run_snpeff.sh         # snpEff annotation (GRCh38.mane.1.0.ensembl) -> data/interim/snpeff/
bash scripts/03_run_ensembl_vep.sh    # offline VEP annotation (GRCh38) -> data/interim/vep/
```

`scripts/py/filter_variants_snpeff.py` and `scripts/py/filter_variants_vep.py` reduce the annotated outputs to compact candidate CSVs under `data/interim/*/filtered_csv/`. All raw and interim genomic artifacts stay outside Git per `.gitignore`; keep every committed file below 8 MB (enforced by the pre-commit hook).

## Project guidance

See [`AGENTS.md`](AGENTS.md) for repository architecture, data-handling rules, VCF analysis conventions, and QA guidance.
