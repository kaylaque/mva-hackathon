# Repository Guidelines

## Project Overview

`mva-hackathon` is a greenfield Python workspace for exploring the SageBio **Rare Disease, Real Kid: MVA Hackathon 2026** case. The repository currently contains no application code, entry point, scripts, notebooks, or tests. The upstream case is a real-child rare-disease research dataset; treat outputs as research exploration, not diagnosis or medical advice.

The upstream Hugging Face dataset is `SageBio/mva-hackathon-2026-data`. It includes a compressed WGS VCF plus tabix index, a clinical-phenotype DOCX, and very large FASTQ files. Do not commit downloaded patient/genomic data to this repository.

## Architecture & Data Flow

No implementation architecture exists yet. When adding code, keep the first version small and explicit:

1. Fetch or receive source data from Hugging Face into ignored local data storage.
2. Read indexed VCF/BCF regions with `cyvcf2` for fast variant iteration or `pysam` for HTSlib/tabix access.
3. Normalize/filter records into `pandas` tables or `scikit-allel` arrays; preserve chromosome, position, alleles, genotype, quality, and filter information.
4. Use `plotly` for exploratory visualizations and `gradio` for an interactive surface if one is needed.
5. Write reproducible summaries to an ignored `results/` directory or small, reviewed fixtures only.

Prefer a thin pipeline: data access, transformation/filtering, and presentation should be separate functions/modules. Keep large-file access region-based; never load the whole WGS VCF into memory by default. Preserve the source build/reference assembly and sample identifiers in outputs.

## Key Directories

The repository currently contains one root exploratory notebook; use this layout as code is introduced:

- `MVA_WGS_VCF_Early_Exploration.ipynb` — Kaggle/local early profiling of VCF, BAM/CRAM, and FASTQ inputs.
- `src/mva_hackathon/` — reusable VCF exploration and domain logic.
- `scripts/` — small command-line workflows and one-off data preparation.
- `notebooks/` — additional exploratory analysis only; move repeatable logic into `src/`.
- `data/raw/` — local downloads and source artifacts; keep ignored.
- `data/processed/` and `results/` — derived local artifacts; keep ignored unless deliberately curated.
- `tests/` — focused behavioral tests with tiny synthetic VCF fixtures.

## Development Commands

The project uses `uv` and Python 3.13. The environment and lockfile are already prepared:

```bash
uv sync                         # create/update .venv from pyproject.toml and uv.lock
source .venv/bin/activate       # optional; uv run does not require activation
uv run python --version
uv run python -c "import cyvcf2, pysam, allel, pandas, plotly, gradio, huggingface_hub, docx; print('VCF stack ready')"
uv run python -m ipykernel install --user --name mva-hackathon --display-name "Python (mva-hackathon)"
uv run jupyter lab MVA_WGS_VCF_Early_Exploration.ipynb
```

The notebook uses `/kaggle/input` and `/kaggle/working` on Kaggle. Locally it falls back to `data/raw/` and `results/mva_early_insights/`. The notebook setup cell installs missing Python-only packages in Kaggle and reports missing HTS binaries instead of running `apt-get` on macOS.

There is currently no application entry point or documented `run` command. Once an app or script exists, run it through `uv run`, for example `uv run python scripts/explore_vcf.py`. Add the command to the README when established. Add dependencies with `uv add <package>` and refresh the lockfile with `uv lock`/`uv sync`; commit `pyproject.toml`, `.python-version`, and `uv.lock`.

## Code Conventions & Common Patterns

No source conventions are established yet. Follow these defaults:

- Use typed, small functions with `snake_case` names; use `PascalCase` only for classes.
- Keep VCF parsing lazy and region-aware. Accept paths/regions as arguments rather than hard-coding patient data.
- Separate I/O from filtering and analysis so filters can be tested with synthetic records or tables.
- Return structured records/dataframes from analysis functions; reserve printing and UI updates for CLI/UI boundaries.
- Fail loudly on missing indexes, incompatible reference builds, malformed records, or ambiguous sample selection. Include the input path/region in errors; do not silently skip variants.
- Prefer synchronous iteration for local indexed VCF scans. Use async only at network/UI boundaries, and do not hide blocking HTSlib work inside async code.
- Preserve missing genotype values and allele numbering; do not coerce missing values to reference calls.
- Keep notebooks and UI callbacks thin. Put reusable behavior in `src/mva_hackathon/`.

## Important Files

- `pyproject.toml` — project metadata and direct dependencies.
- `uv.lock` — reproducible resolved dependency set; keep it under version control.
- `.python-version` — pins the workspace to Python 3.13, matching the upstream Gradio Space runtime.
- `MVA_WGS_VCF_Early_Exploration.ipynb` — early QC/profile workflow; streams large inputs and writes compact outputs under the environment-specific output directory.
- `README.md` — currently a placeholder; update it with the first runnable workflow and data-handling notes.
- `.gitignore` — ignores `.venv`, Python caches, test/coverage caches, and common generated artifacts.
- `.gitattributes` — currently only enables text normalization; it has no genomic binary/LFS policy.

## Runtime/Tooling Preferences

- Required runtime: Python 3.13 (`uv run` or `.venv/bin/python`).
- Package manager/environment: `uv`; do not introduce Poetry, Pipenv, Conda, or a second lockfile without a clear project-level reason.
- VCF/HTS tooling: `cyvcf2` for efficient iteration, `pysam` for HTSlib/tabix operations, and `scikit-allel` for array-oriented analysis.
- Notebook analysis: `numpy`, `pandas`, `matplotlib`, `seaborn`, and `upsetplot`; JupyterLab and `ipykernel` are development dependencies.
- System HTS tooling: `bcftools`, `samtools`, and `tabix`; on macOS install with `brew install bcftools samtools`.
- Data/presentation: `pandas`, `plotly`, `gradio`; `huggingface-hub` for dataset access; `python-docx` for the clinical phenotype document.
- Large inputs (VCF/BCF/BAM/FASTQ) belong outside Git. Use local ignored data paths or an approved artifact store; do not add the upstream FASTQ files or WGS VCF to commits. If genomic binaries ever need versioning, establish Git LFS rules before adding them.

## Testing & QA

No tests, fixtures, test runner, CI, coverage threshold, or lint configuration exists yet. Do not claim a test command is available until it is added. When implementation begins:

- Add `tests/` with tiny synthetic VCF/BCF fixtures, not patient data.
- Test observable behavior: indexed region boundaries, multiallelic records, missing genotypes, sample selection, malformed input errors, and filter precedence.
- Run focused checks with `uv run pytest`; run the full suite before delivery once a suite exists.
- Add `ruff`/`mypy` configuration only when the codebase needs it, and run them through `uv run`.
- At minimum, smoke-test imports and one small local VCF scan after dependency changes:

```bash
uv run python -c "from cyvcf2 import VCF; print('cyvcf2 import ok')"
uv lock --check
```

Do not download the multi-gigabyte upstream FASTQ files for a routine smoke test. Any analysis using the real child’s data must keep provenance, access controls, and the hackathon’s non-diagnostic scope explicit.
