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
│   ├── README.md                     # private-data contract and stage schemas
│   ├── input/
│   │   ├── patient/                  # source patient data (ignored)
│   │   └── refs/                     # large public reference assets (ignored)
│   ├── interim/                      # resumable private stage artifacts (ignored)
│   └── processed/                    # consolidated evidence tables (ignored)
├── envs/                             # conda environments for VCF/VEP tooling
├── scripts/
│   ├── 01_run_qc.sh                  # bcftools QC + filtering
│   ├── 02_run_snpeff.sh              # snpEff annotation
│   ├── 03_run_ensembl_vep.sh         # offline Ensembl VEP annotation
│   ├── download_references.sh        # resumable public-reference downloads
│   ├── run_enrichment.py/.sh         # merge completed enrichment stage outputs
│   └── run_exomiser.py/.sh           # prepare/run/import Exomiser results
├── src/mva_hackathon/                # typed pipeline and Exomiser objects
└── pyproject.toml / uv.lock
```

## WGS variant annotation and enrichment pipeline

The established QC → snpEff → VEP workflow emits private intermediate files
under `data/interim/`. The enrichment layer consumes the latest VEP
high/moderate table at `data/interim/vep/vep_high_moderate.csv` and completed
candidate-level annotation tables for CADD, gnomAD, REVEL, AlphaMissense,
SpliceAI, AlphaGenome, and optionally Exomiser.

```bash
# Public reference downloads only; never downloads patient data.
bash scripts/download_references.sh

# Merge completed stage outputs. This does not launch CADD, gnomAD, or other annotators.
uv run python scripts/run_enrichment.py --vep /private/path/to/current-chunked-output.csv

# Prepare Exomiser candidate VCF/YAML only.
uv run python scripts/run_exomiser.py \
  --jar /path/to/exomiser-cli.jar \
  --hpo HP:0000001

# Execute Exomiser only when intentionally requested.
uv run python scripts/run_exomiser.py \
  --jar /path/to/exomiser-cli.jar \
  --hpo HP:0000001 \
  --execute
```

Read [`data/README.md`](data/README.md) before placing or naming any private
artifact. `variant_evidence_master.csv` and `prioritized_candidates.csv` are
research-prioritization outputs; they are not diagnostic classifications.
### Existing annotation setup

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

## Analysis resources

All source data and reference databases use GRCh38. Downloaded assets remain
local and ignored; exact public download endpoints used by this repository are
in [`scripts/download_references.sh`](scripts/download_references.sh).

| Resource | Role in this repository | Link |
| --- | --- | --- |
| SageBio MVA Hackathon | Case background and source dataset | [Space](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026), [dataset](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data) |
| Ensembl / VEP | Reference FASTA, offline VEP cache, transcript annotation | [Ensembl](https://www.ensembl.org/), [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html), [release 116 downloads](https://ftp.ensembl.org/pub/release-116/) |
| ClinVar | Clinical-variant annotations | [NCBI ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/), [GRCh38 VCF](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/) |
| gnomAD v4 | Population allele frequency and rarity | [Browser](https://gnomad.broadinstitute.org/), [GraphQL API](https://gnomad.broadinstitute.org/api), [v4.1 joint VCF](https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/) |
| CADD v1.7 | Genome-wide deleteriousness score | [CADD](https://cadd.gs.washington.edu/), [GRCh38 API](https://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.7), [precomputed scores](https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/) |
| REVEL v1.3 | Missense-specific ensemble score | [Zenodo archive](https://zenodo.org/records/7072866) |
| AlphaMissense | Missense pathogenicity score and class | [Project](https://deepmind.google/technologies/alphamissense/), [GRCh38 table](https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz) |
| SpliceAI | Predicted splice donor/acceptor disruption | [GitHub](https://github.com/Illumina/SpliceAI), [GRCh38 reference](https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz) |
| AlphaGenome | Regulatory and molecular-effect evidence | [AlphaGenome](https://deepmind.google/science/alphagenome/) |
| Exomiser | HPO and disease-model phenotype ranking | [Exomiser](https://github.com/exomiser/Exomiser), [data downloads](https://data.monarchinitiative.org/exomiser/) |
| bcftools / samtools | VCF QC, filtering, indexing, and reference indexing | [bcftools](https://samtools.github.io/bcftools/), [samtools](https://www.htslib.org/) |
| SnpEff / SnpSift | Secondary functional annotation and field extraction | [SnpEff](https://pcingola.github.io/SnpEff/) |
| Python / uv stack | Environment and VCF/table processing | [uv](https://docs.astral.sh/uv/), [cyvcf2](https://brentp.github.io/cyvcf2/), [pysam](https://pysam.readthedocs.io/) |

Do not treat any computational predictor, population-frequency annotation, or
rank as diagnostic evidence. Manual review, phenotype fit, inheritance, and
the hackathon's research-only scope remain required.

## Project guidance

See [`AGENTS.md`](AGENTS.md) for repository architecture, data-handling rules, VCF analysis conventions, and QA guidance.
