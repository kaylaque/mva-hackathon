# Project Overview — MVA Hackathon 2026

A reference for the repo layout, the analysis methods, and the case behind this work. For quick-start commands see the root [README.md](../README.md).

## The project

**Rare Disease, Real Kid: MVA Hackathon 2026** ([Hugging Face Space](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026), [dataset](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data)) centers on the genome and clinical story of a real child living with **Mosaic Variegated Aneuploidy (MVA)** — an ultra-rare genetic condition affecting fewer than 50 people worldwide, with no established treatment; care today means managing symptoms. The family has opened the case to the research community hoping someone can find an answer.

Participants receive whole-genome sequencing (WGS) data for the child plus a clinical phenotype document, and analyze them to surface candidate variants and genes that could explain the condition. **This is research exploration, not diagnosis or medical advice.**

## Project structure

```
.
├── data/
│   ├── raw/                          # source data (ignored): WGS VCF + tabix index, clinical phenotype DOCX
│   └── interim/                      # pipeline intermediates (ignored by default)
│       ├── qc/                       # bcftools QC summaries            (committed - curated)
│       ├── filtered/                 # PASS-quality filtered VCF        (ignored)
│       ├── snpeff/                   # snpEff annotation outputs        (ignored)
│       │   └── filtered_csv/         # snpEff candidate tables          (committed - curated)
│       └── vep/                      # VEP annotation outputs           (ignored)
│           └── filtered_csv/         # VEP candidate tables             (committed - curated)
├── envs/                             # conda environments for pipeline tools (mva_qc, mva_snpeff, mva_vep)
├── scripts/
│   ├── 01_run_qc.sh                  # bcftools QC + PASS-quality filtering
│   ├── 02_run_snpeff.sh              # snpEff functional annotation
│   ├── 03_run_ensembl_vep.sh         # offline Ensembl VEP annotation
│   └── py/
│       ├── filter_variants_snpeff.py # reduce snpEff output to candidate CSVs
│       └── filter_variants_vep.py    # reduce VEP output to candidate CSVs
├── mutation_analysis.ipynb           # variant exploration notebook
├── quality_check_analysis.ipynb      # QC review notebook
└── pyproject.toml / uv.lock          # Python 3.13 environment for notebooks
```

Only small, deliberately curated outputs are committed. Raw and large derived data (anything over 8 MB) never enter Git history — enforced by `.gitignore` and a pre-commit hook.

## Methods brief

1. **QC and quality filtering** (`scripts/01_run_qc.sh`, env `mva_qc`): bcftools profiles the raw WGS VCF (`data/raw/WGS_EX2312012_HGWCNDSX7.vcf.gz`) and keeps PASS-quality variants, producing `data/interim/filtered/pass_quality.vcf.gz` plus the QC summaries in `data/interim/qc/`.
2. **snpEff annotation** (`scripts/02_run_snpeff.sh`, env `mva_snpeff`): annotates the filtered VCF against `GRCh38.mane.1.0.ensembl`; `scripts/py/filter_variants_snpeff.py` reduces it to `all_high_moderate_variants.csv`, `recessive_candidate_genes.csv` (homozygous / compound heterozygous patterns relevant to recessive disease), and `phenotype_focused_variants.csv` in `data/interim/snpeff/filtered_csv/`.
3. **Ensembl VEP annotation** (`scripts/03_run_ensembl_vep.sh`, env `mva_vep`): offline VEP run on GRCh38; `scripts/py/filter_variants_vep.py` produces high/moderate impact variants, recessive candidates, pathogenic variants, ClinVar-backed calls, phenotype-focused variants, and top clinical rankings in `data/interim/vep/filtered_csv/`.

The two annotators are run in parallel so their candidate lists can be cross-checked; agreement between snpEff and VEP prioritizes the shortlist for manual review against the child's phenotype document.
