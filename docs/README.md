# Project Overview — MVA Hackathon 2026

A reference for the repo layout, the analysis methods, and the case behind this work. For quick-start commands see the root [README.md](../README.md).

## The project

**Rare Disease, Real Kid: MVA Hackathon 2026** ([Hugging Face Space](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026), [dataset](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data)) centers on the genome and clinical story of a real child living with **Mosaic Variegated Aneuploidy (MVA)** — an ultra-rare genetic condition affecting fewer than 50 people worldwide, with no established treatment; care today means managing symptoms. The family has opened the case to the research community hoping someone can find an answer.

Participants receive whole-genome sequencing (WGS) data for the child plus a clinical phenotype document, and analyze them to surface candidate variants and genes that could explain the condition. **This is research exploration, not diagnosis or medical advice.**

## Project structure

```text
.
├── data/
│   ├── README.md                     # private-data contract and stage schemas
│   ├── raw/                          # source patient data (ignored)
│   ├── interim/                      # resumable private stage artifacts (ignored)
│   └── processed/                    # consolidated evidence tables (ignored)
├── database/refs/                    # public reference assets, local only
├── envs/                             # conda environments for VCF/VEP tooling
├── scripts/
│   ├── 01_run_qc.sh                  # bcftools QC + PASS-quality filtering
│   ├── 02_run_snpeff.sh              # snpEff functional annotation
│   ├── 03_run_ensembl_vep.sh         # offline Ensembl VEP annotation
│   ├── download_references.sh        # public reference downloads
│   ├── run_enrichment.py/.sh         # completed-stage evidence consolidation
│   └── run_exomiser.py/.sh           # Exomiser preparation, execution, import
└── src/mva_hackathon/                # reusable evidence-pipeline classes
```

All patient-derived and derived artifacts are ignored. The repository tracks
only code and documentation; no output table is inspected or committed.
## Methods brief

1. **QC and quality filtering** (`scripts/01_run_qc.sh`, env `mva_qc`): bcftools profiles the raw WGS VCF (`data/raw/WGS_EX2312012_HGWCNDSX7.vcf.gz`) and keeps PASS-quality variants, producing `data/interim/filtered/pass_quality.vcf.gz` plus the QC summaries in `data/interim/qc/`.
2. **snpEff annotation** (`scripts/02_run_snpeff.sh`, env `mva_snpeff`): annotates the filtered VCF against `GRCh38.mane.1.0.ensembl`; `scripts/py/filter_variants_snpeff.py` reduces it to `all_high_moderate_variants.csv`, `recessive_candidate_genes.csv` (homozygous / compound heterozygous patterns relevant to recessive disease), and `phenotype_focused_variants.csv` in `data/interim/snpeff/filtered_csv/`.
3. **Ensembl VEP annotation** (`scripts/03_run_ensembl_vep.sh`, env `mva_vep`): offline VEP run on GRCh38; `scripts/py/filter_variants_vep.py` produces high/moderate impact variants, recessive candidates, pathogenic variants, ClinVar-backed calls, phenotype-focused variants, and top clinical rankings in `data/interim/vep/filtered_csv/`.
4. **Evidence consolidation** (`scripts/run_enrichment.py`): normalizes one VEP table and merges only completed candidate-level CADD, gnomAD, REVEL, AlphaMissense, SpliceAI, AlphaGenome, and Exomiser artifacts. The tool writes the variant-level master table and triage shortlist in `data/processed/`. Use `--vep` to pass the existing private chunked-processing output by its current filename; no file must be renamed or inspected.


The two annotators are run in parallel so their candidate lists can be cross-checked; agreement between snpEff and VEP prioritizes the shortlist for manual review against the child's phenotype document.
