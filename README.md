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

## End-to-end analysis flow

All analysis uses GRCh38. Patient data and every derived artifact remain local
under `data/` and are ignored by Git.

```text
data/input/patient/
  WGS VCF + tabix index, clinical phenotype document
        │
        ▼
[1] QC and quality filtering — bcftools / tabix
  scripts/01_run_qc.sh
  └── data/interim/qc/raw_vcf_stats.txt
  └── data/interim/filtered/pass_quality.vcf.gz + .tbi
        │
        ├───────────────────────────────────────┐
        ▼                                       ▼
[2a] Secondary annotation — SnpEff          [2b] Primary annotation — Ensembl VEP
  scripts/02_run_snpeff.sh                     scripts/03_run_ensembl_vep.sh
  scripts/py/filter_variants_snpeff.py         scripts/py/filter_variants_vep.py
  └── data/interim/snpeff/                     └── data/interim/vep/
                                                   high/moderate VEP candidate CSV
                                                         │
                                                         ▼
[3] Variant-level evidence enrichment
  Input: VEP high/moderate candidate table
  ├── baseline QC and descriptive summaries
  ├── population rarity: AF_POPULATION + gnomAD
  ├── clinical evidence: CLIN_SIG + ACMG fields
  ├── protein/coding: SIFT, PolyPhen, CADD, REVEL, AlphaMissense
  ├── splicing: SpliceAI
  ├── regulatory/molecular: AlphaGenome
  ├── inheritance: genotype and multi-hit gene patterns
  └── phenotype: panel context and optional Exomiser/HPO ranking
        │
        ▼
data/interim/{cadd,gnomad,revel,alphamissense,spliceai,alphagenome,exomiser}/
        │
        ▼
[4] Evidence consolidation — scripts/run_enrichment.py
  ├── data/processed/variant_evidence_master.csv
  └── data/processed/prioritized_candidates.csv
        │
        ▼
[5] Research interpretation
  rare/protein-damaging, splice-driven, regulatory, phenotype-focused,
  and recessive/multi-hit candidates → manual MVA mechanism review
```

SnpEff and VEP run from the same filtered VCF and may be compared, but the
evidence master table is VEP-centred. Each enrichment stage is optional: the
consolidator retains VEP evidence when a completed external-stage CSV is absent.
The integration score is research prioritization, not an ACMG classification or
diagnosis; in-silico predictors are correlated and are not independent clinical
evidence.

### Running each stage

```bash
# Download public reference assets only; no patient data.
bash scripts/download_references.sh

# QC → SnpEff → VEP. These are separate long-running workflows.
conda env create -f envs/mva_qc.yml
conda env create -f envs/mva_snpeff.yml
conda env create -f envs/mva_vep.yml
bash scripts/01_run_qc.sh
bash scripts/02_run_snpeff.sh
bash scripts/03_run_ensembl_vep.sh

# Consolidate existing enrichment-stage outputs. This does not launch CADD,
# gnomAD, REVEL, AlphaMissense, SpliceAI, AlphaGenome, or Exomiser.
uv run python scripts/run_enrichment.py \
  --vep /private/path/to/current-high-moderate-or-chunked-output.csv

# Prepare Exomiser candidate VCF/YAML. Add --execute only to launch Exomiser.
uv run python scripts/run_exomiser.py \
  --jar /path/to/exomiser-cli.jar \
  --hpo HP:0000001
```

Read [`data/README.md`](data/README.md) before placing or naming private
artifacts.

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
| REVEL v1.3 | Missense-specific ensemble score | [Official project](https://sites.google.com/site/revelgenomics/), [Zenodo archive](https://zenodo.org/records/7072866), [research paper](https://doi.org/10.1016/j.ajhg.2016.08.016) |
| AlphaMissense | Missense pathogenicity score and class | [Official project](https://deepmind.google/technologies/alphamissense/), [official repository](https://github.com/google-deepmind/alphamissense), [GRCh38 archive](https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz) |
| SpliceAI | Predicted splice donor/acceptor disruption | [GitHub](https://github.com/Illumina/SpliceAI), [GRCh38 reference](https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz) |
| AlphaGenome | Regulatory and molecular-effect evidence | [Official project](https://deepmind.google/science/alphagenome/), [official documentation](https://www.alphagenomedocs.com/), [official repository](https://github.com/google-deepmind/alphagenome), [research paper](https://doi.org/10.1038/s41586-025-10014-0) |
| Exomiser | HPO and disease-model phenotype ranking | [Exomiser](https://github.com/exomiser/Exomiser), [data downloads](https://data.monarchinitiative.org/exomiser/) |
| bcftools / samtools | VCF QC, filtering, indexing, and reference indexing | [bcftools](https://samtools.github.io/bcftools/), [samtools](https://www.htslib.org/) |
| SnpEff / SnpSift | Secondary functional annotation and field extraction | [SnpEff](https://pcingola.github.io/SnpEff/) |
| Python / uv stack | Environment and VCF/table processing | [uv](https://docs.astral.sh/uv/), [cyvcf2](https://brentp.github.io/cyvcf2/), [pysam](https://pysam.readthedocs.io/) |

Do not treat any computational predictor, population-frequency annotation, or
rank as diagnostic evidence. Manual review, phenotype fit, inheritance, and
the hackathon's research-only scope remain required.

## `variant_evidence_master.csv` data dictionary

The master table contains annotation and research-triage evidence. A row may
represent one transcript annotation rather than one unique genomic variant.
No single column, score, or combination of computational predictors establishes
causality or diagnosis.

### Priority and gene context

| Column | Meaning | Interpretation |
| --- | --- | --- |
| `PRIORITY_SCORE` | Internal workflow ranking score. | Higher usually means review sooner; it is not a calibrated pathogenicity probability. |
| `SYMBOL` | Standard gene symbol, such as `BUB1B`. | The gene assigned to this annotation row; one variant can map to multiple genes. |
| `GENE_ID` | Stable gene identifier, usually from Ensembl. | Identifies the gene independently of symbol changes. |
| `INHERITANCE_MODEL` | Imported disease-inheritance label. | Describes reported gene-disease inheritance, not the patient's confirmed inheritance. |
| `PHENOTYPE_PANEL` | Imported panel or phenotype category. | Explains why a gene was included upstream; it is not an observed patient symptom. |

### Genomic identity

| Column | Meaning | Interpretation |
| --- | --- | --- |
| `CHROM` | Source chromosome name. | May be written as `1`, `chr1`, `X`, and so on. |
| `CHROM_NORM` | Normalized chromosome name. | Standardizes forms such as `chr1` → `1` and `23` → `X`. |
| `POS` | 1-based genomic position. | Location of the first reference base on the claimed reference assembly. |
| `REF` | Reference allele. | Sequence expected in the reference genome. |
| `ALT` | Alternate allele. | Sequence observed instead of `REF`. |
| `VARIANT_KEY` | Normalized key such as `1:12345:A:G`. | Primary key used to join records for the same genomic variant. |
| `VARIANT_ID` | Identifier supplied by the source workflow. | May be another coordinate-based or source-specific identifier. |
| `EXISTING_VARIATION` | Known database identifier. | Often a dbSNP `rs` ID; missing does not prove that the variant is novel. |

### Transcript and HGVS annotation

| Column | Meaning | Interpretation |
| --- | --- | --- |
| `FEATURE_TRANSCRIPT` | Transcript used for annotation. | Different transcripts can produce different consequences. |
| `CANONICAL` | Canonical-transcript flag. | Useful for selecting a representative row, but not always the clinically relevant transcript. |
| `HGVS_C` | HGVS coding-DNA description. | Transcript-relative change such as `c.123A>G`. |
| `HGVS_P` | HGVS protein description. | Predicted protein change such as `p.Gly41Arg`. |
| `CDNA_POS` | Position in the complete transcript cDNA. | Locates the variant relative to the transcript. |
| `CDS_POS` | Position in the coding sequence. | Locates the variant within the protein-coding region. |
| `PROTEIN_POS` | Affected amino-acid position. | May be an integer, range, uncertain value, or missing. |
| `AMINO_ACIDS` | Reference and alternate amino acids. | For example, `G/R` means glycine changes to arginine. |
| `CODONS` | Reference and altered codons. | Shows the codon-level nucleotide change. |

### Consequence and clinical annotation

| Column | Meaning | Interpretation |
| --- | --- | --- |
| `CONSEQUENCE` | VEP consequence term or `&`-separated terms. | Examples include missense, frameshift, stop-gain, and splice effects. |
| `IMPACT` | Broad VEP impact label. | `HIGH`, `MODERATE`, `LOW`, or `MODIFIER`; not proof of pathogenicity. |
| `CLIN_SIG` | Imported clinical-significance text, commonly from ClinVar. | May report P/LP, VUS, benign, or conflicting interpretations without full condition/review context. |
| `ACMG_CLASSIFICATION` | Imported ACMG-style classification. | Not an independently verified classification for this case. |
| `ACMG_CRITERIA` | Upstream ACMG evidence codes. | Codes such as `PM2`, `PP3`, or `PVS1`; their presence does not prove correct application. |

### SIFT, PolyPhen, and source population frequency

| Column | Meaning | Interpretation |
| --- | --- | --- |
| `SIFT_PRED` | SIFT category. | Usually `deleterious` or `tolerated`; mainly relevant to missense variants. |
| `SIFT_SCORE` | SIFT score, normally 0–1. | Lower means a stronger predicted effect; missing is unknown, not benign. |
| `POLYPHEN_PRED` | PolyPhen category. | Commonly benign, possibly damaging, or probably damaging. |
| `POLYPHEN_SCORE` | PolyPhen score, normally 0–1. | Higher means a stronger predicted damaging effect; not disease probability. |
| `AF_POPULATION` | Population allele frequency from the source annotation. | Provenance must be checked; `0.001` means about 0.1% of sampled alleles. |

### CADD

| Column | Meaning | Interpretation |
| --- | --- | --- |
| `CADD_RAW` | Raw CADD model score. | Higher generally means more predicted deleteriousness; the raw scale has no simple clinical meaning. |
| `CADD_PHRED` | Rank-scaled CADD score. | Roughly 10 = top 10%, 20 = top 1%, and 30 = top 0.1% of possible substitutions. |
| `CADD_SOURCE` | CADD annotation provenance. | Identifies the API, local file, version, or variant class used. |

### gnomAD

| Column | Meaning | Interpretation |
| --- | --- | --- |
| `GNOMAD_STATUS` | Lookup result. | `not_found` means reported not found with coverage/provenance unverified, not confirmed absence. |
| `GNOMAD_AF` | Overall alternate-allele frequency. | Frequency across the available gnomAD reference cohort. |
| `GNOMAD_AC` | Alternate-allele count. | Number of observed alternate alleles. |
| `GNOMAD_AN` | Allele number. | Number of evaluated alleles and the denominator for `AC`. |
| `GNOMAD_NHOMALT` | Alternate-homozygote count. | Reference-cohort count; it does not describe this patient's genotype. |
| `GNOMAD_GRPMAX_AF` | Highest raw frequency among population groups. | Detects ancestry-group enrichment; it is not filtering AF/FAF. |
| `GNOMAD_RARITY` | Workflow rarity category. | Threshold-dependent label such as rare, common, or unknown. |
| `GNOMAD_RARE_FLAG` | Workflow Boolean rarity flag. | Check the source and threshold before relying on it. |
| `RARE_AND_CADD_HIGH` | Combined rarity/CADD flag. | Research-triage shortcut, not clinical evidence by itself. |

Useful consistency checks are `AC ≤ AN`, `2 × NHOMALT ≤ AC`, and—when
`AN > 0`—`AF ≈ AC / AN`.

### REVEL and AlphaMissense

| Column | Meaning | Interpretation |
| --- | --- | --- |
| `REVEL_SCORE` | REVEL missense score, normally 0–1. | Higher indicates stronger predicted pathogenicity for missense substitutions only. |
| `AM_SCORE` | AlphaMissense score, normally 0–1. | Higher indicates a stronger predicted protein effect; it is not patient disease probability. |
| `AM_CLASS` | AlphaMissense category. | Likely benign, ambiguous, or likely pathogenic; not an ACMG classification. |
| `AM_UNIPROT` | UniProt protein accession. | Identifies the protein sequence used for the prediction. |
| `AM_TRANSCRIPT` | Transcript used by AlphaMissense. | Compare with `FEATURE_TRANSCRIPT`; disagreement may explain differing predictions. |
| `AM_PROTEIN_VARIANT` | Protein substitution used for lookup. | Identifies the exact amino-acid change, such as `G41R`. |

### SpliceAI

| Column | Meaning | Interpretation |
| --- | --- | --- |
| `SPLICEAI_GENE` | Gene associated with the splice prediction. | May differ from `SYMBOL`; preserve both rather than overwriting either. |
| `DS_AG` | Acceptor-gain delta score. | Predicted creation or strengthening of an acceptor site. |
| `DS_AL` | Acceptor-loss delta score. | Predicted weakening or loss of an acceptor site. |
| `DS_DG` | Donor-gain delta score. | Predicted creation or strengthening of a donor site. |
| `DS_DL` | Donor-loss delta score. | Predicted weakening or loss of a donor site. |
| `DP_AG` | Acceptor-gain position offset. | Add to `POS` to estimate the predicted genomic event position. |
| `DP_AL` | Acceptor-loss position offset. | Add to `POS` to estimate the predicted genomic event position. |
| `DP_DG` | Donor-gain position offset. | Add to `POS` to estimate the predicted genomic event position. |
| `DP_DL` | Donor-loss position offset. | Add to `POS` to estimate the predicted genomic event position. |
| `SPLICEAI_MAX_DS` | Largest imported `DS_*` score. | Summarizes the strongest predicted event; inspect the event-specific fields for meaning. |

`DP_AG`, `DP_AL`, `DP_DG`, and `DP_DL` are signed position offsets, **not**
sequencing depth. SpliceAI values around 0.2, 0.5, and 0.8 are commonly used
as high-recall, general, and high-precision research reference points; they are
not proof of splice disruption.

### Reading a row

Review identity first (`CHROM`, `POS`, `REF`, `ALT`, `VARIANT_KEY`), then the
gene/transcript and consequence, population frequency, clinical annotations,
protein predictors, and splice evidence. Genotype, allele balance, phase,
segregation, confirmed phenotype, cytogenetics, CNV/SV evidence, and mosaic
fraction are not supplied by these columns.

## Project guidance

See [`AGENTS.md`](AGENTS.md) for repository architecture, data-handling rules, VCF analysis conventions, and QA guidance.
