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

## Track 2 Workflows & Pipelines

Run full Phase 0 through Phase 3 pipelines using `uv run` or dedicated shell wrappers:

- **Phase 0 Candidate Consolidation**:
  ```bash
  uv run python scripts/consolidate_candidates.py --evidence <csv> --output-candidates data/processed/track2/00_phase0_consolidated_candidates.csv
  ```
- **Phase 1 Mechanism Matrix**:
  ```bash
  uv run python scripts/build_mechanism_matrix.py --candidates data/processed/track2/00_phase0_consolidated_candidates.csv --output-variant-matrix data/processed/track2/01_phase1_variant_mechanism_matrix.csv --output-gene-matrix data/processed/track2/02_phase1_gene_mechanism_matrix.csv
  ```
- **Phase 2 Network Analysis & Topology**:
  ```bash
  uv run python scripts/run_network_analysis.py --gene-mechanism data/processed/track2/02_phase1_gene_mechanism_matrix.csv --output-dir data/processed/track2/
  ```
- **Phase 2 Interactive Network Visualization**:
  ```bash
  uv run python scripts/run_network_visualization.py --gene-context data/processed/track2/07_phase2_network_gene_context.csv --edges data/processed/track2/04_phase2_string_edges.csv --enrichment data/processed/track2/06_phase2_pathway_enrichment.csv --output-dir data/processed/track2/
  ```
- **Phase 3 Query Gene Panel Generation**:
  ```bash
  uv run python scripts/run_query_gene_panel.py --track2-dir data/processed/track2 --output-dir data/processed/track2 --flat-names
  # Or wrapper: bash scripts/run_query_gene_panel.sh
  ```
- **Phase 3 Multi-Source Drug Repurposing Evidence**:
  ```bash
  uv run python scripts/run_drug_evidence.py --query-panel data/processed/track2/18_phase3_query_gene_panel.csv --gene-mechanism data/processed/track2/02_phase1_gene_mechanism_matrix.csv --network-context data/processed/track2/07_phase2_network_gene_context.csv --output-dir data/processed/track2/ --flat-names
  # Or wrapper: bash scripts/run_drug_evidence.sh
  ```
- **Phase 4 Signed Pathway/Reaction Adjudication & Literature Validation**:
  ```bash
  uv run python scripts/run_phase4_adjudication.py --output-dir data/processed/track2/ --flat-names
  # Or wrapper: bash scripts/run_phase4_adjudication.sh
  ```

## Compass: saved evidence to research decisions

Compass (`src/mva_hackathon/compass.py`, `scripts/run_compass.py`) is an evidence-to-decision
layer over already-saved results. It does **not** run enrichment, submit DAVID jobs, query
STRING, or rerun any Phase 1–3 pipeline. It reads verified post-analysis tables and safe
projected network columns, retrieves public membership/protein/literature/pharmacology
records through one content-addressed cache, and writes tables that state explicitly what is
supported, what is only context, and what is unresolved.

```bash
# All destinations must be outside every Git repository.
uv run python scripts/analyze_pathway_results.py \
  --reuse-pathway-dir data/processed/track2/pathway \
  --david-dir /external/immutable/david-import \
  --output-dir /external/postanalysis-$(date -u +%Y%m%d-%H%M%S)

uv run python scripts/run_compass.py \
  --postanalysis-dir /external/postanalysis-<stamp> \
  --network-dir data/processed/track2 \
  --output-dir /external/compass-<stamp> \
  --cache-dir /external/compass-cache

# Supply source-grounded review, then finalize decisions without repeating public requests.
uv run python scripts/run_compass.py \
  --postanalysis-dir /external/postanalysis-<stamp> \
  --network-dir data/processed/track2 \
  --output-dir /external/compass-<stamp> \
  --cache-dir /external/compass-cache \
  --reviewed-evidence /external/reviewed_evidence.csv \
  --resume

# Replay from the cache alone: no HTTP is issued and the stable tables are identical.
uv run python scripts/run_compass.py \
  --postanalysis-dir /external/postanalysis-<stamp> \
  --output-dir /external/compass-verify-<stamp> \
  --cache-dir /external/compass-cache \
  --reviewed-evidence /external/reviewed_evidence.csv \
  --offline
```

Exit codes: `0` completed processing (scientific nulls and unresolved biology included),
`2` a required public source was incomplete, `1` unsafe or invalid inputs. A new run refuses
an existing output directory; `--resume` requires matching input projection hashes and
selection rules, and a changed reviewed-evidence file invalidates only the decisions that
depend on it — verified public cache records are reused.

### Where each step sits, and what differs

| | g:Profiler pathway | DAVID pathway | Network pathway | **Compass** |
| --- | --- | --- | --- | --- |
| Where it sits | step 1: ask a service about a gene list | step 2: ask DAVID the same question by hand | step 3: build PPI context for the same genes | **after all three: reads their saved files** |
| What it reads | the gene list | the gene list | Phase 1 candidates | the saved g:Profiler, DAVID and network outputs |
| What it does | enrichment over many sources | enrichment + clustering + chart | PPI topology, components, neighbour genes | judges each provider on its own numbers, then looks up public protein/literature/drug facts |
| What it writes | terms, p-values, observed intersection | chart and cluster exports with native p-values | edges, components, node context | eligibility, cited rationales, decisions, hold/reject ledger, test plans |
| What it adds | — | a second opinion, not a second test | which genes cluster together | **why a target is or is not worth testing, and what would falsify it** |
| What it cannot claim | that a term is causal | that its numbers match g:Profiler's | that a neighbour is causal | clinical efficacy, patient direction, or experimental proof |

Compass adds no enrichment of its own. It only reads what the three steps above already
wrote, so a null enrichment result stays null: on the last real run all 16 662 significant
candidates failed the provider threshold, 0 pathways qualified, and the work continued in
the mechanism-prior lane instead.

**Gene annotation is not therapeutic evidence.** Enrichment membership, GO annotation,
protein function text and modification annotation describe what a gene or protein is
annotated to do. They do not measure activity, establish a disease direction, or predict
that modulating a target helps. Compass therefore keeps annotation in its own lane and
requires a cited, reviewed assertion — with model context — before any target reaches
`TEST_FIRST`, and requires recorded safety readouts before `ADVANCE_RESEARCH`.

### Gap ledger

Each row separates what is implemented, what actually happened on the last real run, and
what remains biologically uncertain.

| # | Gap | Implementation status | Run status (last real execution) | Remaining biological uncertainty |
| --- | --- | --- | --- | --- |
| 1 | DAVID stable-term recovery | Implements exact source-ID normalisation; ambiguous titles are retained unmapped with `STOP_UNMAPPED` | 111 unmapped annotations retained and reviewed, none silently matched | Unmapped DAVID titles cannot be compared with g:Profiler without manual curation |
| 2 | Shared pathway eligibility | Per-provider gate on the provider's own corrected statistic, stable identity and observed contributors | 16 773 rows evaluated: 16 662 `STOP_NOT_SIGNIFICANT`, 111 `STOP_UNMAPPED`, 0 eligible | Current enrichment is null for every pathway, so no pathway-supported target exists yet |
| 3 | Member retrieval with completeness | Reactome (nested entities), GO BP (paged, `NOT` excluded), KEGG (human links), WikiPathways (typed nodes) behind one cache, with `complete` flags | No eligible identifier required a fetch, so `expanded_members.csv` is header-only | Any future eligible term must re-verify these endpoints; WikiPathways moved to the assets path |
| 4 | Unified target rationale | `build_rationales` joins eligibility, memberships, protein roles, reviewed claims and network membership per target | 34 panel rationales produced, each with cited reason codes | Panel-wide direction still needs reviewed evidence per target |
| 5 | Disease-model direction, never patient direction | Direction comes only from cited reviewed claims; conflicting claims yield `UNKNOWN`, never a patient-specific inference | 4 targets with a cited reduced-function state, 1 altered-function association, 29 left `UNKNOWN` | Human genetics establishes causality for a syndrome, not a modulation direction |
| 6 | Target-action matching | `assess_target_action` maps exact action tokens; partial agonists, binders and mixed labels stay `UNKNOWN` | Pharmacology queried only for the 2 `TEST_FIRST` targets; both returned no matching human target | Absence of a ChEMBL target is not absence of biology |
| 7 | Context-specific worsening and cancer risk | Cited `HAZARD`/`WORSENS` claims force `REJECT`; unknown safety is recorded as unknown, never safe | TTK rejected on a cited induced-instability finding in tumour cells | The hazard is tumour-cell context; a constitutional risk estimate is not established |
| 8 | Mechanistic literature claims with citations | Bounded PubMed search per candidate plus schema-validated reviewed evidence with source section and context | 291 unreviewed search hits retained for audit; 9 reviewed assertions across 7 genes | Review is abstract-only and partial: the remaining panel genes await review |
| 9 | Evidence classes and conflict review | Evidence classes are kept separate; conflicting reviewed direction yields `HOLD` + `EVIDENCE_CONFLICT` | 1 conflict case exercised in unit tests; no conflicting reviewed pair in this run | Class-specific weightings remain a research decision, not a pipeline output |
| 10 | Tissue, developmental, dose and duration context | Context fields default to unreported and are never invented | Reviewed rows record model/context; dose and duration are empty for every abstract-only claim | Full text is needed for dose, duration and tissue-specific effects |
| 11 | Separate mechanism-prior lane | Finite MVA/SAC/cohesin panel researched independently of enrichment, labelled `MECHANISM_PRIOR` | 34 panel targets researched with no pathway claim (29 `MECHANISM_PRIOR`, 5 also `NETWORK_CONTEXT`) | Panel membership is a curated prior, not a discovery result |
| 12 | Cross-source consolidation without double counting | One request cache keyed by source/endpoint/params; repeated statistics preserved; one source is never replication | 34 gene-identity, 34 protein-role and 34 literature requests plus 2 pharmacology requests, each recorded once in `source_status.csv` | A shared PMID/annotation reused by two databases is still one source |
| 13 | Hold/rejection ledger | Every held, unsupported and rejected entity is written with stage, decision and reason codes | 18 789 ledger rows: 16 662 stopped annotations, 111 unmapped, 1 984 context-only network genes, 31 held targets, 1 rejected target | Nothing withheld is invisible, but the ledger is a status record, not evidence |
| 14 | Proposed experimental validation | `build_test_plans` writes readouts, controls and rescue/worsening criteria; unresolved targets get an evidence-resolution plan | 2 `PROPOSED` plans; 32 `EVIDENCE_RESOLUTION_REQUIRED` plans | Concentrations, thresholds and model choice require experimental calibration |
| 15 | Reproducible evidence-to-decision report | Report, tables and hashed provenance are generated together and published atomically | Live and offline-replay runs produced byte-identical stable tables | Reproducibility covers software and public snapshots, not the biological conclusion |

Superseded artifacts: the earlier derived report directory
`~/Documents/mva-pathway-postanalysis-20261013-143202` and the manually edited report are
superseded by this workflow and were left in place rather than deleted. Historical
`HIT_GENES` fields produced before the g:Profiler membership fix are untrusted and are not
consumed anywhere.

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

# ------------------------------------------------------------------------------
# Phase 0 & Phase 1 Downstream Workflows (Lightweight, Non-Rerunning)
# ------------------------------------------------------------------------------
# Phase 0: Candidate Intake, Quality & Evidence Gating, and Candidate Consolidation
# Consumes caller-supplied evidence and optional Exomiser tables via CLI.
# Does NOT rerun QC, VEP, snpEff, enrichment, or Exomiser.
# Can write to flat Track 2 layout directly or hierarchical paths:
uv run python scripts/consolidate_candidates.py \
  --evidence data/processed/variant_evidence_master.csv \
  --exomiser-variants data/interim/exomiser/mva_patient_exomiser.variants.tsv \
  --exomiser-genes data/interim/exomiser/mva_patient_exomiser.genes.tsv \
  --output-candidates data/processed/track2/00_phase0_consolidated_candidates.csv \
  --output-provenance data/processed/track2/00_phase0_provenance.json

# Phase 1: Biological Mechanism & Directionality Matrix.
# Consumes Phase 0 consolidated candidates table.
uv run python scripts/build_mechanism_matrix.py \
  --candidates data/processed/track2/00_phase0_consolidated_candidates.csv \
  --output-variant-matrix data/processed/track2/01_phase1_variant_mechanism_matrix.csv \
  --output-gene-matrix data/processed/track2/02_phase1_gene_mechanism_matrix.csv \
  --output-provenance data/processed/track2/02_phase1_provenance.json

# ------------------------------------------------------------------------------
# Phase 2 Downstream Workflow: PPI & Pathway Network Context Engine
# ------------------------------------------------------------------------------
# Phase 2: Protein-Protein Interaction (PPI) and Pathway Network Context.
# Consumes caller-supplied Phase 1 gene mechanism matrix (and optional variant matrix).
# Queries external biomedical APIs at runtime (STRING v12.5, Reactome, g:Profiler).
# Use --flat-names to write directly to flat numbered Track 2 filenames (03_.. through 08_..):
uv run python scripts/run_network_analysis.py \
  --gene-mechanism data/processed/track2/02_phase1_gene_mechanism_matrix.csv \
  --variant-matrix data/processed/track2/01_phase1_variant_mechanism_matrix.csv \
  --output-dir data/processed/track2/ \
  --flat-names \
  --string-score 700

# ------------------------------------------------------------------------------
# Phase 2 Downstream Workflow: Interactive Network & Pathway Visualizations
# ------------------------------------------------------------------------------
# Generates deterministic spring-layout Plotly HTML network graphs (largest component,
# MVA-focused component, pathway enrichment bubble chart) and node/component CSV summaries:
uv run python scripts/run_network_visualization.py \
  --gene-context data/processed/track2/07_phase2_network_gene_context.csv \
  --edges data/processed/track2/04_phase2_string_edges.csv \
  --enrichment data/processed/track2/06_phase2_pathway_enrichment.csv \
  --output-dir data/processed/track2/ \
  --flat-names

# ------------------------------------------------------------------------------
# Phase 3 Downstream Workflow: Multi-Source Drug Repurposing Evidence Engine
# ------------------------------------------------------------------------------
# Queries runtime pharmacological sources (Open Targets GraphQL, ChEMBL REST, DGIdb GraphQL)
# to build a normalized drug-evidence matrix, evaluate variant mechanism -> pathway ->
# target action -> safety directionality, flag mechanistic hazards, and produce safety rankings:
uv run python scripts/run_drug_evidence.py \
  --gene-mechanism data/processed/track2/02_phase1_gene_mechanism_matrix.csv \
  --network-context data/processed/track2/07_phase2_network_gene_context.csv \
  --output-dir data/processed/track2/ \
  --flat-names

# ------------------------------------------------------------------------------
# Track 2 Artifact Migration & Validation
# ------------------------------------------------------------------------------
# Safely copy hierarchical phase artifacts (phase0/, phase1/, phase2/) into flat
# ordered Track 2 files directly in data/processed/track2/, keeping phase folders:
uv run python scripts/migrate_track2_layout.py --dry-run        # preview copies
uv run python scripts/migrate_track2_layout.py                  # copy & validate
uv run python scripts/migrate_track2_layout.py --validate-only  # verify SHA-256 matches
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

## Independent pathway workflow

The Phase 1 gene matrix also forks into an independent P0–P7 pathway workflow; it does not modify the existing STRING Phase 2 or drug-evidence Phase 3 branches. Run its privacy-safe synthetic fixture first:

```bash
uv run python scripts/run_pathway_workflow.py --gene-mechanism tests/fixtures/pathway_sample/genes.csv --offline-fixture tests/fixtures/pathway_sample --output-dir /tmp/pathway-sample
```

Its default output and cache locations are private, ignored paths. See `docs/JOURNEY.md`, `scripts/README.md`, and `data/README.md` for limitations and output handling.
[data/README.md#A17B]

## Pathway workflow private artifacts

`processed/track2/pathway/` holds P0–P7 tables, provenance, and cumulative reports; `interim/pathway_cache/` holds response fixtures and request caches. Both are ignored and therefore are **not agent-accessible** under this repository's policy. Only normalized gene identifiers may leave the repository; variant identifiers, genotypes, coordinates, clinical text, patient identifiers, and internal scores must never be sent to external services. Reports are private research artifacts, not clinical evidence.
[scripts/README.md#C610]

### Independent pathway workflow

`run_pathway_workflow.py` accepts `--gene-mechanism`, `--offline-fixture`, `--output-dir`, optional `--network-panel`, and `--min-passing`, plus the DAVID addendum flags `--david`, `--reuse-pathway-dir`, `--david-exports`, `--with-david`, and `--david-output-dir`. Run the synthetic fixture before any authorized real input. Offline mode never falls back to live services. The command returns `0` for a complete comparison and `2` when the STRING baseline is unavailable; a partial run is evidence of unavailable sources, not zero evidence. The DAVID leg is registration-free: it prepares gene-only lists for the native web interface and imports the exports downloaded from it, so no API key, account, or `DAVID_EMAIL` session exists and no DAVID API is called.
[docs/JOURNEY.md#B8CF]

---

## Independent pathway workflow (P0–P7)

Phase 1 forks into an independent pathway branch while existing STRING Phase 2 and drug-evidence Phase 3 remain unchanged:

```text
Phase 1 gene mechanism matrix
├── STRING Phase 2 → existing drug-evidence Phase 3
└── Pathway P0 input QC → P1 enrichment → P2 prioritization → P3 mechanisms
    → P4 full-pathway targets → P5 drug evidence → P6 comparison → P7 report
```

The pathway background is Phase 1 evaluated genes, not the whole-exome callable universe. Provider BH-adjusted values remain distinct from raw hypergeometric values; enrichment is not efficacy, and pharmacologic direction is not biological rescue. HP terms are observational gene-annotation enrichment, never observed patient phenotype. The current branch contains synthetic-only verification; private reports and real-input execution remain unavailable under repository access policy.

### Independent pathway branch

The independent P0–P7 pathway workflow consumes only the Phase 1 gene matrix,
leaving the existing STRING Phase 2 / drug-evidence Phase 3 branch unchanged.

```bash
# Synthetic privacy-safe smoke run; no network access.
uv run python scripts/run_pathway_workflow.py \
  --gene-mechanism tests/fixtures/pathway_sample/genes.csv \
  --offline-fixture tests/fixtures/pathway_sample \
  --output-dir /tmp/mva-pathway-sample

# Authorized gene-only live run; private outputs remain ignored.
bash scripts/run_pathway_workflow.sh
```

The cumulative local report is
`data/processed/track2/pathway/p7_pathway_report.md`; stage tables and
provenance are colocated. Read [`docs/JOURNEY.md`](docs/JOURNEY.md) for the
method and source comparison, and [`scripts/README.md`](scripts/README.md) for
options and partial-run semantics.
