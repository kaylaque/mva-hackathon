# Script Directory & Track 2 CLI Workflows

This directory contains standalone CLI scripts and bash wrappers for running the MVA bioinformatics pipeline, upstream annotations, and Track 2 downstream phases (Phases 0, 1, and 2).

All scripts run without hardcoded patient files and support reproducible paths via CLI flags.

---

## Track 2 Downstream Workflow Scripts

Track 2 translates multimodal variant annotations into prioritized mechanism matrices and PPI/pathway networks without rerunning heavy upstream stages.

### 1. Phase 0: Candidate Intake, Gating & Consolidation
- **Script**: `scripts/consolidate_candidates.py` (Wrapper: `scripts/consolidate_candidates.sh`)
- **Purpose**: Consumes upstream `variant_evidence_master.csv` and optional Exomiser TSVs. Applies transparent quality and evidence gates (gnomAD AF < 0.001, CADD >= 20, SpliceAI >= 0.2, REVEL/AlphaMissense >= 0.5) and tags canonical MVA/SAC genes.
- **Direct Flat Output Usage**:
  ```bash
  uv run python scripts/consolidate_candidates.py \
    --evidence data/processed/variant_evidence_master.csv \
    --exomiser-variants data/interim/exomiser/mva_patient_exomiser.variants.tsv \
    --exomiser-genes data/interim/exomiser/mva_patient_exomiser.genes.tsv \
    --output-candidates data/processed/track2/00_phase0_consolidated_candidates.csv \
    --output-provenance data/processed/track2/00_phase0_provenance.json
  ```
  *(Or run `bash scripts/consolidate_candidates.sh` with no arguments to use the above default flat paths).*

### 2. Phase 1: Mechanism & Directionality Matrix
- **Script**: `scripts/build_mechanism_matrix.py` (Wrapper: `scripts/build_mechanism_matrix.sh`)
- **Purpose**: Classifies molecular consequences (Loss-of-Function, Hypomorphic Missense, Splice Disruptor), tags cellular pathway impacts, and evaluates compound heterozygous / multi-hit archetypes.
- **Direct Flat Output Usage**:
  ```bash
  uv run python scripts/build_mechanism_matrix.py \
    --candidates data/processed/track2/00_phase0_consolidated_candidates.csv \
    --output-variant-matrix data/processed/track2/01_phase1_variant_mechanism_matrix.csv \
    --output-gene-matrix data/processed/track2/02_phase1_gene_mechanism_matrix.csv \
    --output-provenance data/processed/track2/02_phase1_provenance.json
  ```
  *(Or run `bash scripts/build_mechanism_matrix.sh` with no arguments to use default flat paths).*

### 3. Phase 2: PPI and Pathway Network Context Engine
- **Script**: `scripts/run_network_analysis.py` (Wrapper: `scripts/run_network_analysis.sh`)
- **Purpose**: Queries live STRING v12.5, Reactome Analysis Service, and g:Profiler (g:Convert and g:GOSt) to map passing candidate seed genes into interactome topologies and enriched biological pathways.
- **Direct Flat Output Usage**:
  ```bash
  uv run python scripts/run_network_analysis.py \
    --gene-mechanism data/processed/track2/02_phase1_gene_mechanism_matrix.csv \
    --variant-matrix data/processed/track2/01_phase1_variant_mechanism_matrix.csv \
    --output-dir data/processed/track2/ \
    --flat-names \
    --string-score 700
  ```
  *(Or run `bash scripts/run_network_analysis.sh` with no arguments to use default flat paths).*

---

## Safe Layout Migration Tool

If phases were previously executed into subdirectories (`phase0/`, `phase1/`, `phase2/`), use `migrate_track2_layout.py` to copy them into the flat numbered Track 2 layout while keeping the subdirectories intact (duplicates preserved):

```bash
# 1. Preview copy actions without modifying any files:
uv run python scripts/migrate_track2_layout.py --dry-run

# 2. Perform copy and verify file byte sizes:
uv run python scripts/migrate_track2_layout.py

# 3. Check consistency of existing flat files against phase subdirectories:
uv run python scripts/migrate_track2_layout.py --validate-only
```

### Artifact Mapping Table
| Source File | Flat Destination File | Description |
|---|---|---|
| `phase0/consolidated_candidates.csv` | `00_phase0_consolidated_candidates.csv` | Filtered candidates table with triage gates |
| `phase0/provenance.json` | `00_phase0_provenance.json` | Intake provenance & filtering counts |
| `phase1/variant_mechanism_matrix.csv` | `01_phase1_variant_mechanism_matrix.csv` | Variant-level consequence & direction |
| `phase1/gene_mechanism_matrix.csv` | `02_phase1_gene_mechanism_matrix.csv` | Gene-level archetype & phase gates |
| `phase1/provenance.json` | `02_phase1_provenance.json` | Mechanism engine provenance |
| `phase2/seed_genes.csv` | `03_phase2_seed_genes.csv` | Seed genes evaluated in network |
| `phase2/string_edges.csv` | `04_phase2_string_edges.csv` | STRING v12.5 high-confidence PPI edges |
| `phase2/network_summary.csv` | `05_phase2_network_summary.csv` | Graph topology & global metrics |
| `phase2/pathway_enrichment.csv` | `06_phase2_pathway_enrichment.csv` | Reactome & g:Profiler enrichment terms |
| `phase2/network_gene_context.csv` | `07_phase2_network_gene_context.csv` | Unified gene-level network context table |
| `phase2/phase2_network_provenance.json` | `08_phase2_network_provenance.json` | API parameters, versions, and caveats |
