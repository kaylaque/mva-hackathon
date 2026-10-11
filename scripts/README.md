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

### 4. Phase 2 Visualization: Interactive Interactome & Pathway Graphs
- **Script**: `scripts/run_network_visualization.py` (Wrapper: `scripts/run_network_visualization.sh`)
- **Purpose**: Generates reproducible interactive Plotly visualizations and structural component tables without hairball dilution, using a deterministic Fruchterman-Reingold layout with fixed seeds.
- **Direct Flat Output Usage**:
  ```bash
  uv run python scripts/run_network_visualization.py \
    --gene-context data/processed/track2/07_phase2_network_gene_context.csv \
    --edges data/processed/track2/04_phase2_string_edges.csv \
    --enrichment data/processed/track2/06_phase2_pathway_enrichment.csv \
    --output-dir data/processed/track2/ \
    --flat-names
  ```
  *(Or run `bash scripts/run_network_visualization.sh` with no arguments to use default flat paths).*

### 5. Phase 3 Query Gene Panel Generation
- **Script**: `scripts/run_query_gene_panel.py` (Wrapper: `scripts/run_query_gene_panel.sh`)
- **Purpose**: Derives an objective, graph-grounded, multi-tiered gene query panel from Phase 1 and Phase 2 artifacts:
  - Gate 0: Symbol hygiene & non-query quarantine (quarantining unmapped symbols, HLA-* complex genes, readthrough transcripts, and pseudogenes).
  - Gate 1: Core MVA causal genes (7 genes: BUB1B, CEP57, TRIP13, CENATAC, MAD1L1, SLF2, SMC5) injected even if absent from Phase 1.
  - Gate 2: Graph-derived interactome neighbors from the programmatically detected MVA-focus component (differentiating one-hop vs component-member tiers).
  - Gate 3: Biological context modifiers (ranked Phase 1/2 candidates + curated literature references TP53, ATM, NF1, MAD2L2).
  - Gate 4: Positive druggable controls and programmatically derived negative sensory receptor controls.
- **Direct Flat Output Usage**:
  ```bash
  uv run python scripts/run_query_gene_panel.py \
    --track2-dir data/processed/track2 \
    --output-dir data/processed/track2 \
    --flat-names
  ```
  *(Or run `bash scripts/run_query_gene_panel.sh` with no arguments to use default flat paths).*

### 6. Phase 3: Multi-Source Drug Repurposing Evidence & Safety Review
- **Script**: `scripts/run_drug_evidence.py` (Wrapper: `scripts/run_drug_evidence.sh`)
- **Purpose**: Queries Open Targets v4, ChEMBL v2.0, and DGIdb GraphQL to map candidate and interactome seed genes into pharmacological space. Evaluates variant mechanism -> pathway -> target action -> safety directionality, detecting directional contradictions (e.g. BUB1B/SAC inhibitors in a child with baseline checkpoint impairment) and strictly prohibiting clinical treatment recommendations.
- **Direct Flat Output Usage**:
  ```bash
  uv run python scripts/run_drug_evidence.py \
    --query-panel data/processed/track2/18_phase3_query_gene_panel.csv \
    --gene-mechanism data/processed/track2/02_phase1_gene_mechanism_matrix.csv \
    --network-context data/processed/track2/07_phase2_network_gene_context.csv \
    --output-dir data/processed/track2/ \
    --flat-names
  ```
  *(Or run `bash scripts/run_drug_evidence.sh` with no arguments to use default flat paths).*

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
| `phase2/network_nodes.csv` | `09_phase2_network_nodes.csv` | Node metrics, component ranks, and degrees |
| `phase2/component_summary.csv` | `10_phase2_component_summary.csv` | Connected component ranking and MVA hits |
| `phase2/largest_component.html` | `11_phase2_largest_component.html` | Interactive graph of largest component (Comp 1) |
| `phase2/mva_focus_network.html` | `12_phase2_mva_focus_network.html` | Interactive graph of MVA component (Comp 2) |
| `phase2/pathway_enrichment.html` | `13_phase2_pathway_enrichment.html` | Interactive pathway enrichment bubble chart |
| `phase3/drug_evidence_matrix.csv` | `14_phase3_drug_evidence_matrix.csv` | Normalized multi-source drug-target evidence matrix |
| `phase3/safety_ranking.csv` | `15_phase3_safety_ranking.csv` | Drug safety, directionality tier, and research ranking |
| `phase3/mechanism_action_map.csv` | `16_phase3_mechanism_action_map.csv` | Gene-level mechanism-to-action summary map |
| `phase3/phase3_repurposing_provenance.json` | `17_phase3_repurposing_provenance.json` | Provenance metadata, API endpoints, and disclaimers |
| `phase3/query_gene_panel.csv` | `18_phase3_query_gene_panel.csv` | Graph-derived multi-tiered query gene panel across 4 evidence gates |
| `phase3/query_panel_provenance.json` | `18_phase3_query_panel_provenance.json` | Query gene panel provenance, gate counts, and component metrics |
| `phase3/query_panel_citations.csv` | `19_phase3_query_panel_citations.csv` | Normalized long-format citations table with verified PMIDs/DOIs |

## Independent pathway workflow (P0–P7)

`run_pathway_workflow.py` is independent from the existing STRING/drug branch.
It accepts only the Phase 1 gene-level allowlist; P0 rejects missing required
counts and never loads `TOP_VARIANT_KEYS` or other variant fields.

```bash
# Required smoke test: fixture-only, no live fallback.
uv run python scripts/run_pathway_workflow.py \
  --gene-mechanism tests/fixtures/pathway_sample/genes.csv \
  --offline-fixture tests/fixtures/pathway_sample \
  --output-dir /tmp/mva-pathway-sample

# Live gene-only workflow defaults.
bash scripts/run_pathway_workflow.sh
```

Implemented options: `--gene-mechanism`, `--network-context`,
`--network-panel`, `--network-drugs`, `--min-passing`, `--q-threshold`,
`--top-pathways`, `--stage`, `--offline-fixture`, `--no-chembl`,
`--no-open-targets`, `--no-dgidb`, `--david`, `--output-dir`, and `--force`.
`--david` requires `DAVID_EMAIL` from the environment; credentials are never a
CLI argument, cache key, report field, or log value.

P0→P7 execute in dependency order and write atomically. A completed report does
not convert source gaps to zero evidence: exit `0` means all selected stages
completed, `2` means partial coverage (for example, missing STRING baseline or
external source), and `1` means an invalid input or required-stage failure.
Offline mode accepts only declared fixture responses; a missing fixture fails
instead of using the network.
