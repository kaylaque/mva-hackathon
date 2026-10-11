"""Unified MVA and chromosome instability gene sets and classifications.

Centralizes gene set definitions across Phase 0, Phase 1, Phase 2, and Phase 3:
- Canonical MVA causal genes (MVA1, MVA2, MVA3)
- Reported / extended MVA and chromosome instability genes (CENATAC, MAD1L1, SLF2, SMC5)
- Spindle Assembly Checkpoint (SAC) core components
- Cohesin and condensin complex components
- Secondary clinical context genes (aneuploidy surveillance, DDR, cancer predisposition)
"""

from __future__ import annotations

# Canonical MVA genes with established autosomal recessive causality
CANONICAL_MVA_GENES: frozenset[str] = frozenset({"BUB1B", "CEP57", "TRIP13"})

# Reported MVA / chromosomal instability causal genes
REPORTED_MVA_CA_GENES: frozenset[str] = frozenset({"CENATAC", "MAD1L1"})

# Atelis syndrome genes (overlapping mosaic variegated hyperploidy & chromosome fragility)
ATELIS_SYNDROME_GENES: frozenset[str] = frozenset({"SLF2", "SMC5"})

# Extended MVA / chromosome instability genes combining reported MVA and Atelis syndrome
EXTENDED_MVA_GENES: frozenset[str] = frozenset(REPORTED_MVA_CA_GENES | ATELIS_SYNDROME_GENES)

# The canonical 7 core MVA causal genes (3 canonical + 4 extended)
CORE_MVA_CAUSAL_GENES: frozenset[str] = frozenset(CANONICAL_MVA_GENES | EXTENDED_MVA_GENES)

# Core Spindle Assembly Checkpoint (SAC) components
SAC_CORE_GENES: frozenset[str] = frozenset({
    "BUB1", "BUB3", "MAD2L1", "MAD2L2", "CDC20", "AURKA", "AURKB",
    "PLK1", "CENPE", "TTK", "MPS1", "ZW10", "ZWILCH", "KNTC1", "KNL1", "NDC80"
})

# Structural maintenance of chromosomes / cohesin / condensin genes
CHROMOSOME_COHESIN_GENES: frozenset[str] = frozenset({
    "STAG2", "SMC1A", "SMC3", "RAD21", "NIPBL", "NCAPD2", "NCAPD3",
    "NCAPH", "NCAPG", "SMC2", "SMC4"
})

# Centrosome and primary cilia genes
CILIA_CENTROSOME_GENES: frozenset[str] = frozenset({
    "CEP290", "OFD1", "PKD1", "PKD2", "TMEM216", "TMEM67", "BBS1", "BBS2",
    "BBS4", "BBS7", "BBS9", "CEP152", "PCNT", "PLK4", "B9D1", "B9D2",
    "NPHP1", "NPHP3", "NPHP4", "INPP5E", "CC2D2A", "RPGRIP1L"
})

# Nephrocalcinosis and tubular ion transport genes
NEPHROCALCINOSIS_ION_GENES: frozenset[str] = frozenset({
    "CLDN16", "CLDN19", "CASR", "SLC12A1", "KCNJ1", "BSND", "SLC4A1",
    "SLC34A1", "SLC34A3", "ATP6V0A4", "ATP6V1B1", "CYP24A1", "AGXT",
    "GRHPR", "HOGA1", "SLC22A12", "SLC2A9", "MAGED2", "CLCN5", "OCRL",
    "VDR", "AP2S1", "GNA11"
})

# Pediatric cancer and rhabdomyosarcoma (RMS) predisposition genes
RMS_CANCER_GENES: frozenset[str] = frozenset({
    "TP53", "DICER1", "WT1", "NF1", "HRAS", "KRAS", "NRAS", "FGFR4",
    "PIK3CA", "FBXW7", "BCOR", "CTNNB1", "ALK", "MYOD1", "PAX3", "PAX7"
})

# Union of all MVA-relevant pathway genes
ALL_MVA_RELEVANT_GENES: frozenset[str] = (
    CANONICAL_MVA_GENES
    | REPORTED_MVA_CA_GENES
    | ATELIS_SYNDROME_GENES
    | SAC_CORE_GENES
    | CHROMOSOME_COHESIN_GENES
)

# Backward-compatible seed lists for Phase 3 drug evidence engine
# Note: CANONICAL_MVA_SEEDS includes all 7 core MVA causal genes
CANONICAL_MVA_SEEDS: list[str] = [
    "BUB1B",
    "CEP57",
    "TRIP13",
    "MAD1L1",
    "CENATAC",
    "SLF2",
    "SMC5",
]

# Hard-coded expert panel of 13 seeds (core MVA + interactome/checkpoint partners)
EXPANDED_MVA_NETWORK_SEEDS: list[str] = [
    "BUB1B",
    "CEP57",
    "TRIP13",
    "MAD1L1",
    "CENATAC",
    "SLF2",
    "SMC5",
    "TP53",
    "ATM",
    "NF1",
    "MAD2L2",
    "AURKA",
    "PLK1",
]

# Aliases matching network_visualization conventions
MVA_CANONICAL_GENES = CANONICAL_MVA_GENES
MVA_EXTENDED_GENES = EXTENDED_MVA_GENES


def classify_mva_gene_class(symbol: str) -> str:
    """Classify a gene symbol into CANONICAL, EXTENDED, or NONE."""
    if symbol in CANONICAL_MVA_GENES:
        return "CANONICAL"
    if symbol in EXTENDED_MVA_GENES:
        return "EXTENDED"
    return "NONE"
