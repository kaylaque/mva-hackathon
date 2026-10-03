#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# SCRIPT 02: FUNCTIONAL ANNOTATION USING SNPEFF & SNPSIFT
# Environment: mva_snpeff (snpeff, snpsift, openjdk)
# ==============================================================================

# Define base paths
BASE_DIR="/home/nauvalrajwaa/proj/mva_hackathon/MVA"
INTERIM_DIR="${BASE_DIR}/data/interim"

FILTERED_VCF="${INTERIM_DIR}/filtered/pass_quality.vcf.gz"
SNPEFF_OUT_DIR="${INTERIM_DIR}/snpeff"
GENOME_BUILD="GRCh38.mane.1.0.ensembl"

# Check prerequisite file
if [ ! -f "${FILTERED_VCF}" ]; then
    echo "Error: Filtered VCF not found at ${FILTERED_VCF}."
    echo "Please execute 01_run_qc.sh first."
    exit 1
fi

mkdir -p "${SNPEFF_OUT_DIR}"

echo "=== [1/3] Checking SnpEff Reference Database ==="
if ! snpEff databases | grep -q "${GENOME_BUILD}"; then
    echo "Database ${GENOME_BUILD} not found locally. Downloading..."
    snpEff download "${GENOME_BUILD}"
else
    echo "Database ${GENOME_BUILD} is ready."
fi

echo "=== [2/3] Running SnpEff Annotation ==="
# Annotate variants using 8GB of memory allocation
snpEff -Xmx8g -v "${GENOME_BUILD}" "${FILTERED_VCF}" > "${SNPEFF_OUT_DIR}/annotated_snpeff.vcf"

echo "=== [3/3] Extracting Relevant Fields to TSV Table ==="
# Extract key columns into a clean tab-delimited file for tabular inspection
SnpSift extractFields "${SNPEFF_OUT_DIR}/annotated_snpeff.vcf" \
    CHROM POS ID REF ALT QUAL FILTER \
    "ANN[*].GENE" \
    "ANN[*].EFFECT" \
    "ANN[*].IMPACT" \
    "ANN[*].HGVS_C" \
    "ANN[*].HGVS_P" \
    "ANN[*].BIOTYPE" > "${SNPEFF_OUT_DIR}/snpeff_summary.tsv"

echo "SnpEff processing finished successfully."
echo "Annotated VCF : ${SNPEFF_OUT_DIR}/annotated_snpeff.vcf"
echo "TSV Summary   : ${SNPEFF_OUT_DIR}/snpeff_summary.tsv"
