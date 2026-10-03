#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# SCRIPT 03: CLINICAL ANNOTATION USING ENSEMBL VEP (OFFLINE ONLY)
# Environment: mva_vep (ensembl-vep, htslib)
# ==============================================================================

# Define base directory paths
BASE_DIR="/home/nauvalrajwaa/proj/mva_hackathon/MVA"
INTERIM_DIR="${BASE_DIR}/data/interim"

FILTERED_VCF="${INTERIM_DIR}/filtered/pass_quality.vcf.gz"
VEP_OUT_DIR="${INTERIM_DIR}/vep"
CACHE_DIR="${HOME}/.vep"
ASSEMBLY="GRCh38"

# 1. Verify input VCF file exists
if [ ! -f "${FILTERED_VCF}" ]; then
    echo "[ERROR] Filtered VCF not found at: ${FILTERED_VCF}"
    echo "Please execute 01_run_qc.sh first."
    exit 1
fi

mkdir -p "${VEP_OUT_DIR}"

# 2. Strict cache verification (No auto-download to prevent hanging/disk overflow)
echo "=== [1/2] Verifying Ensembl VEP Offline Cache Directory ==="
if [ ! -d "${CACHE_DIR}/homo_sapiens" ]; then
    echo "[ERROR] Local cache directory not found at: ${CACHE_DIR}/homo_sapiens"
    echo "Automatic download is permanently disabled in this script."
    echo "Please manually extract the cache archive into ${CACHE_DIR} before running."
    exit 1
fi

echo "Offline VEP cache verified successfully at: ${CACHE_DIR}"

# 3. Execute offline VEP annotation
echo "=== [2/2] Running Ensembl VEP Annotation (4 Threads) ==="
vep -i "${FILTERED_VCF}" \
    -o "${VEP_OUT_DIR}/vep_summary.tsv" \
    --dir_cache "${CACHE_DIR}" \
    --cache \
    --offline \
    --species homo_sapiens \
    --assembly "${ASSEMBLY}" \
    --tab \
    --symbol \
    --canonical \
    --sift b \
    --polyphen b \
    --af \
    --fork 4 \
    --force_overwrite

echo "Ensembl VEP annotation completed successfully."
echo "Output written to: ${VEP_OUT_DIR}/vep_summary.tsv"