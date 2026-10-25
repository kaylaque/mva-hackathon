#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# SCRIPT 01 (REVISED): QUALITY CONTROL (QC) AND VARIANT FILTERING
# Environment: mva_qc (bcftools, htslib)
# ==============================================================================

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW_DATA_DIR="${BASE_DIR}/data/input/patient"
INTERIM_DIR="${BASE_DIR}/data/interim"

# Set the actual VCF filename locally; keep patient sample codes out of Git.
INPUT_VCF="${RAW_DATA_DIR}/WGS_SUBJECT.vcf.gz"
QC_OUT_DIR="${INTERIM_DIR}/qc"
FILTERED_OUT_DIR="${INTERIM_DIR}/filtered"
FILTERED_VCF="${FILTERED_OUT_DIR}/pass_quality.vcf.gz"

mkdir -p "${QC_OUT_DIR}" "${FILTERED_OUT_DIR}"

echo "=== [1/3] Verifying VCF and Index Integrity ==="
if [ ! -f "${INPUT_VCF}.tbi" ]; then
    echo "Tabix index (.tbi) not found. Creating index..."
    tabix -p vcf "${INPUT_VCF}"
else
    echo "Tabix index confirmed: ${INPUT_VCF}.tbi"
fi

echo "=== [2/3] Generating Raw VCF Quality Control Statistics ==="
# Calculate summary metrics (Ts/Tv ratio, depth, variant counts)
bcftools stats "${INPUT_VCF}" > "${QC_OUT_DIR}/raw_vcf_stats.txt"
echo "QC metrics written to: ${QC_OUT_DIR}/raw_vcf_stats.txt"

# Display quick summary statistics in console (Field 3 = Metric Name, Field 4 = Actual Value)
echo "--- Raw VCF Summary Metrics ---"
grep "^SN" "${QC_OUT_DIR}/raw_vcf_stats.txt" | cut -f 3-4

echo "--- Transition / Transversion (Ts/Tv) Ratio ---"
grep "^TSTV" "${QC_OUT_DIR}/raw_vcf_stats.txt" | cut -f 3,5

echo "=== [3/3] Filtering for High-Quality Variants (PASS and QUAL >= 30) ==="
bcftools view -f PASS -i 'QUAL >= 30' "${INPUT_VCF}" -Oz -o "${FILTERED_VCF}"
tabix -p vcf "${FILTERED_VCF}"

# Count total high-quality variants retained after filtering
TOTAL_FILTERED=$(bcftools index -n "${FILTERED_VCF}")
echo "Filtering completed successfully."
echo "Total High-Quality Variants: ${TOTAL_FILTERED}"
echo "Filtered output: ${FILTERED_VCF}"