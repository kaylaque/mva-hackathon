#!/usr/bin/env bash
# Download public GRCh38 resources into data/input/refs. Never stores patient data.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REFS="${ROOT}/data/input/refs"

fetch() {
  local url="$1" destination="$2"
  mkdir -p "$(dirname "${destination}")"
  if [[ -s "${destination}" ]]; then
    printf '[exists] %s\n' "${destination}"
  else
    curl --fail --location --continue-at - --output "${destination}" "${url}"
  fi
}

fetch "https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" \
  "${REFS}/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz" \
  "${REFS}/clinvar/clinvar.vcf.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi" \
  "${REFS}/clinvar/clinvar.vcf.gz.tbi"
fetch "https://zenodo.org/records/7072866/files/revel-v1.3_all_chromosomes.zip" \
  "${REFS}/revel/revel-v1.3_all_chromosomes.zip"
fetch "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz" \
  "${REFS}/alphamissense/AlphaMissense_hg38.tsv.gz"
fetch "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz.tbi" \
  "${REFS}/alphamissense/AlphaMissense_hg38.tsv.gz.tbi"
fetch "https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz" \
  "${REFS}/vep/homo_sapiens_vep_116_GRCh38.tar.gz"
fetch "https://g-879a9f.f5dc97.75bc.dn.glob.us/data/2602_hg38.zip" \
  "${REFS}/exomiser/2602_hg38.zip"
fetch "https://g-879a9f.f5dc97.75bc.dn.glob.us/data/2602_phenotype.zip" \
  "${REFS}/exomiser/2602_phenotype.zip"

printf 'Downloads complete. Decompress/index resources explicitly for the tool that uses them.\n'
