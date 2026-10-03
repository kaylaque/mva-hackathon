import csv
import os
from pathlib import Path
import sys

# ==============================================================================
# KONFIGURASI PATH
# ==============================================================================
BASE_DIR = Path("/home/nauvalrajwaa/proj/mva_hackathon/MVA")
INPUT_VCF = BASE_DIR / "data/interim/snpeff/annotated_snpeff.vcf"
OUTPUT_DIR = BASE_DIR / "data/interim/snpeff/filtered_csv"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

if not INPUT_VCF.exists():
  print(f"[ERROR] File VCF tidak ditemukan di: {INPUT_VCF}")
  sys.exit(1)

# Panel Gen Fokus Fenotipe (MVA, SAC, Ciliopathy, Nefrokalsinosis, Kanker)
TARGET_GENES = {
    # MVA & Mitotic Spindle Checkpoint
    "BUB1B",
    "CEP57",
    "TRIP13",
    "BUB1",
    "BUB3",
    "MAD2L1",
    "MAD2L2",
    "CDC20",
    "AURKA",
    "AURKB",
    "PLK1",
    "CENPE",
    "TTK",
    # Cohesin / Condensin
    "STAG2",
    "SMC1A",
    "SMC3",
    "RAD21",
    "NIPBL",
    "NCAPD2",
    # Sentrosom & Silia (Terkait fenotipe ginjal/ciliopathy)
    "CEP290",
    "OFD1",
    "PKD1",
    "PKD2",
    "TMEM216",
    "BBS1",
    # Nefrokalsinosis & Pengatur Ion Ginjal
    "CLDN16",
    "CLDN19",
    "CASR",
    "SLC12A1",
    "KCNJ1",
    "BSND",
    "SLC4A1",
    # Predisposisi Rhabdomyosarcoma & Tumor Embrional
    "TP53",
    "DICER1",
    "WT1",
    "NF1",
    "HRAS",
    "KRAS",
    "NRAS",
}


def parse_ann_field(info_str):
  """Mengekstrak anotasi HIGH dan MODERATE dari tag ANN di kolom INFO."""
  annotations = []
  for item in info_str.split(";"):
    if item.startswith("ANN="):
      ann_entries = item[4:].split(",")
      for entry in ann_entries:
        parts = entry.split("|")
        # Format ANN SnpEff standar:
        # 0: Allele | 1: Effect | 2: Impact | 3: Gene | 4: GeneID | 7: Biotype | 9: HGVS.c | 10: HGVS.p
        if len(parts) > 10:
          effect = parts[1]
          impact = parts[2]
          gene = parts[3]
          biotype = parts[7]
          hgvsc = parts[9]
          hgvsp = parts[10]

          if impact in ["HIGH", "MODERATE"]:
            annotations.append({
                "GENE": gene,
                "EFFECT": effect,
                "IMPACT": impact,
                "HGVS_C": hgvsc,
                "HGVS_P": hgvsp,
                "BIOTYPE": biotype,
            })
  return annotations


print(f"=== [1/3] Parsing langsung dari: {INPUT_VCF} ===")
all_damaging = []

with open(INPUT_VCF, "r") as f:
  for line in f:
    if line.startswith("#"):
      continue

    # Cek cepat sebelum pembedahan string untuk efisiensi eksekusi
    if "HIGH" not in line and "MODERATE" not in line:
      continue

    cols = line.strip().split("\t")
    if len(cols) < 8:
      continue

    chrom = cols[0]
    pos = cols[1]
    var_id = cols[2]
    ref = cols[3]
    alt = cols[4]
    qual = cols[5]
    filt = cols[6]
    info = cols[7]

    # Ambil data genotipe jika ada kolom sample
    genotype = ""
    if len(cols) >= 10:
      gt_format = cols[8].split(":")
      sample_data = cols[9].split(":")
      if "GT" in gt_format:
        gt_idx = gt_format.index("GT")
        genotype = sample_data[gt_idx]

    damaging_annots = parse_ann_field(info)

    # Ambil representasi unik per varian untuk menghindari duplikasi transkrip
    seen_genes = set()
    for ann in damaging_annots:
      gene_name = ann["GENE"]
      if gene_name not in seen_genes:
        seen_genes.add(gene_name)
        all_damaging.append({
            "CHROM": chrom,
            "POS": pos,
            "ID": var_id,
            "REF": ref,
            "ALT": alt,
            "QUAL": qual,
            "FILTER": filt,
            "GENOTYPE": genotype,
            "GENE": ann["GENE"],
            "EFFECT": ann["EFFECT"],
            "IMPACT": ann["IMPACT"],
            "HGVS_C": ann["HGVS_C"],
            "HGVS_P": ann["HGVS_P"],
            "BIOTYPE": ann["BIOTYPE"],
        })

print(
    f"Total varian HIGH/MODERATE terdeteksi: {len(all_damaging):,} anomali genik."
)

# ==============================================================================
# [2/3] Menyimpan Output CSV
# ==============================================================================
fieldnames = [
    "CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "GENOTYPE",
    "GENE",
    "EFFECT",
    "IMPACT",
    "HGVS_C",
    "HGVS_P",
    "BIOTYPE",
]

# 1. Simpan semua varian HIGH & MODERATE
file_all = OUTPUT_DIR / "all_high_moderate_variants.csv"
with open(file_all, "w", newline="") as f:
  writer = csv.DictWriter(f, fieldnames=fieldnames)
  writer.writeheader()
  writer.writerows(all_damaging)
print(f"-> Tersimpan: {file_all}")

# 2. Simpan varian pada gen target fenotipe pasien
phenotype_variants = [v for v in all_damaging if v["GENE"] in TARGET_GENES]
file_pheno = OUTPUT_DIR / "phenotype_focused_variants.csv"
with open(file_pheno, "w", newline="") as f:
  writer = csv.DictWriter(f, fieldnames=fieldnames)
  writer.writeheader()
  writer.writerows(phenotype_variants)
print(
    f"-> Tersimpan {len(phenotype_variants)} varian fokus fenotipe di:"
    f" {file_pheno}"
)

# 3. Hitung kandidat gen resesif / compound heterozygous (>= 2 varian damaging)
gene_counts = {}
for v in all_damaging:
  g = v["GENE"]
  gene_counts[g] = gene_counts.get(g, 0) + 1

recessive_genes = {g for g, count in gene_counts.items() if count >= 2 and g}
recessive_variants = [v for v in all_damaging if v["GENE"] in recessive_genes]

file_recessive = OUTPUT_DIR / "recessive_candidate_genes.csv"
with open(file_recessive, "w", newline="") as f:
  writer = csv.DictWriter(f, fieldnames=fieldnames)
  writer.writeheader()
  writer.writerows(recessive_variants)
print(
    f"-> Tersimpan {len(recessive_variants)} varian pada gen dengan multiple"
    f" hits di: {file_recessive}"
)

print("\n=== [3/3] SELESAI ===")