import csv
from pathlib import Path
import re
import sys

# ==============================================================================
# KONFIGURASI PATH
# ==============================================================================
BASE_DIR = Path("/home/nauvalrajwaa/proj/mva_hackathon/MVA")
INPUT_TSV = BASE_DIR / "data/interim/vep/vep_summary.tsv"
OUTPUT_DIR = BASE_DIR / "data/interim/vep/filtered_csv"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

if not INPUT_TSV.exists():
  print(f"[ERROR] File input VEP tidak ditemukan di: {INPUT_TSV}")
  sys.exit(1)

# ==============================================================================
# PANEL FENOTIPE KLINIS PASIEN BERDASARKAN KATEGORI JALUR BIOLOGIS
# ==============================================================================
PHENOTYPE_PANELS = {
    "MVA_SAC_CORE": {
        "BUB1B",
        "CEP57",
        "TRIP13",
        "BUB1",
        "BUB3",
        "MAD1L1",
        "MAD2L1",
        "MAD2L2",
        "CDC20",
        "AURKA",
        "AURKB",
        "PLK1",
        "CENPE",
        "TTK",
        "MPS1",
        "ZW10",
        "ZWILCH",
        "KNTC1",
    },
    "CHROMOSOME_COHESIN": {
        "STAG2",
        "SMC1A",
        "SMC3",
        "RAD21",
        "NIPBL",
        "NCAPD2",
        "NCAPD3",
        "NCAPH",
        "NCAPG",
        "SMC2",
        "SMC4",
    },
    "CILIA_CENTROSOME_KIDNEY": {
        "CEP290",
        "OFD1",
        "PKD1",
        "PKD2",
        "TMEM216",
        "TMEM67",
        "BBS1",
        "BBS2",
        "BBS4",
        "BBS7",
        "BBS9",
        "CEP152",
        "PCNT",
        "PLK4",
        "B9D1",
        "B9D2",
        "NPHP1",
        "NPHP3",
        "NPHP4",
        "INPP5E",
        "CC2D2A",
        "RPGRIP1L",
    },
    "NEPHROCALCINOSIS_ION": {
        "CLDN16",
        "CLDN19",
        "CASR",
        "SLC12A1",
        "KCNJ1",
        "BSND",
        "SLC4A1",
        "SLC34A1",
        "SLC34A3",
        "ATP6V0A4",
        "ATP6V1B1",
        "CYP24A1",
        "AGXT",
        "GRHPR",
        "HOGA1",
        "SLC22A12",
        "SLC2A9",
        "MAGED2",
        "CLCN5",
        "OCRL",
        "VDR",
        "AP2S1",
        "GNA11",
    },
    "RMS_CANCER_PREDISPOSITION": {
        "TP53",
        "DICER1",
        "WT1",
        "NF1",
        "HRAS",
        "KRAS",
        "NRAS",
        "FGFR4",
        "PIK3CA",
        "FBXW7",
        "BCOR",
        "CTNNB1",
        "ALK",
        "MYOD1",
        "PAX3",
        "PAX7",
    },
}

ALL_PHENOTYPE_GENES = {
    gene for panel in PHENOTYPE_PANELS.values() for gene in panel
}

# Kamus konversi asam amino 1-letter ke 3-letter code
AA_MAP = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "E": "Glu",
    "Q": "Gln",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
    "*": "*",
}


def format_hgvs_p(amino_acids, protein_pos):
  """Mengonversi format 'L/*' atau 'N/K' menjadi 'p.Leu737*' atau 'p.Asn1002Lys'."""
  if not amino_acids or amino_acids == "-" or not protein_pos or protein_pos == "-":
    return "-"
  pos_clean = protein_pos.split("-")[0]
  if "/" in amino_acids:
    ref_aa, alt_aa = amino_acids.split("/")[:2]
    ref_3 = AA_MAP.get(ref_aa, ref_aa)
    alt_3 = AA_MAP.get(alt_aa, alt_aa)
    return f"p.{ref_3}{pos_clean}{alt_3}"
  return f"p.{amino_acids}{pos_clean}"


def format_hgvs_c(cds_pos, ref, alt):
  """Mengonversi koordinat CDS menjadi notasi c.2210T>G."""
  if not cds_pos or cds_pos == "-" or not ref or not alt or ref == "-" or alt == "-":
    return "-"
  return f"c.{cds_pos}{ref}>{alt}"


def parse_score_prediction(val):
  """Memecah 'probably_damaging(0.997)' menjadi ('probably_damaging', '0.997')."""
  if not val or val == "-":
    return "-", "-"
  match = re.match(r"([a-zA-Z_]+)\(([\d\.]+)\)", val)
  if match:
    return match.group(1), match.group(2)
  return val, "-"


def get_phenotype_category(symbol):
  """Mengidentifikasi kategori fenotipe klinis dari gen."""
  categories = []
  for name, genes in PHENOTYPE_PANELS.items():
    if symbol in genes:
      categories.append(name)
  return "; ".join(categories) if categories else "Genome-wide"


# ==============================================================================
# ACMG / AMP 2015 CLASSIFICATION ENGINE
# ==============================================================================
def evaluate_acmg_and_priority(row, gene_hit_counts):
  """Mengevaluasi kriteria ACMG/AMP 2015 dan menghitung skor prioritas klinis."""
  criteria = []
  symbol = row["SYMBOL"]
  consequence = row["CONSEQUENCE"]
  impact = row["IMPACT"]
  clin_sig = row["CLIN_SIG"].lower()
  sift_pred = row["SIFT_PRED"].lower()
  pp_pred = row["POLYPHEN_PRED"].lower()
  af_str = row["AF_POPULATION"]
  protein_pos = row["PROTEIN_POS"]

  # Parse frekuensi alel (AF)
  af_val = 0.0
  is_rare = True
  is_common_benign = False
  if af_str != "-" and af_str != "":
    try:
      af_val = float(af_str)
      if af_val >= 0.05:  # >= 5% di populasi
        is_common_benign = True
      elif af_val >= 0.01:  # 1% - 5%
        is_rare = False
      else:
        is_rare = True
    except ValueError:
      is_rare = True

  # 1. BA1 & BS1 (Kriteria Benign Frekuensi Tinggi)
  if is_common_benign:
    criteria.append("BA1")
  elif not is_rare:
    criteria.append("BS1")

  # 2. PVS1 (Null variant di gen LoF: stop_gained, frameshift, splice site)
  is_null_variant = any(
      term in consequence
      for term in [
          "stop_gained",
          "frameshift_variant",
          "splice_acceptor_variant",
          "splice_donor_variant",
      ]
  )
  if is_null_variant and impact == "HIGH":
    criteria.append("PVS1")

  # 3. PM2 (Frekuensi ultra-langka atau tidak tercatat di kontrol populasi)
  if is_rare and not is_common_benign:
    criteria.append("PM2")

  # 4. PM3 (Kandidat Resesif / Compound Heterozygous in trans: >= 2 mutasi merusak pada gen yang sama)
  if gene_hit_counts.get(symbol, 0) >= 2:
    criteria.append("PM3")

  # 5. PM1 (Terletak di Hotspot / Domain Fungsional Terkonversi Kritis)
  try:
    pos_num = int(protein_pos.split("-")[0])
    # Domain kinase BUB1B (residu 725 - 1008)
    if symbol == "BUB1B" and (725 <= pos_num <= 1008):
      criteria.append("PM1")
    # DNA-binding domain TP53 (residu 102 - 292)
    elif symbol == "TP53" and (102 <= pos_num <= 292):
      criteria.append("PM1")
  except ValueError:
    pass

  # 6. PP3 vs BP4 (In silico prediction consensus)
  is_sift_damaging = "deleterious" in sift_pred
  is_pp_damaging = (
      "damaging" in pp_pred
      or "probably_damaging" in pp_pred
      or "possibly_damaging" in pp_pred
  )

  if is_sift_damaging and is_pp_damaging:
    criteria.append("PP3")
  elif "tolerated" in sift_pred and "benign" in pp_pred:
    criteria.append("BP4")

  # 7. PP5 vs BP6 (ClinVar Database Assertion)
  if "pathogenic" in clin_sig:
    criteria.append("PP5")
  elif "benign" in clin_sig and "likely" not in clin_sig:
    criteria.append("BP6")

  # --------------------------------------------------------------------------
  # Logika Keputusan ACMG (Richards et al., 2015)
  # --------------------------------------------------------------------------
  has_ba1 = "BA1" in criteria
  has_bs1 = "BS1" in criteria
  has_bp4 = "BP4" in criteria
  has_bp6 = "BP6" in criteria
  has_pvs1 = "PVS1" in criteria
  has_pp5 = "PP5" in criteria

  pm_count = sum(1 for c in criteria if c.startswith("PM"))
  pp_count = sum(1 for c in criteria if c.startswith("PP"))
  bp_count = sum(1 for c in criteria if c.startswith("BP"))

  if has_ba1 or (has_bs1 and bp_count >= 1):
    classification = "Benign"
  elif has_bp6 or bp_count >= 2:
    classification = "Likely Benign"
  elif has_pvs1 and (has_pp5 or pm_count >= 1 or pp_count >= 2):
    classification = "Pathogenic"
  elif has_pp5 and (pm_count >= 1 or pp_count >= 2):
    classification = "Pathogenic"
  elif has_pvs1 or pm_count >= 3 or (pm_count >= 2 and pp_count >= 1):
    classification = "Likely Pathogenic"
  elif has_pp5:
    classification = "Likely Pathogenic"
  elif pm_count >= 1 and pp_count >= 2:
    classification = "Likely Pathogenic"
  else:
    classification = "Uncertain Significance (VUS)"

  # --------------------------------------------------------------------------
  # Clinical Priority Ranking Score (Untuk Menentukan Top Rank Track 1)
  # --------------------------------------------------------------------------
  priority_score = 0
  # Komponen Klasifikasi
  if classification == "Pathogenic":
    priority_score += 60
  elif classification == "Likely Pathogenic":
    priority_score += 40
  elif classification == "Uncertain Significance (VUS)":
    priority_score += 10
  elif classification == "Benign":
    priority_score -= 50

  # Komponen Panel Fenotipe
  if symbol in PHENOTYPE_PANELS["MVA_SAC_CORE"]:
    priority_score += 45
  elif (
      symbol in PHENOTYPE_PANELS["NEPHROCALCINOSIS_ION"]
      or symbol in PHENOTYPE_PANELS["CILIA_CENTROSOME_KIDNEY"]
  ):
    priority_score += 30
  elif symbol in PHENOTYPE_PANELS["RMS_CANCER_PREDISPOSITION"]:
    priority_score += 25
  elif symbol in PHENOTYPE_PANELS["CHROMOSOME_COHESIN"]:
    priority_score += 20

  # Komponen Model Pewarisan Resesif (Double Hit)
  if gene_hit_counts.get(symbol, 0) >= 2:
    priority_score += 30

  # Komponen Keparahan Mutasi & Kelangkaan
  if impact == "HIGH":
    priority_score += 20
  elif impact == "MODERATE":
    priority_score += 10

  if is_rare:
    priority_score += 15
  if is_sift_damaging and is_pp_damaging:
    priority_score += 10

  return "; ".join(criteria), classification, priority_score


# ==============================================================================
# [1/4] MEMBACA DAN MEMILIH TRANSKRIP REPRESENTATIF TERBAIK (PASS 1)
# ==============================================================================
print(f"=== [1/4] Membaca data VEP TSV dari: {INPUT_TSV} ===")

header_line = None
col_idx = {}
variants_best_transcript = {}
gene_hit_counter = {}

with open(INPUT_TSV, "r") as f:
  for line in f:
    if line.startswith("##"):
      continue

    if line.startswith("#Uploaded_variation") or (
        header_line is None and line.startswith("#")
    ):
      header_line = line.lstrip("#").strip().split("\t")
      col_idx = {name: i for i, name in enumerate(header_line)}
      continue

    if "HIGH" not in line and "MODERATE" not in line:
      continue

    cols = line.strip().split("\t")
    if len(cols) < len(col_idx):
      continue

    impact = cols[col_idx["IMPACT"]] if "IMPACT" in col_idx else "-"
    if impact not in ["HIGH", "MODERATE"]:
      continue

    uploaded_var = cols[col_idx["Uploaded_variation"]]
    symbol = cols[col_idx.get("SYMBOL", col_idx["Gene"])]
    if not symbol or symbol == "-":
      symbol = cols[col_idx["Gene"]]

    canonical = cols[col_idx["CANONICAL"]] if "CANONICAL" in col_idx else "-"
    protein_pos = (
        cols[col_idx["Protein_position"]]
        if "Protein_position" in col_idx
        else "-"
    )
    clin_sig = cols[col_idx["CLIN_SIG"]] if "CLIN_SIG" in col_idx else "-"
    sift_raw = cols[col_idx["SIFT"]] if "SIFT" in col_idx else "-"
    polyphen_raw = cols[col_idx["PolyPhen"]] if "PolyPhen" in col_idx else "-"

    # Bobot prioritas pemilihan transkrip representatif terbaik
    weight = 0
    if canonical == "YES":
      weight += 100
    if impact == "HIGH":
      weight += 50
    elif impact == "MODERATE":
      weight += 30
    if clin_sig != "-":
      weight += 20
    if protein_pos != "-":
      weight += 10
    if sift_raw != "-":
      weight += 5
    if polyphen_raw != "-":
      weight += 5

    var_key = (uploaded_var, symbol)

    if (
        var_key not in variants_best_transcript
        or weight > variants_best_transcript[var_key]["_weight"]
    ):
      # Parse koordinat genomik
      chrom, pos, ref, alt = "-", "-", "-", "-"
      var_parts = uploaded_var.split("_")
      if len(var_parts) >= 3:
        chrom = var_parts[0]
        pos = var_parts[1]
        alleles = var_parts[2].split("/")
        ref = alleles[0]
        alt = alleles[1] if len(alleles) > 1 else "-"
      else:
        loc = cols[col_idx["Location"]].split(":")
        chrom = loc[0]
        pos = loc[1].split("-")[0] if len(loc) > 1 else "-"
        alt = cols[col_idx["Allele"]]

      sift_pred, sift_score = parse_score_prediction(sift_raw)
      pp_pred, pp_score = parse_score_prediction(polyphen_raw)

      amino_acids = (
          cols[col_idx["Amino_acids"]] if "Amino_acids" in col_idx else "-"
      )
      codons = cols[col_idx["Codons"]] if "Codons" in col_idx else "-"
      cds_pos = cols[col_idx["CDS_position"]] if "CDS_position" in col_idx else "-"
      cdna_pos = (
          cols[col_idx["cDNA_position"]] if "cDNA_position" in col_idx else "-"
      )
      feature = cols[col_idx["Feature"]] if "Feature" in col_idx else "-"
      consequence = (
          cols[col_idx["Consequence"]] if "Consequence" in col_idx else "-"
      )
      existing_id = (
          cols[col_idx["Existing_variation"]]
          if "Existing_variation" in col_idx
          else "-"
      )
      allele_freq = cols[col_idx["AF"]] if "AF" in col_idx else "-"

      hgvsp = format_hgvs_p(amino_acids, protein_pos)
      hgvsc = format_hgvs_c(cds_pos, ref, alt)

      variants_best_transcript[var_key] = {
          "_weight": weight,
          "CHROM": chrom,
          "POS": pos,
          "REF": ref,
          "ALT": alt,
          "VARIANT_ID": uploaded_var,
          "SYMBOL": symbol,
          "GENE_ID": cols[col_idx["Gene"]],
          "FEATURE_TRANSCRIPT": feature,
          "CANONICAL": canonical,
          "CONSEQUENCE": consequence,
          "IMPACT": impact,
          "PROTEIN_POS": protein_pos,
          "AMINO_ACIDS": amino_acids,
          "CODONS": codons,
          "CDS_POS": cds_pos,
          "CDNA_POS": cdna_pos,
          "HGVS_P": hgvsp,
          "HGVS_C": hgvsc,
          "EXISTING_VARIATION": existing_id,
          "CLIN_SIG": clin_sig,
          "SIFT_PRED": sift_pred,
          "SIFT_SCORE": sift_score,
          "POLYPHEN_PRED": pp_pred,
          "POLYPHEN_SCORE": pp_score,
          "AF_POPULATION": allele_freq,
          "PHENOTYPE_PANEL": get_phenotype_category(symbol),
      }

# Hitung jumlah mutasi merusak per gen di seluruh genom untuk deteksi Compound Heterozygote
for (uploaded_var, symbol), val in variants_best_transcript.items():
  gene_hit_counter[symbol] = gene_hit_counter.get(symbol, 0) + 1

print(
    f"Total varian unik HIGH & MODERATE teridentifikasi:"
    f" {len(variants_best_transcript):,}"
)

# ==============================================================================
# [2/4] EVALUASI ACMG & SCORING PRIORITAS KLINIS (PASS 2)
# ==============================================================================
print("=== [2/4] Menghitung Klasifikasi ACMG dan Skor Prioritas Klinis ===")

all_enriched_rows = []
for var_key, item in variants_best_transcript.items():
  item.pop("_weight", None)
  symbol = item["SYMBOL"]

  # Tambahkan flag model pewarisan
  hit_count = gene_hit_counter.get(symbol, 1)
  item["INHERITANCE_MODEL"] = (
      f"Compound_Het_Candidate ({hit_count}_hits)"
      if hit_count >= 2
      else "Single_Hit"
  )

  # Hitung ACMG dan skor prioritas
  criteria_str, acmg_class, p_score = evaluate_acmg_and_priority(
      item, gene_hit_counter
  )
  item["ACMG_CRITERIA"] = criteria_str
  item["ACMG_CLASSIFICATION"] = acmg_class
  item["PRIORITY_SCORE"] = p_score

  all_enriched_rows.append(item)

# ==============================================================================
# [3/4] MENYUSUN STRUKTUR KOLOM & FILE OUTPUT CSV
# ==============================================================================
fieldnames = [
    "PRIORITY_SCORE",
    "SYMBOL",
    "HGVS_P",
    "HGVS_C",
    "ACMG_CLASSIFICATION",
    "ACMG_CRITERIA",
    "INHERITANCE_MODEL",
    "PHENOTYPE_PANEL",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "CONSEQUENCE",
    "IMPACT",
    "CLIN_SIG",
    "SIFT_PRED",
    "SIFT_SCORE",
    "POLYPHEN_PRED",
    "POLYPHEN_SCORE",
    "AF_POPULATION",
    "EXISTING_VARIATION",
    "FEATURE_TRANSCRIPT",
    "CANONICAL",
    "PROTEIN_POS",
    "AMINO_ACIDS",
    "CODONS",
    "CDS_POS",
    "CDNA_POS",
    "GENE_ID",
    "VARIANT_ID",
]

# 1. Master File: Seluruh Varian HIGH & MODERATE Genom
file_all = OUTPUT_DIR / "vep_all_high_moderate.csv"
all_enriched_rows.sort(
    key=lambda x: (
        -x["PRIORITY_SCORE"],
        (
            x["CHROM"].replace("chr", "").zfill(2)
            if x["CHROM"].isdigit()
            else x["CHROM"]
        ),
        int(x["POS"]) if x["POS"].isdigit() else 0,
    )
)

with open(file_all, "w", newline="") as f:
  writer = csv.DictWriter(f, fieldnames=fieldnames)
  writer.writeheader()
  writer.writerows(all_enriched_rows)
print(
    f"-> 1. Master HIGH/MODERATE   : {file_all} ({len(all_enriched_rows):,}"
    " varian)"
)

# 2. Panel Fenotipe Klinis Pasien (MVA, Silia, Ginjal, Tumor) - Diurutkan dari Skor Tertinggi
file_pheno = OUTPUT_DIR / "vep_phenotype_focused.csv"
pheno_rows = [
    v for v in all_enriched_rows if v["SYMBOL"] in ALL_PHENOTYPE_GENES
]
pheno_rows.sort(key=lambda x: (-x["PRIORITY_SCORE"], x["SYMBOL"]))

with open(file_pheno, "w", newline="") as f:
  writer = csv.DictWriter(f, fieldnames=fieldnames)
  writer.writeheader()
  writer.writerows(pheno_rows)
print(f"-> 2. Phenotype Focused      : {file_pheno} ({len(pheno_rows)} varian)")

# 3. Model Pewarisan Resesif (Gen dengan >= 2 mutasi merusak)
recessive_genes = {g for g, cnt in gene_hit_counter.items() if cnt >= 2}
recessive_rows = [v for v in all_enriched_rows if v["SYMBOL"] in recessive_genes]
recessive_rows.sort(
    key=lambda x: (
        -x["PRIORITY_SCORE"],
        x["SYMBOL"],
        int(x["POS"]) if x["POS"].isdigit() else 0,
    )
)

file_recessive = OUTPUT_DIR / "vep_recessive_candidates.csv"
with open(file_recessive, "w", newline="") as f:
  writer = csv.DictWriter(f, fieldnames=fieldnames)
  writer.writeheader()
  writer.writerows(recessive_rows)
print(
    f"-> 3. Recessive Candidates   : {file_recessive} ({len(recessive_rows):,}"
    f" varian pada {len(recessive_genes)} gen)"
)

# 4. ClinVar & ACMG Pathogenic (Varian dengan bukti klinis definitif)
file_pathogenic = OUTPUT_DIR / "vep_pathogenic_variants.csv"
patho_rows = [
    v
    for v in all_enriched_rows
    if v["ACMG_CLASSIFICATION"] in ["Pathogenic", "Likely Pathogenic"]
    or any(
        term in v["CLIN_SIG"].lower()
        for term in ["pathogenic", "likely_pathogenic"]
    )
]
patho_rows.sort(key=lambda x: (-x["PRIORITY_SCORE"], x["SYMBOL"]))

with open(file_pathogenic, "w", newline="") as f:
  writer = csv.DictWriter(f, fieldnames=fieldnames)
  writer.writeheader()
  writer.writerows(patho_rows)
print(
    f"-> 4. Pathogenic / Likely P  : {file_pathogenic} ({len(patho_rows)}"
    " varian)"
)

# 5. Top Clinical Rankings untuk Submission Track 1
file_top = OUTPUT_DIR / "vep_top_clinical_rankings.csv"
top_rows = pheno_rows[:25]
with open(file_top, "w", newline="") as f:
  writer = csv.DictWriter(f, fieldnames=fieldnames)
  writer.writeheader()
  writer.writerows(top_rows)
print(
    f"-> 5. Top 25 Clinical Ranks  : {file_top} (Kandidat Utama Siap"
    " Dikumpulkan)"
)

print("\n=== [4/4] SELESAI! Seluruh anotasi ACMG dan ranking siap dipakai ===")