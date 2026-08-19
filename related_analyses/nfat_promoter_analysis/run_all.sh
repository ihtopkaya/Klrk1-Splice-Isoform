#!/bin/bash
# ============================================================================
# run_all.sh - Klrk1 promotor NFAT/AP-1 binding-site analizi (TAM PIPELINE)
# Tekrarlanabilirlik icin tum komutlar. Calistirma: bash run_all.sh
# ============================================================================
set -euo pipefail
BASE=~/GVHD_splicing_project/nfat_promoter_analysis
REF=~/GVHD_splicing_project/star_analysis/ref/chr6.fa   # GRCm39 chr6, Ensembl 111
cd "$BASE"
mkdir -p seq pwm results scripts

# --- 1. Promotor sekansi cikar (Klrk1-203 TSS=129,599,647, minus strand) ----
#     Pencere 6:129,599,447-129,601,735 (canonical TSS -2000 .. +200 bp)
samtools faidx "$REF" 6:129599447-129601735 > seq/klrk1_promoter_plus.fa
samtools faidx -i "$REF" 6:129599447-129601735 > seq/klrk1_promoter_sense.fa

# --- 2. JASPAR'dan NFAT + AP-1 PWM indir (CORE vertebrates) -----------------
JURL="https://jaspar.elixir.no/api/v1/matrix"
for ID in MA0624.3 MA0152.3 MA0625.3 MA1525.3 MA0606.3 \
          MA0099.3 MA1634.2 MA0462.2; do
  curl -sL "${JURL}/${ID}/?format=jaspar" -o "pwm/${ID}.jaspar"
done
# MA0624.3=Nfatc1 MA0152.3=Nfatc2 MA0625.3=NFATC3 MA1525.3=NFATC4 MA0606.3=Nfat5
# MA0099.3=FOS::JUN MA1634.2=BATF MA0462.2=BATF::JUN

# --- 3. NFAT taramasi (FIMO mantigi: p<1e-4, her iki strand) -----------------
python3 scripts/scan_nfat.py

# --- 4. NFAT + AP-1 + composite analizi -------------------------------------
python3 scripts/scan_ap1_composite.py

# --- 5. Promotor haritasi figuru --------------------------------------------
python3 scripts/make_figure.py

echo "TAMAM. Sonuclar: results/"
