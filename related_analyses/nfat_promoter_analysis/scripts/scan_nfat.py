#!/usr/bin/env python3
# ============================================================================
# scan_nfat.py
# Klrk1 promotorunde NFAT baglanma bolgesi (binding site) taramasi
#
# Mantik (FIMO ile ayni):
#   1. JASPAR PFM -> PSSM (log-odds), pseudocount + lokal background
#   2. p-value bazli esik (FPR) ile her iki strand taranir
#   3. Hit'ler genomik koordinata + TSS'e uzakliga cevrilir
#
# Gen: Klrk1 (ENSMUSG00000030149), chr6, MINUS strand
# Hedef transkript: Klrk1-203 (ENSMUST00000137660), TSS = 129,599,647
# Canonical:        Klrk1-201 (ENSMUST00000032252), TSS = 129,599,735
#
# Promotor penceresi: 6:129,599,447 - 129,601,735 (plus strand genomik)
#   seq[0] = genomik 129,599,447 (1-based)
# ============================================================================

import sys
from Bio import motifs
from Bio.Seq import Seq

# ---- Parametreler -----------------------------------------------------------
WINDOW_START_GENOMIC = 129599447   # seq[0] 1-based genomik koordinat
TSS_203  = 129599647               # Klrk1-203 TSS (minus strand 5' uc)
TSS_201  = 129599735               # Klrk1-201 canonical TSS
STRAND   = '-'                     # gen minus strand

PSEUDOCOUNT = 0.5
FPR_THRESHOLD = 1e-4               # FIMO varsayilan p-value esigi (1e-4)
REL_SCORE_REPORT = 0.80           # ayrica relative score >= 0.80 raporla

PWM_FILES = {
    "Nfatc1": "pwm/MA0624.3.jaspar",
    "Nfatc2": "pwm/MA0152.3.jaspar",
    "NFATC3": "pwm/MA0625.3.jaspar",
    "NFATC4": "pwm/MA1525.3.jaspar",
    "Nfat5":  "pwm/MA0606.3.jaspar",
}

PROMOTER_FA = "seq/klrk1_promoter_plus.fa"

# ---- Sekansi yukle ----------------------------------------------------------
def read_fasta(path):
    seq = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    return "".join(seq).upper()

promoter = read_fasta(PROMOTER_FA)
L = len(promoter)
print(f"Promotor uzunlugu: {L} bp")
print(f"Genomik pencere: 6:{WINDOW_START_GENOMIC}-{WINDOW_START_GENOMIC+L-1}")

# ---- Lokal background (promotor nukleotid kompozisyonu) --------------------
counts_bg = {b: promoter.count(b) for b in "ACGT"}
total_bg = sum(counts_bg.values())
background = {b: counts_bg[b] / total_bg for b in "ACGT"}
print("\nLokal background (promotor kompozisyonu):")
for b in "ACGT":
    print(f"  {b}: {background[b]:.4f}")

# ---- Genomik koordinat / TSS uzakligi hesabi -------------------------------
def genomic_coord(start0, motif_len):
    """0-based forward-seq start -> 1-based genomik start/end."""
    g_start = WINDOW_START_GENOMIC + start0
    g_end   = g_start + motif_len - 1
    return g_start, g_end

def dist_to_tss(g_start, g_end, tss):
    """Minus strand gen: upstream (yuksek koord) = negatif pozisyon."""
    center = (g_start + g_end) / 2.0
    return tss - center   # upstream -> negatif

# ---- Her PWM icin tara ------------------------------------------------------
all_hits = []
seq_obj = Seq(promoter)

print("\n" + "=" * 78)
for tf_name, pwm_path in PWM_FILES.items():
    with open(pwm_path) as fh:
        m = motifs.read(fh, "jaspar")
    w = m.length
    m.pseudocounts = PSEUDOCOUNT
    m.background = background
    pssm = m.pssm

    # Skor araligi -> relative score hesabi
    max_score = pssm.max
    min_score = pssm.min

    # p-value esigi (FPR) -> skor esigi
    dist = pssm.distribution(background=background, precision=10**4)
    thr_fpr = dist.threshold_fpr(FPR_THRESHOLD)

    consensus = str(m.consensus)
    print(f"\n### {tf_name}  ({m.matrix_id})  consensus={consensus}  len={w}")
    print(f"    Skor araligi: [{min_score:.2f}, {max_score:.2f}] | "
          f"p<{FPR_THRESHOLD:g} esigi = {thr_fpr:.3f} bit")

    # Her iki strand tara (FPR esigiyle)
    hits = []
    for position, score in pssm.search(seq_obj, threshold=thr_fpr, both=True):
        if position >= 0:
            site_strand = '+'
            start0 = position
        else:
            site_strand = '-'
            start0 = position + L
        g_start, g_end = genomic_coord(start0, w)
        rel = (score - min_score) / (max_score - min_score)
        d203 = dist_to_tss(g_start, g_end, TSS_203)
        d201 = dist_to_tss(g_start, g_end, TSS_201)
        # motif dizisi
        site_seq = promoter[start0:start0+w]
        if site_strand == '-':
            site_seq = str(Seq(site_seq).reverse_complement())
        hits.append({
            "tf": tf_name, "matrix": m.matrix_id, "consensus": consensus,
            "g_start": g_start, "g_end": g_end, "site_strand": site_strand,
            "score": score, "rel_score": rel,
            "dist_TSS203": d203, "dist_TSS201": d201,
            "site_seq": site_seq, "thr_fpr": thr_fpr,
        })
    hits.sort(key=lambda h: -h["score"])
    print(f"    p<{FPR_THRESHOLD:g} ile hit sayisi: {len(hits)}")
    for h in hits:
        gene_strand_note = "SENSE(gen yonu)" if h["site_strand"] == STRAND else "antisense"
        print(f"      {h['site_seq']:10s} skor={h['score']:6.2f} "
              f"rel={h['rel_score']:.2f}  DNA-strand={h['site_strand']} ({gene_strand_note})  "
              f"genomik={h['g_start']}-{h['g_end']}  "
              f"TSS203 mesafe={h['dist_TSS203']:+.0f} bp")
    all_hits.extend(hits)

# ---- Ozet tablo (CSV) -------------------------------------------------------
print("\n" + "=" * 78)
print(f"TOPLAM HIT (p<{FPR_THRESHOLD:g}, tum NFAT matrisleri): {len(all_hits)}")

import csv
with open("results/nfat_hits.csv", "w", newline="") as fh:
    wr = csv.writer(fh)
    wr.writerow(["TF","matrix","consensus","site_seq","DNA_strand",
                 "genomic_start","genomic_end","score","rel_score",
                 "dist_to_TSS203_bp","dist_to_TSS201_bp","pval_threshold_FPR"])
    for h in sorted(all_hits, key=lambda x: x["dist_TSS203"]):
        wr.writerow([h["tf"], h["matrix"], h["consensus"], h["site_seq"],
                     h["site_strand"], h["g_start"], h["g_end"],
                     f"{h['score']:.3f}", f"{h['rel_score']:.3f}",
                     f"{h['dist_TSS203']:.0f}", f"{h['dist_TSS201']:.0f}",
                     f"{FPR_THRESHOLD:g}"])
print("CSV yazildi: results/nfat_hits.csv")

# ---- Proksimal promotor (-500..+100) ozeti ---------------------------------
print("\n=== PROKSIMAL PROMOTOR (-500 ile +100 bp, TSS203'e gore) ===")
prox = [h for h in all_hits if -500 <= h["dist_TSS203"] <= 100]
prox.sort(key=lambda x: x["dist_TSS203"])
if prox:
    for h in prox:
        print(f"  {h['tf']:8s} {h['site_seq']:10s} "
              f"TSS mesafe={h['dist_TSS203']:+.0f} bp  rel={h['rel_score']:.2f}  "
              f"strand={h['site_strand']}")
else:
    print("  (proksimal bolgede p<1e-4 hit yok)")
