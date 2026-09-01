#!/usr/bin/env python3
# ============================================================================
# scan_ap1_composite.py
# NFAT + AP-1 ortak taramasi ve NFAT:AP-1 kompozit element tespiti
#
# NFAT cogu zaman AP-1 ile birlikte calisir (klasik ARRE2 / IL-2 promotor
# composite). Fonksiyonel kompozitte AP-1 yarisi, NFAT GGAAA cekirdeginin
# ~3-15 bp yakininda bulunur. Burada NFAT hit'lerine <=20 bp mesafede
# AP-1 hit'i araniyor.
# ============================================================================

from Bio import motifs
from Bio.Seq import Seq

WINDOW_START_GENOMIC = 129599447
TSS_203 = 129599647
PSEUDOCOUNT = 0.5

NFAT_PWMS = {
    "Nfatc1": "pwm/MA0624.3.jaspar",
    "Nfatc2": "pwm/MA0152.3.jaspar",
}
AP1_PWMS = {
    "FOS::JUN":  "pwm/MA0099.3.jaspar",
    "BATF":      "pwm/MA1634.2.jaspar",
    "BATF::JUN": "pwm/MA0462.2.jaspar",
}

NFAT_FPR = 1e-4
AP1_FPR  = 1e-3   # AP-1 yarisi composite'te degenere olabilir -> biraz gevsek
COMPOSITE_MAX_GAP = 20  # bp, NFAT-AP1 kenar mesafesi

def read_fasta(path):
    s = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                s.append(line.strip())
    return "".join(s).upper()

promoter = read_fasta("seq/klrk1_promoter_plus.fa")
L = len(promoter)
seq_obj = Seq(promoter)

bg_counts = {b: promoter.count(b) for b in "ACGT"}
tot = sum(bg_counts.values())
background = {b: bg_counts[b]/tot for b in "ACGT"}

def scan(pwm_files, fpr):
    hits = []
    for tf, path in pwm_files.items():
        with open(path) as fh:
            m = motifs.read(fh, "jaspar")
        m.pseudocounts = PSEUDOCOUNT
        m.background = background
        pssm = m.pssm
        w = m.length
        dist = pssm.distribution(background=background, precision=10**4)
        thr = dist.threshold_fpr(fpr)
        for position, score in pssm.search(seq_obj, threshold=thr, both=True):
            if position >= 0:
                st = '+'; start0 = position
            else:
                st = '-'; start0 = position + L
            g_start = WINDOW_START_GENOMIC + start0
            g_end = g_start + w - 1
            rel = (score - pssm.min)/(pssm.max - pssm.min)
            site = promoter[start0:start0+w]
            if st == '-':
                site = str(Seq(site).reverse_complement())
            center = (g_start+g_end)/2.0
            hits.append({
                "tf": tf, "matrix": m.matrix_id, "g_start": g_start,
                "g_end": g_end, "center": center, "strand": st,
                "score": score, "rel": rel, "site": site,
                "dist_tss": TSS_203 - center,
            })
    return hits

nfat_hits = scan(NFAT_PWMS, NFAT_FPR)
ap1_hits  = scan(AP1_PWMS, AP1_FPR)

print(f"NFAT hit (p<{NFAT_FPR:g}): {len(nfat_hits)}")
print(f"AP-1 hit (p<{AP1_FPR:g}): {len(ap1_hits)}")

print("\n=== AP-1 HIT'LERI (TSS203'e gore) ===")
for h in sorted(ap1_hits, key=lambda x: x["dist_tss"]):
    print(f"  {h['tf']:10s} {h['site']:12s} rel={h['rel']:.2f} "
          f"strand={h['strand']} genomik={h['g_start']}-{h['g_end']} "
          f"TSS mesafe={h['dist_tss']:+.0f} bp")

# ---- Composite tespiti ------------------------------------------------------
print("\n=== NFAT:AP-1 KOMPOZIT ADAYLARI (kenar mesafesi <= %d bp) ===" % COMPOSITE_MAX_GAP)
def edge_gap(a, b):
    # iki motif arasi bosluk (ortismiyorlarsa pozitif)
    if a["g_end"] < b["g_start"]:
        return b["g_start"] - a["g_end"] - 1
    if b["g_end"] < a["g_start"]:
        return a["g_start"] - b["g_end"] - 1
    return 0  # ortusuyor

found = False
seen = set()
for n in nfat_hits:
    for a in ap1_hits:
        gap = edge_gap(n, a)
        if gap <= COMPOSITE_MAX_GAP:
            key = (round(n["center"]), round(a["center"]))
            if key in seen:
                continue
            seen.add(key)
            found = True
            print(f"  NFAT[{n['tf']}] {n['site']} ({n['g_start']}-{n['g_end']}) "
                  f"<-- gap={gap}bp --> AP1[{a['tf']}] {a['site']} "
                  f"({a['g_start']}-{a['g_end']})  | NFAT TSS mesafe={n['dist_tss']:+.0f} bp")
if not found:
    print("  Composite (NFAT-AP1 bitisik) bulunamadi.")
    # En yakin AP-1'i raporla
    print("\n  Bilgi: her NFAT site'a en yakin AP-1:")
    for n in sorted(nfat_hits, key=lambda x: x["dist_tss"]):
        if not ap1_hits:
            break
        nearest = min(ap1_hits, key=lambda a: edge_gap(n, a))
        print(f"    NFAT {n['site']} (TSS{n['dist_tss']:+.0f}) -> en yakin AP-1 "
              f"{nearest['tf']} {nearest['site']} mesafe={edge_gap(n,nearest)}bp")
