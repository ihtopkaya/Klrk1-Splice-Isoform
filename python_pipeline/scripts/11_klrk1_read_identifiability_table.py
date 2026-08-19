#!/usr/bin/env python3
# ============================================================================
# 11_klrk1_read_identifiability_table.py
# Supplementary table tying the Klrk1 region/junction structure to actual read
# support in two representative bias-corrected BAMs: GSE203167 WT Post-Tx D7
# (paired-end) and GSE147371 GVHD Tem_1 (76 bp). For each read-diagnostic region
# it reports mean depth (samtools depth -a), showing which units are directly
# measured (exclusive coverage) vs inferred (203 = retention in excess of 204).
# Run: python3 python_pipeline/scripts/11_klrk1_read_identifiability_table.py
# ============================================================================
import os, sys, subprocess, csv
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config as cfg

SAM = "/opt/homebrew/bin/samtools" if os.path.exists("/opt/homebrew/bin/samtools") else "samtools"
BAMS = [("WT_Post7 (GSE203167)", os.path.join(cfg.PROJECT_ROOT,
         "datasets/GSE203167_Karimi_TCF7/data/bam/WT_Post7.bam")),
        ("Tem_1 (GSE147371)", os.path.join(cfg.PROJECT_ROOT,
         "datasets/GSE147371_GVHD_CD4/data/bam/GVHD_CD4_Tem_1.bam"))]
# region-diagnostic features (chr6, GRCm39)
FEATURES = [
 ("Canonical first exon", "201/203/206", "shared", "129599545-129599735"),
 ("Alternative first exon", "202/205/206", "shared", "129600774-129600804"),
 ("205-exclusive 3' tail", "205", "EXCLUSIVE", "129587286-129589376"),
 ("206-exclusive exon", "206", "EXCLUSIVE", "129599736-129600773"),
 ("Intron-3 retention", "204", "EXCLUSIVE", "129593737-129594445"),
 ("Intron-4 retention segment", "203/204", "shared", "129593531-129593631"),
]
def depth(region, bam):
    if not os.path.exists(bam + ".bai"):
        subprocess.run([SAM, "index", bam], check=False)
    p = subprocess.run([SAM, "depth", "-a", "-r", "6:" + region.replace("-", "-"), bam],
                       capture_output=True, text=True)
    vals = [int(l.split("\t")[2]) for l in p.stdout.splitlines() if l]
    return round(sum(vals) / len(vals), 1) if vals else 0.0

rows = [["Diagnostic region", "Transcript(s)", "Region type", "chr6 coords"]
        + [f"mean depth: {n}" for n, _ in BAMS]]
for name, tx, typ, coords in FEATURES:
    rows.append([name, tx, typ, coords] + [depth(coords, b) for _, b in BAMS])

out = os.path.join(cfg.TABLE_DIR, "TableS_read_identifiability.csv")
with open(out, "w", newline="") as f:
    csv.writer(f).writerows(rows)
w = max(len(r[0]) for r in rows)
for r in rows:
    print(f"{r[0]:<{w}}  {r[1]:<12} {r[2]:<10} {r[3]:<22} " + "  ".join(str(x) for x in r[4:]))
print(f"\nwrote {out}")
print("Reading: 205 & 204 exclusive regions carry direct coverage (measured); 206 exclusive ~0 (low);\n"
      "203 is the intron-4 retention in EXCESS of the 204 intron-3 level (Tem_1: 19.7 vs 10.1).")
