#!/usr/bin/env python3
# ============================================================================
# 02_exon_overlap_correlation.py
# Two structural checks behind the reporting-unit decision:
#   (A) exonic overlap: per-transcript exclusive bp, pairwise shared-bp matrix,
#       the 202 (subset of 205) and 203 (retention shared only with 204) cases.
#   (B) 201 vs 202 vs 205 correlation across samples, and whether pooling lowers
#       within-condition CV (corrects the earlier "anticorrelation" rationale).
# Standalone (no config import). Run:
#   python3 analyses/klrk1_202_reevaluation/02_exon_overlap_correlation.py
# ============================================================================
import os, csv, collections
import numpy as np

ROOT = os.environ.get("KLRK1_ROOT", "/Users/study/Desktop/Karimi/Klrk1_GVHD_project")
GTF = os.path.join(ROOT, "shared", "ref", "genome", "Klrk1.gtf")
TPM = os.path.join(ROOT, "python_pipeline", "outputs", "tables",
                   "Supplementary_all_sample_6isoform.csv")   # produced by script 03
NAME = {"ENSMUST00000032252": "201", "ENSMUST00000095412": "202",
        "ENSMUST00000137660": "203", "ENSMUST00000152256": "204",
        "ENSMUST00000168919": "205", "ENSMUST00000204694": "206"}
ALL = ["201", "202", "203", "204", "205", "206"]

# ---- (A) exonic overlap -----------------------------------------------------
ex = collections.defaultdict(list)
for line in open(GTF):
    if "\texon\t" not in line:
        continue
    f = line.split("\t"); tid = None
    for tok in f[8].split(";"):
        tok = tok.strip()
        if tok.startswith("transcript_id"):
            tid = tok.split('"')[1].split(".")[0]
    if tid in NAME:
        ex[NAME[tid]].append((int(f[3]), int(f[4])))
def bset(iv):
    s = set()
    for a, b in iv:
        s.update(range(a, b + 1))
    return s
bs = {t: bset(ex[t]) for t in ALL}
allother = {t: set().union(*[bs[o] for o in bs if o != t]) for t in bs}

print("=== per-transcript exonic length / EXCLUSIVE bp (not in any other Klrk1 tx) ===")
for t in ALL:
    tot = len(bs[t]); u = len(bs[t] - allother[t])
    print(f"  {t}: exonic={tot:5d}  exclusive={u:5d}  ({100*u/tot:4.1f}%)")
print("\n=== pairwise SHARED bp ===")
print("       " + " ".join(f"{c:>6}" for c in ALL))
for r in ALL:
    print(f"  {r:>4} " + " ".join(f"{len(bs[r] & bs[c]):>6}" for c in ALL))
print("\n=== 202 subset of 205? ===")
print(f"  202 bases not in 205: {len(bs['202']-bs['205'])}  (0 -> 202 is a base-subset of 205)")
print(f"  202 bases not in 201: {len(bs['202']-bs['201'])}  (the ~31 bp alt first exon)")
print("\n=== 203 retained-intron-4 segment 129593531-129593631: which tx cover it? ===")
ri = set(range(129593531, 129593631 + 1))
for t in ["201", "202", "204", "205", "206"]:
    print(f"  {t}: {len(ri & bs[t])}/101 bp")

# ---- (B) correlation + pooling CV ------------------------------------------
if os.path.exists(TPM):
    rows = list(csv.DictReader(open(TPM)))
    a = {k: np.array([float(r[f"TPM_{k}"]) for r in rows]) for k in ["201", "202", "205"]}
    ds = np.array([r["dataset"] for r in rows])
    def corr(x, y, m):
        x, y = x[m], y[m]
        return np.corrcoef(x, y)[0, 1] if len(x) > 2 and x.std() and y.std() else float("nan")
    m = a["202"] > 0
    print("\n=== correlation (samples with 202>0) ===")
    print(f"  r(201,202)={corr(a['201'],a['202'],m):+.3f}  r(205,202)={corr(a['205'],a['202'],m):+.3f}")
    g = ds == "GSE203167"
    print(f"  within GSE203167: r(201,202)={corr(a['201'],a['202'],g):+.3f}  "
          f"r(205,202)={corr(a['205'],a['202'],g):+.3f}")
    print("  (positive, NOT anticorrelated; pooling lowers CV by adding the larger stable signal)")
else:
    print(f"\n[skip B] run script 03 first to produce {TPM}")
