#!/usr/bin/env python3
# ============================================================================
# 08_klrk1_202_reevaluation.py
#
# PURPOSE
#   Re-evaluate Klrk1-202 (ENSMUST00000095412, protein_coding), which ver16
#   excluded globally as a "bias artefact" (its UNCORRECTED TPM inflated 0->52
#   in the discovery data). This script tests that decision against the
#   BIAS-CORRECTED (--gcBias --seqBias) quant.sf files for every dataset and
#   rebuilds the master isoform table + RI decomposition WITH 202 included.
#
#   Key facts established outside this script (see decision note):
#     * All 4 primary datasets were quantified with --gcBias --seqBias (verified
#       from cmd_info.json). 202 is therefore evaluated with bias correction ON.
#     * 202 is protein_coding, and shares 6/8 exons with canonical 201; the two
#       differ only in a ~31 bp alternative first exon and a 46 bp 3' tail. So
#       the 201<->202 split rests on very little unique sequence and is
#       EM-unstable, while the COMBINED 201+202 protein-coding pool is robust.
#     * Because 202 is an expressed, annotated isoform, the defensible total-
#       Klrk1 denominator includes it (6 isoforms). Excluding it undercounts
#       total Klrk1 and inflates every fraction (RI%, 203%).
#
# OUTPUT (python_pipeline/outputs/tables/):
#   Table4b_master_isoform_202included.csv    per-condition, all 6 isoforms
#   Klrk1_decomposition_202included.csv       RI split + 202, per condition
#   Klrk1_201_202_pool_stability.csv          201 vs 202 vs (201+202) CV
#   Supplementary_all_sample_6isoform.csv     per-sample raw TPM, all 6
#
# RUN: python3 python_pipeline/scripts/08_klrk1_202_reevaluation.py
# ============================================================================
import os, sys, csv
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config as cfg

# ---- six-isoform map (config.KLRK1 deliberately omits 202; add it back here)-
TX = {
    "ENSMUST00000032252": "201", "ENSMUST00000095412": "202",
    "ENSMUST00000137660": "203", "ENSMUST00000152256": "204",
    "ENSMUST00000168919": "205", "ENSMUST00000204694": "206",
}
BIOTYPE = {"201": "protein_coding", "202": "protein_coding", "205": "protein_coding",
           "203": "retained_intron", "204": "retained_intron", "206": "retained_intron"}
RI   = ["203", "204", "206"]
SIX  = ["201", "202", "203", "204", "205", "206"]   # correct total denominator
FIVE = ["201", "203", "204", "205", "206"]          # ver16 denominator (no 202)

META = {"GSE147371": "samples_primary_GSE147371_GSE109125.csv",
        "GSE109125": "samples_primary_GSE147371_GSE109125.csv",
        "GSE203167": "samples_GSE203167.csv",
        "GSE119943": "samples_GSE119943.csv"}

# dataset + condition display order (naive -> effector/memory)
DS_ORDER = ["GSE147371", "GSE109125", "GSE203167", "GSE119943"]
COND_ORDER = {
    "GSE147371": ["GVHD_CD4_Tn", "GVHD_CD4_Tem"],
    "GSE109125": ["Healthy_CD4_Naive", "Healthy_CD8_Naive_6wk", "Healthy_CD8_Naive_7wk",
                  "CD8_Naive_alt", "CD8_MPEC", "CD8_Tcm", "CD8_Effector_SLEC"],
    "GSE203167": ["WT_Pre", "WT_Post7", "TCF7cKO_Pre", "TCF7cKO_Post7"],
    "GSE119943": ["Arm_D4.5_EarlyEffector", "Arm_D7_MPEC", "Arm_D7_SLEC",
                  "Cl13_D7_Progenitor", "Cl13_D7_TermExhausted",
                  "pMIG_Cl13_Progenitor", "pMIG_Cl13_TermExhausted"],
}

def read_tpm(ds, salmon_dir):
    q = cfg.salmon_quant(salmon_dir.split("/")[-1], ds)
    if not os.path.exists(q):
        return None
    d = {v: 0.0 for v in TX.values()}
    with open(q) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            t = r["Name"].split(".")[0]
            if t in TX:
                d[TX[t]] = float(r["TPM"])
    return d

def load_conditions(ds):
    m = cfg.metadata_csv(META[ds]); conds = {}
    with open(m) as f:
        for r in csv.DictReader(f):
            if r.get("dataset", ds) != ds:
                continue
            conds.setdefault(r["condition"], []).append(r["salmon_dir"].split("/")[-1])
    return conds

def frac(d, iso_set, denom):
    tot = sum(d[i] for i in denom)
    return 100 * sum(d[i] for i in iso_set) / tot if tot else 0.0

# ---- gather per-sample rows -------------------------------------------------
per_sample = []   # dataset, condition, sample, 201..206, total6
cond_rows  = {}   # (ds,cond) -> list of per-sample TPM dicts
for ds in DS_ORDER:
    for cond, sds in load_conditions(ds).items():
        rows = [(sd, read_tpm(ds, sd)) for sd in sds]
        rows = [(sd, d) for sd, d in rows if d]
        if not rows:
            continue
        cond_rows[(ds, cond)] = [d for _, d in rows]
        for sd, d in rows:
            per_sample.append([ds, cond, sd] + [round(d[i], 4) for i in SIX]
                              + [round(sum(d[i] for i in SIX), 4)])

# ---- Table 4b: master isoform table, 202 included ---------------------------
master = [["dataset", "condition", "n",
           "TPM_201", "TPM_202", "TPM_203", "TPM_204", "TPM_205", "TPM_206",
           "total6_TPM", "total5_TPM",
           "pct201_6", "pct202_6", "pct203_6", "pct204_6", "pct205_6", "pct206_6",
           "RI_pct_6iso", "Klrk1_203_pct_6iso",
           "RI_pct_5iso_ver16", "Klrk1_203_pct_5iso_ver16", "delta_RI", "delta_203"]]
decomp = [["dataset", "condition", "n", "pct202_of6", "pct203_of6", "pct204_of6",
           "pct206_of6", "RI_pct_6iso", "share_203_of_RI", "pct203_of6_ver16", "RI_pct_5iso_ver16"]]
stab   = [["dataset", "condition", "n", "mean_201", "mean_202", "mean_201plus202",
           "CV_201", "CV_202", "CV_201plus202"]]

def cv(xs):
    xs = np.array(xs, float)
    m = xs.mean()
    return round(100 * xs.std(ddof=1) / m, 1) if len(xs) > 1 and m > 0 else np.nan

for ds in DS_ORDER:
    for cond in COND_ORDER[ds]:
        if (ds, cond) not in cond_rows:
            continue
        ds_rows = cond_rows[(ds, cond)]
        n = len(ds_rows)
        mean_tpm = {i: np.mean([d[i] for d in ds_rows]) for i in SIX}
        tot6 = sum(mean_tpm[i] for i in SIX)
        tot5 = sum(mean_tpm[i] for i in FIVE)
        # mean of per-sample fractions (matches manuscript convention)
        def mpct(iso_set, denom): return np.mean([frac(d, iso_set, denom) for d in ds_rows])
        ri6, ri5 = mpct(RI, SIX), mpct(RI, FIVE)
        p203_6, p203_5 = mpct(["203"], SIX), mpct(["203"], FIVE)
        p202 = mpct(["202"], SIX)
        master.append([ds, cond, n,
            *[round(mean_tpm[i], 3) for i in SIX], round(tot6, 3), round(tot5, 3),
            *[round(mpct([i], SIX), 2) for i in SIX],
            round(ri6, 2), round(p203_6, 2), round(ri5, 2), round(p203_5, 2),
            round(ri6 - ri5, 2), round(p203_6 - p203_5, 2)])
        # decomposition (fractions of total6) + within-RI 203 share
        p204_6, p206_6 = mpct(["204"], SIX), mpct(["206"], SIX)
        share = round(100 * p203_6 / ri6) if ri6 > 0 else 0
        decomp.append([ds, cond, n, round(p202, 2), round(p203_6, 2), round(p204_6, 2),
                       round(p206_6, 2), round(ri6, 2), share,
                       round(p203_5, 2), round(ri5, 2)])
        # 201/202 pooling stability
        v201 = [d["201"] for d in ds_rows]; v202 = [d["202"] for d in ds_rows]
        vsum = [a + b for a, b in zip(v201, v202)]
        stab.append([ds, cond, n, round(np.mean(v201), 2), round(np.mean(v202), 2),
                     round(np.mean(vsum), 2), cv(v201), cv(v202), cv(vsum)])

# ---- write ------------------------------------------------------------------
TBL = cfg.TABLE_DIR
def write(name, rows):
    p = os.path.join(TBL, name)
    with open(p, "w", newline="") as f:
        csv.writer(f).writerows(rows)
    return p

write("Table4b_master_isoform_202included.csv", master)
write("Klrk1_decomposition_202included.csv", decomp)
write("Klrk1_201_202_pool_stability.csv", stab)
write("Supplementary_all_sample_6isoform.csv",
      [["dataset", "condition", "sample", "TPM_201", "TPM_202", "TPM_203",
        "TPM_204", "TPM_205", "TPM_206", "total6_TPM"]] + per_sample)

# ---- UNIT-LEVEL tables (config.REPORT_UNITS: 201 | 205/202 | 203 | 204 | 206) ----
# Name reported quantities by the isoform(s) the reads actually support (see
# scripts/13_klrk1_identifiability_map.py). 202 has no read that separates it from
# 205, so it is reported inside the Klrk1-205/202 unit; total Klrk1 = sum of 5 units
# (= all 6 transcripts). These feed the manuscript Table 4/5 + decomposition.
UNITS = list(cfg.REPORT_UNITS)                    # ordered unit labels
umaster = [["dataset", "condition", "n"] + [f"TPM_{u}" for u in UNITS] + ["total_TPM"]
           + [f"pct_{u}" for u in UNITS] + ["RI_pct", "Klrk1_203_pct"]]
udecomp = [["dataset", "condition", "n"] + [f"pct_{u}" for u in UNITS]
           + ["RI_pct", "share_203_of_RI"]]
for ds in DS_ORDER:
    for cond in COND_ORDER[ds]:
        if (ds, cond) not in cond_rows:
            continue
        ds_rows = cond_rows[(ds, cond)]
        n = len(ds_rows)
        urows = [cfg.to_units(d) for d in ds_rows]         # per-sample unit TPM
        mean_u = {u: np.mean([r[u] for r in urows]) for u in UNITS}
        tot = sum(mean_u.values())
        def upct(u):   # mean of per-sample unit fractions (manuscript convention)
            return np.mean([100 * r[u] / sum(r.values()) if sum(r.values()) else 0 for r in urows])
        pcts = {u: upct(u) for u in UNITS}
        ri = sum(pcts[u] for u in cfg.RI_UNITS)
        share = round(100 * pcts["Klrk1-203"] / ri) if ri > 0 else 0
        umaster.append([ds, cond, n] + [round(mean_u[u], 3) for u in UNITS] + [round(tot, 3)]
                       + [round(pcts[u], 2) for u in UNITS] + [round(ri, 2), round(pcts["Klrk1-203"], 2)])
        udecomp.append([ds, cond, n] + [round(pcts[u], 2) for u in UNITS] + [round(ri, 2), share])
write("Table4_units_master.csv", umaster)
write("Klrk1_decomposition_units.csv", udecomp)

# ---- console summary --------------------------------------------------------
print("=" * 104)
print("Klrk1-202 RE-EVALUATION  (bias-corrected quant.sf; 202 = protein_coding, added to denominator)")
print("=" * 104)
print(f"{'dataset':10s} {'condition':24s} {'n':>2} | {'202%':>5} | "
      f"{'203% v16':>8} {'203% new':>8} | {'RI% v16':>7} {'RI% new':>7} | flag")
print("-" * 104)
for row in decomp[1:]:
    ds, cond, n, p202, p203_6, p204_6, p206_6, ri6, share, p203_5, ri5 = row
    flag = "202 MATERIAL (>=10%)" if p202 >= 10 else ("202 minor" if p202 >= 2 else "202~0")
    print(f"{ds:10s} {cond:24s} {n:>2} | {p202:5.1f} | {p203_5:8.2f} {p203_6:8.2f} | "
          f"{ri5:7.2f} {ri6:7.2f} | {flag}")
print("-" * 104)
print(f"Tables written to {TBL}")
print("  - Table4b_master_isoform_202included.csv")
print("  - Klrk1_decomposition_202included.csv")
print("  - Klrk1_201_202_pool_stability.csv")
print("  - Supplementary_all_sample_6isoform.csv")
