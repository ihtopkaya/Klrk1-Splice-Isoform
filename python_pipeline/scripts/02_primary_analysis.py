#!/usr/bin/env python3
# ============================================================================
# 02_primary_analysis.py  —  Python port of 02_klrk1_analysis_primary.R
#
# Primary discovery analysis (GSE147371 GVHD CD4+ + GSE109125 ImmGen CD8):
# reads Salmon quant.sf for 19 samples, computes per-condition Klrk1 isoform
# usage, and writes Fig1-4 + the master/supplementary tables. Same numbers as
# the R version; figures are matplotlib (styled differently, same content).
#
# Run:  python3 scripts/02_primary_analysis.py
# ============================================================================
import os, sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config as cfg

# condition plotting order (naive -> effector/memory; CD4 GVHD then CD8 cascade)
COND_ORDER = ["Healthy_CD4_Naive", "GVHD_CD4_Tn", "GVHD_CD4_Tem",
              "Healthy_CD8_Naive_6wk", "Healthy_CD8_Naive_7wk", "CD8_Naive_alt",
              "CD8_MPEC", "CD8_Tcm", "CD8_Effector_SLEC"]

# ---- 1. load metadata + read Salmon -----------------------------------------
meta = pd.read_csv(cfg.metadata_csv("samples_primary_GSE147371_GSE109125.csv"))

def read_klrk1(salmon_dir, dataset):
    q = cfg.salmon_quant(salmon_dir, dataset)
    if not os.path.exists(q):
        print(f"[warn] missing {q}"); return None
    df = pd.read_csv(q, sep="\t")
    df["tx"] = df["Name"].str.replace(r"\.\d+$", "", regex=True)
    df = df[df["tx"].isin(cfg.KLRK1)]
    return dict(zip(df["tx"], df["TPM"]))

rows = []
for _, s in meta.iterrows():
    tpm = read_klrk1(s["salmon_dir"], s["dataset"])
    if tpm is None: continue
    for tx in cfg.KLRK1:
        rows.append({"sample_id": s["sample_id"], "condition": s["condition"],
                     "cell_type": s["cell_type"], "dataset": s["dataset"],
                     "transcript_id": tx, "isoform": cfg.KLRK1[tx],
                     "biotype": cfg.BIOTYPE[tx], "TPM": tpm.get(tx, 0.0)})
dat = pd.DataFrame(rows)
print(f"Loaded {len(dat)} isoform measurements from {dat.sample_id.nunique()} samples")

# ---- 2. per-sample totals + proportions -------------------------------------
tot = dat.groupby("sample_id")["TPM"].sum().rename("total_klrk1")
dat = dat.merge(tot, on="sample_id")
dat["proportion"] = np.where(dat.total_klrk1 > 0, dat.TPM / dat.total_klrk1 * 100, 0)

ri = (dat[dat.biotype == "retained_intron"].groupby(["sample_id", "condition"])["TPM"].sum()
      .rename("ri_tpm").reset_index().merge(tot.reset_index(), on="sample_id"))
ri["ri_pct"] = np.where(ri.total_klrk1 > 0, ri.ri_tpm / ri.total_klrk1 * 100, 0)

# ---- 3. condition-level summaries -------------------------------------------
def msd(x):  # mean, sd(NA if n<2)
    return pd.Series({"mean": x.mean(), "sd": x.std(ddof=1) if len(x) > 1 else np.nan, "n": len(x)})

cond_total = (dat.groupby(["condition", "sample_id"])["TPM"].sum().groupby("condition")
              .apply(msd).unstack())
cond_ri   = ri.groupby("condition")["ri_pct"].apply(msd).unstack()
d203 = dat[dat.transcript_id == cfg.ISO_203]
cond_203 = d203.groupby("condition").agg(mean_tpm=("TPM", "mean"),
                                         mean_pct=("proportion", "mean"))
cond_iso = dat.groupby(["condition", "isoform"])["TPM"].mean().unstack()

# ---- 4. console report ------------------------------------------------------
print("\n" + "=" * 68 + "\n  KLRK1 ISOFORM ANALYSIS — PRIMARY (Python)\n" + "=" * 68)
for c in COND_ORDER:
    if c not in cond_total.index: continue
    n = int(cond_total.loc[c, "n"]); tt = cond_total.loc[c, "mean"]
    print(f"\n--- {c} (n={n}) ---")
    print(f"  Total Klrk1: {tt:6.1f} TPM   Klrk1-203: {cond_203.loc[c,'mean_tpm']:6.2f} TPM "
          f"({cond_203.loc[c,'mean_pct']:.1f}%)   Total RI: {cond_ri.loc[c,'mean']:.1f}%")

gt = lambda c: cond_total.loc[c, "mean"]
print("\n--- Fold changes ---")
print(f"  GVHD Tn -> Tem: {gt('GVHD_CD4_Tem')/gt('GVHD_CD4_Tn'):.0f}x")
print(f"  CD8 6wk Naive -> Effector: {gt('CD8_Effector_SLEC')/gt('Healthy_CD8_Naive_6wk'):.0f}x")

# ---- 5. figures -------------------------------------------------------------
present = [c for c in COND_ORDER if c in cond_total.index]
x = np.arange(len(present))

def barfig(vals, ylab, title, fname, err=None, color="#3b6ea5"):
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(x, vals, yerr=err, color=color, capsize=3)
    ax.set_xticks(x); ax.set_xticklabels(present, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(ylab); ax.set_title(title, fontweight="bold")
    fig.tight_layout(); fig.savefig(os.path.join(cfg.FIG_DIR, fname + ".pdf"))
    fig.savefig(os.path.join(cfg.FIG_DIR, fname + ".png"), dpi=200); plt.close(fig)

# Fig1: Klrk1-203 TPM by condition
barfig([cond_203.loc[c, "mean_tpm"] for c in present], "Klrk1-203 (TPM)",
       "Fig1 — Klrk1-203 (ENSMUST00000137660) induction", "Fig1_Klrk1-203", color="#c0392b")
# Fig4: total Klrk1 by condition
barfig([cond_total.loc[c, "mean"] for c in present], "Total Klrk1 (TPM)",
       "Fig4 — Total Klrk1 expression", "Fig4_total_Klrk1",
       err=[0 if np.isnan(cond_total.loc[c, "sd"]) else cond_total.loc[c, "sd"] for c in present])
# Fig2: CD8 differentiation kinetics (total + 203) for CD8 conditions only
cd8 = [c for c in present if c.startswith("CD8") or "CD8" in c]
if cd8:
    xi = np.arange(len(cd8))
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(xi, [cond_total.loc[c, "mean"] for c in cd8], "o-", label="Total Klrk1")
    ax.plot(xi, [cond_203.loc[c, "mean_tpm"] for c in cd8], "s--", color="#c0392b", label="Klrk1-203")
    ax.set_xticks(xi); ax.set_xticklabels(cd8, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("TPM"); ax.set_title("Fig2 — CD8 differentiation kinetics", fontweight="bold")
    ax.legend(); fig.tight_layout()
    fig.savefig(os.path.join(cfg.FIG_DIR, "Fig2_CD8_kinetics.pdf"))
    fig.savefig(os.path.join(cfg.FIG_DIR, "Fig2_CD8_kinetics.png"), dpi=200); plt.close(fig)
# Fig3: stacked isoform proportions by condition — five read-identifiability units
iso2unit = {f"Klrk1-{s}": u for u, ss in cfg.REPORT_UNITS.items() for s in ss}
dat["unit"] = dat["isoform"].map(iso2unit)
usamp = dat.groupby(["sample_id", "condition", "unit"])["proportion"].sum().reset_index()
prop = (usamp.groupby(["condition", "unit"])["proportion"].mean().unstack()
        .reindex(present)[list(cfg.REPORT_UNITS)])
fig, ax = plt.subplots(figsize=(10, 5)); bottom = np.zeros(len(present))
for u in cfg.REPORT_UNITS:
    ax.bar(x, prop[u].values, bottom=bottom, label=u); bottom += prop[u].values
ax.set_xticks(x); ax.set_xticklabels(present, rotation=45, ha="right", fontsize=8)
ax.set_ylabel("Isoform proportion (%)")   # no in-figure title (added in manuscript text)
ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8, frameon=False,
          title="reporting unit")
fig.subplots_adjust(right=0.80); fig.tight_layout()
fig.savefig(os.path.join(cfg.FIG_DIR, "Fig3_isoform_proportions.pdf"))
fig.savefig(os.path.join(cfg.FIG_DIR, "Fig3_isoform_proportions.png"), dpi=200); plt.close(fig)

# ---- 6. tables --------------------------------------------------------------
table4 = cond_iso.copy()
table4["total_klrk1"] = cond_total["mean"]
table4["RI_pct"] = cond_ri["mean"]
table4["Klrk1_203_pct"] = cond_203["mean_pct"]
table4 = table4.reindex(present)
table4.to_csv(os.path.join(cfg.TABLE_DIR, "Table4_master_isoform_comparison.csv"))
dat.to_csv(os.path.join(cfg.TABLE_DIR, "Supplementary_all_sample_isoform_data.csv"), index=False)

print(f"\nFigures -> {cfg.FIG_DIR}\nTables  -> {cfg.TABLE_DIR}")
