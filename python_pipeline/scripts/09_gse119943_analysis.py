#!/usr/bin/env python3
# ============================================================================
# 09_gse119943_analysis.py — Python port of 09_klrk1_analysis_GSE119943.R
# Validation dataset: Yao 2019 LCMV CD8 (19 bulk samples, 50 bp single-end).
# Salmon quant.sf for all 19 samples are present under datasets/GSE119943_Yao_LCMV/data/salmon/.
# (A graceful [skip] guard remains for portability if the data is absent.)
# Run: python3 python_pipeline/scripts/09_gse119943_analysis.py
# ============================================================================
import os, sys
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config as cfg

DATASET = "GSE119943"
meta = pd.read_csv(cfg.metadata_csv("samples_GSE119943.csv"))
rows = []
for _, s in meta.iterrows():
    q = cfg.salmon_quant(s["salmon_dir"], DATASET)
    if not os.path.exists(q):
        continue
    df = pd.read_csv(q, sep="\t"); df["tx"] = df["Name"].str.replace(r"\.\d+$", "", regex=True)
    tpm = dict(zip(df["tx"], df["TPM"]))
    for tx in cfg.KLRK1:
        rows.append({"sample_id": s["sample_id"], "condition": s["condition"],
                     "transcript_id": tx, "isoform": cfg.KLRK1[tx],
                     "biotype": cfg.BIOTYPE[tx], "TPM": tpm.get(tx, 0.0)})

if not rows:
    print("[skip] No GSE119943 Salmon data found under "
          f"{cfg.salmon_quant('<sample>', DATASET)}\n"
          "       Re-fetch the data to run this analysis (see dataset README).")
    sys.exit(0)

dat = pd.DataFrame(rows)
tot = dat.groupby("sample_id")["TPM"].sum().rename("total_klrk1")
dat = dat.merge(tot, on="sample_id")
dat["proportion"] = np.where(dat.total_klrk1 > 0, dat.TPM / dat.total_klrk1 * 100, 0)
d203 = dat[dat.transcript_id == cfg.ISO_203].groupby("condition").agg(tpm=("TPM","mean"), pct=("proportion","mean"))
ctot = dat.groupby(["condition","sample_id"])["TPM"].sum().groupby("condition").mean()
order = list(d203.index)
print(f"{DATASET}: per-condition Klrk1-203")
for c in order:
    print(f"  {c:26s} total={ctot[c]:7.1f}  Klrk1-203={d203.loc[c,'tpm']:6.2f} ({d203.loc[c,'pct']:.1f}%)")

x = np.arange(len(order))
fig, ax = plt.subplots(figsize=(10, 5))
ax.bar(x, [d203.loc[c,"tpm"] for c in order], color="#c0392b")
ax.set_xticks(x); ax.set_xticklabels(order, rotation=45, ha="right", fontsize=7)
ax.set_ylabel("Klrk1-203 (TPM)"); ax.set_title("FigV7 — Klrk1-203 (GSE119943)", fontweight="bold")
fig.tight_layout(); fig.savefig(os.path.join(cfg.FIG_DIR, "FigV7_Klrk1-203_GSE119943.pdf"))
fig.savefig(os.path.join(cfg.FIG_DIR, "FigV7_Klrk1-203_GSE119943.png"), dpi=200); plt.close(fig)
dat.to_csv(os.path.join(cfg.TABLE_DIR, "Supplementary_GSE119943_all_samples.csv"), index=False)

# ---- Combined Figure 8 (A total, B Klrk1-203, C RI%, D isoform units); mean ± SD ----
ORDER = ["Arm_D4.5_EarlyEffector", "Arm_D7_SLEC", "Arm_D7_MPEC", "Cl13_D7_Progenitor",
         "Cl13_D7_TermExhausted", "pMIG_Cl13_Progenitor", "pMIG_Cl13_TermExhausted"]
present = [c for c in ORDER if c in ctot.index]; xi = np.arange(len(present))
tot_by = dat.groupby(["condition", "sample_id"])["TPM"].sum()
ctot_sd = tot_by.groupby("condition").std(ddof=1)
d203_sd = dat[dat.transcript_id == cfg.ISO_203].groupby("condition")["TPM"].std(ddof=1)
riS = (dat[dat.biotype == "retained_intron"].groupby(["condition", "sample_id"])["TPM"].sum()
       .reset_index().merge(tot.reset_index(), on="sample_id"))
riS["ri_pct"] = riS.TPM / riS.total_klrk1 * 100
cri = riS.groupby("condition")["ri_pct"].mean(); cri_sd = riS.groupby("condition")["ri_pct"].std(ddof=1)
iso2unit = {f"Klrk1-{s}": u for u, ss in cfg.REPORT_UNITS.items() for s in ss}
dat["unit"] = dat["isoform"].map(iso2unit)
usamp = dat.groupby(["sample_id", "condition", "unit"])["proportion"].sum().reset_index()
prop = usamp.groupby(["condition", "unit"])["proportion"].mean().unstack().reindex(present)[list(cfg.REPORT_UNITS)]
def err(s): return [0 if np.isnan(s.get(c, np.nan)) else s[c] for c in present]
lbl = [c.replace("_", " ").replace("Arm ", "Arm ").replace("Cl13 D7 ", "Cl13 ") for c in present]

fig, axs = plt.subplots(2, 2, figsize=(12, 8.5))
def panel(ax, vals, sd, ylab, color, tag):
    ax.bar(xi, vals, yerr=sd, capsize=3, color=color)
    ax.set_xticks(xi); ax.set_xticklabels(lbl, rotation=30, ha="right", fontsize=7.5)
    ax.set_ylabel(ylab, fontsize=9)
    ax.text(-0.11, 1.04, tag, transform=ax.transAxes, fontsize=14, fontweight="bold")
panel(axs[0, 0], [ctot[c] for c in present], err(ctot_sd), "Total Klrk1 (TPM)", "#3b6ea5", "A")
panel(axs[0, 1], [d203.loc[c, "tpm"] for c in present], err(d203_sd), "Klrk1-203 (TPM)", "#c0392b", "B")
panel(axs[1, 0], [cri[c] for c in present], err(cri_sd), "Retained-intron fraction (%)", "#7d3c98", "C")
axD = axs[1, 1]; bottom = np.zeros(len(present))
for u in cfg.REPORT_UNITS:
    axD.bar(xi, prop[u].values, bottom=bottom, label=u); bottom += prop[u].values
axD.set_xticks(xi); axD.set_xticklabels(lbl, rotation=30, ha="right", fontsize=7.5)
axD.set_ylabel("Isoform proportion (%)", fontsize=9)
axD.text(-0.11, 1.04, "D", transform=axD.transAxes, fontsize=14, fontweight="bold")
axD.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=7.5, frameon=False, title="unit")
fig.tight_layout()
fig.savefig(os.path.join(cfg.FIG_DIR, "Figure8_GSE119943_combined.pdf"))
fig.savefig(os.path.join(cfg.FIG_DIR, "Figure8_GSE119943_combined.png"), dpi=200); plt.close(fig)
print(f"figure -> {cfg.FIG_DIR}  (+ Figure8_GSE119943_combined)")
