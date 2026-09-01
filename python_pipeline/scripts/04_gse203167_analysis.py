#!/usr/bin/env python3
# ============================================================================
# 04_gse203167_analysis.py — Python port of 04_klrk1_analysis_GSE203167.R
# Validation dataset: WT vs TCF-7 cKO CD8, pre/post transplant (12 samples).
# Reproduces per-condition Klrk1 isoform usage + FigV1-4. Run: python3 scripts/04_*.py
# ============================================================================
import os, sys
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config as cfg

DATASET = "GSE203167"
COND_ORDER = ["WT_Pre", "WT_Post7", "TCF7cKO_Pre", "TCF7cKO_Post7"]

meta = pd.read_csv(cfg.metadata_csv("samples_GSE203167.csv"))
rows = []
for _, s in meta.iterrows():
    q = cfg.salmon_quant(s["salmon_dir"], DATASET)
    if not os.path.exists(q):
        print(f"[warn] missing {q}"); continue
    df = pd.read_csv(q, sep="\t"); df["tx"] = df["Name"].str.replace(r"\.\d+$", "", regex=True)
    tpm = dict(zip(df["tx"], df["TPM"]))
    for tx in cfg.KLRK1:
        rows.append({"sample_id": s["sample_id"], "condition": s["condition"],
                     "transcript_id": tx, "isoform": cfg.KLRK1[tx],
                     "biotype": cfg.BIOTYPE[tx], "TPM": tpm.get(tx, 0.0)})
dat = pd.DataFrame(rows)
tot = dat.groupby("sample_id")["TPM"].sum().rename("total_klrk1")
dat = dat.merge(tot, on="sample_id")
dat["proportion"] = np.where(dat.total_klrk1 > 0, dat.TPM / dat.total_klrk1 * 100, 0)

d203 = dat[dat.transcript_id == cfg.ISO_203].groupby("condition").agg(
    tpm=("TPM", "mean"), pct=("proportion", "mean"))
ctot = dat.groupby(["condition", "sample_id"])["TPM"].sum().groupby("condition").mean()
ri = (dat[dat.biotype == "retained_intron"].groupby(["condition", "sample_id"])["TPM"].sum()
      .reset_index().merge(tot.reset_index(), on="sample_id"))
ri["ri_pct"] = ri.TPM / ri.total_klrk1 * 100
cri = ri.groupby("condition")["ri_pct"].mean()

print(f"{DATASET}: per-condition summary")
for c in COND_ORDER:
    if c in ctot.index:
        print(f"  {c:15s} total={ctot[c]:7.1f}  Klrk1-203={d203.loc[c,'tpm']:6.2f} ({d203.loc[c,'pct']:.1f}%)  RI={cri[c]:.1f}%")

present = [c for c in COND_ORDER if c in ctot.index]; x = np.arange(len(present))
def bar(vals, ylab, title, fn, color="#3b6ea5"):
    fig, ax = plt.subplots(figsize=(7, 4.5)); ax.bar(x, vals, color=color)
    ax.set_xticks(x); ax.set_xticklabels(present, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel(ylab); ax.set_title(title, fontweight="bold"); fig.tight_layout()
    fig.savefig(os.path.join(cfg.FIG_DIR, fn + ".pdf")); fig.savefig(os.path.join(cfg.FIG_DIR, fn + ".png"), dpi=200); plt.close(fig)
bar([d203.loc[c, "tpm"] for c in present], "Klrk1-203 (TPM)", "FigV1 — Klrk1-203 (GSE203167)", "FigV1_Klrk1-203_GSE203167", "#c0392b")
bar([ctot[c] for c in present], "Total Klrk1 (TPM)", "FigV2 — Total Klrk1 (GSE203167)", "FigV2_total_Klrk1_GSE203167")
bar([cri[c] for c in present], "Retained-intron (%)", "FigV4 — RI fraction (GSE203167)", "FigV4_RI_pct_GSE203167", "#7d3c98")
iso2unit = {f"Klrk1-{s}": u for u, ss in cfg.REPORT_UNITS.items() for s in ss}
dat["unit"] = dat["isoform"].map(iso2unit)
usamp = dat.groupby(["sample_id", "condition", "unit"])["proportion"].sum().reset_index()
prop = usamp.groupby(["condition", "unit"])["proportion"].mean().unstack().reindex(present)[list(cfg.REPORT_UNITS)]
fig, ax = plt.subplots(figsize=(7.5, 4.5)); bottom = np.zeros(len(present))
for u in cfg.REPORT_UNITS:
    ax.bar(x, prop[u].values, bottom=bottom, label=u); bottom += prop[u].values
ax.set_xticks(x); ax.set_xticklabels(present, rotation=30, ha="right", fontsize=8)
ax.set_ylabel("Proportion (%)"); ax.set_title("FigV3 — Isoform proportions (GSE203167, 5 units)", fontweight="bold")
ax.legend(ncol=5, fontsize=7); fig.tight_layout()
fig.savefig(os.path.join(cfg.FIG_DIR, "FigV3_isoform_proportions_GSE203167.pdf"))
fig.savefig(os.path.join(cfg.FIG_DIR, "FigV3_isoform_proportions_GSE203167.png"), dpi=200); plt.close(fig)

# ---- Combined Figure 4 (A total, B Klrk1-203, C RI%, D isoform units); mean ± SD ----
def sd_by_cond(series_per_sample):
    return series_per_sample.groupby("condition").std(ddof=1)
tot_by = dat.groupby(["condition", "sample_id"])["TPM"].sum()
ctot_sd = tot_by.groupby("condition").std(ddof=1)
d203_sd = dat[dat.transcript_id == cfg.ISO_203].groupby("condition")["TPM"].std(ddof=1)
cri_sd = ri.groupby("condition")["ri_pct"].std(ddof=1)
def err(s): return [0 if np.isnan(s.get(c, np.nan)) else s[c] for c in present]

fig, axs = plt.subplots(2, 2, figsize=(11, 8))
lbl = [c.replace("_", " ") for c in present]
def panel(ax, vals, sd, ylab, color, tag):
    ax.bar(x, vals, yerr=sd, capsize=3, color=color)
    ax.set_xticks(x); ax.set_xticklabels(lbl, rotation=25, ha="right", fontsize=8)
    ax.set_ylabel(ylab, fontsize=9)
    ax.text(-0.12, 1.04, tag, transform=ax.transAxes, fontsize=14, fontweight="bold")
panel(axs[0, 0], [ctot[c] for c in present], err(ctot_sd), "Total Klrk1 (TPM)", "#3b6ea5", "A")
panel(axs[0, 1], [d203.loc[c, "tpm"] for c in present], err(d203_sd), "Klrk1-203 (TPM)", "#c0392b", "B")
panel(axs[1, 0], [cri[c] for c in present], err(cri_sd), "Retained-intron fraction (%)", "#7d3c98", "C")
# D: stacked 5-unit proportions
axD = axs[1, 1]; bottom = np.zeros(len(present))
for u in cfg.REPORT_UNITS:
    axD.bar(x, prop[u].values, bottom=bottom, label=u); bottom += prop[u].values
axD.set_xticks(x); axD.set_xticklabels(lbl, rotation=25, ha="right", fontsize=8)
axD.set_ylabel("Isoform proportion (%)", fontsize=9)
axD.text(-0.12, 1.04, "D", transform=axD.transAxes, fontsize=14, fontweight="bold")
axD.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=7.5, frameon=False, title="unit")
fig.tight_layout()
fig.savefig(os.path.join(cfg.FIG_DIR, "Figure4_GSE203167_combined.pdf"))
fig.savefig(os.path.join(cfg.FIG_DIR, "Figure4_GSE203167_combined.png"), dpi=200); plt.close(fig)
print(f"figures -> {cfg.FIG_DIR}  (+ Figure4_GSE203167_combined)")
