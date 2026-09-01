#!/usr/bin/env python3
# ============================================================================
# 12_manuscript_fig2_fig3.py
# Regenerate manuscript Figure 2 and Figure 3 under the 5-unit convention,
# each from its correct dataset:
#   Figure 2 = ImmGen (GSE109125) antigen-experienced CD8+ subsets — RI% bar
#              (RI = Klrk1-203 + Klrk1-204 + Klrk1-206 over total Klrk1).
#   Figure 3 = allogeneic GVHD CD4+ (GSE147371) — stacked isoform proportions
#              as the five reporting units (Tn, Tem).
# Reads the same primary metadata + bias-corrected quant.sf as script 02; uses
# config.KLRK1 (six transcripts) and config.REPORT_UNITS (five units).
# Run: python3 python_pipeline/scripts/12_manuscript_fig2_fig3.py
# ============================================================================
import os, sys
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config as cfg

meta = pd.read_csv(cfg.metadata_csv("samples_primary_GSE147371_GSE109125.csv"))
rows = []
for _, s in meta.iterrows():
    q = cfg.salmon_quant(s["salmon_dir"], s["dataset"])
    if not os.path.exists(q):
        print(f"[warn] missing {q}"); continue
    df = pd.read_csv(q, sep="\t"); df["tx"] = df["Name"].str.replace(r"\.\d+$", "", regex=True)
    tpm = dict(zip(df["tx"], df["TPM"]))
    for tx in cfg.KLRK1:
        rows.append({"sample_id": s["sample_id"], "condition": s["condition"],
                     "dataset": s["dataset"], "transcript_id": tx, "isoform": cfg.KLRK1[tx],
                     "biotype": cfg.BIOTYPE[tx], "TPM": tpm.get(tx, 0.0)})
dat = pd.DataFrame(rows)
tot = dat.groupby("sample_id")["TPM"].sum().rename("total_klrk1")
dat = dat.merge(tot, on="sample_id")
dat["proportion"] = np.where(dat.total_klrk1 > 0, dat.TPM / dat.total_klrk1 * 100, 0)
iso2unit = {f"Klrk1-{s}": u for u, ss in cfg.REPORT_UNITS.items() for s in ss}
dat["unit"] = dat["isoform"].map(iso2unit)

# ---- Figure 2 — ImmGen (GSE109125) antigen-experienced CD8+ RI% ----
ri = (dat[dat.biotype == "retained_intron"].groupby(["condition", "sample_id"])["TPM"].sum()
      .reset_index().merge(tot.reset_index(), on="sample_id"))
ri["ri_pct"] = ri.TPM / ri.total_klrk1 * 100
F2 = ["CD8_MPEC", "CD8_Tcm", "CD8_Effector_SLEC"]           # antigen-experienced subsets
lab2 = ["MPEC", "Tcm", "SLEC"]
m2 = ri[ri.condition.isin(F2)].groupby("condition")["ri_pct"].agg(["mean", "std", "count"]).reindex(F2)
fig, ax = plt.subplots(figsize=(5.2, 4.4))
xb = np.arange(len(F2))
ax.bar(xb, m2["mean"], yerr=[0 if (c < 2 or np.isnan(sd)) else sd for sd, c in zip(m2["std"], m2["count"])],
       capsize=4, color="#7d3c98")
for i, v in enumerate(m2["mean"]):
    ax.text(i, v + 0.6, f"{v:.1f}%", ha="center", fontsize=9)
ax.set_xticks(xb); ax.set_xticklabels(lab2)
ax.set_ylabel("Retained-intron fraction of total Klrk1 (%)")
ax.set_ylim(0, max(m2["mean"]) * 1.25)
fig.tight_layout()
fig.savefig(os.path.join(cfg.FIG_DIR, "Figure2_ImmGen_RIpct.pdf"))
fig.savefig(os.path.join(cfg.FIG_DIR, "Figure2_ImmGen_RIpct.png"), dpi=200); plt.close(fig)
print("Figure 2 (ImmGen RI%):", {c: round(v, 1) for c, v in zip(lab2, m2["mean"])})

# ---- Figure 3 — GVHD CD4+ (GSE147371) isoform proportions, 5 units ----
F3 = ["GVHD_CD4_Tn", "GVHD_CD4_Tem"]; lab3 = ["Tn", "Tem"]
gv = dat[dat.condition.isin(F3)]
usamp = gv.groupby(["sample_id", "condition", "unit"])["proportion"].sum().reset_index()
prop = usamp.groupby(["condition", "unit"])["proportion"].mean().unstack().reindex(F3)[list(cfg.REPORT_UNITS)]
fig, ax = plt.subplots(figsize=(6.2, 4.6)); xb = np.arange(len(F3)); bottom = np.zeros(len(F3))
for u in cfg.REPORT_UNITS:
    ax.bar(xb, prop[u].values, bottom=bottom, label=u, width=0.6); bottom += prop[u].values
ax.set_xticks(xb); ax.set_xticklabels(lab3)
ax.set_ylabel("Isoform proportion (%)"); ax.set_ylim(0, 100)
ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8, frameon=False, title="reporting unit")
fig.subplots_adjust(right=0.72); fig.tight_layout()
fig.savefig(os.path.join(cfg.FIG_DIR, "Figure3_GVHD_proportions.pdf"))
fig.savefig(os.path.join(cfg.FIG_DIR, "Figure3_GVHD_proportions.png"), dpi=200); plt.close(fig)
tem_ri = ri[ri.condition == "GVHD_CD4_Tem"]["ri_pct"].mean()
tem_203 = dat[(dat.condition == "GVHD_CD4_Tem") & (dat.transcript_id == cfg.ISO_203)]["proportion"].mean()
print(f"Figure 3 (GVHD): Tem RI% = {tem_ri:.1f}, Tem Klrk1-203 = {tem_203:.1f}%")
print(f"figures -> {cfg.FIG_DIR}")
