#!/usr/bin/env python3
# ============================================================================
# 05_housekeeping_control.py — Python port of 05_housekeeping_control.R
# Negative control: is the retained-intron (RI) proportion at Klrk1 far above the
# RI proportion at housekeeping genes Actb/Gapdh? (rules out gDNA contamination /
# pipeline artifact). Reads ALL Salmon outputs across datasets. Run: python3 scripts/05_*.py
# ============================================================================
import os, sys, glob
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config as cfg

# transcript -> (gene, is_retained_intron)   [from GRCm39 Ensembl cDNA headers]
ACTB_RI  = {"ENSMUST00000165629", "ENSMUST00000164765"}
ACTB_ALL = {"ENSMUST00000100497","ENSMUST00000167721","ENSMUST00000171419",
            "ENSMUST00000163829","ENSMUST00000106216","ENSMUST00000196997",
            "ENSMUST00000167386"} | ACTB_RI
GAPDH_RI  = {"ENSMUST00000192506","ENSMUST00000147954","ENSMUST00000182464","ENSMUST00000144588"}
GAPDH_ALL = {"ENSMUST00000117757","ENSMUST00000118875","ENSMUST00000073605",
             "ENSMUST00000183272","ENSMUST00000182052","ENSMUST00000182277",
             "ENSMUST00000144205","ENSMUST00000182115","ENSMUST00000182670"} | GAPDH_RI
KLRK1_ALL = set(cfg.KLRK1)
KLRK1_RI  = {t for t, b in cfg.BIOTYPE.items() if b == "retained_intron"}
GENES = {"Actb": (ACTB_ALL, ACTB_RI), "Gapdh": (GAPDH_ALL, GAPDH_RI),
         "Klrk1": (KLRK1_ALL, KLRK1_RI)}
ALL_TX = ACTB_ALL | GAPDH_ALL | KLRK1_ALL

quants = glob.glob(os.path.join(cfg.PROJECT_ROOT, "datasets", "*", "data", "salmon", "*", "quant.sf"))
print(f"Salmon quant.sf files found: {len(quants)}")

recs = []
for q in quants:
    sample = os.path.basename(os.path.dirname(q))
    df = pd.read_csv(q, sep="\t"); df["tx"] = df["Name"].str.replace(r"\.\d+$", "", regex=True)
    tpm = dict(zip(df["tx"], df["TPM"]))
    for gene, (allset, riset) in GENES.items():
        total = sum(tpm.get(t, 0) for t in allset)
        ri    = sum(tpm.get(t, 0) for t in riset)
        if total > 0:
            recs.append({"sample": sample, "gene": gene, "ri_pct": 100 * ri / total})
res = pd.DataFrame(recs)
summary = res.groupby("gene")["ri_pct"].agg(median="median", mean="mean", n="count").reindex(["Klrk1","Actb","Gapdh"])
print("\n=== Retained-intron proportion per gene (negative control) ===")
print(summary.round(2).to_string())
kl, ac, ga = summary.loc["Klrk1","median"], summary.loc["Actb","median"], summary.loc["Gapdh","median"]
print(f"\nKlrk1 RI median {kl:.1f}%  is  {kl/ac:.0f}x Actb ({ac:.2f}%)  and  {kl/ga:.0f}x Gapdh ({ga:.2f}%)")

summary.to_csv(os.path.join(cfg.TABLE_DIR, "housekeeping_negative_control.csv"))
fig, ax = plt.subplots(figsize=(6, 4.5))
order = ["Klrk1", "Actb", "Gapdh"]
ax.bar(order, [summary.loc[g, "median"] for g in order], color=["#c0392b", "#7f8c8d", "#95a5a6"])
ax.set_ylabel("Retained-intron proportion (%), median across samples")
ax.set_title("FigV5 — Housekeeping negative control", fontweight="bold")
for i, g in enumerate(order):
    ax.text(i, summary.loc[g, "median"], f"{summary.loc[g,'median']:.1f}%", ha="center", va="bottom")
fig.tight_layout()
fig.savefig(os.path.join(cfg.FIG_DIR, "FigV5_housekeeping_negative_control.pdf"))
fig.savefig(os.path.join(cfg.FIG_DIR, "FigV5_housekeeping_negative_control.png"), dpi=200)
print(f"\nfigure -> {cfg.FIG_DIR}/FigV5_housekeeping_negative_control.png")
