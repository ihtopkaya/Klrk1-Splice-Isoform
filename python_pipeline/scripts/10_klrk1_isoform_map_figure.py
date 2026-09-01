#!/usr/bin/env python3
# ============================================================================
# 10_klrk1_isoform_map_figure.py
# Publication figure of the Klrk1 isoform architecture: the six annotated
# transcripts drawn to genomic scale (minus strand, 5'->3' left-to-right),
# with transcript-EXCLUSIVE exonic regions highlighted, the key discriminating
# junctions marked, and the five read-identifiability reporting units bracketed.
# This is the visual behind config.REPORT_UNITS / docs/KLRK1_202_ANALYSIS_NOTE.md.
# Run: python3 python_pipeline/scripts/10_klrk1_isoform_map_figure.py
# ============================================================================
import os, sys, collections
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config as cfg

GTF = os.path.join(cfg.PROJECT_ROOT, "shared", "ref", "genome", "Klrk1.gtf")
NAME = cfg.SHORT
ALL = ["201", "202", "205", "203", "204", "206"]          # PC group then RI group
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
# exclusive base set per transcript
def bset(iv):
    s = set()
    for a, b in iv: s.update(range(a, b + 1))
    return s
bs = {t: bset(ex[t]) for t in ALL}
excl = {t: bs[t] - set().union(*[bs[o] for o in ALL if o != t]) for t in ALL}
def frac_excl(a, b, t):     # fraction of an exon that is exclusive to t
    seg = set(range(a, b + 1)); return len(seg & excl[t]) / len(seg)

BIO = {"201": "protein-coding", "202": "protein-coding", "205": "protein-coding",
       "203": "retained-intron", "204": "retained-intron", "206": "retained-intron"}
GMIN, GMAX = 129586900, 129601200
PC, RI, EXC = "#3b6ea5", "#c0392b", "#f1c40f"

fig, ax = plt.subplots(figsize=(12, 6.2))
yrows = {t: len(ALL) - i for i, t in enumerate(ALL)}
tr = ax.get_yaxis_transform()   # x = axes fraction [0..1], y = data (row)
for t in ALL:
    y = yrows[t]; exs = sorted(ex[t])
    col = PC if BIO[t] == "protein-coding" else RI
    ax.plot([exs[0][0], exs[-1][1]], [y, y], color="#888", lw=1, zorder=1)
    for a, b in exs:
        ax.add_patch(Rectangle((a, y - 0.28), b - a, 0.56, facecolor=col,
                               edgecolor="black", lw=0.5, zorder=2))
        if frac_excl(a, b, t) > 0.5:   # exon segment exclusive to this transcript
            eseg = sorted(excl[t] & set(range(a, b + 1)))
            if eseg:
                ax.add_patch(Rectangle((eseg[0], y - 0.28), eseg[-1] - eseg[0], 0.56,
                             facecolor=EXC, edgecolor="black", lw=0.5, hatch="///", zorder=3))
    # transcript name on the LEFT (axes-fraction x, immune to x-data inversion)
    ax.text(0.008, y, f"Klrk1-{t}", transform=tr, va="center", ha="left", fontsize=10,
            color=col, fontweight="bold" if t == "203" else "normal", clip_on=False)

# key discriminating features (annotate; dy alternates to avoid collisions)
def mark(x, y, txt, dy=0.6):
    ax.annotate(txt, xy=(x, y + (0.28 if dy > 0 else -0.28)), xytext=(x, y + dy),
                fontsize=7.5, ha="center", va="center", color="#222",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="#222"))
mark(129599600, yrows["201"], "canonical first\nexon (201)", dy=0.62)
mark(129600790, yrows["202"], "alt first exon\n(202/205/206)", dy=-0.6)
mark(129588300, yrows["205"], "205-exclusive 3′ tail", dy=0.6)
mark(129593580, yrows["203"], "retained intron-4\n(203/204)", dy=0.6)
mark(129594090, yrows["204"], "intron-3 retained\n(204 only)", dy=-0.62)
mark(129598600, yrows["206"], "206-exclusive exon", dy=0.6)

# reporting-unit brackets on the RIGHT (axes-fraction x)
unit_rows = {"Klrk1-201": ["201"], "Klrk1-205/202": ["205", "202"],
             "Klrk1-203": ["203"], "Klrk1-204": ["204"], "Klrk1-206": ["206"]}
for u, mem in unit_rows.items():
    ys = [yrows[m] for m in mem]
    ax.plot([0.945, 0.945], [min(ys) - 0.3, max(ys) + 0.3], transform=tr, color="black",
            lw=1.4, clip_on=False)
    ax.text(0.955, sum(ys) / len(ys), u, transform=tr, va="center", ha="left", fontsize=9,
            fontweight="bold" if u == "Klrk1-203" else "normal", clip_on=False)
ax.text(0.945, len(ALL) + 0.75, "reporting units", transform=tr, fontsize=8,
        style="italic", color="#555", ha="left")

from matplotlib.patches import Patch
ax.legend(handles=[Patch(facecolor=PC, label="protein-coding exon"),
                   Patch(facecolor=RI, label="retained-intron exon"),
                   Patch(facecolor=EXC, hatch="///", label="transcript-exclusive segment")],
          loc="lower left", bbox_to_anchor=(0.0, -0.02), fontsize=8, frameon=False, ncol=3)
# high coord on the left (minus strand 5'->3'); wide margins for side labels
ax.set_xlim(GMAX + 4300, GMIN - 5200); ax.set_ylim(0.2, len(ALL) + 1.3)
ax.set_yticks([]); ax.set_xlabel("chr6 coordinate (GRCm39, minus strand; 5′→3′ left to right)", fontsize=9)
ax.set_title("Klrk1 isoform architecture and read-identifiability reporting units",
             fontsize=12, fontweight="bold")
for s in ["top", "right", "left"]:
    ax.spines[s].set_visible(False)
fig.tight_layout()
out = os.path.join(cfg.FIG_DIR, "FigX_Klrk1_isoform_map")
fig.savefig(out + ".png", dpi=200); fig.savefig(out + ".pdf")
print("wrote", out + ".png / .pdf")
