# Klrk1 isoform analysis — pure-Python pipeline

This folder reproduces the **entire Klrk1 (NKG2D) retained-intron analysis using
only Python**. If you have never run this before, follow the steps below in
order — you do not need to understand the biology or the code to reproduce the
results.

> What the analysis shows, in one sentence: mouse T cells switch on a
> retained-intron form of the *Klrk1* gene (called **Klrk1-203**) when they
> become activated/effector cells, and this pipeline measures that from public
> RNA-seq data.

---

## What you need (once)

1. **Python 3.10+** (tested on 3.12). Check with:
   ```bash
   python3 --version
   ```
2. **Install the Python packages** (one command):
   ```bash
   cd python_pipeline
   pip install -r requirements.txt
   ```
   (pandas, numpy, matplotlib, and pysam — pysam is only for the junction step.)
3. **Tell the pipeline where the project is.** Everything is found relative to
   one folder. Set it once in your terminal (change the path to where this
   project lives on your computer):
   ```bash
   export KLRK1_ROOT=/Users/study/Desktop/Karimi/Klrk1_GVHD_project
   ```
   That's the only path you ever change.

---

## Where the data lives

The analysis reads **Salmon** output (transcript counts) from:
```
KLRK1_ROOT/datasets/<DATASET>/data/salmon/<sample>/quant.sf
```
Three datasets' Salmon output are already present (GSE147371, GSE109125,
GSE203167). GSE119943's data is **not** on this machine — its script will simply
skip with a message until you re-download it.

If you are starting completely from scratch (no Salmon output yet), first run the
download+quantify shell scripts in `../datasets/*/scripts/01_*.sh` and
`03_*.sh` (these use the `salmon` program). Otherwise, skip straight to the
Python steps below.

---

## Run the analysis (step by step)

From inside `python_pipeline/`, run these in order. Each prints a summary and
saves figures to `outputs/figures/` and tables to `outputs/tables/`.

```bash
# 1) Primary discovery: GVHD CD4 (GSE147371) + ImmGen CD8 (GSE109125)
python3 scripts/02_primary_analysis.py
#    -> Fig1 (Klrk1-203 induction), Fig2 (CD8 kinetics), Fig3 (isoform
#       proportions), Fig4 (total Klrk1) + master/supplementary tables

# 2) Validation: Karimi WT vs TCF-7 cKO CD8 (GSE203167)
python3 scripts/04_gse203167_analysis.py
#    -> FigV1-4

# 3) Negative control: is Klrk1's retained-intron % really higher than
#    housekeeping genes Actb/Gapdh? (rules out artifacts)
python3 scripts/05_housekeeping_control.py
#    -> FigV5 (expect Klrk1 ~25% vs Actb ~0.3% vs Gapdh ~2.6%)

# 4) Split the retained-intron signal into Klrk1-203 vs 204/206 (all datasets)
python3 scripts/07_ri_decomposition.py
#    -> decomposition table

# 5) Alignment-level check: count reads spanning the intron-4 splice junction
#    and intron-4 coverage, compared to the Salmon estimate.
#    (needs genome-aligned BAMs in datasets/<DATASET>/data/bam/ — produced by
#     ../shared/scripts_common/14_align_junction.sh)
python3 scripts/06_junction_counting.py
#    -> per-dataset junction_intron4_counts.csv

# 6) Validation: Yao LCMV CD8 (GSE119943) — only if you have fetched its data
python3 scripts/09_gse119943_analysis.py
```

That's the whole pipeline. Nothing writes outside `python_pipeline/outputs/`.

---

## What each script is (plain English)

| Script | What it answers |
|--------|-----------------|
| `02_primary_analysis.py` | How much Klrk1-203 is in each T-cell state? (main result) |
| `04_gse203167_analysis.py` | Does an independent lab's data show the same? |
| `05_housekeeping_control.py` | Is the signal real, not a contamination artifact? |
| `06_junction_counting.py` | Do individual reads directly confirm intron-4 retention? |
| `07_ri_decomposition.py` | Is it specifically Klrk1-203 (not the other RI forms)? |
| `09_gse119943_analysis.py` | A fourth dataset (infection model) — extra confirmation |
| `config.py` | The one file that holds all paths + coordinates |

## Notes
- **Read length matters (GSE109125 = 25 bp):** at 25 bp, reads are too short to
  span the splice junction, so the direct junction count is ~0 — this is a
  read-length limit, **not** absence of signal. Step 6 also reports an intron-4
  **coverage** ratio, which does work at 25 bp and tracks the Salmon estimate.
- Numbers here match the R pipeline exactly; figures are drawn with matplotlib so
  they look a little different from the R/ggplot manuscript figures but show the
  same data.
- Re-running a script overwrites its own outputs; it never deletes data.
