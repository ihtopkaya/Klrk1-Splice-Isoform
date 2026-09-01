# Klrk1 (NKG2D) retained-intron isoform analysis

Code and derived results for the analysis of activation-dependent intron
retention at the murine *Klrk1* (NKG2D) locus, focused on the retained-intron
isoform **Klrk1-203** (ENSMUST00000137660), across four independent public
RNA-seq datasets.

The analysis combines transcript-level quantification (Salmon) with orthogonal
alignment-level validation: splice-junction and intron-4 coverage counting,
a housekeeping-gene negative control, coding-potential/NMD assessment, and a
read-identifiability map of the *Klrk1* locus.

## Data availability

All primary data are public. Raw reads were obtained from the Gene Expression
Omnibus; no new sequencing data were generated for this study.

| Accession | Description | Cells | Read length | n |
|---|---|---|---|--:|
| [GSE147371](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147371) | Allogeneic GVHD CD4+ Tn / Tem | CD4+ | 76 bp PE | 4 |
| [GSE109125](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109125) | ImmGen CD8 naive → effector → memory | CD8+ | 25 bp PE | 15 |
| [GSE203167](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203167) | WT vs TCF-7 cKO CD8, pre / post transplant | CD8+ | 51 bp | 12 |
| [GSE119943](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119943) | LCMV Armstrong / Clone-13 CD8 | CD8+ | 50 bp SE | 19 |

Sample sheets in `metadata/` map every sample to its GSM and SRR run accession,
condition, and quantification directory, so each dataset can be re-fetched and
re-quantified from the accessions alone.

Reference annotation: GRCm39, Ensembl release 111.

## Repository layout

```
python_pipeline/          self-contained Python pipeline (primary entry point)
  config.py                 all paths, transcript IDs and locus coordinates
  scripts/                  numbered analysis scripts (02-14)
  outputs/tables/           15 result tables (CSV)
  outputs/figures/          17 figures (PNG + PDF)

cross_dataset_analyses/   analyses spanning all four datasets
  scripts/                  housekeeping control (R), CPAT/NMD, RI decomposition
  results/  figures/

analyses/
  klrk1_202_reevaluation/   Klrk1-202 identifiability re-evaluation and the
                            read-identifiability reporting convention

scripts_common/           shared configuration (config.py/.R/.sh) and the
                          junction-counting method (13, 14)

related_analyses/
  nfat_promoter_analysis/   in silico NFAT / AP-1 promoter motif scan

metadata/                 per-dataset sample sheets (sample → GSM → SRR)
logs_provenance/          Salmon run logs for the quantification steps
```

## Reproducing the analysis

`python_pipeline/` reproduces the full isoform analysis in Python alone and is
the recommended entry point.

**1. Install dependencies**

```bash
cd python_pipeline
pip install -r requirements.txt     # pandas, numpy, matplotlib, pysam
```

**2. Point the pipeline at your copy of the project**

Every path is derived from a single environment variable:

```bash
export KLRK1_ROOT=/path/to/this/repository
```

**3. Run the scripts in order**

```bash
python3 scripts/02_primary_analysis.py     # GSE147371 + GSE109125 — main result
python3 scripts/04_gse203167_analysis.py   # validation: WT vs TCF-7 cKO
python3 scripts/05_housekeeping_control.py # negative control (Actb / Gapdh)
python3 scripts/07_ri_decomposition.py     # split RI into 203 vs 204 / 206
python3 scripts/06_junction_counting.py    # alignment-level junction check
python3 scripts/09_gse119943_analysis.py   # validation: LCMV CD8
```

Each script prints a summary and writes to `outputs/tables/` and
`outputs/figures/`. Nothing is written outside `python_pipeline/outputs/`.

The R scripts in `cross_dataset_analyses/scripts/` and the figure script
`python_pipeline/scripts/14_combined_figures_R.R` require R with `tidyverse`;
install via `scripts_common/install_dependencies.R`.

## Inputs not included in this repository

Raw FASTQ files, genome-aligned BAM files, the reference genome and the
Salmon/HISAT2 indices are too large to distribute here. They are regenerated
from the GEO/SRA accessions listed above:

```bash
# quantification input expected by the pipeline
<KLRK1_ROOT>/datasets/<GSE>/data/salmon/<sample>/quant.sf
```

Salmon was run with `--gcBias --seqBias` throughout; the run logs in
`logs_provenance/salmon_run_logs/` record the exact parameters used for each
sample. Genome alignment for junction counting uses
`scripts_common/14_align_junction.sh` (HISAT2 against a chr6-restricted GRCm39
index), which downloads one run at a time and removes the FASTQ afterwards.

## Read-length caveat

Junction-spanning (N-CIGAR) read counting requires reads long enough to anchor
both sides of a splice boundary, in practice ≥50 bp. GSE109125 (ImmGen) is
25 bp, so its junction-spanning counts are ≈0 — a read-length limitation rather
than absence of signal. For that dataset the read-length-independent intron-4
coverage ratio is used instead, and it tracks the transcript-level estimate.
See `scripts_common/README_junction_extension.md`.

## Software

| Tool | Version |
|---|---|
| Python | 3.12 (3.10+ supported) |
| Salmon | 1.10.3 |
| HISAT2 | 2.2.2 |
| samtools | 1.19 |
| pysam | 0.23.3 |
| sra-tools | 3.3.0 |
| R | 4.6.0 |
| tidyverse | 2.0.0 |

R session details are recorded in `scripts_common/SESSION_INFO.txt`.
