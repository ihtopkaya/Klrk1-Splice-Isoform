# Junction-level confirmation for the two remaining datasets (GSE147371, GSE109125)

**Why this exists.** Earlier manuscript drafts stated that the transcript-level
Salmon Klrk1‑203 fractions were checked against direct, alignment-level counts
"at each step". In fact junction / intron‑4 read counting had only been run for
**GSE203167** (script `11_*.py`) and **GSE119943** (script `10_*.sh`). The two
discovery datasets — **GSE147371** (GVHD CD4+) and **GSE109125** (ImmGen) — were
Salmon-only. These scripts close that gap so the alignment-level check truly
covers all four datasets.

## Scripts

| Script | Role |
|---|---|
| `14_align_junction.sh` | Download one sample from SRA, align to the chr6 HISAT2 index, delete FASTQ. Disk-safe (`fastq-dump --gzip`, one run at a time, `.sra` and FASTQ removed after use). |
| `13_junction_primary_datasets.py` | Count E4–E5 junction-spanning (N‑CIGAR) reads, intron‑4 boundary-retention reads, and a read-length-independent coverage retention ratio at intron 4; compare each to the per-sample Salmon Klrk1‑203 fraction. |

Both reuse the intron‑4 coordinates and method of
`11_GSE203167_PIR_junction_reconciliation.py`, so all four datasets are directly
comparable.

## Coordinates (GRCm39, chr6, minus strand)

- intron 4 (spliced out by the canonical E4–E5 junction): `129,593,296–129,593,631`
- 203‑retained portion of intron 4 (terminal exon of Klrk1‑203 = `129,593,531–129,593,736`): `129,593,531–129,593,631`
- constitutive flanking exon 4: `129,593,632–129,593,736`
- `coverage_retention_ratio = mean_cov(203‑retained region) / mean_cov(exon 4)`

## Datasets and read length

| Dataset | Cells | Read length | Junction-spanning counting |
|---|---|---|---|
| GSE147371 | GVHD CD4+ Tn/Tem | 76 bp PE | **valid** (reads long enough to span the junction) |
| GSE109125 (ImmGen) | CD8 naive→effector→memory | **25 bp** PE | **not feasible** — junction-spanning reads ≈ 0 at 25 bp; use the coverage retention ratio instead |

**Key point for the manuscript.** At 25 bp a read cannot provide enough anchor on
both sides of the exon4/intron4 boundary, so the N‑CIGAR method returns ~0
junction reads for ImmGen — this is a read-length artifact, **not** absence of
signal. The isoform-free **coverage retention ratio** is the correct
alignment-level check at short read length, and it tracks the Salmon Klrk1‑203
fraction (see results below).

## How to reproduce

Requires the `star_env` conda env (hisat2, samtools, pysam) and sra-tools
(`prefetch`, `fastq-dump`). The chr6 HISAT2 index built for the GSE203167 sashimi
work is reused from `HOME GVHD_splicing_project/star_analysis/hisat_index/`.

```bash
conda activate star_env
cd Klrk1_NKG2D_isoform_analysis

# --- ImmGen (GSE109125, 25 bp) ---
# Tcm and Naive BAMs already existed from the sashimi work; effector needs aligning:
bash 14_align_junction.sh ImmGen_CD8_Effector 25 SRR6467043 GSE109125

# --- GSE147371 (GVHD CD4+, 76 bp); each GSM = 2 runs, merged ---
bash 14_align_junction.sh GVHD_CD4_Tem_1 76 SRR11389222,SRR11389223 GSE147371
bash 14_align_junction.sh GVHD_CD4_Tem_2 76 SRR11389224,SRR11389225 GSE147371
bash 14_align_junction.sh GVHD_CD4_Tem_3 76 SRR11389226,SRR11389227 GSE147371
bash 14_align_junction.sh GVHD_CD4_Tn 76 SRR11389228,SRR11389229 GSE147371

# --- count junctions / coverage vs Salmon ---
/opt/anaconda3/envs/star_env/bin/python scripts/13_junction_primary_datasets.py
```

Output: `results/tables/junction_primary_GSE147371_GSE109125.csv`.

## Sample → run accession map

| Label | GEO | Runs | Salmon dir |
|---|---|---|---|
| ImmGen_CD8_Effector | GSM2932627 | SRR6467043 | CD8_Effector_1 |
| ImmGen_CD8_Tcm | GSM2932623 | SRR6467039 | CD8_Tcm_1 |
| ImmGen_CD8_Naive | GSM2932629 | SRR6467045 | CD8_Naive_2 |
| GVHD_CD4_Tem_1 | GSM4427964 | SRR11389222+223 | Tem_GSM4427964_22_23 |
| GVHD_CD4_Tem_2 | GSM4427966 | SRR11389224+225 | Tem_GSM4427966_24_25 |
| GVHD_CD4_Tem_3 | GSM4427968 | SRR11389226+227 | Tem_GSM4427968_26_27 |
| GVHD_CD4_Tn | GSM4427970 | SRR11389228+229 | Tn_GSM4427970_28_29 |

## Results (fill in / regenerate from the CSV)

ImmGen (25 bp) — junction-spanning reads = 0 (read-length artifact); coverage
retention ratio vs Salmon Klrk1‑203 fraction:

| Sample | E4E5 junction reads | coverage retention ratio | Salmon 203 % |
|---|---|---|---|
| ImmGen_CD8_Effector | 0 | 17.0 % | 20.9 % |
| ImmGen_CD8_Tcm | 0 | 17.8 % | 16.3 % |
| ImmGen_CD8_Naive | 0 | ~0 (Klrk1 unexpressed) | 0 % |

GSE147371 (76 bp) — junction-spanning counting valid: _pending alignment; regenerate table_.
