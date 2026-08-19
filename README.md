# Klrk1 (NKG2D) retained-intron isoform study — GVHD splicing project

Activation-dependent intron retention at the murine *Klrk1* (NKG2D) locus, most
notably the retained-intron isoform **Klrk1-203** (ENSMUST00000137660), across
four independent RNA-seq datasets. Transcript-level Salmon quantification plus
orthogonal alignment-level validation (junction/PIR, sashimi), housekeeping
control, CPAT/NMD, and structural mapping.

## Layout — each dataset is self-contained

```
datasets/
  GSE147371_GVHD_CD4/     GVHD CD4+ Tn/Tem (76 bp) ........ primary discovery
  GSE109125_ImmGen/       ImmGen CD8 naive→eff→mem (25 bp) . primary discovery
  GSE203167_Karimi_TCF7/  WT vs TCF-7 cKO CD8 (51 bp) ...... validation
  GSE119943_Yao_LCMV/     Yao LCMV CD8 (50 bp SE) .......... validation
    each: data/{fastq,bam,salmon}  scripts/  results/  figures/  README.md
  _extra_datasets/        GSE288143, GSE83978 (exploratory, not in manuscript)

shared/
  ref/                    reference genome + Salmon/HISAT2/STAR indices (not in git)
  scripts_common/         junction method (13,14) + helpers, used by both primary datasets
  metadata/               per-dataset sample sheets
  cross_dataset_analyses/ analyses spanning all datasets:
      scripts/  05_housekeeping_control.R, 07_cpat_nmd_analysis.sh, 12_RI_isoform_decomposition.py
      figures/  FigV5 (housekeeping), FigV7 (NMD)
      results/  Klrk1_203_vs_204_206_decomposition.csv

manuscript/          current ver15 + old_versions/
literature/          reference PDFs
related_analyses/    nfat_promoter_analysis (regulatory; proteomics kept separate, excluded)
```

To review one dataset's analysis, open its folder: `scripts/`, `figures/`, and
`results/` there are that dataset's only. The two primary datasets (GSE147371,
GSE109125) were analyzed together (Fig1–4 plot both), so those joint scripts and
figures are copied into both folders.

## Reproducing
Per-dataset: run the numbered scripts in that dataset's `scripts/` (Salmon → R
analysis → junction). Cross-dataset steps live in `shared/`. Read lengths differ
(25/50/51/76 bp); the junction-spanning method requires ≥50 bp (see the ImmGen
25 bp caveat in `datasets/GSE109125_ImmGen/README.md`).

> Machine-specific absolute paths inside scripts still need to be made relative
> before a clean external re-run; data/ and ref/ are git-ignored (large).
