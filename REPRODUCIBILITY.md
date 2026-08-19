# Reproducibility guide

Everything needed to re-run the analysis from public data. All scripts derive
their paths from one place, so moving to a new machine is a single change.

## 1. Set the project root (the ONE thing to change on a new machine)

Scripts read `PROJECT_ROOT` from the `KLRK1_ROOT` environment variable, falling
back to a hard-coded default. Either:

```bash
export KLRK1_ROOT=/path/to/Klrk1_GVHD_project     # shell scripts + R pick this up
```
or edit the default line at the top of `shared/scripts_common/config.R`.

`config.R` then derives every path: `salmon_quant(sample, dataset)`,
`all_salmon_quants()`, `results_dir(ds)`, `figures_dir(ds)`, `metadata_csv(name)`.

## 2. Software

| Tool | Version | Used for | Where |
|------|---------|----------|-------|
| R | 4.6.0 | isoform analysis, figures | system `Rscript` |
| tidyverse | 2.0.0 | all R scripts (only R dependency) | `install_dependencies.R` |
| Salmon | 1.10.3 | transcript quantification | `salmon` (homebrew) |
| HISAT2 | 2.2.2 | genome alignment for junctions | conda env `star_env` |
| samtools | 1.19+ | BAM handling | `star_env` |
| pysam | 0.23.3 | junction/coverage counting | `star_env` |
| sra-tools | 3.3.0 | download FASTQ from SRA (`prefetch`, `fastq-dump`) | homebrew |

```bash
Rscript shared/scripts_common/install_dependencies.R   # R packages
conda activate star_env                                # alignment/junction tools
```

## 3. Data (not in git — regenerate)

`datasets/*/data/` and `shared/ref/` are git-ignored (large). Each dataset's
`README.md` lists its GEO accession and run (SRR) IDs. Metadata sample sheets in
`shared/metadata/` map every sample to its SRR and Salmon folder name.

> GSE119943 data is not present locally (it lives on the other lab Mac). Re-fetch
> into `datasets/GSE119943_Yao_LCMV/data/salmon/` to reproduce its analysis.

## 4. Run order

Per dataset (from its `scripts/` folder):
```
# primary discovery (GSE147371 + GSE109125, analyzed together)
bash    01_salmon_primary.sh              # quant (downloads data)
Rscript 02_klrk1_analysis_primary.R       # Fig1-4 + primary tables (mirrored to both)
bash    14_align_junction.sh <LABEL> <READLEN> <SRR..>   # HISAT2 chr6 alignment
python  13_junction_primary_datasets.py   # intron-4 junction/coverage vs Salmon

# validation datasets
Rscript 04_klrk1_analysis_GSE203167.R  |  06_sashimi_HISAT2.sh  |  11_..._reconciliation.py
Rscript 09_klrk1_analysis_GSE119943.R  |  10_sashimi_GSE119943.sh
```
Cross-dataset (`shared/cross_dataset_analyses/scripts/`):
```
Rscript 05_housekeeping_control.R        # Actb/Gapdh negative control (FigV5)
bash    07_cpat_nmd_analysis.sh          # coding potential + NMD 50-nt rule
python  12_RI_isoform_decomposition.py   # split RI into 203 vs 204/206
```

## 5. Read-length caveat (important)
Junction-spanning (N-CIGAR) counting needs ≥50 bp reads. GSE109125 (ImmGen) is
25 bp, so its junction-spanning count is ~0 (a read-length artifact); the
read-length-independent intron-4 coverage ratio is used there instead. See
`shared/scripts_common/README_junction_extension.md`.

## 6. Provenance
R session details captured in `shared/scripts_common/SESSION_INFO.txt`.
