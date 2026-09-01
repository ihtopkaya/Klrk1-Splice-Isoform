# Klrk1-202 re-evaluation — reproducibility bundle

All code behind the Klrk1-202 (ENSMUST00000095412) re-evaluation and the read-identifiability reporting
convention, collected in one place. Full write-up + conclusions: [`docs/KLRK1_202_ANALYSIS_NOTE.md`](../../docs/KLRK1_202_ANALYSIS_NOTE.md).

**Question.** ver16 excluded Klrk1-202 globally as a "bias artefact" (uncorrected TPM inflated 0→52).
Re-tested against the **bias-corrected** quant.sf of every dataset: with correction ON, 202 is large in
GSE203167/GSE119943 (11–38%) and rises 0→63% across the GSE288143 activation time-course. 202 is a real
protein-coding alt-promoter signal but has **no read that separates it from Klrk1-205** (202 ⊆ 205), so it
is reported inside a combined **Klrk1-205/202** unit and kept in the total-Klrk1 denominator.

## Run order

| # | File | What it does | Needs |
|---|---|---|---|
| 0 | `00_verify_bias_flags.sh` | audit `--gcBias/--seqBias` per Salmon run (cmd_info.json) | bash |
| 1 | `01_reexamine_202_all_datasets.py` | 202% per condition; 203%/RI% with vs without 202 (5- vs 6-isoform denominator) | python3+numpy |
| 2 | `02_exon_overlap_correlation.py` | exclusive/shared exonic bp; 202⊆205; 203-retention shared only with 204; 201/202/205 correlation | python3+numpy, GTF, (03 output) |
| 3 | `03_klrk1_202_reevaluation.py` | master + decomposition tables incl. 202; **unit-level** tables; per-sample 6-isoform table | python3+numpy, config.py |
| 4 | `04_klrk1_identifiability_map.py` | region partition + junction inventory + `REPORT_UNITS` basis | python3, config.py, GTF |
| 5 | `05_bam_coverage_202_and_203vs204.sh` | alignment-level: alt-first-exon coverage/splice reads (202/205/206); 203-vs-204 retention | samtools, BAMs |

Scripts 03 and 04 are snapshots of the numbered pipeline scripts `python_pipeline/scripts/08_*` and `09_*`
(canonical copies live there and import `python_pipeline/config.py`, which defines the six transcripts and
`REPORT_UNITS`). 00/01/02/05 are standalone (set `KLRK1_ROOT` if the repo is elsewhere).

## Key outputs (in `python_pipeline/outputs/tables/`)

- `Table4_units_master.csv`, `Klrk1_decomposition_units.csv` — **manuscript-ready** (5 reporting units).
- `Table4b_master_isoform_202included.csv`, `Supplementary_all_sample_6isoform.csv` — raw 6-transcript (supplement).
- `Klrk1_region_partition.csv`, `Klrk1_junction_inventory.csv`, `Klrk1_identifiability_units.csv` — the read-identifiability map.
- `Klrk1_201_202_pool_stability.csv` — within-condition CV of 201 / 202 / 201+202.

## The five reporting units (config.REPORT_UNITS)

`Klrk1-201 · Klrk1-205/202 · Klrk1-203 · Klrk1-204 · Klrk1-206`. Total Klrk1 = sum of the five
(= all six transcripts). RI fraction = (203+204+206)/total. 205/206/204 have exclusive read support;
203 is separable via the {203,204} retention excess + the 203-type intron-3 splice junction; 201 by its
canonical first exon; 202 has no separating read and is merged with 205.
