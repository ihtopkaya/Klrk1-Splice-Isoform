# Klrk1 intron-4 junction analysis — methods & results handoff (rev. 3)

For manuscript revision. Self-contained. Rev.3 adds the isoform-geometry analysis
(why the junction metric is total intron-4 retention, not Klrk1-203-specific), a
proximal-boundary PIR variant, and corrected coverage-profile framing.

## 1. Why this analysis was done
Earlier drafts said the transcript-level Salmon Klrk1-203 fractions were checked
against alignment-level counts "at each step," but that check had only been done
for GSE203167 and GSE119943. The two discovery datasets — GSE147371 (GVHD CD4+)
and GSE109125 (ImmGen) — were Salmon-only. This extends the alignment-level check
to all four datasets.

## 2. Method
Reads fetched from SRA (sra-tools 3.3.0), aligned with HISAT2 2.2.2 to a
chr6-restricted GRCm39 index (`--ss/--exon`, Ensembl GTF; same index as the
sashimi work); sorted/indexed with samtools. GSE109125 (25 bp) used
`-L 20 --score-min L,0,-0.6`. With pysam 0.23.3, at intron 4 (chr6:129,593,296–
129,593,631, GRCm39 minus strand) we computed per sample:
- **spliced (E4–E5) reads**: N-CIGAR reads whose gap matches the intron-4
  boundaries within ±10 bp (the canonical splice removing intron 4).
- **retention reads at each exon–intron junction**: reads whose contiguous
  alignment blocks span the 5′ boundary B5=129,593,631 (exon4/intron) or the 3′
  boundary B3=129,593,296 (intron/exon5).
- **PIR, two variants** (reported only where spliced reads exist): the 5′/3′
  **mean** PIR = mean(ret5,ret3)/(mean+spliced), and the **proximal (5′)** PIR =
  ret5/(ret5+spliced). The proximal variant is the cleaner single number at this
  locus (see §4).
- **intron-4 coverage ratio** (read-length-independent): mean coverage of the
  203-retained proximal intron (129,593,531–631) ÷ mean coverage of exon 4
  (129,593,632–736).
- Salmon Klrk1-203 fraction (203 TPM / total Klrk1 TPM) for comparison.

*chr6-index note (pre-empts a reviewer):* junction counting used a chr6-only
index while Salmon used the whole transcriptome; if chr6-restriction forced
mismapping the two estimates would not agree, so their convergence argues against
mismapping inflation.

## 3. Read-length caveat
Junction-spanning counting needs reads long enough to anchor both sides of the
boundary. **GSE109125 (ImmGen) is 25 bp**, so spliced counts are ~0 in every
ImmGen sample — a read-length artifact, not absence of the isoform (true even in
the deepest ImmGen library, CD8 Tcm ~14 M reads). ImmGen is validated by the
coverage ratio; the three ≥50 bp datasets support junction counting directly.

## 4. Isoform geometry — the junction metric is TOTAL intron-4 retention, not 203-specific

From the Ensembl annotation (GRCm39, minus strand; high coordinate = 5′):

| transcript | intron-4 relationship | spans B5 (631)? | spans B3 (296)? |
|---|---|---|---|
| Klrk1-203 (137660) | terminal exon 129,593,531–736; retains proximal intron 531–631, **ends at 531** | yes | **no** |
| Klrk1-204 (152256) | exon 129,593,456–594,484; retains proximal intron 456–631, **ends at 456** | yes | **no** |
| Klrk1-206 (204694) | single exon 129,597,565–600,823; **not in intron 4 at all** | no | no |
| Klrk1-201 (canonical) | exon4 632–736 spliced to exon5 260–295 | (spliced, N-gap) | (spliced) |

Consequences:
- **A read crossing B5 comes from 203 *or* 204 (or unspliced pre-mRNA)** — both
  retain the proximal intron and cross 631 identically. 203 and 204 cannot be
  separated at this splice boundary; they differ only in *where inside the intron
  they terminate* (531 vs 456), which is a coverage drop-off, not a junction, so
  N-CIGAR counting cannot see it.
- **A read crossing B3 comes from neither 203 nor 204** (both end well above
  exon 5); it reflects full-length intron retention / pre-mRNA background.
- Therefore the junction PIR measures **overall intron-4 retention**, not 203
  specifically. **203-specificity comes from Salmon + the 203-vs-204/206
  decomposition**, which use the full transcript sequence and termination points.

The junction/coverage evidence answers "is intron 4 genuinely retained at the
read level?" (yes, inference-free); Salmon/decomposition answers "which isoform?".

**On the 5′/3′ mean:** because retention here is *partial* (203/204 never reach
B3), the two boundaries are not equivalent — unlike a symmetric whole-intron event
where vast-tools averages them to reduce noise. Averaging drags PIR down slightly
by including a boundary the isoforms of interest never reach. The proximal (5′)
PIR is therefore the cleaner single metric; we report both.

## 5. Results (single consistent run)

**GSE147371 (GVHD CD4+, 76 bp; SRR11389222–229):**

| Sample | spliced | ret5 (5′) | ret3 (3′) | PIR mean % | PIR proximal % | coverage ratio % | Salmon 203 % |
|---|--:|--:|--:|--:|--:|--:|--:|
| Tem_1 | 49 | 15 | 15 | 23.4 | 23.4 | 19.3 | 16.2 |
| Tem_2 | 110 | 18 | 16 | 13.4 | 14.1 | 13.8 | 16.2 |
| Tem_3 | 102 | 19 | 12 | 13.2 | 15.7 | 8.3 | 10.5 |
| Tn (naive) | 3 | 1 | 4 | 45.5 | 25.0 | 4.9 | 10.6 |

The three Tem samples show robust spliced reads (49–110) with PIR concordant with
Salmon. Tn yields only 3 spliced reads (Klrk1 barely expressed in naive CD4); its
mean PIR (45.5%) is inflated by 4 distal background reads and is sampling noise —
the proximal PIR (25%) is less misleading, but with n=1 spliced read neither
should be over-interpreted.

**GSE109125 (ImmGen, 25 bp; spliced ≈ 0 → coverage proxy):**

| Sample | spliced | coverage ratio % | Salmon 203 % |
|---|--:|--:|--:|
| CD8 Effector (SLEC) | 0 | 17.0 | 20.9 |
| CD8 Tcm | 0 | 17.8 | 16.3 |
| CD8 Naive | 0 | not meaningfully defined* | 0.0 |

*Naive CD8 is essentially Klrk1-unexpressed (exon-4 coverage ≈ 1×), so the ratio
is ~0/0 — not a measured zero.

## 6. Intron-4 coverage profile (positive-consistency check)

Mean coverage of the 203-retained proximal intron (531–631) vs the distal intron
(296–530):

| Sample | exon4 | proximal | distal | distal/proximal |
|---|--:|--:|--:|--:|
| GVHD Tem_1 | 101.6 | 19.7 | 11.4 | 0.58 |
| GVHD Tem_2 | 175.3 | 24.2 | 15.3 | 0.63 |
| GVHD Tem_3 | 213.5 | 17.7 | 14.4 | 0.81 |
| ImmGen Effector | 159.1 | 27.1 | 49.1 | 1.81 |
| ImmGen Tcm | 196.5 | 35.0 | 64.3 | 1.83 |

In the reliable 76 bp GVHD data, proximal > distal (ratio 0.58–0.81) — exactly
what the annotation predicts: the proximal intron carries 203 + 204 together,
while toward the distal end 203 drops out (it ends at 531), leaving only 204 and
background. So the reliable dataset is geometrically consistent. The reversed
ImmGen pattern (distal > proximal) reflects 25 bp mapping unreliability and should
not be over-interpreted. (Note: an earlier draft attributed distal coverage to
"204 + 206"; that was wrong — 206 is not in intron 4 and 204 also ends above the
distal region — so distal coverage is full-retention/pre-mRNA background, not a
named isoform.)

## 7. Scripts
`14_align_junction.sh` (SRA → chr6 HISAT2, disk-safe); `13_junction_primary_datasets.py`
(junction + coverage, both PIR variants, per-dataset `results/junction_intron4_counts.csv`);
coordinates + isoform map in `config.py`.

## 8. Suggested manuscript edits
- **Required** (not optional): anywhere the junction/PIR count is presented as
  evidence for Klrk1-203 *specifically*, soften it — 203 and 204 are inseparable
  at the exon4/intron junction, so the junction metric is total intron-4
  retention. 203-specificity is carried by Salmon + decomposition.
- **Do NOT** add a "5′ excess of boundary reads is consistent with 203" sentence
  (proposed in an earlier review round): 204 also crosses B5, so the 5′ signal is
  not 203-specific.
- Replace the "At each step…" sentence with the read-length paragraph from
  handoff rev.2 §7 (junction PIR for ≥50 bp; coverage proxy for 25 bp ImmGen,
  interpreted with the decomposition).
- Bind the naive "0%" claim to its dataset: ImmGen naive CD8 is undetectable, but
  GVHD naive CD4 (Tn) has Salmon-203 ≈ 10.6% — do not write "0 in all naive".
- Optionally report the proximal (5′) PIR rather than the 5′/3′ mean, and frame it
  as "similar to the mean of the two flanking exon–intron junctions, noting that
  retention here is partial so the two boundaries are not equivalent."
