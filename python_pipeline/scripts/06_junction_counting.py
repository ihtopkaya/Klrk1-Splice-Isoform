#!/usr/bin/env python3
# ============================================================================
# 13_junction_primary_datasets.py
#
# PURPOSE
#   Extend the alignment-level (junction / PIR) confirmation of the Salmon
#   transcript-level Klrk1-203 fractions to the two remaining datasets that were
#   previously Salmon-only:
#       GSE147371 (GVHD CD4+, 76 bp paired-end)   -> junction counting is valid
#       GSE109125 (ImmGen,     25 bp paired-end)  -> reported with a read-length
#                                                     caveat (25 bp is short for
#                                                     junction-spanning reads)
#
#   Mirrors the intron-4 method of 11_GSE203167_PIR_junction_reconciliation.py so
#   the numbers are directly comparable across all four datasets.
#
# METHOD (per BAM)
#   - E4-E5 junction (spliced) reads: reads with an N (CIGAR op 3) whose skipped
#     block matches intron-4 boundaries within +/-TOL bp.
#   - Intron-4 retention: reads whose contiguous M blocks span the 5' (B5) or 3'
#     (B3) exon/intron boundary; ret = mean(ret5, ret3).
#   - PIR_boundary = 100 * ret / (ret + spliced)   [read-based, isoform-free]
#   - intron-4 mean coverage (samtools-style) for context.
#   - Salmon Klrk1-203 fraction (%) for the same sample, for concordance.
#
# INPUT is a hard-coded table of (label, bam, salmon_dir, read_len, dataset).
# OUTPUT results/tables/junction_primary_GSE147371_GSE109125.csv + console.
#
# REQUIREMENTS: run with star_env python (has pysam):
#   /opt/anaconda3/envs/star_env/bin/python scripts/13_junction_primary_datasets.py
# ============================================================================

import csv, os, sys
import pysam

# ---- paths + coordinates come from the shared config -----------------------
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config as cfg

CHROM, TOL = cfg.CHROM, cfg.TOL
I4_START, I4_END = cfg.I4_START, cfg.I4_END          # intron 4
B5, B3 = cfg.B5, cfg.B3                               # exon/intron boundaries
RI203_START, RI203_END = cfg.RI203_START, cfg.RI203_END  # 203-retained portion
EX4_START, EX4_END = cfg.EX4_START, cfg.EX4_END      # constitutive exon 4
KLRK1 = cfg.KLRK1
# BAMs are read from each dataset's own data/bam/ folder (cfg.bam_dir(ds)); the
# coverage ratio (RI203 region / exon 4) is a read-length-independent PIR-like
# estimate that still works at 25 bp where junction-spanning reads vanish.

# label, bam file (in BAMDIR), salmon_dir, read_len, dataset
# extend this list as more BAMs become available (aligned by 14_align_*.sh)
JOBS = [
    # ImmGen (GSE109125) — 25 bp PE
    ("ImmGen_CD8_Effector", "ImmGen_CD8_Effector.bam", "CD8_Effector_1", 25, "GSE109125"),
    ("ImmGen_CD8_Tcm",      "SRR6467039_Tcm.bam",      "CD8_Tcm_1",      25, "GSE109125"),
    ("ImmGen_CD8_Naive",    "SRR6467045_Naive.bam",    "CD8_Naive_2",    25, "GSE109125"),
    # GSE147371 (GVHD CD4+) — 76 bp PE (added after 14_align_junction.sh runs)
    ("GVHD_CD4_Tem_1",  "GVHD_CD4_Tem_1.bam",  "Tem_GSM4427964_22_23", 76, "GSE147371"),
    ("GVHD_CD4_Tem_2",  "GVHD_CD4_Tem_2.bam",  "Tem_GSM4427966_24_25", 76, "GSE147371"),
    ("GVHD_CD4_Tem_3",  "GVHD_CD4_Tem_3.bam",  "Tem_GSM4427968_26_27", 76, "GSE147371"),
    ("GVHD_CD4_Tn",     "GVHD_CD4_Tn.bam",     "Tn_GSM4427970_28_29",  76, "GSE147371"),
]

def salmon_203_fraction(salmon_dir, ds):
    p = cfg.salmon_quant(salmon_dir, ds)
    if not os.path.exists(p):
        return None
    d = {}
    with open(p) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            t = r["Name"].split(".")[0]
            if t in KLRK1:
                d[KLRK1[t].split("-")[-1]] = float(r["TPM"])   # key by short code (201..206)
    tot = sum(d.values())
    return round(100 * d.get("203", 0) / tot, 2) if tot else 0.0

def analyze_bam(path, read_len):
    bam = pysam.AlignmentFile(path, "rb")
    spliced = ret5 = ret3 = 0
    for r in bam.fetch(CHROM, I4_START - 200, I4_END + 200):
        if r.is_unmapped:
            continue
        ref = r.reference_start
        is_spliced = False
        for op, ln in r.cigartuples:
            if op in (0, 2, 7, 8):        # M/D/=/X consume reference
                ref += ln
            elif op == 3:                 # N (skip) = intron
                n_start, n_end = ref + 1, ref + ln
                if abs(n_start - I4_START) <= TOL and abs(n_end - I4_END) <= TOL:
                    is_spliced = True
                ref += ln
        if is_spliced:
            spliced += 1
        blocks = r.get_blocks()           # contiguous M blocks (0-based ref)
        if any(s < B5 < e for s, e in blocks): ret5 += 1   # spans exon4/intron4 boundary (diagnostic of intron-4 retention)
        if any(s < B3 < e for s, e in blocks): ret3 += 1   # spans intron4/exon5 boundary — NOT reachable by 203/204 (they retain only 129593531+); kept for diagnostics only
    def meancov(a, b):
        cov = bam.count_coverage(CHROM, a - 1, b, quality_threshold=0)
        return sum(sum(c) for c in cov) / (b - a + 1)
    i4_cov   = meancov(I4_START, I4_END)         # whole intron 4 (context)
    ri203_cov = meancov(RI203_START, RI203_END)  # 203-retained portion
    ex4_cov   = meancov(EX4_START, EX4_END)      # constitutive flanking exon
    bam.close()
    # intron-4 retention numerator = ret5 only (exon4/intron4 boundary). ret3 (intron4/exon5
    # boundary) is NOT diagnostic: 203 retains 129593531-631 and 204 129593456+, neither reaching
    # B3=129593296, so ret3 counts Klrk1-205 exon-5-extension / exon3->exon5 spliced reads. [QA fix]
    ret = float(ret5)
    pir_boundary = 100 * ret / (ret + spliced) if (ret + spliced) else 0
    cov_ratio = 100 * ri203_cov / ex4_cov if ex4_cov else 0   # read-length-independent PIR-like
    return (spliced, round(i4_cov, 1), round(ret, 1), round(pir_boundary, 1),
            round(ri203_cov, 1), round(ex4_cov, 1), round(cov_ratio, 1))

HEADER = ["sample", "dataset", "read_len_bp",
          "E4E5_junction_spanning_reads", "intron4_boundary_retention_reads",
          "PIR_junction_based_pct",
          "RI203_region_meancov", "exon4_meancov",
          "coverage_retention_ratio_pct", "Salmon_203_fraction_pct"]

rows = []
for label, bamname, sdir, rlen, ds in JOBS:
    bam = os.path.join(cfg.bam_dir(ds), bamname)
    if not os.path.exists(bam):
        print(f"[warn] missing BAM {bam}")
        continue
    spliced, i4cov, ret, pir_b, ri203cov, ex4cov, covratio = analyze_bam(bam, rlen)
    salmon203 = salmon_203_fraction(sdir, ds)
    # junction-based PIR is only meaningful when junction-spanning reads exist
    # (they vanish at 25 bp); report NA there and rely on the coverage ratio.
    pir_report = pir_b if spliced > 0 else "NA(no junction reads at this read len)"
    rows.append([label, ds, rlen, spliced, ret, pir_report,
                 ri203cov, ex4cov, covratio, salmon203])

# write per-dataset results into each dataset's own results/ folder
for ds in sorted(set(r[1] for r in rows)):
    out = os.path.join(cfg.results_dir(ds), "junction_intron4_counts.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f); w.writerow(HEADER)
        w.writerows([r for r in rows if r[1] == ds])
    print(f"wrote {out}")
print()
hdr = (f"  {'sample':20s} {'ds':10s} {'rl':>3} {'E4E5j':>6} {'retR':>6} "
       f"{'PIRjunc':>8} {'RI203cov':>8} {'ex4cov':>7} {'covRatio%':>9} {'Salmon203%':>10}")
print(hdr)
for r in rows:
    pir = r[5] if isinstance(r[5], (int, float)) else "NA"
    print(f"  {r[0]:20s} {r[1]:10s} {r[2]:>3} {r[3]:>6} {r[4]:>6} "
          f"{str(pir):>8} {r[6]:>8} {r[7]:>7} {r[8]:>9} {r[9]:>10}")
print("\nNote: at 25 bp (ImmGen) junction-spanning reads are ~0 (a read-length")
print("artifact, not biology); the coverage_retention_ratio is the read-length-")
print("independent alignment-level check and should track Salmon_203_fraction.")
