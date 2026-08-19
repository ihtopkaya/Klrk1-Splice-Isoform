#!/bin/bash
# ============================================================================
# 05_bam_coverage_202_and_203vs204.sh
# Model-free (alignment-level) evidence used in the re-evaluation:
#   (1) GSE203167 WT BAMs: real read coverage + spliced reads over the
#       alternative first exon shared by 202/205/206 (chr6:129,600,774-804)
#       -> the alt-promoter protein-coding signal is genuine, not a GC/seq artefact.
#   (2) GSE147371 Tem vs Tn BAMs: intron-4-proximal retention (203+204) exceeds
#       the 204-exclusive intron-3 retention, and is ~0 in naive -> a 203-type
#       contribution is required and is differentiation-induced.
# Requires samtools. Run:
#   bash analyses/klrk1_202_reevaluation/05_bam_coverage_202_and_203vs204.sh
# ============================================================================
set -euo pipefail
ROOT="${KLRK1_ROOT:-/Users/study/Desktop/Karimi/Klrk1_GVHD_project}"
SAM="$(command -v samtools || echo /opt/homebrew/bin/samtools)"
cd "$ROOT"
depth () { $SAM depth -a -r "$1" "$2" 2>/dev/null | awk '{s+=$3;n++} END{printf "%.2f",(n?s/n:0)}'; }

echo "### (1) GSE203167 WT: alt first exon (202/205/206-shared) 6:129600774-129600804 ###"
for b in WT_Pre WT_Post7; do
  BAM="datasets/GSE203167_Karimi_TCF7/data/bam/$b.bam"
  [ -f "$BAM.bai" ] || $SAM index "$BAM"
  d=$(depth 6:129600774-129600804 "$BAM")
  spl=$($SAM view "$BAM" 6:129600774-129600804 2>/dev/null | awk '$6~/N/{c++} END{print c+0}')
  can=$(depth 6:129599500-129599735 "$BAM")   # 201 canonical first exon (positive control)
  echo "  $b: altFirstExon meanDepth=$d  splicedReads=$spl  |  201-canonicalFirstExon meanDepth=$can"
done

echo ""
echo "### (2) GSE147371: 203-vs-204 identifiability (retention vs 204-exclusive intron-3) ###"
for b in GVHD_CD4_Tem_1 GVHD_CD4_Tem_2 GVHD_CD4_Tem_3 GVHD_CD4_Tn; do
  BAM="datasets/GSE147371_GVHD_CD4/data/bam/$b.bam"
  [ -f "$BAM" ] || continue
  [ -f "$BAM.bai" ] || $SAM index "$BAM"
  ret=$(depth 6:129593531-129593631 "$BAM")   # {203,204} retention segment
  i3=$(depth 6:129593737-129594445 "$BAM")     # 204-exclusive intron-3 retention
  e4=$(depth 6:129593632-129593736 "$BAM")     # constitutive exon4
  echo "  $b: retention(203+204)=$ret  intron3(204-only)=$i3  exon4=$e4"
done
echo "  Interp: in Tem, retention > intron3 (excess = 203-type); in Tn both ~0 (induced)."
