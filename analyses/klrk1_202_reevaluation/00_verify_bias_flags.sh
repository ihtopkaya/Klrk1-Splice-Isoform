#!/bin/bash
# ============================================================================
# 00_verify_bias_flags.sh
# Verify, from each Salmon run's cmd_info.json, whether --gcBias and --seqBias
# were applied. This underpins the whole re-evaluation: the ver16 "202 inflates
# 0->52 without correction" claim only describes the UNCORRECTED discovery run;
# every primary dataset here was quantified WITH bias correction ON.
# Run: bash analyses/klrk1_202_reevaluation/00_verify_bias_flags.sh
# ============================================================================
set -euo pipefail
ROOT="${KLRK1_ROOT:-/Users/study/Desktop/Karimi/Klrk1_GVHD_project}"
cd "$ROOT"

printf "%-52s %-8s %-8s %-10s\n" "SAMPLE_DIR" "gcBias" "seqBias" "salmon_v"
printf "%-52s %-8s %-8s %-10s\n" "----------" "------" "-------" "--------"
find datasets -name "cmd_info.json" | sort | while read -r f; do
  d=$(dirname "$f" | sed 's|datasets/||')
  gc=$(grep -c '"gcBias"' "$f" || true); sb=$(grep -c '"seqBias"' "$f" || true)
  ver=$(grep '"salmon_version"' "$f" | sed 's/.*: *"//; s/".*//')
  gctxt="NO"; sbtxt="NO"
  [ "$gc" -gt 0 ] && gctxt="YES"
  [ "$sb" -gt 0 ] && sbtxt="YES"
  printf "%-52s %-8s %-8s %-10s\n" "$d" "$gctxt" "$sbtxt" "$ver"
done
# Expected: all 4 primary datasets (GSE109125/119943/147371/203167) = YES/YES;
# extra GSE288143 = YES/YES; extra GSE83978 = NO/YES (gcBias off; not comparable).
