#!/usr/bin/env bash
# ============================================================================
# 14_align_junction.sh  —  align one sample to the chr6 HISAT2 index for
#   junction / intron-4 counting, then delete the FASTQ to conserve disk.
#
# Extends alignment-level (junction) validation to the two previously
# Salmon-only datasets: GSE147371 (76 bp) and GSE109125/ImmGen (25 bp).
#
# USAGE
#   bash 14_align_junction.sh <LABEL> <READLEN> <SRR1[,SRR2,...]> <DATASET>
#     LABEL    output BAM name (-> datasets/<DATASET>/data/bam/<LABEL>.bam)
#     READLEN  25 | 76  (selects HISAT2 sensitivity params)
#     SRRs     one or more run accessions; multiple runs of the same GSM are
#              concatenated (comma-separated, no spaces)
#     DATASET  GEO accession (GSE147371 | GSE109125 | ...) -> chooses output folder
#
# Requires star_env on PATH (hisat2, samtools) and sra-tools (prefetch,
# fasterq-dump).  Disk-safe: prefetch + fastq land in SCRATCH and are removed
# after the BAM is written.
# ============================================================================
set -euo pipefail

LABEL="${1:?need LABEL}"
READLEN="${2:?need READLEN}"
SRRS="${3:?need SRR list}"
DATASET="${4:?need DATASET (e.g. GSE147371 or GSE109125)}"

# paths from shared config (set KLRK1_ROOT on a new machine)
KLRK1_ROOT="${KLRK1_ROOT:-/Users/study/Desktop/Karimi/Klrk1_GVHD_project}"
source "${KLRK1_ROOT}/shared/scripts_common/config.sh"

BAMDIR="$(klrk1_data_dir "$DATASET" bam)"   # per-dataset: datasets/<DS>/data/bam
SCRATCH="$SCRATCH_DIR"
IDX="$HISAT2_INDEX"                          # shared/ref/hisat2_index/chr6 (space-free)
mkdir -p "$BAMDIR" "$SCRATCH"

# read-length-specific HISAT2 sensitivity (25 bp needs relaxed params, matching
# the provenance of the existing ImmGen BAMs: -L 20 --score-min L,0,-0.6).
# Use the ${arr[@]+"${arr[@]}"} idiom below so an empty array is safe under
# `set -u` on macOS bash 3.2 (a plain "${HS_EXTRA[@]}" errors as "unbound").
if [ "$READLEN" -le 30 ]; then
  HS_EXTRA=(-L 20 --score-min L,0,-0.6)
else
  HS_EXTRA=()
fi

# Disk-safe: process one run at a time (download -> align -> delete fastq),
# producing a per-run BAM, then merge per-GSM. Peak disk = one run's fastqs.
PART_BAMS=()
IFS=',' read -ra RUNS <<< "$SRRS"
for SRR in "${RUNS[@]}"; do
  echo "[$(date +%H:%M:%S)] $LABEL : fetching $SRR"
  if [ ! -f "$SCRATCH/$SRR/$SRR.sra" ]; then
    prefetch "$SRR" --max-size 40g -O "$SCRATCH" >/dev/null 2>&1
  fi
  # Prefer fasterq-dump (multithreaded, ~5-10x faster) when the disk has room for
  # its large temp footprint; otherwise fall back to fastq-dump --gzip, which has
  # a low disk peak but is slow. HISAT2 reads either plain or gzipped FASTQ.
  if ! fasterq-dump "$SCRATCH/$SRR/$SRR.sra" -O "$SCRATCH" -e "$THREADS" \
        --split-3 -f -t "$SCRATCH/fqtmp" >/dev/null 2>&1; then
    echo "[$(date +%H:%M:%S)] $LABEL : fasterq-dump failed (disk?), using fastq-dump --gzip"
    rm -rf "$SCRATCH/fqtmp"; rm -f "$SCRATCH/${SRR}"*.fastq
    fastq-dump --split-3 --gzip "$SCRATCH/$SRR/$SRR.sra" -O "$SCRATCH" >/dev/null 2>&1
  fi
  rm -rf "$SCRATCH/$SRR" "$SCRATCH/fqtmp"   # drop the .sra + temp immediately
  PART="$BAMDIR/${LABEL}__${SRR}.bam"
  echo "[$(date +%H:%M:%S)] $LABEL : aligning $SRR (readlen=$READLEN)"
  # detect paired vs single and plain vs gzipped output from either dumper
  R1=""; R2=""; U=""
  for e in fastq fastq.gz; do
    [ -f "$SCRATCH/${SRR}_1.$e" ] && R1="$SCRATCH/${SRR}_1.$e"
    [ -f "$SCRATCH/${SRR}_2.$e" ] && R2="$SCRATCH/${SRR}_2.$e"
    [ -f "$SCRATCH/${SRR}.$e" ]   && U="$SCRATCH/${SRR}.$e"
  done
  if [ -n "$R1" ] && [ -n "$R2" ]; then
    hisat2 -x "$IDX" -p "$THREADS" --dta ${HS_EXTRA[@]+"${HS_EXTRA[@]}"} \
      -1 "$R1" -2 "$R2" \
      2> "$BAMDIR/${LABEL}__${SRR}_hisat.log" | samtools sort -@ 2 -o "$PART" -
  elif [ -n "$U" ]; then
    hisat2 -x "$IDX" -p "$THREADS" --dta ${HS_EXTRA[@]+"${HS_EXTRA[@]}"} \
      -U "$U" \
      2> "$BAMDIR/${LABEL}__${SRR}_hisat.log" | samtools sort -@ 2 -o "$PART" -
  else
    echo "[err] no fastq produced for $SRR"; exit 1
  fi
  rm -f "$SCRATCH/${SRR}"*.fastq "$SCRATCH/${SRR}"*.fastq.gz
  PART_BAMS+=("$PART")
done

if [ "${#PART_BAMS[@]}" -eq 1 ]; then
  mv "${PART_BAMS[0]}" "$BAMDIR/${LABEL}.bam"
else
  samtools merge -f "$BAMDIR/${LABEL}.bam" "${PART_BAMS[@]}"
  rm -f "${PART_BAMS[@]}" "${PART_BAMS[@]/%.bam/.bam.bai}"
fi
samtools index "$BAMDIR/${LABEL}.bam"

echo "[$(date +%H:%M:%S)] $LABEL : done -> $BAMDIR/${LABEL}.bam"
echo "  per-run alignment rates:"
for L in "$BAMDIR/${LABEL}__"*_hisat.log; do
  [ -f "$L" ] && echo "    $(basename "$L"): $(grep 'overall alignment rate' "$L" | tail -1)"
done
echo "  reads over Klrk1 intron4:"; samtools view -c "$BAMDIR/${LABEL}.bam" 6:129593296-129593631
