#!/usr/bin/env bash
# ============================================================================
# config.sh — single source of truth for paths used by the shell pipeline
# scripts (Salmon quant, HISAT2 alignment, sashimi, CPAT/NMD).
#
# REPRODUCIBILITY: set KLRK1_ROOT on a new machine, or edit PROJECT_ROOT below.
# Source it from a script:   source "$(dirname "$0")/../../shared/scripts_common/config.sh"
# (or by absolute path). Helpers below resolve per-dataset data locations.
# ============================================================================

PROJECT_ROOT="${KLRK1_ROOT:-/Users/study/Desktop/Karimi/Klrk1_GVHD_project}"

# reference / indices (large, git-ignored)
REF_DIR="${PROJECT_ROOT}/shared/ref"
SALMON_INDEX_DIR="${REF_DIR}/salmon_index"
HISAT2_INDEX="${REF_DIR}/hisat2_index/chr6"
TRANSCRIPTOME_FA="${SALMON_INDEX_DIR}/Mus_musculus.GRCm39.cdna.all.fa.gz"
# read-length-matched Salmon indices
INDEX_K31="${SALMON_INDEX_DIR}/mouse_index_light"      # 76 bp (GSE147371)
INDEX_K21="${SALMON_INDEX_DIR}/mouse_index_light_k21"  # 25 bp (GSE109125)
INDEX_K27="${SALMON_INDEX_DIR}/mouse_index_light_k27"  # 50-51 bp (GSE203167, GSE119943)

SCRATCH_DIR="${KLRK1_SCRATCH:-$HOME/klrk1_junction_scratch}"
SALMON="${SALMON:-/opt/anaconda3/envs/salmon_arm64/bin/salmon}"  # homebrew salmon has a broken libtbb
SALMON_FLAGS="--validateMappings --gcBias --seqBias"
THREADS="${THREADS:-6}"

# map a GEO accession to its dataset folder
klrk1_dataset_folder() {
  case "$1" in
    GSE147371) echo "GSE147371_GVHD_CD4" ;;
    GSE109125) echo "GSE109125_ImmGen" ;;
    GSE203167) echo "GSE203167_Karimi_TCF7" ;;
    GSE119943) echo "GSE119943_Yao_LCMV" ;;
    *) echo "UNKNOWN_$1" ;;
  esac
}
# per-dataset data dirs:  klrk1_data_dir <GSE> <fastq|bam|salmon>
klrk1_data_dir() {
  echo "${PROJECT_ROOT}/datasets/$(klrk1_dataset_folder "$1")/data/$2"
}

export SALMON
export PROJECT_ROOT REF_DIR SALMON_INDEX_DIR HISAT2_INDEX TRANSCRIPTOME_FA \
       INDEX_K31 INDEX_K21 INDEX_K27 SCRATCH_DIR SALMON_FLAGS THREADS
