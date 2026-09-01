# ============================================================================
# config.py — single source of truth for paths + Klrk1 coordinates used by the
# Python pipeline scripts (junction counting, PIR reconciliation, RI decomposition).
#
# REPRODUCIBILITY: on a new machine set the KLRK1_ROOT environment variable, or
# edit PROJECT_ROOT below. Everything else is derived.
#
#   import sys, os
#   sys.path.insert(0, os.path.join(os.environ.get("KLRK1_ROOT",
#       "/Users/study/Desktop/Karimi/Klrk1_GVHD_project"), "shared", "scripts_common"))
#   import config as cfg
# ============================================================================

import os

PROJECT_ROOT = os.environ.get(
    "KLRK1_ROOT", "/Users/study/Desktop/Karimi/Klrk1_GVHD_project")

# GEO accession -> dataset folder under datasets/
DATASET_FOLDER = {
    "GSE147371": "GSE147371_GVHD_CD4",
    "GSE109125": "GSE109125_ImmGen",
    "GSE203167": "GSE203167_Karimi_TCF7",
    "GSE119943": "GSE119943_Yao_LCMV",
}

def dataset_dir(ds):   return os.path.join(PROJECT_ROOT, "datasets", DATASET_FOLDER[ds])
def salmon_dir(ds):    return os.path.join(dataset_dir(ds), "data", "salmon")
def salmon_quant(sample, ds): return os.path.join(salmon_dir(ds), sample, "quant.sf")
def bam_dir(ds):
    d = os.path.join(dataset_dir(ds), "data", "bam"); os.makedirs(d, exist_ok=True); return d
def results_dir(ds):
    d = os.path.join(dataset_dir(ds), "results");     os.makedirs(d, exist_ok=True); return d
def figures_dir(ds):
    d = os.path.join(dataset_dir(ds), "figures");     os.makedirs(d, exist_ok=True); return d
def shared_ref(*p):    return os.path.join(PROJECT_ROOT, "shared", "ref", *p)
def metadata_csv(name):return os.path.join(PROJECT_ROOT, "shared", "metadata", name)

# HISAT2 chr6-restricted index prefix (built for the sashimi/junction work)
HISAT2_INDEX = shared_ref("hisat2_index", "chr6")

# ---- Klrk1 intron-4 coordinates (GRCm39, chr6, minus strand) ---------------
CHROM = "6"
I4_START, I4_END = 129593296, 129593631      # intron 4 (1-based, inclusive)
B5, B3 = 129593631, 129593296                # exon/intron boundaries
TOL = 10                                      # boundary tolerance (bp)
RI203_START, RI203_END = 129593531, 129593631 # 203-retained portion of intron 4
EX4_START,  EX4_END    = 129593632, 129593736 # constitutive flanking exon 4

# transcript id -> short Klrk1 isoform name (202/095412 excluded: bias artefact)
KLRK1 = {
    "ENSMUST00000032252": "201", "ENSMUST00000168919": "205",
    "ENSMUST00000152256": "204", "ENSMUST00000137660": "203",
    "ENSMUST00000204694": "206",
}
