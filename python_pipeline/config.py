# ============================================================================
# config.py — paths for the standalone Python pipeline.
# Set KLRK1_ROOT on a new machine, or edit PROJECT_ROOT below. Nothing else
# needs changing.
# ============================================================================
import os

PROJECT_ROOT = os.environ.get(
    "KLRK1_ROOT", "/Users/study/Desktop/Karimi/Klrk1_GVHD_project")

# where this Python pipeline writes its own outputs (kept separate from the R run)
HERE        = os.path.join(PROJECT_ROOT, "python_pipeline")
FIG_DIR     = os.path.join(HERE, "outputs", "figures")
TABLE_DIR   = os.path.join(HERE, "outputs", "tables")
for _d in (FIG_DIR, TABLE_DIR):
    os.makedirs(_d, exist_ok=True)

DATASET_FOLDER = {
    "GSE147371": "GSE147371_GVHD_CD4",
    "GSE109125": "GSE109125_ImmGen",
    "GSE203167": "GSE203167_Karimi_TCF7",
    "GSE119943": "GSE119943_Yao_LCMV",
}

def salmon_quant(sample_dir, dataset):
    return os.path.join(PROJECT_ROOT, "datasets", DATASET_FOLDER[dataset],
                        "data", "salmon", sample_dir, "quant.sf")

def results_dir(dataset):
    d = os.path.join(TABLE_DIR, dataset); os.makedirs(d, exist_ok=True); return d

def bam_dir(dataset):
    return os.path.join(PROJECT_ROOT, "datasets", DATASET_FOLDER[dataset], "data", "bam")

def metadata_csv(name):
    return os.path.join(PROJECT_ROOT, "shared", "metadata", name)

# All SIX annotated Klrk1 transcripts. Klrk1-202 (ENSMUST00000095412) is NO LONGER
# excluded: with bias correction ON it is a real (protein-coding, alt-promoter) signal,
# and it must sit in the total-Klrk1 denominator. See docs/KLRK1_202_ANALYSIS_NOTE.md.
KLRK1 = {
    "ENSMUST00000032252": "Klrk1-201", "ENSMUST00000095412": "Klrk1-202",
    "ENSMUST00000168919": "Klrk1-205", "ENSMUST00000137660": "Klrk1-203",
    "ENSMUST00000152256": "Klrk1-204", "ENSMUST00000204694": "Klrk1-206",
}
BIOTYPE = {  # protein-coding vs retained-intron
    "ENSMUST00000032252": "protein_coding", "ENSMUST00000095412": "protein_coding",
    "ENSMUST00000168919": "protein_coding", "ENSMUST00000137660": "retained_intron",
    "ENSMUST00000152256": "retained_intron", "ENSMUST00000204694": "retained_intron",
}
ISO_203 = "ENSMUST00000137660"   # the primary NKG2D-TR-analog isoform

# ---- Read-identifiability REPORTING UNITS (see scripts/13_klrk1_identifiability_map.py) --
# Name reported quantities by the isoform(s) the reads actually support. From the region+
# junction map: 201/203/204/205/206 each carry a read-supported distinguishing feature and
# are reported individually; Klrk1-202 has NO exclusive exonic base or junction (202 ⊆ 205,
# linked to 205 by the {202,205}-exclusive alt-promoter junction 129,598,206–129,600,773),
# so it cannot be measured on its own and is merged with 205 as "Klrk1-205/202".
REPORT_UNITS = {                       # unit label -> transcript short-codes it sums
    "Klrk1-201":     ["201"],
    "Klrk1-205/202": ["205", "202"],   # merged: 202 not independently identifiable from 205
    "Klrk1-203":     ["203"],
    "Klrk1-204":     ["204"],
    "Klrk1-206":     ["206"],
}
RI_UNITS = ["Klrk1-203", "Klrk1-204", "Klrk1-206"]   # retained-intron numerator units
SHORT = {tx: name.split("-")[-1] for tx, name in KLRK1.items()}   # ENSMUST.. -> "201".."206"

def to_units(tpm_by_short):
    """Collapse a {'201':tpm,...,'206':tpm} dict into the REPORT_UNITS (5 units)."""
    return {u: sum(tpm_by_short.get(s, 0.0) for s in members)
            for u, members in REPORT_UNITS.items()}

# ---- Klrk1 intron-4 coordinates (GRCm39, chr6, minus strand) — junction step -
CHROM = "6"
I4_START, I4_END = 129593296, 129593631
B5, B3 = 129593631, 129593296
TOL = 10
RI203_START, RI203_END = 129593531, 129593631
EX4_START, EX4_END = 129593632, 129593736
