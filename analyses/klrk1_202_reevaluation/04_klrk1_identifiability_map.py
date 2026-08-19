#!/usr/bin/env python3
# ============================================================================
# 13_klrk1_identifiability_map.py
#
# PURPOSE
#   Determine, from the real GTF, which region of the Klrk1 locus belongs to
#   which isoform(s) — the basis for naming reported quantities by the
#   isoform(s) the reads actually support (exclusive region -> that isoform;
#   shared region -> the named set). Produces:
#     (1) the exonic REGION PARTITION (blocks of constant isoform membership),
#     (2) the SPLICE-JUNCTION inventory (introns labelled by which isoforms
#         splice them), and
#     (3) per-isoform IDENTIFIABILITY (exclusive exon bp + exclusive junctions),
#         which yields the REPORT_UNITS convention in config.py.
#
#   Result (why config.REPORT_UNITS is what it is):
#     205 / 206 / 204 have large exclusive exonic regions -> individually named.
#     203 has 0 exclusive bp but its {203,204} retention segment plus the
#       203-type intron-3 splice junction (absent from 204) make it read-
#       separable -> individually named (footnoted).
#     201 has 0 exclusive bp but is separated by its canonical first exon.
#     202 has NEITHER an exclusive exon base NOR an exclusive junction (202 is a
#       base-subset of 205, joined to 205 by the {202,205}-only alt-promoter
#       junction) -> NOT independently measurable -> merged into "Klrk1-205/202".
#
# INPUT  : shared/ref/genome/Klrk1.gtf
# OUTPUT : python_pipeline/outputs/tables/Klrk1_region_partition.csv
#          python_pipeline/outputs/tables/Klrk1_junction_inventory.csv
#          python_pipeline/outputs/tables/Klrk1_identifiability_units.csv
# RUN    : python3 python_pipeline/scripts/13_klrk1_identifiability_map.py
# ============================================================================
import os, sys, csv, collections
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config as cfg

GTF = os.path.join(cfg.PROJECT_ROOT, "shared", "ref", "genome", "Klrk1.gtf")
NAME = cfg.SHORT                       # ENSMUST.. -> "201".."206"
ALL = ["201", "202", "203", "204", "205", "206"]

# ---- parse exons per transcript --------------------------------------------
ex = collections.defaultdict(list)
with open(GTF) as fh:
    for line in fh:
        if "\texon\t" not in line:
            continue
        f = line.split("\t")
        tid = None
        for tok in f[8].split(";"):
            tok = tok.strip()
            if tok.startswith("transcript_id"):
                tid = tok.split('"')[1].split(".")[0]
        if tid in NAME:
            ex[NAME[tid]].append((int(f[3]), int(f[4])))

# ---- (1) exonic region partition (blocks of constant membership) -----------
memb = collections.defaultdict(set)
for t in ALL:
    for a, b in ex[t]:
        for p in range(a, b + 1):
            memb[p].add(t)
blocks = []
for p in sorted(memb):
    sig = frozenset(memb[p])
    if blocks and p == blocks[-1][1] + 1 and blocks[-1][2] == sig:
        blocks[-1][1] = p
    else:
        blocks.append([p, p, sig])

# ---- (2) splice-junction inventory (introns labelled by member isoforms) ----
jm = collections.defaultdict(set)
for t in ALL:
    exs = sorted(ex[t])
    for i in range(len(exs) - 1):
        jm[(exs[i][1] + 1, exs[i + 1][0] - 1)].add(t)

# ---- (3) per-isoform identifiability ---------------------------------------
excl_bp = {t: sum(b - a + 1 for a, b, s in blocks if s == frozenset({t})) for t in ALL}
excl_j = {t: [(a, b) for (a, b), m in jm.items() if m == frozenset({t})] for t in ALL}

# ---- write -----------------------------------------------------------------
TBL = cfg.TABLE_DIR
with open(os.path.join(TBL, "Klrk1_region_partition.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["start", "end", "length_bp", "members", "type"])
    for a, b, sig in blocks:
        typ = "exclusive:" + next(iter(sig)) if len(sig) == 1 else "shared"
        w.writerow([a, b, b - a + 1, "|".join(sorted(sig)), typ])
with open(os.path.join(TBL, "Klrk1_junction_inventory.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["intron_start", "intron_end", "length_bp", "members", "type"])
    for (a, b) in sorted(jm):
        m = jm[(a, b)]
        typ = "exclusive:" + next(iter(m)) if len(m) == 1 else "shared"
        w.writerow([a, b, b - a + 1, "|".join(sorted(m)), typ])
# basis for the reporting decision per isoform (why individual vs merged)
BASIS = {
    "201": "no self-exclusive feature, but separable by its canonical first exon (129599545-735) -> individual",
    "202": "NO self-exclusive exon base or junction; base-subset of 205; no read separates 202 from 205 -> merged into Klrk1-205/202",
    "203": "no self-exclusive feature, but {203,204} retention (129593531-631) exceeds 204 + 203-type intron-3 splice junction absent from 204 -> individual",
    "204": "exclusive intron-3 retention (129593737-129594445, 709bp) -> individual",
    "205": "exclusive 3' tail (2113bp) + exclusive exon5 acceptor junction -> individual",
    "206": "exclusive giant retained-intron exon (2865bp) -> individual",
}
with open(os.path.join(TBL, "Klrk1_identifiability_units.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["isoform", "exclusive_exon_bp", "exclusive_junctions",
                "self_exclusive_feature", "report_unit", "basis"])
    unit_of = {tx: u for u, txs in cfg.REPORT_UNITS.items() for tx in txs}
    for t in ALL:
        js = ";".join(f"{a}-{b}" for a, b in excl_j[t]) or "NONE"
        has_self = "yes" if (excl_bp[t] > 0 or excl_j[t]) else "no"
        w.writerow([f"Klrk1-{t}", excl_bp[t], js, has_self, unit_of.get(t, "?"), BASIS[t]])

# ---- console ---------------------------------------------------------------
print("REGION PARTITION (exonic blocks by isoform membership):")
for a, b, sig in blocks:
    tag = "  <-- EXCLUSIVE " + next(iter(sig)) if len(sig) == 1 else ""
    print(f"  {a}-{b}  {b-a+1:>5}bp  {{{','.join(sorted(sig))}}}{tag}")
print("\nIDENTIFIABILITY -> reporting units:")
for t in ALL:
    js = ";".join(f"{a}-{b}" for a, b in excl_j[t]) or "NONE"
    self_feat = "self-excl" if (excl_bp[t] > 0 or excl_j[t]) else "no self-excl"
    print(f"  Klrk1-{t}: excl_bp={excl_bp[t]:>4}  excl_junc={js:<20}  {self_feat:<12} {BASIS[t]}")
print("\nREPORT_UNITS (config.py):")
for u, txs in cfg.REPORT_UNITS.items():
    print(f"  {u} = {'+'.join('Klrk1-'+x for x in txs)}")
print(f"\nTables -> {TBL}")
