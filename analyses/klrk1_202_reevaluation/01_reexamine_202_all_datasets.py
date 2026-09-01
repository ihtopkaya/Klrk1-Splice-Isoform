#!/usr/bin/env python3
# ============================================================================
# reexamine_202_all_datasets.py
# Klrk1-202 (ENSMUST00000095412, protein-coding) was excluded globally because,
# in the discovery dataset, its UNCORRECTED TPM inflated 0->52 (0 with bias
# correction). This script re-checks that decision across ALL datasets (with bias
# correction) and recomputes the key isoform fractions WITH vs WITHOUT 202, to see
# whether the central claims change.
#
# Run: /opt/anaconda3/envs/proteomics/bin/python reexamine_202_all_datasets.py
# ============================================================================
import os, csv, numpy as np
ROOT = "/Users/study/Desktop/Karimi/Klrk1_GVHD_project"
DSDIR = {"GSE147371":"GSE147371_GVHD_CD4","GSE109125":"GSE109125_ImmGen",
         "GSE203167":"GSE203167_Karimi_TCF7","GSE119943":"GSE119943_Yao_LCMV"}
TX = {"ENSMUST00000032252":"201","ENSMUST00000168919":"205","ENSMUST00000152256":"204",
      "ENSMUST00000137660":"203","ENSMUST00000204694":"206","ENSMUST00000095412":"202"}
RI = {"203","204","206"}                       # retained-intron isoforms (202 & 201 & 205 are PC)
FIVE = ["201","205","204","203","206"]         # manuscript set (no 202)
SIX  = ["201","205","202","204","203","206"]   # with 202

META = {"GSE147371":"samples_primary_GSE147371_GSE109125.csv",
        "GSE109125":"samples_primary_GSE147371_GSE109125.csv",
        "GSE203167":"samples_GSE203167.csv",
        "GSE119943":"samples_GSE119943.csv"}

def read_tpm(ds, salmon_dir):
    q = os.path.join(ROOT,"datasets",DSDIR[ds],"data","salmon",salmon_dir,"quant.sf")
    if not os.path.exists(q): return None
    d = {v:0.0 for v in TX.values()}
    for r in csv.DictReader(open(q),delimiter="\t"):
        t=r["Name"].split(".")[0]
        if t in TX: d[TX[t]]=float(r["TPM"])
    return d

def load_conditions(ds):
    """return dict condition -> list of salmon_dir, for this dataset."""
    m = os.path.join(ROOT,"shared","metadata",META[ds]); conds={}
    for r in csv.DictReader(open(m)):
        # primary metadata has a 'dataset' column; others are single-dataset
        rowds = r.get("dataset", ds)
        if rowds != ds: continue
        sd = r["salmon_dir"].split("/")[-1]
        conds.setdefault(r["condition"], []).append(sd)
    return conds

def pct(d, iso_set, isoforms):
    tot = sum(d[i] for i in isoforms)
    return 100*sum(d[i] for i in iso_set)/tot if tot else 0.0

lines=[]
def out(s): print(s); lines.append(s)

out("="*92)
out("Klrk1-202 re-examination: does INCLUDING 202 change the isoform fractions?")
out("="*92)
out(f"{'dataset':10s} {'condition':26s} {'n':>2} | {'203% (5iso)':>11} {'203% (6iso)':>11} | "
    f"{'RI% (5iso)':>10} {'RI% (6iso)':>10} | {'202%':>6}")
out("-"*92)
for ds in DSDIR:
    conds = load_conditions(ds)
    for c, sds in conds.items():
        rows=[read_tpm(ds,sd) for sd in sds]; rows=[r for r in rows if r]
        if not rows: continue
        p203_5=np.mean([pct(r,{"203"},FIVE) for r in rows])
        p203_6=np.mean([pct(r,{"203"},SIX)  for r in rows])
        pri_5 =np.mean([pct(r,RI,FIVE) for r in rows])
        pri_6 =np.mean([pct(r,RI,SIX)  for r in rows])
        p202  =np.mean([pct(r,{"202"},SIX) for r in rows])
        out(f"{ds:10s} {c:26s} {len(rows):>2} | {p203_5:10.2f}  {p203_6:10.2f}  | "
            f"{pri_5:9.2f}  {pri_6:9.2f}  | {p202:5.1f}")
    out("-"*92)

out("\nKEY QUESTIONS:")
out(" * 202 is protein-coding; including it enlarges the denominator, so 203% and RI% shrink.")
out(" * Does the DIRECTION of the effect (203 up with effector/memory differentiation) survive?")
out(" * In discovery datasets (GSE147371/GSE109125) 202 is tiny -> headline barely moves.")
out(" * In GSE203167/GSE119943 202 is large -> fractions shift materially (esp. GSE203167).")
open(os.path.join(ROOT,"datasets","_extra_datasets","reexamine_202_results.txt"),"w").write("\n".join(lines))
print("\nwrote datasets/_extra_datasets/reexamine_202_results.txt")
