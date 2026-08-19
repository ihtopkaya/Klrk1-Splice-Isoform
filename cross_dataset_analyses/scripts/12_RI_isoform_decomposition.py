#!/usr/bin/env python3
# ============================================================================
# 12_RI_isoform_decomposition.py
#
# PURPOSE
#   Decompose the Klrk1 retained-intron (RI) fraction into its individual
#   isoforms (Klrk1-203 / 204 / 206) for every condition across all four
#   datasets, so the narrative can separate the characterized NKG2D-TR analog
#   (Klrk1-203) from the total RI fraction. Shows that:
#     - naive populations have a HIGH total RI fraction but 0% Klrk1-203
#       (their RI is entirely Klrk1-204/206 background);
#     - Klrk1-203 specifically is the activation/differentiation-induced isoform.
#
# INPUT  : Salmon quant.sf per sample (BASE_DIR/salmon_output/...)
# OUTPUT : results/tables/Klrk1_203_vs_204_206_decomposition.csv + console table
# RUN    : python3 scripts/12_RI_isoform_decomposition.py   (no pysam needed)
# ============================================================================

import csv, os, statistics as st

import sys
sys.path.insert(0, os.path.join(
    os.environ.get("KLRK1_ROOT", "/Users/study/Desktop/Karimi/Klrk1_GVHD_project"),
    "shared", "scripts_common"))
import config as cfg

TBL   = shared_out = os.path.join(cfg.PROJECT_ROOT, "shared", "cross_dataset_analyses", "results")
os.makedirs(TBL, exist_ok=True)
KLRK1 = cfg.KLRK1

# dataset -> condition -> list of salmon_output sample dir names
GROUPS = {
 "GSE109125": {
   "Healthy_CD4_Naive":      ["Healthy_CD4Nve_1","Healthy_CD4Nve_2","Healthy_CD4Nve_3","Healthy_CD4Nve_4"],
   "Healthy_CD8_Naive_6wk":  ["Healthy_CD8Nve_1","Healthy_CD8Nve_2"],
   "Healthy_CD8_Naive_7wk":  ["Healthy_CD8Nve_3","Healthy_CD8Nve_4"],
   "CD8_MPEC":               ["CD8_MPEC_1","CD8_MPEC_2"],
   "CD8_Tcm":                ["CD8_Tcm_1","CD8_Tcm_2"],
   "CD8_Effector_SLEC":      ["CD8_Effector_1"],
 },
 "GSE147371": {
   "GVHD_CD4_Tn":  ["Tn_GSM4427970_28_29"],
   "GVHD_CD4_Tem": ["Tem_GSM4427964_22_23","Tem_GSM4427966_24_25","Tem_GSM4427968_26_27"],
 },
 "GSE203167": {
   "WT_Pre":        ["WT_Pre_1","WT_Pre_2","WT_Pre_3"],
   "WT_Post7":      ["WT_Post7_1","WT_Post7_2","WT_Post7_3"],
   "TCF7cKO_Pre":   ["TCF7cKO_Pre_1","TCF7cKO_Pre_2","TCF7cKO_Pre_3"],
   "TCF7cKO_Post7": ["TCF7cKO_Post7_1","TCF7cKO_Post7_2","TCF7cKO_Post7_3"],
 },
 "GSE119943": {
   "Arm_D4.5_EarlyEffector": ["GSE119943/Arm_D4.5_EarlyEff_1","GSE119943/Arm_D4.5_EarlyEff_2","GSE119943/Arm_D4.5_EarlyEff_3"],
   "Arm_D7_MPEC":            ["GSE119943/Arm_D7_MPEC_1","GSE119943/Arm_D7_MPEC_2","GSE119943/Arm_D7_MPEC_3"],
   "Arm_D7_SLEC":            ["GSE119943/Arm_D7_SLEC_1","GSE119943/Arm_D7_SLEC_2","GSE119943/Arm_D7_SLEC_3"],
   "Cl13_D7_Progenitor":     ["GSE119943/Cl13_D7_Prog_1","GSE119943/Cl13_D7_Prog_2"],
   "Cl13_D7_TermExhausted":  ["GSE119943/Cl13_D7_TermEx_1","GSE119943/Cl13_D7_TermEx_2"],
   "pMIG_Cl13_Progenitor":   ["GSE119943/pMIG_Cl13_Prog_1","GSE119943/pMIG_Cl13_Prog_2","GSE119943/pMIG_Cl13_Prog_3"],
   "pMIG_Cl13_TermExhausted":["GSE119943/pMIG_Cl13_TermEx_1","GSE119943/pMIG_Cl13_TermEx_2","GSE119943/pMIG_Cl13_TermEx_3"],
 },
}

def read_quant(s, ds):
    sample = s.split("/")[-1]          # tolerate legacy "GSE119943/<sample>" names
    p = cfg.salmon_quant(sample, ds)
    if not os.path.exists(p): return None
    d = {}
    with open(p) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            t = r["Name"].split(".")[0]
            if t in KLRK1: d[KLRK1[t]] = float(r["TPM"])
    return d

out = [["dataset","condition","pct_203","pct_204","pct_206","pct_allRI","share_203_of_RI"]]
print(f"{'dataset':10s} {'condition':24s} {'203%':>6} {'204%':>6} {'206%':>6} {'allRI%':>7} {'203/RI':>7}")
for ds, conds in GROUPS.items():
    for cond, samps in conds.items():
        a = {"203":[],"204":[],"206":[]}
        for s in samps:
            d = read_quant(s, ds)
            if not d: continue
            tot = sum(d.values())
            if tot == 0: continue
            for k in a: a[k].append(100*d.get(k,0)/tot)
        if not a["203"]: continue
        p203, p204, p206 = (round(st.mean(a[k]),1) for k in ("203","204","206"))
        pri = round(p203+p204+p206,1)
        share = round(100*p203/pri) if pri>0 else 0
        out.append([ds,cond,p203,p204,p206,pri,share])
        print(f"{ds:10s} {cond:24s} {p203:6} {p204:6} {p206:6} {pri:7} {share:6}%")

with open(os.path.join(TBL,"Klrk1_203_vs_204_206_decomposition.csv"),"w",newline="") as f:
    csv.writer(f).writerows(out)
print("\nKey: naive RI is high but 0% Klrk1-203 (all 204/206); Klrk1-203 specifically is the")
print("activation/differentiation-induced, NKG2D-TR-analog isoform. Saved decomposition CSV.")
