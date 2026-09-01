#!/usr/bin/env python3
# ============================================================================
# make_figure.py
# Klrk1 promotor NFAT/AP-1 binding-site haritasi
# X ekseni: canonical Klrk1 TSS'e (ENSMUST00000032252, 129,599,735) uzaklik
# ============================================================================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Rectangle
from Bio import motifs
from Bio.Seq import Seq

WINDOW_START = 129599447
TSS201 = 129599735   # canonical (referans 0)
TSS203 = 129599647   # Klrk1-203
PSEUDO = 0.5

def read_fa(p):
    s=[]
    for l in open(p):
        if not l.startswith(">"): s.append(l.strip())
    return "".join(s).upper()

prom = read_fa("seq/klrk1_promoter_plus.fa"); L=len(prom)
bg={b:prom.count(b)/len(prom) for b in "ACGT"}
sobj=Seq(prom)

def scan(files, fpr):
    out=[]
    for tf,path in files.items():
        m=motifs.read(open(path),"jaspar"); m.pseudocounts=PSEUDO; m.background=bg
        ps=m.pssm; w=m.length
        thr=ps.distribution(background=bg,precision=10**4).threshold_fpr(fpr)
        for pos,sc in ps.search(sobj,threshold=thr,both=True):
            st='+' if pos>=0 else '-'; s0=pos if pos>=0 else pos+L
            gs=WINDOW_START+s0; ge=gs+w-1
            rel=(sc-ps.min)/(ps.max-ps.min)
            site=prom[s0:s0+w]
            if st=='-': site=str(Seq(site).reverse_complement())
            out.append(dict(tf=tf,gs=gs,ge=ge,center=(gs+ge)/2,strand=st,
                            score=sc,rel=rel,site=site))
    return out

def dist201(c): return TSS201 - c   # upstream negatif

nfat=scan({"Nfatc1":"pwm/MA0624.3.jaspar","Nfatc2":"pwm/MA0152.3.jaspar"},1e-4)
ap1 =scan({"FOS::JUN":"pwm/MA0099.3.jaspar","BATF":"pwm/MA1634.2.jaspar",
           "BATF::JUN":"pwm/MA0462.2.jaspar"},1e-3)

# AP-1 hit'lerini ortusenlere gore birlestir (en yuksek skoru tut)
def dedupe(hits):
    hits=sorted(hits,key=lambda h:-h["rel"]); kept=[]
    for h in hits:
        if all(not (h["gs"]<=k["ge"] and k["gs"]<=h["ge"]) for k in kept):
            kept.append(h)
    return kept
ap1d=dedupe(ap1)
# NFAT'i da merkeze gore birlestir (ayni site farkli matrisle)
nfatd=[]
for h in sorted(nfat,key=lambda x:-x["rel"]):
    if all(abs(h["center"]-k["center"])>3 for k in nfatd): nfatd.append(h)

# ---- Cizim ----
fig,ax=plt.subplots(figsize=(12,3.6))
xmin,xmax=-2050,200
# promotor ekseni
ax.plot([xmin,xmax],[0,0],color="#555555",lw=2,zorder=1)
# TSS oklari
ax.annotate("",xy=(dist201(TSS201)+60,0),xytext=(dist201(TSS201),0),
            arrowprops=dict(arrowstyle="-|>",color="black",lw=2.2))
ax.text(0,0.32,"Klrk1-201 TSS\n(canonical)",ha="center",va="bottom",fontsize=8,fontweight="bold")
d203=dist201((TSS203))
ax.annotate("",xy=(d203+60,-0.0),xytext=(d203,-0.0),
            arrowprops=dict(arrowstyle="-|>",color="#777777",lw=1.6))
ax.text(d203,-0.42,"Klrk1-203 TSS",ha="center",va="top",fontsize=7,color="#555555")

# NFAT site'lar (kirmizi, ust)
for h in nfatd:
    x=dist201(h["center"])
    ax.add_patch(Rectangle((x-12,0.08),24,0.18,color="#D6604D",zorder=3))
    ax.plot([x,x],[0,0.08],color="#D6604D",lw=1,zorder=2)
    lab=f"NFAT\n{h['site']}\n{x:+.0f}"
    ax.text(x,0.30,lab,ha="center",va="bottom",fontsize=7,color="#B2182B",fontweight="bold")

# AP-1 site'lar (mavi, alt) - sadece rel>=0.90 veya perfect
for h in ap1d:
    if h["rel"]<0.90: continue
    x=dist201(h["center"])
    ax.add_patch(Rectangle((x-12,-0.26),24,0.18,color="#4393C3",zorder=3))
    ax.plot([x,x],[-0.08,0],color="#4393C3",lw=1,zorder=2)
    star=" *" if h["rel"]>=0.999 else ""
    # drop the -1150 label one level so its sequence text doesn't collide with the
    # neighbouring -1072 label (they are only ~78 bp apart); a faint guide line keeps
    # the lowered label visually tied to its peak.
    ylbl=-0.72 if -1160<=x<=-1140 else -0.50
    if ylbl<-0.50:
        ax.plot([x,x],[-0.28,ylbl+0.02],color="#4393C3",lw=0.6,ls=":",zorder=1)
    ax.text(x,ylbl,f"AP-1{star}\n{h['site']}\n{x:+.0f}",ha="center",va="top",
            fontsize=7,color="#2166AC",fontweight="bold")

ax.set_xlim(xmin,xmax); ax.set_ylim(-1.18,0.95)
ax.set_yticks([])
ax.set_xlabel("Distance from canonical Klrk1 TSS (bp)",fontsize=10,fontweight="bold")
ax.set_title("Klrk1 proximal promoter: NFAT and AP-1 binding sites (in silico, JASPAR/Biopython, p<1e-4 NFAT / p<1e-3 AP-1)",
             fontsize=10,fontweight="bold")
for s in ["top","right","left"]: ax.spines[s].set_visible(False)
ax.axvspan(-500,100,color="#FFF2CC",alpha=0.5,zorder=0)  # proksimal promotor
ax.text(-200,0.85,"proximal promoter",fontsize=7,color="#999900",style="italic",ha="center")

plt.tight_layout()
import os
os.makedirs("results",exist_ok=True)
fig.savefig("results/FigN1_Klrk1_NFAT_promoter_map.png",dpi=300,bbox_inches="tight")
fig.savefig("results/FigN1_Klrk1_NFAT_promoter_map.pdf",bbox_inches="tight")
print("Figur kaydedildi: results/FigN1_Klrk1_NFAT_promoter_map.png / .pdf")

# Konsol ozeti
print("\nNFAT site (birlesik):")
for h in sorted(nfatd,key=lambda x:dist201(x["center"])):
    print(f"  {h['site']:9s} TSS201 {dist201(h['center']):+.0f} bp  rel={h['rel']:.2f} strand={h['strand']}")
print("\nAP-1 site (rel>=0.90, birlesik):")
for h in sorted([a for a in ap1d if a['rel']>=0.90],key=lambda x:dist201(x['center'])):
    print(f"  {h['tf']:9s} {h['site']:12s} TSS201 {dist201(h['center']):+.0f} bp rel={h['rel']:.2f}")
