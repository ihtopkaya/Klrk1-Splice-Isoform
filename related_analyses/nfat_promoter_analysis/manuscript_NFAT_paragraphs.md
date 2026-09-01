# Klrk1 NFAT Promoter — Manuscript Paragraphs (ready to paste)

> **Notlar:**
> - `[Mullins et al. ref]` bir placeholder'dır — tam künye (yazarlar, yıl, dergi, cilt, sayfa) eklenmelidir.
> - `Figure X` bizim ürettiğimiz `FigN1_Klrk1_NFAT_promoter_map.pdf/.png` figürüne karşılık gelir; manuscript'teki nihai numara ile değiştirilmelidir.
> - Mesafe konvansiyonu: rapor edilen "131 bp", motifin 5′ kenarının (chr6:129,599,866) canonical TSS'e (129,599,735) uzaklığıdır. Motif merkezi −134.5 bp, uzak kenarı −138 bp'dir; yani element −131…−138 bp aralığında yer alır.

---

## Methods

In silico promoter analysis was performed to assess whether the *Klrk1*
retained-intron isoform *Klrk1-203* (ENSMUST00000137660) could be transcriptionally
regulated by activation-induced transcription factors. Because *Klrk1-203* initiates
on the minus strand only 88 bp downstream of the canonical *Klrk1-201* TSS
(ENSMUST00000032252; chr6:129,599,735, GRCm39/Ensembl 111), the two transcripts
share a common promoter; a 2.3-kb window spanning −2000 to +200 bp relative to the
canonical TSS (chr6:129,599,447–129,601,735) was therefore retrieved with `samtools
faidx`. Position weight matrices for NFAT family members (Nfatc1, MA0624.3; Nfatc2,
MA0152.3; NFATC3, MA0625.3; NFATC4, MA1525.3; NFAT5, MA0606.3) and AP-1 family
members (FOS::JUN, MA0099.3; BATF, MA1634.2; BATF::JUN, MA0462.2) were obtained from
JASPAR 2024 CORE (vertebrates) and scored on both strands using Biopython log-odds
PSSMs (`Bio.motifs`) with a pseudocount of 0.5 and the local promoter nucleotide
composition as background. Match significance was assigned from the score
distribution at a false-positive rate of 1×10⁻⁴ for NFAT and 1×10⁻³ for AP-1
(equivalent to the FIMO p-value approach), and NFAT:AP-1 composite elements were
defined as an edge-to-edge spacing of ≤20 bp.

---

## Results

To examine whether *Klrk1* induction could be driven by T-cell–activation
transcription factors, we scanned the shared *Klrk1* proximal promoter for NFAT
binding sites. A perfect-consensus NFAT element (AATGGAAA; relative score 1.00) was
identified 131 bp upstream of the canonical *Klrk1* TSS (motif 5′ edge;
chr6:129,599,866–129,599,873) and was recognized independently by both Nfatc1 and
Nfatc2 matrices at p < 1×10⁻⁴ (Figure X). Two additional NFAT core motifs (GGAAA)
were present further upstream, at −532 bp and −1386 bp. The promoter also contained a
perfect AP-1/TRE element (TGACTCA) at −1072 bp together with two lower-scoring AP-1
motifs; however, no tightly spaced NFAT:AP-1 composite (≤20 bp) was detected, the
nearest AP-1 site lying 48 bp from the proximal NFAT element. The *Klrk1* promoter
thus harbors an independent proximal NFAT element flanked by distal AP-1 sites, and
because *Klrk1-203* shares this promoter with the canonical transcript, NFAT-dependent
activation is positioned to drive coordinate transcription of both isoforms.

---

## Discussion

The presence of a high-affinity, perfect-consensus NFAT element only 131 bp upstream
of the *Klrk1* TSS (motif 5′ edge; chr6:129,599,866–129,599,873) provides a plausible
molecular link between T-cell receptor signaling and this locus: TCR engagement
elevates intracellular Ca²⁺, activates calcineurin, and drives nuclear NFAT
translocation, which could then transactivate *Klrk1* from this proximal site. Because
the retained-intron isoform *Klrk1-203* initiates from the same promoter as the
canonical transcript, this architecture would couple NFAT activation to the
simultaneous production of both isoforms, offering a candidate mechanism for the
differentiation-coupled increase in intron retention observed here, in which a fraction
of rapidly transcribed *Klrk1* pre-mRNA escapes complete splicing. The accompanying
distal AP-1/TRE elements are consistent with the known cooperation between the
Ca²⁺–calcineurin–NFAT and MAPK–AP-1 pathways downstream of TCR signaling, although the
absence of a classical, tightly spaced NFAT:AP-1 composite suggests that NFAT acts here
largely independently rather than within an ARRE2-type cooperative module. In silico
promoter analysis revealed a perfect-consensus NFAT element 131 bp upstream of the
*Klrk1* TSS (motif 5′ edge; chr6:129,599,866–129,599,873), raising the possibility that
TCR-driven NFAT activation contributes to the differentiation-coupled induction of
*Klrk1-203* observed here — a hypothesis consistent with the demonstrated role of NFAT
in regulating KLR family members in human CD8 T cells [Mullins et al. ref] and
warranting experimental validation.

---

## Future perspectives (2 sentences)

Experimental validation — luciferase promoter-reporter assays with and without
calcineurin inhibition (cyclosporin A or FK506) and site-directed mutagenesis of the
−131 bp NFAT element — will be required to establish whether this site is functionally
responsible for activation-induced *Klrk1* transcription. NFAT chromatin
immunoprecipitation (ChIP-qPCR/ChIP-seq) in resting versus activated CD8⁺ T cells
would further confirm direct, inducible NFAT occupancy at the *Klrk1* promoter in vivo.
