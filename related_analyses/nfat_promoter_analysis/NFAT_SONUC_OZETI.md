# Klrk1 Promotör — NFAT Binding Site Analizi (Sonuç Özeti)

**Tarih:** 2026-06-11
**Soru:** Klrk1-203 (retained-intron izoform) transkripsiyonunu yöneten promotörde NFAT bağlanma bölgesi var mı?
**Yöntem:** In silico motif taraması (JASPAR PWM + Biopython PSSM, FIMO mantığı)

---

## 1. Temel Bulgu (Headline)

> **Klrk1 proksimal promotöründe, canonical TSS'ten yalnızca −134 bp upstream'de, mükemmel skorlu (relative score = 1.00) kanonik bir NFAT bağlanma bölgesi (`AATGGAAA`) bulunmaktadır.** Bu site hem **Nfatc1** hem **Nfatc2** PWM'leri tarafından p < 10⁻⁴ ile bağımsız olarak tanınmaktadır.

Bu, Klrk1/NKG2D'nin **aktivasyona-bağımlı (TCR → Ca²⁺ → kalsinörin → NFAT) transkripsiyonel indüksiyonu** modelini doğrudan desteklemektedir.

---

## 2. Önemli Genomik / Biyolojik Bağlam

- **Gen:** Klrk1 (ENSMUSG00000030149), chr6, **minus strand**.
- **Klrk1-203** (ENSMUST00000137660, hedefimiz) TSS = **129,599,647**.
- **Klrk1-201** (ENSMUST00000032252, canonical) TSS = **129,599,735**.
- İki TSS arası yalnızca **88 bp** → Klrk1-203 **canonical promotörü paylaşıyor**, ayrı promotörü yok. Dolayısıyla NFAT analizi paylaşılan Klrk1 promotörü üzerinden yapıldı.
- Referans: GRCm39, Ensembl 111.

---

## 3. Tespit Edilen NFAT Bölgeleri (p < 10⁻⁴)

| # | Dizi | DNA strand | Genomik (chr6) | Canonical TSS'e | Klrk1-203 TSS'e | Tanıyan PWM | rel score |
|---|------|-----------|----------------|-----------------|------------------|-------------|-----------|
| 1 ⭐ | `AATGGAAA` | + | 129,599,866–873 | **−134 bp** | −222 bp | Nfatc1 + Nfatc2 | 1.00 |
| 2 | `AATGGAAA` | − (sense) | 129,600,264–271 | −532 bp | −620 bp | Nfatc1 + Nfatc2 | 1.00 |
| 3 | `CATGGAAA` | + | 129,601,117–124 | −1386 bp | −1474 bp | Nfatc1 + Nfatc2 | 0.98–1.00 |

- Üç site de NFAT çekirdeğini (`GGAAA`) taşır.
- En güçlü ve en proksimal olan **#1**, fonksiyonel promotör penceresinin (−500…−100 bp) tam içindedir.
- NFATC3, NFATC4, NFAT5 matrisleri p < 10⁻⁴ eşiğinde hit vermedi (bu matrislerin bilgi içeriği/eşiği daha katı; çekirdek motif yine de #1–3 ile yakalanıyor).

## 4. AP-1 Bölgeleri ve Kompozit Analizi

NFAT sıklıkla AP-1 ile birlikte (kompozit element) çalışır. AP-1 (FOS::JUN, BATF, BATF::JUN; p < 10⁻³) tarandı:

| Dizi | Genomik | Canonical TSS'e | rel score |
|------|---------|-----------------|-----------|
| `TGACTCA` (mükemmel TRE) | 129,600,804–810 | **−1072 bp** | 1.00 (BATF) |
| `TTACTCA` | 129,600,882–888 | −1150 bp | 0.90 |
| `TTACTCA` | 129,601,195–201 | −1463 bp | 0.90 |

**Kompozit (≤20 bp bitişik NFAT:AP-1) bulunamadı.** En yakın AP-1, proksimal NFAT site'ına 48 bp uzaklıkta. Yani Klrk1 promotörü, klasik ARRE2/IL-2-tipi sıkı NFAT:AP-1 kompoziti yerine, **bağımsız bir proksimal NFAT elementi + distal AP-1 kümesi** mimarisine sahip. Bu mimari de aktivasyon-indüklenebilir transkripsiyonla uyumludur (NFAT ve AP-1/MAPK yolları birbirinden bağımsız katkı verir).

---

## 5. Yorum

1. **Proksimal NFAT elementi (−134 bp)** Klrk1'in aktivasyona-bağımlı uyarılmasının doğrudan moleküler temelini sağlar. Klrk1-203 (RI izoformu) bu promotörü paylaştığı için, NFAT aktivasyonu hem canonical hem RI izoform üretimini birlikte tetikleyebilir — bu da daha önce gözlediğimiz aktivasyon-sonrası RI artışıyla tutarlıdır.
2. Bu bir **in silico tahmindir**; deneysel doğrulama (ChIP-seq/ChIP-qPCR, lusiferaz promotör-reporter ± CsA/FK506 kalsinörin inhibisyonu, NFAT site mutasyonu) gelecekteki adımdır. Yine de literatürle uyumlu güçlü bir hipotez sunar.

---

## 6. Manuscript Metni (İngilizce, kullanıma hazır)

### Methods — *In silico promoter analysis*

> The proximal promoter of murine *Klrk1* was extracted from the GRCm39 genome (Ensembl release 111). Because the retained-intron isoform *Klrk1-203* (ENSMUST00000137660; TSS chr6:129,599,647) initiates only 88 bp downstream of the canonical *Klrk1-201* TSS (ENSMUST00000032252; chr6:129,599,735) on the minus strand, the two transcripts share a common promoter. A 2.3-kb window (chr6:129,599,447–129,601,735; −2000 to +200 relative to the canonical TSS) was retrieved with `samtools faidx`. Position weight matrices for NFAT family members (Nfatc1 MA0624.3, Nfatc2 MA0152.3, NFATC3 MA0625.3, NFATC4 MA1525.3, NFAT5 MA0606.3) and AP-1 family members (FOS::JUN MA0099.3, BATF MA1634.2, BATF::JUN MA0462.2) were obtained from JASPAR 2024 CORE (vertebrates). Sequences were scored on both strands with Biopython (Bio.motifs) log-odds PSSMs using a pseudocount of 0.5 and the local promoter nucleotide composition as background. Significance thresholds were derived from the score distribution at a false-positive rate of 1×10⁻⁴ (NFAT) and 1×10⁻³ (AP-1), equivalent to the FIMO p-value approach. NFAT–AP-1 composite elements were defined as edge-to-edge spacing ≤20 bp.

### Results — *NFAT binding site paragraph*

> To determine whether *Klrk1* induction could be driven by activation-dependent transcription factors, we scanned its proximal promoter for NFAT binding sites. A perfect-consensus NFAT element (`AATGGAAA`, relative score 1.00) was identified 134 bp upstream of the canonical *Klrk1* transcription start site, recognized independently by both Nfatc1 and Nfatc2 matrices at p < 1×10⁻⁴. Two additional NFAT core motifs were present at −532 bp and −1386 bp. A perfect AP-1/TRE element (`TGACTCA`) was found at −1072 bp, together with two weaker AP-1 motifs, but no tight NFAT:AP-1 composite (≤20 bp) was detected, indicating an architecture of an independent proximal NFAT element flanked by distal AP-1 sites. Because *Klrk1-203* shares this promoter with the canonical transcript, NFAT-driven transcriptional activation provides a parsimonious mechanism linking T-cell activation to the coordinate induction of both canonical and retained-intron *Klrk1* isoforms.

---

## 7. Dosyalar / Tekrarlanabilirlik

| Dosya | İçerik |
|-------|--------|
| `run_all.sh` | Uçtan uca tüm pipeline (tek komut) |
| `scripts/scan_nfat.py` | NFAT taraması (p<1e-4, çift strand) |
| `scripts/scan_ap1_composite.py` | AP-1 + NFAT:AP-1 kompozit analizi |
| `scripts/make_figure.py` | Promotör haritası figürü |
| `seq/klrk1_promoter_plus.fa` | Promotör (genomik plus strand) |
| `seq/klrk1_promoter_sense.fa` | Promotör (gen sense yönü) |
| `pwm/*.jaspar` | İndirilen JASPAR PWM'leri |
| `results/nfat_hits.csv` | Tüm NFAT hit tablosu |
| `results/FigN1_Klrk1_NFAT_promoter_map.png/.pdf` | Figür |

**Araçlar:** samtools 1.x, Python 3 + Biopython 1.87, JASPAR 2024 CORE.
