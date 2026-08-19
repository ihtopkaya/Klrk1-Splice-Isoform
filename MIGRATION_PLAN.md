# 📦 Taşıma Planı — Klrk1_GVHD_project

Hedef iskelet `Klrk1_GVHD_project/` altında kuruldu (boş). Aşağıdaki adımlarla
dosyaları oraya taşı. Path'ler tam (absolute). İstersen Finder'da yap, istersen
komutları Terminal'e yapıştır.

**Kısaltmalar (Terminal'de kullanmak istersen bu 3 satırı önce çalıştır):**
```bash
NEW="/Users/study/Desktop/Karimi/Klrk1_GVHD_project"
HOME_D="/Users/study/Desktop/Karimi/HOME GVHD_splicing_project"
CS="/Users/study/Desktop/Karimi/Zaten burda olan dosya GVHD Splicing Isoforms Study /Claude Study area/GVHD Splicing Isoforms Study Study Area"
ZATEN="/Users/study/Desktop/Karimi/Zaten burda olan dosya GVHD Splicing Isoforms Study "
```

---

## ⚠️ ÖNCE OKU
- **Faz 2'yi (star_analysis) junction job'ı bitmeden YAPMA** — job şu an
  `HOME/star_analysis/bam/`'e BAM yazıyor ve `hisat_index`'i okuyor. Taşırsan bozarsın.
- Silme (Faz 3) en son; iCloud klasörü senkronsa cihazdan da siler — sen doğrula.
- Kod reposu (`Klrk1_NKG2D_isoform_analysis/`) yerinde kalıyor; taşıma sonrası
  script path'lerini güncelleriz.

---

## FAZ 1 — ŞİMDİ taşınabilir (job bunlara dokunmuyor)

### 1a) Manuscript'ler → `manuscript/`
```bash
# güncel ver15'ler + ver13/ver14 → manuscript/
mv "$CS"/ver15*.docx "$CS"/ver15*.pdf                 "$NEW/manuscript/"
mv "$CS"/ver13*.docx "$CS"/ver14*.docx "$CS"/ver14*.md "$NEW/manuscript/old_versions/"
# Claude Study area içindeki 'old versions/' (ver10–ver12) → old_versions/
mv "$CS/old versions/"*                                 "$NEW/manuscript/old_versions/"
```

### 1b) Literatür PDF'leri → `literature/`
```bash
mv "$CS/Articles : Kaynaklar/"*                         "$NEW/literature/"
```

### 1c) NFAT promoter analizi → `related_analyses/` (sadece Claude Study area'da var!)
```bash
mv "$CS/nfat_promoter_analysis/"*                       "$NEW/related_analyses/nfat_promoter_analysis/"
```

### 1d) Proteomik analiz (211M) → `related_analyses/`
```bash
mv "$ZATEN/Klrk1-203 Proteomic Analysis/"*             "$NEW/related_analyses/proteomics_Klrk1_203/"
```

### 1e) Salmon çıktıları → dataset başına `data/<GSE>/salmon/`
```bash
cd "$HOME_D/salmon_output"
mv CD8_* Healthy_*        "$NEW/data/GSE109125_ImmGen/salmon/"       # 15 örnek (ImmGen)
mv Tem_* Tn_*             "$NEW/data/GSE147371_GVHD_CD4/salmon/"     # 4 örnek (GVHD CD4)
mv WT_* TCF7cKO_*         "$NEW/data/GSE203167_Karimi_TCF7/salmon/"  # 12 örnek (Karimi)
mv GSE288143_*            "$NEW/data/_extra_datasets/GSE288143/salmon/"  # 8 (makale dışı)
mv GSE83978_*             "$NEW/data/_extra_datasets/GSE83978/salmon/"   # 5 (makale dışı)
```

### 1f) Salmon indexleri → `data/ref/salmon_index/`
```bash
mv "$HOME_D/salmon_index/"*                             "$NEW/data/ref/salmon_index/"
```

### 1g) ImmGen ham FASTQ'ları (HOME kökündeki gevşek dosyalar) → `data/GSE109125_ImmGen/fastq/`
```bash
cd "$HOME_D"
mv SRR6467037_*.fastq.gz SRR6467038_*.fastq.gz SRR6467039_*.fastq.gz \
   SRR6467040_*.fastq.gz SRR6467045_*.fastq.gz          "$NEW/data/GSE109125_ImmGen/fastq/"
# sıkıştırılmamış tekrarları SİL (aynısının .gz'i taşındı — yer kazan):
rm -f SRR6467039_1.fastq SRR6467039_2.fastq SRR6467045_1.fastq SRR6467045_2.fastq
```

### 1h) Eski analiz/provenance klasörleri → `logs_provenance/` (silme, incele)
```bash
mv "$HOME_D/R_analysis"          "$NEW/logs_provenance/"
mv "$HOME_D/cpat_analysis"       "$NEW/logs_provenance/"
mv "$HOME_D/fastqc_results"      "$NEW/logs_provenance/"
mv "$HOME_D/validation_analysis" "$NEW/logs_provenance/"
mv "$HOME_D/cd8_klrk1_data.csv" "$HOME_D/salmon_GSE203167.log" "$NEW/logs_provenance/"
```

---

## FAZ 2 — junction job BİTTİKTEN sonra (star_analysis)

### 2a) Genom-hizalı BAM'ler → dataset başına `bam/`
```bash
cd "$HOME_D/star_analysis/bam"
mv WT_Pre.bam* WT_Post7.bam*                            "$NEW/data/GSE203167_Karimi_TCF7/bam/"
mv SRR6467039_Tcm.bam* SRR6467045_Naive.bam* ImmGen_CD8_Effector.bam* \
                                                        "$NEW/data/GSE109125_ImmGen/bam/"
mv GVHD_CD4_*.bam*                                      "$NEW/data/GSE147371_GVHD_CD4/bam/"
# hisat log'ları da provenance için:
mv *_hisat.log *.tsv                                    "$NEW/logs_provenance/"
```

### 2b) Referans + indexler → `data/ref/`
```bash
mv "$HOME_D/star_analysis/hisat_index/"*  "$NEW/data/ref/hisat2_index/"
mv "$HOME_D/star_analysis/index/"*        "$NEW/data/ref/star_index/"
mv "$HOME_D/star_analysis/ref/"*          "$NEW/data/ref/genome/"
```

---

## FAZ 3 — tekrar kopyaları SİL (her şey taşındıktan ve doğrulandıktan sonra)

Bunlar kanıtlanmış birebir tekrar (md5 eşleşti) — özgün içerik yok:
```bash
# 1) iCloud yedeği (Claude Study area'nın eski alt-kümesi) — ⚠️ iCloud senkronunu doğrula
rm -rf "/Users/study/Desktop/Karimi/Karimi lab icloud directory"

# 2) Claude Study area içindeki VERİ kopyası (HOME'da zaten var) — manuscript/nfat zaten Faz 1'de alındı
rm -rf "$CS/Analysis"

# 3) boşalan HOME kabuğu (her şey taşınınca) — önce içini kontrol et
#    kalan: 6th revision reanalysis, raw_data → incele, gerekirse logs_provenance'a taşı
```

---

## Sonuç: nihai yapı
```
Klrk1_GVHD_project/
├── data/{ref, GSE147371_GVHD_CD4, GSE109125_ImmGen, GSE203167_Karimi_TCF7,
│         GSE119943_Yao_LCMV(boş, diğer Mac'te), _extra_datasets}
├── manuscript/{current ver15, old_versions}
├── related_analyses/{nfat_promoter_analysis, proteomics_Klrk1_203}
├── literature/
└── logs_provenance/
Klrk1_NKG2D_isoform_analysis/   ← kod reposu (yerinde, GitHub)
```
