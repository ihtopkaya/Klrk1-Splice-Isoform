#!/usr/bin/env Rscript
# ============================================================================
# 14_combined_figures_R.R
# R/ggplot versions of the GSE203167 (Fig4) and GSE119943 (Fig8) 2x2 combined
# figures, in the same style as the other R figures. Reads the ground-truth
# per-sample table (Supplementary_all_sample_6isoform.csv) so the numbers are
# GUARANTEED identical to Table4_units_master.csv. Five reporting units
# (Klrk1-205/202 = 205+202); RI% = mean-of-per-sample-ratios; mean +/- SD.
# Run: Rscript python_pipeline/scripts/14_combined_figures_R.R
# ============================================================================
suppressMessages({library(dplyr); library(tidyr); library(ggplot2); library(gridExtra)})
ROOT <- Sys.getenv("KLRK1_ROOT", "/Users/study/Desktop/Karimi/Klrk1_GVHD_project")
tbl  <- file.path(ROOT, "python_pipeline/outputs/tables")
figd <- file.path(ROOT, "python_pipeline/outputs/figures"); dir.create(figd, showWarnings=FALSE, recursive=TRUE)

d <- read.csv(file.path(tbl, "Supplementary_all_sample_6isoform.csv"), check.names=FALSE)
d <- d %>% mutate(
  U_201 = TPM_201, `U_205/202` = TPM_205 + TPM_202, U_203 = TPM_203, U_204 = TPM_204, U_206 = TPM_206,
  RI = 100*(TPM_203+TPM_204+TPM_206)/total6_TPM,
  p201 = 100*TPM_201/total6_TPM, `p205/202` = 100*(TPM_205+TPM_202)/total6_TPM,
  p203 = 100*TPM_203/total6_TPM, p204 = 100*TPM_204/total6_TPM, p206 = 100*TPM_206/total6_TPM)

unit_cols <- c("Klrk1-201"="#E41A1C","Klrk1-205/202"="#FF7F00","Klrk1-203"="#984EA3",
               "Klrk1-204"="#4DAF4A","Klrk1-206"="#377EB8")
th <- theme_classic(base_size=11) + theme(plot.title=element_text(face="bold", size=11),
        axis.text.x=element_text(angle=30, hjust=1, size=8))
msd <- function(x) c(m=mean(x), s=ifelse(length(x)>1, sd(x), 0))

make_panels <- function(ds, ord, labs_short, title){
  s <- d %>% filter(dataset==ds) %>% mutate(condition=factor(condition, levels=ord))
  agg <- s %>% group_by(condition) %>% summarise(
      tot_m=mean(total6_TPM), tot_s=sd(total6_TPM),
      k203_m=mean(TPM_203),  k203_s=sd(TPM_203),
      ri_m=mean(RI),         ri_s=sd(RI), .groups="drop") %>%
    mutate(across(where(is.numeric), ~replace_na(.,0)), lab=labs_short[as.character(condition)])
  agg$lab <- factor(agg$lab, levels=labs_short[ord])
  bar <- function(m,s,ylab,fill,tag) ggplot(agg, aes(lab, .data[[m]])) +
      geom_col(fill=fill, width=0.7) +
      geom_errorbar(aes(ymin=pmax(0,.data[[m]]-.data[[s]]), ymax=.data[[m]]+.data[[s]]), width=0.25, linewidth=0.4) +
      labs(x=NULL, y=ylab, tag=tag) + th
  pA <- bar("tot_m","tot_s","Total Klrk1 (TPM, mean +/- SD)","#2166AC","A")
  pB <- bar("k203_m","k203_s","Klrk1-203 (TPM, mean +/- SD)","#984EA3","B")
  pC <- bar("ri_m","ri_s","Retained-intron fraction (%, mean +/- SD)","#7d3c98","C")
  prop <- s %>% group_by(condition) %>% summarise(
      `Klrk1-201`=mean(p201), `Klrk1-205/202`=mean(`p205/202`), `Klrk1-203`=mean(p203),
      `Klrk1-204`=mean(p204), `Klrk1-206`=mean(p206), .groups="drop") %>%
    pivot_longer(-condition, names_to="unit", values_to="pct") %>%
    mutate(lab=factor(labs_short[as.character(condition)], levels=labs_short[ord]),
           unit=factor(unit, levels=names(unit_cols)))
  pD <- ggplot(prop, aes(lab, pct, fill=unit)) + geom_col(width=0.7) +
      scale_fill_manual(values=unit_cols) +
      labs(x=NULL, y="Isoform proportion (%)", fill="Unit", tag="D") + th +
      theme(legend.position="right", legend.key.size=unit(0.4,"cm"))
  arrangeGrob(pA,pB,pC,pD, ncol=2, top=title)
}

g4 <- make_panels("GSE203167",
  c("WT_Pre","WT_Post7","TCF7cKO_Pre","TCF7cKO_Post7"),
  c(WT_Pre="WT Pre-Tx", WT_Post7="WT Post-Tx D7", TCF7cKO_Pre="TCF-7 cKO Pre-Tx", TCF7cKO_Post7="TCF-7 cKO Post-Tx D7"),
  "Klrk1 in GSE203167 (WT and TCF-7 cKO CD8+ T cells, pre vs day-7 post-transplant)")
ggsave(file.path(figd,"Figure4_GSE203167_R.pdf"), g4, width=10, height=8)
ggsave(file.path(figd,"Figure4_GSE203167_R.png"), g4, width=10, height=8, dpi=200)

g8 <- make_panels("GSE119943",
  c("Arm_D4.5_EarlyEffector","Arm_D7_SLEC","Arm_D7_MPEC","Cl13_D7_Progenitor",
    "Cl13_D7_TermExhausted","pMIG_Cl13_Progenitor","pMIG_Cl13_TermExhausted"),
  c(Arm_D4.5_EarlyEffector="Arm D4.5 Early Effector", Arm_D7_SLEC="Arm D7 SLEC", Arm_D7_MPEC="Arm D7 MPEC",
    Cl13_D7_Progenitor="Cl13 D7 Progenitor", Cl13_D7_TermExhausted="Cl13 D7 Term. Exhausted",
    pMIG_Cl13_Progenitor="pMIG Cl13 Progenitor", pMIG_Cl13_TermExhausted="pMIG Cl13 Term. Exhausted"),
  "Klrk1 across CD8+ differentiation and exhaustion (GSE119943)")
ggsave(file.path(figd,"Figure8_GSE119943_R.pdf"), g8, width=12, height=8.5)
ggsave(file.path(figd,"Figure8_GSE119943_R.png"), g8, width=12, height=8.5, dpi=200)

# ---- verify vs master ----
m <- read.csv(file.path(tbl,"Table4_units_master.csv"), check.names=FALSE)
chk <- d %>% group_by(dataset,condition) %>% summarise(RI=round(mean(RI),2), .groups="drop") %>%
  inner_join(m %>% transmute(dataset,condition,RI_master=RI_pct), by=c("dataset","condition")) %>%
  mutate(diff=abs(RI-RI_master))
cat("max RI diff vs master:", round(max(chk$diff),3), "(should be ~0)\n")
cat("figures -> ", figd, " (Figure4_GSE203167_R, Figure8_GSE119943_R)\n")
