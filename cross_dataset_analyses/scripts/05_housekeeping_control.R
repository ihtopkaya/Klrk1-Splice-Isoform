# ============================================================================
# Task 1: Housekeeping Gene Negative Control
# Tests whether the activation-dependent retained-intron pattern at Klrk1 is
# gene-specific or a pipeline artifact (e.g. gDNA contamination).
#
# Hypothesis: Actb and Gapdh - which are constitutively transcribed but NOT
# subject to activation-dependent splicing regulation - should show
# negligible retained-intron isoform expression across all conditions.
# ============================================================================

library(tidyverse)

# cross-dataset control: paths from shared config
PROJECT_ROOT <- Sys.getenv("KLRK1_ROOT",
                  "/Users/study/Desktop/Karimi/Klrk1_GVHD_project")
source(file.path(PROJECT_ROOT, "shared", "scripts_common", "config.R"))
table_dir <- shared_dir("results")   # shared/cross_dataset_analyses/results
fig_dir   <- shared_dir("figures")   # shared/cross_dataset_analyses/figures

# ---- Transcript annotations (from GRCm39 Ensembl 111 cDNA FASTA headers) ----

actb_meta <- tibble(
  transcript_id = c(
    "ENSMUST00000100497", "ENSMUST00000167721", "ENSMUST00000171419",
    "ENSMUST00000163829", "ENSMUST00000106216",
    "ENSMUST00000196997", "ENSMUST00000167386",
    "ENSMUST00000165629", "ENSMUST00000164765"
  ),
  gene = "Actb",
  biotype = c(rep("protein_coding", 5),
              rep("protein_coding_CDS_not_defined", 2),
              rep("retained_intron", 2))
)

gapdh_meta <- tibble(
  transcript_id = c(
    "ENSMUST00000117757", "ENSMUST00000118875", "ENSMUST00000073605",
    "ENSMUST00000183272", "ENSMUST00000182052", "ENSMUST00000182277",
    "ENSMUST00000144205", "ENSMUST00000182115", "ENSMUST00000182670",
    "ENSMUST00000192506", "ENSMUST00000147954", "ENSMUST00000182464",
    "ENSMUST00000144588"
  ),
  gene = "Gapdh",
  biotype = c(rep("protein_coding", 6),
              rep("protein_coding_CDS_not_defined", 3),
              rep("retained_intron", 4))
)

klrk1_meta <- tibble(
  transcript_id = c(
    "ENSMUST00000032252", "ENSMUST00000095412", "ENSMUST00000168919",
    "ENSMUST00000152256", "ENSMUST00000137660", "ENSMUST00000204694"
  ),                                    # all SIX Klrk1 transcripts (Klrk1-202 now included)
  gene = "Klrk1",
  biotype = c("protein_coding", "protein_coding", "protein_coding",
              "retained_intron", "retained_intron", "retained_intron")
)

all_meta <- bind_rows(actb_meta, gapdh_meta, klrk1_meta)
target_tx <- all_meta$transcript_id

# ---- All Salmon outputs across every dataset (datasets/*/data/salmon/*/quant.sf) ----
qfiles <- all_salmon_quants()
cat("Salmon quant.sf files found:", length(qfiles), "\n")

read_targets <- function(qpath) {
  if (!file.exists(qpath)) return(NULL)
  # sample name = the sample folder name (parent dir of quant.sf)
  sample_name <- basename(dirname(qpath))
  read_tsv(qpath, col_types = "cdddd", show_col_types = FALSE) %>%
    mutate(transcript_id = sub("\\.\\d+$", "", Name)) %>%
    filter(transcript_id %in% target_tx) %>%
    transmute(sample = sample_name, transcript_id, TPM)
}

all_tpm <- map_dfr(qfiles, read_targets) %>%
  left_join(all_meta, by = "transcript_id")

cat("Total measurements:", nrow(all_tpm), "\n")

# ---- Per-gene summary by sample ----
gene_totals <- all_tpm %>%
  group_by(sample, gene) %>%
  summarise(total_gene_tpm = sum(TPM), .groups = "drop")

ri_per_sample <- all_tpm %>%
  filter(biotype == "retained_intron") %>%
  group_by(sample, gene) %>%
  summarise(ri_tpm = sum(TPM), .groups = "drop") %>%
  left_join(gene_totals, by = c("sample", "gene")) %>%
  mutate(ri_pct = ifelse(total_gene_tpm > 0, ri_tpm / total_gene_tpm * 100, 0))

# ---- Master summary table ----
summary_tbl <- ri_per_sample %>%
  group_by(gene) %>%
  summarise(
    n_samples              = n(),
    mean_total_gene_tpm    = mean(total_gene_tpm),
    sd_total_gene_tpm      = sd(total_gene_tpm),
    mean_RI_tpm            = mean(ri_tpm),
    max_RI_tpm             = max(ri_tpm),
    n_samples_RI_above_1   = sum(ri_tpm > 1),
    n_samples_RI_above_5   = sum(ri_tpm > 5),
    mean_RI_pct            = mean(ri_pct),
    sd_RI_pct              = sd(ri_pct),
    .groups = "drop"
  ) %>%
  arrange(factor(gene, levels = c("Actb", "Gapdh", "Klrk1")))

cat("\n=== HOUSEKEEPING NEGATIVE CONTROL - SUMMARY ===\n")
print(summary_tbl, width = Inf)

write_csv(summary_tbl,    file.path(table_dir, "TableV2_housekeeping_summary.csv"))
write_csv(ri_per_sample,  file.path(table_dir, "TableV2_housekeeping_per_sample.csv"))
write_csv(all_tpm,        file.path(table_dir, "TableV2_housekeeping_all_transcripts.csv"))

# ---- Per-RI-transcript max TPM (the strict test) ----
ri_tx_max <- all_tpm %>%
  filter(biotype == "retained_intron") %>%
  group_by(gene, transcript_id) %>%
  summarise(
    max_TPM = max(TPM),
    mean_TPM = mean(TPM),
    n_samples_above_1_TPM = sum(TPM > 1),
    .groups = "drop"
  ) %>%
  arrange(gene, desc(max_TPM))

cat("\n=== PER-RI-TRANSCRIPT MAX EXPRESSION (across all samples) ===\n")
print(ri_tx_max, n = 30, width = Inf)

write_csv(ri_tx_max, file.path(table_dir, "TableV2_RI_transcript_max_per_gene.csv"))

# ---- Visualization ----
plot_data <- ri_per_sample %>%
  mutate(gene = factor(gene, levels = c("Actb", "Gapdh", "Klrk1")))

p <- ggplot(plot_data, aes(x = gene, y = ri_pct, fill = gene)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.4, alpha = 0.65) +
  scale_fill_manual(values = c("Actb" = "#7FBF7B", "Gapdh" = "#7FBF7B", "Klrk1" = "#D6604D")) +
  labs(
    title    = "Retained-Intron Isoform Proportion: Klrk1 vs Housekeeping Genes",
    subtitle = paste0("Each dot = one biological sample; bulk RNA-seq Salmon TPM (n=",
                      dplyr::n_distinct(all_tpm$sample), " samples total)"),
    x        = NULL,
    y        = "Retained-intron isoforms (% of total gene TPM)"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 12))

ggsave(file.path(fig_dir, "FigV5_housekeeping_negative_control.pdf"), p, width = 7, height = 5.5)
ggsave(file.path(fig_dir, "FigV5_housekeeping_negative_control.png"), p, width = 7, height = 5.5, dpi = 300)

cat("\n=== Outputs ===\n")
cat("Tables:", table_dir, "\n")
cat("Figure:", file.path(fig_dir, "FigV5_housekeeping_negative_control.png"), "\n")
