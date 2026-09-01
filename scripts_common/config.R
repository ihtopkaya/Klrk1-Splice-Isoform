# ============================================================================
# config.R  —  single source of truth for all paths in the R analysis.
#
# REPRODUCIBILITY: on a new machine, change ONLY the PROJECT_ROOT default below
# (or set the KLRK1_ROOT environment variable). Everything else is derived.
#
# Each analysis script starts with:
#     PROJECT_ROOT <- Sys.getenv("KLRK1_ROOT",
#                       "/Users/study/Desktop/Karimi/Klrk1_GVHD_project")
#     source(file.path(PROJECT_ROOT, "shared", "scripts_common", "config.R"))
# ============================================================================

if (!exists("PROJECT_ROOT")) {
  PROJECT_ROOT <- Sys.getenv("KLRK1_ROOT",
                    "/Users/study/Desktop/Karimi/Klrk1_GVHD_project")
}
stopifnot(dir.exists(PROJECT_ROOT))

# GEO accession -> dataset folder name under datasets/
DATASET_FOLDER <- c(
  GSE147371 = "GSE147371_GVHD_CD4",
  GSE109125 = "GSE109125_ImmGen",
  GSE203167 = "GSE203167_Karimi_TCF7",
  GSE119943 = "GSE119943_Yao_LCMV"
)

# path to a sample's Salmon quant.sf, given the sample dir name and its dataset
salmon_quant <- function(sample_dir, dataset) {
  folder <- DATASET_FOLDER[[dataset]]
  if (is.null(folder)) stop("unknown dataset: ", dataset)
  file.path(PROJECT_ROOT, "datasets", folder, "data", "salmon",
            sample_dir, "quant.sf")
}

# every sample's quant.sf across all local datasets (for cross-dataset scripts)
all_salmon_quants <- function() {
  Sys.glob(file.path(PROJECT_ROOT, "datasets", "*", "data", "salmon",
                     "*", "quant.sf"))
}

# per-dataset output dirs (created on demand); figures + results live in the
# dataset's own folder so each dataset stays self-contained
results_dir <- function(dataset) {
  d <- file.path(PROJECT_ROOT, "datasets", DATASET_FOLDER[[dataset]], "results")
  dir.create(d, recursive = TRUE, showWarnings = FALSE); d
}
figures_dir <- function(dataset) {
  d <- file.path(PROJECT_ROOT, "datasets", DATASET_FOLDER[[dataset]], "figures")
  dir.create(d, recursive = TRUE, showWarnings = FALSE); d
}
shared_dir <- function(sub) {
  d <- file.path(PROJECT_ROOT, "shared", "cross_dataset_analyses", sub)
  dir.create(d, recursive = TRUE, showWarnings = FALSE); d
}
metadata_csv <- function(name) {
  file.path(PROJECT_ROOT, "shared", "metadata", name)
}

message("config.R loaded — PROJECT_ROOT = ", PROJECT_ROOT)
