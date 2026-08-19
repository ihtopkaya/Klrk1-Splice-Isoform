# install_dependencies.R — one-time setup for the R analysis
# Run once:  Rscript shared/scripts_common/install_dependencies.R
#
# Tested with R 4.6.0. The analysis uses a single CRAN meta-package.
pkgs <- c("tidyverse")   # ggplot2, dplyr, readr, tidyr, purrr, ...
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
} else {
  message("All required R packages already installed.")
}
