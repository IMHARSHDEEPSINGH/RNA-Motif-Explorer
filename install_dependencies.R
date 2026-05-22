#!/usr/bin/env Rscript
# =============================================================================
# RNA Motif Discovery Explorer
# File: install_dependencies.R
# Purpose: Install all required R and Bioconductor packages.
# Run once with: Rscript install_dependencies.R
# =============================================================================

cat("=== RNA Motif Discovery Explorer — Package Installer ===\n\n")
cat("R version:", paste(R.version$major, R.version$minor, sep="."), "\n")

# CRAN packages
cran_pkgs <- c(
  "shiny",
  "shinydashboard",
  "shinyjs",
  "DT",
  "plotly",
  "ggplot2",
  "ggseqlogo",
  "dplyr",
  "tidyr",
  "tibble",
  "stringr",
  "writexl",
  "rmarkdown",
  "knitr",
  "gridExtra",
  "scales",
  "RColorBrewer",
  "plumber",
  "testthat",
  "igraph"
)

bioc_pkgs <- c(
  "Biostrings",
  "universalmotif"
  # Note: RNAstructure is optional - if unavailable, heuristic prediction is used
)

cat("Installing CRAN packages...\n")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing: %s\n", pkg))
    install.packages(pkg,
                     repos   = "https://cloud.r-project.org",
                     quiet   = TRUE,
                     dependencies = TRUE)
  } else {
    cat(sprintf("  Already installed: %s\n", pkg))
  }
}

cat("\nInstalling Bioconductor packages...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install(bioc_pkgs, ask = FALSE, quiet = TRUE, update = FALSE)

# Verify all packages load
cat("\nVerifying installations...\n")
all_pkgs   <- c(cran_pkgs, bioc_pkgs)
failed     <- character(0)

for (pkg in all_pkgs) {
  ok <- tryCatch({
    requireNamespace(pkg, quietly = TRUE)
  }, error = function(e) FALSE)
  if (ok) {
    cat(sprintf("  [OK]  %s\n", pkg))
  } else {
    cat(sprintf("  [FAIL] %s\n", pkg))
    failed <- c(failed, pkg)
  }
}

if (length(failed) == 0) {
  cat("\n=== All packages installed successfully! ===\n")
  if (!requireNamespace("rmarkdown", quietly = TRUE) || !rmarkdown::pandoc_available("1.12.3")) {
    cat("WARNING: pandoc 1.12.3+ is required for PDF report generation.\n")
    cat("If you do not have pandoc, install it from https://pandoc.org or use RStudio.\n")
  }
  cat("Launch the app with:\n")
  cat("  shiny::runApp('shiny_app/')\n\n")
} else {
  cat(sprintf("\n=== WARNING: %d package(s) failed to install: %s ===\n",
              length(failed), paste(failed, collapse = ", ")))
  cat("Try installing manually or check system dependencies.\n\n")
}
