# ============================================================
# Empirical Work — Econometrics I (EESP, 2026)
# 00_master.R — Master script
#
# Instructions:
#   1. Set the working directory below to your local project root
#   2. Run this script: source("code/00_master.R")
#   All outputs will be saved to output/figures/ and output/tables/
# ============================================================

# ---- Set working directory ----
# setwd("/path/to/empirical-work-econometrics1")  # <-- edit this line

# ---- Install packages if needed ----
required_packages <- c("AER", "dplyr", "ggplot2", "lmtest", "sandwich")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(new_packages)) install.packages(new_packages)

# ---- Create output directories if they don't exist ----
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)

# ---- Run all scripts in order ----
cat("========================================\n")
cat("Running 01_descriptive.R...\n")
cat("========================================\n")
source("code/01_descriptive.R")

cat("========================================\n")
cat("Running 02_firststage.R...\n")
cat("========================================\n")
source("code/02_firststage.R")

cat("========================================\n")
cat("Running 03_regressions.R...\n")
cat("========================================\n")
source("code/03_regressions.R")

cat("========================================\n")
cat("Running 04_tests.R...\n")
cat("========================================\n")
source("code/04_tests.R")

cat("\n========================================\n")
cat("All done. Outputs saved to output/\n")
cat("========================================\n")
