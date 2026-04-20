# ============================================================
# Empirical Work — Econometrics I (EESP, 2026)
# 00_master.R — Master Script
#
# Instructions:
#   1. Set the working directory below to your local project root
#   2. Run: source("code/00_master.R")
#   All outputs saved to output/figures/ and output/tables/
#
# Script structure:
#   01_data_and_benchmarks.R  — Descriptive stats, balance check,
#                               first stage, OLS benchmarks
#   02_iv_estimation.R        — Simple IV and CE-TSLS
#   03_tests.R                — Tests for constant marginal effects
#
# Note on memory: ivreg() with vcovHC() on 440k obs exceeds memory.
# IV models are implemented manually as two-stage OLS.
# ============================================================

# ---- Set working directory ----
# setwd("/path/to/empirical-work-econometrics1")  # <-- edit this

# ---- Install packages if needed ----
required <- c("AER", "dplyr", "ggplot2", "lmtest", "sandwich", "xtable")
missing  <- required[!(required %in% installed.packages()[,"Package"])]
if (length(missing)) install.packages(missing)
invisible(lapply(required, library, character.only=TRUE))

# ---- Create output directories ----
dir.create("output/figures", recursive=TRUE, showWarnings=FALSE)
dir.create("output/tables",  recursive=TRUE, showWarnings=FALSE)

# ---- Run scripts ----
run_step <- function(script, name) {
  cat("\n", strrep("=", 55), "\n", sep="")
  cat(" ", name, "\n")
  cat(strrep("=", 55), "\n", sep="")
  source(file.path("code", script))
  cat("[", name, "complete ]\n")
}

run_step("01_data_and_benchmarks.R",
         "01 — Data, Benchmarks & Identification Checks")

run_step("02_iv_estimation.R",
         "02 — Simple IV and CE-TSLS")

run_step("03_tests.R",
         "03 — Tests for Constant Marginal Effects")

# ---- Summary ----
cat("\n", strrep("=", 55), "\n", sep="")
cat(" All done. Outputs:\n")
cat(strrep("=", 55), "\n")

cat("\nFigures:\n")
for (f in list.files("output/figures", pattern="\\.png$"))
  cat("  -", f, "\n")

cat("\nTables (.tex):\n")
for (t in list.files("output/tables", pattern="\\.tex$"))
  cat("  -", t, "\n")
