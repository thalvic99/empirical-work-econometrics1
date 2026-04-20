# ============================================================
# Empirical Work — Econometrics I (EESP, 2026)
# 00_master.R — Master Script
#
# Instructions:
#   1. Open empirical-work-econometrics1.Rproj in RStudio
#   2. If first time: renv::restore()
#   3. Run: source("code/00_master.R")
#
# Script structure:
#   01_data_and_benchmarks.R  — Descriptive stats, balance check,
#                               first stage + reduced form, OLS
#   02_iv_estimation.R        — Simple IV and CE-TSLS (main results)
#   03_tests.R                — All tests: constant effects, separability,
#                               and robustness check (prenatal visits)
#
# Note: IV models implemented manually as two-stage OLS.
# ivreg() with vcovHC() exceeds memory on 440k observations.
# Packages managed by renv — run renv::restore() before first use.
# ============================================================

# ---- Load libraries ----
required <- c("AER","dplyr","ggplot2","lmtest","sandwich","xtable")
invisible(lapply(required, library, character.only=TRUE))

# ---- Create output directories ----
dir.create("output/figures", recursive=TRUE, showWarnings=FALSE)
dir.create("output/tables",  recursive=TRUE, showWarnings=FALSE)

# ---- Helper ----
run_step <- function(script, name) {
  cat("\n", strrep("=", 55), "\n", sep="")
  cat(" ", name, "\n")
  cat(strrep("=", 55), "\n", sep="")
  source(file.path("code", script))
  cat("[", name, "complete ]\n")
}

# ---- Run all scripts ----
run_step("01_data_and_benchmarks.R",
         "01 — Data, Benchmarks & Identification Checks")

run_step("02_iv_estimation.R",
         "02 — Simple IV and CE-TSLS")

run_step("03_tests.R",
         "03 — Tests and Robustness Checks")

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
