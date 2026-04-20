# ============================================================
# Empirical Work — Econometrics I (EESP, 2026)
# 00_master.R — Master Script
#
# Instructions:
#   1. Set the working directory below to your local project root
#   2. Run: source("code/00_master.R")
#   All outputs will be saved to output/figures/ and output/tables/
#
# Scripts run in order:
#   01_descriptive.R  — Descriptive statistics and balance tests
#   02_firststage.R   — First stage and reduced form
#   03_ols.R          — OLS benchmark regressions
#   04_simple_iv.R    — Simple IV (constant marginal effects)
#   05_ce_tsls.R      — CE-TSLS (Caetano & Escanciano 2021)
#   06_tests.R        — Tests for constant marginal effects
#
# Note on memory: ivreg() with vcovHC() on 440k obs exceeds memory
# limits. Steps 4 and 5 implement 2SLS manually for this reason.
# Point estimates are identical to ivreg().
# ============================================================

# ---- Set working directory ----
# setwd("/path/to/empirical-work-econometrics1")  # <-- edit this

# ---- Install packages if needed ----
required <- c("AER", "dplyr", "ggplot2", "lmtest", "sandwich", "xtable")
missing  <- required[!(required %in% installed.packages()[,"Package"])]
if (length(missing)) {
  cat("Installing missing packages:", paste(missing, collapse=", "), "\n")
  install.packages(missing)
}

# ---- Load all libraries upfront ----
invisible(lapply(required, library, character.only = TRUE))

# ---- Create output directories if they don't exist ----
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Run all scripts in order
# ============================================================

run_step <- function(script, step_name) {
  cat("\n", strrep("=", 55), "\n", sep="")
  cat(" Step:", step_name, "\n")
  cat(strrep("=", 55), "\n", sep="")
  source(file.path("code", script))
  cat("\n[", step_name, "complete ]\n")
}

run_step("01_descriptive.R",
         "01 — Descriptive Statistics and Balance Tests")

run_step("02_firststage.R",
         "02 — First Stage and Reduced Form")

run_step("03_ols.R",
         "03 — OLS Benchmark Regressions")

run_step("04_simple_iv.R",
         "04 — Simple IV (Constant Marginal Effects)")

run_step("05_ce_tsls.R",
         "05 — CE-TSLS (Caetano & Escanciano 2021)")

run_step("06_tests.R",
         "06 — Tests for Constant Marginal Effects")

# ============================================================
# Summary of outputs
# ============================================================
cat("\n", strrep("=", 55), "\n", sep="")
cat(" All steps complete. Outputs:\n")
cat(strrep("=", 55), "\n", sep="")

cat("\nFigures (output/figures/):\n")
figs <- list.files("output/figures", pattern="\\.png$")
for (f in figs) cat("  -", f, "\n")

cat("\nTables (output/tables/):\n")
tbls <- list.files("output/tables", pattern="\\.tex$")
for (t in tbls) cat("  -", t, "\n")