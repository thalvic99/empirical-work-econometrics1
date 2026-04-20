# ============================================================
# Empirical Work — Econometrics I (EESP, 2026)
# 03_ols.R — OLS Benchmark Regressions
# ============================================================

library(AER)
library(dplyr)
library(lmtest)
library(sandwich)
library(xtable)

# ---- Load cleaned data ----
df <- readRDS("output/df_clean.rds")

# ---- Helper: export table to LaTeX ----
save_tex <- function(tbl, filename, caption, label, digits = 3) {
  xt <- xtable(tbl,
               caption = caption,
               label   = label,
               digits  = digits)
  print(xt,
        include.rownames  = FALSE,
        caption.placement = "top",
        booktabs          = TRUE,
        file              = file.path("output/tables", filename))
  cat("LaTeX table saved:", filename, "\n")
}

# ---- Helper: extract coefficient + SE from model ----
extract_row <- function(model, varname, se_type = "HC1") {
  ct <- coeftest(model, vcov = vcovHC(model, type = se_type))
  if (!varname %in% rownames(ct)) return(c(NA, NA))
  c(round(ct[varname, "Estimate"], 4),
    round(ct[varname, "Std. Error"], 4))
}

fmt <- function(est, se) paste0(est, " (", se, ")")

# ============================================================
# OLS 1: Baby weight on binary smoking dummy (no controls)
# Mirrors question 1.1 exactly
# ============================================================
cat("\n--- OLS 1: Binary smoking dummy, no controls ---\n")

ols1 <- lm(baby_weight_kg ~ I(smoked_trimesters > 0), data = df)
ct1  <- coeftest(ols1, vcov = vcovHC(ols1, type = "HC1"))
print(ct1)

# ============================================================
# OLS 2: Smoking level (trimesters) as continuous, no controls
# Mirrors question 1.5 — constant marginal effects assumption
# ============================================================
cat("\n--- OLS 2: Smoking level (continuous), no controls ---\n")

ols2 <- lm(baby_weight_kg ~ smoked_trimesters, data = df)
ct2  <- coeftest(ols2, vcov = vcovHC(ols2, type = "HC1"))
print(ct2)

# ============================================================
# OLS 3: Smoking level (trimesters) as continuous, with controls
# ============================================================
cat("\n--- OLS 3: Smoking level (continuous), with controls ---\n")

ols3 <- lm(baby_weight_kg ~ smoked_trimesters + education + age_group + prenatal_high,
           data = df)
ct3  <- coeftest(ols3, vcov = vcovHC(ols3, type = "HC1"))
print(ct3)

# ============================================================
# OLS 4: Three smoking dummies (X1, X2, X3), with controls
# Mirrors question 1.4 — unrestricted model
# ============================================================
cat("\n--- OLS 4: Three smoking dummies, with controls ---\n")

ols4 <- lm(baby_weight_kg ~ X1 + X2 + X3 + education + age_group + prenatal_high,
           data = df)
ct4  <- coeftest(ols4, vcov = vcovHC(ols4, type = "HC1"))
print(ct4)

cat("\nMarginal effects from OLS 4:\n")
cat("  beta1 (0->1 trimester):", round(ct4["X1","Estimate"], 4), "\n")
cat("  beta2 (1->2 trimesters):", round(ct4["X2","Estimate"], 4), "\n")
cat("  beta3 (2->3 trimesters):", round(ct4["X3","Estimate"], 4), "\n")

# ============================================================
# BUILD SUMMARY TABLE: OLS specifications side by side
# ============================================================
cat("\n--- SUMMARY TABLE: OLS Specifications ---\n")

vars_to_show <- c(
  "I(smoked_trimesters > 0)TRUE",
  "smoked_trimesters",
  "X1", "X2", "X3",
  "educationHigh School", "educationCollege", "educationMasters",
  "age_group21-25", "age_group26-30", "age_group31-40",
  "prenatal_high"
)

labels <- c(
  "Smoked (binary)",
  "Smoked trimesters (level)",
  "X1: smoked >= 1 trimester",
  "X2: smoked >= 2 trimesters",
  "X3: smoked >= 3 trimesters",
  "Education: High School", "Education: College", "Education: Masters",
  "Age: 21-25", "Age: 26-30", "Age: 31-40",
  "Prenatal visits (high)"
)

models <- list(ols1, ols2, ols3, ols4)

out_rows <- lapply(seq_along(vars_to_show), function(i) {
  v <- vars_to_show[i]
  row <- data.frame(Variable = labels[i])
  for (j in seq_along(models)) {
    r   <- extract_row(models[[j]], v)
    val <- if (is.na(r[1])) "" else fmt(r[1], r[2])
    row[[paste0("OLS_", j)]] <- val
  }
  row
})

out_tbl <- do.call(rbind, out_rows)

# Add controls indicator
ctrl_row <- data.frame(
  Variable = "Controls",
  OLS_1    = "No",
  OLS_2    = "No",
  OLS_3    = "Yes",
  OLS_4    = "Yes"
)

# Add N and R-squared
n_row <- data.frame(
  Variable = "N",
  OLS_1    = format(nrow(df), big.mark=","),
  OLS_2    = format(nrow(df), big.mark=","),
  OLS_3    = format(nrow(df), big.mark=","),
  OLS_4    = format(nrow(df), big.mark=",")
)

r2_row <- data.frame(
  Variable = "R-squared",
  OLS_1    = round(summary(ols1)$r.squared, 3),
  OLS_2    = round(summary(ols2)$r.squared, 3),
  OLS_3    = round(summary(ols3)$r.squared, 3),
  OLS_4    = round(summary(ols4)$r.squared, 3)
)

colnames(out_tbl) <- colnames(ctrl_row) <- colnames(n_row) <- colnames(r2_row) <-
  c("Variable", "(1)", "(2)", "(3)", "(4)")

out_tbl <- rbind(out_tbl, ctrl_row, n_row, r2_row)
rownames(out_tbl) <- NULL

cat("\nOLS Summary Table:\n")
print(out_tbl)

write.csv(out_tbl, "output/tables/table4_ols.csv", row.names = FALSE)
save_tex(out_tbl,
         filename = "table4_ols.tex",
         caption  = paste(
           "OLS Estimates: Effect of Smoking on Birth Weight.",
           "Dependent variable: birth weight in kilograms.",
           "Column (1): binary smoking dummy, no controls.",
           "Column (2): smoking level (trimesters), no controls.",
           "Column (3): smoking level (trimesters), with controls.",
           "Column (4): three smoking dummies (X1, X2, X3), with controls.",
           "Heteroskedasticity-robust standard errors in parentheses.",
           "All OLS estimates are biased due to endogeneity of smoking;",
           "see Section 1.1 for discussion.",
           "Baseline: Middle School education, age 15--20, low prenatal visits."
         ),
         label = "tab:ols")

cat("\n--- 03_ols.R complete ---\n")
