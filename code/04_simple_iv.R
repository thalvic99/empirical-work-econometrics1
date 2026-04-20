# ============================================================
# Empirical Work — Econometrics I (EESP, 2026)
# 04_simple_iv.R — Simple IV (Constant Marginal Effects)
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

fmt <- function(est, se) paste0(est, " (", se, ")")

extract_iv <- function(model, varname) {
  ct <- coeftest(model, vcov = vcovHC(model, type = "HC1"))
  if (!varname %in% rownames(ct)) return(c(NA, NA))
  c(round(ct[varname, "Estimate"], 4),
    round(ct[varname, "Std. Error"], 4))
}

# ============================================================
# WALD ESTIMATOR: Reduced form / First stage (no controls)
# The simplest possible IV estimate
# ============================================================
cat("\n--- WALD ESTIMATOR (no controls) ---\n")

rf_simple <- lm(baby_weight_kg ~ support_group, data = df)
fs_simple <- lm(smoked_trimesters ~ support_group, data = df)

rf_coef <- coef(rf_simple)["support_group"]
fs_coef <- coef(fs_simple)["support_group"]
wald    <- rf_coef / fs_coef

cat("Reduced form (Z -> Y):", round(rf_coef, 4), "\n")
cat("First stage  (Z -> S):", round(fs_coef, 4), "\n")
cat("Wald estimate (RF/FS):", round(wald, 4), "\n")

# ============================================================
# IV 1: ivreg — smoking level as treatment, no controls
# Corresponds to model from question 1.5
# ============================================================
cat("\n--- IV 1: Smoking level, no controls ---\n")

iv1 <- ivreg(baby_weight_kg ~ smoked_trimesters |
               support_group,
             data = df)
ct_iv1 <- coeftest(iv1, vcov = vcovHC(iv1, type = "HC1"))
print(ct_iv1)

# ============================================================
# IV 2: ivreg — smoking level as treatment, with controls
# Main simple IV specification
# ============================================================
cat("\n--- IV 2: Smoking level, with controls ---\n")

iv2 <- ivreg(baby_weight_kg ~ smoked_trimesters +
               education + age_group + prenatal_high |
               support_group +
               education + age_group + prenatal_high,
             data = df)
ct_iv2 <- coeftest(iv2, vcov = vcovHC(iv2, type = "HC1"))
print(ct_iv2)

cat("\nSimple IV estimate (with controls):",
    round(ct_iv2["smoked_trimesters","Estimate"], 4), "\n")
cat("Compare with OLS (with controls): -0.3635\n")
cat("Compare with Wald (no controls): ", round(wald, 4), "\n")

# ============================================================
# BUILD SUMMARY TABLE: OLS vs Simple IV comparison
# ============================================================
cat("\n--- SUMMARY TABLE: OLS vs Simple IV ---\n")

# Load OLS models from Step 3 by re-estimating
ols3 <- lm(baby_weight_kg ~ smoked_trimesters +
             education + age_group + prenatal_high,
           data = df)

vars_show <- c("smoked_trimesters",
               "educationHigh School", "educationCollege", "educationMasters",
               "age_group21-25", "age_group26-30", "age_group31-40",
               "prenatal_high")

labels <- c("Smoked trimesters (level)",
            "Education: High School", "Education: College", "Education: Masters",
            "Age: 21--25", "Age: 26--30", "Age: 31--40",
            "Prenatal visits (high)")

models    <- list(ols3, iv1, iv2)
mod_names <- c("OLS", "Simple IV", "Simple IV + Controls")

out_rows <- lapply(seq_along(vars_show), function(i) {
  v   <- vars_show[i]
  row <- data.frame(Variable = labels[i])
  for (j in seq_along(models)) {
    ct  <- coeftest(models[[j]], vcov = vcovHC(models[[j]], type = "HC1"))
    val <- if (!v %in% rownames(ct)) "" else
      fmt(round(ct[v,"Estimate"], 4), round(ct[v,"Std. Error"], 4))
    row[[mod_names[j]]] <- val
  }
  row
})

out_tbl <- do.call(rbind, out_rows)

# Add method, N, R-squared rows
method_row <- data.frame(
  Variable             = "Method",
  OLS                  = "OLS",
  "Simple IV"          = "IV",
  "Simple IV + Controls" = "IV",
  check.names          = FALSE
)
ctrl_row <- data.frame(
  Variable             = "Controls",
  OLS                  = "Yes",
  "Simple IV"          = "No",
  "Simple IV + Controls" = "Yes",
  check.names          = FALSE
)
n_row <- data.frame(
  Variable             = "N",
  OLS                  = format(nrow(df), big.mark=","),
  "Simple IV"          = format(nrow(df), big.mark=","),
  "Simple IV + Controls" = format(nrow(df), big.mark=","),
  check.names          = FALSE
)

colnames(out_tbl) <- c("Variable", mod_names)
colnames(method_row) <- colnames(ctrl_row) <- colnames(n_row) <- c("Variable", mod_names)

out_tbl <- rbind(out_tbl, method_row, ctrl_row, n_row)
rownames(out_tbl) <- NULL

cat("\nOLS vs Simple IV:\n")
print(out_tbl)

write.csv(out_tbl, "output/tables/table5_simple_iv.csv", row.names = FALSE)
save_tex(out_tbl,
         filename = "table5_simple_iv.tex",
         caption  = paste(
           "OLS and Simple IV Estimates.",
           "Dependent variable: birth weight in kilograms.",
           "Simple IV uses support group assignment (Z) as instrument",
           "for the level of smoking (trimesters smoked),",
           "imposing constant marginal effects across smoking levels.",
           "Heteroskedasticity-robust standard errors in parentheses.",
           "Baseline: Middle School education, age 15--20,",
           "low prenatal visits (8 visits)."
         ),
         label = "tab:simple_iv")

cat("\n--- 04_simple_iv.R complete ---\n")
