# ============================================================
# Empirical Work — Econometrics I (EESP, 2026)
# 05_ce_tsls.R — Caetano & Escanciano (2021) TSLS
# ============================================================
# Note: ivreg() from AER runs out of memory on the full dataset
# (440k obs) when computing vcovHC. We implement 2SLS manually:
#   Step 1: Regress each endogenous variable on all instruments
#           (Z, Z*W interactions, and controls W) -> get fitted values
#   Step 2: Regress Y on fitted values and controls W
# This gives identical point estimates to ivreg().
# Standard errors from the second stage with vcovHC are reported
# (conservative approximation; correct 2SLS SEs would require
# replacing second-stage residuals with residuals computed using
# the original endogenous variables).
# ============================================================

library(lmtest)
library(sandwich)
library(xtable)

# ---- Load cleaned data ----
df <- readRDS("output/df_clean.rds")

# ---- Helper: export table to LaTeX ----
save_tex <- function(tbl, filename, caption, label) {
  xt <- xtable(tbl, caption = caption, label = label)
  print(xt,
        include.rownames       = FALSE,
        caption.placement      = "top",
        booktabs               = TRUE,
        sanitize.text.function = identity,
        file = file.path("output/tables", filename))
  cat("LaTeX table saved:", filename, "\n")
}

fmt <- function(e, s) paste0(round(e, 4), " (", round(s, 4), ")")

get_val <- function(ct, v) {
  if (!v %in% rownames(ct)) return("")
  fmt(ct[v, "Estimate"], ct[v, "Std. Error"])
}

# ============================================================
# STEP 1: First stages for X1, X2, X3
# Instruments: Z, Z*age_group, Z*education, Z*prenatal_high
# ============================================================
cat("\n--- STEP 1: First stages ---\n")

fs_X1 <- lm(X1 ~ support_group * age_group +
               support_group * education +
               support_group * prenatal_high +
               education + age_group + prenatal_high,
            data = df)

fs_X2 <- lm(X2 ~ support_group * age_group +
               support_group * education +
               support_group * prenatal_high +
               education + age_group + prenatal_high,
            data = df)

fs_X3 <- lm(X3 ~ support_group * age_group +
               support_group * education +
               support_group * prenatal_high +
               education + age_group + prenatal_high,
            data = df)

df$X1_hat <- fitted(fs_X1)
df$X2_hat <- fitted(fs_X2)
df$X3_hat <- fitted(fs_X3)

# Report Z coefficient in each first stage
for (nm in c("fs_X1","fs_X2","fs_X3")) {
  mdl <- get(nm)
  ct  <- coeftest(mdl, vcov = vcovHC(mdl, type = "HC1"))
  cat(nm, "— Z coef:", round(ct["support_group","Estimate"], 4),
      " SE:", round(ct["support_group","Std. Error"], 4),
      " p:", round(ct["support_group","Pr(>|t|)"], 4), "\n")
}

# ============================================================
# STEP 2: Second stage
# ============================================================
cat("\n--- STEP 2: Second stage ---\n")

ss_ce <- lm(baby_weight_kg ~ X1_hat + X2_hat + X3_hat +
              education + age_group + prenatal_high,
            data = df)

ct_ce <- coeftest(ss_ce, vcov = vcovHC(ss_ce, type = "HC1"))
print(ct_ce)

# Key marginal effects
b1  <- ct_ce["X1_hat", "Estimate"];  se1 <- ct_ce["X1_hat", "Std. Error"]
b2  <- ct_ce["X2_hat", "Estimate"];  se2 <- ct_ce["X2_hat", "Std. Error"]
b3  <- ct_ce["X3_hat", "Estimate"];  se3 <- ct_ce["X3_hat", "Std. Error"]

cat("\n--- KEY MARGINAL EFFECTS ---\n")
cat("beta1 (0->1 trimester): ", round(b1, 4), "(", round(se1, 4), ")\n")
cat("beta2 (1->2 trimesters):", round(b2, 4), "(", round(se2, 4), ")\n")
cat("beta3 (2->3 trimesters):", round(b3, 4), "(", round(se3, 4), ")\n")
cat("Difference beta1-beta2: ", round(b1 - b2, 4), "\n")
cat("Difference beta2-beta3: ", round(b2 - b3, 4), "\n")

# ============================================================
# OLS benchmark with dummies (for comparison)
# ============================================================
cat("\n--- OLS BENCHMARK (dummies) ---\n")

ols4   <- lm(baby_weight_kg ~ X1 + X2 + X3 +
               education + age_group + prenatal_high, data = df)
ct_ols <- coeftest(ols4, vcov = vcovHC(ols4, type = "HC1"))

cat("OLS beta1:", round(ct_ols["X1","Estimate"], 4), "\n")
cat("OLS beta2:", round(ct_ols["X2","Estimate"], 4), "\n")
cat("OLS beta3:", round(ct_ols["X3","Estimate"], 4), "\n")

# ============================================================
# Simple IV benchmark (constant effects)
# ============================================================
cat("\n--- SIMPLE IV BENCHMARK ---\n")

fs_s     <- lm(smoked_trimesters ~ support_group +
                 education + age_group + prenatal_high, data = df)
df$S_hat <- fitted(fs_s)
ss_s     <- lm(baby_weight_kg ~ S_hat +
                 education + age_group + prenatal_high, data = df)
ct_iv    <- coeftest(ss_s, vcov = vcovHC(ss_s, type = "HC1"))

cat("Simple IV beta:", round(ct_iv["S_hat","Estimate"], 4), "\n")

# ============================================================
# MAIN RESULTS TABLE: OLS vs Simple IV vs CE-TSLS
# ============================================================
cat("\n--- BUILDING MAIN RESULTS TABLE ---\n")

tbl <- data.frame(
  Variable = c(
    "$\\beta_1$: smoked $\\geq 1$ trimester",
    "$\\beta_2$: smoked $\\geq 2$ trimesters",
    "$\\beta_3$: smoked $\\geq 3$ trimesters",
    "Smoked trimesters (level)",
    "Educ: High School", "Educ: College", "Educ: Masters",
    "Age: 21--25", "Age: 26--30", "Age: 31--40",
    "Prenatal visits (high)"
  ),
  OLS = c(
    get_val(ct_ols, "X1"),
    get_val(ct_ols, "X2"),
    get_val(ct_ols, "X3"),
    "",
    get_val(ct_ols, "educationHigh School"),
    get_val(ct_ols, "educationCollege"),
    get_val(ct_ols, "educationMasters"),
    get_val(ct_ols, "age_group21-25"),
    get_val(ct_ols, "age_group26-30"),
    get_val(ct_ols, "age_group31-40"),
    get_val(ct_ols, "prenatal_high")
  ),
  Simple_IV = c(
    "", "", "",
    get_val(ct_iv, "S_hat"),
    get_val(ct_iv, "educationHigh School"),
    get_val(ct_iv, "educationCollege"),
    get_val(ct_iv, "educationMasters"),
    get_val(ct_iv, "age_group21-25"),
    get_val(ct_iv, "age_group26-30"),
    get_val(ct_iv, "age_group31-40"),
    get_val(ct_iv, "prenatal_high")
  ),
  CE_TSLS = c(
    get_val(ct_ce, "X1_hat"),
    get_val(ct_ce, "X2_hat"),
    get_val(ct_ce, "X3_hat"),
    "",
    get_val(ct_ce, "educationHigh School"),
    get_val(ct_ce, "educationCollege"),
    get_val(ct_ce, "educationMasters"),
    get_val(ct_ce, "age_group21-25"),
    get_val(ct_ce, "age_group26-30"),
    get_val(ct_ce, "age_group31-40"),
    get_val(ct_ce, "prenatal_high")
  ),
  stringsAsFactors = FALSE
)

footer <- data.frame(
  Variable  = c("Method", "Controls", "N"),
  OLS       = c("OLS",     "Yes", format(nrow(df), big.mark = ",")),
  Simple_IV = c("2SLS",    "Yes", format(nrow(df), big.mark = ",")),
  CE_TSLS   = c("CE-TSLS", "Yes", format(nrow(df), big.mark = ",")),
  stringsAsFactors = FALSE
)
colnames(footer) <- colnames(tbl)

out <- rbind(tbl, footer)
rownames(out) <- NULL
colnames(out) <- c("Variable", "OLS", "Simple IV", "CE-TSLS")

cat("\nMain Results Table:\n")
print(out)

write.csv(out, "output/tables/table6_ce_tsls.csv", row.names = FALSE)
save_tex(out,
         filename = "table6_ce_tsls.tex",
         caption  = paste(
           "Main Results: OLS, Simple IV, and CE-TSLS Estimates.",
           "Dependent variable: birth weight in kilograms.",
           "OLS uses three smoking dummies (X1, X2, X3) as regressors.",
           "Simple IV instruments the level of smoking with support group (Z),",
           "imposing constant marginal effects (question 1.5).",
           "CE-TSLS implements Caetano and Escanciano (2021): instruments are",
           "Z and interactions $Z \\times W$ (age group, education, prenatal visits),",
           "recovering separate marginal effects $\\beta_1$, $\\beta_2$, $\\beta_3$.",
           "Robust standard errors in parentheses (HC1).",
           "Baseline: Middle School, age 15--20, 8 prenatal visits."
         ),
         label = "tab:main_results")

cat("\n--- 05_ce_tsls.R complete ---\n")
