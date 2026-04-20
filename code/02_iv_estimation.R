# ============================================================
# Empirical Work — Econometrics I (EESP, 2026)
# 02_iv_estimation.R
#
# Contents:
#   1. Simple IV / 2SLS (constant marginal effects)
#   2. CE-TSLS (Caetano & Escanciano, 2021)
#   3. Main results comparison table (Table 5)
#
# Note: ivreg() with vcovHC() exceeds memory on 440k obs.
# Both models are implemented manually as two-stage OLS.
# Point estimates are identical to ivreg().
# ============================================================

library(lmtest)
library(sandwich)
library(xtable)

# ---- Load cleaned data ----
df <- readRDS("output/df_clean.rds")

# ---- Helpers ----
save_tex <- function(tbl, filename, caption, label) {
  xt <- xtable(tbl, caption=caption, label=label)
  print(xt,
        include.rownames       = FALSE,
        caption.placement      = "top",
        booktabs               = TRUE,
        sanitize.text.function = identity,
        table.placement        = "h!",
        size                   = "small",
        file = file.path("output/tables", filename))
  cat("Saved:", filename, "\n")
}

fmt <- function(e, s) paste0(formatC(e, digits=3, format="f"), 
                             " (", formatC(s, digits=3, format="f"), ")")

get_val <- function(ct, v) {
  if (!v %in% rownames(ct)) return("")
  fmt(ct[v,"Estimate"], ct[v,"Std. Error"])
}

# ============================================================
# 1. SIMPLE IV (Table 5, Column 1)
# Instruments S with Z — imposes constant marginal effects
# ============================================================
cat("\n--- SIMPLE IV ---\n")

# Stage 1
fs_s     <- lm(smoked_trimesters ~ support_group + education +
                 age_group + prenatal_high, data=df)
df$S_hat <- fitted(fs_s)

# Stage 2
ss_s  <- lm(baby_weight_kg ~ S_hat + education +
              age_group + prenatal_high, data=df)
ct_iv <- coeftest(ss_s, vcov=vcovHC(ss_s, type="HC1"))

cat("Simple IV estimate (smoked trimesters):",
    round(ct_iv["S_hat","Estimate"], 4),
    "(SE:", round(ct_iv["S_hat","Std. Error"], 4), ")\n")

# ============================================================
# 2. CE-TSLS (Table 5, Column 2)
# Instruments X1, X2, X3 with Z and Z x W interactions
# Recovers separate beta1, beta2, beta3
# ============================================================
cat("\n--- CE-TSLS ---\n")

# Stage 1: first stage for each Xj
# Instruments: Z, Z*age_group, Z*education, Z*prenatal_high
fs_X1 <- lm(X1 ~ support_group * age_group +
               support_group * education +
               support_group * prenatal_high +
               education + age_group + prenatal_high, data=df)

fs_X2 <- lm(X2 ~ support_group * age_group +
               support_group * education +
               support_group * prenatal_high +
               education + age_group + prenatal_high, data=df)

fs_X3 <- lm(X3 ~ support_group * age_group +
               support_group * education +
               support_group * prenatal_high +
               education + age_group + prenatal_high, data=df)

df$X1_hat <- fitted(fs_X1)
df$X2_hat <- fitted(fs_X2)
df$X3_hat <- fitted(fs_X3)

cat("First stages complete.\n")
for (nm in c("fs_X1","fs_X2","fs_X3")) {
  mdl <- get(nm)
  ct  <- coeftest(mdl, vcov=vcovHC(mdl, type="HC1"))
  cat(nm, "— Z coef:", round(ct["support_group","Estimate"],4),
      " SE:", round(ct["support_group","Std. Error"],4), "\n")
}

# Stage 2
ss_ce <- lm(baby_weight_kg ~ X1_hat + X2_hat + X3_hat +
              education + age_group + prenatal_high, data=df)
ct_ce <- coeftest(ss_ce, vcov=vcovHC(ss_ce, type="HC1"))
V     <- vcovHC(ss_ce, type="HC1")

b1  <- ct_ce["X1_hat","Estimate"]; se1 <- ct_ce["X1_hat","Std. Error"]
b2  <- ct_ce["X2_hat","Estimate"]; se2 <- ct_ce["X2_hat","Std. Error"]
b3  <- ct_ce["X3_hat","Estimate"]; se3 <- ct_ce["X3_hat","Std. Error"]

cat("\nCE-TSLS marginal effects:\n")
cat("  beta1 (0->1 trimester): ", round(b1,4), "(", round(se1,4), ")\n")
cat("  beta2 (1->2 trimesters):", round(b2,4), "(", round(se2,4), ")\n")
cat("  beta3 (2->3 trimesters):", round(b3,4), "(", round(se3,4), ")\n")

# Save CE-TSLS results for Code 3 (tests)
saveRDS(list(ct_ce=ct_ce, V=V, b1=b1, b2=b2, b3=b3, ss_ce=ss_ce),
        "output/ce_tsls_results.rds")
cat("CE-TSLS results saved for Code 3.\n")

# ============================================================
# 3. OLS WITH DUMMIES (benchmark for comparison table)
# ============================================================
ols4   <- lm(baby_weight_kg ~ X1 + X2 + X3 + education +
               age_group + prenatal_high, data=df)
ct_ols <- coeftest(ols4, vcov=vcovHC(ols4, type="HC1"))

# ============================================================
# 4. MAIN RESULTS TABLE (Table 5)
# OLS vs Simple IV vs CE-TSLS
# ============================================================
cat("\n--- TABLE 5: Main Results ---\n")

vars_show   <- c("X1","X2","X3",
                 "smoked_trimesters",
                 "educationHigh School","educationCollege","educationMasters",
                 "age_group21-25","age_group26-30","age_group31-40",
                 "prenatal_high")
vars_iv     <- c("S_hat", vars_show[5:11])
vars_ce     <- c("X1_hat","X2_hat","X3_hat", vars_show[5:11])

labels_show <- c(
  "$\\beta_1$: smoked $\\geq 1$ trimester",
  "$\\beta_2$: smoked $\\geq 2$ trimesters",
  "$\\beta_3$: smoked $\\geq 3$ trimesters",
  "Smoked trimesters (level)",
  "Education: High School","Education: College","Education: Masters",
  "Age: 21--25","Age: 26--30","Age: 31--40",
  "Prenatal visits (high)"
)

tbl5 <- data.frame(
  Variable  = labels_show,
  OLS       = c(get_val(ct_ols,"X1"),
                get_val(ct_ols,"X2"),
                get_val(ct_ols,"X3"),
                "",
                get_val(ct_ols,"educationHigh School"),
                get_val(ct_ols,"educationCollege"),
                get_val(ct_ols,"educationMasters"),
                get_val(ct_ols,"age_group21-25"),
                get_val(ct_ols,"age_group26-30"),
                get_val(ct_ols,"age_group31-40"),
                get_val(ct_ols,"prenatal_high")),
  Simple_IV = c("","","",
                get_val(ct_iv,"S_hat"),
                get_val(ct_iv,"educationHigh School"),
                get_val(ct_iv,"educationCollege"),
                get_val(ct_iv,"educationMasters"),
                get_val(ct_iv,"age_group21-25"),
                get_val(ct_iv,"age_group26-30"),
                get_val(ct_iv,"age_group31-40"),
                get_val(ct_iv,"prenatal_high")),
  CE_TSLS   = c(get_val(ct_ce,"X1_hat"),
                get_val(ct_ce,"X2_hat"),
                get_val(ct_ce,"X3_hat"),
                "",
                get_val(ct_ce,"educationHigh School"),
                get_val(ct_ce,"educationCollege"),
                get_val(ct_ce,"educationMasters"),
                get_val(ct_ce,"age_group21-25"),
                get_val(ct_ce,"age_group26-30"),
                get_val(ct_ce,"age_group31-40"),
                get_val(ct_ce,"prenatal_high")),
  stringsAsFactors = FALSE
)

footer5 <- data.frame(
  Variable  = c("Method","Controls","N"),
  OLS       = c("OLS",     "Yes", format(nrow(df),big.mark=",")),
  Simple_IV = c("2SLS",    "Yes", format(nrow(df),big.mark=",")),
  CE_TSLS   = c("CE-TSLS", "Yes", format(nrow(df),big.mark=",")),
  stringsAsFactors = FALSE
)
colnames(footer5) <- colnames(tbl5)
tbl5 <- rbind(tbl5, footer5)
rownames(tbl5) <- NULL
colnames(tbl5) <- c("Variable","OLS","Simple IV","CE-TSLS")

print(tbl5)
write.csv(tbl5, "output/tables/table5_main_results.csv", row.names=FALSE)
save_tex(tbl5,
         filename = "table5_main_results.tex",
         caption  = paste(
           "Main Results: OLS, Simple IV, and CE-TSLS Estimates.",
           "Dependent variable: birth weight in kilograms.",
           "OLS uses three smoking dummies (X1, X2, X3).",
           "Simple IV instruments the level of smoking with Z,",
           "imposing constant marginal effects (Section~1.5).",
           "CE-TSLS implements Caetano and Escanciano (2021):",
           "instruments are Z and $Z \\times W$ interactions,",
           "recovering separate marginal effects $\\beta_1$, $\\beta_2$, $\\beta_3$.",
           "Robust standard errors (HC1) in parentheses.",
           "Baseline: Middle School, age 15--20, 8 prenatal visits."),
         label = "tab:main_results")

cat("\n--- 02_iv_estimation.R complete ---\n")
