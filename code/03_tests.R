# ============================================================
# Empirical Work — Econometrics I (EESP, 2026)
# 03_tests.R — All Specification Tests
#
# Contents:
#   PART 1: Tests for constant marginal effects
#     1. H0: beta1 = beta2  (t-test)
#     2. H0: beta2 = beta3  (t-test)
#     3. H0: beta1 = beta2 = beta3  (joint Wald F-test)
#     4. Separability test: age_group x X1_hat interaction
#   PART 2: Robustness check
#     5. Simple IV and CE-TSLS with and without prenatal visits
#        Tests sensitivity to bad control concern from balance test
#
# Loads CE-TSLS results saved by 02_iv_estimation.R
# ============================================================

library(lmtest)
library(sandwich)
library(xtable)

# ---- Load data and CE-TSLS results ----
df      <- readRDS("output/df_clean.rds")
results <- readRDS("output/ce_tsls_results.rds")

ct_ce <- results$ct_ce
V     <- results$V
b1    <- results$b1
b2    <- results$b2
b3    <- results$b3
ss_ce <- results$ss_ce
dof   <- nrow(df) - ncol(model.matrix(ss_ce))

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

format_p <- function(p) {
  ifelse(p < 0.001, "$<$0.001", as.character(round(p, 3)))
}

cat("\nCE-TSLS estimates loaded:\n")
cat("  beta1:", round(b1,4), "\n")
cat("  beta2:", round(b2,4), "\n")
cat("  beta3:", round(b3,4), "\n")

# ============================================================
# PART 1: TESTS FOR CONSTANT MARGINAL EFFECTS
# ============================================================

# ---- TEST 1: H0: beta1 = beta2 ----
cat("\n--- TEST 1: H0: beta1 = beta2 ---\n")

var_diff12 <- V["X1_hat","X1_hat"] + V["X2_hat","X2_hat"] -
              2*V["X1_hat","X2_hat"]
t12 <- (b1 - b2) / sqrt(var_diff12)
p12 <- 2 * pt(-abs(t12), df=dof)

cat("Difference (b1-b2):", round(b1-b2, 4), "\n")
cat("SE(b1-b2):         ", round(sqrt(var_diff12), 4), "\n")
cat("t-statistic:       ", round(t12, 3), "\n")
cat("p-value:           ", round(p12, 4), "\n")

# ---- TEST 2: H0: beta2 = beta3 ----
cat("\n--- TEST 2: H0: beta2 = beta3 ---\n")

var_diff23 <- V["X2_hat","X2_hat"] + V["X3_hat","X3_hat"] -
              2*V["X2_hat","X3_hat"]
t23 <- (b2 - b3) / sqrt(var_diff23)
p23 <- 2 * pt(-abs(t23), df=dof)

cat("Difference (b2-b3):", round(b2-b3, 4), "\n")
cat("SE(b2-b3):         ", round(sqrt(var_diff23), 4), "\n")
cat("t-statistic:       ", round(t23, 3), "\n")
cat("p-value:           ", round(p23, 4), "\n")

# ---- TEST 3: Joint H0: beta1 = beta2 = beta3 ----
cat("\n--- TEST 3: Joint H0: beta1 = beta2 = beta3 ---\n")

coef_idx <- c("X1_hat","X2_hat","X3_hat")
R        <- matrix(c(1,-1,0, 0,1,-1), nrow=2, byrow=TRUE)
b_vec    <- c(b1, b2, b3)
V3       <- V[coef_idx, coef_idx]
Rb       <- R %*% b_vec
RVR      <- R %*% V3 %*% t(R)
F_stat   <- as.numeric(t(Rb) %*% solve(RVR) %*% Rb / 2)
p_F      <- pf(F_stat, df1=2, df2=dof, lower.tail=FALSE)

cat("F-statistic:", round(F_stat, 3), "\n")
cat("p-value:    ", round(p_F, 4), "\n")
cat("Conclusion: ", ifelse(p_F < 0.05,
    "Reject H0 — marginal effects are NOT constant.",
    "Fail to reject H0."), "\n")

# ---- TEST 4: Separability test ----
cat("\n--- TEST 4: Separability test (age_group x X1_hat) ---\n")

df$X1_hat <- fitted(lm(X1 ~ support_group*age_group +
                         support_group*education +
                         support_group*prenatal_high +
                         education + age_group + prenatal_high, data=df))
df$X2_hat <- fitted(lm(X2 ~ support_group*age_group +
                         support_group*education +
                         support_group*prenatal_high +
                         education + age_group + prenatal_high, data=df))
df$X3_hat <- fitted(lm(X3 ~ support_group*age_group +
                         support_group*education +
                         support_group*prenatal_high +
                         education + age_group + prenatal_high, data=df))

ss_sep    <- lm(baby_weight_kg ~ X1_hat + X2_hat + X3_hat +
                  age_group + education + prenatal_high +
                  age_group:X1_hat, data=df)
ct_sep    <- coeftest(ss_sep, vcov=vcovHC(ss_sep, type="HC1"))
int_terms <- grep("age_group.*:X1_hat|X1_hat:age_group",
                  rownames(ct_sep), value=TRUE)

cat("Interaction terms (age_group x X1_hat):\n")
for (t in int_terms) {
  cat(" ", t, ":", round(ct_sep[t,"Estimate"],4),
      "(p =", round(ct_sep[t,"Pr(>|t|)"],3), ")\n")
}
cat("Significant interactions suggest separability may not hold.\n")

# ---- SEPARABILITY TABLE ----
sep_tbl <- do.call(rbind, lapply(int_terms, function(t) {
  data.frame(
    "Interaction term" = t,
    Estimate     = round(ct_sep[t,"Estimate"], 4),
    "Std. Error"  = round(ct_sep[t,"Std. Error"], 4),
    "t-stat"      = round(ct_sep[t,"t value"], 3),
    "P-value"     = format_p(ct_sep[t,"Pr(>|t|)"]),
    check.names   = FALSE,
    stringsAsFactors = FALSE
  )
}))
sep_tbl[,"Interaction term"] <- c(
  "Age group 21--25 $\\times$ $\\hat{X}_1$",
  "Age group 26--30 $\\times$ $\\hat{X}_1$",
  "Age group 31--40 $\\times$ $\\hat{X}_1$"
)
rownames(sep_tbl) <- NULL

write.csv(sep_tbl, "output/tables/table_separability.csv", row.names=FALSE)
save_tex(sep_tbl,
         filename = "table_separability.tex",
         caption  = paste(
           "Separability Test: Interactions between Age Group and $\\hat{X}_1$.",
           "The structural equation is augmented with interactions between",
           "age group and the fitted value of $\\hat{X}_1$.",
           "Under the separability assumption of Caetano and Escanciano (2021),",
           "these coefficients should be jointly zero.",
           "Statistically significant interactions suggest the effect of starting",
           "to smoke may vary with the mother's age group.",
           "Robust standard errors (HC1) in parentheses."
         ),
         label = "tab:separability")

# ---- TESTS SUMMARY TABLE (Table 6) ----
cat("\n--- TABLE 6: Tests Summary ---\n")

tbl6 <- data.frame(
  Test         = c("$H_0: \\beta_1 = \\beta_2$",
                   "$H_0: \\beta_2 = \\beta_3$",
                   "$H_0: \\beta_1 = \\beta_2 = \\beta_3$"),
  Statistic    = c(round(t12,3), round(t23,3), round(F_stat,3)),
  Distribution = c("$t$","$t$","$F(2, N-k)$"),
  "P-value"    = c(format_p(p12), format_p(p23), format_p(p_F)),
  Conclusion   = c(
    ifelse(p12<0.05,"Reject $H_0$","Fail to reject"),
    ifelse(p23<0.05,"Reject $H_0$","Fail to reject"),
    ifelse(p_F<0.05,"Reject $H_0$","Fail to reject")
  ),
  check.names      = FALSE,
  stringsAsFactors = FALSE
)

print(tbl6)
write.csv(tbl6, "output/tables/table6_tests.csv", row.names=FALSE)
save_tex(tbl6,
         filename = "table6_tests.tex",
         caption  = paste(
           "Tests for Constant Marginal Effects.",
           "All tests based on CE-TSLS estimates from Table~\\ref{tab:main_results}.",
           "Tests 1 and 2 are t-tests on pairwise differences of marginal effects,",
           "using the joint variance-covariance matrix of the second-stage estimator.",
           "Test 3 is a joint Wald F-test of $\\beta_1 = \\beta_2 = \\beta_3$.",
           "All three tests reject the null of constant marginal effects",
           "at the 5\\% significance level."
         ),
         label = "tab:tests")

# ============================================================
# PART 2: ROBUSTNESS CHECK — PRENATAL VISITS
# ============================================================
# Motivation: balance test shows Z increases prenatal visits
# (p < 0.001), suggesting it may be a mediator (bad control).
# We test sensitivity by running all IV specs with and without
# prenatal visits as a control.
# ============================================================

cat("\n\n============================================================\n")
cat(" PART 2: ROBUSTNESS CHECK — PRENATAL VISITS\n")
cat("============================================================\n")

# Simple IV without prenatal
fs_no    <- lm(smoked_trimesters ~ support_group + education + age_group,
               data=df)
df$S_no  <- fitted(fs_no)
ss_no    <- lm(baby_weight_kg ~ S_no + education + age_group, data=df)
ct_iv_no <- coeftest(ss_no, vcov=vcovHC(ss_no, type="HC1"))
cat("Simple IV WITHOUT prenatal:", round(ct_iv_no["S_no","Estimate"],4), "\n")

# Simple IV with prenatal
fs_yes    <- lm(smoked_trimesters ~ support_group + education +
                  age_group + prenatal_high, data=df)
df$S_yes  <- fitted(fs_yes)
ss_yes    <- lm(baby_weight_kg ~ S_yes + education +
                  age_group + prenatal_high, data=df)
ct_iv_yes <- coeftest(ss_yes, vcov=vcovHC(ss_yes, type="HC1"))
cat("Simple IV WITH prenatal:   ", round(ct_iv_yes["S_yes","Estimate"],4), "\n")

# CE-TSLS without prenatal
fs1n <- lm(X1 ~ support_group*age_group + support_group*education +
             education + age_group, data=df)
fs2n <- lm(X2 ~ support_group*age_group + support_group*education +
             education + age_group, data=df)
fs3n <- lm(X3 ~ support_group*age_group + support_group*education +
             education + age_group, data=df)
df$X1n <- fitted(fs1n); df$X2n <- fitted(fs2n); df$X3n <- fitted(fs3n)
ss_cen <- lm(baby_weight_kg ~ X1n + X2n + X3n + education + age_group,
             data=df)
ct_cen <- coeftest(ss_cen, vcov=vcovHC(ss_cen, type="HC1"))

cat("CE-TSLS WITHOUT prenatal: beta1=", round(ct_cen["X1n","Estimate"],4),
    "beta2=", round(ct_cen["X2n","Estimate"],4),
    "beta3=", round(ct_cen["X3n","Estimate"],4), "\n")

# CE-TSLS with prenatal
fs1y <- lm(X1 ~ support_group*age_group + support_group*education +
             support_group*prenatal_high +
             education + age_group + prenatal_high, data=df)
fs2y <- lm(X2 ~ support_group*age_group + support_group*education +
             support_group*prenatal_high +
             education + age_group + prenatal_high, data=df)
fs3y <- lm(X3 ~ support_group*age_group + support_group*education +
             support_group*prenatal_high +
             education + age_group + prenatal_high, data=df)
df$X1y <- fitted(fs1y); df$X2y <- fitted(fs2y); df$X3y <- fitted(fs3y)
ss_cey <- lm(baby_weight_kg ~ X1y + X2y + X3y +
               education + age_group + prenatal_high, data=df)
ct_cey <- coeftest(ss_cey, vcov=vcovHC(ss_cey, type="HC1"))

cat("CE-TSLS WITH prenatal:    beta1=", round(ct_cey["X1y","Estimate"],4),
    "beta2=", round(ct_cey["X2y","Estimate"],4),
    "beta3=", round(ct_cey["X3y","Estimate"],4), "\n")

# ---- ROBUSTNESS TABLE (Table 7) ----
cat("\n--- TABLE 7: Robustness Check ---\n")

tbl7 <- data.frame(
  Variable = c(
    "$\\beta_1$: smoked $\\geq 1$ trimester",
    "$\\beta_2$: smoked $\\geq 2$ trimesters",
    "$\\beta_3$: smoked $\\geq 3$ trimesters",
    "Smoked trimesters (level)",
    "Education: High School","Education: College","Education: Masters",
    "Age: 21--25","Age: 26--30","Age: 31--40",
    "Prenatal visits (high)"
  ),
  "Simple IV (1)" = c(
    "","","",
    get_val(ct_iv_no,"S_no"),
    get_val(ct_iv_no,"educationHigh School"),
    get_val(ct_iv_no,"educationCollege"),
    get_val(ct_iv_no,"educationMasters"),
    get_val(ct_iv_no,"age_group21-25"),
    get_val(ct_iv_no,"age_group26-30"),
    get_val(ct_iv_no,"age_group31-40"), ""
  ),
  "Simple IV (2)" = c(
    "","","",
    get_val(ct_iv_yes,"S_yes"),
    get_val(ct_iv_yes,"educationHigh School"),
    get_val(ct_iv_yes,"educationCollege"),
    get_val(ct_iv_yes,"educationMasters"),
    get_val(ct_iv_yes,"age_group21-25"),
    get_val(ct_iv_yes,"age_group26-30"),
    get_val(ct_iv_yes,"age_group31-40"),
    get_val(ct_iv_yes,"prenatal_high")
  ),
  "CE-TSLS (3)" = c(
    get_val(ct_cen,"X1n"), get_val(ct_cen,"X2n"), get_val(ct_cen,"X3n"), "",
    get_val(ct_cen,"educationHigh School"),
    get_val(ct_cen,"educationCollege"),
    get_val(ct_cen,"educationMasters"),
    get_val(ct_cen,"age_group21-25"),
    get_val(ct_cen,"age_group26-30"),
    get_val(ct_cen,"age_group31-40"), ""
  ),
  "CE-TSLS (4)" = c(
    get_val(ct_cey,"X1y"), get_val(ct_cey,"X2y"), get_val(ct_cey,"X3y"), "",
    get_val(ct_cey,"educationHigh School"),
    get_val(ct_cey,"educationCollege"),
    get_val(ct_cey,"educationMasters"),
    get_val(ct_cey,"age_group21-25"),
    get_val(ct_cey,"age_group26-30"),
    get_val(ct_cey,"age_group31-40"),
    get_val(ct_cey,"prenatal_high")
  ),
  check.names      = FALSE,
  stringsAsFactors = FALSE
)

footer7 <- data.frame(
  Variable        = c("Prenatal control","N"),
  "Simple IV (1)" = c("No",  format(nrow(df),big.mark=",")),
  "Simple IV (2)" = c("Yes", format(nrow(df),big.mark=",")),
  "CE-TSLS (3)"   = c("No",  format(nrow(df),big.mark=",")),
  "CE-TSLS (4)"   = c("Yes", format(nrow(df),big.mark=",")),
  check.names      = FALSE,
  stringsAsFactors = FALSE
)
colnames(footer7) <- colnames(tbl7)
tbl7 <- rbind(tbl7, footer7)
rownames(tbl7) <- NULL

print(tbl7)
write.csv(tbl7, "output/tables/table7_robustness.csv", row.names=FALSE)
save_tex(tbl7,
         filename = "table7_robustness.tex",
         caption  = paste(
           "Robustness Check: With and Without Prenatal Visits as Control.",
           "Dependent variable: birth weight in kilograms.",
           "Columns (1) and (3) exclude prenatal visits from controls.",
           "Columns (2) and (4) include prenatal visits as control.",
           "The balance test (Table~\\ref{tab:balance}) shows the support group",
           "significantly increases prenatal visits, suggesting it may be a",
           "mediator rather than an exogenous covariate.",
           "The stability of estimates across columns confirms that heterogeneous",
           "marginal effects hold regardless of how prenatal visits are treated.",
           "Robust standard errors (HC1) in parentheses.",
           "Baseline: Middle School, age 15--20."
         ),
         label = "tab:robustness")

cat("\n--- 03_tests.R complete ---\n")
