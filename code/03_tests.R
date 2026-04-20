# ============================================================
# Empirical Work — Econometrics I (EESP, 2026)
# 03_tests.R
#
# Contents:
#   1. H0: beta1 = beta2  (t-test)
#   2. H0: beta2 = beta3  (t-test)
#   3. H0: beta1 = beta2 = beta3  (joint Wald F-test)
#   4. Separability test: age_group x X1_hat interaction
#
# Loads CE-TSLS results saved by 02_iv_estimation.R
# ============================================================

library(lmtest)
library(sandwich)
library(xtable)

# ---- Load cleaned data and CE-TSLS results ----
df      <- readRDS("output/df_clean.rds")
results <- readRDS("output/ce_tsls_results.rds")

ct_ce <- results$ct_ce
V     <- results$V
b1    <- results$b1
b2    <- results$b2
b3    <- results$b3
ss_ce <- results$ss_ce

dof <- nrow(df) - ncol(model.matrix(ss_ce))

# ---- Helper ----
save_tex <- function(tbl, filename, caption, label) {
  xt <- xtable(tbl, caption=caption, label=label)
  print(xt,
        include.rownames       = FALSE,
        caption.placement      = "top",
        booktabs               = TRUE,
        sanitize.text.function = identity,
        file = file.path("output/tables", filename))
  cat("Saved:", filename, "\n")
}

cat("\nCE-TSLS estimates loaded:\n")
cat("  beta1:", round(b1,4), "\n")
cat("  beta2:", round(b2,4), "\n")
cat("  beta3:", round(b3,4), "\n")

# ============================================================
# TEST 1: H0: beta1 = beta2
# ============================================================
cat("\n--- TEST 1: H0: beta1 = beta2 ---\n")

var_diff12 <- V["X1_hat","X1_hat"] + V["X2_hat","X2_hat"] -
              2*V["X1_hat","X2_hat"]
t12 <- (b1 - b2) / sqrt(var_diff12)
p12 <- 2 * pt(-abs(t12), df=dof)

cat("Difference (b1-b2):", round(b1-b2, 4), "\n")
cat("SE(b1-b2):         ", round(sqrt(var_diff12), 4), "\n")
cat("t-statistic:       ", round(t12, 3), "\n")
cat("p-value:           ", round(p12, 4), "\n")

# ============================================================
# TEST 2: H0: beta2 = beta3
# ============================================================
cat("\n--- TEST 2: H0: beta2 = beta3 ---\n")

var_diff23 <- V["X2_hat","X2_hat"] + V["X3_hat","X3_hat"] -
              2*V["X2_hat","X3_hat"]
t23 <- (b2 - b3) / sqrt(var_diff23)
p23 <- 2 * pt(-abs(t23), df=dof)

cat("Difference (b2-b3):", round(b2-b3, 4), "\n")
cat("SE(b2-b3):         ", round(sqrt(var_diff23), 4), "\n")
cat("t-statistic:       ", round(t23, 3), "\n")
cat("p-value:           ", round(p23, 4), "\n")

# ============================================================
# TEST 3: H0: beta1 = beta2 = beta3 (joint Wald F-test)
# R = [1 -1  0]   such that R*[b1,b2,b3]' = 0 under H0
#     [0  1 -1]
# ============================================================
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

# ============================================================
# TEST 4: Separability test (age_group x X1_hat)
# ============================================================
cat("\n--- TEST 4: Separability test ---\n")

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

# ============================================================
# TESTS SUMMARY TABLE (Table 6)
# ============================================================
cat("\n--- TABLE 6: Tests Summary ---\n")

tbl6 <- data.frame(
  Test         = c("$H_0: \\beta_1 = \\beta_2$",
                   "$H_0: \\beta_2 = \\beta_3$",
                   "$H_0: \\beta_1 = \\beta_2 = \\beta_3$"),
  Statistic    = c(round(t12,3), round(t23,3), round(F_stat,3)),
  Distribution = c("$t$","$t$","$F(2, N-k)$"),
  "P-value"    = c(round(p12,4), round(p23,4), round(p_F,4)),
  Conclusion   = c(
    ifelse(p12<0.05,"Reject $H_0$","Fail to reject"),
    ifelse(p23<0.05,"Reject $H_0$","Fail to reject"),
    ifelse(p_F<0.05,"Reject $H_0$","Fail to reject")
  ),
  check.names  = FALSE,
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

cat("\n--- 03_tests.R complete ---\n")
