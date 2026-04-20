# ============================================================
# Empirical Work — Econometrics I (EESP, 2026)
# 02_firststage.R — First Stage and Reduced Form
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

# ---- Helper: extract clean regression table ----
reg_table <- function(model, se_type = "HC1") {
  coefs  <- coeftest(model, vcov = vcovHC(model, type = se_type))
  data.frame(
    Coefficient = rownames(coefs),
    Estimate    = round(coefs[,1], 4),
    Std_Error   = round(coefs[,2], 4),
    T_stat      = round(coefs[,3], 3),
    P_value     = round(coefs[,4], 3)
  )
}

# ============================================================
# FIRST STAGE 1: Effect of Z on smoked_trimesters (level)
# For the simple IV (constant effects) model
# ============================================================
cat("\n--- FIRST STAGE 1: Z on smoked_trimesters ---\n")

fs1 <- lm(smoked_trimesters ~ support_group + education + age_group + prenatal_high,
           data = df)

fs1_tbl <- reg_table(fs1)
cat("\nFirst stage (outcome: smoked trimesters):\n")
print(fs1_tbl)

# F-statistic for instrument strength
fs1_robust <- coeftest(fs1, vcov = vcovHC(fs1, type = "HC1"))
fs1_fstat  <- (fs1_robust["support_group", "Estimate"] /
               fs1_robust["support_group", "Std. Error"])^2
cat("\nF-statistic for Z (support_group):", round(fs1_fstat, 2), "\n")
cat("(Rule of thumb: F > 10 indicates strong instrument)\n")

# ============================================================
# FIRST STAGE 2: Effect of Z on X1, X2, X3 separately
# For the CE-TSLS model
# ============================================================
cat("\n--- FIRST STAGE 2: Z on X1, X2, X3 separately ---\n")

fs_X1 <- lm(X1 ~ support_group + education + age_group + prenatal_high, data = df)
fs_X2 <- lm(X2 ~ support_group + education + age_group + prenatal_high, data = df)
fs_X3 <- lm(X3 ~ support_group + education + age_group + prenatal_high, data = df)

for (nm in c("X1","X2","X3")) {
  mdl <- get(paste0("fs_", nm))
  cat("\nFirst stage for", nm, ":\n")
  ct  <- coeftest(mdl, vcov = vcovHC(mdl, type = "HC1"))
  cat("  Z coefficient:", round(ct["support_group","Estimate"], 4),
      " SE:", round(ct["support_group","Std. Error"], 4),
      " p:", round(ct["support_group","Pr(>|t|)"], 3), "\n")
}

# ============================================================
# FIRST STAGE 3: Heterogeneity — Z x age_group interactions
# This documents the covariance completeness condition
# ============================================================
cat("\n--- FIRST STAGE 3: Heterogeneity in first stage ---\n")

# Run first stage with Z x age_group interaction for smoked_trimesters
fs_het <- lm(smoked_trimesters ~ support_group * age_group + education + prenatal_high,
             data = df)
ct_het <- coeftest(fs_het, vcov = vcovHC(fs_het, type = "HC1"))

# Extract Z effect by age group (main + interaction)
cat("\nEffect of Z by age group (main effect + interaction):\n")
age_levels <- levels(df$age_group)
baseline   <- ct_het["support_group", "Estimate"]

effects <- data.frame(
  Age_group = age_levels,
  Z_effect  = c(
    baseline,
    baseline + ct_het["support_group:age_group21-25", "Estimate"],
    baseline + ct_het["support_group:age_group26-30", "Estimate"],
    baseline + ct_het["support_group:age_group31-40", "Estimate"]
  )
)
effects$Z_effect <- round(effects$Z_effect, 4)
print(effects)

cat("\nNote: Variation in Z effect across age groups supports covariance completeness.\n")

write.csv(effects, "output/tables/table_fs_heterogeneity.csv", row.names = FALSE)
save_tex(effects,
         filename = "table_fs_heterogeneity.tex",
         caption  = paste("First Stage Heterogeneity: Effect of Support Group",
                          "Assignment (Z) on Trimesters Smoked by Age Group.",
                          "Estimates from a regression of smoked trimesters on Z,",
                          "Z $\\times$ age group interactions, education, and prenatal visits.",
                          "Variation across age groups supports the covariance",
                          "completeness condition of Caetano and Escanciano (2021)."),
         label    = "tab:fs_het")

# ============================================================
# REDUCED FORM: Effect of Z on baby_weight_kg
# ============================================================
cat("\n--- REDUCED FORM: Z on baby_weight_kg ---\n")

rf <- lm(baby_weight_kg ~ support_group + education + age_group + prenatal_high,
         data = df)

rf_tbl <- reg_table(rf)
cat("\nReduced form (outcome: baby weight):\n")
print(rf_tbl)

rf_coef <- coeftest(rf, vcov = vcovHC(rf, type = "HC1"))
cat("\nReduced form estimate of Z:", round(rf_coef["support_group","Estimate"], 4))
cat("\n(This is the ITT: total effect of support group on birth weight)\n")

# ============================================================
# BUILD SUMMARY TABLE: First stage + Reduced form side by side
# ============================================================
cat("\n--- SUMMARY TABLE: First Stage and Reduced Form ---\n")

# Extract key rows from both models
vars_to_show <- c("support_group",
                  "educationHigh School", "educationCollege", "educationMasters",
                  "age_group21-25", "age_group26-30", "age_group31-40",
                  "prenatal_high")

labels <- c("Support group (Z)",
            "Education: High School", "Education: College", "Education: Masters",
            "Age: 21-25", "Age: 26-30", "Age: 31-40",
            "Prenatal visits (high)")

extract_row <- function(model, varname, se_type = "HC1") {
  ct <- coeftest(model, vcov = vcovHC(model, type = se_type))
  if (!varname %in% rownames(ct)) return(c(NA, NA))
  c(round(ct[varname, "Estimate"], 4),
    round(ct[varname, "Std. Error"], 4))
}

summary_tbl <- data.frame(
  Variable       = labels,
  FS_Estimate    = sapply(vars_to_show, function(v) extract_row(fs1, v)[1]),
  FS_SE          = sapply(vars_to_show, function(v) extract_row(fs1, v)[2]),
  RF_Estimate    = sapply(vars_to_show, function(v) extract_row(rf,  v)[1]),
  RF_SE          = sapply(vars_to_show, function(v) extract_row(rf,  v)[2])
)
rownames(summary_tbl) <- NULL

# Format as "estimate (SE)"
summary_tbl$`First Stage`  <- paste0(summary_tbl$FS_Estimate,
                                      " (", summary_tbl$FS_SE, ")")
summary_tbl$`Reduced Form` <- paste0(summary_tbl$RF_Estimate,
                                      " (", summary_tbl$RF_SE, ")")

out_tbl <- summary_tbl[, c("Variable","First Stage","Reduced Form")]

# Add N and R-squared rows
n_row  <- data.frame(Variable = "N",
                     "First Stage"  = format(nrow(df), big.mark=","),
                     "Reduced Form" = format(nrow(df), big.mark=","),
                     check.names = FALSE)
r2_row <- data.frame(Variable = "R-squared",
                     "First Stage"  = round(summary(fs1)$r.squared, 3),
                     "Reduced Form" = round(summary(rf)$r.squared,  3),
                     check.names = FALSE)

out_tbl <- rbind(out_tbl, n_row, r2_row)

cat("\nFirst Stage and Reduced Form:\n")
print(out_tbl)
write.csv(out_tbl, "output/tables/table3_fs_rf.csv", row.names = FALSE)
save_tex(out_tbl,
         filename = "table3_fs_rf.tex",
         caption  = paste("First Stage and Reduced Form Estimates.",
                          "Dependent variables: trimesters smoked (first stage)",
                          "and birth weight in kg (reduced form).",
                          "Heteroskedasticity-robust standard errors in parentheses.",
                          "Baseline categories: Middle School education, age group 15--20,",
                          "low prenatal visits (8 visits)."),
         label    = "tab:fs_rf")

cat("\n--- 02_firststage.R complete ---\n")
