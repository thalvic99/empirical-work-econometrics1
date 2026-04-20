# ============================================================
# Empirical Work — Econometrics I (EESP, 2026)
# 01_data_and_benchmarks.R
#
# Contents:
#   1. Descriptive statistics (Table 1)
#   2. Balance check (Table 2)
#   3. First stage on S, X1, X2, X3 + reduced form (Table 3)
#   4. First stage heterogeneity by age group (Figure 2)
#   5. OLS benchmarks (Table 4)
#   6. Figure 1: Birth weight by smoking level
# ============================================================

library(AER)
library(dplyr)
library(ggplot2)
library(lmtest)
library(sandwich)
library(xtable)

# ---- Load data ----
df <- read.csv("data/empirical_work_dataset.csv")

# ---- Recode variables ----
df$education <- factor(df$education,
                       levels = c("Middle School","High School","College","Masters"))
df$age_group  <- factor(df$age_group,
                        levels = c("15-20","21-25","26-30","31-40"))
df$prenatal_high <- as.integer(df$prenatal_visits == 16)
df$X1 <- as.integer(df$smoked_trimesters >= 1)
df$X2 <- as.integer(df$smoked_trimesters >= 2)
df$X3 <- as.integer(df$smoked_trimesters >= 3)

# Save cleaned data for subsequent scripts
saveRDS(df, "output/df_clean.rds")

# ---- Helper: save LaTeX table ----
save_tex <- function(tbl, filename, caption, label, digits=3) {
  xt <- xtable(tbl, caption=caption, label=label, digits=digits)
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

extract <- function(ct, v) {
  if (!v %in% rownames(ct)) return("")
  fmt(ct[v,"Estimate"], ct[v,"Std. Error"])
}

# ============================================================
# 1. DESCRIPTIVE STATISTICS (Table 1)
# ============================================================
cat("\n--- TABLE 1: Descriptive Statistics ---\n")

desc_vars <- list(
  "Baby weight (kg)"         = df$baby_weight_kg,
  "Smoked trimesters"        = df$smoked_trimesters,
  "Support group (Z)"        = df$support_group,
  "Prenatal visits (high)"   = df$prenatal_high,
  "Age: 15--20"              = as.integer(df$age_group == "15-20"),
  "Age: 21--25"              = as.integer(df$age_group == "21-25"),
  "Age: 26--30"              = as.integer(df$age_group == "26-30"),
  "Age: 31--40"              = as.integer(df$age_group == "31-40"),
  "Education: Middle School" = as.integer(df$education == "Middle School"),
  "Education: High School"   = as.integer(df$education == "High School"),
  "Education: College"       = as.integer(df$education == "College"),
  "Education: Masters"       = as.integer(df$education == "Masters")
)

tbl1 <- data.frame(
  Variable = names(desc_vars),
  Mean     = sapply(desc_vars, mean),
  SD       = sapply(desc_vars, sd),
  Min      = sapply(desc_vars, min),
  Max      = sapply(desc_vars, max)
)
rownames(tbl1) <- NULL
tbl1[,-1] <- round(tbl1[,-1], 3)
print(tbl1)

write.csv(tbl1, "output/tables/table1_descriptive.csv", row.names=FALSE)
save_tex(tbl1,
         filename = "table1_descriptive.tex",
         caption  = paste0("Descriptive Statistics (N = ",
                           format(nrow(df), big.mark=","), ")"),
         label    = "tab:descriptive")

# ============================================================
# 2. BALANCE CHECK (Table 2)
# ============================================================
cat("\n--- TABLE 2: Balance Test ---\n")

balance_vars <- list(
  "Prenatal visits (high)"   = df$prenatal_high,
  "Age: 15--20"              = as.integer(df$age_group == "15-20"),
  "Age: 21--25"              = as.integer(df$age_group == "21-25"),
  "Age: 26--30"              = as.integer(df$age_group == "26-30"),
  "Age: 31--40"              = as.integer(df$age_group == "31-40"),
  "Education: Middle School" = as.integer(df$education == "Middle School"),
  "Education: High School"   = as.integer(df$education == "High School"),
  "Education: College"       = as.integer(df$education == "College"),
  "Education: Masters"       = as.integer(df$education == "Masters")
)

bal_results <- lapply(names(balance_vars), function(vname) {
  y  <- balance_vars[[vname]]
  z  <- df$support_group
  m0 <- mean(y[z == 0])
  m1 <- mean(y[z == 1])
  tt <- t.test(y[z == 1], y[z == 0])
  
  # Format p-value
  pval <- tt$p.value
  pval_fmt <- if (pval < 0.001) "$<$0.001" else formatC(pval, digits=3, format="f")
  
  data.frame(
    Variable          = vname,
    "No support (Z=0)" = formatC(m0, digits=3, format="f"),
    "Support (Z=1)"    = formatC(m1, digits=3, format="f"),
    Difference        = formatC(m1 - m0, digits=3, format="f"),
    "P-value"         = pval_fmt,
    check.names       = FALSE
  )
})

tbl2 <- do.call(rbind, bal_results)
rownames(tbl2) <- NULL
print(tbl2)

write.csv(tbl2, "output/tables/table2_balance.csv", row.names=FALSE)
save_tex(tbl2,
         filename = "table2_balance.tex",
         caption  = paste("Balance Test: Pre-determined Characteristics by",
                          "Support Group Assignment.",
                          "P-values from two-sided t-tests."),
         label    = "tab:balance")

# ============================================================
# 3. FIRST STAGE + REDUCED FORM (Table 3)
# ============================================================
cat("\n--- TABLE 3: First Stage and Reduced Form ---\n")

# First stage: effect of Z on smoked_trimesters (level)
fs <- lm(smoked_trimesters ~ support_group + education +
           age_group + prenatal_high, data=df)
ct_fs <- coeftest(fs, vcov=vcovHC(fs, type="HC1"))

# F-statistic for instrument strength
fs_fstat <- (ct_fs["support_group","Estimate"] /
             ct_fs["support_group","Std. Error"])^2
cat("First stage F-statistic:", round(fs_fstat, 1), "\n")

# Reduced form: effect of Z on baby_weight_kg
rf <- lm(baby_weight_kg ~ support_group + education +
           age_group + prenatal_high, data=df)
ct_rf <- coeftest(rf, vcov=vcovHC(rf, type="HC1"))

# Wald estimate
wald <- ct_rf["support_group","Estimate"] / ct_fs["support_group","Estimate"]
cat("Wald (RF/FS):", round(wald, 4), "\n")

# Build table
vars_show <- c("support_group",
               "educationHigh School","educationCollege","educationMasters",
               "age_group21-25","age_group26-30","age_group31-40",
               "prenatal_high")
labels <- c("Support group (Z)",
            "Education: High School","Education: College","Education: Masters",
            "Age: 21--25","Age: 26--30","Age: 31--40",
            "Prenatal visits (high)")

tbl3 <- data.frame(
  Variable      = labels,
  "First Stage" = sapply(vars_show, function(v) extract(ct_fs, v)),
  "Reduced Form"= sapply(vars_show, function(v) extract(ct_rf, v)),
  check.names   = FALSE
)

footer3 <- data.frame(
  Variable       = c("N","R-squared"),
  "First Stage"  = c(format(nrow(df),big.mark=","),
                     round(summary(fs)$r.squared,3)),
  "Reduced Form" = c(format(nrow(df),big.mark=","),
                     round(summary(rf)$r.squared,3)),
  check.names    = FALSE
)
colnames(footer3) <- colnames(tbl3)
tbl3 <- rbind(tbl3, footer3)
rownames(tbl3) <- NULL
print(tbl3)

write.csv(tbl3, "output/tables/table3_fs_rf.csv", row.names=FALSE)
save_tex(tbl3,
         filename = "table3_fs_rf.tex",
         caption  = paste("First Stage and Reduced Form Estimates.",
                          "Dependent variables: trimesters smoked (first stage)",
                          "and birth weight in kg (reduced form).",
                          "Robust standard errors (HC1) in parentheses.",
                          "Baseline: Middle School, age 15--20, 8 prenatal visits."),
         label    = "tab:fs_rf")

# ============================================================
# 4. FIRST STAGE HETEROGENEITY BY AGE GROUP (Figure 2)
# Also documents rank condition for CE-TSLS
# ============================================================
cat("\n--- FIGURE 2: First Stage Heterogeneity ---\n")

# First stages on X1, X2, X3 — checks that Z affects each margin
for (xvar in c("X1","X2","X3")) {
  mdl <- lm(as.formula(paste(xvar, "~ support_group + education +
                              age_group + prenatal_high")), data=df)
  ct  <- coeftest(mdl, vcov=vcovHC(mdl, type="HC1"))
  cat(xvar, "— Z coef:", round(ct["support_group","Estimate"],4),
      " SE:", round(ct["support_group","Std. Error"],4), "\n")
}

# Heterogeneity: Z x age_group interaction
fs_het <- lm(smoked_trimesters ~ support_group * age_group +
               education + prenatal_high, data=df)
ct_het <- coeftest(fs_het, vcov=vcovHC(fs_het, type="HC1"))

baseline <- ct_het["support_group","Estimate"]
effects <- data.frame(
  Age_group = levels(df$age_group),
  Z_effect  = round(c(
    baseline,
    baseline + ct_het["support_group:age_group21-25","Estimate"],
    baseline + ct_het["support_group:age_group26-30","Estimate"],
    baseline + ct_het["support_group:age_group31-40","Estimate"]
  ), 4)
)
cat("First stage effect of Z by age group:\n")
print(effects)

write.csv(effects, "output/tables/table_fs_heterogeneity.csv", row.names=FALSE)
save_tex(effects,
         filename = "table_fs_heterogeneity.tex",
         caption  = paste("First Stage Heterogeneity by Age Group.",
                          "Effect of Z on trimesters smoked by age group,",
                          "from a regression with Z $\\times$ age group interactions.",
                          "Variation across groups supports covariance completeness."),
         label    = "tab:fs_het")

# Figure 2
fs_age <- df %>%
  group_by(age_group, support_group) %>%
  summarise(mean_S = mean(smoked_trimesters), .groups="drop") %>%
  mutate(support_group = ifelse(support_group==1,
                                "Support group (Z=1)",
                                "No support group (Z=0)"))

p2 <- ggplot(fs_age, aes(x=age_group, y=mean_S, fill=support_group)) +
  geom_col(position="dodge", width=0.6) +
  scale_fill_manual(values=c("steelblue","tomato")) +
  labs(
    title    = "Figure 2: First Stage Heterogeneity by Age Group",
    subtitle = "Mean trimesters smoked by support group assignment and age group",
    x        = "Age group",
    y        = "Mean trimesters smoked",
    fill     = "",
    caption  = paste("Note: Gap between bars within each age group is the first",
                     "stage effect of Z.\nHeterogeneity across groups supports",
                     "covariance completeness (Caetano & Escanciano, 2021).")
  ) +
  theme_bw(base_size=13) +
  theme(legend.position="bottom")

print(p2)

ggsave("output/figures/figure2_firststage_heterogeneity.png",
       p2, width=7, height=5, dpi=150)
cat("Figure 2 saved.\n")

# ============================================================
# 5. OLS BENCHMARKS (Table 4)
# ============================================================
cat("\n--- TABLE 4: OLS Benchmarks ---\n")

ols1 <- lm(baby_weight_kg ~ I(smoked_trimesters > 0), data=df)
ols2 <- lm(baby_weight_kg ~ smoked_trimesters, data=df)
ols3 <- lm(baby_weight_kg ~ smoked_trimesters + education +
             age_group + prenatal_high, data=df)
ols4 <- lm(baby_weight_kg ~ X1 + X2 + X3 + education +
             age_group + prenatal_high, data=df)

ct1 <- coeftest(ols1, vcov=vcovHC(ols1, type="HC1"))
ct2 <- coeftest(ols2, vcov=vcovHC(ols2, type="HC1"))
ct3 <- coeftest(ols3, vcov=vcovHC(ols3, type="HC1"))
ct4 <- coeftest(ols4, vcov=vcovHC(ols4, type="HC1"))

cat("OLS marginal effects (Column 4):\n")
cat("  beta1:", round(ct4["X1","Estimate"],4), "\n")
cat("  beta2:", round(ct4["X2","Estimate"],4), "\n")
cat("  beta3:", round(ct4["X3","Estimate"],4), "\n")

vars_ols  <- c("I(smoked_trimesters > 0)TRUE","smoked_trimesters",
               "X1","X2","X3",
               "educationHigh School","educationCollege","educationMasters",
               "age_group21-25","age_group26-30","age_group31-40",
               "prenatal_high")
labels_ols <- c("Smoked (binary)","Smoked trimesters (level)",
                "X1: smoked $\\geq$ 1 trimester",
                "X2: smoked $\\geq$ 2 trimesters",
                "X3: smoked $\\geq$ 3 trimesters",
                "Education: High School","Education: College","Education: Masters",
                "Age: 21--25","Age: 26--30","Age: 31--40",
                "Prenatal visits (high)")

cts <- list(ct1, ct2, ct3, ct4)
out_rows <- lapply(seq_along(vars_ols), function(i) {
  v   <- vars_ols[i]
  row <- data.frame(Variable=labels_ols[i])
  for (j in 1:4) {
    row[[paste0("(",j,")")]] <- extract(cts[[j]], v)
  }
  row
})
tbl4 <- do.call(rbind, out_rows)

footer4 <- data.frame(
  Variable = c("Controls","N","R-squared"),
  "(1)"    = c("No",  format(nrow(df),big.mark=","), formatC(summary(ols1)$r.squared, digits=3, format="f")),
  "(2)"    = c("No",  format(nrow(df),big.mark=","), formatC(summary(ols2)$r.squared, digits=3, format="f")),
  "(3)"    = c("Yes", format(nrow(df),big.mark=","), formatC(summary(ols3)$r.squared, digits=3, format="f")),
  "(4)"    = c("Yes", format(nrow(df),big.mark=","), formatC(summary(ols4)$r.squared, digits=3, format="f")),
  check.names = FALSE
)
colnames(footer4) <- colnames(tbl4)
tbl4 <- rbind(tbl4, footer4)
rownames(tbl4) <- NULL
print(tbl4)

write.csv(tbl4, "output/tables/table4_ols.csv", row.names=FALSE)
save_tex(tbl4,
         filename = "table4_ols.tex",
         caption  = paste("OLS Benchmarks: Effect of Smoking on Birth Weight.",
                          "Dependent variable: birth weight in kilograms.",
                          "(1) Binary smoking dummy, no controls.",
                          "(2) Smoking level, no controls.",
                          "(3) Smoking level, with controls.",
                          "(4) Three smoking dummies X1, X2, X3, with controls.",
                          "All OLS estimates are biased due to endogeneity of smoking.",
                          "Robust standard errors (HC1) in parentheses.",
                          "Baseline: Middle School, age 15--20, 8 prenatal visits."),
         label    = "tab:ols")

# ============================================================
# 6. FIGURE 1: Birth weight by smoking level
# ============================================================
cat("\n--- FIGURE 1: Birth Weight by Smoking Level ---\n")

means_by_s <- df %>%
  group_by(smoked_trimesters) %>%
  summarise(
    mean_weight = mean(baby_weight_kg),
    se          = sd(baby_weight_kg)/sqrt(n()),
    .groups     = "drop"
  )

p1 <- ggplot(means_by_s, aes(x=factor(smoked_trimesters), y=mean_weight)) +
  geom_col(fill="steelblue", width=0.6) +
  geom_errorbar(aes(ymin=mean_weight-1.96*se,
                    ymax=mean_weight+1.96*se), width=0.2) +
  labs(
    title    = "Figure 1: Mean Birth Weight by Trimesters Smoked",
    subtitle = "Error bars represent 95% confidence intervals",
    x        = "Trimesters smoked during pregnancy",
    y        = "Mean birth weight (kg)",
    caption  = "Note: Raw means, not controlling for confounders. N = 440,856."
  ) +
  theme_bw(base_size=13)

print(p1)

ggsave("output/figures/figure1_birthweight_by_smoking.png",
       p1, width=7, height=5, dpi=150)
cat("Figure 1 saved.\n")

cat("\n--- 01_data_and_benchmarks.R complete ---\n")
