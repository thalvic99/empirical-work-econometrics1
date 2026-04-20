# ============================================================
# Empirical Work — Econometrics I (EESP, 2026)
# 01_descriptive.R — Descriptive Statistics and Balance Tests
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

# Dummy endogenous variables (for CE-TSLS later)
df$X1 <- as.integer(df$smoked_trimesters >= 1)
df$X2 <- as.integer(df$smoked_trimesters >= 2)
df$X3 <- as.integer(df$smoked_trimesters >= 3)

# Save cleaned data object for use in subsequent scripts
saveRDS(df, "output/df_clean.rds")

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

# ============================================================
# TABLE 1: DESCRIPTIVE STATISTICS
# ============================================================
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

cat("\n--- TABLE 1: Descriptive Statistics (N =", nrow(df), ") ---\n")
print(tbl1)
write.csv(tbl1, "output/tables/table1_descriptive.csv", row.names = FALSE)
save_tex(tbl1,
         filename = "table1_descriptive.tex",
         caption  = paste0("Descriptive Statistics (N = ",
                           format(nrow(df), big.mark=","), ")"),
         label    = "tab:descriptive")

# ============================================================
# TABLE 2: BALANCE TEST
# ============================================================
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
  data.frame(
    Variable           = vname,
    "No support (Z=0)" = round(m0, 3),
    "Support (Z=1)"    = round(m1, 3),
    Difference         = round(m1 - m0, 3),
    "P-value"          = round(tt$p.value, 3),
    check.names        = FALSE
  )
})

tbl2 <- do.call(rbind, bal_results)
rownames(tbl2) <- NULL

cat("\n--- TABLE 2: Balance Test ---\n")
print(tbl2)
cat("Note: P-values from two-sided t-tests.\n")
write.csv(tbl2, "output/tables/table2_balance.csv", row.names = FALSE)
save_tex(tbl2,
         filename = "table2_balance.tex",
         caption  = paste("Balance Test: Pre-determined Characteristics by",
                          "Support Group Assignment.",
                          "P-values from two-sided t-tests.",
                          "No statistically significant differences support",
                          "random assignment of the instrument."),
         label    = "tab:balance")

# ============================================================
# FIGURE 1: Mean birth weight by smoking level
# ============================================================
means_by_s <- df %>%
  group_by(smoked_trimesters) %>%
  summarise(
    mean_weight = mean(baby_weight_kg),
    se          = sd(baby_weight_kg) / sqrt(n()),
    .groups     = "drop"
  )

p1 <- ggplot(means_by_s, aes(x = factor(smoked_trimesters), y = mean_weight)) +
  geom_col(fill = "steelblue", width = 0.6) +
  geom_errorbar(aes(ymin = mean_weight - 1.96*se,
                    ymax = mean_weight + 1.96*se), width = 0.2) +
  labs(
    title    = "Figure 1: Mean Birth Weight by Trimesters Smoked",
    subtitle = "Error bars represent 95% confidence intervals",
    x        = "Trimesters smoked during pregnancy",
    y        = "Mean birth weight (kg)",
    caption  = "Note: Raw means, not controlling for confounders. N = 440,856."
  ) +
  theme_bw(base_size = 13)

ggsave("output/figures/figure1_birthweight_by_smoking.png",
       p1, width = 7, height = 5, dpi = 150)
cat("Figure 1 saved.\n")

# ============================================================
# FIGURE 2: First stage heterogeneity by age group
# ============================================================
fs_age <- df %>%
  group_by(age_group, support_group) %>%
  summarise(mean_S = mean(smoked_trimesters), .groups = "drop") %>%
  mutate(support_group = ifelse(support_group == 1,
                                "Support group (Z=1)",
                                "No support group (Z=0)"))

p2 <- ggplot(fs_age, aes(x = age_group, y = mean_S, fill = support_group)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("steelblue","tomato")) +
  labs(
    title    = "Figure 2: First Stage Heterogeneity by Age Group",
    subtitle = "Mean trimesters smoked by support group assignment and age group",
    x        = "Age group",
    y        = "Mean trimesters smoked",
    fill     = "",
    caption  = paste("Note: The gap between bars within each age group represents",
                     "the first stage effect of Z.",
                     "\nHeterogeneity across age groups drives identification",
                     "in Caetano & Escanciano (2021).")
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave("output/figures/figure2_firststage_heterogeneity.png",
       p2, width = 7, height = 5, dpi = 150)
cat("Figure 2 saved.\n")

cat("\n--- 01_descriptive.R complete ---\n")
