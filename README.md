# Empirical Work — Econometrics I (EESP, 2026)

## Overview

This repository contains the replication code and outputs for the Empirical Work
assignment of the Econometrics I course (Prof. Bruno Ferman, EESP, Q1 2026).

The analysis estimates the marginal effect of maternal smoking during pregnancy
on infant birth weight, using the method of Caetano and Escanciano (2021) to
recover heterogeneous marginal effects with a single binary instrument.

The main finding is that marginal effects are strongly heterogeneous: smoking
during the first trimester reduces birth weight by approximately 0.57 kg, while
smoking during the third trimester reduces it by only 0.10 kg. The null of
constant marginal effects is rejected at all conventional significance levels
(joint F-statistic = 52.49, p < 0.001).

## Reference

Caetano, C. and J. C. Escanciano (2021). Identifying multiple marginal effects
with a single instrument. *Econometric Theory* 37(3), 464–494.
Joshua D. Angrist and Jorn-Steffen Pischke. Mostly Harmless Econometrics. Princeton Uni-
versity Press, 2009.

## Data

The dataset (`data/empirical_work_dataset.csv`) contains 440,856 birth
records from a random Brazilian city, with the following variables:

| Variable            | Role          | Description                                        |
|---------------------|---------------|----------------------------------------------------|
| `baby_weight_kg`    | Outcome (Y)   | Birth weight in kilograms                          |
| `smoked_trimesters` | Treatment (X) | Number of trimesters smoked during pregnancy (0–3) |
| `support_group`     | Instrument (Z)| Randomly offered participation in support group    |
| `education`         | Control (W)   | Highest education level                            |
| `age_group`         | Control (W)   | Mother's age group                                 |
| `prenatal_visits`   | Control (W)   | Number of prenatal visits: 8 or 16 (binary)        |

## Repository Structure

```
├── empirical-work-econometrics1.Rproj   # RStudio project file
├── README.md
├── renv.lock                            # Package version lockfile
├── .Rprofile                            # Activates renv on project open
├── data/
│   └── empirical_work_dataset.csv
├── code/
│   ├── 00_master.R                      # Master script — runs everything
│   ├── 01_data_and_benchmarks.R         # Descriptive stats, balance check,
│   │                                    # first stage + reduced form, OLS
│   ├── 02_iv_estimation.R               # Simple IV and CE-TSLS (main results)
│   └── 03_tests.R                       # All tests: constant effects,
│                                        # separability, robustness check
├── output/
│   ├── figures/
│   │   ├── figure1_birthweight_by_smoking.png
│   │   └── figure2_firststage_heterogeneity.png
│   ├── tables/
│   │   ├── table1_descriptive.tex       # Descriptive statistics
│   │   ├── table2_balance.tex           # Balance check
│   │   ├── table3_fs_rf.tex             # First stage + reduced form
│   │   ├── table4_ols.tex               # OLS benchmarks
│   │   ├── table5_main_results.tex      # Main results: OLS vs IV vs CE-TSLS
│   │   ├── table6_tests.tex             # Tests for constant marginal effects
│   │   └── table7_robustness.tex        # Robustness: prenatal visits
│   └── ce_tsls_results.rds             # Saved CE-TSLS results (intermediate)
└── report/
    └── empirical_work.pdf               # Final submitted report
```

## How to Reproduce

### First time setup

1. Clone the repository:
   ```bash
   git clone https://github.com/[your-username]/empirical-work-econometrics1
   ```

2. Open `empirical-work-econometrics1.Rproj` in RStudio.
   This automatically sets the working directory to the project root.

3. Restore exact package versions using `renv`:
   ```r
   renv::restore()
   ```

4. Run the master script:
   ```r
   source("code/00_master.R")
   ```

### Without renv (alternative)

```r
install.packages(c("AER", "dplyr", "ggplot2", "lmtest", "sandwich", "xtable"))
source("code/00_master.R")
```

Note that different package versions may produce slightly different results.

## Output Summary

| Script | Tables produced | Figures produced |
|--------|----------------|-----------------|
| `01_data_and_benchmarks.R` | Tables 3–7 | Figures 1–2 |
| `02_iv_estimation.R` | Table 8 | — |
| `03_tests.R` | Tables 9-10 | — |

## Implementation Notes

- **CE-TSLS** is implemented manually as two-stage OLS rather than using
  `ivreg()` from the AER package. `ivreg()` with `vcovHC()` on the full
  dataset (440k obs) exceeds memory limits. The manual approach yields
  identical point estimates. Standard errors are HC1-robust second-stage SEs.

- **Instruments** for CE-TSLS first stages: `support_group` (Z) and all
  interactions `Z × age_group`, `Z × education`, `Z × prenatal_high`.

- **Prenatal visits** takes only two values (8 or 16) and is recoded as a
  binary dummy (`prenatal_high = 1` if 16 visits).

- **CE-TSLS results** are saved to `output/ce_tsls_results.rds` by
  `02_iv_estimation.R` and loaded by `03_tests.R` to avoid recomputation.

- **Robustness check** in `03_tests.R` runs all IV specifications with and
  without prenatal visits as a control, addressing the bad control concern
  raised by the balance test (Table 2).

## Software and Package Versions

Exact versions are recorded in `renv.lock`. Key packages:

| Software  | Version |
|-----------|---------|
| R         | 4.3.3   |
| AER       | 1.2.12  |
| dplyr     | 1.2.1   |
| ggplot2   | 4.0.2   |
| lmtest    | 0.9.40  |
| sandwich  | 3.1.0   |
| xtable    | 1.8.8   |
