# Empirical Work — Econometrics I (EESP, 2026)

## Overview

This repository contains the code and outputs for the Empirical Work assignment
of the Econometrics I course (Bruno Ferman, EESP, Q1 2026).

The analysis estimates the marginal effect of maternal smoking during pregnancy
on infant birth weight, using the method of Caetano and Escanciano (2021) to
recover heterogeneous marginal effects with a single binary instrument.

## Reference

Caetano, C. and J. C. Escanciano (2021). Identifying multiple marginal effects
with a single instrument. *Econometric Theory* 37(3), 464–494.

## Data

The dataset (`data/empirical_work_dataset.csv`) contains 440,856 birth records
from a simulated Brazilian city, with the following variables:

| Variable            | Description                                      |
|---------------------|--------------------------------------------------|
| `baby_weight_kg`    | Birth weight in kilograms (outcome)              |
| `smoked_trimesters` | Number of trimesters smoked (0–3) (treatment)    |
| `support_group`     | Offered participation in support group (instrument) |
| `education`         | Highest education level (control)                |
| `prenatal_visits`   | Number of prenatal visits: 8 or 16 (control)     |
| `age_group`         | Mother's age group (control)                     |

## Repository Structure

```
├── README.md
├── data/
│   └── empirical_work_dataset.csv
├── code/
│   ├── 00_master.R          # Master script — runs everything
│   ├── 01_descriptive.R     # Descriptive statistics and balance tests
│   ├── 02_firststage.R      # First stage and reduced form
│   ├── 03_regressions.R     # OLS, simple IV, and CE-TSLS
│   └── 04_tests.R           # Constant effects test, overidentification
├── output/
│   ├── figures/             # All figures (PNG)
│   └── tables/              # All tables (CSV/TXT)
└── report/
    └── empirical_work.pdf   # Final submitted report
```

## How to Reproduce

1. Clone the repository:
   ```
   git clone https://github.com/[your-username]/empirical-work-econometrics1
   ```

2. Open R and set the working directory to the project root:
   ```r
   setwd("/path/to/empirical-work-econometrics1")
   ```

3. Install required packages (first time only):
   ```r
   install.packages(c("AER", "dplyr", "ggplot2", "lmtest", "sandwich", "xtable"))
   ```

4. Run the master script:
   ```r
   source("code/00_master.R")
   ```

All figures and tables will be saved automatically to the `output/` folder.

## Software and Package Versions

| Software  | Version    |
|-----------|------------|
| R         | 4.3.3      |
| AER       | 1.2.12     |
| dplyr     | 1.1.4      |
| ggplot2   | 3.4.4      |
| lmtest    | 0.9.40     |
| sandwich  | 3.1.0      |
