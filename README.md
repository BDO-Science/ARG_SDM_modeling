# SalmonCountR Analysis & Shiny Application

This repository contains end-to-end R scripts and a Shiny application to simulate American River Fall‑run Chinook salmon spawner dynamics under alternative power bypass scenarios. From precomputing temperature‑dependent mortality (TDM) and life‑cycle calibrations to interactive exploration in Shiny, the workflow comprises four key scripts:

```
precompute.R    # Data loading, TDM simulations, calibration & forecast (run_analysis.R)
functions.R     # Core TDM and life‑cycle modeling functions
global.R        # Shiny data preload, constants, and helper lookups
app.R           # Shiny UI and server logic
```

---

## Prerequisites

* **R** ≥ 4.0
* **Packages**:

  ```r
  tidyverse
  lubridate
  furrr
  data.table
  compiler
  shiny
  DT
  ```
* **Directory structure**:

  ```
  SalmonCountR/
  ├── precompute.R
  ├── functions.R
  ├── global.R
  ├── app.R
  └── app_data/
      ├── env_ext_list.rds
      ├── df_all.rds
      ├── carcassdet_*.csv
      ├── grandtab_*.csv
      ├── egg_summary.rds
      ├── surv_lookup_full.rds
      ├── base_P_list.rds
      ├── S_seed.rds
      ├── stoch_SAR_opts.rds
      ├── spawn_dates_vec.rds
      └── american_river_instream.rda
  ```

---

## Script: precompute.R

**Purpose:**

1. Clear workspace and load libraries/helper functions (`functions.R`).
2. Read input data: environmental time series, carcass surveys, escapement estimates.
3. Define TDM variants (`exp_WF`, `exp_SM`, `lin_Martin`).
4. Build redd simulation records for observed years (2011–2024) and bootstrap future years (2025–2060).
5. Compile and run TDM per redd × environment in parallel (`furrr`).
6. Calibrate the life‑cycle model against observed escapement (2011–2024) via `optim()`.
7. Specify stochastic SAR options and generate SAR time series.
8. Forecast full life‑cycle to 2060 under stochastic SAR.
9. Save intermediate and final objects to `app_data/` for Shiny.

**Usage:**

```r
# In R at project root
tidyverse::walk(c("tidyverse","lubridate","furrr","data.table","compiler"), library, character.only = TRUE)
source("precompute.R")
```

---

## Script: functions.R

**Purpose:** Defines all modeling functions sourced by both `precompute.R` and Shiny:

* **egg\_model(T):** days to hatch (958 ATU).
* **hatch\_model(T):** days to emergence (417 ATU).
* **tdm\_exp(temps, calib):** exponential TDM (WaterForum2020 or SALMOD2006).
* **tdm\_lin\_martin(temps, α, β):** linear TDM (Martin et al. 2017).
* **compute\_surv(rdr, ...):** wrapper to extract temperature slice and compute survival.
* **surv\_adult\_prespawn(deg\_day, intercept, beta):** logistic adult pre‑spawn survival.
* **compute\_deg\_day\_adult(...):** cumulative degree‑day calculation per brood year.
* **simulate\_variant(...):** full age‑structured life‑cycle simulation (egg → adult).
* **generate\_SAR\_vec(n\_years, opts):** SAR time series generator with timing options.

**Usage:**

```r
# In any script
eval(parse("functions.R"))
```

---

## Script: global.R

**Purpose:** Preload data and define constants for the Shiny app:

1. Load required libraries and `functions.R`.
2. Read all precomputed `.rds` and `.rda` data (`env_ext_list`, `egg_summary`, etc.).
3. Define simulation constants (`n_sim`, `real_years`, `sim_years_full`).
4. Compute instream habitat lookup and expose `get_K_spawners(flow)`.
5. Expose UI constants: `env_levels`, default `n_sim`.

**Usage:**

```r
# Called at top of app.R
source("global.R")
```

---

## Script: app.R

**Purpose:** Interactive Shiny application allowing users to explore and compare scenarios:

* **UI (navbarPage):** Tabs for About, Single Alternative, Compare Alternatives, and Settings.
* **Inputs:** Simulation length, mode (Deterministic/Stochastic), TDM variant, alternative(s), flow, SAR distribution parameters.
* **Server logic:**

  * `observeEvent(run_sim)`: runs single-scenario simulation and stores in `sim_data`.
  * `observeEvent(run_cmp)`: runs multi-alternative comparison and stores in `cmp_data`.
  * `renderDT` tables with formatted outputs and calibration shading.
  * `renderPlot` for time series, distributions, heatmaps, boxplots, and SAR series using ggplot2 & Viridis.
* **Launch:** `shiny::runApp('app.R')` or via RStudio “Run App.”

---

## Notes & Conventions

* **Parallelization:** TDM simulations use up to 6 workers via `furrr`.
* **Calibration window:** brood years 2011–2024 (first 14 years).
* **File organization:** All R scripts at project root; precomputed data in `app_data/`.
* **Dependencies:** ensure R packages installed before running any script.

---

*End of consolidated README*
