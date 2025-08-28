# SalmonCountR Analysis & Shiny Application

This repository contains end-to-end R scripts and a Shiny application to simulate American River Fall-run Chinook salmon spawner dynamics under alternative (env in the code) power bypass scenarios. From precomputing temperature-dependent mortality (TDM) and life-cycle calibrations to interactive exploration in Shiny, the workflow comprises four key scripts:

```
precompute.R    # Data loading, TDM simulations, calibration & forecast
functions.R     # Core TDM and life-cycle modeling functions
global.R        # Shiny data preload, constants, and helper lookups
app.R           # Shiny UI and server logic
```

---

## 🚀 Quickstart

```r
# From the SalmonCountR project root in R:
source("precompute.R")       # Precompute survivals, calibration, forecasts
shiny::runApp("app.R")       # Launch the Shiny app
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
      ├── calib_results.rds
      ├── calib_pred_by_variant.rds
      ├── results_full.rds
      ├── S_seed_calib.rds
      ├── S_seed_fore_list.rds
      ├── spawn_dates_vec.rds
      ├── spawn_dates_by_env.rds
      ├── stoch_SAR_opts.rds
      └── american_river_instream.rda
  ```

---

## Script: precompute.R

**Purpose:**

1. Clear workspace and load libraries/helper functions (`functions.R`).
2. Read input data: environmental time series, carcass surveys, escapement estimates.
3. Define TDM variants (`exp_WF`, `exp_SM`, `lin_Martin`).
4. Build redd simulation records for observed years (2011–2024) and bootstrap future years.
5. Compile and run TDM per redd × alternative in parallel (`furrr`).
6. Calibrate the life-cycle model against observed escapement (2011–2024).
7. Specify stochastic SAR options and generate SAR time series.
8. Forecast full life-cycle to 2060 under stochastic SAR.
9. Save intermediate and final objects to `app_data/` for Shiny.

---

## Script: functions.R

Defines all modeling functions sourced by both `precompute.R` and Shiny:

* **hatch\_model(T):** days to hatch (958 ATU).
* **emergence\_model(T):** days to emergence (417 ATU).
* **tdm\_exp(temps, calib):** exponential TDM (WaterForum2020 or SALMOD2006; stage-split).
* **tdm\_lin\_martin(temps, α, β):** linear threshold TDM (Martin et al. 2017).
* **compute\_surv / compute\_surv\_by\_atu:** wrappers to extract temp slices and compute cumulative survival.
* **surv\_adult\_prespawn(deg\_day):** logistic adult pre-spawn survival.
* **compute\_deg\_day\_adult(...):** cumulative degree-day calculation per brood year.
* **simulate\_variant(...):** full age-structured life-cycle simulation (egg → adult).
* **generate\_SAR\_vec(...):** SAR time series generator with distribution/timing options.
* **sim\_forecast\_fn(...):** forecast wrapper (variant × env).

---

## Script: global.R

Preloads data and constants for the Shiny app:

1. Loads libraries and `functions.R`.
2. Reads all precomputed `.rds` and `.rda` data (`env_ext_list`, `egg_summary`, etc.).
3. Defines simulation constants (`n_sim`, `real_years`, `sim_years`).
4. Computes instream habitat lookup (`get_K_spawners(flow)`).
5. Exposes UI constants (`env_levels`, default `n_sim`).

---

## Script: app.R

Interactive Shiny application:

* **UI:** Tabs for About, Single Alternative, Compare Alternatives, Calibration, and Settings.
* **Inputs:** Simulation length, mode (Deterministic/Stochastic), TDM variant, alternative(s), flow, SAR distribution parameters.
* **Server:**

  * Runs calibration & forecasts.
  * `renderDT` tables with formatted outputs.
  * `renderPlot` for time series, survival curves, distributions, boxplots, SAR series.

Launch with:

```r
shiny::runApp("app.R")
```

---

## 📦 Key Saved Objects (in `app_data/`)

| File                             | Contents                                                                               | Used By                   |
| -------------------------------- | -------------------------------------------------------------------------------------- | ------------------------- |
| **egg\_summary.rds**             | Per env × variant × year mean egg→fry survival from TDM                                | Shiny plots & tables      |
| **base\_P\_list.rds**            | Calibrated life-cycle parameters (SAR\_mean, rear\_surv, K\_spawners) by variant × env | Forecast engine           |
| **results\_full.rds**            | Full life-cycle forecast outputs (spawners, survival, DD, etc.)                        | Shiny Compare/Single tabs |
| **calib\_results.rds**           | Optimized calibration fits (SAR\_mean, rear\_surv, SSE)                                | Calibration tab           |
| **calib\_pred\_by\_variant.rds** | Observed vs predicted spawner time series (2011–2024)                                  | Calibration tab           |
| **spawn\_dates\_by\_env.rds**    | Median spawn date series per env × year (LOCF/backfill)                                | Pre-spawn survival calc   |
| **S\_seed\_fore\_list.rds**      | Final 3 years of calibrated run used to seed forecasts                                 | Forecast engine           |
| **stoch\_SAR\_opts.rds**         | Distribution/timing settings for stochastic SAR draws                                  | Forecast engine           |
| **env\_ext\_list.rds**           | Daily temp inputs per env × site                                                       | TDM, degree-days          |
| **df\_all.rds**                  | Reference data table for sites & metadata                                              | Input joiners             |






