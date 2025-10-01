# SalmonCountR — Lower American River Power Bypass Decision Support

End-to-end R workflow and Shiny app to simulate American River fall-run Chinook spawner dynamics under alternative power-bypass (water temperature) scenarios. It precomputes temperature-dependent mortality (TDM), calibrates a simple life cycle, and serves interactive comparisons across alternatives.

---

## Quickstart

```r
# From the repo root:
source("precompute.R")           # (run first time or whenever inputs change)
shiny::runApp("SalmonCountR")   # Launch the app
```

---

## Requirements

* **R** ≥ 4.0
* Packages:

  ```r
  tidyverse
  lubridate
  here
  data.table
  furrr
  future         # parallel backend for furrr
  compiler
  shiny
  DT
  # optional: ordinal (for CLM helpers in functions.R)
  ```

---

## Repository layout

```
SalmonCountR/
├── app.R               # Shiny UI & server
├── global.R            # Loads app_data and exposes helpers/constants
├── functions.R         # Core modeling library (TDM, survival, lifecycle, utilities)
├── precompute.R        # Builds inputs, runs calibration, writes app_data
└── app_data/
    ├── env_ext_list.rds
    ├── df_all.rds
    ├── egg_summary.rds
    ├── surv_lookup_full.rds
    ├── base_P_list.rds
    ├── calib_results.rds
    ├── calib_pred_by_variant.rds
    ├── S_seed_calib.rds
    ├── S_seed_fore_list.rds
    ├── stoch_SAR_opts.rds
    ├── sim_years.rds
    ├── spawn_dates_by_env.rds
    └── american_river_instream.rds   # WUA (m²) & flow → K_spawners lookup
```

> **Terminology:** In code and data, **`env`** means “power-bypass alternative.”

---

## What the scripts do

### `functions.R` (modeling library)

* **Egg/alevin development:** `hatch_model()`, `emergence_model()`
* **TDM models:** `tdm_exp()` (WF2020/SALMOD2006; egg/alevin stages), `tdm_lin_martin()`
* **ATU utilities:** `.stage_indices_by_atu()`, `.slice_by_atu()`, `compute_surv_by_atu()`
* **Adult holding:** `surv_adult_prespawn()` (logistic vs degree-days)
* **Season helpers:** `season_posix()`, `season_year()`
* **Spawn timing vector:** `build_spawn_vec_for_env()`
* **Degree-days:** `compute_deg_day_adult()`, `deg_day_cal_for()`
* **CLM helper:** `predict_clm_probs()` (category probabilities)
* **Forecast temps:** `build_forecast_temps()` (Oct/Nov, standardized)
* **Calibration objective:** `modular_sse()`
* **Life cycle:** `simulate_variant()` (age 3–5 returns, BH density dependence)
* **Forecast factory:** `sim_forecast_fn()` (encapsulates variant×env run)
* **SAR generator:** `generate_SAR_vec()` (Normal/Lognormal/Beta/Gamma; all/block/pulse)

### `precompute.R`

1. Loads inputs, builds Oct/Nov summaries, standardization, etc.
2. Constructs spawn date vectors and survival lookups: `surv_lookup_by_variant`, `surv_lookup_full`.
3. Calibrates (`SAR_mean`, `rear_surv`) by minimizing SSE via `modular_sse()`.
4. Saves: `calib_results.rds`, `calib_pred_by_variant.rds`, `base_P_list.rds`, `S_seed_calib.rds`, `S_seed_fore_list.rds`, `stoch_SAR_opts.rds`, `sim_years.rds`, `spawn_dates_by_env.rds`, `surv_lookup_full.rds`.

### `global.R`

* Reads all `app_data/*.rds`.
* Computes **`K_spawners`** from **`american_river_instream.rds`**: `K_spawners = FR_spawn_wua / 9.29` (m² per redd per SIT DSM).
* Exposes helper `get_K_spawners(flow_cfs)`.

### `app.R`

* Tabs: **About**, **Calibration**, **Single Alternative**, **Compare Alternatives**.
* Deterministic vs **stochastic SAR** (Normal/Lognormal/Beta/Gamma) with **all/block/pulse** timing.
* Outputs: tables (CSV export), time series, distributions, fry×DD, heatmaps, and boxplots.

---

## Data inputs (high level)

* **`env_ext_list.rds`**: per-env, per-site daily temps (Date, temp); used for ATU, Oct/Nov summaries, and degree-days.
* **`american_river_instream.rds`**: flow (cfs), **FR\_spawn\_wua** (m²). `global.R` computes `K_spawners = FR_spawn_wua / 9.29` and interpolates to arbitrary flows.

---

## Citation

Please cite this repo (tagged release) in reports/manuscripts.
