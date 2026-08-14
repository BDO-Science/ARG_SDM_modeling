# SalmonCountR — Lower American River Power Bypass Decision Support

End-to-end R workflow and Shiny app to simulate American River fall-run Chinook spawner dynamics under alternative power-bypass (water temperature) scenarios. It precomputes temperature-dependent mortality (TDM), calibrates a simple life cycle, and serves interactive comparisons across alternatives.

---

## Quickstart

```r
# From the repo root:
source("SalmonCountR/precompute.R")   # rebuild app_data (first run, or when inputs change)
shiny::runApp("SalmonCountR")         # launch the app
```

**Run `precompute.R` start to finish in a clean session.** The simulated redd set
is drawn from `set.seed(123)` at the top of the script; re-executing a chunk
interactively shifts the RNG stream and changes the draw. The script saves
`sim_redds` and `sim_future` so the redd set behind a given set of results can be
recovered without a full re-run.

Reproducing the manuscript figures and tables does **not** require re-running
`precompute.R` — the scripts in `analysis/` read the saved `app_data/*.rds` and
finish in seconds.

---

## Requirements

* **R** ≥ 4.0

Install everything in one go:

```r
install.packages(c(
  # app  (app.R + global.R)
  "shiny", "shinyjs", "shinyWidgets", "DT", "tidyverse", "scales", "ggrepel", "here",
  # precompute.R additionally
  "furrr", "future", "data.table", "MASS", "ordinal", "ggridges", "readxl",
  # manuscript figure/table scripts in analysis/
  "patchwork", "viridis", "janitor"
))
```

`compiler` ships with R. **`shinyWidgets` is required for the app to start and
`furrr` for `precompute.R` to run** — check both before a season kickoff.

---

## Repository layout

```
analysis/               # standalone analysis & figure scripts (see below)
│   └── additions/      # exploratory one-offs
data_raw/               # raw inputs consumed by analysis/ scripts
│   ├── american_river_data/
│   ├── juvenile_data/          # rotary screw trap catch, efficiency, operations
│   └── *.xlsx, *.csv           # temperature modelling, swing weighting, CWT, HCI
figures/                # publication figures (written by analysis/)
output/                 # tables and reports (written by analysis/)
archive/                # superseded, kept for provenance — nothing here is live
│   ├── calibration/    # BT-SPAS-X MCMC output (CODA chains, TSPDE results, .stan)
│   ├── old_scripts/    # earlier life-cycle and ETF survival implementations
│   ├── review/         # external review copies of the app and precompute
│   └── legacy_outputs/ # outputs from earlier versions of analysis/age_dist.R
SalmonCountR/
├── app.R               # Shiny UI & server
├── global.R            # Loads app_data and exposes helpers/constants
├── functions.R         # Core modeling library (TDM, survival, lifecycle, utilities)
├── precompute.R        # Builds inputs, runs calibration, writes app_data
├── scenario_engine.R   # Evaluates an uploaded deliverable (the "Add a Year" tab)
└── app_data/           # written by precompute.R, read by global.R
    ├── env_ext_list.rds            # per-alternative daily temps (Date, site, temp, alt)
    ├── df_all.rds                  # same, long form with `env`
    ├── egg_summary.rds             # cumulative egg survival by env × variant × year
    ├── surv_lookup_full.rds
    ├── base_P.rds, base_P_list.rds # life-cycle parameters (calibrated)
    ├── calib_results.rds           # calibrated SAR_mean and rear_surv
    ├── S_seed_calib.rds            # 2011-2013 observed escapement (calibration seed)
    ├── S_seed_fore_list.rds        # 2022-2024 observed escapement (forecast seed)
    ├── stoch_SAR_opts.rds
    ├── sim_years.rds
    ├── spawn_dates_by_alt.rds
    ├── results_full.rds            # the 114-year projection, all 36 env × 3 variants
    ├── steelhead_metrics.rds       # days < 18.3 °C in Oct/Nov, by env
    ├── swing_scenario_results.rds  # Chinook metric per alternative (default weights)
    ├── steelhead_scenario_results.rds
    ├── swing_ranges.rds            # worst/best case per objective, for swing weighting
    ├── american_river_instream.rds # WUA (m²) & flow → K_spawners lookup
    ├── spawn_timing_model.rds      # CLM behind "Add a Year" (analysis/build_spawn_timing_model.R)
    └── data_vintage.rds            # provenance stamp (analysis/refresh_data_year.R --apply)
```

> **Terminology:** In code and data, **`env`** means “power-bypass alternative.”
> There are 36 of them: 9 alternatives (NB, PB1, PB2, PB2b, PB2c, PB3–PB6) ×
> 4 meteorological years (2011, 2014, 2017, 2020). `env` 1–9 are the 2011 met
> year, 10–18 are 2014, 19–27 are 2017, 28–36 are 2020; within each block the
> order is NB, PB1, PB2, PB2b, PB2c, PB3, PB4, PB5, PB6.

> **Year labels in `results_full`:** this object is one continuous 114-year
> projection seeded from the observed 2022–2024 escapement, but its `year`
> column is labelled 2011–2124. Treat it as projection year, not calendar year.

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
* Loads `data_vintage.rds` and exposes `data_vintage_label()`, `data_vintage_sources()`
  and `data_vintage_lines()` — see **Data provenance** below.

### `app.R`

* Tabs: **About**, **Temperature Explorer**, **Compare Alternatives**, **Swing Weighting**, **Decision Support**, **Add a Year**.
* Deterministic vs **stochastic SAR** (Normal/Lognormal/Beta/Gamma) with **all/block/pulse** timing.
* Outputs: tables (CSV export), time series, distributions, fry×DD, heatmaps, and boxplots.
* Default objective weights are 0.40 Chinook, 0.10 steelhead, 0.50 hydropower — the
  elicited set used in the manuscript.

---

## Data provenance

`analysis/refresh_data_year.R --apply` writes `app_data/data_vintage.rds`: when the
refresh ran, the repo commit at the time, and the size and modification date of every
snapshot it found. The app surfaces it in three places, so a number on screen or in a
downloaded file can be traced back to the data behind it:

* a one-line stamp in the footer of **every** tab;
* a **Data Provenance** section on the About tab, listing each snapshot;
* `#` comment lines at the head of every CSV export, which `read.csv(..., comment.char = "#")`
  and `readr::read_csv(..., comment = "#")` both skip.

The file is absent until the first `--apply` run. That is not an error — the app reports
"not recorded" and runs normally, but nothing then ties `app_data/` to a dated snapshot
(`OUTSTANDING_ITEMS.md` G4).

---

## Data inputs (high level)

* **`env_ext_list.rds`**: per-env, per-site daily temps (Date, temp); used for ATU, Oct/Nov summaries, and degree-days.
* **`american_river_instream.rds`**: flow (cfs), **FR\_spawn\_wua** (m²). `global.R` computes `K_spawners = FR_spawn_wua / 9.29` and interpolates to arbitrary flows.

---

## Manuscript figures and tables

Standalone scripts in `analysis/`, each reading `SalmonCountR/app_data/`
and writing to `figures/` and `output/`:

| Script | Produces |
|---|---|
| `figures.R` | Figure 4 — adult population index by alternative, faceted climate year × TDM model; also the baseline spawner forecast barchart |
| `mcda.R` | Figure 5 — composite MCDA scores, stacked by objective, with numeric bar labels |
| `figure3_tdm_curves.R` | Figure 3 — TDM daily survival and cumulative egg-to-fry survival, 10–18 °C |
| `tdm_weight_sensitivity.R` | TDM weight sensitivity of the composite score, plus the Martin-weight sweep |
| `elicitation_tables.R` | SI Tables S2-7 and S2-8 from the TDM elicitation scoresheet |
| `evpi.R` | Expected value of perfect information under both objective weight sets |
| `frontloading_cohort_decomposition.R` | Front-loading mechanism: crossover dates, hazard split, spawn-cohort decomposition |
| `calibration_fit_statistics.R` | Calibration predictions and fit statistics; repopulates `calib_pred_by_variant.rds` |
| `sar_from_cwt.R` | SAR from American River CWT release groups — the provenance for every SAR figure quoted in SI §S2.6, with the text checked against the data |

`elicitation_tables.R` reads a scoresheet that lives with the manuscript rather
than in this repo; point it there with the `ARG_SCORESHEET` environment variable.

All of these read precomputed `.rds` files, so they run in seconds without
re-running `precompute.R`.

The remaining scripts in `analysis/` derive model parameters from raw data in
`data_raw/` (age structure, CWT returns, flow, spawning habitat, juvenile
abundance, temperature scenarios). They document provenance and are not part of
the app pipeline.

### Adding a new year of temperature modelling

The app has an **Add a Year** tab: upload the temperature deliverable exactly as
the modelling team sends it, and get consequence-table and MCDA results back. No
R, no `precompute.R`, nothing saved server-side. A full run over four
meteorological years and nine scenarios takes under two seconds.

| Script | Purpose |
|---|---|
| `SalmonCountR/scenario_engine.R` | The engine. Pure functions, no global state. |
| `analysis/build_spawn_timing_model.R` | Refits and saves the spawn-timing CLM. **Required** — the tab is disabled without it. |
| `analysis/test_scenario_engine.R` | Regression test against the published deliverable. |
| `analysis/refresh_data_year.R` | Annual data refresh. Safe by default; `--apply` to commit. |
| `analysis/data_sources.R` | Register of every external input and how to fetch it. |

Full rationale and the options that were considered:
`APP_DATA_UPDATE_OPTIONS.md`.

> **Comparing runs:** the engine weights every spawn date by its probability
> where `precompute.R` draws a finite sample of redds, so its composite scores
> are not interchangeable with the published ones. Compare an upload against the
> baseline **re-run through the same engine** — which is what the app's "Against
> the baseline" tab does.

### Revision reports

| File | Purpose |
|---|---|
| `MANUSCRIPT_REVISION_HANDOFF.md` | Findings organised as manuscript/SI/response edits — the one to hand to someone making revisions |
| `REVISION_FINDINGS.md` | The same work as an analysis record, with method detail |

---

## Citation

Please cite this repo (tagged release) in reports/manuscripts.
