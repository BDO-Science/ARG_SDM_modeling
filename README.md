# SalmonCountR — Lower American River Power Bypass Decision Support

A life-cycle model of American River fall-run Chinook salmon, and a Shiny app for
comparing Folsom Dam power-bypass alternatives with it. Given daily water
temperatures for each alternative, the model projects spawner abundance 100 years
ahead through spawn timing, pre-spawn survival, temperature-dependent egg
mortality, density dependence and ocean returns. The app combines that with
steelhead and hydropower objectives for multi-criteria decision support.

This repository holds the model, the app, and the scripts behind every figure
and table in the accompanying manuscript. The published results are the 2025
analysis in `SalmonCountR/app_data/`.

| I want to… | Go to |
|---|---|
| Open the app | [Quickstart](#quickstart) |
| Run the model on my own temperature results | [`docs/new-temperature-scenario.md`](docs/new-temperature-scenario.md) |
| Understand what the model does | [`analysis/equations.qmd`](analysis/equations.qmd) |
| Reproduce a manuscript figure or number | [Manuscript figures and tables](#manuscript-figures-and-tables) |
| Deploy the app to shinyapps.io | [`SalmonCountR/DEPLOY.md`](SalmonCountR/DEPLOY.md) |

---

## Quickstart

Open `american_river_SDM.Rproj` in RStudio or Positron, so the working directory
is the repo root. Install the packages below, then:

```r
shiny::runApp("SalmonCountR")
```

The app reads precomputed results, so it starts in a few seconds and runs no
simulations itself. Everything it needs is committed in `SalmonCountR/app_data/`.

## Requirements

**R 4.1 or later** (the code uses the native `|>` pipe).

```r
install.packages(c(
  # the app (app.R, global.R, years.R, functions.R)
  "shiny", "shinyjs", "shinyWidgets", "DT", "tidyverse", "scales", "ggrepel", "data.table",
  # rebuilding the temperature series (analysis/temperature_data.R)
  "here", "readxl", "writexl", "dataRetrieval",
  # rerunning the model (SalmonCountR/precompute.R)
  "furrr", "future", "ordinal", "MASS", "ggridges",
  # manuscript figure and table scripts in analysis/
  "patchwork", "viridisLite", "janitor", "ggh4x",
  # deploying (SalmonCountR/deploy.R)
  "rsconnect"
))
```

The scripts that derive model parameters from raw data (see
[Other analysis scripts](#other-analysis-scripts)) use further packages, among
them `lme4`, `bbmle`, `coda`, `sf`, `cowplot`, `magick`, `rvest` and `pdftools`.
Install those as a script asks for them.

---

## How the pieces fit

```
data_raw/SDM Power Bypass Temperature Modeling Results.xlsx     CE-QUAL-W2 deliverable
  + USGS gauge data (downloaded)
        │  analysis/temperature_data.R
        ▼
SalmonCountR/app_data/env_ext_list.rds, df_all.rds             daily temperature, 36 alternatives
  + carcass survey, GrandTab escapement, spawning habitat       (fixed inputs, in app_data/)
        │  SalmonCountR/precompute.R            ~15 min
        ▼
SalmonCountR/app_data/results_full.rds and others               projections, objectives
        │  SalmonCountR/global.R + years.R
        ▼
SalmonCountR/app.R                                              the app
        │
        └──► analysis/*.R                                       manuscript figures and tables
```

**Terminology.** In code and data, **`env`** means "power-bypass alternative", not
"environment". There are 36: 9 alternatives (NB, PB1, PB2, PB2b, PB2c, PB3, PB4,
PB5, PB6) × 4 meteorological years (2011, 2014, 2017, 2020). `env` 1–9 are the 2011
met year, 10–18 are 2014, 19–27 are 2017, 28–36 are 2020; within each block the
order is NB, PB1, PB2, PB2b, PB2c, PB3, PB4, PB5, PB6.

**Weights.** The three egg-mortality (TDM) models are combined with elicited
weights: 0.51 Bratovich et al. 2020 ("Water Forum", `exp_WF`), 0.24 Bartholow &
Heasley 2006 ("SALMOD", `exp_SM`), 0.25 Martin et al. 2017 (`lin_Martin`). The four
met years are weighted equally. Default objective weights are 0.40 Chinook, 0.10
steelhead, 0.50 hydropower, the elicited set used in the manuscript.
`tdm_weight_sensitivity.R` and `evpi.R` also use a post-hoc set (0.73 / 0.05 / 0.22)
for the value-of-information analysis.

## Running the model on new temperatures

Step by step, including the workbook layout and what can and cannot change:
**[`docs/new-temperature-scenario.md`](docs/new-temperature-scenario.md)**.
In short:

```r
# 1. build the daily series from your workbook, into its own folder
Sys.setenv(ARG_TEMP_FILE    = "data_raw/my_temperatures.xlsx",
           ARG_APP_DATA_DIR = "SalmonCountR/app_data/my_scenario")
source("analysis/temperature_data.R")

# 2. in a fresh R session: run the model into the same folder
Sys.setenv(ARG_APP_DATA_DIR = "SalmonCountR/app_data/my_scenario")
source("SalmonCountR/precompute.R")

# 3. add a "my_scenario" entry to ARG_YEARS in SalmonCountR/years.R, then
shiny::runApp("SalmonCountR")   # and pick it under "Analysis year"
```

Without `ARG_APP_DATA_DIR`, both scripts write to the flat `app_data/` and
overwrite the published results.

---

## The app

A navbar app. The **Analysis year** selector at the top switches every tab
between result sets; see [Analysis years](#analysis-years).

| Tab | What it does |
|---|---|
| **About** | Background, the nine alternatives and their volumes and costs, model components, references. |
| **Temperature Explorer** | Daily temperature by alternative at Watt or Hazel Avenue, for October–December or the full first projection year, weighted across met years. Summary statistics table. |
| **Compare Alternatives** | Reweight the precomputed projections by TDM model and met year, and rescale for a downstream flow (500–5,000 cfs, default 1,000). Nothing runs until **Run Comparison** is clicked. Time series, boxplots of the last *N* years, and a summary table with CSV export. |
| **Swing Weighting** | Rank and score three hypothetical extreme alternatives to derive objective weights, then send them to Decision Support. |
| **Decision Support** | Equal or manual objective weights → overall weighted score, per-objective contributions, the consequence table, and trade-off plots. Uses the Compare Alternatives settings once a comparison has been run. |

How flow works in the app: the projections are run at a fixed 1,000 cfs spawning
capacity (K = 33,185 redds). The flow slider rescales those results for the change in
capacity after the fact; it does not rerun the model.

### Analysis years

An *analysis year* is one round of modelling: a temperature deliverable, the
precompute run built from it, and the design inputs (objective weights,
hydropower costs) for that round. It is not a year inside the projection. Years
are registered in `ARG_YEARS` in `SalmonCountR/years.R`, each pointing at a folder:

| Entry | Folder | State |
|---|---|---|
| `2025` (default) | `SalmonCountR/app_data/` | The published analysis |
| `2026` | `SalmonCountR/app_data/2026/` | Placeholder; shows *(data not loaded)* until filled. See its [README](SalmonCountR/app_data/2026/README.md) |

A folder counts as loaded once it holds the seven files in `ARG_YEAR_FILES`. The
app reads the list of years once at startup, so restart it after adding one. The
2025 results load at startup whatever year is selected; `functions.R` is not
year-aware, which is fine because the app does not call its simulation functions.

---

## Rerunning precompute

You do not need to rerun anything to use the app or reproduce the manuscript.
The scripts in `analysis/` read the committed `app_data/*.rds` and finish in
seconds (except the multi-seed scripts below).

To rerun, from the repo root:

```r
source("SalmonCountR/precompute.R")
```

- **Run it start to finish in a fresh R session.** The simulated redd set is drawn
  from `set.seed(123)` at the top; re-executing a chunk interactively shifts the
  RNG stream and changes the draw. The script saves `sim_redds` and `sim_future`,
  so the redd set behind a given set of results can be recovered without a re-run.
- **Do not edit `precompute.R` (or `functions.R`) while a run is in flight.** R
  reads the file incrementally, so an edit part-way through shifts byte offsets
  under the running process and it dies with a syntax error that is not in the
  file. A run takes about fifteen minutes.

### Settings

All optional. Set them with `Sys.setenv()` before sourcing; `precompute.R` clears
the workspace but not environment variables.

| Variable | Read by | Default | Effect |
|---|---|---|---|
| `ARG_APP_DATA_DIR` | `temperature_data.R`, `precompute.R` | `SalmonCountR/app_data` | Folder the temperature series is written to and read from, and precompute's outputs are written to. Point it at a subfolder for a new scenario. |
| `ARG_TEMP_FILE` | `temperature_data.R` | `data_raw/SDM Power Bypass Temperature Modeling Results.xlsx` | The temperature deliverable. |
| `ARG_OBS_END` | `temperature_data.R` | the day before the deliverable starts (2025-09-21) | The **decision date**: observed USGS temperatures are used through it, scenario temperatures after. |
| `ARG_OBS_REFRESH` | `temperature_data.R` | unset | `1` downloads the USGS record again instead of reusing the frozen copy (`observed_temps.rds`) in the scenario folder. |
| `ARG_SCENARIO_START` | `temperature_data.R` | the deliverable's first date | First day (`MM-DD`) the scenario temperatures are used. The published 2025 results used `10-18`; set that to reproduce them. |
| `ARG_SEED` | `precompute.R` | `123` | RNG seed, for replicating across seeds. |
| `ARG_SPAWN_TIMING` | `precompute.R` | `alternative` | `alternative`: each alternative is evaluated against the redds simulated under its own temperatures. `pooled`: the **superseded** behaviour, in which redds were pooled across all 36 alternatives before survival was computed. |

Pooling severed the temperature → spawn timing → thermal mortality channel that
the model's ordinal spawn-timing regression exists to represent, leaving
alternatives to differ only through incubation exposure. Changed to
alternative-specific on 2026-08-18; `pooled` is retained so earlier results can
be regenerated, and the repository state before the change is tagged
`pre-g1-correction`. What it moves, and by how much: **`docs/spawn-timing.md`**.

---

## Repository layout

```
SalmonCountR/                 the app, and the model run that feeds it
├── app.R                     Shiny UI and server
├── global.R                  loads the 2025 results at startup
├── years.R                   analysis-year registry (ARG_YEARS) and app path resolution
├── functions.R               model library: TDM, survival, spawn timing, life cycle
├── precompute.R              runs the model; writes app_data/  (not deployed)
├── deploy.R, DEPLOY.md       one-command, verified deploy to shinyapps.io
├── todo.md                   open modelling questions
└── app_data/                 inputs and results (below)
    └── 2026/                 placeholder for the next analysis year
analysis/                     standalone scripts: manuscript exhibits, parameter derivation,
│                             data refresh; model write-ups (equations.qmd, math.qmd)
└── additions/                exploratory one-offs
data_raw/                     raw inputs: temperature deliverables, swing weighting, CWT, HCI,
│                             american_river_data/, juvenile_data/ (screw traps)
docs/                         model change records and guides (see below)
figures/                      figures written by analysis/ (figures/bdsc/: conference talk)
output/                       tables and reports written by analysis/
presentations/                conference slides (Quarto)
archive/                      superseded code and outputs, kept for provenance; nothing here is live
```

### `SalmonCountR/app_data/`

| File | Written by | Read by |
|---|---|---|
| `env_ext_list.rds` | `analysis/temperature_data.R` | `precompute.R`, `global.R` — daily temperature per alternative (`Date`, `site`, `temp` °C, `alt`), 2011-09-01 to 2151-08-31 |
| `df_all.rds` | `analysis/temperature_data.R` | same data in one long table with `env`; app, steelhead metric |
| `observed_temps.rds` | `analysis/temperature_data.R` | the USGS record as downloaded for the decision date, reused on later runs so the analysis stays frozen (written on the next run; not deployed) |
| `carcassdet_*.csv` | SacPAS carcass survey download | `precompute.R` (spawn timing) |
| `grandtab_*.csv` | CDFW GrandTab download | `precompute.R` (calibration) |
| `american_river_instream.rds` | `analysis/spawn_habitat.R` | `precompute.R`, app — flow (cfs) → spawning WUA (m²); K = WUA ÷ 9.29 m² per redd |
| `results_full.rds` | `precompute.R` | app — the projection for every alternative × TDM model (114 labelled years; see below) |
| `swing_scenario_results.rds`, `steelhead_scenario_results.rds`, `swing_ranges.rds`, `steelhead_metrics.rds` | `precompute.R` | app — objective values and ranges |
| `egg_summary.rds`, `surv_lookup_full.rds`, `spawn_dates_by_alt.rds`, `base_P.rds`, `base_P_list.rds`, `calib_results.rds`, `S_seed_calib.rds`, `S_seed_fore_list.rds`, `sim_years.rds`, `stoch_SAR_opts.rds` | `precompute.R` | egg survival, calibrated parameters, seeds; read by `global.R` and `analysis/` |
| `sim_redds.rds`, `sim_future.rds` | `precompute.R` | the simulated redd set, kept for reproducibility |
| `calib_pred_by_variant.rds` | `analysis/calibration_fit_statistics.R` | calibration predictions |
| `spawn_timing_model.rds` | `analysis/build_spawn_timing_model.R` | the fitted spawn-timing model, standalone |
| `data_vintage.rds` | `analysis/refresh_data_year.R --apply` | app banner — see [Data provenance](#data-provenance) |
| `SAR LAR Releases.xlsx` | CWT release data | `analysis/sar_from_cwt.R` |

Also in the folder, and read by nothing: `S_seed.rds`, `american_river_instream.rda`,
`generate_SAR_vec.rds`, `simulate_variant.rds`, `rear_surv_lookup.rds`,
`spawn_dates.rds`, `spawn_dates_by_env.rds`, `spawn_dates_vec.rds`,
`swing_extreme_combos.rds`, `nonsalmon_objectives.csv`, `steelhead_objective.csv`.
Leftovers from earlier pipeline versions; `deploy.R` keeps them out of the bundle.

**Year labels in `results_full`.** One continuous 114-year run: the 2011–2024
calibration years followed by a 100-year projection seeded from observed 2022–2024
escapement. Its `year` column runs 2011–2124. Every reported metric uses the
projection years only (`year > 2024`) and takes the median of the final 20.

### `functions.R`

| Area | Functions |
|---|---|
| Egg and alevin development | `hatch_model()`, `emergence_model()`; thresholds `egg_ATU` = 400, `total_ATU` = 958 |
| Egg mortality (TDM) | `tdm_exp()` (Water Forum 2020 / SALMOD 2006, separate egg and alevin parameters), `tdm_lin_martin()` |
| Incubation window | `.stage_indices_by_atu()`, `.slice_by_atu()`, `compute_surv_by_atu()`, `compute_surv_vec()`, `memo_surv()` |
| Adult holding | `surv_adult_prespawn()`, `compute_deg_day_adult()`, `deg_day_cal_for()` |
| Spawn timing | `predict_clm_probs()`, `build_spawn_vec_for_env()`, `sample_dates_fast()`, `assign_period()` |
| Forecast temperatures | `build_forecast_temps()` |
| Survival by year | `eval_year()`, `pairs_for_env_year()` |
| Calibration | `combined_sse()` — one SAR and rearing survival fitted jointly across the three TDM models |
| Life cycle | `simulate_variant()` (ages 3–5, Beverton–Holt), `sim_forecast_fn()` |
| Stochastic SAR | `generate_SAR_vec()` (not used by the published run or the app) |
| Utilities | `get_scenario_alternatives()`, `season_posix()`, `season_year()`, `safe_range()`, `trim_trailing_text()` |

---

## Manuscript figures and tables

Standalone scripts in `analysis/`, each reading `SalmonCountR/app_data/` and
writing to `figures/` and `output/`. All run in seconds.

| Script | Produces |
|---|---|
| `reporting_values.R` | **The numbers reported in the manuscript**: adult index, composite score and rank, volume-normalised benefit, and the volume-vs-benefit rank correlation, all from the committed run |
| `figures.R` | Figure 4 — adult population index by alternative, faceted climate year × TDM model; also the baseline spawner forecast barchart |
| `mcda.R` | Figure 5 — composite MCDA scores, stacked by objective, with numeric bar labels |
| `figure3_tdm_curves.R` | Figure 3 — TDM daily survival and cumulative egg-to-fry survival, 10–18 °C |
| `tdm_weight_sensitivity.R` | TDM weight sensitivity of the composite score, plus the Martin-weight sweep |
| `elicitation_tables.R` | SI Tables S2-7 and S2-8 from the TDM elicitation scoresheet |
| `evpi.R` | Expected value of perfect information under both objective weight sets |
| `frontloading_cohort_decomposition.R` | Front-loading mechanism: crossover dates, hazard split, and the two-channel spawn-cohort decomposition |
| `calibration_fit_statistics.R` | Calibration predictions and fit statistics; writes `calib_pred_by_variant.rds` |
| `sar_from_cwt.R` | SAR from American River CWT release groups — provenance for every SAR figure in SI §S2.6 |

`elicitation_tables.R` reads a scoresheet that lives with the manuscript rather
than in this repo; point it there with the `ARG_SCORESHEET` environment variable.
Shared styling is in `figure_theme.R`; composite and normalisation helpers in
`composite_helpers.R`.

### Multi-seed scripts

These need `app_data` snapshots from several `precompute.R` runs rather than the
single committed copy, so they do **not** run in seconds. Point
`SPAWN_TIMING_SNAPROOT` at a directory of `<mode>_seed<n>/` folders, built by
running `precompute.R` with `ARG_SPAWN_TIMING` in `{alternative,pooled}` and
`ARG_SEED` in `{123,456,789,1011,1213}` (about 15 minutes per run; with
`ARG_APP_DATA_DIR` each run can write straight into its snapshot folder).

| Script | Produces |
|---|---|
| `spawn_timing_effect.R` | Every objective value the spawn-timing change moves, as means over seeds with the run-to-run range |
| `figure4_seed_uncertainty.R` | Figure 4 variant drawing both uncertainty tiers — year-to-year IQR and run-to-run spread across seeds |
| `compare_spawn_timing.R` | One pooled-vs-alternative pair, for a quick check |
| `compare_spawn_timing_seeds.R` | The full multi-seed replication, separating the effect from Monte Carlo noise |

**Why seeds matter here.** Each alternative is evaluated against its own redd
sample rather than a pooled sample ~36× larger, so per-alternative estimates carry
±173–217 fish of run-to-run noise. Differences under roughly 400 fish are not
resolvable from a single run; quote them as paired within-seed contrasts. See
`docs/spawn-timing.md`.

### Other analysis scripts

| Script | Purpose |
|---|---|
| `temperature_data.R` | Builds the daily temperature series from a deliverable plus USGS gauge data |
| `refresh_data_year.R` | Annual refresh of carcass, GrandTab and HCI data. Dry run by default; `--apply` to commit |
| `data_sources.R` | Register of every external input and how to fetch it |
| `build_spawn_timing_model.R` | Refits and saves the spawn-timing model on its own |
| `age_dist.R`, `cwt_data.R` | Age at return from coded-wire tag recoveries (RMIS) |
| `parameter_data.R` | Averages of published parameter tables |
| `spawn_habitat.R` | Builds `american_river_instream.rds`, the flow → spawning habitat (WUA) lookup |
| `flow_data.R` | Daily discharge and October–January monthly means |
| `2011_2024_spawners.R` | Modelled against observed escapement, 2011–2024 |
| `fall_run_spawn_dist_american.R` | Spatial and temporal spawning distribution from carcass data |
| `juvenile_rst_data.R`, `juvenile_rst_abundance.R` | Rotary screw trap juvenile abundance |
| `presentation_bdsc.R`, `presentation_plots.R` | Conference figures (`figures/bdsc/`) |
| `garage.R`, `additions/` | Scratch and exploratory work; not maintained |

---

## Data provenance

`analysis/refresh_data_year.R --apply` writes `app_data/data_vintage.rds`: when
the refresh ran, the repo commit at the time, and the size and modification date
of every data snapshot it found. The app shows it in the banner under the
**Analysis year** selector ("Data refreshed …"). A year folder without one shows
"Data vintage not recorded" and otherwise works normally. `deploy.R` requires it
for the published year.

## Known issues

| Issue | State |
|---|---|
| **The published 2025 results use scenario temperatures only from October 18.** A fixed day-291 start, left over from an earlier deliverable, replaced the 2025 deliverable's September 22–October 17 values with the long-term average, hiding the early cooling of PB1, PB2, PB2b and PB6. The adult index moves by less than the run-to-run noise when corrected; the steelhead metric (days below 18.3 °C in October–November) moves by up to 5.5 days. | **Fixed for future deliverables** (2026-09): the window now starts on the deliverable's first date. The 2025 analysis is kept as published; `ARG_SCENARIO_START = "10-18"` reproduces its temperature series exactly |
| **Forecast temperatures misalign in leap years.** The series is built by day of year, so leap years shift a day against non-leap years — up to 1.08 °C on a given calendar date. | Not fixed. Immaterial to the conclusions |
| **`sar_percent` in `app_data/SAR LAR Releases.xlsx` equals `sar`**, never multiplied by 100. | No live consumer is affected: `sar_from_cwt.R` recomputes it and `data_sources.R` records the defect. Fix when the workbook is next regenerated; do not edit the snapshot |
| **The "Add a Year" upload tab is not on `main`.** It was built on `revision-2026-08` (upload a deliverable in the app, get results back without R) and held back when app development moved to a contractor (commit `6b574ac`). `analysis/test_scenario_engine.R` tests that engine and fails on `main`. | Restore with `git revert 6b574ac` if wanted; design notes in `APP_DATA_UPDATE_OPTIONS.md` |

## Model documentation

| File | Contents |
|---|---|
| `analysis/equations.qmd` | The model, equation by equation, with parameter values |
| `analysis/math.qmd` | Derivation of the cumulative survival formulas |
| `docs/new-temperature-scenario.md` | Running the model on new temperatures |
| `docs/spawn-timing.md` | The alternative-specific spawn-timing correction: what changed, why, what it moves, and how to reproduce either behaviour |
| `APP_DATA_UPDATE_OPTIONS.md` | Options considered for annual updates, including the upload tab |

---

## Citation and license

Code: <https://github.com/BDO-Science/ARG_SDM_modeling>. Please cite a tagged
release (listed under the repository's tags; `v1.3.0` for the 2025 analysis as presented in 2026-09) together with the accompanying manuscript.

Licensed under the Apache License 2.0; see `LICENSE`.
