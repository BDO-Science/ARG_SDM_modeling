# app_data/2026

The app's **2026 (draft)** analysis year reads from this folder. It holds the
2026 temperature deliverable run through the 2025 model.

## What is here

| File | From |
|---|---|
| `alt_key.rds`, `env_ext_list.rds`, `df_all.rds`, `observed_temps.rds`, `temperature_alternatives.xlsx` | `analysis/temperature_data.R` on `data_raw/TemperatureModelingResults_9-30-26.xlsx`, run on 2026-09-30 with `ARG_OBS_END = "2025-09-21"` and `ARG_ATSP_BASE = "41"` |
| `results_full.rds` and the other model outputs | `SalmonCountR/precompute.R` with `ARG_APP_DATA_DIR` pointed here |

The deliverable (30 Sep 2026) has ten scenarios. It replaces the 24 Sep file,
which had eight and is kept in `archive/legacy_inputs/` with the 23 Sep one
(six). It adds `PB5` and `PB6` and revises the temperatures of the other eight
(by up to 1.1 °C on single days at Watt Avenue for `PB3`, under 0.5 °C for the
rest). The schedules are in `data_raw/2026BypassModelingScenarios.xlsx`:
descriptions for Scenarios 1–4 on `Sheet1`, daily bypass flow for Scenarios 1–6
on `Timeseries`. The 24 Sep ARG ad hoc deck's schedule chart, which draws
`PB1`–`PB4` only, is in `SalmonCountR/www/2026_bypass_schedules_draft.png` and
shown on the About tab.

| Code | Workbook label | Scenario |
|---|---|---|
| `NB` | `NBP` | No Bypass, ATSP 41 |
| `NB-38` | `NBP - 38` | No Bypass, ATSP 38 |
| `PB1` … `PB6` | `PB1` … `PB6` | Scenario 1 … Scenario 6, ATSP 41 |
| `PB1-38`, `PB2-38` | `PB1-38`, `PB2-38` | Scenario 1, Scenario 2, ATSP 38 |

This workbook labels its scenarios by code, where the earlier ones wrote them
out (`ATSP 38 - Scenario 1`); the reader takes either. A code-style label does
not name the base ATSP schedule, so `ARG_ATSP_BASE = "41"` records it in the
key.

## When an updated deliverable arrives

Drop the new workbook in `data_raw/`, then from the repo root, in a fresh R
session:

```r
Sys.setenv(ARG_TEMP_FILE    = "data_raw/<new file>.xlsx",
           ARG_APP_DATA_DIR = "SalmonCountR/app_data/2026",
           ARG_OBS_END      = "2025-09-21",     # see below
           ARG_ATSP_BASE    = "41")
source("analysis/temperature_data.R")
```

Check the printed code mapping (`NB-38`, `PB1-38`, `PB2-38` for the ATSP 38
scenarios). Then, in another fresh session:

```r
Sys.setenv(ARG_APP_DATA_DIR = "SalmonCountR/app_data/2026")
source("SalmonCountR/precompute.R")             # about 15 minutes
```

`ARG_HYDRO_COST_2026` in `SalmonCountR/years.R` lists all ten codes. If the new
file uses other labels, add their codes there or the year shows as
*(data not loaded)* and the banner says which alternative has no cost. Restart
the app and pick **2026 (draft)** under Analysis year. Commit the folder.

## What is placeholder

- **Hydropower costs and specifications** (`ARG_ALT_SPECS_2026` in `years.R`)
  are Reclamation's 2026 valuation of `PB1`–`PB6`
  (`docs/2026_hydropower_costs_draft.png`, received 30 Sep 2026), with
  operations data from the `september_update.xlsx` forecast. It replaces the
  23 Sep draft valuation of `PB1`–`PB4` (90% exceedance forecast, v1.4.1).
  ATSP variants take their base scenario's values, as the table itself lists
  them, and No Bypass costs nothing under either schedule.
  **The 2026 schedules are not the 2025 ones with the same code**: 2026 PB1 is
  the 2025 PB2 schedule, 2026 PB2 the 2025 PB4 schedule and 2026 PB6 the 2025
  PB2b schedule; PB3 and PB4 (500 cfs from Oct 28 and from Oct 15, to Nov 30)
  and PB5 (250 cfs from Oct 21, 500 cfs from Oct 28, 250 cfs from Nov 11, to
  Nov 30) are new.
- **Objective weights** are the 2025 elicited set. They were elicited against
  2025's local swings, so they do not yet describe the global ranges below;
  re-elicit them in the Swing Weighting tab, which now presents those ranges.
- **Objective ranges** (`ARG_OBJECTIVE_RANGES_2026` in `years.R`) are the
  global-scaling ranges B. Mahardja proposed in September 2026 as a starting
  point: Chinook 0–25,000, steelhead 0–61 days, hydropower $0–3 million. The
  2025 tab is unaffected and keeps its local scaling.
- **Bypass volumes** are not in the temperature workbook (no `Scenario Summary`
  sheet); the About tab takes them from the valuation table.

## Why the decision date is 2025-09-21

The run uses the **2025 model**: calibrated on 2011–2024 escapement, with the
projection starting in 2025. The first projection year has to carry the
scenario temperatures, so the observed record is cut at 2025-09-21 and the 2026
scenario pattern is used from the 2025 fall onward. The temperatures are
matched by day of year, so the calendar year on the workbook does not matter.
`first_projection_year` for the `"2026"` entry in `years.R` is therefore 2025.

A true 2026 start, adding the 2025 escapement and carcass data, needs code
changes first:

- `real_years` in `SalmonCountR/precompute.R` and `SalmonCountR/global.R`, and the
  literal `2011, 2024` / `2011:2024` filters on GrandTab and carcass data in
  `precompute.R`;
- `year > 2024` in `arg_prepare_bundle()` in `SalmonCountR/years.R`;
- `ARG_OBS_END` set to the 2026 decision date when running `temperature_data.R`;
- the carcass and GrandTab snapshots, via `analysis/refresh_data_year.R`;
- `first_projection_year = 2026` in `years.R`.

Changing the model inputs themselves is a Reclamation decision, not something
this scaffolding assumes.
