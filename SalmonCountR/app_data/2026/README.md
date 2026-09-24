# app_data/2026

The app's **2026 (draft)** analysis year reads from this folder. It holds the
draft 2026 temperature deliverable run through the 2025 model.

## What is here

| File | From |
|---|---|
| `alt_key.rds`, `env_ext_list.rds`, `df_all.rds`, `observed_temps.rds`, `temperature_alternatives.xlsx` | `analysis/temperature_data.R` on `data_raw/TemperatureModelingResults_9-24-26.xlsx`, run on 2026-09-24 with `ARG_OBS_END = "2025-09-21"` |
| `results_full.rds` and the other model outputs | `SalmonCountR/precompute.R` with `ARG_APP_DATA_DIR` pointed here |

The deliverable (V. Martinez, 24 Sep 2026, with the ARG ad hoc deck
`DRAFT_WaterTemperatureModeling_9_24_26_ARGAdhoc.pptx`) is a draft with eight
scenarios. It replaces the 23 Sep file, which had six and is kept in
`archive/legacy_inputs/`. The deck's schedule chart is in
`SalmonCountR/www/2026_bypass_schedules_draft.png` and shown on the About tab.

| Code | Workbook label |
|---|---|
| `NB` | ATSP 41 - No Bypass |
| `NB-38` | ATSP 38 - No Bypass |
| `PB1` … `PB4` | ATSP 41 - Scenario 1 … Scenario 4 |
| `PB1-38`, `PB2-38` | ATSP 38 - Scenario 1, Scenario 2 |

## When the updated deliverable arrives

Drop the new workbook in `data_raw/`, then from the repo root, in a fresh R
session:

```r
Sys.setenv(ARG_TEMP_FILE    = "data_raw/<new file>.xlsx",
           ARG_APP_DATA_DIR = "SalmonCountR/app_data/2026",
           ARG_OBS_END      = "2025-09-21")     # see below
source("analysis/temperature_data.R")
```

Check the printed code mapping (`NB-38`, `PB1-38`, `PB2-38` for the ATSP 38
scenarios). Then, in another fresh session:

```r
Sys.setenv(ARG_APP_DATA_DIR = "SalmonCountR/app_data/2026")
source("SalmonCountR/precompute.R")             # about 15 minutes
```

`ARG_HYDRO_COST_2026` in `SalmonCountR/years.R` already lists all eight codes.
If the new file uses other labels, add their codes there or the year shows as
*(data not loaded)* and the banner says which alternative has no cost. Restart
the app and pick **2026 (draft)** under Analysis year. Commit the folder.

## What is placeholder

- **Hydropower costs and specifications** (`ARG_ALT_SPECS_2026` in `years.R`)
  are Reclamation's draft 2026 valuation of Scenarios 1–4
  (`docs/2026_hydropower_costs_draft.png`, 23 Sep 2026), on the 90% exceedance (per K. Thielen)
  WY2026 operations forecast; the 50% case is recorded in the comments there.
  S1–S4 in that table are `PB1`–`PB4` here (the cooling at Hazel Avenue starts
  on each schedule's first bypass day). ATSP variants take their base
  scenario's values, and No Bypass costs nothing under either schedule.
  **The 2026 schedules are not the 2025 ones with the same code**: 2026 PB1 is
  the 2025 PB2 schedule, 2026 PB2 is the 2025 PB4 schedule, and PB3 and PB4 are
  new (500 cfs from Oct 28 and from Oct 15, to Nov 30). S5–S8 are TBD.
- **Objective weights** are the 2025 elicited set. They were elicited against
  2025's local swings, so they do not yet describe the global ranges below;
  re-elicit them in the Swing Weighting tab, which now presents those ranges.
- **Objective ranges** (`ARG_OBJECTIVE_RANGES_2026` in `years.R`) are the
  global-scaling ranges B. Mahardja proposed in September 2026 as a starting
  point: Chinook 0–25,000, steelhead 0–61 days, hydropower $0–3 million. The
  2025 tab is unaffected and keeps its local scaling.
- **Bypass volumes** are not in the deliverable (no `Scenario Summary` sheet), so
  nothing here reports volume per alternative.

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
