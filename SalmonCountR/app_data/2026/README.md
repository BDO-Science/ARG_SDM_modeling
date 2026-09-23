# app_data/2026

The app's **2026 (draft)** analysis year reads from this folder. It holds the
draft 2026 temperature deliverable run through the 2025 model.

## What is here

| File | From |
|---|---|
| `alt_key.rds`, `env_ext_list.rds`, `df_all.rds`, `observed_temps.rds`, `temperature_alternatives.xlsx` | `analysis/temperature_data.R` on `data_raw/TemperatureModelingResults_9-23-26.xlsx`, run on 2026-09-23 with `ARG_OBS_END = "2025-09-21"` |
| `results_full.rds` and the other model outputs | `SalmonCountR/precompute.R` with `ARG_APP_DATA_DIR` pointed here |

The deliverable (V. Martinez, 23 Sep 2026) is a draft. It has six scenarios,
and two more are expected: Scenarios 1 and 2 on the ATSP 38 schedule.

| Code | Workbook label |
|---|---|
| `NB` | ATSP 41 - No Bypass |
| `NB-38` | ATSP 38 - No Bypass |
| `PB1` … `PB4` | ATSP 41 - Scenario 1 … Scenario 4 |
| `PB1-38`, `PB2-38` | ATSP 38 - Scenario 1, Scenario 2 (expected) |

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

- **Hydropower costs** (`ARG_HYDRO_COST_2026` in `years.R`) carry each scenario's
  2025 value, on the assumption that the bypass schedules are unchanged; ATSP
  variants take their base scenario's cost, and No Bypass costs nothing under
  either schedule. Replace when Reclamation values the 2026 alternatives. Confirm
  with the temperature modeller that Scenarios 1–4 are the 2025 definitions.
- **Objective weights** are the 2025 elicited set.
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
