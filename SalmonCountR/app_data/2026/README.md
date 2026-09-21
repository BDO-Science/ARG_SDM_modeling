# app_data/2026

The app's **2026** analysis year reads from this folder. Until every file below
exists, the year selector shows **2026 (data not loaded)** and the banner names
what is missing. It will not error, and it will not silently fall back to 2025.

Required (see `ARG_YEAR_FILES` in `SalmonCountR/years.R`):

```
results_full.rds
steelhead_metrics.rds
swing_ranges.rds
american_river_instream.rds
df_all.rds
swing_scenario_results.rds
steelhead_scenario_results.rds
```

Optional: `data_vintage.rds`, written by `analysis/refresh_data_year.R`. If it
is present the banner shows the refresh date; if not it says the vintage was not
recorded.

## Filling it

Follow [`docs/new-temperature-scenario.md`](../../../docs/new-temperature-scenario.md)
with `ARG_APP_DATA_DIR = "SalmonCountR/app_data/2026"`.
`analysis/temperature_data.R` writes `df_all.rds` here, and `precompute.R` writes
the rest (including a copy of `american_river_instream.rds`). The 2025 results in
the parent folder are not touched.

## Before calling it a 2026 analysis

That procedure runs the 2026 temperatures through the **2025 model**: calibrated
on 2011–2024 escapement, with the projection starting in 2025. A 2026 analysis
that also adds the 2025 escapement and carcass data needs code changes first:

- `real_years` in `SalmonCountR/precompute.R` and `SalmonCountR/global.R`, and the
  literal `2011, 2024` / `2011:2024` filters on GrandTab and carcass data in
  `precompute.R`;
- `year > 2024` in `arg_prepare_bundle()` in `SalmonCountR/years.R`;
- `ARG_OBS_END` when running `temperature_data.R`: set it to the 2026 decision
  date, so observed temperatures run up to the decision and the scenarios take
  over after it (the run then freezes that record in `observed_temps.rds`);
- the carcass and GrandTab snapshots, via `analysis/refresh_data_year.R`.

Then, in the `"2026"` entry of `ARG_YEARS` in `years.R`:

- `first_projection_year` is set to 2026. It must match the first year
  precompute actually projects: **2025 unless `real_years` was extended.**
- `default_weights` and `hydro_cost` are copied from 2025 as placeholders. Replace
  them with the 2026 elicitation and valuation.

Changing the model inputs themselves is a Reclamation decision, not something
this scaffolding assumes.
