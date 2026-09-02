# app_data/2026

Drop the 2026 precompute output here. Until every file below exists, the app's
year selector shows **2026 (data not loaded)** and the banner names what is
missing — it will not error, and it will not silently fall back to 2025.

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

Optional: `data_vintage.rds` — written by `analysis/refresh_data_year.R`. If it
is present the banner shows the refresh date; if not it says the vintage was not
recorded.

## Also update, in `SalmonCountR/years.R`

The `"2026"` entry currently carries placeholders:

- `default_weights` — copied from 2025. Replace with the 2026 elicitation.
- `hydro_cost` — copied from 2025. Replace with the 2026 valuation.
- `first_projection_year` — set to 2026. This drives the forecast filter and the
  Temperature Explorer labels; check it matches what precompute actually produced.

## What generates these files

`SalmonCountR/precompute.R`, run against the 2026 CE-QUAL-W2 temperature
deliverable. It currently writes to the flat `app_data/` directory — point its
output here, or run it and move the results, so the 2025 bundle stays intact.

Changing the model inputs themselves is a Reclamation decision, not something
this scaffolding assumes.
