# Annual data refresh -- 29 July 2026

**Mode: dry run.** Nothing will be written. Re-run with `--apply` to commit.

## Sources

| Source | Automatic | Snapshot |
|---|---|---|
| SacPAS Chinook Carcass Survey Detail | no -- manual | `SalmonCountR/app_data/carcassdet_1752789274_15.csv` |
| SacPAS CDFW GrandTab Adult Escapement | no -- manual | `SalmonCountR/app_data/grandtab_1752793045_337.csv` |
| SacPAS Water Year Hydrologic Classification Indices | no -- manual | `data_raw/hci_1753379372_198.csv` |
| USGS NWIS instantaneous water temperature | no -- manual | `embedded in SalmonCountR/app_data/env_ext_list.rds` |
| CE-QUAL-W2 power bypass scenario temperatures | no -- manual | `data_raw/SDM Power Bypass Temperature Modeling Results.xlsx` |
| PHABSIM weighted usable area | no -- manual | `SalmonCountR/app_data/american_river_instream.rds` |
| Coded wire tag releases and returns | no -- manual | `SalmonCountR/app_data/SAR LAR Releases.xlsx` |
| Hydropower revenue loss by alternative | no -- manual | `hard-coded in analysis/mcda.R and SalmonCountR/app.R` |

> No SacPAS query links are configured yet, so nothing can be downloaded
> automatically. See the instructions at the top of `analysis/data_sources.R`:
> tick **Generate Query Result Link Only** on each SacPAS query page, then
> paste the link it produces into `DATA_SOURCES`. Until then this script
> reports on the bundled snapshots only.

## Current snapshots

| Source | File | Modified | Size |
|---|---|---|---|
| SacPAS Chinook Carcass Survey Detail | `carcassdet_1752789274_15.csv` | 2026-07-28 | 957 KB |
| SacPAS CDFW GrandTab Adult Escapement | `grandtab_1752793045_337.csv` | 2026-07-28 | 2 KB |
| SacPAS Water Year Hydrologic Classification Indices | `hci_1753379372_198.csv` | 2026-07-28 | 13 KB |
| USGS NWIS instantaneous water temperature | env_ext_list.rds | -- | not on disk |
| CE-QUAL-W2 power bypass scenario temperatures | `SDM Power Bypass Temperature Modeling Results.xlsx` | 2026-07-28 | 106 KB |
| PHABSIM weighted usable area | `american_river_instream.rds` | 2026-07-28 | 0 KB |
| Coded wire tag releases and returns | `SAR LAR Releases.xlsx` | 2026-07-28 | 22 KB |
| Hydropower revenue loss by alternative | app.R | -- | not on disk |

## What changed

- **Carcass survey detail**: no new download, nothing to compare.
- **GrandTab escapement**: no new download, nothing to compare.
- **Hydrologic index**: no new download, nothing to compare.

## Calibration window

- Escapement available through **2024**.
- The model calibrates on **2011-2024**.
- No new escapement years since the last calibration.

## Spawn-timing model

- Fitted 2026-07-28 on 10087 carcass observations, brood years 2011-2024.
- Coefficients: Oct_std +0.1860, Nov_std -0.0698.
- Refit with `Rscript analysis/build_spawn_timing_model.R` after new carcass data lands.

## Actions

- Nothing to apply.

### Still to do by hand

1. Re-run `SalmonCountR/precompute.R` start to finish in a clean session.
   The redd draw is seeded but only reproducible on a full run.
2. Re-run `analysis/build_spawn_timing_model.R`.
3. Re-run `analysis/calibration_fit_statistics.R` and **look at the fit** before
   accepting the new calibration. As of 2026-07-28 it was poor (R2 0.13,
   Nash-Sutcliffe -0.59); tracked in the project revision notes.
4. Re-run `analysis/test_scenario_engine.R` and confirm it still passes.
5. Regenerate figures if any published number moved.
6. Drop the new CE-QUAL-W2 deliverable in `data_raw/` so the app's baseline
   comparison uses the current year.

