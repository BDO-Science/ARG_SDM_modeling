# archive/

Nothing in this folder is live. No script outside it reads anything here, and
nothing here is deployed. It is kept so that earlier results can be traced.

| Folder | What it holds |
|---|---|
| `app_data_snapshot/` | `SalmonCountR/app_data/` as it stood before the 2026 revision (pre-`g1` spawn timing, old ATU thresholds, fixed K) |
| `calibration/` | BTSPAS run for juvenile abundance: JAGS model, inits and the three raw CODA chains (~107 MB). Produced by `analysis/juvenile_rst_abundance.R` |
| `legacy_inputs/` | Raw inputs no current script reads: the October 2024 temperature deliverable and its placeholders, an earlier HCI download, NOAA stoplight indicators, and old American River lookups |
| `legacy_outputs/` | Outputs with no current writer: age-composition tables from an earlier `age_dist.R`, superseded forecast summaries, orphan figures, and `app_data_leftovers/` (files earlier pipeline versions wrote into `app_data/`) |
| `old_scripts/` | Superseded model code (`old_life_cycle/`, `etf_survival_old/`, `misc_old/`) and unmaintained scratch (`scratch/`: `garage.R`, `early_spawners.R`, the BTSPAS `model.txt`) |
| `review/` | The reviewer's edited copy of the repository (`ARG_SDM_Ayers_edits-master/`), as received |
| `stale_calibration_figures/` | Calibration figures from before the 2026 recalibration |
| `code_review.txt` | A 2025 code review of the app and model. `code_review.md` is the same text with markdown escaping |
