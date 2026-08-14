# Branch `revision-2026-08` — app and modeling changes

This branch carries 18 commits off `main` (tip of `main` at `26d61b8`). It
reorganises the repository, restructures the analysis scripts, and adds a
scenario engine and a new Shiny tab. This file covers the app and modeling
side; the manuscript-side documents in the repo root are a separate track and
are not needed to work on the code.

---

## Read the diff with rename detection on

GitHub will announce **283 files changed, ~9,000 insertions, 5,933,646
deletions**. That deletion count is an artefact. With rename detection:

| Change | Files |
|----|----|
| Renamed / moved (176 of them byte-identical) | 187 |
| Added | 43 |
| Deleted | 40 |
| Modified | 13 |

Nearly all 5.9M "deleted" lines are the BT-SPAS-X MCMC output — `CODAchain1-3.txt`
at ~2M lines each, plus `TSPDE--results.txt` and the `inits*.txt` files — which
**moved** from `scripts_data_misc/calibration/` to `archive/calibration/`. They
are still tracked and unchanged. GitHub's default view gives up on rename
detection at that file size, so it renders them as delete + add. `git diff -M40%
origin/main...HEAD` shows the truth.

The 40 real deletions are all rendered Quarto build artefacts —
`math_files/`, `scripts_data_misc/equations_files/`, `equations.html`,
`math.html` — plus a stray `.RDataTmp`. The sources are retained as
`analysis/equations.qmd` and `analysis/math.qmd`, and `.gitignore` now covers
the render output.

**No file that any script or the app depends on was deleted on this branch.**

---

## What changed

**Repository layout.** `scripts_data_misc/` is gone as a catch-all, replaced
by `analysis/` (live scripts), `data_raw/` (raw inputs), and `archive/`
(superseded work, kept for provenance — nothing there is live, including the
old life-cycle and ETF survival implementations). `README.md` documents the
resulting layout and every `app_data/*.rds` object in full.

**Analysis scripts.** Each derived quantity now has a standalone script in
`analysis/` that reads the saved `app_data/*.rds` and runs in seconds — TDM
weight sensitivity, EVPI, calibration fit statistics, SAR from CWT release
groups, the front-loading cohort decomposition, and the elicitation tables.
None of them require re-running `precompute.R`. Outputs go to `figures/` and
`output/`.

**Figure styling.** A shared `analysis/figure_theme.R` — viridis throughout,
larger bold type. Every plot script was refactored onto it and had its output
path corrected.

**Scenario engine and the "Add a Year" tab.** `SalmonCountR/scenario_engine.R`
evaluates a temperature deliverable in the form the modelling team sends it and
returns consequence-table and MCDA results. Pure functions, no global state, no
`precompute.R`, nothing saved server-side; a full run over four meteorological
years and nine scenarios takes under two seconds. Surfaced as a new app tab.

| File | Role |
|----|----|
| `SalmonCountR/scenario_engine.R` | The engine |
| `analysis/build_spawn_timing_model.R` | Refits and saves the spawn-timing CLM. **Required** — the tab is disabled without `app_data/spawn_timing_model.rds` |
| `analysis/test_scenario_engine.R` | Regression test against the published deliverable |
| `analysis/refresh_data_year.R` | Annual data refresh; safe by default, `--apply` to commit |
| `analysis/data_sources.R` | Register of every external input and how to fetch it |

Rationale and the options that were rejected: `APP_DATA_UPDATE_OPTIONS.md`.

> **Comparing runs:** the engine weights every spawn date by its probability
> where `precompute.R` draws a finite sample of redds, so its composite scores
> are **not** interchangeable with the published ones. Compare an upload against
> the baseline re-run through the same engine — which is what the tab's "Against
> the baseline" view does.

**Data provenance.** `refresh_data_year.R --apply` writes
`app_data/data_vintage.rds` — when the refresh ran, the repo commit at the time,
and the size and modification date of every snapshot found. The app surfaces it
in every tab footer, in a Data Provenance section on About, and as `#` comment
lines atop every CSV export (both `read.csv(comment.char = "#")` and
`readr::read_csv(comment = "#")` skip them). The file is absent until the first
`--apply` run; the app reports "not recorded" and runs normally.

---

## Known code issues — decide before this goes public

The manuscript points at this repository, so these are publication decisions
rather than housekeeping. Full text in `OUTSTANDING_ITEMS.md` §G.

| # | Issue | State |
|----|----|----|
| G1 | **Alternative-specific spawn timing is computed and then discarded.** `precompute.R:659` splits redds by alternative and year; line 773 overwrites that with a split by year only, pooling all 36 alternatives. The CLM's temperature-driven shift in spawn timing therefore never reaches any result. | Flagged in the code, deliberately **not** changed — fixing it moves every downstream number. Decide whether to disclose, fix, or defend it as holding spawn timing constant so differences are attributable to incubation exposure alone |
| G2 | **Forecast temperatures misalign in leap years.** The series is built by day-of-year, so leap years shift a day against non-leap years — up to 1.08 °C on a given calendar date. | Not fixed. Immaterial to the conclusions, but a real defect anyone reading the code will find |
| G3 | **`sar_percent` in `app_data/SAR LAR Releases.xlsx` is identical to `sar`** — never multiplied by 100, so reading it as a percentage gives values 100× too small. | Not fixed |
| G4 | **`app_data/*.rds` provenance is unstamped.** Everything is built and wired; it needs one applied refresh to populate. | Run `analysis/refresh_data_year.R --apply` |

---

## Picking it up

**G4 is the smallest useful first commit** — one applied refresh run, and the
provenance stamp stops reading "not recorded" everywhere in the app. G1 is the
one that needs a real decision before anyone quotes a number from this repo.

Before running anything: `precompute.R` must be run start to finish in a clean
session. The simulated redd set comes from `set.seed(123)` at the top, so
re-executing a chunk interactively shifts the RNG stream and changes the draw.
The script saves `sim_redds` and `sim_future` so the redd set behind a given
result can be recovered without a full re-run. Reproducing anything in
`figures/` or `output/` does not need `precompute.R` at all.

Neither this file nor the branch name is load-bearing — both can go when the
branch merges.
