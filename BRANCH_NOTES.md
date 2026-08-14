# Branch `revision-2026-08` — what changed and where to pick up

This branch carries the revision-cycle work off `main` (17 commits, tip of
`main` at `26d61b8`). It is a handoff, not a finished PR: the analysis and the
repo reorganisation are done, the manuscript edits they support are not. Open
the PR when you are ready to take it on.

---

## Read the diff with rename detection on

GitHub will announce **282 files changed, 8,887 insertions, 5,933,646
deletions**. That deletion count is an artefact. With rename detection:

| Change | Files |
|----|----|
| Renamed / moved (176 of them byte-identical) | 187 |
| Added | 42 |
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
`analysis/equations.qmd` and `analysis/math.qmd`, so nothing is lost. They are
now covered by `.gitignore`.

**No file that anything depends on was deleted on this branch.**

---

## What changed

**Repository reorganisation.** `scripts_data_misc/` is gone as a catch-all,
replaced by `analysis/` (live scripts), `data_raw/` (raw inputs), and
`archive/` (superseded work, kept for provenance — nothing there is live).
`README.md` documents the resulting layout in full.

**Revision analysis** — the front-loading mechanism, TDM weight sensitivity,
EVPI, calibration fit statistics, SAR from CWT, and the elicitation tables.
Each has a script in `analysis/` that reads saved `.rds` files and runs in
seconds; none require re-running `precompute.R`. Outputs land in `figures/` and
`output/`.

**Figure restyling.** A shared `analysis/figure_theme.R`, viridis throughout,
larger bold type. All plot scripts were refactored onto it and their output
paths corrected.

**Scenario engine and the "Add a Year" tab.** `SalmonCountR/scenario_engine.R`
evaluates a temperature deliverable as the modelling team sends it — no
`precompute.R`, nothing saved server-side, under two seconds for four
meteorological years by nine scenarios. Surfaced as a new app tab, with
`analysis/refresh_data_year.R` for the annual data refresh and
`analysis/test_scenario_engine.R` as a regression test against the published
deliverable. Rationale and rejected options: `APP_DATA_UPDATE_OPTIONS.md`.

**Data provenance.** `refresh_data_year.R --apply` stamps `app_data/data_vintage.rds`;
the app shows it in every tab footer, on the About tab, and as `#` comment lines
atop every CSV export.

**Documentation.** Four new top-level documents, described below.

---

## State of the work

**Done and verified.** The analysis scripts all run and their outputs are
current. The mixed-run problem is closed — no value from a superseded run
survives in the manuscript or SI. The SAR provenance reproduces from CWT data.
The scenario engine passes its regression test.

**Not done.** The manuscript, SI and response have not received the figures and
tables this analysis produced. Four exhibits were promised to the reviewers and
none are in the documents yet; every one of them already exists in this repo at
a named path. That, plus seven internal inconsistencies a reviewer would catch,
is the open work.

Read these in this order:

| File | What it is |
|----|----|
| `OUTSTANDING_ITEMS.md` | **Start here.** What is still open, with a suggested order at the end |
| `MANUSCRIPT_REVISION_HANDOFF.md` | The findings organised by where the edit goes — for whoever edits the documents |
| `REVISION_FINDINGS.md` | The same work as an analysis record, with method detail |
| `APP_DATA_UPDATE_OPTIONS.md` | The annual app-data update question, options and recommendation |

Two things worth knowing before you rely on a number:

- **Four code items need a publication decision** (`OUTSTANDING_ITEMS.md` §G),
  because the manuscript points at this repository. G1 is the substantive one:
  alternative-specific spawn timing is computed and then overwritten, so the
  CLM's temperature-driven timing shift never reaches any result. It is flagged
  in the code and deliberately **not** changed — fixing it moves every
  downstream number.
- **`MANUSCRIPT_REVISION_HANDOFF.md` names commit `26d61b8` as its source of
  truth.** That is the tip of `main` and it is correct as a record of where the
  analysis was run, but this branch is 17 commits past it. For anything about
  the current state of the code, this branch is authoritative.

---

## Picking it up

Nothing here is load-bearing on the branch name or on this file — both can go
when the PR merges. `OUTSTANDING_ITEMS.md` §"Suggested order" is the actual
work queue, and G4 (one applied `refresh_data_year.R --apply` run to populate
the provenance stamp) is the smallest useful first commit.
