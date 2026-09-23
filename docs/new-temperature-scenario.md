# Running the model on your own temperatures

This walks through taking a new set of CE-QUAL-W2 temperature results (a new
deliverable, a revised set of bypass schedules, a sensitivity run), running the
salmon model on them, and looking at the results in the Shiny app next to the
published 2025 analysis.

You do not need to understand the model to do this. You do need R, about half an
hour, and an internet connection. For what the model does, see
[`analysis/equations.qmd`](../analysis/equations.qmd).

```
your workbook (.xlsx)
      │  analysis/temperature_data.R        a few minutes; downloads USGS gauge data
      ▼
env_ext_list.rds, df_all.rds   ──┐
      │  SalmonCountR/precompute.R  ~15 min
      ▼                          │  all in one scenario folder:
results_full.rds and friends   ──┘  SalmonCountR/app_data/<name>/
      │  one entry in SalmonCountR/years.R
      ▼
the app's "Analysis year" selector
```

Nothing here touches the published results, as long as you do step 2.

---

## 0. Before you start

- **R 4.1 or later**, run from the repository root (open `american_river_SDM.Rproj`
  in RStudio or Positron and you are there).
- Packages:

  ```r
  install.packages(c(
    "tidyverse", "here", "readxl", "writexl", "dataRetrieval",   # temperature_data.R
    "furrr", "future", "data.table", "ordinal", "MASS", "ggrepel", "ggridges",  # precompute.R
    "shiny", "shinyjs", "shinyWidgets", "DT", "scales"           # the app
  ))
  ```

## 1. Put your workbook in the expected layout

The reader **discovers the scenarios from the labels in row 1**, so a workbook
can carry any number of scenarios, in any order, and adding one is a data
change. Both deliverables so far read this way: the published 2025 file,
`data_raw/SDM Power Bypass Temperature Modeling Results.xlsx`, and the 2026
draft, `data_raw/TemperatureModelingResults_9-23-26.xlsx`.

**One sheet per meteorological year**, named exactly `2011`, `2014`, `2017`,
`2020`. Other sheets (`Scenario Summary`, `Flow`, `Averaged All Years`) are
ignored.

On each of those sheets:

| Row | Column A | Column B | Scenario block | Scenario block | … |
|---|---|---|---|---|---|
| 1 | | | `No Bypass` | `ATSP 38 - Scenario 1` | … |
| 2 | `Date` | `JDAY` | site names, e.g. `AveFol`, `AveWatt`, `AveHazel` | `AveWatt`, `AveHazel` | … |
| 3+ | a date | day of year | °C per site | °C per site | … |

- **A scenario block starts at every non-empty cell in row 1** from column C
  onward and runs to the column before the next label. The order of blocks is
  the order the alternatives appear in the app.
- **Row 2 names the sites.** Each block must contain an `AveWatt` and an
  `AveHazel` column; any other site column in the block (the 2026 file has
  `AveFol`) is ignored, as is anything after the last block (the 2025 file has
  `Target Temp` on the 2011 sheet). The columns can be in any order.
- **Every met-year sheet must list the same scenarios in the same order.**
- **Daily mean water temperature in degrees Celsius.** A blank day in a scenario
  falls back to the climatology for that day.
- **Dates stored as Excel dates**, not text.

**How labels become codes.** `No Bypass` becomes `NB` and `Scenario 2b` becomes
`PB2b`. A label may carry an ATSP schedule, as in `ATSP 41 - Scenario 1`. When
the workbook mixes schedules, the base schedule keeps the plain codes so they
line up with earlier years, and the others get a suffix: `ATSP 38 - No Bypass`
becomes `NB-38`. The base is the most common schedule in the file; set
`ARG_ATSP_BASE` (for example `"41"`) to choose it explicitly. The mapping is
printed when the script runs and saved as `alt_key.rds` (see step 3).

`temperature_data.R` checks all of this before it does anything else and stops
with a plain message if something is off: a missing sheet, a block without both
site columns, sheets that disagree about the scenarios, values that look like
Fahrenheit, dates stored as text, or two labels that map to the same code.

### What the model actually uses from your file

This is the part most likely to surprise you.

- **From your file's first date to December 31.** The window starts on the
  earliest date in the workbook, whatever it is, so a deliverable that begins
  earlier or later is used from its own start. Every day with a value is used.
  Days outside the file (in the published deliverable, which runs September 22
  to November 30: everything before September 22 and all of December) use the
  2011–2025 average for that day at that site, **the same for every scenario**.
  If your bypass schedules differ in December, include December in the file.
  (The published 2025 analysis used a fixed October 18 start instead; see the
  README's known issues. To reproduce it, set `ARG_SCENARIO_START = "10-18"`.)
- **The calendar year in your file is ignored.** Dates are matched by day of year,
  and that seasonal pattern repeats in every year of the 100-year projection. A
  file dated 2026 behaves exactly like one dated 2025.
- **Before the projection starts, temperatures are observed gauge data,** not your
  file. Those come from USGS (Watt Avenue 11446980, Hazel Avenue 11446500) and are
  identical for every scenario.

## 2. Make a scenario folder

Pick a short name with no spaces and create the folder:

```
SalmonCountR/app_data/<name>/        e.g. SalmonCountR/app_data/2025_revised
```

Everything for your run goes here. **Do not skip this step.** The default folder
is the flat `SalmonCountR/app_data/`, which holds the published results, and both
scripts overwrite what they find there.

## 3. Build the daily temperature series

In a fresh R session at the repo root:

```r
Sys.setenv(
  ARG_TEMP_FILE    = "data_raw/my_temperatures.xlsx",       # your workbook
  ARG_APP_DATA_DIR = "SalmonCountR/app_data/2025_revised"   # your folder
)
source("analysis/temperature_data.R")
```

This checks the workbook, downloads the observed gauge record, and writes
`env_ext_list.rds`, `df_all.rds`, `alt_key.rds` and `temperature_alternatives.xlsx`
(the workbook split into one sheet per run, scenario × met year) into your
folder. It finishes with two plots of the series; look at them.

**`alt_key.rds` is the map from run number to alternative.** Runs are numbered
met-year-major (every scenario for 2011, then every scenario for 2014, and so
on), and the key records `env`, `alt` (the code), `label` (the workbook's own
label), `met_year` and `atsp` for each. `precompute.R` and the app read it, so
they never assume how many alternatives there are.

**The decision date.** Observed USGS temperatures are used through the decision
date and the scenario temperatures after it. Before the decision every
alternative is the same river, so observed data is the right input up to that
day. Set it with `ARG_OBS_END`, for example `ARG_OBS_END = "2026-10-20"`. If you
don't, it defaults to the day before your workbook starts. If the gauge record
stops short of the date you give, the script says so and uses observed data
only as far as it goes.

**Freezing it.** The first run saves the downloaded record as
`observed_temps.rds` in your folder, and every later run reuses it, so the
inputs stay as they were at the decision. USGS revises provisional data for
months afterwards, and a fresh download would otherwise change them quietly.
To download again on purpose, set `ARG_OBS_REFRESH = "1"`.

**If the decision date is in a new water year,** extend the calibration years
as well; see [`SalmonCountR/app_data/2026/README.md`](../SalmonCountR/app_data/2026/README.md).
`precompute.R` stops with an explanation if you forget, because the first
projection year would otherwise be observed data for every alternative.

## 4. Run the model

Start **a new R session** (Session → Restart R), then:

```r
Sys.setenv(ARG_APP_DATA_DIR = "SalmonCountR/app_data/2025_revised")
source("SalmonCountR/precompute.R")
```

It prints `Scenario data folder: …` near the top. Check that it is your folder.
The run takes about fifteen minutes. **Do not edit `precompute.R` or
`functions.R` while it runs.** R reads the script as it goes, and an edit
mid-run kills it with a misleading syntax error.

When it finishes it prints `Outputs written to …`, and your folder holds
`results_full.rds` and the other files the app reads, plus a copy of the habitat
lookup (`american_river_instream.rds`).

What stays fixed: the model is recalibrated on each run, but against the same
observed escapement (2011–2024), carcass data and habitat curve. Calibration is
fitted to the observed temperature record, which your workbook does not touch,
so the calibrated parameters should not move. The random redd draw uses the same seed (`ARG_SEED`, default 123)
as the published run. Differences between your run and the published one
therefore come from your temperatures, up to the ±173–217 adults of run-to-run
noise described in [`spawn-timing.md`](spawn-timing.md).

## 5. Register it with the app

Add an entry to `ARG_YEARS` in `SalmonCountR/years.R`. Copy the `"2025"` entry
and change it:

```r
  "2025_revised" = list(
    label                 = "2025 (revised bypass schedules)",  # what the selector shows
    dir                   = file.path("app_data", "2025_revised"),
    default_weights       = c(chinook = 0.40, steelhead = 0.10, hydro = 0.50),
    hydro_cost            = ARG_HYDRO_COST_2025,
    first_projection_year = 2025,
    note                  = "Revised schedules from <source>, <date>."
  ),
```

- **`first_projection_year` must be `2025`.** The projection starts in the year
  after the calibration period, which is fixed at 2011–2024. It does not follow
  the date in your workbook.
- **`temperature_year`** is the calendar year your workbook models (its dates).
  The Temperature Explorer shows and labels that year. Leave it out to use the
  first projection year.
- **`hydro_cost`** is the hydropower replacement cost for each alternative, named
  by the codes in your `alt_key.rds`. It is a design input, not a model output.
  **Every alternative in the key needs a cost**, or the year shows as *(data not
  loaded)* and the banner names the alternative; an alternative with no cost
  would otherwise look free. If your scenarios change the bypass volumes, supply
  new costs; otherwise the hydropower objective in Decision Support is wrong for
  your scenarios.
- **`default_weights`** are the objective weights the Decision Support sliders
  start at.
- **`objective_ranges`** (optional) switches the year to global scaling: fixed
  0–1 ranges for each objective, as in the `"2026"` entry. Leave it out to scale
  locally over the alternatives in the run, as the published 2025 analysis does.
  See OBJECTIVE SCALING in `years.R`.

## 6. Look at it

```r
shiny::runApp("SalmonCountR")
```

Choose your entry in **Analysis year** at the top of any tab. If it shows
*(data not loaded)*, the banner lists the files it could not find in your folder.
The app reads its list of years when it starts, so restart it after editing
`years.R`.

To put it on shinyapps.io, see [`SalmonCountR/DEPLOY.md`](../SalmonCountR/DEPLOY.md).
Only the files the app reads are uploaded from a scenario folder.

---

## What you cannot change without editing code

The pipeline assumes the shape of the 2025 decision. Each of these is hard-coded
in several places, so changing one is a code change, not a data change:

| Assumption | Where |
|---|---|
| Four met years 2011, 2014, 2017, 2020, weighted equally | `HYDRO_YEARS` in `temperature_data.R`; the met-year weight sliders in `app.R`; the About tab |
| Two sites, `AveWatt` and `AveHazel` | `temperature_data.R`; carcass section → site map in `precompute.R` (§9) |
| Calibration period 2011–2024, projection starting 2025 | `real_years` in `precompute.R` and `global.R`; GrandTab and carcass filters in `precompute.R`; `year > 2024` in `arg_prepare_bundle()` in `years.R` |
| Scenario temperatures used from the file's first date to Dec 31 | `threshold_start` / `threshold_end` in `temperature_data.R` |
| The About tab's table of alternative specifications (volumes, MWh, cost) | `app.R`; it describes the 2025 alternatives |

The number and names of the scenarios are **not** on this list any more: they
come from the workbook, through `alt_key.rds`. The one thing a new scenario
needs from you is a hydropower cost in `years.R` (step 5).

The manuscript scripts in `analysis/` are the exception. They document the
published 2025 analysis, read the flat `SalmonCountR/app_data/`, and carry that
run's nine alternatives, costs and bypass volumes as literals. They are not
meant to run on a scenario folder.

## Known limitations

- **Leap years.** Matching by day of year shifts the seasonal pattern by one day
  in leap years of the projection.
- **Capacity is fixed at 1,000 cfs.** Spawning capacity (K = 33,185 redds) is held
  at the 1,000 cfs reference flow for every scenario. The app's flow slider
  rescales the results afterwards; it does not rerun the model.

## If something goes wrong

| Symptom | Cause |
|---|---|
| `Sheet(s) 2017 not found …` | Sheet names must be exactly the four met years. |
| `row 1 has no scenario labels from column C onward` | Scenario names are not in row 1, or the sheet has extra header rows. |
| `scenario '…' has no AveWatt or AveHazel column` | A site column is missing or misspelt under that label. |
| `Sheet 2014 lists different scenarios … from sheet 2011` | Met-year sheets disagree; make the row-1 labels identical on every sheet. |
| `Scenario labels map to duplicate codes` | Two labels reduce to the same code (e.g. two `Scenario 1` on the same ATSP). |
| `hydro_cost in years.R for …` in the app banner | An alternative in `alt_key.rds` has no cost in `years.R` (step 5). |
| `temperatures run 54.0 to 68.1 …` | Values are in Fahrenheit. |
| `column A has dates in year 1 …` | Dates stored as text. Reformat the column as dates. |
| `ARG_APP_DATA_DIR does not exist` | Create the folder first (step 2), or fix the path. |
| `readNWISdata` errors | No internet, or USGS is down. Try again later. |
| `unexpected symbol` part-way through precompute | The script was edited while it ran. Restart the run. |
| App shows *(data not loaded)* | Files missing from your folder; the banner names them. Did precompute finish? |
| Your entry is not in the selector | Restart the app after editing `years.R`. |
