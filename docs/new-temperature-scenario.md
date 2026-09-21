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

The model reads the workbook **by column position**, so the layout matters more
than the labels. It must match the published deliverable,
`data_raw/SDM Power Bypass Temperature Modeling Results.xlsx`. The easiest route
is to copy that file and paste your numbers over it.

**One sheet per meteorological year**, named exactly `2011`, `2014`, `2017`,
`2020`. Other sheets (such as `Scenario Summary` or `Flow`) are copied along but
not used.

On each of those sheets:

| Row | Column A | Column B | Columns C–D | Columns E–F | … | Columns S–T |
|---|---|---|---|---|---|---|
| 1 | | | `No Bypass` | `Scenario 1` | … | `Scenario 6` |
| 2 | `Date` | `JDAY` | `AveWatt`, `AveHazel` | `AveWatt`, `AveHazel` | … | `AveWatt`, `AveHazel` |
| 3+ | a date | day of year | °C, °C | °C, °C | … | °C, °C |

- **Scenario order on row 1** (in C, E, G, … S): `No Bypass`, `Scenario 1`,
  `Scenario 2`, `Scenario 2b`, `Scenario 2c`, `Scenario 3`, `Scenario 4`,
  `Scenario 5`, `Scenario 6`. These become NB, PB1, PB2, PB2b, PB2c, PB3, PB4, PB5,
  PB6 in the app.
- **Sites on row 2:** Watt Avenue first, then Hazel Avenue, for every scenario.
- **Daily mean water temperature in degrees Celsius.**
- **Dates stored as Excel dates**, not text.
- Columns after T (the published file has a `Target Temp` column on the 2011
  sheet) are ignored.

`temperature_data.R` checks all of this before it does anything else and stops
with a plain message if something is off: a missing sheet, labels in the wrong
order, Watt and Hazel swapped, values that look like Fahrenheit, or dates stored
as text.

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
`env_ext_list.rds`, `df_all.rds` and `temperature_alternatives.xlsx` (the workbook
split into 36 sheets, one per scenario × met year) into your folder. It finishes
with two plots of the series; look at them.

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
- **`hydro_cost`** is the hydropower replacement cost for each alternative. It is a
  design input, not a model output. If your scenarios change the bypass volumes,
  supply new costs; otherwise the hydropower objective in Decision Support is
  wrong for your scenarios.
- **`default_weights`** are the objective weights the Decision Support sliders
  start at.

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
| Nine scenarios, in the order above, named NB … PB6 | `analysis/temperature_data.R`; `get_scenario_alternatives()` in `SalmonCountR/functions.R` **and** `SalmonCountR/app.R`; `precompute.R` (§37); `years.R` (`ARG_HYDRO_COST_*`, `arg_prepare_bundle()`); app dropdowns |
| Four met years 2011, 2014, 2017, 2020, weighted equally | the same places |
| Two sites, `AveWatt` and `AveHazel` | `temperature_data.R`; carcass section → site map in `precompute.R` (§9) |
| Calibration period 2011–2024, projection starting 2025 | `real_years` in `precompute.R` and `global.R`; GrandTab and carcass filters in `precompute.R`; `year > 2024` in `arg_prepare_bundle()` in `years.R` |
| Scenario temperatures used from the file's first date to Dec 31 | `threshold_start` / `threshold_end` in `temperature_data.R` |

Scenario numbering matters here: internally the 36 runs are numbered
met-year-major (1–9 = 2011 NB…PB6, 10–18 = 2014, and so on). A workbook with a
different number of scenarios would shift that numbering, which is why the check
in step 3 refuses one.

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
| `row 1: scenario labels … must be, in order …` | Scenario columns reordered or renamed. |
| `row 2: columns C to T must alternate AveWatt, AveHazel` | Site columns swapped or renamed. |
| `temperatures run 54.0 to 68.1 …` | Values are in Fahrenheit. |
| `column A has dates in year 1 …` | Dates stored as text. Reformat the column as dates. |
| `ARG_APP_DATA_DIR does not exist` | Create the folder first (step 2), or fix the path. |
| `readNWISdata` errors | No internet, or USGS is down. Try again later. |
| `unexpected symbol` part-way through precompute | The script was edited while it ran. Restart the run. |
| App shows *(data not loaded)* | Files missing from your folder; the banner names them. Did precompute finish? |
| Your entry is not in the selector | Restart the app after editing `years.R`. |
