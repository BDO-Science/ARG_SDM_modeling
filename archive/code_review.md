Awesome—this script is doing a lot. Here’s a structured, detailed walkthrough of what each section is doing, why, and where to watch for edge cases.



---



\# 1) Setup \& data ingest



\*\*What it does\*\*



\* Wipes the R session with `rm(list=ls())`.

\* Loads core packages: tidyverse, lubridate, furrr (parallel map), data.table (fast tabular ops), compiler (to JIT-compile hot functions), here (paths), MASS/ordinal (for CLMM), ggridges (plots).

\* Sources your modeling functions from `SalmonCountR/functions.R`.



\*\*Inputs read\*\*



\* `env\_ext\_list.rds`: daily temps/flows (per \*environment\*; looks like a site × date × temp long table per env).

\* `df\_all.rds`: broader data (not used later in this file).

\* `carcassdet\_\*.csv`: carcass survey data (for timing \& site info).

\* `grandtab\_\*.csv`: observed escapement (2011–2024).



\*\*Parallel / threading\*\*



\* Detects cores, picks `n\_workers` up to 6, sets `plan(multisession, workers = n\_workers)` (future) and `setDTthreads(n\_workers)` (data.table). Good balance between R-level and C-level threading.



\*\*Why\*\*



\* Everything that follows is compute-heavy; this keeps things fast and reproducible across machines.



---



\# 2) TDM model menu \& site/date indices



\*\*What it does\*\*



\* Declares the three TDM “variants” you want to evaluate:



&nbsp; \* `exp\_WF` (WaterForum 2020)

&nbsp; \* `exp\_SM` (SALMOD 2006)

&nbsp; \* `lin\_Martin` (Martin 2017).

\* Builds per-env “lookups”:



&nbsp; \* `site\_temps\_list`: `site -> numeric temp vector` (fast daily extraction).

&nbsp; \* `site\_dates\_list`: `site -> Date vector` (aligned to temps).

&nbsp; \* `date\_idx\_list`: `site -> named integer index` (names = ISO date strings), so you can turn a date into a fast vector position and slice temps directly.



\*\*Why\*\*



\* `compute\_surv()` later needs to turn `(site, spawn date)` into a daily temperature \*slice\* quickly.



---



\# 3) Carcass data → usable spawn-timing covariates



\*\*What it does\*\*



\* Creates `obs\_df` with: `spawn\_dt` (7 days before carcass date), `brood\_year` (Sept-based), `site` (recode from sections), `fork\_length`.

\* Joins \*\*October/November mean temps\*\* by (site, brood\\\_year); then standardizes to z-scores (`Oct\_std`, `Nov\_std`).



\*\*Why\*\*



\* The CLMM you fit uses standardized Oct/Nov temps and a year RE to model \*which 10-day bin fish spawn in\*.



---



\# 4) Define 10-day “spawn bins” and assign them to observations



\*\*What it does\*\*



\* Defines 11 fixed bins from Oct 10 → Jan 27, with \*\*dummy years\*\* 2000–2001 to keep Oct–Jan contiguous.

\* Maps each `spawn\_dt` into `spawn\_bin` via a 2000/2001 dummy date comparison.



\*\*Why\*\*



\* The CLMM is ordinal over these 11 bins; dummy year simplifies bin boundary logic.



---



\# 5) CLMM fit (spawn timing model)



\*\*What it does\*\*



\* Filters to complete rows and fits:



&nbsp; ```

&nbsp; clmm(spawn\_bin ~ Oct\_std + Nov\_std + (1 | brood\_year), link = "logit")

&nbsp; ```

\* Extracts:



&nbsp; \* fixed slopes `beta = {Oct\_std, Nov\_std}`

&nbsp; \* cumulative thresholds `zeta` (K−1 cutpoints)

&nbsp; \* per-year random effects (BLUPs).



\*\*Why\*\*



\* Gives a temperature-driven (Oct/Nov) estimate of \*spawn-bin probabilities\*, with year-to-year departures captured as a random effect.



---



\# 6) Year random effects → AR(1) for forecasting



\*\*What it does\*\*



\* Centers the historical BLUPs `u\_hat`, estimates AR(1) `phi` from lag-1 ACF, and innovation sd `sigma`.

\* Safeguards if estimates are degenerate; clamps `phi` to ±0.95.

\* `simulate\_u()` creates a forward series of year effects `u\_fc` for future years.



\*\*Why\*\*



\* Keeps “personality” of years: warm/late or cool/early years persist somewhat, but regress to mean.



---



\# 7) Sanity plots of CLMM predictions (optional diagnostics)



\*\*What it does\*\*



\* Builds a small `predict\_clmm()` that turns `(beta, zeta)` and new covariates into per-bin probabilities via logit CDF differences.

\* Produces two diagnostic plots:



&nbsp; \* Probabilities vs \*\*Oct\\\_std\*\* (holding Nov at mean).

&nbsp; \* Probabilities vs \*\*Nov\\\_std\*\* (holding Oct at mean).



\*\*Why\*\*



\* Quick visual QA that slopes go in the expected direction and thresholds look sensible.



---



\# 8) Forecast-year spawn-date simulation (per environment)



\*\*Preparation\*\*



\* `forecast\_temps`: for each env × future year, computes \*\*mean Oct \& Nov temps\*\*, then z-scores \*\*using the historical means/sds\*\* you computed above. (Great: keeps scale consistent with the CLMM fit.)



\*\*Bin probabilities\*\*



\* Determine `present\_bins` (those observed while fitting) and `all\_bins` (the full 11).

\* For each env × year:



&nbsp; \* Compute linear predictor `xβ + u\_y` (using `u\_fc` for future).

&nbsp; \* Turn thresholds into cumulative probs with `plogis(ζ\_k − xβ − u\_y)`.

&nbsp; \* Difference to category probs; set missing bins to 0; renormalize.



\*\*Sampling redds \& dates\*\*



\* For each env × year:



&nbsp; \* Determine `N\_redd`: rounded mean redds/year from carcass observations (2011–2024). (Acts as a \*count\* to sample spawn dates; doesn’t constrain final escapement.)

&nbsp; \* Draw `N\_redd` bins from that year’s probability vector.

&nbsp; \* Convert bins to dates via `sample\_dates\_fast()`:



&nbsp;   \* Vectorized: builds start/end Date per redd from precomputed bin month/day and the brood year offset (2000↔2001 logic carried across).

&nbsp;   \* Picks a uniform day within the 10-day window.



\* Samples `section`, `fork\_length`, `site` from historical pools for that env.



&nbsp; \* \*\*Watch out:\*\* if `site\_pool` is missing, you fall back to \*\*`env\_sites`\*\*, but that object isn’t defined in this script. You probably want `site\_env\_map %>% filter(env == env\_i) %>% pull(site)` (which you already computed) instead of `env\_sites`.



\*\*Outputs\*\*



\* `sim\_future`: a big table of simulated redds per (env, year), with dates and some attributes for realism.

\* Several QA plots (ridges, heatmaps) over \*\*dummy years\*\* to compare distributions year-to-year.



\*\*Why\*\*



\* This builds realistic future \*\*spawn-date distributions\*\* coherent with the CLMM + AR(1) and env-specific fall temps.



---



\# 9) Combine historical \& future redds → median spawn date per year



\*\*What it does\*\*



\* `sim\_actual`: observed redds in 2011–2024.

\* `sim\_redds <- bind\_rows(sim\_actual, sim\_future)`.

\* `spawn\_dates\_real`: \*\*median spawn date by year\*\* across \*all envs combined\*.

\* Builds `sim\_years\_vec = 2011:(2011+n\_sim−1)`.

\* Builds `spawn\_dates\_vec`: forces the per-year median \*\*month/day\*\* into the desired `sim\_years\_vec` (i.e., you copy month/day to target year).



\*\*Why\*\*



\* You use these per-year medians as \*\*end dates\*\* for pre-spawn degree-day accumulation (Oct 1 → spawn date) in `compute\_deg\_day\_adult()` during calibration \& forecasts.



**\*\*Note / gotcha\*\***



**\* `spawn\_dates\_vec` is \*\*global per year\*\* (not per env). If you want \*\*env-specific\*\* degree-day end dates, compute medians by `(env, sim\_year)` and pass env-specific vectors into `compute\_deg\_day\_adult()`. As written, \*\*all envs share the same spawn date per year\*\*.**



---



\# 10) Compile hot functions \& build caches



\*\*What it does\*\*



\* JIT-compiles `hatch\_model`, `emergence\_model`, `tdm\_exp`, `tdm\_lin\_martin`.

\* Converts `sim\_redds` to data.table and splits by `sim\_year` for fast access.

\* Builds `env\_cache` with, for each env:



&nbsp; \* `date\_idx\_env` (site → date index)

&nbsp; \* `temps\_env` (site → temp vector)

&nbsp; \* min/max dates available (to guard bounds)



\*\*Why\*\*



\* These are in the inner loop of survival computations; compilation + caching cuts runtime a lot.



---



\# 11) TDM survival per env × variant × year



\*\*What it does\*\*



\* Loops over \*\*all simulation years\*\* (`sim\_years` = calibration + forecast horizon).



\* For each year:



&nbsp; \* Pulls redds (site, spawn date) from `sim\_redds\_split`.

&nbsp; \* Converts each spawn date to a \*\*cohort calendar date\*\* (`yy` logic moves Jan spawns into next calendar), builds a vector `rdr`.

&nbsp; \* Drops redds outside the env’s data date range.

&nbsp; \* Collapses duplicate `(site, rdr)` pairs into unique rows with a \*\*weight N\*\* (how many redds share that pair) → big speedup.



\* For each TDM variant:



&nbsp; \* Calls `compute\_surv()` for every unique `(site, rdr)`:



&nbsp;   \* Turns `(site, date)` into a temp slice via `date\_idx\_env` + `temps\_env`.

&nbsp;   \* \*\*Exp\*\*: `tdm\_exp()` uses stage-specific hazards with ATU-based egg/alevin split; survival `S = exp(-sum(h))`.

&nbsp;   \* \*\*Linear\*\*: `tdm\_lin\_martin()` uses `h = α·max(T−β,0)` then cumulates.

&nbsp; \* Returns a \*\*weighted mean\*\* survival for the year using weights = number of redds sharing the pair.



\* Produces `results\_obs\_fast` with `mean\_cum\_surv` per `env, variant, sim\_year`.



\*\*Why\*\*



\* This is the core egg→fry survival time series you’ll feed into the life-cycle model.



**\*\*Note / nuance\*\***



**\* `compute\_surv()` selects the window length from \*\*spawn-day temperature T0\*\* using `hatch\_model(T0) + emergence\_model(T0)`. That’s a \*constant-T\* approximation to determine the number of days to slice from the temps. Inside `tdm\_exp()`, the \*\*egg/alevin switch\*\* is done by \*\*accumulated ATU over the slice\*\*. If you want the slice itself to \*end\* exactly when total ATU is reached under the \*actual\* temps, you could replace the length heuristic with a “grow the window until ATU ≥ 1375” logic. The current approach is fast and usually close; it can misalign when temps drift far from T0.**



---



\# 12) Summaries \& survival lookups



\*\*What it does\*\*



\* `egg\_summary`: mean of `mean\_cum\_surv` by env/variant/year (already one per year).

\* `surv\_lookup\_by\_variant`: for \*\*calibration env only\*\* (`ref\_env` = first env name), pulls the vector (2011–2024) per variant.

\* `surv\_lookup\_full`: for \*\*all\*\* env×variant, stores the full-year vectors keyed by `"env\_variant"`.



\*\*Why\*\*



\* These lookups feed the life-cycle model in both \*\*calibration\*\* (ref env) and \*\*forecasts\*\* (all envs).



---



\# 13) Life-cycle base parameters \& calibration targets



\*\*What it does\*\*



\* Defines `base\_P` with fecundity, female fraction, density-dependence, K, SAR mean/sd (mean set later), age-return probs, rear survival (set later).

\* Sets calibration targets:



&nbsp; \* `obs\_spawners` (2011–2024),

&nbsp; \* `S\_seed\_calib` (first 3 obs years),

&nbsp; \* index of fit years (4–14).



\*\*Why\*\*



\* Spawner recursion needs these; calibration picks `SAR\_mean` and `rear\_surv` that best fit escapement.



---



\# 14) Calibrate SAR\\\_mean and rear\\\_surv (per variant)



\*\*What it does\*\*



\* Minimizes SSE on years 4–14 between predicted and observed spawners:



&nbsp; \* Runs `simulate\_variant()` using:



&nbsp;   \* `surv\_vec` from `surv\_lookup\_by\_variant\[\[variant]]`,

&nbsp;   \* `deg\_day\_adult = 0` (i.e., pre-spawn survival = plogis(θ0) only),

&nbsp;   \* K fixed,

&nbsp;   \* SAR constant = parameter being optimized.

\* L-BFGS-B over `\[0,1]×\[0,1]` for `(SAR\_mean, rear\_surv)`.



\*\*Why\*\*



\* Fits the aggregate pieces of the life cycle to the observed spawner series, \*holding egg→fry survival fixed\* to the TDM outputs for the calibration env.



**\*\*Note\*\***



**\* Pre-spawn thermal effects are \*\*not\*\* used in calibration (`deg\_day\_adult = 0`). That’s fine if intentional; otherwise pass real degree-days for 2011–2024 (per env) for a tighter match.**



---



\# 15) Build env-specific parameter lists with calibrated values



\*\*What it does\*\*



\* `base\_P\_list`: for each variant, clones `base\_P` across envs and fills in the \*\*calibrated\*\* `(SAR\_mean, rear\_surv)` for that variant.



\*\*Why\*\*



\* Lets every env share the same calibrated life-cycle rates per variant, differing only by their survival vector and degree-days.



---



\# 16) Precompute calibration-tab predictions (for Shiny)



\*\*What it does\*\*



\* For each variant, simulates 2011–2024 with the \*\*calibrated\*\* parameters to produce a table with `year, observed, predicted, SAR\_mean, rear\_surv, sse`.



\*\*Why\*\*



\* The Shiny Calibration tab can show this instantly (no re-sim).



---



\# 17) Build seeds for forecast runs



\*\*What it does\*\*



\* For each variant, re-simulates 2011–2024 with:



&nbsp; \* `surv\_vec` from `surv\_lookup\_full\[\[paste(ref\_env, v)]]` (same numbers),

&nbsp; \* \*\*degree-days now computed\*\* via `compute\_deg\_day\_adult(ref\_env, real\_years, spawn\_dates\_vec, env\_ext\_list)`,

\* Takes the \*\*last 3 spawner values\*\* as forecast seeds for that variant (`S\_seed\_fore\_list`).



\*\*Why\*\*



\* Forecasts need initial S in years 1–3; using the calibrated trajectory’s tail is standard.



---



\# 18) Forecasts for every env × variant



\*\*What it does\*\*



\* `keys <- names(surv\_lookup\_full)` gives you all `"env\_variant"` pairs.

\* For each key:



&nbsp; \* creates a simulator `sim\_forecast\_fn(var, env, flow\_cfs=NULL, S\_seed = S\_seed\_fore\_list\[\[var]])`,

&nbsp; \* runs it to produce forward series using:



&nbsp;   \* `surv\_vec` from `surv\_lookup\_full`,

&nbsp;   \* `deg\_day\_adult` from `compute\_deg\_day\_adult(env, sim\_years, spawn\_dates\_vec, env\_ext\_list)`,

&nbsp;   \* `SAR\_vec` \*\*stochastic\*\* (because `use\_stochastic\_SAR <- TRUE`); mean is set internally to that env/variant’s `SAR\_mean` when the simulator builds the options.



\*\*Why\*\*



\* Produces the master “precomputed” forecast data frame for Shiny.



**\*\*Note\*\***



**\* Because `use\_stochastic\_SAR` is `TRUE`, your saved `results\_full` is a stochastic draw. If you want deterministic defaults in the app, set it to `FALSE` here and let the UI control stochasticity later.**



---



\# 19) Save all artifacts for Shiny



\*\*What it does\*\*



\* Saves RDS files: calibration results and predictions, full forecast results, survival summaries/lookups, parameter lists, seeds, SAR options, `sim\_years`, and `spawn\_dates\_vec`.



\*\*Why\*\*



\* Shiny loads these so it can render instantly and only recompute what the user changes interactively.



---



\## Key “does it do what I want?” checks \& tweaks



1\. \*\*Env-specific spawn dates in degree-day calc\*\*



\* Current: `spawn\_dates\_vec` is \*one date per year\*, applied to \*\*all envs\*\* in `compute\_deg\_day\_adult()`.

\* If you want env-specific degree-day end dates, compute `spawn\_dates\_vec` \*\*by env\*\* and pass the env’s vector when calling `compute\_deg\_day\_adult()`.



2\. \*\*Fallback `env\_sites` is undefined\*\*



\* In the sim loop for future redds, if a site pool is missing you call `env\_sites` (not defined). Replace with:



&nbsp; ```r

&nbsp; st\_pool <- ft\_aug$site\_pool\[\[i]]

&nbsp; if (is.null(st\_pool)) {

&nbsp;   st\_pool <- site\_env\_map %>% filter(env == env\_i) %>% pull(site) %>% unique()

&nbsp; }

&nbsp; ```



3\. \*\*Window length for `compute\_surv()`\*\*



\* Uses \*\*T0-based\*\* length `ceil(958/T0 + 417/T0)`. Fast, but could miss true ATU end if temps drift. If you need exactness:



&nbsp; \* Build the slice by extending days until `cumsum(temp) >= 1375`.

&nbsp; \* Then let `tdm\_exp()` split egg/alevin internally (it already does).



4\. \*\*Calibration and pre-spawn thermal stress\*\*



\* You set `deg\_day\_adult = 0` in calibration. If you want pre-spawn effects “in the fit,” compute and pass degree-days for 2011–2024 (per env).



5\. \*\*Reproducibility\*\*



\* You set `seed=TRUE` for `furrr`, but bin sampling \& `runif()` in `sample\_dates\_fast()` and SAR draws aren’t globally seeded here. Add `set.seed()` near the top if you need bitwise repeatability.



6\. \*\*CLMM standardization\*\*



\* Forecast z-scores use `oct\_mean/sd` and `nov\_mean/sd` from your \*historical\* carcass-joined data—good. Ensure this is the same centering/scaling used at fit time (it is, since those came from the same data). If the fit used different scaling (e.g., within env), align them.



7\. \*\*Reference environment\*\*



\* `ref\_env <- names(env\_ext\_list)\[1]`. If your “calibration env” should be a specific alternative, set it explicitly (e.g., `"NoAction"`).



8\. \*\*Unused objects\*\*



\* `site\_props` is computed but never used. Safe to remove unless you plan to weight site sampling by it.



9\. \*\*Results variability\*\*



\* `results\_full` depends on `use\_stochastic\_SAR`. Decide whether the app should boot with deterministic (set FALSE here) and let users switch on stochastic in the UI.



---



If you want, I can sketch the minimal code diffs to (a) fix `env\_sites`, (b) make env-specific `spawn\_dates\_vec`, and (c) add an ATU-precise `compute\_surv\_by\_atu()`—just say the word and I’ll drop in patch-ready snippets.



