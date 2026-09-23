# years.R ---------------------------------------------------------------------
# Analysis-year registry for the SalmonCountR app.
#
# One place that declares what changes between analysis years, so adding a year
# is a data drop plus one entry in ARG_YEARS below -- not an edit to app.R.
#
# WHAT AN ANALYSIS YEAR IS. The vintage of the modelling deliverable the app is
# showing: the CE-QUAL-W2 scenario temperatures, the precomputed forecast built
# from them, and the design inputs (hydropower replacement cost, elicited
# objective weights) that belong to that round of the decision process. It is
# NOT the simulation year inside a projection -- see first_projection_year.
#
# HOW TO ADD A YEAR.
#   1. Put that year's precompute output in the app's app_data/<year>/.
#      The 2025 bundle is the flat app_data/ directory for historical reasons;
#      every later year gets its own subdirectory. Nothing has to move.
#   2. Add an entry to ARG_YEARS with its weights, hydropower costs and first
#      projection year.
#   3. That is all. The selector picks it up, and a year whose files are not
#      present yet shows in the selector as unavailable rather than erroring.
#
# WHAT IS NOT WIRED. functions.R reads several objects from the global
# environment (env_ext_list, sim_years, spawn_dates_by_alt and others). The app
# does not call those functions -- it reads precomputed results and calls only
# get_scenario_alternatives(), which is pure -- so the year switch is complete
# for the app as it stands. Anyone who later makes the app run simulations live
# must parameterise functions.R first, or those paths will silently use whatever
# global.R loaded at startup.
# -----------------------------------------------------------------------------

# ---- App directory resolution -----------------------------------------------
# One rule that holds locally and on the Shiny server, so the same code deploys
# without being restructured first.
#
# WHAT WAS WRONG. This used here::here("SalmonCountR", ...). here() ignores the
# working directory on purpose and walks up looking for a project marker
# (.git / .Rproj). Locally that lands on the repo root, so prefixing
# "SalmonCountR" is correct. In a deployed bundle there is no marker, so here()
# falls back to the working directory -- which Shiny has ALREADY set to the app
# directory -- and the "SalmonCountR" prefix pointed one level too deep. Hence
# the manual reshuffle before every deploy.
#
# WHAT IT DOES NOW. Anchor on the directory that actually contains app_data/,
# probing the two places it can be. Resolved once and cached.
#
#   deployed / runApp()      wd is the app dir          -> "."
#   sourced from repo root   wd is the repo root        -> "SalmonCountR"
#
# Add a candidate here if a third layout ever appears; do not reintroduce here().
.arg_app_dir_cache <- NULL

arg_app_dir <- function() {
  if (!is.null(.arg_app_dir_cache)) return(.arg_app_dir_cache)
  candidates <- c(".", "SalmonCountR")
  hit <- Find(function(p) dir.exists(file.path(p, "app_data")), candidates)
  if (is.null(hit)) {
    stop("Cannot locate the SalmonCountR app directory: no app_data/ found ",
         "relative to '", getwd(), "'. Looked in: ",
         paste(candidates, collapse = ", "), call. = FALSE)
  }
  .arg_app_dir_cache <<- normalizePath(hit, winslash = "/", mustWork = TRUE)
  .arg_app_dir_cache
}

arg_app_path <- function(...) file.path(arg_app_dir(), ...)

# The objects the app actually reads, and the file each comes from.
ARG_YEAR_FILES <- c(
  results_full               = "results_full.rds",
  steelhead_metrics          = "steelhead_metrics.rds",
  swing_ranges               = "swing_ranges.rds",
  instream                   = "american_river_instream.rds",
  df_all_orig                = "df_all.rds",
  swing_scenario_results     = "swing_scenario_results.rds",
  steelhead_scenario_results = "steelhead_scenario_results.rds",
  alt_key                    = "alt_key.rds"   # which run is which alternative; see functions.R
)

# Hydropower replacement cost ($) by alternative. A design input, not a model
# output, so it is declared per year rather than derived. Every alternative in
# the year's alt_key must have a cost here, or the year shows as not loaded
# (arg_year_missing) -- an alternative without a cost would otherwise look free.
ARG_HYDRO_COST_2025 <- c(
  NB   = 0,      PB1 = 111422, PB2 = 376671,
  PB2b = 470090, PB2c = 433215, PB3 = 201552,
  PB4  = 241590, PB5 = 199382, PB6 = 348806
)

# 2026 draft deliverable (TemperatureModelingResults_9-23-26.xlsx): Scenarios
# 1-4 on the ATSP 41 schedule, No Bypass on ATSP 41 and ATSP 38, with Scenarios
# 1 and 2 on ATSP 38 still to come. PLACEHOLDERS: no 2026 valuation exists yet.
# Each scenario carries its 2025 cost on the assumption that the bypass
# schedule, and so the foregone generation, is unchanged; an ATSP variant costs
# the same as its base scenario because ATSP is a shutter schedule, not a
# bypass. No Bypass costs nothing under either schedule. Replace when Reclamation
# values the 2026 alternatives.
ARG_HYDRO_COST_2026 <- c(
  NB      = 0,
  `NB-38` = 0,
  PB1     = ARG_HYDRO_COST_2025[["PB1"]],
  PB2     = ARG_HYDRO_COST_2025[["PB2"]],
  PB3     = ARG_HYDRO_COST_2025[["PB3"]],
  PB4     = ARG_HYDRO_COST_2025[["PB4"]],
  `PB1-38` = ARG_HYDRO_COST_2025[["PB1"]],
  `PB2-38` = ARG_HYDRO_COST_2025[["PB2"]]
)

# OBJECTIVE SCALING. Decision Support and Swing Weighting put each objective on
# a 0-1 scale before weighting. Two ways to do that:
#
#   local   0 = the worst of the alternatives in the analysis, 1 = the best.
#           The 2025 analysis did this (min-max over the nine), and its
#           published weights were elicited against those swings. Adding an
#           alternative rescales every other one.
#   global  0 and 1 are FIXED ends of a plausible range for the objective,
#           declared here per year. Scores stay comparable across years and
#           across additions to the alternative set, and the swing weights are
#           elicited against these same ranges. Values outside the range are
#           clamped.
#
# A year with `objective_ranges` scales globally; a year without scales locally.
# Ranges are c(lo, hi) in the objective's raw units: Chinook adult index,
# steelhead days below 18.3 C in Oct-Nov (at most 61), hydropower replacement
# cost in $ (lower is better; the scale is inverted for it).
ARG_OBJECTIVE_RANGES_2026 <- list(
  chinook   = c(0, 25000),
  steelhead = c(0, 61),
  hydro     = c(0, 3e6)
)

ARG_YEARS <- list(
  "2025" = list(
    label                 = "2025",
    dir                   = "app_data",
    default_weights       = c(chinook = 0.40, steelhead = 0.10, hydro = 0.50),
    hydro_cost            = ARG_HYDRO_COST_2025,
    first_projection_year = 2025,
    temperature_year      = 2025,
    # No objective_ranges: the published analysis scales locally, and stays so.
    note                  = "Published analysis. Elicited weights from the 2025 SDM workshop."
  ),
  "2026" = list(
    label                 = "2026 (draft)",
    dir                   = file.path("app_data", "2026"),
    # Weights carried forward from 2025 as a starting point. Replace when the
    # 2026 elicitation is done -- these are placeholders, not results.
    default_weights       = c(chinook = 0.40, steelhead = 0.10, hydro = 0.50),
    hydro_cost            = ARG_HYDRO_COST_2026,
    # Global scaling on the ranges B. Mahardja proposed in September 2026 as a
    # starting point (salmon 0-25,000; hydro $0-3M; steelhead the full Oct-Nov
    # window). The 2025 weights above were elicited against LOCAL swings, so
    # they are only placeholders here until re-elicited against these ranges.
    objective_ranges      = ARG_OBJECTIVE_RANGES_2026,
    # The 2026 temperatures are run through the 2025 model: calibration ends in
    # 2024 and the projection starts in 2025, with the decision date held at
    # 2025-09-21 so that the first projection year carries the scenario
    # temperatures. See app_data/2026/README.md for what a true 2026 start needs.
    first_projection_year = 2025,
    # The deliverable models the fall of 2026, and that is what the Temperature
    # Explorer shows, even though the population model projects from 2025.
    temperature_year      = 2026,
    hydro_cost_note       = "Placeholders: each scenario carries its 2025 cost, ATSP variants their base scenario's, No Bypass $0. Replace when the 2026 alternatives are valued.",
    note                  = "Draft 2026 temperature deliverable (23 Sep 2026) run through the 2025 model. Objective weights and hydropower costs are placeholders carried over from 2025."
  )
)

#' The calendar year a year's temperature deliverable models, for the
#' Temperature Explorer. Falls back to the first projection year, which is what
#' it is whenever the model has been recalibrated up to the decision.
arg_temperature_year <- function(cfg) {
  if (!is.null(cfg$temperature_year)) cfg$temperature_year else cfg$first_projection_year
}

ARG_DEFAULT_YEAR <- "2025"

arg_year_ids <- function() names(ARG_YEARS)

arg_year_cfg <- function(year) {
  year <- as.character(year)
  if (!year %in% names(ARG_YEARS)) {
    stop("Unknown analysis year: ", year, call. = FALSE)
  }
  ARG_YEARS[[year]]
}

arg_year_dir <- function(year) arg_app_path(arg_year_cfg(year)$dir)

#' A year's objective ranges, validated, or NULL for local scaling.
arg_objective_ranges <- function(cfg) {
  r <- cfg$objective_ranges
  if (is.null(r)) return(NULL)
  need <- c("chinook", "steelhead", "hydro")
  if (!all(need %in% names(r))) {
    stop("objective_ranges for ", cfg$label, " must name ", paste(need, collapse = ", "),
         call. = FALSE)
  }
  for (n in need) {
    if (length(r[[n]]) != 2 || !all(is.finite(r[[n]])) || r[[n]][2] <= r[[n]][1]) {
      stop("objective_ranges$", n, " for ", cfg$label, " must be c(lo, hi) with hi > lo",
           call. = FALSE)
    }
  }
  r[need]
}

#' 0-1 score for one objective under a year's scaling. `range` NULL means
#' local min-max over `x`; otherwise the fixed range, clamped. `lower_better`
#' inverts the scale (hydropower cost).
arg_scale_objective <- function(x, range = NULL, lower_better = FALSE) {
  if (is.null(range)) {
    lo <- min(x, na.rm = TRUE); hi <- max(x, na.rm = TRUE)
    if (!is.finite(lo) || !is.finite(hi) || hi == lo) return(rep(0.5, length(x)))
  } else {
    lo <- range[1]; hi <- range[2]
  }
  s <- (x - lo) / (hi - lo)
  if (lower_better) s <- 1 - s
  pmin(pmax(s, 0), 1)
}

#' One line saying how a year's objectives are scaled, for the app.
arg_scaling_note <- function(cfg) {
  r <- arg_objective_ranges(cfg)
  if (is.null(r)) {
    return(paste0("Objectives are scaled locally: 0 is the worst and 1 the best of the ",
                  "alternatives in the ", cfg$label, " analysis."))
  }
  fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  paste0("Objectives are scaled on fixed (global) ranges for ", cfg$label, ": Chinook ",
         fmt(r$chinook[1]), "-", fmt(r$chinook[2]),
         " adults; steelhead ", r$steelhead[1], "-", r$steelhead[2], " days below 18.3 °C; ",
         "hydropower $", fmt(r$hydro[1]), "-$", fmt(r$hydro[2]), " (lower is better). ",
         "Values outside a range are clamped.")
}

#' Which required files a year is missing. character(0) when the year is ready.
#' A year whose alt_key names an alternative with no hydropower cost in
#' years.R counts as not ready too, and says so, rather than loading with that
#' alternative looking free in Decision Support.
arg_year_missing <- function(year) {
  d <- arg_year_dir(year)
  if (!dir.exists(d)) return(unname(ARG_YEAR_FILES))
  missing <- ARG_YEAR_FILES[!file.exists(file.path(d, ARG_YEAR_FILES))] |> unname()
  if (!"alt_key.rds" %in% missing) {
    key <- readRDS(file.path(d, "alt_key.rds"))
    no_cost <- setdiff(unique(key$alt), names(arg_year_cfg(year)$hydro_cost))
    if (length(no_cost)) {
      missing <- c(missing, paste0("hydro_cost in years.R for ", paste(no_cost, collapse = ", ")))
    }
  }
  missing
}

arg_year_available <- function(year) length(arg_year_missing(year)) == 0

arg_years_available <- function() {
  Filter(arg_year_available, arg_year_ids())
}

#' Choices for the year selector. Unavailable years stay visible and are
#' labelled, rather than being hidden -- someone looking for 2026 should be able
#' to see that the app knows about it and is waiting on data.
arg_year_choices <- function() {
  ids <- arg_year_ids()
  labs <- vapply(ids, function(y) {
    lab <- arg_year_cfg(y)$label
    if (arg_year_available(y)) lab else paste0(lab, " (data not loaded)")
  }, character(1))
  stats::setNames(ids, labs)
}

#' Provenance stamp for a year, or NULL. Written by analysis/refresh_data_year.R.
arg_year_vintage <- function(year) {
  p <- file.path(arg_year_dir(year), "data_vintage.rds")
  if (file.exists(p)) readRDS(p) else NULL
}

# ---- Bundle loading ---------------------------------------------------------

# Loaded bundles are cached per year. Reading them is cheap, but the year
# reactive can invalidate often and there is no reason to hit disk each time.
.arg_bundle_cache <- new.env(parent = emptyenv())

#' Post-load preparation shared by every year.
#'
#' Kept here rather than in global.R so a year switched into at runtime gets
#' exactly the same treatment the startup year got.
arg_prepare_bundle <- function(raw, cfg) {
  # The alternative key drives everything that used to count to nine: the
  # env -> alternative -> met year mapping is data, not arithmetic.
  key <- as.data.frame(raw$alt_key, stringsAsFactors = FALSE)
  key$env <- as.character(key$env)
  raw$alt_key   <- key
  raw$alt_codes <- unique(key$alt)
  raw$met_years <- unique(key$met_year)
  raw$alt_labels <- stats::setNames(
    vapply(raw$alt_codes, function(a) unique(key$label[key$alt == a])[1], character(1)),
    raw$alt_codes)

  # Fixed normalisation bounds for the Chinook objective, computed across all
  # alternatives AND all three TDM models.
  #
  # WHY THIS EXISTS. The MCDA tab used to min-max the Chinook scores within
  # whatever the current TDM weighting produced. Push the TDM weights far enough
  # and the scale moves underneath the composite, which breaks the scoring:
  # swing weights are defined relative to each objective's range, so a moving
  # range means the elicited weights no longer describe the trade-off that was
  # elicited. Same defect B. Mahardja found in the weight-sensitivity analysis;
  # same fix, and the same bounds analysis/evpi.R and
  # analysis/tdm_weight_sensitivity.R use.
  #
  # NOTE this changes the composite values the app displays relative to the
  # published Table/Figure 5, which min-max within the nine alternatives at the
  # elicited TDM weighting. That is a deliberate divergence: the app has to stay
  # correct while a user moves the weights, the paper does not.
  raw$salmon_bounds <- local({
    rf <- raw$results_full
    if (is.null(rf) || !all(c("env", "variant", "year", "spawners") %in% names(rf))) {
      return(NULL)
    }
    met_w <- 1 / length(raw$met_years)
    per_state <- rf |>
      dplyr::filter(year > 2024) |>
      dplyr::group_by(env, variant) |>
      dplyr::slice_tail(n = 20) |>
      dplyr::summarise(med = stats::median(spawners, na.rm = TRUE), .groups = "drop") |>
      dplyr::mutate(env = as.character(env)) |>
      dplyr::inner_join(key[, c("env", "alt")], by = "env") |>
      dplyr::group_by(alt, variant) |>
      dplyr::summarise(value = sum(med * met_w), .groups = "drop")
    c(lo = min(per_state$value), hi = max(per_state$value))
  })

  raw$get_K_spawners <- local({
    tbl <- dplyr::mutate(raw$instream, K_spawners = FR_spawn_wua / 9.29)
    interp <- stats::approxfun(tbl$flow_cfs, tbl$K_spawners, rule = 2)
    function(flow_vec) interp(flow_vec)
  })

  # Temperature Explorer data, precomputed once. The filter is on the year the
  # deliverable models (temperature_year), which is the first projection year
  # unless the model is lagging the deliverable -- it used to be a literal
  # 2025, which would have returned zero rows for any later deliverable.
  raw$df_temp_first_year <- if (is.data.frame(raw$df_all_orig)) {
    raw$df_all_orig |>
      dplyr::filter(lubridate::year(Date) == arg_temperature_year(cfg)) |>
      dplyr::mutate(month_num = lubridate::month(Date), env = as.character(env)) |>
      dplyr::left_join(
        stats::setNames(key[, c("env", "met_year")], c("env", "climate")), by = "env")
  } else NULL

  raw$cfg     <- cfg
  raw$vintage <- NULL   # filled by load_year_bundle, which knows the year
  raw
}

#' Load one analysis year. Returns NULL if the year's files are not present,
#' so callers can render a message instead of failing.
load_year_bundle <- function(year, refresh = FALSE) {
  year <- as.character(year)
  if (!refresh && !is.null(.arg_bundle_cache[[year]])) return(.arg_bundle_cache[[year]])
  if (!arg_year_available(year)) return(NULL)

  d   <- arg_year_dir(year)
  cfg <- arg_year_cfg(year)

  raw <- lapply(ARG_YEAR_FILES, function(f) readRDS(file.path(d, f)))
  names(raw) <- names(ARG_YEAR_FILES)

  bundle <- arg_prepare_bundle(raw, cfg)
  bundle$year    <- year
  bundle$vintage <- arg_year_vintage(year)

  .arg_bundle_cache[[year]] <- bundle
  bundle
}
