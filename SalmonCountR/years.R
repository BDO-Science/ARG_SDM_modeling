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
#   1. Put that year's precompute output in SalmonCountR/app_data/<year>/.
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

# Paths are built exactly the way global.R has always built them, so deployment
# behaviour is unchanged. If the here() convention is ever fixed, fix it here.
arg_app_path <- function(...) here::here("SalmonCountR", ...)

# The objects the app actually reads, and the file each comes from.
ARG_YEAR_FILES <- c(
  results_full               = "results_full.rds",
  steelhead_metrics          = "steelhead_metrics.rds",
  swing_ranges               = "swing_ranges.rds",
  instream                   = "american_river_instream.rds",
  df_all_orig                = "df_all.rds",
  swing_scenario_results     = "swing_scenario_results.rds",
  steelhead_scenario_results = "steelhead_scenario_results.rds"
)

# Hydropower replacement cost ($) by alternative. A design input, not a model
# output, so it is declared per year rather than derived. The 2025 values are
# the ones app.R carried hard-coded.
ARG_HYDRO_COST_2025 <- c(
  NB   = 0,      PB1 = 111422, PB2 = 370826,
  PB2b = 470090, PB2c = 433215, PB3 = 201552,
  PB4  = 241590, PB5 = 199382, PB6 = 348806
)

ARG_YEARS <- list(
  "2025" = list(
    label                 = "2025",
    dir                   = "app_data",
    default_weights       = c(chinook = 0.40, steelhead = 0.10, hydro = 0.50),
    hydro_cost            = ARG_HYDRO_COST_2025,
    first_projection_year = 2025,
    note                  = "Published analysis. Elicited weights from the 2025 SDM workshop."
  ),
  "2026" = list(
    label                 = "2026",
    dir                   = file.path("app_data", "2026"),
    # Carried forward from 2025 as a starting point. Replace when the 2026
    # elicitation and valuation are done -- these are placeholders, not results.
    default_weights       = c(chinook = 0.40, steelhead = 0.10, hydro = 0.50),
    hydro_cost            = ARG_HYDRO_COST_2025,
    first_projection_year = 2026,
    note                  = "Awaiting the 2026 temperature deliverable and precompute run. Weights and hydropower costs are placeholders carried over from 2025."
  )
)

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

#' Which required files a year is missing. character(0) when the year is ready.
arg_year_missing <- function(year) {
  d <- arg_year_dir(year)
  if (!dir.exists(d)) return(unname(ARG_YEAR_FILES))
  ARG_YEAR_FILES[!file.exists(file.path(d, ARG_YEAR_FILES))] |> unname()
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
  # Fixed normalisation bounds for the Chinook objective, computed across all
  # nine alternatives AND all three TDM models.
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
    per_state <- rf |>
      dplyr::filter(year > 2024) |>
      dplyr::group_by(env, variant) |>
      dplyr::slice_tail(n = 20) |>
      dplyr::summarise(med = stats::median(spawners, na.rm = TRUE), .groups = "drop") |>
      dplyr::mutate(env = as.integer(env),
                    scenario = ((env - 1) %% 9) + 1) |>
      dplyr::group_by(scenario, variant) |>
      dplyr::summarise(value = sum(med * 0.25), .groups = "drop")
    c(lo = min(per_state$value), hi = max(per_state$value))
  })

  raw$get_K_spawners <- local({
    tbl <- dplyr::mutate(raw$instream, K_spawners = FR_spawn_wua / 9.29)
    interp <- stats::approxfun(tbl$flow_cfs, tbl$K_spawners, rule = 2)
    function(flow_vec) interp(flow_vec)
  })

  # Temperature Explorer data, precomputed once. The filter is on the first
  # projection year of THIS bundle -- it used to be a literal 2025, which would
  # have returned zero rows for any later deliverable.
  raw$df_temp_first_year <- if (is.data.frame(raw$df_all_orig)) {
    raw$df_all_orig |>
      dplyr::filter(lubridate::year(Date) == cfg$first_projection_year) |>
      dplyr::mutate(
        month_num = lubridate::month(Date),
        climate = dplyr::case_when(
          env %in% as.character(1:9)   ~ "2011",
          env %in% as.character(10:18) ~ "2014",
          env %in% as.character(19:27) ~ "2017",
          env %in% as.character(28:36) ~ "2020"
        )
      )
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
