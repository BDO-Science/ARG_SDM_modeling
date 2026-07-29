# ============================================================================
# Scenario engine -- evaluate a new temperature deliverable without re-running
# the pipeline
# ============================================================================
# Job A of output/APP_DATA_UPDATE_OPTIONS.md: let someone drop in the
# temperature modelling deliverable exactly as the modelling team sends it, and
# get consequence-table and MCDA results back, with no R and no precompute.
#
# Why this is cheap. Three facts established 2026-07-28:
#
#   1. Calibration does not depend on scenario temperatures. SAR and rear_surv
#      are fitted against observed escapement using only the reference
#      alternative, so new scenarios change the forecast and nothing else.
#   2. The forecast temperature series is exactly periodic from 2026 onward --
#      every projection year has an identical daily profile (verified: 0.000 degC
#      difference between years 2030, 2031 and 2075). So egg-to-fry survival
#      need only be computed for ONE season per scenario, not 114.
#   3. Survival over an incubation window is a difference of cumulative sums,
#      so it vectorises. precompute.R is slow because it loops over 114 years x
#      36 alternatives x ~400 individual redds with a memoised vapply.
#
# Together those make a nine-scenario run sub-second rather than minutes.
#
# DESIGN CONSTRAINT: every function here is pure. Nothing assigns into the
# global environment. The Shiny app is deployed on shinyapps.io, where several
# sessions share one R process, so mutating globals would leak one user's
# uploaded scenario into another user's results.
#
# Deterministic by construction. precompute.R samples individual redds from the
# spawn-timing CLM; this weights every spawn date by its CLM probability
# instead. Same model, expected value rather than one draw, so repeat runs of
# the same file give identical answers.
#
# Requires: SalmonCountR/functions.R, and app_data/spawn_timing_model.rds
#           (built by analysis/build_spawn_timing_model.R).
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(readxl); library(lubridate)
})

# Scenario labels as they appear in row 1 of a year sheet, and the short codes
# used everywhere else in the project.
SCENARIO_LABELS <- c("No Bypass", "Scenario 1", "Scenario 2", "Scenario 2b",
                     "Scenario 2c", "Scenario 3", "Scenario 4", "Scenario 5",
                     "Scenario 6")
SCENARIO_CODES  <- c("NB", "PB1", "PB2", "PB2b", "PB2c", "PB3", "PB4", "PB5", "PB6")
ENGINE_SITES    <- c("AveWatt", "AveHazel")

# Window in which the deliverable overrides climatology, from temperature_data.R
OVERRIDE_DOY <- c(start = 291L, end = 365L)   # 18 Oct - 31 Dec
TEMP_FLOOR   <- c(AveHazel = 7, AveWatt = 8)  # gauge floors applied in temperature_data.R

# Plausible daily water temperature for the lower American River in autumn.
# Outside this a file is almost certainly in Fahrenheit or mis-parsed.
TEMP_PLAUSIBLE <- c(4, 30)


# ---- 1. Read -----------------------------------------------------------------

#' Read a temperature modelling deliverable
#'
#' Handles the layout the modelling team ships: one sheet per meteorological
#' year, scenario names merged across row 1, column headers on row 2, then
#' Date / JDAY followed by AveWatt-AveHazel pairs per scenario. Trailing
#' columns such as "Target Temp" are ignored.
#'
#' @param path Path to the .xlsx.
#' @return list(temps, met_years, scenarios, sheets, path, read_errors)
read_deliverable <- function(path) {
  empty <- function(errs) list(temps = tibble(), met_years = character(),
                               scenarios = character(), sheets = character(),
                               path = path, read_errors = errs)

  if (!file.exists(path)) return(empty("The file could not be found."))

  # A user can upload anything. Fail with a message rather than an R error.
  sheets <- tryCatch(readxl::excel_sheets(path), error = function(e) NULL)
  if (is.null(sheets)) {
    return(empty(paste("This is not a readable Excel workbook. Save the deliverable",
                       "as .xlsx and upload it again.")))
  }

  # Year sheets are those named as a four-digit year
  year_sheets <- sheets[grepl("^\\d{4}$", sheets)]
  if (!length(year_sheets)) {
    return(empty(sprintf(paste("No sheet is named as a four-digit year. Sheets found:",
                               "%s. Each meteorological year needs its own sheet,",
                               "named for example '2026'."),
                         paste(sheets, collapse = ", "))))
  }
  read_errors <- character()

  temps <- purrr::map_dfr(year_sheets, function(sh) {
    hdr <- tryCatch(
      readxl::read_excel(path, sheet = sh, n_max = 2, col_names = FALSE,
                         .name_repair = "minimal"),
      error = function(e) NULL)
    dat <- tryCatch(
      readxl::read_excel(path, sheet = sh, skip = 1, .name_repair = "minimal"),
      error = function(e) NULL)
    if (is.null(hdr) || is.null(dat) || nrow(dat) == 0) {
      read_errors <<- c(read_errors, sprintf("Sheet '%s' could not be read.", sh))
      return(tibble())
    }

    # Row 1 carries the scenario name above the first column of each pair, with
    # NA in the second. Carry it forward to label both columns.
    lab <- as.character(unlist(hdr[1, ], use.names = FALSE))
    lab <- vctrs::vec_fill_missing(lab, direction = "down")
    site_row <- as.character(unlist(hdr[2, ], use.names = FALSE))

    keep <- which(site_row %in% ENGINE_SITES & !is.na(lab))
    if (!length(keep)) {
      read_errors <<- c(read_errors,
                        sprintf("Sheet '%s' has no AveWatt/AveHazel columns.", sh))
      return(tibble())
    }

    dates <- suppressWarnings(as.Date(dat[[1]]))
    purrr::map_dfr(keep, function(j) {
      if (j > ncol(dat)) return(tibble())
      tibble(met_year = sh,
             scenario_label = lab[j],
             site = site_row[j],
             Date = dates,
             temp = suppressWarnings(as.numeric(dat[[j]])))
    })
  })

  if (nrow(temps)) {
    temps <- temps %>%
      filter(!is.na(Date)) %>%
      mutate(scenario = SCENARIO_CODES[match(scenario_label, SCENARIO_LABELS)],
             doy = lubridate::yday(Date))
  }

  list(temps       = temps,
       met_years   = sort(unique(temps$met_year)),
       scenarios   = unique(temps$scenario_label),
       sheets      = sheets,
       path        = path,
       read_errors = read_errors)
}


# ---- 2. Validate -------------------------------------------------------------

#' Check a parsed deliverable and report problems in plain language
#'
#' @return tibble(check, status in {ok, warning, error}, message)
validate_deliverable <- function(parsed) {
  add <- function(acc, check, status, message) {
    dplyr::bind_rows(acc, tibble(check = check, status = status, message = message))
  }
  out <- tibble(check = character(), status = character(), message = character())
  tm  <- parsed$temps

  for (e in parsed$read_errors) out <- add(out, "File structure", "error", e)

  if (!nrow(tm)) {
    return(add(out, "File structure", "error",
               paste("No temperature data could be read. Expected one sheet per",
                     "meteorological year, named as a four-digit year (for example",
                     "'2026'), with scenario names on the first row and",
                     "AveWatt / AveHazel column headers on the second.")))
  }

  out <- add(out, "Meteorological years", "ok",
             sprintf("Found %d: %s.", length(parsed$met_years),
                     paste(parsed$met_years, collapse = ", ")))

  # Scenario names
  unknown <- setdiff(parsed$scenarios, SCENARIO_LABELS)
  missing <- setdiff(SCENARIO_LABELS, parsed$scenarios)
  if (length(unknown)) {
    out <- add(out, "Scenario names", "error",
               sprintf(paste("These scenario names were not recognised: %s.",
                             "Expected exactly: %s."),
                       paste(unknown, collapse = ", "),
                       paste(SCENARIO_LABELS, collapse = ", ")))
  }
  if (length(missing)) {
    out <- add(out, "Scenario names", "warning",
               sprintf(paste("These scenarios are absent and will be skipped: %s.",
                             "Results will cover only the scenarios present."),
                       paste(missing, collapse = ", ")))
  }
  if (!length(unknown) && !length(missing)) {
    out <- add(out, "Scenario names", "ok", "All nine scenarios present and named as expected.")
  }

  # Sites
  per <- tm %>% count(met_year, scenario, site) %>% count(met_year, scenario, name = "n_site")
  if (any(per$n_site < 2)) {
    bad <- per %>% filter(n_site < 2)
    out <- add(out, "Monitoring sites", "error",
               sprintf("Missing an AveWatt or AveHazel column for: %s.",
                       paste(sprintf("%s %s", bad$met_year, bad$scenario), collapse = "; ")))
  } else {
    out <- add(out, "Monitoring sites", "ok", "Both AveWatt and AveHazel present throughout.")
  }

  # Date coverage against the window that actually drives incubation
  cov <- tm %>%
    group_by(met_year) %>%
    summarise(from = min(Date), to = max(Date),
              n_days = n_distinct(Date), .groups = "drop")
  out <- add(out, "Date coverage", "ok",
             paste(sprintf("%s: %s to %s (%d days)", cov$met_year,
                           format(cov$from, "%d %b"), format(cov$to, "%d %b"),
                           cov$n_days), collapse = "; "))

  overlap <- tm %>%
    filter(doy >= OVERRIDE_DOY[["start"]], doy <= OVERRIDE_DOY[["end"]]) %>%
    count(met_year, name = "n")
  if (!nrow(overlap) || any(overlap$n == 0)) {
    out <- add(out, "Incubation window", "error",
               paste("The file covers no dates between 18 October and 31 December.",
                     "That is the window in which the scenario temperatures are used,",
                     "so nothing would change."))
  } else {
    out <- add(out, "Incubation window", "ok",
               "Covers the 18 October - 31 December window used for incubation.")
  }

  # Gaps within each series
  gaps <- tm %>%
    group_by(met_year, scenario, site) %>%
    arrange(Date, .by_group = TRUE) %>%
    summarise(gap = sum(as.integer(diff(Date)) > 1), .groups = "drop") %>%
    filter(gap > 0)
  if (nrow(gaps)) {
    out <- add(out, "Missing days", "warning",
               sprintf(paste("%d series have gaps in their dates. Missing days fall",
                             "back to the long-term average, which may understate",
                             "differences between scenarios."), nrow(gaps)))
  } else {
    out <- add(out, "Missing days", "ok", "No gaps in any daily series.")
  }

  # Values
  nas <- sum(is.na(tm$temp))
  if (nas > 0) {
    out <- add(out, "Missing values", "warning",
               sprintf("%d temperature cells are blank and will fall back to the long-term average.", nas))
  } else {
    out <- add(out, "Missing values", "ok", "No blank temperature cells.")
  }

  rng <- range(tm$temp, na.rm = TRUE)
  if (!all(is.finite(rng))) {
    out <- add(out, "Temperature values", "error", "No usable numeric temperatures found.")
  } else if (rng[1] < TEMP_PLAUSIBLE[1] || rng[2] > TEMP_PLAUSIBLE[2]) {
    out <- add(out, "Temperature values", "error",
               sprintf(paste("Temperatures run %.1f to %.1f. Expected roughly %g-%g degrees",
                             "Celsius. Values above about 30 usually mean the file is in",
                             "Fahrenheit."),
                       rng[1], rng[2], TEMP_PLAUSIBLE[1], TEMP_PLAUSIBLE[2]))
  } else {
    out <- add(out, "Temperature values", "ok",
               sprintf("All values between %.1f and %.1f degrees Celsius.", rng[1], rng[2]))
  }

  out
}

#' TRUE if nothing blocks a run
deliverable_ok <- function(report) !any(report$status == "error")


# ---- 3. Build daily series ---------------------------------------------------

#' Build one representative forecast season per scenario
#'
#' Takes the shape of an existing forecast season -- which already blends
#' observed climatology with the modelled window exactly as temperature_data.R
#' built it -- and overwrites the 18 Oct - 31 Dec window with the uploaded
#' scenario. Anything the deliverable does not cover keeps the climatology,
#' which is the same behaviour the pipeline has.
#'
#' @param parsed From read_deliverable().
#' @param env_ext_list The published temperature series.
#' @param base_env Which published alternative supplies the climatological
#'   backbone. Any will do outside the override window; they are identical there.
#' @param season Calendar year of the representative season to borrow.
scenario_series <- function(parsed, env_ext_list, base_env = "1", season = 2050) {
  base <- env_ext_list[[base_env]]
  stopifnot(!is.null(base))

  skeleton <- base %>%
    filter(lubridate::year(Date) == season) %>%
    transmute(site = as.character(site), Date, doy = lubridate::yday(Date), clim = temp)

  ov <- parsed$temps %>%
    filter(!is.na(temp), !is.na(scenario),
           doy >= OVERRIDE_DOY[["start"]], doy <= OVERRIDE_DOY[["end"]]) %>%
    select(met_year, scenario, site, doy, temp_alt = temp)

  grid <- tidyr::expand_grid(
    distinct(ov, met_year, scenario),
    skeleton
  )

  grid %>%
    left_join(ov, by = c("met_year", "scenario", "site", "doy")) %>%
    mutate(temp = dplyr::coalesce(temp_alt, clim),
           temp = pmax(temp, unname(TEMP_FLOOR[site]))) %>%
    arrange(met_year, scenario, site, Date) %>%
    select(met_year, scenario, site, Date, doy, temp)
}


# ---- 4. Spawn-date weights from the CLM --------------------------------------

#' Redd weight by spawn date, driven by the scenario's own Oct/Nov temperatures
#'
#' The pipeline samples redds from this model; here every candidate spawn date
#' carries its CLM probability instead, which removes the sampling noise while
#' keeping the temperature response.
scenario_spawn_weights <- function(series, spawn_model) {
  sm <- spawn_model
  s  <- sm$standardisation

  oct_nov <- series %>%
    filter(lubridate::month(Date) %in% c(10, 11)) %>%
    group_by(met_year, scenario, site, month = lubridate::month(Date)) %>%
    summarise(mt = mean(temp, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = month, values_from = mt, names_prefix = "m") %>%
    rename(Oct = m10, Nov = m11) %>%
    mutate(Oct_std = (Oct - s$oct_mean) / s$oct_sd,
           Nov_std = (Nov - s$nov_mean) / s$nov_sd)

  probs <- predict_clm_probs(sm$beta, sm$zeta,
                             as.data.frame(oct_nov[, c("Oct_std", "Nov_std")]),
                             offset = 0)
  colnames(probs) <- sm$spawn_bins_model
  probs <- probs[, sm$present_bins, drop = FALSE]
  probs <- probs / rowSums(probs)

  # Spread each bin's probability evenly over the days it contains, mapped onto
  # the representative season by month-day.
  bins <- sm$bin_defs %>%
    filter(period %in% sm$present_bins) %>%
    mutate(days = purrr::map2(start, end, ~ format(seq(.x, .y, by = "day"), "%m-%d")))

  site_w <- setNames(sm$site_props$prop, sm$site_props$site)

  # One long table of (bin, day, share-within-bin), joined to the per-series bin
  # probabilities. Building it this way avoids a nested row-wise loop.
  bin_days <- tibble(period = rep(bins$period, lengths(bins$days)),
                     md     = unlist(bins$days, use.names = FALSE)) %>%
    group_by(period) %>% mutate(share = 1 / n()) %>% ungroup()

  probs_long <- as_tibble(probs) %>%
    mutate(.row = dplyr::row_number()) %>%
    tidyr::pivot_longer(-.row, names_to = "period", values_to = "p") %>%
    left_join(oct_nov %>% mutate(.row = dplyr::row_number()) %>%
                select(.row, met_year, scenario, site), by = ".row") %>%
    select(-.row)

  # Deliberately many-to-many: each bin spans several days, and every series
  # carries a probability for every bin.
  probs_long %>%
    inner_join(bin_days, by = "period", relationship = "many-to-many") %>%
    mutate(w = p * share) %>%
    select(met_year, scenario, site, md, w) %>%
    mutate(w = w * unname(site_w[site])) %>%
    group_by(met_year, scenario) %>%
    mutate(w = w / sum(w)) %>%
    ungroup()
}


# ---- 5. Egg-to-fry survival --------------------------------------------------

EXP_PARS <- list(
  exp_WF = list(egg = c(3.408486e-11, 1.21122), alv = c(1.017550e-10, 1.24092)),
  exp_SM = list(egg = c(1.475e-11,    1.392),   alv = c(2.521e-12,    1.461))
)
MARTIN <- c(alpha = 0.026, beta = 12.14)

#' Egg-to-fry survival by scenario and TDM variant, one season
#'
#' Vectorised: cumulative ATU locates emergence with findInterval(), and the
#' hazard over the window is one subtraction. Uses the shipped ATU boundaries
#' (egg_ATU, total_ATU) so it stays in step with functions.R.
scenario_survival <- function(series, weights) {
  # Split once. Filtering inside the loop is quadratic in the number of
  # scenario-site combinations and dominated the runtime.
  chunks <- split(series[order(series$Date), ],
                  list(series$met_year[order(series$Date)],
                       series$scenario[order(series$Date)],
                       series$site[order(series$Date)]),
                  drop = TRUE)

  per <- purrr::map_dfr(chunks, function(x) {
    met_year <- x$met_year[1]; scenario <- x$scenario[1]; site <- x$site[1]
    tt <- x$temp
    n  <- length(tt)
    cum_atu <- c(0, cumsum(pmax(tt, 0)))
    cum_mar <- c(0, cumsum(MARTIN[["alpha"]] * pmax(tt - MARTIN[["beta"]], 0)))
    cum_wfe <- c(0, cumsum(EXP_PARS$exp_WF$egg[1] * exp(EXP_PARS$exp_WF$egg[2] * tt)))
    cum_wfa <- c(0, cumsum(EXP_PARS$exp_WF$alv[1] * exp(EXP_PARS$exp_WF$alv[2] * tt)))
    cum_sme <- c(0, cumsum(EXP_PARS$exp_SM$egg[1] * exp(EXP_PARS$exp_SM$egg[2] * tt)))
    cum_sma <- c(0, cumsum(EXP_PARS$exp_SM$alv[1] * exp(EXP_PARS$exp_SM$alv[2] * tt)))

    i0   <- seq_len(n)
    i_em <- pmax(pmin(findInterval(cum_atu[i0] + total_ATU, cum_atu), n), i0)
    i_h  <- pmax(pmin(findInterval(cum_atu[i0] + egg_ATU,   cum_atu), i_em), i0)

    tibble(met_year = met_year, scenario = scenario, site = site,
           md = format(x$Date, "%m-%d"),
           lin_Martin = exp(-(cum_mar[i_em + 1] - cum_mar[i0])),
           exp_WF = exp(-((cum_wfe[i_h + 1] - cum_wfe[i0]) +
                          (cum_wfa[i_em + 1] - cum_wfa[i_h + 1]))),
           exp_SM = exp(-((cum_sme[i_h + 1] - cum_sme[i0]) +
                          (cum_sma[i_em + 1] - cum_sma[i_h + 1]))))
  })

  per %>%
    tidyr::pivot_longer(c(lin_Martin, exp_WF, exp_SM),
                        names_to = "variant", values_to = "S") %>%
    inner_join(weights, by = c("met_year", "scenario", "site", "md")) %>%
    group_by(met_year, scenario, variant) %>%
    summarise(egg_surv = sum(S * w) / sum(w), .groups = "drop")
}


# ---- 6. Project and score ----------------------------------------------------

#' Project each scenario forward using the published calibration
#'
#' Egg survival is constant across forecast years (the temperature series is
#' periodic), so one value per scenario-variant drives the whole projection.
scenario_project <- function(survival, base_P_list, S_seed_fore_list,
                             sim_years, deg_day = NULL, flow_cfs = NULL,
                             K_fn = NULL) {
  n <- length(sim_years)
  purrr::pmap_dfr(survival, function(met_year, scenario, variant, egg_surv) {
    P <- base_P_list[[variant]][[1]]
    if (!is.null(flow_cfs) && !is.null(K_fn)) P$K_spawners <- K_fn(flow_cfs)
    dd <- if (is.null(deg_day)) rep(0, n) else deg_day
    out <- simulate_variant(
      surv_vec       = rep(egg_surv, n),
      P              = P,
      years          = n,
      S_init         = S_seed_fore_list[[variant]],
      SAR_vec        = rep(P$SAR_mean, n),
      K_spawners_vec = rep(P$K_spawners, n),
      deg_day_adult  = dd,
      sim_years_vec  = sim_years
    )
    tibble(met_year = met_year, scenario = scenario, variant = variant,
           adult_index = median(tail(out$spawners, 20), na.rm = TRUE))
  })
}

#' Steelhead metric: days below 18.3 degC at Hazel Avenue in Oct-Nov
scenario_steelhead <- function(series) {
  series %>%
    filter(site == "AveHazel", lubridate::month(Date) %in% c(10, 11)) %>%
    group_by(met_year, scenario) %>%
    summarise(steelhead_score = sum(temp < 18.3), .groups = "drop")
}

#' Composite MCDA score, matching app.R and analysis/mcda.R
#'
#' @param hydro_loss Named numeric, revenue loss per scenario code.
scenario_scores <- function(projection, steelhead, hydro_loss,
                            tdm_w = c(exp_WF = .51, exp_SM = .24, lin_Martin = .25),
                            met_w = NULL,
                            w_salmon = .40, w_hydro = .50, w_steel = .10) {
  mets <- sort(unique(projection$met_year))
  if (is.null(met_w)) met_w <- setNames(rep(1 / length(mets), length(mets)), mets)

  norm_hi <- function(x) if (diff(range(x)) == 0) rep(0.5, length(x)) else (x - min(x)) / diff(range(x))
  norm_lo <- function(x) if (diff(range(x)) == 0) rep(0.5, length(x)) else (max(x) - x) / diff(range(x))

  chin <- projection %>%
    mutate(w = tdm_w[variant] * met_w[met_year]) %>%
    group_by(scenario) %>%
    summarise(adult_index = sum(adult_index * w), .groups = "drop")

  steel <- steelhead %>%
    mutate(w = met_w[met_year]) %>%
    group_by(scenario) %>%
    summarise(steelhead_score = sum(steelhead_score * w), .groups = "drop")

  chin %>%
    left_join(steel, by = "scenario") %>%
    mutate(hydro_raw = unname(hydro_loss[scenario]),
           chinook_norm   = norm_hi(adult_index),
           steelhead_norm = norm_hi(steelhead_score),
           hydro_norm     = norm_lo(hydro_raw),
           composite = w_salmon * chinook_norm + w_steel * steelhead_norm +
                       w_hydro * hydro_norm) %>%
    arrange(desc(composite)) %>%
    mutate(rank = row_number())
}


# ---- 7. One-call wrapper -----------------------------------------------------

#' Run a deliverable end to end
#'
#' @return list(report, ok, series, survival, projection, steelhead, scores)
run_scenario_deliverable <- function(path, env_ext_list, spawn_model,
                                     base_P_list, S_seed_fore_list, sim_years,
                                     hydro_loss, flow_cfs = NULL, K_fn = NULL,
                                     progress = function(frac, msg) NULL) {
  progress(0.05, "Reading the file")
  parsed <- read_deliverable(path)

  progress(0.15, "Checking the file")
  report <- validate_deliverable(parsed)
  if (!deliverable_ok(report)) {
    return(list(report = report, ok = FALSE))
  }

  progress(0.30, "Building daily temperature series")
  series <- scenario_series(parsed, env_ext_list)

  progress(0.50, "Predicting spawn timing")
  weights <- scenario_spawn_weights(series, spawn_model)

  progress(0.65, "Computing egg-to-fry survival")
  survival <- scenario_survival(series, weights)

  progress(0.80, "Projecting the population")
  projection <- scenario_project(survival, base_P_list, S_seed_fore_list, sim_years,
                                 flow_cfs = flow_cfs, K_fn = K_fn)

  progress(0.92, "Scoring alternatives")
  steel  <- scenario_steelhead(series)
  scores <- scenario_scores(projection, steel, hydro_loss)

  progress(1, "Done")
  list(report = report, ok = TRUE, parsed = parsed, series = series,
       survival = survival, projection = projection, steelhead = steel,
       scores = scores)
}
