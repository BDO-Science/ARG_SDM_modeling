# functions.R
# ═══════════════════════════════════════════════════════════════════════════════
# SALMON POPULATION MODEL - CORE FUNCTIONS LIBRARY
# ═══════════════════════════════════════════════════════════════════════════════
# This file contains all core functions for the salmon population model including:
# - Temperature-Dependent Mortality (TDM) models
# - Life cycle simulation functions  
# - Statistical helper functions
# - Date/time utilities for spawn timing
# ═══════════════════════════════════════════════════════════════════════════════

#' @importFrom stats plogis rnorm rbeta rgamma rlnorm
#' @importFrom lubridate month day mday make_date year hour minute second make_datetime
#' @importFrom dplyr group_by summarize mutate transmute filter arrange desc slice
#' @importFrom tibble tibble
#' @importFrom purrr map imap_dfr map_dfr
#' @importFrom tidyr pivot_wider crossing
#' @importFrom utils modifyList tail
NULL

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1: TEMPERATURE-DEPENDENT MORTALITY (TDM) MODELS
# ═══════════════════════════════════════════════════════════════════════════════
# These functions model egg and alevin mortality as a function of water temperature

#' Estimate egg development time (days) from temperature
#'
#' Computes the number of days required for salmon eggs to reach hatching,
#' assuming a constant requirement of 958 accumulated thermal units (ATU).
#' This is based on the principle that development rate is proportional to
#' temperature above a base threshold (here assumed to be 0°C).
#'
#' @param T Numeric vector of daily mean water temperatures (°C).
#'
#' @return Numeric vector of the same length as `T` with days to hatching.
#' @examples
#' hatch_model(10)             # ~95.8 days at constant 10°C
#' hatch_model(c(8, 12))       # Different times for different temperatures
#' @export
hatch_model <- function(T) {
  # ATU requirement / daily temperature = days to accumulate required ATUs
  958 / T
}


#' Estimate fry emergence time (days) from temperature
#'
#' Computes days until fry emergence from the alevin stage (post-hatch),
#' assuming 417 accumulated thermal units (ATU) are required for this
#' developmental transition.
#'
#' @param T Numeric vector of daily mean water temperatures (°C).
#'
#' @return Numeric vector of the same length as `T` with days to emergence.
#' @examples
#' emergence_model(10)           # ~41.7 days at constant 10°C
#' emergence_model(c(8, 12))     # Faster at warmer temperatures
#' @export
emergence_model <- function(T) {
  # ATU requirement / daily temperature = days to accumulate required ATUs
  417 / T
}


# ───────────────────────────────────────────────────────────────────────────────
# ATU stage boundaries (°C·days) used for egg→hatch→emergence transitions
# Based on empirical salmon development studies
# ───────────────────────────────────────────────────────────────────────────────
egg_ATU   <- 958     # ATUs required from fertilization to hatch
alev_ATU  <- 417     # ATUs required from hatch to emergence
total_ATU <- egg_ATU + alev_ATU  # Total ATUs for complete incubation (1375)


#' Locate hatch and emergence day indices from a temperature series
#'
#' Uses cumulative ATU (accumulated thermal units = cumsum(temps)) to find 
#' the first day at which key developmental milestones are reached:
#'  - Hatch is reached at >= 958 °C·days
#'  - Emergence is reached at >= 1375 °C·days (958 + 417)
#'
#' This function is used internally to split temperature series into
#' egg and alevin stages for stage-specific mortality calculations.
#'
#' @param temps Numeric vector of daily mean temperatures (°C).
#' @param egg_atu Numeric, ATU threshold for hatching (default 958).
#' @param total_atu Numeric, ATU threshold for emergence (default 1375).
#' 
#' @return List with integer elements:
#'   \item{hatch}{Index of hatch day (1-based)}
#'   \item{emerge}{Index of emergence day (1-based)}
#'   Returns last day index if threshold never reached.
#' @keywords internal
.stage_indices_by_atu <- function(temps, egg_atu = egg_ATU, total_atu = total_ATU) {
  # Handle empty or invalid input
  if (length(temps) == 0L || all(!is.finite(temps))) {
    return(list(hatch = 0L, emerge = 0L))
  }
  
  # Calculate cumulative ATUs (negative temperatures contribute 0)
  atu <- cumsum(pmax(temps, 0))
  
  # Find first day each threshold is exceeded
  i_h <- which(atu >= egg_atu)[1]
  i_e <- which(atu >= total_atu)[1]
  
  # If threshold never reached, use last day
  if (is.na(i_h)) i_h <- length(temps)
  if (is.na(i_e)) i_e <- length(temps)
  
  list(hatch = as.integer(i_h), emerge = as.integer(i_e))
}


#' Exponential Temperature-Dependent Mortality model with stage-specific parameters
#'
#' Computes cumulative egg-to-fry survival under an exponential
#' temperature-dependent mortality model. This model uses different
#' mortality parameters for egg vs. alevin stages, with transitions
#' determined by accumulated thermal units (ATU).
#'
#' @details
#' The model implements daily hazard as: h(T) = α * exp(β * T)
#' where α and β differ between egg and alevin stages.
#' 
#' Cumulative survival is: S = exp(-Σ h(T_i))
#' 
#' Two calibration parameter sets are available:
#' \itemize{
#'   \item \strong{WaterForum2020}: Recent calibration from Water Forum studies
#'         - Egg: α=3.408486e-11, β=1.21122
#'         - Alevin: α=1.01755e-10, β=1.24092
#'   \item \strong{SALMOD2006}: Historical SALMOD model parameters
#'         - Egg: α=1.475e-11, β=1.392
#'         - Alevin: α=2.521e-12, β=1.461
#' }
#' 
#' When use_stages = TRUE (default), the model:
#' 1. Calculates cumulative ATUs to find hatch and emergence days
#' 2. Applies egg parameters up to hatch
#' 3. Applies alevin parameters from hatch to emergence
#' 4. Stops accumulating mortality after emergence
#' 
#' Set use_stages = FALSE for backward compatibility (applies egg
#' parameters to all days regardless of developmental stage).
#'
#' @param temps Numeric vector of daily incubation temperatures (°C).
#' @param calib Character, calibration set: "WaterForum2020" or "SALMOD2006".
#' @param use_stages Logical, if TRUE (default) use stage-specific parameters.
#'
#' @return Scalar survival probability in [0, 1].
#' @examples
#' # Constant 11°C for 70 days with Water Forum 2020 calibration:
#' tdm_exp(rep(11, 70), "WaterForum2020")
#'
#' # SALMOD 2006 without stage splitting (legacy mode):
#' tdm_exp(rep(12, 60), "SALMOD2006", use_stages = FALSE)
#' @export
tdm_exp <- function(temps,
                    calib = c("WaterForum2020", "SALMOD2006"),
                    use_stages = TRUE) {
  # Validate calibration choice
  calib <- match.arg(calib)
  
  # Handle empty input
  if (length(temps) == 0L) return(NA_real_)
  
  # Define parameter sets for each calibration and life stage
  pars <- switch(
    calib,
    "WaterForum2020" = list(
      egg    = list(alpha = 3.408488e-11, beta = 1.21122),
      alevin = list(alpha = 1.017554e-10, beta = 1.24092)
    ),
    "SALMOD2006" = list(
      egg    = list(alpha = 1.475e-11, beta = 1.392),
      alevin = list(alpha = 2.521e-12, beta = 1.461)
    )
  )
  
  # Hazard function: daily mortality rate
  haz <- function(T, a, b) a * exp(b * T)
  
  # Legacy mode: apply egg parameters to entire period
  if (!use_stages) {
    H <- sum(haz(temps, pars$egg$alpha, pars$egg$beta), na.rm = TRUE)
    return(exp(-H))
  }
  
  # Stage-aware mode: split incubation period by developmental stage
  idx <- .stage_indices_by_atu(temps)
  i_h <- max(0L, min(idx$hatch,  length(temps)))
  i_e <- max(i_h, min(idx$emerge, length(temps)))
  
  # Extract temperature slices for each stage
  egg_slice    <- if (i_h > 0L)     temps[1:i_h]        else numeric(0)
  alevin_slice <- if (i_e > i_h)    temps[(i_h+1):i_e]  else numeric(0)
  
  # Calculate stage-specific cumulative hazards
  H_egg <- if (length(egg_slice))
    sum(haz(egg_slice, pars$egg$alpha, pars$egg$beta), na.rm = TRUE) else 0
  H_alv <- if (length(alevin_slice))
    sum(haz(alevin_slice, pars$alevin$alpha, pars$alevin$beta), na.rm = TRUE) else 0
  
  # Return cumulative survival
  exp(-(H_egg + H_alv))
}


#' Linear threshold Temperature-Dependent Mortality model
#'
#' Implements the Martin et al. (2017) linear threshold model where
#' mortality only occurs above a threshold temperature. This model
#' assumes zero mortality below the threshold and linear increase above it.
#'
#' @details
#' Cumulative survival formula:
#' S = exp(-α * Σ max(T_i - β, 0))
#' 
#' where:
#' - T_i is daily temperature
#' - β is the threshold temperature (no mortality below this)
#' - α is the mortality coefficient (rate of mortality increase above threshold)
#' 
#' Default parameters (α=0.026, β=12.14) are from Martin et al. (2017)
#' calibrated for Sacramento River Chinook salmon.
#'
#' @param temps Numeric vector of daily incubation temperatures (°C).
#' @param α Numeric, mortality coefficient (default 0.026).
#' @param β Numeric, threshold temperature in °C (default 12.14).
#'
#' @return Scalar survival probability in [0, 1].
#' @references Martin, B.T. et al. (2017). Phenomenological vs. biophysical 
#'   models of thermal stress in aquatic eggs. *Ecology Letters* 20:50–59.
#' @examples
#' tdm_lin_martin(c(10, 11, 12))  # Only 12°C contributes to mortality
#' tdm_lin_martin(c(10, 11, 12), α = 0.03, β = 11)  # Custom parameters
#' @export
tdm_lin_martin <- function(temps, α = 0.026, β = 12.14) {
  # Calculate excess temperature above threshold
  # Sum only positive exceedances (temperatures below threshold contribute 0)
  exp(-α * sum(pmax(temps - β, 0)))
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2: ATU-BASED SURVIVAL CALCULATIONS
# ═══════════════════════════════════════════════════════════════════════════════

#' Helper function to slice temperature series by ATU thresholds
#' 
#' Marches forward through temperatures until total ATU reaches emergence.
#' Used internally to determine the exact incubation window.
#' 
#' @param temps Numeric vector of daily temperatures
#' @param egg_atu ATU threshold for hatch (default 958)
#' @param total_atu ATU threshold for emergence (default 1375)
#' @return Integer vector of indices to include in incubation period
#' @keywords internal
.slice_by_atu <- function(temps, egg_atu = egg_ATU, total_atu = total_ATU) {
  if (!length(temps)) return(integer())
  
  # Calculate cumulative ATUs
  atu <- cumsum(pmax(temps, 0))
  
  # Find emergence day
  i_end <- which(atu >= total_atu)[1]
  if (is.na(i_end)) i_end <- length(temps)
  
  # Return indices from start to emergence
  seq_len(i_end)
}


#' Compute egg-to-fry survival using ATU-precise incubation window
#'
#' This function provides the most accurate survival calculation by:
#' 1. Starting from the actual redd (spawning) date
#' 2. Marching forward through realized daily temperatures
#' 3. Accumulating ATUs until reaching hatch (958) and emergence (1375)
#' 4. Applying the TDM model to exactly that temperature window
#' 
#' This approach accounts for variable development rates due to
#' temperature fluctuations, unlike fixed-window approximations.
#'
#' @param rdr Date or IDate, redd (spawning) date.
#' @param site Character, site identifier.
#' @param date_idx_env Named list mapping dates to array indices for each site.
#' @param temps_env Named list of temperature vectors for each site.
#' @param model Character, "exp" for exponential or "lin" for linear TDM.
#' @param calib Character, calibration for exponential model (ignored if linear).
#' @param max_days Integer, maximum look-ahead window in days (default 300).
#'
#' @return Numeric survival in [0,1]; NA_real_ if data unavailable.
#' @examples
#' # compute_surv_by_atu(as.Date("2020-11-15"), "AveHazel", 
#' #                     date_idx_env, temps_env, "exp", "WaterForum2020")
#' @export
compute_surv_by_atu <- function(rdr, site, date_idx_env, temps_env, 
                                model, calib, max_days = 300L) {
  # Look up position of redd date in temperature array
  pos <- date_idx_env[[site]][as.character(rdr)]
  if (is.na(pos) || pos < 1) return(NA_real_)
  
  # Extract temperature slice (up to max_days for safety)
  end_idx <- min(length(temps_env[[site]]), pos + max_days)
  Tvec <- temps_env[[site]][pos:end_idx]
  
  # Handle missing data
  if (!length(Tvec) || all(!is.finite(Tvec))) return(NA_real_)
  
  # Determine incubation period based on ATU accumulation
  idx <- .slice_by_atu(Tvec, egg_ATU, total_ATU)
  slice <- Tvec[idx]
  
  # Apply appropriate TDM model
  if (model == "exp") {
    tdm_exp(slice, calib)
  } else {
    tdm_lin_martin(slice)
  }
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3: ADULT PRE-SPAWN SURVIVAL
# ═══════════════════════════════════════════════════════════════════════════════

#' Adult pre-spawn survival from degree-days (logistic model)
#'
#' Models the probability that adult salmon survive from river entry
#' to spawning as a logistic function of accumulated degree-days.
#' Higher degree-days (warmer water × longer holding time) reduce survival.
#'
#' @details
#' Survival probability: S = 1 / (1 + exp(-(intercept + beta * deg_day)))
#' 
#' Default parameters (intercept=3.0, beta=-0.00067) are from 
#' Colvin et al. (2018) for Columbia River Chinook salmon.
#' 
#' The negative beta coefficient means survival decreases with
#' increasing degree-days, reflecting thermal stress during holding.
#'
#' @param deg_day Numeric vector of cumulative pre-spawn degree-days (°C·days).
#' @param intercept Numeric, logit intercept (default 3.0).
#' @param beta Numeric, logit slope per °C·day (default -0.00067).
#'
#' @return Numeric vector of survival probabilities in (0, 1).
#' @references Colvin et al. (2018). Identifying optimal water temperature 
#'   and flow regimes for anadromous fish. *River Research and Applications* 
#'   34(6):621–632.
#' @examples
#' surv_adult_prespawn(0)         # ~95% survival at 0 degree-days
#' surv_adult_prespawn(seq(0, 2000, 500))  # Declining survival curve
#' @export
surv_adult_prespawn <- function(deg_day,
                                intercept = 3.0,
                                beta      = -0.00067) {
  # Logistic function: plogis(x) == 1/(1+exp(-x))
  plogis(intercept + beta * deg_day)
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4: DATE AND SEASON UTILITIES
# ═══════════════════════════════════════════════════════════════════════════════

#' Convert dates to season-anchored timestamps
#'
#' Transforms dates to a common seasonal reference frame, useful for
#' aligning multi-year data to a single seasonal cycle (e.g., water year
#' or biological year starting October 1).
#'
#' @details
#' This function determines which "season year" each date belongs to based
#' on an anchor date (default Oct 1). Dates before the anchor in a calendar
#' year belong to the previous season year. The output preserves month/day/time
#' but replaces the year with the season year.
#' 
#' Example: With Oct 1 anchor:
#' - 2020-09-15 → season year 2019 → output 2019-09-15
#' - 2020-10-15 → season year 2020 → output 2020-10-15
#'
#' @param x Date or POSIXct vector to transform.
#' @param anchor_mmdd Anchor date as "MM-DD", c(MM,DD), or Date (default Oct 1).
#' @param tz Character, output timezone for POSIXct results (default "UTC").
#'
#' @return POSIXct vector with dates adjusted to season years.
#' @examples
#' dates <- as.Date(c("2020-09-15", "2020-10-15"))
#' season_posix(dates)  # Different season years
#' @export
season_posix <- function(x, anchor_mmdd = c(10, 1), tz = "UTC") {
  # Input validation
  stopifnot(inherits(x, "Date") || inherits(x, "POSIXct") || inherits(x, "POSIXt"))
  
  # Parse anchor date from various formats
  if (is.character(anchor_mmdd) && length(anchor_mmdd) == 1L) {
    mmdd <- strsplit(anchor_mmdd, "[-/\\.]")[[1]]
    mm <- as.integer(mmdd[1])
    dd <- as.integer(mmdd[2])
  } else if (inherits(anchor_mmdd, "Date")) {
    mm <- lubridate::month(anchor_mmdd)
    dd <- lubridate::day(anchor_mmdd)
  } else if (is.numeric(anchor_mmdd) && length(anchor_mmdd) >= 2L) {
    mm <- as.integer(anchor_mmdd[1])
    dd <- as.integer(anchor_mmdd[2])
  } else {
    stop("anchor_mmdd must be 'MM-DD', c(MM,DD), or a Date.")
  }
  
  if (anyNA(c(mm, dd))) stop("Could not parse anchor_mmdd.")
  
  # Convert to POSIXct for uniform processing
  is_date   <- inherits(x, "Date")
  is_posix  <- inherits(x, "POSIXct") || inherits(x, "POSIXt")
  x_posix   <- if (is_date) as.POSIXct(x) else as.POSIXct(x)
  yrs       <- lubridate::year(x_posix)
  
  # Determine season year for each date
  anchor_this_year <- lubridate::make_datetime(year = yrs, month = mm, day = dd, tz = tz)
  season_year <- ifelse(x_posix >= anchor_this_year, yrs, yrs - 1L)
  
  # Rebuild dates with season year
  out <- lubridate::make_datetime(
    year  = season_year,
    month = lubridate::month(x_posix),
    day   = lubridate::mday(x_posix),
    hour  = lubridate::hour(x_posix),
    min   = lubridate::minute(x_posix),
    sec   = lubridate::second(x_posix),
    tz    = tz
  )
  
  # Return as POSIXct
  out
}


#' Extract season year as integer
#'
#' Helper function that returns just the season year number.
#'
#' @param x Date or POSIXct vector.
#' @param anchor_mmdd Anchor date (see season_posix).
#' @return Integer vector of season years.
#' @export
season_year <- function(x, anchor_mmdd = c(10, 1)) {
  as.integer(lubridate::year(season_posix(x, anchor_mmdd)))
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5: SPAWN DATE VECTOR CONSTRUCTION
# ═══════════════════════════════════════════════════════════════════════════════

#' Build environment-specific spawn date vector aligned to simulation years
#'
#' Creates a complete spawn date vector for all simulation years, filling
#' gaps using Last Observation Carried Forward (LOCF) and backfill strategies.
#' This ensures every year has a spawn date for degree-day calculations.
#'
#' @details
#' Gap-filling strategy:
#' 1. Match observed dates to simulation years
#' 2. Forward-fill gaps (LOCF)
#' 3. Backward-fill remaining gaps from first known date
#' 4. Use fallback date for any remaining NAs
#' 
#' The function handles various input formats (Date, POSIXct, numeric, character)
#' and ensures output dates have the correct year for degree-day accumulation.
#'
#' @param df_env Data frame with columns:
#'   \item{sim_year}{Integer year}
#'   \item{spawn_dt}{Spawn date (Date/POSIXct/numeric/character)}
#' @param sim_years Integer vector of all simulation years to generate.
#' @param fallback_md Integer c(month, day) for missing data (default Nov 20).
#'
#' @return Date vector aligned to sim_years with no gaps.
#' @export
build_spawn_vec_for_env <- function(df_env, sim_years, fallback_md = c(11, 20)) {
  # Validate required columns
  stopifnot(all(c("sim_year","spawn_dt") %in% names(df_env)))
  
  # Coerce spawn_dt to Date regardless of input type
  v_raw <- df_env$spawn_dt
  v <- if (inherits(v_raw, "Date")) {
    v_raw
  } else if (inherits(v_raw, c("POSIXct","POSIXt"))) {
    as.Date(v_raw)
  } else if (is.numeric(v_raw)) {
    as.Date(v_raw, origin = "1970-01-01")
  } else {
    # Character or other format - attempt conversion
    suppressWarnings(as.Date(v_raw))
  }
  
  df_env2 <- df_env
  df_env2$spawn_dt <- v
  
  # Align to simulation years
  m <- match(sim_years, df_env2$sim_year)
  v <- df_env2$spawn_dt[m]  # May contain NAs
  
  # If completely empty, use fallback for all years
  if (length(v) == 0L || all(is.na(v))) {
    return(lubridate::make_date(
      year = sim_years,
      month = fallback_md[1],
      day = fallback_md[2]
    ))
  }
  
  # Last Observation Carried Forward (LOCF)
  for (i in seq_along(v)) {
    if (is.na(v[i]) && i > 1) v[i] <- v[i - 1]
  }
  
  # Backfill from first known date
  if (anyNA(v)) {
    first_ok <- which(!is.na(v))[1]
    if (length(first_ok) == 1L && !is.na(first_ok)) {
      v[seq_len(first_ok - 1)] <- v[first_ok]
    }
  }
  
  # Apply fallback to any remaining NAs
  if (anyNA(v)) {
    v[is.na(v)] <- lubridate::make_date(
      year  = sim_years[is.na(v)],
      month = fallback_md[1],
      day   = fallback_md[2]
    )
  }
  
  # Rebuild dates with correct years and month/day from filled vector
  lubridate::make_date(
    year  = sim_years,
    month = lubridate::month(v),
    day   = lubridate::day(v)
  )
}


#' Calculate cumulative adult pre-spawn degree-days for each brood year
#'
#' Computes accumulated thermal units experienced by adult salmon from
#' river entry (fixed date) to spawning (variable date). Higher accumulation
#' indicates more thermal stress and reduced pre-spawn survival.
#'
#' @details
#' For each brood year:
#' 1. Start accumulation on start_month/start_day
#' 2. Sum daily temperatures above base_temp through spawn date
#' 3. Return total degree-days
#' 
#' Formula: Σ max(T - base_temp, 0) for all days in holding period
#' 
#' If spawn_dates length differs from sim_years, the vector is
#' padded or truncated to match.
#'
#' @param env_nm Character, environment name (key in env_ext_list).
#' @param sim_years Integer vector of brood years to calculate.
#' @param spawn_dates Date vector of spawn dates aligned to sim_years.
#' @param env_ext_list Named list of data.frames with Date and temp columns.
#' @param start_month Integer, month of river entry (default 10 = October).
#' @param start_day Integer, day of river entry (default 1).
#' @param base_temp Numeric, temperature threshold in °C (default 0).
#'
#' @return Numeric vector of degree-days aligned to sim_years.
#' @examples
#' # compute_deg_day_adult("Alt1", 2011:2015, spawn_dates, env_ext_list)
#' @export
compute_deg_day_adult <- function(
    env_nm,
    sim_years,
    spawn_dates,
    env_ext_list,
    start_month = 10,
    start_day   = 1,
    base_temp   = 0
) {
  # Validate environment name
  if (!is.character(env_nm) || length(env_nm) != 1 ||
      !env_nm %in% names(env_ext_list)) {
    stop("`env_nm` must be a single name in names(env_ext_list).")
  }
  
  # Align spawn_dates length to sim_years
  if (length(spawn_dates) < length(sim_years)) {
    # Pad with last value
    spawn_dates <- c(
      spawn_dates,
      rep(tail(spawn_dates, 1), length(sim_years) - length(spawn_dates))
    )
  } else if (length(spawn_dates) > length(sim_years)) {
    # Truncate to match
    spawn_dates <- spawn_dates[seq_len(length(sim_years))]
  }
  
  # Build daily mean temperature series
  daily_df <- env_ext_list[[env_nm]] %>%
    dplyr::group_by(Date) %>%
    dplyr::summarize(Tmean = mean(temp, na.rm = TRUE), .groups = "drop")
  
  # Calculate degree-days for each brood year
  vapply(seq_along(sim_years), function(i) {
    # Define holding period
    start_date <- as.Date(sprintf(
      "%d-%02d-%02d",
      sim_years[i], start_month, start_day
    ))
    end_date <- spawn_dates[i]
    
    # Extract temperatures for holding period
    temps <- daily_df$Tmean[
      daily_df$Date >= start_date &
        daily_df$Date <= end_date
    ]
    
    # Return accumulated degree-days
    if (length(temps) == 0) return(NA_real_)
    sum(pmax(temps - base_temp, 0), na.rm = TRUE)
  }, numeric(1))
}


#' Helper: Calculate degree-days for calibration period
#'
#' Convenience function for computing degree-days during the
#' historical calibration period for a given environment.
#'
#' @param env_nm Character, environment name.
#' @return Numeric vector of degree-days for real_years
#' @keywords internal
deg_day_cal_for <- function(env_nm) {
  # Get spawn dates for this environment
  sd_env <- spawn_dates_by_env[[env_nm]]
  stopifnot(length(sd_env) >= length(sim_years))
  
  # Align to calibration years
  sd_cal <- sd_env[match(real_years, sim_years)]
  
  # Compute degree-days
  compute_deg_day_adult(
    env_nm       = env_nm,
    sim_years    = real_years,
    spawn_dates  = sd_cal,
    env_ext_list = env_ext_list
  )
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 6: CUMULATIVE LINK MODEL (CLM) UTILITIES
# ═══════════════════════════════════════════════════════════════════════════════

#' Predict probabilities from Cumulative Link Model
#'
#' Computes category probabilities from a fitted CLM given coefficients,
#' thresholds, and new predictor data. Optionally adds random offset for
#' stochastic predictions.
#'
#' @details
#' The CLM models ordered categories using cumulative logits:
#' - Cumulative prob: P(Y ≤ j) = logit^(-1)(ζ_j - x'β)
#' - Category prob: P(Y = j) = P(Y ≤ j) - P(Y ≤ j-1)
#' 
#' Where ζ_j are threshold parameters and β are covariate effects.
#' Adding an offset introduces individual-level variation.
#'
#' @param beta Named numeric vector of regression coefficients.
#' @param zeta Numeric vector of threshold parameters (K-1 for K categories).
#' @param newdata Data frame with columns matching names(beta).
#' @param offset Numeric scalar or vector, added to linear predictor (default 0).
#'
#' @return Matrix with rows = nrow(newdata), columns = K categories.
#'   Each row sums to 1 (probability distribution over categories).
#' @keywords internal
predict_clm_probs <- function(beta, zeta, newdata, offset = 0) {
  # Extract predictor matrix
  X  <- as.matrix(newdata[, names(beta), drop = FALSE])
  
  # Calculate linear predictor with optional offset
  xb <- as.vector(X %*% beta) + as.numeric(offset)
  
  # Convert to probabilities for each observation
  t(vapply(xb, function(xi) {
    # Cumulative probabilities
    cp <- plogis(zeta - xi)
    
    # Category probabilities (differences)
    p  <- c(cp, 1) - c(0, cp)
    
    # Normalize to ensure sum = 1
    p / sum(p)
  }, numeric(length(zeta) + 1)))
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 7: TEMPERATURE DATA PROCESSING
# ═══════════════════════════════════════════════════════════════════════════════

#' Build forecast temperature data for future years
#'
#' Extracts October and November mean temperatures for forecast years,
#' carrying forward the last available year's data when future projections
#' are incomplete. Standardizes using calibration period statistics.
#'
#' @details
#' Process:
#' 1. Extract Oct/Nov temperatures from each environment
#' 2. Calculate monthly means by year
#' 3. Fill missing forecast years with last available data
#' 4. Standardize using calibration mean/SD
#'
#' @param env_ext_list Named list of temperature data frames.
#' @param yrs_forecast Integer vector of forecast years.
#' @param sc Data frame with standardization parameters (o_m, o_s, n_m, n_s).
#'
#' @return Data frame with columns:
#'   env, sim_year, Oct, Nov, Oct_std, Nov_std
#' @keywords internal
build_forecast_temps <- function(env_ext_list, yrs_forecast, sc) {
  # Extract Oct/Nov temperatures for all environments
  raw <- purrr::imap_dfr(env_ext_list, function(df_env, env_nm) {
    df_env %>%
      dplyr::filter(lubridate::month(Date) %in% c(10, 11)) %>%
      dplyr::mutate(
        env   = as.character(env_nm),
        year  = lubridate::year(Date),
        month = lubridate::month(Date)
      ) %>%
      dplyr::group_by(env, year, month) %>%
      dplyr::summarise(mean_temp = mean(temp, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = month, values_from = mean_temp, names_prefix = "m_") %>%
      dplyr::rename(Oct = m_10, Nov = m_11)
  })
  
  # Process each environment separately
  purrr::map_dfr(split(raw, raw$env), function(d) {
    if (!nrow(d)) {
      return(tibble(env = unique(d$env), sim_year = yrs_forecast, 
                    Oct = NA_real_, Nov = NA_real_))
    }
    
    # Identify available and missing years
    have <- d %>% 
      dplyr::rename(sim_year = year) %>% 
      dplyr::select(env, sim_year, Oct, Nov)
    
    missing_years <- setdiff(yrs_forecast, have$sim_year)
    
    # Carry forward last year's data for missing years
    cf_rows <- if (length(missing_years)) {
      last_row <- d %>% 
        dplyr::arrange(dplyr::desc(year)) %>% 
        dplyr::slice(1) %>% 
        dplyr::transmute(env = unique(env), Oct, Nov)
      tidyr::crossing(last_row, sim_year = missing_years)
    } else {
      tibble()
    }
    
    dplyr::bind_rows(have, cf_rows) %>%
      dplyr::arrange(sim_year)
  }) %>%
    dplyr::mutate(
      # Standardize using calibration statistics
      Oct_std = (Oct - sc$o_m) / sc$o_s,
      Nov_std = (Nov - sc$n_m) / sc$n_s
    )
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 8: MODEL CALIBRATION
# ═══════════════════════════════════════════════════════════════════════════════

#' Sum of Squared Errors objective for life-cycle calibration
#'
#' Computes SSE between predicted and observed spawner abundances for
#' calibration years. Used to optimize SAR_mean and rear_surv parameters.
#'
#' @details
#' This function:
#' 1. Sets SAR and rearing survival from parameter vector
#' 2. Runs life-cycle simulation for calibration period
#' 3. Calculates SSE for fit years (typically years 4-14)
#' 
#' The optimization finds parameters that minimize discrepancy between
#' model predictions and observed escapement.
#'
#' @param par Numeric length-2: c(SAR_mean, rear_surv).
#' @param variant Character, TDM variant name (e.g., "exp_WF").
#'
#' @return Numeric scalar SSE; returns large value if simulation fails.
#' @keywords internal
modular_sse <- function(par, variant) {
  # Update parameters
  P_tmp <- base_P
  P_tmp$SAR_mean  <- par[1]
  P_tmp$rear_surv <- par[2]
  
  # Run simulation with trial parameters
  out <- simulate_variant(
    surv_vec       = surv_lookup_by_variant[[variant]][1:n_calib],
    P              = P_tmp,
    years          = n_calib,
    S_init         = S_seed_calib,
    SAR_vec        = rep(P_tmp$SAR_mean,  n_calib),
    K_spawners_vec = rep(P_tmp$K_spawners, n_calib),
    deg_day_adult  = deg_day_cal_ref,  # Use actual degree-days
    sim_years_vec  = real_years
  )
  
  # Calculate SSE for fit years only
  preds <- out$spawners[fit_idx]
  if (!all(is.finite(preds))) return(.Machine$double.xmax)
  
  sum((preds - obs_spawners[fit_idx])^2)
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 9: POPULATION FORECASTING
# ═══════════════════════════════════════════════════════════════════════════════

#' Build a forecast simulator function for one variant-environment pair
#'
#' Creates a closure (function factory) that encapsulates all parameters
#' and data needed to simulate future population dynamics for a specific
#' TDM variant and environmental scenario.
#'
#' @details
#' The returned function:
#' 1. Uses environment-specific temperatures and spawn timing
#' 2. Calculates pre-spawn degree-days for adult survival
#' 3. Applies calibrated life-cycle parameters
#' 4. Optionally includes stochastic SAR variation
#' 
#' This design allows parallel execution across variant-environment
#' combinations without parameter conflicts.
#'
#' @param var_nm Character, TDM variant name (e.g., "exp_WF").
#' @param env_nm Character, environment/alternative name.
#' @param flow_cfs Numeric or NULL. If provided, adjusts carrying capacity.
#' @param S_seed Numeric vector, initial spawner abundances.
#' @param spawn_dates_by_env Named list, env -> Date vector aligned to sim_years.
#'
#' @return Zero-argument function that returns a tibble with columns:
#'   year, spawners, deg_day, pre_spawn, dd, fry_dd, egg_surv, eff_surv,
#'   rear_surv, SAR_used, K_spawners, env, variant
#' @export
sim_forecast_fn <- function(var_nm,
                            env_nm,
                            flow_cfs = NULL,
                            S_seed,
                            spawn_dates_by_env) {
  # Force evaluation of arguments in parent scope
  force(var_nm); force(env_nm); force(flow_cfs); force(S_seed); force(spawn_dates_by_env)
  
  # Return closure that captures these values
  function() {
    # Get parameters for this variant-environment
    P_tmp <- base_P_list[[var_nm]][[env_nm]]
    if (!is.null(flow_cfs)) P_tmp$K_spawners <- get_K_spawners(flow_cfs)
    
    # Get survival vector for this combination
    surv_vec <- surv_lookup_full[[paste(env_nm, var_nm, sep = "_")]]
    
    # Get environment-specific spawn dates
    if (is.null(spawn_dates_by_env[[env_nm]])) {
      stop("spawn_dates_by_env[[", env_nm, "]] is NULL.")
    }
    spawn_dates_env <- spawn_dates_by_env[[env_nm]]
    if (length(spawn_dates_env) < length(sim_years)) {
      stop("spawn_dates_by_env[[", env_nm, "]] length < sim_years length.")
    }
    
    # Calculate pre-spawn degree-days
    deg_day_vec <- compute_deg_day_adult(
      env_nm       = env_nm,
      sim_years    = sim_years,
      spawn_dates  = spawn_dates_env,
      env_ext_list = env_ext_list
    )
    
    # Set up carrying capacity and SAR vectors
    K_vec   <- rep(P_tmp$K_spawners, length(sim_years))
    SAR_vec <- if (use_stochastic_SAR) {
      generate_SAR_vec(length(sim_years), 
                       modifyList(stoch_SAR_opts, list(mean = P_tmp$SAR_mean)))
    } else {
      rep(P_tmp$SAR_mean, length(sim_years))
    }
    
    # Run simulation
    sim_out <- simulate_variant(
      surv_vec       = surv_vec,
      P              = P_tmp,
      years          = length(sim_years),
      S_init         = S_seed,
      SAR_vec        = SAR_vec,
      K_spawners_vec = K_vec,
      deg_day_adult  = deg_day_vec,
      sim_years_vec  = sim_years
    )
    
    # Add environment and variant labels
    dplyr::mutate(sim_out, env = env_nm, variant = var_nm)
  }
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 10: LIFE CYCLE SIMULATION
# ═══════════════════════════════════════════════════════════════════════════════

#' Simulate full salmon life cycle with age structure
#'
#' Core population dynamics model that tracks cohorts through:
#' egg → fry → smolt → ocean → adult returns at ages 3-5.
#' Includes density dependence, stage-specific mortality, and
#' temperature effects on adults and eggs.
#'
#' @details
#' Life cycle stages and processes:
#' 1. Adult pre-spawn survival (temperature-dependent via degree-days)
#' 2. Spawning (females × fecundity)
#' 3. Egg-to-fry survival (TDM from surv_vec)
#' 4. Density dependence (Beverton-Holt on fry)
#' 5. Freshwater rearing survival
#' 6. Ocean survival (SAR)
#' 7. Age-structured returns (3-5 years later)
#' 
#' The model tracks brood years and properly accounts for the
#' multi-year lag between spawning and adult returns.
#'
#' @param surv_vec Numeric vector, egg-to-fry survival by year (from TDM).
#' @param P List of biological parameters:
#'   \item{female_fraction}{Proportion of spawners that are female}
#'   \item{fec}{Eggs per female}
#'   \item{S0}{Maximum fry survival at low density}
#'   \item{rear_surv}{Freshwater rearing survival}
#'   \item{K_spawners}{Spawner carrying capacity}
#' @param years Integer, number of years to simulate.
#' @param S_init Numeric vector, initial spawner abundances.
#' @param SAR_vec Numeric vector, smolt-to-adult return rates by brood year.
#' @param K_spawners_vec Numeric vector, carrying capacity by year.
#' @param deg_day_adult Numeric vector, pre-spawn degree-days by year.
#' @param sim_years_vec Integer vector of calendar years.
#' @param pre_int Numeric, logistic intercept for pre-spawn survival (default 3.0).
#' @param pre_beta Numeric, logistic slope for pre-spawn survival (default -0.00067).
#'
#' @return Tibble with annual metrics:
#'   \item{year}{Calendar year}
#'   \item{spawners}{Adult spawner abundance}
#'   \item{deg_day}{Pre-spawn degree-days}
#'   \item{pre_spawn}{Pre-spawn survival rate}
#'   \item{dd}{Density-dependent survival multiplier}
#'   \item{fry_dd}{Fry production after density dependence}
#'   \item{egg_surv}{Egg-to-fry survival (TDM)}
#'   \item{eff_surv}{Effective survival (egg × density)}
#'   \item{rear_surv}{Freshwater rearing survival}
#'   \item{SAR_used}{Smolt-to-adult return rate}
#'   \item{K_spawners}{Carrying capacity}
#' @export
simulate_variant <- function(
    surv_vec, P, years, S_init, SAR_vec, K_spawners_vec,
    deg_day_adult = NULL,
    sim_years_vec,
    pre_int = 3.0, pre_beta = -0.00067
) {
  # Extend inputs to full simulation length
  surv_vec       <- rep_len(surv_vec, years)
  K_spawners_vec <- rep_len(K_spawners_vec, years)
  SAR_vec        <- rep_len(SAR_vec, years)
  seed_len       <- min(length(S_init), years)
  
  # Initialize storage arrays
  S                 <- numeric(years)  # Spawner abundance
  S_pre             <- numeric(years)  # Pre-spawn survival
  deg_day_adult_vec <- if (is.null(deg_day_adult)) rep(0, years) else rep_len(deg_day_adult, years)
  dd_vec            <- numeric(years)  # Density dependence
  fry_dd            <- numeric(years)  # Fry after density dependence
  rear_surv_vec     <- numeric(years)  # Rearing survival
  reared_vec        <- numeric(years)  # Smolts produced
  
  # Seed initial years with observed/specified abundances
  if (seed_len > 0) S[1:seed_len] <- S_init[1:seed_len]
  
  # Main simulation loop
  for (t in seq_len(years)) {
    # Use seeded value if available
    if (t <= seed_len && !is.na(S_init[t])) S[t] <- S_init[t]
    
    # Pre-spawn survival (temperature effect on holding adults)
    S_pre[t] <- surv_adult_prespawn(deg_day_adult_vec[t], 
                                    intercept = pre_int, beta = pre_beta)
    
    # Calculate egg production
    redds <- S[t] * P$female_fraction * S_pre[t]
    eggs  <- redds * P$fec
    
    # Apply density dependence (Beverton-Holt)
    dd_vec[t]        <- P$S0 / (1 + redds / K_spawners_vec[t])
    
    # Calculate fry production
    fry_dd[t]        <- eggs * surv_vec[t] * dd_vec[t]
    
    # Freshwater rearing
    rear_surv_vec[t] <- P$rear_surv
    reared_vec[t]    <- fry_dd[t] * P$rear_surv
    
    # Distribute returns across future years (age structure)
    for (age in 3:5) {
      return_year <- t + age
      if (return_year <= years) {
        # Add returns from this brood year at specified age
        prob_age <- P$lag_probs[as.character(age)]
        S[return_year] <- S[return_year] + 
          reared_vec[t] * SAR_vec[t] * prob_age
      }
    }
  }
  
  # Return results as tibble
  tibble::tibble(
    year       = sim_years_vec,
    spawners   = S,
    deg_day    = deg_day_adult_vec,
    pre_spawn  = S_pre,
    dd         = dd_vec,
    fry_dd     = fry_dd,
    egg_surv   = surv_vec,
    eff_surv   = surv_vec * dd_vec,
    rear_surv  = rear_surv_vec,
    SAR_used   = SAR_vec,
    K_spawners = K_spawners_vec
  )
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 11: STOCHASTIC SAR GENERATION
# ═══════════════════════════════════════════════════════════════════════════════

#' Generate stochastic Smolt-to-Adult Return (SAR) time series
#'
#' Creates a time series of SAR values with optional temporal patterns.
#' Supports multiple probability distributions and timing options for
#' incorporating environmental variability into population projections.
#'
#' @details
#' Distribution options:
#' - Normal: Standard normal with specified mean/SD
#' - Lognormal: Log-normal to ensure positive values
#' - Beta: Bounded between 0 and 1
#' - Gamma: Positive values with flexible shape
#' 
#' Timing options:
#' - "all": Stochastic values for all years
#' - "block": Stochastic only for specified year ranges
#' - "pulse": Baseline with occasional perturbations
#' 
#' All values are truncated at 0 (no negative SARs allowed).
#'
#' @param n_years Integer, length of time series to generate.
#' @param opts List of options:
#'   \describe{
#'     \item{model}{"normal", "lognormal", "beta", or "gamma"}
#'     \item{mean}{Mean of the SAR distribution}
#'     \item{sd}{Standard deviation (Normal/Lognormal)}
#'     \item{shape1, shape2}{Shape parameters (Beta/Gamma)}
#'     \item{timing}{"all", "block", or "pulse"}
#'     \item{block_years}{Years to apply stochasticity (if timing="block")}
#'     \item{pulse_years}{Years for pulses (if timing="pulse")}
#'     \item{pulse_sd}{SD for pulse years (if timing="pulse")}
#'   }
#'
#' @return Non-negative numeric vector of length n_years.
#' @examples
#' # Constant SAR with occasional bad years
#' generate_SAR_vec(30, list(model="normal", mean=0.003, sd=0.001, 
#'                          timing="pulse", pulse_years=c(5,15,25), 
#'                          pulse_sd=0.002))
#' @export
generate_SAR_vec <- function(n_years, opts) {
  # Validate input
  n_years <- as.integer(n_years)[1]
  if (is.na(n_years) || n_years < 1) {
    stop("generate_SAR_vec(): 'n_years' must be a positive integer")
  }
  
  # Generate base series according to distribution
  vec <- switch(opts$model,
                normal = rnorm(n_years, opts$mean, opts$sd),
                
                lognormal = {
                  # Convert mean/sd to lognormal parameters
                  mu    <- log(opts$mean^2 / sqrt(opts$sd^2 + opts$mean^2))
                  sigma <- sqrt(log(1 + opts$sd^2 / opts$mean^2))
                  rlnorm(n_years, mu, sigma)
                },
                
                beta = rbeta(n_years, opts$shape1, opts$shape2),
                
                gamma = rgamma(n_years, 
                               shape = opts$shape1,
                               scale = opts$mean / opts$shape1),
                
                stop("Unknown model: ", opts$model)
  )
  
  # Apply temporal pattern: block timing
  if (opts$timing == "block") {
    # Parse block years
    blk <- as.integer(opts$block_years)
    blk <- blk[!is.na(blk) & blk >= 1 & blk <= n_years]
    
    # Reset non-block years to mean
    if (length(blk) > 0) vec[-blk] <- opts$mean
  }
  
  # Apply temporal pattern: pulse timing
  if (opts$timing == "pulse") {
    # Start with constant baseline
    vec[] <- opts$mean
    
    # Parse pulse years
    yrs_idx <- as.numeric(opts$pulse_years)
    yrs_idx <- yrs_idx[!is.na(yrs_idx)]
    
    # Clamp to valid range
    yrs_idx <- yrs_idx[yrs_idx >= 1 & yrs_idx <= n_years]
    
    # Apply pulses
    if (length(yrs_idx) > 0) {
      vec[yrs_idx] <- rnorm(length(yrs_idx), opts$mean, opts$pulse_sd)
    }
  }
  
  # Ensure non-negative (SARs cannot be negative)
  pmax(vec, 0)
}

# ═══════════════════════════════════════════════════════════════════════════════
# END OF FUNCTIONS.R
# ═══════════════════════════════════════════════════════════════════════════════