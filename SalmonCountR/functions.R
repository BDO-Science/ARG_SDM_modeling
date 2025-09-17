# functions.R

#' @importFrom stats plogis rnorm
#' @importFrom lubridate month day mday make_date
#' @importFrom dplyr group_by summarize mutate transmute
#' @importFrom tibble tibble
#' @importFrom purrr map
#' @importFrom tidyr pivot_wider
NULL

##########################
#######TDM MODELS#########
##########################

#' Estimate egg development time (days) from temperature
#'
#' Computes the number of days required for salmon eggs to reach hatching,
#' assuming a constant requirement of 958 accumulated thermal units (ATU).
#'
#' @param T Numeric vector of daily mean water temperatures (°C).
#'
#' @return Numeric vector of the same length as `T` with days to hatching.
#' @examples
#' hatch_model(10)             # ~95.8 days at 10 °C
#' hatch_model(c(8, 12))
#' @export
hatch_model <- function(T) 958 / T


#' Estimate fry emergence time (days) from temperature
#'
#' Computes days until fry emergence, assuming 417 accumulated thermal units (ATU).
#'
#' @param T Numeric vector of daily mean water temperatures (°C).
#'
#' @return Numeric vector of the same length as `T` with days to emergence.
#' @examples
#' emergence_model(10)           # ~41.7 days at 10 °C
#' emergence_model(c(8, 12))
#' @export
emergence_model <- function(T) 417 / T


# -------------------------------------------------------------------
# ATU stage boundaries (°C·days) used for egg→hatch→emergence splits
# -------------------------------------------------------------------
egg_ATU   <- 958     # to hatch
alev_ATU  <- 417     # hatch -> emergence
total_ATU <- egg_ATU + alev_ATU

#' Locate hatch and emergence day indices from a temperature series
#'
#' Uses cumulative ATU (= cumsum(temps)) to find the first day at which
#'  - hatch is reached (>= 958 °C·days), and
#'  - emergence is reached (>= 958 + 417 °C·days).
#'
#' Returns integer indices in 1..length(temps); if a threshold is never reached,
#' the corresponding index is set to the last day.
#'
#' @param temps Numeric vector of daily mean temperatures (°C).
#' @param egg_atu Numeric, ATU to hatch (default 958).
#' @param total_atu Numeric, ATU to emergence (default 1375).
#' @return List with integer elements `hatch`, `emerge`.
#' @keywords internal
.stage_indices_by_atu <- function(temps, egg_atu = egg_ATU, total_atu = total_ATU) {
  if (length(temps) == 0L || all(!is.finite(temps))) {
    return(list(hatch = 0L, emerge = 0L))
  }
  atu <- cumsum(pmax(temps, 0))                  # accumulated °C·days
  i_h <- which(atu >= egg_atu)[1]
  i_e <- which(atu >= total_atu)[1]
  if (is.na(i_h)) i_h <- length(temps)
  if (is.na(i_e)) i_e <- length(temps)
  list(hatch = as.integer(i_h), emerge = as.integer(i_e))
}

#' Exponential TDM with stage-specific hazards (egg vs alevin)
#'
#' Computes cumulative egg→fry survival under an exponential
#' temperature-dependent mortality model with **stage-specific**
#' parameters (egg vs. alevin), switching stages at **hatch** and
#' stopping at **emergence** based on accumulated thermal units (ATU).
#'
#'#' @details
#' When \code{use_stages = TRUE}, hazards switch from egg to alevin parameters
#' at the ATU-inferred hatch day, and stop at emergence (total 1375 °C·days).
#' Setting \code{use_stages = FALSE} applies egg parameters to all days for
#' backward compatibility.
#'
#' Daily hazard is \eqn{h(T) = \alpha \exp(\beta T)} and cumulative
#' survival is \eqn{S = \exp\!\big(-\sum h(T_i)\big)}.
#'
#' Calibrations implemented:
#' \itemize{
#'   \item \strong{WaterForum2020}:
#'         egg \eqn{\alpha=3.408486\times 10^{-11}}, \eqn{\beta=1.21122};
#'         alevin \eqn{\alpha=1.01755\times 10^{-10}}, \eqn{\beta=1.24092}.
#'   \item \strong{SALMOD2006}:
#'         egg \eqn{\alpha=1.475\times 10^{-11}}, \eqn{\beta=1.392};
#'         alevin \eqn{\alpha=2.521\times 10^{-12}}, \eqn{\beta=1.461}.
#' }
#'
#' Set \code{use_stages = FALSE} to apply the egg parameters to all days
#' (backward-compatible single-set behavior).
#'
#' @param temps Numeric vector of daily incubation temperatures (°C).
#' @param calib Character(1), one of \code{"WaterForum2020"} or \code{"SALMOD2006"}.
#' @param use_stages Logical, if \code{TRUE} (default) split hazards at hatch/emergence.
#'
#' @return Scalar survival in [0, 1].
#' @examples
#' # Constant 11 °C for 70 days, Water Forum 2020:
#' tdm_exp(rep(11, 70), "WaterForum2020")
#'
#' # SALMOD 2006 with stage split disabled (legacy behavior):
#' tdm_exp(rep(12, 60), "SALMOD2006", use_stages = FALSE)
#' @export
tdm_exp <- function(temps,
                    calib = c("WaterForum2020", "SALMOD2006"),
                    use_stages = TRUE) {
  calib <- match.arg(calib)
  if (length(temps) == 0L) return(NA_real_)
  
  # Parameter sets (egg vs alevin) for each calibration
  pars <- switch(
    calib,
    "WaterForum2020" = list(
      egg    = list(alpha = 3.408488e-11, beta = 1.21122),
      alevin = list(alpha = 1.017554e-10,  beta = 1.24092)
    ),
    "SALMOD2006" = list(
      egg    = list(alpha = 1.475e-11,  beta = 1.392),
      alevin = list(alpha = 2.521e-12,  beta = 1.461)
    )
  )
  
  haz <- function(T, a, b) a * exp(b * T)
  
  if (!use_stages) {
    # Back-compat: apply egg params across all days
    H <- sum(haz(temps, pars$egg$alpha, pars$egg$beta), na.rm = TRUE)
    return(exp(-H))
  }
  
  # Split the incubation by ATU-inferred hatch and emergence
  idx <- .stage_indices_by_atu(temps)
  i_h <- max(0L, min(idx$hatch,  length(temps)))
  i_e <- max(i_h, min(idx$emerge, length(temps)))
  
  egg_slice    <- if (i_h > 0L)      temps[1:i_h]        else numeric(0)
  alevin_slice <- if (i_e > i_h)     temps[(i_h+1):i_e]   else numeric(0)
  
  H_egg <- if (length(egg_slice))
    sum(haz(egg_slice,    pars$egg$alpha,    pars$egg$beta),    na.rm = TRUE) else 0
  H_alv <- if (length(alevin_slice))
    sum(haz(alevin_slice, pars$alevin$alpha, pars$alevin$beta), na.rm = TRUE) else 0
  
  exp(-(H_egg + H_alv))
}


#' Linear threshold TDM (Martin et al. 2017)
#'
#' Cumulative survival:
#' \deqn{S = \exp\left(-\alpha \sum \max(T_i - \beta, 0)\right)}
#'
#' @param temps Numeric vector of daily incubation temperatures (°C).
#' @param α Numeric, mortality coefficient (default 0.026).
#' @param β Numeric, threshold temperature (°C; default 12.14).
#'
#' @return Scalar survival in [0, 1].
#' @references Martin, B.T. et al. (2017) *Ecology Letters* 20:50–59.
#' @examples
#' tdm_lin_martin(c(10, 11, 12))
#' tdm_lin_martin(c(10, 11, 12), α = 0.03, β = 11)
#' @export

tdm_lin_martin <- function(temps, α = 0.026, β = 12.14) {
  exp(-α * sum(pmax(temps - β, 0)))
}

#' Compute egg→fry survival using ATU-precise incubation window
#'
#' Unlike \code{compute_surv()}, this function marches forward through realized
#' daily temperatures starting on the redd date until accumulated thermal units
#' (ATU) reach \eqn{958} (hatch) and then \eqn{958+417=1375} (emergence), and
#' evaluates the TDM over exactly that window.
#'
#' @inheritParams compute_surv
#' @param max_days Integer, safety cap for the look-ahead window (default 300).
#'
#' @return Numeric survival in [0,1]; \code{NA_real_} if inputs are missing or
#'   the series cannot be evaluated.
#' @seealso \code{\link{compute_surv}} for the constant-T window approximation.
#' @export

# March forward through realized daily temps until crossing ATU thresholds
.slice_by_atu <- function(temps, egg_atu = egg_ATU, total_atu = total_ATU) {
  if (!length(temps)) return(integer())
  atu <- cumsum(pmax(temps, 0))
  i_end <- which(atu >= total_atu)[1]
  if (is.na(i_end)) i_end <- length(temps)
  seq_len(i_end)
}

compute_surv_by_atu <- function(rdr, site, date_idx_env, temps_env, model, calib, max_days = 300L) {
  pos <- date_idx_env[[site]][as.character(rdr)]
  if (is.na(pos) || pos < 1) return(NA_real_)
  # take a long-enough look-ahead; 300 is conservative
  end_idx <- min(length(temps_env[[site]]), pos + max_days)
  Tvec <- temps_env[[site]][pos:end_idx]
  if (!length(Tvec) || all(!is.finite(Tvec))) return(NA_real_)
  idx <- .slice_by_atu(Tvec, egg_ATU, total_ATU)
  slice <- Tvec[idx]
  if (model == "exp") tdm_exp(slice, calib) else tdm_lin_martin(slice)
}

#' Adult pre-spawn survival from degree-days (logistic form)
#'
#' @param deg_day Numeric vector of cumulative pre-spawn degree-days (°C·days).
#' @param intercept Numeric, logit intercept (default 3.0).
#' @param beta Numeric, logit slope per °C·day (default -0.00067).
#'
#' @return Numeric vector in (0, 1).
#' @references Colvin et al. (2018) *RRA* 34(6):621–632.
#' @examples
#' surv_adult_prespawn(0)
#' surv_adult_prespawn(seq(0, 2000, 500))
#' @export

surv_adult_prespawn <- function(deg_day,
                                intercept = 3.0,
                                beta      = -0.00067) {
  # plogis(x) == exp(x)/(1+exp(x))
  plogis(intercept + beta * deg_day)
}

# Season-anchored timestamp (e.g., hydrologic year anchored at Oct-01)
# - x: Date or POSIXct
# - anchor_mmdd: "MM-DD", c(MM,DD), or a Date (month/day used; year ignored)
# - tz: output timezone if returning POSIXct
season_posix <- function(x, anchor_mmdd = c(10, 1), tz = "UTC") {
  stopifnot(inherits(x, "Date") || inherits(x, "POSIXct") || inherits(x, "POSIXt"))
  # parse anchor
  if (is.character(anchor_mmdd) && length(anchor_mmdd) == 1L) {
    mmdd <- strsplit(anchor_mmdd, "[-/\\.]")[[1]]
    mm <- as.integer(mmdd[1]); dd <- as.integer(mmdd[2])
  } else if (inherits(anchor_mmdd, "Date")) {
    mm <- lubridate::month(anchor_mmdd); dd <- lubridate::day(anchor_mmdd)
  } else if (is.numeric(anchor_mmdd) && length(anchor_mmdd) >= 2L) {
    mm <- as.integer(anchor_mmdd[1]); dd <- as.integer(anchor_mmdd[2])
  } else {
    stop("anchor_mmdd must be 'MM-DD', c(MM,DD), or a Date.")
  }
  if (anyNA(c(mm, dd))) stop("Could not parse anchor_mmdd.")
  
  # normalize input to POSIXct for comparison; remember original class
  is_date   <- inherits(x, "Date")
  is_posix  <- inherits(x, "POSIXct") || inherits(x, "POSIXt")
  x_posix   <- if (is_date) as.POSIXct(x) else as.POSIXct(x)
  yrs       <- lubridate::year(x_posix)
  
  # anchor date in the SAME calendar year as the timestamp
  anchor_this_year <- lubridate::make_datetime(year = yrs, month = mm, day = dd, tz = tz)
  
  # season year: if timestamp is before the anchor in that calendar year, season_year = year - 1
  season_year <- ifelse(x_posix >= anchor_this_year, yrs, yrs - 1L)
  
  # rebuild a POSIXct keeping the original month/day/time but replacing the year with season_year
  out <- lubridate::make_datetime(
    year  = season_year,
    month = lubridate::month(x_posix),
    day   = lubridate::mday(x_posix),
    hour  = lubridate::hour(x_posix),
    min   = lubridate::minute(x_posix),
    sec   = lubridate::second(x_posix),
    tz    = tz
  )
  
  # always return POSIXct (as the name suggests)
  out
}

# (Optional) Just the season year as an integer, same anchor rules
season_year <- function(x, anchor_mmdd = c(10, 1)) {
  as.integer(lubridate::year(season_posix(x, anchor_mmdd)))
}

#' Build an env-specific spawn date vector aligned to simulation years
#'
#' Uses LOCF then backfill to fill gaps in env×year median dates, and coerces
#' types safely (Date/POSIXct/numeric-date/character). Final dates are rebuilt
#' with the target \code{sim_years} and the month/day from the filled series,
#' so prespawn degree-days end on the correct year’s date.
#'
#' @param df_env Data frame with columns \code{sim_year} (integer) and
#'   \code{spawn_dt} (Date/POSIXct/numeric-date/character).
#' @param sim_years Integer vector of simulation years.
#' @param fallback_md Integer length-2 \code{c(month, day)} used when an env
#'   has no dates at all (default \code{c(11, 20)}).
#' @return Date vector of length \code{sim_years}.
#' @export

# Robust env-specific spawn date vector (LOCF/backfill + fallback)
build_spawn_vec_for_env <- function(df_env, sim_years, fallback_md = c(11, 20)) {
  stopifnot(all(c("sim_year","spawn_dt") %in% names(df_env)))
  # Coerce spawn_dt to Date regardless of incoming type
  v_raw <- df_env$spawn_dt
  v <- if (inherits(v_raw, "Date")) {
    v_raw
  } else if (inherits(v_raw, c("POSIXct","POSIXt"))) {
    as.Date(v_raw)
  } else if (is.numeric(v_raw)) {
    as.Date(v_raw, origin = "1970-01-01")
  } else {
    # character (YYYY-mm-dd) or other → best effort
    suppressWarnings(as.Date(v_raw))
  }
  
  df_env2 <- df_env
  df_env2$spawn_dt <- v
  
  # Align to sim_years
  m <- match(sim_years, df_env2$sim_year)
  v <- df_env2$spawn_dt[m]  # Date vector (may be all NA)
  
  # If empty or all NA → fill with fallback month/day for all years
  if (length(v) == 0L || all(is.na(v))) {
    return(lubridate::make_date(year = sim_years,
                                month = fallback_md[1],
                                day   = fallback_md[2]))
  }
  
  # LOCF
  for (i in seq_along(v)) if (is.na(v[i]) && i > 1) v[i] <- v[i - 1]
  # Backfill from first known
  if (anyNA(v)) {
    first_ok <- which(!is.na(v))[1]
    if (length(first_ok) == 1L && !is.na(first_ok)) v[seq_len(first_ok - 1)] <- v[first_ok]
  }
  # Any residual NA → fallback on those positions
  if (anyNA(v)) {
    v[is.na(v)] <- lubridate::make_date(
      year  = sim_years[is.na(v)],
      month = fallback_md[1],
      day   = fallback_md[2]
    )
  }
  
  # Rebuild dates with target year + month/day from v
  lubridate::make_date(
    year  = sim_years,
    month = lubridate::month(v),
    day   = lubridate::day(v)
  )
}


#' Cumulative adult pre-spawn degree-days for each brood year
#'
#'#' @details
#' Degree-days are summed from \code{start_month/start_day} in each \code{sim_year}
#' up to (and including) that year's spawn date. Supply env-specific spawn dates
#' (e.g., from \code{build_spawn_vec_for_env}) to ensure alternatives use the
#' correct timing.
#' 
#' Sums \code{max(T - base_temp, 0)} from a fixed start date (month/day)
#' in each brood year through that year's spawn date.
#' If `spawn_dates` length differs from `sim_years`, it is padded/truncated.
#'
#' @param env_nm Character(1), name found in `env_ext_list`.
#' @param sim_years Integer vector of brood years.
#' @param spawn_dates Date vector of spawn dates aligned to `sim_years`
#'   (padding/truncation performed if needed).
#' @param env_ext_list Named list of data.frames with `Date` (Date) and `temp` (°C).
#' @param start_month,start_day Integers for accumulation start (default 10/1).
#' @param base_temp Numeric, contribution threshold (default 0).
#'
#' @return Numeric vector of degree-days, length `sim_years`.
#' @examples
#' # compute_deg_day_adult("Alt1", 2011:2015, spawns, env_ext_list)
#' @export

compute_deg_day_adult <- function(
    env_nm,
    sim_years,
    spawn_dates,    # must already be Date
    env_ext_list,
    start_month = 10,
    start_day   = 1,
    base_temp   = 0
) {
  # 1) validate env_nm
  if (!is.character(env_nm) || length(env_nm) != 1 ||
      !env_nm %in% names(env_ext_list)) {
    stop("`env_nm` must be a single name in names(env_ext_list).")
  }
  
  # 2) pad or truncate spawn_dates to match sim_years
  if (length(spawn_dates) < length(sim_years)) {
    spawn_dates <- c(
      spawn_dates,
      rep(tail(spawn_dates, 1), length(sim_years) - length(spawn_dates))
    )
  } else if (length(spawn_dates) > length(sim_years)) {
    spawn_dates <- spawn_dates[seq_len(length(sim_years))]
  }
  
  # 3) build daily‐mean temperature series
  daily_df <- env_ext_list[[env_nm]] %>%
    dplyr::group_by(Date) %>%
    dplyr::summarize(Tmean = mean(temp, na.rm = TRUE), .groups = "drop")
  
  # 4) compute degree‐days for each brood year
  vapply(seq_along(sim_years), function(i) {
    start_date <- as.Date(sprintf(
      "%d-%02d-%02d",
      sim_years[i], start_month, start_day
    ))
    end_date <- spawn_dates[i]
    
    temps <- daily_df$Tmean[
      daily_df$Date >= start_date &
        daily_df$Date <= end_date
    ]
    if (length(temps) == 0) return(NA_real_)
    sum(pmax(temps - base_temp, 0), na.rm = TRUE)
  }, numeric(1))
}

# helper: adult degree-days for a given env across calibration years
deg_day_cal_for <- function(env_nm) {
  # spawn_dates_by_env[[env_nm]] is a Date vector indexed by sim_years
  sd_env <- spawn_dates_by_env[[env_nm]]
  stopifnot(length(sd_env) >= length(sim_years))
  # align to real calibration years
  sd_cal <- sd_env[match(real_years, sim_years)]
  compute_deg_day_adult(
    env_nm       = env_nm,
    sim_years    = real_years,
    spawn_dates  = sd_cal,
    env_ext_list = env_ext_list
  )
}

# ---- 2) Helper to get probabilities from CLM (optionally with i.i.d. offset)
predict_clm_probs <- function(beta, zeta, newdata, offset = 0) {
  X  <- as.matrix(newdata[, names(beta), drop = FALSE])
  xb <- as.vector(X %*% beta) + as.numeric(offset)
  t(vapply(xb, function(xi) {
    cp <- plogis(zeta - xi)      # cumulative logits → cum probs
    p  <- c(cp, 1) - c(0, cp)    # category probs
    p / sum(p)
  }, numeric(length(zeta) + 1)))
}

# 1) Build forecast_temps (env × sim_year = forecast_years)
build_forecast_temps <- function(env_ext_list, yrs_forecast, sc) {
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
  
  # Carry forward last available Oct/Nov if a forecast year isn't in the series
  purrr::map_dfr(split(raw, raw$env), function(d) {
    if (!nrow(d)) return(tibble(env = unique(d$env), sim_year = yrs_forecast, Oct = NA_real_, Nov = NA_real_))
    want <- tibble(sim_year = yrs_forecast)
    have <- d %>% dplyr::rename(sim_year = year) %>% dplyr::select(env, sim_year, Oct, Nov)
    
    # for any missing (env, year), fill with last available year for that env
    missing_years <- setdiff(yrs_forecast, have$sim_year)
    cf_rows <- if (length(missing_years)) {
      last_row <- d %>% dplyr::arrange(dplyr::desc(year)) %>% dplyr::slice(1) %>% dplyr::transmute(env = unique(env), Oct, Nov)
      tidyr::crossing(last_row, sim_year = missing_years)
    } else tibble()
    
    dplyr::bind_rows(have, cf_rows) %>%
      dplyr::arrange(sim_year)
  }) %>%
    dplyr::mutate(
      Oct_std = (Oct - sc$o_m) / sc$o_s,
      Nov_std = (Nov - sc$n_m) / sc$n_s
    )
}

#' SSE objective for joint calibration of SAR_mean and rear_surv (modular)
#'
#' Uses TDM egg→fry survivals for the first 14 calibration years in a single
#' alternative/environment and variant, then computes SSE against observed
#' spawners for fit years (4:14).
#'
#' @param par Numeric length-2: `c(SAR_mean, rear_surv)`.
#' @param alt Character(1), environment/alternative key.
#' @param variant Character(1), TDM variant key.
#'
#' @return Numeric scalar SSE.
#' @export

# objective: SSE between predicted and observed spawners
modular_sse <- function(par, variant) {
  P_tmp <- base_P
  P_tmp$SAR_mean  <- par[1]
  P_tmp$rear_surv <- par[2]
  
  out <- simulate_variant(
    surv_vec       = surv_lookup_by_variant[[variant]][1:n_calib],
    P              = P_tmp,
    years          = n_calib,
    S_init         = S_seed_calib,
    SAR_vec        = rep(P_tmp$SAR_mean,  n_calib),
    K_spawners_vec = rep(P_tmp$K_spawners, n_calib),
    deg_day_adult  = deg_day_cal_ref,      # <-- KEY FIX: real degree-days, not zeros
    sim_years_vec  = real_years
  )
  
  preds <- out$spawners[fit_idx]
  if (!all(is.finite(preds))) return(.Machine$double.xmax)
  sum((preds - obs_spawners[fit_idx])^2)
}

#' Build a forecast simulator for one (variant, environment)
#'
#' Returns a zero-argument function that simulates the life cycle over the next
#' `n_sim` brood years using TDM survivals and age-structured returns.
#'
#' @param var_nm Character(1), TDM variant name (e.g., "exp_WF").
#' @param env_nm Character(1), environment/alternative key.
#' @param flow_cfs Numeric or `NULL`. If non-NULL, used to set `K_spawners`.
#' @param S_seed Numeric vector of initial spawner abundances to seed the first
#'   `length(S_seed)` years of the run.
#' @param spawn_dates_by_env Optional named list env -> Date vector (length = length(sim_years)).
#'   If NULL, falls back to global `spawn_dates_vec` for backward compatibility.
#' @param spawn_dates_by_env Optional named list \code{env -> Date vector}
#'   (each vector aligned to \code{sim_years}). If provided, prespawn degree-days
#'   are ended at the env-specific date for each year; otherwise the function
#'   errors. (Backward-compat behavior using a global \code{spawn_dates_vec}
#'   has been removed to prevent cross-env timing bugs.)
#' @return A function with no arguments that returns a tibble:
#'   `year, spawners, deg_day, pre_spawn, dd, fry_dd, egg_surv, eff_surv,
#'    rear_surv, SAR_used, K_spawners, env, variant`.
#' @export
sim_forecast_fn <- function(var_nm,
                            env_nm,
                            flow_cfs = NULL,
                            S_seed,
                            spawn_dates_by_env) {  # <-- NEW arg (a named list: env -> Date vector)
  force(var_nm); force(env_nm); force(flow_cfs); force(S_seed); force(spawn_dates_by_env)
  
  function() {
    P_tmp <- base_P_list[[var_nm]][[env_nm]]
    if (!is.null(flow_cfs)) P_tmp$K_spawners <- get_K_spawners(flow_cfs)
    
    # 1) survivals for this env × variant (already env-specific)
    surv_vec <- surv_lookup_full[[paste(env_nm, var_nm, sep = "_")]]
    
    # 2) env-specific spawn dates (MUST exist & match sim_years)
    if (is.null(spawn_dates_by_env[[env_nm]])) {
      stop("spawn_dates_by_env[[", env_nm, "]] is NULL.")
    }
    spawn_dates_env <- spawn_dates_by_env[[env_nm]]
    if (length(spawn_dates_env) < length(sim_years)) {
      stop("spawn_dates_by_env[[", env_nm, "]] length < sim_years length.")
    }
    
    # 3) prespawn DD ends at (env-specific) spawn date
    deg_day_vec <- compute_deg_day_adult(
      env_nm       = env_nm,
      sim_years    = sim_years,
      spawn_dates  = spawn_dates_env,      # <-- key change
      env_ext_list = env_ext_list
    )
    
    # 4) K and SAR
    K_vec   <- rep(P_tmp$K_spawners, length(sim_years))
    SAR_vec <- if (use_stochastic_SAR) {
      generate_SAR_vec(length(sim_years), modifyList(stoch_SAR_opts, list(mean = P_tmp$SAR_mean)))
    } else {
      rep(P_tmp$SAR_mean, length(sim_years))
    }
    
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
    
    dplyr::mutate(sim_out, env = env_nm, variant = var_nm)
  }
}

#' Simulate full salmon life cycle (egg → fry → age-structured returns)
#'
#' @param surv_vec Numeric vector, egg→fry survival per year (TDM output).
#' @param P List of biological parameters; must include
#'   `female_fraction`, `fec`, `S0`, `rear_surv`, `K_spawners`.
#' @param years Integer, number of years to simulate.
#' @param S_init Numeric vector, initial spawners for the first years (seed).
#' @param SAR_vec Numeric vector, smolt-to-adult return ratios by brood year.
#' @param K_spawners_vec Numeric vector, carrying capacity by brood year.
#' @param deg_day_adult Numeric vector or `NULL`, prespawn degree-days by year
#'   (used by `surv_adult_prespawn`); if `NULL`, zeros are used.
#' @param sim_years_vec Integer vector of brood years.
#' @param pre_int Numeric, logistic intercept for prespawn survival (default 3.0).
#' @param pre_beta Numeric, logistic slope for prespawn survival (default -0.00067).
#'
#' @return Tibble with columns:
#'   `year, spawners, deg_day, pre_spawn, dd, fry_dd, egg_surv, eff_surv,
#'    rear_surv, SAR_used, K_spawners`.
#' @examples
#' # simulate_variant(surv_vec, P, years, S_init, SAR_vec, K_vec, NULL, 2011:2060)
#' @export

simulate_variant <- function(
    surv_vec, P, years, S_init, SAR_vec, K_spawners_vec,
    deg_day_adult = NULL,
    sim_years_vec,
    pre_int = 3.0, pre_beta = -0.00067
) {
  # repeat inputs to length 'years'
  surv_vec       <- rep_len(surv_vec, years)
  K_spawners_vec <- rep_len(K_spawners_vec, years)
  SAR_vec        <- rep_len(SAR_vec, years)
  seed_len       <- min(length(S_init), years)
  
  # storage
  S                 <- numeric(years)
  S_pre             <- numeric(years)
  deg_day_adult_vec <- if (is.null(deg_day_adult)) rep(0, years) else rep_len(deg_day_adult, years)
  dd_vec            <- numeric(years)
  fry_dd            <- numeric(years)
  rear_surv_vec     <- numeric(years)
  reared_vec        <- numeric(years)
  
  # seed initial years
  if (seed_len > 0) S[1:seed_len] <- S_init[1:seed_len]
  
  for (t in seq_len(years)) {
    if (t <= seed_len && !is.na(S_init[t])) S[t] <- S_init[t]
    
    # pre-spawn survival (logistic function of deg-day)
    S_pre[t] <- surv_adult_prespawn(deg_day_adult_vec[t], intercept = pre_int, beta = pre_beta)
    
    # redds and eggs
    redds <- S[t] * P$female_fraction * S_pre[t]
    eggs  <- redds * P$fec
    
    # density dependence + fry production
    dd_vec[t]        <- P$S0 / (1 + redds / K_spawners_vec[t])
    fry_dd[t]        <- eggs * surv_vec[t] * dd_vec[t]
    rear_surv_vec[t] <- P$rear_surv
    reared_vec[t]    <- fry_dd[t] * P$rear_surv
    
    # age-structured returns (3–5)
    for (age in 3:5) {
      ry <- t + age
      if (ry <= years) {
        S[ry] <- S[ry] + reared_vec[t] * SAR_vec[t] * P$lag_probs[as.character(age)]
      }
    }
  }
  
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

#' Generate a SAR (smolt-to-adult) time series with optional timing
#'
#' Draws SAR values from one of Normal, Lognormal, Beta, or Gamma distributions,
#' then optionally applies block- or pulse-year patterns.
#'
#'#' @details
#' Values are truncated at 0 (no negatives). When \code{timing = "block"}, only
#' indices listed in \code{block_years} are stochastic; all others are set to
#' the mean. When \code{timing = "pulse"}, all years are set to the mean except
#' \code{pulse_years}, which are drawn with \code{sd = pulse_sd}.
#'
#' @param n_years Integer(1), length of the series (≥ 1).
#' @param opts List of options:
#'   \describe{
#'     \item{model}{`"normal"`, `"lognormal"`, `"beta"`, or `"gamma"`.}
#'     \item{mean}{Numeric mean of the SAR distribution.}
#'     \item{sd}{Numeric standard deviation (Normal/Lognormal).}
#'     \item{shape1, shape2}{Shape parameters (Beta/Gamma).}
#'     \item{timing}{`"all"`, `"block"`, or `"pulse"`.}
#'     \item{block_years}{Integer indices to keep stochastic when `timing="block"`.}
#'     \item{pulse_years}{Integer indices to pulse when `timing="pulse"`.}
#'     \item{pulse_sd}{Numeric sd for pulses (when `timing="pulse"`).}
#'   }
#'
#' @return Non-negative numeric vector of length `n_years`.
#' @examples
#' generate_SAR_vec(30, list(model="normal", mean=0.003, sd=0.001, timing="all"))
#' @export

generate_SAR_vec <- function(n_years, opts) {
  # ensure n_years is a single integer
  n_years <- as.integer(n_years)[1]
  if (is.na(n_years) || n_years < 1) {
    stop("generate_SAR_vec(): 'n_years' must be a positive integer")
  }
  
  # draw the base series
  vec <- switch(opts$model,
                normal    = rnorm(n_years, opts$mean, opts$sd),
                lognormal = {
                  mu    <- log(opts$mean^2 / sqrt(opts$sd^2 + opts$mean^2))
                  sigma <- sqrt(log(1 + opts$sd^2 / opts$mean^2))
                  rlnorm(n_years, mu, sigma)
                },
                beta      = rbeta(n_years, opts$shape1, opts$shape2),
                gamma     = rgamma(n_years, shape = opts$shape1,
                                   scale = opts$mean / opts$shape1),
                stop("Unknown model: ", opts$model)
  )
  
  # block timing
  if (opts$timing == "block") {
    # reset everything outside the block
    blk <- as.integer(opts$block_years)
    blk <- blk[!is.na(blk) & blk >= 1 & blk <= n_years]
    # only reset if we have valid block indices
    if (length(blk) > 0) vec[-blk] <- opts$mean
  }
  
  # pulse timing
  if (opts$timing == "pulse") {
    # step 1: reset all to the mean
    vec[] <- opts$mean
    
    # step 2: parse the user’s pulse years
    yrs_idx <- as.numeric(opts$pulse_years)
    yrs_idx <- yrs_idx[!is.na(yrs_idx)]
    
    # step 3: clamp to [1, n_years]
    yrs_idx <- yrs_idx[yrs_idx >= 1 & yrs_idx <= n_years]
    
    # step 4: only pulse valid indices
    if (length(yrs_idx) > 0) {
      vec[yrs_idx] <- rnorm(length(yrs_idx), opts$mean, opts$pulse_sd)
    }
  }
  
  # never negative
  pmax(vec, 0)
}



