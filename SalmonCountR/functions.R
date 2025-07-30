# functions.R

library(lubridate)
library(tibble)

##########################
#######TDM MODELS#########
##########################

#’ -----------------------------------------------------------------------------
#’ Estimate egg development time in days from daily temperature
#’
#’ Calculates the number of days required for salmon eggs to hatch, assuming
#’ a constant requirement of 958 accumulated thermal units (ATU).
#’
#’ @param T  Numeric vector of daily water temperatures (°C).
#’
#’ @return Numeric vector of the same length as `T`, where each element is
#’         the estimated days to hatch for that day’s temperature.
#’
#’ @examples
#’ egg_model(10)       # 95.8 days at 10 °C
#’ egg_model(c(8,12))  # two values
#’ -----------------------------------------------------------------------------
egg_model <- function(T) 958 / T


#’ -----------------------------------------------------------------------------
#’ Estimate fry emergence time in days from daily temperature
#’
#’ Calculates the number of days required for salmon fry to emerge, assuming
#’ a constant requirement of 417 accumulated thermal units (ATU).
#’
#’ @param T  Numeric vector of daily water temperatures (°C).
#’
#’ @return Numeric vector of the same length as `T`, where each element is
#’         the estimated days to emergence for that day’s temperature.
#’
#’ @examples
#’ hatch_model(10)       # 41.7 days at 10 °C
#’ hatch_model(c(8,12))  # two values
#’ -----------------------------------------------------------------------------
hatch_model <- function(T) 417 / T


#’ -----------------------------------------------------------------------------
#’ Exponential temperature‑dependent mortality (TDM) model
#’
#’ Computes cumulative egg‑to‑fry survival using an exponential form:
#’   S = exp( - Σ α · exp(β · Tᵢ) )
#’ where (α, β) are calibration parameters.
#’
#’ @param temps  Numeric vector of daily temperatures (°C) during incubation.
#’ @param calib  Character; one of `"WaterForum2020"` or `"SALMOD2006"`,
#’               selecting the α and β calibration parameters.
#’
#’ @return Single numeric survival value (0–1) for the entire `temps` series.
#’
#’ @examples
#’ # Water Forum 2020 calibration
#’ tdm_exp(c(10,11,12), "WaterForum2020")
#’ # SALMOD 2006 calibration
#’ tdm_exp(c(10,11,12), "SALMOD2006")
#’ -----------------------------------------------------------------------------
tdm_exp <- function(temps, calib) {
  p <- list(
    WaterForum2020 = list(α = 3.40848e-11, β = 1.21122),
    SALMOD2006     = list(α = 1.475e-11,   β = 1.392)
  )[[calib]]
  exp(-sum(p$α * exp(p$β * temps)))
}



#’ -----------------------------------------------------------------------------
#’ Linear temperature‑dependent mortality (TDM) model
#’
#’ Computes cumulative egg‑to‑fry survival using a linear threshold form:
#’   S = exp( - α · Σ max(Tᵢ − β, 0) )
#’ where α and β are calibration parameters.
#’
#’ @param temps  Numeric vector of daily temperatures (°C) during incubation.
#’ @param α      Numeric mortality coefficient (default 0.026).
#’ @param β      Numeric threshold temperature (°C; default 12.14).
#’
#’ @return Single numeric survival value (0–1) for the entire `temps` series.
#’
#’ @examples
#’ # default Martin et al. 2017 parameters
#’ tdm_lin_martin(c(10,11,12))
#’ # custom parameters
#’ tdm_lin_martin(c(10,11,12), α = 0.03, β = 11)
#’ -----------------------------------------------------------------------------
tdm_lin_martin <- function(temps, α = 0.026, β = 12.14) {
  exp(-α * sum(pmax(temps - β, 0)))
}

#’ -----------------------------------------------------------------------------
#’ Compute egg‑to‑fry survival for a single redd (off‑line use in precompute)
#’
#’ Given a redd date, site, and pre‑indexed daily temperatures, computes the
#’ cumulative survival probability over the egg and fry incubation period.
#’
#’ @param rdr            Date of redd emergence (as a Date or character convertible to Date).
#’ @param site           Character; site name matching names(temps_env) and date_idx_env.
#’ @param date_idx_env   Named list mapping each site to a named integer vector:
#’                       names = dates (as “YYYY‑mm‑dd”), values = row indices in temps_env[[site]].
#’ @param temps_env      Named list of numeric temperature vectors, one per site.
#’ @param model          Character; `"exp"` or `"lin_martin"`, choosing the TDM variant.
#’ @param calib          Character or NA; for exponential variants, one of `"WaterForum2020"` or `"SALMOD2006"`.
#’
#’ @return Numeric survival probability (0–1) for that redd, or `NA` if data are missing.
#’
#’ @examples
#’ # assume date_idx_env and temps_env already built for “SiteA”
#’ compute_surv(as.Date("2015-10-01"), "SiteA", date_idx_env, temps_env, "exp", "WaterForum2020")
#’ -----------------------------------------------------------------------------
compute_surv <- function(rdr, site, date_idx_env, temps_env, model, calib) {
  pos <- date_idx_env[[site]][as.character(rdr)]
  if (is.na(pos) || pos < 1) return(NA_real_)
  T0 <- temps_env[[site]][pos]
  if (is.na(T0) || T0 <= 0) return(NA_real_)
  td <- ceiling(egg_model(T0) + hatch_model(T0))
  slice <- temps_env[[site]][ pos + seq_len(td) - 1L ]
  if (model == "exp") tdm_exp(slice, calib)
  else                tdm_lin_martin(slice)
}

#’ -----------------------------------------------------------------------------
#’ Compute adult pre‑spawn survival probability
#’
#’ @description
#’ A logistic (“inv.logit”) function that converts cumulative pre‑spawn degree‑days
#’ into a survival probability for adult salmon prior to spawning.
#’
#’ @param deg_day   Numeric vector of cumulative pre‑spawn degree‑days (°C·days).
#’ @param intercept Numeric. The logit intercept (default = 3.0).
#’ @param beta      Numeric. The logit slope per °C·day (default = –0.00067).
#’
#’ @return
#’ A numeric vector of the same length as `deg_day`, each entry in (0,1), giving
#’ the adult survival probability up to spawning.
#’
#’ @examples
#’ # No thermal stress (deg_day = 0) → high survival
#’ surv_adult_prespawn(0)
#’
#’ # As degree‑days accumulate, survival declines
#’ degs <- seq(0, 2000, by = 500)
#’ surv_adult_prespawn(degs, intercept = 3.0, beta = -0.00067)
#’
#’ @export
surv_adult_prespawn <- function(deg_day,
                                intercept = 3.0,
                                beta      = -0.00067) {
  # plogis(x) == exp(x)/(1+exp(x))
  plogis(intercept + beta * deg_day)
}


#' Compute cumulative adult pre‑spawn degree‑days (°C·days)
#'
#' For each brood year, this function sums daily‐mean water temperatures from a
#' specified start date (e.g. September 1 of the brood year) up to that year’s
#' observed spawn date.  Returns one °C·day total per year.
#'
#' @param env_nm        Character(1). Name of the scenario/alternative; must match
#'                      a component of `env_ext_list`.
#' @param sim_years     Integer vector. Brood years to simulate (e.g. `2011:2060`).
#' @param spawn_dates   Date vector. Observed spawn dates, same length as
#'                      `sim_years`.
#' @param env_ext_list  Named list of data.frames.  Each element must contain
#'                      columns `Date` (class Date) and `temp` (daily mean °C).
#' @param start_month   Integer in `[1,12]`. Month to begin accumulation
#'                      (default `9` = September).
#' @param start_day     Integer in `[1,31]`.  Day of month to begin
#'                      accumulation (default `1`).
#' @param base_temp     Numeric.  Temperature threshold (°C); only
#'                      `temp > base_temp` contribute to the sum
#'                      (default `0`).
#'
#' @return A numeric vector of length `length(sim_years)`.  Each entry is the
#'         sum over the selected date window of
#'         `pmax(daily_mean_temp - base_temp, 0)`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Suppose env_ext_list[["Alt1"]] has Date & temp columns:
#' years <- 2011:2015
#' spawns <- as.Date(c("2011-10-15","2012-10-20","2013-10-18",
#'                     "2014-10-22","2015-10-19"))
#'
#' compute_deg_day_adult(
#'   env_nm      = "Alt1",
#'   sim_years   = years,
#'   spawn_dates = spawns,
#'   env_ext_list= env_ext_list,
#'   start_month = 9,
#'   start_day   = 1,
#'   base_temp   = 0
#' )
#' }
compute_deg_day_adult <- function(
    env_nm,
    sim_years,
    spawn_dates,    # must already be Date
    env_ext_list,
    start_month = 10,
    start_day   = 1,
    base_temp   = 0
) {
  # Precompute one daily‐mean series per alternative
  daily_df <- env_ext_list[[env_nm]] %>%
    group_by(Date) %>%
    summarize(
      Tmean = mean(temp, na.rm = TRUE),
      .groups = "drop"
    )
  
  vapply(seq_along(sim_years), function(i) {
    season_start <- as.Date(sprintf(
      "%d-%02d-%02d",
      sim_years[i], start_month, start_day
    ))
    season_end <- spawn_dates[i]  # already Date
    
    temps <- daily_df$Tmean[
      daily_df$Date >= season_start &
        daily_df$Date <= season_end
    ]
    
    sum(pmax(temps - base_temp, 0), na.rm = TRUE)
  }, numeric(1))
}
#’ -----------------------------------------------------------------------------
#’ Simulate full life‑cycle including pre‑spawn and post‑spawn stages
#’
#’ Runs an age‑structured salmon life‑cycle model over a specified number of
#’ years, incorporating egg‑to‑fry survival, density dependence, rearing survival,
#’ smolt‑to‑adult ratio (SAR) variability, and adult pre‑spawn mortality.
#’
#’ @param surv_vec        Numeric vector of egg‑to‑fry survival rates (length ≥ years).
#’ @param P               List of life‑cycle parameters:
#’                        - female_fraction, fec, S0, K_spawners, S_rear, SAR_mean, lag_probs.
#’ @param years           Integer; number of years to simulate.
#’ @param S_init          Numeric vector of initial adult spawner counts (seed years).
#’ @param SAR_vec         Numeric vector of SAR values (length years; stochastic or constant).
#’ @param K_spawners_vec  Numeric vector of carrying capacities (redds) per year.
#’ @param deg_day_adult   Numeric vector of cumulative pre‑spawn °C·days per year.
#’ @param pre_int         Numeric intercept for logistic pre‑spawn survival (default 3.0).
#’ @param pre_beta        Numeric slope for logistic pre‑spawn survival (default –0.00067).
#’
#’ @return A tibble with columns:
#’   - year: 1:years  
#’   - spawners: predicted adult spawners each year  
#’   - dd: density‑dependence multiplier each year  
#’   - fry_dd: fry production after TDM & density dependence  
#’   - eff_surv: combined egg & density survival  
#’   - SAR_used: smolt‑to‑adult ratio used each year  
#’   - K_spawners: capacity each year  
#’   - pre_spawn: adult pre‑spawn survival each year  
#’
#’ @examples
#’ simulate_variant(
#’   surv_vec  = rep(0.5, 20),
#’   P         = list(female_fraction=0.5, fec=5000, S0=0.3,
#’                    K_spawners=10000, S_rear=0.8, SAR_mean=0.003,
#’                    lag_probs=c(`3`=0.75,`4`=0.249,`5`=0.001)),
#’   years     = 20,
#’   S_init    = c(100,120,110),
#’   SAR_vec   = rep(0.003, 20),
#’   K_spawners_vec = rep(10000, 20),
#’   deg_day_adult  = rep(50, 20)
#’ )
#’ -----------------------------------------------------------------------------
simulate_variant <- function(
    surv_vec,
    P,
    years,
    S_init,
    SAR_vec        = rep(P$SAR_mean,    years),
    K_spawners_vec = rep(P$K_spawners,  years),
    deg_day_adult  = NULL,               # vector of length ≥ years
    pre_int        = 3.0,                # default logistic intercept
    pre_beta       = -0.00067            # default logistic slope
) {
  # — pad or truncate egg‐survival series to 'years' —
  if (length(surv_vec) == 1) {
    surv_vec <- rep(surv_vec, years)
  } else if (length(surv_vec) < years) {
    surv_vec <- c(surv_vec, rep(tail(surv_vec,1), years - length(surv_vec)))
  } else {
    surv_vec <- surv_vec[1:years]
  }
  
  # — ensure deg_day_adult is numeric of length 'years' —
  if (is.null(deg_day_adult)) {
    deg_day_adult <- rep(0, years)
  } else {
    deg_day_adult <- as.numeric(deg_day_adult)[1:years]
    deg_day_adult[is.na(deg_day_adult)] <- 0
  }
  
  # — seed guard & allocation —
  seed_len <- min(length(S_init), years)
  S        <- numeric(years)
  reared   <- numeric(years)
  dd_vec   <- numeric(years)
  fry_dd   <- numeric(years)
  SAR_used <- numeric(years)
  S_pre    <- numeric(years)  # store pre‑spawn survival
  
  # initialize
  S[1:seed_len] <- S_init[1:seed_len]
  
  for (t in seq_len(years)) {
    # 1) adult pre‑spawn survival
    S_pre[t] <- surv_adult_prespawn(deg_day_adult[t], intercept = pre_int, beta = pre_beta)
    
    # 2) redds after pre‑spawn mortality
    redds    <- S[t] * P$female_fraction * S_pre[t]
    
    # 3) density‐dependence using K_spawners_vec[t]
    eggs     <- redds * P$fec
    Kt       <- K_spawners_vec[t]
    dd_vec[t]<- P$S0 / (1 + redds / Kt)
    
    # 4) fry production
    fry_dd[t]<- eggs * surv_vec[t] * dd_vec[t]
    
    # 5) rearing to smolt
    reared_t <- fry_dd[t] * P$S_rear
    reared[t]<- reared_t
    
    # 6) smolt‑to‑adult
    SAR_t       <- SAR_vec[t]
    SAR_used[t] <- SAR_t
    
    # 7) age‑structured returns (ages 3–5)
    for (age in 3:5) {
      ry <- t + age
      if (ry <= years) {
        S[ry] <- S[ry] + reared_t * SAR_t * P$lag_probs[as.character(age)]
      }
    }
  }
  
  tibble(
    year        = seq_len(years),
    spawners    = S,
    dd          = dd_vec,
    fry_dd      = fry_dd,
    eff_surv    = surv_vec * dd_vec,
    SAR_used    = SAR_used,
    K_spawners  = K_spawners_vec,
    pre_spawn   = S_pre
  )
}


#’ -----------------------------------------------------------------------------
#’ Generate a smolt‑to‑adult ratio (SAR) time series with optional stochastic timing
#’
#’ Draws a SAR vector of length `n_years` from one of Normal, Lognormal, Beta, or Gamma
#’ distributions, then optionally applies “block” or “pulse” timing adjustments.
#’
#’ @param n_years  Integer (length 1).  Number of years (length of the output vector).
#’ @param opts     List of parameters controlling the draw and timing.  Required elements:
#’                 
#’                 * `model` – one of `"normal"`, `"lognormal"`, `"beta"`, `"gamma"`.  
#’                 * `mean`  – numeric mean of the SAR distribution.  
#’                 * `sd`    – numeric standard deviation (for Normal or to match moments in Lognormal).  
#’                 * `shape1`, `shape2` – numeric shape parameters (for Beta/Gamma).  
#’                 * `timing`        – one of `"all"`, `"block"`, or `"pulse"`.  
#’                 * `block_years`   – integer vector of years to _keep_ stochastic (others reset to `mean`).  
#’                 * `pulse_years`   – integer vector of years to _pulse_ only (all others reset to `mean`).  
#’                 * `pulse_sd`      – numeric sd used for pulses (when `timing == "pulse"`).  
#’
#’ @return Numeric vector of length `n_years`.  Values are ≥ 0; when `timing == "all"`, it’s
#’         just the raw draw, otherwise non‑selected years are reset to `opts$mean` and
#’         selected “block” or “pulse” years draw fresh noise.  
#’
#’ @examples
#’ # 50-year Normal SAR series with base mean=0.003, sd=0.001, no timing filter
#’ opts <- list(
#’   model       = "normal",
#’   mean        = 0.003,
#’   sd          = 0.001,
#’   timing      = "all",
#’   shape1      = NA,  # unused for "all"
#’   shape2      = NA,
#’   block_years = integer(0),
#’   pulse_years = integer(0),
#’   pulse_sd    = NA
#’ )
#’ sar_all <- generate_SAR_vec(50, opts)
#’
#’ # 50-year series with a 10–20 block of stochastic noise only
#’ opts$timing      <- "block"
#’ opts$block_years <- 10:20
#’ sar_block <- generate_SAR_vec(50, opts)
#’
#’ # 50-year series with pulses in years 5, 15, 25
#’ opts$timing      <- "pulse"
#’ opts$pulse_years <- c(5, 15, 25)
#’ opts$pulse_sd    <- 0.002
#’ sar_pulse <- generate_SAR_vec(50, opts)
#’ -----------------------------------------------------------------------------
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



