# functions.R

library(lubridate)
library(tibble)

##########################
#######TDM MODELS#########
##########################

# — TDM helper functions —
egg_model <- function(T) 958 / T
hatch_model <- function(T) 417 / T

tdm_exp <- function(temps, calib) {
  p <- list(
    WaterForum2020 = list(α = 3.40848e-11, β = 1.21122),
    SALMOD2006     = list(α = 1.475e-11,   β = 1.392)
  )[[calib]]
  exp(-sum(p$α * exp(p$β * temps)))
}

tdm_lin_martin <- function(temps, α = 0.026, β = 12.14) {
  exp(-α * sum(pmax(temps - β, 0)))
}

# — Off‑line only: compute_surv used in precompute_data.R —
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

# — Life‑cycle model (used in both precompute & app) —
simulate_variant <- function(
    surv_vec, P, years, S_init, SAR_vec = rep(P$SAR_mean, years)
) {
  if (length(surv_vec) == 1) surv_vec <- rep(surv_vec, years)
  if (length(surv_vec) <  years) surv_vec <- c(surv_vec, rep(tail(surv_vec,1), years - length(surv_vec)))
  if (length(surv_vec) >  years) surv_vec <- surv_vec[1:years]
  
  n_init <- length(S_init)
  if (n_init > years) stop("…")
  seed_len <- min(n_init, years)
  
  S        <- numeric(years)
  reared   <- numeric(years)
  dd_vec   <- numeric(years)
  fry_dd   <- numeric(years)
  SAR_used <- numeric(years)
  
  S[1:seed_len] <- S_init[1:seed_len]
  
  for (t in seq_len(years)) {
    redds      <- S[t] * P$female_fraction
    eggs       <- redds * P$fec
    dd         <- P$S0 / (1 + redds / P$K_spawners)
    fry_dd[t]  <- eggs * surv_vec[t] * dd
    dd_vec[t]  <- dd
    reared[t]  <- fry_dd[t] * P$S_rear
    SAR_t      <- SAR_vec[t]
    SAR_used[t]<- SAR_t
    for (age in 3:5) {
      ry <- t + age
      if (ry <= years) {
        S[ry] <- S[ry] + reared[t] * SAR_t * P$lag_probs[as.character(age)]
      }
    }
  }
  
  tibble(
    year      = seq_len(years),
    spawners  = S,
    dd        = dd_vec,
    fry_dd    = fry_dd,
    eff_surv  = surv_vec * dd_vec,
    SAR_used  = SAR_used
  )
}

# — SAR‐generator (used in app and optionally in precompute) —
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
