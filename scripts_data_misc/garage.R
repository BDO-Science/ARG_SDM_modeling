# future brood‐years (2025–20??) by sampling from the real 14 years
future_years <- setdiff(sim_years, real_years)
sim_future <- map_dfr(future_years, function(y) {
  pick <- sample(real_years, 1)
  obs_df %>%
    filter(brood_year == pick) %>%
    transmute(
      sim_year     = y,
      spawn_dt,
      section,
      fork_length 
    )
})

# combine
sim_redds <- bind_rows(sim_actual, sim_future) %>%
  arrange(sim_year, spawn_dt)

# ───────────────────────────────────────────────────────────────────────────────
# Map sections → sites
# ───────────────────────────────────────────────────────────────────────────────
# picking the sections in the carcass surveys that apply to the sites we have forecasted/real water temperatures for
sim_redds <- sim_redds %>%
  mutate(site = case_when(
    section %in% c("NB", "W", "1a","1b", "1a/1b", "2") ~ "AveHazel",
    section %in% c("3")    ~ "AveWatt"
  ))

# -----------------------------------------------------------------------------
# OUTMIGRANT LOGISTIC MODEL COVARIATE CALCULATIONS
# -----------------------------------------------------------------------------
# This section prepares per-redd covariates for the outmigration logistic model.
# We compute temperature-dependent hatch times, spawn dates, and habitat treatment
# indicators for each redd, organized by environment.

# Helper: pad a vector of length 14 (calibration years) up to n_sim years by
# repeating its last value. This ensures static covariates extend through forecasts.
pad <- function(x) c(x, rep(tail(x,1), n_sim - length(x)))

# Define a static habitat treatment vector: treatments applied only in real calibration years.
# Here, real_years 2015–2018 get treatment=1, others 0; then pad to n_sim.
hab_treat_real <- if_else(real_years %in% 2015:2018, 1L, 0L)
hab_treat_pad  <- pad(hab_treat_real)

# Ensure sim_redds is a data.table for efficient joins by site and date
setDT(sim_redds)   

# --- 0) Ensure sim_redds has the columns we need ------------------------------
if (!is.data.table(sim_redds)) setDT(sim_redds)
needed <- c("sim_year","spawn_dt","site","fork_length")
missing_cols <- setdiff(needed, names(sim_redds))
if (length(missing_cols)) {
  # add missing columns with NAs
  for (cc in missing_cols) sim_redds[, (cc) := if (cc == "fork_length") NA_real_ else NA]
}

# --- 1) Build env-specific covariates with proper interval averaging ----------
redd_covars_by_env <- lapply(names(env_ext_list), function(env_nm) {
  df_env <- env_ext_list[[env_nm]]
  
  # daily temps
  dt_env <- as.data.table(df_env)[, .(Date, site, temp)]
  setkey(dt_env, site, Date)
  
  # redds on sites covered by this env
  reds_env <- as.data.table(sim_redds)[site %in% unique(dt_env$site),
                                       .(sim_year, site, spawn_dt)]
  if (nrow(reds_env) == 0L) {
    return(data.table(sim_year = integer(), mean_temp = numeric(),
                      spawn_DOY = integer(), hab_treat = integer()))
  }
  
  # join (explicit projection to avoid i.* names)
  reds_w_temp0 <- dt_env[
    reds_env, on = .(site, Date = spawn_dt), nomatch = 0L,
    .(sim_year = i.sim_year, site, spawn_dt = i.spawn_dt, temp)
  ]
  
  # hatch days from spawn-day temp
  reds_w_temp0[, hatch_days := ifelse(temp > 0, ceiling(958/temp), NA_integer_)]
  
  # intervals
  reds_iv <- reds_w_temp0[!is.na(hatch_days),
                          .(sim_year, site,
                            start = spawn_dt,
                            end   = spawn_dt + pmax(hatch_days - 1L, 0L))]
  
  # if no valid intervals, return NAs for mean_temp and finish
  if (nrow(reds_iv) == 0L) {
    out0 <- reds_w_temp0[, .(sim_year, spawn_dt)]
    out0[, raw_j := yday(spawn_dt)]
    out0[, ord_j := ifelse(month(spawn_dt) >= 8L, raw_j, raw_j + 365L)]
    out0[, spawn_DOY := ord_j - min(ord_j, na.rm = TRUE) + 1L, by = .(site, sim_year)]
    out0[, hab_treat := hab_treat_pad[sim_year - real_years[1] + 1L]]
    out0[, mean_temp := NA_real_]
    return(out0[, .(sim_year, mean_temp, spawn_DOY, hab_treat)])
  }
  
  # daily temps as intervals for overlap join
  dt_env_iv <- copy(dt_env)[, `:=`(start = Date, end = Date)]
  setkey(dt_env_iv, site, start, end)
  setkey(reds_iv,   site, start, end)
  
  # overlap and average
  ov <- foverlaps(dt_env_iv, reds_iv, nomatch = 0L)
  mean_by_redd <- ov[, .(mean_temp = mean(temp, na.rm = TRUE)),
                     by = .(sim_year, site, start, end)]
  mean_by_redd[, spawn_dt := start][]
  
  # season-aligned DOY (Aug–Jul)
  mean_by_redd[, raw_j := yday(spawn_dt)]
  mean_by_redd[, ord_j := ifelse(month(spawn_dt) >= 8L, raw_j, raw_j + 365L)]
  mean_by_redd[, spawn_DOY := ord_j - min(ord_j, na.rm = TRUE) + 1L, by = .(site, sim_year)]
  mean_by_redd[, hab_treat := hab_treat_pad[sim_year - real_years[1] + 1L]]
  
  mean_by_redd[, .(sim_year, mean_temp, spawn_DOY, hab_treat)]
})
names(redd_covars_by_env) <- names(env_ext_list)

# reference env & calibration years
ref_env <- names(redd_covars_by_env)[1]

surv_lookup_by_variant <- results_obs_fast %>%
  filter(env == ref_env, sim_year <= max(real_years)) %>%
  arrange(variant, sim_year) %>%
  group_by(variant) %>%
  summarise(surv_vec = list(mean_cum_surv), .groups="drop") %>%
  deframe()

#  4b) full horizon by env_variant key
surv_lookup_full <- egg_summary %>%
  arrange(env,variant,sim_year) %>%
  group_by(env,variant) %>%
  summarise(surv_vec=list(mean_cum_surv), .groups="drop") %>%
  mutate(key=paste(env,variant,sep="_")) %>%
  select(key,surv_vec) %>%
  deframe()

# Subset covariates to calibration years only (2011–2024)
dt_ref  <- redd_covars_by_env[[ref_env]][sim_year %in% real_years]

#  4d) compute `scales` for `predict_outmig()`
# Compute z-scaling parameters for continuous predictors from reference env
# These scales feed into predict_outmig() for covariate standardization.
scales <- list(
  fork_length_mean = mean(dt_ref$fork_length, na.rm=TRUE),
  fork_length_sd   = sd(  dt_ref$fork_length, na.rm=TRUE),
  spawn_DOY_mean   = mean(dt_ref$spawn_DOY,   na.rm=TRUE),
  spawn_DOY_sd     = sd(  dt_ref$spawn_DOY,   na.rm=TRUE),
  mean_temp_mean   = mean(dt_ref$mean_temp,   na.rm=TRUE),
  mean_temp_sd     = sd(  dt_ref$mean_temp,   na.rm=TRUE)
)

beta_outmig <- c("(Intercept)"=0,hab_treat=0.51, N_adults=-0.83,
                 fork_length=0.16, spawn_DOY=-0.56, mean_temp=-0.80)

#’ -----------------------------------------------------------------------------
#’ Predict fry-to-smolt out-migration probability
#’
#’ @description
#’ Computes the probability that fry from a given redd out-migrate (survive to smolt),
#’ using a logistic model with both categorical (habitat treatment) and scaled continuous covariates.
#’
#’ @param hab_treat    Numeric (0/1).  Habitat treatment indicator (0 = control, 1 = treated).
#’ @param N_adults     Numeric.  Number of adult female spawners this brood year.
#’ @param fork_length  Numeric vector.  Female fork lengths (mm) for each redd.
#’ @param spawn_DOY    Numeric vector.  Day-of-year of spawning for each redd.
#’ @param mean_temp    Numeric vector.  Mean incubation temperature (°C) for each redd.
#’ @param scales       Named list with elements  
#’                     - `N_adults_mean`, `N_adults_sd`  
#’                     - `fork_length_mean`, `fork_length_sd`  
#’                     - `spawn_DOY_mean`, `spawn_DOY_sd`  
#’                     - `mean_temp_mean`, `mean_temp_sd`  
#’
#’ @details
#’ 1. Each continuous predictor \(X\) is z-transformed:  
#’    \deqn{Z = \frac{X - \mu_X}{\sigma_X},}{Z = (X - mean_X) / sd_X,}  
#’    using means/SDs from `scales`.  
#’ 2. The linear predictor is  
#’    \deqn{ 
#’      \eta = \beta_0 
#’        + \beta_1 \times \text{hab\_treat}
#’        + \beta_2 \times Z(N\_adults)
#’        + \beta_3 \times Z(\text{fork\_length})
#’        + \beta_4 \times Z(\text{spawn\_DOY})
#’        + \beta_5 \times Z(\text{mean\_temp})
#’    }{ 
#’      lp = beta_outmig["(Intercept)"]
#’        + beta_outmig["hab_treat"]   * hab_treat 
#’        + beta_outmig["N_adults"]    * z_adults
#’        + beta_outmig["fork_length"] * z_FL
#’        + beta_outmig["spawn_DOY"]   * z_DOY
#’        + beta_outmig["mean_temp"]   * z_temp
#’    }  
#’ 3. The probability is \(\text{plogis}(\eta)\).  
#’ 4. Any NA in the inputs yields NA in the output.
#’
#’ @return
#’ A numeric vector of out-migration probabilities in (0,1), same length as `fork_length` etc.
#’
#’ @references
#’ Blankenship, L. et al. (2024). Logistic “zero” component for fry-to-smolt out-migration.  
#’ *Aquatic Ecology*, Table 3.
#’
#’ @examples
#’ scales <- list(
#’   N_adults_mean    = 1000, N_adults_sd    = 200,
#’   fork_length_mean = 800,  fork_length_sd = 100,
#’   spawn_DOY_mean   = 45,   spawn_DOY_sd   = 10,
#’   mean_temp_mean   = 14,   mean_temp_sd   = 1.2
#’ )
#’ # Single‐redd example
#’ predict_outmig(
#’   hab_treat   = 1,
#’   N_adults    = 500,
#’   fork_length = 750,
#’   spawn_DOY   = 50,
#’   mean_temp   = 15,
#’   scales      = scales
#’ )
#’
#’ @export
predict_outmig <- function(hab_treat,
                           N_adults, fork_length,
                           spawn_DOY, mean_temp,
                           scales) {
  
  print("Inside predict_outmig")
  print(beta_outmig)
  print(str(beta_outmig))
  print(list(
    hab_treat=hab_treat, N_adults=N_adults,
    fork_length=fork_length, spawn_DOY=spawn_DOY,
    mean_temp=mean_temp
  ))
  
  # z-scale each continuous predictor
  z_adults   <- (N_adults    - scales$N_adults_mean)    / scales$N_adults_sd
  z_FL       <- (fork_length - scales$fork_length_mean) / scales$fork_length_sd
  z_DOY      <- (spawn_DOY   - scales$spawn_DOY_mean)   / scales$spawn_DOY_sd
  z_temp     <- (mean_temp   - scales$mean_temp_mean)   / scales$mean_temp_sd
  
  # build linear predictor
  lp <- beta_outmig["(Intercept)"] +
    beta_outmig["hab_treat"]   * hab_treat +
    beta_outmig["N_adults"]    * z_adults +
    beta_outmig["fork_length"] * z_FL +
    beta_outmig["spawn_DOY"]   * z_DOY +
    beta_outmig["mean_temp"]   * z_temp
  
  print(lp)
  plogis(lp)
}

#' SSE objective for integrated calibration (SAR only)
#'
#' Fits a single SAR_mean to observed spawners assuming egg→smolt survival
#' is entirely captured outside the TDM (i.e., survival vector of ones).
#'
#' @param par Numeric length-1: trial `SAR_mean`.
#'
#' @return Numeric scalar SSE.
#' @export

integrated_sse <- function(par) {
  # par[1] = SAR_mean (the only parameter being optimized)
  P_tmp <- base_P
  P_tmp$SAR_mean <- par[1]
  
  # For integrated, no TDM: survival vector all 1s
  out <- simulate_variant(
    surv_vec       = rep(1, 14),                # No TDM in this calibration
    P              = P_tmp,
    years          = 14,
    S_init         = S_seed_calib,
    SAR_vec        = rep(par[1], 14),
    K_spawners_vec = rep(NA, 14),               # Not used for integrated method
    deg_day_adult  = rep(0, 14),
    sim_years_vec  = real_years,
    surv_method    = "integrated"
  )
  
  preds <- out$spawners[fit_idx]                # Model predictions for years being fit
  
  # Troubleshooting output
  cat(sprintf("par: %f, preds: %s\n", par[1], paste(round(preds, 1), collapse = ", ")))
  
  # If any predictions are not finite, return huge value (so optimizer avoids this)
  if (!all(is.finite(preds))) return(.Machine$double.xmax)
  
  # Calculate sum of squared errors
  sse <- sum((preds - obs_spawners[fit_idx])^2)
  if (!is.finite(sse)) sse <- .Machine$double.xmax
  return(sse)
}

#' Calibrate SAR_mean (integrated method) via L-BFGS-B
#'
#' @param init_SAR Numeric, initial SAR_mean (default 0.0025).
#'
#' @return Tibble with `year, observed, predicted, SAR_mean, rear_surv, sse`.
#' @export

run_integrated_calibration <- function(init_SAR = 0.0025) {
  opt <- optim(
    par    = init_SAR,
    fn     = integrated_sse,
    method = "L-BFGS-B",
    lower  = 0,
    upper  = 1
  )
  SAR_fit <- opt$par[1]
  P_tmp   <- base_P
  P_tmp$SAR_mean <- SAR_fit
  sim_out <- simulate_variant(
    surv_vec       = rep(1, 14),
    P              = P_tmp,
    years          = 14,
    S_init         = S_seed_calib,
    SAR_vec        = rep(SAR_fit, 14),
    K_spawners_vec = rep(NA, 14),
    deg_day_adult  = rep(0, 14),
    sim_years_vec  = real_years,
    surv_method    = "integrated"
  )
  tibble(
    year      = real_years,
    observed  = obs_spawners,
    predicted = sim_out$spawners,
    SAR_mean  = SAR_fit,
    rear_surv = NA_real_,
    sse       = sum((sim_out$spawners[fit_idx] - obs_spawners[fit_idx])^2)
  )
}