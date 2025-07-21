rm(list=ls())
# precompute_data.R

# ─── Libraries ────────────────────────────────────────────────────────────────
library(tidyverse); library(lubridate); library(furrr); library(data.table)
source("SalmonCountR/functions.R")
# (plus any others needed for the heavy loops)

# ─── 1) Read in raw inputs ────────────────────────────────────────────────────
env_ext_list     <- readRDS("env_ext_list.rds")   # or read from CSV, NWIS, etc.
df_all           <- readRDS("df_all.rds")
carcass_raw      <- read.csv("SalmonCountR/carcassdet_1752789274_15.csv")

# LOAD OBSERVED ESCAPEMENT --------------------------------------------------
esc_obs <- read_csv(
  "SalmonCountR/grandtab_1752793045_337.csv",
  col_types = cols(`End Year of Monitoring Period` = col_character(),
                   `Population Estimate`           = col_double())
) %>%
  # 1) pull the raw text, strip any “*” (or other non‑digits), parse as integer
  mutate(
    year = parse_number(`End Year of Monitoring Period`)
  ) %>%
  # 2) now you can safely filter by numeric year
  filter(year >= 2011, year <= 2024) %>%
  rename(spawners = `Population Estimate`) %>%
  arrange(year)

obs_spawners <- esc_obs$spawners           # length = 14
S_seed       <- obs_spawners[1:3]         # brood 2011–2013
n_calib      <- length(obs_spawners)      # 14
fit_idx      <- (length(S_seed) + 1):n_calib  # 4:14

# detect how many physical cores you have (on Windows this will report logical threads,
# so we subtract 1 to leave one for the OS)
n_threads <- parallel::detectCores(logical = TRUE) - 1

# but cap at your physical cores (6) to avoid oversubscription
n_workers <- min(6, n_threads)

plan(multisession, workers = n_workers)

# and for data.table’s own internal C threading:
setDTthreads(n_workers)

# ─── 2) Build sim_redds & run the TDM simulations (results_obs_fast, egg_summary)  
#     *This is your ~2 min chunk*  
#     (all of the code up through building egg_summary)

# Define TDM variants to run: exponential with two calibrations, and linear
tdm_defs <- tribble(
  ~model,       ~calib,           ~variant,
  "exp",        "WaterForum2020", "exp_WF",
  "exp",        "SALMOD2006",     "exp_SM",
  "lin_martin", NA,                "lin_Martin"
)

# 1) site → temp vectors
site_temps_list <- map(env_ext_list, ~ split(.x$temp, .x$site))

# 2) site → Date vectors
site_dates_list <- map(env_ext_list, ~ split(.x$Date, .x$site))

# 3) site → named‐index lookup
date_idx_list <- map(site_dates_list, function(dlist) {
  map(dlist, ~ setNames(seq_along(.x), as.character(.x)))
})

#prepping carcass data
obs_df <- carcass_raw %>%
  #filter(clip == "No Clip") %>%
  mutate(
    Date      = as.Date(surveydate),
    spawn_dt  = Date - days(7),
    # brood year runs from Sep 15 of brood_year to Apr 30 of brood_year+1
    brood_year = if_else(
      month(spawn_dt) >= 9, 
      year(spawn_dt), 
      year(spawn_dt) - 1
    ),
    section   = section
  ) %>%
  select(brood_year, spawn_dt, section)

# Which brood years are “real”?
real_years <- 2011:2024
n_real     <- length(real_years)   # 14
n_sim      <- 100
sim_years  <- real_years[1] + seq(0, n_sim - 1)  # e.g. 2011–2060

# ───────────────────────────────────────────────────────────────────────────────
# 1) Build sim_redds: actual brood years 2011–2024, then bootstrap beyond
# ───────────────────────────────────────────────────────────────────────────────
# A) real brood years
sim_actual <- obs_df %>%
  filter(brood_year %in% real_years) %>%
  rename(sim_year = brood_year)

# B) future brood‐years (2025–2060) by sampling from the real 14 years
future_years <- setdiff(sim_years, real_years)
sim_future <- map_dfr(future_years, function(y) {
  pick <- sample(real_years, 1)
  obs_df %>%
    filter(brood_year == pick) %>%
    transmute(
      sim_year = y,
      spawn_dt,
      section
    )
})

# C) combine
sim_redds <- bind_rows(sim_actual, sim_future) %>%
  arrange(sim_year, spawn_dt)

# ───────────────────────────────────────────────────────────────────────────────
# 2) Map sections → sites
# ───────────────────────────────────────────────────────────────────────────────
sim_redds <- sim_redds %>%
  mutate(site = case_when(
    section %in% c("NB", "W", "1a","1b", "1a/1b", "2") ~ "AveHazel",
    section %in% c("3")    ~ "AveWatt"
  ))

# ───────────────────────────────────────────────────────────────────────────────
# 3) Plug into your TDM loop exactly as before, but now sim_year is brood_year
# ───────────────────────────────────────────────────────────────────────────────

# 1) compile the “hot” functions
egg_model_c   <- cmpfun(egg_model)
hatch_model_c <- cmpfun(hatch_model)
tdm_exp_c     <- cmpfun(tdm_exp)
tdm_lin_c     <- cmpfun(tdm_lin_martin)

# convert sim_redds in place to a data.table
setDT(sim_redds)

# 3) the main loop
results_obs_fast <- future_map_dfr(sim_years, function(sim_yr) {
  # pull only that year — returns a data.table, so .() will work
  red_this <- sim_redds[sim_year == sim_yr]
  
  env_res <- lapply(names(env_ext_list), function(env_nm) {
    date_idx_env <- date_idx_list[[env_nm]]
    temps_env     <- site_temps_list[[env_nm]]
    env_dates     <- env_ext_list[[env_nm]]$Date
    
    # compute each redd’s calendar date in water‑year
    rdrs_sites <- red_this[
      , .(
        rdr = make_date(
          year  = if_else(month(spawn_dt) >= 9L, sim_yr, sim_yr+1L),
          month = month(spawn_dt),
          day   = day(spawn_dt)
        ),
        site
      )
    ][
      # filter to your env window
      rdr >= min(env_dates) & rdr <= max(env_dates)
    ]
    
    if (nrow(rdrs_sites) == 0L) return(NULL)
    
    # loop over TDM definitions
    variant_res <- lapply(seq_len(nrow(tdm_defs)), function(i) {
      model   <- tdm_defs$model[i]
      calib   <- tdm_defs$calib[i]
      variant <- tdm_defs$variant[i]
      
      survs <- sapply(
        seq_len(nrow(rdrs_sites)),
        function(j) {
          compute_surv(
            rdr          = rdrs_sites$rdr[j],
            site         = rdrs_sites$site[j],
            date_idx_env = date_idx_env,
            temps_env    = temps_env,
            model        = model,
            calib        = calib
          )
        },
        simplify = TRUE, USE.NAMES = FALSE
      )
      
      # build one small data.table
      data.table(
        sim_year      = sim_yr,
        env           = env_nm,
        variant       = variant,
        method        = "observed",
        mean_cum_surv = mean(survs, na.rm = TRUE)
      )
    })
    
    rbindlist(variant_res)
  })
  
  rbindlist(env_res)
}, .options = furrr_options(seed = TRUE))

# LIFE CYCLE MODEL: CALIBRATION & FORECAST 

# LOAD TDM OUTPUTS ----------------------------------------------------------
egg_summary <- results_obs_fast


# ─── 3) Calibrate the life‐cycle model and make base_P_list, S_seed, stoch_SAR_opts  
#     (everything up through saving your helper functions)

# SET PARAMETER LIST --------------------------------------------------------
base_P <- list(
  female_fraction = 0.5,
  fec             = 5522,
  S0              = 0.347,
  K_spawners      = 10000,
  S_rear          = 0.8,
  SAR_mean        = 0.00308,
  SAR_sd          = 0.001,
  lag_probs       = c(`3` = 0.75, `4` = 0.249, `5` = 0.001)
)

# 1) Build a surv‐lookup *by variant* (2011–2024 only)
surv_lookup_by_variant <- results_obs_fast %>%
  filter(sim_year <= 2024) %>%       # calibration window
  arrange(variant, sim_year) %>%
  group_by(variant) %>%
  summarise(
    surv_vec = list(mean_cum_surv),
    .groups = "drop"
  ) %>%
  deframe()                          # named list: $exp_WF, $exp_SM, $lin_Martin

# egg_summary has env, variant, sim_year, mean_cum_surv
surv_lookup_full <- egg_summary %>%
  arrange(env, variant, sim_year) %>%
  group_by(env, variant) %>%
  summarise(surv_vec = list(mean_cum_surv), .groups="drop") %>%
  mutate(key = paste(env, variant, sep = "_")) %>%
  select(key, surv_vec) %>%
  deframe()



# factory to make an optim‑ready function for each variant
make_calib_fn <- function(variant_name) {
  force(variant_name)
  function(par) {
    P_tmp     <- base_P
    P_tmp$S_rear <- par[1]
    SAR_vec   <- rep(par[2], n_calib)
    
    # --- here’s the slice to exactly 14 years ---
    full_ser <- surv_lookup_by_variant[[variant_name]]
    surv_ser <- full_ser[1:n_calib]  
    
    sim_out <- simulate_variant(
      surv_vec = surv_ser,
      P        = P_tmp,
      years    = n_calib,
      S_init   = S_seed,
      SAR_vec  = SAR_vec
    )
    pred <- sim_out$spawners[fit_idx]
    obs  <- obs_spawners[fit_idx]
    sum((pred - obs)^2)
  }
}

variant_names <- names(surv_lookup_by_variant)
init_par      <- c(base_P$S_rear, base_P$SAR_mean)

calib_results <- map_dfr(variant_names, function(varnm) {
  opt <- optim(
    par    = init_par,
    fn     = make_calib_fn(varnm),
    method = "L-BFGS-B",
    lower  = c(0, 0),
    upper  = c(1.0, 1)
  )
  tibble(
    variant   = varnm,
    S_rear    = opt$par[1],
    SAR_mean  = opt$par[2],
    obj_value = opt$value
  )
})

print(calib_results)


# a named list of base_P for each variant
base_P_list <- calib_results %>%
  split(.$variant) %>%
  map(~ {
    P <- base_P
    P$S_rear   <- .x$S_rear
    P$SAR_mean <- .x$SAR_mean
    P
  })

# ────────────────────────────────────────────────────────────────────────────
# STOCHASTIC SAR OPTIONS  ----------------------------------------------------
# ────────────────────────────────────────────────────────────────────────────

# 1) Pick your model and its parameters:
stoch_SAR_opts <- list(
  model       = "normal",       # "normal" | "lognormal" | "beta" | "gamma"
  mean        = base_P$SAR_mean,
  sd          = base_P$SAR_sd,
  
  # only used for beta/gamma:
  shape1      = 2,              
  shape2      = 5,              
  
  # timing:
  timing      = "pulse",        # "all" | "block" | "pulse"
  block_years = 20:30,          # if timing == "block", these years get special treatment
  pulse_years = c(10, 15, 20, 25, 30, 35, 40),  # if timing == "pulse", these are the only stochastic years
  pulse_sd    = 0.002           # magnitude of pulses (if desired)
)


# ────────────────────────────────────────────────────────────────────────────
# FORECAST (2011–2060) with stochastic SAR -----------------------------------
# ────────────────────────────────────────────────────────────────────────────

use_stochastic_SAR <- TRUE

SAR_by_year <- if (use_stochastic_SAR) {
  generate_SAR_vec(n_sim, stoch_SAR_opts)
} else {
  rep(base_P$SAR_mean, n_sim)
}

# 9) POSTPROCESS RESULTS ------------------------------------------------------
results_full_altvar <- map_dfr(names(surv_lookup_full), function(key) {
  message("Processing key: ", key)
  
  parts <- strsplit(key, "_")[[1]]
  if (length(parts) < 2) stop("Unexpected key format: ", key)
  
  env_nm <- parts[1]
  var_nm <- paste(parts[-1], collapse = "_")
  
  surv_vec <- surv_lookup_full[[key]]
  P        <- base_P_list[[var_nm]]
  
  if (is.null(surv_vec)) stop("surv_vec is NULL for key: ", key)
  if (is.null(P)) stop("P is NULL for variant: ", var_nm)
  
  result <- tryCatch({
    sim <- simulate_variant(
      surv_vec = surv_vec,
      P        = P,
      years    = n_sim,
      S_init   = S_seed,
      SAR_vec  = SAR_by_year
    )
    
    lag_df <- tibble(
      lag_p3 = P$lag_probs[1],
      lag_p4 = P$lag_probs[2],
      lag_p5 = P$lag_probs[3]
    )
    
    sim %>%
      mutate(
        variant_full = key,
        env = env_nm,
        variant = var_nm,
        !!!as.list(P[map_lgl(P, ~ length(.x) == 1)])
      ) %>%
      bind_cols(lag_df)
  }, error = function(e) {
    message("Error in key: ", key, " — ", conditionMessage(e))
    tibble()
  })
  
  result
})


# ─── 4) Save the results for Shiny to load
saveRDS(egg_summary,      "SalmonCountR/egg_summary.rds")
saveRDS(surv_lookup_full, "SalmonCountR/surv_lookup_full.rds")
saveRDS(base_P_list,      "SalmonCountR/base_P_list.rds")
saveRDS(S_seed,           "SalmonCountR/S_seed.rds")
saveRDS(stoch_SAR_opts,   "SalmonCountR/stoch_SAR_opts.rds")
