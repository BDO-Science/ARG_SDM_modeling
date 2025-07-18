# 0) CLEAN & LIBRARIES
rm(list=ls())
library(readxl)      # for excel
library(dataRetrieval) # for readNWISdata()
library(tidyverse)
library(lubridate)   # dates
library(progressr)   # progress bars
library(furrr)       # parallel map if you want
library(scales)
library(ggridges)
handlers("txtprogressbar")

# 1) PARAMETERS
obs_start   <- as.Date("2011-09-01")
obs_end     <- as.Date("2024-05-31")
last_wy     <- 2060                       # final water‐year to simulate
sim_end     <- as.Date(sprintf("%04d-08-31", last_wy + 1))

# 2) GET OBSERVED GAUGE TEMPS (Daily avg)
stations   <- c("11446980","11446500")
param_cd   <- "00010"

amer_obs <- map_df(stations, function(stn) {
  readNWISdata(
    sites       = stn,
    parameterCd = param_cd,
    service     = "uv",
    startDate   = obs_start,
    endDate     = obs_end
  ) %>%
    transmute(
      Date  = as.Date(dateTime),
      site  = recode(site_no,
                     "11446980" = "AveWatt",
                     "11446500" = "AveHazel"),
      temp  = X_00010_00000
    )
}) %>%
  group_by(Date, site) %>%
  summarise(temp = mean(temp, na.rm=TRUE), .groups="drop")

amer_obs <- amer_obs %>%
  group_by(Date, site) %>%
  summarise(temp = mean(temp, na.rm=TRUE), .groups="drop") %>%
  mutate(
    # Hazel floored at 7°C, all others at 8°C
    temp = if_else(site == "AveHazel",
                   pmax(if_else(is.na(temp)|is.nan(temp), 7, temp), 7),
                   pmax(if_else(is.na(temp)|is.nan(temp), 8, temp), 8))
  )

# 3) READ ALT FORECASTS & BUILD A DOY→temp TABLE FOR EACH ALT
xlsx_path <- "ARG_LAR_TempModeling_10-21-24.xlsx"
alts      <- excel_sheets(xlsx_path)[-1]  # assume first sheet is metadata

pred_by_doy <- map_df(alts, function(alt) {
  read_excel(xlsx_path, sheet=alt) %>%
    mutate(Date = as.Date(Date)) %>%
    pivot_longer(
      cols      = starts_with("Ave"),
      names_to  = "site",
      values_to = "temp_alt"
    ) %>%
    mutate(
      doy = yday(Date),
      alt = alt
    ) %>%
    select(alt, site, doy, temp_alt)
})

# 4) MAKE FUTURE DATE SKELETON
future_dates <- tibble(
  Date = seq(obs_end + 1, sim_end, by="day")
) %>% mutate(doy = yday(Date))

# ───────────────────────────────────────────────────────────────────────────────
# 1) Build a 14‑year DOY climatology of amer_obs
# ───────────────────────────────────────────────────────────────────────────────
clim14 <- amer_obs %>%
  mutate(doy = yday(Date)) %>%
  group_by(site, doy) %>%
  summarize(clim_temp = mean(temp, na.rm = TRUE), .groups = "drop")

# DOY thresholds for Oct 18 and Dec 31 (use a non‑leap year for consistency)
threshold_start <- yday(as.Date("2021-10-18"))  # 291
threshold_end   <- yday(as.Date("2021-12-31"))  # 365

# ───────────────────────────────────────────────────────────────────────────────
# 2) Rebuild env_ext_list with climatology‑fill
# ───────────────────────────────────────────────────────────────────────────────
env_ext_list <- map(alts, function(alt_nm) {
  # a) Observed block unchanged
  obs_block <- amer_obs %>%
    filter(Date <= obs_end) %>%
    mutate(alt = alt_nm)
  
  # b) Prepare the raw DOY→temp_alt pattern
  dedup_pattern <- pred_by_doy %>%
    filter(alt == alt_nm) %>%
    group_by(site, doy) %>%
    summarize(temp_alt = mean(temp_alt, na.rm = TRUE), .groups = "drop")
  
  # c) Build the full future skeleton
  future_skel <- future_dates %>%
    expand_grid(site = unique(dedup_pattern$site))
  
  # d) Join forecast + climatology, then apply the fill rule
  pred_block <- future_skel %>%
    left_join(dedup_pattern, by = c("doy","site")) %>%
    left_join(clim14,        by = c("doy","site")) %>%
    mutate(
      # fill logic as before
      temp_raw = if_else(
        doy >= threshold_start & doy <= threshold_end & !is.na(temp_alt),
        temp_alt,
        clim_temp
      ),
      # now floor by site
      temp = if_else(site == "AveHazel",
                     pmax(temp_raw, 7),
                     pmax(temp_raw, 8))
    ) %>%
    select(Date, site, temp) %>%
    mutate(alt = alt_nm)
  
  bind_rows(obs_block, pred_block) %>%
    arrange(Date, site)
}) %>% set_names(alts)

# 6) Now rebuild your site/Date lookup lists *after* the full horizons exist
site_temps_list <- map(env_ext_list, ~ split(.x$temp, .x$site))
site_dates_list <- map(env_ext_list, ~ split(.x$Date, .x$site))
date_idx_list   <- map(site_dates_list, function(dlist) {
  map(dlist, ~ setNames(seq_along(.x), as.character(.x)))
})


df_all <- bind_rows(env_ext_list, .id = "env") 

ggplot(df_all, aes(Date, temp, color = site)) +
  geom_line(size = 0.5, alpha = 0.8) +
  facet_wrap(~ env, ncol = 2, scales = "free_y") +
  labs(
    title = "Observed + Predicted Temp by Alternative",
    x     = "Date",
    y     = "Temperature (°C)",
    color = "Site"
  ) +
  theme_minimal(base_size = 12)

test_date <- obs_end + 2500  # 2024‑06‑01
 
  map_df(alts, function(alt_nm) {
    env_ext_list[[alt_nm]] %>%
      filter(Date == test_date) %>%
      mutate(alt = alt_nm)
    }) %>%
  select(alt, site, temp) %>%
  arrange(alt, site)

##########################
#######TDM MODELS#########
##########################

# Egg development model: estimates days to egg hatch as 958 ATU / daily temp
egg_model <- function(T) 958 / T
# Fry development model: estimates days to emergence as 417 ATU / daily temp
hatch_model <- function(T) 417 / T

# TDM exponential model (two calibrations available)
tdm_exp <- function(temps, calib) {
  # Calibration parameters for each curve
  p <- list(
    WaterForum2020 = list(α = 3.40848e-11, β = 1.21122),
    SALMOD2006     = list(α = 1.475e-11,   β = 1.392)
  )[[calib]]
  # Exponential cumulative survival: exp(-Σ α * exp(β * T_i))
  exp(-sum(p$α * exp(p$β * temps)))
}

# Linear TDM model from Martin et al. 2017
# Survival = exp(-α × Σ max(T - β, 0))
tdm_lin_martin <- function(temps, α = 0.026, β = 12.14) {
  exp(-α * sum(pmax(temps - β, 0)))
}

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


# ───────────────────────────────────────────────────────────────────────────────
# 0) Read & prep your observed carcass data, assign brood_year by spawn_dt
# ───────────────────────────────────────────────────────────────────────────────
carcass_raw <- read.csv("carcassdet_1752789274_15.csv")

obs_df <- carcass_raw %>%
  filter(clip == "No Clip") %>%
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
n_sim      <- 50
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
plan(multisession, workers = 4)               # or multisession/multicore as you like

results_obs_fast <- 
  future_map_dfr(sim_years, function(sim_yr) {
    # all observed redds for this brood‐year
    red_this <- filter(sim_redds, sim_year == sim_yr)
    
    map_dfr(names(env_ext_list), function(env_nm) {
      date_idx   <- date_idx_list[[env_nm]]
      temps_list <- site_temps_list[[env_nm]]
      env_dates  <- env_ext_list[[env_nm]]$Date
      
      # 2) pull out month/day of each spawn
      spawn_m <- month(red_this$spawn_dt)
      spawn_d <- day( red_this$spawn_dt)
      
      rdrs <- make_date(
        year  = if_else(month(red_this$spawn_dt) >= 9L, sim_yr, sim_yr + 1L),
        month = month(red_this$spawn_dt),
        day   = day(  red_this$spawn_dt)
      )
      
      # 4) optionally filter only those that fall into our spawn window
      in_window <- (rdrs >= min(env_dates)) & (rdrs <= max(env_dates))
      rdrs <- rdrs[in_window]
      sites <- red_this$site[in_window]
      
      # 5) now do the exact same vapply loop you had:
      map_dfr(seq_len(nrow(tdm_defs)), function(i) {
        model   <- tdm_defs$model[i]
        calib   <- tdm_defs$calib[i]
        variant <- tdm_defs$variant[i]
        
        survs <- map_dbl(seq_along(rdrs), function(j) {
          site_j <- sites[j]
          rdr_j  <- rdrs[j]
          pos    <- date_idx_list[[env_nm]][[site_j]][as.character(rdr_j)]
          
          # bail if anything is off
          if (is.na(pos) || pos < 1)            return(NA_real_)
          temp0 <- site_temps_list[[env_nm]][[site_j]][pos]
          if (is.na(temp0) || temp0 <= 0)       return(NA_real_)
          
          d_egg   <- egg_model(temp0)
          d_hatch <- hatch_model(temp0)
          td      <- ceiling(d_egg + d_hatch)
          if (is.na(td) || td < 1)              return(NA_real_)
          
          slice   <- temps_list[[site_j]][pos + seq_len(td) - 1]
          if (model == "exp") tdm_exp(slice, calib) else tdm_lin_martin(slice)
        })
        
        tibble(
          sim_year      = sim_yr,
          env           = env_nm,
          variant       = variant,
          method        = "observed",
          mean_cum_surv = mean(survs, na.rm = TRUE)
        )
      })
    })
  },
  .options = furrr_options(seed = TRUE)
  )

print(results_obs_fast, n = Inf)

ggplot(results_obs_fast, aes(x = sim_year, y = mean_cum_surv, color = variant)) +
  geom_line(size = 1) +
  facet_wrap(~ env, ncol = 2) +
  scale_x_continuous(breaks = pretty(results_obs_fast$sim_year, 10)) +
  labs(
    title = "Observed‐based TDM: Mean Egg→Emergence Survival by Brood Year",
    x     = "Brood Year",
    y     = "Mean cumulative survival",
    color = "Variant"
  ) +
  theme_minimal(base_size = 14)

ggplot(results_obs_fast, aes(x = variant, y = mean_cum_surv)) +
  geom_boxplot(outlier.size = 1) +
  facet_wrap(~ env) +
  labs(
    title = "Distribution of Brood‐Year Survival by Variant & Env",
    x     = "TDM Variant",
    y     = "Mean cumulative survival"
  ) +
  theme_minimal(base_size = 14)


####################################
#FISH SPAWNING BEFORE 11-15 EACH YEAR
#####################################
results_obs_early <- future_map_dfr(sim_years, function(sim_yr) {
  # 1) pull only the redds for this brood_year AND before Nov 15
  red_this <- sim_redds %>%
    filter(sim_year == sim_yr) %>%
    # keep spawn_dt ≤ yyyy-11-15
    filter(spawn_dt <= as.Date(sprintf("%04d-11-15", sim_yr)))
  
  # if none, skip
  if (nrow(red_this)==0) {
    return(tibble(
      sim_year      = sim_yr,
      env           = NA_character_,
      variant       = NA_character_,
      method        = "early",
      mean_cum_surv = NA_real_
    ))
  }
  
  # 2) proceed exactly as in results_obs_fast, but on this red_this subset
  map_dfr(names(env_ext_list), function(env_nm) {
    date_idx   <- date_idx_list[[env_nm]]
    temps_list <- site_temps_list[[env_nm]]
    env_dates  <- env_ext_list[[env_nm]]$Date
    
    rdrs <- make_date(
      year  = if_else(month(red_this$spawn_dt) >= 9L, sim_yr, sim_yr + 1L),
      month = month(red_this$spawn_dt),
      day   = day(  red_this$spawn_dt)
    )
    
    in_window <- (rdrs >= min(env_dates)) & (rdrs <= max(env_dates))
    rdrs  <- rdrs[in_window]
    sites <- red_this$site[in_window]
    
    map_dfr(seq_len(nrow(tdm_defs)), function(i) {
      model   <- tdm_defs$model[i]
      calib   <- tdm_defs$calib[i]
      variant <- tdm_defs$variant[i]
      
      survs <- map_dbl(seq_along(rdrs), function(j) {
        site_j <- sites[j]
        rdr_j  <- rdrs[j]
        pos    <- date_idx[[site_j]][as.character(rdr_j)]
        if (!isTRUE(pos >= 1))            return(NA_real_)
        temp0 <- temps_list[[site_j]][pos]
        if (!isTRUE(temp0 > 0))           return(NA_real_)
        d_egg   <- egg_model(temp0)
        d_hatch <- hatch_model(temp0)
        td      <- ceiling(d_egg + d_hatch)
        if (!isTRUE(td >= 1))             return(NA_real_)
        slice   <- temps_list[[site_j]][pos + seq_len(td) - 1]
        if (model=="exp") tdm_exp(slice, calib) else tdm_lin_martin(slice)
      })
      
      tibble(
        sim_year      = sim_yr,
        env           = env_nm,
        variant       = variant,
        method        = "early",
        mean_cum_surv = mean(survs, na.rm=TRUE)
      )
    })
  })
}, .options = furrr_options(seed = TRUE))

print(results_obs_early, n = Inf)

ggplot(results_obs_early, aes(x = sim_year, y = mean_cum_surv, color = variant)) +
  geom_line(size = 1) +
  facet_wrap(~ env, ncol = 2) +
  scale_x_continuous(breaks = pretty(results_obs_early$sim_year, 10)) +
  labs(
    title = "Observed‐based TDM: Mean Egg→Emergence Survival by Brood Year",
    x     = "Brood Year",
    y     = "Mean cumulative survival",
    color = "Variant"
  ) +
  theme_minimal(base_size = 14)

ggplot(results_obs_early, aes(x = variant, y = mean_cum_surv)) +
  geom_boxplot(outlier.size = 1) +
  facet_wrap(~ env) +
  labs(
    title = "Distribution of Brood‐Year Survival by Variant & Env",
    x     = "TDM Variant",
    y     = "Mean cumulative survival"
  ) +
  theme_minimal(base_size = 14)

# LIFE CYCLE MODEL: CALIBRATION & FORECAST 

# LOAD TDM OUTPUTS ----------------------------------------------------------
egg_summary <- results_obs_fast
# egg_summary must have: env, variant, mean_cum_surv, year

simulate_variant <- function(
    surv_vec,
    P,
    years,
    S_init,
    SAR_vec = rep(P$SAR_mean, years)
) {
  # — ensure surv_vec is exactly `years` long —
  if (length(surv_vec) == 1) {
    surv_vec <- rep(surv_vec, years)
  } else if (length(surv_vec) < years) {
    surv_vec <- c(surv_vec, rep(tail(surv_vec,1), years - length(surv_vec)))
  } else if (length(surv_vec) > years) {
    surv_vec <- surv_vec[1:years]
  }
  
  # seed years guard
  n_init <- length(S_init)
  if (n_init > years) {
    stop("length(S_init) = ", n_init,
         " > years = ", years,
         ". Please shorten S_init or increase years.")
  }
  seed_len <- min(n_init, years)
  
  # preallocate
  S        <- numeric(years)
  reared   <- numeric(years)
  dd_vec   <- numeric(years)
  fry_dd   <- numeric(years)
  SAR_used <- numeric(years)
  
  # seed
  S[1:seed_len] <- S_init[1:seed_len]
  
  # loop
  for (t in seq_len(years)) {
    redds      <- S[t] * P$female_fraction
    eggs       <- redds * P$fec
    dd         <- P$S0 / (1 + redds / P$K_spawners)
    fry_dd[t]  <- eggs * surv_vec[t] * dd
    dd_vec[t]  <- dd
    reared[t]  <- fry_dd[t] * P$S_rear
    SAR_t      <- SAR_vec[t]
    SAR_used[t]<- SAR_t
    # age‐structured returns
    for (age in 3:5) {
      ry <- t + age
      if (ry <= years) {
        S[ry] <- S[ry] + reared[t] * SAR_t * P$lag_probs[as.character(age)]
      }
    }
  }
  
  # return tibble
  tibble(
    year      = seq_len(years),
    spawners  = S,
    dd        = dd_vec,
    fry_dd    = fry_dd,
    eff_surv  = surv_vec * dd_vec,
    SAR_used  = SAR_used
  )
}

# 5) SET PARAMETER LIST --------------------------------------------------------
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

# 6) LOAD OBSERVED ESCAPEMENT --------------------------------------------------
esc_obs <- read_csv(
  "grandtab_1752793045_337.csv",
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
n_total      <- 50

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

# 2) A generator function
generate_SAR_vec <- function(n_years, opts) {
  # draw a full-length series
  vec <- switch(opts$model,
                normal    = rnorm(n_years, opts$mean, opts$sd),
                lognormal = {
                  # match lognormal moments to mean/sd
                  mu    <- log(opts$mean^2 / sqrt(opts$sd^2 + opts$mean^2))
                  sigma <- sqrt(log(1 + opts$sd^2 / opts$mean^2))
                  rlnorm(n_years, mu, sigma)
                },
                beta      = rbeta(n_years, opts$shape1, opts$shape2),
                gamma     = rgamma(n_years, shape = opts$shape1,
                                   scale = opts$mean / opts$shape1)
  )
  # now apply timing adjustments
  if (opts$timing == "block") {
    # zero‑out (or reset to mean) outside the block
    vec[-opts$block_years] <- opts$mean
  }
  if (opts$timing == "pulse") {
    # pulses in just those years
    vec[] <- opts$mean
    vec[opts$pulse_years] <- rnorm(length(opts$pulse_years),
                                   opts$mean, opts$pulse_sd)
  }
  # no timing filter for "all"
  pmax(vec, 0)  
}

# ────────────────────────────────────────────────────────────────────────────
# FORECAST (2011–2060) with stochastic SAR -----------------------------------
# ────────────────────────────────────────────────────────────────────────────

n_total <- 50
use_stochastic_SAR <- TRUE

SAR_by_year <- if (use_stochastic_SAR) {
  generate_SAR_vec(n_total, stoch_SAR_opts)
} else {
  rep(base_P$SAR_mean, n_total)
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
      years    = n_total,
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

pct_results <- results_full_altvar %>%
  group_by(env, year) %>%
  mutate(
    max_spawners = max(spawners, na.rm = TRUE),
    pct_diff     = (spawners - max_spawners) / max_spawners * 100
  ) %>%
  ungroup() %>%
  select(env, alt_variant = variant_full, variant, year, spawners, pct_diff)



last_year <- max(pct_results$year)
cat("% differences from best variant in year", last_year, ":\n")
pct_results %>%
  filter(year == last_year) %>%
  arrange(env, desc(pct_diff)) %>%
  print(n = Inf)

# 10) PLOTTING -----------------------------------------------------------------
# define the exact order you want
alt_levels <- paste0("Alt ", 1:10)
tdm_levels <- c("exp_SM","exp_WF","lin_Martin") 

results2 <- results_full_altvar %>%
  mutate(
    # 1) env → fill aesthetic
    env = factor(env,
                 levels = as.character(1:10),    # original codes 1–10
                 labels = alt_levels),           # prettied-up labels
    
    # 2) variant → x-axis aesthetic
    variant = factor(variant,
                     levels = tdm_levels)      # in the order you like
  )

# now re-calc max
max_spawn <- max(results2$spawners, na.rm = TRUE)
last_year <- max(results2$year)

# A) Time series
p1 <- ggplot(results2, aes(x = year, y = spawners, color = variant)) +
  geom_line(size = 1) +
  facet_wrap(~ env) +
  scale_y_continuous(limits = c(0, max_spawn), breaks = pretty_breaks(10), labels = comma) +
  labs(x = "Simulation Year", y = "Spawners", color = "TDM Variant") +
  theme_minimal() + theme(legend.position = "bottom")

print(p1)


ggplot(results2, aes(x = year, y = variant, fill = spawners)) +
  geom_tile() +
  facet_wrap(~env, ncol = 2) +
  scale_fill_viridis_c(option = "magma", name = "Spawners") +
  labs(x = "Year", y = "TDM Variant") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7))

final_decade <- results2 %>% filter(year >= (max(year)-9))

# 1) p2: end‑of‑simulation barplot
p2 <- ggplot(
  filter(results2, year == last_year),
  aes(x = variant, y = spawners, fill = env)
) +
  geom_col(position = position_dodge(width = 0.8)) +
  scale_y_continuous(
    limits = c(0, 24000),                     # <-- clamp y to 0–40k
    breaks  = pretty_breaks(10),
    labels  = comma
  ) +
  scale_fill_discrete(breaks = alt_levels) +
  labs(
    title = paste("End‑of‑Simulation Spawners (Year", last_year, ")"),
    x     = "TDM Variant",
    y     = "Spawners",
    fill  = "Alternative"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p2)

# 2) facet_grid line‐plot
ggplot(results2, aes(x = year, y = spawners)) +
  geom_line(color = "steelblue", size = 0.8) +
  facet_grid(variant ~ env, scales = "free_y") +
  scale_y_continuous(
    limits = c(0, 40000),                     # <-- clamp here as well
    labels  = comma
  ) +
  labs(x = "Year", y = "Spawners") +
  theme_minimal() +
  theme(
    strip.text   = element_text(size = 6),
    axis.text.y  = element_text(size = 6)
  )

# 3) boxplot of final decade
ggplot(final_decade, aes(x = variant, y = spawners)) +
  geom_boxplot() +
  facet_wrap(~env, nrow = 2) +
  scale_y_continuous(
    limits = c(0, 22000),                     # <-- and here
    labels  = comma
  ) +
  labs(
    title = "Distribution of Spawners in Final 10 Years by Variant",
    x     = "TDM Variant",
    y     = "Spawners"
  ) +
  theme_minimal()

# 4) density‑ridge: spawners is on the x‑axis, so clamp x instead
ggplot(results2, aes(x = spawners, y = variant, fill = variant)) +
  geom_density_ridges(alpha = 0.7, scale = 1.2) +
  facet_wrap(~env, ncol = 2) +
  scale_x_continuous(
    limits = c(0, 40000),                     # <-- clamp the x‑axis
    labels  = comma
  ) +
  labs(x = "Spawners", y = "TDM Variant") +
  theme_minimal() +
  theme(legend.position = "none")
