# ───────────────────────────────────────────────────────────────────────────────
# spatial_tdm_multi_year.R - Thoroughly commented version
#   Reads multiple sheets of daily temperature data,
#   pads missing dates, generates redd distributions,
#   assigns redds to sites, runs TDM variants,
#   and outputs mean survival by year × environment × variant.
# ───────────────────────────────────────────────────────────────────────────────

# 0) CLEAN WORKSPACE & LOAD REQUIRED LIBRARIES
# Remove all objects from the global environment to start fresh
rm(list = ls())

# Load libraries for reading Excel, data manipulation, dates, plotting, and parallel processing
library(readxl)    # read Excel files
library(dplyr)     # data manipulation (filter, mutate, etc.)
library(tidyr)     # data reshaping (pivot_longer)
library(purrr)     # functional mapping (map, map_dfr)
library(lubridate) # date handling (ymd, days)
library(ggplot2)   # plotting
library(progressr) # progress bars for long operations
library(furrr)     # parallel mapping
handlers("txtprogressbar")  # use text progress bar in Rscript

# 1) USER PARAMETERS
# Path to Excel file containing 10 sheets of daily temperature time series
xlsx_path      <- "ARG_LAR_TempModeling_10-21-24.xlsx"
# Starting water year for simulation
watershed_year <- 2024
# Number of years to simulate (e.g., 50 years)
n_years        <- 50
# Number of redds (spawning sites) to generate per simulation year
n_redds        <- 1000
# Standard deviation for redd timing (in days)
sigma_days     <- 20
# Standard deviation for peak spawn date shift (in days)
shift_sd       <- 14
# Valid site names matching columns in Excel (e.g., "AveHazel", "AveWatt")
valid_sites    <- c("AveHazel", "AveWatt")
# Probability of assigning a redd to each site (70% AveHazel, 30% AveWatt)
site_probs     <- c(0.7, 0.3)

# 2) READ & TIDY ALL ENVIRONMENT SHEETS
# Get sheet names (excluding the first sheet if it's metadata)
data_sheets <- excel_sheets(xlsx_path)[-1]

# Read each sheet, reshape to long format, pad missing dates, and floor temperatures
env_temps_ext <- map_dfr(data_sheets, function(sh) {
  # Read data: Date + columns starting with "Ave"
  df <- read_excel(xlsx_path, sheet = sh) %>%
    select(Date, starts_with("Ave")) %>%
    mutate(Date = as.Date(Date)) %>%
    pivot_longer(cols = starts_with("Ave"),
                 names_to = "site", values_to = "temp")
  
  # Determine first and last dates in this sheet
  fd   <- min(df$Date)
  ld   <- max(df$Date)
  yr0  <- year(fd)
  
  # Create padding before and after data window:
  #   Pre: Sept 1 of year0 up to day before first date
  #   Post: day after last date up to May 31 of year0+1
  pre  <- seq.Date(ymd(paste0(yr0, "-09-01")), fd - 1, by = "day")
  post <- seq.Date(ld + 1, ymd(paste0(yr0 + 1, "-05-31")), by = "day")
  
  # Get all site names present
  all_sites <- unique(df$site)
  
  # Combine padded rows with original data,
  # setting default temps: 17°C pre, 9°C post,
  # then enforce minimum temp of 8°C across entire period.
  bind_rows(
    expand_grid(Date = pre,  site = all_sites) %>% mutate(temp = 17),
    df,
    expand_grid(Date = post, site = all_sites) %>% mutate(temp = 9)
  ) %>%
    mutate(
      temp = pmax(temp, 8),  # floor temperature at 8°C
      env  = sh             # tag environment with sheet name
    )
})

# Split combined data by environment for easy lookup
env_ext_list <- split(env_temps_ext, env_temps_ext$env)

# 3) PRE-CALCULATIONS: Spawn window & TDM FUNCTIONS
# Define the sequence of potential spawn dates: Sept 15 (WY) to Apr 30 (WY+1)
dates     <- seq(ymd(sprintf("%d-09-15", watershed_year)),
                 ymd(sprintf("%d-04-30", watershed_year + 1)), by = "day")
# Compute base peak spawn day offset relative to first date
base_peak <- as.numeric(ymd(sprintf("%d-12-01", watershed_year)) - dates[1])

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

# 4) PREPARE LOOKUP STRUCTURES FOR EFFICIENT INDEXING
# Split full temperature series by site for each environment
site_temps_list <- map(env_ext_list, ~ split(.x$temp, .x$site))
# Also split dates by site for indexing
site_dates_list <- map(env_ext_list, ~ split(.x$Date,  .x$site))

# Build a lookup: for each env/site, a named vector mapping date strings to positions
date_idx_list <- map(site_dates_list, function(dlist) {
  map(dlist, ~ setNames(seq_along(.x), as.character(.x)))
})

# 5) MULTI-YEAR SIMULATION LOOP (PARALLEL)
# Set up parallel plan
plan(multisession)

# Run simulations for each of n_years, returning a data frame of mean survival
results_multi_fast <- future_map_dfr(
  seq_len(n_years),       # loop index = simulation year
  .progress = TRUE,
  .options  = furrr_options(seed = TRUE),
  function(yr_idx) {
    # Random shift of spawn peak around base_peak
    peak_shift <- rnorm(1, 0, shift_sd)
    spawn_mean <- base_peak + peak_shift
    
    # Generate redd spawn dates: normally distributed around spawn_mean
    offs    <- round(rnorm(n_redds, spawn_mean, sigma_days))
    # Clamp offsets to valid range
    offs    <- pmin(pmax(offs, 0), length(dates) - 1)
    redds_j <- dates[1] + offs  # actual spawn dates
    
    # For each environment, compute variant survivals
    map_dfr(names(env_ext_list), function(env_nm) {
      # Lookup structures for this environment
      date_idx   <- date_idx_list[[env_nm]]
      temps_list <- site_temps_list[[env_nm]]
      
      # Clamp redds to this env's available date range
      minD <- min(env_ext_list[[env_nm]]$Date)
      maxD <- max(env_ext_list[[env_nm]]$Date)
      rdrs <- pmin(pmax(redds_j, minD), maxD)
      
      # Assign redds to sites based on probability 70/30
      site_assign <- sample(valid_sites, n_redds, TRUE, site_probs)
      
      # For each TDM variant, compute survival for all redds
      map_dfr(seq_len(nrow(tdm_defs)), function(i) {
        model   <- tdm_defs$model[i]
        calib   <- tdm_defs$calib[i]
        variant <- tdm_defs$variant[i]
        
        # Preallocate vector for survival values
        survs <- numeric(n_redds)
        
        # Loop through each redd
        for (j in seq_len(n_redds)) {
          this_site <- site_assign[j]
          idx_map   <- date_idx[[this_site]]
          # Find position of spawn date in site-specific date vector
          pos       <- idx_map[as.character(rdrs[j])]
          # Days to egg and fry development based on temperature at spawn
          d_egg   <- egg_model(temps_list[[this_site]][pos])
          d_hatch <- hatch_model(temps_list[[this_site]][pos])
          td      <- ceiling(d_egg + d_hatch)
          # Temperature slice covering the entire incubation period
          slice   <- temps_list[[this_site]][pos + seq_len(td) - 1]
          # Compute survival depending on model
          survs[j] <- if (model == "exp") {
            tdm_exp(slice, calib)
          } else {
            tdm_lin_martin(slice)
          }
        }
        
        # Return tibble with mean survival for this variant and environment
        tibble(
          year          = yr_idx,
          env           = env_nm,
          variant       = variant,
          peak_shift    = peak_shift,
          mean_cum_surv = mean(survs, na.rm = TRUE)
        )
      })
    })
  }
)

# Inspect first few rows of simulation results
print(head(results_multi_fast, 10))

# Optional: plot mean survival across years by variant and environment
ggplot(results_multi_fast, aes(x = year, y = mean_cum_surv, color = variant)) +
  geom_line() +
  facet_wrap(~ env) +
  labs(
    title = "Mean Egg → Emergence Survival over 50 Years",
    x     = "Simulation Year",
    y     = "Mean survival"
  ) +
  theme_minimal()

library(dplyr)
library(lubridate)
library(purrr)

# ───────────────────────────────────────────────────────────────────────────────
# 0) Read & prep your observed carcass data, assign brood_year by spawn_dt
# ───────────────────────────────────────────────────────────────────────────────
carcass_raw <- read.csv("carcassdet_1752789274_15.csv")

obs_df <- carcass_raw %>%
  #filter(ConditionCD == "F") %>%
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
    section %in% c("1a","1b") ~ "AveHazel",
    section %in% c("2","3")    ~ "AveWatt",
    TRUE                       ~ sample(valid_sites, n(), TRUE, site_probs)
  ))

# ───────────────────────────────────────────────────────────────────────────────
# 3) Plug into your TDM loop exactly as before, but now sim_year is brood_year
# ───────────────────────────────────────────────────────────────────────────────
library(furrr)
plan(multisession, workers = 4)               # or multisession/multicore as you like

results_obs_fast <- 
  future_map_dfr(sim_years, function(sim_yr) {
    # all observed redds for this brood‐year
    red_this <- filter(sim_redds, sim_year == sim_yr)
    
    map_dfr(names(env_ext_list), function(env_nm) {
      date_idx   <- date_idx_list[[env_nm]]
      temps_list <- site_temps_list[[env_nm]]
      env_dates  <- env_ext_list[[env_nm]]$Date
      
      # 1) figure out which water‐year we're using:
      start_yr <- year(min(env_dates))
      
      # 2) pull out month/day of each spawn
      spawn_m <- month(red_this$spawn_dt)
      spawn_d <- day( red_this$spawn_dt)
      
      # 3) re‐construct each spawn date *in* our water‐year:
      #    if month >= 9, it belongs to start_yr; else start_yr+1
      rdrs <- make_date(start_yr + (spawn_m < 9L), spawn_m, spawn_d)
      
      # 4) optionally filter only those that fall into our spawn window
      in_window <- (rdrs >= min(env_dates)) & (rdrs <= max(env_dates))
      rdrs <- rdrs[in_window]
      sites <- red_this$site[in_window]
      
      # 5) now do the exact same vapply loop you had:
      map_dfr(seq_len(nrow(tdm_defs)), function(i) {
        model   <- tdm_defs$model[i]
        calib   <- tdm_defs$calib[i]
        variant <- tdm_defs$variant[i]
        
        survs <- vapply(seq_along(rdrs), function(j) {
          pos     <- date_idx[[sites[j]]][ as.character(rdrs[j]) ]
          d_egg   <- egg_model(temps_list[[sites[j]]][pos])
          d_hatch <- hatch_model(temps_list[[sites[j]]][pos])
          td      <- ceiling(d_egg + d_hatch)
          slice   <- temps_list[[sites[j]]][ pos + seq_len(td) - 1 ]
          if (model == "exp") tdm_exp(slice, calib) else tdm_lin_martin(slice)
        }, numeric(1))
        
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
