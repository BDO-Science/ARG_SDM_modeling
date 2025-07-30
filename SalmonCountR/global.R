# global.R

library(tidyverse); library(lubridate)
source("functions.R")

# load precomputed data (fast)
env_ext_list     <- readRDS("app_data/env_ext_list.rds")
df_all           <- readRDS("app_data/df_all.rds")
egg_summary      <- readRDS("app_data/egg_summary.rds")
surv_lookup_full <- readRDS("app_data/surv_lookup_full.rds")
base_P_list      <- readRDS("app_data/base_P_list.rds")
S_seed           <- readRDS("app_data/S_seed.rds")
stoch_SAR_opts   <- readRDS("app_data/stoch_SAR_opts.rds")
load("app_data/american_river_instream.rda")
spawn_dates_vec    <- readRDS("app_data/spawn_dates_vec.rds")

instream <- american_river_instream
n_sim      <- 100
real_years <- 2011:2024  # or whatever years you actually need
sim_years_full <- real_years[1] + seq(0, n_sim - 1)   # e.g. 2011:2060

# 2) compute redd capacity per flow
instream <- instream %>%
  mutate(
    K_spawners = FR_spawn_wua / 9.29  # m² per redd
  )

# 3) a helper to interpolate your K_spawners for any daily flow series
get_K_spawners <- function(flow_vec, lookup = instream) {
  # flow_vec: numeric vector of daily (or annual) CFS you want to match
  # returns a numeric vector same length as flow_vec
  approx(
    x    = lookup$flow_cfs,
    y    = lookup$K_spawners,
    xout = flow_vec,
    rule = 2
  )$y
}

# any constants UI needs
n_sim      <- 50  # or hard‐code 50
env_levels <- as.character(sort(as.numeric(unique(egg_summary$env))))
