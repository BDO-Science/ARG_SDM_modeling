# global.R

library(tidyverse); library(lubridate); library(here)
source("functions.R")

# load precomputed data (fast)
env_ext_list       <- readRDS(here("SalmonCountR","app_data","env_ext_list.rds"))
df_all             <- readRDS(here("SalmonCountR","app_data", "df_all.rds"))
egg_summary        <- readRDS(here("SalmonCountR","app_data", "egg_summary.rds"))
surv_lookup_full   <- readRDS(here("SalmonCountR","app_data", "surv_lookup_full.rds"))
base_P_list        <- readRDS(here("SalmonCountR","app_data", "base_P_list.rds"))
#rear_surv_lookup <- readRDS(here("SalmonCountR","app_data", "rear_surv_lookup.rds"))

calib_results         <- readRDS(here("SalmonCountR","app_data","calib_results.rds"))
calib_pred_by_variant <- readRDS(here("SalmonCountR","app_data","calib_pred_by_variant.rds"))

# calibration seed
S_seed_calib       <- readRDS(here("SalmonCountR","app_data", "S_seed_calib.rds"))
# forecast seeds
S_seed_fore_list   <- readRDS(here("SalmonCountR","app_data", "S_seed_fore_list.rds"))

stoch_SAR_opts     <- readRDS(here("SalmonCountR","app_data", "stoch_SAR_opts.rds"))

# full‐horizon year + date vectors
sim_years          <- readRDS(here("SalmonCountR","app_data", "sim_years.rds"))

spawn_dates_by_env <- readRDS(here("SalmonCountR","app_data","spawn_dates_by_env.rds"))

# legacy .rda data
load(here("SalmonCountR","app_data", "american_river_instream.rda"))

instream <- american_river_instream
real_years <- 2011:2024 
n_sim      <- 114
sim_years_full <- real_years[1] + seq(0, n_sim - 1)   # e.g. 2011:2060

# 2) compute redd capacity per flow
instream <- instream %>%
  mutate(
    K_spawners = FR_spawn_wua / 9.29  # m² per redd (SIT DSM)
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
env_levels <- as.character(sort(as.numeric(unique(egg_summary$env))))
