# global.R

library(tidyverse); library(lubridate); library(here)
source("functions.R")

# ---- Scenario engine (the "Add a Year" tab) ---------------------------------
# Loaded defensively: if either the engine or the spawn-timing model is absent
# the rest of the app must still start, and the tab reports the problem rather
# than erroring on launch.
scenario_engine_loaded <- FALSE
spawn_timing_model     <- NULL
try({
  source(here("SalmonCountR", "scenario_engine.R"))
  stm <- here("SalmonCountR", "app_data", "spawn_timing_model.rds")
  if (file.exists(stm)) {
    spawn_timing_model     <- readRDS(stm)
    scenario_engine_loaded <- TRUE
  } else {
    warning("app_data/spawn_timing_model.rds not found; the 'Add a Year' tab will be disabled. ",
            "Run analysis/build_spawn_timing_model.R to create it.")
  }
}, silent = TRUE)

# load precomputed data (fast)
# NOTE: This script now assumes that the 'env_ext_list.rds' file in the
# 'app_data' directory contains the real temperature data for all alternatives.
env_ext_list       <- readRDS(here("SalmonCountR", "app_data","env_ext_list.rds"))
df_all_orig        <- readRDS(here("SalmonCountR", "app_data", "df_all.rds"))
egg_summary_orig   <- readRDS(here("SalmonCountR", "app_data", "egg_summary.rds"))
surv_lookup_full   <- readRDS(here("SalmonCountR", "app_data", "surv_lookup_full.rds"))
base_P_list        <- readRDS(here("SalmonCountR", "app_data", "base_P_list.rds"))
base_P             <- readRDS(here("SalmonCountR", "app_data", "base_P.rds"))
S_seed_calib       <- readRDS(here("SalmonCountR", "app_data", "S_seed_calib.rds"))
S_seed_fore_list   <- readRDS(here("SalmonCountR", "app_data", "S_seed_fore_list.rds"))
stoch_SAR_opts     <- readRDS(here("SalmonCountR", "app_data", "stoch_SAR_opts.rds"))
sim_years          <- readRDS(here("SalmonCountR", "app_data", "sim_years.rds"))
spawn_dates_by_alt <- readRDS(here("SalmonCountR", "app_data","spawn_dates_by_alt.rds"))
instream <- readRDS(here("SalmonCountR", "app_data","american_river_instream.rds"))
steelhead_metrics <- readRDS(here("SalmonCountR", "app_data", "steelhead_metrics.rds"))
swing_ranges <- readRDS(here("SalmonCountR", "app_data", "swing_ranges.rds"))
results_full <- readRDS(here("SalmonCountR", "app_data", "results_full.rds"))

real_years <- 2011:2024 
n_sim      <- 114
sim_years_full <- real_years[1] + seq(0, n_sim - 1)   # e.g. 2011:2060

# ---- Redd Capacity vs Flow (American River) ----------------------------------
instream <- instream %>%
  mutate(K_spawners = FR_spawn_wua / 9.29)

get_K_spawners <- function(flow_vec) {
  interp_fun <- approxfun(instream$flow_cfs, instream$K_spawners, rule = 2)
  interp_fun(flow_vec)
}

# ---- Working copies of the loaded data --------------------------------------
# Historical note: this file used to synthesise placeholder rows for a set of
# "new alternatives" whose model output did not exist yet, then merge them into
# df_all / egg_summary. That block was removed because it had stopped doing
# anything useful and had become actively misleading:
#   * df_all_orig has columns env/Date/site/temp/alt -- no year, variant or
#     sar_name -- so the df_all placeholder branch never fired and df_all was
#     left as an empty tibble.
#   * egg_summary_orig DOES have a `variant` column but no `alternative`
#     column, so the merge kept only the synthetic rows and silently dropped
#     the real data. `egg_summary` was therefore 100% runif() noise.
#   * The factor-ordering pass keyed off an `alternative` column that none of
#     these objects has, so it was a no-op on every object including
#     env_ext_list.
# Nothing in app.R referenced df_all or egg_summary (only df_temp_2025, built
# from df_all_orig below), so no app behaviour changes. The 36 alternatives are
# real model output now and are addressed by `env` throughout.
df_all      <- df_all_orig
egg_summary <- egg_summary_orig

# ---- Pre-computed temperature explorer data (2025 only, climate-labelled) ----
# Doing year(Date) / month(Date) on the full df_all_orig on every render is
# the main cause of the 20-second lag in the Temperature Explorer tab.
# We do it once here at startup instead.
if (exists("df_all_orig") && is.data.frame(df_all_orig)) {
  df_temp_2025 <- df_all_orig %>%
    filter(year(Date) == 2025) %>%
    mutate(
      month_num = month(Date),
      climate = case_when(
        env %in% as.character(1:9)   ~ "2011",
        env %in% as.character(10:18) ~ "2014",
        env %in% as.character(19:27) ~ "2017",
        env %in% as.character(28:36) ~ "2020"
      )
    )
}