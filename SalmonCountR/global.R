# global.R

library(tidyverse); library(lubridate)

# Bootstrap only: locate the app directory well enough to source years.R, which
# then owns path resolution for everything else (see arg_app_dir there for why
# here() was removed). This probe is deliberately the same test years.R uses --
# the one place the duplication is unavoidable, because we cannot call
# arg_app_path() before the file that defines it has been sourced.
.arg_boot_dir <- if (dir.exists("app_data")) {
  "."
} else if (dir.exists(file.path("SalmonCountR", "app_data"))) {
  "SalmonCountR"
} else {
  stop("Cannot locate the SalmonCountR app directory from '", getwd(), "'.",
       call. = FALSE)
}

source(file.path(.arg_boot_dir, "functions.R"))
source(file.path(.arg_boot_dir, "years.R"))   # registry -- read its header first

# Startup loads below are the DEFAULT analysis year, unchanged from before the
# year selector existed. They are what functions.R reads out of the global
# environment, so they stay. Anything the app itself reads now comes from a
# per-session bundle (see load_year_bundle) and follows the selector instead.

# load precomputed data (fast)
# NOTE: This script now assumes that the 'env_ext_list.rds' file in the
# 'app_data' directory contains the real temperature data for all alternatives.
env_ext_list       <- readRDS(arg_app_path("app_data", "env_ext_list.rds"))
df_all_orig        <- readRDS(arg_app_path("app_data", "df_all.rds"))
egg_summary_orig   <- readRDS(arg_app_path("app_data", "egg_summary.rds"))
surv_lookup_full   <- readRDS(arg_app_path("app_data", "surv_lookup_full.rds"))
base_P_list        <- readRDS(arg_app_path("app_data", "base_P_list.rds"))
base_P             <- readRDS(arg_app_path("app_data", "base_P.rds"))
S_seed_calib       <- readRDS(arg_app_path("app_data", "S_seed_calib.rds"))
S_seed_fore_list   <- readRDS(arg_app_path("app_data", "S_seed_fore_list.rds"))
stoch_SAR_opts     <- readRDS(arg_app_path("app_data", "stoch_SAR_opts.rds"))
sim_years          <- readRDS(arg_app_path("app_data", "sim_years.rds"))
spawn_dates_by_alt <- readRDS(arg_app_path("app_data", "spawn_dates_by_alt.rds"))
instream <- readRDS(arg_app_path("app_data", "american_river_instream.rds"))
steelhead_metrics <- readRDS(arg_app_path("app_data", "steelhead_metrics.rds"))
swing_ranges <- readRDS(arg_app_path("app_data", "swing_ranges.rds"))
results_full <- readRDS(arg_app_path("app_data", "results_full.rds"))

# precompute.R writes these two and app.R's Decision Support tab reads them, but
# nothing ever loaded them -- see BUGFIX note in app.R's performance_data_full().
swing_scenario_results     <- readRDS(arg_app_path("app_data", "swing_scenario_results.rds"))
steelhead_scenario_results <- readRDS(arg_app_path("app_data", "steelhead_scenario_results.rds"))

real_years <- 2011:2024
n_sim      <- 114
sim_years_full <- real_years[1] + seq(0, n_sim - 1)   # 2011:2124

# ---- Redd Capacity vs Flow (American River) ----------------------------------
instream <- instream %>%
  mutate(K_spawners = FR_spawn_wua / 9.29)

get_K_spawners <- function(flow_vec) {
  interp_fun <- approxfun(instream$flow_cfs, instream$K_spawners, rule = 2)
  interp_fun(flow_vec)
}

# ---- Default analysis-year bundle -------------------------------------------
# Everything the app reads, for the startup year. Loading it here means a cold
# start pays the cost once; the selector swaps to another year's bundle without
# touching these globals, so concurrent sessions on shinyapps.io cannot tread on
# each other.
#
# This also builds the Temperature Explorer's pre-filtered frame (previously
# df_temp_2025, computed here at startup because doing year(Date)/month(Date)
# over the full df_all_orig on every render was the main cause of the 20-second
# lag on that tab). It now filters on the bundle's first projection year rather
# than a literal 2025, which would return zero rows for any later deliverable.
ARG_BUNDLE_DEFAULT <- load_year_bundle(ARG_DEFAULT_YEAR)

if (is.null(ARG_BUNDLE_DEFAULT)) {
  stop(
    "Default analysis year '", ARG_DEFAULT_YEAR, "' is missing required files: ",
    paste(arg_year_missing(ARG_DEFAULT_YEAR), collapse = ", "),
    call. = FALSE
  )
}