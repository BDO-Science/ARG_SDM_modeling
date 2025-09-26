# global.R

library(tidyverse); library(lubridate); library(here)
source("functions.R")

# load precomputed data (fast)
# NOTE: This script now assumes that the 'env_ext_list.rds' file in the
# 'app_data' directory contains the real temperature data for all alternatives.
env_ext_list       <- readRDS(here("SalmonCountR", "app_data","env_ext_list.rds"))
df_all_orig        <- readRDS(here("SalmonCountR", "app_data", "df_all.rds"))
egg_summary_orig   <- readRDS(here("SalmonCountR", "app_data", "egg_summary.rds"))
surv_lookup_full   <- readRDS(here("SalmonCountR", "app_data", "surv_lookup_full.rds"))
base_P_list        <- readRDS(here("SalmonCountR", "app_data", "base_P_list.rds"))
base_P             <- readRDS(here("SalmonCountR", "app_data", "base_P.rds"))
calib_results         <- readRDS(here("SalmonCountR", "app_data","calib_results.rds"))
calib_pred_by_variant <- readRDS(here("SalmonCountR", "app_data","calib_pred_by_variant.rds"))
S_seed_calib       <- readRDS(here("SalmonCountR", "app_data", "S_seed_calib.rds"))
S_seed_fore_list   <- readRDS(here("SalmonCountR", "app_data", "S_seed_fore_list.rds"))
stoch_SAR_opts     <- readRDS(here("SalmonCountR", "app_data", "stoch_SAR_opts.rds"))
sim_years          <- readRDS(here("SalmonCountR", "app_data", "sim_years.rds"))
spawn_dates_by_alt <- readRDS(here("SalmonCountR", "app_data","spawn_dates_by_alt.rds"))
instream <- readRDS(here("SalmonCountR", "app_data","american_river_instream.rds"))

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

# ---- ADD PLACEHOLDERS FOR NEW ALTERNATIVES ----
new_alt_names <- c(
  "NB: No Bypass",
  "PB1: 125-250cfs Oct15-Nov14",
  "PB2: 250-500cfs Oct15-Nov30",
  "PB3: 250-500cfs Oct21-Nov14",
  "PB4: 250-500cfs Oct21-Nov21",
  "PB5: 500-250cfs Oct28-Nov21",
  "PB6: 100-500cfs Step-up Oct1-Nov14"
)

hydro_years_placeholder <- c(2011, 2014, 2017, 2020)
placeholder_df_all <- tibble()
placeholder_egg_summary <- tibble()

# --- Safely generate placeholder data ---
# Only create placeholders if the original data has the required structure
required_cols_df <- c("year", "variant", "sar_name")
if (exists("df_all_orig") && is.data.frame(df_all_orig) && all(required_cols_df %in% names(df_all_orig))) {
  all_sim_years <- unique(df_all_orig$year)
  variants <- unique(df_all_orig$variant)
  sar_options <- unique(df_all_orig$sar_name)
  
  placeholder_df_all_list <- list()
  if (length(all_sim_years) > 0 && length(variants) > 0 && length(sar_options) > 0) {
    for (alt in new_alt_names) {
      for (yr in all_sim_years) {
        for (var in variants) {
          for (sar in sar_options) {
            set.seed(which(new_alt_names == alt) * yr + which(variants == var))
            placeholder_df_all_list[[length(placeholder_df_all_list) + 1]] <-
              tibble(
                alternative = alt,
                year = yr,
                spawners = round(runif(1, 15000, 40000)),
                variant = var,
                sar_name = sar
              )
          }
        }
      }
    }
    placeholder_df_all <- bind_rows(placeholder_df_all_list)
  }
}

required_cols_egg <- c("variant")
if (exists("egg_summary_orig") && is.data.frame(egg_summary_orig) && all(required_cols_egg %in% names(egg_summary_orig))) {
  variants <- unique(egg_summary_orig$variant)
  placeholder_egg_summary_list <- list()
  if (length(variants) > 0) {
    for (alt in new_alt_names) {
      for (yr in hydro_years_placeholder) {
        for (var in variants) {
          set.seed(which(new_alt_names == alt) * yr * 2 + which(variants == var))
          placeholder_egg_summary_list[[length(placeholder_egg_summary_list) + 1]] <-
            tibble(
              alternative = alt,
              year = yr,
              days_lt_18.3C = round(runif(1, 50, 61)),
              variant = var
            )
        }
      }
    }
    placeholder_egg_summary <- bind_rows(placeholder_egg_summary_list)
  }
}

# ---- DATA COMBINATION AND FACTORING (ROBUST METHOD) ----

# 1. Safely combine original and placeholder data.
# Start with placeholders, which are known to be well-formed.
df_all <- placeholder_df_all
egg_summary <- placeholder_egg_summary

# If original data exists and has the 'alternative' column, add it.
if (exists("df_all_orig") && is.data.frame(df_all_orig) && "alternative" %in% names(df_all_orig)) {
  df_all <- bind_rows(
    ungroup(df_all_orig) %>% mutate(alternative = as.character(alternative)), 
    df_all
  )
}
if (exists("egg_summary_orig") && is.data.frame(egg_summary_orig) && "alternative" %in% names(egg_summary_orig)) {
  egg_summary <- bind_rows(
    ungroup(egg_summary_orig) %>% mutate(alternative = as.character(alternative)), 
    egg_summary
  )
}

# 2. Proceed with factoring only if we have data with an 'alternative' column.
if (nrow(df_all) > 0 && "alternative" %in% names(df_all)) {
  all_found_alts <- unique(as.character(df_all$alternative))
  
  # Separate numeric-like from text-like alternatives
  suppressWarnings({
    numeric_alts <- all_found_alts[!is.na(as.numeric(all_found_alts))]
    text_alts <- all_found_alts[is.na(as.numeric(all_found_alts))]
  })
  
  # Sort the numeric ones numerically, and the text ones alphabetically
  sorted_numeric_alts <- if(length(numeric_alts) > 0) as.character(sort(as.numeric(numeric_alts))) else character(0)
  sorted_text_alts <- sort(text_alts)
  
  # The final, robust order
  new_order <- c(sorted_numeric_alts, sorted_text_alts)
  
  # 3. Apply this final factor ordering to all relevant data frames.
  df_all$alternative <- factor(df_all$alternative, levels = new_order)
  
  if (nrow(egg_summary) > 0 && "alternative" %in% names(egg_summary)) {
    egg_summary$alternative <- factor(egg_summary$alternative, levels = new_order)
  }
  
  # Also update the temperature data list
  for(yr_name in names(env_ext_list)) {
    if (!is.null(env_ext_list[[yr_name]]) && nrow(env_ext_list[[yr_name]]) > 0 && "alternative" %in% names(env_ext_list[[yr_name]])) {
      env_ext_list[[yr_name]]$alternative <- factor(as.character(env_ext_list[[yr_name]]$alternative), levels = new_order)
    }
  }
}
