# ═══════════════════════════════════════════════════════════════════════════════
# SALMON POPULATION MODELING AND FORECASTING SYSTEM
# ═══════════════════════════════════════════════════════════════════════════════
# This script implements a comprehensive salmon population model that:
# 1. Analyzes historical spawning patterns based on water temperature
# 2. Calibrates survival models using Temperature-Dependent Mortality (TDM)
# 3. Forecasts future population dynamics under different environmental scenarios
# ═══════════════════════════════════════════════════════════════════════════════

# Clear the workspace to ensure a clean environment
rm(list=ls())

# ─── Libraries ────────────────────────────────────────────────────────────────
# Load required packages for data manipulation, parallel processing, and visualization
library(tidyverse)      # Core tidyverse packages for data manipulation
library(lubridate)      # Date/time manipulation
library(furrr)          # Parallel processing with future + purrr
library(data.table)     # High-performance data manipulation
library(compiler)       # Byte-compile R functions for speed
library(here)           # Relative file paths
library(ggrepel)        # Text/label positioning in ggplot2
library(MASS)           # Statistical functions
library(ordinal)        # Ordinal regression models (CLM)
library(ggridges)       # Ridge plots for distributions

# Load custom functions specific to this salmon counting project
source(here("SalmonCountR", "functions.R"))

# Set random seed for reproducibility across all stochastic processes
set.seed(123)

# ═══════════════════════════════════════════════════════════════════════════════

# ---- SECTION 1: INITIALIZE TIME PARAMETERS AND SIMULATION SETTINGS ----

# ═══════════════════════════════════════════════════════════════════════════════
# .----

#Alex, calibration for what, spawner to escapement model?

# Define the historical calibration period (14 years of observed data)
real_years <- 2011:2024 
n_calib <- length(real_years)  # Number of calibration years (14)

# Set up the forecast horizon
max_forecast <- 100  # Forecast 100 years into the future
n_sim <- n_calib + max_forecast  # Total simulation years (114)

# Create year sequences for different purposes
sim_years <- real_years[1] + seq(0, n_sim-1)  # All simulation years starting from 2011
forecast_years <- (max(real_years)+1):(max(sim_years))  # Future years only (2025-2124)

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 2: DATA INPUT DOCUMENTATION AND LOADING ----
# ═══════════════════════════════════════════════════════════════════════════════
# This section loads all required input data files with detailed documentation
# of expected formats and contents for each file type

# ── Data inputs overview ──────────────────────────────
# env_ext_list.rds:
#   • A master list of tibbles/data.frames, one per environment/alternative.
#   • Each element contains water temperature extracted by site and time.
#   • Expected columns include:
#       - site: Location identifier
#       - Date: Time index
#       - temp: Water temperature (numeric)
#   • Used to map environment → sites and to assemble df_all for analysis.
#

#Alex, dumb question but is this temperature data from the same source as env_ext_list.rds? In both cases is this empirical field temperature or modeled temperature based on various management alternatives

# df_all.rds:
#   • A consolidated table containing observed water temperature across
#     all alternatives, sites, and dates for analysis
#
# carcassdet_*.csv:
#   • Raw carcass detection records from field surveys
#   • Contains individual fish information including:
#     - Survey date, reach/section, sex, fork length, tag status
#   • Used to determine actual spawning timing and fish characteristics
#

#Alex, is the grandtab.csv filtered for particular run(s) of Chinook salmon somewhere in this code?
# grandtab_*.csv:
#   • California Department of Fish and Wildlife (CDFW) GrandTab escapement estimates
#   • Key columns: 
#     - End Year of Monitoring Period
#     - Population Estimate (renamed to 'spawners')
#   • Filtered to years 2011-2024 for this project

# ── Read data inputs ──────────────────────────────────
# Load the master list of environmental data (water temps by site/time for each alternative)
env_ext_list <- readRDS(
  here("SalmonCountR", "app_data", "env_ext_list.rds")
)

#Alex, minor point but the shiny app refers to alternative as "management alternatives and not "environmental alternatives"
# NOTE: 'env' throughout this code refers to environmental alternatives/scenarios
# Create a mapping between environmental scenarios and their associated sites
env_sites <- purrr::imap_dfr(env_ext_list, ~ {
  tibble(env = as.character(.y), site = unique(.x$site))
}) %>% distinct(env, site)                                              

# Load consolidated temperature data across all scenarios
df_all <- readRDS(
  here("SalmonCountR", "app_data", "df_all.rds")
)

# ═══════════════════════════════════════════════════════════════════════════════
# UTILITY FUNCTION: Clean up CSV files with trailing empty rows
# ═══════════════════════════════════════════════════════════════════════════════

#Alex, minor point but if you relocated the trim_trailing_text and other functions to the functions.R script it wouldn't muddy up your document outline window in R
# Many CSV exports contain trailing empty rows that need to be removed
trim_trailing_text <- function(df) {
  # Count non-NA numeric values in each row
  nnum <- df %>% mutate(.nn = rowSums(across(where(is.numeric), ~ !is.na(.)))) %>% pull(.nn)
  # Find the last row with actual data
  last <- suppressWarnings(max(which(nnum > 0), na.rm = TRUE))
  if (!is.finite(last)) return(df)
  # Return only rows up to the last row with data
  df[seq_len(last), , drop = FALSE]
}
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 3: LOAD AND PROCESS CARCASS SURVEY DATA ----
# ═══════════════════════════════════════════════════════════════════════════════
# Load raw carcass detection data from field surveys
# TN: Note - removing excess NA rows here to avoid downstream issues
carcass_raw <- readr::read_csv(
  here::here("SalmonCountR", "app_data", "carcassdet_1752789274_15.csv")
) %>%
  trim_trailing_text()  # Remove trailing empty rows
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 4: LOAD AND PROCESS ESCAPEMENT (SPAWNER) DATA ----
# ═══════════════════════════════════════════════════════════════════════════════
# Load historical escapement estimates from CDFW GrandTab
# TN: There are parsing issues (6 total, 1 did parse but incorrectly)
#Alex, same question as before. Do you need to omit some Chinook salmon runs from grandtab?
esc_obs <- read_csv(
  here("SalmonCountR", "app_data", "grandtab_1752793045_337.csv"),
  col_types = cols(
    `End Year of Monitoring Period` = col_character(),
    `Population Estimate`           = col_double()
  ),
  show_col_types = FALSE
) %>%
  trim_trailing_text() %>%  # Remove footer notes
  mutate(year = parse_number(`End Year of Monitoring Period`)) %>%
  filter(between(year, 2011, 2024)) %>%  # Keep only calibration period
  rename(spawners = `Population Estimate`) %>%
  arrange(year)
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 5: CONFIGURE PARALLEL PROCESSING ----
# ═══════════════════════════════════════════════════════════════════════════════
# Set up parallel processing to speed up computations

# Detect available CPU cores (on Windows this reports logical threads)
# Subtract 1 to leave one core for the operating system
n_threads <- parallel::detectCores(logical = TRUE) - 1

# Cap at 6 physical cores to avoid oversubscription
n_workers <- min(6, n_threads)

# Initialize parallel session for furrr package
plan(multisession, workers = n_workers)

# Configure data.table's internal C-level threading
setDTthreads(n_workers)
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 6: DEFINE TEMPERATURE-DEPENDENT MORTALITY (TDM) MODEL VARIANTS ----
# ═══════════════════════════════════════════════════════════════════════════════
# Define three TDM variants with different calibrations:
# 1. Exponential model with Water Forum 2020 calibration
# 2. Exponential model with SALMOD 2006 calibration  
# 3. Linear Martin model (no calibration needed)
tdm_defs <- tribble(
  ~model,       ~calib,           ~variant,
  "exp",        "WaterForum2020", "exp_WF",
  "exp",        "SALMOD2006",     "exp_SM",
  "lin_martin", NA,               "lin_Martin"
)
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 7: PREPARE TEMPERATURE AND DATE LOOKUPS FOR FAST ACCESS ----
# ═══════════════════════════════════════════════════════════════════════════════
# Create optimized lookup structures for temperature and date data

#Alex, Global comment: I'm not quite sure how to read the arrow within the context of the code description.
# Create site → temperature vector mappings for each environment
site_temps_list <- map(env_ext_list, ~ split(.x$temp, .x$site))

# Create site → Date vector mappings for each environment
site_dates_list <- map(env_ext_list, ~ split(.x$Date, .x$site))

# Create site → named-index lookup for fast date-to-index conversion
# This allows O(1) lookup of array position for any date
date_idx_list <- map(site_dates_list, function(dlist) {
  map(dlist, ~ setNames(seq_along(.x), as.character(.x)))
})
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 8: PROCESS CARCASS DATA AND ASSIGN TO SITES ----
# ═══════════════════════════════════════════════════════════════════════════════
# Transform raw carcass data into analysis-ready format
#Alex, it seems like "carcass_df" would be a more intuitive object name than "obs_df"
obs_df <- carcass_raw %>%
  mutate(
    # Convert survey date to Date object
    Date = as.Date(surveydate),
    # Estimate spawn date as 7 days before survey (decomposition time)
    spawn_dt = Date - days(7),
    # Determine brood year (Sept-Aug cycle: Sept+ = current year, <Sept = previous year)
    brood_year = if_else(month(spawn_dt) >= 9, year(spawn_dt), year(spawn_dt) - 1), 
#Alex, "survey section" is fairly intuitive to me and likely means sections of the American River.
#Alex, "analysis sites" is less intuitive to me. What does this mean?
    # Map survey sections to analysis sites
    site = case_when(
      section %in% c("NB","W","1a","1b","1a/1b","2") ~ "AveHazel",
      section %in% c("3")                            ~ "AveWatt",
      TRUE                                           ~ NA_character_
    )
  ) %>%
  # Select and rename columns for analysis (transmute avoids dplyr::select mask)
  dplyr::transmute(
    brood_year, 
    spawn_dt, 
    section,
    fork_length = .data[["forklength"]],  # Rename forklength to fork_length
    site
  )

# Calculate the proportion of observations at each site
# Used later for distributing simulated fish across sites
site_props <- obs_df %>%
  count(site) %>%
  mutate(prop = n / sum(n))
site_props
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 9: CREATE SPAWNING PERIOD BINS ----
# ═══════════════════════════════════════════════════════════════════════════════
# Divide the spawning season into temporal bins for analysis

# Ensure spawn_dt is properly formatted as Date
if (!inherits(obs_df$spawn_dt, "Date")) {
  obs_df <- obs_df %>% mutate(spawn_dt = as.Date(spawn_dt))
}

# Step 1: Find the anchor date (earliest spawning date in Aug-Dec)
# Prefer Aug-Dec dates; fallback to any date if none found
anchor_mmdd <- obs_df %>%
  filter(!is.na(spawn_dt)) %>%
  mutate(mmdd = format(spawn_dt, "%m-%d"), m = month(spawn_dt)) %>%
  filter(m >= 8) %>%  # August or later
  summarise(anchor = min(mmdd, na.rm = TRUE)) %>%
  pull(anchor)

# Fallback if no Aug-Dec dates found
if (length(anchor_mmdd) == 0 || is.na(anchor_mmdd)) {
  anchor_mmdd <- obs_df %>%
    filter(!is.na(spawn_dt)) %>%
    mutate(mmdd = format(spawn_dt, "%m-%d")) %>%
    summarise(anchor = min(mmdd, na.rm = TRUE)) %>%
    pull(anchor)
}

# Step 2: Define bin parameters
bin_width <- 10L  # Each bin is 10 days wide

# Calculate season day relative to anchor and determine number of bins needed
md <- format(obs_df$spawn_dt, "%m-%d")
season_year  <- year(obs_df$spawn_dt) - (md < anchor_mmdd)  # Adjust year if before anchor
season_start <- as.Date(paste0(season_year, "-", anchor_mmdd))
season_day   <- as.integer(obs_df$spawn_dt - season_start)
needed_bins  <- ceiling((max(season_day, na.rm = TRUE) + 1) / bin_width)
n_bins <- max(12L, needed_bins)  # Ensure at least 12 bins

# Step 3: Create function to assign dates to period bins
assign_period <- function(dates, anchor_mmdd, bin_width, n_bins) {
  md <- format(dates, "%m-%d")
  season_year  <- year(dates) - (md < anchor_mmdd)
  season_start <- as.Date(paste0(season_year, "-", anchor_mmdd))
  season_day   <- as.integer(dates - season_start)
  # Calculate bin index (1-based, capped at n_bins)
  bin_idx <- pmax(1L, pmin(floor(season_day / bin_width) + 1L, n_bins))
  # Return as ordered factor
  factor(paste0("p", bin_idx), levels = paste0("p", seq_len(n_bins)), ordered = TRUE)
}

# Assign period to each observation
obs_df <- obs_df %>%
  mutate(period = assign_period(spawn_dt, anchor_mmdd, bin_width, n_bins))

# Step 4: Create bin definitions for documentation and plotting
bin_defs <- tibble(period = paste0("p", seq_len(n_bins))) %>%
  mutate(
    # Use year 2000 as reference for bin boundaries
    start = as.Date(paste0("2000-", anchor_mmdd)) + (row_number() - 1) * bin_width,
    end   = start + (bin_width - 1)
  )

# Assign observations to spawn bins with proper year handling
obs_df$spawn_bin <- sapply(obs_df$spawn_dt, function(d) {
  if (is.na(d)) return(NA_character_)  # Skip NAs safely
  
  # TN: Account for January being in the next year
  dummy_date <- if_else(month(d) >= 9, 
                        as.Date(paste0("2000-", format(d, "%m-%d"))),
                        as.Date(paste0("2001-", format(d, "%m-%d"))))
  
  # Find matching bin
  match_idx <- which(dummy_date >= bin_defs$start & dummy_date <= bin_defs$end)
  if (length(match_idx) == 1) bin_defs$period[match_idx] else NA
})

# Convert to ordered factor
obs_df$spawn_bin <- factor(obs_df$spawn_bin, levels = bin_defs$period, ordered = TRUE)

# Extract unique years from observed data
yrs <- sort(unique(as.integer(obs_df$brood_year)))
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 10: VERIFY TEMPERATURE CONSISTENCY ACROSS ALTERNATIVES ----
# ═══════════════════════════════════════════════════════════════════════════════
# Check that observed temperatures are identical across alternatives
# (They should be the same since alternatives only differ in future projections)
#Alex, the previous comment line is unclear to me. 
env_ext_list %>%
  imap_dfr(~ mutate(.x, alt = .y)) %>%
  filter(year(Date) %in% yrs) %>%
  group_by(site, Date) %>%
  summarise(n_temp = n_distinct(temp), .groups = "drop") %>%
  summarise(all_equal = all(n_temp == 1))
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 11: CALCULATE OCTOBER/NOVEMBER MEAN TEMPERATURES ----
# ═══════════════════════════════════════════════════════════════════════════════
# These months are critical for spawning timing

# Collapse temperature data across alternatives (they're identical for historical period)
# Calculate mean temperatures for October and November by site and year
oct_nov_temps <- env_ext_list %>%
  imap_dfr(~ mutate(.x, alt = as.character(.y))) %>%  # Tag with alternative name
  filter(year(Date) %in% yrs, month(Date) %in% c(10, 11)) %>%  # Oct & Nov only
  group_by(site, year = year(Date), month = month(Date)) %>%  # No 'alt' - collapses across alternatives
  summarise(mean_temp = mean(temp, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = month, values_from = mean_temp, names_prefix = "m_") %>%
  rename(Oct = m_10, Nov = m_11)

# Sanity check: ensure unique site-year combinations (prevents many-to-many joins)
stopifnot(!any(duplicated(oct_nov_temps[c("site", "year")])))

# Join historical temperatures to observation data
#Alex, this is an example where it would be easier to undersand what data is being left joined if "obs_df" was called "carcass_df"
obs_df2 <- obs_df %>%
  mutate(brood_year = as.integer(brood_year)) %>%
  left_join(oct_nov_temps, by = c("site", "brood_year" = "year"))
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 12: STANDARDIZE TEMPERATURE VARIABLES ----
# ═══════════════════════════════════════════════════════════════════════════════
# Standardize temperatures using training data statistics (z-scores)

# Calculate mean and standard deviation for standardization
oct_mean <- mean(obs_df2$Oct, na.rm = TRUE)
oct_sd <- sd(obs_df2$Oct, na.rm = TRUE)
nov_mean <- mean(obs_df2$Nov, na.rm = TRUE)
nov_sd <- sd(obs_df2$Nov, na.rm = TRUE)

# Create analysis dataset with standardized temperatures
obs_fit <- obs_df2 %>%
  mutate(
    # Standardize temperatures (z-scores)
    Oct_std = (Oct - oct_mean) / oct_sd,
    Nov_std = (Nov - nov_mean) / nov_sd,
    # Ensure factors are properly set
    spawn_bin  = factor(spawn_bin, levels = bin_defs$period, ordered = TRUE),
    brood_year = factor(brood_year)
  ) %>%
  # Remove incomplete cases
  filter(!is.na(spawn_bin), !is.na(Oct_std), !is.na(Nov_std), !is.na(brood_year))
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 13: FIT CUMULATIVE LINK MODEL (CLM) FOR SPAWN TIMING ----
# ═══════════════════════════════════════════════════════════════════════════════
# Model spawning period as function of October and November temperatures

# Fit CLM without random effects
# This models the probability of spawning in each time bin based on temperature
spawn_clm <- clm(spawn_bin ~ Oct_std + Nov_std,
                 data = obs_fit, link = "logit")
summary(spawn_clm)

# Extract model coefficients
cf   <- coef(spawn_clm)
beta <- cf[c("Oct_std","Nov_std")]  # Temperature effects
zeta <- unname(cf[!(names(cf) %in% names(beta))])  # Threshold parameters
K    <- length(zeta) + 1  # Number of categories

# Get the bin levels from the model
#Alex, what do bin levels refer to again? 10 day spawning windows?
bins_model <- levels(spawn_clm$model$spawn_bin)
if (length(bins_model) != K) bins_model <- bin_defs$period[seq_len(K)]
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 14: BUILD DETERMINISTIC FORECAST ----
# ═══════════════════════════════════════════════════════════════════════════════
# Generate future spawning patterns based on projected temperatures

# Create/verify scaling parameters if not already made
#Alex, can you explain the month specific scaling a bit more
if (!exists("sc")) {
  sc <- obs_df2 %>%
    filter(brood_year %in% yrs) %>%
    summarise(
      o_m = mean(Oct, na.rm = TRUE), o_s = sd(Oct, na.rm = TRUE),
      n_m = mean(Nov, na.rm = TRUE), n_s = sd(Nov, na.rm = TRUE)
    )
}

# Build forecast temperatures for future years
yrs_forecast   <- forecast_years
forecast_temps <- build_forecast_temps(env_ext_list, yrs_forecast, sc)
stopifnot(all(c("env","sim_year","Oct","Nov","Oct_std","Nov_std") %in% names(forecast_temps)))

# Calculate CLM probabilities for forecast period (deterministic, no noise)
probs_all <- predict_clm_probs(beta, zeta, forecast_temps[, c("Oct_std","Nov_std")], offset = 0)
colnames(probs_all) <- bins_model

# Normalize probabilities to sum to 1
rs <- rowSums(probs_all)
rs[!is.finite(rs) | rs <= 0] <- 1
probs_all <- probs_all / rs

# Ensure columns match the bins actually present in observed data
present_bins <- levels(droplevels(obs_df$spawn_bin))
if (!identical(colnames(probs_all), present_bins)) {
  keep <- match(present_bins, colnames(probs_all))
  probs_all <- probs_all[, keep, drop = FALSE]
  colnames(probs_all) <- present_bins
}
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 15: JOIN ENVIRONMENTAL DATA AND CREATE SAMPLING POOLS ----
# ═══════════════════════════════════════════════════════════════════════════════
# Prepare data structures for simulation

# Join observations with environment mapping
#Alex can you elaborate on what observations and environment mean here?
#Alex, I think environment means management alternative but maybe not
#Alex, Nitpicky comment but "env" is pretty vague compared another description like "mgt_alt
obs_df_env <- obs_df %>% left_join(
  purrr::imap_dfr(env_ext_list, ~ tibble(env = as.character(.y), site = unique(.x$site))) %>%
    distinct(site, env),
  by = "site"
)

# Calculate average number of redds per environment
env_nredd <- obs_df_env %>%
  filter(brood_year %in% as.integer(2011:2024), !is.na(env)) %>%
  count(env, brood_year, name = "n") %>%
  group_by(env) %>%
  summarise(N_redd = max(1L, round(mean(n, na.rm = TRUE))), .groups = "drop")

# Create sampling pools for each environment

# These pools are used to randomly sample characteristics for simulated fish
env_pools <- obs_df_env %>%
  filter(brood_year %in% as.integer(2011:2024), !is.na(env)) %>%
  group_by(env) %>%
  summarise(
#Alex, I'd suggest saying something like "carcass survey sections" instead of just "sections" which is a fairly generic descriptor.
    section_pool     = list(section),      # Pool of observed sections
    fork_length_pool = list(fork_length),  # Pool of observed fork lengths
    site_pool        = list(unique(site)), # Pool of unique sites
    .groups = "drop"
  )

# Augment forecast temperatures with redd counts and sampling pools
#Alex, "augment" seems like perhaps a ambiguous or misleading word here.
ft_aug <- forecast_temps %>%
  left_join(env_nredd, by = "env") %>%
  left_join(env_pools, by = "env")
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 16: CREATE FAST DATE SAMPLING FUNCTION ----
# ═══════════════════════════════════════════════════════════════════════════════
# Efficiently sample spawn dates within bins

if (!exists("sample_dates_fast")) {
  # Pre-compute bin date components for fast sampling
#Alex, can you elaborate on what you mean by bin here. Is this the 10 day spawning period bin?
  bin_tbl <- bin_defs %>%
    transmute(
      period,
      start_m = lubridate::month(start), 
      start_d = lubridate::mday(start), 
      start_yoff = ifelse(start_m >= 10, 0L, 1L),  # Year offset for Oct+ vs Jan
      end_m   = lubridate::month(end),   
      end_d   = lubridate::mday(end),   
      end_yoff   = ifelse(end_m >= 10, 0L, 1L)
    )
  
  # Function to sample random dates within specified bins
  sample_dates_fast <- function(bins_chr, year_int) {
    # Match bins to pre-computed table
    idx <- match(bins_chr, bin_tbl$period)
    
    # Extract date components
    sm  <- bin_tbl$start_m[idx]; sd <- bin_tbl$start_d[idx]; sy <- year_int + bin_tbl$start_yoff[idx]
    em  <- bin_tbl$end_m[idx];   ed <- bin_tbl$end_d[idx];   ey <- year_int + bin_tbl$end_yoff[idx]
    
    # Construct start and end dates
    start <- as.Date(sprintf("%04d-%02d-%02d", sy, sm, sd))
    end   <- as.Date(sprintf("%04d-%02d-%02d", ey, em, ed))
    
    # Calculate bin lengths and sample random offsets
#Alex, aren't bin lengths always 10 days? Or is this referring to some other type of bin?
    len <- as.integer(end - start) + 1L
    len[len <= 0 | is.na(len)] <- 1L
    off <- if (length(len)) floor(runif(length(len), min = 0, max = pmax(1L, len))) else integer()
    
    # Return sampled dates
    start + off
  }
}
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 17: BUILD SIMULATED FUTURE SPAWNING DATA ----
# ═══════════════════════════════════════════════════════════════════════════════
# Generate synthetic spawning data for forecast period
out_list <- vector("list", nrow(ft_aug))

for (i in seq_len(nrow(ft_aug))) {
  # Extract parameters for this environment-year
  env_i <- ft_aug$env[i]
  yr_i  <- ft_aug$sim_year[i]
  n_i   <- ft_aug$N_redd[i]
  if (is.na(n_i) || n_i <= 0) n_i <- 1L  # Ensure at least 1 redd
  
  # Sample bins based on predicted probabilities
  #Alex, can you elaborate on this?
  p_i <- probs_all[i, ]  # Already normalized to sum to 1
  bins_i  <- sample(present_bins, n_i, replace = TRUE, prob = p_i)
  
  # Sample dates within selected bins
  dates_i <- sample_dates_fast(bins_i, yr_i)
  
  # Get sampling pools for this environment
  #Alex, is this generating sampling pools by carcass survey section?
  sec_pool <- ft_aug$section_pool[[i]]
  if (is.null(sec_pool)) sec_pool <- NA
  
  #Alex, is this generating sampling pools by fish fork length?
  fl_pool  <- ft_aug$fork_length_pool[[i]]
  if (is.null(fl_pool)) fl_pool <- NA_real_
  
  #Alex, is this generating sampling pools by site within the american river?
  st_pool  <- ft_aug$site_pool[[i]]
  if (is.null(st_pool)) {
    st_pool <- obs_df_env %>% filter(env == env_i) %>% pull(site) %>% unique()
  }
  
  # Create simulated observations
  #Alex, can you provide a little more context on what is meant by simulated observations and what they're used for in the big picture?
  out_list[[i]] <- tibble(
    env         = env_i,
    sim_year    = yr_i,
    spawn_dt    = dates_i,
    section     = sample(sec_pool, n_i, replace = TRUE),
    fork_length = sample(fl_pool,  n_i, replace = TRUE),
    site        = sample(st_pool,  n_i, replace = TRUE)
  )
}

# Combine all simulated future data
#Alex, how about "Combine all simulated data for future water years" Is that correct?
sim_future <- dplyr::bind_rows(out_list)
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 18: VISUALIZE CLM PREDICTIONS ----
# ═══════════════════════════════════════════════════════════════════════════════
# Create plots showing how temperature affects spawning timing probabilities

# Create pretty labels for bins (e.g., "p1: Oct 05–Oct 14")
#Alex, global comment. I'd replace all instances of "bins" with something like "10-day spawning period bins"
idx <- match(bins_model, bin_defs$period)
pretty_labels <- paste0(
  bins_model, ": ",
  format(bin_defs$start[idx], "%b %d"), "–", format(bin_defs$end[idx], "%b %d")
)

# Define safe ranges from training data
#Alex, this description is pretty vague. Can you provide more detail for what's happening here?
safe_range <- function(x, fallback = c(-2, 2)) {
  x <- x[is.finite(x)]
  if (length(x)) range(x) else fallback
}
rng_oct <- safe_range(obs_fit$Oct_std)
rng_nov <- safe_range(obs_fit$Nov_std)

# Create sequences for prediction
#Alex, do you mean predictions of spawning time? If so, I'd indicate that.
oct_range <- seq(rng_oct[1], rng_oct[2], length.out = 100)
nov_range <- seq(rng_nov[1], rng_nov[2], length.out = 200)

# Predict probabilities across October temperatures (November at mean)
#Alex, do you mean spawn timing probabilities?
newdata_oct <- data.frame(
  Oct_std = oct_range,
  Nov_std = mean(obs_fit$Nov_std, na.rm = TRUE)
)
pred_oct <- predict_clm_probs(beta, zeta, newdata_oct, offset = 0)
colnames(pred_oct) <- bins_model

# Convert to long format for plotting
pred_oct_long <- as.data.frame(pred_oct) %>%
  mutate(Oct_std = oct_range) %>%
  pivot_longer(all_of(bins_model), names_to = "spawn_bin", values_to = "probability") %>%
  mutate(spawn_bin = factor(spawn_bin, levels = bins_model, labels = pretty_labels))

# Predict probabilities across November temperatures (October at mean)
#Alex, do you mean spawn timing probabilities?
newdata_nov <- data.frame(
  Oct_std = mean(obs_fit$Oct_std, na.rm = TRUE),
  Nov_std = nov_range
)
pred_nov <- predict_clm_probs(beta, zeta, newdata_nov, offset = 0)
colnames(pred_nov) <- bins_model

# Convert to long format for plotting
pred_nov_long <- as.data.frame(pred_nov) %>%
  mutate(Nov_std = nov_range) %>%
  pivot_longer(all_of(bins_model), names_to = "spawn_bin", values_to = "probability") %>%
  mutate(spawn_bin = factor(spawn_bin, levels = bins_model, labels = pretty_labels))
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 19: CREATE VISUALIZATION PLOTS ----
# ═══════════════════════════════════════════════════════════════════════════════

# Plot 1: October temperature effects with facets (one panel per spawn bin)
p_oct_facet <- ggplot(pred_oct_long, aes(x = Oct_std, y = probability, color = spawn_bin)) +
  geom_line(size = 1) +
  facet_wrap(~ spawn_bin, ncol = 3) +  # Arrange in 3 columns
  scale_color_viridis_d(guide = "none") +  # Hide legend (redundant with facets)
  labs(
    title = "Predicted spawn-bin probability vs standardized October temperature",
    x = "Standardized October Temperature (Oct_std)",
    y = "Probability"
  ) +
  theme_minimal(base_size = 14)

# Plot 2: November temperature effects with facets
p_nov_facet <- ggplot(pred_nov_long, aes(x = Nov_std, y = probability, color = spawn_bin)) +
  geom_line(size = 1) +
  facet_wrap(~ spawn_bin, ncol = 3) +
  scale_color_viridis_d(guide = "none") +
  labs(
    title = "Predicted spawn-bin probability vs standardized November temperature",
    x = "Standardized November Temperature (Nov_std)",
    y = "Probability"
  ) +
  theme_minimal(base_size = 14)

# Plot 3: October temperature effects with legend (all bins on one plot)
p_oct_legend <- ggplot(pred_oct_long, aes(x = Oct_std, y = probability, color = spawn_bin)) +
  geom_line(size = 1) +
  scale_color_viridis_d(name = "Spawn Bin") +
  labs(
    title = "Predicted spawn-bin probability vs standardized October temperature",
    x = "Standardized October Temperature (Oct_std)",
    y = "Probability"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

# Plot 4: November temperature effects with legend
p_nov_legend <- ggplot(pred_nov_long, aes(x = Nov_std, y = probability, color = spawn_bin)) +
  geom_line(size = 1) +
  scale_color_viridis_d(name = "Spawn Bin") +
  labs(
    title = "Predicted spawn-bin probability vs standardized November temperature",
    x = "Standardized November Temperature (Nov_std)",
    y = "Probability"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

# Display the plots (choose which ones to show)
p_oct_facet
p_oct_legend
p_nov_facet
p_nov_legend
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 20: PREPARE FORECAST DATA FOR COMPARISON ----
# ═══════════════════════════════════════════════════════════════════════════════

# Add season-aligned dates and source label to forecast data
#Alex, I believe this is the first time you mention forecast data. Can you explain in a little more detail what you mean by that?
fc_df <- sim_future %>%
  mutate(
    season_date = season_posix(spawn_dt, anchor_mmdd),
    source = "Forecast"
  )
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 21: COMPARE OBSERVED VS MODELED SPAWNING PATTERNS ----
# ═══════════════════════════════════════════════════════════════════════════════
# This section validates the model by comparing actual vs predicted spawn timing

# Step 1: Standardize temperatures using training data statistics
sc <- obs_df2 %>%
  filter(brood_year %in% yrs) %>%
  summarise(
    o_m = mean(Oct, na.rm=TRUE), o_s = sd(Oct, na.rm=TRUE),
    n_m = mean(Nov, na.rm=TRUE), n_s = sd(Nov, na.rm=TRUE)
  )

# Calculate yearly mean temperatures and standardize
hist_temps <- obs_df2 %>%
  filter(brood_year %in% yrs) %>%
  group_by(brood_year) %>%
  summarise(Oct = mean(Oct, na.rm=TRUE), Nov = mean(Nov, na.rm=TRUE), .groups="drop") %>%
  mutate(
    Oct_std = (Oct - sc$o_m)/sc$o_s,
    Nov_std = (Nov - sc$n_m)/sc$n_s
  ) %>%
  rename(sim_year = brood_year)

# Step 2: Generate CLM probabilities for each historical year
#Alex, would you want to say the following instead? "Generate CLM probabilities for each year of empirical data"
cf   <- coef(spawn_clm)
beta <- cf[c("Oct_std","Nov_std")]
zeta <- unname(cf[!(names(cf) %in% names(beta))])
bins <- levels(spawn_clm$model$spawn_bin)

# Define prediction function for CLM probabilities
predict_clm_probs <- function(beta, zeta, newdata) {
  xb <- as.matrix(newdata[, names(beta), drop=FALSE]) %*% beta
  t(vapply(as.vector(xb), function(xi) {
    cp <- plogis(zeta - xi)
    p <- c(cp,1) - c(0,cp)
    p/sum(p)
  }, numeric(length(zeta)+1)))
}

# Calculate probabilities for historical years
#Alex, would you want to say the following instead? "Calculate CLM probabilities for each year of empirical data"
probs_hist <- predict_clm_probs(beta, zeta, hist_temps[, c("Oct_std","Nov_std")])
colnames(probs_hist) <- bins

# Step 3: Sample modeled dates matching observed counts per year
#Alex, can you provide broader context to all the coding comments for this step?
n_by_year <- obs_df %>% 
  filter(brood_year %in% yrs) %>% 
  count(brood_year, name="N")

set.seed(42)  # For reproducibility of sampling
sim_pred <- lapply(seq_len(nrow(hist_temps)), function(i){
  yr <- hist_temps$sim_year[i]
  # Get actual count for this year
  n  <- n_by_year$N[n_by_year$brood_year==yr]
  if (length(n)==0 || is.na(n) || n<=0) n <- 1L
  
  # Sample bins based on model probabilities
  p  <- probs_hist[i,]
  p <- p/sum(p)  # Ensure normalized
  b  <- sample(bins, n, replace=TRUE, prob=p)
  
  # Sample dates within selected bins
  tibble(
    sim_year = yr, 
    spawn_dt = sample_dates_fast(b, yr), 
    source = "Modeled"
  )
}) %>% bind_rows()

# Step 4: Extract observed carcass dates
sim_actual <- obs_df %>% 
  filter(brood_year %in% yrs) %>%
  transmute(
    sim_year = brood_year, 
    spawn_dt, 
    source = "Observed"
  )
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 22: CREATE COMPARISON VISUALIZATIONS ----
# ═══════════════════════════════════════════════════════════════════════════════

# (A) Season-aligned ridge plot comparing observed vs modeled distributions

# Find anchor date (earliest Aug-Dec spawning date)
anchor_mmdd <- sim_actual %>%
  mutate(m=month(spawn_dt), mmdd=format(spawn_dt,"%m-%d")) %>%
  filter(m>=8) %>% 
  summarise(a=min(mmdd, na.rm=TRUE)) %>% 
  pull(a)

# Fallback if no Aug-Dec dates
if (!length(anchor_mmdd) || is.na(anchor_mmdd)) anchor_mmdd <- "10-05"

anchor_base <- as.Date(paste0("2000-", anchor_mmdd))

# Function to convert dates to season-aligned dates
#Alex, it might be worth mentioning what is meant by "season-aligned dates" here. I'm assuming it refers to the water or spawning year convention but making that more explicit here and perhaps elsewhere would be helpful to the unacquainted.
season_posix <- function(d, mmdd){
  md <- format(d, "%m-%d")
  y0 <- year(d) - (md < mmdd)  # Adjust year if before anchor
  anchor <- as.Date(paste0(y0, "-", mmdd))
  anchor_base + as.integer(d - anchor)
}

# Combine observed and modeled data with season alignment
comp_df <- bind_rows(sim_actual, sim_pred) %>%
  mutate(season_date = season_posix(spawn_dt, anchor_mmdd))

# Set up x-axis breaks
rng <- range(comp_df$season_date, na.rm=TRUE)
x_breaks <- seq(rng[1], rng[2], by = "10 days")

# Create ridge plot showing density distributions by year
p_ridge <- ggplot(comp_df, aes(x=season_date, y=factor(sim_year))) +
  geom_density_ridges_gradient(
    aes(height=after_stat(density), fill=source),
    scale=2,                    # Height scaling
    rel_min_height=0.001,      # Minimum height to show
    linewidth=0.2,             # Line thickness
    bandwidth=10,              # Smoothing bandwidth
    trim=TRUE                  # Trim tails
  ) +
  scale_x_date(breaks=x_breaks, date_labels="%b %d", expand=c(0.01,0)) +
  scale_fill_manual(values=c(Observed="steelblue", Modeled="tomato")) +
  labs(
    title="Observed vs Modeled spawn timing by year (2011–2024)",
    x=paste0("Season day (anchor ", anchor_mmdd, ")"), 
    y="Year", 
    fill=NULL
  ) +
  theme_minimal(base_size=13)
p_ridge

# (B) Calendar-date boxplot comparing medians and distributions
p_box <- ggplot(bind_rows(sim_actual, sim_pred),
                aes(x=factor(sim_year), y=spawn_dt, fill=source)) +
  geom_boxplot(
    alpha=0.65,                              # Transparency
    position=position_dodge(width=0.7),      # Side-by-side boxes
    outlier.alpha=0.3                        # Outlier transparency
  ) +
  scale_fill_manual(values=c(Observed="steelblue", Modeled="tomato")) +
  labs(
    title="Observed vs Modeled spawning dates by year (2011–2024)",
    x="Brood year", 
    y="Spawn date", 
    fill=NULL
  ) +
  theme_minimal(base_size=13)
p_box
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 23: BUILD FINAL SIMULATION DATASETS ----
# ═══════════════════════════════════════════════════════════════════════════════
# Combine observed and simulated data for Temperature-Dependent Mortality modeling

# Combine observed (historical) and simulated (future) spawning data
sim_redds <- bind_rows(sim_actual, sim_future) %>%
  arrange(env, sim_year, spawn_dt) %>%
  dplyr::select(env, sim_year, spawn_dt, site)

# Convert to data.table for efficient processing
if (!data.table::is.data.table(sim_redds)) data.table::setDT(sim_redds)
data.table::setkey(sim_redds, env, sim_year)

# Split by environment-year for parallel processing
# drop=TRUE removes the key columns from each chunk
sim_redds_split <- split(sim_redds, by = c("env","sim_year"), drop = TRUE, keep.by = FALSE)

# Keep only essential columns in each chunk
sim_redds_split <- lapply(sim_redds_split, function(dt) dt[, .(spawn_dt, site)])

# Calculate median spawn dates per environment-year
#Alex, another question about "environment". Does this mean management alternative? I think I mentioned this earlier by environment is kind of a vague descriptor that prevents an intuitive understanding of what's going on.
spawn_medians_env_year <- sim_redds %>%
  as_tibble() %>%
  group_by(env, sim_year) %>%
  summarise(
    spawn_dt = as.Date(median(as.numeric(spawn_dt), na.rm = TRUE), origin = "1970-01-01"),
    .groups = "drop"
  )

# Get unique environments
envs <- sort(unique(spawn_medians_env_year$env))

# Build spawn date vectors for each environment
spawn_dates_by_env <- map(
  .x = set_names(envs),
  .f = ~ build_spawn_vec_for_env(
    df_env    = filter(spawn_medians_env_year, env == .x),
    sim_years = sim_years
  )
)
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 24: VISUALIZE FUTURE SPAWNING DISTRIBUTIONS ----
# ═══════════════════════════════════════════════════════════════════════════════
# Create multiple visualizations of projected spawning patterns

# Add dummy dates for seasonal alignment
sim_future_comp <- sim_future %>% 
  mutate(dummy_date = as.Date(paste0("2000-", format(spawn_dt, "%m-%d"))))

# Recode Oct-Dec as year 2000 and Jan as 2001 for visualization
sim_future_comp <- sim_future_comp %>%
  mutate(
    m = month(dummy_date),
    d = day(dummy_date),
    dummy_date_season = if_else(
      m >= 10,
      make_date(2000, m, d),  # Oct-Dec → 2000
      make_date(2001, m, d)   # Jan → 2001
    )
  )

# Plot 1: Density distributions for all years overlaid
ggplot(sim_future_comp, aes(x = dummy_date, color = factor(sim_year), group = sim_year)) +
  geom_density(alpha = 0.3) +
  labs(
    title = "Spawning Date Distributions (all years aligned to same season)",
    x = "Spawning Date (Oct–Jan, dummy year)",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)

# Plot 2: Ridge plot showing distributions by year
ggplot(sim_future_comp, aes(x = dummy_date_season, y = factor(sim_year), fill = factor(sim_year))) +
  geom_density_ridges(
    scale = 3,               # Height scaling
    rel_min_height = 0.01,   # Minimum height
    alpha = 0.6              # Transparency
  ) +
  scale_x_date(
    limits = c(as.Date("2000-10-01"), as.Date("2001-01-30")),
    date_labels = "%b %d"
  ) +
  labs(
    title = "Spawning Date Distributions by Year",
    x = "Spawning Date (Oct–Jan, dummy year)",
    y = "Simulation Year"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# Calculate summary statistics for spawn dates
summary_dates <- sim_future_comp %>%
  dplyr::group_by(sim_year) %>%
  dplyr::summarize(
    median_spawn = as.Date(median(as.numeric(spawn_dt), na.rm = TRUE), origin = "1970-01-01"),
    p10 = as.Date(quantile(as.numeric(spawn_dt), 0.10, na.rm = TRUE), origin = "1970-01-01"),
    p90 = as.Date(quantile(as.numeric(spawn_dt), 0.90, na.rm = TRUE), origin = "1970-01-01"),
    .groups = "drop"
  )

# Plot 3: Time series of median spawn date with uncertainty bands
ggplot(summary_dates, aes(x = sim_year)) +
  geom_line(aes(y = median_spawn), color = "blue") +
  geom_ribbon(aes(ymin = p10, ymax = p90), alpha = 0.2, fill = "blue") +
  labs(
    title = "Median and Spread of Spawn Dates by Year",
    x = "Simulation Year",
    y = "Spawn Date"
  ) +
  theme_minimal(base_size = 14)

# Plot 4: Heatmap of spawning intensity
sim_future_comp %>%
  count(sim_year, dummy_date) %>%
  ggplot(aes(x = dummy_date, y = sim_year, fill = n)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(
    title = "Heatmap of Spawning Intensity by Date and Year",
    x = "Spawning Date (Oct–Jan, dummy year)",
    y = "Simulation Year"
  ) +
  theme_minimal(base_size = 14)
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 25: CONFIGURE PARALLEL PROCESSING FOR TDM CALCULATIONS ----
# ═══════════════════════════════════════════════════════════════════════════════

# Increase memory limit for parallel workers (4 GiB)
options(future.globals.maxSize = 4 * 1024^3)

# Reset parallel plan with appropriate worker count
plan(multisession, workers = min(4, n_workers))

# Extract real median spawn dates from observed data
spawn_dates_real <- sim_redds %>%
  arrange(sim_year) %>%
  group_by(sim_year) %>%
  summarize(spawn_dt = median(spawn_dt), .groups="drop") %>%
  pull(spawn_dt)

# Repeat the month-day patterns to fill the entire simulation period
spawn_dm <- rep(spawn_dates_real, length.out = n_sim)

# Generate simulation year sequence
sim_years_vec <- (real_years[1] + seq_len(n_sim) - 1)
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 26: TEMPERATURE-DEPENDENT MORTALITY (TDM) CALCULATIONS ----
# ═══════════════════════════════════════════════════════════════════════════════
# This is the core TDM processing loop that calculates egg-to-fry survival

# Step 1: Compile performance-critical functions for speed
hatch_model_c     <- cmpfun(hatch_model)
emergence_model_c <- cmpfun(emergence_model)
tdm_exp_c        <- cmpfun(tdm_exp)
tdm_lin_c        <- cmpfun(tdm_lin_martin)

# Convert sim_redds to data.table if not already
setDT(sim_redds)

# Split data by simulation year for parallel processing
sim_redds_split <- split(sim_redds[, .(spawn_dt, site)], sim_redds$sim_year)

# Step 2: Cache environment assets for fast access
env_cache <- lapply(names(env_ext_list), function(env_nm) {
  list(
    env_nm       = env_nm,
    date_idx_env = date_idx_list[[env_nm]],
    temps_env    = site_temps_list[[env_nm]],
    date_min     = min(env_ext_list[[env_nm]]$Date),
    date_max     = max(env_ext_list[[env_nm]]$Date)
  )
})
names(env_cache) <- names(env_ext_list)

# Ensure TDM definitions use character types (not factors)
tdm_defs <- tdm_defs %>%
  dplyr::mutate(
    model   = as.character(model),
    calib   = as.character(calib),
    variant = as.character(variant)
  )
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 27: DATA QUALITY CONTROL - REMOVE INVALID SITES ----
# ═══════════════════════════════════════════════════════════════════════════════
# Remove sites that don't have temperature data across all environments

# Find sites that have temperature data in ALL environments
sites_with_temps <- Reduce(
  intersect,
  lapply(env_cache, function(ec) names(ec$date_idx_env))
)

# Filter sim_redds_split to only include valid sites
sim_redds_split <- lapply(
  sim_redds_split,
  function(dt) {
    if (is.null(dt) || !nrow(dt)) return(dt)
    # Ensure required columns exist
    stopifnot(all(c("site", "spawn_dt") %in% names(dt)))
    # Keep only sites with temperature data
    dt[dt$site %in% sites_with_temps, , drop = FALSE]
  }
)
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 28: BUILD LOOKUP TABLES FOR FAST DATE-POSITION MAPPING ----
# ═══════════════════════════════════════════════════════════════════════════════
# Create indexed lookup tables for efficient date-to-array-position conversion

env_lookup <- lapply(env_cache, function(ec) {
  # Build data.table with site, redd date (rdr), and position
  DT <- rbindlist(
    lapply(names(ec$date_idx_env), function(s) {
      v <- ec$date_idx_env[[s]]
      if (is.null(v) || !length(v)) {
        # Empty table if no data
        data.table(
          site = factor(character(0)), 
          rdr = as.IDate(integer(0)), 
          pos = integer(0)
        )
      } else {
        data.table(
          site = s,
          rdr  = as.IDate(names(v)),  # Convert date strings to IDate
          pos  = as.integer(v)         # Array position
        )
      }
    }),
    use.names = TRUE, 
    fill = TRUE
  )
  
  # Convert site to factor and set key for fast lookups
  DT[, site := as.factor(site)]
  setkey(DT, site, rdr)
  DT
})
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 29: IMPLEMENT MEMOIZATION FOR SURVIVAL CALCULATIONS ----
# ═══════════════════════════════════════════════════════════════════════════════
# Cache survival calculations to avoid redundant computations

# Create environment for memoization cache
.mem_env <- new.env(parent = emptyenv())

# Memoized survival calculation function
memo_surv <- function(env_nm, site, rdr, model, calib) {
  # Create unique key for this calculation
  k <- paste(env_nm, site, as.integer(rdr), model, calib, sep = "|")
  
  # Check cache
  hit <- .mem_env[[k]]
  if (!is.null(hit)) return(hit)
  
  # Calculate if not cached
  ec <- env_cache[[env_nm]]
  out <- compute_surv_by_atu(
    rdr          = rdr,
    site         = site,
    date_idx_env = ec$date_idx_env,
    temps_env    = ec$temps_env,
    model        = model,
    calib        = calib
  )
  
  # Store in cache
  .mem_env[[k]] <- out
  out
}

# Vectorized wrapper for survival calculations
compute_surv_vec <- function(env_nm, site, rdr, model, calib) {
  vapply(
    seq_along(site),
    function(i) memo_surv(env_nm, site[i], rdr[i], model, calib),
    numeric(1)
  )
}
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 30: PREPARE REDD DATA FOR ONE ENVIRONMENT-YEAR ----
# ═══════════════════════════════════════════════════════════════════════════════

pairs_for_env_year <- function(red_this, sim_yr, env_nm) {
  # Get lookup table for this environment
  lk <- env_lookup[[env_nm]]
  if (is.null(lk) || !nrow(lk)) return(NULL)
  
  # Calculate cohort redd date from spawn date
  # Account for brood year vs calendar year difference
  mm  <- data.table::month(red_this$spawn_dt)
  dd  <- lubridate::day(red_this$spawn_dt)
  yy  <- fifelse(mm >= 9L, sim_yr, sim_yr + 1L)  # Sept+ stays in brood year
  rdr <- as.IDate(lubridate::make_date(yy, mm, dd))
  
  # Create pairs and collapse duplicates
  pairs <- data.table(site = red_this$site, rdr = rdr)[
    , .N, by = .(site, rdr)  # Count occurrences
  ]
  
  if (!nrow(pairs)) return(NULL)
  
  # Join to ensure (site,rdr) exists in the environment index
  setkeyv(pairs, c("site", "rdr"))
  pairs <- lk[pairs, nomatch = 0L]  # Inner join - keep only valid pairs
  if (!nrow(pairs)) return(NULL)
  
  # Return columns: site, rdr, pos, N
  pairs[]
}
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 31: EVALUATE TDM FOR ONE YEAR ACROSS ALL VARIANTS ----
# ═══════════════════════════════════════════════════════════════════════════════

eval_year <- function(sim_yr, sim_redds_split, env_cache, tdm_defs) {
  # Get redd data for this year
  red_this <- sim_redds_split[[as.character(sim_yr)]]
  if (is.null(red_this) || !nrow(red_this)) return(data.table())
  
  env_names <- names(env_cache)
  
  # Process each environment
  env_tables <- lapply(env_names, function(env_nm) {
    # Get site-redd pairs for this environment-year
    pairs <- pairs_for_env_year(red_this, sim_yr, env_nm)
    if (is.null(pairs) || !nrow(pairs)) return(NULL)
    
    # Calculate survival for each TDM variant
    rbindlist(lapply(seq_len(nrow(tdm_defs)), function(i) {
      # Compute survivals for all redds
      survs <- compute_surv_vec(
        env_nm = env_nm,
        site   = pairs$site,
        rdr    = pairs$rdr,
        model  = tdm_defs$model[i],
        calib  = tdm_defs$calib[i]
      )
      
      # Calculate weighted mean survival
      data.table(
        sim_year      = sim_yr,
        env           = env_nm,
        variant       = tdm_defs$variant[i],
        method        = "observed",
        mean_cum_surv = {
          wsum <- sum(pairs$N)
          if (!is.finite(wsum) || wsum <= 0) {
            NA_real_
          } else {
            sum(survs * pairs$N, na.rm = TRUE) / wsum
          }
        }
      )
    }), use.names = TRUE, fill = TRUE)
  })
  
  # Filter out NULLs and combine
  env_tables <- Filter(Negate(is.null), env_tables)
  if (!length(env_tables)) return(data.table())
  rbindlist(env_tables, use.names = TRUE, fill = TRUE)
}
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 32: OPTIMIZE AND EXECUTE TDM CALCULATIONS IN PARALLEL ----
# ═══════════════════════════════════════════════════════════════════════════════

# Compile hot functions for additional speed
if (requireNamespace("compiler", quietly = TRUE)) {
  compute_surv_vec    <- compiler::cmpfun(compute_surv_vec)
  pairs_for_env_year  <- compiler::cmpfun(pairs_for_env_year)
  eval_year           <- compiler::cmpfun(eval_year)
}

# Set data.table threads to cooperate with parallel futures
data.table::setDTthreads(2L)

# Set up parallel execution plan
plan(multisession, workers = max(1L, parallel::detectCores() - 1L))

# Execute TDM calculations in parallel across all years
results_obs_fast <- furrr::future_map_dfr(
  sim_years,
  ~eval_year(.x, sim_redds_split, env_cache, tdm_defs),
  .options = furrr::furrr_options(seed = TRUE, scheduling = 20)
)

# results_obs_fast now contains: sim_year, env, variant, method="observed", mean_cum_surv
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 33: LIFE CYCLE MODEL CALIBRATION AND FORECASTING ----
# ═══════════════════════════════════════════════════════════════════════════════
# Calibrate life-cycle parameters and run population projections

# Select reference environment for calibration (using first available)
ref_env <- names(env_ext_list)[1]

# Step 1: Summarize TDM egg-to-fry survival results
egg_summary <- results_obs_fast %>%
  arrange(env, variant, sim_year) %>%
  group_by(env, variant, sim_year) %>%
  summarise(mean_cum_surv = mean(mean_cum_surv, na.rm = TRUE), .groups = "drop")

# Get unique variant names for processing
variant_names <- egg_summary %>% pull(variant) %>% unique() %>% sort()

# Step 2: Create survival lookups for calibration and forecasting

# Survival vectors for reference environment (calibration only)
surv_lookup_by_variant <- results_obs_fast %>%
  filter(env == ref_env, sim_year %in% real_years) %>%
  arrange(variant, sim_year) %>%
  group_by(variant) %>%
  summarise(surv_vec = list(mean_cum_surv), .groups = "drop") %>%
  tibble::deframe()  # Convert to named list: variant -> numeric vector

# Full survival lookup for all environment-variant combinations
surv_lookup_full <- egg_summary %>%
  arrange(env, variant, sim_year) %>%
  group_by(env, variant) %>%
  summarise(surv_vec = list(mean_cum_surv), .groups = "drop") %>%
  mutate(key = paste(env, variant, sep = "_")) %>%
  dplyr::select(key, surv_vec) %>%
  tibble::deframe()  # Named list: "env_variant" -> survival vector

# Step 3: Define base life-cycle parameters (without site-specific covariates)
base_P <- list(
  female_fraction = 0.5,    # Proportion of spawners that are female
  fec = 5522,               # Fecundity: eggs per female
  S0 = 0.347,               # Egg survival rate
  K_spawners = 12493,       # Carrying capacity for spawners
  SAR_mean = NA_real_,      # Smolt-to-Adult Return rate (to be calibrated)
  SAR_sd = 0.00237,         # Standard deviation of SAR
  # Age structure: probability of returning at age 3, 4, or 5
  lag_probs = c(`3` = 0.828982, `4` = 0.168885, `5` = 0.002105),  # From CWT data 2011-2024
  rear_surv = NA_real_      # Rearing survival (to be calibrated)
)

# Step 4: Set up calibration parameters
obs_spawners  <- esc_obs$spawners         # Observed spawners 2011-2024
S_seed_calib  <- obs_spawners[1:3]        # Seed population (first 3 years)
n_calib       <- length(obs_spawners)     # Number of calibration years (14)
fit_idx       <- (length(S_seed_calib) + 1):n_calib  # Years to fit (4-14)
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 34: CALIBRATE LIFE-CYCLE PARAMETERS ----
# ═══════════════════════════════════════════════════════════════════════════════
# Use optimization to find SAR and rearing survival that best match observations

# Select reference environment for calibration
ref_env <- names(env_ext_list)[1]

# Calculate degree-days for adult spawners in reference environment
deg_day_cal_ref <- deg_day_cal_for(ref_env)
stopifnot(
  length(deg_day_cal_ref) == length(real_years),
  all(is.finite(deg_day_cal_ref))
)

# Run calibration for each TDM variant in parallel
calib_results <- furrr::future_map_dfr(
  variant_names,
  function(v) {
    # Optimize SAR_mean and rear_surv to minimize sum of squared errors
    opt <- optim(
      par    = c(0.0025, 0.5419),  # Initial values from O'Farrell et al. (2018)
      fn     = modular_sse,         # Objective function (sum of squared errors)
      variant= v,                   # TDM variant
      method = "L-BFGS-B",          # Bounded optimization
      lower  = c(0, 0),             # Lower bounds
      upper  = c(1, 1)              # Upper bounds (both are probabilities)
    )
    
    # Return calibrated parameters
    tibble::tibble(
      variant   = v,
      SAR_mean  = opt$par[1],
      rear_surv = opt$par[2],
      sse       = opt$value
    )
  },
  .options = furrr::furrr_options(seed = TRUE)
)

# Create parameter lists for each variant-environment combination
base_P_list <- calib_results %>%
  split(.$variant) %>%
  purrr::map(function(df_v) {
    # Extract calibrated values for this variant
    SARv  <- df_v$SAR_mean[1]
    rearv <- df_v$rear_surv[1]
    
    # Create parameter set for each environment
    rlang::set_names(
      lapply(names(env_ext_list), function(env_nm) {
        P <- base_P
        P$SAR_mean  <- SARv
        P$rear_surv <- rearv
        P
      }),
      nm = names(env_ext_list)
    )
  })
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 35: GENERATE CALIBRATION PREDICTIONS FOR VALIDATION ----
# ═══════════════════════════════════════════════════════════════════════════════
# Compare model predictions with observed data during calibration period

calib_pred_by_variant <- rlang::set_names(
  lapply(variant_names, function(v) {
    # Get parameters for this variant
    P0 <- base_P_list[[v]][[ref_env]]
    
    # Run simulation for calibration period
    out <- simulate_variant(
      surv_vec       = surv_lookup_by_variant[[v]][1:n_calib],
      P              = P0,
      years          = n_calib,
      S_init         = S_seed_calib,
      SAR_vec        = rep(P0$SAR_mean, n_calib),
      K_spawners_vec = rep(P0$K_spawners, n_calib),
      deg_day_adult  = deg_day_cal_ref,  # Use consistent degree-days
      sim_years_vec  = real_years
    )
    
    # Compare predictions with observations
    tibble(
      year      = real_years,
      observed  = esc_obs$spawners,
      predicted = out$spawners,
      SAR_mean  = P0$SAR_mean,
      rear_surv = P0$rear_surv,
      sse       = sum((out$spawners[fit_idx] - esc_obs$spawners[fit_idx])^2)
    )
  }),
  variant_names
)
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 36: PREPARE SEED POPULATIONS FOR FORECASTING ----
# ═══════════════════════════════════════════════════════════════════════════════
# Use the last 3 years of calibrated run as initial conditions for forecasts

env_seed <- ref_env
k <- length(S_seed_calib)  # Number of seed years (3)

# Build forecast seed populations for each variant
S_seed_fore_list <- rlang::set_names(
  purrr::map(calib_results$variant, function(v) {
    years_cal <- length(real_years)
    
    # Get survival vector for this variant-environment
    surv_vec  <- surv_lookup_full[[paste(ref_env, v, sep = "_")]][1:years_cal]
    
    # Get parameters for this variant
    Ptmp <- base_P_list[[v]][[ref_env]]
    
    # Calculate degree-days for calibration period
    deg_day_cal <- compute_deg_day_adult(
      env_nm       = ref_env,
      sim_years    = real_years,
      spawn_dates  = spawn_dates_by_env[[ref_env]][match(real_years, sim_years)],
      env_ext_list = env_ext_list
    )
    
    # Run full calibration simulation
    out <- simulate_variant(
      surv_vec       = surv_vec,
      P              = Ptmp,
      years          = years_cal,
      S_init         = S_seed_calib,
      SAR_vec        = rep(Ptmp$SAR_mean, years_cal),
      K_spawners_vec = rep(Ptmp$K_spawners, years_cal),
      deg_day_adult  = deg_day_cal,
      sim_years_vec  = real_years
    )
    
    # Return last 3 years as seed for forecasts
    tail(out$spawners, length(S_seed_calib))
  }),
  calib_results$variant
)
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 37: RUN POPULATION FORECASTS ----
# ═══════════════════════════════════════════════════════════════════════════════
# Project populations forward under different scenarios

# Get all environment-variant combinations
keys <- names(surv_lookup_full)

# Configure stochasticity options (default: deterministic)
use_stochastic_SAR <- FALSE  # Can be changed in Shiny app

# Define stochastic SAR options for potential use in UI
stoch_SAR_opts <- list(
  model       = "normal",      # Distribution type: "normal", "lognormal", "beta", "gamma"
  mean        = base_P$SAR_mean,
  sd          = base_P$SAR_sd,
  shape1      = 2,             # Beta/gamma shape parameter 1
  shape2      = 5,             # Beta/gamma shape parameter 2
  timing      = "pulse",       # When to apply: "all", "block", "pulse"
  block_years = 20:30,         # Years for block application
  pulse_years = c(10, 15, 20, 25, 30, 35, 40),  # Years for pulse application
  pulse_sd    = 0.002          # Additional SD for pulse years
)

# Run forecasts for all environment-variant combinations
results_full <- purrr::map_dfr(keys, function(key) {
  # Parse key to get environment and variant
  parts  <- strsplit(key, "_")[[1]]
  env_nm <- parts[1]
  var_nm <- paste(parts[-1], collapse = "_")
  
  # Get seed population for this variant
  seed_vec <- S_seed_fore_list[[var_nm]]
  
  # Run forecast simulation
  sim_forecast_fn(
    var_nm,
    env_nm,
    flow_cfs = NULL,              # No flow covariates in this version
    S_seed   = seed_vec,
    spawn_dates_by_env = spawn_dates_by_env  # Pass spawn dates
  )()
})
# .----

# ═══════════════════════════════════════════════════════════════════════════════
# ---- SECTION 38: SAVE ALL OUTPUTS FOR SHINY APPLICATION ----
# ═══════════════════════════════════════════════════════════════════════════════
# Save all processed data and results for use in the Shiny dashboard

# Save calibration results (SAR and rearing survival parameters)
saveRDS(calib_results,         here("SalmonCountR","app_data","calib_results.rds"))

# Save calibration predictions for validation plots
saveRDS(calib_pred_by_variant, here("SalmonCountR","app_data","calib_pred_by_variant.rds"))

# Save full forecast results (main output)
saveRDS(results_full,       here("SalmonCountR","app_data","results_full.rds"))

# Save egg-to-fry survival summaries
saveRDS(egg_summary,        here("SalmonCountR","app_data","egg_summary.rds"))

# Save complete survival lookup table
saveRDS(surv_lookup_full,   here("SalmonCountR","app_data","surv_lookup_full.rds"))

# Save calibrated parameter sets
saveRDS(base_P_list,        here("SalmonCountR","app_data","base_P_list.rds"))

# Save initial seed population for calibration
saveRDS(S_seed_calib,       here("SalmonCountR","app_data","S_seed_calib.rds"))

# Save seed populations for forecasting
saveRDS(S_seed_fore_list,   here("SalmonCountR","app_data","S_seed_fore_list.rds"))

# Save stochastic SAR options
saveRDS(stoch_SAR_opts,     here("SalmonCountR","app_data","stoch_SAR_opts.rds"))

# Save simulation years vector
saveRDS(sim_years,          here("SalmonCountR","app_data","sim_years.rds"))

# Save spawn dates by environment
saveRDS(spawn_dates_by_env, here("SalmonCountR","app_data","spawn_dates_by_env.rds"))

# ═══════════════════════════════════════════════════════════════════════════════
# END OF SCRIPT
# ═══════════════════════════════════════════════════════════════════════════════
# This completes the salmon population modeling pipeline.
# All results are saved and ready for visualization in the Shiny application.