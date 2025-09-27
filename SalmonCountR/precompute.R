# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║         FALL-RUN CHINOOK SALMON POPULATION MODELING AND FORECASTING           ║
# ╔══════════════════════════════════════════════════════════════════════════════╗
# This script implements a comprehensive fall-run Chinook salmon population model that:
# 1. Analyzes historical spawning patterns based on water temperature
# 2. Calibrates survival models using Temperature-Dependent Mortality (TDM)
# 3. Forecasts future population dynamics under different management alternatives
# ╚══════════════════════════════════════════════════════════════════════════════╝

# Clear the workspace to ensure a clean environment
rm(list=ls())

# ---- 1. LIBRARIES ----
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

# ---- 2. TIME PARAMETERS ----
# Define the historical calibration period for the spawner-escapement model
# This period contains observed spawner abundance data from CDFW GrandTab
real_years <- 2011:2024 
n_calib <- length(real_years)  # Number of calibration years (14)

# Set up the forecast horizon
max_forecast <- 100  # Forecast 100 years into the future
n_sim <- n_calib + max_forecast  # Total simulation years (114)

# Create year sequences for different purposes
sim_years <- real_years[1] + seq(0, n_sim-1)  # All simulation years starting from 2011
forecast_years <- (max(real_years)+1):(max(sim_years))  # Future years only (2025-2124)

# ---- 3. DATA INPUT ----
# ---- 3.1 Environmental Temperature Data ----
# Load the master list of water temperature data organized by management alternative
# Each list element represents a different water management alternative/scenario
# Contains modeled daily water temperatures for different operational scenarios
env_ext_list <- readRDS(
  here("SalmonCountR", "app_data", "env_ext_list.rds")
)

# NOTE: Throughout this code 'env' refers to management alternatives (e.g., different 
# operational scenarios for water releases, not "environmental" in the general sense)
# Create a mapping between management alternatives and their associated river sites
mgt_alt_sites <- purrr::imap_dfr(env_ext_list, ~ {
  tibble(mgt_alt = as.character(.y), site = unique(.x$site))
}) %>% distinct(mgt_alt, site)                                              

# ---- 3.2 Consolidated Temperature Data ----
# Load temperature data consolidated across all management alternatives
# This contains the same modeled temperature data as env_ext_list but in a single dataframe
df_all <- readRDS(
  here("SalmonCountR", "app_data", "df_all.rds")
)



# ---- 4. CARCASS SURVEY DATA ----
# Load raw carcass detection data from American River field surveys
# These are empirical observations of dead salmon used to estimate spawn timing
carcass_raw <- readr::read_csv(
  here::here("SalmonCountR", "app_data", "carcassdet_1752789274_15.csv"),
  show_col_types = FALSE
) %>%
  trim_trailing_text()  # Remove trailing empty rows that cause parsing issues

# ---- 5. ESCAPEMENT DATA ----
# Load historical escapement estimates from CDFW GrandTab
# This file contains fall-run Chinook salmon population estimates for the American River
# Note: GrandTab data is pre-filtered for fall-run only (other runs removed upstream)
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

# ---- 6. PARALLEL PROCESSING SETUP ----
# Configure parallel processing to speed up computations
n_threads <- parallel::detectCores(logical = TRUE) - 1
n_workers <- min(6, n_threads)  # Cap at 6 physical cores
plan(multisession, workers = n_workers)
setDTthreads(n_workers)

# ---- 7. TDM MODEL VARIANTS ----
# Define three Temperature-Dependent Mortality model variants:
tdm_defs <- tribble(
  ~model,       ~calib,           ~variant,
  "exp",        "WaterForum2020", "exp_WF",      # Exponential with 2020 parameters
  "exp",        "SALMOD2006",     "exp_SM",      # Exponential with 2006 parameters  
  "lin_martin", NA,               "lin_Martin"   # Linear threshold model
)

# ---- 7.1 TDM MODEL WEIGHTS ----
# Define weights for combining the TDM model variants into an ensemble average
tdm_weights <- tribble(
  ~variant,     ~weight,
  "exp_WF",     0.51,  # Water Forum
  "exp_SM",     0.24,  # SALMOD
  "lin_Martin", 0.25   # Martin
)
stopifnot(sum(tdm_weights$weight) == 1.0) # Ensure weights sum to 1

# ---- 8. TEMPERATURE LOOKUPS ----
# Create optimized lookup structures for fast temperature data access
# These mappings allow O(1) lookup of temperatures by site and date

# Map each site to its temperature time series for each management alternative
site_temps_list <- map(env_ext_list, ~ split(.x$temp, .x$site))

# Map each site to its date sequence for each management alternative
site_dates_list <- map(env_ext_list, ~ split(.x$Date, .x$site))

# Create indexed lookups: date string -> array position for fast access
date_idx_list <- map(site_dates_list, function(dlist) {
  map(dlist, ~ setNames(seq_along(.x), as.character(.x)))
})

# ---- 9. PROCESS CARCASS DATA ----
# Transform raw carcass survey data into analysis-ready format
carcass_df <- carcass_raw %>%
  mutate(
    # Convert survey date to Date object
    Date = as.Date(surveydate),
    # Estimate spawn date as 7 days before survey (decomposition time)
    spawn_dt = Date - days(7),
    # Determine brood year (Sept-Aug cycle: Sept+ = current year, <Sept = previous year)
    brood_year = if_else(month(spawn_dt) >= 9, year(spawn_dt), year(spawn_dt) - 1), 
    # Map American River survey sections to temperature monitoring sites
    # Sections NB,W,1a,1b,2 are near the Hazel Ave gauge
    # Section 3 is near the Watt Ave gauge
    site = case_when(
      section %in% c("NB","W","1a","1b","1a/1b","2") ~ "AveHazel",
      section %in% c("3")                            ~ "AveWatt",
      TRUE                                           ~ NA_character_
    )
  ) %>%
  # Select and rename columns for analysis
  dplyr::transmute(
    brood_year, 
    spawn_dt, 
    section,      # American River section from carcass survey
    fork_length = .data[["forklength"]],
    site          # Temperature monitoring site
  )

# Calculate the proportion of observations at each site
# Used later for distributing simulated fish across sites
site_props <- carcass_df %>%
  count(site) %>%
  mutate(prop = n / sum(n))
print("Site proportions:")
print(site_props)

# ---- 10. CREATE 10-DAY SPAWNING PERIOD BINS ----
# Divide the spawning season into 10-day temporal bins for CLM analysis

# ---- 10.1 Find Anchor Date ----
# Ensure spawn_dt is properly formatted as Date
if (!inherits(carcass_df$spawn_dt, "Date")) {
  carcass_df <- carcass_df %>% mutate(spawn_dt = as.Date(spawn_dt))
}

# Find the earliest spawning date in Aug-Dec to use as season anchor
anchor_mmdd <- carcass_df %>%
  filter(!is.na(spawn_dt)) %>%
  mutate(mmdd = format(spawn_dt, "%m-%d"), m = month(spawn_dt)) %>%
  filter(m >= 8) %>%  # August or later
  summarise(anchor = min(mmdd, na.rm = TRUE)) %>%
  pull(anchor)

# Fallback if no Aug-Dec dates found
if (length(anchor_mmdd) == 0 || is.na(anchor_mmdd)) {
  anchor_mmdd <- carcass_df %>%
    filter(!is.na(spawn_dt)) %>%
    mutate(mmdd = format(spawn_dt, "%m-%d")) %>%
    summarise(anchor = min(mmdd, na.rm = TRUE)) %>%
    pull(anchor)
}

# ---- 10.2 Define Spawning Period Bins ----
bin_width <- 10L  # Each spawning period bin is 10 days wide

# Calculate season day relative to anchor and determine number of bins needed
md <- format(carcass_df$spawn_dt, "%m-%d")
season_year  <- year(carcass_df$spawn_dt) - (md < anchor_mmdd)  
season_start <- as.Date(paste0(season_year, "-", anchor_mmdd))
season_day   <- as.integer(carcass_df$spawn_dt - season_start)
needed_bins  <- ceiling((max(season_day, na.rm = TRUE) + 1) / bin_width)
n_bins <- max(12L, needed_bins)  # Ensure at least 12 bins

# Assign period to each observation
carcass_df <- carcass_df %>%
  mutate(period = assign_period(spawn_dt, anchor_mmdd, bin_width, n_bins))

# ---- 10.3 Create Bin Definitions ----
# Define the 10-day spawning period bins for the season
bin_defs <- tibble(period = paste0("p", seq_len(n_bins))) %>%
  mutate(
    # Use year 2000 as reference for bin boundaries
    start = as.Date(paste0("2000-", anchor_mmdd)) + (row_number() - 1) * bin_width,
    end   = start + (bin_width - 1)
  )

# Assign observations to spawn bins with proper year handling
carcass_df$spawn_bin <- sapply(carcass_df$spawn_dt, function(d) {
  if (is.na(d)) return(NA_character_)
  
  # Account for spawning year convention (January spawning belongs to previous year's run)
  dummy_date <- if_else(month(d) >= 9, 
                        as.Date(paste0("2000-", format(d, "%m-%d"))),
                        as.Date(paste0("2001-", format(d, "%m-%d"))))
  
  # Find matching 10-day bin
  match_idx <- which(dummy_date >= bin_defs$start & dummy_date <= bin_defs$end)
  if (length(match_idx) == 1) bin_defs$period[match_idx] else NA
})

# Convert to ordered factor
carcass_df$spawn_bin <- factor(carcass_df$spawn_bin, levels = bin_defs$period, ordered = TRUE)

# Extract unique years from observed data
yrs <- sort(unique(as.integer(carcass_df$brood_year)))

# ---- 11. VERIFY TEMPERATURE CONSISTENCY ----
# Check that historical temperatures are identical across management alternatives
# (They should be the same for 2011-2024 since alternatives only differ in future projections)
temp_consistency_check <- env_ext_list %>%
  imap_dfr(~ mutate(.x, alt = .y)) %>%
  filter(year(Date) %in% yrs) %>%
  group_by(site, Date) %>%
  summarise(n_temp = n_distinct(temp), .groups = "drop") %>%
  summarise(all_equal = all(n_temp == 1))

print(paste("Historical temperatures consistent across alternatives:", temp_consistency_check$all_equal))

# ---- 12. CALCULATE OCTOBER/NOVEMBER MEAN TEMPERATURES ----
# These months are critical for predicting spawning timing

# Collapse temperature data across alternatives (identical for historical period)
# Calculate mean temperatures for October and November by site and year
oct_nov_temps <- env_ext_list %>%
  imap_dfr(~ mutate(.x, alt = as.character(.y))) %>%
  filter(year(Date) %in% yrs, month(Date) %in% c(10, 11)) %>%
  group_by(site, year = year(Date), month = month(Date)) %>%
  summarise(mean_temp = mean(temp, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = month, values_from = mean_temp, names_prefix = "m_") %>%
  rename(Oct = m_10, Nov = m_11)

# Sanity check: ensure unique site-year combinations
stopifnot(!any(duplicated(oct_nov_temps[c("site", "year")])))

# Join historical temperatures to carcass observation data
carcass_df2 <- carcass_df %>%
  mutate(brood_year = as.integer(brood_year)) %>%
  left_join(oct_nov_temps, by = c("site", "brood_year" = "year"))

# ---- 13. STANDARDIZE TEMPERATURES ----
# Standardize Oct/Nov temperatures to z-scores for CLM model fitting
# This scaling improves model convergence and coefficient interpretability

# Calculate mean and standard deviation for standardization
oct_mean <- mean(carcass_df2$Oct, na.rm = TRUE)
oct_sd <- sd(carcass_df2$Oct, na.rm = TRUE)
nov_mean <- mean(carcass_df2$Nov, na.rm = TRUE)
nov_sd <- sd(carcass_df2$Nov, na.rm = TRUE)

# Create analysis dataset with standardized temperatures
obs_fit <- carcass_df2 %>%
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

# ---- 14. FIT CLM MODEL FOR SPAWN TIMING ----
# Fit Cumulative Link Model to predict spawning period probabilities
# based on October and November water temperatures

spawn_clm <- clm(spawn_bin ~ Oct_std + Nov_std,
                 data = obs_fit, link = "logit")
summary(spawn_clm)

# Extract model coefficients
cf   <- coef(spawn_clm)
beta <- cf[c("Oct_std","Nov_std")]  # Temperature effects
zeta <- unname(cf[!(names(cf) %in% names(beta))])  # Threshold parameters
K    <- length(zeta) + 1  # Number of 10-day spawning period categories

# Get the spawning period bin levels from the model
spawn_bins_model <- levels(spawn_clm$model$spawn_bin)
if (length(spawn_bins_model) != K) spawn_bins_model <- bin_defs$period[seq_len(K)]

# ---- 15. BUILD DETERMINISTIC FORECAST ----
# Generate future spawning patterns based on projected temperatures

# Create temperature scaling parameters using historical data
if (!exists("sc")) {
  sc <- carcass_df2 %>%
    filter(brood_year %in% yrs) %>%
    summarise(
      o_m = mean(Oct, na.rm=TRUE), o_s = sd(Oct, na.rm=TRUE),
      n_m = mean(Nov, na.rm=TRUE), n_s = sd(Nov, na.rm=TRUE)
    )
}

# Build forecast temperatures for future years
yrs_forecast   <- forecast_years
forecast_temps <- build_forecast_temps(env_ext_list, yrs_forecast, sc)
stopifnot(all(c("env","sim_year","Oct","Nov","Oct_std","Nov_std") %in% names(forecast_temps)))

# Calculate CLM spawn timing probabilities for forecast period
probs_all <- predict_clm_probs(beta, zeta, forecast_temps[, c("Oct_std","Nov_std")], offset = 0)
colnames(probs_all) <- spawn_bins_model

# Normalize probabilities to sum to 1
rs <- rowSums(probs_all)
rs[!is.finite(rs) | rs <= 0] <- 1
probs_all <- probs_all / rs

# Ensure columns match the bins actually present in observed data
present_bins <- levels(droplevels(carcass_df$spawn_bin))
if (!identical(colnames(probs_all), present_bins)) {
  keep <- match(present_bins, colnames(probs_all))
  probs_all <- probs_all[, keep, drop = FALSE]
  colnames(probs_all) <- present_bins
}

# ---- 16. JOIN CARCASS DATA WITH MANAGEMENT ALTERNATIVES ----
# Map observed carcass data to management alternatives based on site

# Join carcass observations with management alternative mapping
carcass_df_mgt <- carcass_df %>% left_join(
  purrr::imap_dfr(env_ext_list, ~ tibble(mgt_alt = as.character(.y), site = unique(.x$site))) %>%
    distinct(site, mgt_alt),
  by = "site"
)

# Calculate average number of redds per management alternative
mgt_alt_nredd <- carcass_df_mgt %>%
  filter(brood_year %in% as.integer(2011:2024), !is.na(mgt_alt)) %>%
  count(mgt_alt, brood_year, name = "n") %>%
  group_by(mgt_alt) %>%
  summarise(N_redd = max(1L, round(mean(n, na.rm = TRUE))), .groups = "drop")

# Create sampling pools for each management alternative
# These pools are used to randomly sample characteristics for simulated fish
mgt_alt_pools <- carcass_df_mgt %>%
  filter(brood_year %in% as.integer(2011:2024), !is.na(mgt_alt)) %>%
  group_by(mgt_alt) %>%
  summarise(
    section_pool     = list(section),      # Pool of carcass survey sections
    fork_length_pool = list(fork_length),  # Pool of observed fork lengths
    site_pool        = list(unique(site)), # Pool of temperature monitoring sites
    .groups = "drop"
  )

# Join forecast temperatures with redd counts and sampling pools
ft_joined <- forecast_temps %>%
  left_join(mgt_alt_nredd, by = c("env" = "mgt_alt")) %>%
  left_join(mgt_alt_pools, by = c("env" = "mgt_alt"))

# ---- 17. SIMULATE FUTURE SPAWNING DATA ----
# Generate synthetic spawning observations for future water years
# These simulated observations are used to calculate TDM survival rates

out_list <- vector("list", nrow(ft_joined))

for (i in seq_len(nrow(ft_joined))) {
  # Extract parameters for this management alternative and year
  mgt_alt_i <- ft_joined$env[i]
  yr_i  <- ft_joined$sim_year[i]
  n_i   <- ft_joined$N_redd[i]
  if (is.na(n_i) || n_i <= 0) n_i <- 1L  # Ensure at least 1 redd
  
  # Sample 10-day spawning bins based on temperature-driven probabilities
  p_i <- probs_all[i, ]  # Spawn timing probabilities from CLM model
  bins_i  <- sample(present_bins, n_i, replace = TRUE, prob = p_i)
  
  # Sample specific dates within selected 10-day bins
  # FIX: Pass the bin_defs data frame to the function
  dates_i <- sample_dates_fast(bins_i, yr_i, bin_defs = bin_defs)
  
  # Get sampling pools for this management alternative
  # Pool of carcass survey sections
  sec_pool <- ft_joined$section_pool[[i]]
  if (is.null(sec_pool)) sec_pool <- NA
  
  # Pool of observed fork lengths
  fl_pool  <- ft_joined$fork_length_pool[[i]]
  if (is.null(fl_pool)) fl_pool <- NA_real_
  
  # Pool of temperature monitoring sites (AveHazel, AveWatt)
  st_pool  <- ft_joined$site_pool[[i]]
  if (is.null(st_pool)) {
    st_pool <- carcass_df_mgt %>% filter(mgt_alt == mgt_alt_i) %>% pull(site) %>% unique()
  }
  
  # Create simulated spawning observations
  # These represent individual redds with spawn dates and locations
  out_list[[i]] <- tibble(
    mgt_alt     = mgt_alt_i,
    sim_year    = yr_i,
    spawn_dt    = dates_i,
    section     = sample(sec_pool, n_i, replace = TRUE),
    fork_length = sample(fl_pool,  n_i, replace = TRUE),
    site        = sample(st_pool,  n_i, replace = TRUE)
  )
}

# Combine all simulated data for future water years
sim_future <- dplyr::bind_rows(out_list)

# ---- 18. VISUALIZE CLM PREDICTIONS ----
# Create plots showing how temperature affects spawning timing probabilities

# ---- 19.1 Prepare Prediction Data ----
# Create pretty labels for bins (e.g., "p1: Oct 05–Oct 14")
idx <- match(spawn_bins_model, bin_defs$period)
pretty_labels <- paste0(
  spawn_bins_model, ": ",
  format(bin_defs$start[idx], "%b %d"), "–", format(bin_defs$end[idx], "%b %d")
)

rng_oct <- safe_range(obs_fit$Oct_std)
rng_nov <- safe_range(obs_fit$Nov_std)

# Create sequences for spawn timing predictions
oct_range <- seq(rng_oct[1], rng_oct[2], length.out = 100)
nov_range <- seq(rng_nov[1], rng_nov[2], length.out = 200)

# ---- 19.2 October Temperature Effects on Spawn Timing ----
newdata_oct <- data.frame(
  Oct_std = oct_range,
  Nov_std = mean(obs_fit$Nov_std, na.rm = TRUE)
)
pred_oct <- predict_clm_probs(beta, zeta, newdata_oct, offset = 0)
colnames(pred_oct) <- spawn_bins_model

# Convert to long format for plotting
pred_oct_long <- as.data.frame(pred_oct) %>%
  mutate(Oct_std = oct_range) %>%
  pivot_longer(all_of(spawn_bins_model), names_to = "spawn_bin", values_to = "probability") %>%
  mutate(spawn_bin = factor(spawn_bin, levels = spawn_bins_model, labels = pretty_labels))

# ---- 19.3 November Temperature Effects on Spawn Timing ----
newdata_nov <- data.frame(
  Oct_std = mean(obs_fit$Oct_std, na.rm = TRUE),
  Nov_std = nov_range
)
pred_nov <- predict_clm_probs(beta, zeta, newdata_nov, offset = 0)
colnames(pred_nov) <- spawn_bins_model

# Convert to long format for plotting
pred_nov_long <- as.data.frame(pred_nov) %>%
  mutate(Nov_std = nov_range) %>%
  pivot_longer(all_of(spawn_bins_model), names_to = "spawn_bin", values_to = "probability") %>%
  mutate(spawn_bin = factor(spawn_bin, levels = spawn_bins_model, labels = pretty_labels))

# ---- 19.4 Create Plots ----
# Plot 1: October temperature effects with facets
p_oct_facet <- ggplot(pred_oct_long, aes(x = Oct_std, y = probability, color = spawn_bin)) +
  geom_line(size = 1) +
  facet_wrap(~ spawn_bin, ncol = 3) +
  scale_color_viridis_d(guide = "none") +
  labs(
    title = "Spawn timing probability vs standardized October temperature",
    subtitle = "Each panel shows probability for a different 10-day spawning period",
    x = "Standardized October Temperature",
    y = "Probability"
  ) +
  theme_minimal(base_size = 14)

# Plot 2: November temperature effects with facets
p_nov_facet <- ggplot(pred_nov_long, aes(x = Nov_std, y = probability, color = spawn_bin)) +
  geom_line(size = 1) +
  facet_wrap(~ spawn_bin, ncol = 3) +
  scale_color_viridis_d(guide = "none") +
  labs(
    title = "Spawn timing probability vs standardized November temperature",
    subtitle = "Each panel shows probability for a different 10-day spawning period",
    x = "Standardized November Temperature",
    y = "Probability"
  ) +
  theme_minimal(base_size = 14)

# Display the plots
print(p_oct_facet)
print(p_nov_facet)

# ---- 20. PREPARE FORECAST DATA FOR COMPARISON ----
# Add spawning year-aligned dates and source label to forecast data
# Spawning year runs Sept-Aug, so dates are aligned to this biological cycle
fc_df <- sim_future %>%
  mutate(
    season_date = season_posix(spawn_dt, anchor_mmdd),
    source = "Forecast"
  )

# ---- 21. COMPARE OBSERVED VS MODELED SPAWNING PATTERNS ----
# This section validates the CLM model by comparing actual vs predicted spawn timing

# Step 1: Standardize temperatures using empirical data statistics
sc <- carcass_df2 %>%
  filter(brood_year %in% yrs) %>%
  summarise(
    o_m = mean(Oct, na.rm=TRUE), o_s = sd(Oct, na.rm=TRUE),
    n_m = mean(Nov, na.rm=TRUE), n_s = sd(Nov, na.rm=TRUE)
  )

# Calculate yearly mean temperatures and standardize
hist_temps <- carcass_df2 %>%
  filter(brood_year %in% yrs) %>%
  group_by(brood_year) %>%
  summarise(Oct = mean(Oct, na.rm=TRUE), Nov = mean(Nov, na.rm=TRUE), .groups="drop") %>%
  mutate(
    Oct_std = (Oct - sc$o_m)/sc$o_s,
    Nov_std = (Nov - sc$n_m)/sc$n_s
  ) %>%
  rename(sim_year = brood_year)

# Step 2: Generate CLM probabilities for empirical years
cf   <- coef(spawn_clm)
beta <- cf[c("Oct_std","Nov_std")]
zeta <- unname(cf[!(names(cf) %in% names(beta))])
bins <- levels(spawn_clm$model$spawn_bin)

# Calculate spawn timing probabilities for empirical years
probs_hist <- predict_clm_probs(beta, zeta, hist_temps[, c("Oct_std","Nov_std")])
colnames(probs_hist) <- bins

# Step 3: Sample modeled dates matching observed counts per year
# This creates synthetic spawn dates based on CLM predictions to validate the model
n_by_year <- carcass_df %>% 
  filter(brood_year %in% yrs) %>% 
  count(brood_year, name="N")

set.seed(42)  # For reproducibility of sampling
sim_pred <- lapply(seq_len(nrow(hist_temps)), function(i){
  yr <- hist_temps$sim_year[i]
  # Get actual count for this year
  n  <- n_by_year$N[n_by_year$brood_year==yr]
  if (length(n)==0 || is.na(n) || n<=0) n <- 1L
  
  # Sample 10-day bins based on model probabilities
  p  <- probs_hist[i,]
  p <- p/sum(p)  # Ensure normalized
  b  <- sample(bins, n, replace=TRUE, prob=p)
  
  # Sample dates within selected bins
  tibble(
    sim_year = yr, 
    spawn_dt = sample_dates_fast(b, yr, bin_defs = bin_defs), 
    source = "Modeled"
  )
}) %>% bind_rows()

# Step 4: Extract observed carcass dates
sim_actual <- carcass_df %>% 
  filter(brood_year %in% yrs) %>%
  transmute(
    sim_year = brood_year, 
    spawn_dt, 
    source = "Observed"
  )

# ---- 22. CREATE COMPARISON VISUALIZATIONS ----

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

# Combine observed and modeled data with season alignment
comp_df <- bind_rows(sim_actual, sim_pred) %>%
  mutate(season_date = as.Date(season_posix(spawn_dt, anchor_mmdd)))

# Set up x-axis breaks
rng <- range(comp_df$season_date, na.rm=TRUE)
x_breaks <- seq(rng[1], rng[2], by = "10 days")

# Create ridge plot showing density distributions by year
p_ridge <- ggplot(comp_df, aes(x=season_date, y=factor(sim_year))) +
  geom_density_ridges_gradient(
    aes(height=after_stat(density), fill=source),
    scale=2, rel_min_height=0.001, linewidth=0.2, bandwidth=10, trim=TRUE
  ) +
  scale_x_date(breaks=x_breaks, date_labels="%b %d", expand=c(0.01,0)) +
  scale_fill_manual(values=c(Observed="steelblue", Modeled="tomato")) +
  labs(
    title="Observed vs Modeled spawn timing by year (2011–2024)",
    x=paste0("Spawning year day (anchor ", anchor_mmdd, ")"), 
    y="Year", 
    fill=NULL
  ) +
  theme_minimal(base_size=13)
print(p_ridge)

# (B) Calendar-date boxplot comparing medians and distributions
p_box <- ggplot(bind_rows(sim_actual, sim_pred),
                aes(x=factor(sim_year), y=spawn_dt, fill=source)) +
  geom_boxplot(alpha=0.65, position=position_dodge(width=0.7), outlier.alpha=0.3) +
  scale_fill_manual(values=c(Observed="steelblue", Modeled="tomato")) +
  labs(
    title="Observed vs Modeled spawning dates by year (2011–2024)",
    x="Brood year", y="Spawn date", fill=NULL
  ) +
  theme_minimal(base_size=13)
print(p_box)

# ---- 23. BUILD FINAL SIMULATION DATASETS ----
# Combine observed and simulated data for Temperature-Dependent Mortality modeling

# Combine observed (historical) and simulated (future) spawning data
sim_redds <- bind_rows(
  sim_actual %>% mutate(mgt_alt = NA_character_, site = NA_character_),
  sim_future %>% dplyr::select(mgt_alt, sim_year, spawn_dt, site)
) %>%
  arrange(mgt_alt, sim_year, spawn_dt)

# Convert to data.table for efficient processing
if (!data.table::is.data.table(sim_redds)) data.table::setDT(sim_redds)
data.table::setkey(sim_redds, mgt_alt, sim_year)

# Split by management alternative and year for parallel processing
sim_redds_split <- split(sim_redds, by = c("mgt_alt","sim_year"), drop = TRUE, keep.by = FALSE)

# Keep only essential columns in each chunk
sim_redds_split <- lapply(sim_redds_split, function(dt) dt[, .(spawn_dt, site)])

# Calculate median spawn dates per management alternative and year
spawn_medians_alt_year <- sim_redds %>%
  as_tibble() %>%
  filter(!is.na(mgt_alt)) %>%
  group_by(mgt_alt, sim_year) %>%
  summarise(
    spawn_dt = as.Date(median(as.numeric(spawn_dt), na.rm = TRUE), origin = "1970-01-01"),
    .groups = "drop"
  )

# Get unique management alternatives
mgt_alts <- sort(unique(spawn_medians_alt_year$mgt_alt))

# Build spawn date vectors for each management alternative
spawn_dates_by_alt <- map(
  .x = set_names(mgt_alts),
  .f = ~ build_spawn_vec_for_env(
    df_env    = filter(spawn_medians_alt_year, mgt_alt == .x),
    sim_years = sim_years
  )
)

# ---- 24. VISUALIZE FUTURE SPAWNING DISTRIBUTIONS ----
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
p1 <- ggplot(sim_future_comp, aes(x = dummy_date, color = factor(sim_year), group = sim_year)) +
  geom_density(alpha = 0.3) +
  labs(
    title = "Projected Spawning Date Distributions (all years aligned to same season)",
    x = "Spawning Date (Oct–Jan, reference year)",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)
print(p1)

# Calculate summary statistics for spawn dates
summary_dates <- sim_future_comp %>%
  dplyr::group_by(sim_year) %>%
  dplyr::summarize(
    median_spawn = as.Date(median(as.numeric(spawn_dt), na.rm = TRUE), origin = "1970-01-01"),
    p10 = as.Date(quantile(as.numeric(spawn_dt), 0.10, na.rm = TRUE), origin = "1970-01-01"),
    p90 = as.Date(quantile(as.numeric(spawn_dt), 0.90, na.rm = TRUE), origin = "1970-01-01"),
    .groups = "drop"
  )

# Plot 2: Time series of median spawn date with uncertainty bands
p2 <- ggplot(summary_dates, aes(x = sim_year)) +
  geom_line(aes(y = median_spawn), color = "blue") +
  geom_ribbon(aes(ymin = p10, ymax = p90), alpha = 0.2, fill = "blue") +
  labs(
    title = "Projected Median Spawn Date and 80% Confidence Interval",
    x = "Simulation Year",
    y = "Spawn Date"
  ) +
  theme_minimal(base_size = 14)
print(p2)

# ---- 25. CONFIGURE PARALLEL PROCESSING FOR TDM CALCULATIONS ----

# Increase memory limit for parallel workers (4 GiB)
options(future.globals.maxSize = 4 * 1024^3)

# Reset parallel plan with appropriate worker count
plan(multisession, workers = min(4, n_workers))

# Extract real median spawn dates from observed data
spawn_dates_real <- sim_redds %>%
  filter(!is.na(sim_year)) %>%
  arrange(sim_year) %>%
  group_by(sim_year) %>%
  summarize(spawn_dt = median(spawn_dt, na.rm=TRUE), .groups="drop") %>%
  pull(spawn_dt)

# Repeat the month-day patterns to fill the entire simulation period
spawn_dm <- rep(spawn_dates_real, length.out = n_sim)

# Generate simulation year sequence
sim_years_vec <- (real_years[1] + seq_len(n_sim) - 1)

# ---- 26. TEMPERATURE-DEPENDENT MORTALITY (TDM) CALCULATIONS ----
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

# Step 2: Cache management alternative assets for fast access
alt_cache <- lapply(names(env_ext_list), function(alt_nm) {
  list(
    alt_nm       = alt_nm,
    date_idx_alt = date_idx_list[[alt_nm]],
    temps_alt    = site_temps_list[[alt_nm]],
    date_min     = min(env_ext_list[[alt_nm]]$Date),
    date_max     = max(env_ext_list[[alt_nm]]$Date)
  )
})
names(alt_cache) <- names(env_ext_list)

# Ensure TDM definitions use character types (not factors)
tdm_defs <- tdm_defs %>%
  dplyr::mutate(
    model   = as.character(model),
    calib   = as.character(calib),
    variant = as.character(variant)
  )

# Continue with remaining sections...
print("TDM calculations configured. Ready for parallel processing.")

# ---- 27. DATA QUALITY CONTROL - REMOVE INVALID SITES ----
# Ensure we only use sites that have temperature data across all management alternatives.
# This prevents errors during the TDM calculations.

# Create a cache for environment-specific assets to speed up access.
# This avoids repeatedly accessing the large env_ext_list.
env_cache <- alt_cache

# Find sites that have temperature data in ALL management alternatives
sites_with_temps <- Reduce(
  intersect,
  lapply(env_cache, function(ec) names(ec$date_idx_alt))
)

# Filter the split redd data to only include redds from these common, valid sites.
sim_redds_split <- lapply(
  sim_redds_split,
  function(dt) {
    if (is.null(dt) || !nrow(dt)) return(dt)
    # Ensure required columns exist before filtering
    stopifnot(all(c("site", "spawn_dt") %in% names(dt)))
    # Keep only rows where the site is in our list of valid sites
    dt[dt$site %in% sites_with_temps, , drop = FALSE]
  }
)

# ---- 28. BUILD LOOKUP TABLES FOR FAST DATE-POSITION MAPPING ----
# Create indexed lookup tables for efficient date-to-array-position conversion.
# Using a keyed data.table provides near-instantaneous (O(1)) lookups, which is a
# major performance improvement over repeated searching.

env_lookup <- lapply(env_cache, function(ec) {
  # Build a data.table with site, redd date (rdr), and position
  DT <- rbindlist(
    lapply(names(ec$date_idx_alt), function(s) {
      v <- ec$date_idx_alt[[s]]
      if (is.null(v) || !length(v)) {
        # Return an empty table with the correct structure if no data exists
        data.table(
          site = factor(character(0)),
          rdr = as.IDate(integer(0)),
          pos = integer(0)
        )
      } else {
        data.table(
          site = s,
          rdr  = as.IDate(names(v)),  # Use fast IDate for integer-based date
          pos  = as.integer(v)         # The position in the temperature array
        )
      }
    }),
    use.names = TRUE,
    fill = TRUE
  )
  
  # Convert site to factor and set a two-column key for rapid lookups
  DT[, site := as.factor(site)]
  setkey(DT, site, rdr)
  DT
})

# ---- 29. IMPLEMENT MEMOIZATION FOR SURVIVAL CALCULATIONS ----
# Cache survival calculations (memoization) to avoid redundant computations, as many
# redds share the same site, spawn date, and model variant.

# Create a dedicated environment to act as a cache
.mem_env <- new.env(parent = emptyenv())

# ---- 30. OPTIMIZE AND EXECUTE TDM CALCULATIONS IN PARALLEL ----
# Run the main TDM calculation loop in parallel for maximum speed.

# Byte-compile hot-path functions for a potential speed boost
if (requireNamespace("compiler", quietly = TRUE)) {
  compute_surv_vec    <- compiler::cmpfun(compute_surv_vec)
  pairs_for_env_year  <- compiler::cmpfun(pairs_for_env_year)
  eval_year           <- compiler::cmpfun(eval_year)
}

# Set data.table threads low to avoid contention with parallel futures
data.table::setDTthreads(2L)

# Set up the parallel execution plan using 'future'
plan(multisession, workers = max(1L, parallel::detectCores() - 1L))

# Execute TDM calculations in parallel across all simulation years
results_obs_fast <- furrr::future_map_dfr(
  sim_years,
  ~eval_year(.x, sim_redds_split, env_cache, tdm_defs),
  .options = furrr::furrr_options(seed = TRUE, scheduling = 20)
)
# Result contains: sim_year, mgt_alt, variant, method, mean_cum_surv

# ---- 31. LIFE CYCLE MODEL CALIBRATION AND FORECASTING ----
# Prepare for and run the full life-cycle model, using the TDM survival rates
# calculated above as a key input.

# Select a reference management alternative for calibration
ref_env <- names(env_ext_list)[1]

# ---- 31.1 Summarize TDM Survival Results ----
egg_summary <- results_obs_fast %>%
  rename(env = mgt_alt) %>% # Rename for consistency with downstream code
  arrange(env, variant, sim_year) %>%
  group_by(env, variant, sim_year) %>%
  summarise(mean_cum_surv = mean(mean_cum_surv, na.rm = TRUE), .groups = "drop")

# Calculate weighted average survival across variants
egg_summary_weighted <- egg_summary %>%
  left_join(tdm_weights, by = "variant") %>%
  # Ensure only variants with weights are included
  filter(!is.na(weight)) %>%
  group_by(env, sim_year) %>%
  summarise(
    mean_cum_surv = sum(mean_cum_surv * weight, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(variant = "weighted_avg") # Assign new variant name

# Combine original and weighted results
egg_summary <- bind_rows(egg_summary, egg_summary_weighted)

variant_names <- egg_summary %>% pull(variant) %>% unique() %>% sort()

# ---- 31.2 Create Survival Lookup Structures ----
# For reference environment (used for calibration)
surv_lookup_by_variant <- egg_summary %>%
  filter(env == ref_env, sim_year %in% real_years) %>%
  arrange(variant, sim_year) %>%
  group_by(variant) %>%
  summarise(surv_vec = list(mean_cum_surv), .groups = "drop") %>%
  tibble::deframe() # Converts to a named list: variant -> numeric vector

# Full lookup for all environment-variant combinations (for forecasting)
surv_lookup_full <- egg_summary %>%
  arrange(env, variant, sim_year) %>%
  group_by(env, variant) %>%
  summarise(surv_vec = list(mean_cum_surv), .groups = "drop") %>%
  mutate(key = paste(env, variant, sep = "_")) %>%
  dplyr::select(key, surv_vec) %>%
  tibble::deframe() # Converts to a named list: "env_variant" -> survival vector

# ---- 31.3 Define Base Life-Cycle Parameters ----
base_P <- list(
  female_fraction = 0.5,
  fec = 5522,               # Fecundity: eggs per female
  S0 = 0.347,               # Base egg survival (pre-TDM)
  K_spawners = 12493,       # Spawner carrying capacity
  SAR_mean = NA_real_,      # Smolt-to-Adult Return rate (to be calibrated)
  SAR_sd = 0.00237,         # Std dev of SAR
  lag_probs = c(`3` = 0.828982, `4` = 0.168885, `5` = 0.002105),
  rear_surv = NA_real_      # Rearing survival (to be calibrated)
)

# ---- 31.4 Set Up Calibration Data ----
obs_spawners  <- esc_obs$spawners
S_seed_calib  <- obs_spawners[1:3]        # Seed population (first 3 years)
n_calib       <- length(obs_spawners)
fit_idx       <- (length(S_seed_calib) + 1):n_calib # Years to fit model (yr 4-14)

# ---- 32. CALIBRATE LIFE-CYCLE PARAMETERS ----
# Calculate degree-days for adult spawners for calibration
ref_env <- names(env_ext_list)[1]
deg_day_cal_ref <- deg_day_cal_for(ref_env)
stopifnot(length(deg_day_cal_ref) == length(real_years))

# Run calibration in parallel for each TDM variant
calib_results <- furrr::future_map_dfr(
  variant_names, # This will be c("exp_SM", "exp_WF", "lin_Martin")
  function(v) {
    opt <- optim(
      par    = c(0.0025, 0.5419),
      fn     = modular_sse,
      variant= v,
      method = "L-BFGS-B",
      lower  = c(0, 0), upper  = c(1, 1)
    )
    tibble::tibble(
      variant   = v,
      SAR_mean  = opt$par[1],
      rear_surv = opt$par[2],
      sse       = opt$value
    )
  },
  .options = furrr::furrr_options(seed = TRUE)
)

# Create a nested list of parameter sets for every variant-environment combination
base_P_list <- calib_results %>%
  split(.$variant) %>%
  purrr::map(function(df_v) {
    SARv  <- df_v$SAR_mean[1]
    rearv <- df_v$rear_surv[1]
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

# ---- 33. GENERATE CALIBRATION PREDICTIONS FOR VALIDATION ----
# Run the model over the historical period using the calibrated parameters
# to visualize how well the model fits the observed data.

calib_pred_by_variant <- rlang::set_names(
  lapply(variant_names, function(v) {
    P0 <- base_P_list[[v]][[ref_env]]
    # Get the correct survival vector (now includes weighted_avg)
    surv_vec <- surv_lookup_by_variant[[v]][1:n_calib]
    # Run simulation for the calibration period
    out <- simulate_variant(
      surv_vec       = surv_vec,
      P              = P0,
      years          = n_calib,
      S_init         = S_seed_calib,
      SAR_vec        = rep(P0$SAR_mean, n_calib),
      K_spawners_vec = rep(P0$K_spawners, n_calib),
      deg_day_adult  = deg_day_cal_ref,
      sim_years_vec  = real_years
    )
    # Combine observations and predictions for plotting
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

# ---- 34. PREPARE SEED POPULATIONS FOR FORECASTING ----
# Use the last 3 years of the calibrated simulation run as the initial "seed"
# population for the forecasts to ensure a smooth transition.

k <- length(S_seed_calib)  # Number of seed years (3)

S_seed_fore_list <- rlang::set_names(
  purrr::map(calib_results$variant, function(v) {
    years_cal <- length(real_years)
    surv_vec  <- surv_lookup_full[[paste(ref_env, v, sep = "_")]][1:years_cal]
    Ptmp <- base_P_list[[v]][[ref_env]]
    deg_day_cal <- compute_deg_day_adult(
      env_nm       = ref_env,
      sim_years    = real_years,
      spawn_dates  = spawn_dates_by_alt[[ref_env]][match(real_years, sim_years)],
      env_ext_list = env_ext_list
    )
    # Run the full calibration simulation one last time
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
    # Return the last 3 years of the simulation as the seed
    tail(out$spawners, k)
  }),
  calib_results$variant
)

# ---- 35. RUN POPULATION FORECASTS ----
# Project populations 100 years forward for every combination of management
# alternative and TDM model variant.

keys <- names(surv_lookup_full)
use_stochastic_SAR <- FALSE # Can be toggled in the Shiny app

# Define a list of stochastic SAR options for potential UI control
stoch_SAR_opts <- list(
  model       = "normal",
  mean        = base_P$SAR_mean,
  sd          = base_P$SAR_sd,
  shape1      = 2,
  shape2      = 5,
  timing      = "pulse",
  block_years = 20:30,
  pulse_years = c(10, 15, 20, 25, 30, 35, 40),
  pulse_sd    = 0.002
)

# Run forecasts for all environment-variant combinations
results_full <- purrr::map_dfr(keys, function(key) {
  # Parse key to get environment and variant names
  parts  <- strsplit(key, "_")[[1]]
  env_nm <- parts[1]
  var_nm <- paste(parts[-1], collapse = "_")
  seed_vec <- S_seed_fore_list[[var_nm]]
  # `sim_forecast_fn` is a helper defined in functions.R
  sim_forecast_fn(
    var_nm,
    env_nm,
    flow_cfs = NULL, # Flow covariates not used in this version
    S_seed   = seed_vec,
    spawn_dates_by_alt = spawn_dates_by_alt
  )()
})

# ---- 36.STEELHEAD PERFORMANCE METRIC CALCULATION ----
# This metric is the number of days below 18.3°C in Oct/Nov for each alternative.

# Calculate the metric for each of the 28 alternatives
steelhead_metrics <- df_all_orig %>%
  # 1. Filter for the relevant time period first
  filter(
    month(Date) %in% c(10, 11),      # October and November
    year(Date) == 2025               # Standard forecast year
  ) %>%
  # 2. For each day in each alternative, calculate the average temperature across sites
  group_by(env, Date) %>%
  summarise(avg_temp = mean(temp, na.rm = TRUE), .groups = "drop") %>%
  # 3. For each alternative, count how many days had an average temp below the threshold
  group_by(env) %>%
  summarise(
    steelhead_score = sum(avg_temp < 18.3, na.rm = TRUE)
  ) %>%
  ungroup()

# ---- 37. SAVE ALL OUTPUTS FOR SHINY APPLICATION ----
# Save all processed data frames and model objects as .rds files. These files
# will be loaded directly by the Shiny dashboard for fast startup.

saveRDS(calib_results,         here("SalmonCountR","app_data","calib_results.rds"))
saveRDS(calib_pred_by_variant, here("SalmonCountR","app_data","calib_pred_by_variant.rds"))
saveRDS(results_full,          here("SalmonCountR","app_data","results_full.rds"))
saveRDS(egg_summary,           here("SalmonCountR","app_data","egg_summary.rds"))
saveRDS(surv_lookup_full,      here("SalmonCountR","app_data","surv_lookup_full.rds"))
saveRDS(base_P_list,           here("SalmonCountR","app_data","base_P_list.rds"))
saveRDS(base_P,              here("SalmonCountR","app_data","base_P.rds"))
saveRDS(S_seed_calib,          here("SalmonCountR","app_data","S_seed_calib.rds"))
saveRDS(S_seed_fore_list,      here("SalmonCountR","app_data","S_seed_fore_list.rds"))
saveRDS(stoch_SAR_opts,        here("SalmonCountR","app_data","stoch_SAR_opts.rds"))
saveRDS(sim_years,             here("SalmonCountR","app_data","sim_years.rds"))
saveRDS(spawn_dates_by_alt,    here("SalmonCountR","app_data","spawn_dates_by_alt.rds"))
saveRDS(steelhead_metrics, here("SalmonCountR","app_data", "steelhead_metrics.rds"))

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                             END OF PRECOMPUTE SCRIPT                          ║
# ╚══════════════════════════════════════════════════════════════════════════════╝
