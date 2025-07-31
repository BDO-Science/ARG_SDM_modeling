#!/usr/bin/env Rscript

# ───────────────────────────────────────────────────────────────────────────────
# scalar_by_env_variant.R
#
# 1. Load per-env, per-variant egg→fry survival from `etf_survival_actual_temps.R`
# 2. Build a surv0 lookup for every (env, variant)
# 3. Define and run a density-dependent scalar population model
# 4. Compute %-difference from the best variant within each env
# 5. Plot trajectories and end-of-simulation spawners
# ───────────────────────────────────────────────────────────────────────────────

# 0) CLEAN WORKSPACE & LOAD LIBRARIES
#    - clear all variables
rm(list = ls())

#    - load tidyverse for data manipulation & ggplot2, scales for formatting
library(tidyverse)
library(scales)

# 1) LOAD YOUR TDM OUTPUTS ------------------------------------------------------
#    We assume `etf_survival_actual_temps.R` sets up an object `results_multi_fast`
egg_env     <- new.env()
sys.source("etf_survival_actual_temps.R", envir = egg_env)
egg_summary <- egg_env$results_obs_fast
# egg_summary should contain columns: env, variant, mean_cum_surv, year

# 2) BUILD surv0 LOOKUP (one value per env×variant) -----------------------------
#    For each env/variant, take the maximum mean_cum_surv as your base survival
surv_lookup_env <- egg_summary %>%
  group_by(env, variant) %>%
  summarise(
    surv0 = max(mean_cum_surv, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # combine env and variant into a single key for mapping
  mutate(alt_variant = paste(env, variant, sep = "_"))

#    Turn that into a named vector: names() are "env_variant", values are surv0
surv_lookup <- deframe(surv_lookup_env %>% select(alt_variant, surv0))

# -------------------------------------------------------------------------------
# 3) SET PARAMETERS FOR THE SCALAR POPULATION MODEL
# -------------------------------------------------------------------------------
base_P <- list(
  female_fraction = 0.5,      # fraction of spawners that are female
  fec             = 5522,     # eggs per female
  S0              = 0.347,    # Beverton–Holt asymptotic survival at low density
  K_spawners      = 10000,    # half‐saturation constant for redd density
  S_rear          = 0.8,      # rearing survival rate (no overflow modeled)
  SAR_mean        = 0.00308,  # mean smolt→adult return rate
  SAR_sd          = 0.002,    # SD for interannual variability in SAR
  lag_probs       = c(`3` = 0.75, `4` = 0.249, `5` = 0.001)  # age structure
)

# -------------------------------------------------------------------------------
# 4) DEFINE simulate_variant(): density‐dependent life‐cycle with external SAR
# -------------------------------------------------------------------------------
simulate_variant <- function(
    surv_vec,    # vector of egg→fry survival rates per year
    P,           # parameter list
    years,       # total simulation years
    S_init,      # initial seed spawners for the first few years
    SAR_vec = rep(P$SAR_mean, years)  # external SAR time series (default constant)
) {
  # replicate surv_vec if only one value is provided
  if (length(surv_vec) == 1) surv_vec <- rep(surv_vec, years)
  
  n_init <- length(S_init)       # how many years are pre‐seeded
  
  # Preallocate diagnostic and state vectors
  S        <- numeric(years)     # adult spawners each year
  reared   <- numeric(years)     # smolts leaving rearing each year
  dd_vec   <- numeric(years)     # density‐dependence factors
  fry_dd   <- numeric(years)     # fry after TDM and DD
  SAR_used <- numeric(years)     # SAR actually used each year
  
  # Seed the first n_init years with provided spawner counts
  S[1:n_init] <- S_init
  
  # Main loop over years
  for (t in seq_len(years)) {
    # 1) Spawn → Eggs
    redds <- S[t] * P$female_fraction      # number of nests
    eggs  <- redds * P$fec                  # total eggs
    
    # 2) Beverton–Holt density dependence
    dd        <- P$S0 / (1 + redds / P$K_spawners)
    fry_dd[t] <- eggs * surv_vec[t] * dd    # fry after TDM×DD
    dd_vec[t] <- dd
    
    # 3) Rearing → Smolts (no explicit overflow)
    reared[t] <- fry_dd[t] * P$S_rear
    
    # 4) Marine survival & age‐structured returns
    #    Use the external SAR_vec[t] rather than a fixed SAR
    SAR_t       <- SAR_vec[t]
    SAR_used[t] <- SAR_t
    
    # distribute returning adults at ages 3–5
    for (age in 3:5) {
      ret_year <- t + age
      if (ret_year <= years) {
        S[ret_year] <- S[ret_year] +
          reared[t] * SAR_t * P$lag_probs[as.character(age)]
      }
    }
  }
  
  # Return diagnostics as a tibble
  tibble(
    year      = seq_len(years),
    spawners  = S,
    dd        = dd_vec,
    fry_dd    = fry_dd,
    eff_surv  = surv_vec * dd_vec,  # combined egg→fry and DD survival
    SAR_used  = SAR_used            # actual SAR applied each year
  )
}

# ------------------------------------------------------------------------------
# 5) CONTROL: deterministic vs. stochastic SAR
#    - If TRUE, draw one SAR per year; if FALSE, use a constant SAR each year
# ------------------------------------------------------------------------------
use_stochastic_SAR <- TRUE

# spin‑up parameters
years   <- 50
S_seed  <- c(16383, 37252, 45541)

# Build shared SAR_by_year ------------------------------------------------------
if (use_stochastic_SAR) {
  # stochastic: one draw from N(mean, sd) per year, floored at zero
  SAR_by_year <- pmax(
    rnorm(n = years, mean = base_P$SAR_mean, sd = base_P$SAR_sd),
    0
  )
} else {
  # deterministic: constant SAR every year
  SAR_by_year <- rep(base_P$SAR_mean, years)
}

# Optional: inspect the SAR time series
print(tibble(year = 1:years, SAR = SAR_by_year))

# ------------------------------------------------------------------------------
# 6) RUN THE LIFE‑CYCLE SIM FOR EACH ENV×VARIANT
#    - uses shared SAR_by_year vector
# ------------------------------------------------------------------------------
results_full <- map_dfr(names(surv_lookup), function(key) {
  simulate_variant(
    surv_vec = surv_lookup[[key]],  # TDM egg→fry survival
    P        = base_P,
    years    = years,
    S_init   = S_seed,
    SAR_vec  = SAR_by_year           # external SAR schedule
  ) %>%
    mutate(variant_full = key)      # preserve which alt
})

# ------------------------------------------------------------------------------
# 7) SPLIT variant_full into env & variant, select core columns
# ------------------------------------------------------------------------------
results <- results_full %>%
  tidyr::extract(
    variant_full,
    into  = c("env","variant"),
    regex = "^([^_]+)_(.+)$"
  ) %>%
  mutate(alt_variant = variant) %>%
  select(env, variant, alt_variant, year, spawners)

# ------------------------------------------------------------------------------
# 8) COMPUTE PERCENT DIFFERENCES WITHIN EACH ENV × YEAR
# ------------------------------------------------------------------------------
pct_results <- results %>%
  group_by(env, year) %>%
  mutate(
    max_spawners = max(spawners),
    pct_diff     = (spawners - max_spawners) / max_spawners * 100
  ) %>%
  ungroup() %>%
  select(env, alt_variant, variant, year, spawners, pct_diff)

# Print end-of-simulation comparisons
last_year <- max(pct_results$year)
cat("\n% differences from best variant in year", last_year, ":\n")
pct_results %>%
  filter(year == last_year) %>%
  arrange(env, desc(pct_diff)) %>%
  print(n = Inf)

# ------------------------------------------------------------------------------
# 9) PLOTTING
#    A) Time series, faceted by env
#    B) Bar chart of final-year spawners
# ------------------------------------------------------------------------------
# Re-factor env codes for prettier labels if desired
results2 <- results %>%
  mutate(env = factor(env, levels = as.character(1:10), labels = paste("Alt",1:10)))

# Determine y-axis max
max_spawn <- max(results2$spawners, na.rm = TRUE)

# 9a) Time series
ggplot(results2, aes(x = year, y = spawners, color = variant)) +
  geom_line(size = 1) +
  facet_wrap(~ env, scales = "fixed") +
  scale_y_continuous(
    limits       = c(0, max_spawn),
    breaks       = pretty_breaks(n = 10),
    minor_breaks = seq(0, max_spawn, length.out = 41),
    labels       = comma
  ) +
  labs(x = "Simulation Year", y = "Escapement (Spawners)", color = "TDM Variant") +
  theme_minimal() +
  theme(legend.position = "bottom")

# 9b) End-of-simulation bar chart
ggplot(
  results2 %>% filter(year == last_year),
  aes(x = variant, y = spawners, fill = env)
) +
  geom_col(position = "dodge") +
  scale_y_continuous(
    limits       = c(0, max_spawn),
    expand       = c(0, 0),
    breaks       = pretty_breaks(n = 30),
    minor_breaks = seq(0, max_spawn, length.out = 30),
    labels       = comma
  ) +
  labs(
    title = paste("End-of-Simulation Spawners (Year", last_year, ")"),
    x = "TDM Variant",
    y = "Spawners",
    fill = "Alternative"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
