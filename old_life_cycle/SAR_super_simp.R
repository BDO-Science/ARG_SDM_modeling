#!/usr/bin/env Rscript

# ───────────────────────────────────────────────────────────────────────────────
# scalar_by_env_variant.R
#
# 1. Load per-env, per-variant egg→fry survival from `etf_survival_actual_temps.R`
# 2. Build a surv0 lookup for every (env, model, calib)
# 3. Run your density-dependent scalar model for each env×variant
# 4. Compute %-difference from the best variant *within each env*
# 5. Plot trajectories and end-of-simulation spawners
# ───────────────────────────────────────────────────────────────────────────────

# 0) CLEAN & LIBRARIES  
rm(list=ls())
library(tidyverse)
library(scales)

# 1) LOAD YOUR TDM OUTPUTS  
egg_env     <- new.env()
sys.source("etf_survival_actual_temps.R", envir = egg_env)
egg_summary <- egg_env$results_multi_fast
# egg_summary must contain: env, model, calib, mean_cum_surv, year
egg_summary <- results_obs
# 2) BUILD surv0 LOOKUP (one value per env×model×calib)  
surv_lookup_env <- egg_summary %>%
  group_by(env, variant) %>%
  summarise(
    surv0 = max(mean_cum_surv, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    alt_variant = paste(env, variant, sep = "_")
  )

surv_lookup <- deframe( surv_lookup_env %>% select(alt_variant, surv0) )

# 3) PARAMETERS & SCALAR MODEL  
base_P <- list(
  female_fraction     = 0.5,
  fec                  = 5522,
  S0                   = 0.347,
  K_spawners           = 50000,
  #K_rear               = 5e6,
  S_rear               = 0.15,
  #S_overflow           = 0.1,
  SAR                  = 0.0155,
  lag_probs            = c(`3`=0.75, `4`=0.249, `5`=0.001)
)

simulate_variant <- function(surv_vec, P, years, S_init) {
  if (length(surv_vec) == 1) surv_vec <- rep(surv_vec, years)
  n_init <- length(S_init)
  S      <- numeric(years)
  reared <- numeric(years)
  
  dd_vec     <- numeric(years)  
  fry_dd_vec <- numeric(years)
  
  S[1:n_init] <- S_init
  
  for (t in seq_len(years)) {
    redds  <- S[t] * P$female_fraction
    eggs   <- redds * P$fec
    dd     <- P$S0 / (1 + redds / P$K_spawners)
    fry_dd <- eggs * surv_vec[t] * dd
    # record for output
    dd_vec[t]     <- dd
    fry_dd_vec[t] <- fry_dd
    
    in_cap   <- pmin(fry_dd, P$K_rear)
    overflow <- pmax(fry_dd - P$K_rear, 0)
    reared[t] <- in_cap * P$S_rear + overflow * P$S_overflow
    
    for (age in 3:5) {
      ret_year <- t + age
      if (ret_year <= years) {
        S[ret_year] <- S[ret_year] +
          reared[t] * P$SAR * P$lag_probs[as.character(age)]
      }
    }
  }
  
  tibble(
     year      = seq_len(years),
     spawners  = S,
     dd        = dd_vec,
     fry_dd    = fry_dd_vec,
     eff_surv  = surv_vec * dd_vec     # TDM × density‐dependence
    )
}

# spin-up settings  
years  <- 50
S_seed <- c(16383, 37252, 45541)

# 4) RUN THE MODEL FOR EVERY ENV×VARIANT  
results_full <- map_dfr(names(surv_lookup), function(key) {
  surv0 <- surv_lookup[[key]]
  simulate_variant(surv0, base_P, years, S_seed) %>%
    mutate(variant_full = key)
  })

# 5) SPLIT BACK INTO ENV, VARIANT 
results <- results_full %>%
  # split on the first “_” only
  tidyr::extract(
    variant_full,
    into = c("env","variant"),
    regex = "^([^_]+)_(.+)$"
  ) %>%
  # if you still need an alt_variant for pct‐diff grouping, just duplicate:
  mutate(alt_variant = variant) %>%
  select(env, variant, alt_variant, year, spawners)

# 6) PERCENT-DIFFERENCES WITHIN EACH ENV & YEAR  
pct_results <- results %>%
  group_by(env, year) %>%
  mutate(
    max_spawners = max(spawners),
    pct_diff     = (spawners - max_spawners) / max_spawners * 100
  ) %>%
  ungroup() %>%
  select(env, alt_variant, variant, year, spawners, pct_diff)

# 7) PRINT END-OF-SIMULATION COMPARISONS  
last_year <- max(pct_results$year)
cat("\n% differences from best variant in year", last_year, ":\n")
pct_results %>%
  filter(year == last_year) %>%
  arrange(env, desc(pct_diff)) %>%
  print(n = Inf)

# 8) PLOTS  

# find the overall max so all panels share the same top
max_spawn <- max(results$spawners, na.rm = TRUE)

ggplot(results, aes(year, spawners, color = variant)) +
  geom_line(size = 1) +
  facet_wrap(~ env, scales = "fixed") +    # same y axis on every panel
  scale_y_continuous(
    limits      = c(0, max_spawn),
    breaks      = pretty_breaks(n = 20),          # ≈20 major ticks
    minor_breaks = seq(0, max_spawn, length.out = 41),  # add minor ticks
    labels      = comma
  ) +
  labs(
    x = "Simulation Year",
    y = "Escapement (Spawners)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# 8b) End-of-simulation spawners by env & variant  
ggplot(results %>% filter(year == last_year),
       aes(x = variant, y = spawners, fill = env)) +
  geom_col(position = "dodge") +
  scale_y_continuous(
    limits      = c(0, max_spawn),
    breaks      = pretty_breaks(n = 30),          # ≈20 major ticks
    minor_breaks = seq(0, max_spawn, length.out = 30),  # add minor ticks
    labels      = comma
  ) +
  labs(
    title = paste("End-of-Simulation Spawners (Year", last_year, ")"),
    x = "TDM Variant",
    y = "Spawners"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
