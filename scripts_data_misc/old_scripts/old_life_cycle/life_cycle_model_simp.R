# R_model_scalar.R
rm(list = ls())
# 0) packages
library(tidyverse)

# 1) source egg_growth_survival.R into its own environment
egg_env <- new.env()
sys.source("egg_growth_survival.R", envir = egg_env)
egg_summary <- egg_env$summary_tbl

# 2) scenario grid from egg_summary
scenarios <- egg_summary %>%
  filter(egg_m   == "mech",
         hatch_m == "linear",
         dd_type == "none") %>%
  distinct(model, calib) %>%
  arrange(model, calib)

# 3) lookup table of max spawn→fry survival per variant
surv_lookup <- egg_summary %>%
  filter(egg_m   == "mech",
         hatch_m == "linear",
         dd_type == "none") %>%
  group_by(model, calib) %>%
  summarise(
    surv0   = max(mean_cum_surv, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(key = paste(model, calib, sep = "_")) %>%
  select(key, surv0) %>%
  deframe()


# 2) parameters
base_P <- list(
  surv_pre_spawn_hold = 0.85,
  female_fraction     = 0.5,
  fec                  = 5522,
  S0                   = 0.347,
  K_spawners           = 200,
  K_rear               = 5e6,
  S_rear               = 0.5,
  S_overflow           = 0.1,
  surv_smolts_out      = 0.11,
  surv_marine          = 0.5
)
ds_mult <- with(base_P, surv_smolts_out * surv_marine)

# 3) scalar step (DD on redds only)
step_scalar <- function(S_prev, surv0, P, ds) {
  redds    <- S_prev * P$surv_pre_spawn_hold * P$female_fraction
  eggs     <- redds * P$fec                # make sure `P$fec` exists
  dd       <- P$S0 / (1 + redds / P$K_spawners)
  fry_dd   <- eggs * surv0 * dd            # <— dd, not density
  in_cap   <- pmin(fry_dd,      P$K_rear)
  overflow <- pmax(fry_dd - P$K_rear, 0)
  reared   <- in_cap * P$S_rear + overflow * P$S_overflow
  reared * ds
}

# 4) spin‐up constants
years     <- 5
total_escapement <- 45541
S0_year1  <- total_escapement * base_P$surv_pre_spawn_hold

# 5) loop over each variant
results <- map_dfr(names(surv_lookup), function(key) {
  surv0 <- surv_lookup[[key]]
  S     <- numeric(years); S[1] <- S0_year1
  for (t in 2:years) {
    S[t] <- step_scalar(S[t-1], surv0, base_P, ds_mult)
  }
  tibble(year = 1:years, spawners = S, variant = key)
})

# 6) plot
ggplot(results, aes(year, spawners, color=variant)) +
  geom_line(size=1) +
  theme_minimal() +
  labs(
    title="Scalar Model (DD on Redds) Across TDM Variants",
    x="Year", y="Spawners"
  )

# pull out only the last year
last_year   <- max(results$year)
last_spawns <- results %>%
  filter(year == last_year)

ggplot(last_spawns, aes(x = variant, y = spawners, fill = variant)) +
  geom_col() +
  theme_minimal() +
  labs(
    title = paste("End-of-Simulation Spawner Counts (Year", last_year, ")"),
    x     = "TDM Variant",
    y     = "Spawners"
  ) +
  theme(legend.position = "none")
