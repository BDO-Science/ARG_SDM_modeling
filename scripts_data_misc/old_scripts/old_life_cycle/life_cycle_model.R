# R_model.R

# 0) packages
library(tidyverse)

# 1) source egg_growth_survival.R into its own environment
egg_env <- new.env()
sys.source("egg_growth_survival.R", envir = egg_env)

# now egg_env$summary_tbl is your 330×8 tibble
egg_summary <- egg_env$summary_tbl

# 1) build your scenario grid directly from egg_summary
scenarios <- egg_summary %>%
  filter(
    egg_m   == "mech",
    hatch_m == "linear",
    dd_type == "none"
  ) %>%
  distinct(model, calib) %>%
  arrange(model, calib)

# inspect:
print(scenarios)

# 2) pull out all the base survivals into a named lookup
surv_lookup <- egg_summary %>%
  filter(
    egg_m   == "mech",
    hatch_m == "linear",
    dd_type == "none"
  ) %>%
  select(model, calib, mean_cum_surv) %>%
  # paste model & calib into a single key
  mutate(key = paste(model, calib, sep = "_")) %>%
  select(key, mean_cum_surv) %>%
  deframe()

# ───────────────────────────────────────────────────────────────────────────────
# 2) BUILD PARAMETER LIST -------------------------------------------------------
base_P <- list(
  # life‐history
  spawn_at_age         = c(0, 0, 0.7, 0.2, 0.1),  # ages 1–5; only 3:5 spawn
  female_fraction      = 0.5,                    # fraction of spawners that are female
  fecundity_per_female = 5522,
  surv_pre_spawn      = 0.85,                   # pre‐spawn holding survival
  
  # egg→fry
  surv_egg_to_fry      = surv_base,              # from egg_summary  
  S0                   = 0.347,                  # base density survival intercept
  K_spawners           = 50000,                    # capacity in redd‐units
  
  # fry rearing
  prop_stream_fry      = 0.8,
  prop_fp_fry          = 0.2,
  cap_stream_fry       = 5e5,
  cap_floodplain_fry   = 1e5,
  surv_stream_rear     = 0.3,
  surv_fp_surv         = 0.4,
  overflow_fry_surv    = 0.1,
  
  # smolts
  cap_smolts           = 2e5,
  surv_smolts_instream = 0.3,
  surv_smolt_overflow  = 0.1,
  
  # delta‐inflow survival to ocean
  surv_delta_instream  = 0.1,
  
  # ocean survival per year
  surv_ocean_annual    = 0.8
)

# ───────────────────────────────────────────────────────────────────────────────
# 3) STEP‐COHORT: ONE YEAR OF the LC LOOP ---------------------------------------
step_cohort <- function(cohort, P) {
  # 1) adult spawners & holding
  total_spawners <- sum(cohort * P$spawn_at_age)
  spawners       <- total_spawners * P$surv_pre_spawn
  
  # 2) redd count
  redds <- spawners * P$female_fraction
  
  # 3) eggs produced
  eggs  <- redds * P$fecundity_per_female
  
  # 4) temperature‐driven survival (constant here)
  S_temp <- P$surv_egg_to_fry
  fry_T  <- eggs * S_temp
  
  # 5) density‐dependence: S0 / (1 + A / K)
  S_density <- P$S0 / (1 + redds / P$K_spawners)
  fry_D     <- fry_T * S_density
  
  # 6) fry rearing + capacity constraints
  fry_raw         <- fry_D
  fry_stream_prop <- fry_raw * P$prop_stream_fry
  fry_fp_prop     <- fry_raw * P$prop_fp_fry
  
  fry_stream   <- min(fry_stream_prop, P$cap_stream_fry)
  spill_stream <- max(fry_stream_prop - P$cap_stream_fry, 0)
  
  fry_fp       <- min(fry_fp_prop + spill_stream, P$cap_floodplain_fry)
  spill_fp     <- max((fry_fp_prop + spill_stream) - P$cap_floodplain_fry, 0)
  fry_overflow <- spill_fp
  
  survivors_stream   <- fry_stream * P$surv_stream_rear
  survivors_fp       <- fry_fp   * P$surv_fp_surv
  survivors_overflow <- fry_overflow * P$overflow_fry_surv
  reared_fry         <- survivors_stream + survivors_fp + survivors_overflow
  
  # 7) smolt pipeline
  smolts_raw        <- reared_fry
  smolts_instream   <- min(smolts_raw, P$cap_smolts)
  smolts_overflow   <- max(smolts_raw - P$cap_smolts, 0)
  
  surv_sm_instream  <- smolts_instream * P$surv_smolts_instream
  surv_sm_overflow  <- smolts_overflow * P$surv_smolt_overflow
  smolts_total      <- surv_sm_instream + surv_sm_overflow
  
  delta_instream    <- smolts_total * P$surv_delta_instream
  recruits_age1     <- delta_instream * P$surv_ocean_annual
  
  # 8) age‐up into next year’s cohort
  new_cohort <- numeric(length(P$spawn_at_age))
  new_cohort[1] <- recruits_age1
  for (a in 2:length(new_cohort)) {
    new_cohort[a] <- cohort[a-1] * P$surv_ocean_annual
  }
  # plus‐group
  new_cohort[length(new_cohort)] <-
    new_cohort[length(new_cohort)] +
    cohort[length(new_cohort)] * P$surv_ocean_annual
  
  return(new_cohort)
}

# ──────────────────────────────────────────────────────────────────────────────
# 3) Escapement‐only spin‐up (as you had it)
# ──────────────────────────────────────────────────────────────────────────────
total_escapement <- 45541
age3_frac <- 0.7; age4_frac <- 0.2; age5_frac <- 0.1
esc <- c(0,0,
         total_escapement*age3_frac,
         total_escapement*age4_frac,
         total_escapement*age5_frac)

# ───────────────────────────────────────────────────────────────────────────────
# 4) RUN‐YEARS & SCENARIOS  (same as before) ------------------------------------
run_years <- function(years, init_cohort, P) {
  mat <- matrix(0, nrow = years, ncol = length(P$spawn_at_age))
  mat[1,] <- init_cohort
  for (t in seq_len(years-1)) {
    mat[t+1,] <- step_cohort(mat[t,], P)
  }
  as_tibble(mat) %>%
    mutate(year = 1:years) %>%
    pivot_longer(-year, names_to="age", values_to="N")
}

# define init_cohort exactly as you did before…
# define scenarios: a tribble of tdm_model/tdm_calib if needed…

# ───────────────────────────────────────────────────────────────────────────────
# 5) SIMULATE & PLOT ------------------------------------------------------------

# 3) run your loop over that dynamic grid
results <- scenarios %>%
  mutate(sim = pmap(
    list(model, calib),
    function(tm, tc) {
      key   <- paste(tm, tc, sep = "_")
      surv0 <- surv_lookup[[key]]
      
      # two param‐sets: spin‐up (no‐dd) vs run (with‐dd)
      Pspin <- modifyList(base_P, list(
        surv_egg_to_fry = surv0,
        S0               = 1
      ))
      Prun  <- modifyList(base_P, list(
        surv_egg_to_fry = surv0,
        S0               = base_P$S0
      ))
      
      # escapement spin‐up (as before)
      cohort_y1 <- esc
      cohort_y2 <- step_cohort(cohort_y1, Pspin)
      cohort_y3 <- step_cohort(cohort_y2, Pspin)
      init_cohort <- c(cohort_y2[1], cohort_y3[2], esc[3:5])
      
      # full 50-yr run
      run_years(50, init_cohort, Prun) %>%
        filter(age %in% paste0("V", 3:5)) %>%
        group_by(year) %>%
        summarise(
          spawners = sum(
            N *
              Prun$spawn_at_age[as.integer(sub("V","", age))]
          ),
          variant = key,
          .groups = "drop"
        )
    }
  )) %>%
  unnest(sim)

# 4) plot
ggplot(results, aes(year, spawners, color = variant)) +
  geom_line() +
  theme_minimal() +
  labs(
    title = "Spawner Trajectories for Every TDM Variant",
    x     = "Year",
    y     = "Spawners"
  )
