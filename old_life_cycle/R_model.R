# R_model.R
# 1) pull in egg_growth_survival.R and capture its return value
egg_summary <- source("egg_growth_survival.R", local = TRUE)$value

# 2) pick the bit you actually want
surv_base <- egg_summary %>%
  filter(model=="exp",
         calib=="WaterForum2020",
         egg_m=="mech",
         hatch_m=="linear",
         dd_type=="none") %>%
  pull(mean_cum_surv)

# 3) build your parameter list
base_P <- list(
  # … everything else …,
  surv_egg_base = surv_base,
  dd_type       = "beverton-holt"
)

# ──────────────────────────────────────────────────────────────────────────────
# 0) Dependencies
# ──────────────────────────────────────────────────────────────────────────────
library(tidyverse)

# ──────────────────────────────────────────────────────────────────────────────
# 1) One‐year step and multi‐year runner (exactly as you had)
# ──────────────────────────────────────────────────────────────────────────────
step_cohort <- function(cohort, P) {
  # 1) spawners & holding
  total_spawners <- sum(cohort * P$spawn_at_age)
  spawners       <- total_spawners * P$surv_pre_spawn_hold
  
  # 2) female spawners = redd count
  redds <- spawners * P$female_fraction
  
  # 3) eggs produced
  eggs  <- redds * P$fecundity_per_female
  
  # 4) temperature‐driven survival (from SacPAS or fixed TDM)
  S_temp <- P$surv_egg_to_fry
  fry_T  <- eggs * S_temp
  
  # 5) density‐dependent survival (Essington et al. 2000)
  #    S0 / (1 + A/K)  where A = redd count, K = capacity in same units
  S_density <- P$S0 / (1 + redds / P$K_spawners)
  fry_D     <- fry_T * S_density
  
  # now fry_D is your egg→fry cohort, already down‐weighted
  # by both temperature AND density effects
  
  # 6) allocate fry_D into rearing bins and proceed exactly as before…
  fry_raw <- fry_D
  
  # proportional + capacity‐constrained allocation
  fry_stream_prop <- fry_raw * P$prop_stream_fry
  fry_fp_prop     <- fry_raw * P$prop_fp_fry
  
  fry_stream   <- min(fry_stream_prop, P$cap_stream_fry)
  spill_stream <- max(fry_stream_prop - P$cap_stream_fry, 0)
  
  fry_fp       <- min(fry_fp_prop + spill_stream, P$cap_floodplain_fry)
  spill_fp     <- max((fry_fp_prop + spill_stream) - P$cap_floodplain_fry, 0)
  
  fry_overflow <- spill_fp
  
  # bin‐specific survival
  survivors_stream   <- fry_stream   * P$surv_stream_rear
  survivors_fp       <- fry_fp       * P$surv_fp_surv
  survivors_fry_over <- fry_overflow * P$overflow_fry_surv
  
  reared_fry <- survivors_stream + survivors_fp + survivors_fry_over
  
  # 7) smolt capacity, delta, ocean, aging… (same as before)
  smolts_raw      <- reared_fry
  smolts_instream <- min(smolts_raw, P$cap_smolts)
  smolts_overflow <- max(smolts_raw - P$cap_smolts, 0)
  
  survivors_smolts_in   <- smolts_instream   * P$surv_smolts_instream
  survivors_smolts_over <- smolts_overflow   * P$surv_smolt_overflow
  smolts_total_surviving<- survivors_smolts_in + survivors_smolts_over
  
  delta_instream <- smolts_total_surviving * P$surv_delta_instream
  # if you want split‐delta, do the same two‐pool logic here…
  
  recruits_age1 <- delta_instream * P$surv_ocean_annual
  
  # 8) age up
  new_cohort <- numeric(length(P$spawn_at_age))
  new_cohort[1] <- recruits_age1
  for (a in 2:length(new_cohort)) {
    new_cohort[a] <- cohort[a-1] * P$surv_ocean_annual
  }
  new_cohort[length(new_cohort)] <- new_cohort[length(new_cohort)] +
    cohort[length(new_cohort)] * P$surv_ocean_annual
  
  return(new_cohort)
}

run_years <- function(years, init_cohort, P) {
  mat <- matrix(0, nrow = years, ncol = length(P$spawn_at_age))
  mat[1,] <- init_cohort
  for (t in seq_len(years-1)) {
    mat[t+1,] <- step_cohort(mat[t,], P)
  }
  as_tibble(mat) %>%
    mutate(year = 1:years) %>%
    pivot_longer(-year, names_to = "age", values_to = "N")
}

# ──────────────────────────────────────────────────────────────────────────────
# 4) Wrapper to run entire sim at a given TDM
# ──────────────────────────────────────────────────────────────────────────────
run_with_tdm <- function(tdm, years, init_cohort, base_P){
  P <- base_P
  P$surv_egg_to_fry <- tdm
  df <- run_years(years, init_cohort, P)
  df %>% 
    mutate(age_num = as.integer(parse_number(age)),
           spawners = if_else(age_num >= 3,
                              N * P$spawn_at_age[age_num], 0)) %>%
    group_by(year) %>%
    summarise(spawners = sum(spawners), .groups="drop") %>%
    mutate(tdm = tdm)
}

# ──────────────────────────────────────────────────────────────────────────────
# 2) Base parameter template (tdm left blank)
# ──────────────────────────────────────────────────────────────────────────────
base_P <- list(
  # spawning
  fecundity_per_female   = 2000,
  female_fraction        = 0.5,
  surv_pre_spawn_hold    = 0.85,
  
  # egg→fry (we’ll overwrite this each run)
  surv_egg_to_fry        = NA,
  # Essington BH density‐dependence (S₀/(1 + A/K))
  S0                = 0.347,    # S₀: baseline egg→fry survival
  K_spawners        = 10000,  # K: capacity in # of redds
  
  # fry‐stage: proportions + capacities + bin‐specific survival
  prop_stream_fry        = 0.8,     # 60% of fry choose main‐stem first
  prop_fp_fry            = 0.2,     # 40% choose floodplain first
  cap_stream_fry         = 5e5,     # max stream fry
  cap_floodplain_fry     = 1e5,     # max floodplain fry
  surv_stream_rear       = 0.3,     # stream rearing survival
  surv_fp_surv           = 0.4,     # floodplain rearing survival
  overflow_fry_surv      = 0.1,     # survival for overflow fry
  
  # smolt stage: capacity + overflow survival
  cap_smolts             = 2e5,
  surv_smolts_instream   = 0.3,
  surv_smolt_overflow    = 0.1,
  
  # delta and ocean
  surv_delta_instream    = 0.1,
  surv_delta_overflow    = 0.02,
  surv_ocean_annual      = 0.8,
  
  # adult return
  prob_enter_spawn_river = 0.9,
  surv_adult_immigration = 0.9,
  
  # age‐at‐return (only ages 3–5 spawn)
  spawn_at_age           = c(0, 0.2, 0.9, 1.0, 1.0)
)

# ──────────────────────────────────────────────────────────────────────────────
# 3) Escapement‐only spin‐up (as you had it)
# ──────────────────────────────────────────────────────────────────────────────
total_escapement <- 15000
age3_frac <- 0.7; age4_frac <- 0.2; age5_frac <- 0.1
esc <- c(0,0,
         total_escapement*age3_frac,
         total_escapement*age4_frac,
         total_escapement*age5_frac)

# two‐year spin‐up
cohort_y1 <- esc
cohort_y2 <- step_cohort(cohort_y1, modifyList(base_P, list(surv_egg_to_fry=0.5)))
cohort_y3 <- step_cohort(cohort_y2, modifyList(base_P, list(surv_egg_to_fry=0.5)))

init_cohort <- c(cohort_y2[1], cohort_y3[2], esc[3:5])


# ──────────────────────────────────────────────────────────────────────────────
# 5) Run for whatever tdms your heart desires
# ──────────────────────────────────────────────────────────────────────────────

###if multiple tdm values
tdm_vals <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
all_res  <- map_dfr(tdm_vals, run_with_tdm, 
                    years=50, init_cohort=init_cohort, base_P=base_P)

ggplot(all_res, aes(year, spawners, color = factor(tdm))) +
  geom_line(size = 1) +
  coord_cartesian(ylim = c(0, 30000)) +
  labs(
    x     = "Year",
    y     = "Spawners",
    color = "TDM"
  ) +
  theme_minimal()
  


###############purrrrring#####################
# keep everything else the same, then:
res <- map_dfr(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), 
               ~ run_with_tdm(.x, 
                              years=100, init_cohort=init_cohort, base_P=base_P),
               .id="scenario")

ggplot(res, aes(year, spawners, color = factor(tdm))) +
  geom_line(size = 1) +
  scale_y_continuous(
    limits = c(0, 30000),  # set y-axis from 0 to 30k
    expand = c(0, 0)       # drop any padding at the ends
  ) +
  labs(
    x     = "Year",
    y     = "Spawners",
    color = "TDM"
  ) +
  theme_minimal()













