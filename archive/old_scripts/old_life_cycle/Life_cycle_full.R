# ───────────────────────────────────────────────────────────────────────────────
# 0) PACKAGES & SETUP -----------------------------------------------------------
# ───────────────────────────────────────────────────────────────────────────────
library(tidyverse)
library(lubridate)

# 0.1) simulate daily water temperatures & spawn‐window
dates        <- seq(as.Date("2024-10-15"), as.Date("2025-02-13"), by="day")
day0         <- as.numeric(dates - dates[1])
temps        <- 15 - 0.03 * day0 + rnorm(length(day0), sd = 0.5)
spawn_window <- seq(as.Date("2024-10-15"), as.Date("2024-11-29"), by="day")

# ───────────────────────────────────────────────────────────────────────────────
# 1) TDM FUNCTIONS ---------------------------------------------------------------
# ───────────────────────────────────────────────────────────────────────────────
tdm_lin_martin <- function(T, α=0.026, β=12.14) exp(-α * pmax(T - β, 0))
tdm_lin_anderson <- function(T, day, α=0.026, β=12.14) {
  exp(-ifelse(day <= 4, α * pmax(T - β, 0), 0))
}
tdm_exp <- function(T, calib, stage) {
  p <- list(
    WaterForum2020 = list(egg=list(α=3.40848e-11,β=1.21122),
                          alevin=list(α=1.017554e-10,β=1.24092)),
    ZeugCramer2012 = list(egg=list(α=1.35e-8,β=0.9054),
                          alevin=list(α=1.35e-8,β=0.9054)),
    SALMOD2006     = list(egg=list(α=1.475e-11,β=1.392),
                          alevin=list(α=2.521e-12,β=1.461))
  )[[calib]][[stage]]
  exp(- (p$α * exp(p$β * T)))
}
tdm_weibull <- function(T, kL=25, kU=29) {
  M <- 1 - (1 - exp(-T/kL))^(T/kU); exp(-M)
}
tdm_none  <- function(T) rep(1, length(T))
tdm_const <- function(T, S0=0.26) rep(S0, length(T))

# ───────────────────────────────────────────────────────────────────────────────
# 2) EGG & ALEVIN DEVELOPMENT ---------------------------------------------------
# ───────────────────────────────────────────────────────────────────────────────
egg_models <- list(
  linear  = function(T)  958 / T,
  mech    = function(T)  1e3*200^0.2/(pmax(T+3,0.1)^1.5),
  jensen  = function(T)  495.28 * T^-0.136,
  beacham = function(T)  1/exp(10.404 - 2.043*log(T+7.575)),
  usgs    = function(T)  1/exp(6.872  -    log(T+0.332))
)
hatch_models <- list(
  linear = function(T) 417 / T,
  mech   = function(T) 2.5 * (pmax(T+2,0.1)^(-0.5))
)

# ───────────────────────────────────────────────────────────────────────────────
# 3) DENSITY‐DEPENDENCE ----------------------------------------------------------
# ───────────────────────────────────────────────────────────────────────────────
density_fun <- function(dd, N, S0=0.347, mu1=-0.0000188, K=10000) {
  switch(dd,
         none            = 1,
         linear          = S0 + mu1 * N,
         `beverton-holt` = S0 / (1 + N/K)
  )
}

# ───────────────────────────────────────────────────────────────────────────────
# 4) PER‐REDD SURVIVAL SIMULATOR ------------------------------------------------
# ───────────────────────────────────────────────────────────────────────────────
simulate_combo <- function(model, calib,
                           egg_fn, hatch_fn,
                           spawn_dates, dates, temps,
                           dd_type) {
  # returns survival per redd
  sapply(spawn_dates, function(sp) {
    i0         <- match(sp, dates)
    days_egg   <- egg_fn(   temps[i0] )
    days_alev  <- hatch_fn( temps[i0] )
    total_days <- ceiling(days_egg + days_alev)
    idx        <- pmin(i0 + seq_len(total_days) - 1, length(dates))
    survs      <- sapply(seq_along(idx), function(k) {
      Tday    <- temps[idx[k]]
      S_tdm   <- switch(model,
                        lin_martin   = tdm_lin_martin(Tday),
                        lin_anderson = tdm_lin_anderson(Tday, k),
                        exp          = {
                          stage <- if(k <= ceiling(days_egg)) "egg" else "alevin"
                          tdm_exp(Tday, calib, stage)
                        },
                        weibull  = tdm_weibull(Tday),
                        none     = tdm_none(Tday),
                        const    = tdm_const(Tday)
      )
      N_active <- length(spawn_dates)
      S_tdm * density_fun(dd_type, N_active)
    })
    prod(survs)
  })
}

# ───────────────────────────────────────────────────────────────────────────────
# 5) WRAPPER FOR EGG→FRY SURVIVAL ------------------------------------------------
# ───────────────────────────────────────────────────────────────────────────────
compute_egg_to_fry_surv <- function(n_redds,
                                    model, calib,
                                    egg_m, hatch_m,
                                    dd_type,
                                    dates, temps, spawn_window) {
  # — convert to a clean integer number of redds —
  n_redds_int <- as.integer(floor(n_redds))
  
  # — if it's NA or zero, just return 100% survival (eggs*1 = 0 fry when eggs=0) —
  if (is.na(n_redds_int) || n_redds_int <= 0) {
    return(1)
  }
  
  # — now safe to sample that many spawn dates —
  spawn_dates <- sample(spawn_window, n_redds_int, replace = TRUE)
  
  # — the rest is untouched —
  surv_vec <- simulate_combo(
    model       = model,
    calib       = calib,
    egg_fn      = egg_models[[egg_m]],
    hatch_fn    = hatch_models[[hatch_m]],
    spawn_dates = spawn_dates,
    dates       = dates,
    temps       = temps,
    dd_type     = dd_type
  )
  mean(surv_vec, na.rm = TRUE)
}

# ───────────────────────────────────────────────────────────────────────────────
# 6) LIFE‐CYCLE FUNCTIONS -------------------------------------------------------
# ───────────────────────────────────────────────────────────────────────────────
step_cohort <- function(cohort, P) {
  # adults → spawners
  total_spawners <- sum(cohort * P$spawn_at_age)
  spawners       <- total_spawners * P$surv_pre_spawn_hold
  
  # redds & eggs
  redds <- spawners * P$female_fraction
  eggs  <- redds * P$fecundity_per_female
  
  # dynamic egg→fry survival
  S_egg   <- compute_egg_to_fry_surv(
    n_redds      = redds,
    model        = P$tdm_model,
    calib        = P$tdm_calib,
    egg_m        = P$egg_model,
    hatch_m      = P$hatch_model,
    dd_type      = P$dd_type,
    dates        = P$dates,
    temps        = P$temps,
    spawn_window = P$spawn_window
  )
  fry_raw <- eggs * S_egg
  
  # fry rearing & smolt pipeline (unchanged)...
  fry_stream_prop <- fry_raw * P$prop_stream_fry
  fry_fp_prop     <- fry_raw * P$prop_fp_fry
  
  fry_stream   <- min(fry_stream_prop, P$cap_stream_fry)
  spill_stream <- max(fry_stream_prop - P$cap_stream_fry, 0)
  
  fry_fp       <- min(fry_fp_prop + spill_stream, P$cap_fp_fry)
  spill_fp     <- max((fry_fp_prop + spill_stream) - P$cap_fp_fry, 0)
  fry_overflow <- spill_fp
  
  survivors_stream   <- fry_stream   * P$surv_stream_rear
  survivors_fp       <- fry_fp       * P$surv_fp_surv
  survivors_overflow <- fry_overflow * P$overflow_fry_surv
  reared_fry         <- survivors_stream + survivors_fp + survivors_overflow
  
  smolts_raw      <- reared_fry
  smolts_instream <- min(smolts_raw, P$cap_smolts)
  smolts_overflow <- max(smolts_raw - P$cap_smolts, 0)
  
  survivors_smin <- smolts_instream * P$surv_smolts_instream
  survivors_smov <- smolts_overflow * P$surv_smolt_overflow
  smolts_total   <- survivors_smin + survivors_smov
  
  delta_instream <- smolts_total * P$surv_delta_instream
  recruits       <- delta_instream * P$surv_ocean_annual
  
  # age‐up
  new_cohort   <- numeric(length(P$spawn_at_age))
  new_cohort[1] <- recruits
  for(a in 2:length(new_cohort)) {
    new_cohort[a] <- cohort[a-1] * P$surv_ocean_annual
  }
  # plus‐group
  new_cohort[length(new_cohort)] <-
    new_cohort[length(new_cohort)] +
    cohort[length(cohort)] * P$surv_ocean_annual
  
  new_cohort
}

run_years <- function(years, init_cohort, P) {
  mat <- matrix(0, nrow=years, ncol=length(P$spawn_at_age))
  mat[1,] <- init_cohort
  for(t in 1:(years-1)) {
    mat[t+1,] <- step_cohort(mat[t,], P)
  }
  as_tibble(mat) %>%
    mutate(year=1:years) %>%
    pivot_longer(-year, names_to="age", values_to="N")
}

# ───────────────────────────────────────────────────────────────────────────────
# 7) SPIN‐UP & SCENARIOS --------------------------------------------------------
# ───────────────────────────────────────────────────────────────────────────────
# base parameters
base_P <- list(
  fecundity_per_female = 2000,
  female_fraction      = 0.5,
  surv_pre_spawn_hold  = 0.85,
  # placeholder for egg→fry
  dates        = dates,
  temps        = temps,
  spawn_window = spawn_window,
  # default TDM settings (will be overridden per scenario)
  tdm_model    = "exp",
  tdm_calib    = "WaterForum2020",
  egg_model    = "linear",
  hatch_model  = "linear",
  dd_type      = "beverton-holt",
  # fry rearing
  prop_stream_fry      = 0.8,
  prop_fp_fry          = 0.2,
  cap_stream_fry       = 5e5,
  cap_fp_fry           = 1e5,
  surv_stream_rear     = 0.3,
  surv_fp_surv         = 0.4,
  overflow_fry_surv    = 0.1,
  # smolts
  cap_smolts           = 2e5,
  surv_smolts_instream = 0.3,
  surv_smolt_overflow  = 0.1,
  surv_delta_instream  = 0.1,
  surv_ocean_annual    = 0.8,
  # adult spawning
  spawn_at_age         = c(0,0.2,0.9,1.0,1.0)  # ages 3–5 spawn
)

# spin‐up escapement‐only to get init_cohort
esc_total <- 15000
esc      <- c(0,0, esc_total*c(0.7,0.2,0.1))
co1      <- esc
co2      <- step_cohort(co1, modifyList(base_P, list(tdm_model="none", tdm_calib=NA)))
co3      <- step_cohort(co2, modifyList(base_P, list(tdm_model="none", tdm_calib=NA)))
init_cohort <- c(co2[1], co3[2], esc[3:5])

# define multiple TDM scenarios
scenarios <- tribble(
  ~tdm_model,    ~tdm_calib,
  "lin_martin",  NA,
  "lin_anderson",NA,
  "weibull",     NA,
  "none",        NA,
  "const",       NA,
  "exp",      "WaterForum2020",
  "exp",      "ZeugCramer2012",
  "exp",      "SALMOD2006"
)

# run each for 50 years, track age‐5 spawners
results <- scenarios %>%
  mutate(sim = pmap(
    list(tdm_model, tdm_calib),
    function(tm, tc) {
      P <- modifyList(base_P, list(tdm_model=tm, tdm_calib=tc))
      run_years(50, init_cohort, P) %>%
        filter(age == paste0("V",5)) %>%           # age‐5 survivors
        transmute(year, spawners=N, variant=paste(tm,tc,sep="_"))
    }
  )) %>%
  unnest(sim)

# ───────────────────────────────────────────────────────────────────────────────
# 8) PLOT RESULTS ---------------------------------------------------------------
# ───────────────────────────────────────────────────────────────────────────────
ggplot(results, aes(year, spawners, color=variant)) +
  geom_line(size=1) +
  labs(
    title = "Spawning Cohort Dynamics under Different TDM Scenarios",
    x     = "Year",
    y     = "Number of Age-5 Spawners",
    color = "TDM Scenario"
  ) +
  theme_minimal()
