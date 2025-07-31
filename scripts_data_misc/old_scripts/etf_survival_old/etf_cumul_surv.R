# 0) PACKAGES -------------------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(lubridate)
library(future.apply)
library(progressr)
handlers("txtprogressbar")   # simple console bar

# 1) CUMULATIVE TDM & DEVELOPMENT FUNCTIONS --------------------------------------

# -- Cumulative Linear Martin (daily excess T over threshold) -------------------
cum_tdm_lin_martin <- function(temps, α = 0.026, β = 12.14) {
  # temps: vector of daily temperatures
  H <- α * sum(pmax(temps - β, 0))
  exp(-H)
}

# -- Cumulative Linear Anderson (4 days pre-hatch organogenesis) ---------------
cum_tdm_lin_anderson <- function(temps, days_egg, α = 0.026, β = 12.14) {
  # temps: full slice spawn→emergence
  # days_egg: egg-phase length (days)
  egg_temps <- head(temps, ceiling(days_egg))
  organo   <- tail(egg_temps, 4)
  H        <- α * sum(pmax(organo - β, 0))
  exp(-H)
}

# -- Cumulative Exponential (Jensen‐type) ---------------------------------------
cum_tdm_exp <- function(temps, calib = c("WaterForum2020","ZeugCramer2012","SALMOD2006"), stage = c("egg","alevin")) {
  calib <- match.arg(calib)
  stage <- match.arg(stage)
  ps <- list(
    WaterForum2020 = list(egg   = list(α=3.40848e-11, β=1.21122),
                          alevin= list(α=1.017554e-10, β=1.24092)),
    ZeugCramer2012 = list(egg   = list(α=1.35e-8,     β=0.9054),
                          alevin= list(α=1.35e-8,     β=0.9054)),
    SALMOD2006     = list(egg   = list(α=1.475e-11,   β=1.392),
                          alevin= list(α=2.521e-12,   β=1.461))
  )[[calib]][[stage]]
  H <- sum(ps$α * exp(ps$β * temps))
  exp(-H)
}

# -- Cumulative Weibull‐type (Jager 2011) direct cumulative over thermal units ----------
cum_tdm_weibull <- function(temps,
                            TL = 3,     # lower‐temp scaler (°C)
                            TU = 15.4,  # upper‐temp scaler (°C)
                            kL = 25,    # lower exponent
                            kU = 29) {  # upper exponent
  # temps: vector of daily temperatures (°C)
  # T_cum: cumulative thermal units or degree-days
  T_cum <- sum(temps)
  # cumulative survival as per Jager (2011)
  S <- (1 - exp(-T_cum / TL))^kL * (exp(-T_cum / TU))^kU
  return(S)
}

# -- Cumulative no‐T effect & constant survival ---------------------------------
cum_tdm_none  <- function(temps) rep(1, length(temps)) |> prod()
cum_tdm_const <- function(temps, S0 = 0.26) rep(S0, length(temps)) |> prod()

# -- Egg & Alevin development models --------------------------------------------
egg_models <- list(
  linear  = function(T)  958      / T,
  mech    = function(T)  1e3 * 200^0.2 / (pmax(T + 3, 0.1)^1.5),
  jensen  = function(T)  495.28   * T^-0.136,
  beacham = function(T)  1/exp(10.404 - 2.043*log(T + 7.575)),
  usgs    = function(T)  1/exp(6.872  -    log(T + 0.332))
)
hatch_models <- list(
  linear = function(T) 417 / T,
  mech   = function(T) 2.5 * (pmax(T + 2, 0.1)^(-0.5))
)

# -- Density‐dependence ---------------------------------------------------------
density_fun <- function(dd_type, N, S0 = 0.347, mu1 = -1.88e-05, K = 50000) {
  switch(dd_type,
         none           = 1,
         linear         = S0 + mu1 * N,
         `beverton-holt` = S0 / (1 + N / K)
  )
}

# 2) SETUP AND TEMPERATURE/HATCH DATA ------------------------------------------
set.seed(2025)
years      <- 2024:(2024+20-1)
n_redds    <- 1000
peak_date  <- "12-01"
sigma_days <- 20
start_temp <- 19; min_temp <- 8; end_temp <- 14; noise_sd <- 0.5

year_list <- map(years, function(yy) {
  start_date <- ymd(sprintf("%d-10-01", yy))
  feb28      <- ymd(sprintf("%d-02-28", yy+1))
  end_date   <- ymd(sprintf("%d-04-30", yy+1))
  dates      <- seq(start_date, end_date, by = "day")
  # seasonal baseline + noise
  d   <- as.numeric(dates - start_date)
  n1  <- as.numeric(feb28 - start_date)
  n2  <- as.numeric(end_date - feb28)
  base_temp <- ifelse(
    dates <= feb28,
    start_temp + (min_temp - start_temp) * (d / n1),
    min_temp   + (end_temp   - min_temp)   * ((d - n1) / n2)
  )
  temps <- base_temp + rnorm(length(dates), sd = noise_sd)
  # spawn dates
  peak_idx <- as.numeric(ymd(sprintf("%d-%s", yy, peak_date)) - start_date)
  offs     <- round(rnorm(n_redds, mean = peak_idx, sd = sigma_days))
  offs     <- pmin(pmax(offs, 0), length(dates) - 1)
  redds    <- start_date + offs
  list(year = yy, dates = dates, temps = temps, redds = redds)
})
names(year_list) <- as.character(years)

# TDM definitions for variants
library(tibble)
tdm_defs <- tribble(
  ~model,        ~calib,           ~stage,   ~variant,
  "exp",        "WaterForum2020","egg",    "exp_WF_egg",
  "exp",        "WaterForum2020","alevin", "exp_WF_alevin",
  "exp",        "ZeugCramer2012","egg",    "exp_ZC_egg",
  "exp",        "ZeugCramer2012","alevin", "exp_ZC_alevin",
  "exp",        "SALMOD2006",    "egg",    "exp_SM_egg",
  "exp",        "SALMOD2006",    "alevin", "exp_SM_alevin",
  "lin_martin", NA,               NA,        "lin_martin",
  "lin_anderson",NA,              NA,        "lin_anderson",
  "weibull",    NA,               NA,        "weibull",
  "none",       NA,               NA,        "none",
  "const",      NA,               NA,        "const"
)

###########################
#MULTI-YEAR VERSION
#########################
# 3) RUN SIMULATION & AGGREGATE ------------------------------------------------
plan(multisession, workers = 4)
handlers("txtprogressbar")
all_years_summary <- with_progress({
  p <- progressor(along = year_list)
  map_dfr(year_list, function(yr) {
    p()
    dates    <- yr$dates
    t_series <- yr$temps
    redds    <- yr$redds
    
    # cross‐variations
    tidyr::crossing(
      tdm_defs,
      egg_m   = names(egg_models),
      hatch_m = names(hatch_models),
      dd_type = c("none","linear","beverton-holt")
    ) %>% rowwise() %>% mutate(
      survs = list(vapply(
        redds, function(sp) {
          i0    <- match(sp, dates)
          days_egg  <- egg_models[[egg_m]](t_series[i0])
          days_alev <- hatch_models[[hatch_m]](t_series[i0])
          td         <- ceiling(days_egg + days_alev)
          temps_slice <- t_series[i0 + seq_len(td) - 1]
          
          # cumulative TDM for this variant
          base_surv <- switch(model,
                              exp          = cum_tdm_exp(temps_slice, calib, stage),
                              lin_martin   = cum_tdm_lin_martin(temps_slice),
                              lin_anderson = cum_tdm_lin_anderson(temps_slice, days_egg),
                              weibull      = cum_tdm_weibull(temps_slice),
                              none         = cum_tdm_none(temps_slice),
                              const        = cum_tdm_const(temps_slice)
          )
          base_surv
        }, numeric(1)
      )),
      base_surv     = mean(unlist(survs), na.rm = TRUE),
      dd_factor     = density_fun(dd_type, N = length(redds)),
      mean_cum_surv = base_surv * dd_factor,
      # development summary
      dev0 = {
        Te  <- t_series[match(redds, dates)]
        de  <- egg_models[[egg_m]](Te) + hatch_models[[hatch_m]](Te)
        atu <- de * Te
        tibble(mean_dev_days = mean(de), mean_ATU = mean(atu))
      }
    ) %>% unnest_wider(dev0) %>% ungroup() %>%
      mutate(across(c(mean_dev_days, mean_ATU, mean_cum_surv), ~ round(.x,3)), year = yr$year) %>%
      select(year, everything())
  })
})

# 4) PLOT RESULTS ---------------------------------------------------------------
ggplot(all_years_summary, aes(x = year, y = mean_cum_surv, color = variant)) +
  geom_line(size = 1) +
  facet_grid(egg_m ~ dd_type) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(x = "Water Year", y = "Mean Egg→Emergence Survival", color = "Variant") +
  theme_minimal() + theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))


#######################
#SINGLE YEAR VERSION
#########################
single_year <- year_list["2024"]

one_year_summary <- map_dfr(single_year, function(yr) {
    p()
    dates    <- yr$dates
    t_series <- yr$temps
    redds    <- yr$redds
    
    # cross‐variations
    tidyr::crossing(
      tdm_defs,
      egg_m   = names(egg_models),
      hatch_m = names(hatch_models),
      dd_type = c("none","linear","beverton-holt")
    ) %>% rowwise() %>% mutate(
      survs = list(vapply(
        redds, function(sp) {
          i0    <- match(sp, dates)
          days_egg  <- egg_models[[egg_m]](t_series[i0])
          days_alev <- hatch_models[[hatch_m]](t_series[i0])
          td         <- ceiling(days_egg + days_alev)
          temps_slice <- t_series[i0 + seq_len(td) - 1]
          
          # cumulative TDM for this variant
          base_surv <- switch(model,
                              exp          = cum_tdm_exp(temps_slice, calib, stage),
                              lin_martin   = cum_tdm_lin_martin(temps_slice),
                              lin_anderson = cum_tdm_lin_anderson(temps_slice, days_egg),
                              weibull      = cum_tdm_weibull(temps_slice),
                              none         = cum_tdm_none(temps_slice),
                              const        = cum_tdm_const(temps_slice)
          )
          base_surv
        }, numeric(1)
      )),
      base_surv     = mean(unlist(survs), na.rm = TRUE),
      dd_factor     = density_fun(dd_type, N = length(redds)),
      mean_cum_surv = base_surv * dd_factor,
      # development summary
      dev0 = {
        Te  <- t_series[match(redds, dates)]
        de  <- egg_models[[egg_m]](Te) + hatch_models[[hatch_m]](Te)
        atu <- de * Te
        tibble(mean_dev_days = mean(de), mean_ATU = mean(atu))
      }
    ) %>% unnest_wider(dev0) %>% ungroup() %>%
      mutate(across(c(mean_dev_days, mean_ATU, mean_cum_surv), ~ round(.x,3)), year = yr$year) %>%
      select(year, everything())
  })

# 4) PLOT RESULTS ---------------------------------------------------------------
ggplot(one_year_summary, aes(x = year, y = mean_cum_surv, color = variant)) +
  geom_line(size = 1) +
  facet_grid(egg_m ~ dd_type) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(x = "Water Year", y = "Mean Egg→Emergence Survival", color = "Variant") +
  theme_minimal() + theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

