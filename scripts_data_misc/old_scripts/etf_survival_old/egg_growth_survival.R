rm(list=ls())
# ───────────────────────────────────────────────────────────────────────────────
# 0) PACKAGES -------------------------------------------------------------------
# ───────────────────────────────────────────────────────────────────────────────
library(tidyverse)
library(lubridate)
library(future.apply)
library(progressr)
handlers("txtprogressbar")   # simple console bar

# ───────────────────────────────────────────────────────────────────────────────
# 1) TDM & DEV FUNCS -----------------------------------------------------------
# ───────────────────────────────────────────────────────────────────────────────

# -- TDM variants --
tdm_lin_martin <- function(temp, α=0.026, β=12.14)      exp(-α*pmax(temp-β,0))
tdm_lin_anderson <- function(temp, day, α=0.026, β=12.14) {
  exp(-ifelse(day<=4, α*pmax(temp-β,0), 0))
}
tdm_exp <- function(temp, calib, stage) {
  ps <- list(
    WaterForum2020 = list(egg   = list(α=3.40848e-11,β=1.21122),
                          alevin= list(α=1.017554e-10,β=1.24092)),
    ZeugCramer2012 = list(egg   = list(α=1.35e-8,    β=0.9054),
                          alevin= list(α=1.35e-8,    β=0.9054)),
    SALMOD2006     = list(egg   = list(α=1.475e-11, β=1.392),
                          alevin= list(α=2.521e-12, β=1.461))
  )[[calib]][[stage]]
  exp(- (ps$α * exp(ps$β * temp)) )
}
# Weibull‐type (Jager 2011) for the egg phase defaults:
tdm_weibull <- function(temp,
                        TL = 3,     # lower‐temp scaler (°C)
                        TU = 15.4,  # upper‐temp scaler (°C)
                        kL = 25,    # lower exponent
                        kU = 29) {  # upper exponent
  term1 <- (1 - exp(-temp / TL))^kL
  term2 <- (exp(-temp / TU))^kU
  term1 * term2
}
tdm_none     <- function(temp) rep(1, length(temp))
tdm_const    <- function(temp, S0=0.26) rep(S0, length(temp))

# -- Egg & Alevin development --
egg_models <- list(
  linear  = function(T)  958  / T,
  mech    = function(T)  1e3 * 200^0.2 / (pmax(T+3,0.1)^1.5),
  jensen  = function(T)  495.28 * T^-0.136,
  beacham = function(T)  1/exp(10.404 - 2.043*log(T+7.575)),
  usgs    = function(T)  1/exp(6.872  -    log(T+0.332))
)
hatch_models <- list(
  linear  = function(T)  417  / T,
  mech    = function(T)  2.5  * (pmax(T+2,0.1)^(-0.5))
)

# ───────────────────────────────────────────────────────────────────────────────
# 2) STATIC DAILY SURVIVAL CURVES (Figure) -------------------------------------
# ───────────────────────────────────────────────────────────────────────────────

temps <- seq(2,18, by=0.1)
plot_df <- expand_grid(
  temp  = temps,
  model = c("lin_martin","lin_anderson",
            "exp_WF_egg","exp_WF_alevin",
            "exp_ZC_egg","exp_ZC_alevin",
            "exp_SM_egg","exp_SM_alevin",
            "weibull","none","const")
) %>%
  mutate(
    surv = case_when(
      model=="lin_martin"   ~ tdm_lin_martin(temp),
      model=="lin_anderson" ~ tdm_lin_anderson(temp, day=1),  # any day<=4
      str_detect(model,"^exp_WF_egg")   ~ tdm_exp(temp,"WaterForum2020","egg"),
      str_detect(model,"^exp_WF_alevin")~ tdm_exp(temp,"WaterForum2020","alevin"),
      str_detect(model,"^exp_ZC_egg")   ~ tdm_exp(temp,"ZeugCramer2012","egg"),
      str_detect(model,"^exp_ZC_alevin")~ tdm_exp(temp,"ZeugCramer2012","alevin"),
      str_detect(model,"^exp_SM_egg")   ~ tdm_exp(temp,"SALMOD2006","egg"),
      str_detect(model,"^exp_SM_alevin")~ tdm_exp(temp,"SALMOD2006","alevin"),
      model=="weibull" ~ tdm_weibull(temp),
      model=="none"    ~ tdm_none(temp),
      model=="const"   ~ tdm_const(temp)
    )
  )

ggplot(plot_df, aes(temp, surv, color=model)) +
  geom_line(size=1) +
  scale_y_continuous(lim=c(0,1)) +
  labs(title="Parameterized Daily Survival vs Temp",
       x="Temp (°C)", y="Daily Survival") +
  theme_minimal()

# ───────────────────────────────────────────────────────────────────────────────
# 3) SIMULATION: COMBINED EGG→ALEVIN SURVIVAL & SUMMARY TABLE -------------------
# ───────────────────────────────────────────────────────────────────────────────

set.seed(2025)

# ───────────────────────────────────────────────────────────────────────────────
# Parameters
# ───────────────────────────────────────────────────────────────────────────────
years      <- 2024:(2024 + 2 - 1)   # e.g. 2024–2073
n_redds    <- 1000
peak_date  <- "12-01"                # mm-dd for the normal spawn peak
sigma_days <- 20                     # ~3-week std-dev in spawning
start_temp <- 19;   min_temp <- 8;   end_temp <- 14
noise_sd   <- 0.5

# ───────────────────────────────────────────────────────────────────────────────
# Build a list of 50 years of data
# ───────────────────────────────────────────────────────────────────────────────
year_list <- map(years, function(yy){
  # 1) define the Oct1→Apr30 window for this water-year
  start_date <- ymd(sprintf("%d-10-01", yy))
  feb28      <- ymd(sprintf("%d-02-28", yy+1))
  end_date   <- ymd(sprintf("%d-04-30", yy+1))
  dates      <- seq(start_date, end_date, by = "day")
  n1 <- as.numeric(feb28 - start_date)
  n2 <- as.numeric(end_date - feb28)
  d  <- as.numeric(dates - start_date)
  
  # 2) piecewise linear seasonal baseline
  base_temp <- ifelse(
    dates <= feb28,
    start_temp + (min_temp - start_temp) * (d / n1),
    min_temp   + (end_temp - min_temp)   * ((d - n1) / n2)
  )
  # add some day-to-day noise
  temps <- base_temp + rnorm(length(dates), sd = noise_sd)
  
  # 3) sample redd dates from a normal bell around Dec 1
  peak_idx <- as.numeric(ymd(sprintf("%d-%s", yy, peak_date)) - start_date)
  offs     <- round(rnorm(n_redds, mean = peak_idx, sd = sigma_days))
  offs     <- pmin(pmax(offs, 0), length(dates)-1)  # clamp into [0, N-1]
  redds    <- start_date + offs
  
  # 4) return a small list for this year
  list(
    year      = yy,
    dates     = dates,
    temps     = temps,
    redds     = redds
  )
})

# convert to a named list for easy lookup
names(year_list) <- as.character(years)

# ───────────────────────────────────────────────────────────────────────────────
# Example: plot one year’s temperature & redd histogram
# ───────────────────────────────────────────────────────────────────────────────
yr2024 <- year_list[["2024"]]

# temp series
ggplot(tibble(date = yr2024$dates, temp = yr2024$temps),
       aes(date, temp)) +
  geom_line() +
  labs(title="Example Oct→Apr Temp (2024–2025)",
       x="", y="Temperature (°C)") +
  theme_minimal()

# redd histogram
ggplot(tibble(redd=yr2024$redds), aes(redd)) +
  geom_histogram(binwidth=7, fill="steelblue", color="white") +
  labs(title="Spawn‐date Distribution (2024–2025)",
       x="Spawn Date", y="Count") +
  theme_minimal()

simulate_combo <- function(variant, dates, redds, daily_surv) {
  survs <- vapply(
    redds, 
    function(sp) {
      i0 <- match(sp, dates)
      Te <- daily_surv[[variant]][i0]    # spawn‐day survival
      # compute durations once per redd:
      days_egg  <- egg_fn(Te); days_alev <- hatch_fn(Te)
      total_days<- ceiling(days_egg + days_alev)
      
      # just take the product of the precomputed curve slice:
      prod(daily_surv[[variant]][i0 + seq_len(total_days)-1], na.rm=TRUE)
    },
    numeric(1)
  )
  survs
}

# ───────────────────────────────────────────────────────────────────────────────
# 4) DENSITY‐DEPENDENCE FUNCTIONS ------------------------------------------------
# ───────────────────────────────────────────────────────────────────────────────

density_fun <- function(dd_type,
                        N,
                        S0  = 0.347,       # base survival (Sacramento)
                        mu1 = -0.0000188,  # linear redd effect
                        K   = 50000       # Beverton–Holt K
) {
  switch(dd_type,
         none           = 1,
         linear         = S0 + mu1 * N,
         `beverton-holt` = S0 / (1 + N / K)
  )
}

dd_tab <- tibble(
  dd_type   = c("none","linear","beverton-holt"),
  dd_factor = map_dbl(dd_type, ~ density_fun(.x, N=50000))
)

# ───────────────────────────────────────────────────────────────────────────────
# 5) FULL SUMMARY: TDM × Egg × Hatch × Density ----------------------------------
# ───────────────────────────────────────────────────────────────────────────────

# 2) list out each TDM variant with its metadata:
tdm_defs <- tribble(
  ~model,        ~calib,             ~stage,   ~variant,
  "exp",         "WaterForum2020",   "egg",    "exp_WF_egg",
  "exp",         "WaterForum2020",   "alevin", "exp_WF_alevin",
  "exp",         "ZeugCramer2012",   "egg",    "exp_ZC_egg",
  "exp",         "ZeugCramer2012",   "alevin", "exp_ZC_alevin",
  "exp",         "SALMOD2006",       "egg",    "exp_SM_egg",
  "exp",         "SALMOD2006",       "alevin", "exp_SM_alevin",
  "lin_martin",  NA,                 NA,       "lin_martin",
  "lin_anderson",NA,                 NA,       "lin_anderson",
  "weibull",     NA,                 NA,       "weibull",
  "none",        NA,                 NA,       "none",
  "const",       NA,                 NA,       "const"
)


plan(multisession, workers = 4)

all_years_summary <- with_progress({
  p <- progressor(along = year_list)
  map_dfr(year_list, function(yr) {
    p()   # advance the progress bar by 1
    
  # unpack the year’s data
  dates    <- yr$dates
  t_series <- yr$temps
  redds    <- yr$redds
  
  # ── 1) Precompute daily survival curves for each TDM variant ──────────────
  daily_surv <- map(tdm_defs$variant, function(v) {
    row   <- filter(tdm_defs, variant == v)
    model <- row$model[[1]]
    calib <- row$calib[[1]]
    stage <- row$stage[[1]]
    temps <- t_series
    
    switch(model,
           exp          = tdm_exp(temps, calib, stage),
           lin_martin   = tdm_lin_martin(temps),
           lin_anderson = tdm_lin_anderson(temps, day = seq_along(temps)),
           weibull      = tdm_weibull(temps),
           none         = tdm_none(temps),
           const        = tdm_const(temps)
    )
  }) %>% set_names(tdm_defs$variant)
  
  # ── 2) Now simulate each variant × egg_m × hatch_m × dd_type ─────────────
  summary_tbl_year <- tidyr::crossing(
    tdm_defs,
    egg_m   = names(egg_models),
    hatch_m = names(hatch_models),
    dd_type = c("none","linear","beverton-holt")
  ) %>%
    rowwise() %>%
    mutate(
      survs = list({
        # use the precomputed daily_surv vector for this variant
        ds   <- daily_surv[[variant]]
        # vectorized redd loop with vapply
        vapply(redds, function(sp) {
          i0 <- match(sp, dates)
          # days to hatch+emerge at spawn‐day temp
          days_egg  <- egg_models[[egg_m]](t_series[i0])
          days_alev <- hatch_models[[hatch_m]](t_series[i0])
          td        <- ceiling(days_egg + days_alev)
          slice     <- ds[i0 + seq_len(td) - 1]
          prod(slice, na.rm = TRUE)
        }, numeric(1))
      }),
      base_surv     = mean(survs, na.rm = TRUE),
      dd_factor     = density_fun(dd_type, N = length(redds)),
      mean_cum_surv = base_surv * dd_factor,
      dev0 = {
        Te  <- t_series[ match(redds, dates) ]
        de  <- egg_models[[egg_m]](Te) + hatch_models[[hatch_m]](Te)
        atu <- de * Te
        tibble(
          mean_dev_days = mean(de),
          mean_ATU      = mean(atu)
        )
      }
    ) %>%
    unnest_wider(dev0) %>%
    ungroup() %>%
    mutate(
      across(c(mean_dev_days, mean_ATU, mean_cum_surv), ~ round(.x,3)),
      year = yr$year
    ) %>%
    select(year, everything())
  
  summary_tbl_year
  })
})

# now you’ve got a big 50×N table
print(all_years_summary)

single_year <- year_list["2024"]

one_year_summary <- with_progress({
  p <- progressor(steps = length(single_year))
  map_dfr(single_year, function(yr) {
    p()
      
      # unpack the year’s data
      dates    <- yr$dates
      t_series <- yr$temps
      redds    <- yr$redds
      
      # ── 1) Precompute daily survival curves for each TDM variant ──────────────
      daily_surv <- map(tdm_defs$variant, function(v) {
        row   <- filter(tdm_defs, variant == v)
        model <- row$model[[1]]
        calib <- row$calib[[1]]
        stage <- row$stage[[1]]
        temps <- t_series
        
        switch(model,
               exp          = tdm_exp(temps, calib, stage),
               lin_martin   = tdm_lin_martin(temps),
               lin_anderson = tdm_lin_anderson(temps, day = seq_along(temps)),
               weibull      = tdm_weibull(temps),
               none         = tdm_none(temps),
               const        = tdm_const(temps)
        )
      }) %>% set_names(tdm_defs$variant)
      
      # ── 2) Now simulate each variant × egg_m × hatch_m × dd_type ─────────────
      summary_tbl_year <- tidyr::crossing(
        tdm_defs,
        egg_m   = names(egg_models),
        hatch_m = names(hatch_models),
        dd_type = c("none","linear","beverton-holt")
      ) %>%
        rowwise() %>%
        mutate(
          survs = list({
            # use the precomputed daily_surv vector for this variant
            ds   <- daily_surv[[variant]]
            # vectorized redd loop with vapply
            vapply(redds, function(sp) {
              i0 <- match(sp, dates)
              # days to hatch+emerge at spawn‐day temp
              days_egg  <- egg_models[[egg_m]](t_series[i0])
              days_alev <- hatch_models[[hatch_m]](t_series[i0])
              td        <- ceiling(days_egg + days_alev)
              slice     <- ds[i0 + seq_len(td) - 1]
              prod(slice, na.rm = TRUE)
            }, numeric(1))
          }),
          base_surv     = mean(survs, na.rm = TRUE),
          dd_factor     = density_fun(dd_type, N = length(redds)),
          mean_cum_surv = base_surv * dd_factor,
          dev0 = {
            Te  <- t_series[ match(redds, dates) ]
            de  <- egg_models[[egg_m]](Te) + hatch_models[[hatch_m]](Te)
            atu <- de * Te
            tibble(
              mean_dev_days = mean(de),
              mean_ATU      = mean(atu)
            )
          }
        ) %>%
        unnest_wider(dev0) %>%
        ungroup() %>%
        mutate(
          across(c(mean_dev_days, mean_ATU, mean_cum_surv), ~ round(.x,3)),
          year = yr$year
        ) %>%
        select(year, everything())
      
      summary_tbl_year
  })
})

# 1) Plot mean cumulative survival vs. year, colored by TDM variant
ggplot(one_year_summary, aes(x = year, y = mean_cum_surv, color = variant)) +
  geom_line(size = 1) +
  facet_grid(egg_m ~ dd_type) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(
    x     = "Water Year",
    y     = "Mean Egg→Emergence Survival",
    color = "TDM Variant"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(angle = 45, hjust = 1)
  )

