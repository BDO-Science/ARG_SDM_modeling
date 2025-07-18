rm(list = ls())
library(tidyverse)
library(lubridate)

set.seed(2025)
years      <- 2024:(2024+20-1)
n_redds    <- 1000
peak_date  <- "12-01"
sigma_days <- 20
start_temp <- 19; min_temp <- 8; end_temp <- 14; noise_sd <- 0.5

year_list <- map(years, function(yy) {
  start_date <- ymd(sprintf("%d-09-01", yy))
  feb28      <- ymd(sprintf("%d-02-28", yy+1))
  end_date   <- ymd(sprintf("%d-05-31", yy+1))
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


single_year <- year_list["2024"]
# ─────────────────────────────────────────────────
# make sure these exist in your session:
yr       <- single_year[[1]]
dates    <- yr$dates
temps <- yr$temps
peak_idx <- as.numeric(ymd(sprintf("%d-12-01", yr$year)) - min(dates))
sigma_days <- 20
n_redds  <- length(yr$redds)

# ─────────────────────────────────────────────────
# 2) Define only the TDM variants you actually want
tdm_defs <- tribble(
  ~model,       ~calib,           ~variant,
  "exp",       "WaterForum2020",  "exp_WF",
  "exp",       "SALMOD2006",      "exp_SM",
  "lin_martin", NA,               "lin_Martin"
)

# single egg‐only development + ATU hatch
egg_model   <- function(T)  958 / T      # Zeug et al. 2012
hatch_model <- function(T) 417 / T       # ATU hatching

# exponential TDM (egg‐only)
tdm_exp <- function(temps, calib) {
  params <- list(
    WaterForum2020 = list(α=3.40848e-11, β=1.21122),
    SALMOD2006     = list(α=1.475e-11,   β=1.392)
  )[[calib]]
  H <- sum(params$α * exp(params$β * temps))
  exp(-H)
}

# linear Martin TDM (egg‐only)
tdm_lin_martin <- function(temps, α=0.026, β=12.14) {
  H <- α * sum(pmax(temps - β, 0))
  exp(-H)
}

# ─────────────────────────────────────────────────
# 3) One‐year summary function
one_year_summary <- function(dates, temps, redds, tdm_defs) {
  tdm_defs %>% rowwise() %>% mutate(
    # loop over redd dates:
    survs = list(vapply(redds, function(sp) {
      i0       <- match(sp, dates)
      d_egg    <- egg_model(temps[i0])
      d_hatch  <- hatch_model(temps[i0])
      td       <- ceiling(d_egg + d_hatch)
      temp_sl  <- temps[i0 + seq_len(td) - 1]
      
      switch(model,
             exp        = tdm_exp(temp_sl, calib),
             lin_martin = tdm_lin_martin(temp_sl)
      )
    }, numeric(1))),
    
    mean_cum_surv = mean(unlist(survs), na.rm=TRUE),
    year          = yr$year
  ) %>% ungroup()
}

# ─────────────────────────────────────────────────
# 4) Run N replicates with spawn‐date jitter
N <- 20
set.seed(42)

# parameters for shifting the peak, in days
shift_sd   <- 14    # how much the population mean can move
sigma_days <- 20    # within-replicate spread (as before)

results_shifted <- map_dfr(year_list, function(yr) {
  dates   <- yr$dates
  temps   <- yr$temps
  n_redds <- length(yr$redds)
  
  # compute original peak index
  peak_idx <- as.numeric(ymd(sprintf("%d-12-01", yr$year)) - min(dates))
  
  # 1) draw a replicate‐level shift of the spawn peak
  peak_shift <- rnorm(1, mean = 0, sd = shift_sd)
  spawn_mean <- peak_idx + peak_shift
  
  # 2) sample individual redds around that shifted mean
  offs   <- round(rnorm(n_redds, mean = spawn_mean, sd = sigma_days))
  offs   <- pmin(pmax(offs, 0), length(dates)-1)
  redds_j <- min(dates) + offs
  
  # 3) run your summary (no‐DD version)
  one_year_summary(dates, temps, redds_j, tdm_defs) %>%
    mutate(
      year       = yr$year,
      peak_shift = peak_shift
    )
})
# ─────────────────────────────────────────────────
# 1) Prepare data frames
temps_df <- tibble(date = dates, temp = temps)

yr <- single_year[[1]]   # first (and only) element of that list
spawn_df <- tibble(
  date = yr$redds
)

redds_rep_df <- map_dfr(seq_len(N), function(r) {
  offs <- round(rnorm(length(yr$redds), mean = peak_idx, sd = sigma_days))
  offs <- pmin(
    pmax(offs, peak_idx - 2 * sigma_days),
    peak_idx + 2 * sigma_days
  )
  tibble(rep = r, date = min(dates) + offs)
})

# 2) Plot daily temperature time series
p1 <- ggplot(temps_df, aes(date, temp)) +
  geom_line(size=0.8, color="steelblue") +
  labs(title="Daily Water Temperatures",
       x="Date", y="Temperature (°C)") +
  theme_minimal()
p1

# 3) Plot observed redd distribution
p2 <- ggplot(spawn_df, aes(date)) +
  geom_histogram(binwidth=7, fill="tomato", color="white", alpha=0.8) +
  labs(title="Observed Spawning Dates",
       x="Spawn Date", y="Count (per week)") +
  theme_minimal()
p2

# 4) Plot jittered spawn‐date distributions across replicates (density)
p3 <- ggplot(redds_rep_df, aes(date, ..density.., group=rep)) +
  geom_density(alpha=0.2, color=NA, fill="darkgreen") +
  labs(title=sprintf("Spawn‐Date Variability (N = %d replicates)", N),
       x="Date", y="Density") +
  theme_minimal()
p3

# 1) Build a data-frame of all redd dates + year
spawn_df <- map_dfr(year_list, function(yr) {
  tibble(
    year = yr$year,
    date = yr$redds
  )
})

# 2) Optionally convert to day-of-year for easier alignment
spawn_df <- spawn_df %>%
  mutate(doy = yday(date))

# 3a) Histogram of weekly spawn counts, faceted by year
ggplot(spawn_df, aes(x = date)) +
  geom_histogram(binwidth = 7, fill = "tomato", color = "white", alpha = 0.8) +
  facet_wrap(~ year, ncol = 5, scales = "free_y") +
  labs(
    title = "Weekly Spawning Distributions by Year",
    x     = "Spawn Date",
    y     = "Number of Redds (per week)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 3b) Alternatively: histogram of day-of-year, so all years share the same x-axis
ggplot(spawn_df, aes(x = doy)) +
  geom_histogram(binwidth = 7, fill = "darkgreen", color = "white", alpha = 0.8) +
  facet_wrap(~ year, ncol = 5, scales = "free_y") +
  labs(
    title = "Spawning Distributions (Day-of-Year) by Year",
    x     = "Day of Year",
    y     = "Number of Redds (per week)"
  ) +
  theme_minimal()

# 3c) Density curves of spawn timing for each year, all on one plot
ggplot(spawn_df, aes(x = doy, color = factor(year), group = year)) +
  geom_density(size = 1) +
  labs(
    title = "Spawn-Date Density by Year",
    x     = "Day of Year",
    y     = "Density",
    color = "Year"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
