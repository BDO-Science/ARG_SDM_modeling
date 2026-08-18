library(tidyverse)
library(pdftools)
library(stringr)
library(dplyr)
library(readxl)
library(lubridate)
library(dataRetrieval)

rm(list=ls())

inriver_df <- tibble(
  Year = 1952:2024,
  Escapement = c(
    25000, 28000, 29000, 9000,  4900,  6832,  17300, 17900, 25200, 11200,
    14400, 37810, 38500, 25000, 18600, 18000, 26100, 44200, 28680, 41680,
    17459, 82242, 53596, 32132, 23159, 41605, 12929, 37315, 34259, 43462,
    33000, 26400, 27447, 56120, 49372, 39885, 24889, 19183,  5339, 17683,
    5911, 31027, 33598, 70618, 69745, 47195, 50457, 55339,100852,135384,
    124252,163742, 99230, 62679, 24540, 10120,  2514,  5211, 14689, 25626,
    38328, 58228, 26475, 15739, 14473,  9663, 21092, 27030, 22456, 11232,
    16383, 37321, 45541
  )
)

inriver_df

# ───────────────────────────────────────────────────────────────────────────────
# 1) READ IN YOUR FOUR SHEETS
# ───────────────────────────────────────────────────────────────────────────────
raw_catch_fp   <- "Raw Catch - Chinook.csv"
trap_eff_fp    <- "Trap Efficiency Summary.csv"
trap_ops_fp    <- "Trap Operations.csv"
env_fp         <- "Environmentals.csv"

raw_catch      <- read_csv(raw_catch_fp)
trap_eff       <- read_csv(trap_eff_fp)
trap_ops       <- read_csv(trap_ops_fp)
environmentals <- read_csv(env_fp)


# ───────────────────────────────────────────────────────────────────────────────
# 2) PROCESS TRAP OPERATIONS → daily effort proxy
# ───────────────────────────────────────────────────────────────────────────────
trap_ops2 <- trap_ops %>%
  mutate(
    date       = as_date(Date),
    mean_rpm   = (rpmRevolutionsAtStart + rpmRevolutionsAtEnd)/2,
    effort_rev = mean_rpm * counterAtEnd
  ) %>%
  select(date, effort_rev)

# ───────────────────────────────────────────────────────────────────────────────
# 1. Determine start and end date from trap_eff2
# ───────────────────────────────────────────────────────────────────────────────
start_date <- min(trap_ops2$date, na.rm = TRUE)
end_date   <- max(trap_ops2$date, na.rm = TRUE)

# ───────────────────────────────────────────────────────────────────────────────
# 3) FIT TRAP‐EFFICIENCY MODEL USING COVARIATES
#     (we use glm to model p_detect from flow — extend as needed)
# ───────────────────────────────────────────────────────────────────────────────

# filter and prep trap_eff data
trap_eff2 <- trap_eff %>%
  filter(
    markedTaxon == "Chinook salmon",
    is.na(includeTestComments) | includeTestComments == ""
  ) %>%
  mutate(
    date = as_date(releaseTime)
  ) 

# ───────────────────────────────────────────────────────────────────────────────
# 2. Download USGS flow data (discharge, parameter 00060)
#    Site 11446500 = American River near 
# ───────────────────────────────────────────────────────────────────────────────
site_number   <- "11446500"
parameter_cd  <- c("00060", "00010")  # discharge (cfs), water temperature

cov_data <- readNWISdata(
  sites = site_number,
  parameterCd = parameter_cd,
  startDate   = start_date,
  endDate     = end_date,
  service = "uv"
)

# ───────────────────────────────────────────────────────────────────────────────
# 3. Process flow data — average flow per day
# ───────────────────────────────────────────────────────────────────────────────
# Rename columns from NWIS
cov_daily <- cov_data %>%
  rename(
    flow_cfs = X_00060_00000,
    temp_C   = X_00010_00000
  ) %>%
  mutate(date = as_date(dateTime)) %>%
  group_by(date) %>%
  summarise(
    flow = mean(flow_cfs, na.rm = TRUE),
    temp = mean(temp_C,   na.rm = TRUE),
    .groups = "drop"
  )

# ───────────────────────────────────────────────────────────────────────────────
# 4. Join daily flow into trap_eff2
# ───────────────────────────────────────────────────────────────────────────────
trap_eff2 <- trap_eff2 %>%
  left_join(cov_daily, by = "date") %>%
  select(date, nRecaptured, nReleased, flow, temp) %>%
  mutate(
    year  = lubridate::year(date),
    month = lubridate::floor_date(date, "month"),
    flow_scaled = scale(flow)[, 1],
    temp_scaled = scale(temp)[, 1]
  )

library(lme4)
# fit trap efficiency model (can add more covariates like rpm, temp, etc.)
# Candidate model 1: intercept only
# Model 0: Intercept only (random effect for month)
m0 <- glmer(
  cbind(nRecaptured, nReleased - nRecaptured) ~ 1 + (1 | year),
  family = binomial,
  data = trap_eff2
)

# Model 1: Scaled flow only
m1 <- glmer(
  cbind(nRecaptured, nReleased - nRecaptured) ~ flow_scaled + (1 | year),
  family = binomial,
  data = trap_eff2
)

# Model 2: Scaled temp only
m2 <- glmer(
  cbind(nRecaptured, nReleased - nRecaptured) ~ temp_scaled + (1 | year),
  family = binomial,
  data = trap_eff2
)

# Model 3: Scaled flow + temp
m3 <- glmer(
  cbind(nRecaptured, nReleased - nRecaptured) ~ flow_scaled + temp_scaled + (1 | year),
  family = binomial,
  data = trap_eff2
)

library(bbmle)

AICtab(m0, m1, m2, m3, weights = TRUE)

# Compute scaling params from trap_eff2
flow_mean <- mean(trap_eff2$flow, na.rm = TRUE)
flow_sd   <- sd(trap_eff2$flow, na.rm = TRUE)

temp_mean <- mean(trap_eff2$temp, na.rm = TRUE)
temp_sd   <- sd(trap_eff2$temp, na.rm = TRUE)

# Apply same scaling to cov_daily
cov_predict <- cov_daily %>%
  mutate(
    flow_scaled = (flow - flow_mean) / flow_sd,
    temp_scaled = (temp - temp_mean) / temp_sd
  )

# ───────────────────────────────────────────────────────────────────────────────
# 4) PREDICT DAILY DETECTION PROBABILITY FROM ENVIRONMENTAL DATA
# ───────────────────────────────────────────────────────────────────────────────
p_detect_daily <- cov_predict %>%
  mutate(p_detect = predict(m3, newdata = ., type = "response", re.form = NA))

# ───────────────────────────────────────────────────────────────────────────────
# 5) AGGREGATE RAW CATCH TO DAILY AND MERGE WITH p_detect
# ───────────────────────────────────────────────────────────────────────────────
raw2 <- raw_catch %>%
  filter(finalRun == "Fall" & fishOrigin == "Natural") %>%
  mutate(date = as_date(Date)) %>%
  count(date, name = "catch_raw")  # daily catch count

# merge all: raw catch, detection probabilities, and trap effort (if needed)
daily <- raw2 %>%
  left_join(p_detect_daily, by = "date") %>%
  left_join(trap_ops2, by = "date") %>%
  distinct(date, .keep_all = TRUE) %>%  # ✅ remove duplicates by date
  mutate(
    juv_hat = catch_raw / p_detect
  )
  

# ───────────────────────────────────────────────────────────────────────────────
# 6) AGGREGATE TO ANNUAL TOTALS
# ───────────────────────────────────────────────────────────────────────────────
obs <- daily %>%
  mutate(year = year(date)) %>%
  group_by(year) %>%
  summarise(
    C_juv     = sum(catch_raw, na.rm = TRUE),
    N_juv_hat = sum(juv_hat, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  left_join(inriver_df, by = c("year" = "Year")) %>%
  rename(escapement = Escapement)

# ───────────────────────────────────────────────────────────────────────────────
# 7) BUILD DATA FOR STOCK-RECRUITMENT MODEL
# ───────────────────────────────────────────────────────────────────────────────
# Step 1: Shift juvenile year back to match brood year (spawning year)
recruits_df <- obs %>%
  transmute(
    brood_year = year - 1,  # juveniles were produced by spawners one year prior
    Recruits   = N_juv_hat
  )

# Step 2: Spawners as brood year
spawners_df <- obs %>%
  transmute(
    brood_year = year,         # year fish spawned
    Spawners   = escapement
  )

# Step 3: Join brood_year to align spawner → outmigrant cohort
d <- recruits_df %>%
  inner_join(spawners_df, by = "brood_year") %>%
  select(brood_year, Spawners, Recruits)

# Inspect
glimpse(d)



# ───────────────────────────────────────────────────────────────────────────────
# B) Negative‐log‐likelihood for Beverton–Holt (Mike’s version)
#    parameterized so that:
#      mu = a * S / (b + S)
#    where a = asymptotic max recruits,
#          b = half‐saturation constant,
#    and we fit on the log‐scale with lognormal error.
# ───────────────────────────────────────────────────────────────────────────────
nll_BVH <- function(loga, logb, logsigma) {
  # exponentiate parameters
  a     <- exp(loga)
  b     <- exp(logb)
  sigma <- exp(logsigma)
  
  # predicted log‐mean recruitment for each row of d
  #   log(mu) = log(a * S/(b + S))
  mu_log <- log(a * d$Spawners / (b + d$Spawners))
  
  # negative sum of log‐density of observed log‐recruits
  -sum(dnorm(log(d$Recruits),
             mean = mu_log,
             sd   = sigma,
             log  = TRUE))
}

# ───────────────────────────────────────────────────────────────────────────────
# C) Fit via MLE with bbmle
# ───────────────────────────────────────────────────────────────────────────────
library(bbmle)

# Random starting guesses on the log‐scale
start_list <- list(loga     = log(50000),  # start a ~ 50k
                   logb     = log(20000),  # start b ~ 20k
                   logsigma = log(0.2))    # start σ ~ 0.2

BVH_mle <- mle2(nll_BVH,
                start = start_list,
                method = "L-BFGS-B",
                # enforce positivity if you want
                lower = c(loga = -Inf, logb = -Inf, logsigma = -Inf),
                upper = c(loga =  Inf, logb =  Inf, logsigma =  Inf)
)

# View results
summary(BVH_mle)

# Extract MLEs on original scale
loga     <- coef(BVH_mle)["loga"]
logb     <- coef(BVH_mle)["logb"]
logsigma <- coef(BVH_mle)["logsigma"]

a_hat     <- exp(loga)       # asymptotic recruit capacity
b_hat     <- exp(logb)       # half‐saturation
sigma_hat <- exp(logsigma)   # ln‐error SD

# 95% CI (on original scale)
ci <- exp(confint(BVH_mle))

# ───────────────────────────────────────────────────────────────────────────────
# D) Build a handy BH function
# ───────────────────────────────────────────────────────────────────────────────
BVH_fun <- function(S) {
  # note: returns *mean* recruits (no error)
  a_hat * S / (b_hat + S)
}

# ───────────────────────────────────────────────────────────────────────────────
# E) Quick plot of fit
# ───────────────────────────────────────────────────────────────────────────────
plot(d$Spawners, d$Recruits,
     pch   = 21, bg = "steelblue", cex = 1.2,
     xlab  = "Escapement (Spawners)",
     ylab  = "Estimated Juvenile Outmigrants",
     main  = "Beverton–Holt Fit: R = a·S/(b + S)")

# overlay MLE curve
S_seq <- seq(min(d$Spawners), max(d$Spawners), length.out = 200)
lines(S_seq, BVH_fun(S_seq), col = "purple", lwd = 2)

# ───────────────────────────────────────────────────────────────────────────────
# Now you can plug ‘a_hat’, ‘b_hat’, and ‘sigma_hat’ back into your integrated
# Monte Carlo chain or use them to inform your Stan priors below.
# ───────────────────────────────────────────────────────────────────────────────

# 1. Build a function for mean predicted recruits
BH_predict <- function(S, a = a_hat, b = b_hat){
  a * S / (b + S)
}

# 2. Apply it to your data
obs <- obs %>%
  mutate(
    R_hat = BH_predict(N_juv_hat)        # mean predicted escapement
  )

# 3. (Optional) Add log‐normal process noise to get a predictive distribution
set.seed(42)
obs <- obs %>%
  rowwise() %>%
  mutate(
    R_pred = exp(rnorm(1,
                       mean = log(R_hat),
                       sd   = sigma_hat))
  ) %>%
  ungroup()

# 4. Plot observed vs. predicted
library(ggplot2)
ggplot(obs, aes(x = N_juv_hat)) +
  geom_point(aes(y = escapement), color = "steelblue", size = 3) +
  geom_point(aes(y = R_pred),   color = "purple",   size = 2, alpha = 0.6) +
  geom_line(aes(y = R_hat),      color = "darkred",  size = 1) +
  labs(
    x = "Estimated juv Out-migrants (N_juv_hat)",
    y = "Adult Escapement",
    title = "Observed (blue) vs. Predicted BH Mean (red) & One Random Draw (purple)"
  ) +
  theme_minimal()








# Adult escapement in year t
esc_df <- obs %>%
  transmute(
    year = year,
    Escapement = escapement
  )

# Juvenile outmigrants in year t-3
juv_df <- obs %>%
  transmute(
    year = year + 3,       # shift forward so they align with return year
    Outmigrants = N_juv_hat
  )

# Join to calculate survival
adult_surv_df <- esc_df %>%
  left_join(juv_df, by = "year") %>%
  mutate(
    S_adult = Escapement / Outmigrants
  ) %>%
  filter(!is.na(S_adult) & is.finite(S_adult))

  summary(adult_surv_df$S_adult)

# Inspect
glimpse(adult_surv_df)
