library(tidyverse)
library(pdftools)
library(stringr)
library(dplyr)
library(readxl)
library(lubridate)
library(dataRetrieval)

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