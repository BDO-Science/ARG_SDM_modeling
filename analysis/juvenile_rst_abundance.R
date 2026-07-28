library(tidyverse)
library(pdftools)
library(stringr)
library(dplyr)
library(readxl)
library(lubridate)
library(dataRetrieval)
library(BTSPAS)

# ───────────────────────────────────────────────────────────────────────────────
# 1) READ IN YOUR FOUR SHEETS
# ───────────────────────────────────────────────────────────────────────────────
raw_catch <- read.csv(
  here("data_raw", "juvenile_data", "Raw Catch - Chinook.csv"),
  stringsAsFactors = FALSE
)

trap_eff  <- read.csv(
  here("data_raw", "juvenile_data", "Trap Efficiency Summary.csv"),
  stringsAsFactors = FALSE
)

trap_ops  <- read.csv(
  here("data_raw", "juvenile_data", "Trap Operations.csv"),
  stringsAsFactors = FALSE
)

environmentals      <- read.csv(
  here("data_raw", "juvenile_data", "Environmentals.csv"),
  stringsAsFactors = FALSE
)


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
    # parse the datetime string…
    datetime = ymd_hm(releaseTime, tz = "UTC"),
    # …then extract just the date
    date     = as_date(datetime)
  )

# verify
trap_eff2 %>% select(releaseTime, datetime, date) %>% head() 

#########################################
#Environmental data (could use later)
#########################################
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

# ───────────────────────────────────────────────────────────────────────────────
# 5) AGGREGATE RAW CATCH TO DAILY AND MERGE WITH EFFICIENCY TRIAL DATA
# ───────────────────────────────────────────────────────────────────────────────
raw2 <- raw_catch %>%
  filter(finalRun == "Fall" & fishOrigin == "Natural") %>%
  mutate(date = as_date(Date)) %>%
  count(date, name = "catch_raw")  # daily catch count

# merge all: raw catch, detection probabilities, and trap effort (if needed)
daily <- raw2 %>%
  left_join(trap_ops2, by = "date") %>%
  distinct(date, .keep_all = TRUE) %>%
  mutate(juv_hat = catch_raw / (effort_rev/100)) %>%
  # join your daily flow/temp
  left_join(cov_daily, by = "date")

# 1) Pick your week‐start.  Here I’ll use ISO weeks (Mon–Sun):
to_week <- function(d) floor_date(d, "week", week_start = 1)

# 2) Build a single weekly table:
weekly_tbl <- daily %>%             # your daily catch/effectiveness table
  mutate(week = to_week(date)) %>%
  group_by(week) %>%
  summarise(
    # un‐marked fish: sum of daily catch, effort‐corrected
    u2 = sum(catch_raw, na.rm = TRUE),
    # any additional covariates: mean flow, temp, etc
    flow  = mean(flow,  na.rm = TRUE),
    temp  = mean(temp,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Now join on your trap‐efficiency trials, also rolled up to weeks:
  left_join(
    trap_eff2 %>%
      mutate(week = to_week(date)) %>%
      group_by(week) %>%
      summarise(
        m2 = sum(nRecaptured, na.rm = TRUE),
        M1 = sum(nReleased,   na.rm = TRUE),
        .groups = "drop"
      ),
    by = "week"
  ) %>%
  # turn flow/temp into a covariate matrix:
  mutate(
    flow_s  = scale(flow)[,1],
    temp_s  = scale(temp)[,1]
  ) %>%
  arrange(week) %>%
  mutate(time = row_number())    # 1,2,3,… each weekly strata

# 1) Filter out weeks with no trial
weekly_fit_tbl <- weekly_tbl %>%
  filter(!is.na(m2) & !is.na(M1)) %>%
  arrange(week)

# 2) Pull your inputs
time_vec <- weekly_fit_tbl$time       # 1, 2, 3, … up to J
n1_vec   <- weekly_fit_tbl$M1
m2_vec   <- weekly_fit_tbl$m2
u2_vec   <- weekly_fit_tbl$u2
covp_mat <- as.matrix(weekly_fit_tbl %>% select(flow_s, temp_s))

# (Optional sanity check)
length(time_vec)  # should equal length(n1_vec), length(m2_vec), length(u2_vec), nrow(covp_mat)

# 3) Fit the diagonal‐error TP estimator
fit <- TimeStratPetersenDiagError_fit(
  time       = time_vec,
  n1         = n1_vec,
  m2         = m2_vec,
  u2         = u2_vec,
  logitP.cov = covp_mat,
  n.chains   = 3,
  n.iter     = 2000,
  n.burnin   = 500
)

# 4) Extract
U_post    <- fit$samples$U
run_total <- fit$samples$Utot

# Posterior density + 95% CI
plot(fit, parms = c("Utot"))

library(coda)

# extract the coda::mcmc.list
jags_mcmc <- as.mcmc.list(fit$jags)  

# Gelman–Rubin R̂
gelman.diag(jags_mcmc)

# Effective sample size
effectiveSize(jags_mcmc)

# Pairwise trace‐ and density‐plots
plot(jags_mcmc)

weekly_summ <- tibble(
  week = weekly_fit_tbl$week
) %>% 
  mutate(
    median = apply(U_post, 2, median),
    lwr    = apply(U_post, 2, quantile, .025),
    upr    = apply(U_post, 2, quantile, .975)
  )

# Plot
ggplot(weekly_summ, aes(week, median)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  labs(y = "Unmarked abundance (U_t)")

# 0) Inspect what U_post actually is
str(U_post)
dim(U_post)   # should show something like (n_iter, J), or NULL if it’s a vector

# 1) Coerce to a 2D matrix if needed
U_post_mat <- if (is.null(dim(U_post))) {
  # it’s a vector of length n_iter
  matrix(U_post, ncol = 1)
} else {
  # it already has dimensions (n_iter × J)
  U_post
}

# confirm
dim(U_post_mat)  # should be c(n_iter, J)

# a) extract the calendar year for each weekly stratum
years     <- year(weekly_fit_tbl$week)
yr_levels <- sort(unique(years))

# b) sum the weekly U's by year, for each iteration
post_yearly <- sapply(yr_levels, function(Y) {
  cols <- which(years == Y)
  # this always returns a length‐n_iter vector
  rowSums(U_post_mat[, cols, drop = FALSE])
})
colnames(post_yearly) <- yr_levels

# c) now post_yearly is n_iter × n_years
dim(post_yearly)  # should be c(n_iter, length(yr_levels))

# 4) Plot
ggplot(annual_summ, aes(year, median)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  labs(y = "Total run size per year", x = "Year")
