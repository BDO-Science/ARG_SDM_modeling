rm(list=ls())

# ─── Libraries ────────────────────────────────────────────────────────────────
# library(tidyverse); 
library(dplyr); library(ggplot2); library(tidyr)
library(lubridate); library(furrr); library(data.table); 
# library(compiler); 
library(here); library(ggrepel); library(MASS); library(ordinal); library(ggridges)
source(here("SalmonCountR", "functions.R"))

# ── Read data inputs ──────────────────────────────────
env_ext_list <- readRDS(
  here("SalmonCountR", "app_data", "env_ext_list.rds")
)

#NOTE: env refers to the alternatives, ChatGPT stuck that name on and I have rolled with it
# --- env -> site map ---------------------------------------------------- 
env_sites <- purrr::imap_dfr(env_ext_list, ~ {
  tibble(env = as.character(.y), site = unique(.x$site))
}) %>% distinct(env, site)                                              

df_all <- readRDS(
  here("SalmonCountR", "app_data", "df_all.rds")
)

# TN: I think it's best to remove the excess rows here. I see that these NA rows
# have to be accounted for later, which can cleanly be remedied if we ignore them here
carcass_raw <- read.csv(
  here("SalmonCountR", "app_data", "carcassdet_1752789274_15.csv")
)

# TN: problem with parsing OK, 6 issues, 1 did parse but its wrong
esc_obs <- readr::read_csv(
  here("SalmonCountR", "app_data", "grandtab_1752793045_337.csv"),
  col_types = cols(
    `End Year of Monitoring Period` = col_character(),
    `Population Estimate`           = col_double()
  )
) %>%
  mutate(year = parse_number(`End Year of Monitoring Period`)) %>%
  filter(year >= 2011, year <= 2024) %>%
  rename(spawners = `Population Estimate`) %>%
  arrange(year)

# detect how many physical cores you have (on Windows this will report logical threads,
# so we subtract 1 to leave one for the OS)
n_threads <- parallel::detectCores(logical = TRUE) - 1

# but cap at your physical cores (6) to avoid oversubscription
n_workers <- min(6, n_threads)

plan(multisession, workers = n_workers)

# and for data.table’s own internal C threading:
setDTthreads(n_workers)

# Define TDM variants to run: exponential with two calibrations, and linear
tdm_defs <- tribble(
  ~model,       ~calib,           ~variant,
  "exp",        "WaterForum2020", "exp_WF",
  "exp",        "SALMOD2006",     "exp_SM",
  "lin_martin", NA,                "lin_Martin"
)

# site → temp vectors
site_temps_list <- map(env_ext_list, ~ split(.x$temp, .x$site))

# site → Date vectors
site_dates_list <- map(env_ext_list, ~ split(.x$Date, .x$site))

# site → named‐index lookup
date_idx_list <- map(site_dates_list, function(dlist) {
  map(dlist, ~ setNames(seq_along(.x), as.character(.x)))
})

# prepping carcass data (keep fork_length_mm)
obs_df <- carcass_raw %>%
  mutate(
    Date       = as.Date(surveydate),
    spawn_dt   = Date - days(7),
    brood_year = if_else(month(spawn_dt) >= 9, year(spawn_dt), year(spawn_dt) - 1), 
    site = case_when(
      section %in% c("NB","W","1a","1b","1a/1b","2") ~ "AveHazel",
      section %in% c("3")                            ~ "AveWatt",
      TRUE                                           ~ NA_character_
    )
  ) %>%
  # transmute avoids the masked select(); make forklength → fork_length
  dplyr::transmute(
    brood_year, spawn_dt, section,
    fork_length = .data[["forklength"]],
    site
  )

site_props <- obs_df %>%
  count(site) %>%
  mutate(prop = n / sum(n))
site_props

# ---- Assign 10-day spawn bins ----
# # need dummy date because fall-run spawn across calendar years, this simplifies things
# bin_defs <- data.frame(
#   period = paste0("p", 1:11),
#   start = as.Date(c("2000-10-10", "2000-10-20", "2000-10-30", "2000-11-09",
#                     "2000-11-19", "2000-11-29", "2000-12-09", "2000-12-19",
#                     "2000-12-29", "2001-01-08", "2001-01-18")),
#   end = as.Date(c("2000-10-19", "2000-10-29", "2000-11-08", "2000-11-18",
#                   "2000-11-28", "2000-12-08", "2000-12-18", "2000-12-28",
#                   "2001-01-07", "2001-01-17", "2001-01-27"))
# )

# TN: There are dates in obs_df that are earlier than 10-10 that do not get assigned.
# Going to assume "period" is just an arbitrary label
bin_defs <- data.frame(period = paste0("p", 1:12)) %>%
  mutate(
    start = as.Date("2000-10-05") + (row_number() - 1) * 10,
    end = start + (10 - 1)
  )

obs_df$spawn_bin <- sapply(obs_df$spawn_dt, function(d) {
  if (is.na(d)) return(NA_character_)  # <---- skip NAs safely
  # dummy_date <- as.Date(paste0("2000-", format(d, "%m-%d"))) 
  # TN: Need to account for Jan being the next year as specified in bin_defs
  dummy_date <- if_else(month(d) >= 9, 
                        as.Date(paste0("2000-", format(d, "%m-%d"))),
                        as.Date(paste0("2001-", format(d, "%m-%d"))))
  match_idx <- which(dummy_date >= bin_defs$start & dummy_date <= bin_defs$end)
  if (length(match_idx) == 1) bin_defs$period[match_idx] else NA
})
obs_df$spawn_bin <- factor(obs_df$spawn_bin, levels = bin_defs$period, ordered = TRUE)

# Aggregate Oct/Nov mean temps by year/site
# TN: This collapses all alternatives into a singular data frame, giving you duplications
# Is this intended or did you want to keep the "alt"?
oct_nov_temps <- map_dfr(names(env_ext_list), function(site) {
  dat <- env_ext_list[[site]]
  dat %>%
    filter(month(Date) %in% c(10, 11)) %>%
    mutate(month = month(Date), year = year(Date)) %>%
    group_by(site, year, month) %>%
    summarize(mean_temp = mean(temp, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = month, values_from = mean_temp, names_prefix = "month_") %>%
    rename(Oct = month_10, Nov = month_11)
})

# Join to obs_df by year and site
# TN: This gives a warning about many-to-many relationship as there is duplications in oct_nov_temps
# TN: obs_df is of the historical data right? oct_nov_temps currently has data from the past AND
# data that extends into the future. Is this what you were intending? I would assume having
# projected data would misrepresent what happened in the observed data?
obs_df <- left_join(obs_df, oct_nov_temps, by = c("brood_year" = "year", "site"))

oct_mean <- mean(obs_df$Oct, na.rm = TRUE)
oct_sd   <- sd(obs_df$Oct, na.rm = TRUE)
nov_mean <- mean(obs_df$Nov, na.rm = TRUE)
nov_sd   <- sd(obs_df$Nov, na.rm = TRUE)

obs_df <- obs_df %>%
  mutate(
    Oct_std = (Oct - oct_mean) / oct_sd,
    Nov_std = (Nov - nov_mean) / nov_sd
  )

# Ensure clean rows & factor coding
obs_fit <- obs_df %>%
  # TN: I generally like to know where my NAs are coming from. Think these are from 
  # carcass_raw and should be resolved if the last few lines of that file is removed from the read
  filter(!is.na(spawn_bin), !is.na(Oct_std), !is.na(Nov_std), !is.na(brood_year)) %>% 
  mutate(
    spawn_bin  = factor(spawn_bin, levels = bin_defs$period, ordered = TRUE),
    brood_year = factor(brood_year)
  )

# TN: Potential red flag. Again, oct_nov_temps from above mixes temp across the 10 alternatives
# and also uses projected data in the future to calculate mean/sd. Might not be
# giving accurate coefficients for the model below? Is there an alternative that is
# a "No alternative" or "baseline" option?
# clmm model with standardized mean monthly temps and a year random effect
spawn_clmm <- clmm(spawn_bin ~ Oct_std + Nov_std + (1 | brood_year),
                   data = obs_fit, link = "logit")
# TN: Some warnings here
summary(spawn_clmm)

# Extract fixed effects (β) and thresholds (ζ)
coef_all <- coef(spawn_clmm)                 # named vector
beta     <- coef_all[c("Oct_std","Nov_std")] # fixed slopes
zeta     <- unname(coef_all[!names(coef_all) %in% names(beta)])  # thresholds (length K-1)

# --- robust BLUP extraction (Best Linear Unbiased Predictor)
re_by  <- ranef(spawn_clmm)$brood_year              # a data.frame with rows = years
u_blup <- setNames(as.numeric(re_by[,1]), rownames(re_by))  # vector with names=years

yrs_hist <- as.character(2011:2024)
avail    <- intersect(names(u_blup), yrs_hist)
u_hat    <- if (length(avail)) u_blup[avail] else numeric(0)
u_hat    <- u_hat - mean(u_hat, na.rm = TRUE)

# AR(1) phi and sigma as before (with guards)
phi_raw <- tryCatch(acf(u_hat, plot = FALSE, lag.max = 1)$acf[2], error = function(e) NA_real_)
phi     <- ifelse(is.finite(phi_raw), phi_raw, 0)
phi     <- max(min(phi, 0.95), -0.95)

eps   <- diff(u_hat) - phi * head(u_hat, -1)
sigma <- sd(eps, na.rm = TRUE)
if (!is.finite(sigma) || sigma == 0) sigma <- sd(u_hat, na.rm = TRUE)
if (!is.finite(sigma) || sigma == 0) sigma <- 0.2   # final fallback

# pick u0 safely
u0 <- if (length(u_hat) && is.finite(tail(u_hat, 1))) tail(u_hat, 1) else 0

# --- creates a forward series of year effects `u` for future years
simulate_u <- function(years_vec, phi, sigma, u0 = 0) {
  n <- length(years_vec)
  if (n == 0) return(numeric())
  u <- numeric(n)
  u[1] <- ifelse(length(u0) == 1 && is.finite(u0), u0, 0)
  if (n > 1) {
    for (i in 2:n) u[i] <- phi * u[i-1] + rnorm(1, 0, sigma)
  }
  setNames(u, years_vec)
}

real_years <- 2011:2024 
n_calib <- length(real_years) # 14 
max_forecast <- 100 
n_sim <- n_calib + max_forecast # 114 
sim_years <- real_years[1] + seq(0, n_sim-1)
forecast_years <- (max(real_years)+1):(max(sim_years))
set.seed(42)
u_fc <- simulate_u(forecast_years, phi, sigma, u0 = u0)

# If you want to use pretty date labels for each bin in legend/labels:
pretty_labels <- paste0(
  bin_defs$period, ": ",
  format(bin_defs$start, "%b %d"), "–", format(bin_defs$end, "%b %d")
)
label_map <- setNames(pretty_labels, bin_defs$period)

# 1. Generate prediction data for plotting
oct_range <- seq(min(obs_df$Oct_std, na.rm=TRUE), max(obs_df$Oct_std, na.rm=TRUE), length.out = 100)
mean_Nov_std <- mean(obs_df$Nov_std, na.rm = TRUE)
newdata <- data.frame(Oct_std = oct_range, Nov_std = mean_Nov_std)
# Make sure to use predict from the ordinal package
# Function to predict probabilities from clmm model
predict_clmm <- function(model, newdata, include_re = FALSE, re_values = NULL) {
  
  # Extract fixed effect coefficients and thresholds
  all_coef <- coef(model)
  beta_names <- c("Oct_std", "Nov_std")  # your fixed effects
  beta <- all_coef[beta_names]
  thresholds <- all_coef[!names(all_coef) %in% beta_names]
  
  # Calculate linear predictor from fixed effects
  X <- as.matrix(newdata[, beta_names, drop = FALSE])
  linear_pred <- as.vector(X %*% beta)
  
  # Add random effects if specified
  if (include_re && !is.null(re_values)) {
    # re_values should be a vector matching rows of newdata
    linear_pred <- linear_pred + re_values
  }
  
  # Calculate probabilities for each spawn bin
  n_obs <- nrow(newdata)
  n_categories <- length(thresholds) + 1
  
  # Initialize probability matrix
  probs <- matrix(0, nrow = n_obs, ncol = n_categories)
  
  # Calculate cumulative probabilities
  cum_probs <- matrix(0, nrow = n_obs, ncol = length(thresholds))
  for (i in seq_along(thresholds)) {
    cum_probs[, i] <- plogis(thresholds[i] - linear_pred)
  }
  
  # Convert cumulative to individual category probabilities
  probs[, 1] <- cum_probs[, 1]  # First category
  
  if (length(thresholds) > 1) {
    for (j in 2:length(thresholds)) {
      probs[, j] <- cum_probs[, j] - cum_probs[, j-1]
    }
  }
  
  probs[, n_categories] <- 1 - cum_probs[, length(thresholds)]  # Last category
  
  # Add column names (spawn bins)
  colnames(probs) <- levels(model$model$spawn_bin)
  
  return(probs)
}

# Use the function
pred_probs <- predict_clmm(spawn_clmm, newdata = newdata)
pred_df <- as.data.frame(pred_probs)
pred_df$Oct_std <- oct_range

# 2. Pivot to long
pred_long <- pred_df %>%
  pivot_longer(-Oct_std, names_to = "spawn_bin", values_to = "probability") %>%
  mutate(
    spawn_bin = factor(spawn_bin, levels = bin_defs$period, labels = pretty_labels)
  )

# 3. For direct labels: get the last (rightmost) point for each line
label_df <- pred_long %>%
  group_by(spawn_bin) %>%
  slice_max(Oct_std, n = 1)

# 4. Plot: direct labels AND legend (color matches)
ggplot(pred_long, aes(x = Oct_std, y = probability, color = spawn_bin)) +
  geom_line(size = 1) +
  facet_wrap(~ spawn_bin, ncol = 3) +
  scale_color_viridis_d(name = "Spawn Bin", option = "C", guide = "none") +
  labs(
    title = "Predicted spawn bin probability by standardized October temperature",
    x = "Standardized October Temperature (Oct_std)",
    y = "Probability"
  ) +
  theme_minimal(base_size = 14)

ggplot(pred_long, aes(x = Oct_std, y = probability, color = spawn_bin)) +
  geom_line(size = 1) +
  scale_color_viridis_d(name = "Spawn Bin", option = "C") +
  labs(
    title = "Predicted spawn bin probability by standardized October temperature",
    x = "Standardized October Temperature (Oct_std)",
    y = "Probability"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

# Step 1: Create a new data frame for predictions over Nov_std
nov_grid <- tibble(
  Oct_std = 0,                              # Keep October temp at mean (0)
  Nov_std = seq(-2, 2, length.out = 200)    # Vary November
)

# Step 2: Get predicted probabilities for each bin from your model
nov_probs <- predict_clmm(spawn_clmm, newdata = nov_grid)

# Step 3: Convert to long format for ggplot
pred_nov_long <- nov_probs %>%
  as.data.frame() %>%
  mutate(Nov_std = nov_grid$Nov_std) %>%
  pivot_longer(
    cols = starts_with("p"),   # assuming bin names are like p1, p2...
    names_to = "spawn_bin",
    values_to = "probability"
  )

# Step 4: Plot, faceted by spawn_bin
ggplot(pred_nov_long, aes(x = Nov_std, y = probability, color = spawn_bin)) +
  geom_line(size = 1) +
  scale_color_viridis_d(name = "Spawn Bin", option = "C") +
  labs(
    title = "Predicted spawn bin probability by standardized November temperature",
    x = "Standardized November Temperature (Nov_std)",
    y = "Probability"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

# ── 0) Site↔env map (unique) ────────────────────────────────────────────────
site_env_map <- purrr::imap_dfr(env_ext_list, ~ tibble(env = as.character(.y),
                                                       site = unique(.x$site))) %>%
  distinct(site, env)

# ── 1) Historical redds with env (sim_actual) ───────────────────────────────
sim_actual <- obs_df %>%
  filter(brood_year %in% real_years) %>%
  left_join(site_env_map, by = "site") %>%
  transmute(env, sim_year = brood_year, spawn_dt, site)

# ── 2) Forecast temps (env × year Oct/Nov means + z-scores) ────────────────
forecast_years <- (max(real_years)+1):(max(real_years)+max_forecast)  # or your existing vector

forecast_temps <- map_dfr(names(env_ext_list), function(env_nm) {
  env_ext_list[[env_nm]] %>%
    filter(year(Date) %in% forecast_years, month(Date) %in% c(10,11)) %>%
    mutate(sim_year = year(Date), month = month(Date)) %>%
    group_by(sim_year, month) %>%
    summarise(mean_temp = mean(temp, na.rm=TRUE), .groups="drop") %>%
    pivot_wider(names_from = month, values_from = mean_temp, names_prefix="m") %>%
    transmute(
      env      = env_nm,
      sim_year = sim_year,
      Oct      = m10,
      Nov      = m11
    ) %>%
    mutate(
      Oct_std = (Oct - oct_mean) / oct_sd,
      Nov_std = (Nov - nov_mean) / nov_sd
    )
}) %>% mutate(env = as.character(env))

# ── 3) Pools, counts, AR(1) year effect, bin probabilities ─────────────────
obs_df_env <- obs_df %>% left_join(site_env_map, by = "site")

env_nredd <- obs_df_env %>%
  filter(brood_year %in% real_years, !is.na(env)) %>%
  count(env, brood_year, name = "n") %>%
  group_by(env) %>%
  summarise(N_redd = max(1L, round(mean(n, na.rm = TRUE))), .groups = "drop")

env_pools <- obs_df_env %>%
  filter(brood_year %in% real_years, !is.na(env)) %>%
  group_by(env) %>%
  summarise(
    section_pool     = list(section),
    fork_length_pool = list(fork_length),
    site_pool        = list(unique(site)),
    .groups = "drop"
  )

# AR(1) forward year effects for forecast years
u_fc <- simulate_u(sort(unique(forecast_temps$sim_year)), phi, sigma)

present_bins <- levels(droplevels(obs_df$spawn_bin))
all_bins     <- bin_defs$period

x_beta <- function(Oct_std, Nov_std, beta) beta["Oct_std"]*Oct_std + beta["Nov_std"]*Nov_std
cumlogit_to_probs <- function(eta_row) {
  cp <- plogis(eta_row)
  p  <- c(cp, 1) - c(0, cp)
  p / sum(p)
}

probs_list <- vector("list", nrow(forecast_temps))
for (i in seq_len(nrow(forecast_temps))) {
  row <- forecast_temps[i, ]
  xb  <- x_beta(row$Oct_std, row$Nov_std, beta)
  uy  <- u_fc[as.character(row$sim_year)]
  eta <- zeta - xb - uy
  p_i <- cumlogit_to_probs(eta)
  p_i[!is.finite(p_i) | p_i < 0] <- 0
  s <- sum(p_i)
  probs_list[[i]] <- if (s > 0) p_i / s else rep(1/length(present_bins), length(present_bins))
}
probs_present <- do.call(rbind, probs_list)
colnames(probs_present) <- present_bins

probs_all <- matrix(0, nrow = nrow(forecast_temps), ncol = length(all_bins))
colnames(probs_all) <- all_bins
probs_all[, present_bins] <- probs_present
rs <- rowSums(probs_all); rs[!is.finite(rs) | rs <= 0] <- 1
probs_all <- probs_all / rs

# ── 4) Fast spawn‐date sampler pieces ───────────────────────────────────────
bin_tbl <- bin_defs %>%
  transmute(
    period,
    start_m = month(start), start_d = mday(start), start_yoff = ifelse(start_m >= 10, 0L, 1L),
    end_m   = month(end),   end_d   = mday(end),   end_yoff   = ifelse(end_m   >= 10, 0L, 1L)
  )
sample_dates_fast <- function(bins_chr, year_int) {
  idx <- match(bins_chr, bin_tbl$period)
  sm  <- bin_tbl$start_m[idx]; sd <- bin_tbl$start_d[idx]; sy <- year_int + bin_tbl$start_yoff[idx]
  em  <- bin_tbl$end_m[idx];   ed <- bin_tbl$end_d[idx];   ey <- year_int + bin_tbl$end_yoff[idx]
  start <- as.Date(sprintf("%04d-%02d-%02d", sy, sm, sd))
  end   <- as.Date(sprintf("%04d-%02d-%02d", ey, em, ed))
  len <- as.integer(end - start) + 1L; len[len <= 0 | is.na(len)] <- 1L
  off <- if (length(len)) floor(runif(length(len), min = 0, max = pmax(1L, len))) else integer()
  start + off
}

# ── 5) Build sim_future (this MUST happen before sim_redds!) ────────────────
ft_aug <- forecast_temps %>%
  left_join(env_nredd, by = "env") %>%
  left_join(env_pools, by = "env")

out_list <- vector("list", nrow(ft_aug))
for (i in seq_len(nrow(ft_aug))) {
  env_i <- ft_aug$env[i]
  yr_i  <- ft_aug$sim_year[i]
  n_i   <- ft_aug$N_redd[i]; if (is.na(n_i) || n_i <= 0) n_i <- 1L
  p_i   <- probs_all[i, ]; p_i <- if (sum(p_i) > 0) p_i/sum(p_i) else rep(1/length(p_i), length(p_i))
  bins  <- sample(bin_defs$period, n_i, replace = TRUE, prob = p_i)
  dates <- sample_dates_fast(bins, yr_i)
  sec_pool <- ft_aug$section_pool[[i]]; if (is.null(sec_pool)) sec_pool <- NA
  fl_pool  <- ft_aug$fork_length_pool[[i]]; if (is.null(fl_pool)) fl_pool <- NA_real_
  st_pool  <- ft_aug$site_pool[[i]]
  if (is.null(st_pool)) {
    st_pool <- env_sites %>% filter(env == env_i) %>% pull(site) %>% unique()
  }
  out_list[[i]] <- tibble(
    env         = env_i,
    sim_year    = yr_i,
    spawn_dt    = dates,
    section     = sample(sec_pool, n_i, replace = TRUE),
    fork_length = sample(fl_pool,  n_i, replace = TRUE),
    site        = sample(st_pool,  n_i, replace = TRUE)
  )
}
sim_future <- bind_rows(out_list)

# ── 6) NOW build sim_redds and any splits/medians ───────────────────────────
sim_redds <- bind_rows(sim_actual, sim_future) %>%
  arrange(env, sim_year, spawn_dt) %>%
  dplyr::select(env, sim_year, spawn_dt, site)

if (!data.table::is.data.table(sim_redds)) data.table::setDT(sim_redds)
data.table::setkey(sim_redds, env, sim_year)

# split by env-year (data.table method) and drop the key cols inside each chunk
sim_redds_split <- split(sim_redds, by = c("env","sim_year"), drop = TRUE, keep.by = FALSE)
# keep only what's needed inside each chunk
sim_redds_split <- lapply(sim_redds_split, function(dt) dt[, .(spawn_dt, site)])

# medians per env-year (for spawn_dates_by_env)
spawn_medians_env_year <- sim_redds %>%
  as_tibble() %>%
  group_by(env, sim_year) %>%
  summarise(
    spawn_dt = as.Date(median(as.numeric(spawn_dt), na.rm = TRUE), origin = "1970-01-01"),
    .groups = "drop"
  )


# helper: produce a Date vector of length(sim_years) with LOCF/backfill
build_spawn_vec_for_env <- function(df_env, sim_years) {                 # NEW
  m <- match(sim_years, df_env$sim_year)
  v <- df_env$spawn_dt[m]
  # carry forward missing
  for (i in seq_along(v)) if (is.na(v[i]) && i > 1) v[i] <- v[i-1]
  # backfill from first known
  if (anyNA(v)) {
    first_ok <- which(!is.na(v))[1]; if (!is.na(first_ok)) v[seq_len(first_ok-1)] <- v[first_ok]
  }
  # enforce month/day into target years
  as.Date(sprintf("%04d-%02d-%02d", sim_years, lubridate::month(v), lubridate::day(v)))
}

# env -> spawn_dates vector (aligned to sim_years)
spawn_dates_by_env <- purrr::map(
  setNames(nm = unique(spawn_medians_env_year$env)),
  ~build_spawn_vec_for_env(
    df_env = dplyr::filter(spawn_medians_env_year, env == .x),
    sim_years = sim_years
  )
)

##############################
#plotting spawning dates
##############################
sim_future_comp <- sim_future %>% mutate( dummy_date = as.Date(paste0("2000-", format(spawn_dt, "%m-%d"))) )

# 1) Recode Oct–Dec as year 2000 and Jan as 2001
sim_future_comp <- sim_future_comp %>%
  mutate(
    m = month(dummy_date),            # <-- use the column that is already a Date
    d = day(dummy_date),
    dummy_date_season = if_else(
      m >= 10,
      make_date(2000, m, d),         # Oct–Dec → 2000
      make_date(2001, m, d)          # Jan → 2001
    )
  )

ggplot(sim_future_comp, aes(x = dummy_date, color = factor(sim_year), group = sim_year)) +
  geom_density(alpha = 0.3) +
  labs(
    title = "Spawning Date Distributions (all years aligned to same season)",
    x = "Spawning Date (Oct–Jan, dummy year)",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)

ggplot(sim_future_comp, aes(x = dummy_date_season, y = factor(sim_year), fill = factor(sim_year))) +
  geom_density_ridges(scale = 3, rel_min_height = 0.01, alpha = 0.6) +
  scale_x_date(
    limits = c(as.Date("2000-10-01"), as.Date("2001-01-30")),
    date_labels = "%b %d"
  ) +
  labs(
    title = "Spawning Date Distributions by Year",
    x = "Spawning Date (Oct–Jan, dummy year)",
    y = "Simulation Year"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

summary_dates <- sim_future_comp %>%
  dplyr::group_by(sim_year) %>%
  dplyr::summarize(
    median_spawn = as.Date(median(as.numeric(spawn_dt), na.rm = TRUE), origin = "1970-01-01"),
    p10          = as.Date(quantile(as.numeric(spawn_dt), 0.10, na.rm = TRUE), origin = "1970-01-01"),
    p90          = as.Date(quantile(as.numeric(spawn_dt), 0.90, na.rm = TRUE), origin = "1970-01-01"),
    .groups = "drop"
  )

ggplot(summary_dates, aes(x = sim_year)) +
  geom_line(aes(y = median_spawn), color = "blue") +
  geom_ribbon(aes(ymin = p10, ymax = p90), alpha = 0.2, fill = "blue") +
  labs(
    title = "Median and Spread of Spawn Dates by Year",
    x = "Simulation Year",
    y = "Spawn Date"
  ) +
  theme_minimal(base_size = 14)

sim_future_comp %>%
  count(sim_year, dummy_date) %>%
  ggplot(aes(x = dummy_date, y = sim_year, fill = n)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(
    title = "Heatmap of Spawning Intensity by Date and Year",
    x = "Spawning Date (Oct–Jan, dummy year)",
    y = "Simulation Year"
  ) +
  theme_minimal(base_size = 14)

options(future.globals.maxSize = 4 * 1024^3)  # 4 GiB; pick a number safely above 0.75 GiB

plan(multisession, workers = min(4, n_workers))

# after building sim_years_vec <- 2011:20?? and real median dates:
spawn_dates_real <- sim_redds %>%
  arrange(sim_year) %>%
  group_by(sim_year) %>%
  summarize(spawn_dt = median(spawn_dt), .groups="drop") %>%
  pull(spawn_dt)

# repeat those month‐day pairs out to length n_sim:
spawn_dm <- rep(spawn_dates_real, length.out = n_sim)

# Generate a sequence of simulation years starting from the first real year:
sim_years_vec   <- (real_years[1] + seq_len(n_sim) - 1)


# ───────────────────────────────────────────────────────────────────────────────
# Plug into TDM loop exactly as before, but now sim_year is brood_year
# ───────────────────────────────────────────────────────────────────────────────

# 2) Compile “hot” functions ---------------------------------------------------
# compile the “hot” functions
hatch_model_c   <- cmpfun(hatch_model)
emergence_model_c <- cmpfun(emergence_model)
tdm_exp_c     <- cmpfun(tdm_exp)
tdm_lin_c     <- cmpfun(tdm_lin_martin)

# convert sim_redds in place to a data.table
setDT(sim_redds)

# split sim_redds once (fast subsetting later)
sim_redds_split <- split(sim_redds[, .(spawn_dt, site)], sim_redds$sim_year)

# cache env assets once
env_cache <- lapply(names(env_ext_list), function(env_nm) {
  list(
    env_nm  = env_nm,
    date_idx_env = date_idx_list[[env_nm]],
    temps_env    = site_temps_list[[env_nm]],
    date_min     = min(env_ext_list[[env_nm]]$Date),
    date_max     = max(env_ext_list[[env_nm]]$Date)
  )
})
names(env_cache) <- names(env_ext_list)

# ensure character (not factors)
tdm_defs <- tdm_defs %>%
  dplyr::mutate(
    model   = as.character(model),
    calib   = as.character(calib),
    variant = as.character(variant)
  )

results_obs_fast <- furrr::future_map_dfr(
  sim_years,
  function(sim_yr) {
    red_this <- sim_redds_split[[as.character(sim_yr)]]
    if (is.null(red_this) || nrow(red_this) == 0L) return(data.table::data.table())
    
    env_res <- lapply(names(env_cache), function(env_nm) {
      ec <- env_cache[[env_nm]]
      
      # cohort dates (vectorized)
      mm  <- data.table::month(red_this$spawn_dt)
      dd  <- lubridate::day(red_this$spawn_dt)
      yy  <- ifelse(mm >= 9L, sim_yr, sim_yr + 1L)
      rdr <- lubridate::make_date(year = yy, month = mm, day = dd)
      
      keep <- (rdr >= ec$date_min) & (rdr <= ec$date_max)
      if (!any(keep)) return(NULL)
      
      # ---- KEY SPEEDUP: collapse to unique (site, rdr) with weights ----
      pairs <- data.table::data.table(
        site = red_this$site[keep],
        rdr  = rdr[keep]
      )[, .N, by = .(site, rdr)]  # N = multiplicity (weight)
      
      # small helper: weighted mean survival for one variant
      variant_mean <- function(model_i, calib_i) {
        USE_ATU_WINDOW <- TRUE  # flip to FALSE to go back to constant-T window
        
        f_surv <- if (USE_ATU_WINDOW) compute_surv_by_atu else compute_surv
        survs <- mapply(
          FUN = f_surv,
          rdr  = pairs$rdr,
          site = pairs$site,
          MoreArgs = list(date_idx_env = ec$date_idx_env, temps_env = ec$temps_env, model = model_i, calib = calib_i),
          USE.NAMES = FALSE
        )
        # weight by how many redds shared that (site, rdr)
        sum(survs * pairs$N, na.rm = TRUE) / sum(pairs$N)
      }
      
      data.table::rbindlist(lapply(seq_len(nrow(tdm_defs)), function(i) {
        data.table::data.table(
          sim_year      = sim_yr,
          env           = env_nm,
          variant       = tdm_defs$variant[i],
          method        = "observed",
          mean_cum_surv = variant_mean(tdm_defs$model[i], tdm_defs$calib[i])
        )
      }), use.names = TRUE, fill = TRUE)
    })
    
    data.table::rbindlist(env_res, use.names = TRUE, fill = TRUE)
  },
  .options = furrr::furrr_options(seed = TRUE, scheduling = 20)  # coarse scheduling reduces overhead
)

# LIFE CYCLE MODEL: CALIBRATION & FORECAST  (no outmigration covariates)

# ---- Inputs assumed available ----
# results_obs_fast  : TDM output (env, variant, sim_year, mean_cum_surv)
# env_ext_list      : list of environments (just used to pick names)
# esc_obs           : observed spawners (2011–2024)
# real_years        : 2011:2024
# simulate_variant  : function using surv_vec and life-cycle params (no redd_covars)
# sim_forecast_fn   : function for forward sims (should not require redd_covars)

ref_env <- names(env_ext_list)[1]   # choose any env as your calibration reference

# 1) Summarize TDM survivals
egg_summary <- results_obs_fast %>%
  arrange(env, variant, sim_year) %>%
  group_by(env, variant, sim_year) %>%
  summarise(mean_cum_surv = mean(mean_cum_surv, na.rm = TRUE), .groups = "drop")

variant_names <- egg_summary %>% pull(variant) %>% unique() %>% sort()

# 2) Survival lookups (reference env for calibration; all envs for forecasts)
surv_lookup_by_variant <- results_obs_fast %>%
  filter(env == ref_env, sim_year %in% real_years) %>%
  arrange(variant, sim_year) %>%
  group_by(variant) %>%
  summarise(surv_vec = list(mean_cum_surv), .groups = "drop") %>%
  tibble::deframe()  # named list: variant -> numeric vector (length = length(real_years))

surv_lookup_full <- egg_summary %>%
  arrange(env, variant, sim_year) %>%
  group_by(env, variant) %>%
  summarise(surv_vec = list(mean_cum_surv), .groups = "drop") %>%
  mutate(key = paste(env, variant, sep = "_")) %>%
  dplyr::select(key, surv_vec) %>%
  tibble::deframe()

# 3) Life-cycle base parameters (no covariates)
base_P <- list(
  female_fraction = 0.5, fec = 5522, S0 = 0.347,
  K_spawners = 12493,
  SAR_mean = NA_real_, SAR_sd = 0.00237,
  lag_probs = c(`3` = 0.75, `4` = 0.249, `5` = 0.001),
  rear_surv = NA_real_
)

# 4) Calibration targets
obs_spawners  <- esc_obs$spawners                 # 2011–2024
S_seed_calib  <- obs_spawners[1:3]                # seeds for brood 2011–2013
n_calib       <- length(obs_spawners)             # 14
fit_idx       <- (length(S_seed_calib) + 1):n_calib

modular_sse_once <- function(par, variant) {
  P_tmp <- base_P
  P_tmp$SAR_mean  <- par[1]
  P_tmp$rear_surv <- par[2]
  
  out <- simulate_variant(
    surv_vec       = surv_lookup_by_variant[[variant]][1:n_calib],
    P              = P_tmp,
    years          = n_calib,
    S_init         = S_seed_calib,
    SAR_vec        = rep(par[1], n_calib),
    K_spawners_vec = rep(P_tmp$K_spawners, n_calib),
    deg_day_adult  = rep(0, n_calib),   # zeros in calibration
    sim_years_vec  = real_years
  )
  
  preds <- out$spawners[fit_idx]
  if (!all(is.finite(preds))) return(.Machine$double.xmax)
  sum((preds - obs_spawners[fit_idx])^2)
}

calib_results <- furrr::future_map_dfr(variant_names, function(v) {
  opt <- optim(
    par    = c(0.0025, 0.8),
    fn     = modular_sse_once,
    variant= v,
    method = "L-BFGS-B",
    lower  = c(0, 0),
    upper  = c(1, 1)
  )
  tibble::tibble(variant = v, SAR_mean = opt$par[1], rear_surv = opt$par[2], sse = opt$value)
}, .options = furrr::furrr_options(seed = TRUE))

# 5) Build parameter lists (copy calibrated values across envs)
base_P_list <- calib_results %>%
  split(.$variant) %>%
  purrr::map(function(df_v) {
    SARv  <- df_v$SAR_mean[1]
    rearv <- df_v$rear_surv[1]
    rlang::set_names(
      lapply(names(env_ext_list), function(env_nm) {
        P <- base_P
        P$SAR_mean  <- SARv
        P$rear_surv <- rearv
        P
      }),
      nm = names(env_ext_list)
    )
  })


# 6) Precompute Calibration-tab predictions (so Shiny doesn’t simulate) ----
calib_pred_by_variant <- rlang::set_names(
  lapply(variant_names, function(v) {
    P0 <- base_P_list[[v]][[ref_env]]
    out <- simulate_variant(
      surv_vec       = surv_lookup_by_variant[[v]][1:n_calib],
      P              = P0,
      years          = n_calib,
      S_init         = S_seed_calib,
      SAR_vec        = rep(P0$SAR_mean, n_calib),
      K_spawners_vec = rep(P0$K_spawners, n_calib),
      deg_day_adult  = rep(0, n_calib),
      sim_years_vec  = real_years
    )
    tibble(
      year      = real_years,
      observed  = esc_obs$spawners,
      predicted = out$spawners,
      SAR_mean  = P0$SAR_mean,
      rear_surv = P0$rear_surv,
      sse       = sum((out$spawners[fit_idx] - esc_obs$spawners[fit_idx])^2)
    )
  }),
  variant_names
)


# 7) Build forecast seeds (last 3 years of calibrated run) for each variant
env_seed <- ref_env
k <- length(S_seed_calib)

S_seed_fore_list <- rlang::set_names(
  purrr::map(calib_results$variant, function(v) {
    years_cal <- length(real_years)
    surv_vec  <- surv_lookup_full[[paste(ref_env, v, sep = "_")]][1:years_cal]
    Ptmp      <- base_P_list[[v]][[ref_env]]
    
    deg_day_cal <- compute_deg_day_adult(
      env_nm       = ref_env,
      sim_years    = real_years,
      spawn_dates  = spawn_dates_by_env[[ref_env]][match(real_years, sim_years)],
      env_ext_list = env_ext_list
    )
    
    out <- simulate_variant(
      surv_vec       = surv_vec,
      P              = Ptmp,
      years          = years_cal,
      S_init         = S_seed_calib,
      SAR_vec        = rep(Ptmp$SAR_mean, years_cal),
      K_spawners_vec = rep(Ptmp$K_spawners, years_cal),
      deg_day_adult  = deg_day_cal,
      sim_years_vec  = real_years
    )
    tail(out$spawners, length(S_seed_calib))
  }),
  calib_results$variant
)

# 8) Forecasts (modular only)
keys <- names(surv_lookup_full)

# default to false, can add stochasticity in shiny app
use_stochastic_SAR <- FALSE

# 9) Stochastic SAR options (kept for UI; not used in code above)
stoch_SAR_opts <- list(
  model       = "normal",      # or "lognormal", "beta", "gamma"
  mean        = base_P$SAR_mean,
  sd          = base_P$SAR_sd,
  shape1      = 2,             # for beta/gamma only
  shape2      = 5,
  timing      = "pulse",       # "all", "block", "pulse"
  block_years = 20:30,
  pulse_years = c(10, 15, 20, 25, 30, 35, 40),
  pulse_sd    = 0.002
)

# Now try running your code
results_full <- purrr::map_dfr(keys, function(key) {
  parts  <- strsplit(key, "_")[[1]]
  env_nm <- parts[1]
  var_nm <- paste(parts[-1], collapse = "_")
  
  seed_vec <- S_seed_fore_list[[var_nm]]
  
  sim_forecast_fn(
    var_nm,
    env_nm,
    flow_cfs = NULL,
    S_seed   = seed_vec,
    spawn_dates_by_env = spawn_dates_by_env   # <-- NEW
  )()
})


# 10) Save outputs used by Shiny
saveRDS(calib_results,         here("SalmonCountR","app_data","calib_results.rds"))
saveRDS(calib_pred_by_variant, here("SalmonCountR","app_data","calib_pred_by_variant.rds"))
saveRDS(results_full,       here("SalmonCountR","app_data","results_full.rds"))
saveRDS(egg_summary,        here("SalmonCountR","app_data","egg_summary.rds"))
saveRDS(surv_lookup_full,   here("SalmonCountR","app_data","surv_lookup_full.rds"))
saveRDS(base_P_list,        here("SalmonCountR","app_data","base_P_list.rds"))
saveRDS(S_seed_calib,       here("SalmonCountR","app_data","S_seed_calib.rds"))
saveRDS(S_seed_fore_list,   here("SalmonCountR","app_data","S_seed_fore_list.rds"))
saveRDS(stoch_SAR_opts,     here("SalmonCountR","app_data","stoch_SAR_opts.rds"))
saveRDS(sim_years,          here("SalmonCountR","app_data","sim_years.rds"))
saveRDS(spawn_dates_by_env, here("SalmonCountR","app_data","spawn_dates_by_env.rds"))  # NEW
