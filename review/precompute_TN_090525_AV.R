rm(list=ls())

# ─── Libraries ────────────────────────────────────────────────────────────────
library(tidyverse); library(lubridate); library(furrr); library(data.table); library(compiler); 
library(here); library(ggrepel); library(MASS); library(ordinal); library(ggridges)
source(here("SalmonCountR", "functions.R"))

set.seed(123)

real_years <- 2011:2024 
n_calib <- length(real_years) # 14 
max_forecast <- 100 
n_sim <- n_calib + max_forecast # 114 
sim_years <- real_years[1] + seq(0, n_sim-1)
forecast_years <- (max(real_years)+1):(max(sim_years))

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

trim_trailing_text <- function(df) {
  nnum <- df %>% mutate(.nn = rowSums(across(where(is.numeric), ~ !is.na(.)))) %>% pull(.nn)
  last <- suppressWarnings(max(which(nnum > 0), na.rm = TRUE))
  if (!is.finite(last)) return(df)
  df[seq_len(last), , drop = FALSE]
}

# TN: I think it's best to remove the excess rows here. I see that these NA rows
# have to be accounted for later, which can cleanly be remedied if we ignore them here
carcass_raw <- readr::read_csv(
  here::here("SalmonCountR", "app_data", "carcassdet_1752789274_15.csv")
) %>%
  trim_trailing_text()

# TN: problem with parsing OK, 6 issues, 1 did parse but its wrong
esc_obs <- read_csv(
  here("SalmonCountR", "app_data", "grandtab_1752793045_337.csv"),
  col_types = cols(
    `End Year of Monitoring Period` = col_character(),
    `Population Estimate`           = col_double()
  ),
  show_col_types = FALSE
) %>%
  trim_trailing_text() %>%                     # ← drop the footer notes
  mutate(year = parse_number(`End Year of Monitoring Period`)) %>%
  filter(between(year, 2011, 2024)) %>%
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

# ensure spawn_dt is Date
if (!inherits(obs_df$spawn_dt, "Date")) {
  obs_df <- obs_df %>% mutate(spawn_dt = as.Date(spawn_dt))
}

# 1) Anchor = earliest month-day seen (prefer Aug–Dec; fallback to all)
anchor_mmdd <- obs_df %>%
  filter(!is.na(spawn_dt)) %>%
  mutate(mmdd = format(spawn_dt, "%m-%d"), m = month(spawn_dt)) %>%
  filter(m >= 8) %>%
  summarise(anchor = min(mmdd, na.rm = TRUE)) %>%
  pull(anchor)

if (length(anchor_mmdd) == 0 || is.na(anchor_mmdd)) {
  anchor_mmdd <- obs_df %>%
    filter(!is.na(spawn_dt)) %>%
    mutate(mmdd = format(spawn_dt, "%m-%d")) %>%
    summarise(anchor = min(mmdd, na.rm = TRUE)) %>%
    pull(anchor)
}

# 2) Choose bin width, compute how many bins are needed (dynamic)
bin_width <- 10L
md <- format(obs_df$spawn_dt, "%m-%d")
season_year  <- year(obs_df$spawn_dt) - (md < anchor_mmdd)
season_start <- as.Date(paste0(season_year, "-", anchor_mmdd))
season_day   <- as.integer(obs_df$spawn_dt - season_start)
needed_bins  <- ceiling((max(season_day, na.rm = TRUE) + 1) / bin_width)
n_bins <- max(12L, needed_bins)   # at least 12, expand if season is longer

# 3) Assign ordered period factor p1..pN
assign_period <- function(dates, anchor_mmdd, bin_width, n_bins) {
  md <- format(dates, "%m-%d")
  season_year  <- year(dates) - (md < anchor_mmdd)
  season_start <- as.Date(paste0(season_year, "-", anchor_mmdd))
  season_day   <- as.integer(dates - season_start)
  bin_idx <- pmax(1L, pmin(floor(season_day / bin_width) + 1L, n_bins))
  factor(paste0("p", bin_idx), levels = paste0("p", seq_len(n_bins)), ordered = TRUE)
}

obs_df <- obs_df %>%
  mutate(period = assign_period(spawn_dt, anchor_mmdd, bin_width, n_bins))

# 4) (optional) Bin labels for docs/plots
bin_defs <- tibble(period = paste0("p", seq_len(n_bins))) %>%
  mutate(
    start = as.Date(paste0("2000-", anchor_mmdd)) + (row_number() - 1) * bin_width,
    end   = start + (bin_width - 1)
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

# years present in your observed data
yrs <- sort(unique(as.integer(obs_df$brood_year)))

#verifying that these observed temps are the same across alternatives like they should be
env_ext_list %>%
  imap_dfr(~ mutate(.x, alt = .y)) %>%
  filter(year(Date) %in% yrs) %>%
  group_by(site, Date) %>%
  summarise(n_temp = n_distinct(temp), .groups = "drop") %>%
  summarise(all_equal = all(n_temp == 1))


# Collapse env_ext_list across alts → one row per site–year–(Oct/Nov)
oct_nov_temps <- env_ext_list %>%
  imap_dfr(~ mutate(.x, alt = as.character(.y))) %>%          # tag alt (not used further)
  filter(year(Date) %in% yrs, month(Date) %in% c(10, 11)) %>%
  group_by(site, year = year(Date), month = month(Date)) %>%   # <— no 'alt' here: collapses alts
  summarise(mean_temp = mean(temp, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = month, values_from = mean_temp, names_prefix = "m_") %>%
  rename(Oct = m_10, Nov = m_11)

# sanity: ensure unique site–year (prevents many-to-many on join)
stopifnot(!any(duplicated(oct_nov_temps[c("site", "year")])))

# Join ONLY historical temps to obs_df
obs_df2 <- obs_df %>%
  mutate(brood_year = as.integer(brood_year)) %>%
  left_join(oct_nov_temps, by = c("site", "brood_year" = "year"))

# Standardize with training values
oct_mean <- mean(obs_df2$Oct, na.rm = TRUE);  oct_sd <- sd(obs_df2$Oct, na.rm = TRUE)
nov_mean <- mean(obs_df2$Nov, na.rm = TRUE);  nov_sd <- sd(obs_df2$Nov, na.rm = TRUE)

obs_fit <- obs_df2 %>%
  mutate(
    Oct_std = (Oct - oct_mean) / oct_sd,
    Nov_std = (Nov - nov_mean) / nov_sd,
    spawn_bin  = factor(spawn_bin, levels = bin_defs$period, ordered = TRUE),
    brood_year = factor(brood_year)
  ) %>%
  filter(!is.na(spawn_bin), !is.na(Oct_std), !is.na(Nov_std), !is.na(brood_year))

# ---- 0) Fit CLM (no random effects)
spawn_clm <- clm(spawn_bin ~ Oct_std + Nov_std,
                 data = obs_fit, link = "logit")
summary(spawn_clm)

# ---- 1) Pull coefficients and categories from model
cf   <- coef(spawn_clm)
beta <- cf[c("Oct_std","Nov_std")]
zeta <- unname(cf[!(names(cf) %in% names(beta))])
K    <- length(zeta) + 1

bins_model <- levels(spawn_clm$model$spawn_bin)
if (length(bins_model) != K) bins_model <- bin_defs$period[seq_len(K)]

# ---------- Minimal deterministic forecast to build sim_future ----------

# 0) Training scalers (if not already made)
if (!exists("sc")) {
  sc <- obs_df2 %>%
    filter(brood_year %in% yrs) %>%
    summarise(
      o_m = mean(Oct, na.rm = TRUE), o_s = sd(Oct, na.rm = TRUE),
      n_m = mean(Nov, na.rm = TRUE), n_s = sd(Nov, na.rm = TRUE)
    )
}

yrs_forecast   <- forecast_years
forecast_temps <- build_forecast_temps(env_ext_list, yrs_forecast, sc)
stopifnot(all(c("env","sim_year","Oct","Nov","Oct_std","Nov_std") %in% names(forecast_temps)))

# 2) Deterministic CLM probabilities (no noise/jitter)
probs_all <- predict_clm_probs(beta, zeta, forecast_temps[, c("Oct_std","Nov_std")], offset = 0)
colnames(probs_all) <- bins_model
rs <- rowSums(probs_all); rs[!is.finite(rs) | rs <= 0] <- 1
probs_all <- probs_all / rs

# Ensure columns match the bins present in observed data
present_bins <- levels(droplevels(obs_df$spawn_bin))
if (!identical(colnames(probs_all), present_bins)) {
  keep <- match(present_bins, colnames(probs_all))
  probs_all <- probs_all[, keep, drop = FALSE]
  colnames(probs_all) <- present_bins
}

# 3) Join counts and sampling pools
obs_df_env <- obs_df %>% left_join(
  purrr::imap_dfr(env_ext_list, ~ tibble(env = as.character(.y), site = unique(.x$site))) %>%
    distinct(site, env),
  by = "site"
)

env_nredd <- obs_df_env %>%
  filter(brood_year %in% as.integer(2011:2024), !is.na(env)) %>%
  count(env, brood_year, name = "n") %>%
  group_by(env) %>%
  summarise(N_redd = max(1L, round(mean(n, na.rm = TRUE))), .groups = "drop")

env_pools <- obs_df_env %>%
  filter(brood_year %in% as.integer(2011:2024), !is.na(env)) %>%
  group_by(env) %>%
  summarise(
    section_pool     = list(section),
    fork_length_pool = list(fork_length),
    site_pool        = list(unique(site)),
    .groups = "drop"
  )

ft_aug <- forecast_temps %>%
  left_join(env_nredd, by = "env") %>%
  left_join(env_pools, by = "env")

# 4) Date sampler per bin (define once)
if (!exists("sample_dates_fast")) {
  bin_tbl <- bin_defs %>%
    transmute(
      period,
      start_m = lubridate::month(start), start_d = lubridate::mday(start), start_yoff = ifelse(start_m >= 10, 0L, 1L),
      end_m   = lubridate::month(end),   end_d   = lubridate::mday(end),   end_yoff   = ifelse(end_m   >= 10, 0L, 1L)
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
}

# 5) Build sim_future (deterministic except for within-bin date draw)
out_list <- vector("list", nrow(ft_aug))
for (i in seq_len(nrow(ft_aug))) {
  env_i <- ft_aug$env[i]
  yr_i  <- ft_aug$sim_year[i]
  n_i   <- ft_aug$N_redd[i]; if (is.na(n_i) || n_i <= 0) n_i <- 1L
  
  p_i <- probs_all[i, ]  # already sums to 1
  bins_i  <- sample(present_bins, n_i, replace = TRUE, prob = p_i)
  dates_i <- sample_dates_fast(bins_i, yr_i)
  
  sec_pool <- ft_aug$section_pool[[i]]; if (is.null(sec_pool)) sec_pool <- NA
  fl_pool  <- ft_aug$fork_length_pool[[i]]; if (is.null(fl_pool)) fl_pool <- NA_real_
  st_pool  <- ft_aug$site_pool[[i]]
  if (is.null(st_pool)) {
    st_pool <- obs_df_env %>% filter(env == env_i) %>% pull(site) %>% unique()
  }
  
  out_list[[i]] <- tibble(
    env         = env_i,
    sim_year    = yr_i,
    spawn_dt    = dates_i,
    section     = sample(sec_pool, n_i, replace = TRUE),
    fork_length = sample(fl_pool,  n_i, replace = TRUE),
    site        = sample(st_pool,  n_i, replace = TRUE)
  )
}
sim_future <- dplyr::bind_rows(out_list)

# -------------------------------------------------------------------


# 0) Pretty labels p1..pK → "p#: Mon DD–Mon DD"
idx <- match(bins_model, bin_defs$period)
pretty_labels <- paste0(
  bins_model, ": ",
  format(bin_defs$start[idx], "%b %d"), "–", format(bin_defs$end[idx], "%b %d")
)

# 1) Safe ranges from TRAINING data used to fit CLM
safe_range <- function(x, fallback = c(-2, 2)) {
  x <- x[is.finite(x)]
  if (length(x)) range(x) else fallback
}
rng_oct <- safe_range(obs_fit$Oct_std)
rng_nov <- safe_range(obs_fit$Nov_std)

oct_range <- seq(rng_oct[1], rng_oct[2], length.out = 100)
nov_range <- seq(rng_nov[1], rng_nov[2], length.out = 200)

# 2) Predicted probabilities across Oct_std (Nov at its training mean)
newdata_oct <- data.frame(
  Oct_std = oct_range,
  Nov_std = mean(obs_fit$Nov_std, na.rm = TRUE)
)
pred_oct <- predict_clm_probs(beta, zeta, newdata_oct, offset = 0)
colnames(pred_oct) <- bins_model

pred_oct_long <- as.data.frame(pred_oct) %>%
  mutate(Oct_std = oct_range) %>%
  pivot_longer(all_of(bins_model), names_to = "spawn_bin", values_to = "probability") %>%
  mutate(spawn_bin = factor(spawn_bin, levels = bins_model, labels = pretty_labels))

# 3) Predicted probabilities across Nov_std (Oct at its training mean)
newdata_nov <- data.frame(
  Oct_std = mean(obs_fit$Oct_std, na.rm = TRUE),
  Nov_std = nov_range
)
pred_nov <- predict_clm_probs(beta, zeta, newdata_nov, offset = 0)
colnames(pred_nov) <- bins_model

pred_nov_long <- as.data.frame(pred_nov) %>%
  mutate(Nov_std = nov_range) %>%
  pivot_longer(all_of(bins_model), names_to = "spawn_bin", values_to = "probability") %>%
  mutate(spawn_bin = factor(spawn_bin, levels = bins_model, labels = pretty_labels))

# 4) PLOTS — facet like before
p_oct_facet <- ggplot(pred_oct_long, aes(x = Oct_std, y = probability, color = spawn_bin)) +
  geom_line(size = 1) +
  facet_wrap(~ spawn_bin, ncol = 3) +
  scale_color_viridis_d(guide = "none") +
  labs(
    title = "Predicted spawn-bin probability vs standardized October temperature",
    x = "Standardized October Temperature (Oct_std)",
    y = "Probability"
  ) +
  theme_minimal(base_size = 14)

p_nov_facet <- ggplot(pred_nov_long, aes(x = Nov_std, y = probability, color = spawn_bin)) +
  geom_line(size = 1) +
  facet_wrap(~ spawn_bin, ncol = 3) +
  scale_color_viridis_d(guide = "none") +
  labs(
    title = "Predicted spawn-bin probability vs standardized November temperature",
    x = "Standardized November Temperature (Nov_std)",
    y = "Probability"
  ) +
  theme_minimal(base_size = 14)

# 5) PLOTS — single panel with legend (also like you had)
p_oct_legend <- ggplot(pred_oct_long, aes(x = Oct_std, y = probability, color = spawn_bin)) +
  geom_line(size = 1) +
  scale_color_viridis_d(name = "Spawn Bin") +
  labs(
    title = "Predicted spawn-bin probability vs standardized October temperature",
    x = "Standardized October Temperature (Oct_std)",
    y = "Probability"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

p_nov_legend <- ggplot(pred_nov_long, aes(x = Nov_std, y = probability, color = spawn_bin)) +
  geom_line(size = 1) +
  scale_color_viridis_d(name = "Spawn Bin") +
  labs(
    title = "Predicted spawn-bin probability vs standardized November temperature",
    x = "Standardized November Temperature (Nov_std)",
    y = "Probability"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

# Print any of them:
p_oct_facet
p_oct_legend
p_nov_facet
p_nov_legend

# Forecast (you already built sim_future: env, sim_year, spawn_dt, ...)
fc_df <- sim_future %>%
  mutate(
    season_date = season_posix(spawn_dt, anchor_mmdd),
    source = "Forecast"
  )

# 1) Standardize temps (train scalers on 2011–2024), then get per-year means
sc <- obs_df2 %>%
  filter(brood_year %in% yrs) %>%
  summarise(o_m = mean(Oct, na.rm=TRUE), o_s = sd(Oct, na.rm=TRUE),
            n_m = mean(Nov, na.rm=TRUE), n_s = sd(Nov, na.rm=TRUE))
hist_temps <- obs_df2 %>%
  filter(brood_year %in% yrs) %>%
  group_by(brood_year) %>%
  summarise(Oct = mean(Oct, na.rm=TRUE), Nov = mean(Nov, na.rm=TRUE), .groups="drop") %>%
  mutate(Oct_std = (Oct - sc$o_m)/sc$o_s,
         Nov_std = (Nov - sc$n_m)/sc$n_s) %>%
  rename(sim_year = brood_year)

# 2) CLM → probabilities per year
cf   <- coef(spawn_clm)
beta <- cf[c("Oct_std","Nov_std")]
zeta <- unname(cf[!(names(cf) %in% names(beta))])
bins <- levels(spawn_clm$model$spawn_bin)

predict_clm_probs <- function(beta, zeta, newdata) {
  xb <- as.matrix(newdata[, names(beta), drop=FALSE]) %*% beta
  t(vapply(as.vector(xb), function(xi) {
    cp <- plogis(zeta - xi); p <- c(cp,1) - c(0,cp); p/sum(p)
  }, numeric(length(zeta)+1)))
}
probs_hist <- predict_clm_probs(beta, zeta, hist_temps[, c("Oct_std","Nov_std")])
colnames(probs_hist) <- bins

# 3) Sample modeled dates (match carcass counts per year)
n_by_year <- obs_df %>% filter(brood_year %in% yrs) %>% count(brood_year, name="N")
set.seed(42)
sim_pred <- lapply(seq_len(nrow(hist_temps)), function(i){
  yr <- hist_temps$sim_year[i]
  n  <- n_by_year$N[n_by_year$brood_year==yr]; if (length(n)==0 || is.na(n) || n<=0) n <- 1L
  p  <- probs_hist[i,]; p <- p/sum(p)
  b  <- sample(bins, n, replace=TRUE, prob=p)
  tibble(sim_year=yr, spawn_dt=sample_dates_fast(b, yr), source="Modeled")
}) %>% bind_rows()

# 4) Observed carcass dates
sim_actual <- obs_df %>% filter(brood_year %in% yrs) %>%
  transmute(sim_year = brood_year, spawn_dt, source="Observed")

# 5) Plots
# (A) Season-aligned ridgelines (anchor = earliest Aug–Dec date in observed)
anchor_mmdd <- sim_actual %>%
  mutate(m=month(spawn_dt), mmdd=format(spawn_dt,"%m-%d")) %>%
  filter(m>=8) %>% summarise(a=min(mmdd, na.rm=TRUE)) %>% pull(a)
if (!length(anchor_mmdd) || is.na(anchor_mmdd)) anchor_mmdd <- "10-05"
anchor_base <- as.Date(paste0("2000-", anchor_mmdd))
season_posix <- function(d, mmdd){
  md <- format(d, "%m-%d"); y0 <- year(d) - (md < mmdd)
  anchor <- as.Date(paste0(y0, "-", mmdd))
  anchor_base + as.integer(d - anchor)
}
comp_df <- bind_rows(sim_actual, sim_pred) %>%
  mutate(season_date = season_posix(spawn_dt, anchor_mmdd))

rng <- range(comp_df$season_date, na.rm=TRUE)
x_breaks <- seq(rng[1], rng[2], by = "10 days")

p_ridge <- ggplot(comp_df, aes(x=season_date, y=factor(sim_year))) +
  geom_density_ridges_gradient(aes(height=after_stat(density), fill=source),
                               scale=2, rel_min_height=0.001, linewidth=0.2,
                               bandwidth=10, trim=TRUE) +
  scale_x_date(breaks=x_breaks, date_labels="%b %d", expand=c(0.01,0)) +
  scale_fill_manual(values=c(Observed="steelblue", Modeled="tomato")) +
  labs(title="Observed vs Modeled spawn timing by year (2011–2024)",
       x=paste0("Season day (anchor ", anchor_mmdd, ")"), y="Year", fill=NULL) +
  theme_minimal(base_size=13)
p_ridge

# (B) Calendar-date boxplots (clear medians)
p_box <- ggplot(bind_rows(sim_actual, sim_pred),
                aes(x=factor(sim_year), y=spawn_dt, fill=source)) +
  geom_boxplot(alpha=0.65, position=position_dodge(width=0.7), outlier.alpha=0.3) +
  scale_fill_manual(values=c(Observed="steelblue", Modeled="tomato")) +
  labs(title="Observed vs Modeled spawning dates by year (2011–2024)",
       x="Brood year", y="Spawn date", fill=NULL) +
  theme_minimal(base_size=13)
p_box


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

envs <- sort(unique(spawn_medians_env_year$env))

spawn_dates_by_env <- map(
  .x = set_names(envs),
  .f = ~ build_spawn_vec_for_env(
    df_env    = filter(spawn_medians_env_year, env == .x),
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

# ==========================================
# 0) DROP invalid sites globally (no forecast temps)
# ==========================================
sites_with_temps <- Reduce(
  intersect,
  lapply(env_cache, function(ec) names(ec$date_idx_env))
)

sim_redds_split <- lapply(
  sim_redds_split,
  function(dt) {
    if (is.null(dt) || !nrow(dt)) return(dt)
    # Ensure required columns exist
    stopifnot(all(c("site", "spawn_dt") %in% names(dt)))
    dt[dt$site %in% sites_with_temps, , drop = FALSE]
  }
)

# ==========================================
# 1) BUILD per-env (site, rdr) -> pos lookup tables
# ==========================================
env_lookup <- lapply(env_cache, function(ec) {
  DT <- rbindlist(
    lapply(names(ec$date_idx_env), function(s) {
      v <- ec$date_idx_env[[s]]
      if (is.null(v) || !length(v)) {
        data.table(site = factor(character(0)), rdr = as.IDate(integer(0)), pos = integer(0))
      } else {
        data.table(
          site = s,
          rdr  = as.IDate(names(v)),   # names like "YYYY-MM-DD"
          pos  = as.integer(v)
        )
      }
    }),
    use.names = TRUE, fill = TRUE
  )
  
  DT[, site := as.factor(site)]
  setkey(DT, site, rdr)   # <- set key on the object, not .SD
  DT
})

# ==========================================
# 2) LIGHTWEIGHT MEMOIZER for survival calls
#    (Keyed by env|site|rdr|model|calib)
# ==========================================
.mem_env <- new.env(parent = emptyenv())

memo_surv <- function(env_nm, site, rdr, model, calib) {
  # rdr is IDate (integer), site is factor/character
  k <- paste(env_nm, site, as.integer(rdr), model, calib, sep = "|")
  hit <- .mem_env[[k]]
  if (!is.null(hit)) return(hit)
  
  ec <- env_cache[[env_nm]]
  out <- compute_surv_by_atu(
    rdr   = rdr,
    site  = site,
    date_idx_env = ec$date_idx_env,
    temps_env    = ec$temps_env,
    model = model,
    calib = calib
  )
  .mem_env[[k]] <- out
  out
}

# Vectorized facade (vapply) if compute_surv_by_atu isn't natively vectorized
compute_surv_vec <- function(env_nm, site, rdr, model, calib) {
  vapply(
    seq_along(site),
    function(i) memo_surv(env_nm, site[i], rdr[i], model, calib),
    numeric(1)
  )
}

# ==========================================
# 3) Build (site, rdr, N) for one env-year (with join)
# ==========================================
pairs_for_env_year <- function(red_this, sim_yr, env_nm) {
  lk <- env_lookup[[env_nm]]
  if (is.null(lk) || !nrow(lk)) return(NULL)
  
  # cohort rdr from spawn_dt (brood year → calendar)
  mm  <- data.table::month(red_this$spawn_dt)
  dd  <- lubridate::day(red_this$spawn_dt)
  yy  <- fifelse(mm >= 9L, sim_yr, sim_yr + 1L)
  rdr <- as.IDate(lubridate::make_date(yy, mm, dd))
  
  # collapse duplicates first
  pairs <- data.table(site = red_this$site, rdr = rdr)[
    , .N, by = .(site, rdr)
  ]
  
  if (!nrow(pairs)) return(NULL)
  
  # join to ensure (site,rdr) exists in the env index
  setkeyv(pairs, c("site", "rdr"))
  pairs <- lk[pairs, nomatch = 0L]  # keeps only valid pairs with 'pos'
  if (!nrow(pairs)) return(NULL)
  
  # columns: site, rdr, pos, N
  pairs[]
}

# ==========================================
# 4) eval_year(): compute all variants per env for one year
# ==========================================
eval_year <- function(sim_yr, sim_redds_split, env_cache, tdm_defs) {
  red_this <- sim_redds_split[[as.character(sim_yr)]]
  if (is.null(red_this) || !nrow(red_this)) return(data.table())
  
  env_names <- names(env_cache)
  
  env_tables <- lapply(env_names, function(env_nm) {
    pairs <- pairs_for_env_year(red_this, sim_yr, env_nm)
    if (is.null(pairs) || !nrow(pairs)) return(NULL)
    
    rbindlist(lapply(seq_len(nrow(tdm_defs)), function(i) {
      survs <- compute_surv_vec(
        env_nm = env_nm,
        site   = pairs$site,
        rdr    = pairs$rdr,
        model  = tdm_defs$model[i],
        calib  = tdm_defs$calib[i]
      )
      data.table(
        sim_year      = sim_yr,
        env           = env_nm,
        variant       = tdm_defs$variant[i],
        method        = "observed",
        mean_cum_surv = {
          wsum <- sum(pairs$N)
          if (!is.finite(wsum) || wsum <= 0) NA_real_ else sum(survs * pairs$N, na.rm = TRUE) / wsum
        }
      )
    }), use.names = TRUE, fill = TRUE)
  })
  
  env_tables <- Filter(Negate(is.null), env_tables)
  if (!length(env_tables)) return(data.table())
  rbindlist(env_tables, use.names = TRUE, fill = TRUE)
}

# ==========================================
# 5) Compile hot paths (optional) & set threads
# ==========================================
if (requireNamespace("compiler", quietly = TRUE)) {
  compute_surv_vec    <- compiler::cmpfun(compute_surv_vec)
  pairs_for_env_year  <- compiler::cmpfun(pairs_for_env_year)
  eval_year           <- compiler::cmpfun(eval_year)
}

# cooperate with futures (avoid oversubscription)
data.table::setDTthreads(2L)

# ==========================================
# 6) PARALLEL EXECUTION
# ==========================================
plan(multisession, workers = max(1L, parallel::detectCores() - 1L))

results_obs_fast <- furrr::future_map_dfr(
  sim_years,
  ~eval_year(.x, sim_redds_split, env_cache, tdm_defs),
  .options = furrr::furrr_options(seed = TRUE, scheduling = 20)
)

# results_obs_fast has columns: sim_year, env, variant, method="observed", mean_cum_surv
# ==========================================

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

# ----------------------------
# CALIBRATION (uses real deg-days)
# ----------------------------

# pick a reference env for calibration (or set explicitly, e.g. "Alt1")
ref_env <- names(env_ext_list)[1]

# degree-day vector for the reference env (length == n_calib)
deg_day_cal_ref <- deg_day_cal_for(ref_env)
stopifnot(length(deg_day_cal_ref) == length(real_years),
          all(is.finite(deg_day_cal_ref)))


# run per-variant calibration
calib_results <- furrr::future_map_dfr(
  variant_names,
  function(v) {
    opt <- optim(
      par    = c(0.0025, 0.8),
      fn     = modular_sse,
      variant= v,
      method = "L-BFGS-B",
      lower  = c(0, 0),
      upper  = c(1, 1)
    )
    tibble::tibble(
      variant  = v,
      SAR_mean = opt$par[1],
      rear_surv= opt$par[2],
      sse      = opt$value
    )
  },
  .options = furrr::furrr_options(seed = TRUE)
)

# parameter lists copied across envs (as before)
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

# precompute calibration-tab predictions (uses the SAME deg-day vector)
calib_pred_by_variant <- rlang::set_names(
  lapply(variant_names, function(v) {
    P0 <- base_P_list[[v]][[ref_env]]
    out <- simulate_variant(
      surv_vec       = surv_lookup_by_variant[[v]][1:n_calib],
      P              = P0,
      years          = n_calib,
      S_init         = S_seed_calib,
      SAR_vec        = rep(P0$SAR_mean,  n_calib),
      K_spawners_vec = rep(P0$K_spawners, n_calib),
      deg_day_adult  = deg_day_cal_ref,    # <-- keep consistent with calibration
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
