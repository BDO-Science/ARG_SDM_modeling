# ============================================================================
# Build and save the spawn-timing model
# ============================================================================
# precompute.R fits a cumulative link model predicting the 10-day spawn period
# from standardised October and November water temperature, uses it, and then
# throws it away -- it is never saved. That is why nothing outside precompute.R
# can work out when fish spawn under a new temperature scenario.
#
# This script refits exactly the same model from the same inputs and saves it,
# along with everything needed to apply it: the bin definitions, the bins
# actually observed, the standardisation constants, and the site split. It
# touches nothing else in app_data, so the published results are unaffected.
#
# Why it matters: Martin egg-to-fry survival is near-vertical with respect to
# spawn date in November, so a fixed spawn-date distribution cannot reproduce
# the model's behaviour. Without this object the scenario engine disagrees with
# the pipeline (Martin correlation 0.82, mean difference -0.08).
#
# Inputs  : SalmonCountR/app_data/carcassdet_1752789274_15.csv
#           SalmonCountR/app_data/env_ext_list.rds  (historical portion only)
# Output  : SalmonCountR/app_data/spawn_timing_model.rds
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(lubridate); library(readr); library(ordinal); library(here)
})

source(here("SalmonCountR", "functions.R"))

app <- function(...) here("SalmonCountR", "app_data", ...)

carcass_raw  <- read_csv(app("carcassdet_1752789274_15.csv"), show_col_types = FALSE)
env_ext_list <- readRDS(app("env_ext_list.rds"))

# ---- 1. Carcass observations -> spawn dates, sites, brood years -------------
# Mirrors precompute.R section 9. Spawn date is the survey date minus seven days
# (decomposition time); sections NB/W/1a/1b/2 sit near Hazel, section 3 near Watt.
carcass_df <- carcass_raw %>%
  mutate(
    Date       = as.Date(surveydate),
    spawn_dt   = Date - days(7),
    brood_year = if_else(month(spawn_dt) >= 9, year(spawn_dt), year(spawn_dt) - 1),
    site = case_when(
      section %in% c("NB", "W", "1a", "1b", "1a/1b", "2") ~ "AveHazel",
      section %in% c("3")                                 ~ "AveWatt",
      TRUE                                                ~ NA_character_
    )
  ) %>%
  filter(!is.na(spawn_dt), !is.na(site))

# ---- 2. Season anchor and 10-day bins (precompute.R section 10) -------------
anchor_mmdd <- carcass_df %>%
  mutate(mmdd = format(spawn_dt, "%m-%d"), m = month(spawn_dt)) %>%
  filter(m >= 8) %>%
  summarise(anchor = min(mmdd, na.rm = TRUE)) %>%
  pull(anchor)
if (!length(anchor_mmdd) || is.na(anchor_mmdd)) anchor_mmdd <- "10-05"

bin_width   <- 10L
md          <- format(carcass_df$spawn_dt, "%m-%d")
season_year <- year(carcass_df$spawn_dt) - (md < anchor_mmdd)
season_day  <- as.integer(carcass_df$spawn_dt - as.Date(paste0(season_year, "-", anchor_mmdd)))
n_bins      <- max(12L, ceiling((max(season_day, na.rm = TRUE) + 1) / bin_width))

bin_defs <- tibble(period = paste0("p", seq_len(n_bins))) %>%
  mutate(start = as.Date(paste0("2000-", anchor_mmdd)) + (row_number() - 1) * bin_width,
         end   = start + (bin_width - 1))

# January spawning belongs to the previous year's run, so it maps onto 2001 in
# the reference season that bin_defs is built on.
carcass_df$spawn_bin <- vapply(carcass_df$spawn_dt, function(d) {
  if (is.na(d)) return(NA_character_)
  dd  <- if (month(d) >= 9) as.Date(paste0("2000-", format(d, "%m-%d")))
         else                as.Date(paste0("2001-", format(d, "%m-%d")))
  i <- which(dd >= bin_defs$start & dd <= bin_defs$end)
  if (length(i) == 1) bin_defs$period[i] else NA_character_
}, character(1))

carcass_df$spawn_bin <- factor(carcass_df$spawn_bin, levels = bin_defs$period, ordered = TRUE)
yrs <- sort(unique(as.integer(carcass_df$brood_year)))

# ---- 3. Observed October/November temperatures (precompute.R section 12) ----
# Historical temperatures are identical across alternatives, so any of the 36
# series gives the same answer; this collapses them and checks that assumption.
hist_temps <- env_ext_list %>%
  imap_dfr(~ mutate(.x, alt = as.character(.y))) %>%
  filter(year(Date) %in% yrs, month(Date) %in% c(10, 11))

# The check must compare alternatives on the SAME day -- grouping by month first
# would count variation across days and always trip.
alt_spread <- hist_temps %>%
  group_by(site, Date) %>%
  summarise(spread = max(temp) - min(temp), .groups = "drop") %>%
  summarise(max_spread = max(spread), n_differing = sum(spread > 1e-9))

if (alt_spread$n_differing > 0) {
  warning(sprintf(paste0("Historical temperatures differ across alternatives on %d ",
                         "site-dates (max %.4f degC); the CLM inputs are no longer ",
                         "alternative-independent."),
                  alt_spread$n_differing, alt_spread$max_spread))
} else {
  cat(sprintf("Historical temperatures identical across all %d alternatives (%d site-dates checked).\n",
              length(env_ext_list), nrow(hist_temps) / length(env_ext_list)))
}

oct_nov_temps <- hist_temps %>%
  group_by(site, year = year(Date), month = month(Date)) %>%
  summarise(mean_temp = mean(temp, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = month, values_from = mean_temp, names_prefix = "m_") %>%
  rename(Oct = m_10, Nov = m_11)

stopifnot(!any(duplicated(oct_nov_temps[c("site", "year")])))

carcass_df2 <- carcass_df %>%
  mutate(brood_year = as.integer(brood_year)) %>%
  left_join(oct_nov_temps, by = c("site", "brood_year" = "year"))

# ---- 4. Standardise and fit (precompute.R sections 13-14) ------------------
oct_mean <- mean(carcass_df2$Oct, na.rm = TRUE); oct_sd <- sd(carcass_df2$Oct, na.rm = TRUE)
nov_mean <- mean(carcass_df2$Nov, na.rm = TRUE); nov_sd <- sd(carcass_df2$Nov, na.rm = TRUE)

obs_fit <- carcass_df2 %>%
  mutate(Oct_std    = (Oct - oct_mean) / oct_sd,
         Nov_std    = (Nov - nov_mean) / nov_sd,
         spawn_bin  = factor(spawn_bin, levels = bin_defs$period, ordered = TRUE),
         brood_year = factor(brood_year)) %>%
  filter(!is.na(spawn_bin), !is.na(Oct_std), !is.na(Nov_std)) %>%
  droplevels()

spawn_clm <- clm(spawn_bin ~ Oct_std + Nov_std, data = obs_fit, link = "logit")
print(summary(spawn_clm))

beta  <- coef(spawn_clm)[c("Oct_std", "Nov_std")]
zeta  <- spawn_clm$Theta
if (is.matrix(zeta)) zeta <- zeta[1, ]
spawn_bins_model <- levels(spawn_clm$model$spawn_bin)
present_bins     <- levels(droplevels(carcass_df$spawn_bin))

# Redd share by site, so a scenario can be weighted the way the surveys observed
site_props <- carcass_df %>% count(site) %>% mutate(prop = n / sum(n))

# ---- 5. Save ----------------------------------------------------------------
spawn_timing_model <- list(
  beta             = beta,
  zeta             = zeta,
  bin_defs         = bin_defs,
  bin_width        = bin_width,
  anchor_mmdd      = anchor_mmdd,
  spawn_bins_model = spawn_bins_model,
  present_bins     = present_bins,
  standardisation  = list(oct_mean = oct_mean, oct_sd = oct_sd,
                          nov_mean = nov_mean, nov_sd = nov_sd),
  site_props       = site_props,
  n_obs            = nrow(obs_fit),
  brood_years      = yrs,
  fitted           = Sys.time(),
  source           = "analysis/build_spawn_timing_model.R; refit of the CLM in precompute.R section 14"
)

saveRDS(spawn_timing_model, app("spawn_timing_model.rds"))

cat(sprintf("\nFitted on %d carcass observations, brood years %d-%d\n",
            nrow(obs_fit), min(yrs), max(yrs)))
cat(sprintf("Oct_std %+0.4f   Nov_std %+0.4f\n", beta[["Oct_std"]], beta[["Nov_std"]]))
cat(sprintf("%d bins in the model, %d observed: %s\n",
            length(spawn_bins_model), length(present_bins),
            paste(present_bins, collapse = ", ")))
cat(sprintf("Standardisation: Oct %.2f +/- %.2f, Nov %.2f +/- %.2f degC\n",
            oct_mean, oct_sd, nov_mean, nov_sd))
cat(sprintf("Site split: %s\n",
            paste(sprintf("%s %.1f%%", site_props$site, 100 * site_props$prop),
                  collapse = ", ")))
cat("\nWrote SalmonCountR/app_data/spawn_timing_model.rds\n")
