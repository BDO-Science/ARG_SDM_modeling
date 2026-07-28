# ============================================================================
# Calibration fit statistics, regenerated
# ============================================================================
# Restores the calibration-prediction step that precompute.R section 33 stopped
# producing. That section now reads
#
#     # Since we're not calibrating, create empty calib_pred_by_variant for
#     # compatibility
#     calib_pred_by_variant <- list()
#
# and writes a 0-byte calib_pred_by_variant.rds, so the fit statistics quoted in
# SI Section S2.8 (R2 = 0.72, RMSE = 8,240 spawners, MAPE = 24%) could not be
# reproduced from the repository. This script reproduces them from saved
# artifacts alone, without re-running precompute.R and therefore without
# disturbing any published number.
#
# Method, matching functions.R::combined_sse exactly:
#   * survival input  : egg_summary for the reference alternative (env 1) over
#                       the calibration years, one vector per TDM variant
#   * parameters      : calib_results (calibrated SAR and rear_surv) on top of
#                       base_P
#   * seed            : S_seed_calib, the observed 2011-2013 escapement
#   * degree-days     : deg_day_cal_for(ref_env)
#   * simulation      : simulate_variant over real_years
#   * fit window      : fit_idx, i.e. 2014-2024, the 11 years not used as seed
#
# "TDM-weighted predictions" in S2.8 is taken to mean the elicited panel weights
# (0.51 Bratovich / 0.24 Bartholow / 0.25 Martin) applied to the three variant
# predictions. Per-variant statistics are reported as well so the weighting
# choice is transparent.
#
# Inputs  : SalmonCountR/app_data/{egg_summary,calib_results,base_P,
#           S_seed_calib,spawn_dates_by_alt,sim_years,env_ext_list,
#           american_river_instream,grandtab_*}.
# Outputs : output/calibration_fit_statistics.csv
#           output/calibration_predictions.csv
#           figures/calibration_observed_vs_predicted.png
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr)
  library(ggplot2); library(here)
})

source(here("SalmonCountR", "functions.R"))
source(here("analysis", "figure_theme.R"))

app <- function(...) here("SalmonCountR", "app_data", ...)

egg_summary        <- readRDS(app("egg_summary.rds"))
calib_results      <- readRDS(app("calib_results.rds"))
base_P             <- readRDS(app("base_P.rds"))
S_seed_calib       <- readRDS(app("S_seed_calib.rds"))
spawn_dates_by_alt <- readRDS(app("spawn_dates_by_alt.rds"))
sim_years          <- readRDS(app("sim_years.rds"))
env_ext_list       <- readRDS(app("env_ext_list.rds"))

real_years <- 2011:2024
ref_env    <- names(env_ext_list)[1]
variants   <- c("exp_WF", "exp_SM", "lin_Martin")
tdm_w      <- c(exp_WF = 0.51, exp_SM = 0.24, lin_Martin = 0.25)

# ---- Observed escapement ----------------------------------------------------
esc_obs <- read_csv(app("grandtab_1752793045_337.csv"),
                    col_types = cols(`End Year of Monitoring Period` = col_character(),
                                     `Population Estimate` = col_double())) %>%
  mutate(year = parse_number(`End Year of Monitoring Period`)) %>%
  filter(between(year, min(real_years), max(real_years))) %>%
  transmute(year, observed = `Population Estimate`) %>%
  arrange(year)

stopifnot(nrow(esc_obs) == length(real_years))
obs_spawners <- esc_obs$observed

n_calib <- length(obs_spawners)
fit_idx <- (length(S_seed_calib) + 1):n_calib      # years 4-14 => 2014-2024

# ---- Survival vectors for the reference alternative -------------------------
surv_lookup_by_variant <- egg_summary %>%
  filter(env == ref_env, sim_year %in% real_years, variant %in% variants) %>%
  arrange(variant, sim_year) %>%
  group_by(variant) %>%
  summarise(surv_vec = list(mean_cum_surv), .groups = "drop") %>%
  tibble::deframe()

stopifnot(all(lengths(surv_lookup_by_variant) == n_calib))

# ---- Pre-spawn degree-days for the reference alternative --------------------
deg_day_cal_ref <- deg_day_cal_for(ref_env)
stopifnot(length(deg_day_cal_ref) == n_calib)

# ---- Reproduce the calibration predictions ---------------------------------
best_SAR       <- calib_results$SAR_mean[1]
best_rear_surv <- calib_results$rear_surv[1]

cat(sprintf("Calibrated parameters in calib_results.rds: SAR = %.6f, rear_surv = %.6f\n",
            best_SAR, best_rear_surv))

calib_pred_by_variant <- map(setNames(variants, variants), function(v) {
  P <- base_P
  P$SAR_mean  <- best_SAR
  P$rear_surv <- best_rear_surv
  simulate_variant(
    surv_vec       = surv_lookup_by_variant[[v]],
    P              = P,
    years          = n_calib,
    S_init         = S_seed_calib,
    SAR_vec        = rep(best_SAR, n_calib),
    K_spawners_vec = rep(P$K_spawners, n_calib),
    deg_day_adult  = deg_day_cal_ref,
    sim_years_vec  = real_years
  )
})

pred <- imap_dfr(calib_pred_by_variant,
                 ~ tibble(year = real_years, variant = .y, predicted = .x$spawners)) %>%
  left_join(esc_obs, by = "year")

# TDM-weighted prediction
pred_w <- pred %>%
  mutate(w = tdm_w[variant]) %>%
  group_by(year, observed) %>%
  summarise(predicted = sum(predicted * w) / sum(w), .groups = "drop") %>%
  mutate(variant = "TDM-weighted") %>%
  select(year, variant, predicted, observed)

all_pred <- bind_rows(pred, pred_w) %>%
  mutate(in_fit_window = year %in% real_years[fit_idx]) %>%
  arrange(variant, year)

# ---- Fit statistics ---------------------------------------------------------
# R2 is reported two ways because they are not interchangeable and the SI does
# not say which was used:
#   r2_cor   = squared Pearson correlation, the "variance explained" reading
#   r2_nash  = 1 - SSE/SST (Nash-Sutcliffe), the goodness-of-fit reading, which
#              can be negative and is the stricter of the two
fit_stats <- function(obs, prd) {
  ok  <- is.finite(obs) & is.finite(prd)
  obs <- obs[ok]; prd <- prd[ok]
  resid <- prd - obs
  tibble(
    n        = length(obs),
    r2_cor   = suppressWarnings(cor(obs, prd))^2,
    r2_nash  = 1 - sum(resid^2) / sum((obs - mean(obs))^2),
    rmse     = sqrt(mean(resid^2)),
    mae      = mean(abs(resid)),
    mape_pct = 100 * mean(abs(resid / obs)),
    bias     = mean(resid)
  )
}

stats_fit <- all_pred %>%
  filter(in_fit_window) %>%
  group_by(variant) %>%
  group_modify(~ fit_stats(.x$observed, .x$predicted)) %>%
  ungroup() %>%
  mutate(window = "2014-2024 (fit years)")

stats_all <- all_pred %>%
  group_by(variant) %>%
  group_modify(~ fit_stats(.x$observed, .x$predicted)) %>%
  ungroup() %>%
  mutate(window = "2011-2024 (all years, incl. seed)")

stats <- bind_rows(stats_fit, stats_all) %>%
  select(window, variant, n, r2_cor, r2_nash, rmse, mae, mape_pct, bias)

cat("\n=== Calibration fit statistics ===\n")
print(as.data.frame(stats %>% mutate(across(where(is.numeric), ~round(., 4)))),
      row.names = FALSE)

published <- stats %>%
  filter(window == "2014-2024 (fit years)", variant == "TDM-weighted")

cat("\n=== Against the values quoted in SI Section S2.8 ===\n")
cat(sprintf("  SI quotes    : R2 = 0.72, RMSE = 8,240 spawners, MAPE = 24%%\n"))
cat(sprintf("  Regenerated  : R2 = %.2f (squared correlation) / %.2f (Nash-Sutcliffe),\n",
            published$r2_cor, published$r2_nash))
cat(sprintf("                 RMSE = %s spawners, MAPE = %.0f%%\n",
            format(round(published$rmse), big.mark = ","), published$mape_pct))

# ---- Outputs ----------------------------------------------------------------
dir.create(here("output"),  showWarnings = FALSE)
dir.create(here("figures"), showWarnings = FALSE)

write_csv(stats,    here("output", "calibration_fit_statistics.csv"))
write_csv(all_pred, here("output", "calibration_predictions.csv"))
saveRDS(calib_pred_by_variant, app("calib_pred_by_variant.rds"))

plot_df <- all_pred %>%
  mutate(variant = factor(variant,
                          levels = c("TDM-weighted", "exp_WF", "exp_SM", "lin_Martin"),
                          labels = c("TDM-weighted", "Bratovich (2020)",
                                     "Bartholow & Heasley (2006)", "Martin et al. (2017)")))

series_cols <- setNames(PAIR_COLS, c("Observed escapement", "Model prediction"))

plot_long <- plot_df %>%
  select(year, variant, Observed = observed, Predicted = predicted) %>%
  pivot_longer(c(Observed, Predicted), names_to = "series", values_to = "spawners") %>%
  mutate(series = factor(series, levels = c("Observed", "Predicted"),
                         labels = names(series_cols)))

p <- ggplot(plot_long, aes(year, spawners, colour = series)) +
  annotate("rect", xmin = min(real_years) - 0.5,
           xmax = real_years[length(S_seed_calib)] + 0.5,
           ymin = -Inf, ymax = Inf, fill = "grey88") +
  annotate("text", x = mean(real_years[1:length(S_seed_calib)]),
           y = Inf, label = "seed years", vjust = 1.5, size = 3.6,
           fontface = "bold", colour = "black") +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.1) +
  facet_wrap(~ variant, ncol = 2) +
  scale_colour_manual(values = series_cols, name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(breaks = seq(2011, 2024, 3)) +
  labs(subtitle = "Calibration fit, 2011-2024: observed escapement vs predicted spawners",
       x = NULL, y = "Spawners") +
  theme_arg(base_size = 14, legend = "top")

ggsave(here("figures", "calibration_observed_vs_predicted.png"), p,
       width = 12, height = 8.5, dpi = 300, bg = "white")

cat("\nWrote output/calibration_fit_statistics.csv,",
    "output/calibration_predictions.csv,",
    "figures/calibration_observed_vs_predicted.png,",
    "and repopulated app_data/calib_pred_by_variant.rds\n")
