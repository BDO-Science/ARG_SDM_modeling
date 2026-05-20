# ============================================================
# Presentation Plots — fully corrected
# Updated 2026-05-07
#
# Fixes from prior version:
#   - ATU constants updated: egg_ATU=400, alev_ATU=558, total_ATU=958
#     (was 958/417/1375 — those are the legacy SacPAS values)
#   - egg_summary plotting uses correct columns: sim_year (not year),
#     mean_cum_surv (not days_lt_18.3C), env (not alternative)
#   - Days-below-18.3°C (steelhead metric) plotted from its proper
#     source: steelhead_metrics.rds, not egg_summary
# ============================================================

library(here)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(scales)
library(tibble)

# ============================================================
# 1) DENSITY DEPENDENCE CURVES (flow-dependent K)
# ============================================================

instream <- readRDS(here("SalmonCountR", "app_data", "american_river_instream.rds")) %>%
  mutate(K_spawners = FR_spawn_wua / 9.29)  # 9.29 m^2 per redd (Burner 1951)

get_K_spawners <- function(flow_vec, lookup = instream) {
  approx(
    x    = lookup$flow_cfs,
    y    = lookup$K_spawners,
    xout = flow_vec,
    rule = 2
  )$y
}

flows <- seq(1000, 5000, length.out = 5)
redds <- seq(0, 100000, by = 100)
S0    <- 0.347

dd_plot_df <- map_dfr(flows, function(flow) {
  K <- get_K_spawners(flow)
  surv <- S0 / (1 + redds / K)
  data.frame(redds = redds, survival = surv, flow_cfs = factor(flow))
})

p_dd <- ggplot(dd_plot_df, aes(x = redds, y = survival, color = flow_cfs)) +
  geom_line(linewidth = 1.2) +
  labs(
    x = "Redds (Female Spawner Abundance)",
    y = "Density-Dependent Survival Rate",
    color = "Flow (cfs)"
  ) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.85) +
  scale_y_continuous(limits = c(0, 0.35), expand = c(0, 0)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title   = element_text(face = "bold", size = 14),
    axis.text    = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 13),
    legend.text  = element_text(size = 11)
  )

print(p_dd)

# ============================================================
# 2) EGG-TO-FRY (ETF) CUMULATIVE SURVIVAL — TIME SERIES
# ============================================================
# egg_summary columns: env, variant, sim_year, mean_cum_surv

egg_summary <- readRDS(here("SalmonCountR", "app_data", "egg_summary.rds"))

egg_summary_2011_2024 <- egg_summary %>%
  filter(env == "1",
         sim_year >= 2011, sim_year <= 2024,
         variant != "weighted_avg")

p_etf_ts <- ggplot(egg_summary_2011_2024,
                   aes(x = sim_year, y = mean_cum_surv, color = variant)) +
  geom_line(linewidth = 1) +
  scale_color_viridis_d() +
  scale_x_continuous(
    breaks = seq(2011, 2024, by = 1),
    name = "Simulation Year"
  ) +
  scale_y_continuous(
    name = "Mean Cumulative Egg-to-Fry Survival",
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(color = "TDM Variant") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title       = element_text(face = "bold"),
    axis.text        = element_text(face = "bold"),
    legend.title     = element_text(face = "bold"),
    legend.text      = element_text(size = 11),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  )

print(p_etf_ts)

# ============================================================
# 3) ETF SURVIVAL BY ALTERNATIVE x TDM VARIANT (forecast year 2025)
# ============================================================

alt_names <- c("NB", "PB1", "PB2", "PB2b", "PB2c", "PB3", "PB4", "PB5", "PB6")

egg_summary_2025 <- egg_summary %>%
  filter(sim_year == 2025,
         variant != "weighted_avg",
         env %in% as.character(1:9)) %>%
  mutate(alternative = factor(alt_names[as.integer(env)], levels = alt_names))

p_etf_alt <- ggplot(egg_summary_2025,
                    aes(x = variant, y = mean_cum_surv, fill = alternative)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_viridis_d(name = "Alternative") +
  scale_y_continuous(
    name = "Mean Cumulative Egg-to-Fry Survival",
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(x = "TDM Variant") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title       = element_text(face = "bold"),
    axis.text        = element_text(face = "bold"),
    legend.title     = element_text(face = "bold"),
    legend.text      = element_text(size = 11),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  )

print(p_etf_alt)

# ============================================================
# 4) STEELHEAD METRIC: DAYS BELOW 18.3 deg C
# ============================================================
# This is the "days_lt_18.3C" data the original code was reaching for.
# It lives in steelhead_metrics.rds, with column `steelhead_score`.

steelhead_metrics <- readRDS(here("SalmonCountR", "app_data", "steelhead_metrics.rds"))

sh_by_alt <- steelhead_metrics %>%
  mutate(
    env_num     = as.integer(env),
    alt_idx     = ((env_num - 1) %% 9) + 1,
    alternative = factor(alt_names[alt_idx], levels = alt_names),
    climate     = c("2011","2014","2017","2020")[((env_num - 1) %/% 9) + 1]
  )

p_sh <- ggplot(sh_by_alt, aes(x = alternative, y = steelhead_score, fill = climate)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_viridis_d(name = "Climate Year") +
  scale_y_continuous(name = expression("Days Below 18.3" * degree * "C (Oct-Nov)")) +
  labs(x = "Management Alternative") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title       = element_text(face = "bold"),
    axis.text        = element_text(face = "bold"),
    legend.title     = element_text(face = "bold"),
    legend.text      = element_text(size = 11),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  )

print(p_sh)

# ============================================================
# 5) DAILY TDM SURVIVAL CURVES (illustrative, no data dependency)
# ============================================================

s_day_exp    <- function(T, alpha, beta) exp(-alpha * exp(beta * T))
s_day_martin <- function(T, alpha = 0.026, beta = 12.14) exp(-alpha * pmax(T - beta, 0))

SM_egg    <- list(alpha = 1.475e-11,    beta = 1.392)
SM_alevin <- list(alpha = 2.521e-12,    beta = 1.461)
WF_egg    <- list(alpha = 3.408488e-11, beta = 1.21122)
WF_alevin <- list(alpha = 1.017554e-10, beta = 1.24092)
Martin    <- list(alpha = 0.026,        beta = 12.14)

T_seq <- seq(6, 22, by = 0.1)

curves_daily <- bind_rows(
  tibble(T = T_seq, family = "Bratovich et al. (2020)",   stage = "Egg",
         S_day = s_day_exp(T_seq, WF_egg$alpha,    WF_egg$beta)),
  tibble(T = T_seq, family = "Bratovich et al. (2020)",   stage = "Alevin",
         S_day = s_day_exp(T_seq, WF_alevin$alpha, WF_alevin$beta)),
  tibble(T = T_seq, family = "Bartholow & Heasley (2006)", stage = "Egg",
         S_day = s_day_exp(T_seq, SM_egg$alpha,    SM_egg$beta)),
  tibble(T = T_seq, family = "Bartholow & Heasley (2006)", stage = "Alevin",
         S_day = s_day_exp(T_seq, SM_alevin$alpha, SM_alevin$beta)),
  tibble(T = T_seq, family = "Martin et al. (2017)",      stage = "Incubation",
         S_day = s_day_martin(T_seq, Martin$alpha, Martin$beta))
) %>%
  mutate(S_day = pmin(pmax(S_day, 0), 1))

stage_types <- c("Egg" = "solid", "Alevin" = "dashed", "Incubation" = "dotted")

p_daily <- ggplot(curves_daily, aes(T, S_day, color = family, linetype = stage)) +
  geom_line(linewidth = 1.3) +
  scale_color_viridis_d(name = "TDM Model") +
  scale_linetype_manual(values = stage_types, name = "Stage") +
  coord_cartesian(xlim = c(6, 22), ylim = c(0, 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Temperature (deg C)", y = "Daily survival rate") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title       = element_text(face = "bold"),
    axis.text        = element_text(face = "bold", size = 11),
    legend.position  = "right",
    legend.title     = element_text(face = "bold"),
    legend.text      = element_text(size = 11),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  )

print(p_daily)
ggsave("tdm_daily_survival.png", p_daily, width = 11, height = 7, dpi = 300, bg = "white")

# ============================================================
# 6) CUMULATIVE TDM SURVIVAL OVER TIME (illustrative scenarios)
# ============================================================

# ---- ATU stage boundaries (UPDATED to SacPAS Fish Model v3 / Zeug et al. 2012) ----
egg_ATU   <- 400    # was 958
alev_ATU  <- 558    # was 417
total_ATU <- egg_ATU + alev_ATU   # 958 (was 1375)

cum_surv_series_exp <- function(T, egg, alevin) {
  n <- length(T)
  if (n == 0L) return(tibble(day = integer(), H = numeric(), S = numeric()))
  
  atu <- cumsum(pmax(T, 0))
  i_h <- which(atu >= egg_ATU)[1]
  if (is.na(i_h)) i_h <- n
  i_e <- which(atu >= total_ATU)[1]
  if (is.na(i_e)) i_e <- n
  i_h <- max(0L, min(i_h, n))
  i_e <- max(i_h, min(i_e, n))
  
  h <- numeric(n)
  if (i_h >= 1L) {
    h[1:i_h] <- egg$alpha * exp(egg$beta * T[1:i_h])
  }
  if (i_e > i_h) {
    idx <- (i_h + 1L):i_e
    h[idx] <- alevin$alpha * exp(alevin$beta * T[idx])
  }
  
  H <- cumsum(h)
  tibble(day = seq_len(n), H = H, S = exp(-H))
}

cum_surv_series_martin <- function(T, pars) {
  n <- length(T)
  if (n == 0L) return(tibble(day = integer(), H = numeric(), S = numeric()))
  h <- pars$alpha * pmax(T - pars$beta, 0)
  h[!is.finite(h)] <- 0
  H <- cumsum(h)
  tibble(day = seq_len(n), H = H, S = exp(-H))
}

n_days <- 60
T_cool  <- rep(10.5, n_days)
T_ramp  <- seq(10.5, 14.5, length.out = n_days)
T_pulse <- {x <- rep(11.5, n_days); x[26:32] <- 16.5; x}

scens <- list(
  "Cool steady (10.5 deg C)"            = T_cool,
  "Warming ramp (10.5 -> 14.5 deg C)"   = T_ramp,
  "Heatwave pulse (11.5 -> 16.5 deg C)" = T_pulse
)

cum_panel <- purrr::imap_dfr(scens, function(Ts, nm) {
  bind_rows(
    cum_surv_series_exp(Ts, WF_egg, WF_alevin) %>% mutate(model = "Bratovich et al. (2020)"),
    cum_surv_series_exp(Ts, SM_egg, SM_alevin) %>% mutate(model = "Bartholow & Heasley (2006)"),
    cum_surv_series_martin(Ts, Martin)         %>% mutate(model = "Martin et al. (2017)")
  ) %>% mutate(scenario = nm)
})

p_cum_time <- ggplot(cum_panel, aes(day, S, color = model)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ scenario, ncol = 1) +
  scale_color_viridis_d(name = "TDM Model") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Incubation day", y = "Cumulative survival") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title       = element_text(face = "bold"),
    axis.text        = element_text(face = "bold", size = 11),
    strip.text       = element_text(face = "bold", size = 11),
    legend.title     = element_text(face = "bold"),
    legend.text      = element_text(size = 11),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  )

print(p_cum_time)
ggsave("tdm_cumulative_over_time.png", p_cum_time, width = 12, height = 7, dpi = 300, bg = "white")