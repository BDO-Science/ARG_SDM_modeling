# Load required libraries
library(ggplot2)
library(dplyr)
library(purrr)

sim_year <- 114

# Load instream flow vs K_spawners data
load("SalmonCountR/app_data/american_river_instream.rda")  # should load object `instream`

# 2) Compute redd capacity for each flow (if not already done)
instream <- instream %>%
  mutate(K_spawners = FR_spawn_wua / 9.29)  # 9.29 m² per redd

# 3) Function to get K_spawners at any flow
get_K_spawners <- function(flow_vec, lookup = instream) {
  approx(
    x    = lookup$flow_cfs,
    y    = lookup$K_spawners,
    xout = flow_vec,
    rule = 2
  )$y
}

# Define flows and female spawner values
flows <- seq(1000, 5000, length.out = 5)
redds <- seq(0, 100000, by = 100)
S0    <- 0.347  # max survival

# Create plotting data
dd_plot_df <- map_dfr(flows, function(flow) {
  K <- get_K_spawners(flow)
  surv <- S0 / (1 + redds / K)
  data.frame(redds = redds, survival = surv, flow_cfs = factor(flow))
})
# 5) Plot
ggplot(dd_plot_df, aes(x = redds, y = survival, color = flow_cfs)) +
  geom_line(size = 1.2) +
  labs(
    title    = NULL,
    subtitle = NULL,
    x = "Redds (Female Spawner Abundance)",
    y = "Density-Dependent Survival Rate",
    color = "Flow (cfs)"
  ) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.85) +
  scale_y_continuous(limits = c(0, 0.35), expand = c(0, 0)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title  = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title  = element_text(face = "bold", size = 14),
    axis.text   = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 13),
    legend.text  = element_text(size = 11)
  )

#####################
#ETF survival plotting
#######################
# Filter years 2011–2024
egg_summary_2011_2024 <- egg_summary %>%
  filter(year >= 2011, year <= 2024)

ggplot(egg_summary_2011_2024, aes(x = year, y = days_lt_18.3C, color = variant)) +
  geom_line(linewidth = 1) +
  # geom_point(size = 2) +
  scale_color_viridis_d() +
  scale_x_continuous(
    breaks = seq(2011, 2024, by = 1),  # more x-axis ticks
    name = "Simulation Year"
  ) +
  scale_y_continuous(
    name = expression("Days Below 18.3" * degree * "C")
  ) +
  labs(color = "TDM Variant") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title  = element_text(face = "bold"),
    axis.text   = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(size = 11),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.border     = element_rect(colour = "black", fill = NA, size = 0.5)
  )


egg_summary_2025 <- egg_summary %>%
  filter(year == 2025) %>%
  mutate(alternative = factor(as.integer(alternative), levels = sort(unique(as.integer(alternative)))))

ggplot(egg_summary_2025, aes(x = variant, y = days_lt_18.3C, fill = alternative)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_viridis_d(name = "Alternative") +
  scale_y_continuous(
    name = expression("Days Below 18.3" * degree * "C")
  ) +
  labs(x = "TDM Variant") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title  = element_text(face = "bold"),
    axis.text   = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(size = 11),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.border     = element_rect(colour = "black", fill = NA, size = 0.5)
  )


# ---- Libraries ----
library(tidyverse)
library(scales)

# ---- Daily survival (hazard → per-day survival) ----
# Exponential: h(T) = α * exp(β T);  S_day = exp(-h)
s_day_exp    <- function(T, alpha, beta) exp(-alpha * exp(beta * T))
# Martin 2017: h(T) = α * max(T − β, 0); S_day = exp(-h)
s_day_martin <- function(T, alpha = 0.026, beta = 12.14) exp(-alpha * pmax(T - beta, 0))

# ---- Parameter sets (UPDATED) ----
# SALMOD 2006 / USBR 2008 / HCI 1996
SM_egg    <- list(alpha = 1.475e-11,  beta = 1.392)
SM_alevin <- list(alpha = 2.521e-12,  beta = 1.461)

# Water Forum 2020
WF_egg    <- list(alpha = 3.408488e-11, beta = 1.21122)
WF_alevin <- list(alpha = 1.017554e-10,  beta = 1.24092)

# Martin 2017 (linear-threshold — single curve for incubation)
Martin <- list(alpha = 0.026, beta = 12.14)

# ---- Temperature grid ----
T_seq <- seq(6, 22, by = 0.1)

# ---- Build tidy curves (color = family, linetype = stage) ----
curves_daily <- bind_rows(
  tibble(T = T_seq, family = "Bratovich et al. (2020)", stage = "Egg",
         S_day = s_day_exp(T_seq, WF_egg$alpha,    WF_egg$beta)),
  tibble(T = T_seq, family = "Bratovich et al. (2020)", stage = "Alevin",
         S_day = s_day_exp(T_seq, WF_alevin$alpha, WF_alevin$beta)),
  tibble(T = T_seq, family = "Bartholow & Heasley (2006)", stage = "Egg",
         S_day = s_day_exp(T_seq, SM_egg$alpha,    SM_egg$beta)),
  tibble(T = T_seq, family = "Bartholow & Heasley (2006)", stage = "Alevin",
         S_day = s_day_exp(T_seq, SM_alevin$alpha, SM_alevin$beta)),
  tibble(T = T_seq, family = "Martin et al. (2017)", stage = "Incubation",
         S_day = s_day_martin(T_seq, Martin$alpha, Martin$beta))
) %>%
  mutate(S_day = pmin(pmax(S_day, 0), 1))

# ---- Scales & styling ----
stage_types <- c("Egg" = "solid", "Alevin" = "dashed", "Incubation" = "dotted")

p <- ggplot(curves_daily, aes(T, S_day, color = family, linetype = stage)) +
  geom_line(linewidth = 1.3) +
  scale_color_viridis_d(name = "TDM Model") +
  scale_linetype_manual(values = stage_types, name = "Stage") +
  coord_cartesian(xlim = c(6, 22), ylim = c(0, 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Temperature (°C)",
    y = "Daily survival rate"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 12, margin = margin(b = 10)),
    axis.title    = element_text(face = "bold"),
    axis.text     = element_text(face = "bold", size = 11),
    legend.position = "right",
    legend.title    = element_text(face = "bold"),
    legend.text     = element_text(size = 11),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.border     = element_rect(colour = "black", fill = NA, size = 0.5)
  )

print(p)

# Slide-ready export (prevents black background in PPT)
ggsave("tdm_daily_survival.png", p, width = 11, height = 7, dpi = 300, bg = "white")

# ---- ATU stage boundaries ----
egg_ATU   <- 958
alev_ATU  <- 417
total_ATU <- egg_ATU + alev_ATU

# ---- Robust cumulative survival helpers (fixed length) -----------------------

cum_surv_series_exp <- function(T, egg, alevin) {
  n  <- length(T)
  if (n == 0L) return(tibble(day = integer(), H = numeric(), S = numeric()))
  
  # ATU thresholds -> stage indices (clamped to [0, n])
  atu <- cumsum(pmax(T, 0))
  i_h <- which(atu >= 958)[1];         if (is.na(i_h)) i_h <- n
  i_e <- which(atu >= 958 + 417)[1];   if (is.na(i_e)) i_e <- n
  i_h <- max(0L, min(i_h, n))
  i_e <- max(i_h, min(i_e, n))
  
  # build hazard vector of length n
  h <- numeric(n)
  if (i_h >= 1L) {
    h[1:i_h] <- egg$alpha    * exp(egg$beta    * T[1:i_h])
  }
  if (i_e > i_h) {
    idx <- (i_h + 1L):i_e
    h[idx] <- alevin$alpha * exp(alevin$beta * T[idx])
  }
  # (i_e+1 .. n) remain zero (no hazard after emergence)
  
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


# ---- Mock scenarios (same as before) ----
n_days <- 60
T_cool  <- rep(10.5, n_days)
T_ramp  <- seq(10.5, 14.5, length.out = n_days)
T_pulse <- {x <- rep(11.5, n_days); x[26:32] <- 16.5; x}

scens <- list(
  "Cool steady (10.5°C)"         = T_cool,
  "Warming ramp (10.5→14.5°C)"   = T_ramp,
  "Heatwave pulse (11.5→16.5°C)" = T_pulse
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
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Incubation day", y = "Cumulative survival"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    axis.title    = element_text(face = "bold"),
    axis.text     = element_text(face = "bold", size = 11),
    strip.text    = element_text(face = "bold", size = 11),
    legend.title  = element_text(face = "bold"),
    legend.text   = element_text(size = 11),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.border     = element_rect(colour = "black", fill = NA, size = 0.5)
  )

print(p_cum_time)
ggsave("tdm_cumulative_over_time.png", p_cum_time, width = 12, height = 7, dpi = 300, bg = "white")