library(tidyverse)
library(here)

# Load data
source(here("SalmonCountR", "functions.R"))
env_ext_list <- readRDS(here("SalmonCountR", "app_data", "env_ext_list.rds"))
surv_lookup_full <- readRDS(here("SalmonCountR", "app_data", "surv_lookup_full.rds"))
base_P_list <- readRDS(here("SalmonCountR", "app_data", "base_P_list.rds"))
spawn_dates_by_alt <- readRDS(here("SalmonCountR", "app_data", "spawn_dates_by_alt.rds"))
sim_years <- readRDS(here("SalmonCountR", "app_data", "sim_years.rds"))

# Load observed data
esc_obs <- read_csv(
  here("SalmonCountR", "app_data", "grandtab_1752793045_337.csv"),
  col_types = cols(
    `End Year of Monitoring Period` = col_character(),
    `Population Estimate` = col_double()
  )
) %>%
  mutate(year = parse_number(`End Year of Monitoring Period`)) %>%
  filter(between(year, 2011, 2024)) %>%
  rename(spawners = `Population Estimate`)

# Define historical period
real_years <- 2011:2024
n_calib <- length(real_years)
S_seed_calib <- esc_obs$spawners[1:3]

# Reference environment
ref_env <- names(env_ext_list)[1]

# Calculate degree-days for historical period
deg_day_cal <- compute_deg_day_adult(
  env_nm = ref_env,
  sim_years = real_years,
  spawn_dates = spawn_dates_by_alt[[ref_env]][match(real_years, sim_years)],
  env_ext_list = env_ext_list
)

# Run simulation for each TDM variant
variant_names <- c("exp_WF", "exp_SM", "lin_Martin")

results_by_variant <- map_dfr(variant_names, function(v) {
  P0 <- base_P_list[[v]][[ref_env]]
  surv_vec <- surv_lookup_full[[paste(ref_env, v, sep = "_")]][1:n_calib]
  
  out <- simulate_variant(
    surv_vec = surv_vec,
    P = P0,
    years = n_calib,
    S_init = S_seed_calib,
    SAR_vec = rep(P0$SAR_mean, n_calib),
    K_spawners_vec = rep(P0$K_spawners, n_calib),
    deg_day_adult = deg_day_cal,
    sim_years_vec = real_years
  )
  
  tibble(
    year = real_years,
    variant = v,
    predicted = out$spawners
  )
})

# Add observed data
plot_data <- results_by_variant %>%
  left_join(esc_obs %>% dplyr::select(year, observed = spawners), by = "year")

# Calculate fit statistics
fit_stats <- plot_data %>%
  filter(year >= 2014) %>%  # Exclude seed years
  group_by(variant) %>%
  summarise(
    RMSE = sqrt(mean((predicted - observed)^2)),
    MAE = mean(abs(predicted - observed)),
    R2 = cor(predicted, observed)^2,
    .groups = "drop"
  )

print("Fit Statistics (2014-2024):")
print(fit_stats)

# Plot 1: Time series comparison
p1 <- ggplot(plot_data, aes(x = year)) +
  geom_line(aes(y = observed), color = "black", size = 1.5, alpha = 0.8) +
  geom_point(aes(y = observed), color = "black", size = 3) +
  geom_line(aes(y = predicted, color = variant), size = 1) +
  geom_point(aes(y = predicted, color = variant), size = 2) +
  geom_vline(xintercept = 2013.5, linetype = "dashed", alpha = 0.5) +
  annotate("text", x = 2012, y = max(plot_data$observed) * 0.95, 
           label = "Seed Years", hjust = 0.5, size = 3.5) +
  scale_color_viridis_d(name = "TDM Variant",
                        labels = c("Water Forum 2020", "SALMOD 2006", "Martin 2017")) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Fixed Parameters: Predicted vs Observed Spawners (2011-2024)",
    subtitle = "SAR = 0.0025, Rearing Survival = 0.5419",
    x = "Year",
    y = "Spawner Abundance",
    caption = "Black line = Observed | Colored lines = Predicted with fixed parameters"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(p1)
ggsave(here("output", "fixed_params_vs_observed.png"), p1, width = 12, height = 7, dpi = 300)

# Plot 2: Predicted vs Observed scatter
p2 <- plot_data %>%
  filter(year >= 2014) %>%
  ggplot(aes(x = observed, y = predicted, color = variant)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
  facet_wrap(~ variant, labeller = labeller(
    variant = c("exp_WF" = "Water Forum 2020", 
                "exp_SM" = "SALMOD 2006", 
                "lin_Martin" = "Martin 2017")
  )) +
  scale_color_viridis_d(guide = "none") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Predicted vs Observed Spawners (2014-2024)",
    subtitle = "Dashed line = perfect fit",
    x = "Observed Spawners",
    y = "Predicted Spawners"
  ) +
  theme_minimal(base_size = 14)

print(p2)
ggsave(here("output", "fixed_params_scatter.png"), p2, width = 12, height = 5, dpi = 300)

# Plot 3: Residuals over time
p3 <- plot_data %>%
  filter(year >= 2014) %>%
  mutate(residual = predicted - observed) %>%
  ggplot(aes(x = year, y = residual, color = variant)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_viridis_d(name = "TDM Variant",
                        labels = c("Water Forum 2020", "SALMOD 2006", "Martin 2017")) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Model Residuals Over Time (2014-2024)",
    subtitle = "Positive = Over-prediction, Negative = Under-prediction",
    x = "Year",
    y = "Residual (Predicted - Observed)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(p3)
ggsave(here("output", "fixed_params_residuals.png"), p3, width = 12, height = 6, dpi = 300)

cat("\nPlots saved to output/ directory\n")