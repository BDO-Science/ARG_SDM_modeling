# Load required libraries
library(ggplot2)
library(dplyr)
library(purrr)

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
  filter(sim_year >= 2011, sim_year <= 2024)

ggplot(egg_summary_2011_2024, aes(x = sim_year, y = mean_cum_surv, color = variant)) +
  geom_line(linewidth = 1) +
  # geom_point(size = 2) +
  scale_color_viridis_d() +
  scale_x_continuous(
    breaks = seq(2011, 2024, by = 1),  # more x-axis ticks
    name = "Simulation Year"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1),      # more y-axis ticks
    name = "Mean Cumulative Egg Survival"
  ) +
  labs(color = "TDM Variant") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text    = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )


egg_summary_2025 <- egg_summary %>%
  filter(sim_year == 2025) %>%
  mutate(env = factor(as.integer(env), levels = 1:10))  # forces numeric order 1–10

ggplot(egg_summary_2025, aes(x = variant, y = mean_cum_surv, fill = env)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_viridis_d(name = "Alternative") +
  scale_y_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.1),
                     name = "Mean Cumulative Egg-to-fry Survival") +
  labs(x = "TDM Variant") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text    = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )


