# Load necessary libraries
library(tidyverse)
library(scales)
library(here)
library(ggh4x)

# --- 1. Load Pre-computed Model Results ---
# This file contains the full 100-year forecasts for all 36 environment combinations
# and all 3 TDM model variants.
results_full <- readRDS(here("SalmonCountR", "app_data", "results_full.rds"))
if (!exists("results_full") || !is.data.frame(results_full)) {
  stop("Error: 'results_full.rds' not found or is not a dataframe. Please ensure the file is in the 'app_data' directory.")
}

# --- 2. Helper Function to Get Alternative Mapping ---
get_scenario_alternatives <- function(scenario, hydro_year) {
  scenario_map <- c("NB"=1, "PB1"=2, "PB2"=3, "PB2b"=4, "PB2c"=5, "PB3"=6, "PB4"=7, "PB5"=8, "PB6"=9)
  hydro_map <- c("2011"=0, "2014"=9, "2017"=18, "2020"=27)
  if (hydro_year == "all") {
    base_idx <- scenario_map[scenario]
    return(c(base_idx, base_idx+9, base_idx+18, base_idx+27))
  } else {
    return(hydro_map[hydro_year] + scenario_map[scenario])
  }
}

# --- 3. Helper Function to Get Alternative Mapping ---
# Maps a scenario name (e.g., "PB1") to its corresponding numerical
# alternative IDs for each of the four climate years.
get_scenario_alternatives <- function(scenario, hydro_year) {
  scenario_map <- c("NB"=1, "PB1"=2, "PB2"=3, "PB2b"=4, "PB2c"=5, "PB3"=6, "PB4"=7, "PB5"=8, "PB6"=9)
  hydro_map <- c("2011"=0, "2014"=9, "2017"=18, "2020"=27)
  if (hydro_year == "all") {
    base_idx <- scenario_map[scenario]
    return(c(base_idx, base_idx+9, base_idx+18, base_idx+27))
  } else {
    return(hydro_map[hydro_year] + scenario_map[scenario])
  }
}

# --- 4. Core Function to Calculate Weighted Spawner Abundance ---
# Replicates the Shiny app logic to calculate a weighted average spawner forecast
# based on user-defined weights for climate and TDM models.
calculate_weighted_spawners <- function(scenario, hydro_weights, tdm_weights) {
  alts <- get_scenario_alternatives(scenario, "all")
  hydro_years <- c("2011", "2014", "2017", "2020")
  final_series <- tibble(year = 2025:(2025 + 99), spawners = 0)
  for (i in seq_along(alts)) {
    alt_id <- as.character(alts[i])
    hydro_year <- hydro_years[i]
    alt_data <- results_full %>% filter(env == alt_id, year >= 2025)
    if (nrow(alt_data) == 0) next
    tdm_weighted_spawners <- alt_data %>%
      group_by(year) %>%
      summarise(spawners = sum(case_when(variant == "exp_WF" ~ spawners * tdm_weights["exp_WF"], variant == "exp_SM" ~ spawners * tdm_weights["exp_SM"], variant == "lin_Martin" ~ spawners * tdm_weights["lin_Martin"], TRUE ~ 0), na.rm = TRUE), .groups = "drop") %>%
      arrange(year) %>% pull(spawners)
    final_series$spawners <- final_series$spawners + (tdm_weighted_spawners * hydro_weights[hydro_year])
  }
  return(final_series)
}

# --- 5. Define All Weighting Combinations ---
scenarios <- c("NB", "PB1", "PB2", "PB2b", "PB2c", "PB3", "PB4", "PB5", "PB6")
tdm_full_names <- c("exp_WF" = "Bratovich et al. (2020)", "exp_SM" = "Bartholow & Heasley (2006)", "lin_Martin" = "Martin et al. (2017)")

# Baseline (Status Quo) Weights
baseline_weights <- list(list(
  type = "Baseline",
  hydro_weights = setNames(c(0.25, 0.25, 0.25, 0.25), c("2011", "2014", "2017", "2020")),
  tdm_weights = setNames(c(0.51, 0.24, 0.25), c("exp_WF", "exp_SM", "lin_Martin"))
))

# Generate the 12 Extreme Weight Combinations
extreme_weights_list <- list()
hydro_options <- c("2011", "2014", "2017", "2020")
tdm_options <- c("exp_WF", "exp_SM", "lin_Martin")
for (h_year in hydro_options) {
  for (t_model in tdm_options) {
    hydro_w <- setNames(rep(0, 4), hydro_options); hydro_w[h_year] <- 1.0
    tdm_w <- setNames(rep(0, 3), tdm_options); tdm_w[t_model] <- 1.0
    extreme_weights_list[[length(extreme_weights_list) + 1]] <- list(
      type = "Extreme", climate_year = h_year, tdm_model = tdm_full_names[t_model],
      hydro_weights = hydro_w, tdm_weights = tdm_w
    )
  }
}

# Combine baseline and extreme definitions into one list
all_weight_combos <- c(baseline_weights, extreme_weights_list)

# --- 6. Run Simulations for ALL Combinations (Baseline and Extreme) ---
all_results <- map_dfr(all_weight_combos, function(weights) {
  map_dfr(scenarios, function(scen) {
    calculate_weighted_spawners(
      scenario = scen,
      hydro_weights = weights$hydro_weights,
      tdm_weights = weights$tdm_weights
    ) %>%
      mutate(
        scenario = scen,
        type = weights$type,
        climate_year = if (weights$type == "Extreme") weights$climate_year else "Default Mix",
        tdm_model = if (weights$type == "Extreme") weights$tdm_model else "Default Mix"
      )
  })
})

# --- 7. Generate and Save the BASELINE PLOT ---

# Filter and summarize data for the baseline plot
summary_data_baseline <- all_results %>%
  filter(type == "Baseline", year >= max(year) - 19) %>%
  group_by(scenario) %>%
  summarise(
    median_spawners = median(spawners, na.rm = TRUE),
    p25 = quantile(spawners, 0.25, na.rm = TRUE),
    p75 = quantile(spawners, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# Create the baseline plot object
baseline_barchart <- ggplot(summary_data_baseline, aes(x = factor(scenario, levels = scenarios), y = median_spawners, fill = scenario)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_errorbar(aes(ymin = p25, ymax = p75), width = 0.25, linewidth = 0.7) +
  scale_y_continuous(labels = comma, limits = c(0, NA)) +
  scale_fill_viridis_d(guide = "none") +
  labs(
    title = NULL,
    x = "Management Alternative",
    y = "Adult Population Index"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text   = element_text(face = "bold", size = 11),
    axis.title  = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(size = 11),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5)
  )

# Save the baseline plot
ggsave(
  filename = "spawner_forecast_baseline_barchart.png",
  plot = baseline_barchart,
  width = 6.5, height = 7, dpi = 300
)

# --- 8. Generate and Save the EXTREME SCENARIOS PLOT ---

# Filter and summarize data for the extreme scenarios plot
summary_data_extreme <- all_results %>%
  filter(type == "Extreme", year >= max(year) - 19) %>%
  group_by(scenario, climate_year, tdm_model) %>%
  summarise(
    median_spawners = median(spawners, na.rm = TRUE),
    p25 = quantile(spawners, 0.25, na.rm = TRUE),
    p75 = quantile(spawners, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# Set factor levels for correct grid ordering
summary_data_extreme$climate_year <- factor(summary_data_extreme$climate_year, levels = hydro_options)
summary_data_extreme$tdm_model <- factor(summary_data_extreme$tdm_model, levels = tdm_full_names)

# Create the extreme scenarios plot object
extreme_barchart <- ggplot(summary_data_extreme, aes(x = factor(scenario, levels = scenarios), y = median_spawners, fill = scenario)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_errorbar(aes(ymin = p25, ymax = p75), width = 0.25, linewidth = 0.7) +
  scale_y_continuous(labels = comma, limits = c(0, NA)) +
  scale_fill_viridis_d(guide = "none") +
  facet_grid2(climate_year ~ tdm_model, scales = "free_y", independent = "y") +
  labs(
    title = NULL,
    x = "Management Alternative",
    y = "Adult Population Index"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    strip.text = element_text(face = "bold", size = 9, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.title  = element_text(face = "bold", color = "black"),
    legend.title = element_text(face = "bold", color = "black"),
    legend.text  = element_text(size = 11, color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  )

print(extreme_barchart)

# Save the extreme scenarios plot
ggsave(
  filename = "spawner_forecast_extremes_barchart.png",
  plot = extreme_barchart,
  width = 9, height = 7, dpi = 300
)

# --- End of Script ---
cat("Script finished. Two files have been saved:\n")
cat("- spawner_forecast_baseline_barchart.png\n")
cat("- spawner_forecast_extremes_barchart.png\n")