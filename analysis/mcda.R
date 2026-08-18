# ============================================================================
# Figure 5. MCDA Composite Scores — Standalone Extraction from Shiny App
# ============================================================================
# Replicates the stacked-bar "Contribution to Score by Objective" plot from
# app.R (output$plot_score_contribution) using pre-computed swing weighting
# results saved by precompute.R.
#
# Required files (in SalmonCountR/app_data/):
#   - swing_scenario_results.rds   (chinook spawner metric per scenario)
#   - steelhead_scenario_results.rds (steelhead score per scenario)
# ============================================================================

library(tidyverse)
library(scales)
library(here)

# --- 1. Load pre-computed scenario results ---
swing_scenario_results    <- readRDS(here("SalmonCountR", "app_data", "swing_scenario_results.rds"))
steelhead_scenario_results <- readRDS(here("SalmonCountR", "app_data", "steelhead_scenario_results.rds"))

# --- 2. Hard-coded hydropower revenue-loss scores ($ from SDM alternatives) ---
hardcoded_hydro_scores <- c(
  "NB"  = 0,      "PB1" = 111422, "PB2" = 370826,
  "PB2b" = 470090, "PB2c" = 433215, "PB3" = 201552,
  "PB4" = 241590, "PB5" = 199382,  "PB6" = 348806
)

# --- 3. Assemble performance matrix ---
perf_data <- swing_scenario_results %>%
  rename(chinook_raw = spawner_metric) %>%
  left_join(
    steelhead_scenario_results %>% rename(steelhead_raw = steelhead_score),
    by = "scenario"
  ) %>%
  left_join(
    tibble(scenario = names(hardcoded_hydro_scores),
           hydro_raw = hardcoded_hydro_scores),
    by = "scenario"
  )

# --- 4. Min-max normalization (matches app.R logic) ---
normalize_minmax <- function(x) (x - min(x)) / (max(x) - min(x))
normalize_invert <- function(x) (max(x) - x) / (max(x) - min(x))  # lower = better

perf_data <- perf_data %>%
  mutate(
    chinook_norm   = normalize_minmax(chinook_raw),
    steelhead_norm = normalize_minmax(steelhead_raw),
    hydro_norm     = normalize_invert(hydro_raw)
  )

# --- 5. Reclamation decision-maker weights ---
w_hydro     <- 0.5
w_chinook   <- 0.4
w_steelhead <- 0.1

# --- 6. Weighted scores for stacked bar ---
scenarios_ordered <- c("NB", "PB1", "PB2", "PB2b", "PB2c", "PB3", "PB4", "PB5", "PB6")

plot_df <- perf_data %>%
  mutate(
    Chinook    = chinook_norm   * w_chinook,
    Steelhead  = steelhead_norm * w_steelhead,
    Hydropower = hydro_norm     * w_hydro
  ) %>%
  select(scenario, Chinook, Steelhead, Hydropower) %>%
  pivot_longer(cols = c(Chinook, Steelhead, Hydropower),
               names_to = "Objective", values_to = "Contribution") %>%
  mutate(
    scenario  = factor(scenario, levels = scenarios_ordered),
    Objective = factor(Objective, levels = c("Hydropower", "Steelhead", "Chinook"))
  )

# --- 6b. Composite totals for the bar labels ---
# The revised Figure 5 caption promises numeric labels: the top-five alternatives
# are separated by only a few thousandths, which the bar heights cannot resolve.
composite_totals <- plot_df %>%
  group_by(scenario) %>%
  summarise(Composite = sum(Contribution), .groups = "drop")

# --- 7. Plot ---
p_mcda <- ggplot(plot_df, aes(x = scenario, y = Contribution, fill = Objective)) +
  geom_col(position = "stack") +
  geom_text(
    data = composite_totals,
    aes(x = scenario, y = Composite, label = sprintf("%.3f", Composite)),
    inherit.aes = FALSE,
    vjust = -0.6, size = 4, fontface = "bold", colour = "black"
  ) +
  scale_fill_viridis_d(name = "Objective") +
  scale_y_continuous(
    limits = c(0, max(composite_totals$Composite) * 1.12),
    breaks = seq(0, 0.7, 0.1)
  ) +
  labs(
    title = NULL,
    x = "Management Alternative",
    y = "Composite Score"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 18, color = "black"),
    plot.subtitle = element_text(size = 12, margin = margin(b = 10), color = "black"),
    axis.title    = element_text(face = "bold", color = "black"),
    axis.text     = element_text(face = "bold", size = 11, color = "black"),
    legend.position = "right",
    legend.title    = element_text(face = "bold", color = "black"),
    legend.text     = element_text(size = 11, color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  )

print(p_mcda)

# Write into figures/ rather than the working directory, so reruns do not leave
# stray copies at the repo root.
dir.create(here("figures"), showWarnings = FALSE)
ggsave(here("figures", "mcda_composite_scores.png"), p_mcda,
       width = 9, height = 7, dpi = 300, bg = "white")

# --- 8. Companion table (values behind the bars) ---
dir.create(here("output"), showWarnings = FALSE)
perf_data %>%
  left_join(composite_totals, by = "scenario") %>%
  arrange(desc(Composite)) %>%
  write.csv(here("output", "mcda_composite_scores_v2.csv"), row.names = FALSE)
