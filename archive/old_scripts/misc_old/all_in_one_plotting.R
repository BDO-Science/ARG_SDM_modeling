library(ggridges)
library(tidyverse)

#####
#TDM#
#####
###############################
print(results_obs_fast, n = Inf)

ggplot(results_obs_fast, aes(x = sim_year, y = mean_cum_surv, color = variant)) +
  geom_line(size = 1) +
  facet_wrap(~ env, ncol = 2) +
  scale_x_continuous(breaks = pretty(results_obs_fast$sim_year, 10)) +
  labs(
    title = "Observed‐based TDM: Mean Egg→Emergence Survival by Brood Year",
    x     = "Brood Year",
    y     = "Mean cumulative survival",
    color = "Variant"
  ) +
  theme_minimal(base_size = 14)

ggplot(results_obs_fast, aes(x = variant, y = mean_cum_surv)) +
  geom_boxplot(outlier.size = 1) +
  facet_wrap(~ env) +
  labs(
    title = "Distribution of Brood‐Year Survival by Variant & Env",
    x     = "TDM Variant",
    y     = "Mean cumulative survival"
  ) +
  theme_minimal(base_size = 14)

###################
#LIFE CYCLE MODELING
###################
pct_results <- results_full_altvar %>%
  group_by(env, year) %>%
  mutate(
    max_spawners = max(spawners, na.rm = TRUE),
    pct_diff     = (spawners - max_spawners) / max_spawners * 100
  ) %>%
  ungroup() %>%
  select(env, alt_variant = variant_full, variant, year, spawners, pct_diff)



last_year <- max(pct_results$year)
cat("% differences from best variant in year", last_year, ":\n")
pct_results %>%
  filter(year == last_year) %>%
  arrange(env, desc(pct_diff)) %>%
  print(n = Inf)

# 10) PLOTTING -----------------------------------------------------------------
# define the exact order you want
alt_levels <- paste0("Alt ", 1:10)
tdm_levels <- c("exp_SM","exp_WF","lin_Martin") 

results2 <- results_full_altvar %>%
  mutate(
    # 1) env → fill aesthetic
    env = factor(env,
                 levels = as.character(1:10),    # original codes 1–10
                 labels = alt_levels),           # prettied-up labels
    
    # 2) variant → x-axis aesthetic
    variant = factor(variant,
                     levels = tdm_levels)      # in the order you like
  )

# now re-calc max
max_spawn <- max(results2$spawners, na.rm = TRUE)
last_year <- max(results2$year)

# A) Time series
p1 <- ggplot(results2, aes(x = year, y = spawners, color = variant)) +
  geom_line(size = 1) +
  facet_wrap(~ env) +
  scale_y_continuous(limits = c(0, max_spawn), breaks = pretty_breaks(10), labels = comma) +
  labs(x = "Simulation Year", y = "Spawners", color = "TDM Variant") +
  theme_minimal() + theme(legend.position = "bottom")

print(p1)


ggplot(results2, aes(x = year, y = variant, fill = spawners)) +
  geom_tile() +
  facet_wrap(~env, ncol = 2) +
  scale_fill_viridis_c(option = "magma", name = "Spawners") +
  labs(x = "Year", y = "TDM Variant") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7))

final_decade <- results2 %>% filter(year >= (max(year)-9))

# 1) p2: end‑of‑simulation barplot
p2 <- ggplot(
  filter(results2, year == last_year),
  aes(x = variant, y = spawners, fill = env)
) +
  geom_col(position = position_dodge(width = 0.8)) +
  scale_y_continuous(
    limits = c(0, 24000),                     # <-- clamp y to 0–40k
    breaks  = pretty_breaks(10),
    labels  = comma
  ) +
  scale_fill_discrete(breaks = alt_levels) +
  labs(
    title = paste("End‑of‑Simulation Spawners (Year", last_year, ")"),
    x     = "TDM Variant",
    y     = "Spawners",
    fill  = "Alternative"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p2)

# 2) facet_grid line‐plot
ggplot(results2, aes(x = year, y = spawners)) +
  geom_line(color = "steelblue", size = 0.8) +
  facet_grid(variant ~ env, scales = "free_y") +
  scale_y_continuous(
    limits = c(0, 40000),                     # <-- clamp here as well
    labels  = comma
  ) +
  labs(x = "Year", y = "Spawners") +
  theme_minimal() +
  theme(
    strip.text   = element_text(size = 6),
    axis.text.y  = element_text(size = 6)
  )

# 3) boxplot of final decade
ggplot(final_decade, aes(x = variant, y = spawners)) +
  geom_boxplot() +
  facet_wrap(~env, nrow = 2) +
  scale_y_continuous(
    limits = c(0, 22000),                     # <-- and here
    labels  = comma
  ) +
  labs(
    title = "Distribution of Spawners in Final 10 Years by Variant",
    x     = "TDM Variant",
    y     = "Spawners"
  ) +
  theme_minimal()

# 4) density‑ridge: spawners is on the x‑axis, so clamp x instead
ggplot(results2, aes(x = spawners, y = variant, fill = variant)) +
  geom_density_ridges(alpha = 0.7, scale = 1.2) +
  facet_wrap(~env, ncol = 2) +
  scale_x_continuous(
    limits = c(0, 40000),                     # <-- clamp the x‑axis
    labels  = comma
  ) +
  labs(x = "Spawners", y = "TDM Variant") +
  theme_minimal() +
  theme(legend.position = "none")
