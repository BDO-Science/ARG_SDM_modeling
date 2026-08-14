# ============================================================================
# Figure 4 with run-to-run uncertainty
# ============================================================================
# The published Figure 4 (analysis/figures.R) draws one error bar: the
# interquartile range of the adult population index across the 114 projection
# years WITHIN a single run. That is year-to-year spread of the projection.
#
# It shows nothing about how much the bar itself would move on a different
# random draw of redds. That second quantity — run-to-run uncertainty in the
# estimate — is what determines whether two alternatives are actually
# distinguishable, and pooling (G1) was suppressing it by a factor of ~14.
#
# This script draws both, from the multi-seed snapshots:
#
#   thin black whisker  = year-to-year IQR   (as published, mean across seeds)
#   thick orange bar    = run-to-run range   (min-max of the per-seed medians)
#
# Usage (from the repo root):
#   G1_SNAPROOT=<dir> Rscript analysis/figure4_seed_uncertainty.R
#   G1_ARM=fixed|legacy   (default fixed)
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(scales); library(here); library(ggh4x)
})
source(here("analysis", "figure_theme.R"))

snaproot <- Sys.getenv("G1_SNAPROOT")
arm      <- Sys.getenv("G1_ARM", "fixed")
stopifnot(dir.exists(snaproot))

snaps <- list.dirs(snaproot, recursive = FALSE)
snaps <- snaps[grepl(paste0("^", arm, "_seed[0-9]+$"), basename(snaps))]
if (!length(snaps)) stop("no ", arm, "_seed* snapshots in ", snaproot)
cat(sprintf("arm = %s, %d seed snapshots\n", arm, length(snaps)))

scenarios      <- c("NB", "PB1", "PB2", "PB2b", "PB2c", "PB3", "PB4", "PB5", "PB6")
hydro_options  <- c("2011", "2014", "2017", "2020")
tdm_options    <- c("exp_WF", "exp_SM", "lin_Martin")
tdm_full_names <- c("exp_WF"     = "Bratovich et al. (2020)",
                    "exp_SM"     = "Bartholow & Heasley (2006)",
                    "lin_Martin" = "Martin et al. (2017)")

scenario_map <- c("NB"=1, "PB1"=2, "PB2"=3, "PB2b"=4, "PB2c"=5,
                  "PB3"=6, "PB4"=7, "PB5"=8, "PB6"=9)
hydro_map    <- c("2011"=0, "2014"=9, "2017"=18, "2020"=27)

# Same weighting logic as analysis/figures.R, but taking results_full as an
# argument instead of reading it from the global environment.
weighted_spawners <- function(results_full, scenario, hydro_weights, tdm_w) {
  base <- scenario_map[[scenario]]
  alts <- base + unname(hydro_map)
  out  <- tibble(year = 2025:(2025 + 99), spawners = 0)
  for (i in seq_along(alts)) {
    alt_data <- results_full %>% filter(env == as.character(alts[i]), year >= 2025)
    if (!nrow(alt_data)) next
    s <- alt_data %>%
      group_by(year) %>%
      summarise(spawners = sum(case_when(
        variant == "exp_WF"     ~ spawners * tdm_w["exp_WF"],
        variant == "exp_SM"     ~ spawners * tdm_w["exp_SM"],
        variant == "lin_Martin" ~ spawners * tdm_w["lin_Martin"],
        TRUE ~ 0), na.rm = TRUE), .groups = "drop") %>%
      arrange(year) %>% pull(spawners)
    out$spawners <- out$spawners + s * hydro_weights[[hydro_options[i]]]
  }
  out
}

# Per-seed summary over the 12 extreme (single climate year x single TDM) cells
per_seed <- map_dfr(snaps, function(dir) {
  seed <- as.integer(sub("^.*_seed", "", basename(dir)))
  rf   <- readRDS(file.path(dir, "results_full.rds"))
  map_dfr(hydro_options, function(h) {
    map_dfr(tdm_options, function(t) {
      hw <- setNames(rep(0, 4), hydro_options); hw[h] <- 1
      tw <- setNames(rep(0, 3), tdm_options);   tw[t] <- 1
      map_dfr(scenarios, function(scen) {
        s <- weighted_spawners(rf, scen, hw, tw) %>%
          filter(year >= max(year) - 19)
        tibble(seed = seed, scenario = scen, climate_year = h,
               tdm_model = tdm_full_names[[t]],
               median_spawners = median(s$spawners, na.rm = TRUE),
               p25 = quantile(s$spawners, 0.25, na.rm = TRUE),
               p75 = quantile(s$spawners, 0.75, na.rm = TRUE))
      })
    })
  })
})

# Collapse across seeds: the bar is the mean of the per-seed medians; the
# run-to-run range is the min-max of those medians.
summary_data <- per_seed %>%
  group_by(scenario, climate_year, tdm_model) %>%
  # NOTE: compute the across-seed spread BEFORE collapsing median_spawners.
  # summarise() evaluates sequentially, so assigning median_spawners first
  # would leave min()/max() reading the already-collapsed scalar and report a
  # spread of exactly zero.
  summarise(
    seed_lo = min(median_spawners),
    seed_hi = max(median_spawners),
    seed_sd = sd(median_spawners),
    mean_median = mean(median_spawners),
    p25 = mean(p25), p75 = mean(p75),
    .groups = "drop"
  ) %>%
  rename(median_spawners = mean_median) %>%
  mutate(
    climate_year = factor(climate_year, levels = hydro_options),
    tdm_model    = factor(tdm_model, levels = unname(tdm_full_names))
  )

# Colour carries no identity here — the x axis already names each alternative —
# so the two uncertainty tiers are separated by width AND lightness AND a
# legend, never by hue alone.
col_seed <- "#D55E00"  # Okabe-Ito vermillion: distinct from black and viridis
col_iqr  <- "black"

p <- ggplot(summary_data,
            aes(x = factor(scenario, levels = scenarios), y = median_spawners)) +
  geom_bar(aes(fill = scenario), stat = "identity", alpha = 0.8) +
  geom_errorbar(aes(ymin = p25, ymax = p75, colour = "Year-to-year (IQR)"),
                width = 0.25, linewidth = 0.6) +
  geom_linerange(aes(ymin = seed_lo, ymax = seed_hi,
                     colour = "Run-to-run (5 seeds)"),
                 linewidth = 2.1, alpha = 0.95) +
  scale_y_continuous(labels = comma, limits = c(0, NA)) +
  scale_fill_viridis_d(guide = "none") +
  scale_colour_manual(
    name   = NULL,
    values = c("Year-to-year (IQR)" = col_iqr,
               "Run-to-run (5 seeds)" = col_seed),
    breaks = c("Year-to-year (IQR)", "Run-to-run (5 seeds)")
  ) +
  facet_grid2(climate_year ~ tdm_model, scales = "free_y", independent = "y") +
  labs(x = "Management Alternative", y = "Adult Population Index") +
  theme_arg(base_size = 14) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, face = "bold", colour = "black"),
    strip.text      = element_text(size = 10),
    legend.position = "top"
  ) +
  guides(colour = guide_legend(override.aes = list(linewidth = c(0.8, 2.1))))

ggsave(here("figures", "figure4_with_seed_uncertainty.png"),
       plot = p, width = 11, height = 9, dpi = 300, bg = "white")

# How the two tiers compare, in numbers
cmp <- summary_data %>%
  mutate(seed_width = seed_hi - seed_lo, iqr_width = p75 - p25,
         # Martin projects near-extinction in some cells, where the projection
         # is flat and the year-to-year IQR is exactly 0.
         seed_pct_of_iqr = ifelse(iqr_width > 0, 100 * seed_width / iqr_width,
                                  NA_real_)) %>%
  group_by(tdm_model) %>%
  summarise(mean_median      = mean(median_spawners),
            mean_seed_width  = mean(seed_width),
            mean_seed_sd     = mean(seed_sd),
            mean_iqr_width   = mean(iqr_width),
            mean_seed_pct_of_iqr = mean(seed_pct_of_iqr), .groups = "drop")
cat("\nRun-to-run range vs year-to-year IQR, by TDM model:\n")
print(as.data.frame(cmp), digits = 4, row.names = FALSE)

write_csv(summary_data, here("output", "figure4_seed_uncertainty.csv"))
cat("\nWrote figures/figure4_with_seed_uncertainty.png and",
    "output/figure4_seed_uncertainty.csv\n")
