# ============================================================================
# 2026 bypass schedules, for the app's About tab
# ============================================================================
# Draws the daily bypass flow of the six 2026 scenarios from Reclamation's
# schedule workbook (data_raw/2026BypassModelingScenarios.xlsx, Timeseries
# sheet, received 30 Sep 2026) as one panel per alternative, and writes
# SalmonCountR/www/2026_bypass_schedules.png. years.R names that file as the
# 2026 schedule_image.
#
# Run from the repo root:  source("analysis/bypass_schedules_2026.R")
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(here)
})
source(here("analysis", "figure_theme.R"))

ts <- read_excel(here("data_raw", "2026BypassModelingScenarios.xlsx"), sheet = "Timeseries")

# One row per scenario and day. A day's flow holds from that date to the next,
# so each day is drawn as a block one day wide. The sheet's last row (Dec 1) is
# the end of the bypass: zero, or blank for Scenarios 5 and 6.
sched <- ts %>%
  pivot_longer(starts_with("Scenario"), names_to = "scenario", values_to = "cfs") %>%
  mutate(Date = as.Date(Date),
         cfs  = replace_na(as.numeric(cfs), 0),
         alt  = factor(str_replace(scenario, "Scenario ", "PB"), levels = paste0("PB", 1:6)))

# Bypass volume per alternative, for the panel titles (cfs-days to acre-feet).
vol <- sched %>%
  group_by(alt) %>%
  summarise(af = sum(cfs) * 86400 / 43560, .groups = "drop") %>%
  mutate(strip = sprintf("%s  (%s AF)", alt, format(round(af, -2), big.mark = ",")))
sched <- sched %>%
  left_join(vol, by = "alt") %>%
  mutate(strip = factor(strip, levels = vol$strip))

p <- ggplot(filter(sched, cfs > 0)) +
  geom_rect(aes(xmin = Date, xmax = Date + 1, ymin = 0, ymax = cfs),
            fill = arg_pal(1, begin = 0.35)) +
  facet_wrap(~ strip, ncol = 3, axes = "all") +   # x and y axis on every panel
  scale_x_date(breaks = as.Date(c("2026-10-15", "2026-11-01", "2026-11-15", "2026-12-01")),
               labels = function(d) sub(" 0", " ", format(d, "%b %d")),
               limits = as.Date(c("2026-10-15", "2026-12-01")),
               expand = expansion(mult = 0.02)) +
  scale_y_continuous(breaks = c(0, 250, 500), limits = c(0, 520),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Date (2026)", y = "Bypass flow (cfs)") +
  theme_arg(base_size = 13, border = FALSE) +
  theme(axis.line         = element_line(colour = "black", linewidth = 0.6),
        axis.ticks        = element_line(colour = "black", linewidth = 0.6),
        axis.ticks.length = unit(5, "pt"),
        panel.grid.major.x = element_blank(),
        panel.spacing.x = unit(34, "pt"),
        panel.spacing.y = unit(16, "pt"),
        plot.margin = margin(6, 22, 6, 6))

ggsave(here("SalmonCountR", "www", "2026_bypass_schedules.png"), p,
       width = 10, height = 5.6, dpi = 180)
