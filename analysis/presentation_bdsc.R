# ============================================================================
# Bay-Delta Science Conference 2026 - presentation figures
# ============================================================================
# Dark-surface versions of the talk's figures. Separate from the manuscript
# figure scripts on purpose: nothing here touches analysis/figure_theme.R or any
# published exhibit, so the paper's styling cannot drift because of a talk.
#
# DARK MODE IS RE-STEPPED, NOT FLIPPED. The manuscript palette runs viridis
# 0.10-0.80, whose dark end disappears against a near-black slide. The talk
# takes its hues from the same ramp at 0.50-0.97, so identity is consistent with
# the paper while every series stays legible on the dark surface. Colours still
# follow the entity (a TDM model, an objective), never its rank.
#
# NOT VALIDATED BY SCRIPT: the palette checker needs node, which is not
# installed on this machine. Viridis is perceptually uniform and CVD-safe by
# construction, and three widely spaced steps are a conservative use of it, but
# this was reasoned rather than measured. Re-run the checker before this is used
# anywhere that matters.
#
# Inputs  : output/{figure3_etf_survival_by_temp,reporting_values,
#           consequence_table}.csv
# Outputs : figures/bdsc/*.png
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(ggrepel); library(here); library(viridisLite)
})

OUT <- here("figures", "bdsc")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- Dark surface and ink ---------------------------------------------------
# Measured against SURFACE (WCAG contrast ratio):
#   INK       #FFFFFF  18.2:1   AAA
#   INK_SOFT  #D8DEE6  13.0:1   AAA   (captions only)
#   AXIS      #8D97A6   6.0:1   AAA   (rules and ticks, not text)
#   GRID      #333B47   1.6:1   deliberately below text thresholds - it is a
#                               background rule, and must stay recessive
# All chart text is white. The earlier muted grey measured 7.2:1, which passes
# on a monitor, but a conference projector in a lit room loses most of that.
SURFACE   <- "#12161C"   # slide background
INK       <- "#FFFFFF"   # all chart text
INK_SOFT  <- "#D8DEE6"   # captions
AXIS      <- "#8D97A6"   # axis lines and ticks
GRID      <- "#333B47"   # recessive grid

# Same ramp as the manuscript, stepped for the dark surface. Starting at 0.50
# rather than 0.40 lifts the darkest series from 3.61:1 to 4.69:1 against the
# slide, which matters far more on a projector than the slightly wider hue
# spread the lower start would have given.
pal_dark <- function(n) viridisLite::viridis(n, begin = 0.50, end = 0.97)

TDM_DARK <- setNames(pal_dark(3),
                     c("Bratovich et al. (2020)",
                       "Bartholow & Heasley (2006)",
                       "Martin et al. (2017)"))

OBJ_DARK <- setNames(pal_dark(3), c("Chinook salmon", "Steelhead", "Hydropower"))

theme_talk <- function(base_size = 26, legend = "top") {
  theme_minimal(base_size = base_size) +
    theme(
      text             = element_text(colour = INK),
      plot.title       = element_text(face = "bold", size = base_size + 4, colour = INK),
      plot.subtitle    = element_text(size = base_size - 4, colour = INK,
                                      margin = margin(b = 10)),
      plot.caption     = element_text(size = base_size - 10, colour = INK_SOFT),
      axis.title       = element_text(colour = INK, size = base_size - 3),
      axis.text        = element_text(colour = INK, size = base_size - 4),
      strip.text       = element_text(face = "bold", colour = INK, size = base_size - 1),
      legend.title     = element_blank(),
      legend.text      = element_text(colour = INK, size = base_size - 4),
      legend.position  = legend,
      # Axes are drawn, not implied. A floating panel reads as unanchored from
      # the back of a conference room.
      axis.line        = element_line(colour = AXIS, linewidth = 0.6),
      axis.ticks       = element_line(colour = AXIS, linewidth = 0.6),
      axis.ticks.length = unit(5, "pt"),
      panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = SURFACE, colour = NA),
      plot.background  = element_rect(fill = SURFACE, colour = NA),
      legend.background = element_rect(fill = SURFACE, colour = NA),
      legend.key       = element_rect(fill = SURFACE, colour = NA),
      plot.margin      = margin(16, 24, 12, 16)
    )
}

# Wider and shorter than before: the slide gives a 16:9 well under the heading,
# and a taller image gets letterboxed rather than filling it.
save_talk <- function(plot, file, width = 14, height = 7.0) {
  ggsave(file.path(OUT, file), plot, width = width, height = height,
         dpi = 200, bg = SURFACE)
  cat("wrote", file, "\n")
}

# ============================================================================
# 1. The disagreement -- same temperature, three different answers
# ============================================================================
tdm <- read_csv(here("output", "figure3_etf_survival_by_temp.csv"),
                show_col_types = FALSE) %>%
  pivot_longer(-T_C, names_to = "model", values_to = "surv") %>%
  mutate(model = recode(model,
                        Bratovich_pct = "Bratovich et al. (2020)",
                        Bartholow_pct = "Bartholow & Heasley (2006)",
                        Martin_pct    = "Martin et al. (2017)"))

# Observed October-November range at Hazel Avenue across all alternatives and
# meteorological years, from the manuscript's Figure 3 caption.
OBS_LO <- 14.1; OBS_HI <- 18.2

ends <- tdm %>% group_by(model) %>% slice_max(T_C, n = 1) %>% ungroup()

p1 <- ggplot(tdm, aes(T_C, surv, colour = model)) +
  annotate("rect", xmin = OBS_LO, xmax = OBS_HI, ymin = -Inf, ymax = Inf,
           fill = "#FFFFFF", alpha = 0.06) +
  annotate("text", x = (OBS_LO + OBS_HI) / 2, y = 103,
           label = "observed Oct-Nov range", colour = INK, size = 6.5) +
  geom_line(linewidth = 1.6) +
  geom_text_repel(data = ends, aes(label = model), hjust = 0, nudge_x = 0.35,
                  direction = "y", size = 6.5, segment.colour = AXIS,
                  show.legend = FALSE) +
  scale_colour_manual(values = TDM_DARK, guide = "none") +
  scale_x_continuous(limits = c(10, 26), breaks = seq(10, 22, 2)) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25)) +
  labs(x = "Water temperature (°C)", y = "Egg-to-fry survival (%)") +
  theme_talk(legend = "none")

save_talk(p1, "01_tdm_disagreement.png")

# ============================================================================
# 2. Headline -- projected adult population index by alternative
# ============================================================================
rv <- read_csv(here("output", "reporting_values.csv"), show_col_types = FALSE) %>%
  mutate(
    scenario = factor(scenario, levels = scenario[order(adult_index)]),
    role = case_when(scenario == "NB"  ~ "No bypass (reference)",
                     scenario == "PB4" ~ "Selected by the team",
                     TRUE              ~ "Other alternatives")
  )

REF_FILL <- "#9AA4B2"   # neutral fill for the reference bar (a mark, not text)
ROLE_COLS <- c("No bypass (reference)" = REF_FILL,
               "Selected by the team"  = unname(pal_dark(3)[3]),
               "Other alternatives"    = unname(pal_dark(3)[1]))

p2 <- ggplot(rv, aes(scenario, adult_index, fill = role)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = format(round(adult_index), big.mark = ",")),
            hjust = -0.12, colour = INK, size = 6.5) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = ROLE_COLS, breaks = names(ROLE_COLS)) +
  scale_y_continuous(limits = c(0, 13200), expand = c(0, 0),
                     labels = function(x) format(x, big.mark = ",", trim = TRUE)) +
  labs(x = NULL, y = "Projected adult population index") +
  theme_talk() +
  theme(panel.grid.major.y = element_blank())

save_talk(p2, "02_adult_index.png")

# ============================================================================
# 3. Three objectives, not one trade-off
# ============================================================================
# Deliberately small multiples rather than a 27-bar grouped chart: the point is
# that the alternatives rank differently under each objective, which grouping
# buries. Framing follows the operational context -- several competing
# objectives being balanced, not fish against power.
ct <- read_csv(here("output", "consequence_table.csv"), show_col_types = FALSE) %>%
  select(Alternative,
         `Chinook salmon` = `Chinook (0-1)`,
         Steelhead        = `Steelhead (0-1)`,
         Hydropower       = `Hydro (0-1)`) %>%
  pivot_longer(-Alternative, names_to = "objective", values_to = "score") %>%
  mutate(objective   = factor(objective, levels = names(OBJ_DARK)),
         Alternative = factor(Alternative,
                              levels = c("NB","PB1","PB2","PB2b","PB2c","PB3","PB4","PB5","PB6")))

p3 <- ggplot(ct, aes(Alternative, score, fill = objective)) +
  geom_col(width = 0.74) +
  facet_wrap(~ objective, nrow = 1) +
  scale_fill_manual(values = OBJ_DARK, guide = "none") +
  scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Normalised performance (0-1, higher is better)",
       caption = "0 = worst of the nine by construction (min-max), not zero benefit") +
  theme_talk(base_size = 24, legend = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                   size = 15, colour = INK),
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(1.4, "lines"))

save_talk(p3, "03_three_objectives.png")

# ============================================================================
# 4. Volume is not the whole story -- schedule matters
# ============================================================================
# PB3 and PB5 bypass the identical volume on different schedules. Under the
# corrected model they differ by 261 (+/- 36) adults, a paired within-seed
# contrast rather than a subtraction of the two levels shown here.
vb <- read_csv(here("output", "reporting_values.csv"), show_col_types = FALSE) %>%
  filter(scenario != "NB") %>%
  mutate(pair = scenario %in% c("PB3", "PB5"))

p4 <- ggplot(vb, aes(volume_Mm3, gain_vs_NB)) +
  geom_smooth(method = "lm", se = FALSE, colour = GRID, linewidth = 0.9,
              linetype = "dashed", formula = y ~ x) +
  geom_point(aes(colour = pair, size = pair)) +
  geom_text_repel(aes(label = scenario, colour = pair), size = 6.8,
                  box.padding = 0.55, segment.colour = GRID, show.legend = FALSE) +
  annotate("text", x = 21.4, y = 2450,
           label = "same volume,\ndifferent schedule:\n261 ± 36 adults apart",
           colour = INK, size = 6.5, lineheight = 1.05, hjust = 0.5) +
  scale_colour_manual(values = c(`TRUE` = unname(pal_dark(3)[3]),
                                 `FALSE` = unname(pal_dark(3)[1])), guide = "none") +
  scale_size_manual(values = c(`TRUE` = 6.5, `FALSE` = 4), guide = "none") +
  scale_x_continuous(limits = c(8, 54)) +
  scale_y_continuous(limits = c(0, 4300), labels = function(x) format(x, big.mark = ",")) +
  labs(x = "Bypass volume (million m³)",
       y = "Additional adults vs no bypass") +
  theme_talk(legend = "none")

save_talk(p4, "04_volume_vs_schedule.png")

cat("\nAll BDSC figures written to figures/bdsc/\n")
