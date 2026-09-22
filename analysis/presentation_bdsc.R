# ============================================================================
# Bay-Delta Science Conference 2026 - presentation figures
# ============================================================================
# Dark-surface versions of the talk's figures. Separate from the manuscript
# figure scripts on purpose: nothing here touches analysis/figure_theme.R or any
# published exhibit, so the paper's styling cannot drift because of a talk.
#
# PALETTE. Categorical series use a four-colour Okabe-Ito subset (CAT_SAFE,
# below) chosen for the dark surface; yellow is reserved for one meaning only,
# "this is where water temperature acts". The manuscript's viridis ramp is not
# used here - see the palette block for why it was tried and abandoned. Colours
# follow the entity (a TDM model, an objective), never its rank.
#
# Inputs  : output/{figure3_etf_survival_by_temp,reporting_values,
#           consequence_table,calibration_predictions}.csv
#           SalmonCountR/app_data/{results_full,env_ext_list,spawn_timing_model,
#           base_P,calib_pred_by_variant,american_river_instream}.rds
#           data_raw/study_area_spatial.rds        (clipped CDFW map cache)
#           figures/bdsc/lifecycle/0[1-5]_*.png    (life-stage illustrations)
# Outputs : figures/bdsc/*.png
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(ggrepel); library(here)
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

# Same ramp as the manuscript, stepped for the dark surface, and TRUNCATED SHORT
# OF YELLOW. Yellow carries one meaning in this deck - this is where water
# temperature acts - so it cannot also be series three of a categorical scale.
# It was previously doing both, which meant the selected alternative, the
# hydropower objective and a capacity level all came out the same colour as the
# temperature annotations.
#
# Truncating viridis was tried first and abandoned. Viridis below about 0.45 is
# too dark to survive a projector (0.25 measures 2.15:1 against this surface),
# so cutting the yellow end leaves only 0.45-0.80 - teal to green - and three
# steps inside that band are not separable by hue. On the TDM figure the
# Bratovich and Bartholow curves became untraceable where they cross.
#
# So the categorical scale is Okabe-Ito instead: the standard qualitative set
# built for colour-vision deficiency, which gives real hue separation without
# relying on lightness alone. Members that fail here are dropped - #0072B2 blue
# measures 2.9:1 against the slide, and the amber and yellow entries are held
# back so nothing competes with the temperature yellow.
#
# Contrast against SURFACE: sky blue 8.2:1, green 5.4:1, pink 5.5:1,
# vermillion 4.7:1. All clear 4.5:1.
#
# This does mean the talk no longer shares hues with the manuscript's viridis
# figures. That consistency was already partly gone - the dark surface forced a
# re-step - and legibility in a dark room matters more than matching a printed
# page nobody is holding.
CAT_SAFE <- c("#56B4E9",  # sky blue
              "#009E73",  # bluish green
              "#CC79A7",  # reddish purple
              "#D55E00")  # vermillion

pal_dark <- function(n) {
  if (n <= length(CAT_SAFE)) return(CAT_SAFE[seq_len(n)])
  # More categories than the safe set holds. Interpolating would break the CVD
  # guarantee, so say so rather than silently returning unsafe colours.
  warning("pal_dark(): ", n, " categories exceeds the ", length(CAT_SAFE),
          " CVD-safe steps; falling back to interpolation, which is NOT verified")
  grDevices::colorRampPalette(CAT_SAFE)(n)
}

# The one place yellow is legitimate as a data colour: a warm-versus-cool
# contrast IS the temperature meaning, not a competing one.
TEMP_WARM <- "#F1E51D"

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
      # Bold axes: the presenter opens every figure by reading the axes aloud,
      # and from the back of the room regular-weight axis text is the first
      # thing to go.
      axis.title       = element_text(face = "bold", colour = INK, size = base_size - 3),
      axis.text        = element_text(face = "bold", colour = INK, size = base_size - 4),
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

# October-November range at Hazel Avenue across all 36 alternative x met-year
# runs: the 5th and 95th percentiles that analysis/figure3_tdm_curves.R prints
# from env_ext_list. Mostly CE-QUAL-W2 scenario output (the 2011-2024 Octobers
# and Novembers in that series are the USGS record), so it is NOT an observed
# range and must not be labelled as one.
OBS_LO <- 14.1; OBS_HI <- 18.2

ends <- tdm %>% group_by(model) %>% slice_max(T_C, n = 1) %>% ungroup()

p1 <- ggplot(tdm, aes(T_C, surv, colour = model)) +
  annotate("rect", xmin = OBS_LO, xmax = OBS_HI, ymin = -Inf, ymax = Inf,
           fill = "#FFFFFF", alpha = 0.06) +
  annotate("text", x = (OBS_LO + OBS_HI) / 2, y = 103,
           label = "Oct-Nov range across the model runs", colour = INK, size = 6.5) +
  geom_line(linewidth = 1.6) +
  geom_text_repel(data = ends, aes(label = model), hjust = 0, nudge_x = 0.35,
                  direction = "y", size = 6.5, segment.colour = AXIS,
                  show.legend = FALSE) +
  scale_colour_manual(values = TDM_DARK, guide = "none") +
  scale_x_continuous(limits = c(10, 25), breaks = seq(10, 18, 2)) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25)) +
  # These are S_cum from analysis/figure3_tdm_curves.R: cumulative egg-to-fry
  # survival at a CONSTANT temperature, not a daily rate. The constant-
  # temperature qualifier matters - real incubation sees a falling autumn curve.
  # "Cumulative" needs a period attached or the axis is unanswerable. At a
  # constant temperature the incubation length is set by the ATU thresholds
  # (958 ATU to emergence), so it runs 96 days at 10 C down to 54 days at 18 C
  # (18 x 53 = 954 < 958, so emergence lands on day 54; .slice_by_atu agrees).
  labs(x = "Constant water temperature (°C)",
       y = "Egg-to-emergence survival (%)",
       subtitle = "Cumulative over the whole incubation: 96 days at 10 °C, 54 days at 18 °C") +
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
  # Sits in the empty upper-left rather than over the trend line: at the old
  # x = 21.4 the right edge of the label ran into PB4.
  annotate("text", x = 14.5, y = 2850,
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

# ============================================================================
# 5. Does the TDM choice change the ranking?
# ============================================================================
# A bump chart answers this directly: flat lines at the top mean the preferred
# alternatives are stable across formulations; crossing lines at the bottom mean
# the disagreement is about which alternative is WORST. Bars would need three
# free y-scales (Bratovich is in the thousands, Martin in the tens), which buries
# the ranking question under a scale problem.
rf <- readRDS(here("SalmonCountR", "app_data", "results_full.rds"))
scen_map <- c(NB=1,PB1=2,PB2=3,PB2b=4,PB2c=5,PB3=6,PB4=7,PB5=8,PB6=9)
TDM_LAB  <- c(exp_WF = "Bratovich", exp_SM = "Bartholow", lin_Martin = "Martin")

rank_by_tdm <- rf %>%
  filter(year > 2024) %>%
  group_by(env, variant) %>% slice_tail(n = 20) %>%
  summarise(med = median(spawners, na.rm = TRUE), .groups = "drop") %>%
  mutate(env = as.integer(env),
         scenario = names(scen_map)[match((env - 1) %% 9 + 1, scen_map)]) %>%
  group_by(scenario, variant) %>%
  summarise(value = sum(med * 0.25), .groups = "drop") %>%
  group_by(variant) %>%
  mutate(rank = rank(-value)) %>%
  ungroup() %>%
  mutate(model = factor(TDM_LAB[variant], levels = unname(TDM_LAB)))

# Colour only the alternatives whose rank actually moves; the rest stay neutral
movers <- rank_by_tdm %>% group_by(scenario) %>%
  summarise(span = max(rank) - min(rank), .groups = "drop") %>%
  filter(span >= 3) %>% pull(scenario)

rank_by_tdm <- rank_by_tdm %>%
  mutate(grp = ifelse(scenario %in% movers, scenario, "stable"))
mover_cols <- setNames(pal_dark(length(movers)), movers)
mover_cols["stable"] <- "#7A8492"

ends_r <- rank_by_tdm %>% filter(model == levels(model)[nlevels(model)])
starts_r <- rank_by_tdm %>% filter(model == levels(model)[1])

p5 <- ggplot(rank_by_tdm, aes(model, rank, group = scenario, colour = grp)) +
  geom_line(aes(linewidth = grp %in% movers)) +
  geom_point(size = 4) +
  geom_text(data = starts_r, aes(label = scenario), hjust = 1.35, size = 6, show.legend = FALSE) +
  geom_text(data = ends_r,   aes(label = scenario), hjust = -0.35, size = 6, show.legend = FALSE) +
  scale_y_reverse(breaks = 1:9, expand = expansion(add = 0.6)) +
  scale_x_discrete(expand = expansion(add = c(0.55, 0.55))) +
  scale_colour_manual(values = mover_cols, guide = "none") +
  scale_linewidth_manual(values = c(`TRUE` = 1.8, `FALSE` = 0.8), guide = "none") +
  labs(x = NULL, y = "Rank by adult population index  (1 = best)") +
  theme_talk(legend = "none") +
  theme(panel.grid.major.x = element_blank())

save_talk(p5, "05_rank_by_tdm.png", width = 12, height = 7.2)


# ============================================================================
# 6. Study area -- context and reach
# ============================================================================
# Replaces the manuscript's Figure 1 for the talk. Two co-equal panels rather
# than a reach map with a near-empty state inset: a Bay-Delta audience knows the
# estuary and needs the American River placed against it.
#
# Geometry comes from data_raw/study_area_spatial.rds, a small clipped cache
# built once from the CDFW 100k stream layer (cdfg_100k_2003_6) and California
# Lakes, which live in the sibling Reclamation repos rather than here. Caching
# it keeps this script runnable on a machine that does not have those.
#
# deltamapr's WW_Watershed was tried first and dropped: it is Delta-focused, so
# the Sacramento River essentially vanishes north of the confluence, and the
# American River highlight had a hole in it where Lake Natoma sits.
suppressPackageStartupMessages({library(sf); library(patchwork)})
old_s2 <- sf::sf_use_s2(); sf::sf_use_s2(FALSE)
on.exit(sf::sf_use_s2(old_s2), add = TRUE)

geo <- readRDS(here("data_raw", "study_area_spatial.rds"))
RIVER    <- "#4A7FB5"                 # named rivers, light enough to read on dark
RIVER_HI <- unname(pal_dark(3)[1])    # the American, same teal as elsewhere
DELTA_F  <- "#1B3A5B"

amer   <- geo$rivers %>% filter(NAME == "American River")
others <- geo$rivers %>% filter(NAME != "American River")

stations <- tibble::tribble(
  ~name,          ~lon,       ~lat,     ~kind,
  "Folsom Dam",   -121.1558,  38.7075,  "dam",
  "Nimbus Dam",   -121.2211,  38.6353,  "dam",
  "Hazel Ave",    -121.2266,  38.6355,  "gauge",
  "Watt Ave",     -121.3899,  38.5666,  "gauge"
) %>% st_as_sf(coords = c("lon","lat"), crs = 4326, remove = FALSE)

reach_box <- st_as_sfc(st_bbox(c(xmin = -121.52, ymin = 38.50,
                                 xmax = -121.10, ymax = 38.76), crs = 4326))

map_theme <- function() {
  theme_talk(base_size = 20, legend = "none") +
    theme(panel.grid.major = element_line(colour = "#232B36", linewidth = 0.25),
          axis.title = element_blank(),
          axis.text  = element_text(size = 12, colour = INK_SOFT),
          axis.line = element_blank(), axis.ticks = element_blank(),
          plot.subtitle = element_text(size = 17, colour = INK))
}

pA_map <- ggplot() +
  geom_sf(data = geo$delta, fill = DELTA_F, colour = NA) +
  geom_sf(data = others, colour = RIVER, linewidth = 0.7) +
  geom_sf(data = amer,   colour = RIVER_HI, linewidth = 1.7) +
  geom_sf(data = reach_box, fill = NA, colour = "#F1E51D", linewidth = 0.9) +
  annotate("text", x = -121.02, y = 38.80, label = "American R.",
           colour = RIVER_HI, size = 6.2, fontface = "bold", hjust = 1) +
  annotate("text", x = -121.90, y = 39.60, label = "Sacramento R.",
           colour = INK, size = 5, hjust = 0.5) +
  annotate("text", x = -121.60, y = 38.15, label = "Delta",
           colour = INK, size = 5.2, hjust = 0.5) +
  annotate("text", x = -122.02, y = 37.68, label = "San Francisco\nBay",
           colour = INK, size = 4.8, hjust = 0, lineheight = 0.95) +
  coord_sf(xlim = c(-122.6, -120.8), ylim = c(37.3, 40.1), expand = FALSE) +
  labs(subtitle = "Sacramento Valley to the Bay") +
  map_theme()

pB_map <- ggplot() +
  geom_sf(data = others, colour = RIVER, linewidth = 1.0) +
  geom_sf(data = geo$lakes, fill = DELTA_F, colour = RIVER, linewidth = 0.3) +
  geom_sf(data = amer, colour = RIVER_HI, linewidth = 2.0) +
  geom_sf(data = filter(stations, kind == "dam"),   colour = "#FFFFFF", size = 4.2, shape = 15) +
  geom_sf(data = filter(stations, kind == "gauge"), colour = "#F1E51D", size = 4.6, shape = 19) +
  ggrepel::geom_text_repel(data = stations, aes(lon, lat, label = name),
                           colour = INK, size = 5.4, box.padding = 0.7,
                           segment.colour = AXIS, seed = 1) +
  annotate("text", x = -121.135, y = 38.745, label = "Folsom Lake",
           colour = INK_SOFT, size = 4.4, hjust = 1, fontface = "italic") +
  annotate("text", x = -121.512, y = 38.535, label = "Sacramento River,
to the Bay",
           colour = INK, size = 4.8, hjust = 0, lineheight = 0.95) +
  coord_sf(xlim = c(-121.52, -121.10), ylim = c(38.50, 38.76), expand = FALSE) +
  labs(subtitle = "Lower American River spawning reach") +
  map_theme()

fig_map <- pA_map + pB_map + plot_layout(widths = c(1, 1.25))
save_talk(fig_map, "06_study_area.png", width = 14, height = 6.6)

# ============================================================================
# 7. One figure: the life cycle laid onto the corridor
# ============================================================================
# Merges "where this happens" and "where the model touches the life cycle".
# The argument the layout makes: every place temperature enters the model sits
# in one short reach over about two months, while the life it governs spans the
# estuary and the ocean over several years. That asymmetry is why a two-month
# operating decision can move adult returns at all.
suppressPackageStartupMessages({library(cowplot); library(magick)})

LC <- here("figures", "bdsc", "lifecycle")
stage_img <- function(f) image_read(file.path(LC, f))

# Corridor: narrow river on the left, opening to the ocean on the right.
X0 <- 0.035; X1 <- 0.980; YC <- 0.495
band <- tibble::tibble(
  x = c(X0, 0.32, 0.60, 0.80, X1,  X1,  0.80, 0.60, 0.32, X0),
  y = c(YC + 0.010, YC + 0.014, YC + 0.026, YC + 0.048, YC + 0.085,
        YC - 0.085, YC - 0.048, YC - 0.026, YC - 0.014, YC - 0.010)
)

corridor <- ggplot() +
  geom_polygon(data = band, aes(x, y), fill = "#1B3A5B") +
  geom_segment(aes(x = X0, xend = X1, y = YC, yend = YC),
               colour = pal_dark(3)[1], linewidth = 1.4, alpha = 0.9) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void() +
  theme(plot.background  = element_rect(fill = SURFACE, colour = NA),
        panel.background = element_rect(fill = SURFACE, colour = NA))

XS <- c(spawn = 0.095, egg = 0.275, parr = 0.470, smolt = 0.665, ocean = 0.880)

g <- ggdraw(corridor) +
  draw_image(stage_img("01_spawning_adult.png"), x = XS[["spawn"]] - 0.080, y = 0.575, width = 0.160, height = 0.235) +
  draw_image(stage_img("02_eggs.png"),           x = XS[["egg"]]   - 0.058, y = 0.592, width = 0.116, height = 0.200) +
  draw_image(stage_img("03_parr.png"),           x = XS[["parr"]]  - 0.085, y = 0.585, width = 0.170, height = 0.215) +
  draw_image(stage_img("04_smolt.png"),          x = XS[["smolt"]] - 0.078, y = 0.592, width = 0.156, height = 0.200) +
  draw_image(stage_img("05_ocean_adult.png"),    x = XS[["ocean"]] - 0.098, y = 0.575, width = 0.196, height = 0.235)

labels <- tibble::tribble(
  ~x,             ~top,             ~sub,
  XS[["spawn"]],  "Spawners",       "Nimbus to Watt",
  XS[["egg"]],    "Eggs",           "in the gravel",
  XS[["parr"]],   "Fry / parr",     "lower river",
  XS[["smolt"]],  "Smolts",         "Delta and Bay",
  XS[["ocean"]],  "Ocean returns",  "age 3 to 5"
)
for (i in seq_len(nrow(labels))) {
  g <- g +
    draw_label(labels$top[i], x = labels$x[i], y = 0.385, size = 16.5,
               colour = INK, fontface = "bold", hjust = 0.5) +
    draw_label(labels$sub[i], x = labels$x[i], y = 0.328, size = 13,
               colour = INK_SOFT, hjust = 0.5)
}

g <- g +
  # where temperature enters, bracketed over the first two stages only
  draw_line(x = c(X0, 0.335), y = c(0.870, 0.870), colour = "#F1E51D", linewidth = 1) +
  draw_line(x = c(X0, X0),    y = c(0.870, 0.838), colour = "#F1E51D", linewidth = 1) +
  draw_line(x = c(0.335, 0.335), y = c(0.870, 0.838), colour = "#F1E51D", linewidth = 1) +
  draw_label("Temperature enters here. Only here.",
             x = 0.185, y = 0.935, size = 19, colour = "#F1E51D",
             fontface = "bold", hjust = 0.5) +
  draw_label("spawn timing   pre-spawn survival   egg-to-fry survival",
             x = 0.185, y = 0.822, size = 13.5, colour = "#F1E51D", hjust = 0.5) +
  # anchors sitting on the corridor itself
  draw_label("Folsom Dam", x = X0, y = YC + 0.048, size = 12.5,
             colour = INK_SOFT, hjust = 0) +
  # time, the other half of the asymmetry
  draw_line(x = c(X0, X1), y = c(0.235, 0.235), colour = AXIS, linewidth = 0.7) +
  draw_line(x = c(X0, X0),       y = c(0.235, 0.212), colour = AXIS, linewidth = 0.7) +
  draw_line(x = c(0.335, 0.335), y = c(0.235, 0.212), colour = AXIS, linewidth = 0.7) +
  draw_line(x = c(X1, X1),       y = c(0.235, 0.212), colour = AXIS, linewidth = 0.7) +
  draw_label("October to January", x = 0.185, y = 0.155, size = 14.5, colour = INK,      hjust = 0.5) +
  draw_label("spring",             x = 0.560, y = 0.155, size = 14.5, colour = INK_SOFT, hjust = 0.5) +
  draw_label("2 to 4 years",       x = 0.900, y = 0.155, size = 14.5, colour = INK_SOFT, hjust = 0.5)

ggsave(file.path(OUT, "07_corridor_lifecycle.png"), g,
       width = 14, height = 5.9, dpi = 200, bg = SURFACE)
cat("wrote 07_corridor_lifecycle.png\n")




# ============================================================================
# 8. The life cycle placed on the map
# ============================================================================
# Two columns. Left is real geography, right is the reach blown up with the
# stages that occur in it. Yellow means one thing throughout: temperature. It is
# spent on four things only — the headline, the reach box and its connector, the
# inset border, and the two temperature-sensitive stages. Folsom Lake is drawn
# as an ordinary reservoir: it was yellow to mark the cold-water source, but
# that put a second, competing yellow object next to the box and diluted the
# one meaning the colour is supposed to carry.
OCEAN_F <- "#16324D"
RES_F   <- "#24507A"

# Map labels sit over water, land and grid lines at once. A filled label box
# would occlude the geography it names, so each label gets a dark halo instead:
# the text is drawn 16 times in the surface colour on a small circle, then once
# in its real colour on top. `r` is in degrees and must be set per panel, since
# the wide map and the inset are at very different scales; `yf` compensates for
# the non-square aspect so the halo is even on all sides.
halo_text <- function(x, y, label, colour = INK, size = 5.6, hjust = 0.5,
                      fontface = "plain", r = 0.010, yf = 0.55,
                      halo = "#0B0E13") {
  ang <- head(seq(0, 2 * pi, length.out = 17), -1)
  c(lapply(ang, function(a)
      annotate("text", x = x + r * cos(a), y = y + r * yf * sin(a),
               label = label, colour = halo, size = size, hjust = hjust,
               fontface = fontface)),
    list(annotate("text", x = x, y = y, label = label, colour = colour,
                  size = size, hjust = hjust, fontface = fontface)))
}

# ---- Layout geometry, derived rather than nudged -----------------------------
# The connector kept missing the inset's corners because two separate offsets
# sit between a draw_plot rectangle and the yellow border you can actually see:
#   1. `reach` carries a 1pt plot.margin, so the panel is inset from the
#      rectangle by that much on every side.
#   2. coord_sf holds a true geographic aspect. If the rectangle's aspect does
#      not match the data's, ggplot letterboxes the PANEL inside it and centres
#      the remainder - and panel.border is drawn on the panel, not the
#      rectangle. That gap was about 0.005 of figure width, which is exactly
#      the "hair off" that eyeballing could not close.
# So: pick the inset's longitude window FROM the rectangle's aspect, which
# makes the letterbox zero by construction, then derive the corners.
FIG_W <- 16.5   # the slide well is wider than it is tall once the heading is
FIG_H <- 8.0    # taken out, so a wider figure uses room a 15x8 left empty
MAP_W <- 0.532
INS_X <- 0.539; INS_W <- 0.456
INS_Y <- 0.300; INS_H <- 0.655
INS_M <- 1 / 72                      # `reach` plot.margin, inches

WIDE_LON <- c(-123.75, -120.90)      # extended west: more ocean to label into
WIDE_LAT <- c( 37.25,   39.05)

REACH_LAT <- c(38.50, 38.75)
REACH_CTR <- -121.325
ins_rect_ar <- (INS_W * FIG_W - 2 * INS_M) / (INS_H * FIG_H - 2 * INS_M)
reach_lon_span <- ins_rect_ar * diff(REACH_LAT) /
                  cos(mean(REACH_LAT) * pi / 180)
REACH_LON <- REACH_CTR + c(-1, 1) * reach_lon_span / 2

# Visible border of the inset. No aspect padding now, so only the margin.
INS_L <- INS_X + INS_M / FIG_W
INS_T <- INS_Y + INS_H - INS_M / FIG_H
INS_B <- INS_Y + INS_M / FIG_H

# The wide map stays width-filling and letterboxed in y; same derivation as
# before, now in terms of the constants.
map_ar <- (diff(WIDE_LON) * cos(mean(WIDE_LAT) * pi / 180)) / diff(WIDE_LAT)
draw_h <- (MAP_W * FIG_W / map_ar) / FIG_H
pad_y  <- (1 - draw_h) / 2
fx <- function(lon) ((lon - WIDE_LON[1]) / diff(WIDE_LON)) * MAP_W
fy <- function(lat) pad_y + ((lat - WIDE_LAT[1]) / diff(WIDE_LAT)) * draw_h

resv   <- geo$reservoirs
folsom <- resv[resv$NAME == "Folsom Lake", ]
# The box and the inset now show exactly the same ground, so "this rectangle
# opens into that one" is literally true rather than approximately.
reach_box <- st_as_sfc(st_bbox(c(xmin = REACH_LON[1], ymin = REACH_LAT[1],
                                 xmax = REACH_LON[2], ymax = REACH_LAT[2]),
                               crs = 4326))

wide <- ggplot() +
  geom_sf(data = geo$ocean, fill = OCEAN_F, colour = NA) +
  geom_sf(data = geo$delta, fill = DELTA_F, colour = NA) +
  # rivers first, reservoirs over them: several channels run straight through a
  # reservoir and read as a line drawn on top of the water body otherwise
  geom_sf(data = geo$rivers, colour = RIVER, linewidth = 0.85) +
  geom_sf(data = resv, fill = RES_F, colour = RES_F, linewidth = 0.3) +
  geom_sf(data = reach_box, fill = NA, colour = "#F1E51D", linewidth = 1.2) +
  # The audience knows salmon are anadromous; the migration arrow was telling
  # them so. Removed, which also frees the left half for the ocean labels.
  # Place names run about 12% smaller than they did: they are orientation, not
  # content, and at the previous size they were competing with the inset.
  halo_text(-123.72, 38.02, "Pacific Ocean", size = 5.6, hjust = 0,
            fontface = "italic") +
  halo_text(-122.80, 37.78, "Gulf of the Farallones", colour = INK_SOFT,
            size = 4.4, hjust = 1, fontface = "italic") +
  halo_text(-121.95, 38.94, "Sacramento River", size = 5.1) +
  halo_text(-121.62, 38.06, "Sacramento–San Joaquin Delta", size = 4.7) +
  halo_text(-120.97, 38.92, "Folsom Lake", size = 4.7, hjust = 1) +
  halo_text(-121.37, 38.78, "American River", size = 4.9) +
  coord_sf(xlim = WIDE_LON, ylim = WIDE_LAT, expand = FALSE) +
  # This is a locator map, not a measuring instrument. The graticule labels cost
  # real estate on all four sides and nobody reads a coordinate off a slide, so
  # dropping them is what actually makes the map bigger.
  theme_void() +
  theme(plot.background  = element_rect(fill = SURFACE, colour = NA),
        panel.background = element_rect(fill = SURFACE, colour = NA),
        plot.margin = margin(0, 0, 0, 0))

reach <- ggplot() +
  geom_sf(data = geo$rivers, colour = RIVER, linewidth = 1.6) +
  geom_sf(data = geo$lakes, fill = RES_F, colour = RIVER, linewidth = 0.3) +
  geom_sf(data = filter(stations, kind == "dam"),   colour = "#FFFFFF", size = 5.2, shape = 15) +
  geom_sf(data = filter(stations, kind == "gauge"), colour = "#FFFFFF", size = 5.6, shape = 21, stroke = 1.4, fill = SURFACE) +
  halo_text(-121.172, 38.712, "Folsom Dam", size = 6.3, hjust = 1,
            r = 0.0014, yf = 0.80) +
  halo_text(-121.238, 38.652, "Nimbus Dam", size = 6.3, hjust = 1,
            r = 0.0014, yf = 0.80) +
  halo_text(-121.224, 38.606, "Hazel Ave", size = 6.3, hjust = 1,
            r = 0.0014, yf = 0.80) +
  halo_text(-121.384, 38.543, "Watt Ave", size = 6.3, hjust = 0.5,
            r = 0.0014, yf = 0.80) +
  coord_sf(xlim = REACH_LON, ylim = REACH_LAT, expand = FALSE) +
  theme_void() +
  theme(plot.background  = element_rect(fill = "#171D26", colour = NA),
        panel.background = element_rect(fill = "#171D26", colour = NA),
        panel.border = element_rect(fill = NA, colour = "#F1E51D", linewidth = 1.8),
        plot.margin = margin(1, 1, 1, 1))

# Where the reach box lands on the assembled figure. coord_sf holds a true
# geographic aspect, so with the graticule gone the panel fills MAP_W in x and
# is letterboxed in y: at 38.15 degrees N the 2.55 deg of longitude are worth
# 2.55 * cos(38.15) = 2.004 deg of latitude against 1.80 deg of height, so the
# drawn map is 1.113 times wider than tall and the spare height splits evenly
# above and below. Deriving the corners rather than eyeballing them is what lets
# the connector actually meet the box.

fig8 <- ggdraw() +
  draw_plot(wide, x = 0.000, y = 0.000, width = MAP_W, height = 1.000) +
  # A magnifier frustum, not a pointer: one line from each right-hand corner of
  # the reach box to the nearest corner of the inset. It reads as "this rectangle
  # opens into that one" without needing an arrowhead to say so.
  # Both ends are computed: box corner from fx/fy, inset corner from INS_L/T/B.
  draw_line(x = c(fx(REACH_LON[2]), INS_L), y = c(fy(REACH_LAT[2]), INS_T),
            colour = "#F1E51D", linewidth = 0.7) +
  draw_line(x = c(fx(REACH_LON[2]), INS_L), y = c(fy(REACH_LAT[1]), INS_B),
            colour = "#F1E51D", linewidth = 0.7) +
  draw_label("Temperature enters here. Only here.",
             x = INS_X, y = 0.972, size = 21, colour = "#F1E51D",
             fontface = "bold", hjust = 0) +
  draw_plot(reach, x = INS_X, y = INS_Y, width = INS_W, height = INS_H) +
  # Stage row: tightened about 15% horizontally and each label pulled up toward
  # its own image, so the three read as one group rather than three islands.
  draw_image(stage_img("01_spawning_adult.png"), x = 0.572, y = 0.078, width = 0.152, height = 0.192) +
  draw_image(stage_img("02_eggs.png"),           x = 0.736, y = 0.088, width = 0.097, height = 0.170) +
  draw_image(stage_img("03_parr.png"),           x = 0.845, y = 0.108, width = 0.103, height = 0.135) +
  draw_label("Spawners",   x = 0.648, y = 0.050, size = 17, colour = "#F1E51D", fontface = "bold", hjust = 0.5) +
  draw_label("Eggs",       x = 0.785, y = 0.050, size = 17, colour = "#F1E51D", fontface = "bold", hjust = 0.5) +
  draw_label("Fry / parr", x = 0.897, y = 0.050, size = 17, colour = INK,      fontface = "bold", hjust = 0.5) +
  draw_image(stage_img("04_smolt.png"), x = 0.200, y = 0.470, width = 0.114, height = 0.090) +
  draw_label("Smolts", x = 0.257, y = 0.454, size = 15, colour = INK, fontface = "bold", hjust = 0.5) +
  draw_image(stage_img("05_ocean_adult.png"), x = 0.038, y = 0.150, width = 0.146, height = 0.110) +
  draw_label("Ocean returns", x = 0.110, y = 0.136, size = 15, colour = INK, fontface = "bold", hjust = 0.5) +
  theme(plot.background = element_rect(fill = SURFACE, colour = NA))

ggsave(file.path(OUT, "08_map_lifecycle.png"), fig8,
       width = FIG_W, height = FIG_H, dpi = 200, bg = SURFACE)
cat("wrote 08_map_lifecycle.png\n")


# ============================================================================
# 9. What the alternatives actually do to the river
# ============================================================================
# The design slide used to assert that the alternatives differ in volume and
# schedule. This shows it. Each panel is one bypass alternative; the filled
# curve is how much cooler Hazel Avenue runs than it would with no bypass, on
# the same day, under the same weather.
#
# Meteorology is averaged across the four met years at equal weight, which is
# exactly how the model combines them, so the panels show the quantity the
# decision is actually made on rather than one arbitrary year.
app <- function(...) here("SalmonCountR", "app_data", ...)

ALT9 <- c("NB", "PB1", "PB2", "PB2b", "PB2c", "PB3", "PB4", "PB5", "PB6")
MET4 <- c(2011, 2014, 2017, 2020)
alt_idx <- tibble(i   = 1:36,
                  alt = ALT9[((0:35) %% 9) + 1],
                  met = MET4[((0:35) %/% 9) + 1])

env_ext <- readRDS(app("env_ext_list.rds"))

# One projected fall. The forecast series is periodic from 2026, so any fall in
# the projection carries each met year's signature; 2030 is simply inside it.
WIN <- as.Date(c("2030-09-20", "2030-12-20"))

temps <- purrr::map_dfr(1:36, function(i) {
  env_ext[[as.character(i)]] %>%
    filter(site == "AveHazel", Date >= WIN[1], Date <= WIN[2]) %>%
    transmute(Date, temp, i = i)
}) %>% left_join(alt_idx, by = "i")

by_alt <- temps %>%
  group_by(alt, Date) %>%
  summarise(temp = mean(temp), .groups = "drop")

nb_ref <- by_alt %>% filter(alt == "NB") %>% select(Date, nb = temp)

cooling <- by_alt %>%
  filter(alt != "NB") %>%
  left_join(nb_ref, by = "Date") %>%
  mutate(delta = nb - temp,
         alt = factor(alt, levels = c("PB1", "PB2", "PB2b", "PB2c",
                                      "PB3", "PB4", "PB5", "PB6")))

g9 <- ggplot(cooling, aes(Date, delta)) +
  geom_hline(yintercept = 0, colour = AXIS, linewidth = 0.5) +
  # Two steps of the same ramp, not the yellow end: yellow is spent on pointing
  # in this deck, and an eight-panel yellow outline would outrank the map.
  geom_area(fill = pal_dark(3)[1], alpha = 0.85) +
  geom_line(colour = pal_dark(3)[2], linewidth = 0.7) +
  facet_wrap(~ alt, nrow = 2) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  labs(x = NULL, y = "°C cooler than no bypass") +
  theme_talk(base_size = 22, legend = "none") +
  theme(panel.spacing = unit(14, "pt"))

save_talk(g9, "09_alternatives_temperature.png", width = 14, height = 6.6)


# ============================================================================
# 10. Spawn timing -- the component the correction turned back on
# ============================================================================
# An ordinal (cumulative-link) model on October and November temperature gives
# the distribution of spawn dates. Warmer Octobers push spawning later, which
# moves eggs deeper into the warm tail. Holding November at its mean and moving
# October one standard deviation either side of its mean isolates the effect.
st <- readRDS(app("spawn_timing_model.rds"))
sd_ <- st$standardisation

bin_probs <- function(oct_C, nov_C) {
  eta <- unname(st$beta["Oct_std"]) * (oct_C - sd_$oct_mean) / sd_$oct_sd +
         unname(st$beta["Nov_std"]) * (nov_C - sd_$nov_mean) / sd_$nov_sd
  cum <- plogis(st$zeta - eta)              # P(bin <= k)
  diff(c(0, cum, 1))                        # per-bin probability
}

bin_start <- as.Date(st$bin_defs$start[seq_along(st$present_bins)])
oct_lo <- sd_$oct_mean - sd_$oct_sd; oct_hi <- sd_$oct_mean + sd_$oct_sd
nov_lo <- sd_$nov_mean - sd_$nov_sd; nov_hi <- sd_$nov_mean + sd_$nov_sd

# ---- No uncertainty band, deliberately --------------------------------------
# An earlier version drew 90% intervals from the CLM's coefficient covariance.
# They are gone, for two reasons.
#
# 1. The model downstream is deterministic. Coefficient uncertainty from this
#    component is not propagated into spawn dates, egg survival, or the adult
#    index. Putting a band on one link while every other link is a point
#    estimate advertises a rigour the projection does not have.
# 2. Those intervals were too narrow anyway. October and November temperature
#    are YEAR-level covariates, so all ~10,000 redds within a brood year share
#    one value: there are 14 independent temperature observations, not 10,000,
#    and a redd-level CLM understates the standard errors accordingly.
#
# The covariance is kept in output/spawn_timing_uncertainty.rds if this is ever
# done properly - which means a year-level random effect, and propagation
# through the whole chain rather than a ribbon on one figure.
band <- function(oct_C, nov_C) {
  tibble(date = bin_start, p = bin_probs(oct_C, nov_C))
}

# Each panel moves one month one standard deviation either side of its own mean
# and holds the other at its mean, so the two effects are not confounded.
# "Cool" and "Warm" mean nothing on their own, so each panel's strip carries the
# two actual temperatures. They also make the asymmetry visible: October swings
# 1.6 °C between cool and warm, November only 0.7 °C.
lab_oct <- sprintf("October
cool %.1f °C  ·  warm %.1f °C", oct_lo, oct_hi)
lab_nov <- sprintf("November
cool %.1f °C  ·  warm %.1f °C", nov_lo, nov_hi)

spawn_dist <- bind_rows(
  band(oct_lo,        sd_$nov_mean) %>% mutate(panel = lab_oct, lvl = "Cool"),
  band(oct_hi,        sd_$nov_mean) %>% mutate(panel = lab_oct, lvl = "Warm"),
  band(sd_$oct_mean,  nov_lo)       %>% mutate(panel = lab_nov, lvl = "Cool"),
  band(sd_$oct_mean,  nov_hi)       %>% mutate(panel = lab_nov, lvl = "Warm")
) %>%
  mutate(panel = factor(panel, levels = c(lab_oct, lab_nov)))

# Median spawn date under each scenario, so the shift is given in days rather
# than left to the eye to judge. Bins are 10 days wide; the median is
# interpolated inside whichever bin the cumulative distribution crosses 0.5.
med_date <- function(p) {
  cum <- cumsum(p); k <- which(cum >= 0.5)[1]
  lo  <- if (k == 1) 0 else cum[k - 1]
  as.numeric(bin_start[k]) + st$bin_width * (0.5 - lo) / p[k]
}
shift_oct <- med_date(bin_probs(oct_hi, sd_$nov_mean)) -
             med_date(bin_probs(oct_lo, sd_$nov_mean))
shift_nov <- med_date(bin_probs(sd_$oct_mean, nov_hi)) -
             med_date(bin_probs(sd_$oct_mean, nov_lo))

# The two months pull in opposite directions: Oct_std is +0.186 and Nov_std is
# -0.070, so a warm October delays spawning and a warm November advances it,
# with October about two and a half times the stronger of the two.
shift_lab <- tibble(
  panel = factor(c(lab_oct, lab_nov), levels = c(lab_oct, lab_nov)),
  txt   = sprintf("warm %s: median %+.1f days", c("Oct", "Nov"),
                  c(shift_oct, shift_nov)))

g10 <- ggplot(spawn_dist, aes(date, p, colour = lvl)) +
  geom_line(linewidth = 1.9) +
  geom_point(size = 3.2) +
  # Top left. Both distributions start at zero on 5 October and the peak sits
  # right of centre, so the left shoulder is the clear corner in both panels.
  # The y-scale carries extra headroom above the peak so this cannot collide
  # with the warm curve however the numbers move.
  geom_text(data = shift_lab, aes(x = min(bin_start), y = Inf, label = txt),
            inherit.aes = FALSE, vjust = 1.5, hjust = 0, colour = INK_SOFT,
            size = 6.4, fontface = "italic") +
  facet_wrap(~ panel) +
  # Warm keeps the deck's temperature yellow; cool takes the ramp's dark end.
  scale_colour_manual(values = c(Cool = unname(pal_dark(2)[1]),
                                 Warm = TEMP_WARM)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%d %b") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.02, 0.20))) +
  # No methodology subtitle: the standard-deviation framing is something the
  # presenter says, and dropping it gives the panels the height instead.
  labs(x = NULL, y = "Share of redds") +
  # Legend to the right rather than on top: two series need very little width,
  # and moving it off the top hands the whole slide height to the panels.
  theme_talk(base_size = 24, legend = "right") +
  theme(panel.spacing = unit(22, "pt"),
        strip.text = element_text(size = 18, colour = INK, face = "bold", lineheight = 1.15))

save_talk(g10, "10_spawn_timing.png", width = 15.5, height = 7.6)


# ---- Parameters for the two structural figures ------------------------------
# Read from the model's own base parameter set rather than retyped, so these
# figures cannot drift away from what the model actually runs. The pre-spawn
# defaults are restated here instead of sourcing functions.R, which would drag
# in the whole simulation stack for one two-parameter logistic.
P_BASE     <- readRDS(app("base_P.rds"))
P_BASE_S0  <- P_BASE$S0            # 0.347
P_BASE_K   <- P_BASE$K_spawners    # 33,185 redds
surv_adult_prespawn <- function(dd, intercept = 3.0, beta = -0.00067) {
  plogis(intercept + beta * dd)
}

# ============================================================================
# 11. Pre-spawn survival
# ============================================================================
# Logistic in accumulated pre-spawn degree-days. The parameters are Colvin et
# al. (2018), not fitted here, so the honest thing to show alongside the curve
# is how little of it the American River actually occupies: the observed range
# sits high on the curve where it is nearly flat.
dd_obs <- readRDS(app("calib_pred_by_variant.rds"))[[1]]$deg_day
dd_grid <- seq(0, 4000, by = 10)

g12 <- ggplot(tibble(dd = dd_grid, s = surv_adult_prespawn(dd_grid)),
              aes(dd, s)) +
  annotate("rect", xmin = min(dd_obs), xmax = max(dd_obs),
           ymin = -Inf, ymax = Inf, fill = "#FFFFFF", alpha = 0.07) +
  annotate("text", x = mean(range(dd_obs)), y = 0.06,
           label = sprintf("observed 2011-2024\n%.0f-%.0f °C·days",
                           min(dd_obs), max(dd_obs)),
           colour = INK, size = 6, lineheight = 0.95) +
  geom_line(colour = pal_dark(3)[2], linewidth = 1.8) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(x = "Accumulated pre-spawn degree-days (°C·days)",
       y = "Pre-spawn survival",
       subtitle = "logit S = 3.0 − 0.00067 · degree-days   (Colvin et al. 2018, not refitted here)") +
  theme_talk(base_size = 24, legend = "none") +
  theme(plot.subtitle = element_text(size = 15, colour = INK_SOFT))

save_talk(g12, "12_prespawn.png", width = 14, height = 6.4)


# ============================================================================
# 12. Density dependence, and where its capacity comes from
# ============================================================================
# Two panels because the Beverton-Holt story is two steps: flow sets spawning
# capacity through weighted usable area, and capacity then sets how hard the
# density-dependent brake bites. The second panel is the reason a large fry
# gain does not become a proportional adult gain.
instream <- readRDS(app("american_river_instream.rds")) %>%
  mutate(K_spawners = FR_spawn_wua / 9.29)

g13a <- ggplot(instream, aes(flow_cfs, K_spawners)) +
  geom_line(colour = pal_dark(3)[1], linewidth = 1.8) +
  geom_point(colour = pal_dark(3)[1], size = 2.6) +
  geom_hline(yintercept = P_BASE_K, colour = AXIS, linetype = "22") +
  annotate("text", x = max(instream$flow_cfs), y = P_BASE_K,
           label = sprintf("base K = %s", scales::comma(round(P_BASE_K))),
           colour = INK_SOFT, size = 5.6, hjust = 1, vjust = -0.7) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Flow (cfs)", y = "Spawning capacity K (redds)",
       subtitle = "Weighted usable area ÷ 9.29 m² per redd") +
  theme_talk(base_size = 21, legend = "none") +
  theme(plot.subtitle = element_text(size = 14, colour = INK_SOFT))

bh <- expand_grid(redds = seq(0, 60000, by = 500),
                  K = c(20000, P_BASE_K, 45000)) %>%
  mutate(s = P_BASE_S0 / (1 + redds / K),
         K_lab = factor(sprintf("K = %s", scales::comma(round(K)))))

g13b <- ggplot(bh, aes(redds, s, colour = K_lab)) +
  geom_line(linewidth = 1.7) +
  scale_colour_manual(values = unname(pal_dark(3))) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  # NOT "egg-to-fry survival": this term multiplies the temperature-dependent
  # egg survival rather than replacing it (fry = eggs x S_temp x S_dd), so
  # naming it egg-to-fry would claim the whole step for one of its two factors.
  labs(x = "Redds", y = "Density-dependent survival",
       subtitle = sprintf("S = %.3f / (1 + redds / K), applied on top of temperature survival",
                          P_BASE_S0)) +
  theme_talk(base_size = 21) +
  theme(plot.subtitle = element_text(size = 14, colour = INK_SOFT))

# patchwork composes onto its own background, which defaults to theme_grey's
# white and shows as a frame around both panels and a seam between them.
g13 <- (g13a + g13b) +
  patchwork::plot_annotation(
    theme = theme(plot.background = element_rect(fill = SURFACE, colour = NA),
                  panel.background = element_rect(fill = SURFACE, colour = NA)))

save_talk(g13, "13_density_dependence.png", width = 15, height = 6.8)


# ============================================================================
# 13. Calibration against observed escapement
# ============================================================================
# Observed GrandTab escapement against the model's prediction, TDM-weighted.
# 2011-2013 seed the run and are therefore reproduced by construction; only
# 2014 onward is a test. The fit over that window is weak, and the figure is
# drawn so that it says so rather than hiding it -- see the speaker notes on
# the slide, and output/calibration_fit_statistics.csv for the numbers.
calib <- read_csv(here("output", "calibration_predictions.csv"),
                  show_col_types = FALSE) %>%
  filter(variant == "TDM-weighted")

calib_long <- calib %>%
  select(year, Observed = observed, Predicted = predicted) %>%
  pivot_longer(-year, names_to = "series", values_to = "n")

seed_band <- data.frame(xmin = 2010.5, xmax = 2013.5, ymin = -Inf, ymax = Inf)

g11 <- ggplot() +
  geom_rect(data = seed_band,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "#FFFFFF", alpha = 0.05) +
  annotate("text", x = 2012, y = -Inf, label = "seeded, not a test",
           colour = INK_SOFT, size = 6, vjust = -1.1, fontface = "italic") +
  geom_line(data = calib_long, aes(year, n, colour = series), linewidth = 1.6) +
  geom_point(data = calib_long, aes(year, n, colour = series), size = 3.4) +
  scale_colour_manual(values = c(Observed = INK, Predicted = pal_dark(3)[2])) +
  scale_x_continuous(breaks = seq(2011, 2024, 2)) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Spawner escapement") +
  theme_talk(base_size = 24)

save_talk(g11, "11_calibration.png", width = 14, height = 6.4)

cat("\nAll BDSC figures written to figures/bdsc/\n")
