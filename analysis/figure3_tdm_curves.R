# ============================================================================
# Figure 3. TDM model comparison — daily hazard and cumulative egg-to-fry survival
# ============================================================================
# Revision notes:
#  * Martin threshold is 12.14 deg C (Martin et al. 2017), matching
#    SalmonCountR/functions.R::tdm_lin_martin. The old caption's 12.8 was wrong.
#  * x range extended to 10-18 deg C. The previous version topped out near
#    15 deg C, where the exponential models still look almost flat — the most
#    plausible reason Reviewer 2 read TDM.1 as having no temperature response.
#  * Panel (b) added: cumulative egg-to-fry survival, which is what the life
#    cycle model actually consumes.
#  * The crossover rules at 16.6/17.3 deg C were dropped (2026-09-03). They were
#    arithmetically right but unreadable: every curve is under 0.5% survival
#    there, so the crossing is a sub-pixel event and one label was describing
#    two different crossings. See output/figure3_etf_survival_by_temp.csv if the
#    numbers are ever wanted.
#  * No annotation layers (2026-09-03). The shaded operational range, the
#    12.14 deg C threshold rule and their labels were all removed: the curves
#    carry the comparison on their own, and the furniture was competing with
#    them. The Hazel Avenue range is still computed and printed below, because
#    the manuscript text quotes it — it is just no longer drawn on the figure.
#
# Outputs: figures/figure3_tdm_curves.png
#          output/figure3_etf_survival_by_temp.csv
# ============================================================================

library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
library(scales); library(here); library(lubridate)

source(here("SalmonCountR", "functions.R"))
source(here("analysis", "figure_theme.R"))

# ---- 1. Operational temperature range, Oct-Nov at Hazel Avenue --------------
# Reported to the console only. Nothing here is drawn on the figure any more;
# it is kept so the range quoted in the text stays reproducible from this script.
env_ext_list <- readRDS(here("SalmonCountR", "app_data", "env_ext_list.rds"))

oct_nov <- bind_rows(lapply(names(env_ext_list), function(nm) {
  d <- env_ext_list[[nm]]
  if (is.null(d) || !nrow(d)) return(NULL)
  d %>% mutate(env = nm)
})) %>%
  filter(site == "AveHazel", month(Date) %in% c(10, 11)) %>%
  filter(is.finite(temp))

op_range <- quantile(oct_nov$temp, c(0.05, 0.95), na.rm = TRUE)
cat(sprintf("Hazel Ave Oct-Nov temperatures across all alternatives (n = %d):\n",
            nrow(oct_nov)))
cat(sprintf("  min %.2f  5%% %.2f  median %.2f  95%% %.2f  max %.2f  (deg C)\n",
            min(oct_nov$temp), op_range[1], median(oct_nov$temp),
            op_range[2], max(oct_nov$temp)))

# ---- 2. Panel (a): daily survival ------------------------------------------
s_day_exp    <- function(T, alpha, beta) exp(-alpha * exp(beta * T))
s_day_martin <- function(T, alpha = 0.026, beta = 12.14) exp(-alpha * pmax(T - beta, 0))

WF_egg    <- list(alpha = 3.408488e-11, beta = 1.21122)
WF_alevin <- list(alpha = 1.017554e-10, beta = 1.24092)
SM_egg    <- list(alpha = 1.475e-11,    beta = 1.392)
SM_alevin <- list(alpha = 2.521e-12,    beta = 1.461)

T_seq <- seq(10, 18, by = 0.02)

# Viridis, matching every other figure in the manuscript (see analysis/figure_theme.R)
fam_cols <- TDM_COLS
stage_labels <- c("Egg" = "Egg", "Alevin" = "Alevin",
                  "Incubation" = "Not stage-specific")
stage_types  <- c("Egg" = "solid", "Alevin" = "dashed",
                  "Not stage-specific" = "dotted")

curves_daily <- bind_rows(
  tibble(T = T_seq, family = "Bratovich et al. (2020)",    stage = "Egg",
         S_day = s_day_exp(T_seq, WF_egg$alpha,    WF_egg$beta)),
  tibble(T = T_seq, family = "Bratovich et al. (2020)",    stage = "Alevin",
         S_day = s_day_exp(T_seq, WF_alevin$alpha, WF_alevin$beta)),
  tibble(T = T_seq, family = "Bartholow & Heasley (2006)", stage = "Egg",
         S_day = s_day_exp(T_seq, SM_egg$alpha,    SM_egg$beta)),
  tibble(T = T_seq, family = "Bartholow & Heasley (2006)", stage = "Alevin",
         S_day = s_day_exp(T_seq, SM_alevin$alpha, SM_alevin$beta)),
  tibble(T = T_seq, family = "Martin et al. (2017)",       stage = "Incubation",
         S_day = s_day_martin(T_seq))
) %>%
  mutate(S_day = pmin(pmax(S_day, 0), 1),
         # Order the legend TDM.1, TDM.2, TDM.3 as the text numbers them,
         # rather than letting ggplot sort it alphabetically.
         family = factor(family, levels = names(TDM_COLS)),
         # Developmental order, not alphabetical. Martin is not stage-specific,
         # so its single curve is named as such rather than sitting in the list
         # as if it were a third stage the other two models also have.
         stage = factor(stage, levels = names(stage_labels),
                        labels = unname(stage_labels)))

pA <- ggplot(curves_daily, aes(T, S_day, colour = family, linetype = stage)) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = fam_cols, name = "TDM model") +
  scale_linetype_manual(values = stage_types, name = "Parameter set") +
  # Grey the linetype keys so they do not read as a fourth coloured series,
  # and draw the colour keys solid so they match the curves they name.
  guides(colour   = guide_legend(order = 1, override.aes = list(linetype = "solid")),
         linetype = guide_legend(order = 2, override.aes = list(colour = "grey30"))) +
  scale_x_continuous(limits = c(10, 18), breaks = seq(10, 18, 1), expand = c(0, 0)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(subtitle = "(a) Daily survival rate", x = NULL, y = "Daily survival") +
  theme_arg(base_size = 14)

# ---- 3. Panel (b): cumulative egg-to-fry survival ---------------------------
# Uses the shipped implementation so the panel matches what the life cycle runs:
# ATU-paced stage boundaries (hatch 400, emergence 958) at constant temperature.
S_cum <- function(Temp, model) {
  temps <- rep(Temp, 400)
  if (model == "martin") tdm_lin_martin(temps[.slice_by_atu(temps)])
  else tdm_exp(temps, calib = model, use_stages = TRUE)
}

curves_cum <- tibble(T = T_seq) %>%
  mutate(`Bratovich et al. (2020)`    = vapply(T, S_cum, 0, "WaterForum2020"),
         `Bartholow & Heasley (2006)` = vapply(T, S_cum, 0, "SALMOD2006"),
         `Martin et al. (2017)`       = vapply(T, S_cum, 0, "martin")) %>%
  pivot_longer(-T, names_to = "family", values_to = "S") %>%
  mutate(family = factor(family, levels = names(TDM_COLS)))

pB <- ggplot(curves_cum, aes(T, S, colour = family)) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = fam_cols, name = "TDM model") +
  scale_x_continuous(limits = c(10, 18), breaks = seq(10, 18, 1), expand = c(0, 0)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(subtitle = "(b) Cumulative egg-to-fry survival at constant temperature",
       x = "Temperature (°C)", y = "Egg-to-fry survival") +
  theme_arg(base_size = 14)

# Pin the collected guides to the top. Left centred, patchwork splits them
# across the full figure height and the linetype guide ends up alongside panel
# (b), which has no linetypes at all — it reads as a legend for the wrong panel.
fig <- pA / pB + plot_layout(guides = "collect") &
  theme(legend.position = "right",
        legend.justification = "top",
        legend.box.just = "left",
        legend.spacing.y = unit(10, "pt"))

dir.create(here("figures"), showWarnings = FALSE)
dir.create(here("output"),  showWarnings = FALSE)
ggsave(here("figures", "figure3_tdm_curves.png"), fig,
       width = 11, height = 9, dpi = 300, bg = "white")

# ---- 4. Companion table -----------------------------------------------------
tbl <- tibble(T_C = seq(10, 18, 0.5)) %>%
  mutate(Bratovich_pct = 100 * vapply(T_C, S_cum, 0, "WaterForum2020"),
         Bartholow_pct = 100 * vapply(T_C, S_cum, 0, "SALMOD2006"),
         Martin_pct    = 100 * vapply(T_C, S_cum, 0, "martin"))
write.csv(tbl %>% mutate(across(-T_C, ~round(., 2))),
          here("output", "figure3_etf_survival_by_temp.csv"), row.names = FALSE)

cat("\nEgg-to-fry survival (%) at constant temperature:\n")
print(as.data.frame(tbl %>% mutate(across(-T_C, ~round(., 1)))), row.names = FALSE)
cat("\nWrote figures/figure3_tdm_curves.png and output/figure3_etf_survival_by_temp.csv\n")
