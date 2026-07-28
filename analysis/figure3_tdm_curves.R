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
#    cycle model actually consumes. The exponential models overtake Martin at
#    16.6-17.3 deg C, which the daily-hazard panel alone does not show.
#  * The Lower American River October-November operational range is shaded,
#    computed from the CE-QUAL-W2 scenario temperatures at Hazel Avenue.
#
# Outputs: figures/figure3_tdm_curves.png
#          output/figure3_etf_survival_by_temp.csv
# ============================================================================

library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
library(scales); library(here); library(lubridate)

source(here("SalmonCountR", "functions.R"))

# ---- 1. Operational temperature range, Oct-Nov at Hazel Avenue --------------
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
  mutate(S_day = pmin(pmax(S_day, 0), 1))

# purple / teal / vermillion — CVD-safe trio, all >=3:1 on white
fam_cols <- c("Bratovich et al. (2020)"    = "#440154",
              "Bartholow & Heasley (2006)" = "#1F9E89",
              "Martin et al. (2017)"       = "#D55E00")
stage_types <- c("Egg" = "solid", "Alevin" = "dashed", "Incubation" = "dotted")

# Build a fresh layer per panel — a single layer object cannot be shared
# between two ggplot builds. Clamp to the plotted x range: a rect whose xmax
# falls outside scale_x_continuous(limits = ...) is censored and dropped whole.
X_LO <- 10; X_HI <- 18
shade <- function() {
  annotate("rect",
           xmin = max(as.numeric(op_range[1]), X_LO),
           xmax = min(as.numeric(op_range[2]), X_HI),
           ymin = -Inf, ymax = Inf, fill = "#B0AFA8", alpha = 0.32)
}

pA <- ggplot(curves_daily, aes(T, S_day, colour = family, linetype = stage)) +
  shade() +
  annotate("text", x = mean(as.numeric(op_range)), y = 0.30,
           label = sprintf("Lower American River\nOct-Nov range\n(%.1f-%.1f °C)",
                           op_range[1], op_range[2]),
           size = 3, colour = "grey20", lineheight = 0.95) +
  geom_vline(xintercept = 12.14, linetype = "dotdash",
             colour = "grey35", linewidth = 0.4) +
  annotate("text", x = 12.14, y = 0.06, label = "Martin threshold 12.14",
           hjust = -0.05, size = 3, colour = "grey25") +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = fam_cols, name = "TDM model") +
  scale_linetype_manual(values = stage_types, name = "Stage") +
  scale_x_continuous(limits = c(10, 18), breaks = seq(10, 18, 1), expand = c(0, 0)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(subtitle = "(a) Daily survival rate", x = NULL, y = "Daily survival") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.border  = element_rect(colour = "grey40", fill = NA, linewidth = 0.5),
        plot.subtitle = element_text(face = "bold"),
        legend.position = "right")

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
  pivot_longer(-T, names_to = "family", values_to = "S")

# Crossovers: where each exponential model becomes harsher than Martin
cross <- sapply(c("WaterForum2020", "SALMOD2006"), function(m) {
  f  <- function(x) S_cum(x, m) - S_cum(x, "martin")
  gr <- seq(13, 20, by = 0.001); v <- sapply(gr, f)
  i  <- which(diff(sign(v)) != 0)
  if (length(i)) gr[i[length(i)]] else NA_real_
})
cat(sprintf("\nCrossover (exponential model drops below Martin): Bratovich %.2f C, Bartholow %.2f C\n",
            cross[["WaterForum2020"]], cross[["SALMOD2006"]]))

pB <- ggplot(curves_cum, aes(T, S, colour = family)) +
  shade() +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = as.numeric(cross), linetype = "dotted",
             colour = "grey30", linewidth = 0.45) +
  annotate("text", x = cross[["SALMOD2006"]], y = 0.62,
           label = sprintf("exponential models fall\nbelow Martin: %.1f / %.1f °C",
                           cross[["SALMOD2006"]], cross[["WaterForum2020"]]),
           hjust = 1.06, size = 3, colour = "grey20", lineheight = 0.95) +
  scale_colour_manual(values = fam_cols, name = "TDM model") +
  scale_x_continuous(limits = c(10, 18), breaks = seq(10, 18, 1), expand = c(0, 0)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(subtitle = "(b) Cumulative egg-to-fry survival at constant temperature",
       x = "Temperature (°C)", y = "Egg-to-fry survival") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.border  = element_rect(colour = "grey40", fill = NA, linewidth = 0.5),
        plot.subtitle = element_text(face = "bold"),
        legend.position = "right")

fig <- pA / pB + plot_layout(guides = "collect") &
  theme(legend.position = "right")

dir.create(here("figures"), showWarnings = FALSE)
dir.create(here("output"),  showWarnings = FALSE)
ggsave(here("figures", "figure3_tdm_curves.png"), fig,
       width = 9.5, height = 8, dpi = 300, bg = "white")

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
