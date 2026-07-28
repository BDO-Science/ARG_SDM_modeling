# ============================================================================
# TDM weight sensitivity of the composite MCDA score
# ============================================================================
# Responds to Reviewer 2's request for a robustness-under-structural-
# uncertainty analysis: how much does the TDM panel weighting drive the
# preferred alternative?
#
# Holds fixed:
#   objective weights   0.40 Chinook / 0.50 hydropower / 0.10 steelhead
#   meteorological year weights   0.25 each (2011, 2014, 2017, 2020)
# Varies:
#   TDM model weights across the five weightings below, plus a continuous
#   sweep of the Martin weight from 0 to 1.
#
# Inputs  : SalmonCountR/app_data/{results_full,steelhead_scenario_results}.rds
# Outputs : output/tdm_weight_sensitivity.csv
#           output/tdm_weight_sensitivity_martin_sweep.csv
#           figures/tdm_weight_sensitivity.png
# ============================================================================

library(dplyr); library(tidyr); library(purrr); library(tibble)
library(ggplot2); library(patchwork); library(here)

# ---- 1. Inputs --------------------------------------------------------------
results_full               <- readRDS(here("SalmonCountR","app_data","results_full.rds"))
steelhead_scenario_results <- readRDS(here("SalmonCountR","app_data","steelhead_scenario_results.rds"))

real_years  <- 2011:2024
scenarios   <- c("NB","PB1","PB2","PB2b","PB2c","PB3","PB4","PB5","PB6")
scen_map    <- c(NB=1, PB1=2, PB2=3, PB2b=4, PB2c=5, PB3=6, PB4=7, PB5=8, PB6=9)
hydro_years <- c("2011","2014","2017","2020")
wyt_default <- c("2011"=.25, "2014"=.25, "2017"=.25, "2020"=.25)

# Hydropower revenue loss ($) per alternative — same source as mcda.R
hydro_loss <- c(NB=0, PB1=111422, PB2=370826, PB2b=470090, PB2c=433215,
                PB3=201552, PB4=241590, PB5=199382, PB6=348806)

# ---- 2. Adult index: median of the final 20 forecast years ------------------
# env 1-9 = 2011 met year, 10-18 = 2014, 19-27 = 2017, 28-36 = 2020;
# within each block 1=NB, 2=PB1, 3=PB2, 4=PB2b, 5=PB2c, 6=PB3, 7=PB4, 8=PB5, 9=PB6
med20 <- results_full %>%
  filter(year > max(real_years)) %>%
  group_by(env, variant) %>%
  slice_tail(n = 20) %>%
  summarise(med = median(spawners, na.rm = TRUE), .groups = "drop") %>%
  mutate(env      = as.integer(env),
         scenario = names(scen_map)[match((env - 1) %% 9 + 1, scen_map)],
         wyt      = hydro_years[(env - 1) %/% 9 + 1])

spawner_metric <- function(w) {
  med20 %>%
    mutate(wt = w[variant] * wyt_default[wyt]) %>%
    group_by(scenario) %>%
    summarise(metric = sum(med * wt), .groups = "drop")
}

# ---- 3. Composite score -----------------------------------------------------
# Min-max normalisation across the nine alternatives, matching app.R and mcda.R.
norm_hi <- function(x) (x - min(x)) / (max(x) - min(x))   # higher is better
norm_lo <- function(x) (max(x) - x) / (max(x) - min(x))   # lower  is better

composite <- function(w, w_salmon = .40, w_hydro = .50, w_steel = .10) {
  spawner_metric(w) %>%
    left_join(steelhead_scenario_results, by = "scenario") %>%
    mutate(hydro_raw = hydro_loss[scenario],
           chinook_n = norm_hi(metric),
           steel_n   = norm_hi(steelhead_score),
           hydro_n   = norm_lo(hydro_raw),
           score     = w_salmon*chinook_n + w_steel*steel_n + w_hydro*hydro_n)
}

W <- list(
  "Bratovich only\n(1, 0, 0)"        = c(exp_WF=1,   exp_SM=0,   lin_Martin=0),
  "Bartholow only\n(0, 1, 0)"        = c(exp_WF=0,   exp_SM=1,   lin_Martin=0),
  "Martin only\n(0, 0, 1)"           = c(exp_WF=0,   exp_SM=0,   lin_Martin=1),
  "Equal\n(1/3, 1/3, 1/3)"           = c(exp_WF=1/3, exp_SM=1/3, lin_Martin=1/3),
  "Elicited\n(0.51, 0.24, 0.25)"     = c(exp_WF=.51, exp_SM=.24, lin_Martin=.25)
)

all_scores <- imap_dfr(W, ~ composite(.x) %>% mutate(weighting = .y)) %>%
  group_by(weighting) %>%
  mutate(rank = rank(-score), is_top = rank == 1) %>%
  ungroup() %>%
  mutate(weighting = factor(weighting, levels = names(W)),
         scenario  = factor(scenario, levels = scenarios))

# ---- 4. Continuous Martin-weight sweep --------------------------------------
# Martin weight from 0 to 1; the other two hold the elicited 0.51 : 0.24 ratio.
ratio <- c(exp_WF = .51/.75, exp_SM = .24/.75)
sweep <- map_dfr(seq(0, 1, by = 0.002), function(m) {
  w <- c(exp_WF = ratio[["exp_WF"]]*(1-m), exp_SM = ratio[["exp_SM"]]*(1-m), lin_Martin = m)
  composite(w) %>% mutate(w_martin = m)
})

sweep_top <- sweep %>% group_by(w_martin) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>% ungroup()

runs <- rle(sweep_top$scenario)
invariance <- tibble(
  top    = runs$values,
  from   = sweep_top$w_martin[cumsum(c(1, head(runs$lengths, -1)))],
  to     = sweep_top$w_martin[cumsum(runs$lengths)]
)

# ---- 5. Console report ------------------------------------------------------
cat("\n=== Composite score by TDM weighting ===\n")
print(as.data.frame(
  all_scores %>% select(scenario, weighting, score) %>%
    mutate(weighting = gsub("\n", " ", weighting)) %>%
    pivot_wider(names_from = weighting, values_from = score) %>%
    mutate(across(where(is.numeric), ~round(., 4)))), row.names = FALSE)

cat("\n=== Top-ranked alternative under each weighting ===\n")
print(as.data.frame(
  all_scores %>% filter(is_top) %>%
    select(weighting, scenario, score) %>%
    mutate(weighting = gsub("\n", " ", weighting), score = round(score, 4))), row.names = FALSE)

cat("\n=== Martin-weight invariance intervals ===\n")
print(as.data.frame(invariance), row.names = FALSE)

# ---- 6. Figure --------------------------------------------------------------
# Panel A: composite score by alternative, one facet per TDM weighting.
# Top-ranked alternative in each facet is filled; the rest are outlined.
pA <- ggplot(all_scores, aes(x = scenario, y = score)) +
  geom_col(aes(fill = is_top), colour = "grey25", linewidth = 0.3, width = 0.75) +
  # Labels are set vertically: at nine bars per facet, horizontal labels collide.
  geom_text(aes(label = sprintf("%.3f", score)),
            hjust = -0.12, size = 2.4, colour = "grey15", angle = 90) +
  facet_wrap(~ weighting, nrow = 1) +
  scale_fill_manual(values = c(`TRUE` = "#31688E", `FALSE` = "grey88"),
                    breaks = "TRUE", labels = "Top-ranked", name = NULL) +
  scale_y_continuous(limits = c(0, 0.82), breaks = seq(0, 0.8, 0.2), expand = c(0, 0)) +
  labs(subtitle = "(a) Composite score under alternative TDM model weightings",
       x = NULL, y = "Composite score") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.border  = element_rect(colour = "grey70", fill = NA, linewidth = 0.4),
        strip.text    = element_text(face = "bold", size = 8.5, lineheight = 1.05),
        axis.text.x   = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7.5),
        plot.subtitle = element_text(face = "bold", size = 11),
        legend.position = "top")

# Panel B: continuous Martin-weight sweep. Only the alternatives that ever take
# the top rank somewhere on the sweep are coloured; the rest are grey. Three
# colours keeps every pair well separated for colour-vision deficiency, and the
# lines are direct-labelled so identity never rests on colour alone.
focal <- sweep %>% group_by(w_martin) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>% ungroup() %>%
  distinct(scenario) %>% pull(scenario)
# purple / teal / vermillion — a CVD-safe trio, all >=3:1 on a white surface
focal_cols <- setNames(c("#440154", "#1F9E89", "#D55E00")[seq_along(focal)], focal)

sweep_plot <- sweep %>%
  mutate(grp = ifelse(scenario %in% focal, scenario, "Other"),
         grp = factor(grp, levels = c(focal, "Other")))

ends <- sweep_plot %>% filter(w_martin == max(w_martin), scenario %in% focal)

pB <- ggplot() +
  geom_line(data = filter(sweep_plot, grp == "Other"),
            aes(w_martin, score, group = scenario),
            colour = "grey82", linewidth = 0.5) +
  geom_line(data = filter(sweep_plot, grp != "Other"),
            aes(w_martin, score, colour = grp), linewidth = 0.9) +
  geom_vline(xintercept = 0.25, linetype = "dashed", colour = "grey35", linewidth = 0.4) +
  annotate("text", x = 0.25, y = 0.66, label = "elicited\nMartin weight = 0.25",
           hjust = -0.08, size = 2.8, colour = "grey25", lineheight = 0.95) +
  ggrepel::geom_text_repel(data = ends, aes(w_martin, score, label = scenario, colour = grp),
                           direction = "y", hjust = 0, nudge_x = 0.02,
                           segment.size = 0.2, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = focal_cols, name = NULL) +
  scale_x_continuous(limits = c(0, 1.09), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(subtitle = "(b) Composite score across the Martin weight (Bratovich : Bartholow held at 0.51 : 0.24)",
       x = "Weight on Martin et al. (2017) TDM model", y = "Composite score") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.border  = element_rect(colour = "grey70", fill = NA, linewidth = 0.4),
        plot.subtitle = element_text(face = "bold", size = 11),
        legend.position = "top")

fig <- pA / pB + plot_layout(heights = c(1, 0.95))

dir.create(here("figures"), showWarnings = FALSE)
dir.create(here("output"),  showWarnings = FALSE)
ggsave(here("figures","tdm_weight_sensitivity.png"), fig,
       width = 10, height = 8.5, dpi = 300, bg = "white")

# ---- 7. Tables --------------------------------------------------------------
write.csv(
  all_scores %>%
    mutate(weighting = gsub("\n", " ", weighting)) %>%
    select(weighting, scenario, adult_index = metric, steelhead_score,
           hydro_revenue_loss = hydro_raw, chinook_norm = chinook_n,
           steelhead_norm = steel_n, hydro_norm = hydro_n, composite = score, rank),
  here("output","tdm_weight_sensitivity.csv"), row.names = FALSE)

write.csv(sweep_top %>% select(w_martin, top_alternative = scenario, top_score = score),
          here("output","tdm_weight_sensitivity_martin_sweep.csv"), row.names = FALSE)

cat("\nWrote figures/tdm_weight_sensitivity.png,",
    "output/tdm_weight_sensitivity.csv,",
    "output/tdm_weight_sensitivity_martin_sweep.csv\n")
