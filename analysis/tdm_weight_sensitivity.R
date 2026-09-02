# ============================================================================
# TDM weight sensitivity of the composite MCDA score
# ============================================================================
# Responds to Reviewer 2's request for a robustness-under-structural-
# uncertainty analysis: how much does the TDM panel weighting drive the
# preferred alternative?
#
# METHOD FROM B. MAHARDJA (tdm_weight_sensitivity_BM.R, 2026-09-02), reworked
# here to read the committed app_data instead of a hand-built spreadsheet, so
# the analysis regenerates from the model rather than from a transcription.
# Verified to reproduce his numbers: composites agree to the fourth decimal and
# the Martin-sweep boundaries to ~0.01, the residual being integer rounding in
# ConsequenceTable.xlsx (e.g. NB under Bratovich 12,838 there, 12,827.5 here).
# His spreadsheet is kept at data_raw/ConsequenceTable.xlsx for provenance.
#
# WHY THE NORMALISATION IS GLOBAL. An earlier version re-derived the min-max
# bounds inside every weighting. That is invalid, and for a sharper reason than
# "it rescales things": swing weights are *defined relative to the range of each
# objective*. Move the range and the elicited weights no longer describe the
# trade-off that was elicited - you would have to redo the swing weighting at
# every point of the sweep. Bounds are therefore computed once, across all nine
# alternatives AND all three TDM models, before any weighting is applied. This
# matches the pooled convention in analysis/evpi.R.
#
# OBJECTIVE WEIGHTS are the post-hoc VoI set (0.73 salmon / 0.22 hydropower /
# 0.05 steelhead), not the elicited set. Those weights were reverse-engineered
# to reproduce the observed ranking, which is exactly why this analysis and the
# VoI analysis belong in the response to reviewers rather than the main text.
#
# Holds fixed:
#   objective weights   0.73 salmon / 0.22 hydropower / 0.05 steelhead
#   normalisation scale global, across alternatives x TDM models
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

source(here("analysis", "figure_theme.R"))

# ---- 1. Inputs --------------------------------------------------------------
results_full               <- readRDS(here("SalmonCountR","app_data","results_full.rds"))
steelhead_scenario_results <- readRDS(here("SalmonCountR","app_data","steelhead_scenario_results.rds"))

real_years  <- 2011:2024
scenarios   <- c("NB","PB1","PB2","PB2b","PB2c","PB3","PB4","PB5","PB6")
scen_map    <- c(NB=1, PB1=2, PB2=3, PB2b=4, PB2c=5, PB3=6, PB4=7, PB5=8, PB6=9)
wyt_default <- 0.25

# Hydropower revenue loss ($) per alternative - same source as mcda.R
hydro_loss <- c(NB=0, PB1=111422, PB2=370826, PB2b=470090, PB2c=433215,
                PB3=201552, PB4=241590, PB5=199382, PB6=348806)

W_SALMON <- 0.73
W_HYDRO  <- 0.22
W_STEEL  <- 0.05

# ---- 2. Consequence table: alternative x TDM model --------------------------
# Adult index, median of the final 20 forecast years, averaged over the four
# meteorological years. env 1-9 = 2011 met year, 10-18 = 2014, 19-27 = 2017,
# 28-36 = 2020; within each block 1=NB ... 9=PB6.
salmon <- results_full %>%
  filter(year > max(real_years)) %>%
  group_by(env, variant) %>%
  slice_tail(n = 20) %>%
  summarise(med = median(spawners, na.rm = TRUE), .groups = "drop") %>%
  mutate(env      = as.integer(env),
         scenario = names(scen_map)[match((env - 1) %% 9 + 1, scen_map)]) %>%
  group_by(scenario, variant) %>%
  summarise(value = sum(med * wyt_default), .groups = "drop")

# ---- 3. Global normalisation, before any weighting --------------------------
SALMON_LO <- min(salmon$value); SALMON_HI <- max(salmon$value)

salmon_n <- salmon %>% mutate(n = (value - SALMON_LO) / (SALMON_HI - SALMON_LO))

steel_n <- steelhead_scenario_results %>%
  transmute(scenario,
            n = (steelhead_score - min(steelhead_score)) /
                (max(steelhead_score) - min(steelhead_score)))

hydro_n <- tibble(scenario = names(hydro_loss), raw = as.numeric(hydro_loss)) %>%
  mutate(n = (max(raw) - raw) / (max(raw) - min(raw)))   # lower cost is better

# ---- 4. Composite score -----------------------------------------------------
# Salmon enters as the TDM-weighted mean of its NORMALISED values. Because the
# scale is fixed, that is identical to normalising the weighted mean - the two
# only diverge when the bounds move, which is the defect this version removes.
composite <- function(w, w_salmon = W_SALMON, w_hydro = W_HYDRO, w_steel = W_STEEL) {
  salmon_n %>%
    mutate(wt = w[variant]) %>%
    group_by(scenario) %>%
    summarise(chinook_n = sum(n * wt), .groups = "drop") %>%
    left_join(steel_n %>% rename(steel_n = n), by = "scenario") %>%
    left_join(hydro_n %>% select(scenario, hydro_n = n, hydro_raw = raw), by = "scenario") %>%
    mutate(score = w_salmon * chinook_n + w_steel * steel_n + w_hydro * hydro_n)
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

# ---- 5. Continuous Martin-weight sweep --------------------------------------
# Martin weight from 0 to 1; the other two hold the elicited 0.51 : 0.24 ratio.
ratio <- c(exp_WF = .51/.75, exp_SM = .24/.75)
sweep <- map_dfr(seq(0, 1, by = 0.002), function(m) {
  w <- c(exp_WF = ratio[["exp_WF"]]*(1-m), exp_SM = ratio[["exp_SM"]]*(1-m), lin_Martin = m)
  composite(w) %>% mutate(w_martin = m)
})

sweep_top <- sweep %>% group_by(w_martin) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>% ungroup()

runs <- rle(as.character(sweep_top$scenario))
invariance <- tibble(
  top    = runs$values,
  from   = sweep_top$w_martin[cumsum(c(1, head(runs$lengths, -1)))],
  to     = sweep_top$w_martin[cumsum(runs$lengths)]
)

# ---- 6. Console report ------------------------------------------------------
cat(sprintf("\nGlobal salmon scale (all alternatives x all TDM models): %.1f to %.1f\n",
            SALMON_LO, SALMON_HI))
cat(sprintf("Objective weights: %.2f salmon / %.2f hydropower / %.2f steelhead\n",
            W_SALMON, W_HYDRO, W_STEEL))

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

# ---- 7. Figure --------------------------------------------------------------
pA <- ggplot(all_scores, aes(x = scenario, y = score)) +
  geom_col(aes(fill = is_top), colour = "grey25", linewidth = 0.3, width = 0.75) +
  geom_text(aes(label = sprintf("%.3f", score)),
            hjust = -0.10, size = 3.1, fontface = "bold",
            colour = "black", angle = 90) +
  facet_wrap(~ weighting, nrow = 1) +
  scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "grey85"),
                    breaks = "TRUE", labels = "Top-ranked", name = NULL) +
  scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(subtitle = "(a) Composite score under alternative TDM model weightings",
       x = NULL, y = "Composite score") +
  theme_arg(base_size = 14, legend = "top") +
  theme(panel.grid.major.x = element_blank(),
        strip.text  = element_text(face = "bold", size = 11, lineheight = 1.05,
                                   colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                   face = "bold", size = 10, colour = "black"))

focal <- sweep %>% group_by(w_martin) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>% ungroup() %>%
  distinct(scenario) %>% pull(scenario) %>% as.character()

focal_cols  <- setNames(arg_pal(length(focal), begin = 0, end = 0.70), focal)
focal_types <- setNames(rep(c("solid", "longdash", "dotted", "dotdash"),
                            length.out = length(focal)), focal)

sweep_plot <- sweep %>%
  mutate(grp = ifelse(as.character(scenario) %in% focal, as.character(scenario), "Other"),
         grp = factor(grp, levels = c(focal, "Other")))

ends <- sweep_plot %>% filter(w_martin == max(w_martin), grp != "Other")

y_hi <- max(sweep$score); y_lo <- min(sweep$score); y_pad <- 0.06 * (y_hi - y_lo)

pB <- ggplot() +
  geom_line(data = filter(sweep_plot, grp == "Other"),
            aes(w_martin, score, group = scenario),
            colour = GREY_CTX, linewidth = 0.6) +
  geom_line(data = filter(sweep_plot, grp != "Other"),
            aes(w_martin, score, colour = grp, linetype = grp), linewidth = 1.2) +
  geom_vline(xintercept = 0.25, linetype = "dashed", colour = "grey25", linewidth = 0.5) +
  annotate("text", x = 0.25, y = y_hi - y_pad, label = "elicited\nMartin weight = 0.25",
           hjust = -0.08, size = 3.8, fontface = "bold", colour = "black",
           lineheight = 0.95) +
  annotate("text", x = 0.02, y = y_lo + y_pad,
           label = "scale fixed across all TDM models; objective weights 0.73 / 0.22 / 0.05",
           hjust = 0, size = 3.4, fontface = "italic", colour = "grey25") +
  ggrepel::geom_text_repel(data = ends, aes(w_martin, score, label = scenario, colour = grp),
                           direction = "y", hjust = 0, nudge_x = 0.02,
                           segment.size = 0.2, size = 4, fontface = "bold",
                           show.legend = FALSE) +
  scale_colour_manual(values = focal_cols, name = NULL) +
  scale_linetype_manual(values = focal_types, name = NULL) +
  scale_x_continuous(limits = c(0, 1.09), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(subtitle = "(b) Composite score across the Martin weight (Bratovich : Bartholow held at 0.51 : 0.24)",
       x = "Weight on Martin et al. (2017) TDM model", y = "Composite score") +
  theme_arg(base_size = 14, legend = "top")

fig <- pA / pB + plot_layout(heights = c(1, 0.95))

dir.create(here("figures"), showWarnings = FALSE)
dir.create(here("output"),  showWarnings = FALSE)
ggsave(here("figures","tdm_weight_sensitivity.png"), fig,
       width = 13, height = 10, dpi = 300, bg = "white")

# ---- 8. Tables --------------------------------------------------------------
write.csv(
  all_scores %>%
    mutate(weighting = gsub("\n", " ", weighting)) %>%
    select(weighting, scenario, chinook_norm = chinook_n, steelhead_norm = steel_n,
           hydro_revenue_loss = hydro_raw, hydro_norm = hydro_n,
           composite = score, rank),
  here("output","tdm_weight_sensitivity.csv"), row.names = FALSE)

write.csv(sweep_top %>% select(w_martin, top_alternative = scenario, top_score = score),
          here("output","tdm_weight_sensitivity_martin_sweep.csv"), row.names = FALSE)

cat("\nWrote figures/tdm_weight_sensitivity.png,",
    "output/tdm_weight_sensitivity.csv,",
    "output/tdm_weight_sensitivity_martin_sweep.csv\n")
