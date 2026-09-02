# ============================================================================
# EVPI on the TDM model choice, under both objective weight sets
# ============================================================================
# Reviewer 1 asked for the rationale behind the objective weights used for the
# EVPI calculation. The published 0.026 (4.6%) uses weights re-derived to
# reproduce the observed ranking (0.73 salmon / 0.22 hydropower / 0.05
# steelhead) rather than the elicited weights (0.40 / 0.50 / 0.10). Computing
# it under both lets the response report a range instead of a single figure.
#
# States of nature : the three TDM models, with the elicited panel weights
#                    as their probabilities (0.51 / 0.24 / 0.25).
# EVPI = E_m[ max_a score(a, m) ] - max_a E_m[ score(a, m) ]
#
# Two normalisations are reported, because the choice matters here:
#   "per-state"  min-max within each TDM model (what app.R/mcda.R do for a
#                single weighting) -- convenient, but it rescales each state
#                separately, which inflates the value of information.
#   "pooled"     one min-max over all alternatives AND all states, so scores
#                are on a single comparable scale. This is the defensible one
#                for EVPI.
#
# Output: output/evpi.csv
# ============================================================================

library(dplyr); library(tidyr); library(purrr); library(tibble); library(here)

results_full               <- readRDS(here("SalmonCountR","app_data","results_full.rds"))
steelhead_scenario_results <- readRDS(here("SalmonCountR","app_data","steelhead_scenario_results.rds"))

real_years  <- 2011:2024
scen_map    <- c(NB=1, PB1=2, PB2=3, PB2b=4, PB2c=5, PB3=6, PB4=7, PB5=8, PB6=9)
hydro_years <- c("2011","2014","2017","2020")
wyt_w       <- c("2011"=.25, "2014"=.25, "2017"=.25, "2020"=.25)
tdm_prob    <- c(exp_WF = .51, exp_SM = .24, lin_Martin = .25)

hydro_loss <- c(NB=0, PB1=111422, PB2=370826, PB2b=470090, PB2c=433215,
                PB3=201552, PB4=241590, PB5=199382, PB6=348806)

# Adult index by alternative × TDM model, averaged over met years
med20 <- results_full %>%
  filter(year > max(real_years)) %>%
  group_by(env, variant) %>% slice_tail(n = 20) %>%
  summarise(med = median(spawners, na.rm = TRUE), .groups = "drop") %>%
  mutate(env = as.integer(env),
         scenario = names(scen_map)[match((env - 1) %% 9 + 1, scen_map)],
         wyt      = hydro_years[(env - 1) %/% 9 + 1])

by_state <- med20 %>%
  mutate(wt = wyt_w[wyt]) %>%
  group_by(scenario, variant) %>%
  summarise(metric = sum(med * wt), .groups = "drop") %>%
  left_join(steelhead_scenario_results, by = "scenario") %>%
  mutate(hydro_raw = hydro_loss[scenario])

norm_hi <- function(x, lo = min(x), hi = max(x)) (x - lo) / (hi - lo)
norm_lo <- function(x, lo = min(x), hi = max(x)) (hi - x) / (hi - lo)

score_table <- function(w_salmon, w_hydro, w_steel, pooled) {
  d <- by_state
  if (pooled) {
    lo <- min(d$metric); hi <- max(d$metric)          # one scale across states
    d <- d %>% mutate(chinook_n = norm_hi(metric, lo, hi))
  } else {
    d <- d %>% group_by(variant) %>%                   # rescale within state
      mutate(chinook_n = norm_hi(metric)) %>% ungroup()
  }
  d %>% mutate(steel_n = norm_hi(steelhead_score),
               hydro_n = norm_lo(hydro_raw),
               score   = w_salmon*chinook_n + w_steel*steel_n + w_hydro*hydro_n)
}

evpi <- function(w_salmon, w_hydro, w_steel, pooled) {
  s <- score_table(w_salmon, w_hydro, w_steel, pooled)
  # max over alternatives of the expected score (best action under uncertainty)
  exp_score <- s %>% mutate(p = tdm_prob[variant]) %>%
    group_by(scenario) %>% summarise(es = sum(score * p), .groups = "drop")
  best_under_uncertainty <- exp_score %>% slice_max(es, n = 1, with_ties = FALSE)
  # expected max: best alternative in each state, weighted by state probability
  per_state_best <- s %>% group_by(variant) %>%
    slice_max(score, n = 1, with_ties = FALSE) %>% ungroup() %>%
    mutate(p = tdm_prob[variant])
  e_max <- sum(per_state_best$score * per_state_best$p)
  tibble(
    baseline_alt   = best_under_uncertainty$scenario,
    max_expected   = best_under_uncertainty$es,
    expected_max   = e_max,
    evpi           = e_max - best_under_uncertainty$es,
    evpi_pct       = 100 * (e_max - best_under_uncertainty$es) / best_under_uncertainty$es,
    best_by_state  = paste(sprintf("%s->%s", per_state_best$variant,
                                   per_state_best$scenario), collapse = "; ")
  )
}

grid <- tribble(
  ~weight_set,                        ~w_salmon, ~w_hydro, ~w_steel,
  "elicited (0.40/0.50/0.10)",              .40,      .50,      .10,
  "re-derived (0.73/0.22/0.05)",            .73,      .22,      .05
) %>%
  crossing(normalisation = c("pooled", "per-state")) %>%
  mutate(res = pmap(list(w_salmon, w_hydro, w_steel, normalisation == "pooled"), evpi)) %>%
  unnest(res)

out <- grid %>% select(weight_set, normalisation, baseline_alt, max_expected,
                       expected_max, evpi, evpi_pct, best_by_state)

cat("=== EVPI on the TDM model choice ===\n")
print(as.data.frame(out %>% mutate(across(where(is.numeric), ~round(., 4)))),
      row.names = FALSE)

dir.create(here("output"), showWarnings = FALSE)
write.csv(out, here("output", "evpi.csv"), row.names = FALSE)
cat("\nWrote output/evpi.csv\n")
