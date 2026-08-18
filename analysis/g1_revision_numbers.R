# ============================================================================
# G1 effect on the reported objective values
# ============================================================================
# Quantifies what the G1 correction (alternative-specific spawn timing) changes
# in the model's output, for every quantity the project reports: the adult
# population index, the MCDA composite, volume-normalised benefit, and the
# volume-vs-benefit rank correlation.
#
# WHY MEANS ACROSS SEEDS. The correction reduces the redd sample behind each
# survival estimate from ~36 x N_redd to N_redd, so per-alternative estimates
# carry +/-173-217 fish of Monte Carlo noise; pooling suppressed this to
# +/-15-21 and made single runs look far more precise than they were. Any
# difference smaller than roughly 400 fish is not resolvable from one run, so
# every value here is a mean over seeds with the run-to-run range alongside it.
#
# Contrasts are paired WITHIN seed. Averaging each arm separately and
# subtracting would fold the common seed shift into the contrast.
#
# Inputs  : snapshots from running precompute.R with ARG_G1_ALT_SPAWN in {0,1}
#           and ARG_SEED in {123,456,789,1011,1213}. Point G1_SNAPROOT at the
#           directory holding the <arm>_seed<n>/ folders.
# Outputs : output/g1_revision_numbers.csv    per-alternative, both arms
#           output/g1_revision_contrasts.csv  paired contrasts
#           output/g1_revision_efficiency.csv volume-normalised benefit + rho
#
# Usage:
#   G1_SNAPROOT=<dir> Rscript analysis/g1_revision_numbers.R
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(readr); library(here); library(stringr)
})
source(here("analysis", "g1_helpers.R"))

snaproot <- Sys.getenv("G1_SNAPROOT")
stopifnot(nzchar(snaproot), dir.exists(snaproot))

SCEN  <- c("NB","PB1","PB2","PB2b","PB2c","PB3","PB4","PB5","PB6")
SEEDS <- c(123, 456, 789, 1011, 1213)
ARMS  <- c("legacy", "fixed")

# ---------------------------------------------------------------------------
# Bypass volume (Mm3) per alternative - a design input, not a model output.
# ---------------------------------------------------------------------------
# Do NOT derive these from the hydropower revenue-loss scores. Cost scales
# roughly with volume but not exactly - PB3 and PB5 bypass identical volumes and
# differ in cost by about 1% - and the rank correlation below is meant to relate
# volume to fish benefit independently of the cost model.
G1_VOLUME_MM3 <- c(
  "NB"   =  0.0, "PB1" = 12.2, "PB2" = 42.2,
  "PB2b" = 49.5, "PB2c" = 45.9, "PB3" = 21.4,
  "PB4"  = 25.7, "PB5" = 21.4,  "PB6" = 37.2
)

# Override with G1_VOLUMES="NB=0,PB1=12.2,..." if the operational definitions
# of the alternatives change.
parse_volumes <- function(s) {
  if (!nzchar(s)) return(G1_VOLUME_MM3[SCEN])
  kv <- str_split(str_split(s, ",")[[1]], "=")
  v  <- setNames(as.numeric(map_chr(kv, 2)), str_trim(map_chr(kv, 1)))
  if (!all(SCEN %in% names(v))) {
    warning("G1_VOLUMES is missing: ", paste(setdiff(SCEN, names(v)), collapse = ", "),
            " - falling back to the built-in values")
    return(G1_VOLUME_MM3[SCEN])
  }
  v[SCEN]
}
VOL  <- parse_volumes(Sys.getenv("G1_VOLUMES", ""))
volt <- tibble(scenario = SCEN, volume_Mm3 = as.numeric(VOL))

# ---------------------------------------------------------------------------
# Load every arm x seed
# ---------------------------------------------------------------------------
load_one <- function(arm, seed) {
  d <- file.path(snaproot, sprintf("%s_seed%d", arm, seed))
  if (!dir.exists(d)) { warning("missing snapshot: ", d); return(NULL) }
  sw <- g1_read(d, "swing_scenario_results")
  cm <- g1_composite(d)
  if (is.null(sw) || is.null(cm)) return(NULL)
  sw %>%
    rename(adult_index = spawner_metric) %>%
    left_join(cm %>% select(scenario, composite, rank), by = "scenario") %>%
    mutate(arm = arm, seed = seed)
}

raw <- expand_grid(arm = ARMS, seed = SEEDS) %>%
  pmap_dfr(function(arm, seed) load_one(arm, seed))
stopifnot(nrow(raw) > 0)
cat(sprintf("loaded %d arm x seed combinations\n", nrow(distinct(raw, arm, seed))))

per_alt <- raw %>%
  group_by(arm, scenario) %>%
  summarise(n_seed     = n(),
            adult_mean = mean(adult_index), adult_sd = sd(adult_index),
            adult_lo   = min(adult_index),  adult_hi = max(adult_index),
            comp_mean  = mean(composite),   comp_sd  = sd(composite),
            .groups = "drop") %>%
  mutate(scenario = factor(scenario, levels = SCEN)) %>%
  arrange(arm, scenario)

# ---------------------------------------------------------------------------
# Paired contrasts
# ---------------------------------------------------------------------------
# `metric` matters: G1 moves the adult index and the composite differently,
# because the composite's min-max normalisation cancels the shift common to all
# alternatives. The one reproducible ranking change (PB2b vs PB5) is a COMPOSITE
# change - both arms keep PB2b well ahead on raw abundance - so contrasting that
# pair on abundance alone would miss it.
contrast <- function(a, b, metric = "adult_index") {
  raw %>%
    filter(scenario %in% c(a, b)) %>%
    select(arm, seed, scenario, all_of(metric)) %>%
    pivot_wider(names_from = scenario, values_from = all_of(metric)) %>%
    mutate(diff = .data[[a]] - .data[[b]]) %>%
    group_by(arm) %>%
    summarise(pair = paste(a, "-", b), metric = metric,
              mean_diff = mean(diff), sd_diff = sd(diff),
              lo = min(diff), hi = max(diff),
              same_sign = n_distinct(sign(diff)) == 1,
              .groups = "drop")
}

contrasts <- bind_rows(
  contrast("PB3",  "PB5"),
  contrast("PB6",  "PB4"),
  contrast("PB2b", "PB5"),
  contrast("PB2b", "PB5", metric = "composite")
)

# ---------------------------------------------------------------------------
# Volume-normalised benefit, and the volume-vs-benefit rank correlation
# ---------------------------------------------------------------------------
# Benefit is measured against no-bypass, the only alternative with zero volume.
nb <- per_alt %>% filter(scenario == "NB") %>% select(arm, nb_adult = adult_mean)

eff <- per_alt %>%
  left_join(volt, by = "scenario") %>%
  left_join(nb, by = "arm") %>%
  filter(volume_Mm3 > 0) %>%
  mutate(gain = adult_mean - nb_adult, adults_per_Mm3 = gain / volume_Mm3) %>%
  select(arm, scenario, volume_Mm3, gain, adults_per_Mm3) %>%
  arrange(arm, scenario)

# rho is computed per seed then summarised, so it carries a spread like
# everything else. In practice the spread is zero: rho is rank-based, and the
# run-to-run noise shifts all nine alternatives together without reordering
# them. This is the one quantity here that is stable from a single run.
spear <- raw %>%
  left_join(volt, by = "scenario") %>%
  group_by(arm, seed) %>%
  summarise(rho = cor(volume_Mm3, adult_index, method = "spearman"),
            .groups = "drop") %>%
  group_by(arm) %>%
  summarise(rho_mean = mean(rho), rho_lo = min(rho), rho_hi = max(rho),
            .groups = "drop")

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
g1_rule("Per-alternative adult population index and composite")
print(as.data.frame(per_alt), digits = 5, row.names = FALSE)

g1_rule("Paired contrasts (within-seed differences)")
print(as.data.frame(contrasts), digits = 4, row.names = FALSE)

g1_rule("Volume-normalised benefit: additional adults per million m3 vs no-bypass")
print(as.data.frame(
  eff %>% select(arm, scenario, volume_Mm3, adults_per_Mm3) %>%
    pivot_wider(names_from = arm, values_from = adults_per_Mm3)),
  digits = 4, row.names = FALSE)

g1_rule("Rank correlation: bypass volume vs adult index")
print(as.data.frame(spear), digits = 4, row.names = FALSE)

dir.create(here("output"), showWarnings = FALSE)
write_csv(per_alt,   here("output", "g1_revision_numbers.csv"))
write_csv(contrasts, here("output", "g1_revision_contrasts.csv"))
write_csv(eff %>% left_join(spear, by = "arm"),
          here("output", "g1_revision_efficiency.csv"))

cat("\nWrote output/g1_revision_{numbers,contrasts,efficiency}.csv\n")
cat("Interpretation and the reporting write-up live in G1_FINDINGS.md and",
    "G1_DISCLOSURE.md, not here.\n")
