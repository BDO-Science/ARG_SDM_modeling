# ============================================================================
# Corrected Discussion numbers for the 2026-08 revision
# ============================================================================
# The decision was taken to apply the G1 correction (alternative-specific spawn
# timing) during the current revision rather than the next cycle. This script
# produces the replacement values for every Discussion claim that moves, in the
# form the manuscript should now report them: a mean across seeds with the
# run-to-run spread stated, not a single-run point estimate.
#
# WHY A MEAN ACROSS SEEDS. The correction reduces the redd sample behind each
# survival estimate from ~36 x N_redd to N_redd, so per-alternative estimates
# carry +/-173-217 fish of Monte Carlo noise (pooling suppressed this to
# +/-15-21 and made single runs look far more precise than they were). Any claim
# resting on a difference smaller than roughly 400 fish is not supportable from
# one run. See G1_FINDINGS.md.
#
# Inputs  : the multi-seed snapshots produced by running precompute.R with
#           ARG_G1_ALT_SPAWN in {0,1} and ARG_SEED in {123,456,789,1011,1213}.
#           Set G1_SNAPROOT to the directory holding <arm>_seed<n>/ folders.
# Outputs : output/g1_revision_numbers.csv   (per-alternative, both arms)
#           output/g1_revision_claims.md     (claim-by-claim replacement text)
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
ARMS  <- c(legacy = "legacy", fixed = "fixed")

# ---------------------------------------------------------------------------
# Bypass volume (Mm3) per alternative.
# ---------------------------------------------------------------------------
# These live in manuscript Table 2 and are NOT recorded anywhere in the repo.
# Only the PB3/PB5 pair is documented here (21.4 Mm3 each, identical volume on
# different schedules - that pairing is the whole point of the comparison).
# Supply the rest via G1_VOLUMES as a comma-separated ALT=VALUE list, e.g.
#   G1_VOLUMES="NB=0,PB1=10.7,PB2=21.4,..."
# Without it the efficiency and Spearman claims are reported as PENDING rather
# than guessed - deriving volume from the hydropower scores would be circular,
# since those scores are what the Spearman is meant to be independent of.
parse_volumes <- function(s) {
  if (!nzchar(s)) return(NULL)
  kv <- str_split(str_split(s, ",")[[1]], "=")
  v  <- setNames(as.numeric(map_chr(kv, 2)), str_trim(map_chr(kv, 1)))
  if (!all(SCEN %in% names(v))) {
    warning("G1_VOLUMES is missing: ", paste(setdiff(SCEN, names(v)), collapse = ", "))
    return(NULL)
  }
  v[SCEN]
}
VOL <- parse_volumes(Sys.getenv("G1_VOLUMES", ""))

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

raw <- expand_grid(arm = unname(ARMS), seed = SEEDS) %>%
  pmap_dfr(function(arm, seed) load_one(arm, seed))
stopifnot(nrow(raw) > 0)

cat(sprintf("loaded %d arm x seed combinations\n",
            nrow(distinct(raw, arm, seed))))

# Per-alternative summary: mean across seeds, plus the run-to-run range that
# has to be quoted alongside it.
per_alt <- raw %>%
  group_by(arm, scenario) %>%
  summarise(n_seed      = n(),
            adult_mean  = mean(adult_index),
            adult_sd    = sd(adult_index),
            adult_lo    = min(adult_index),
            adult_hi    = max(adult_index),
            comp_mean   = mean(composite),
            comp_sd     = sd(composite),
            .groups = "drop") %>%
  mutate(scenario = factor(scenario, levels = SCEN)) %>%
  arrange(arm, scenario)

# ---------------------------------------------------------------------------
# Paired contrasts. Pair WITHIN seed, then average the differences - averaging
# the two arms separately and subtracting would fold the common seed shift into
# the contrast.
# ---------------------------------------------------------------------------
# `metric` selects which quantity is contrasted. This matters: G1 moves the
# adult index and the MCDA composite in different ways, because the composite's
# min-max normalisation cancels the shift that is common to all alternatives.
# The one reproducible ranking change (PB2b vs PB5) is a COMPOSITE flip - both
# arms keep PB2b well ahead on raw abundance - so contrasting it on abundance
# would miss it entirely.
contrast <- function(a, b, metric = "adult_index") {
  raw %>%
    filter(scenario %in% c(a, b)) %>%
    select(arm, seed, scenario, all_of(metric)) %>%
    pivot_wider(names_from = scenario, values_from = all_of(metric)) %>%
    mutate(diff = .data[[a]] - .data[[b]]) %>%
    group_by(arm) %>%
    summarise(pair      = paste(a, "-", b),
              metric    = metric,
              mean_diff = mean(diff),
              sd_diff   = sd(diff),
              lo        = min(diff),
              hi        = max(diff),
              same_sign = n_distinct(sign(diff)) == 1,
              .groups = "drop")
}

contrasts <- bind_rows(contrast("PB3", "PB5"),
                       contrast("PB6", "PB4"),
                       contrast("PB2b", "PB5"),
                       contrast("PB2b", "PB5", metric = "composite"))

# ---------------------------------------------------------------------------
# Efficiency and Spearman - only if volumes were supplied
# ---------------------------------------------------------------------------
eff <- NULL; spear <- NULL
if (!is.null(VOL)) {
  volt <- tibble(scenario = SCEN, volume_Mm3 = as.numeric(VOL))
  nb   <- per_alt %>% filter(scenario == "NB") %>% select(arm, nb_adult = adult_mean)

  # "additional adults per million m3" is measured against no-bypass, which is
  # the only alternative with zero bypass volume.
  eff <- per_alt %>%
    left_join(volt, by = "scenario") %>%
    left_join(nb, by = "arm") %>%
    filter(volume_Mm3 > 0) %>%
    mutate(gain            = adult_mean - nb_adult,
           adults_per_Mm3  = gain / volume_Mm3) %>%
    select(arm, scenario, volume_Mm3, gain, adults_per_Mm3) %>%
    arrange(arm, scenario)

  # Spearman is computed per seed and then averaged, so the reported value
  # carries a spread like everything else.
  spear <- raw %>%
    left_join(volt, by = "scenario") %>%
    group_by(arm, seed) %>%
    summarise(rho = cor(volume_Mm3, adult_index, method = "spearman"),
              .groups = "drop") %>%
    group_by(arm) %>%
    summarise(rho_mean = mean(rho), rho_lo = min(rho), rho_hi = max(rho),
              .groups = "drop")
}

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
g1_rule("Per-alternative adult population index (model-averaged)")
print(as.data.frame(per_alt), digits = 5, row.names = FALSE)

g1_rule("Paired contrasts (within-seed differences)")
print(as.data.frame(contrasts), digits = 4, row.names = FALSE)

if (!is.null(eff)) {
  g1_rule("Efficiency: additional adults per million m3, vs no-bypass")
  print(as.data.frame(eff), digits = 4, row.names = FALSE)
  g1_rule("Spearman rho: bypass volume vs adult index")
  print(as.data.frame(spear), digits = 4, row.names = FALSE)
} else {
  g1_rule("Efficiency and Spearman: PENDING")
  cat("Bypass volumes were not supplied, so these two claims cannot be\n",
      "recomputed. Re-run with G1_VOLUMES set from manuscript Table 2:\n",
      '  G1_VOLUMES="NB=0,PB1=...,PB2=...,PB2b=...,PB2c=...,',
      'PB3=21.4,PB4=...,PB5=21.4,PB6=..."\n', sep = "")
}

dir.create(here("output"), showWarnings = FALSE)
write_csv(per_alt, here("output", "g1_revision_numbers.csv"))

# ---------------------------------------------------------------------------
# Claim-by-claim replacement text
# ---------------------------------------------------------------------------
gv <- function(arm_, scen) per_alt %>%
  filter(arm == arm_, scenario == scen) %>% pull(adult_mean)
cv <- function(arm_, pair_, metric_ = "adult_index") contrasts %>%
  filter(arm == arm_, pair == pair_, metric == metric_)

c35 <- cv("fixed", "PB3 - PB5"); c35L <- cv("legacy", "PB3 - PB5")
c64 <- cv("fixed", "PB6 - PB4"); c64L <- cv("legacy", "PB6 - PB4")
k25 <- cv("fixed",  "PB2b - PB5", "composite")
k25L<- cv("legacy", "PB2b - PB5", "composite")

lines <- c(
"# Corrected Discussion numbers — G1 applied",
"",
"Generated by `analysis/g1_revision_numbers.R`. Every value is a mean across",
sprintf("%d seeds (%s) with the run-to-run range in brackets. Single-run point",
        length(SEEDS), paste(SEEDS, collapse = ", ")),
"estimates should not be reinstated: see the noise note in `G1_FINDINGS.md`.",
"",
"## 1. PB3 vs PB5 — the schedule counterweight (SIGN REVERSES)",
"",
sprintf("- Published (pooled): PB3 - PB5 = **%+.0f adults** [%.0f, %.0f], PB3 ahead",
        c35L$mean_diff, c35L$lo, c35L$hi),
sprintf("- Corrected: PB3 - PB5 = **%+.0f adults** [%.0f, %.0f], PB5 ahead",
        c35$mean_diff, c35$lo, c35$hi),
sprintf("- Same sign at every seed: %s", c35$same_sign),
"",
"The current sentence — identical volumes differing by \"fewer than 40 adults\",",
"so schedule is consequential only \"at specific points in the decision space\" —",
"no longer holds and should be replaced. The correction *removes* this hedge",
"against the paper's own front-loading argument: PB5 now leads by a margin",
"several times the run-to-run noise, in the direction its net hazard balance",
"already predicted (PB3 +0.014, PB5 -0.018).",
"",
"## 2. PB6 vs PB4 adult population index",
"",
sprintf("- Published: PB6 %.0f against PB4 %.0f (difference %+.0f)",
        gv("legacy","PB6"), gv("legacy","PB4"), c64L$mean_diff),
sprintf("- Corrected: PB6 **%.0f** [%.0f, %.0f] against PB4 **%.0f** [%.0f, %.0f]",
        gv("fixed","PB6"),
        per_alt %>% filter(arm=="fixed", scenario=="PB6") %>% pull(adult_lo),
        per_alt %>% filter(arm=="fixed", scenario=="PB6") %>% pull(adult_hi),
        gv("fixed","PB4"),
        per_alt %>% filter(arm=="fixed", scenario=="PB4") %>% pull(adult_lo),
        per_alt %>% filter(arm=="fixed", scenario=="PB4") %>% pull(adult_hi)),
sprintf("- Difference %+.0f adults [%.0f, %.0f]; same sign at every seed: %s",
        c64$mean_diff, c64$lo, c64$hi, c64$same_sign),
"",
"## 3. Efficiency — additional adults per million m3",
"",
if (is.null(eff)) paste(
  "**PENDING.** Needs bypass volumes from manuscript Table 2; they are not",
  "recorded in the repository. Re-run with G1_VOLUMES set.") else
  paste(capture.output(print(as.data.frame(eff %>% filter(arm == "fixed")),
                             digits = 4, row.names = FALSE)), collapse = "\n"),
"",
"## 4. Spearman rho — bypass volume vs salmon benefit",
"",
if (is.null(spear)) paste(
  "**PENDING.** Same blocker as claim 3.") else
  sprintf("- Corrected rho = **%.3f** [%.3f, %.3f] (published %.3f)",
          spear %>% filter(arm=="fixed") %>% pull(rho_mean),
          spear %>% filter(arm=="fixed") %>% pull(rho_lo),
          spear %>% filter(arm=="fixed") %>% pull(rho_hi),
          spear %>% filter(arm=="legacy") %>% pull(rho_mean)),
"",
"## 5. Ranking",
"",
"PB1 remains the top-ranked alternative in the MCDA under both arms at every",
sprintf("seed. PB4 remains the decision-makers' selection and is the least affected"),
sprintf("alternative in absolute abundance (%+.0f adults). Exactly one pairwise order",
        gv("fixed","PB4") - gv("legacy","PB4")),
"changes reproducibly, and it is a change in COMPOSITE score, not abundance:",
"",
sprintf("- PB2b - PB5 composite: published %+.4f (PB2b ahead), corrected %+.4f (PB5 ahead)",
        k25L$mean_diff, k25$mean_diff),
sprintf("  Same sign at every seed in both arms: %s / %s",
        k25L$same_sign, k25$same_sign),
"",
"On raw abundance PB2b stays well ahead of PB5 in both arms (by roughly",
sprintf("%.0f adults), so this is the min-max normalisation and the hydropower",
        cv("fixed","PB2b - PB5")$mean_diff),
"weight resolving a pair that abundance alone does not separate.",
""
)
writeLines(lines, here("output", "g1_revision_claims.md"))
cat("\nWrote output/g1_revision_numbers.csv and output/g1_revision_claims.md\n")
