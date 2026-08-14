# ============================================================================
# G1 seed replication — is the effect larger than Monte Carlo noise?
# ============================================================================
# The alternative-specific fix reduces the redd sample behind each survival
# estimate from ~36 x N_redd (all alternatives pooled) to N_redd, so the
# corrected pipeline is inherently noisier. A raw legacy-vs-fixed difference
# therefore confounds the real signal with sampling noise.
#
# This script separates them by running both arms across several RNG seeds:
#
#   noise  = spread WITHIN an arm across seeds
#   effect = paired (fixed - legacy) difference at the SAME seed
#
# The effect is only credible where it is consistent in sign across seeds and
# large relative to the within-arm spread.
#
# Usage (from the repo root):
#   G1_SNAPROOT=<dir> Rscript analysis/compare_g1_seeds.R
#
# Expects subdirectories named  legacy_seed<N>  and  fixed_seed<N>.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(here)
})
source(here("analysis", "g1_helpers.R"))

snaproot <- Sys.getenv("G1_SNAPROOT")
stopifnot(dir.exists(snaproot))

dirs <- list.dirs(snaproot, recursive = FALSE)
meta <- tibble(path = dirs, base = basename(dirs)) %>%
  filter(grepl("^(legacy|fixed)_seed[0-9]+$", base)) %>%
  mutate(arm  = sub("_seed[0-9]+$", "", base),
         seed = as.integer(sub("^.*_seed", "", base)))

if (!nrow(meta)) stop("no legacy_seed*/fixed_seed* snapshots found in ", snaproot)

# Keep only seeds present in BOTH arms so every comparison is paired
paired_seeds <- meta %>% count(seed) %>% filter(n == 2) %>% pull(seed) %>% sort()
meta <- meta %>% filter(seed %in% paired_seeds) %>% arrange(arm, seed)

cat("paired seeds:", paste(paired_seeds, collapse = ", "),
    sprintf("  (%d replicate%s per arm)\n", length(paired_seeds),
            if (length(paired_seeds) == 1) "" else "s"))
if (length(paired_seeds) < 3) {
  cat("NOTE: fewer than 3 paired seeds — treat the noise estimate as indicative only.\n")
}

scores <- meta %>%
  mutate(comp = map(path, g1_composite)) %>%
  filter(!map_lgl(comp, is.null)) %>%
  select(arm, seed, comp) %>%
  unnest(comp)

# ---------------------------------------------------------------------------
g1_rule("1. RANKING BY ARM AND SEED")
# ---------------------------------------------------------------------------
rankings <- scores %>%
  arrange(arm, seed, rank) %>%
  group_by(arm, seed) %>%
  summarise(ranking = paste(scenario, collapse = " > "),
            top     = scenario[which.min(rank)], .groups = "drop")
print(as.data.frame(rankings), row.names = FALSE)

within_legacy <- n_distinct(rankings$ranking[rankings$arm == "legacy"])
within_fixed  <- n_distinct(rankings$ranking[rankings$arm == "fixed"])
cat(sprintf("\ndistinct rankings within legacy across seeds: %d\n", within_legacy))
cat(sprintf("distinct rankings within fixed  across seeds: %d\n", within_fixed))
cat(sprintf("legacy ranking stable across seeds: %s\n", within_legacy == 1))
cat(sprintf("fixed  ranking stable across seeds: %s\n", within_fixed  == 1))

paired_rank_same <- rankings %>%
  select(arm, seed, ranking) %>%
  pivot_wider(names_from = arm, values_from = ranking) %>%
  mutate(same = legacy == fixed)
cat("\nlegacy vs fixed ranking identical, by seed:\n")
print(as.data.frame(paired_rank_same %>% select(seed, same)), row.names = FALSE)

# ---------------------------------------------------------------------------
g1_rule("2. CHINOOK SPAWNER METRIC — NOISE vs EFFECT")
# ---------------------------------------------------------------------------
noise_effect <- function(df, value_col, digits = 4) {
  v <- rlang::sym(value_col)

  noise <- df %>%
    group_by(scenario, arm) %>%
    summarise(sd_across_seeds = sd(!!v), .groups = "drop") %>%
    pivot_wider(names_from = arm, values_from = sd_across_seeds,
                names_prefix = "noise_sd_")

  effect <- df %>%
    select(scenario, arm, seed, val = !!v) %>%
    pivot_wider(names_from = arm, values_from = val) %>%
    mutate(diff = fixed - legacy) %>%
    group_by(scenario) %>%
    summarise(mean_effect = mean(diff),
              sd_effect   = sd(diff),
              n_pos       = sum(diff > 0),
              n_neg       = sum(diff < 0),
              .groups     = "drop")

  noise %>%
    left_join(effect, by = "scenario") %>%
    mutate(
      consistent_sign = pmax(n_pos, n_neg) == (n_pos + n_neg),
      # effect relative to the noisier arm's own seed-to-seed spread
      effect_vs_noise = abs(mean_effect) / pmax(noise_sd_fixed, noise_sd_legacy)
    )
}

chin <- noise_effect(scores %>% select(scenario, arm, seed, chinook_raw),
                     "chinook_raw")
print(as.data.frame(chin), digits = 5, row.names = FALSE)
cat("\nnoise_sd_* = SD across seeds within that arm (Monte Carlo noise)\n")
cat("mean_effect = mean paired (fixed - legacy) difference\n")
cat("effect_vs_noise = |mean_effect| / larger within-arm SD; >2 suggests signal\n")

# ---------------------------------------------------------------------------
g1_rule("3. MCDA COMPOSITE SCORE — NOISE vs EFFECT")
# ---------------------------------------------------------------------------
comp <- noise_effect(scores %>% select(scenario, arm, seed, composite),
                     "composite")
print(as.data.frame(comp), digits = 5, row.names = FALSE)

# ---------------------------------------------------------------------------
g1_rule("4. VERDICT")
# ---------------------------------------------------------------------------
n_signal <- sum(chin$effect_vs_noise > 2 & chin$consistent_sign, na.rm = TRUE)
cat(sprintf("scenarios where the Chinook effect exceeds 2x noise AND is\n"))
cat(sprintf("  sign-consistent across all seeds: %d of %d\n",
            n_signal, nrow(chin)))
cat(sprintf("max |mean effect| on the composite score: %.5f\n",
            max(abs(comp$mean_effect), na.rm = TRUE)))
cat(sprintf("top choice by arm: legacy = %s | fixed = %s\n",
            paste(unique(rankings$top[rankings$arm == "legacy"]), collapse = "/"),
            paste(unique(rankings$top[rankings$arm == "fixed"]),  collapse = "/")))
cat(sprintf("\nranking changes attributable to G1 (not noise): %s\n",
            if (all(paired_rank_same$same)) "NONE — ranking identical at every seed"
            else if (within_legacy > 1 || within_fixed > 1)
              "AMBIGUOUS — ranking is not even stable across seeds within an arm"
            else "YES — ranking differs between arms and is stable within each arm"))
