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
g1_rule("4. PAIRWISE ORDER — WHICH RANK CHANGES ARE REAL")
# ---------------------------------------------------------------------------
# A whole-ranking comparison is too coarse: one noisy pair makes the entire
# ranking look unstable and hides the pairs that flip reproducibly. Order each
# pair of alternatives separately, in each arm, at every seed.
scs    <- sort(unique(scores$scenario))
combos <- utils::combn(scs, 2)

pairwise <- map_dfr(seq_len(ncol(combos)), function(i) {
  a <- combos[1, i]; b <- combos[2, i]
  d <- scores %>%
    filter(scenario %in% c(a, b)) %>%
    select(arm, seed, scenario, composite) %>%
    pivot_wider(names_from = scenario, values_from = composite) %>%
    mutate(a_ahead = .data[[a]] > .data[[b]])
  leg <- d$a_ahead[d$arm == "legacy"]
  fix <- d$a_ahead[d$arm == "fixed"]
  tibble(
    pair          = paste(a, "vs", b),
    legacy_stable = n_distinct(leg) == 1,
    fixed_stable  = n_distinct(fix) == 1,
    legacy_order  = if (n_distinct(leg) == 1)
                      (if (leg[1]) paste(a, ">", b) else paste(b, ">", a)) else "UNSTABLE",
    fixed_order   = if (n_distinct(fix) == 1)
                      (if (fix[1]) paste(a, ">", b) else paste(b, ">", a)) else "UNSTABLE",
    real_flip     = n_distinct(leg) == 1 && n_distinct(fix) == 1 && leg[1] != fix[1]
  )
})

real <- pairwise %>% filter(real_flip)
cat("Pairs that flip CONSISTENTLY between arms (stable within each arm) —\n")
cat("these are attributable to G1, not to noise:\n")
if (nrow(real)) {
  print(as.data.frame(real %>% select(pair, legacy_order, fixed_order)), row.names = FALSE)
} else {
  cat("  none\n")
}

unstable <- pairwise %>% filter(!fixed_stable | !legacy_stable)
cat("\nPairs whose order is NOT determinable from a single run (seed-dependent\n")
cat("in at least one arm) — a single-seed ranking overstates precision here:\n")
if (nrow(unstable)) {
  print(as.data.frame(unstable %>% select(pair, legacy_order, fixed_order)), row.names = FALSE)
} else {
  cat("  none\n")
}

# ---------------------------------------------------------------------------
g1_rule("5. VERDICT")
# ---------------------------------------------------------------------------
cat(sprintf("pairs flipping consistently (real G1 effect): %d of %d\n",
            nrow(real), nrow(pairwise)))
cat(sprintf("pairs not determinable from one run (noise):  %d of %d\n",
            nrow(unstable), nrow(pairwise)))
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
