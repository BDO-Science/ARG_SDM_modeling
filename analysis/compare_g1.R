# ============================================================================
# G1 comparison — alternative-specific vs pooled spawn timing
# ============================================================================
# Diffs two snapshots of SalmonCountR/app_data and reports what changes in the
# quantities the decision rests on: the Chinook spawner metric, the MCDA
# composite scores, and the resulting ranking of alternatives.
#
# Usage:
#   G1_A=<dir> G1_A_LABEL=<name> G1_B=<dir> G1_B_LABEL=<name> \
#     Rscript analysis/compare_g1.R
#
# Both directories must contain the .rds files written by precompute.R.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(here)
})
source(here("analysis", "g1_helpers.R"))

dir_a   <- Sys.getenv("G1_A")
dir_b   <- Sys.getenv("G1_B")
label_a <- Sys.getenv("G1_A_LABEL", "A")
label_b <- Sys.getenv("G1_B_LABEL", "B")
stopifnot(dir.exists(dir_a), dir.exists(dir_b))

rd   <- g1_read
rule <- g1_rule

pct <- function(new, old) ifelse(old == 0, NA_real_, 100 * (new - old) / old)

ca <- g1_composite(dir_a)
cb <- g1_composite(dir_b)

rule(sprintf("1. CHINOOK SPAWNER METRIC  (%s -> %s)", label_a, label_b))
chin <- ca %>% select(scenario, a = chinook_raw) %>%
  left_join(cb %>% select(scenario, b = chinook_raw), by = "scenario") %>%
  mutate(abs_change = b - a, pct_change = pct(b, a)) %>%
  arrange(scenario)
print(as.data.frame(chin), digits = 6, row.names = FALSE)
cat(sprintf("\nmedian |%% change| = %.3f%%   max |%% change| = %.3f%%\n",
            median(abs(chin$pct_change), na.rm = TRUE),
            max(abs(chin$pct_change), na.rm = TRUE)))

rule(sprintf("2. MCDA COMPOSITE SCORE AND RANKING  (%s -> %s)", label_a, label_b))
comp <- ca %>% select(scenario, comp_a = composite, rank_a = rank) %>%
  left_join(cb %>% select(scenario, comp_b = composite, rank_b = rank),
            by = "scenario") %>%
  mutate(comp_change = comp_b - comp_a,
         rank_change = rank_a - rank_b) %>%
  arrange(rank_a)
print(as.data.frame(comp), digits = 6, row.names = FALSE)

order_a <- ca$scenario
order_b <- cb$scenario
cat("\nranking ", label_a, ": ", paste(order_a, collapse = " > "), "\n", sep = "")
cat(  "ranking ", label_b, ": ", paste(order_b, collapse = " > "), "\n", sep = "")

if (identical(order_a, order_b)) {
  cat("\n*** RANKING UNCHANGED ***\n")
} else {
  moved <- comp %>% filter(rank_change != 0)
  cat("\n*** RANKING CHANGED —", nrow(moved), "alternative(s) moved ***\n")
  print(as.data.frame(moved %>% select(scenario, rank_a, rank_b, rank_change)),
        row.names = FALSE)
}
cat(sprintf("\ntop choice: %s = %s | %s = %s%s\n",
            label_a, order_a[1], label_b, order_b[1],
            if (identical(order_a[1], order_b[1])) "  (same)" else "  <-- TOP CHOICE CHANGED"))

rule("3. STEELHEAD METRIC (expected identical — does not depend on spawn timing)")
sa <- rd(dir_a, "steelhead_scenario_results")
sb <- rd(dir_b, "steelhead_scenario_results")
if (isTRUE(all.equal(sa, sb))) {
  cat("identical.\n")
} else {
  cat("DIFFERENT — investigate:\n"); print(all.equal(sa, sb))
}

rule("4. EGG-TO-FRY SURVIVAL (egg_summary) BY VARIANT")
ea <- rd(dir_a, "egg_summary"); eb <- rd(dir_b, "egg_summary")
if (!is.null(ea) && !is.null(eb)) {
  num_col <- intersect(names(ea), c("mean_cum_surv", "cum_surv", "egg_surv"))[1]
  key     <- intersect(names(ea), c("env", "mgt_alt", "variant", "sim_year", "year"))
  if (!is.na(num_col)) {
    j <- ea %>% select(all_of(key), a = all_of(num_col)) %>%
      inner_join(eb %>% select(all_of(key), b = all_of(num_col)), by = key) %>%
      mutate(d = b - a)
    cat(sprintf("rows compared: %d (of %d / %d)\n", nrow(j), nrow(ea), nrow(eb)))
    by_variant <- j %>% group_by(variant) %>%
      summarise(mean_a = mean(a, na.rm = TRUE),
                mean_b = mean(b, na.rm = TRUE),
                mean_abs_diff = mean(abs(d), na.rm = TRUE),
                max_abs_diff  = max(abs(d), na.rm = TRUE),
                .groups = "drop")
    print(as.data.frame(by_variant), digits = 6, row.names = FALSE)
  } else {
    cat("could not identify survival column; names:", paste(names(ea), collapse = ", "), "\n")
  }
}

rule("5. CALIBRATED PARAMETERS (calib_results)")
ka <- rd(dir_a, "calib_results"); kb <- rd(dir_b, "calib_results")
cat(label_a, ":\n"); print(ka)
cat(label_b, ":\n"); print(kb)

rule("6. SWING RANGES")
cat(label_a, ":\n"); print(rd(dir_a, "swing_ranges"))
cat(label_b, ":\n"); print(rd(dir_b, "swing_ranges"))

rule("SUMMARY")
cat(sprintf("Chinook metric: median %.3f%%, max %.3f%% absolute change\n",
            median(abs(chin$pct_change), na.rm = TRUE),
            max(abs(chin$pct_change), na.rm = TRUE)))
cat(sprintf("Composite score: max absolute change %.5f\n",
            max(abs(comp$comp_change), na.rm = TRUE)))
cat(sprintf("Ranking: %s\n",
            if (identical(order_a, order_b)) "UNCHANGED" else "CHANGED"))
