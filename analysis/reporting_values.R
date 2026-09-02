# ============================================================================
# Headline objective values
# ============================================================================
# Single source for the summary values quoted from this model - adult population
# index, MCDA composite and rank, volume-normalised benefit, and the
# volume-vs-benefit rank correlation. Reads the committed app_data, so a clone
# reproduces them in seconds without re-running precompute.R.
#
# REPORTING CONVENTION (see docs/spawn-timing.md)
#
#   Levels     - reported from this single committed run. Every figure, the app
#                and this table therefore trace to one artefact.
#   Spread     - the run-to-run range from the five-seed replication is quoted
#                alongside, so the precision of a level is never overstated.
#   Contrasts  - differences between alternatives are reported as PAIRED
#                within-seed means from the five-seed set, not as a subtraction
#                of two levels here.
#
# Why contrasts are handled differently: the seed shift is largely common across
# alternatives, so it cancels in a paired difference. Levels carry +/-173-217
# fish of run-to-run noise while the contrasts that matter carry only +/-21-36.
# Subtracting two levels from this file would discard that cancellation and
# understate what the model can actually resolve. Contrasts come from
# analysis/spawn_timing_effect.R (needs SPAWN_TIMING_SNAPROOT).
#
# Outputs: output/reporting_values.csv
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(readr); library(here)
})
source(here("analysis", "composite_helpers.R"))

APP <- here("SalmonCountR", "app_data")

# Bypass volume (Mm3) per alternative - a design input. Kept in step with
# analysis/spawn_timing_effect.R; see the note there on why these must not be
# derived from the hydropower revenue-loss scores.
VOLUME_MM3 <- c("NB" = 0.0, "PB1" = 12.2, "PB2" = 42.2, "PB2b" = 49.5,
                "PB2c" = 45.9, "PB3" = 21.4, "PB4" = 25.7, "PB5" = 21.4,
                "PB6" = 37.2)

comp <- composite_scores(APP)
stopifnot(!is.null(comp))

nb_index <- comp$chinook_raw[comp$scenario == "NB"]

vals <- comp %>%
  transmute(scenario,
            adult_index = chinook_raw,
            composite,
            rank,
            volume_Mm3 = unname(VOLUME_MM3[scenario]),
            gain_vs_NB = adult_index - nb_index,
            adults_per_Mm3 = ifelse(volume_Mm3 > 0, gain_vs_NB / volume_Mm3, NA_real_))

rho <- cor(vals$volume_Mm3, vals$adult_index, method = "spearman")

rule("Reported values - committed run")
print(as.data.frame(
  vals %>% arrange(desc(composite)) %>%
    mutate(across(c(adult_index, gain_vs_NB), ~round(.x, 0)),
           across(c(composite, adults_per_Mm3), ~round(.x, 3)))),
  row.names = FALSE)

cat(sprintf("\nMCDA ordering: %s\n", paste(vals$scenario[order(-vals$composite)],
                                           collapse = " > ")))
cat(sprintf("Adult index range: %.0f to %.0f\n",
            min(vals$adult_index), max(vals$adult_index)))
cat(sprintf("Spearman rho, bypass volume vs adult index: %.4f\n", rho))

cat("\nQuote alongside these levels: run-to-run spread across five seeds is\n",
    "+/-173-217 fish per alternative. Differences smaller than roughly 400 fish\n",
    "are not resolvable from one run - report those as paired contrasts instead.\n",
    sep = "")

dir.create(here("output"), showWarnings = FALSE)
write_csv(vals, here("output", "reporting_values.csv"))
cat("\nWrote output/reporting_values.csv\n")
