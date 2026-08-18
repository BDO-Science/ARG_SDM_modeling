# ============================================================================
# Shared MCDA composite helpers for the spawn-timing comparison scripts
# ============================================================================
# Sourced by analysis/compare_spawn_timing.R (single pair) and
# analysis/compare_spawn_timing_seeds.R (multi-seed replication).
#
# The composite score here mirrors analysis/mcda.R exactly: min-max
# normalisation within the nine scenarios, hydropower inverted, and the
# elicited weights 0.40 Chinook / 0.10 steelhead / 0.50 hydropower.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

# Hydropower revenue loss ($) per alternative — hard-coded in app.R and mcda.R
HYDRO_SCORES <- c(
  "NB"   = 0,      "PB1" = 111422, "PB2" = 370826,
  "PB2b" = 470090, "PB2c" = 433215, "PB3" = 201552,
  "PB4"  = 241590, "PB5" = 199382,  "PB6" = 348806
)

W_CHINOOK   <- 0.4
W_STEELHEAD <- 0.1
W_HYDRO     <- 0.5

read_snap <- function(dir, name) {
  p <- file.path(dir, paste0(name, ".rds"))
  if (!file.exists(p)) return(NULL)
  readRDS(p)
}

norm_minmax <- function(x) (x - min(x)) / (max(x) - min(x))
norm_invert <- function(x) (max(x) - x) / (max(x) - min(x))  # lower = better

#' Composite MCDA score for one app_data snapshot
#'
#' @param dir Directory holding precompute.R's .rds output.
#' @return Tibble ordered best-to-worst with `scenario`, raw and normalised
#'   objective values, `composite`, and `rank`. NULL if inputs are absent.
composite_scores <- function(dir) {
  sw <- read_snap(dir, "swing_scenario_results")
  st <- read_snap(dir, "steelhead_scenario_results")
  if (is.null(sw) || is.null(st)) return(NULL)

  sw %>%
    rename(chinook_raw = spawner_metric) %>%
    left_join(st %>% rename(steelhead_raw = steelhead_score), by = "scenario") %>%
    left_join(tibble(scenario  = names(HYDRO_SCORES),
                     hydro_raw = as.numeric(HYDRO_SCORES)),
              by = "scenario") %>%
    mutate(
      chinook_norm   = norm_minmax(chinook_raw),
      steelhead_norm = norm_minmax(steelhead_raw),
      hydro_norm     = norm_invert(hydro_raw),
      composite      = W_CHINOOK   * chinook_norm +
                       W_STEELHEAD * steelhead_norm +
                       W_HYDRO     * hydro_norm
    ) %>%
    arrange(desc(composite)) %>%
    mutate(rank = row_number())
}

#' Ranking as a single comparable string
ranking_string <- function(comp) paste(comp$scenario, collapse = " > ")

rule <- function(txt) {
  cat("\n", strrep("=", 78), "\n", txt, "\n", strrep("=", 78), "\n", sep = "")
}
