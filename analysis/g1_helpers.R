# ============================================================================
# Shared helpers for the G1 comparison scripts
# ============================================================================
# Sourced by analysis/compare_g1.R (single pair) and
# analysis/compare_g1_seeds.R (multi-seed replication).
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
G1_HYDRO_SCORES <- c(
  "NB"   = 0,      "PB1" = 111422, "PB2" = 370826,
  "PB2b" = 470090, "PB2c" = 433215, "PB3" = 201552,
  "PB4"  = 241590, "PB5" = 199382,  "PB6" = 348806
)

G1_W_CHINOOK   <- 0.4
G1_W_STEELHEAD <- 0.1
G1_W_HYDRO     <- 0.5

g1_read <- function(dir, name) {
  p <- file.path(dir, paste0(name, ".rds"))
  if (!file.exists(p)) return(NULL)
  readRDS(p)
}

g1_norm_minmax <- function(x) (x - min(x)) / (max(x) - min(x))
g1_norm_invert <- function(x) (max(x) - x) / (max(x) - min(x))  # lower = better

#' Composite MCDA score for one app_data snapshot
#'
#' @param dir Directory holding precompute.R's .rds output.
#' @return Tibble ordered best-to-worst with `scenario`, raw and normalised
#'   objective values, `composite`, and `rank`. NULL if inputs are absent.
g1_composite <- function(dir) {
  sw <- g1_read(dir, "swing_scenario_results")
  st <- g1_read(dir, "steelhead_scenario_results")
  if (is.null(sw) || is.null(st)) return(NULL)

  sw %>%
    rename(chinook_raw = spawner_metric) %>%
    left_join(st %>% rename(steelhead_raw = steelhead_score), by = "scenario") %>%
    left_join(tibble(scenario  = names(G1_HYDRO_SCORES),
                     hydro_raw = as.numeric(G1_HYDRO_SCORES)),
              by = "scenario") %>%
    mutate(
      chinook_norm   = g1_norm_minmax(chinook_raw),
      steelhead_norm = g1_norm_minmax(steelhead_raw),
      hydro_norm     = g1_norm_invert(hydro_raw),
      composite      = G1_W_CHINOOK   * chinook_norm +
                       G1_W_STEELHEAD * steelhead_norm +
                       G1_W_HYDRO     * hydro_norm
    ) %>%
    arrange(desc(composite)) %>%
    mutate(rank = row_number())
}

#' Ranking as a single comparable string
g1_ranking <- function(comp) paste(comp$scenario, collapse = " > ")

g1_rule <- function(txt) {
  cat("\n", strrep("=", 78), "\n", txt, "\n", strrep("=", 78), "\n", sep = "")
}
