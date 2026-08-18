# ============================================================================
# Regression test for the scenario engine
# ============================================================================
# Feeds the engine the SAME deliverable the published run was built from, and
# checks it lands in the same place. This is the guard that says the upload path
# in the app is telling the truth.
#
# Exact agreement is not expected and would be suspicious: the engine weights
# every spawn date by its CLM probability, where precompute.R draws a finite
# sample of redds. What must hold is that the ranking and the relative spacing
# of alternatives come out the same.
# ============================================================================

suppressPackageStartupMessages({library(dplyr); library(tibble); library(here)})

source(here("SalmonCountR", "functions.R"))
source(here("SalmonCountR", "scenario_engine.R"))

app <- function(...) here("SalmonCountR", "app_data", ...)

env_ext_list     <- readRDS(app("env_ext_list.rds"))
spawn_model      <- readRDS(app("spawn_timing_model.rds"))
base_P_list      <- readRDS(app("base_P_list.rds"))
S_seed_fore_list <- readRDS(app("S_seed_fore_list.rds"))
sim_years        <- readRDS(app("sim_years.rds"))
results_full     <- readRDS(app("results_full.rds"))
instream         <- readRDS(app("american_river_instream.rds")) %>%
  mutate(K_spawners = FR_spawn_wua / 9.29)

K_fn <- function(f) approx(instream$flow_cfs, instream$K_spawners, xout = f, rule = 2)$y

hydro_loss <- c(NB = 0, PB1 = 111422, PB2 = 370826, PB2b = 470090, PB2c = 433215,
                PB3 = 201552, PB4 = 241590, PB5 = 199382, PB6 = 348806)

deliverable <- here("data_raw", "SDM Power Bypass Temperature Modeling Results.xlsx")

t0 <- Sys.time()
res <- run_scenario_deliverable(
  path = deliverable, env_ext_list = env_ext_list, spawn_model = spawn_model,
  base_P_list = base_P_list, S_seed_fore_list = S_seed_fore_list,
  sim_years = sim_years, hydro_loss = hydro_loss,
  progress = function(f, m) cat(sprintf("  [%3.0f%%] %s\n", 100 * f, m))
)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cat("\n=== Validation report ===\n")
print(as.data.frame(res$report), row.names = FALSE)
stopifnot(res$ok)

cat(sprintf("\n=== Ran in %.2f seconds ===\n", elapsed))

cat("\n=== Engine scores ===\n")
print(as.data.frame(res$scores %>%
  transmute(rank, scenario, adult_index = round(adult_index),
            steelhead = round(steelhead_score, 1),
            composite = round(composite, 4))), row.names = FALSE)

# ---- Compare against the published run --------------------------------------
# Only the Chinook metric is directly comparable. swing_scenario_results.rds is
# built the same way the engine builds adult_index -- TDM-weighted median of the
# final 20 forecast years, averaged over met years -- so it is a like-for-like
# test. The steelhead metric and the composite are NOT comparable; see below.

published_chinook   <- readRDS(app("swing_scenario_results.rds"))
published_steelhead <- readRDS(app("steelhead_scenario_results.rds"))

cmp <- res$scores %>%
  select(scenario, engine = adult_index) %>%
  left_join(published_chinook %>% rename(published = spawner_metric), by = "scenario") %>%
  mutate(ratio = engine / published) %>%
  arrange(desc(published)) %>%
  mutate(published_rank = row_number()) %>%
  arrange(desc(engine)) %>%
  mutate(engine_rank = row_number()) %>%
  arrange(published_rank)

cat("\n=== Chinook metric: engine vs published (like for like) ===\n")
print(as.data.frame(cmp %>% transmute(
  scenario, engine = round(engine), published = round(published),
  ratio = round(ratio, 3), engine_rank, published_rank)), row.names = FALSE)

rho <- cor(cmp$engine, cmp$published, method = "spearman")
cat(sprintf("\n  Spearman rank correlation : %.3f\n", rho))
cat(sprintf("  Same rank position        : %d of %d\n",
            sum(cmp$engine_rank == cmp$published_rank), nrow(cmp)))
cat(sprintf("  Ratio range               : %.3f to %.3f\n",
            min(cmp$ratio), max(cmp$ratio)))

# Steelhead: ordering should agree, levels should not
sh <- res$steelhead %>%
  group_by(scenario) %>% summarise(engine = mean(steelhead_score), .groups = "drop") %>%
  left_join(published_steelhead %>% rename(published = steelhead_score), by = "scenario")
sh_rho <- cor(sh$engine, sh$published, method = "spearman")

cat("\n=== Steelhead metric: ordering agrees, level does not (expected) ===\n")
print(as.data.frame(sh %>% transmute(scenario, engine = round(engine, 1),
                                     published = round(published, 2),
                                     offset = round(engine - published, 1))),
      row.names = FALSE)
cat(sprintf("\n  Spearman rank correlation : %.3f\n", sh_rho))
cat("  The engine counts more days below 18.3 degC because its representative\n")
cat("  season uses climatology for 1-17 October, which the deliverable does not\n")
cat("  cover. The offset is near-constant, so the ordering is preserved.\n")

# ---- Verdict ----------------------------------------------------------------
cat("\n=== Verdict ===\n")
pass_chinook <- rho >= 0.80 && all(abs(cmp$ratio - 1) < 0.15)
pass_steel   <- sh_rho >= 0.80

cat(sprintf("  Chinook metric tracks published : %s\n", if (pass_chinook) "PASS" else "FAIL"))
cat(sprintf("  Steelhead ordering preserved    : %s\n", if (pass_steel) "PASS" else "FAIL"))

cat("\nIMPORTANT: the engine's composite score is NOT comparable with the published\n")
cat("composite. Both the steelhead level and the min-max normalisation differ, so\n")
cat("a new year must be compared against a BASELINE RUN THROUGH THIS SAME ENGINE,\n")
cat("never against the published table. The app does it that way.\n")

if (!pass_chinook || !pass_steel) {
  stop("Scenario engine regression test FAILED -- do not trust the upload path.")
}
cat("\nRegression test passed.\n")
