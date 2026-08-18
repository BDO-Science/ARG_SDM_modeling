# ============================================================================
# Annual data refresh
# ============================================================================
# Job B of output/APP_DATA_UPDATE_OPTIONS.md. Run this ONCE A YEAR, before the
# season, by someone who can judge whether the new model fit is acceptable.
#
# Why a script and not a button. Adding a year of carcass and escapement data
# changes the calibration -- the fitted smolt-to-adult return rate and rearing
# survival both move. Deciding that a new fit is good enough to base a
# recommendation on is a scientific judgement. This script does the mechanical
# work and lays out the evidence; a person makes the call.
#
# It is deliberately SAFE BY DEFAULT: nothing in app_data/ is overwritten unless
# you pass --apply. A dry run tells you exactly what would change.
#
#   Rscript analysis/refresh_data_year.R              # dry run, report only
#   Rscript analysis/refresh_data_year.R --apply      # write the new snapshots
#
# After --apply you must still:
#   1. re-run SalmonCountR/precompute.R start to finish in a clean session
#   2. re-run analysis/calibration_fit_statistics.R and look at the fit
#   3. re-run analysis/test_scenario_engine.R and confirm it still passes
#   4. regenerate the manuscript figures if any published number moved
#
# Output: output/data_refresh_report_<date>.md
#         SalmonCountR/app_data/data_vintage.rds
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr)
  library(lubridate); library(here)
})

source(here("analysis", "data_sources.R"))

args  <- commandArgs(trailingOnly = TRUE)
APPLY <- "--apply" %in% args
today <- Sys.Date()

log_lines <- character()
say <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "")
  log_lines <<- c(log_lines, msg)
}

say("# Annual data refresh -- %s", format(today, "%d %B %Y"))
say("")
say(if (APPLY) "**Mode: APPLY.** Snapshots will be replaced."
    else "**Mode: dry run.** Nothing will be written. Re-run with `--apply` to commit.")
say("")

# ---- 1. What is configured for automatic download ---------------------------
auto <- configured_sources()
say("## Sources")
say("")
say("| Source | Automatic | Snapshot |")
say("|---|---|---|")
for (k in names(DATA_SOURCES)) {
  s <- DATA_SOURCES[[k]]
  say("| %s | %s | `%s` |", s$name,
      if (k %in% auto) "yes" else "no -- manual",
      if (is.null(s$snapshot)) "-" else s$snapshot)
}
say("")

if (!length(auto)) {
  say("> No SacPAS query links are configured yet, so nothing can be downloaded")
  say("> automatically. See the instructions at the top of `analysis/data_sources.R`:")
  say("> tick **Generate Query Result Link Only** on each SacPAS query page, then")
  say("> paste the link it produces into `DATA_SOURCES`. Until then this script")
  say("> reports on the bundled snapshots only.")
  say("")
}

# ---- 2. Current vintage of each snapshot ------------------------------------
say("## Current snapshots")
say("")
say("| Source | File | Modified | Size |")
say("|---|---|---|---|")
vintage_rows <- purrr::map_dfr(names(DATA_SOURCES), function(k) {
  p <- source_snapshot_path(k)
  if (is.na(p) || !file.exists(p)) {
    say("| %s | %s | -- | not on disk |", DATA_SOURCES[[k]]$name,
        ifelse(is.na(p), "n/a", basename(p)))
    return(tibble(source = k, file = NA_character_, modified = as.POSIXct(NA), bytes = NA_integer_))
  }
  fi <- file.info(p)
  say("| %s | `%s` | %s | %.0f KB |", DATA_SOURCES[[k]]$name, basename(p),
      format(fi$mtime, "%Y-%m-%d"), fi$size / 1024)
  tibble(source = k, file = basename(p), modified = fi$mtime, bytes = as.integer(fi$size))
})
say("")

# ---- 3. Download, if configured ---------------------------------------------
downloaded <- list()
if (length(auto)) {
  say("## Downloads")
  say("")
  for (k in auto) {
    s   <- DATA_SOURCES[[k]]
    dst <- file.path(tempdir(), sprintf("%s_%s.csv", k, format(today, "%Y%m%d")))
    ok  <- tryCatch({
      utils::download.file(s$url, dst, quiet = TRUE, mode = "wb"); TRUE
    }, error = function(e) { say("- **%s**: download failed -- %s", s$name, conditionMessage(e)); FALSE })
    if (ok && file.exists(dst) && file.info(dst)$size > 0) {
      downloaded[[k]] <- dst
      say("- **%s**: downloaded, %.0f KB", s$name, file.info(dst)$size / 1024)
    } else if (ok) {
      say("- **%s**: download produced an empty file; keeping the snapshot", s$name)
    }
  }
  say("")
}

# ---- 4. What changed --------------------------------------------------------
say("## What changed")
say("")

compare_csv <- function(key, label, id_col = NULL) {
  old_p <- source_snapshot_path(key)
  new_p <- downloaded[[key]]
  if (is.null(new_p)) {
    say("- **%s**: no new download, nothing to compare.", label); return(invisible(NULL))
  }
  old <- suppressWarnings(try(read_csv(old_p, show_col_types = FALSE), silent = TRUE))
  new <- suppressWarnings(try(read_csv(new_p, show_col_types = FALSE), silent = TRUE))
  if (inherits(old, "try-error") || inherits(new, "try-error")) {
    say("- **%s**: could not be compared (one of the files would not parse).", label); return(invisible(NULL))
  }
  say("- **%s**: %d rows -> %d rows (%+d).", label, nrow(old), nrow(new), nrow(new) - nrow(old))
  if (!is.null(id_col) && id_col %in% names(old) && id_col %in% names(new)) {
    added <- setdiff(unique(new[[id_col]]), unique(old[[id_col]]))
    if (length(added)) say("  New values of `%s`: %s.", id_col, paste(sort(added), collapse = ", "))
  }
}

compare_csv("carcass",    "Carcass survey detail", "surveydate")
compare_csv("escapement", "GrandTab escapement",   "End Year of Monitoring Period")
compare_csv("hci",        "Hydrologic index")
say("")

# ---- 5. Escapement years available vs used ----------------------------------
esc_p <- source_snapshot_path("escapement")
if (file.exists(esc_p)) {
  esc <- suppressWarnings(read_csv(esc_p, show_col_types = FALSE)) %>%
    mutate(year = suppressWarnings(parse_number(`End Year of Monitoring Period`))) %>%
    filter(!is.na(year), !is.na(`Population Estimate`))
  used_to <- 2024   # the calibration window baked into precompute.R
  say("## Calibration window")
  say("")
  say("- Escapement available through **%d**.", max(esc$year))
  say("- The model calibrates on **2011-%d**.", used_to)
  if (max(esc$year) > used_to) {
    say("- **%d additional year(s) are available and unused.** Extending the window",
        max(esc$year) - used_to)
    say("  means editing `real_years` in `SalmonCountR/precompute.R` and re-running it.")
    say("  Expect the calibrated SAR and `rear_surv` to move.")
  } else {
    say("- No new escapement years since the last calibration.")
  }
  say("")
}

# ---- 6. Spawn-timing model --------------------------------------------------
say("## Spawn-timing model")
say("")
stm_p <- here("SalmonCountR", "app_data", "spawn_timing_model.rds")
if (file.exists(stm_p)) {
  stm <- readRDS(stm_p)
  say("- Fitted %s on %d carcass observations, brood years %d-%d.",
      format(stm$fitted, "%Y-%m-%d"), stm$n_obs, min(stm$brood_years), max(stm$brood_years))
  say("- Coefficients: Oct_std %+0.4f, Nov_std %+0.4f.",
      stm$beta[["Oct_std"]], stm$beta[["Nov_std"]])
  say("- Refit with `Rscript analysis/build_spawn_timing_model.R` after new carcass data lands.")
} else {
  say("- **Missing.** Run `Rscript analysis/build_spawn_timing_model.R`.")
  say("  Without it the app's 'Add a Year' tab is disabled.")
}
say("")

# ---- 7. Apply ---------------------------------------------------------------
say("## Actions")
say("")
if (APPLY && length(downloaded)) {
  for (k in names(downloaded)) {
    dst <- source_snapshot_path(k)
    bak <- paste0(dst, ".bak-", format(today, "%Y%m%d"))
    file.copy(dst, bak, overwrite = TRUE)
    file.copy(downloaded[[k]], dst, overwrite = TRUE)
    say("- Replaced `%s` (previous version kept as `%s`).", basename(dst), basename(bak))
  }
} else if (length(downloaded)) {
  say("- Dry run: %d snapshot(s) would be replaced. Re-run with `--apply`.", length(downloaded))
} else {
  say("- Nothing to apply.")
}
say("")

say("### Still to do by hand")
say("")
say("1. Re-run `SalmonCountR/precompute.R` start to finish in a clean session.")
say("   The redd draw is seeded but only reproducible on a full run.")
say("2. Re-run `analysis/build_spawn_timing_model.R`.")
say("3. Re-run `analysis/calibration_fit_statistics.R` and **look at the fit** before")
say("   accepting the new calibration. As of 2026-07-28 it was poor (R2 0.13,")
say("   Nash-Sutcliffe -0.59); see docs/spawn-timing.md D1.")
say("4. Re-run `analysis/test_scenario_engine.R` and confirm it still passes.")
say("5. Regenerate figures if any published number moved.")
say("6. Drop the new CE-QUAL-W2 deliverable in `data_raw/` so the app's baseline")
say("   comparison uses the current year.")
say("")

# ---- 8. Vintage stamp -------------------------------------------------------
#' Current commit, read straight from .git
#'
#' Deliberately does NOT shell out to git. This repository lives on a OneDrive
#' path where git subprocesses can block on a lock held by the sync client, and
#' a hung child process would take the whole refresh with it.
current_commit <- function() {
  tryCatch({
    gd <- here(".git")
    if (!dir.exists(gd)) return(NA_character_)
    head_txt <- readLines(file.path(gd, "HEAD"), warn = FALSE)[1]
    if (grepl("^ref: ", head_txt)) {
      ref <- sub("^ref: ", "", head_txt)
      rp  <- file.path(gd, ref)
      if (file.exists(rp)) return(substr(readLines(rp, warn = FALSE)[1], 1, 7))
      pk <- file.path(gd, "packed-refs")          # ref may only be in packed-refs
      if (file.exists(pk)) {
        hit <- grep(paste0(" ", ref, "$"), readLines(pk, warn = FALSE), value = TRUE)
        if (length(hit)) return(substr(sub(" .*$", "", hit[1]), 1, 7))
      }
      return(NA_character_)
    }
    substr(head_txt, 1, 7)                         # detached HEAD
  }, error = function(e) NA_character_)
}

vintage <- list(
  refreshed   = Sys.time(),
  applied     = APPLY,
  sources     = vintage_rows,
  downloaded  = names(downloaded),
  git_commit  = current_commit(),
  note        = "Written by analysis/refresh_data_year.R"
)

if (APPLY) {
  saveRDS(vintage, here("SalmonCountR", "app_data", "data_vintage.rds"))
  say("Wrote `SalmonCountR/app_data/data_vintage.rds`.")
}

dir.create(here("output"), showWarnings = FALSE)
report_path <- here("output", sprintf("data_refresh_report_%s.md", format(today, "%Y-%m-%d")))
writeLines(log_lines, report_path)
cat("\nReport written to", report_path, "\n")
