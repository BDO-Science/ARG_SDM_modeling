# deploy.R --------------------------------------------------------------------
# One-command deploy for SalmonCountR, plus a check that proves the bundle works
# BEFORE it is uploaded.
#
# USAGE. From the repo root or from SalmonCountR/, in a fresh R session:
#
#     source("SalmonCountR/deploy.R")
#     arg_verify_bundle()                              # optional, see below
#     arg_deploy(account = "reclamation-bdo-science")
#
# arg_deploy() runs arg_verify_bundle() itself and refuses to upload if it
# fails, so the middle line is only for checking without deploying.
#
# WHY THIS FILE EXISTS. Three things used to make deploying this app a manual
# chore. The first is fixed in the code; the other two are handled here.
#
#   1. FIXED IN CODE, NOT HERE. global.R used here::here("SalmonCountR", ...).
#      here() ignores the working directory and walks up for a .git/.Rproj
#      marker. A deployed bundle has no marker, so here() fell back to the
#      working directory -- already the app directory -- and the "SalmonCountR"
#      prefix pointed one level too deep. That is why the folder had to be
#      rearranged before every deploy. See arg_app_dir() in years.R.
#
#   2. precompute.R sits in the app directory. rsconnect scans EVERY .R file in
#      the bundle for library() calls, so shipping precompute.R made the server
#      install furrr, ordinal, MASS, ggridges, readxl and here. The app never
#      runs any of them, and two of them compile from source -- slow, and a
#      failure mode that has nothing to do with the app.
#
#   3. 16 of app_data's 25 MB are precompute intermediates (sim_future.rds
#      alone is 13 MB) that no runtime code reads.
#
# The exclusion below is a BLACKLIST on purpose. A whitelist would silently drop
# any new data file someone adds; a blacklist ships new files by default and
# only omits the ones named here as known-unused.
#
# IF THE APP EVER STARTS READING ONE OF THESE FILES, delete it from the list.
# arg_verify_bundle() will catch the mistake either way -- it builds the app
# from a copy containing only what would actually be uploaded.
# -----------------------------------------------------------------------------

# Locate the app directory the same way the app itself does, so this script can
# be sourced from the repo root or from inside SalmonCountR/.
local({
  probe <- if (dir.exists("app_data")) "." else
           if (dir.exists(file.path("SalmonCountR", "app_data"))) "SalmonCountR" else
           stop("Cannot find the SalmonCountR app directory from '", getwd(), "'.",
                call. = FALSE)
  source(file.path(probe, "years.R"))
})

#' Files that exist in the app directory but that no runtime code reads.
#'
#' Verified against every readRDS/read.csv/read_excel call in app.R, global.R,
#' years.R and functions.R, plus the ARG_YEAR_FILES registry.
ARG_DEPLOY_EXCLUDE <- c(
  # Not part of the running app.
  "precompute.R",   # see note 2 above -- this is the important one
  "deploy.R",       # would make rsconnect a server dependency
  "todo.md",
  ".rscignore",

  # precompute.R intermediates. Sizes are why they are listed.
  file.path("app_data", c(
    "sim_future.rds",                  # 13 MB
    "sim_redds.rds",                   # 2.1 MB
    "carcassdet_1752789274_15.csv",    # 960 KB
    "simulate_variant.rds",            # 56 KB
    "generate_SAR_vec.rds",            # 52 KB
    "SAR LAR Releases.xlsx",           # 24 KB
    "grandtab_1752793045_337.csv",
    "american_river_instream.rda",     # superseded by the .rds
    "swing_extreme_combos.rds",
    "spawn_timing_model.rds",
    "spawn_dates_vec.rds",
    "spawn_dates_by_env.rds",
    "spawn_dates.rds",
    "calib_pred_by_variant.rds",
    "calib_results.rds",
    "rear_surv_lookup.rds",
    "nonsalmon_objectives.csv",
    "steelhead_objective.csv",
    "S_seed.rds"
  ))
)

#' The exact file list that will be uploaded, relative to the app directory.
arg_bundle_files <- function(app_dir = arg_app_dir()) {
  all <- list.files(app_dir, recursive = TRUE, all.files = FALSE, no.. = TRUE)
  # rsconnect's own deployment records: never part of the app.
  all <- all[!startsWith(all, "rsconnect/")]
  setdiff(all, ARG_DEPLOY_EXCLUDE)
}

#' Build the app from a copy of the bundle, in a separate R process.
#'
#' This is the check that matters. It copies ONLY the files that would be
#' uploaded into a temporary directory with no .git/.Rproj above it -- the same
#' conditions as the server -- and confirms app.R evaluates and the Shiny app
#' object constructs. A missing data file or a path that only resolves locally
#' fails here, on the desk, instead of after a five-minute upload.
#'
#' @return TRUE invisibly on success; stops with the child process output on
#'   failure.
arg_verify_bundle <- function(app_dir = arg_app_dir(), quiet = FALSE) {
  files <- arg_bundle_files(app_dir)

  # Cheap check first: everything the app declares it needs must be present.
  need <- c(unname(ARG_YEAR_FILES), "data_vintage.rds")
  missing_year <- need[!file.exists(file.path(app_dir, "app_data", need))]
  if (length(missing_year)) {
    stop("app_data is missing required files: ",
         paste(missing_year, collapse = ", "), call. = FALSE)
  }

  sim <- file.path(tempdir(), paste0("argbundle_", as.integer(Sys.time())))
  on.exit(unlink(sim, recursive = TRUE, force = TRUE), add = TRUE)

  for (f in files) {
    dest <- file.path(sim, f)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(app_dir, f), dest, overwrite = TRUE)
  }

  if (!quiet) {
    mb <- sum(file.size(file.path(app_dir, files)), na.rm = TRUE) / 1024^2
    message(sprintf("Verifying %d files (%.1f MB) in %s", length(files), mb, sim))
  }

  rscript <- file.path(R.home("bin"),
                       if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

  # The child must run with the bundle as its working directory -- that is the
  # whole point, since it is what Shiny sets on the server. Have the CHILD
  # setwd() rather than doing it around the call: leaving this process parked
  # inside a directory that on.exit() then tries to unlink() is a good way to
  # crash R on Windows. (system2()'s env= is ignored on Windows and is not
  # needed anyway -- a child of the same R.home() inherits the library paths.)
  sim_r <- normalizePath(sim, winslash = "/", mustWork = TRUE)

  out <- suppressWarnings(system2(
    rscript,
    c("-e", shQuote(paste(
      sprintf("setwd('%s')", sim_r),
      "app <- shiny::shinyAppFile('app.R')",
      "stopifnot(inherits(app, 'shiny.appobj'))",
      "cat('BUNDLE_OK\n')", sep = "; "))),
    stdout = TRUE, stderr = TRUE
  ))

  if (!any(grepl("BUNDLE_OK", out, fixed = TRUE))) {
    stop("Bundle verification FAILED. The app does not build from the files ",
         "that would be uploaded.\n\n", paste(out, collapse = "\n"), call. = FALSE)
  }
  if (!quiet) message("Bundle OK -- the app builds from the deploy set alone.")
  invisible(TRUE)
}

#' Regenerate .rscignore from ARG_DEPLOY_EXCLUDE.
#'
#' arg_deploy() passes the file list explicitly, so it does not need .rscignore.
#' That file is the safety net for the OTHER way people deploy: the Publish
#' button in RStudio, which calls deployApp() without appFiles. Generating it
#' from the same vector means the two cannot drift apart.
#'
#' Re-run this after editing ARG_DEPLOY_EXCLUDE. Needs rsconnect >= 0.8.17 to
#' have any effect; on older versions the Publish button simply ships more than
#' it needs to, which is slower but still correct.
arg_write_rscignore <- function(app_dir = arg_app_dir()) {
  path <- file.path(app_dir, ".rscignore")
  writeLines(c(
    "# Generated by deploy.R -- edit ARG_DEPLOY_EXCLUDE there, then run",
    "# arg_write_rscignore(). Do not hand-edit this file.",
    "#",
    "# Keeps precompute.R and its intermediates out of the deployed bundle.",
    "# precompute.R matters most: rsconnect scans every .R file in the bundle",
    "# for library() calls, so shipping it makes the server install furrr,",
    "# ordinal, MASS, ggridges, readxl and here, none of which the app runs.",
    "",
    ARG_DEPLOY_EXCLUDE
  ), path)
  message("Wrote ", path)
  invisible(path)
}

#' Verify, then deploy.
#'
#' @param account shinyapps.io account name, e.g. "reclamation-bdo-science".
#'   Must already be configured -- see the setup note in DEPLOY.md.
#' @param app_name Name of the app on the server.
#' @param verify Set FALSE only to skip the pre-upload check deliberately.
arg_deploy <- function(account,
                       app_name = "SalmonCountR",
                       app_dir  = arg_app_dir(),
                       verify   = TRUE) {
  if (!requireNamespace("rsconnect", quietly = TRUE)) {
    stop("Install rsconnect first: install.packages('rsconnect')", call. = FALSE)
  }
  if (missing(account)) {
    stop("Give the account explicitly, e.g. ",
         'arg_deploy(account = "reclamation-bdo-science")', call. = FALSE)
  }
  if (verify) arg_verify_bundle(app_dir)

  rsconnect::deployApp(
    appDir      = app_dir,
    appFiles    = arg_bundle_files(app_dir),
    appName     = app_name,
    account     = account,
    forceUpdate = TRUE
  )
}
