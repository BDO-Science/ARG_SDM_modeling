# ============================================================================
# Where the project's data comes from
# ============================================================================
# One place recording every external input, so the annual refresh is a matter of
# running a script rather than remembering a procedure.
#
# ---------------------------------------------------------------------------
# HOW TO FILL IN THE SacPAS URLs (one-time, takes about a minute each)
# ---------------------------------------------------------------------------
# SacPAS supports scripted access officially. On each query page below:
#
#   1. Set the query the way you want it (all years, Fall-run, etc.)
#   2. Tick "Generate Query Result Link Only"
#   3. Click "Submit Query"
#   4. Copy the link it produces and paste it in as the `url` here
#
# That link is a stable data endpoint. Leave a url as NA and the refresh script
# will skip the download and keep using the bundled snapshot.
#
# NOTE ON DATA VINTAGE: every SacPAS export carries the notice "data presented
# here are preliminary and subject to revision". That is why this project pins
# dated snapshots and refreshes them deliberately, rather than fetching live at
# runtime. Two people running the same analysis months apart must be able to get
# the same answer.
# ============================================================================

DATA_SOURCES <- list(

  carcass = list(
    name        = "SacPAS Chinook Carcass Survey Detail",
    page        = "https://www.cbr.washington.edu/sacramento/data/query_carcass_detail.html",
    url         = NA_character_,      # <-- paste the generated query link here
    snapshot    = "SalmonCountR/app_data/carcassdet_1752789274_15.csv",
    provider    = "Columbia Basin Research, University of Washington (data courtesy CDFW via CalFish)",
    used_for    = "spawn timing model, redd site assignment, spawn-date distribution",
    preliminary = TRUE
  ),

  escapement = list(
    name        = "SacPAS CDFW GrandTab Adult Escapement",
    page        = "https://www.cbr.washington.edu/sacramento/data/query_adult_grandtab.html",
    url         = NA_character_,
    snapshot    = "SalmonCountR/app_data/grandtab_1752793045_337.csv",
    provider    = "Columbia Basin Research, University of Washington (CDFW GrandTab)",
    used_for    = "calibration target, seed population, fit statistics",
    preliminary = TRUE
  ),

  hci = list(
    name        = "SacPAS Water Year Hydrologic Classification Indices",
    page        = "https://www.cbr.washington.edu/sacramento/data/query_hci.html",
    url         = NA_character_,
    snapshot    = "data_raw/hci_1753379372_198.csv",
    provider    = "Columbia Basin Research, University of Washington",
    used_for    = "water year type labelling",
    preliminary = FALSE
  ),

  water_temp = list(
    name        = "USGS NWIS instantaneous water temperature",
    page        = "https://waterdata.usgs.gov/nwis",
    stations    = c(AveHazel = "11446500", AveWatt = "11446980"),
    parameter   = "00010",
    url         = "fetched via dataRetrieval::readNWISdata()",
    snapshot    = "embedded in SalmonCountR/app_data/env_ext_list.rds",
    provider    = "U.S. Geological Survey",
    used_for    = "observed temperature history, climatology outside the modelled window",
    preliminary = FALSE
  ),

  scenario_temps = list(
    name        = "CE-QUAL-W2 power bypass scenario temperatures",
    page        = NA_character_,
    url         = NA_character_,
    snapshot    = "data_raw/SDM Power Bypass Temperature Modeling Results.xlsx",
    provider    = "Reclamation temperature modelling team (emailed deliverable)",
    used_for    = "the nine bypass alternatives under each meteorological year",
    preliminary = FALSE,
    note        = "Not machine-retrievable. Arrives as a spreadsheet; upload it through the app's 'Add a Year' tab, or drop it in data_raw/ for the annual refresh."
  ),

  wua = list(
    name        = "PHABSIM weighted usable area",
    url         = NA_character_,
    snapshot    = "SalmonCountR/app_data/american_river_instream.rds",
    provider    = "Reclamation",
    used_for    = "spawning carrying capacity as a function of flow",
    preliminary = FALSE,
    note        = "Static. Changes only if the habitat study is redone."
  ),

  cwt = list(
    name        = "Coded wire tag releases and returns",
    snapshot    = "SalmonCountR/app_data/SAR LAR Releases.xlsx",
    provider    = "RMIS, hand-assembled",
    used_for    = "smolt-to-adult return ratio",
    preliminary = FALSE,
    note        = paste("Hand-built. The sar_percent column is a copy of sar, NOT a",
                        "percentage -- see OUTSTANDING_ITEMS.md G3.")
  ),

  hydropower = list(
    name        = "Hydropower revenue loss by alternative",
    snapshot    = "hard-coded in analysis/mcda.R and SalmonCountR/app.R",
    provider    = "internal valuation",
    used_for    = "the hydropower objective in the MCDA",
    preliminary = FALSE,
    note        = "Static per scenario. If a new year changes the valuation, update both places."
  )
)

#' Path to a source's bundled snapshot, resolved from the repo root
source_snapshot_path <- function(key) {
  s <- DATA_SOURCES[[key]]$snapshot
  if (is.null(s) || !grepl("/", s)) return(NA_character_)
  here::here(s)
}

#' Which sources are configured for automatic download
configured_sources <- function() {
  names(Filter(function(s) !is.null(s$url) && !is.na(s$url) &&
                 grepl("^https?://", s$url), DATA_SOURCES))
}
