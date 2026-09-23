# Builds the daily temperature series that precompute.R and the app read, from a
# CE-QUAL-W2 temperature deliverable plus observed USGS gauge data.
#
#   deliverable .xlsx  ->  env_ext_list.rds, df_all.rds, alt_key.rds  (scenario folder)
#
# Settings (environment variables):
#   ARG_TEMP_FILE     the deliverable workbook
#                     default data_raw/SDM Power Bypass Temperature Modeling Results.xlsx
#   ARG_APP_DATA_DIR  where to write the .rds files -- use the same folder for
#                     precompute.R. default SalmonCountR/app_data (the published
#                     results; point it at a subfolder for a new scenario)
#   ARG_OBS_END       the decision date: observed gauge temperatures are used
#                     through it, scenario temperatures after. default: the day
#                     before the deliverable starts (2025-09-21 for the 2025 file)
#   ARG_OBS_REFRESH   "1" to download the observed record again instead of
#                     reusing the frozen copy (observed_temps.rds) in the folder
#   ARG_SCENARIO_START  first day (MM-DD) the scenario temperatures are used.
#                     default: the deliverable's first date. The published 2025
#                     results used "10-18"; set that to reproduce them.
#   ARG_ATSP_BASE     when the workbook mixes ATSP schedules, the one that keeps
#                     the plain codes (NB, PB1, ...); the others get a suffix
#                     (NB-38, PB1-38, ...). default: the most common one.
#
# To rebuild the published 2025 temperature series exactly:
#   Sys.setenv(ARG_SCENARIO_START = "10-18"); source("analysis/temperature_data.R")
#
# The workbook layout, and what the model does with it, is described in
# docs/new-temperature-scenario.md. Needs an internet connection (USGS NWIS).

# Load necessary libraries
library(dataRetrieval) # for readNWISdata()
library(tidyverse)
library(readxl)
library(writexl)
library(here)

.env_path <- function(var, default) {
  p <- Sys.getenv(var, default)
  if (!grepl("^([A-Za-z]:)?[/\\\\]", p)) p <- here(p)
  p
}
temp_file    <- .env_path("ARG_TEMP_FILE",
                          file.path("data_raw", "SDM Power Bypass Temperature Modeling Results.xlsx"))
app_data_dir <- .env_path("ARG_APP_DATA_DIR", file.path("SalmonCountR", "app_data"))
if (!dir.exists(app_data_dir)) dir.create(app_data_dir, recursive = TRUE)
is_published_dir <- normalizePath(app_data_dir, winslash = "/") ==
  normalizePath(here("SalmonCountR", "app_data"), winslash = "/")

# The meteorological years, one sheet each. The rest of the pipeline (met-year
# weight sliders in the app, the equal 0.25 weighting in precompute.R) assumes
# these four; a deliverable with different years is a code change, not a data
# change.
HYDRO_YEARS <- c("2011", "2014", "2017", "2020")
SITES       <- c("AveWatt", "AveHazel")

# ---- Workbook layout ---------------------------------------------------------
# Scenarios are DISCOVERED from row 1 rather than counted: every non-empty cell
# in row 1 from column C onward starts a scenario block, and the block runs to
# the column before the next label. Within a block, row 2 names the sites; the
# AveWatt and AveHazel columns are used and anything else (the 2026 deliverable
# adds an AveFol column per scenario; the 2025 one has a Target Temp column at
# the end) is ignored. So a workbook with six, nine or eleven scenarios reads
# the same way, and appending a scenario is a data change.
read_sheet_layout <- function(input_path, sheet) {
  hdr <- read_excel(input_path, sheet = sheet, col_names = FALSE, n_max = 2,
                    .name_repair = "minimal")
  r1 <- trimws(as.character(unlist(hdr[1, ])))
  r2 <- trimws(as.character(unlist(hdr[2, ])))
  pos <- which(!is.na(r1) & nzchar(r1) & seq_along(r1) >= 3)
  if (!length(pos)) {
    stop("Sheet ", sheet, ", row 1 has no scenario labels from column C onward.")
  }
  ends <- c(pos[-1] - 1L, length(r1))
  cols <- map2(pos, ends, function(a, b) {
    rng <- a:b
    vapply(SITES, function(s) {
      hit <- rng[which(!is.na(r2[rng]) & r2[rng] == s)]
      if (length(hit)) hit[1] else NA_integer_
    }, integer(1))
  })
  tibble(label = r1[pos], AveWatt = map_int(cols, "AveWatt"), AveHazel = map_int(cols, "AveHazel"))
}

# Short code for each scenario label, e.g.
#   "No Bypass"            -> NB        "Scenario 2b"           -> PB2b
#   "ATSP 41 - No Bypass"  -> NB        "ATSP 38 - Scenario 1"  -> PB1-38
# The ATSP suffix appears only when the workbook mixes schedules; the base
# schedule (ARG_ATSP_BASE, default the most common) keeps the plain codes so
# they line up with earlier years.
alt_codes_from_labels <- function(labels) {
  atsp <- str_match(labels, regex("ATSP\\s*(\\d+)", ignore_case = TRUE))[, 2]
  core <- str_trim(str_remove(labels, regex("^\\s*ATSP\\s*\\d+\\s*[-:–]\\s*", ignore_case = TRUE)))
  scen <- str_match(core, regex("^scenario\\s*(\\S+)$", ignore_case = TRUE))[, 2]
  code <- case_when(
    str_detect(core, regex("^no\\s*bypass$", ignore_case = TRUE)) ~ "NB",
    !is.na(scen)                                                  ~ paste0("PB", scen),
    TRUE                                                          ~ make.names(core)
  )
  schedules <- unique(na.omit(atsp))
  if (length(schedules) > 1) {
    base <- Sys.getenv("ARG_ATSP_BASE", "")
    if (!nzchar(base)) {
      counts <- table(atsp)
      top    <- names(counts)[counts == max(counts)]
      base   <- if (length(top) == 1) top else atsp[!is.na(atsp)][1]
    }
    if (!base %in% schedules) {
      stop("ARG_ATSP_BASE = ", base, " but the workbook has ATSP ",
           paste(schedules, collapse = ", "))
    }
    code <- ifelse(!is.na(atsp) & atsp != base, paste0(code, "-", atsp), code)
    message("ATSP schedules ", paste(schedules, collapse = ", "), " in the workbook; ",
            base, " keeps the plain codes.")
  }
  if (anyDuplicated(code)) {
    stop("Scenario labels map to duplicate codes: ",
         paste(code[duplicated(code)], collapse = ", "), ". Labels: ",
         paste(labels, collapse = " | "))
  }
  list(code = code, atsp = atsp)
}

# Checks the deliverable before anything is read from it, and returns the
# layout (scenario labels and the columns each reads from). Columns are read
# by position, so a relabelled or reordered file would otherwise run to
# completion and give wrong answers without any error.
check_temperature_file <- function(input_path, hydro_years) {
  if (!file.exists(input_path)) stop("Temperature file not found: ", input_path)
  sheets  <- excel_sheets(input_path)
  missing <- setdiff(hydro_years, sheets)
  if (length(missing)) {
    stop("Sheet(s) ", paste(missing, collapse = ", "), " not found in ", basename(input_path),
         ". Expected one sheet per meteorological year, named ",
         paste(hydro_years, collapse = ", "), ". Sheets found: ",
         paste(sheets, collapse = ", "))
  }
  layout <- NULL
  for (hy in hydro_years) {
    lay <- read_sheet_layout(input_path, hy)
    bad <- lay %>% filter(is.na(AveWatt) | is.na(AveHazel))
    if (nrow(bad)) {
      stop("Sheet ", hy, ": scenario '", bad$label[1], "' has no ",
           paste(SITES[is.na(c(bad$AveWatt[1], bad$AveHazel[1]))], collapse = " or "),
           " column. Row 2 must name the site columns AveWatt and AveHazel under each ",
           "row-1 scenario label.")
    }
    if (is.null(layout)) {
      layout <- lay
    } else if (!identical(lay$label, layout$label)) {
      stop("Sheet ", hy, " lists different scenarios (", paste(lay$label, collapse = ", "),
           ") from sheet ", hydro_years[1], " (", paste(layout$label, collapse = ", "),
           "). Every met-year sheet must have the same scenarios in the same order.")
    }
    dat   <- suppressMessages(read_excel(input_path, sheet = hy, skip = 2, col_names = FALSE,
                                         .name_repair = "minimal"))
    dates <- tryCatch(suppressWarnings(as.Date(dat[[1]])),
                      error = function(e) rep(as.Date(NA), nrow(dat)))
    if (all(is.na(dates))) {
      stop("Sheet ", hy, ": column A (from row 3) does not contain dates. ",
           "Store them as Excel dates, not text.")
    }
    yrs <- as.integer(format(dates[!is.na(dates)], "%Y"))
    if (any(yrs < 1950 | yrs > 2200)) {
      stop("Sheet ", hy, ": column A has dates in year ", min(yrs), ", which usually means ",
           "they were stored as text (for example day/month/year) and misread. ",
           "Store them as Excel dates.")
    }
    used  <- c(lay$AveWatt, lay$AveHazel)
    temps <- suppressWarnings(as.numeric(unlist(dat[, used])))
    if (all(is.na(temps))) stop("Sheet ", hy, ": temperature columns are not numeric.")
    rng <- range(temps, na.rm = TRUE)
    if (rng[1] < 0 || rng[2] > 30) {
      stop(sprintf(paste0("Sheet %s: temperatures run %.1f to %.1f. Expected degrees ",
                          "Celsius (roughly 4-30); values above 30 usually mean Fahrenheit."),
                   hy, rng[1], rng[2]))
    }
    message(sprintf("Sheet %s: %s to %s, %d days, %d scenarios.",
                    hy, format(min(dates, na.rm = TRUE)), format(max(dates, na.rm = TRUE)),
                    sum(!is.na(dates)), nrow(lay)))
    # Spawning, incubation and the Oct/Nov means that drive spawn timing all
    # fall in October-December; a file that misses them has nothing to act on.
    if (!any(month(dates) %in% 10:12, na.rm = TRUE)) {
      stop("Sheet ", hy, " has no dates in October-December, when the model ",
           "uses the scenario temperatures.")
    }
  }
  layout
}

# Reads the deliverable, reformats it as one sheet per run (scenario x met
# year, numbered met-year-major), and writes the alternative key.
prepare_temperature_file <- function() {
  input_path <- temp_file
  sheets     <- excel_sheets(input_path)
  layout     <- check_temperature_file(input_path, HYDRO_YEARS)
  codes      <- alt_codes_from_labels(layout$label)
  message("Scenarios: ", paste(sprintf("%s = %s", codes$code, layout$label), collapse = "; "))

  all_sheets <- list()

  # Copy metadata sheets, when the deliverable has them
  for (meta in c("Scenario Summary", "Flow")) {
    if (meta %in% sheets) all_sheets[[meta]] <- read_excel(input_path, sheet = meta)
  }

  # The alternative key: one row per run. Saved as alt_key.rds for precompute.R
  # and the app, and as the intermediate workbook's metadata sheet.
  alt_key <- expand_grid(met_year = HYDRO_YEARS, k = seq_len(nrow(layout))) %>%
    mutate(env   = row_number(),
           alt   = codes$code[k],
           label = layout$label[k],
           atsp  = codes$atsp[k]) %>%
    select(env, alt, label, met_year, atsp)
  all_sheets[["metadata"]] <- alt_key %>%
    transmute(Alternative = env, Code = alt, Scenario = label, Hydro_Year = met_year, ATSP = atsp) %>%
    as.data.frame()

  # Extract data for each run
  for (i in seq_len(nrow(alt_key))) {
    hy  <- alt_key$met_year[i]
    k   <- alt_key$env[i] - (match(hy, HYDRO_YEARS) - 1L) * nrow(layout)
    if (is.null(all_sheets[[hy]])) {
      all_sheets[[hy]] <- suppressMessages(
        read_excel(input_path, sheet = hy, skip = 2, col_names = FALSE, .name_repair = "minimal"))
    }
    df_hy <- all_sheets[[hy]]
    all_sheets[[as.character(alt_key$env[i])]] <- data.frame(
      Date     = as.Date(df_hy[[1]]),
      AveWatt  = suppressWarnings(as.numeric(df_hy[[layout$AveWatt[k]]])),
      AveHazel = suppressWarnings(as.numeric(df_hy[[layout$AveHazel[k]]]))
    )
  }
  all_sheets <- all_sheets[setdiff(names(all_sheets), HYDRO_YEARS)]

  # The published run keeps its intermediate workbook in data_raw/; a scenario
  # run keeps it in its own folder so the published one is not overwritten.
  output_path <- if (is_published_dir) {
    here("data_raw", "ARG_LAR_TempModeling_alternatives.xlsx")
  } else {
    file.path(app_data_dir, "temperature_alternatives.xlsx")
  }
  write_xlsx(all_sheets, output_path)
  saveRDS(alt_key, file.path(app_data_dir, "alt_key.rds"))
  print(paste("Created", nrow(alt_key), "runs (", nrow(layout), "scenarios x",
              length(HYDRO_YEARS), "met years ) in", output_path))
  print(paste("Saved alt_key.rds to", app_data_dir))
  list(path = output_path, key = alt_key)
}

# --- SCRIPT EXECUTION STARTS HERE ---

# 1) PREPARE DATA: Run the function to reformat the alternatives file
prepared  <- prepare_temperature_file()
xlsx_path <- prepared$path
alt_key   <- prepared$key

# 2) SET PARAMETERS
# obs_end is the DECISION DATE: observed gauge temperatures are used through
# it, the scenario temperatures after it. Before the decision every alternative
# is the same river, so observed data is the right input up to that day.
# Default: the day before the deliverable starts (2025-09-21 for the 2025 file,
# which is what the published run used).
deliverable_start <- min(do.call(c, lapply(HYDRO_YEARS, function(s)
  as.Date(suppressMessages(read_excel(temp_file, sheet = s, skip = 2, col_names = FALSE,
                                      .name_repair = "minimal"))[[1]]))), na.rm = TRUE)
obs_start <- as.Date("2011-09-01")
obs_end   <- as.Date(Sys.getenv("ARG_OBS_END", format(deliverable_start - 1)))
message(sprintf("Observed temperatures through %s (decision date); scenario temperatures after.",
                format(obs_end)))
last_wy   <- 2150 # Final water-year to simulate
sim_end   <- as.Date(sprintf("%04d-08-31", last_wy + 1))

# 3) GET OBSERVED GAUGE DATA from NWIS -- frozen at the decision date
# The download is saved as observed_temps.rds in the scenario folder and reused
# on every later run, so the analysis stays frozen as it stood at the decision:
# USGS revises provisional data for months afterwards, and a fresh download
# would quietly change the inputs. Set ARG_OBS_REFRESH = "1" to download again.
stations <- c("11446980", "11446500") # Watt Avenue and Hazel Avenue
param_cd <- "00010" # Temperature parameter code
obs_cache <- file.path(app_data_dir, "observed_temps.rds")
# Reused only if it was downloaded for this same decision date. (A download
# ending on a date also carries a few readings from the next day, because USGS
# timestamps are UTC; those feed the climatology, so a record fetched for a
# later date is not equivalent to one fetched for this date.)
use_cache <- file.exists(obs_cache) && Sys.getenv("ARG_OBS_REFRESH") != "1" && {
  cached <- readRDS(obs_cache)
  identical(attr(cached, "obs_start"), obs_start) && identical(attr(cached, "obs_end"), obs_end)
}

amer_obs <- if (use_cache) {
  message(sprintf("Using the frozen observed record in %s (downloaded %s, decision date %s).",
                  obs_cache, attr(cached, "downloaded"), format(attr(cached, "obs_end"))))
  attr(cached, "obs_start") <- attr(cached, "obs_end") <- attr(cached, "downloaded") <- NULL
  cached
} else map_df(stations, function(stn) {
  readNWISdata(
    sites       = stn,
    parameterCd = param_cd,
    service     = "uv", # instantaneous data
    startDate   = obs_start,
    endDate     = obs_end
  ) %>%
    transmute(
      Date  = as.Date(dateTime),
      site  = recode(site_no, "11446980" = "AveWatt", "11446500" = "AveHazel"),
      temp  = X_00010_00000
    )
}) %>%
  group_by(Date, site) %>%
  summarise(temp = mean(temp, na.rm=TRUE), .groups="drop") %>%
  mutate(
    temp = if_else(site == "AveHazel",
                   pmax(if_else(is.na(temp)|is.nan(temp), 7, temp), 7),
                   pmax(if_else(is.na(temp)|is.nan(temp), 8, temp), 8))
  )
if (!use_cache) {
  frozen <- amer_obs
  attr(frozen, "obs_start")  <- obs_start
  attr(frozen, "obs_end")    <- obs_end
  attr(frozen, "downloaded") <- format(Sys.time(), "%Y-%m-%d %H:%M")
  saveRDS(frozen, obs_cache)
  message("Observed record saved to ", obs_cache)
}
# If the gauge record stops short of the decision date, end the observed block
# where the data ends, so the series has no missing days; scenario temperatures
# (or climatology outside the deliverable) start the day after.
obs_last <- max(amer_obs$Date)
if (obs_last < obs_end) {
  warning(sprintf(paste0("USGS data ends %s, before the decision date %s. Observed ",
                         "temperatures are used through %s and scenario temperatures after."),
                  format(obs_last), format(obs_end), format(obs_last)))
  obs_end <- obs_last
}

# 4) BUILD CLIMATOLOGY: Create a 14-year daily average temperature
clim14 <- amer_obs %>%
  mutate(doy = yday(Date)) %>%
  group_by(site, doy) %>%
  summarize(clim_temp = mean(temp, na.rm = TRUE), .groups = "drop")

# 5) READ ALTERNATIVES: Load the generated runs, in key order
alts <- as.character(alt_key$env)

pred_by_doy <- map_df(alts, function(alt) {
  read_excel(xlsx_path, sheet=alt) %>%
    mutate(Date = as.Date(Date)) %>%
    pivot_longer(
      cols      = starts_with("Ave"),
      names_to  = "site",
      values_to = "temp_alt"
    ) %>%
    mutate(doy = yday(Date), alt = alt) %>%
    select(alt, site, doy, temp_alt)
})

# 6) CREATE FUTURE TIMESERIES: Generate dates from end of observed to end of simulation
future_dates <- tibble(
  Date = seq(obs_end + 1, sim_end, by="day")
) %>% mutate(doy = yday(Date))

# Define the period where forecast data is used: from the deliverable's first
# day to Dec 31. Days after the deliverable ends (Dec 1-31 for the 2025 file)
# have no scenario value and fall back to climatology below.
#
# THE PUBLISHED 2025 RESULTS USED A FIXED OCT 18 START (day 291), left over from
# an earlier deliverable. It replaced the 2025 deliverable's Sep 22 - Oct 17
# values with climatology, hiding the early cooling of PB1, PB2, PB2b and PB6.
# The 2025 analysis is kept as published; the corrected start applies from the
# next deliverable. To reproduce the published env_ext_list.rds / df_all.rds,
# set ARG_SCENARIO_START = "10-18".
.scen_start <- Sys.getenv("ARG_SCENARIO_START", "")
threshold_start <- if (nzchar(.scen_start)) {
  yday(as.Date(paste0("2021-", .scen_start)))   # MM-DD, non-leap reference year
} else {
  min(pred_by_doy$doy[!is.na(pred_by_doy$temp_alt)], na.rm = TRUE)
}
threshold_end   <- yday(as.Date("2021-12-31")) # 365
message(sprintf("Scenario temperatures used from day of year %d (%s) to Dec 31.",
                threshold_start, format(as.Date(threshold_start - 1, origin = "2021-01-01"), "%b %d")))

# 7) COMBINE OBSERVED AND FUTURE DATA
# For each alternative, combine the historical data with a future projection.
# The future uses the alternative's forecast from the deliverable's first day
# to Dec 31, and the 14-year climatology for all other times.
env_ext_list <- map(alts, function(alt_nm) {
  obs_block <- amer_obs %>%
    filter(Date <= obs_end) %>%
    mutate(alt = alt_nm)

  dedup_pattern <- pred_by_doy %>%
    filter(alt == alt_nm) %>%
    group_by(site, doy) %>%
    summarize(temp_alt = mean(temp_alt, na.rm = TRUE), .groups = "drop")

  future_skel <- future_dates %>%
    expand_grid(site = unique(dedup_pattern$site))

  pred_block <- future_skel %>%
    left_join(dedup_pattern, by = c("doy","site")) %>%
    left_join(clim14,        by = c("doy","site")) %>%
    mutate(
      temp_raw = if_else(
        doy >= threshold_start & doy <= threshold_end & !is.na(temp_alt),
        temp_alt,
        clim_temp
      ),
      temp = if_else(site == "AveHazel", pmax(temp_raw, 7), pmax(temp_raw, 8))
    ) %>%
    select(Date, site, temp) %>%
    mutate(alt = alt_nm)

  bind_rows(obs_block, pred_block) %>%
    arrange(Date, site)
}) %>% set_names(alts)

# 8) FINAL DATA PREPARATION & SAVING
# Combine all alternatives into a single data frame
df_all <- bind_rows(env_ext_list, .id = "env")

# Save the final data objects
saveRDS(env_ext_list, file.path(app_data_dir, "env_ext_list.rds"))
saveRDS(df_all, file.path(app_data_dir, "df_all.rds"))

print(paste("Saved env_ext_list.rds with", length(env_ext_list), "runs to", app_data_dir))
print(paste("Saved df_all.rds with", nrow(df_all), "rows"))

# 9) VISUALIZATIONS
# Plot of all runs
ggplot(df_all, aes(Date, temp, color = site)) +
  geom_line(size = 0.5, alpha = 0.8) +
  facet_wrap(~ env, ncol = 6, scales = "free_y") +
  labs(
    title = sprintf("Observed + Predicted Temp by run (%d total)", length(alts)),
    x     = "Date",
    y     = "Temperature (°C)",
    color = "Site"
  ) +
  theme_minimal(base_size = 10)

# Plot focusing on the first scenario window after the decision date
win_start <- as.Date(sprintf("%d-10-01", year(obs_end) + (month(obs_end) >= 10)))
future_temp <- df_all %>%
  filter(Date >= win_start & Date <= win_start + 91) %>%
  mutate(env = factor(env, levels = alts))

ggplot(future_temp, aes(x = Date, y = temp, color = site)) +
  geom_line(size = 1) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.9) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  scale_y_continuous(name = "Temperature (°C)", breaks = seq(0, 25, 2)) +
  labs(x = NULL, y = "Temperature (°C)", color = "Site") +
  facet_wrap(~ env, ncol = 6) +
  theme_minimal(base_size = 10) +
  theme(
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text    = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold")
  )

print("Temperature data processing complete!")
print(sprintf("Created %d runs (%d scenarios x %d met years): %s",
              length(alts), n_distinct(alt_key$alt), length(HYDRO_YEARS),
              paste(unique(alt_key$alt), collapse = ", ")))
