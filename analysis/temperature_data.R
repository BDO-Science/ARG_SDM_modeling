# Builds the daily temperature series that precompute.R and the app read, from a
# CE-QUAL-W2 temperature deliverable plus observed USGS gauge data.
#
#   deliverable .xlsx  ->  env_ext_list.rds, df_all.rds  (in the scenario folder)
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

# Checks the deliverable against the layout the rest of this script assumes.
# Columns are read by POSITION, so a reordered or relabelled file would otherwise
# run to completion and give wrong answers without any error.
check_temperature_file <- function(input_path, hydro_years, scenarios) {
  if (!file.exists(input_path)) stop("Temperature file not found: ", input_path)
  sheets  <- excel_sheets(input_path)
  missing <- setdiff(hydro_years, sheets)
  if (length(missing)) {
    stop("Sheet(s) ", paste(missing, collapse = ", "), " not found in ", basename(input_path),
         ". Expected one sheet per meteorological year, named ",
         paste(hydro_years, collapse = ", "), ". Sheets found: ",
         paste(sheets, collapse = ", "))
  }
  for (hy in hydro_years) {
    hdr <- read_excel(input_path, sheet = hy, col_names = FALSE, n_max = 2,
                      .name_repair = "minimal")
    if (ncol(hdr) < 20) {
      stop("Sheet ", hy, " has ", ncol(hdr), " columns; expected at least 20 ",
           "(Date, JDAY, then an AveWatt/AveHazel pair for each of 9 scenarios).")
    }
    got_scen <- trimws(unlist(hdr[1, seq(3, 19, 2)]))
    if (!identical(unname(got_scen), scenarios)) {
      stop("Sheet ", hy, ", row 1: scenario labels in columns C, E, G, ... S must be, in order: ",
           paste(scenarios, collapse = ", "), ". Found: ",
           paste(got_scen, collapse = ", "))
    }
    got_site <- trimws(unlist(hdr[2, 3:20]))
    want_site <- rep(c("AveWatt", "AveHazel"), 9)
    if (!identical(unname(got_site), want_site)) {
      stop("Sheet ", hy, ", row 2: columns C to T must alternate AveWatt, AveHazel. Found: ",
           paste(got_site, collapse = ", "))
    }
    dat   <- suppressMessages(read_excel(input_path, sheet = hy, skip = 1))
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
    temps <- unlist(dat[, 3:20])
    if (!is.numeric(temps)) stop("Sheet ", hy, ": temperature columns are not numeric.")
    rng <- range(temps, na.rm = TRUE)
    if (rng[1] < 0 || rng[2] > 30) {
      stop(sprintf(paste0("Sheet %s: temperatures run %.1f to %.1f. Expected degrees ",
                          "Celsius (roughly 4-30); values above 30 usually mean Fahrenheit."),
                   hy, rng[1], rng[2]))
    }
    message(sprintf("Sheet %s: %s to %s, %d days.",
                    hy, format(min(dates, na.rm = TRUE)), format(max(dates, na.rm = TRUE)),
                    sum(!is.na(dates))))
    # Spawning, incubation and the Oct/Nov means that drive spawn timing all
    # fall in October-December; a file that misses them has nothing to act on.
    if (!any(month(dates) %in% 10:12, na.rm = TRUE)) {
      stop("Sheet ", hy, " has no dates in October-December, when the model ",
           "uses the scenario temperatures.")
    }
  }
  invisible(TRUE)
}

# This function reads the temperature modeling results and reformats into alternatives
prepare_temperature_file <- function() {
  input_path <- temp_file

  sheets <- excel_sheets(input_path)
  hydro_years <- c("2011", "2014", "2017", "2020")
  check_temperature_file(
    input_path, hydro_years,
    c("No Bypass", "Scenario 1", "Scenario 2", "Scenario 2b", "Scenario 2c",
      "Scenario 3", "Scenario 4", "Scenario 5", "Scenario 6")
  )

  all_sheets <- list()
  
  # Copy metadata sheets
  if ("Scenario Summary" %in% sheets) {
    all_sheets[["Scenario Summary"]] <- read_excel(input_path, sheet = "Scenario Summary")
  }
  if ("Flow" %in% sheets) {
    all_sheets[["Flow"]] <- read_excel(input_path, sheet = "Flow")
  }
  
  # Read each hydro year sheet
  for (hydro_year in hydro_years) {
    if (hydro_year %in% sheets) {
      # Skip the first row to use the second row as headers
      df <- read_excel(input_path, sheet = hydro_year, skip = 1)
      all_sheets[[hydro_year]] <- df
    }
  }
  
  # Create alternatives format
  alternatives_list <- list()
  
  # Updated scenarios list to match the actual data
  scenarios <- c("No Bypass", "Scenario 1", "Scenario 2", "Scenario 2b", 
                 "Scenario 2c", "Scenario 3", "Scenario 4", "Scenario 5", "Scenario 6")
  
  # Metadata
  metadata_rows <- list()
  alt_num <- 1
  for (hydro_year in hydro_years) {
    for (scenario_idx in 1:length(scenarios)) {
      metadata_rows[[alt_num]] <- data.frame(
        Alternative = alt_num,
        Scenario = scenarios[scenario_idx],
        Hydro_Year = hydro_year
      )
      alt_num <- alt_num + 1
    }
  }
  alternatives_list[["metadata"]] <- do.call(rbind, metadata_rows)
  
  # Extract data for each alternative
  alt_counter <- 1
  for (hydro_year in hydro_years) {
    if (hydro_year %in% names(all_sheets)) {
      df_hydro <- all_sheets[[hydro_year]]
      
      # Column pairs for each scenario (columns 3-4, 5-6, 7-8, etc.)
      # No Bypass: cols 3-4
      # Scenario 1: cols 5-6
      # Scenario 2: cols 7-8
      # Scenario 2b: cols 9-10
      # Scenario 2c: cols 11-12
      # Scenario 3: cols 13-14
      # Scenario 4: cols 15-16
      # Scenario 5: cols 17-18
      # Scenario 6: cols 19-20
      scenario_cols <- list(c(3,4), c(5,6), c(7,8), c(9,10), c(11,12), 
                            c(13,14), c(15,16), c(17,18), c(19,20))
      
      for (cols in scenario_cols) {
        alt_data <- data.frame(
          Date = as.Date(df_hydro[[1]]),
          AveWatt = df_hydro[[cols[1]]],
          AveHazel = df_hydro[[cols[2]]]
        )
        alternatives_list[[as.character(alt_counter)]] <- alt_data
        alt_counter <- alt_counter + 1
      }
    }
  }
  
  # The published run keeps its intermediate workbook in data_raw/; a scenario
  # run keeps it in its own folder so the published one is not overwritten.
  output_path <- if (is_published_dir) {
    here("data_raw", "ARG_LAR_TempModeling_alternatives.xlsx")
  } else {
    file.path(app_data_dir, "temperature_alternatives.xlsx")
  }
  write_xlsx(alternatives_list, output_path)
  print(paste("Created", alt_counter - 1, "alternatives in", output_path))
  return(output_path)
}

# --- SCRIPT EXECUTION STARTS HERE ---

# 1) PREPARE DATA: Run the function to reformat the alternatives file
xlsx_path <- prepare_temperature_file()

# 2) SET PARAMETERS
# obs_end is the DECISION DATE: observed gauge temperatures are used through
# it, the scenario temperatures after it. Before the decision every alternative
# is the same river, so observed data is the right input up to that day.
# Default: the day before the deliverable starts (2025-09-21 for the 2025 file,
# which is what the published run used).
deliverable_start <- min(do.call(c, lapply(c("2011", "2014", "2017", "2020"), function(s)
  as.Date(suppressMessages(read_excel(temp_file, sheet = s, skip = 1))[[1]]))), na.rm = TRUE)
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

# 5) READ ALTERNATIVES: Load the generated alternatives
alts <- excel_sheets(xlsx_path)
alts <- alts[alts != "metadata"] # Remove the metadata sheet

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

print(paste("Saved env_ext_list.rds with", length(env_ext_list), "alternatives to", app_data_dir))
print(paste("Saved df_all.rds with", nrow(df_all), "rows"))

# 9) VISUALIZATIONS
# Plot of all alternatives
ggplot(df_all, aes(Date, temp, color = site)) +
  geom_line(size = 0.5, alpha = 0.8) +
  facet_wrap(~ env, ncol = 6, scales = "free_y") +
  labs(
    title = "Observed + Predicted Temp by Alternative (36 total)",
    x     = "Date",
    y     = "Temperature (°C)",
    color = "Site"
  ) +
  theme_minimal(base_size = 10)

# Plot focusing on the 2024 forecast window
future_temp <- df_all %>%
  filter(site != "AveFol") %>%
  filter(Date >= as.Date("2024-10-18") & Date <= as.Date("2024-12-31")) %>%
  mutate(env = factor(env, levels = as.character(1:36)))

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
print(paste("Created", length(alts), "alternatives (9 scenarios × 4 hydro years)"))