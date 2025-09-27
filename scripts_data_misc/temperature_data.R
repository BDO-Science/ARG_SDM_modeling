# Load necessary libraries
library(dataRetrieval) # for readNWISdata()
library(tidyverse)
library(readxl)
library(writexl)

# This function reads the raw temperature modeling results,
# fills in empty placeholder scenarios with simulated data,
# and then reformats the data into 28 distinct alternatives (7 scenarios x 4 hydro years).
prepare_temperature_file <- function() {
  input_path <- "scripts_data_misc/SDM Power Bypass Temperature Modeling Results.xlsx"
  
  sheets <- excel_sheets(input_path)
  hydro_years <- c("2011", "2014", "2017", "2020")
  
  all_sheets_filled <- list()
  
  # Copy metadata sheets
  if ("Scenario Summary" %in% sheets) {
    all_sheets_filled[["Scenario Summary"]] <- read_excel(input_path, sheet = "Scenario Summary")
  }
  if ("Flow" %in% sheets) {
    all_sheets_filled[["Flow"]] <- read_excel(input_path, sheet = "Flow")
  }
  
  # Fill empty scenarios in each hydro year
  for (hydro_year in hydro_years) {
    if (hydro_year %in% sheets) {
      # Skip the first row to use the second row as headers
      df <- read_excel(input_path, sheet = hydro_year, skip = 1)
      
      # Explicitly convert columns to numeric to avoid errors
      base_watt <- as.numeric(df[[3]])
      base_hazel <- as.numeric(df[[4]])
      jday <- as.numeric(df[[2]])
      
      set.seed(123 + as.numeric(hydro_year))
      
      # ====================================================================
      # FIX APPLIED HERE:
      # Removed the `if (ncol(df) >= ...)` checks. This forces R to
      # create columns 5 through 16, which is the desired behavior.
      # ====================================================================
      
      # Create and fill columns for Scenario 1
      adj <- runif(nrow(df), -0.5, -0.2)
      df[[5]] <- round(base_watt + adj, 2)
      df[[6]] <- round(base_hazel + adj, 2)
      
      # Create and fill columns for Scenario 2
      adj <- runif(nrow(df), -1.0, -0.5)
      df[[7]] <- round(base_watt + adj, 2)
      df[[8]] <- round(base_hazel + adj, 2)
      
      # Create and fill columns for Scenario 3
      adj <- runif(nrow(df), -1.5, -0.8)
      df[[9]] <- round(base_watt + adj, 2)
      df[[10]] <- round(base_hazel + adj, 2)
      
      # Create and fill columns for Scenario 4
      time_factor <- (jday - min(jday, na.rm = TRUE)) / (max(jday, na.rm = TRUE) - min(jday, na.rm = TRUE))
      adj <- -0.3 - time_factor * 1.2
      df[[11]] <- round(base_watt + adj, 2)
      df[[12]] <- round(base_hazel + adj, 2)
      
      # Create and fill columns for Scenario 5
      adj <- runif(nrow(df), -2.0, -1.0)
      df[[13]] <- round(base_watt + adj, 2)
      df[[14]] <- round(base_hazel + adj, 2)
      
      # Create and fill columns for Scenario 6
      adj <- runif(nrow(df), -0.3, 0.3)
      df[[15]] <- round(base_watt + adj, 2)
      df[[16]] <- round(base_hazel + adj, 2)
      
      # Apply temperature floors
      for (col in seq(3, ncol(df), 2)) {
        if (col <= ncol(df)) df[[col]] <- pmax(df[[col]], 8, na.rm = TRUE)
      }
      for (col in seq(4, ncol(df), 2)) {
        if (col <= ncol(df)) df[[col]] <- pmax(df[[col]], 7, na.rm = TRUE)
      }
      
      all_sheets_filled[[hydro_year]] <- df
    }
  }
  
  # Create 28 alternatives format
  alternatives_list <- list()
  scenarios <- c("No Bypass", "Scenario 1", "Scenario 2", "Scenario 3", 
                 "Scenario 4", "Scenario 5", "Scenario 6")
  
  # Metadata
  metadata_rows <- list()
  alt_num <- 1
  for (hydro_year in hydro_years) {
    for (scenario_idx in 1:7) {
      metadata_rows[[alt_num]] <- data.frame(
        Alternative = alt_num,
        Scenario = scenarios[scenario_idx],
        Hydro_Year = hydro_year
      )
      alt_num <- alt_num + 1
    }
  }
  alternatives_list[["metadata"]] <- do.call(rbind, metadata_rows)
  
  # Extract data for each alternative (this part will now work)
  alt_counter <- 1
  for (hydro_year in hydro_years) {
    if (hydro_year %in% names(all_sheets_filled)) {
      df_hydro <- all_sheets_filled[[hydro_year]]
      
      scenario_cols <- list(c(3,4), c(5,6), c(7,8), c(9,10), c(11,12), c(13,14), c(15,16))
      
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
  
  output_path <- "scripts_data_misc/ARG_LAR_TempModeling_placeholders.xlsx"
  write_xlsx(alternatives_list, output_path)
  print(paste("Created", alt_counter - 1, "alternatives in", output_path))
  return(output_path)
}

# --- SCRIPT EXECUTION STARTS HERE ---

# 1) PREPARE DATA: Run the function to generate the alternatives file
xlsx_path <- prepare_temperature_file()

# 2) SET PARAMETERS
obs_start <- as.Date("2011-09-01")
obs_end   <- as.Date("2025-09-21")
last_wy   <- 2150 # Final water-year to simulate
sim_end   <- as.Date(sprintf("%04d-08-31", last_wy + 1))

# 3) GET OBSERVED GAUGE DATA from NWIS
stations <- c("11446980", "11446500") # Watt Avenue and Hazel Avenue
param_cd <- "00010" # Temperature parameter code

amer_obs <- map_df(stations, function(stn) {
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

# 4) BUILD CLIMATOLOGY: Create a 14-year daily average temperature
clim14 <- amer_obs %>%
  mutate(doy = yday(Date)) %>%
  group_by(site, doy) %>%
  summarize(clim_temp = mean(temp, na.rm = TRUE), .groups = "drop")

# 5) READ ALTERNATIVES: Load the 28 generated alternatives
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

# Define the period where forecast data is used (Oct 18 - Dec 31)
threshold_start <- yday(as.Date("2021-10-18")) # 291
threshold_end   <- yday(as.Date("2021-12-31")) # 365

# 7) COMBINE OBSERVED AND FUTURE DATA
# For each alternative, combine the historical data with a future projection.
# The future uses the alternative's forecast for the Oct-Dec window,
# and the 14-year climatology for all other times.
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
saveRDS(env_ext_list, "SalmonCountR/app_data/env_ext_list.rds")
saveRDS(df_all, "SalmonCountR/app_data/df_all.rds")

print(paste("Saved env_ext_list.rds with", length(env_ext_list), "alternatives"))
print(paste("Saved df_all.rds with", nrow(df_all), "rows"))

# 9) VISUALIZATIONS
# Plot of all 28 alternatives
ggplot(df_all, aes(Date, temp, color = site)) +
  geom_line(size = 0.5, alpha = 0.8) +
  facet_wrap(~ env, ncol = 4, scales = "free_y") +
  labs(
    title = "Observed + Predicted Temp by Alternative (28 total)",
    x     = "Date",
    y     = "Temperature (°C)",
    color = "Site"
  ) +
  theme_minimal(base_size = 10)

# Plot focusing on the 2024 forecast window
future_temp <- df_all %>%
  filter(site != "AveFol") %>%
  filter(Date >= as.Date("2024-10-18") & Date <= as.Date("2024-12-31")) %>%
  mutate(env = factor(env, levels = as.character(1:28)))

ggplot(future_temp, aes(x = Date, y = temp, color = site)) +
  geom_line(size = 1) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.9) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  scale_y_continuous(name = "Temperature (°C)", breaks = seq(0, 25, 2)) +
  labs(x = NULL, y = "Temperature (°C)", color = "Site") +
  facet_wrap(~ env, ncol = 4) +
  theme_minimal(base_size = 10) +
  theme(
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text    = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold")
  )

print("Temperature data processing complete!")
print(paste("Created", length(alts), "alternatives (7 scenarios × 4 hydro years)"))