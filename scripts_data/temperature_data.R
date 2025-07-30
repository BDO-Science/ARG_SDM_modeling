library(dataRetrieval) # for readNWISdata()
library(tidyverse)
library(readxl)

# 1) PARAMETERS
obs_start   <- as.Date("2011-09-01")
obs_end     <- as.Date("2024-05-31")
last_wy     <- 2150                    # final water‐year to simulate
sim_end     <- as.Date(sprintf("%04d-08-31", last_wy + 1))

# 2) GET OBSERVED GAUGE TEMPS (Daily avg)
stations   <- c("11446980","11446500")
param_cd   <- c("00010")

amer_obs <- map_df(stations, function(stn) {
  readNWISdata(
    sites       = stn,
    parameterCd = param_cd,
    service     = "uv",
    startDate   = obs_start,
    endDate     = obs_end
  ) %>%
    transmute(
      Date  = as.Date(dateTime),
      site  = recode(site_no,
                     "11446980" = "AveWatt",
                     "11446500" = "AveHazel"),
      temp  = X_00010_00000
    )
}) %>%
  group_by(Date, site) %>%
  summarise(temp = mean(temp, na.rm=TRUE), .groups="drop")

amer_obs <- amer_obs %>%
  group_by(Date, site) %>%
  summarise(temp = mean(temp, na.rm=TRUE), .groups="drop") %>%
  mutate(
    # Hazel floored at 7°C, all others at 8°C
    temp = if_else(site == "AveHazel",
                   pmax(if_else(is.na(temp)|is.nan(temp), 7, temp), 7),
                   pmax(if_else(is.na(temp)|is.nan(temp), 8, temp), 8))
  )

# ───────────────────────────────────────────────────────────────────────────────
# 1) Build a 14‑year DOY climatology of amer_obs
# ───────────────────────────────────────────────────────────────────────────────
clim14 <- amer_obs %>%
  mutate(doy = yday(Date)) %>%
  group_by(site, doy) %>%
  summarize(clim_temp = mean(temp, na.rm = TRUE), .groups = "drop")

# 3) READ ALT FORECASTS & BUILD A DOY→temp TABLE FOR EACH ALT
xlsx_path <- "ARG_LAR_TempModeling_10-21-24.xlsx"
alts      <- excel_sheets(xlsx_path)[-1]  # assume first sheet is metadata

pred_by_doy <- map_df(alts, function(alt) {
  read_excel(xlsx_path, sheet=alt) %>%
    mutate(Date = as.Date(Date)) %>%
    pivot_longer(
      cols      = starts_with("Ave"),
      names_to  = "site",
      values_to = "temp_alt"
    ) %>%
    mutate(
      doy = yday(Date),
      alt = alt
    ) %>%
    select(alt, site, doy, temp_alt)
})

# 4) MAKE FUTURE DATE SKELETON
future_dates <- tibble(
  Date = seq(obs_end + 1, sim_end, by="day")
) %>% mutate(doy = yday(Date))

# DOY thresholds for Oct 18 and Dec 31 (use a non‑leap year for consistency)
threshold_start <- yday(as.Date("2021-10-18"))  # 291
threshold_end   <- yday(as.Date("2021-12-31"))  # 365

# ───────────────────────────────────────────────────────────────────────────────
# 2) Rebuild env_ext_list with climatology‑fill
# ───────────────────────────────────────────────────────────────────────────────
env_ext_list <- map(alts, function(alt_nm) {
  # a) Observed block unchanged
  obs_block <- amer_obs %>%
    filter(Date <= obs_end) %>%
    mutate(alt = alt_nm)
  
  # b) Prepare the raw DOY→temp_alt pattern
  dedup_pattern <- pred_by_doy %>%
    filter(alt == alt_nm) %>%
    group_by(site, doy) %>%
    summarize(temp_alt = mean(temp_alt, na.rm = TRUE), .groups = "drop")
  
  # c) Build the full future skeleton
  future_skel <- future_dates %>%
    expand_grid(site = unique(dedup_pattern$site))
  
  # d) Join forecast + climatology, then apply the fill rule
  pred_block <- future_skel %>%
    left_join(dedup_pattern, by = c("doy","site")) %>%
    left_join(clim14,        by = c("doy","site")) %>%
    mutate(
      # fill logic as before
      temp_raw = if_else(
        doy >= threshold_start & doy <= threshold_end & !is.na(temp_alt),
        temp_alt,
        clim_temp
      ),
      # now floor by site
      temp = if_else(site == "AveHazel",
                     pmax(temp_raw, 7),
                     pmax(temp_raw, 8))
    ) %>%
    select(Date, site, temp) %>%
    mutate(alt = alt_nm)
  
  bind_rows(obs_block, pred_block) %>%
    arrange(Date, site)
}) %>% set_names(alts)

# 6) Now rebuild your site/Date lookup lists *after* the full horizons exist
site_temps_list <- map(env_ext_list, ~ split(.x$temp, .x$site))
site_dates_list <- map(env_ext_list, ~ split(.x$Date, .x$site))
date_idx_list   <- map(site_dates_list, function(dlist) {
  map(dlist, ~ setNames(seq_along(.x), as.character(.x)))
})

df_all <- bind_rows(env_ext_list, .id = "env") 

ggplot(df_all, aes(Date, temp, color = site)) +
  geom_line(size = 0.5, alpha = 0.8) +
  facet_wrap(~ env, ncol = 2, scales = "free_y") +
  labs(
    title = "Observed + Predicted Temp by Alternative",
    x     = "Date",
    y     = "Temperature (°C)",
    color = "Site"
  ) +
  theme_minimal(base_size = 12)

test_date <- obs_end + 2500  # 2024‑06‑01

map_df(alts, function(alt_nm) {
  env_ext_list[[alt_nm]] %>%
    filter(Date == test_date) %>%
    mutate(alt = alt_nm)
}) %>%
  select(alt, site, temp) %>%
  arrange(alt, site)

saveRDS(env_ext_list, "env_ext_list.rds")
saveRDS(df_all, "df_all.rds")

# 1. Filter to Oct–Feb and observed period
obs_temp <- amer_obs %>%
  filter(Date >= as.Date("2011-10-01") & Date <= as.Date("2024-02-29")) %>%
  filter(month(Date) %in% c(10, 11, 12, 1, 2)) %>%
  mutate(month_day = format(Date, "%m-%d"))

# 2. Create dummy reference date from 2023-10-01 to 2024-02-29
ref_dates <- tibble(
  dummy_date = seq(as.Date("2023-10-01"), as.Date("2024-02-29"), by = "1 day")
) %>%
  mutate(month_day = format(dummy_date, "%m-%d"))

obs_temp <- obs_temp %>%
  left_join(ref_dates, by = "month_day")

# 3. Summarize percentiles
obs_summary <- obs_temp %>%
  group_by(site, dummy_date) %>%
  summarise(
    p10    = quantile(temp, 0.10, na.rm = TRUE),
    median = median(temp, na.rm = TRUE),
    p90    = quantile(temp, 0.90, na.rm = TRUE),
    .groups = "drop"
  )

# 4. Plot
ggplot(obs_summary, aes(x = dummy_date)) +
  geom_ribbon(aes(ymin = p10, ymax = p90, fill = site), alpha = 0.4) +
  geom_line(aes(y = median, color = site), size = 0.9) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.9) +
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.9) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  scale_y_continuous(name = "Temperature (°C)", breaks = seq(0, 25, 2)) +
  labs(x = NULL, y = "Temperature (°C)", color = "Site", fill = "Site") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_blank(),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text    = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold")
  )

# Filter df_all for future prediction window: Oct 18 – Dec 31, 2024
future_temp <- df_all %>%
  filter(site != "AveFol") %>%
  filter(Date >= as.Date("2024-10-18") & Date <= as.Date("2024-12-31"))

future_temp <- future_temp %>%
  mutate(env = factor(env, levels = as.character(1:10)))  # ensures 10 is last

# Plot with lines per site, faceted by alternative
ggplot(future_temp, aes(x = Date, y = temp, color = site)) +
  geom_line(size = 1) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.9) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  scale_y_continuous(name = "Temperature (°C)", breaks = seq(0, 25, 2)) +
  labs(x = NULL, y = "Temperature (°C)", color = "Site") +
  facet_wrap(~ env, ncol = 2) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_blank(),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text    = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold")
  )

# Assume clim14 is the climatology table by doy and site
future_anomaly <- future_temp %>%
  mutate(doy = yday(Date)) %>%
  left_join(clim14, by = c("site", "doy")) %>%
  mutate(temp_anomaly = temp - clim_temp)

ggplot(future_anomaly, aes(x = Date, y = temp_anomaly, color = site)) +
  geom_line() +
  facet_wrap(~ env, ncol = 2) +
  labs(x = "Date", y = "Temp Anomaly (°C)", color = "Site") +
  scale_color_viridis_d() +
  theme_minimal(base_size = 14)

ddays <- future_temp %>%
  mutate(thresh_temp = pmax(temp - 0, 0)) %>%
  group_by(env, site) %>%
  summarise(degree_days = sum(thresh_temp), .groups = "drop")

ggplot(ddays, aes(x = site, y = degree_days, fill = env)) +
  geom_col(position = position_dodge()) +
  labs(y = "Cumulative Degree Days > 0°C", x = "Site", fill = "Alt") +
  scale_fill_viridis_d() +
  theme_minimal(base_size = 14)

ggplot(future_temp, aes(x = env, y = temp, fill = site)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "Alternative", y = "Daily Temperature (°C)", fill = "Site") +
  scale_fill_viridis_d() +
  theme_minimal(base_size = 14)

threshold_summary <- future_temp %>%
  mutate(thresh_exceed = temp > 12.8) %>%
  group_by(env, site) %>%
  summarise(n_exceed_days = sum(thresh_exceed), .groups = "drop")

ggplot(threshold_summary, aes(x = site, y = n_exceed_days, fill = env)) +
  geom_col(position = "dodge") +
  labs(y = "# Days > 12.8°C", x = "Site", fill = "Alt") +
  scale_fill_viridis_d() +
  theme_minimal(base_size = 14)

