library(dataRetrieval)
library(tidyverse)
library(lubridate)

# PARAMETERS
obs_start <- as.Date("2011-09-01")
obs_end   <- as.Date("2024-05-31")
stations  <- "11446500"    # only AveHazel
param_cd  <- "00060"

# 1) Pull daily mean discharge and rename columns
flow_dv <- readNWISdv(
  siteNumbers = stations,
  parameterCd = param_cd,
  startDate   = obs_start,
  endDate     = obs_end
) %>%
  renameNWISColumns()  # X_00060_00003 → Flow

# 2) Compute water‐year Oct–Jan monthly means
amer_flow_monthly <- flow_dv %>%
  mutate(
    site       = "AveHazel",
    month      = month(Date),
    water_year = if_else(month >= 10, year(Date) + 1, year(Date))
  ) %>%
  filter(month %in% c(10, 11, 12, 1)) %>%
  group_by(site, month) %>%
  summarise(
    mean_flow_cfs = mean(Flow, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(month = month.name[month])

# View the result
print(amer_flow_monthly)



# 2) Filter to Oct 1 – Feb 29
flow_season <- flow_dv %>%
  mutate(
    year        = year(Date),
    month       = month(Date),
    water_year  = if_else(month >= 10, year + 1, year),
    doy         = yday(Date)
  ) %>%
  filter(
    (month %in% c(10, 11, 12)) | (month %in% c(1, 2))
  ) %>%
  filter(!(month == 2 & day(Date) > 29)) %>%
  filter(!is.na(Flow))

# 3) Assign a "dummy date" to align Oct 1 – Feb 29 as continuous timeline
# We'll use 2023-10-01 to 2024-02-29 as a fake year sequence
reference_dates <- tibble(
  date_dummy = seq(as.Date("2023-10-01"), as.Date("2024-02-29"), by = "1 day")
) %>%
  mutate(day_of_season = row_number())

# Create daily summary across years
flow_summary <- flow_season %>%
  group_by(date_dummy) %>%
  summarise(
    p10   = quantile(Flow, 0.10, na.rm = TRUE),
    median = median(Flow, na.rm = TRUE),
    p90   = quantile(Flow, 0.90, na.rm = TRUE),
    .groups = "drop"
  )

# Plot percentiles with ribbon
ggplot(flow_summary, aes(x = date_dummy)) +
  geom_ribbon(aes(ymin = p10, ymax = p90), fill = "lightgrey", alpha = 0.5) +
  geom_line(aes(y = median), color = "black", size = 1) +
  scale_y_continuous(
    labels = scales::comma,
    breaks = seq(0, 25000, by = 1500)  # add ticks every 5,000 cfs
  ) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  labs(
    x = NULL,
    y = "Flow (cfs)",
    title = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text    = element_text(face = "bold", size = 12)
  )
