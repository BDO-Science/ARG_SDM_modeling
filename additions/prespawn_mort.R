library(tidyverse)
library(lubridate)
library(janitor)
library(ggridges)
library(Hmisc)  # for wtd.quantile()

#exploring yuba river fall-run adult upstream passage to get an idea of when they get into the river
# Load your Excel file (update path if needed)
df <- read.csv("yuba_daily_corrected_passage.csv") %>%
  clean_names() %>%
  mutate(
    date      = as_date(date),
    # July–Dec belong to that calendar year, Jan–Jun to previous calendar year
    bio_year  = if_else(month(date) >= 7, year(date), year(date) - 1),
    # Map July–Dec → 2000, Jan–Jun → 2001 (for plotting on a common axis)
    season_date = case_when(
      month(date) >= 7  ~ as_date(sprintf("2000-%02d-%02d", month(date), day(date))),
      TRUE              ~ as_date(sprintf("2001-%02d-%02d", month(date), day(date)))
    )
  )

# 2. Aggregate & filter
fall_counts <- df %>%
  filter(run == "fall") %>%
  group_by(bio_year, season_date) %>%
  summarise(daily_count = sum(count, na.rm = TRUE), .groups="drop") %>%
  filter(
    !is.na(season_date),
    is.finite(daily_count),
    season_date >= as_date("2000-07-01"),
    season_date <= as_date("2001-06-30")
  )

xlims <- as_date(c("2000-07-01", "2001-06-30"))

# A) Faceted weighted histogram by bio_year
p_hist <- ggplot(fall_counts, aes(x = season_date, weight = daily_count)) +
  geom_histogram(binwidth = 14, fill = "steelblue", alpha = 0.7, color = "white") +
  facet_wrap(~ bio_year, ncol = 4, scales = "free_y") +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b",
    limits      = xlims
  ) +
  labs(
    title = "Seasonal Histogram of Fall‑Run Counts by Bio‑Year (Jul 1–Jun 30)",
    x     = "Season Date (Jul 1 – Jun 30)",
    y     = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# B) Faceted weighted density by bio_year
p_den <- ggplot(fall_counts, aes(x = season_date, weight = daily_count)) +
  geom_density(fill = "coral", alpha = 0.6) +
  facet_wrap(~ bio_year, ncol = 4, scales = "free_y") +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b",
    limits      = xlims
  ) +
  labs(
    title = "Seasonal Density of Fall‑Run Counts by Bio‑Year (Jul 1–Jun 30)",
    x     = "Season Date (Jul 1 – Jun 30)",
    y     = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# C) Overall seasonal histogram & density
p_all_hist <- ggplot(fall_counts, aes(x = season_date, weight = daily_count)) +
  geom_histogram(binwidth = 14, fill = "darkgreen", alpha = 0.7, color = "white") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b", limits = xlims) +
  labs(
    title = "Overall Seasonal Histogram of Fall‑Run Counts (Jul 1–Jun 30)",
    x     = "Season Date",
    y     = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Calculate the weighted median
median_date <- as.Date(wtd.quantile(fall_counts$season_date, weights = fall_counts$daily_count, probs = 0.5))

# Plot with vertical line at the weighted median
p_all_den <- ggplot(fall_counts, aes(x = season_date, weight = daily_count)) +
  geom_density(fill = "lightblue", alpha = 0.6) +
  geom_vline(xintercept = as.numeric(median_date), color = "black", linetype = "dashed", linewidth = 1) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d", limits = xlims) +
  labs(
    title = NULL,
    x     = NULL,
    y     = "Density"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y  = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

# 3. Print
print(p_hist)
print(p_den)
print(p_all_hist)
print(p_all_den)

# 1) A simple logistic “pre‐spawn” survival function
#    deg_day   = cumulative °C‐days adults experience before spawning
#    intercept = logit‐intercept (you’ll fill in)
#    beta      = slope (per °C‐day; you’ll fill in)
# 1) Logistic pre‑spawn survival
surv_adult_prespawn <- function(deg_day,
                                intercept = 3.0,
                                beta      = -0.00067) {
  # plogis(x) == exp(x)/(1+exp(x))
  plogis(intercept + beta * deg_day)
}

# generate a sequence of cumulative degree-days (°C·days)
deg_days <- seq(0, 10000, length.out = 1000)
surv <- surv_adult_prespawn(deg_days)

ggplot(data.frame(deg_days, surv), aes(deg_days, surv)) +
  geom_line() +
  scale_y_continuous(
    limits      = c(0, 1),
    breaks      = seq(0, 1, by = 0.05),
    minor_breaks= seq(0, 1, by = 0.05)
  ) +
  labs(
    x     = "Cumulative Degree-Days (°C·days)",
    y     = "Survival Probability",
    title = "Logistic Pre-Spawn Survival Curve\n(intercept = 3.0, beta = -0.00067)"
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title       = element_text(face = "bold", size = 14),
    axis.text        = element_text(face = "bold", size = 12),
    axis.ticks.length= unit(0.25, "cm")
  )
