library(tidyverse)
library(janitor)

# 1) Read the data in (you can paste in your block of text)
dat <- read_csv(
  "year,metric,value
  2018,Survival to Tower Bridge,0.681
  2018,SE to Tower Bridge,0.038
  2018,Survival from TB to Benicia,NA
  2018,SE from TB to Benicia,NA
  2018,Survival overall,0.02
  2018,SE overall,0.011
  2023,Survival to Tower Bridge,0.703
  2023,SE to Tower Bridge,0.026
  2023,Survival from TB to Benicia,0.459
  2023,SE from TB to Benicia,0.28
  2023,Survival overall,0.34
  2023,SE overall,0.021
  2024,Survival to Tower Bridge,0.198
  2024,SE to Tower Bridge,0.025
  2024,Survival from TB to Benicia,0.078
  2024,SE from TB to Benicia,0.038
  2024,Survival overall,0.016
  2024,SE overall,0.008
  2025,Survival to Tower Bridge,0.791
  2025,SE to Tower Bridge,0.027
  2025,Survival from TB to Benicia,0.098
  2025,SE from TB to Benicia,0.022
  2025,Survival overall,0.077
  2025,SE overall,0.017",
  show_col_types = FALSE
)

# 2) Compute the average for each metric
summary_outmigrants <- dat %>%
  group_by(metric) %>%
  summarize(
    mean_value = mean(value, na.rm = TRUE),
    .groups = "drop"
  )


# 2) copy–paste your table (with single spaces or tabs) into a string—
#    I’ve replaced the “Brood Year” header with “Brood” so read_table2 will parse it cleanly.
txt <- "
Brood  SAR      Standard_Error  Estimated_Returns  Number_Released  Tag_Codes  Singles
2021   0        0               0.00               1162794          9          N
2020   0.000533 0.000139        540.46             1013809          9          N
2019   0.007814 0.002417        9767.33            1249975          9          N
2018   0.006473 0.001837        8008.69            1237183          9          Y
2017   0.011337 0.004006       11323.14             998809          6          Y
2016   0.026738 0.004509       23232.18             868883          6          Y
2015   0.005602 0.002546        5833.67            1041278          6          Y
2014   0.007632 0.001713        7493.41             981822          6          Y
2013   0.008113 0.001412        7272.82             896419          4          N
2012   0.004322 0.000859        5227.49            1209525          4          N
2011   0.007487 0.002409       10528.53            1406264          5          Y
2010   0.020044 0.001845       27741.39            1384035          6          Y
2009   0.027241 0.003729       36723.21            1348086          5          Y
2008   0.006971 0.001765        8674.53            1244405          5          Y
2007   0.001464 0.000277        1784.50            1218755          3          Y
2006   0.000305 0.000019         380.69            1248450          4          N
2003   0.000144 0.000095          13.44             93591          5          N
2001   0.013149 0.001544        9383.59             713619          3          Y
2000   0.038315 0.005142       25542.24             666646          9          N
1989   0.001845 0.000505         342.22            185466          4          N
1988   0.014186 0.001891        2675.21            188580          4          N
1987   0.015561 0.003131        1968.48            126498          3          N
1986   0.021215 0.018381        2087.67             98407          2          N
1985   0.092756 0.007340        9460.36            101992          2          Y
1984   0.049160 0.038030        2427.47             49379          2          N
1983   0.028675 0.002922        4252.63            148304          3          N
1982   0.000609 0.000169          86.11            141285          3          N
"

# 3) parse it
sar_df <- read_table2(txt)

# 4) compute the mean SAR
sar_df %>%
  summarise(
    mean_SAR = mean(SAR, na.rm = TRUE),
    Q1_SAR   = quantile(SAR, 0.25, na.rm = TRUE),
    Q3_SAR   = quantile(SAR, 0.75, na.rm = TRUE),
    max_SAR = max(SAR, na.rm = TRUE),
    min_SAR = min(SAR, na.rm = TRUE)
  )


# Filter SAR data to 2011–2021
sar_filtered <- sar_df %>%
  filter(Brood >= 2011 & Brood <= 2021)

# Calculate mean, Q1, and Q3 SAR for this range
sar_stats <- sar_filtered %>%
  summarise(
    mean_SAR = mean(SAR, na.rm = TRUE),
    Q1_SAR   = quantile(SAR, 0.25, na.rm = TRUE),
    Q3_SAR   = quantile(SAR, 0.75, na.rm = TRUE)
  )

index_df <- read_csv("scripts_data/hci_1753379372_198.csv") %>%
  clean_names() %>%   # standardizes to e.g., wy, code
  mutate(brood = as.numeric(wy)) %>%
  select(brood, code) %>%
  distinct(brood, .keep_all = TRUE)  # ensure only one row per brood year

sar_labeled <- sar_df %>%
  filter(Brood >= 2011, Brood <= 2021) %>%
  left_join(index_df, by = c("Brood" = "brood"))

# Plot
ggplot(sar_labeled, aes(x = Brood, y = SAR)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 2) +
  geom_errorbar(aes(ymin = SAR - Standard_Error, ymax = SAR + Standard_Error),
                width = 0.5, color = "gray30") +
  geom_hline(yintercept = mean(sar_labeled$SAR), color = "black", linetype = "solid", linewidth = 0.8) +
  geom_hline(yintercept = quantile(sar_labeled$SAR, 0.25), color = "gray40", linetype = "dashed") +
  geom_hline(yintercept = quantile(sar_labeled$SAR, 0.75), color = "gray40", linetype = "dashed") +
  geom_text(aes(label = code), vjust = -1.2, size = 3.5, fontface = "bold") +
  scale_x_continuous(breaks = 2011:2021) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(x = "Brood Year", y = "SAR (%)") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(face = "bold")
  )



