# 1. load libraries
library(tidyverse)
library(here)

# 2. point to the folder containing your CSVs
data_dir <- here("data_raw", "american_river_data")

# 1. list any CSV with both “inland” and “taf” in the name (case-insensitive)
csvs <- list.files(
  data_dir,
  pattern    = "inland.*taf.*\\.csv$",
  full.names = TRUE,
  ignore.case = TRUE
)

# 2. read & row-bind
all_data <- csvs %>% 
  set_names() %>% 
  map_dfr(read_csv, .id = "source")

# 3. filter to AMN
amn_only <- all_data %>% 
  filter(location == "AMN")

# 4. preview
print(amn_only)