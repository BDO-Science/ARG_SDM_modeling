# scripts/make_american_river_instream.R
# Build american_river_instream.rds from a 2-col table (Flow cfs, Habitat acres)

library(tidyverse)
library(here)
# If your source is an Excel file, enable one of these:
# install.packages("readxl")
# library(readxl)

# --- A) READ THE SOURCE TABLE -------------------------------------------------
# Option 1: from Excel (uncomment and set your file + sheet + range as needed)
raw_tbl <- readxl::read_excel(
   path  = here("scripts_data_misc", "american_river_data", "LAR_2023_Chinook_Spawning_Habitat.xlsx"),
   sheet = 1,     # or the named sheet
   range = "A2:B21"  # adjust to the actual block with numbers
 )

# --- B) CLEAN + CONVERT UNITS -------------------------------------------------
# 1 acre = 4046.8564224 m^2
ACRE_TO_M2 <- 4046.8564224

instream <- raw_tbl %>%
  janitor::clean_names() %>%
  rename(
    flow_cfs = flow_cfs,                # will exist after clean_names()
    habitat_acres = suitable_chinook_salmon_spawning_habitat_acres       # will exist after clean_names()
  ) %>%
  mutate(
    FR_spawn_wua = habitat_acres * ACRE_TO_M2,  # WUA in m^2
    watershed    = "American River"
  ) %>%
  dplyr::select(flow_cfs, FR_spawn_wua) %>%
  arrange(flow_cfs)

# Quick sanity checks
stopifnot(all(is.finite(instream$flow_cfs)))
stopifnot(all(instream$FR_spawn_wua > 0))
stopifnot(identical(names(instream),
                    c("flow_cfs","FR_spawn_wua")))

# --- C) SAVE AS RDS -----------------------------------------------------------
saveRDS(instream, here("SalmonCountR","app_data","american_river_instream.rds"))

# --- D) OPTIONAL: preview K_spawners calc works ------------------------------
instream %>%
  mutate(K_spawners = FR_spawn_wua / 9.29) 

