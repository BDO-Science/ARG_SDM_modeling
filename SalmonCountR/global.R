# global.R

library(tidyverse); library(lubridate)
source("functions.R")

# load precomputed data (fast)
env_ext_list     <- readRDS("env_ext_list.rds")
df_all           <- readRDS("df_all.rds")
egg_summary      <- readRDS("egg_summary.rds")
surv_lookup_full <- readRDS("surv_lookup_full.rds")
base_P_list      <- readRDS("base_P_list.rds")
S_seed           <- readRDS("S_seed.rds")
stoch_SAR_opts   <- readRDS("stoch_SAR_opts.rds")

# any constants UI needs
n_sim      <- 50  # or hard‐code 50
env_levels <- as.character(sort(as.numeric(unique(egg_summary$env))))





