rm(list = ls())

# ================================
# RMIS -> Escapement age summaries with release/recovery filters
# ================================
suppressPackageStartupMessages({
  library(rvest)
  library(xml2)
  library(dplyr)
  library(stringr)
  library(janitor)
  library(readr)
  library(here)
  library(ggplot2)
})

# -------- helpers --------
to_int <- function(x) as.integer(suppressWarnings(as.numeric(gsub("[^0-9.-]", "", as.character(x)))))
to_num <- function(x)  as.numeric(suppressWarnings(as.numeric(gsub("[^0-9.-]", "", as.character(x)))))

# -------- parser (includes release_site) --------
rmis_read_html <- function(path) {
  stopifnot(file.exists(path))
  doc <- read_html(path)
  
  rec_tables <- xml_find_all(
    doc,
    "//table[@border='1' and .//tr[1]/td[normalize-space()='Return Year']]"
  )
  if (!length(rec_tables)) return(tibble())
  
  bind_rows(lapply(seq_along(rec_tables), function(i) {
    rec_tbl <- rec_tables[[i]]
    
    meta_tbl <- xml_find_first(
      rec_tbl,
      "preceding::table[.//td[normalize-space()='Brood Year:']][1]"
    )
    if (inherits(meta_tbl, "xml_missing")) return(tibble())
    
    brood_year   <- xml_text(xml_find_first(meta_tbl, ".//td[normalize-space()='Brood Year:']/following-sibling::td[1]")) |> to_int()
    tag_code     <- xml_text(xml_find_first(meta_tbl, ".//td[normalize-space()='Tag Code:']/following-sibling::td[1]"))
    release_site <- xml_text(xml_find_first(meta_tbl, ".//td[normalize-space()='Release Site:']/following-sibling::td[1]"))
    
    raw_tb <- html_table(rec_tbl, header = TRUE, trim = TRUE, fill = TRUE) |> as_tibble()
    if (nrow(raw_tb) == 0) return(tibble())
    rec <- clean_names(raw_tb)
    
    return_col <- names(rec)[str_detect(names(rec), "^return[_ ]?year$|^year$")][1]
    est_col    <- names(rec)[str_detect(names(rec), "estimated|x_estimated|number_estimated|_est$")][1]
    obs_col    <- names(rec)[str_detect(names(rec), "observed|x_observed|number_observed|_obs$")][1]
    fish_col   <- names(rec)[str_detect(names(rec), "^psc.*fishery$|^fishery$")][1]
    site_col   <- names(rec)[str_detect(names(rec), "^site$|^site_name$")][1]
    if (is.na(return_col)) return(tibble())
    
    rec |>
      transmute(
        tag_code     = tag_code,
        brood_year   = brood_year,
        release_site = release_site,
        return_year  = to_int(.data[[return_col]]),
        age          = return_year - brood_year,
        psc_fishery  = if (!is.na(fish_col)) as.character(.data[[fish_col]]) else NA_character_,
        site         = if (!is.na(site_col))  as.character(.data[[site_col]])  else NA_character_,
        n_est        = dplyr::coalesce(
          if (!is.na(est_col)) to_num(.data[[est_col]]) else NA_real_,
          if (!is.na(obs_col)) to_num(.data[[obs_col]]) else NA_real_
        )
      ) |>
      filter(!is.na(return_year), !is.na(age), !is.na(n_est), n_est > 0)
  }))
}

# ================= RUN =================
# Path to your saved page
path <- here("data_raw","american_river_data","rmis_lar.html")

# Parse full page
raw <- rmis_read_html(path)

# Escapement only (no harvest fisheries)
raw_esc <- raw %>% filter(psc_fishery %in% c("Hatchery", "Spawning Ground"))

# ---- Core analysis subset: recovery at AMERICAN or NIMBUS sites, ages 2+ (all years) ----
raw_esc_sites <- raw_esc %>%
  filter(str_detect(site, regex("AMERICAN|NIMBUS", ignore_case = TRUE)),
         !age %in% c(1))

# Output dir
out_dir <- here("output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ================= PLOTS & TABLES =================

# 1) AMERICAN RELEASES ONLY — histogram by age, faceted by YEAR (ALL YEARS)
counts_rel_AMrel_byYear <- raw_esc_sites %>%
  filter(str_detect(release_site, regex("AMERICAN", ignore_case = TRUE))) %>%
  mutate(release_group = "ReleaseSite: AMERICAN") %>%
  group_by(release_group, return_year, age) %>%
  summarise(n = sum(n_est, na.rm = TRUE), .groups = "drop")

write_csv(counts_rel_AMrel_byYear,
          file.path(out_dir, "esc_counts_by_year_age_AMERICAN_RELEASESITE_only_RECOVERY_ANsites_age3plus.csv"))

ggplot(counts_rel_AMrel_byYear,
       aes(x = factor(age), y = n, fill = release_group)) +
  geom_col(position = "dodge") +
  facet_wrap(~return_year, scales = "free_y") +
  labs(x = "Age", y = "# of fish (estimated)", fill = "Release group",
       title = "Escapement age counts (ages 2+) — AMERICAN release site only\n(Recovery sites: AMERICAN/NIMBUS)") +
  theme_minimal(base_size = 11)

# 2) ALL RELEASE GROUPS — histogram by age, faceted by YEAR (ALL YEARS)
counts_rel_groups_byYear <- raw_esc_sites %>%
  mutate(release_group = if_else(str_detect(release_site, regex("AMERICAN", TRUE)),
                                 "ReleaseSite: AMERICAN", "ReleaseSite: OTHER")) %>%
  group_by(release_group, return_year, age) %>%
  summarise(n = sum(n_est, na.rm = TRUE), .groups = "drop")

write_csv(counts_rel_groups_byYear,
          file.path(out_dir, "esc_counts_by_year_age_releaseGroups_RECOVERY_ANsites_age2plus.csv"))

ggplot(counts_rel_groups_byYear,
       aes(x = factor(age), y = n, fill = release_group)) +
  geom_col(position = "dodge") +
  facet_wrap(~return_year, scales = "free_y") +
  labs(x = "Age", y = "# of fish (estimated)", fill = "Release group",
       title = "Escapement age counts (ages 2+) — ALL release groups by year\n(Recovery sites: AMERICAN/NIMBUS)") +
  theme_minimal(base_size = 11)

# 3) ALL RELEASE GROUPS — one figure, ALL YEARS POOLED (counts by age)
counts_rel_groups_allYears <- raw_esc_sites %>%
  mutate(release_group = if_else(str_detect(release_site, regex("AMERICAN", TRUE)),
                                 "ReleaseSite: AMERICAN", "ReleaseSite: OTHER")) %>%
  group_by(release_group, age) %>%
  summarise(n = sum(n_est, na.rm = TRUE), .groups = "drop")

write_csv(counts_rel_groups_allYears,
          file.path(out_dir, "esc_counts_by_age_releaseGroups_POOLED_RECOVERY_ANsites_age3plus.csv"))

ggplot(counts_rel_groups_allYears,
       aes(x = factor(age), y = n, fill = release_group)) +
  geom_col(position = "dodge") +
  labs(x = "Age", y = "# of fish (estimated)", fill = "Release group",
       title = "Escapement age counts (ages 3+) — all years pooled\n(Recovery sites: AMERICAN/NIMBUS)") +
  theme_minimal(base_size = 11)

prop_by_age_all <- counts_rel_groups_allYears %>%
  group_by(age) %>%
  summarise(n_age = sum(n)) %>%
  mutate(prop = n_age / sum(n_age))

# ============ Optional sanity checks ============
# See which years remain after filters (helps explain missing years)
# raw_esc_sites3 %>% distinct(return_year) %>% arrange(return_year)

# Quick totals by year after all filters
# raw_esc_sites3 %>% group_by(return_year) %>% summarise(total = sum(n_est), .groups = "drop") %>% arrange(return_year)

# =========================================
# Escapement ONLY (Hatchery + Spawning Ground), ages 3+,
# NO recovery-site filter; compare release groups
# =========================================

# Base data: escapement, ages >=2, all years (no AMERICAN/NIMBUS recovery-site filter)
raw_esc_age2_allSites <- raw_esc %>%
  filter(!age %in% c(1)) %>%
  mutate(release_group = if_else(
    str_detect(release_site, regex("AMERICAN", ignore_case = TRUE)),
    "ReleaseSite: AMERICAN", "ReleaseSite: OTHER"
  ))

# ---- By-year counts: release group x year x age ----
counts_relGroups_byYear_allSites <- raw_esc_age2_allSites %>%
  group_by(release_group, return_year, age) %>%
  summarise(n = sum(n_est, na.rm = TRUE), .groups = "drop") %>%
  arrange(return_year, release_group, age)

# Plot: histogram-style counts by age, faceted by return year
ggplot(counts_relGroups_byYear_allSites,
       aes(x = factor(age), y = n, fill = release_group)) +
  geom_col(position = "dodge") +
  facet_wrap(~ return_year, scales = "free_y") +
  labs(
    x = "Age", y = "# of fish (estimated)", fill = "Release group",
    title = "Escapement age counts (ages 3+) — by release site group (all recovery sites)"
  ) +
  theme_minimal(base_size = 11)

# Save table
readr::write_csv(
  counts_relGroups_byYear_allSites,
  file.path(out_dir, "esc_counts_by_year_age_releaseGroups_ALLSITES_age3plus.csv")
)

# ---- All-years pooled counts by release group x age ----
counts_relGroups_pooled_allSites <- raw_esc_age3_allSites %>%
  group_by(release_group, age) %>%
  summarise(n = sum(n_est, na.rm = TRUE), .groups = "drop") %>%
  arrange(release_group, age)

# Plot: pooled histogram (no facets)
ggplot(counts_relGroups_pooled_allSites,
       aes(x = factor(age), y = n, fill = release_group)) +
  geom_col(position = "dodge") +
  labs(
    x = "Age", y = "# of fish (estimated)", fill = "Release group",
    title = "Escapement age counts (ages 3+) — all recovery sites, all years pooled"
  ) +
  theme_minimal(base_size = 11)

# Save table
readr::write_csv(
  counts_relGroups_pooled_allSites,
  file.path(out_dir, "esc_counts_by_age_releaseGroups_POOLED_ALLSITES_age3plus.csv")
)
