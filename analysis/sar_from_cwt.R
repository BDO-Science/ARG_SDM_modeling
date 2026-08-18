# ============================================================================
# Smolt-to-adult return ratio (SAR) from American River coded wire tag returns
# ============================================================================
# Derives every SAR statistic quoted in SI Section S2.6, so the numbers in the
# text can be checked against the data instead of taken on trust. This is the
# citation the SI needs (OUTSTANDING_ITEMS.md D3): the source is the hand-built
# `app_data/SAR LAR Releases.xlsx`, filtered to American River release groups,
# where SAR = expanded returns / number released.
#
# Runs in seconds. Reads only saved data; does not touch precompute.R.
#
# Outputs: output/sar_by_brood.csv, output/sar_statistics.csv
#
# NOTE the two windows. The mean, standard deviation and group count are over
# ALL American release groups (brood years 2008-2019, n = 27), because 0.25% is
# the optimizer starting value the SI ties to that mean. The pooled brood-year
# RANGE must be quoted over the same window: 0.02% (BY2017) to 0.95% (BY2010).
# Quoting 0.56% (BY2016) as the maximum silently drops brood years 2008-2010 --
# see the check at the bottom of this script.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readxl); library(janitor); library(stringr); library(here)
})

xlsx <- here("SalmonCountR", "app_data", "SAR LAR Releases.xlsx")

# `sar_percent` in this workbook is a copy of `sar`, NOT a percentage
# (OUTSTANDING_ITEMS.md G3), so it is dropped here rather than used.
am <- read_excel(xlsx) %>%
  clean_names() %>%
  select(-any_of("sar_percent")) %>%
  filter(str_detect(toupper(release_location), "AMERICAN"))

stopifnot(nrow(am) > 0, all(c("by", "sar", "number_released") %in% names(am)))

# ---- Pooled SAR by brood year -----------------------------------------------
# Pooled, not the mean of the release-group ratios: returns and releases are
# summed within a brood year first, so large release groups carry their weight.
by_brood <- am %>%
  group_by(brood_year = by) %>%
  summarise(
    release_groups = n(),
    number_released = sum(number_released, na.rm = TRUE),
    expanded_returns = sum(expanded_fsamp_and_fprod, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(sar_pooled = expanded_returns / number_released,
         sar_pooled_percent = 100 * sar_pooled)

# ---- Statistics over release groups -----------------------------------------
set.seed(123)                                  # bootstrap CI is seed-dependent
boot_mean <- replicate(10000, mean(sample(am$sar, replace = TRUE), na.rm = TRUE))

stats <- tibble::tibble(
  statistic = c("release groups", "first brood year", "last brood year",
                "mean SAR", "sd of SAR", "median SAR",
                "bootstrap 95% CI lower", "bootstrap 95% CI upper",
                "pooled min", "pooled min brood year",
                "pooled max", "pooled max brood year"),
  value = c(nrow(am), min(am$by), max(am$by),
            mean(am$sar, na.rm = TRUE), sd(am$sar, na.rm = TRUE), median(am$sar, na.rm = TRUE),
            unname(quantile(boot_mean, 0.025)), unname(quantile(boot_mean, 0.975)),
            min(by_brood$sar_pooled), by_brood$brood_year[which.min(by_brood$sar_pooled)],
            max(by_brood$sar_pooled), by_brood$brood_year[which.max(by_brood$sar_pooled)])
)

dir.create(here("output"), showWarnings = FALSE)
write.csv(by_brood, here("output", "sar_by_brood.csv"), row.names = FALSE)
write.csv(stats,    here("output", "sar_statistics.csv"), row.names = FALSE)

cat("\n=== Pooled SAR by brood year, American River release groups ===\n")
print(as.data.frame(by_brood %>%
  transmute(brood_year, release_groups, number_released,
            sar_percent = round(sar_pooled_percent, 3))))

cat(sprintf("\n=== Release-group statistics (brood years %d-%d, n = %d) ===\n",
            min(am$by), max(am$by), nrow(am)))
cat(sprintf("  mean SAR   : %.4f  (%.3f%%)\n", mean(am$sar), 100 * mean(am$sar)))
cat(sprintf("  sd         : %.5f\n", sd(am$sar)))
cat(sprintf("  bootstrap  : %.4f to %.4f  (95%% percentile CI on the mean)\n",
            quantile(boot_mean, 0.025), quantile(boot_mean, 0.975)))

# ---- Check the SI text against the data -------------------------------------
cat("\n=== SI Section S2.6 text, checked against the data ===\n")
check <- function(claim, got, want, tol) {
  ok <- abs(got - want) <= tol
  cat(sprintf("  %-34s SI says %-9s data give %-9s %s\n", claim, format(want),
              format(round(got, 5)), if (ok) "ok" else "*** MISMATCH ***"))
  ok
}
hi <- by_brood %>% slice_max(sar_pooled)
lo <- by_brood %>% slice_min(sar_pooled)
ok <- c(
  check("release groups",        nrow(am), 27, 0),
  check("first brood year",      min(am$by), 2008, 0),
  check("last brood year",       max(am$by), 2019, 0),
  check("mean SAR, percent",     100 * mean(am$sar), 0.25, 0.005),
  check("sd of SAR",             sd(am$sar), 0.0024, 0.00005),
  check("pooled minimum, percent", 100 * lo$sar_pooled, 0.02, 0.005),
  check("pooled min brood year", lo$brood_year, 2017, 0),
  check("pooled maximum, percent", 100 * hi$sar_pooled, 0.56, 0.005),
  check("pooled max brood year", hi$brood_year, 2016, 0)
)

if (!all(ok)) {
  cat("\n  The last two are the known error. Over brood years 2008-2019 the\n")
  cat(sprintf("  maximum pooled SAR is %.2f%% in brood year %d; %s is the maximum only\n",
              100 * hi$sar_pooled, hi$brood_year, "0.56% (BY2016)"))
  cat("  once brood years 2008-2010 are excluded, which the stated window keeps.\n")
  cat("  Fix the SI to read: ranged from 0.02% (brood year 2017) to 0.95% (brood\n")
  cat("  year 2010). Restricting the window instead would change n to 24 and the\n")
  cat("  mean to 0.22%, breaking the tie to the 0.0025 optimizer starting value.\n")
} else {
  cat("\n  Every claim in the SI paragraph reproduces.\n")
}

cat("\nWrote output/sar_by_brood.csv and output/sar_statistics.csv\n")
