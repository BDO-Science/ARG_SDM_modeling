# ============================================================================
# SI Tables S2-7 and S2-8: TDM model weight elicitation
# ============================================================================
# Builds the two supporting-information tables promised in the response to
# reviewers from the Round 1 elicitation scoresheet:
#
#   Table S2-7  individual TDM model weightings by panelist (anonymized)
#   Table S2-8  written justifications, verbatim
#
# The scoresheet lives with the manuscript, not in this repo. Point
# SCORESHEET at it (or override via the ARG_SCORESHEET environment variable).
#
# Outputs: output/table_S2-7_panelist_weights.csv
#          output/table_S2-8_justifications.md
# ============================================================================

library(readxl); library(dplyr); library(here)

SCORESHEET <- Sys.getenv(
  "ARG_SCORESHEET",
  unset = file.path(
    "C:/Users/avaisvil/OneDrive - California Department of Water Resources",
    "Documents/projects/american_power_bypass_sdm",
    "2026-07-20 - First Review and Revision/ScoreSheet_ARG_TDM_Round1.xlsx")
)

if (!file.exists(SCORESHEET)) {
  stop("Scoresheet not found: ", SCORESHEET,
       "\nSet ARG_SCORESHEET to its location.")
}

raw <- suppressMessages(read_excel(SCORESHEET, sheet = 1, col_names = FALSE))

# Row 1 is the header; rows 2-6 are the five panelists; row 13 is the aggregate.
panel <- raw[2:6, ] |>
  setNames(c("participant", "name", "martin", "bratovich", "bartholow", "justification")) |>
  mutate(across(c(participant, martin, bratovich, bartholow), as.numeric))

# ---- Table S2-7: anonymized individual weights ------------------------------
# Names are withheld pending consent; keep `name` out of the written output.
tbl_s27 <- panel |>
  transmute(
    Panelist              = paste("Panelist", participant),
    `Martin et al. (2017)`                     = martin,
    `Bratovich et al. (2020) / Water Forum`    = bratovich,
    `Bartholow & Heasley (2006) / SALMOD`      = bartholow,
    Sum                                        = martin + bratovich + bartholow
  )

aggregate_row <- tibble(
  Panelist = "Panel mean (weights used)",
  `Martin et al. (2017)`                  = mean(panel$martin),
  `Bratovich et al. (2020) / Water Forum` = mean(panel$bratovich),
  `Bartholow & Heasley (2006) / SALMOD`   = mean(panel$bartholow),
  Sum = mean(panel$martin) + mean(panel$bratovich) + mean(panel$bartholow)
)

tbl_s27 <- bind_rows(tbl_s27, aggregate_row)

cat("=== Table S2-7: individual TDM weightings ===\n")
print(as.data.frame(tbl_s27), row.names = FALSE)

# Each panelist's allocation should sum to 1
bad <- panel |> filter(abs(martin + bratovich + bartholow - 1) > 1e-8)
if (nrow(bad)) {
  cat("\nWARNING: allocations not summing to 1 for participant(s): ",
      paste(bad$participant, collapse = ", "), "\n")
} else {
  cat("\nAll five allocations sum to 1.\n")
}

# ---- Cross-check: weights stated in prose vs weights recorded ---------------
# Panelist 5's justification quotes 45% / 20% / 35%, which is not what the
# numeric columns record. Flag it rather than silently choosing one.
cat("\n=== Cross-check of weights quoted in the justification text ===\n")
stated <- tibble(
  participant = 5L,
  martin_stated = 0.45, bratovich_stated = 0.20, bartholow_stated = 0.35
)
chk <- panel |> select(participant, martin, bratovich, bartholow) |>
  inner_join(stated, by = "participant")
print(as.data.frame(chk), row.names = FALSE)

alt_means <- c(
  martin    = mean(c(panel$martin[1:4],    0.45)),
  bratovich = mean(c(panel$bratovich[1:4], 0.20)),
  bartholow = mean(c(panel$bartholow[1:4], 0.35))
)
cat(sprintf(
  "\nPanel mean as recorded : Martin %.2f  Bratovich %.2f  Bartholow %.2f\n",
  mean(panel$martin), mean(panel$bratovich), mean(panel$bartholow)))
cat(sprintf(
  "Panel mean if P5 prose : Martin %.2f  Bratovich %.2f  Bartholow %.2f\n",
  alt_means[["martin"]], alt_means[["bratovich"]], alt_means[["bartholow"]]))

# ---- Table S2-8: verbatim justifications ------------------------------------
clean <- function(x) {
  x <- gsub("\\\\r\\\\n", "\n", x)   # literal \r\n sequences stored as text
  x <- gsub("\r\n", "\n", x)
  trimws(x)
}

dir.create(here("output"), showWarnings = FALSE)

md <- c(
  "# Table S2-8. Panelist justifications for TDM model weightings (verbatim)",
  "",
  paste0("Source: `", basename(SCORESHEET), "`, Round 1. ",
         "Panelist numbering matches Table S2-7. Text is reproduced verbatim; ",
         "only line breaks have been normalised."),
  ""
)
for (i in seq_len(nrow(panel))) {
  md <- c(md,
          sprintf("## Panelist %d", panel$participant[i]),
          "",
          sprintf("*Weights: Martin %.2f, Bratovich %.2f, Bartholow %.2f*",
                  panel$martin[i], panel$bratovich[i], panel$bartholow[i]),
          "",
          clean(panel$justification[i]),
          "")
}
writeLines(md, here("output", "table_S2-8_justifications.md"))

write.csv(tbl_s27, here("output", "table_S2-7_panelist_weights.csv"), row.names = FALSE)

cat("\nWrote output/table_S2-7_panelist_weights.csv and",
    "output/table_S2-8_justifications.md\n")
