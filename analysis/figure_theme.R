# ============================================================================
# Shared figure styling
# ============================================================================
# House style, taken from the figures that were already in the manuscript
# (analysis/mcda.R, analysis/figures.R, analysis/presentation_plots.R):
#
#   theme_minimal(base_size = 14)
#   bold black axis titles and axis text
#   bold legend and strip titles
#   white panel and plot background, thin black panel border
#   viridis for every categorical colour or fill scale
#
# Source this at the top of a figure script and use theme_arg() in place of a
# hand-rolled theme() block, so the manuscript figures stay consistent with one
# another and a change here propagates to all of them.
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(viridisLite)
})

#' House theme for manuscript figures
#'
#' @param base_size Base font size. 14 matches the existing figures; drop to 12
#'   only for multi-panel figures where 14 would collide.
#' @param legend Legend position.
#' @param border Draw the black panel border (TRUE in the existing figures).
theme_arg <- function(base_size = 14, legend = "right", border = TRUE) {
  th <- theme_minimal(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", size = base_size + 4, colour = "black"),
      plot.subtitle = element_text(face = "bold", size = base_size - 1,
                                   margin = margin(b = 6), colour = "black"),
      axis.title    = element_text(face = "bold", colour = "black"),
      axis.text     = element_text(face = "bold", size = base_size - 3, colour = "black"),
      strip.text    = element_text(face = "bold", size = base_size - 2, colour = "black"),
      legend.title  = element_text(face = "bold", colour = "black"),
      legend.text   = element_text(size = base_size - 3, colour = "black"),
      legend.position  = legend,
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA)
    )
  if (border) {
    th <- th + theme(panel.border = element_rect(colour = "black", fill = NA,
                                                 linewidth = 0.5))
  }
  th
}

#' n evenly spaced viridis colours, in the range the existing figures use
arg_pal <- function(n, begin = 0.15, end = 0.88, option = "D") {
  viridisLite::viridis(n, begin = begin, end = end, option = option)
}

# ---- Fixed colours for recurring entities ----------------------------------
# Assigned once so the same thing is the same colour in every figure.

# The three TDM formulations, dark to light in the order they appear in the text
TDM_COLS <- setNames(arg_pal(3, begin = 0.10, end = 0.80),
                     c("Bratovich et al. (2020)",
                       "Bartholow & Heasley (2006)",
                       "Martin et al. (2017)"))

TDM_COLS_VARIANT <- setNames(unname(TDM_COLS), c("exp_WF", "exp_SM", "lin_Martin"))

# Alternatives that get singled out in the front-loading figures
ALT_COLS <- setNames(arg_pal(3, begin = 0.10, end = 0.80), c("NB", "PB4", "PB6"))

# Two-level contrasts (before/after, observed/predicted). Kept far apart on the
# viridis ramp so they separate in greyscale as well as colour.
PAIR_COLS <- arg_pal(2, begin = 0.20, end = 0.75)

# Neutral grey for context lines and de-emphasised series
GREY_CTX <- "grey70"
