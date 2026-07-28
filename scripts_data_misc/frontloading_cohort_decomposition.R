# ============================================================================
# Front-loading mechanism, sub-checks 2 and 3
# ============================================================================
# Task 12 of the revision list. Sub-checks 1 and 4 were done earlier; 1 is
# reproduced here (Part 1) so the whole check lives in one script. This script
# adds the two that were outstanding:
#
#   2. Cross the crossover dates against the spawn-date distribution to get the
#      fraction of redds whose incubation window overlaps the period when a
#      front-loaded alternative is warmer than no-bypass.
#   3. Decompose the deficit by spawn cohort for PB6 vs PB4 under the Martin
#      et al. (2017) formulation, to test whether it sits in late-spawning
#      cohorts rather than being spread evenly.
#
# ---------------------------------------------------------------------------
# Why this does not simply re-run the pipeline's own redd set
# ---------------------------------------------------------------------------
# precompute.R draws the simulated redds with sample() and no seed (section 17),
# never saves the result, and app_data/sim_redds.rds is a stale artifact from an
# older pipeline version. egg_summary.rds therefore cannot be reproduced exactly
# from the repo as it stands - Part 0 quantifies the gap.
#
# So this script computes the survival response deterministically instead: for
# every possible spawn date, under every alternative, from the actual CE-QUAL-W2
# temperatures. That is the mechanism itself, free of any sampling noise. The
# cohort weights needed to turn it into a decomposition are then applied from
# two independent spawn-date distributions (observed carcass surveys 2011-2024,
# and the stale simulated set) and the conclusion is checked against both.
#
# Method note. Both TDM families are additive in the log,
#     -log S = sum_d h(T_d)   over the incubation window,
# so the difference in log-survival between two alternatives decomposes *exactly*
# over calendar days, and the difference in mean survival decomposes exactly over
# spawn cohorts. Nothing in the attribution below is an approximation.
#
# Inputs  : SalmonCountR/app_data/{env_ext_list,sim_redds,egg_summary,results_full}.rds
#           SalmonCountR/app_data/carcassdet_1752789274_15.csv
#           SalmonCountR/functions.R
# Outputs : output/frontloading_crossover_dates.csv
#           output/frontloading_survival_by_spawn_date.csv
#           output/frontloading_incubation_overlap.csv
#           output/frontloading_cohort_decomposition.csv
#           figures/frontloading_cohort_decomposition.png
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(data.table); library(lubridate); library(ggplot2)
  library(patchwork); library(here)
})

source(here("SalmonCountR", "functions.R"))

env_ext_list <- readRDS(here("SalmonCountR", "app_data", "env_ext_list.rds"))
egg_summary  <- readRDS(here("SalmonCountR", "app_data", "egg_summary.rds"))
results_full <- readRDS(here("SalmonCountR", "app_data", "results_full.rds"))
sim_redds    <- as.data.table(readRDS(here("SalmonCountR", "app_data", "sim_redds.rds")))

# env 1-9 = 2011 met year, 10-18 = 2014, 19-27 = 2017, 28-36 = 2020;
# within each block 1=NB, 2=PB1, 3=PB2, 4=PB2b, 5=PB2c, 6=PB3, 7=PB4, 8=PB5, 9=PB6
scen_map  <- c(NB = 1, PB1 = 2, PB2 = 3, PB2b = 4, PB2c = 5,
               PB3 = 6, PB4 = 7, PB5 = 8, PB6 = 9)
scenarios <- names(scen_map)
met_years <- c("2011", "2014", "2017", "2020")

env_key <- data.table(env = as.character(1:36))[
  , i := as.integer(env)][
  , `:=`(scenario = names(scen_map)[match((i - 1) %% 9 + 1, scen_map)],
         met_year = met_years[(i - 1) %/% 9 + 1])][]

env_of <- function(scen, met) env_key[scenario == scen & met_year == met, env]

MARTIN_ALPHA <- 0.026
MARTIN_BETA  <- 12.14
SEASONS      <- 2025:2124          # forecast seasons, matching the adult index
CUT_END      <- "11-30"

exp_par <- list(egg = c(a = 3.408486e-11, b = 1.21122),
                alv = c(a = 1.017550e-10, b = 1.24092))   # WaterForum2020 / Bratovich

# ---- Cumulative-sum series, one per alternative and site --------------------
# A redd starting at index p emerges at the first index where cumulative ATU has
# risen by 958, which findInterval() locates in O(log n); the hazard over that
# window is then a single subtraction. cum vectors carry a leading 0 so any
# window i..j is cum[j + 1] - cum[i].

build_series <- function(env_nm) {
  d <- as.data.table(env_ext_list[[env_nm]])[order(site, Date)]
  lapply(split(d, d$site), function(x) {
    tt <- x$temp; tt[!is.finite(tt)] <- NA_real_
    list(date       = x$Date,
         temp       = tt,
         cum_atu    = c(0, cumsum(pmax(tt, 0))),
         cum_mar    = c(0, cumsum(MARTIN_ALPHA * pmax(tt - MARTIN_BETA, 0))),
         cum_wf_egg = c(0, cumsum(exp_par$egg["a"] * exp(exp_par$egg["b"] * tt))),
         cum_wf_alv = c(0, cumsum(exp_par$alv["a"] * exp(exp_par$alv["b"] * tt))))
  })
}

message("Building temperature series for 36 alternatives ...")
series <- lapply(setNames(names(env_ext_list), names(env_ext_list)), build_series)

#' Survival for a vector of start indices under one alternative/site.
surv_at <- function(ss, i0) {
  n    <- length(ss$temp)
  i_em <- pmax(pmin(findInterval(ss$cum_atu[i0] + total_ATU, ss$cum_atu), n), i0)
  i_h  <- pmax(pmin(findInterval(ss$cum_atu[i0] + egg_ATU,   ss$cum_atu), i_em), i0)
  list(i_emerge   = i_em,
       i_hatch    = i_h,
       nll_martin = ss$cum_mar[i_em + 1] - ss$cum_mar[i0],
       nll_wf     = (ss$cum_wf_egg[i_h + 1]  - ss$cum_wf_egg[i0]) +
                    (ss$cum_wf_alv[i_em + 1] - ss$cum_wf_alv[i_h + 1]))
}

# ---- 0. Reproducibility check ----------------------------------------------
# Recompute egg_summary from the stale sim_redds.rds along the pipeline's own
# path, to document how far apart the two are.

setDT(sim_redds)
redd_pairs <- sim_redds[, .(N = .N), by = .(sim_year, site, spawn_dt)][
  , rdr := as.Date(make_date(fifelse(month(spawn_dt) >= 9L, sim_year, sim_year + 1L),
                             month(spawn_dt), day(spawn_dt)))][
  , site := as.character(site)][]

recomputed <- rbindlist(lapply(names(series), function(e) {
  rbindlist(lapply(unique(redd_pairs$site), function(st) {
    ss <- series[[e]][[st]]; if (is.null(ss)) return(NULL)
    p  <- redd_pairs[site == st]
    i0 <- match(p$rdr, ss$date); keep <- !is.na(i0)
    p  <- p[keep]; i0 <- i0[keep]
    if (!length(i0)) return(NULL)
    s <- surv_at(ss, i0)
    data.table(env = e, sim_year = p$sim_year, N = p$N,
               S_martin = exp(-s$nll_martin))
  }))
}))[, .(recomputed = sum(S_martin * N) / sum(N)), by = .(env, sim_year)]

cmp <- merge(recomputed,
             as.data.table(egg_summary)[variant == "lin_Martin",
                                        .(env = as.character(env), sim_year,
                                          pipeline = mean_cum_surv)],
             by = c("env", "sim_year"))

cat("\n=== Part 0. Can egg_summary.rds be reproduced from the repo? ===\n")
cat(sprintf("  n = %d env-years compared (Martin variant)\n", nrow(cmp)))
cat(sprintf("  mean pipeline %.4f vs mean recomputed %.4f\n",
            mean(cmp$pipeline), mean(cmp$recomputed)))
cat(sprintf("  max |difference| %.4f, correlation %.3f\n",
            max(abs(cmp$pipeline - cmp$recomputed)),
            cor(cmp$pipeline, cmp$recomputed)))
cat("  -> NO. precompute.R draws the redds with sample() and no seed, never\n")
cat("     saves them, and app_data/sim_redds.rds predates the current pipeline.\n")
cat("     Everything below is therefore computed deterministically instead.\n")

# ---- 1. Crossover dates (sub-check 1) --------------------------------------
hazel <- rbindlist(lapply(names(env_ext_list), function(e) {
  d <- as.data.table(env_ext_list[[e]])[site == "AveHazel" & month(Date) %in% c(10, 11)]
  d[, .(env = e, Date, temp)]
}))
hazel <- merge(hazel, env_key[, .(env, scenario, met_year)], by = "env")
hazel[, md := format(Date, "%m-%d")]

daily <- hazel[is.finite(temp), .(temp = mean(temp)), by = .(scenario, md)]
daily <- merge(daily, daily[scenario == "NB", .(md, temp_nb = temp)], by = "md")
daily[, `:=`(delta = temp - temp_nb,
             doy   = as.integer(as.Date(paste0("2001-", md)) - as.Date("2001-10-01")) + 1L)]
setorder(daily, scenario, doy)

crossover <- daily[scenario != "NB", {
  run <- rev(cumprod(rev(as.integer(delta > 0))))   # 1 while warm holds to Nov 30
  i   <- which(run == 1L)[1]
  .(crossover     = if (is.na(i)) NA_character_ else md[i],
    d_oct16_nov15 = mean(delta[md >= "10-16" & md <= "11-15"]),
    d_nov16_30    = mean(delta[md >= "11-16"]),
    max_delta     = max(delta))
}, by = scenario]
setorder(crossover, scenario)

cat("\n=== Part 1. Hazel Avenue crossover vs no-bypass (sub-check 1) ===\n")
print(as.data.frame(crossover[, lapply(.SD, function(x)
  if (is.numeric(x)) round(x, 3) else x)]), row.names = FALSE)
write.csv(crossover, here("output", "frontloading_crossover_dates.csv"), row.names = FALSE)

cross_md <- setNames(crossover$crossover, crossover$scenario)

# ---- 2. Deterministic survival by spawn date -------------------------------
# Every spawn date Oct 1 - Jan 31 of every forecast season, every alternative,
# every site. Also splits the Martin hazard at the alternative's crossover date.

spawn_grid <- function(env_nm) {
  rbindlist(lapply(names(series[[env_nm]]), function(st) {
    ss  <- series[[env_nm]][[st]]
    md  <- format(ss$date, "%m-%d")
    mo  <- month(ss$date)
    sy  <- year(ss$date) - (mo < 9L)                # season year (Sep-Aug)
    sel <- which(mo %in% c(10L, 11L, 12L, 1L) & sy %in% SEASONS)
    if (!length(sel)) return(NULL)
    s <- surv_at(ss, sel)
    data.table(env = env_nm, site = st, season = sy[sel], md = md[sel],
               i0 = sel, i_emerge = s$i_emerge,
               d_emerge = ss$date[s$i_emerge],
               inc_days = s$i_emerge - sel + 1L,
               S_martin = exp(-s$nll_martin),
               S_wf     = exp(-s$nll_wf))
  }))
}

message("Computing survival for every spawn date under 36 alternatives ...")
grid <- rbindlist(lapply(names(series), spawn_grid))
grid <- merge(grid, env_key[, .(env, scenario, met_year)], by = "env")

# Season-averaged survival by scenario / site / calendar spawn date
surv_md <- grid[, .(S_martin = mean(S_martin), S_wf = mean(S_wf),
                    inc_days = mean(inc_days)),
                by = .(scenario, site, md)]
write.csv(surv_md[order(scenario, site, md)],
          here("output", "frontloading_survival_by_spawn_date.csv"), row.names = FALSE)

# ---- 2b. Sub-check 2: overlap with the warm-anomaly window ------------------
# For each spawn date, split the change in Martin -log(survival) relative to NB
# into the part accumulated before the alternative's crossover date and the part
# after it. The split is exact.

split_at_crossover <- function(sc) {
  cut_md <- cross_md[[sc]]
  rbindlist(lapply(met_years, function(m) {
    ea <- env_of(sc, m); en <- env_of("NB", m)
    rbindlist(lapply(names(series[[ea]]), function(st) {
      sa <- series[[ea]][[st]]; sn <- series[[en]][[st]]
      if (is.null(sn)) return(NULL)
      g  <- grid[env == ea & site == st]
      i0 <- g$i0
      sA <- surv_at(sa, i0); sN <- surv_at(sn, i0)
      cutd  <- as.Date(paste0(g$season, "-", cut_md))
      endd  <- as.Date(paste0(g$season, "-", CUT_END))
      i_cut <- pmax(pmin(findInterval(as.integer(cutd), as.integer(sa$date)), sA$i_emerge), i0 - 1L)
      i_cn  <- pmax(pmin(findInterval(as.integer(cutd), as.integer(sn$date)), sN$i_emerge), i0 - 1L)
      data.table(
        scenario = sc, met_year = m, site = st, season = g$season, md = g$md,
        # still in the gravel when the alternative turns warmer than NB?
        overlaps = sa$date[sA$i_emerge] >= cutd & g$i0 <= as.integer(endd),
        pre  = (sa$cum_mar[i_cut + 1] - sa$cum_mar[i0]) -
               (sn$cum_mar[i_cn  + 1] - sn$cum_mar[i0]),
        post = (sa$cum_mar[sA$i_emerge + 1] - sa$cum_mar[i_cut + 1]) -
               (sn$cum_mar[sN$i_emerge + 1] - sn$cum_mar[i_cn  + 1]),
        S_alt = exp(-sA$nll_martin), S_nb = exp(-sN$nll_martin))
    }))
  }))
}

message("Splitting the Martin hazard at each alternative's crossover date ...")
splits <- rbindlist(lapply(setdiff(scenarios, "NB"), split_at_crossover))

# ---- 2c. Spawn-cohort weights ----------------------------------------------
# Two independent distributions. Neither is the pipeline's own draw (which is
# unrecoverable), so the decomposition is reported under both.

carcass_raw <- fread(here("SalmonCountR", "app_data", "carcassdet_1752789274_15.csv"))
carcass <- carcass_raw[, .(Date = as.Date(surveydate), section)][
  , spawn_dt := Date - 7L][
  , site := fifelse(section %in% c("NB", "W", "1a", "1b", "1a/1b", "2"), "AveHazel",
            fifelse(section == "3", "AveWatt", NA_character_))][
  !is.na(site) & !is.na(spawn_dt) & month(spawn_dt) %in% c(10L, 11L, 12L, 1L)]

w_obs <- carcass[, .(N = .N), by = .(site, md = format(spawn_dt, "%m-%d"))]
w_obs[, w := N / sum(N)]

w_sim <- sim_redds[month(spawn_dt) %in% c(10L, 11L, 12L, 1L),
                   .(N = .N), by = .(site = as.character(site),
                                     md = format(spawn_dt, "%m-%d"))]
w_sim[, w := N / sum(N)]

weightings <- list(observed = w_obs, simulated = w_sim)

cohort_bin <- function(md) {
  b <- fifelse(md < "10-16", "Oct 1-15",
       fifelse(md < "11-01", "Oct 16-31",
       fifelse(md < "11-16", "Nov 1-15",
       fifelse(md < "12-01", "Nov 16-30",
       fifelse(md < "12-16", "Dec 1-15", "Dec 16 - Jan")))))
  factor(b, levels = c("Oct 1-15", "Oct 16-31", "Nov 1-15", "Nov 16-30",
                       "Dec 1-15", "Dec 16 - Jan"))
}

overlap_summary <- rbindlist(lapply(names(weightings), function(wn) {
  W <- weightings[[wn]][, .(site, md, w)]
  x <- merge(splits[, .(overlaps = mean(overlaps), pre = mean(pre), post = mean(post)),
                    by = .(scenario, site, md)], W, by = c("site", "md"))
  x[, .(weighting = wn,
        frac_overlapping = sum(w * overlaps) / sum(w),
        d_nll_pre        = sum(w * pre)  / sum(w),
        d_nll_post       = sum(w * post) / sum(w)),
    by = scenario]
}))[, `:=`(d_nll_net = d_nll_pre + d_nll_post,
           crossover = cross_md[scenario])][
  , pct_offset := fifelse(d_nll_pre < 0, 100 * -d_nll_post / d_nll_pre, NA_real_)][]
setorder(overlap_summary, weighting, scenario)

cat("\n=== Part 2. Sub-check 2: incubation overlap with the warm anomaly ===\n")
cat("d_nll = change in Martin -log(survival) relative to no-bypass.\n")
cat("Negative before the crossover = the October cooling benefit.\n")
cat("Positive after = the late-November penalty.\n")
cat("pct_offset = how much of the benefit is handed back after the crossover.\n\n")
print(as.data.frame(overlap_summary[, .(weighting, scenario, crossover,
  frac_overlapping = round(frac_overlapping, 4),
  d_nll_pre = round(d_nll_pre, 5), d_nll_post = round(d_nll_post, 5),
  d_nll_net = round(d_nll_net, 5), pct_offset = round(pct_offset, 1))]),
  row.names = FALSE)
write.csv(overlap_summary, here("output", "frontloading_incubation_overlap.csv"),
          row.names = FALSE)

# ---- 3. Sub-check 3: cohort decomposition ----------------------------------
# The alternatives share one redd distribution, so the difference in mean
# egg-to-fry survival decomposes exactly by cohort:
#   dS = sum_c w_c * (S_c,A - S_c,B),  w_c = redd share of cohort c.

surv_by_md <- grid[, .(S = mean(S_martin)), by = .(scenario, site, md)]

cohort_table <- function(sc_a, sc_b, wn) {
  W <- weightings[[wn]][, .(site, md, w)]
  a <- surv_by_md[scenario == sc_a, .(site, md, S_a = S)]
  b <- surv_by_md[scenario == sc_b, .(site, md, S_b = S)]
  x <- merge(merge(a, b, by = c("site", "md")), W, by = c("site", "md"))
  x[, cohort := cohort_bin(md)]
  wt <- sum(x$w)
  out <- x[, .(redd_share = sum(w) / wt,
               S_a        = sum(w * S_a) / sum(w),
               S_b        = sum(w * S_b) / sum(w),
               contrib    = sum(w * (S_a - S_b)) / wt),
           by = cohort][order(cohort)]
  out[, `:=`(share_of_gap = contrib / sum(contrib),
             comparison   = paste(sc_a, "-", sc_b), weighting = wn)][]
}

comparisons <- list(c("PB6", "PB4"), c("PB6", "NB"), c("PB4", "NB"))
cohorts <- rbindlist(lapply(names(weightings), function(wn)
  rbindlist(lapply(comparisons, function(p) cohort_table(p[1], p[2], wn)))))

cat("\n=== Part 3. Sub-check 3: cohort decomposition under Martin ===\n")
for (wn in names(weightings)) for (cmp in unique(cohorts$comparison)) {
  x <- cohorts[weighting == wn & comparison == cmp]
  cat(sprintf("\n-- %s, %s spawn-date weights -- total difference in mean egg-to-fry survival: %+0.5f\n",
              cmp, wn, sum(x$contrib)))
  print(as.data.frame(x[, .(cohort, redd_share = round(redd_share, 4),
                            S_first = round(S_a, 4), S_second = round(S_b, 4),
                            contribution = round(contrib, 5),
                            share_of_gap = round(share_of_gap, 3))]), row.names = FALSE)
}
write.csv(cohorts, here("output", "frontloading_cohort_decomposition.csv"), row.names = FALSE)

# ---- 4. Is egg survival really the channel? --------------------------------
chan <- as.data.table(results_full)[variant == "lin_Martin"][
  , env := as.character(env)]
chan <- merge(chan, env_key[, .(env, scenario, met_year)], by = "env")
chan <- chan[year > 2024][order(env, year)][, tail(.SD, 20), by = env][
  , .(egg_surv = mean(egg_surv, na.rm = TRUE),
      pre_spawn = mean(pre_spawn, na.rm = TRUE),
      spawners  = median(spawners, na.rm = TRUE)), by = scenario][
  scenario %in% c("NB", "PB4", "PB6")]
setorder(chan, scenario)

cat("\n=== Part 4. Channels through which the alternative acts (Martin only) ===\n")
print(as.data.frame(chan[, lapply(.SD, function(x)
  if (is.numeric(x)) round(x, 5) else x)]), row.names = FALSE)
cat(sprintf("\n  egg survival PB6/PB4 = %.4f   pre-spawn PB6/PB4 = %.4f\n",
            chan[scenario == "PB6", egg_surv] / chan[scenario == "PB4", egg_surv],
            chan[scenario == "PB6", pre_spawn] / chan[scenario == "PB4", pre_spawn]))

# ---- 5. What the Martin-only corner actually contains -----------------------
# Task 1 found the top-ranked alternative flips to NB at 100% Martin weight.
# This is the context for that result: under Martin alone every alternative's
# projected population is functionally extinct, so the Chinook axis is a
# min-max normalisation over a range of a few dozen fish.

med20 <- as.data.table(results_full)[, env := as.integer(env)][
  , `:=`(scenario = names(scen_map)[match((env - 1) %% 9 + 1, scen_map)],
         met_year = met_years[(env - 1) %/% 9 + 1])][
  year > 2024][order(env, variant, year)][, tail(.SD, 20), by = .(env, variant)][
  , .(med = median(spawners, na.rm = TRUE)), by = .(scenario, met_year, variant)]

tdm_w <- c(exp_WF = .51, exp_SM = .24, lin_Martin = .25)

by_variant <- dcast(med20[, .(idx = sum(med * 0.25)), by = .(scenario, variant)],
                    scenario ~ variant, value.var = "idx")
by_variant <- merge(by_variant,
                    med20[, .(model_averaged = sum(med * tdm_w[variant] * 0.25)),
                          by = scenario], by = "scenario")
setorder(by_variant, -model_averaged)

cat("\n=== Part 5. Adult index by TDM variant (median of final 20 forecast years) ===\n")
print(as.data.frame(by_variant[, lapply(.SD, function(x)
  if (is.numeric(x)) round(x, 1) else x)]), row.names = FALSE)

mm <- med20[variant == "lin_Martin", .(idx = sum(med * 0.25)), by = scenario][
  , chinook_norm := (idx - min(idx)) / (max(idx) - min(idx))][order(-idx)]
cat("\n  Martin-only: the whole nine-alternative range is",
    sprintf("%.1f to %.1f spawners.\n", min(mm$idx), max(mm$idx)))
print(as.data.frame(mm[, .(scenario, adult_index = round(idx, 1),
                           chinook_norm = round(chinook_norm, 3))]), row.names = FALSE)
write.csv(by_variant, here("output", "frontloading_index_by_variant.csv"), row.names = FALSE)

# ---- 6. Figure --------------------------------------------------------------
pal <- c(NB = "#440154", PB4 = "#1F9E89", PB6 = "#D55E00")

pA <- ggplot(daily[scenario != "NB"], aes(doy, delta, group = scenario)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_line(colour = "grey82", linewidth = 0.5) +
  geom_line(data = daily[scenario %in% c("PB4", "PB6")],
            aes(colour = scenario), linewidth = 0.9) +
  scale_colour_manual(values = pal[c("PB4", "PB6")], name = NULL) +
  scale_x_continuous(breaks = c(1, 16, 32, 47, 61),
                     labels = c("Oct 1", "Oct 16", "Nov 1", "Nov 16", "Nov 30"),
                     expand = c(0.01, 0)) +
  labs(subtitle = "(a) Hazel Avenue daily temperature minus no-bypass (grey: other alternatives)",
       x = NULL, y = expression(Delta*" temperature ("*degree*"C)")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.4),
        plot.subtitle = element_text(face = "bold", size = 10.5),
        legend.position = "top")

md_lv <- sort(unique(surv_by_md$md))
md_ord <- md_lv[order(ifelse(substr(md_lv, 1, 2) == "01", paste0("13", substr(md_lv, 4, 5)),
                             paste0(substr(md_lv, 1, 2), substr(md_lv, 4, 5))))]
sm <- surv_by_md[site == "AveHazel" & scenario %in% c("NB", "PB4", "PB6")]
sm[, x := match(md, md_ord)]

x_breaks <- match(c("10-01", "11-01", "12-01", "01-01"), md_ord)
x_labs   <- c("Oct 1", "Nov 1", "Dec 1", "Jan 1")

# Levels first: the season splits into a lethal Oct-Nov window and a safe
# December onwards. Only NB is drawn, because at this scale the three
# alternatives are visually identical - which is the point panel (c) makes.
pB <- ggplot(sm[scenario == "NB"], aes(x, S)) +
  geom_line(linewidth = 0.9, colour = pal[["NB"]]) +
  scale_x_continuous(breaks = x_breaks, labels = x_labs, expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(subtitle = "(b) Martin egg-to-fry survival by spawn date under no-bypass, Hazel Avenue (mean of 100 forecast seasons)",
       x = NULL, y = "Egg-to-fry survival") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.4),
        plot.subtitle = element_text(face = "bold", size = 10.5))

# Differences, which are what the decomposition acts on and are far too small
# to read off panel (b).
sm_d <- merge(sm[scenario != "NB", .(scenario, x, S)],
              sm[scenario == "NB", .(x, S_nb = S)], by = "x")
sm_d[, dS := S - S_nb]

pC <- ggplot(sm_d, aes(x, dS, colour = scenario)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = pal[c("PB4", "PB6")], name = NULL) +
  scale_x_continuous(breaks = x_breaks, labels = x_labs, expand = c(0.01, 0)) +
  labs(subtitle = "(c) Same curves minus no-bypass: cohorts spawning in November lose, cohorts spawning in October gain",
       x = NULL, y = expression(Delta*" survival vs NB")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.4),
        plot.subtitle = element_text(face = "bold", size = 10.5),
        legend.position = "top")

pD <- overlap_summary[weighting == "observed",
                      .(scenario, Before = d_nll_pre, After = d_nll_post)] |>
  melt(id.vars = "scenario", variable.name = "period", value.name = "d_nll") |>
  ggplot(aes(scenario, d_nll, fill = period)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_col(colour = "grey25", linewidth = 0.3, width = 0.72) +
  scale_fill_manual(values = c(Before = "#1F9E89", After = "#D55E00"),
                    name = "relative to crossover") +
  labs(subtitle = "(d) Change in Martin -log(survival) vs no-bypass, split at the crossover date",
       x = NULL, y = expression(Delta*" -log S")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.4),
        plot.subtitle = element_text(face = "bold", size = 10.5),
        legend.position = "top")

coh_plot <- cohorts[weighting == "observed" & comparison == "PB6 - PB4"]
pad <- 0.35 * diff(range(c(0, coh_plot$contrib)))

pE <- ggplot(coh_plot, aes(cohort, contrib)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_col(fill = "#D55E00", colour = "grey25", linewidth = 0.3, width = 0.72) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * share_of_gap)),
            vjust = -0.6, size = 3, colour = "grey15") +
  scale_y_continuous(limits = c(min(coh_plot$contrib) - pad, pad)) +
  labs(subtitle = "(e) PB6 minus PB4 egg-to-fry survival by spawn cohort (Martin); labels are each cohort's share of the gap",
       x = "Spawn cohort", y = "Contribution to the gap") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.4),
        plot.subtitle = element_text(face = "bold", size = 10.5))

fig <- pA / pB / pC / pD / pE + plot_layout(heights = c(1, 0.8, 0.9, 0.9, 0.9))

dir.create(here("figures"), showWarnings = FALSE)
dir.create(here("output"),  showWarnings = FALSE)
ggsave(here("figures", "frontloading_cohort_decomposition.png"), fig,
       width = 10, height = 15, dpi = 300, bg = "white")

cat("\nWrote output/frontloading_crossover_dates.csv,",
    "output/frontloading_survival_by_spawn_date.csv,",
    "output/frontloading_incubation_overlap.csv,",
    "output/frontloading_cohort_decomposition.csv,",
    "output/frontloading_index_by_variant.csv,",
    "figures/frontloading_cohort_decomposition.png\n")
