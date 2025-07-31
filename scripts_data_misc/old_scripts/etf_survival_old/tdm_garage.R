# 3) cross your defs with the time series, compute daily & cumulative survival:
tdm_daily <-
  crossing(
    date = dates,
    temp = t_series,
    tdm_defs
  ) %>%
  arrange(variant, date) %>%
  # 1) give each row within a variant its day‐index for lin_anderson:
  group_by(variant) %>%
  mutate(day = row_number()) %>%
  ungroup() %>%
  # 2) compute surv_daily with pmap_dbl + switch (one row at a time):
  mutate(
    surv_daily = pmap_dbl(
      list(model, temp, calib, stage, day),
      function(model, temp, calib, stage, day) {
        switch(model,
               exp           = tdm_exp(temp, calib, stage),
               lin_martin    = tdm_lin_martin(temp),
               lin_anderson  = tdm_lin_anderson(temp, day = day),
               weibull       = tdm_weibull(temp),
               none          = tdm_none(temp),
               const         = tdm_const(temp),
               NA_real_
        )
      }
    )
  ) %>%
  # 3) cumulative product within each variant:
  group_by(variant) %>%
  mutate(cum_surv = cumprod(surv_daily)) %>%
  ungroup() %>%
  select(model, calib, stage, surv_daily, date, cum_surv, variant)

# 2) Expand into every combination
# … inside your summary_tbl block, replace the `dev0 = compute_combo(…)` bit with:
summary_tbl <- tidyr::crossing(
  tdm_defs,
  egg_m   = names(egg_models),
  hatch_m = names(hatch_models),
  dd_type = c("none","linear","beverton-holt")
) %>%
  rowwise() %>%
  mutate(
    # previous bits…
    survs        = list(simulate_combo(model, calib,
                                       egg_models[[egg_m]],
                                       hatch_models[[hatch_m]])),
    base_surv    = mean(survs, na.rm=TRUE),
    dd_factor    = density_fun(dd_type, N = length(redds)),
    mean_cum_surv= base_surv * dd_factor,
    
    # ←— NEW: compute development stats inline instead of compute_combo()
    dev0 = {
      # get the spawn temps for all redds
      Te  <- t_series[ match(redds, dates) ]
      # days to emergence
      de  <- egg_models[[egg_m]](Te) + hatch_models[[hatch_m]](Te)
      # accumulate thermal units = temp * days
      atu <- de * Te
      tibble(
        mean_dev_days = mean(de),
        mean_ATU      = mean(atu)
      )
    }
  ) %>%
  tidyr::unnest_wider(dev0) %>%
  ungroup() %>%
  # rounding and select as before
  mutate(across(c(mean_dev_days, mean_ATU, mean_cum_surv), ~ round(.x, 3))) %>%
  select(model, calib, egg_m, hatch_m, dd_type,
         mean_dev_days, mean_ATU, mean_cum_surv)


# Print the full table
print(summary_tbl, n = Inf)
invisible(summary_tbl)










# — assume you already have:
#    dates     = seq(Date1, Date2, by="day")
#    t_series  = your real‐world temp series vector of the same length
#    redds     = vector of spawn dates for each redd
#    tdm_defs  = tibble of (model, calib, variant) for each TDM variant
#    egg_models, hatch_models, density_fun, and your tdm_*() functions

# 1) helper for a single redd’s time series
simulate_redd <- function(spawn_date, model, calib, egg_fn, hatch_fn) {
  i0 <- match(spawn_date, dates)
  days_egg   <- egg_fn(   t_series[i0] )
  days_alev  <- hatch_fn( t_series[i0] )
  total_days <- ceiling(days_egg + days_alev)
  idx        <- seq(i0, min(i0 + total_days - 1, length(dates)))
  # compute daily survival
  surv_daily <- map2_dbl(seq_along(idx), idx, ~{
    Tday <- t_series[.y]
    if (model == "exp") {
      stage <- if (.x <= ceiling(days_egg)) "egg" else "alevin"
      tdm_exp(Tday, calib, stage)
    } else if (model=="lin_martin") {
      tdm_lin_martin(Tday)
    } else if (model=="lin_anderson") {
      tdm_lin_anderson(Tday, day = .x)
    } else if (model=="weibull") {
      tdm_weibull(Tday)
    } else if (model=="none") {
      tdm_none(Tday)
    } else if (model=="const") {
      tdm_const(Tday)
    }
  })
  cum_surv <- cumprod(surv_daily)
  tibble(
    redd       = as.character(spawn_date),
    spawn_date = spawn_date,
    date       = dates[idx],
    cum_surv   = cum_surv
  )
}

# 2) simulate for every combo & every redd, then bind
all_ts <- tidyr::crossing(
  tdm_defs,
  egg_m   = names(egg_models),
  hatch_m = names(hatch_models),
  dd_type = c("none","linear","beverton-holt")
) %>%
  rowwise() %>%
  mutate(
    # run each redd through simulate_redd() and tag on the density factor
    ts = list(
      map_dfr(
        redds,
        ~ simulate_redd(.x, model, calib,
                        egg_models[[egg_m]],
                        hatch_models[[hatch_m]])
      ) %>%
        mutate(
          dd_factor = density_fun(dd_type, N = length(redds)),
          cum_surv_dd = cum_surv * dd_factor
        )
    )
  ) %>%
  unnest(ts) %>%
  ungroup()

# 3) now compute the MEAN cumulative survival across redds at each date
mean_ts <- all_ts %>%
  group_by(model, calib, variant, egg_m, dd_type, date) %>%
  summarize(
    mean_surv    = mean(cum_surv),
    mean_surv_dd = mean(cum_surv_dd),
    .groups="drop"
  )

# 4) finally, plot it
ggplot(mean_ts, aes(date, mean_surv, color=variant, group=variant)) +
  geom_line(size=1) +
  facet_grid(egg_m ~ dd_type) +
  scale_y_continuous(limits = c(0,1)) +
  labs(
    title    = "Mean Spawn-to-Emergence Survival Over Time",
    subtitle = "Rows = egg model; cols = density dependence",
    x        = "Calendar Date",
    y        = "Mean Cumulative Survival",
    color    = "TDM Variant"
  ) +
  theme_minimal() +
  theme(legend.position="bottom")


# e.g. pick one TDM variant & density‐dependence to inspect:
sel_variant <- "exp_WF_egg"
sel_dd      <- "beverton-holt"



# First build a daily‐survival + density table per year, date & variant
daily_dd <- map_dfr(year_list, function(yr){
  tibble(
    year = yr$year,
    date = yr$dates,
    temp = yr$temps
  ) %>%
    crossing(variant = tdm_defs$variant, dd_type = dd_tab$dd_type) %>%
    left_join(tdm_defs, by = "variant") %>%
    mutate(
      surv_daily = pmap_dbl(list(model, temp, calib, stage, row_number()),
                            ~ switch(..1,
                                     exp          = tdm_exp(..2, ..3, ..4),
                                     lin_martin   = tdm_lin_martin(..2),
                                     lin_anderson = tdm_lin_anderson(..2, day=..5),
                                     weibull      = tdm_weibull(..2),
                                     none         = tdm_none(..2),
                                     const        = tdm_const(..2),
                                     NA_real_)),
      surv_dd = surv_daily * density_fun(dd_type, N = n_redds)
    )
})

# ───────────────────────────────────────────────────────────────────────────────
# Build mean_ts: average (daily × density) survival by year, variant & date
# ───────────────────────────────────────────────────────────────────────────────
mean_ts <- daily_dd %>%
  group_by(year, variant, dd_type, date) %>%
  summarise(
    mean_surv_dd = mean(surv_dd, na.rm = TRUE),
    .groups = "drop"
  )

mean_ts %>%
  filter(variant == sel_variant,
         dd_type == sel_dd) %>%
  ggplot(aes(date, mean_surv_dd, color = variant)) +
  geom_line(size = 1) +
  facet_wrap(~ year, ncol = 5) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(
    title = paste("Cumulative Survival:", sel_variant, "+", sel_dd),
    x     = "Date",
    y     = "Mean Cumulative Survival"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Then heatmap of surv_dd
ggplot(daily_dd, aes(x = format(date, "%m-%d"), y = factor(year), fill = surv_dd)) +
  geom_raster() +
  scale_fill_viridis_c(limits = c(0,1)) +
  facet_wrap(~ variant + dd_type, ncol = length(unique(daily_dd$dd_type))) +
  labs(x = "Day of Water‐Year", y = "Year", fill = "Daily × Density") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
    panel.spacing = unit(0.1, "lines")
  )



