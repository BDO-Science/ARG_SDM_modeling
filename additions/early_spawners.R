####################################
#FISH SPAWNING BEFORE 11-15 EACH YEAR
#####################################
results_obs_early <- future_map_dfr(sim_years, function(sim_yr) {
  # 1) pull only the redds for this brood_year AND before Nov 15
  red_this <- sim_redds %>%
    filter(sim_year == sim_yr) %>%
    # keep spawn_dt ≤ yyyy-11-15
    filter(spawn_dt <= as.Date(sprintf("%04d-11-15", sim_yr)))
  
  # if none, skip
  if (nrow(red_this)==0) {
    return(tibble(
      sim_year      = sim_yr,
      env           = NA_character_,
      variant       = NA_character_,
      method        = "early",
      mean_cum_surv = NA_real_
    ))
  }
  
  # 2) proceed exactly as in results_obs_fast, but on this red_this subset
  map_dfr(names(env_ext_list), function(env_nm) {
    date_idx   <- date_idx_list[[env_nm]]
    temps_list <- site_temps_list[[env_nm]]
    env_dates  <- env_ext_list[[env_nm]]$Date
    
    rdrs <- make_date(
      year  = if_else(month(red_this$spawn_dt) >= 9L, sim_yr, sim_yr + 1L),
      month = month(red_this$spawn_dt),
      day   = day(  red_this$spawn_dt)
    )
    
    in_window <- (rdrs >= min(env_dates)) & (rdrs <= max(env_dates))
    rdrs  <- rdrs[in_window]
    sites <- red_this$site[in_window]
    
    map_dfr(seq_len(nrow(tdm_defs)), function(i) {
      model   <- tdm_defs$model[i]
      calib   <- tdm_defs$calib[i]
      variant <- tdm_defs$variant[i]
      
      survs <- map_dbl(seq_along(rdrs), function(j) {
        site_j <- sites[j]
        rdr_j  <- rdrs[j]
        pos    <- date_idx[[site_j]][as.character(rdr_j)]
        if (!isTRUE(pos >= 1))            return(NA_real_)
        temp0 <- temps_list[[site_j]][pos]
        if (!isTRUE(temp0 > 0))           return(NA_real_)
        d_egg   <- egg_model(temp0)
        d_hatch <- hatch_model(temp0)
        td      <- ceiling(d_egg + d_hatch)
        if (!isTRUE(td >= 1))             return(NA_real_)
        slice   <- temps_list[[site_j]][pos + seq_len(td) - 1]
        if (model=="exp") tdm_exp(slice, calib) else tdm_lin_martin(slice)
      })
      
      tibble(
        sim_year      = sim_yr,
        env           = env_nm,
        variant       = variant,
        method        = "early",
        mean_cum_surv = mean(survs, na.rm=TRUE)
      )
    })
  })
}, .options = furrr_options(seed = TRUE))

print(results_obs_early, n = Inf)

ggplot(results_obs_early, aes(x = sim_year, y = mean_cum_surv, color = variant)) +
  geom_line(size = 1) +
  facet_wrap(~ env, ncol = 2) +
  scale_x_continuous(breaks = pretty(results_obs_early$sim_year, 10)) +
  labs(
    title = "Observed‐based TDM: Mean Egg→Emergence Survival by Brood Year",
    x     = "Brood Year",
    y     = "Mean cumulative survival",
    color = "Variant"
  ) +
  theme_minimal(base_size = 14)

ggplot(results_obs_early, aes(x = variant, y = mean_cum_surv)) +
  geom_boxplot(outlier.size = 1) +
  facet_wrap(~ env) +
  labs(
    title = "Distribution of Brood‐Year Survival by Variant & Env",
    x     = "TDM Variant",
    y     = "Mean cumulative survival"
  ) +
  theme_minimal(base_size = 14)
