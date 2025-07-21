pre_spawn_days <- 7              # or make this an input

# pseudo‐code inside your compute_surv() or a new compute_prespawn()
temps_before  <- temps_env[[site_j]][pos - seq_len(pre_spawn_days) + 1L]
cum_dd_presp  <- sum(pmax(temps_before - β_prespawn, 0))

α_prespawn <- 0     # or log(0.95) if you assume a 95% baseline
β_prespawn <- -0.0007

S_prespawn <- exp(α_prespawn + β_prespawn * cum_dd_prespawn)

# old:
redds <- S[t] * P$female_fraction

# new:
# compute S_prespawn[t] for that year/timing/env
redds <- S[t] * P$female_fraction * S_prespawn[t]

prespawn_lookup <- map_dfr(sim_years, function(sim_yr) {
  # for each env & site & spawn date, compute cum_dd_prespawn
  # then S_prespawn as above
  tibble(
    sim_year = sim_yr,
    env      = env_nm,
    variant  = variant,
    S_prespawn = mean(S_prespawn_values)
  )
})

key <- paste(env, variant, sep="_")
S_pre   <- prespawn_lookup_full[[key]][t]    # year‐specific pre‐spawn surv
redds   <- S[t] * P$female_fraction * S_pre
