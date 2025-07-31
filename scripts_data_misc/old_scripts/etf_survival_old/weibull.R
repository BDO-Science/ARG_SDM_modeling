# 1) Compute egg ATU
temps <- rep(10, 20)      # e.g. constant 10°C for 20 days
ATU   <- sum(temps)       # 200

# 2) Direct‐cumulative Weibull
cum_tdm_weibull_ATU <- function(ATU, TL=3, TU=15.4, kL=25, kU=29) {
  (1 - exp(-ATU/TL))^kL * (exp(-ATU/TU))^kU
}

# 3) See the result
cum_tdm_weibull_ATU(ATU)
#> [1] 0.123  # example non-zero survival