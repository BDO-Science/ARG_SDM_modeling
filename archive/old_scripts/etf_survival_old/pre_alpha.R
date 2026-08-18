library(dplyr)
library(purrr)

# 1) Calibrations pulled from your HTML snippet:
exp_vals <- c("waterforum", "zeug", "salmod")
lin_vals <- c("martin", "anderson")

# 2) Build your test tibble
tests <- tibble(
  mortality = c(rep("exp", length(exp_vals)), rep("lin", length(lin_vals))),
  value     = c(exp_vals, lin_vals)
)

# 3) Fetch survival for each combo
test_results <- tests %>%
  mutate(
    surv = map2_dbl(mortality, value, ~ {
      field <- if (.x == "exp") "paramdefault" else "linparamdefault"
      fetch_one(.x, field, .y)
    })
  )

print(test_results)
