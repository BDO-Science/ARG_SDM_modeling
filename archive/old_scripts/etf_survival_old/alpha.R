# ------------------------------------------------------------------------------
# Full combined egg-development + temperature-dependent mortality (TDM) code
# ------------------------------------------------------------------------------

# 0) Libraries
library(httr)
library(rvest)
library(xml2)
library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

# 1) Fetch function
fetch_tdm_eggcombo <- function(
    spawn_date,
    redds,
    eggperredd,
    temp_point,
    
    # Spawning method
    spawn_method       = c("text", "database"),
    db_source          = "American River spawned female carcasses",
    db_year            = format(Sys.Date(), "%Y"),
    carcasses_per_redd = 1,
    timing_offset      = 0,
    
    # TDM + Egg model
    tdm_model    = c("exp","lin","weibull","none","constant"),
    eggmodel     = c("zeug", "Mechanistic", "JensenC", "BeachamMurray", "USGS"),
    
    # Calibration and parameters
    tdm_calib    = NULL,
    pre_alpha    = NULL, pre_beta  = NULL,
    post_alpha   = NULL, post_beta = NULL,
    lin_alpha    = NULL, lin_beta  = NULL,
    weib_TL      = NULL, weib_TU   = NULL,
    weib_kL      = NULL, weib_kU   = NULL,
    
    # Egg model parameters
    atus           = 958,
    eggparameter   = 200,
    hatchmechanism = "hatchatu",
    atushatch      = 417,
    
    # Output options
    densdeptype        = "none",
    rtnCSV             = "1",
    submit_button_text = "Run Emergence (to an output tab)"
) {
  spawn_method <- match.arg(spawn_method)
  tdm_model    <- match.arg(tdm_model)
  eggmodel     <- match.arg(eggmodel)
  year_str     <- format(spawn_date, "%Y")
  doy          <- as.integer(format(spawn_date, "%j"))
  
  body <- list(
    version        = "2.1.5",
    tempsource     = "dbtemp",
    tempyear       = year_str,
    tempprojpoint  = temp_point,
    units          = "centigrade",
    densdeptype    = densdeptype,
    rtnCSV         = rtnCSV,
    userfish       = "sacpas",
    whereoutput    = "Alt.Passage Results",
    modeltemptype  = "degC",
    showdetails    = "yes",
    showcumulative = "no",
    go             = submit_button_text,
    eggmodel       = eggmodel
  )
  
  # Spawning
  if (spawn_method == "text") {
    body$spawnsource     <- "spawntext"
    body$spawntextstring <- paste0("Day,Redds\n", doy, ",", redds)
    body$redddatacolumn  <- "2"
    body$reddyear        <- year_str
    body$eggperredd      <- as.character(eggperredd)
    body$adultperredd    <- "1"
  } else {
    body$spawnsource   <- "spawndb"
    body$spawndbsource <- db_source
    body$reddyear      <- as.character(db_year)
    body$adultperredd  <- as.character(carcasses_per_redd)
    body$timingoffset  <- as.character(timing_offset)
  }
  
  # TDM model
  body$mortality <- tdm_model
  
  if (tdm_model == "exp") {
    if (!is.null(tdm_calib)) {
      body$paramdefault <- case_when(
        tdm_calib == "Water Forum 2020"        ~ "waterforum",
        tdm_calib == "Zeug 2012, Cramer 2010"  ~ "zeug",
        tdm_calib == "OSALMOD 2006"            ~ "osalmod",
        TRUE                                   ~ "waterforum"
      )
    } else {
      if (!is.null(pre_alpha))  body$pre_hatch_alpha  <- as.character(pre_alpha)
      if (!is.null(pre_beta))   body$pre_hatch_beta   <- as.character(pre_beta)
      if (!is.null(post_alpha)) body$post_hatch_alpha <- as.character(post_alpha)
      if (!is.null(post_beta))  body$post_hatch_beta  <- as.character(post_beta)
    }
  }
  if (tdm_model == "lin") {
    if (!is.null(tdm_calib)) {
      body$linparamdefault <- case_when(
        tdm_calib == "Martin (2016, 2017)"     ~ "martin",
        tdm_calib == "Anderson (2022)"         ~ "anderson",
        TRUE                                   ~ "martin"
      )
    } else {
      if (!is.null(lin_alpha)) body$lin_alpha <- as.character(lin_alpha)
      if (!is.null(lin_beta))  body$lin_beta  <- as.character(lin_beta)
    }
  }
  if (tdm_model == "weibull") {
    if (!is.null(weib_TL)) body$weib_eggs_TL <- as.character(weib_TL)
    if (!is.null(weib_TU)) body$weib_eggs_TU <- as.character(weib_TU)
    if (!is.null(weib_kL)) body$weib_eggs_kL <- as.character(weib_kL)
    if (!is.null(weib_kU)) body$weib_eggs_kU <- as.character(weib_kU)
  }
  
  # Egg model settings
  body$atus           <- as.character(atus)
  body$eggparameter   <- as.character(eggparameter)
  body$hatchmechanism <- hatchmechanism
  body$atushatch      <- as.character(atushatch)
  
  # Submit
  res <- POST("https://cbr.washington.edu/sac-bin/grow/emergecontrols.pl",
              encode = "form", body = body)
  stop_for_status(res)
  
  # Parse CSV
  doc      <- read_html(content(res, "text", encoding = "UTF-8"))
  csv_href <- doc %>% html_node(xpath = "//a[contains(., 'download csv')]") %>% html_attr("href")
  csv_url  <- url_absolute(csv_href, "https://cbr.washington.edu")
  txt      <- content(GET(csv_url), "text", encoding = "UTF-8")
  lines    <- strsplit(txt, "\r?\n")[[1]]
  hdr      <- grep("^SpawnDay,", lines)[1]
  df       <- read_csv(I(lines[hdr:length(lines)]), show_col_types = FALSE)
  
  list(
    eggmodel = eggmodel,
    tdm_model = tdm_model,
    cum_surv = df$Survival[1],
    details = df
  )
}

# 2) Calibration presets
tdm_calibrations <- tribble(
  ~tdm_model, ~tdm_calib,                    ~params,
  "exp",      "Water Forum 2020",            list(pre_alpha = 3.40849e-11, pre_beta = 1.2112,
                                                  post_alpha = 1.01755e-10, post_beta = 1.2409),
  "exp",      "Zeug 2012, Cramer 2010",      list(pre_alpha = 5.4e-10, pre_beta = 1.13,
                                                  post_alpha = 6.1e-10, post_beta = 1.31),
  "exp",      "OSALMOD 2006",                list(pre_alpha = 4.77e-10, pre_beta = 1.04,
                                                  post_alpha = 5.88e-10, post_beta = 1.22),
  "lin",      "Martin (2016, 2017)",         list(lin_alpha = 0.0239, lin_beta = 12.07),
  "lin",      "Anderson (2022)",             list(lin_alpha = 0.03, lin_beta = 14.5)
)

# 3) Create grid of combinations
egg_models <- c("zeug", "Mechanistic", "JensenC", "BeachamMurray", "USGS")
other_tdm_models <- tibble(tdm_model = c("weibull", "none", "constant"),
                           tdm_calib = "default", params = list(NULL, NULL, NULL))

full_combos <- expand_grid(
  eggmodel = egg_models,
  join_tdm = bind_rows(tdm_calibrations, other_tdm_models)
) %>% unnest_wider(join_tdm)

# 4) Run all combinations
results <- full_combos %>%
  mutate(
    result = pmap(list(eggmodel, tdm_model, tdm_calib, params), function(egg, tdm, calib, p) {
      message("Running: ", egg, " × ", tdm, " [", calib, "]")
      safe <- safely(fetch_tdm_eggcombo)
      res <- safe(
        spawn_date   = as.Date("2024-10-15"),
        redds        = 50,
        eggperredd   = 3000,
        temp_point   = "AWB:DailyAvg",
        eggmodel     = egg,
        tdm_model    = tdm,
        tdm_calib    = if (calib != "default") calib else NULL,
        pre_alpha    = p$pre_alpha, pre_beta  = p$pre_beta,
        post_alpha   = p$post_alpha, post_beta = p$post_beta,
        lin_alpha    = p$lin_alpha, lin_beta  = p$lin_beta,
        weib_TL      = 13, weib_TU = 15, weib_kL = 2.7, weib_kU = 3.2
      )
      return(res)
    }),
    cum_surv = map_dbl(result, ~ .x$result$cum_surv %||% NA_real_)
  ) %>% 
  select(eggmodel, tdm_model, tdm_calib, cum_surv)

# 5) Print and plot
print(results)

ggplot(results, aes(x = tdm_model, y = cum_surv, fill = eggmodel)) +
  geom_col(position = "dodge") +
  facet_wrap(~ tdm_calib, scales = "free_x") +
  labs(title = "Cumulative Survival by TDM & Egg Development Models",
       x = "TDM Model", y = "Cumulative Survival", fill = "Egg Model") +
  theme_minimal()
