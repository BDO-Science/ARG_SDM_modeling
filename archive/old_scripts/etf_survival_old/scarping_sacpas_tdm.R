fetch_tdm <- function(
    spawn_date,
    redds,
    eggperredd,
    temp_point,
    spawn_method       = c("text", "database"),
    db_source          = "American River spawned female carcasses",
    db_year            = format(Sys.Date(), "%Y"),
    carcasses_per_redd = 1,
    timing_offset      = 0,
    
    ## 1) New TDM args
    tdm_model = c("exp","lin","weibull","none","constant"),
    tdm_calib    = NULL,    # e.g. "Water Forum 2020" or "Zeug 2012, Cramer 2010"
    pre_alpha    = NULL,    # for exponential
    pre_beta     = NULL,
    post_alpha   = NULL,
    post_beta    = NULL,
    lin_alpha    = NULL,    # for linear
    lin_beta     = NULL,
    weib_TL      = NULL,    # for Weibull
    weib_TU      = NULL,
    weib_kL      = NULL,
    weib_kU      = NULL,
    
    densdeptype        = "none",
    rtnCSV             = "1",
    submit_button_text = "Run Emergence (to an output tab)"
) {
  spawn_method <- match.arg(spawn_method)
  tdm_model    <- match.arg(tdm_model)
  year_str     <- format(spawn_date, "%Y")
  doy          <- as.integer(format(spawn_date, "%j"))
  
  ## build the core body
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
    go             = submit_button_text
  )
  
  # spawn text vs database...
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
  
  ## 2) Wire in TDM selection
  ##    Note: replace 'mortalitytype' and other field names with the actual names from the form
  body$mortality <- tdm_model
  
  if (tdm_model == "exponential") {
    if (!is.null(tdm_calib)) body$exp_calibration <- tdm_calib
    if (!is.null(pre_alpha))  body$pre_hatch_alpha   <- as.character(pre_alpha)
    if (!is.null(pre_beta))   body$pre_hatch_beta    <- as.character(pre_beta)
    if (!is.null(post_alpha)) body$post_hatch_alpha  <- as.character(post_alpha)
    if (!is.null(post_beta))  body$post_hatch_beta   <- as.character(post_beta)
  }
  
  if (tdm_model == "linear") {
    if (!is.null(lin_alpha)) body$lin_alpha <- as.character(lin_alpha)
    if (!is.null(lin_beta))  body$lin_beta  <- as.character(lin_beta)
  }
  
  if (tdm_model == "weibull") {
    if (!is.null(weib_TL)) body$weib_eggs_TL <- as.character(weib_TL)
    if (!is.null(weib_TU)) body$weib_eggs_TU <- as.character(weib_TU)
    if (!is.null(weib_kL)) body$weib_eggs_kL <- as.character(weib_kL)
    if (!is.null(weib_kU)) body$weib_eggs_kU <- as.character(weib_kU)
    # …and similarly for the alevin inputs if needed
  }
  
  ## now POST and scrape as before
  res <- POST(
    "https://cbr.washington.edu/sac-bin/grow/emergecontrols.pl",
    encode = "form",
    body   = body
  )
  stop_for_status(res)
  
  doc      <- read_html(content(res, as = "text", encoding = "UTF-8"))
  csv_href <- doc %>% 
    html_node(xpath = "//a[contains(., 'download csv')]") %>% 
    html_attr("href")
  csv_url  <- url_absolute(csv_href, "https://cbr.washington.edu")
  
  txt    <- content(GET(csv_url), as = "text", encoding = "UTF-8")
  lines  <- strsplit(txt, "\r?\n")[[1]]
  hdr    <- grep("^SpawnDay,", lines)[1]
  df     <- read_csv(I(lines[hdr:length(lines)]), show_col_types = FALSE)
  
  list(details = df, cum_surv = df$Survival[1])
}

# 1) load libraries
library(httr)
library(rvest)
library(xml2)
library(readr)

# 2) source your fetch_tdm() definition
#    (or paste the function definition above here)
# source("path/to/your/fetch_tdm.R")

models <- c("exp","lin","weibull","none","constant")
for(m in models){
  cat("=== TDM model:", m, "===\n")
  out <- try(fetch_tdm(
    spawn_date   = as.Date("2024-10-15"),
    redds        = 50,
    eggperredd   = 3000,
    temp_point   = "AWB:DailyAvg",
    spawn_method = "text",
    tdm_model    = m,
    tdm_calib    = "Water Forum 2020",
    pre_alpha    = 3.40848,
    pre_beta     = 1.2112,
    post_alpha   = 1.01755,
    post_beta    = 1.24092,
    lin_alpha    = 0.0239,
    lin_beta     = 12.07,
    weib_TL      = 13,
    weib_TU      = 15,
    weib_kL      = 2.7,
    weib_kU      = 3.2
  ), silent=TRUE)
  
  if(inherits(out,"try-error")){
    cat("  ❌ error with", m, "\n\n")
  } else {
    cat("  ✅ cum_surv =", out$cum_surv, "\n\n")
  }}

