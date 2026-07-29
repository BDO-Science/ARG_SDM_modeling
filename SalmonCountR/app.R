library(shiny)
library(tidyverse)
library(DT)
library(scales)
library(shinyjs) # Using shinyjs for smoother slider updates
library(shinyWidgets) # For styled buttons
library(ggrepel)

# Add this line to fix the select() conflict
select <- dplyr::select

# The global.R file should load all pre-computed .rds files
source("global.R")

# Helper function to normalize a named vector of weights to sum to 1
normalize_weights <- function(weights) {
  total <- sum(weights, na.rm = TRUE)
  if (total > 0) {
    return(weights / total)
  }
  setNames(rep(1 / length(weights), length(weights)), names(weights))
}

# Helper for min-max normalization using actual data ranges
normalize_scores_chinook <- function(scores) {
  # Get the actual min and max from the current scores
  min_s <- min(scores, na.rm = TRUE)
  max_s <- max(scores, na.rm = TRUE)
  if (max_s == min_s) return(rep(0.5, length(scores)))
  (scores - min_s) / (max_s - min_s)
}

normalize_scores_steelhead <- function(scores) {
  # Get the actual min and max from the current scores
  min_s <- min(scores, na.rm = TRUE)
  max_s <- max(scores, na.rm = TRUE)
  if (max_s == min_s) return(rep(0.5, length(scores)))
  (scores - min_s) / (max_s - min_s)
}

normalize_scores_hydro <- function(scores) {
  # For hydropower, lower cost is better, so we invert
  min_s <- min(scores, na.rm = TRUE)
  max_s <- max(scores, na.rm = TRUE)
  if (max_s == min_s) return(rep(0.5, length(scores)))
  (max_s - scores) / (max_s - min_s)
}

get_scenario_alternatives <- function(scenario, hydro_year) {
  # Mapping: Alt 1-9 (2011), 10-18 (2014), 19-27 (2017), 28-36 (2020)
  # Within each hydro year: 1=NB, 2=PB1, 3=PB2, 4=PB2b, 5=PB2c, 6=PB3, 7=PB4, 8=PB5, 9=PB6
  scenario_map <- c("NB"=1, "PB1"=2, "PB2"=3, "PB2b"=4, "PB2c"=5, "PB3"=6, "PB4"=7, "PB5"=8, "PB6"=9)
  hydro_map <- c("2011"=0, "2014"=9, "2017"=18, "2020"=27)
  
  if (hydro_year == "all") {
    base_idx <- scenario_map[scenario]
    return(c(base_idx, base_idx+9, base_idx+18, base_idx+27))
  } else {
    return(hydro_map[hydro_year] + scenario_map[scenario])
  }
}

ui <- navbarPage("Lower American River Power Bypass Decision Support",
                 useShinyjs(), 
                 
                 # ---- About Tab ----
                 tabPanel("About",
                          fluidRow(
                            column(12,
                                   h3("Author & Contact Information"),
                                   tags$ul(
                                     tags$li(strong("Author:"), "Alexander Vaisvil"),
                                     tags$li(strong("Email:"), tags$a(href = "mailto:avaisvil@usbr.gov", "avaisvil@usbr.gov")),
                                     tags$li(strong("Institution/Organization:"), "U.S. Bureau of Reclamation"),
                                     tags$li(strong("GitHub Repository:"), tags$a(href = "https://github.com/BDO-Science/ARG_SDM_modeling", "https://github.com/BDO-Science/ARG_SDM_modeling")),
                                     tags$li(strong("Date Last Updated:"), format(file.info("app.R")$mtime, "%B %d, %Y"))
                                   ),
                                   hr(),
                                   
                                   h3("Application Overview"),
                                   p("This interactive application simulates fall-run Chinook salmon population dynamics on the American River under different Folsom Dam power bypass flow and temperature management alternatives. The model integrates temperature-dependent mortality, spawn timing, and comprehensive life-cycle processes to project spawner abundance through 2125."),
                                   
                                   h4("Management Structure:"),
                                   tags$ul(
                                     tags$li(strong("9 Management Alternatives:"), "No Bypass (NB) and 8 Power Bypass configurations (PB1-PB6, including PB2b and PB2c variants) with varying flow rates and timing"),
                                     tags$li(strong("4 Climate Years:"), "2011 (Cool), 2014 (Warm), 2017 (Warm), 2020 (Cool)"),
                                     tags$li(strong("36 Pre-computed Alternatives:"), "Each alternative modeled under all 4 climate year conditions, allowing dynamic weighting"),
                                     tags$li(strong("Temperature Data:"), "Sept 22 - Nov 30 forecast window from SDM Power Bypass modeling, combined with 14-year USGS gauge climatology (2011-2025) for full annual cycles"),
                                     tags$li(strong("Simulation Period:"), "2025-2125 with user-adjustable weighting of climatological conditions and TDM models")
                                   ),
                                   
                                   h4("Power Bypass Alternative Specifications:"),
                                   tags$div(style = "margin-left: 20px;",
                                            tags$table(class = "table table-striped table-condensed",
                                                       tags$thead(
                                                         tags$tr(
                                                           tags$th("Alternative"),
                                                           tags$th("Bypass (AF)"),
                                                           tags$th("Bypass (MWh)"),
                                                           tags$th("Loss ($)"),
                                                           tags$th("Increase (mTCO2)"),
                                                           tags$th("Description")
                                                         )
                                                       ),
                                                       tags$tbody(
                                                         tags$tr(
                                                           tags$td(strong("NB")),
                                                           tags$td("0"),
                                                           tags$td("0"),
                                                           tags$td("$0"),
                                                           tags$td("0"),
                                                           tags$td("No bypass - baseline operations")
                                                         ),
                                                         tags$tr(
                                                           tags$td(strong("PB1")),
                                                           tags$td("10,163"),
                                                           tags$td("2,424"),
                                                           tags$td("$111,422"),
                                                           tags$td("1,149"),
                                                           tags$td("125 cfs starting Oct 15, 250 cfs on Oct 28, 125 cfs on Nov 7, end bypass on Nov 14")
                                                         ),
                                                         tags$tr(
                                                           tags$td(strong("PB2")),
                                                           tags$td("32,224"),
                                                           tags$td("7,674"),
                                                           tags$td("$376,671"),
                                                           tags$td("3,650"),
                                                           tags$td("250 cfs starting Oct 15, 500 cfs on Oct 28, 250 cfs on Nov 14, end bypass on Nov 30")
                                                         ),
                                                         tags$tr(
                                                           tags$td(strong("PB2b")),
                                                           tags$td("40,156"),
                                                           tags$td("9,558"),
                                                           tags$td("$470,090"),
                                                           tags$td("4,522"),
                                                           tags$td("250 cfs starting Oct 15, 500 cfs on Oct 28, end bypass on Nov 30")
                                                         ),
                                                         tags$tr(
                                                           tags$td(strong("PB2c")),
                                                           tags$td("37,181"),
                                                           tags$td("8,846"),
                                                           tags$td("$433,215"),
                                                           tags$td("4,195"),
                                                           tags$td("250 cfs starting Oct 21, 500 cfs on Oct 28, end bypass on Nov 30")
                                                         ),
                                                         tags$tr(
                                                           tags$td(strong("PB3")),
                                                           tags$td("17,351"),
                                                           tags$td("4,135"),
                                                           tags$td("$201,552"),
                                                           tags$td("1,932"),
                                                           tags$td("250 cfs starting Oct 21, 500 cfs on Oct 28, 250 cfs on Nov 7, end bypass on Nov 14")
                                                         ),
                                                         tags$tr(
                                                           tags$td(strong("PB4")),
                                                           tags$td("20,822"),
                                                           tags$td("4,959"),
                                                           tags$td("$241,590"),
                                                           tags$td("2,350"),
                                                           tags$td("250 cfs starting Oct 21, 500 cfs on Oct 28, 250 cfs on Nov 7, end bypass on Nov 21")
                                                         ),
                                                         tags$tr(
                                                           tags$td(strong("PB5")),
                                                           tags$td("17,351"),
                                                           tags$td("4,130"),
                                                           tags$td("$199,382"),
                                                           tags$td("1,974"),
                                                           tags$td("500 cfs bypass starting Oct 28, reduce to 250 on Nov 7, end bypass on Nov 21")
                                                         ),
                                                         tags$tr(
                                                           tags$td(strong("PB6")),
                                                           tags$td("30,141"),
                                                           tags$td("7,100"),
                                                           tags$td("$348,806"),
                                                           tags$td("3,321"),
                                                           tags$td("100 cfs Oct 1, 200 cfs Oct 8, 300 cfs Oct 15, 400 cfs Oct 22, 500 cfs Nov 1, ending Nov 14")
                                                         )
                                                       )
                                            )
                                   ),
                                   
                                   h3("Model Components"),
                                   
                                   h4("1. Temperature Data Integration"),
                                   p("The model combines multiple temperature data sources:"),
                                   tags$ul(
                                     tags$li(strong("USGS Gauge Data:"), "Daily observations from stations 11446980 (Watt Ave) and 11446500 (Hazel Ave) for 2011-2025"),
                                     tags$li(strong("Power Bypass Forecasts:"), "Modeled temperatures for Oct 18 - Dec 31 under different bypass scenarios"),
                                     tags$li(strong("14-Year DOY Climatology:"), "Day-of-year averages from observed data used to fill gaps outside forecast window"),
                                     tags$li(strong("Temperature Floors:"), "7°C minimum for Hazel Ave, 8°C for Watt Ave to reflect thermal refugia")
                                   ),
                                   
                                   h4("2. Spawn Timing Model (CLM)"),
                                   p("A cumulative link model predicts spawning probability across 10-day temporal bins based on October and November water temperatures:"),
                                   tags$ul(
                                     tags$li("Spawning season divided into 12+ periods from October through January"),
                                     tags$li("Temperature coefficients: β_Oct and β_Nov control timing shifts with temperature"),
                                     tags$li("Threshold parameters (ζ) define cumulative probability boundaries between periods"),
                                     tags$li("Model formula: P(spawn in period j) = logit⁻¹(ζⱼ - β'X) - logit⁻¹(ζⱼ₋₁ - β'X)"),
                                     tags$li("Calibrated on 2011-2024 American River carcass survey data")
                                   ),
                                   
                                   h4("3. Temperature-Dependent Mortality (TDM)"),
                                   p("Three mortality models calculate egg-to-fry survival from spawn to emergence (958 ATU):"),
                                   
                                   tags$div(style = "margin-left: 20px;",
                                            tags$h5("Exponential Models:"),
                                            tags$ul(
                                              tags$li(HTML("<strong>Water Forum 2020:</strong> Recent calibration with stage-specific parameters")),
                                              tags$li(HTML("&nbsp;&nbsp;• Egg stage (0-400 ATU): α = 3.408×10⁻¹¹, β = 1.211")),
                                              tags$li(HTML("&nbsp;&nbsp;• Alevin stage (400-958 ATU): α = 1.018×10⁻¹⁰, β = 1.241")),
                                              tags$li(HTML("<strong>SALMOD 2006:</strong> Historical parameters from USFWS model")),
                                              tags$li(HTML("&nbsp;&nbsp;• Egg: α = 1.475×10⁻¹¹, β = 1.392")),
                                              tags$li(HTML("&nbsp;&nbsp;• Alevin: α = 2.521×10⁻¹², β = 1.461"))
                                            ),
                                            
                                            tags$h5("Linear Threshold Model:"),
                                            tags$ul(
                                              tags$li(HTML("<strong>Martin et al. 2017:</strong> h(T) = 0.026 × max(T - 12.14°C, 0)")),
                                              tags$li("Zero mortality below 12.14°C threshold, linear increase above")
                                            )
                                   ),
                                   
                                   h4("4. Adult Pre-spawn Survival"),
                                   p("Thermal stress accumulates as degree-days from October 1 through spawning:"),
                                   tags$ul(
                                     tags$li(HTML("Survival function: S<sub>pre</sub> = 1 / (1 + exp(-(3.0 - 0.00067 × DD)))")),
                                     tags$li("Degree-days: DD = Σ max(T - 0°C, 0) from Oct 1 to spawn date"),
                                     tags$li("Higher accumulated thermal units reduce pre-spawn survival"),
                                     tags$li("Based on Columbia River Chinook studies (Colvin et al. 2018)")
                                   ),
                                   
                                   h4("5. Population Dynamics"),
                                   p("Age-structured life-cycle model with density dependence:"),
                                   tags$ul(
                                     tags$li(HTML("<strong>Egg Production:</strong> Eggs = Spawners × 0.5 (female fraction) × S<sub>pre</sub> × 5,522 (fecundity)")),
                                     tags$li(HTML("<strong>Density Dependence:</strong> Beverton-Holt at fry stage: dd = 0.347 / (1 + Redds / K)")),
                                     tags$li(HTML("<strong>Carrying Capacity:</strong> K = 33,185 redds at the 1,000 cfs baseline (flow-dependent via slider)")),
                                     tags$li(HTML("<strong>Fry Production:</strong> Fry = Eggs × S<sub>TDM</sub> × dd")),
                                     tags$li(HTML("<strong>Smolt Production:</strong> Smolts = Fry × rear_surv (calibrated at 0.5428)")),
                                     tags$li(HTML("<strong>Non-American River Survival:</strong> SAR = 0.00269 (calibrated value representing 0.27% survival from the Sacramento River to returning adults)")),
                                     tags$li(HTML("<strong>Age Structure:</strong> 82.9% age-3, 16.9% age-4, 0.2% age-5 returns (CWT data)")),
                                     tags$li(HTML("<strong>Returns:</strong> Spawners<sub>t+a</sub> = Smolts<sub>t</sub> × SAR × P(age=a)"))
                                   ),
                                   
                                   h3("Model Parameters"),
                                   
                                   h4("Calibrated Life-Cycle Parameters:"),
                                   p("The model uses pre-specified biological parameters rather than statistical calibration:"),
                                   tags$ul(
                                     tags$li(HTML("<strong>Smolt-to-Adult Return (SAR):</strong> 0.00269 (0.27%) - calibrated value representing ocean survival")),
                                     tags$li(HTML("<strong>Rearing Survival:</strong> 0.5428 (54.28%) - calibrated freshwater survival from fry to smolt")),
                                     tags$li(HTML("<strong>Initial Population:</strong> Years 2011-2013 seeded from CDFW GrandTab observed escapement")),
                                     tags$li(HTML("<strong>Forecast Starting Point:</strong> 2022-2024 observed escapement used as initial conditions for 2025+ projections"))
                                   ),
                                   
                                   p("These calibrated parameters are applied uniformly across all TDM variants, with forecasts driven primarily by temperature-dependent egg-to-fry survival differences between scenarios."),
                                   
                                   h4("Data Sources:"),
                                   tags$ul(
                                     tags$li(strong("Escapement:"), "CDFW GrandTab American River fall-run Chinook (2011-2024)"),
                                     tags$li(strong("Carcass Surveys:"), "American River sections NB, W, 1a, 1b, 2, 3 (2011-2024)"),
                                     tags$li(strong("Temperature:"), "USGS stations 11446980, 11446500 (15-minute intervals aggregated to daily)"),
                                     tags$li(strong("Age Structure:"), "Coded wire tag recoveries from Central Valley database"),
                                     tags$li(strong("Power Bypass Modeling:"), "SDM temperature projections for 9 operational scenarios")
                                   ),
                                   
                                   h3("Key Equations Reference"),
                                   tags$div(style = "background-color: #f5f5f5; padding: 15px; border-radius: 5px;",
                                            tags$ol(
                                              tags$li(HTML("<strong>ATU Accumulation:</strong> ATU = Σ max(T<sub>daily</sub>, 0)"), 
                                                      br(), tags$small(em("Accumulated thermal units drive developmental timing"))),
                                              
                                              tags$li(HTML("<strong>Egg Development:</strong> Days to hatch = 400 / T̄"), 
                                                      br(), tags$small(em("Temperature-dependent development rate to reach alevin stage"))),
                                              
                                              tags$li(HTML("<strong>Alevin Emergence:</strong> Days to emerge = 558 / T̄"), 
                                                      br(), tags$small(em("Additional thermal units from hatch to fry emergence"))),
                                              
                                              tags$li(HTML("<strong>Exponential TDM:</strong> S = exp(-Σ α·exp(β·T<sub>i</sub>))"), 
                                                      br(), tags$small(em("Daily hazard accumulated over incubation period"))),
                                              
                                              tags$li(HTML("<strong>Linear TDM:</strong> S = exp(-α·Σ max(T<sub>i</sub> - 12.14, 0))"), 
                                                      br(), tags$small(em("Threshold model with mortality above critical temperature"))),
                                              
                                              tags$li(HTML("<strong>Pre-spawn Survival:</strong> S<sub>pre</sub> = logit⁻¹(3.0 - 0.00067·DD)"), 
                                                      br(), tags$small(em("Logistic decline with accumulated degree-days"))),
                                              
                                              tags$li(HTML("<strong>Beverton-Holt DD:</strong> S<sub>dd</sub> = S₀ / (1 + Redds/K)"), 
                                                      br(), tags$small(em("Density-dependent survival at maximum S₀ = 0.347"))),
                                              
                                              tags$li(HTML("<strong>Effective Survival:</strong> S<sub>eff</sub> = S<sub>TDM</sub> × S<sub>dd</sub>"), 
                                                      br(), tags$small(em("Combined temperature and density effects on egg-to-fry"))),
                                              
                                              tags$li(HTML("<strong>Cohort Returns:</strong> N<sub>t+a</sub> = Σ Smolts<sub>t</sub> × SAR<sub>t</sub> × P(age=a)"), 
                                                      br(), tags$small(em("Age-structured returns with stochastic ocean survival option")))
                                            )
                                   ),
                                   
                                   h3("Application Features"),
                                   
                                   h4("Interactive Components:"),
                                   tags$ul(
                                     tags$li(strong("Temperature Explorer:"), "Visualize and compare weighted temperature patterns across multiple alternatives and climatological years"),
                                     tags$li(strong("Alternative Comparison:"), "Side-by-side evaluation of multiple alternatives with customizable climatology and TDM weighting"),
                                     tags$li(strong("Swing Weighting:"), "Interactive preference elicitation tool to determine objective importance through hypothetical alternative rankings"),
                                     tags$li(strong("Decision Support:"), "Multi-objective analysis with consequence tables, trade-off plots, and weighted performance scores")
                                   ),
                                   
                                   h4("Weighting Options:"),
                                   tags$ul(
                                     tags$li(strong("Climatology Weights:"), "Adjust relative importance of 4 climate year types (sum to 1.0)"),
                                     tags$li(strong("TDM Model Weights:"), "Combine 3 mortality models with custom weights (sum to 1.0)"),
                                     tags$li(strong("Objective Weights:"), "Three methods: (1) Equal weights (33.3% each), (2) Manual slider adjustment, or (3) Derive from Swing Weighting tab"),
                                     tags$li(strong("Default Weights:"), "TDM: 51% Water Forum, 24% SALMOD, 25% Martin; Climatology: 25% each; Objectives: 40% Chinook, 30% Steelhead, 30% Hydropower")
                                   ),
                                   
                                   h3("Technical Implementation"),
                                   
                                   h4("Computational Approach:"),
                                   tags$ul(
                                     tags$li("Parallel processing using furrr package for TDM calculations"),
                                     tags$li("Memoization of survival calculations to avoid redundant computation"),
                                     tags$li("Data.table for efficient large-scale data manipulation"),
                                     tags$li("Byte-compilation of performance-critical functions"),
                                     tags$li("Pre-computed lookup tables for temperature-date mapping")
                                   ),
                                   
                                   h4("File Structure:"),
                                   tags$ul(
                                     tags$li(strong("temperature_data.R:"), "Processes USGS gauge data and SDM forecasts"),
                                     tags$li(strong("precompute.R:"), "Calibrates models and generates population forecasts"),
                                     tags$li(strong("functions.R:"), "Core TDM and life-cycle simulation functions"),
                                     tags$li(strong("global.R:"), "Loads pre-computed RDS files for app startup"),
                                     tags$li(strong("app.R:"), "Shiny interface and reactive logic")
                                   ),
                                   
                                   h3("References"),
                                   tags$ul(
                                     tags$li("Bartholow, J.M. & Heasley, J. (2006). Evaluation of Shasta Dam Alternatives Using a Salmon Production Model. USGS Open-File Report 2004-1351."),
                                     tags$li("Bratovich, P., Neal, M., Ransom, A., et al. (2020). Chinook Salmon Early Lifestage Survival and Folsom Dam Power Bypass Considerations. Water Forum Technical Memorandum."),
                                     tags$li("Colvin, R., Falke, J.A., Henson, S. (2018). Identifying optimal water temperature and flow regimes for anadromous fish. River Research and Applications 34(6):621-632."),
                                     tags$li("Martin, B.T., Pike, A., John, S.N., Hamda, N., Roberts, J., Lindley, S.T., Danner, E.M. (2017). Phenomenological vs. biophysical models of thermal stress in aquatic eggs. Ecology Letters 20:50-59."),
                                     tags$li("USFWS (2006). SALMOD: Salmon Population Model Version 3.0. Sacramento Fish and Wildlife Office."),
                                     tags$li("CDFW GrandTab (2024). California Central Valley Chinook Salmon Escapement Database. Available at: https://wildlife.ca.gov/Conservation/Fishes/Chinook-Salmon/Anadromous-Assessment")
                                   )
                            )
                          )
                 ),
                 
                 # ---- Temperature Explorer Tab (FIXED) ----
                 tabPanel("Temperature Explorer",
                          sidebarLayout(
                            sidebarPanel(
                              h4("Climatology Weights"),  # Changed
                              sliderInput("temp_w_2011", "2011 (Cool)", value = 0.25, min = 0, max = 1, step = 0.01),
                              sliderInput("temp_w_2014", "2014 (Warm)", value = 0.25, min = 0, max = 1, step = 0.01),
                              sliderInput("temp_w_2017", "2017 (Warm)", value = 0.25, min = 0, max = 1, step = 0.01),
                              sliderInput("temp_w_2020", "2020 (Cool)", value = 0.25, min = 0, max = 1, step = 0.01),
                              hr(),
                              # In the UI, replace the selectInput with:
                              checkboxGroupInput("temp_alternatives", "Alternatives to Compare:",
                                                 choices = c("No Bypass"="NB", "Power Bypass 1"="PB1", 
                                                             "Power Bypass 2"="PB2", "Power Bypass 2b"="PB2b", 
                                                             "Power Bypass 2c"="PB2c", "Power Bypass 3"="PB3",
                                                             "Power Bypass 4"="PB4", "Power Bypass 5"="PB5", 
                                                             "Power Bypass 6"="PB6"),
                                                 selected = c("NB", "PB1", "PB2", "PB2b", "PB2c", "PB3", "PB4", "PB5", "PB6")),  # Default to comparing NB and PB1
                              radioButtons("temp_site", "Site:", choices = c("Ave Watt"="AveWatt", "Ave Hazel"="AveHazel")),
                              radioButtons("temp_period", "Time Period:",
                                           choices = c("Oct-Dec 2025"="oct_dec", "Full Year 2025"="full")),
                              width = 3
                            ),
                            mainPanel(
                              plotOutput("temp_plot", height = "500px"),
                              hr(),
                              h4("Temperature Statistics"),
                              tableOutput("temp_stats"),
                              width = 9
                            )
                          )
                 ),
                 
                 
                 # Compare Alternatives Tab → Compare Alternatives Tab
                 tabPanel("Compare Alternatives",  # Changed from "Compare Alternatives"
                          sidebarLayout(
                            sidebarPanel(
                              h4("Alternatives to Compare"),  # Changed from "Alternatives to Compare"
                              checkboxGroupInput("cmp_scenarios", "Select:",  # Keep the ID for compatibility
                                                 choices = c("No Bypass"="NB", "Power Bypass 1"="PB1", "Power Bypass 2"="PB2",
                                                             "Power Bypass 2b"="PB2b", "Power Bypass 2c"="PB2c",
                                                             "Power Bypass 3"="PB3", "Power Bypass 4"="PB4",
                                                             "Power Bypass 5"="PB5", "Power Bypass 6"="PB6"),
                                                 selected = c("NB", "PB1", "PB2", "PB2b", "PB2c", "PB3", "PB4", "PB5", "PB6")),
                              hr(),
                              h4("Climatology Weights"),  # Changed from "Hydrology Weights"
                              sliderInput("cmp_w_2011", "2011 (Cool)", value = 0.25, min = 0, max = 1, step = 0.01),  # Changed from "Wet"
                              sliderInput("cmp_w_2014", "2014 (Warm)", value = 0.25, min = 0, max = 1, step = 0.01),  # Changed from "Critical"
                              sliderInput("cmp_w_2017", "2017 (Warm)", value = 0.25, min = 0, max = 1, step = 0.01),  # Changed from "Wet"
                              sliderInput("cmp_w_2020", "2020 (Cool)", value = 0.25, min = 0, max = 1, step = 0.01),  # Changed from "Dry"
                              hr(),
                              h4("TDM Weights"),
                              sliderInput("cmp_tdm_wf", "Water Forum", value = 0.51, min = 0, max = 1, step = 0.01),
                              sliderInput("cmp_tdm_sm", "SALMOD", value = 0.24, min = 0, max = 1, step = 0.01),
                              sliderInput("cmp_tdm_martin", "Martin", value = 0.25, min = 0, max = 1, step = 0.01),
                              hr(),
                              # ★★★★★★★★★★★
                              h4("Flow and Capacity"),
                              sliderInput("cmp_flow", "Set Downstream Flow (cfs)",
                                          min = 500, max = 5000, value = 1000, step = 100),
                              hr(),
                              # ★★★★★★★★★★
                              
                              sliderInput("cmp_years", "Forecast Years", value = 100, min = 10, max = 100),
                              actionButton("run_cmp", "Run Comparison"),
                              width = 3
                            ),
                            mainPanel(
                              tabsetPanel(
                                tabPanel("Time Series", plotOutput("cmp_ts_plot")),
                                tabPanel("Boxplot",
                                         sliderInput("last_n", "Last N Years:", min = 5, max = 100, value = 20, step = 1),
                                         plotOutput("cmp_box_plot"),
                                         hr(),
                                         h4("Boxplot Summary"),
                                         tableOutput("cmp_boxplot_stats")),
                                tabPanel("Summary Table", DTOutput("cmp_summary"))
                              ),
                              width = 9
                            )
                          )
                 ),
                 
                 # ---- Swing Weighting Tab (NEW TOP-LEVEL TAB) ----
                 tabPanel("Swing Weighting",
                          fluidRow(
                            column(12,
                                   h3("Swing Weighting for Objective Importance"),
                                   p("Use this tool to determine objective weights by comparing hypothetical alternatives. Once you've calculated weights, you can apply them in the Decision Support tab."),
                                   hr(),
                                   
                                   wellPanel(
                                     h5(icon("info-circle"), "How Swing Weighting Works:"),
                                     tags$ol(
                                       tags$li("Review the three hypothetical alternatives below. Each performs BEST on one objective and WORST on all others."),
                                       tags$li("The 'Worst Alternative' scores worst on all objectives (rank 4, score 0)."),
                                       tags$li(strong("Rank"), "alternatives 1-3 based on which improvement matters most to you (1 = most valuable improvement)."),
                                       tags$li(strong("Score the alternatives:"), 
                                               tags$ul(
                                                 tags$li(strong("Rank 1 gets 100 points"), " (this is the reference)"),
                                                 tags$li("Rank 2: Score between 0-100 based on how valuable that improvement is compared to rank 1"),
                                                 tags$li("Rank 3: Score between 0-100 based on how valuable that improvement is compared to rank 1"),
                                                 tags$li("Example: If rank 2 is half as important as rank 1, give it 50 points")
                                               )),
                                       tags$li("Weights are calculated automatically: weight = score / sum(all scores)."),
                                       tags$li("Click 'Apply to Decision Support' to use these weights in your analysis.")
                                     )
                                   ),
                                   
                                   hr(),
                                   h4("Objective Ranges"),
                                   tableOutput("swing_ranges_table"),
                                   
                                   hr(),
                                   h4("Hypothetical Alternatives to Rank"),
                                   p(em("Each alternative performs BEST on one objective and worst on all others.")),
                                   tableOutput("swing_alternatives_table"),
                                   
                                   hr(),
                                   h4("Your Rankings and Scores"),
                                   fluidRow(
                                     column(4,
                                            wellPanel(
                                              h5("Alternative 1: Best Chinook"),
                                              numericInput("rank_chinook", "Rank (1-3):", value = 1, min = 1, max = 3, step = 1),
                                              numericInput("score_chinook", "Score (0-100):", value = 80, min = 0, max = 100, step = 1)
                                            )
                                     ),
                                     column(4,
                                            wellPanel(
                                              h5("Alternative 2: Best Steelhead"),
                                              numericInput("rank_steelhead", "Rank (1-3):", value = 3, min = 1, max = 3, step = 1),
                                              numericInput("score_steelhead", "Score (0-100):", value = 20, min = 0, max = 100, step = 1)
                                            )
                                     ),
                                     column(4,
                                            wellPanel(
                                              h5("Alternative 3: Best Hydropower"),
                                              numericInput("rank_hydropower", "Rank (1-3):", value = 2, min = 1, max = 3, step = 1),
                                              numericInput("score_hydropower", "Score (0-100):", value = 100, min = 0, max = 100, step = 1)
                                            )
                                     )
                                   ),
                                   
                                   uiOutput("swing_validation"),
                                   
                                   hr(),
                                   h4("Calculated Weights"),
                                   tableOutput("swing_weights_table"),
                                   plotOutput("swing_weights_plot", height = "300px"),
                                   
                                   hr(),
                                   actionButton("apply_swing_weights", "Apply to Decision Support Tab", 
                                                class = "btn-primary btn-lg", icon = icon("arrow-right"))
                            )
                          )
                 ),
                 
                 tabPanel("Decision Support",
                          sidebarLayout(
                            sidebarPanel(
                              h4("Objective Weighting"),
                              radioButtons("weight_method", "Choose Weighting Method:",
                                           choices = c("Equal Weights" = "equal",
                                                       "Manual Weights" = "manual"),
                                           selected = "manual"),
                              hr(),
                              conditionalPanel(
                                condition = "input.weight_method == 'manual'",
                                p("Adjust sliders to reflect importance. Weights sum to 1."),
                                sliderInput("w_chinook", "Fall-run Chinook", min = 0, max = 1, value = 0.4, step = 0.01),
                                sliderInput("w_steelhead", "Steelhead", min = 0, max = 1, value = 0.1, step = 0.01),
                                sliderInput("w_hydro", "Hydropower", min = 0, max = 1, value = 0.5, step = 0.01),
                                hr(),
                                p(em("Tip: Use the Swing Weighting tab to determine these weights systematically."))
                              ),
                              width = 3
                            ),
                            mainPanel(
                              tabsetPanel(
                                tabPanel("Overall Performance",
                                         h4("Overall Weighted Scores"),
                                         plotOutput("overall_scores_plot"),
                                         hr(),
                                         h4("Score Contribution by Objective"),
                                         plotOutput("plot_score_contribution"),
                                         hr(),
                                         h4("Current Objective Weights"),
                                         tableOutput("weights_table")),
                                tabPanel("Consequence Table",
                                         h4("Consequence Table"),
                                         p("Raw and normalized (0-1) scores for each objective by scenario."),
                                         tableOutput("performance_matrix")),
                                tabPanel("Trade-offs",
                                         h4("Multi-Objective Performance"),
                                         plotOutput("performance_plot"),
                                         hr(),
                                         h4("Trade-off Analysis (Chinook vs. Hydropower)"),
                                         plotOutput("tradeoff_plot"))
                              ),
                              width = 9
                            )
                          )
                 ),

                 # ---- Add a Year Tab -------------------------------------------
                 # Lets a non-R user drop in a new temperature modelling
                 # deliverable and get results, without touching precompute.R.
                 tabPanel("Add a Year",
                          sidebarLayout(
                            sidebarPanel(
                              width = 3,
                              h4("1. Upload"),
                              helpText(
                                "Upload the temperature modelling spreadsheet exactly as it ",
                                "arrives from the modelling team. No reformatting needed."
                              ),
                              fileInput("newyear_file", NULL,
                                        accept = c(".xlsx", ".xls"),
                                        buttonLabel = "Choose file...",
                                        placeholder = "No file selected"),
                              tags$hr(),
                              h4("2. Settings"),
                              sliderInput("newyear_flow", "Reference flow (cfs)",
                                          min = 500, max = 3000, value = 1000, step = 100),
                              helpText("Sets spawning carrying capacity. 1,000 cfs is the value used throughout the analysis."),
                              tags$hr(),
                              h4("3. Run"),
                              actionButton("newyear_run", "Run analysis",
                                           class = "btn-primary", width = "100%",
                                           icon = icon("play")),
                              tags$br(), tags$br(),
                              uiOutput("newyear_status"),
                              tags$hr(),
                              h4("4. Download"),
                              downloadButton("newyear_dl_scores", "Scores (CSV)",
                                             style = "width:100%; margin-bottom:6px;"),
                              downloadButton("newyear_dl_detail", "Full detail (CSV)",
                                             style = "width:100%;")
                            ),
                            mainPanel(
                              width = 9,
                              tabsetPanel(
                                tabPanel(
                                  "File check",
                                  br(),
                                  p(strong("Every upload is checked before anything is computed."),
                                    " Errors stop the run; warnings do not."),
                                  DTOutput("newyear_report"),
                                  br(),
                                  uiOutput("newyear_help")
                                ),
                                tabPanel(
                                  "Results",
                                  br(),
                                  uiOutput("newyear_headline"),
                                  br(),
                                  h4("Composite scores"),
                                  plotOutput("newyear_score_plot", height = "420px"),
                                  br(),
                                  h4("Consequence table"),
                                  DTOutput("newyear_scores")
                                ),
                                tabPanel(
                                  "Against the baseline",
                                  br(),
                                  div(
                                    style = "background:#eef4fa; border-left:4px solid #31688E; padding:10px 14px; margin-bottom:12px;",
                                    strong("Like-for-like comparison. "),
                                    "The baseline below is the published four meteorological years ",
                                    "re-run through this same engine, not the numbers printed in the ",
                                    "manuscript. Comparing an upload against published figures would ",
                                    "mix two different calculations."
                                  ),
                                  DTOutput("newyear_vs_base"),
                                  br(),
                                  plotOutput("newyear_vs_base_plot", height = "400px")
                                ),
                                tabPanel(
                                  "Temperatures",
                                  br(),
                                  p("Daily temperature at Hazel Avenue, each uploaded scenario minus no-bypass."),
                                  plotOutput("newyear_temp_plot", height = "480px")
                                ),
                                tabPanel(
                                  "What this does",
                                  br(),
                                  uiOutput("newyear_doc")
                                )
                              )
                            )
                          )
                 )
)

server <- function(input, output, session) {

  # Hard-coded Hydropower Scores
  hardcoded_hydro_scores <- c(
    "NB"  = 0, "PB1" = 111422, "PB2" = 370826, "PB2b" = 470090, "PB2c" = 433215,
    "PB3" = 201552, "PB4" = 241590, "PB5" = 199382, "PB6" = 348806
  )

  # ===========================================================================
  # "Add a Year" tab
  # ===========================================================================
  # Session-scoped: everything lives in this reactiveValues and dies with the
  # session. Nothing is written to disk and no global is modified, which matters
  # because shinyapps.io shares one R process across sessions.

  ny <- reactiveValues(res = NULL, base = NULL, error = NULL, ran = FALSE)

  ny_ready <- reactive({
    isTRUE(exists("scenario_engine_loaded")) && isTRUE(exists("spawn_timing_model")) &&
      !is.null(spawn_timing_model)
  })

  # Baseline: the published deliverable through the SAME engine, computed once
  # per session and cached, so uploads are compared like for like.
  ny_baseline <- function() {
    if (!is.null(ny$base)) return(ny$base)
    p <- here::here("data_raw", "SDM Power Bypass Temperature Modeling Results.xlsx")
    if (!file.exists(p)) return(NULL)
    out <- try(run_scenario_deliverable(
      path = p, env_ext_list = env_ext_list, spawn_model = spawn_timing_model,
      base_P_list = base_P_list, S_seed_fore_list = S_seed_fore_list,
      sim_years = sim_years, hydro_loss = hardcoded_hydro_scores,
      flow_cfs = input$newyear_flow, K_fn = get_K_spawners), silent = TRUE)
    if (inherits(out, "try-error") || !isTRUE(out$ok)) return(NULL)
    ny$base <- out
    out
  }

  # -- Validate as soon as a file arrives, before any Run click ---------------
  ny_parsed <- reactive({
    req(input$newyear_file)
    read_deliverable(input$newyear_file$datapath)
  })

  ny_report <- reactive({
    p <- ny_parsed()
    validate_deliverable(p)
  })

  output$newyear_report <- renderDT({
    rep <- ny_report()
    icon_for <- function(s) switch(s,
      ok      = '<span style="color:#1a7f37;font-weight:bold;">&#10003; OK</span>',
      warning = '<span style="color:#9a6700;font-weight:bold;">&#9888; Check</span>',
      error   = '<span style="color:#b42318;font-weight:bold;">&#10007; Problem</span>')
    df <- data.frame(
      Status = vapply(rep$status, icon_for, character(1)),
      Check  = rep$check,
      Detail = rep$message,
      stringsAsFactors = FALSE
    )
    datatable(df, escape = FALSE, rownames = FALSE,
              options = list(dom = "t", paging = FALSE, ordering = FALSE,
                             columnDefs = list(list(width = "110px", targets = 0),
                                               list(width = "170px", targets = 1))))
  })

  output$newyear_help <- renderUI({
    if (is.null(input$newyear_file)) {
      return(div(style = "color:#555;",
                 p(strong("No file uploaded yet.")),
                 p("The spreadsheet should have one sheet per meteorological year, ",
                   "named as a four-digit year. Scenario names go on the first row, ",
                   "and ", code("AveWatt"), " / ", code("AveHazel"),
                   " column headers on the second.")))
    }
    rep <- ny_report()
    if (any(rep$status == "error")) {
      div(style = "background:#fdf2f0; border-left:4px solid #b42318; padding:10px 14px;",
          strong("This file cannot be used yet. "),
          "Fix the items marked Problem above and upload again. ",
          "Items marked Check are safe to ignore if you expected them.")
    } else {
      div(style = "background:#eefaf0; border-left:4px solid #1a7f37; padding:10px 14px;",
          strong("File looks good. "), "Click ", strong("Run analysis"), " to compute results.")
    }
  })

  # -- Run --------------------------------------------------------------------
  observeEvent(input$newyear_run, {
    req(input$newyear_file)
    ny$error <- NULL; ny$res <- NULL; ny$ran <- FALSE

    if (!ny_ready()) {
      ny$error <- paste("The scenario engine is not available in this deployment.",
                        "app_data/spawn_timing_model.rds is missing -- run",
                        "analysis/build_spawn_timing_model.R and redeploy.")
      return(invisible(NULL))
    }
    if (any(ny_report()$status == "error")) {
      ny$error <- "The uploaded file has problems that must be fixed first. See the File check tab."
      return(invisible(NULL))
    }

    withProgress(message = "Running", value = 0, {
      out <- try(run_scenario_deliverable(
        path = input$newyear_file$datapath,
        env_ext_list = env_ext_list, spawn_model = spawn_timing_model,
        base_P_list = base_P_list, S_seed_fore_list = S_seed_fore_list,
        sim_years = sim_years, hydro_loss = hardcoded_hydro_scores,
        flow_cfs = input$newyear_flow, K_fn = get_K_spawners,
        progress = function(f, m) setProgress(value = f, detail = m)
      ), silent = TRUE)

      if (inherits(out, "try-error")) {
        ny$error <- paste("The analysis failed:", conditionMessage(attr(out, "condition")))
      } else if (!isTRUE(out$ok)) {
        ny$error <- "The file did not pass validation."
      } else {
        ny$res <- out; ny$ran <- TRUE
        setProgress(value = 0.97, detail = "Preparing the baseline comparison")
        ny_baseline()
      }
    })
  })

  output$newyear_status <- renderUI({
    if (!is.null(ny$error))
      return(div(style = "color:#b42318;", strong("Could not run. "), ny$error))
    if (!ny$ran)
      return(div(style = "color:#555;", "Upload a file, then click Run analysis."))
    div(style = "color:#1a7f37;", strong("Done. "), "See the Results tab.")
  })

  # -- Headline ---------------------------------------------------------------
  output$newyear_headline <- renderUI({
    req(ny$res)
    s <- ny$res$scores
    top <- s$scenario[1]
    gap <- if (nrow(s) > 1) s$composite[1] - s$composite[2] else NA_real_
    yrs <- paste(ny$res$parsed$met_years, collapse = ", ")
    div(style = "background:#eef4fa; border-left:4px solid #31688E; padding:12px 16px;",
        h4(style = "margin-top:0;",
           sprintf("Top-ranked alternative: %s", top)),
        p(sprintf("Composite score %.3f, ahead of %s by %.3f. Meteorological years in this file: %s.",
                  s$composite[1], s$scenario[2], gap, yrs)),
        p(style = "margin-bottom:0; color:#444;",
          "Scores are computed with the elicited weights (0.40 Chinook, 0.50 hydropower, ",
          "0.10 steelhead) and the published calibration. Hydropower revenue loss uses the ",
          "published per-scenario values."))
  })

  output$newyear_scores <- renderDT({
    req(ny$res)
    df <- ny$res$scores %>%
      transmute(Rank = rank, Alternative = scenario,
                `Adult index` = round(adult_index),
                `Steelhead (days < 18.3C)` = round(steelhead_score, 1),
                `Revenue loss ($)` = round(hydro_raw),
                `Composite` = round(composite, 4))
    datatable(df, rownames = FALSE, options = list(dom = "t", paging = FALSE))
  })

  output$newyear_score_plot <- renderPlot({
    req(ny$res)
    d <- ny$res$scores %>% mutate(scenario = factor(scenario, levels = scenario))
    ggplot(d, aes(scenario, composite, fill = rank == 1)) +
      geom_col(colour = "grey20", linewidth = 0.35, width = 0.75) +
      geom_text(aes(label = sprintf("%.3f", composite)), vjust = -0.5,
                fontface = "bold", size = 4.5) +
      scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "grey85"),
                        guide = "none") +
      scale_y_continuous(limits = c(0, max(d$composite) * 1.15), expand = c(0, 0)) +
      labs(x = "Management alternative", y = "Composite score") +
      theme_minimal(base_size = 14) +
      theme(axis.title = element_text(face = "bold", colour = "black"),
            axis.text  = element_text(face = "bold", size = 11, colour = "black"),
            panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
  })

  # -- Against the baseline ---------------------------------------------------
  ny_compare <- reactive({
    req(ny$res)
    b <- ny_baseline()
    validate(need(!is.null(b),
                  paste("The baseline could not be computed in this deployment;",
                        "data_raw/SDM Power Bypass Temperature Modeling Results.xlsx",
                        "is missing.")))
    ny$res$scores %>%
      select(scenario, new_index = adult_index, new_comp = composite, new_rank = rank) %>%
      left_join(b$scores %>% select(scenario, base_index = adult_index,
                                    base_comp = composite, base_rank = rank),
                by = "scenario") %>%
      mutate(d_index = new_index - base_index,
             d_comp  = new_comp - base_comp,
             d_rank  = base_rank - new_rank) %>%
      arrange(new_rank)
  })

  output$newyear_vs_base <- renderDT({
    d <- ny_compare()
    df <- d %>% transmute(
      Alternative = scenario,
      `Rank (new)` = new_rank, `Rank (baseline)` = base_rank,
      `Rank change` = ifelse(d_rank > 0, paste0("+", d_rank), as.character(d_rank)),
      `Adult index (new)` = round(new_index),
      `Adult index (baseline)` = round(base_index),
      `Change` = sprintf("%+.0f", d_index),
      `Composite change` = sprintf("%+.4f", d_comp))
    datatable(df, rownames = FALSE, options = list(dom = "t", paging = FALSE))
  })

  output$newyear_vs_base_plot <- renderPlot({
    d <- ny_compare()
    dd <- d %>%
      select(scenario, New = new_comp, Baseline = base_comp) %>%
      pivot_longer(-scenario, names_to = "run", values_to = "composite") %>%
      mutate(run = factor(run, levels = c("Baseline", "New")),
             scenario = factor(scenario, levels = d$scenario))
    ggplot(dd, aes(scenario, composite, fill = run)) +
      geom_col(position = position_dodge(width = 0.78), width = 0.7,
               colour = "grey20", linewidth = 0.3) +
      scale_fill_viridis_d(begin = 0.2, end = 0.75, name = NULL) +
      labs(x = "Management alternative", y = "Composite score") +
      theme_minimal(base_size = 14) +
      theme(axis.title = element_text(face = "bold", colour = "black"),
            axis.text  = element_text(face = "bold", size = 11, colour = "black"),
            legend.text = element_text(size = 11, colour = "black"),
            legend.position = "top",
            panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
  })

  # -- Temperatures -----------------------------------------------------------
  output$newyear_temp_plot <- renderPlot({
    req(ny$res)
    s <- ny$res$series %>%
      filter(site == "AveHazel", lubridate::month(Date) %in% c(10, 11))
    nb <- s %>% filter(scenario == "NB") %>% select(met_year, Date, nb_temp = temp)
    d  <- s %>% filter(scenario != "NB") %>%
      left_join(nb, by = c("met_year", "Date")) %>%
      mutate(delta = temp - nb_temp)
    ggplot(d, aes(Date, delta, colour = scenario)) +
      geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.5) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~ met_year) +
      scale_colour_viridis_d(begin = 0.05, end = 0.9, name = NULL) +
      scale_x_date(date_labels = "%b %d") +
      labs(x = NULL, y = "Difference from no-bypass (deg C)") +
      theme_minimal(base_size = 14) +
      theme(axis.title = element_text(face = "bold", colour = "black"),
            axis.text  = element_text(face = "bold", size = 10, colour = "black"),
            strip.text = element_text(face = "bold", size = 12, colour = "black"),
            legend.position = "top", panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
  })

  # -- Downloads --------------------------------------------------------------
  output$newyear_dl_scores <- downloadHandler(
    filename = function() sprintf("scenario_scores_%s.csv", Sys.Date()),
    content = function(file) {
      req(ny$res)
      write.csv(ny$res$scores, file, row.names = FALSE)
    })

  output$newyear_dl_detail <- downloadHandler(
    filename = function() sprintf("scenario_detail_%s.csv", Sys.Date()),
    content = function(file) {
      req(ny$res)
      out <- ny$res$projection %>%
        left_join(ny$res$survival, by = c("met_year", "scenario", "variant")) %>%
        left_join(ny$res$steelhead, by = c("met_year", "scenario"))
      write.csv(out, file, row.names = FALSE)
    })

  # -- Documentation ----------------------------------------------------------
  output$newyear_doc <- renderUI({
    tagList(
      h4("What happens when you press Run"),
      tags$ol(
        tags$li(strong("Your file is read and checked."), " Nothing runs until it passes."),
        tags$li(strong("Daily temperature series are built."),
                " Your scenario temperatures replace the modelled window (18 October to 31 December); ",
                "everything outside it keeps the long-term average, exactly as the published model does."),
        tags$li(strong("Spawn timing is predicted"),
                " from your October and November temperatures, using the same model the published analysis uses."),
        tags$li(strong("Egg-to-fry survival is computed"),
                " under all three temperature-dependent mortality formulations."),
        tags$li(strong("The population is projected"),
                " using the published calibration, and alternatives are scored.")
      ),
      h4("What this does not do"),
      tags$ul(
        tags$li("It does ", strong("not"), " re-calibrate the model. Survival and return rates ",
                "stay at their published values, which were fitted to observed escapement. ",
                "Adding new carcass or escapement years is a separate annual job."),
        tags$li("It does ", strong("not"), " save anything. Results disappear when you close ",
                "this tab, so download anything you want to keep."),
        tags$li("It does ", strong("not"), " change what anyone else sees.")
      ),
      h4("Why the numbers differ slightly from the manuscript"),
      p("The published run draws a random sample of redds each year; this engine weights every ",
        "spawn date by its probability instead. That removes run-to-run noise but shifts values ",
        "by a few percent. For that reason the ", strong("Against the baseline"),
        " tab compares your upload against the published years re-run through this same engine, ",
        "not against the printed tables."),
      h4("Rough guide to reading the result"),
      p("The composite score combines three objectives with the elicited weights: 0.40 Chinook, ",
        "0.50 hydropower, 0.10 steelhead. Each objective is rescaled 0 to 1 across the ",
        "alternatives in your file, so scores are ", strong("relative"),
        " and only comparable within a single run.")
    )
  })
  
  # Reactive Values
  values <- reactiveValues(
    calib_data = NULL,
    single_data = NULL,
    cmp_data = NULL,
    performance_auto = NULL
  )
  
  # Helper function for 3-slider auto-adjustment with debouncing
  make_3_slider_observers <- function(id1, id2, id3, lock) {
    # Use debounce to prevent rapid-fire updates
    observeEvent(input[[id1]], {
      if(lock()) return()
      lock(TRUE)
      
      # Add a small delay to allow the UI to settle
      Sys.sleep(0.01)
      
      remainder <- 1 - input[[id1]]
      sum_others <- isolate(input[[id2]]) + isolate(input[[id3]])
      
      if (sum_others > 0) {
        new_val2 <- remainder * (isolate(input[[id2]]) / sum_others)
        new_val3 <- remainder * (isolate(input[[id3]]) / sum_others)
      } else {
        new_val2 <- remainder / 2
        new_val3 <- remainder / 2
      }
      
      # Only update if the change is significant (more than 0.001)
      if (abs(isolate(input[[id2]]) - new_val2) > 0.001) {
        updateSliderInput(session, id2, value = new_val2)
      }
      if (abs(isolate(input[[id3]]) - new_val3) > 0.001) {
        updateSliderInput(session, id3, value = new_val3)
      }
      
      lock(FALSE)
    }, ignoreInit = TRUE, priority = 1)
    
    observeEvent(input[[id2]], {
      if(lock()) return()
      lock(TRUE)
      
      Sys.sleep(0.01)
      
      remainder <- 1 - input[[id2]]
      sum_others <- isolate(input[[id1]]) + isolate(input[[id3]])
      
      if (sum_others > 0) {
        new_val1 <- remainder * (isolate(input[[id1]]) / sum_others)
        new_val3 <- remainder * (isolate(input[[id3]]) / sum_others)
      } else {
        new_val1 <- remainder / 2
        new_val3 <- remainder / 2
      }
      
      if (abs(isolate(input[[id1]]) - new_val1) > 0.001) {
        updateSliderInput(session, id1, value = new_val1)
      }
      if (abs(isolate(input[[id3]]) - new_val3) > 0.001) {
        updateSliderInput(session, id3, value = new_val3)
      }
      
      lock(FALSE)
    }, ignoreInit = TRUE, priority = 1)
    
    observeEvent(input[[id3]], {
      if(lock()) return()
      lock(TRUE)
      
      Sys.sleep(0.01)
      
      remainder <- 1 - input[[id3]]
      sum_others <- isolate(input[[id1]]) + isolate(input[[id2]])
      
      if (sum_others > 0) {
        new_val1 <- remainder * (isolate(input[[id1]]) / sum_others)
        new_val2 <- remainder * (isolate(input[[id2]]) / sum_others)
      } else {
        new_val1 <- remainder / 2
        new_val2 <- remainder / 2
      }
      
      if (abs(isolate(input[[id1]]) - new_val1) > 0.001) {
        updateSliderInput(session, id1, value = new_val1)
      }
      if (abs(isolate(input[[id2]]) - new_val2) > 0.001) {
        updateSliderInput(session, id2, value = new_val2)
      }
      
      lock(FALSE)
    }, ignoreInit = TRUE, priority = 1)
  }
  
  # Create lock and apply observer to objective weight sliders
  objective_lock <- reactiveVal(FALSE)
  make_3_slider_observers("w_chinook", "w_steelhead", "w_hydro", objective_lock)
  
  # Run simulation for a scenario with hydro and TDM weighting
  run_scenario_simulation <- function(scenario, hydro_weights, tdm_weights, n_years, flow_val) {
    alts <- get_scenario_alternatives(scenario, "all")
    hydro_years <- c("2011", "2014", "2017", "2020")
    
    # Normalize weights
    hydro_w <- normalize_weights(hydro_weights)
    tdm_w <- normalize_weights(tdm_weights)
    
    # Use the actual flow-habitat relationship from instream data
    K_spawners <- get_K_spawners(flow_val)
    
    # Extract forecast years from results_full
    forecast_data <- results_full %>%
      filter(year >= 2025) %>%
      slice_head(n = n_years)
    
    if (nrow(forecast_data) == 0) {
      stop("No forecast data available in results_full")
    }
    
    # Initialize accumulator for final results
    final_spawners <- rep(0, n_years)
    years_vec <- sort(unique(forecast_data$year))[1:n_years]
    
    # Apply density-dependent scaling based on K change
    # The default K at 1500 cfs (from the slider default)
    base_K <- get_K_spawners(1000)   # was 1500 — now matches both slider default AND precompute
    K_scalar <- K_spawners / base_K
    
    for (i in seq_along(alts)) {
      alt_id <- as.character(alts[i])
      
      # Get data for all TDM variants for this alternative
      alt_data <- results_full %>%
        filter(env == alt_id, year %in% years_vec)
      
      if (nrow(alt_data) == 0) next
      
      # Calculate TDM-weighted spawners for this alternative
      tdm_weighted <- alt_data %>%
        filter(variant %in% c("exp_WF", "exp_SM", "lin_Martin")) %>%
        group_by(year) %>%
        summarise(
          spawners = sum(
            case_when(
              variant == "exp_WF" ~ spawners * tdm_w["wf"],
              variant == "exp_SM" ~ spawners * tdm_w["sm"],
              variant == "lin_Martin" ~ spawners * tdm_w["martin"],
              TRUE ~ 0
            ), na.rm = TRUE
          ),
          .groups = "drop"
        ) %>%
        arrange(year)
      
      # Apply carrying capacity scaling
      if (nrow(tdm_weighted) > 0) {
        tdm_spawners <- tdm_weighted$spawners[1:min(n_years, length(tdm_weighted$spawners))]
        
        # Apply density-dependent adjustment
        # Beverton-Holt style: population approaches K asymptotically
        # This better represents actual density dependence than simple scaling
        tdm_spawners <- tdm_spawners * K_spawners / (base_K + (K_spawners - base_K) * (1 - tdm_spawners/base_K))
        
        final_spawners[1:length(tdm_spawners)] <- final_spawners[1:length(tdm_spawners)] + 
          tdm_spawners * hydro_w[i]
      }
    }
    
    # Build final result dataframe
    result <- tibble(
      year = years_vec,
      spawners = final_spawners,
      scenario = scenario,
      K_spawners = K_spawners  # Include K in output
    )
    
    # Add other columns from template for compatibility
    template_cols <- results_full %>%
      filter(env == as.character(alts[1]), variant == "exp_WF") %>%
      slice_head(n = 1) %>%
      select(-year, -spawners, -env, -variant)
    
    # Don't overwrite K_spawners if it exists in template
    if (ncol(template_cols) > 0) {
      for (col in names(template_cols)) {
        if (col != "K_spawners") {
          result[[col]] <- template_cols[[col]][1]
        }
      }
    }
    
    return(result)
  }
  
  # Temperature Explorer (FIXED)
  # Debounced temp weight inputs - prevents re-render on every slider tick
  temp_w_2011_d <- debounce(reactive(input$temp_w_2011), 300)
  temp_w_2014_d <- debounce(reactive(input$temp_w_2014), 300)
  temp_w_2017_d <- debounce(reactive(input$temp_w_2017), 300)
  temp_w_2020_d <- debounce(reactive(input$temp_w_2020), 300)
  
  # Cache the filtered base data - uses df_temp_2025 (pre-filtered in global.R)
  # so no year(Date)/month(Date) calls happen at render time
  temp_base_data <- reactive({
    req(df_temp_2025, input$temp_alternatives, length(input$temp_alternatives) > 0)
    all_envs <- unlist(lapply(input$temp_alternatives,
                              function(alt) as.character(get_scenario_alternatives(alt, "all"))))
    df <- df_temp_2025 %>%
      filter(env %in% all_envs, site == input$temp_site)
    if (input$temp_period == "oct_dec") {
      df <- df %>% filter(month_num %in% c(10, 11, 12))
    }
    df
  })
  
  output$temp_plot <- renderPlot({
    req(temp_base_data(), length(input$temp_alternatives) > 0)
    weights <- normalize_weights(c("2011" = temp_w_2011_d(), "2014" = temp_w_2014_d(),
                                   "2017" = temp_w_2017_d(), "2020" = temp_w_2020_d()))
    plot_data <- map_dfr(input$temp_alternatives, function(alt) {
      alts <- as.character(get_scenario_alternatives(alt, "all"))
      temp_base_data() %>%
        filter(env %in% alts) %>%
        group_by(Date) %>%
        summarise(temp = sum(temp * weights[climate], na.rm = TRUE), .groups = "drop") %>%
        mutate(Alternative = alt)
    })
    ggplot(plot_data, aes(x = Date, y = temp, color = Alternative)) +
      geom_line(size = 1.2, alpha = 0.8) +
      scale_color_viridis_d(name = "Alternative") +
      scale_y_continuous(limits = c(0, 20)) +
      labs(title = paste("Temperature Comparison at", input$temp_site),
           subtitle = "Weighted averages based on climatology weights",
           y = "Temperature (°C)", x = "Date") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "bottom")
  })
  
  # Temperature statistics table - reuses temp_base_data() so no redundant scan
  output$temp_stats <- renderTable({
    req(temp_base_data(), length(input$temp_alternatives) > 0)
    weights <- normalize_weights(c("2011" = temp_w_2011_d(), "2014" = temp_w_2014_d(),
                                   "2017" = temp_w_2017_d(), "2020" = temp_w_2020_d()))
    map_dfr(input$temp_alternatives, function(alt) {
      alts <- as.character(get_scenario_alternatives(alt, "all"))
      weighted_temps <- temp_base_data() %>%
        filter(env %in% alts) %>%
        group_by(Date) %>%
        summarise(temp = sum(temp * weights[climate], na.rm = TRUE), .groups = "drop")
      tibble(
        Alternative = alt,
        `Mean Temp`   = round(mean(weighted_temps$temp,   na.rm = TRUE), 2),
        `Median Temp` = round(median(weighted_temps$temp, na.rm = TRUE), 2),
        `Min Temp`    = round(min(weighted_temps$temp,    na.rm = TRUE), 2),
        `Max Temp`    = round(max(weighted_temps$temp,    na.rm = TRUE), 2),
        `Std Dev`     = round(sd(weighted_temps$temp,     na.rm = TRUE), 2)
      )
    })
  }, rownames = FALSE)
  
  # Compare Alternatives
  observeEvent(input$run_cmp, {
    req(input$cmp_scenarios)
    showNotification("Running comparison...", duration = NULL, id = "notify")
    
    hydro_w_raw <- c("2011" = input$cmp_w_2011, "2014" = input$cmp_w_2014,
                     "2017" = input$cmp_w_2017, "2020" = input$cmp_w_2020)
    hydro_w <- normalize_weights(hydro_w_raw)
    tdm_w_raw <- c(wf = input$cmp_tdm_wf, sm = input$cmp_tdm_sm, martin = input$cmp_tdm_martin)
    
    cmp_results <- map_dfr(input$cmp_scenarios, ~run_scenario_simulation(., hydro_w_raw, tdm_w_raw, input$cmp_years, input$cmp_flow))
    values$cmp_data <- cmp_results
    
    # Calculate performance for Chinook using LAST 20 YEARS to match other tabs
    perf_data_chinook <- cmp_results %>%
      group_by(scenario) %>%
      slice_tail(n = 20) %>%  # Use last 20 years like boxplot
      summarise(chinook_raw = median(spawners), .groups = "drop")
    
    # Calculate weighted steelhead scores if steelhead_metrics exists
    # Calculate weighted steelhead scores if steelhead_metrics exists
    if(exists("steelhead_metrics")) {
      steelhead_weighted <- map_dfr(input$cmp_scenarios, function(scen) {
        alts <- get_scenario_alternatives(scen, "all")
        hydro_years <- c("2011", "2014", "2017", "2020")
        
        combined_steelhead <- 0
        for (j in seq_along(alts)) {
          alt_id <- as.character(alts[j])
          hydro_year <- hydro_years[j]
          
          steelhead_score <- steelhead_metrics %>% 
            filter(env == alt_id) %>% 
            pull(steelhead_score)
          
          if (length(steelhead_score) > 0) {
            combined_steelhead <- combined_steelhead + (steelhead_score * hydro_w[hydro_year])
          }
        }
        
        tibble(scenario = scen, steelhead_raw = combined_steelhead)
      })
      
      values$performance_auto <- perf_data_chinook %>%
        left_join(steelhead_weighted, by = "scenario")
    } else {
      # If no steelhead data, use placeholder values
      values$performance_auto <- perf_data_chinook %>%
        mutate(steelhead_raw = runif(n(), 40, 60))
    }
    
    removeNotification("notify")
    showNotification("Comparison complete!", duration = 2)
  })
  
  output$cmp_summary <- renderDT({
    req(values$cmp_data) %>%
      mutate(
        spawners   = round(spawners),
        pre_spawn  = round(pre_spawn, 2),
        egg_surv   = round(egg_surv, 2),
        eff_surv   = round(eff_surv, 2),
        rear_surv  = round(rear_surv, 2),
        SAR_used   = round(SAR_used, 4),
        K_spawners = round(K_spawners, 0)
      ) %>%
      dplyr::select(
        Alternative = scenario,
        Year = year,
        `Forecasted Spawners` = spawners,
        `Pre-Spawn Survival` = pre_spawn,
        `Egg Survival` = egg_surv,
        `Effective Survival` = eff_surv,
        `Rearing Survival` = rear_surv,
        `Smolt-to-Adult Ratio (SAR)` = SAR_used,
        `Spawning Capacity (K)` = K_spawners
      ) %>%
      DT::datatable(
        rownames = FALSE,
        extensions = 'Buttons',
        options = list(dom = 'Bfrtip', buttons = c('csv'))
      )
  })
  
  output$cmp_ts_plot <- renderPlot({
    df <- req(values$cmp_data)
    
    # Get the actual range of years in the data
    year_range <- range(df$year, na.rm = TRUE)
    
    ggplot(df, aes(x = year, y = spawners, color = factor(scenario), group = scenario)) +
      geom_line(size = 1.2) +
      expand_limits(y = 0) +
      scale_x_continuous(limits = year_range, 
                         breaks = seq(year_range[1], year_range[2], by = 25)) +
      scale_y_continuous(labels = comma) +
      scale_color_viridis_d(name = "Alternative") +
      labs(title = "Comparison of Weighted Alternatives",
           x = "Year", y = "Forecasted Spawner Abundance") +
      theme_minimal(base_size = 16) +
      theme(legend.position = "bottom")
  })
  
  output$cmp_box_plot <- renderPlot({
    req(values$cmp_data)
    last_yr <- max(values$cmp_data$year)
    df_filt <- values$cmp_data %>%
      filter(year >= (last_yr - input$last_n + 1))
    
    ggplot(df_filt, aes(x = factor(scenario), y = spawners, fill = factor(scenario))) +
      geom_boxplot() +
      scale_fill_viridis_d(guide = "none") +
      scale_y_continuous(labels = comma) +
      labs(title = paste0("Spawner Distribution: Last ", input$last_n, " Years"),
           x = "Alternative", y = "Forecasted Spawner Abundance") +
      theme_minimal(base_size = 16)
  })
  
  output$cmp_boxplot_stats <- renderTable({
    req(values$cmp_data)
    last_yr <- max(values$cmp_data$year)
    values$cmp_data %>%
      filter(year >= (last_yr - input$last_n + 1)) %>%
      group_by(scenario) %>%
      summarise(
        Minimum = min(spawners, na.rm = TRUE),
        `1st Qu.` = quantile(spawners, 0.25, na.rm = TRUE),
        Median = median(spawners, na.rm = TRUE),
        `3rd Qu.` = quantile(spawners, 0.75, na.rm = TRUE),
        Maximum = max(spawners, na.rm = TRUE),
        .groups = "drop") %>%
      mutate(across(where(is.numeric), round, 0))
  }, rownames = FALSE)
  
  # Decision Support
  performance_data_full <- reactive({
    # Check if user has run Compare Alternatives with custom weights
    if (!is.null(values$performance_auto)) {
      # Use the data from Compare Alternatives for consistency
      perf_data <- values$performance_auto
    } else {
      # Use pre-computed swing results (based on default weights)
      perf_data <- swing_scenario_results %>%
        rename(chinook_raw = spawner_metric) %>%
        left_join(
          steelhead_scenario_results %>% rename(steelhead_raw = steelhead_score),
          by = "scenario"
        )
    }
    
    hydro_df <- tibble(
      scenario = names(hardcoded_hydro_scores),
      hydro_raw = hardcoded_hydro_scores
    )
    
    perf_data %>%
      left_join(hydro_df, by = "scenario") %>%
      mutate(hydro_raw = ifelse(is.na(hydro_raw), 50, hydro_raw)) %>%
      mutate(
        chinook_norm = normalize_scores_chinook(chinook_raw),
        steelhead_norm = normalize_scores_steelhead(steelhead_raw),
        hydro_norm = normalize_scores_hydro(hydro_raw)
      )
  })
  
  objective_weights <- reactive({
    if (input$weight_method == "equal") {
      weights <- c(chinook = 1/3, steelhead = 1/3, hydro = 1/3)
    } else {
      weights <- c(chinook = input$w_chinook, steelhead = input$w_steelhead, hydro = input$w_hydro)
    }
    return(weights)
  })
  
  output$performance_matrix <- renderTable({
    req(performance_data_full()) %>%
      select(
        Alternative = scenario,
        `Chinook (Raw)` = chinook_raw,
        `Chinook (0-1)` = chinook_norm,
        `Steelhead (Raw)` = steelhead_raw,
        `Steelhead (0-1)` = steelhead_norm,
        `Hydro (Raw)` = hydro_raw,
        `Hydro (0-1)` = hydro_norm
      ) %>%
      mutate(
        `Chinook (Raw)` = round(`Chinook (Raw)`, 0),
        `Steelhead (Raw)` = round(`Steelhead (Raw)`, 2),
        `Hydro (Raw)` = round(`Hydro (Raw)`, 0)
      )
  })
  
  output$performance_plot <- renderPlot({
    df <- req(performance_data_full())
    df %>%
      select(scenario, chinook_norm, steelhead_norm, hydro_norm) %>%
      pivot_longer(cols = -scenario, names_to = "objective", values_to = "score") %>%
      mutate(objective = str_remove(objective, "_norm")) %>%
      ggplot(aes(x = scenario, y = score, fill = objective)) +
      geom_col(position = "dodge") +
      scale_y_continuous(limits = c(0, 1)) +
      scale_fill_viridis_d(name = "Objective") +
      labs(title = "Normalized Performance Scores by Objective",
           x = "Alternative", y = "Normalized Score (0-1)") +
      theme_minimal(base_size = 16)
  })
  
  output$tradeoff_plot <- renderPlot({
    df <- req(performance_data_full())
    ggplot(df, aes(x = hydro_norm, y = chinook_norm)) +
      geom_point(aes(color = scenario), size = 5) +
      geom_text_repel(aes(label = scenario), size = 5, box.padding = 0.5) +
      scale_x_continuous(limits = c(0, 1)) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_color_viridis_d(guide = "none") +
      labs(title = "Trade-off: Hydropower vs. Chinook Performance",
           x = "Hydropower Score (Normalized)",
           y = "Chinook Score (Normalized)") +
      theme_minimal(base_size = 16)
  })
  
  output$weights_table <- renderTable({
    weights <- objective_weights()
    tibble(
      Objective = c("Fall-run Chinook", "Steelhead", "Hydropower"),
      Weight = scales::percent(weights, accuracy = 0.1)
    )
  })
  
  output$overall_scores_plot <- renderPlot({
    req(performance_data_full())
    weights <- objective_weights()
    
    scores <- performance_data_full() %>%
      mutate(overall_score = (chinook_norm * weights["chinook"]) +
               (steelhead_norm * weights["steelhead"]) +
               (hydro_norm * weights["hydro"]))
    
    ggplot(scores, aes(x = scenario, y = overall_score, fill = scenario)) +
      geom_col() +
      scale_fill_viridis_d(guide = "none") +
      scale_y_continuous(limits = c(0, 1), labels = percent) +
      labs(title = "Overall Performance Scores", y = "Total Weighted Score") +
      theme_minimal(base_size = 14)
  })
  
  # ADDED: Stacked bar plot for score contribution
  output$plot_score_contribution <- renderPlot({
    req(performance_data_full())
    weights <- objective_weights()
    
    df <- performance_data_full() %>%
      mutate(
        w_Chinook = chinook_norm * weights["chinook"],
        w_Steelhead = steelhead_norm * weights["steelhead"],
        w_Hydropower = hydro_norm * weights["hydro"]
      ) %>%
      select(scenario, w_Chinook, w_Steelhead, w_Hydropower) %>%
      pivot_longer(cols = starts_with("w_"), names_to = "objective", values_to = "contribution") %>%
      mutate(objective = str_remove(objective, "w_"))
    
    ggplot(df, aes(x = scenario, y = contribution, fill = objective)) +
      geom_col(position = "stack") +
      scale_fill_viridis_d(name = "Objective") +
      labs(title = "Contribution to Score by Objective",
           x = "Management Alternative",
           y = "Weighted Score Contribution") +
      theme_minimal(base_size = 16) +
      theme(legend.position = "bottom")
  })
  # Validation for rankings
  output$swing_validation <- renderUI({
    ranks <- c(input$rank_chinook, input$rank_steelhead, input$rank_hydropower)
    scores <- c(input$score_chinook, input$score_steelhead, input$score_hydropower)
    
    messages <- list()
    
    # Check for duplicate ranks
    if (length(unique(ranks)) != 3) {
      messages <- c(messages, list(
        tags$div(class = "alert alert-warning", 
                 icon("exclamation-triangle"), 
                 " Each alternative must have a unique rank (1, 2, or 3).")
      ))
    }
    
    # Check that scores align with ranks
    rank_order <- order(ranks)
    score_order <- order(-scores)  # Negative for descending
    
    if (!identical(rank_order, score_order)) {
      messages <- c(messages, list(
        tags$div(class = "alert alert-warning",
                 icon("exclamation-triangle"),
                 " Scores should match rankings: Rank 1 should have the highest score, rank 3 the lowest.")
      ))
    }
    
    if (length(messages) == 0) {
      tags$div(class = "alert alert-success",
               icon("check-circle"),
               " Rankings and scores are valid!")
    } else {
      do.call(tagList, messages)
    }
  })
  
  # Display objective ranges
  output$swing_ranges_table <- renderTable({
    # Get the current performance data to find actual ranges
    perf_data <- performance_data_full()
    
    tibble(
      Objective = c("Fall-run Chinook", "Steelhead", "Hydropower"),
      Direction = c("Maximize", "Maximize", "Minimize"),
      `Worst Case` = c(
        min(perf_data$chinook_raw, na.rm = TRUE),
        min(perf_data$steelhead_raw, na.rm = TRUE),
        max(perf_data$hydro_raw, na.rm = TRUE)
      ),
      `Best Case` = c(
        max(perf_data$chinook_raw, na.rm = TRUE),
        max(perf_data$steelhead_raw, na.rm = TRUE),
        min(perf_data$hydro_raw, na.rm = TRUE)
      )
    ) %>%
      mutate(
        `Worst Case` = if_else(Objective == "Steelhead", round(`Worst Case`, 2), round(`Worst Case`, 0)),
        `Best Case` = if_else(Objective == "Steelhead", round(`Best Case`, 2), round(`Best Case`, 0))
      )
  }, align = 'lccc')
  
  # Display hypothetical alternatives
  output$swing_alternatives_table <- renderTable({
    perf_data <- performance_data_full()
    
    tibble(
      Alternative = c("Worst Alternative", "Alt 1: Best Chinook", "Alt 2: Best Steelhead", "Alt 3: Best Hydropower"),
      `Chinook Abundance` = c(
        min(perf_data$chinook_raw, na.rm = TRUE),
        max(perf_data$chinook_raw, na.rm = TRUE),
        min(perf_data$chinook_raw, na.rm = TRUE),
        min(perf_data$chinook_raw, na.rm = TRUE)
      ),
      `Steelhead Score` = c(
        swing_ranges$worst_case[swing_ranges$objective == "Steelhead"],
        swing_ranges$worst_case[swing_ranges$objective == "Steelhead"],
        swing_ranges$best_case[swing_ranges$objective == "Steelhead"],
        swing_ranges$worst_case[swing_ranges$objective == "Steelhead"]
      ),
      `Hydropower Cost` = c(
        max(hardcoded_hydro_scores),
        max(hardcoded_hydro_scores),
        max(hardcoded_hydro_scores),
        min(hardcoded_hydro_scores)
      )
    ) %>%
      mutate(
        `Chinook Abundance` = round(`Chinook Abundance`, 0),
        `Steelhead Score` = round(`Steelhead Score`, 2),
        `Hydropower Cost` = round(`Hydropower Cost`, 0)
      )
  }, align = 'lccc')
  
  # Calculate and display weights
  swing_weights_calculated <- reactive({
    scores <- c(
      Chinook = input$score_chinook,
      Steelhead = input$score_steelhead,
      Hydropower = input$score_hydropower
    )
    
    total <- sum(scores)
    if (total == 0) {
      weights <- c(Chinook = 0.33, Steelhead = 0.33, Hydropower = 0.34)
    } else {
      weights <- scores / total
    }
    
    # Debug: print to console
    print(paste("Calculated weights:", weights))
    
    tibble(
      Objective = c("Fall-run Chinook", "Steelhead", "Hydropower"),
      Score = scores,
      Weight = weights,
      `Weight %` = scales::percent(weights, accuracy = 0.1)
    )
  })
  
  output$swing_weights_table <- renderTable({
    swing_weights_calculated()
  }, digits = 3)
  
  # Plot weights
  output$swing_weights_plot <- renderPlot({
    df <- swing_weights_calculated()
    
    ggplot(df, aes(x = Objective, y = Weight, fill = Objective)) +
      geom_col() +
      geom_text(aes(label = `Weight %`), vjust = -0.5, size = 5) +
      scale_fill_viridis_d(guide = "none") +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
      labs(title = "Calculated Objective Weights from Swing Weighting",
           y = "Weight", x = NULL) +
      theme_minimal(base_size = 14)
  })
  
  # Apply weights to the main analysis
  observeEvent(input$apply_swing_weights, {
    # Get the calculated weights
    weights_df <- swing_weights_calculated()
    
    # Extract individual weight values
    w_chinook <- weights_df$Weight[weights_df$Objective == "Fall-run Chinook"]
    w_steelhead <- weights_df$Weight[weights_df$Objective == "Steelhead"]
    w_hydro <- weights_df$Weight[weights_df$Objective == "Hydropower"]
    
    # Lock the observers to prevent cascading updates
    objective_lock(TRUE)
    
    # Update the weight sliders with explicit values
    updateSliderInput(session, "w_chinook", value = as.numeric(w_chinook))
    updateSliderInput(session, "w_steelhead", value = as.numeric(w_steelhead))
    updateSliderInput(session, "w_hydro", value = as.numeric(w_hydro))
    
    # Switch to manual weights mode
    updateRadioButtons(session, "weight_method", selected = "manual")
    
    # Small delay to let updates complete
    Sys.sleep(0.1)
    
    # Unlock the observers
    objective_lock(FALSE)
    
    # Show confirmation
    showNotification(
      paste0("Swing weights applied! Chinook: ", round(w_chinook, 2), 
             ", Steelhead: ", round(w_steelhead, 2), 
             ", Hydropower: ", round(w_hydro, 2)),
      type = "message", 
      duration = 5
    )
  })
}

shinyApp(ui, server)