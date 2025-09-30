library(shiny)
library(tidyverse)
library(DT)
library(scales)
library(shinyjs) # Using shinyjs for smoother slider updates
library(shinyWidgets) # For styled buttons

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

# Helper for min-max normalization (scales results to 0-100)
normalize_scores <- function(scores) {
  min_s <- min(scores, na.rm = TRUE)
  max_s <- max(scores, na.rm = TRUE)
  if (max_s == min_s) return(rep(100, length(scores))) # Return 100 if all values are the same
  (scores - min_s) / (max_s - min_s) * 100 
}

# Helper for min-max normalization (scales results to 0-1), one for objectives we're maximizing and another for ones we're minimizing
normalize_scores_max <- function(scores) {
  min_s <- min(scores, na.rm = TRUE)
  max_s <- max(scores, na.rm = TRUE)
  if (max_s == min_s) return(rep(100, length(scores))) # Return 1 if all values are the same
  (scores - min_s) / (max_s - min_s)
}

normalize_scores_min <- function(scores) {
  min_s <- min(scores, na.rm = TRUE)
  max_s <- max(scores, na.rm = TRUE)
  if (max_s == min_s) return(rep(100, length(scores))) # Return 1 if all values are the same
  (max_s - scores) / (max_s - min_s)
}

# Map 28 alternatives to scenarios and hydro years
get_scenario_alternatives <- function(scenario, hydro_year) {
  # Mapping: Alt 1-7 (2011), 8-14 (2014), 15-21 (2017), 22-28 (2020)
  # Within each hydro year: 1=NB, 2=PB1, 3=PB2, 4=PB3, 5=PB4, 6=PB5, 7=PB6
  scenario_map <- c("NB"=1, "PB1"=2, "PB2"=3, "PB3"=4, "PB4"=5, "PB5"=6, "PB6"=7)
  hydro_map <- c("2011"=0, "2014"=7, "2017"=14, "2020"=21)
  
  if (hydro_year == "all") {
    # Return alternatives for this scenario across all hydro years
    base_idx <- scenario_map[scenario]
    return(c(base_idx, base_idx+7, base_idx+14, base_idx+21))
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
                                   p("This interactive application simulates fall-run Chinook salmon population dynamics on the American River under different Folsom Dam power bypass flow and temperature management alternatives. The model integrates temperature-dependent mortality, spawn timing, and comprehensive life-cycle processes to project spawner abundance through 2150."),
                                   
                                   h4("Management Structure:"),
                                   tags$ul(
                                     tags$li(strong("28 Total Alternatives:"), "7 scenarios × 4 hydrological year types"),
                                     tags$li(strong("Scenarios:"), "No Bypass (NB) and 6 Power Bypass configurations (PB1-PB6) with varying flow rates and timing"),
                                     tags$li(strong("Hydrological Years:"), "2011 (Dry), 2014 (Critical), 2017 (Wet), 2020 (Below Normal)"),
                                     tags$li(strong("Temperature Data:"), "Oct 18 - Dec 31 forecast window from SDM Power Bypass modeling, combined with 14-year USGS gauge climatology (2011-2025) for full annual cycles"),
                                     tags$li(strong("Simulation Period:"), "2025-2150 with dynamic weighting of hydrological conditions and TDM models")
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
                                     tags$li("Calibrated on 2011-2024 American River carcass survey data (n = observations)")
                                   ),
                                   
                                   h4("3. Temperature-Dependent Mortality (TDM)"),
                                   p("Three mortality models calculate egg-to-fry survival from spawn to emergence (1375 ATU):"),
                                   
                                   tags$div(style = "margin-left: 20px;",
                                            tags$h5("Exponential Models:"),
                                            tags$ul(
                                              tags$li(HTML("<strong>Water Forum 2020:</strong> Recent calibration with stage-specific parameters")),
                                              tags$li(HTML("&nbsp;&nbsp;• Egg stage (0-958 ATU): α = 3.408×10⁻¹¹, β = 1.211")),
                                              tags$li(HTML("&nbsp;&nbsp;• Alevin stage (958-1375 ATU): α = 1.018×10⁻¹⁰, β = 1.241")),
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
                                     tags$li(HTML("<strong>Carrying Capacity:</strong> K = 12,493 spawners (or flow-dependent if specified)")),
                                     tags$li(HTML("<strong>Fry Production:</strong> Fry = Eggs × S<sub>TDM</sub> × dd")),
                                     tags$li(HTML("<strong>Smolt Production:</strong> Smolts = Fry × rear_surv (calibrated)")),
                                     tags$li(HTML("<strong>Ocean Survival:</strong> SAR (smolt-to-adult return) calibrated to 2011-2024 escapement")),
                                     tags$li(HTML("<strong>Age Structure:</strong> 82.9% age-3, 16.9% age-4, 0.2% age-5 returns (CWT data)")),
                                     tags$li(HTML("<strong>Returns:</strong> Spawners<sub>t+a</sub> = Smolts<sub>t</sub> × SAR × P(age=a)"))
                                   ),
                                   
                                   h3("Calibration & Validation"),
                                   
                                   h4("Calibration Procedure:"),
                                   tags$ol(
                                     tags$li("Years 2011-2013 serve as initial population seed"),
                                     tags$li("Model jointly optimizes SAR_mean and rear_surv parameters"),
                                     tags$li("Minimizes SSE between predicted and observed spawners for 2014-2024"),
                                     tags$li("Observed data: CDFW GrandTab escapement estimates"),
                                     tags$li("Each TDM variant calibrated independently, then weighted ensemble created")
                                   ),
                                   
                                   h4("Data Sources:"),
                                   tags$ul(
                                     tags$li(strong("Escapement:"), "CDFW GrandTab American River fall-run Chinook (2011-2024)"),
                                     tags$li(strong("Carcass Surveys:"), "American River sections NB, W, 1a, 1b, 2, 3 (2011-2024)"),
                                     tags$li(strong("Temperature:"), "USGS stations 11446980, 11446500 (15-minute intervals aggregated to daily)"),
                                     tags$li(strong("Age Structure:"), "Coded wire tag recoveries from Central Valley database"),
                                     tags$li(strong("Power Bypass Modeling:"), "SDM temperature projections for 7 operational scenarios")
                                   ),
                                   
                                   h3("Key Equations Reference"),
                                   tags$div(style = "background-color: #f5f5f5; padding: 15px; border-radius: 5px;",
                                            tags$ol(
                                              tags$li(HTML("<strong>ATU Accumulation:</strong> ATU = Σ max(T<sub>daily</sub>, 0)"), 
                                                      br(), tags$small(em("Accumulated thermal units drive developmental timing"))),
                                              
                                              tags$li(HTML("<strong>Egg Development:</strong> Days to hatch = 958 / T̄"), 
                                                      br(), tags$small(em("Temperature-dependent development rate to reach alevin stage"))),
                                              
                                              tags$li(HTML("<strong>Alevin Emergence:</strong> Days to emerge = 417 / T̄"), 
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
                                     tags$li(strong("Temperature Explorer:"), "Visualize weighted temperature patterns across scenarios and hydrological years"),
                                     tags$li(strong("Calibration Tab:"), "Adjust TDM model weights and examine historical fit"),
                                     tags$li(strong("Single Scenario Analysis:"), "Detailed forecasts with custom hydrology and TDM weighting"),
                                     tags$li(strong("Scenario Comparison:"), "Side-by-side evaluation of multiple alternatives"),
                                     tags$li(strong("Swing Weighting:"), "Multi-objective decision support with importance ratings"),
                                     tags$li(strong("Consequence Table:"), "Raw and normalized (0-1) scores for all objectives")
                                   ),
                                   
                                   h4("Weighting Options:"),
                                   tags$ul(
                                     tags$li(strong("Hydrology Weights:"), "Adjust relative importance of 4 water year types (sum to 1.0)"),
                                     tags$li(strong("TDM Model Weights:"), "Combine 3 mortality models with custom weights (sum to 1.0)"),
                                     tags$li(strong("Default Weights:"), "TDM: 51% Water Forum, 24% SALMOD, 25% Martin; Hydrology: 25% each")
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
                                     tags$li("Bartholow, J.M. & Heasley, J. (2006). Evaluation of Shasta Dam Scenarios Using a Salmon Production Model. USGS Open-File Report 2004-1351."),
                                     tags$li("Bratovich, P., Neal, M., Ransom, A., et al. (2020). Chinook Salmon Early Lifestage Survival and Folsom Dam Power Bypass Considerations. Water Forum Technical Memorandum."),
                                     tags$li("Colvin, R., Falke, J.A., Henson, S. (2018). Identifying optimal water temperature and flow regimes for anadromous fish. River Research and Applications 34(6):621-632."),
                                     tags$li("Martin, B.T., Pike, A., John, S.N., Hamda, N., Roberts, J., Lindley, S.T., Danner, E.M. (2017). Phenomenological vs. biophysical models of thermal stress in aquatic eggs. Ecology Letters 20:50-59."),
                                     tags$li("USFWS (2006). SALMOD: Salmon Population Model Version 3.0. Sacramento Fish and Wildlife Office."),
                                     tags$li("CDFW GrandTab (2024). California Central Valley Chinook Salmon Escapement Database. Available at: https://wildlife.ca.gov/Conservation/Fishes/Chinook-Salmon/Anadromous-Assessment")
                                   ),
                                   
                                   h3("Version History"),
                                   tags$ul(
                                     tags$li(strong("v3.0 (Current):"), "28-alternative structure with SDM Power Bypass temperature integration"),
                                     tags$li(strong("v2.0:"), "10-alternative version with preliminary temperature scenarios"),
                                     tags$li(strong("v1.0:"), "Initial 5-alternative proof of concept")
                                   )
                            )
                          )
                 ),
                 
                 # ---- Temperature Explorer Tab (FIXED) ----
                 tabPanel("Temperature Explorer",
                          sidebarLayout(
                            sidebarPanel(
                              h4("Hydrology Weights"),
                              sliderInput("temp_w_2011", "2011 (Wet)", value = 0.25, min = 0, max = 1, step = 0.01),
                              sliderInput("temp_w_2014", "2014 (Critical)", value = 0.25, min = 0, max = 1, step = 0.01),
                              sliderInput("temp_w_2017", "2017 (Wet)", value = 0.25, min = 0, max = 1, step = 0.01),
                              sliderInput("temp_w_2020", "2020 (Dry)", value = 0.25, min = 0, max = 1, step = 0.01),
                              hr(),
                              selectInput("temp_scenario", "Scenario:",
                                          choices = c("No Bypass (NB)"="NB", "Power Bypass 1"="PB1", "Power Bypass 2"="PB2",
                                                      "Power Bypass 3"="PB3", "Power Bypass 4"="PB4", "Power Bypass 5"="PB5",
                                                      "Power Bypass 6"="PB6")),
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
                 
                 
                 # ---- Compare Scenarios Tab ----
                 tabPanel("Compare Scenarios",
                          sidebarLayout(
                            sidebarPanel(
                              h4("Scenarios to Compare"),
                              checkboxGroupInput("cmp_scenarios", "Select:",
                                                 choices = c("No Bypass"="NB", "Power Bypass 1"="PB1", "Power Bypass 2"="PB2",
                                                             "Power Bypass 3"="PB3", "Power Bypass 4"="PB4",
                                                             "Power Bypass 5"="PB5", "Power Bypass 6"="PB6"),
                                                 selected = c("NB", "PB1", "PB2", "PB3", "PB4", "PB5", "PB6")),
                              hr(),
                              h4("Hydrology Weights"),
                              sliderInput("cmp_w_2011", "2011 (Dry)", value = 0.25, min = 0, max = 1, step = 0.01),
                              sliderInput("cmp_w_2014", "2014 (Critical)", value = 0.25, min = 0, max = 1, step = 0.01),
                              sliderInput("cmp_w_2017", "2017 (Wet)", value = 0.25, min = 0, max = 1, step = 0.01),
                              sliderInput("cmp_w_2020", "2020 (Below Normal)", value = 0.25, min = 0, max = 1, step = 0.01),
                              hr(),
                              h4("TDM Weights"),
                              sliderInput("cmp_tdm_wf", "Water Forum", value = 0.51, min = 0, max = 1, step = 0.01),
                              sliderInput("cmp_tdm_sm", "SALMOD", value = 0.24, min = 0, max = 1, step = 0.01),
                              sliderInput("cmp_tdm_martin", "Martin", value = 0.25, min = 0, max = 1, step = 0.01),
                              hr(),
                              sliderInput("cmp_years", "Forecast Years", value = 50, min = 10, max = 100),
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
                 
                 # ---- Decision Support Tab (WITH STACKED BAR PLOT) ----
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
                                sliderInput("w_chinook", "Fall-run Chinook", min = 0, max = 1, value = 0.4, step = 0.05),
                                sliderInput("w_steelhead", "Steelhead", min = 0, max = 1, value = 0.3, step = 0.05),
                                sliderInput("w_hydro", "Hydropower", min = 0, max = 1, value = 0.3, step = 0.05)
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
                 )
)

server <- function(input, output, session) {
  
  # Hard-coded Hydropower Scores
  hardcoded_hydro_scores <- c(
    "NB"  = 0, "PB1" = 111422, "PB2" = 370826, "PB3" = 201552,
    "PB4" = 241590, "PB5" = 199382, "PB6" = 348806
  )
  
  # Reactive Values
  values <- reactiveValues(
    calib_data = NULL,
    single_data = NULL,
    cmp_data = NULL,
    performance_auto = NULL
  )
  
  # Helper function for 3-slider auto-adjustment
  make_3_slider_observers <- function(id1, id2, id3, lock) {
    observeEvent(input[[id1]], {
      if(lock()) return()
      lock(TRUE)
      remainder <- 1 - input[[id1]]
      sum_others <- isolate(input[[id2]]) + isolate(input[[id3]])
      if (sum_others > 0) {
        updateSliderInput(session, id2, value = remainder * (isolate(input[[id2]]) / sum_others))
        updateSliderInput(session, id3, value = remainder * (isolate(input[[id3]]) / sum_others))
      } else {
        updateSliderInput(session, id2, value = remainder / 2)
        updateSliderInput(session, id3, value = remainder / 2)
      }
      lock(FALSE)
    }, ignoreInit = TRUE)
    
    observeEvent(input[[id2]], {
      if(lock()) return()
      lock(TRUE)
      remainder <- 1 - input[[id2]]
      sum_others <- isolate(input[[id1]]) + isolate(input[[id3]])
      if (sum_others > 0) {
        updateSliderInput(session, id1, value = remainder * (isolate(input[[id1]]) / sum_others))
        updateSliderInput(session, id3, value = remainder * (isolate(input[[id3]]) / sum_others))
      } else {
        updateSliderInput(session, id1, value = remainder / 2)
        updateSliderInput(session, id3, value = remainder / 2)
      }
      lock(FALSE)
    }, ignoreInit = TRUE)
    
    observeEvent(input[[id3]], {
      if(lock()) return()
      lock(TRUE)
      remainder <- 1 - input[[id3]]
      sum_others <- isolate(input[[id1]]) + isolate(input[[id2]])
      if (sum_others > 0) {
        updateSliderInput(session, id1, value = remainder * (isolate(input[[id1]]) / sum_others))
        updateSliderInput(session, id2, value = remainder * (isolate(input[[id2]]) / sum_others))
      } else {
        updateSliderInput(session, id1, value = remainder / 2)
        updateSliderInput(session, id2, value = remainder / 2)
      }
      lock(FALSE)
    }, ignoreInit = TRUE)
  }
  
  # Create lock and apply observer to objective weight sliders
  objective_lock <- reactiveVal(FALSE)
  make_3_slider_observers("w_chinook", "w_steelhead", "w_hydro", objective_lock)
  
  # Run simulation for a scenario with hydro weighting
  run_scenario_simulation <- function(scenario, hydro_weights, tdm_weights, n_years) {
    alts <- get_scenario_alternatives(scenario, "all")
    hydro_years <- c("2011", "2014", "2017", "2020")
    hydro_w <- normalize_weights(hydro_weights)
    tdm_w <- normalize_weights(tdm_weights)
    combined_result <- NULL
    
    for (i in seq_along(alts)) {
      alt_id <- as.character(alts[i])
      hydro_year <- hydro_years[i]
      
      surv_weighted <- tdm_w["wf"] * surv_lookup_full[[paste(alt_id, "exp_WF", sep="_")]] +
        tdm_w["sm"] * surv_lookup_full[[paste(alt_id, "exp_SM", sep="_")]] +
        tdm_w["martin"] * surv_lookup_full[[paste(alt_id, "lin_Martin", sep="_")]]
      
      P_weighted <- base_P
      P_weighted$SAR_mean <- tdm_w["wf"] * base_P_list$exp_WF[[alt_id]]$SAR_mean +
        tdm_w["sm"] * base_P_list$exp_SM[[alt_id]]$SAR_mean +
        tdm_w["martin"] * base_P_list$lin_Martin[[alt_id]]$SAR_mean
      P_weighted$rear_surv <- tdm_w["wf"] * base_P_list$exp_WF[[alt_id]]$rear_surv +
        tdm_w["sm"] * base_P_list$exp_SM[[alt_id]]$rear_surv +
        tdm_w["martin"] * base_P_list$lin_Martin[[alt_id]]$rear_surv
      
      seed_weighted <- tdm_w["wf"] * S_seed_fore_list$exp_WF +
        tdm_w["sm"] * S_seed_fore_list$exp_SM +
        tdm_w["martin"] * S_seed_fore_list$lin_Martin
      
      surv_lookup_full[[paste(alt_id, "weighted", sep="_")]] <<- surv_weighted
      if (!"weighted" %in% names(base_P_list)) base_P_list[["weighted"]] <<- list()
      base_P_list[["weighted"]][[alt_id]] <<- P_weighted
      
      fc_fn <- sim_forecast_fn("weighted", alt_id, NULL, seed_weighted, spawn_dates_by_alt)
      result <- fc_fn() %>% filter(year > max(real_years)) %>% slice_head(n = n_years)
      
      if (is.null(combined_result)) {
        combined_result <- result
        combined_result$spawners <- combined_result$spawners * hydro_w[i]
      } else {
        combined_result$spawners <- combined_result$spawners + (result$spawners * hydro_w[i])
      }
    }
    combined_result$scenario <- scenario
    return(combined_result)
  }
  
  # Temperature Explorer (FIXED)
  output$temp_plot <- renderPlot({
    req(df_all_orig)
    alts <- get_scenario_alternatives(input$temp_scenario, "all")
    
    # Filter for selected period
    if (input$temp_period == "oct_dec") {
      temp_data <- df_all_orig %>%
        filter(env %in% as.character(alts),
               site == input$temp_site,
               month(Date) %in% c(10, 11, 12),
               year(Date) == 2025)
    } else {
      temp_data <- df_all_orig %>%
        filter(env %in% as.character(alts),
               site == input$temp_site,
               year(Date) == 2025)
    }
    
    # Add hydro year labels
    temp_data <- temp_data %>%
      mutate(hydro = case_when(
        env %in% as.character(1:7) ~ "2011",
        env %in% as.character(8:14) ~ "2014",
        env %in% as.character(15:21) ~ "2017",
        env %in% as.character(22:28) ~ "2020"
      ))
    
    # Calculate weighted average
    weights <- c("2011" = input$temp_w_2011, "2014" = input$temp_w_2014,
                 "2017" = input$temp_w_2017, "2020" = input$temp_w_2020)
    weights <- normalize_weights(weights)
    
    # FIXED: Calculate average only from the filtered data
    avg_temp <- temp_data %>%
      group_by(Date) %>%
      summarise(temp = sum(temp * weights[hydro]), .groups = "drop")
    
    ggplot() +
      geom_line(data = temp_data, aes(Date, temp, color = hydro), alpha = 0.5, size = 1) +
      geom_line(data = avg_temp, aes(Date, temp), color = "black", size = 1.5) +
      scale_color_viridis_d(name = "Hydro Year") +
      labs(title = paste("Temperature:", input$temp_scenario, "at", input$temp_site),
           subtitle = "Black line = weighted average",
           y = "Temperature (°C)", x = "Date") +
      theme_minimal(base_size = 14)
  })
  
  # FIXED: Temperature stats without steelhead metric
  output$temp_stats <- renderTable({
    req(df_all_orig)
    alts <- get_scenario_alternatives(input$temp_scenario, "all")
    
    if (input$temp_period == "oct_dec") {
      temp_data <- df_all_orig %>%
        filter(env %in% as.character(alts),
               site == input$temp_site,
               month(Date) %in% c(10, 11, 12),
               year(Date) == 2025)
    } else {
      temp_data <- df_all_orig %>%
        filter(env %in% as.character(alts),
               site == input$temp_site,
               year(Date) == 2025)
    }
    
    temp_data <- temp_data %>%
      mutate(hydro = case_when(
        env %in% as.character(1:7) ~ "2011",
        env %in% as.character(8:14) ~ "2014",
        env %in% as.character(15:21) ~ "2017",
        env %in% as.character(22:28) ~ "2020"
      ))
    
    weights <- c("2011" = input$temp_w_2011, "2014" = input$temp_w_2014,
                 "2017" = input$temp_w_2017, "2020" = input$temp_w_2020)
    weights <- normalize_weights(weights)
    
    weighted_temps <- temp_data %>%
      group_by(Date) %>%
      summarise(temp = sum(temp * weights[hydro]), .groups = "drop")
    
    tibble(
      Metric = c("Mean Temp", "Median Temp", "Min Temp", "Max Temp", "Std Dev", "Days > 12.8°C"),
      Value = c(
        round(mean(weighted_temps$temp), 2),
        round(median(weighted_temps$temp), 2),
        round(min(weighted_temps$temp), 2),
        round(max(weighted_temps$temp), 2),
        round(sd(weighted_temps$temp), 2),
        sum(weighted_temps$temp > 12.8)
      )
    )
  })
  
  # Compare Scenarios
  observeEvent(input$run_cmp, {
    req(input$cmp_scenarios)
    showNotification("Running comparison...", duration = NULL, id = "notify")
    
    hydro_w_raw <- c("2011" = input$cmp_w_2011, "2014" = input$cmp_w_2014,
                     "2017" = input$cmp_w_2017, "2020" = input$cmp_w_2020)
    hydro_w <- normalize_weights(hydro_w_raw)
    tdm_w_raw <- c(wf = input$cmp_tdm_wf, sm = input$cmp_tdm_sm, martin = input$cmp_tdm_martin)
    
    cmp_results <- map_dfr(input$cmp_scenarios, ~run_scenario_simulation(., hydro_w_raw, tdm_w_raw, input$cmp_years))
    values$cmp_data <- cmp_results
    
    # Calculate performance for Chinook
    perf_data_chinook <- cmp_results %>%
      group_by(scenario) %>%
      summarise(chinook_raw = median(spawners), .groups = "drop")
    
    # Calculate weighted steelhead scores if steelhead_metrics exists
    if(exists("steelhead_metrics")) {
      steelhead_weighted <- map_dfr(input$cmp_scenarios, function(scen) {
        alts <- get_scenario_alternatives(scen, "all")
        alt_scores <- steelhead_metrics %>% filter(env %in% as.character(alts))
        alt_scores$hydro_year <- c("2011", "2014", "2017", "2020")
        weighted_score <- sum(alt_scores$steelhead_score * hydro_w[alt_scores$hydro_year])
        tibble(scenario = scen, steelhead_raw = weighted_score)
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
        Scenario = scenario,
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
    ggplot(df, aes(x = year, y = spawners, color = factor(scenario), group = scenario)) +
      geom_line(size = 1.2) +
      expand_limits(y = 0) +
      scale_y_continuous(labels = comma) +
      scale_color_viridis_d(name = "Scenario") +
      labs(title = "Comparison of Weighted Scenarios",
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
           x = "Scenario", y = "Forecasted Spawner Abundance") +
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
    req(values$performance_auto)
    
    hydro_df <- tibble(
      scenario = names(hardcoded_hydro_scores),
      hydro_raw = hardcoded_hydro_scores
    )
    
    values$performance_auto %>%
      left_join(hydro_df, by = "scenario") %>%
      mutate(hydro_raw = ifelse(is.na(hydro_raw), 50, hydro_raw)) %>%
      mutate(
        chinook_norm = normalize_scores_max(chinook_raw),
        steelhead_norm = normalize_scores_max(steelhead_raw),
        hydro_norm = normalize_scores_min(hydro_raw)  # Lower cost is better
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
        Scenario = scenario,
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
           x = "Scenario", y = "Normalized Score (0-1)") +
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
}

shinyApp(ui, server)