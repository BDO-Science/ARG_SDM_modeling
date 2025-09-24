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


ui <- navbarPage("Fall-run Chinook Power Bypass Simulator",
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
                                   h3("Summary"),
                                   p("This interactive application simulates fall-run Chinook salmon population dynamics on the American River under different power-bypass flow and temperature management alternatives. The model integrates temperature-dependent mortality, spawn timing, and life-cycle processes to project spawner abundance."),
                                   h4("Key Features:"),
                                   tags$ul(
                                     tags$li("Analyzes daily water temperature projections across multiple management alternatives"),
                                     tags$li("Uses 2011-2024 carcass survey data and CDFW escapement estimates for calibration"),
                                     tags$li("Models spawn timing using a cumulative link model based on October-November temperatures"),
                                     tags$li("Calculates temperature-dependent egg-to-fry survival using three different models"),
                                     tags$li("Incorporates adult pre-spawn mortality from thermal stress"),
                                     tags$li("Applies Beverton-Holt density dependence at the fry stage"),
                                     tags$li("Projects age-structured returns (ages 3-5) with calibrated smolt-to-adult return (SAR) rates"),
                                     tags$li("Supports dynamic weighting of both hydrological conditions and TDM models.")
                                   ),
                                   h3("Model Components"),
                                   h4("1. Spawn Timing Model"),
                                   p("A cumulative link model (CLM) predicts the probability of spawning in each 10-day bin based on standardized October and November water temperatures. The model divides the spawning season (October-January) into temporal bins, with spawn dates randomly sampled within bins according to temperature-driven probabilities. This 10-day resolution was chosen to balance model complexity with the temporal patterns observed in carcass survey data."),
                                   h4("2. Temperature-Dependent Mortality (TDM)"),
                                   p("Egg-to-fry survival is calculated using temperature-dependent mortality models applied from spawn date until emergence (1375 accumulated thermal units):"),
                                   tags$ul(
                                     tags$li(HTML("<strong>Exponential models</strong> (Water Forum 2020, SALMOD 2006): Daily hazard h(T) = α·exp(β·T) with stage-specific parameters for eggs (0-958 ATU) and alevins (958-1375 ATU)")),
                                     tags$li(HTML("<strong>Linear threshold model</strong> (Martin et al. 2017): Daily hazard h(T) = α·max(T-12.14, 0) applied above threshold temperature"))
                                   ),
                                   h4("3. Adult Pre-spawn Survival"),
                                   p("Adults accumulate thermal stress as degree-days from October 1 through spawning. Survival decreases logistically with accumulated degree-days:"),
                                   p(HTML("<code>S_pre = 1 / (1 + exp(-(3.0 - 0.00067 × degree_days)))</code>")),
                                   h4("4. Population Dynamics"),
                                   p("The life-cycle model tracks cohorts through multiple stages:"),
                                   tags$ul(
                                     tags$li("Female spawners (50% of total) produce 5,522 eggs each"),
                                     tags$li("Base egg survival S0 = 0.347 modified by density dependence"),
                                     tags$li("Carrying capacity K = 12,493 spawners (or flow-based if specified)"),
                                     tags$li("Age structure: 82.9% return at age 3, 16.9% at age 4, 0.2% at age 5 (2011-2024 CWT data)"),
                                     tags$li("Smolt-to-adult return (SAR) and rearing survival calibrated to minimize the sum of squared errors (SSE)")
                                   ),
                                   h3("Calibration & Validation"),
                                   p("The model jointly optimizes SAR_mean and rear_surv parameters by minimizing sum of squared errors between modeled spawners and observed escapement for 2014–2024 of the 2011–2024 calibration period. The years 2011–2013 serve as initial conditions."),
                                   h3("Key Equations"),
                                   tags$div(style = "background-color: #f5f5f5; padding: 15px; border-radius: 5px;",
                                            tags$ol(
                                              tags$li(HTML("<strong>Egg development:</strong> Days = 958 / T (°C)"), br(), tags$small(em("Time required for eggs to accumulate 958 thermal units to reach the alevin stage."))),
                                              tags$li(HTML("<strong>Alevin emergence:</strong> Days = 417 / T (°C)"), br(), tags$small(em("Time required for alevins to accumulate an additional 417 thermal units to emerge as fry."))),
                                              tags$li(HTML("<strong>Cumulative survival:</strong> S = exp(-Σ daily_hazard)"), br(), tags$small(em("Calculates the survival rate over a period by compounding daily mortality hazards."))),
                                              tags$li(HTML("<strong>Degree-days:</strong> DD = Σ max(T - 0, 0) from Oct 1 to spawn"), br(), tags$small(em("Accumulated thermal stress experienced by pre-spawn adults."))),
                                              tags$li(HTML("<strong>Redds:</strong> R = Spawners × 0.5 × S_pre"), br(), tags$small(em("Total number of nests (redds), assuming a 50% female population adjusted for pre-spawn survival."))),
                                              tags$li(HTML("<strong>Density dependence:</strong> dd = 0.347 / (1 + R / K)"), br(), tags$small(em("Adjusts egg survival based on the number of redds relative to habitat carrying capacity (K)."))),
                                              tags$li(HTML("<strong>Fry production:</strong> Fry = Eggs × TDM_survival × dd"), br(), tags$small(em("Total number of fry produced after accounting for temperature-dependent mortality and density effects."))),
                                              tags$li(HTML("<strong>Effective survival:</strong> S_eff = TDM_survival × dd"), br(), tags$small(em("The combined egg-to-fry survival rate from both temperature and density-dependent factors."))),
                                              tags$li(HTML("<strong>Returns at age a:</strong> S[y+a] += Smolts[y] × SAR[y] × P(age=a)"), br(), tags$small(em("Calculates the number of returning adult spawners in a future year based on smolt production, ocean survival (SAR), and age-at-return probability.")))
                                            )
                                   ),
                                   h3("References"),
                                   tags$ul(
                                     tags$li("Bartholow, J.M. & Heasley, J. 2006. Evaluation of Shasta Dam Scenarios Using a Salmon Production Model. USGS."),
                                     tags$li("Bratovich, P., Neal, M., Ransom, A., et al. 2020. Chinook Salmon Early Lifestage Survival and Folsom Dam Power Bypass Considerations. Water Forum Tech. Memo."),
                                     tags$li("Martin, B.T., et al. 2017. Phenomenological vs. biophysical models of thermal stress in aquatic eggs. Ecology Letters 20:50-59.")
                                   )
                            )
                          )
                 ),
                 
                 # ---- Calibration Tab ----
                 tabPanel("Calibration",
                          sidebarLayout(
                            sidebarPanel(
                              h4("TDM Model Weights"),
                              p("Adjust weights to see how the combined model fits historical data. Sliders will auto-adjust to sum to 1."),
                              sliderInput("calib_weight_wf", "Water Forum (Exponential)", value = 0.51, min = 0, max = 1, step = 0.01),
                              sliderInput("calib_weight_sm", "SALMOD (Exponential)", value = 0.24, min = 0, max = 1, step = 0.01),
                              sliderInput("calib_weight_martin", "Martin (Linear)", value = 0.25, min = 0, max = 1, step = 0.01),
                              actionButton("run_calib", "Run Calibration")
                            ),
                            mainPanel(
                              DTOutput("calib_table"),
                              plotOutput("calib_ts_plot")
                            )
                          )
                 ),
                 
                 # ---- Single Scenario Tab ----
                 tabPanel("Single Scenario Analysis",
                          sidebarLayout(
                            sidebarPanel(
                              h4("1. Management Alternative"),
                              selectInput("alternative", "Select Alternative:", choices = 1:5),
                              hr(),
                              h4("2. Hydrology Year Weights"),
                              sliderInput("weight_hydro_2017", "2017 Hydrology", value = 1.0, min = 0, max = 1, step = 0.05),
                              sliderInput("weight_hydro_2020", "2020 Hydrology", value = 0.0, min = 0, max = 1, step = 0.05),
                              hr(),
                              h4("3. TDM Model Weights"),
                              sliderInput("weight_wf", "Water Forum (Exponential)", value = 0.51, min = 0, max = 1, step = 0.01),
                              sliderInput("weight_sm", "SALMOD (Exponential)", value = 0.24, min = 0, max = 1, step = 0.01),
                              sliderInput("weight_martin", "Martin (Linear)", value = 0.25, min = 0, max = 1, step = 0.01),
                              hr(),
                              h4("4. Simulation Settings"),
                              sliderInput("sim_years", "Forecast Horizon (years)", min = 10, max = 100, value = 50, step = 1),
                              actionButton("run_sim", "Run Simulation"),
                              width = 3
                            ),
                            mainPanel(
                              tabsetPanel(
                                tabPanel("Table", DTOutput("results_table")),
                                tabPanel("Time Series", plotOutput("ts_plot")),
                                tabPanel("Distribution", plotOutput("dist_plot")),
                                tabPanel("Heatmap", plotOutput("heatmap_plot")),
                                tabPanel("Fry × DD", plotOutput("fry_dd_plot"))
                              ),
                              width = 9
                            )
                          )
                 ),
                 
                 # ---- Compare Scenarios Tab ----
                 tabPanel("Compare Scenarios",
                          sidebarLayout(
                            sidebarPanel(
                              h4("1. Management Alternatives"),
                              selectizeInput("alts_cmp", "Select alternatives to compare:", choices = 1:5, multiple = TRUE, selected = 1:5),
                              hr(),
                              h4("2. Hydrology Year Weights"),
                              sliderInput("cmp_weight_hydro_2017", "2017 Hydrology", value = 0.5, min = 0, max = 1, step = 0.05),
                              sliderInput("cmp_weight_hydro_2020", "2020 Hydrology", value = 0.5, min = 0, max = 1, step = 0.05),
                              hr(),
                              h4("3. TDM Model Weights"),
                              sliderInput("cmp_weight_wf", "Water Forum (Exponential)", value = 0.51, min = 0, max = 1, step = 0.01),
                              sliderInput("cmp_weight_sm", "SALMOD (Exponential)", value = 0.24, min = 0, max = 1, step = 0.01),
                              sliderInput("cmp_weight_martin", "Martin (Linear)", value = 0.25, min = 0, max = 1, step = 0.01),
                              hr(),
                              h4("4. Simulation Settings"),
                              sliderInput("cmp_years", "Forecast Horizon (years)", min = 10, max = 100, value = 50, step = 1),
                              actionButton("run_cmp", "Run Comparison"),
                              width = 3
                            ),
                            mainPanel(
                              tabsetPanel(
                                tabPanel("Table", DTOutput("cmp_table")),
                                tabPanel("Time Series", plotOutput("cmp_ts_plot")),
                                tabPanel("Boxplot (Last N Years)", 
                                         sliderInput("last_n", "Last N Years:", min = 5, max = 100, value = 10, step = 1),
                                         plotOutput("cmp_box_plot"),
                                         hr(),
                                         h4("Boxplot Summary"),
                                         tableOutput("cmp_boxplot_stats")
                                ),
                                tabPanel("Fry × DD Comparison", plotOutput("cmp_fry_dd_plot")),
                                tabPanel("Heatmap", plotOutput("cmp_heatmap"))
                              )
                            )
                          )
                 ),
                 
                 # ---- Objectives & Decision Support Tab ----
                 tabPanel("Objectives & Decision Support",
                          sidebarLayout(
                            sidebarPanel(
                              h4("Objective Weights"),
                              p("Adjust sliders to reflect the relative importance of each objective. Weights will always sum to 100%."),
                              sliderInput("w_chinook", "Fall-run Chinook Salmon", min = 0, max = 1, value = 0.4, step = 0.05),
                              sliderInput("w_steelhead", "American River Steelhead", min = 0, max = 1, value = 0.3, step = 0.05),
                              sliderInput("w_hydro", "Hydropower Generation", min = 0, max = 1, value = 0.15, step = 0.05),
                              sliderInput("w_hatchery", "Nimbus Hatchery Impact", min = 0, max = 1, value = 0.15, step = 0.05),
                              width = 3
                            ),
                            mainPanel(
                              h3("Performance Metrics"),
                              p("Enter performance scores for each alternative (0-100, higher is better). Chinook scores are auto-populated from the 'Compare Scenarios' results."),
                              uiOutput("performance_metric_inputs"),
                              hr(),
                              h3("Decision Support Plots"),
                              plotOutput("plot_total_score"),
                              plotOutput("plot_score_contribution")
                            )
                          )
                 )
)

server <- function(input, output, session) {
  
  # --- Reactive Values ---
  calib_data <- reactiveVal(NULL); sim_data <- reactiveVal(NULL); cmp_data <- reactiveVal(NULL)
  performance_scores <- reactiveVal(tibble(alternative = as.character(1:5), Chinook = rep(NA_real_, 5), Steelhead = rep(50, 5), Hydropower = rep(50, 5), Hatchery = rep(50, 5)))
  
  # --- Helper Functions ---
  run_weighted_simulation <- function(selected_alt, hydro_weights, tdm_weights, sim_params) {
    hydro_weights_norm <- normalize_weights(hydro_weights); tdm_weights_norm <- normalize_weights(tdm_weights)
    original_alt_names <- list(`2017` = as.character(selected_alt), `2020` = as.character(selected_alt + 5))
    hydro_results <- list()
    for (hydro_name in names(original_alt_names)) {
      alt_name <- original_alt_names[[hydro_name]]
      base_surv_vectors <- list(exp_WF = surv_lookup_full[[paste(alt_name, "exp_WF", sep = "_")]], exp_SM = surv_lookup_full[[paste(alt_name, "exp_SM", sep = "_")]], lin_Martin = surv_lookup_full[[paste(alt_name, "lin_Martin", sep = "_")]])
      weighted_surv_vec <- (base_surv_vectors$exp_WF*tdm_weights_norm["exp_WF"])+(base_surv_vectors$exp_SM*tdm_weights_norm["exp_SM"])+(base_surv_vectors$lin_Martin*tdm_weights_norm["lin_Martin"])
      P_weighted <- base_P
      for(param in c("SAR_mean","rear_surv")){P_weighted[[param]]<-(base_P_list$exp_WF[[alt_name]][[param]]*tdm_weights_norm["exp_WF"])+(base_P_list$exp_SM[[alt_name]][[param]]*tdm_weights_norm["exp_SM"])+(base_P_list$lin_Martin[[alt_name]][[param]]*tdm_weights_norm["lin_Martin"])}
      seed_weighted <- (S_seed_fore_list$exp_WF * tdm_weights_norm["exp_WF"]) + (S_seed_fore_list$exp_SM * tdm_weights_norm["exp_SM"]) + (S_seed_fore_list$lin_Martin * tdm_weights_norm["lin_Martin"])
      original_surv_lookup <- surv_lookup_full; original_base_P <- base_P_list
      temp_key <- paste(alt_name, "weighted_avg", sep = "_"); surv_lookup_full[[temp_key]] <<- weighted_surv_vec
      base_P_list[["weighted_avg"]] <<- list(); base_P_list[["weighted_avg"]][[alt_name]] <<- P_weighted
      fc_fn <- sim_forecast_fn(var_nm = "weighted_avg", env_nm = alt_name, flow_cfs = NULL, S_seed = seed_weighted, spawn_dates_by_alt = spawn_dates_by_alt)
      full_out <- fc_fn()
      surv_lookup_full <<- original_surv_lookup; base_P_list <<- original_base_P
      hydro_results[[hydro_name]] <- full_out %>% dplyr::filter(year > max(real_years)) %>% dplyr::slice_head(n = sim_params$yrs)
    }
    final_result <- hydro_results[[1]]; final_result$spawners <- final_result$spawners * hydro_weights_norm[1]
    if (length(hydro_results) > 1) { for (i in 2:length(hydro_results)) { final_result$spawners <- final_result$spawners + (hydro_results[[i]]$spawners * hydro_weights_norm[i]) } }
    return(final_result)
  }
  
  # --- Dynamic Slider Weighting Logic ---
  make_3_slider_observers <- function(id1, id2, id3, lock) {
    observeEvent(input[[id1]], { if(lock()) return(); lock(TRUE); remainder <- 1 - input[[id1]]; sum_others <- isolate(input[[id2]]) + isolate(input[[id3]]); if (sum_others > 0) { updateSliderInput(session, id2, value = remainder * (isolate(input[[id2]]) / sum_others)); updateSliderInput(session, id3, value = remainder * (isolate(input[[id3]]) / sum_others)) } else { updateSliderInput(session, id2, value = remainder / 2); updateSliderInput(session, id3, value = remainder / 2) }; lock(FALSE) }, ignoreInit = TRUE)
    observeEvent(input[[id2]], { if(lock()) return(); lock(TRUE); remainder <- 1 - input[[id2]]; sum_others <- isolate(input[[id1]]) + isolate(input[[id3]]); if (sum_others > 0) { updateSliderInput(session, id1, value = remainder * (isolate(input[[id1]]) / sum_others)); updateSliderInput(session, id3, value = remainder * (isolate(input[[id3]]) / sum_others)) } else { updateSliderInput(session, id1, value = remainder / 2); updateSliderInput(session, id3, value = remainder / 2) }; lock(FALSE) }, ignoreInit = TRUE)
    observeEvent(input[[id3]], { if(lock()) return(); lock(TRUE); remainder <- 1 - input[[id3]]; sum_others <- isolate(input[[id1]]) + isolate(input[[id2]]); if (sum_others > 0) { updateSliderInput(session, id1, value = remainder * (isolate(input[[id1]]) / sum_others)); updateSliderInput(session, id2, value = remainder * (isolate(input[[id2]]) / sum_others)) } else { updateSliderInput(session, id1, value = remainder / 2); updateSliderInput(session, id2, value = remainder / 2) }; lock(FALSE) }, ignoreInit = TRUE)
  }
  make_2_slider_observers <- function(id1, id2, lock) {
    observeEvent(input[[id1]], { if(lock()) return(); lock(TRUE); updateSliderInput(session, id2, value = 1 - input[[id1]]); lock(FALSE) }, ignoreInit = TRUE)
    observeEvent(input[[id2]], { if(lock()) return(); lock(TRUE); updateSliderInput(session, id1, value = 1 - input[[id2]]); lock(FALSE) }, ignoreInit = TRUE)
  }
  make_4_slider_observers <- function(id1, id2, id3, id4, lock) {
    observeEvent(input[[id1]], { if(lock()) return(); lock(TRUE); remainder <- 1 - input[[id1]]; sum_others <- isolate(input[[id2]])+isolate(input[[id3]])+isolate(input[[id4]]); if (sum_others > 0) { updateSliderInput(session,id2,value=remainder*(isolate(input[[id2]])/sum_others)); updateSliderInput(session,id3,value=remainder*(isolate(input[[id3]])/sum_others)); updateSliderInput(session,id4,value=remainder*(isolate(input[[id4]])/sum_others)) } else { updateSliderInput(session,id2,value=remainder/3); updateSliderInput(session,id3,value=remainder/3); updateSliderInput(session,id4,value=remainder/3) }; lock(FALSE) }, ignoreInit = TRUE)
    observeEvent(input[[id2]], { if(lock()) return(); lock(TRUE); remainder <- 1 - input[[id2]]; sum_others <- isolate(input[[id1]])+isolate(input[[id3]])+isolate(input[[id4]]); if (sum_others > 0) { updateSliderInput(session,id1,value=remainder*(isolate(input[[id1]])/sum_others)); updateSliderInput(session,id3,value=remainder*(isolate(input[[id3]])/sum_others)); updateSliderInput(session,id4,value=remainder*(isolate(input[[id4]])/sum_others)) } else { updateSliderInput(session,id1,value=remainder/3); updateSliderInput(session,id3,value=remainder/3); updateSliderInput(session,id4,value=remainder/3) }; lock(FALSE) }, ignoreInit = TRUE)
    observeEvent(input[[id3]], { if(lock()) return(); lock(TRUE); remainder <- 1 - input[[id3]]; sum_others <- isolate(input[[id1]])+isolate(input[[id2]])+isolate(input[[id4]]); if (sum_others > 0) { updateSliderInput(session,id1,value=remainder*(isolate(input[[id1]])/sum_others)); updateSliderInput(session,id2,value=remainder*(isolate(input[[id2]])/sum_others)); updateSliderInput(session,id4,value=remainder*(isolate(input[[id4]])/sum_others)) } else { updateSliderInput(session,id1,value=remainder/3); updateSliderInput(session,id2,value=remainder/3); updateSliderInput(session,id4,value=remainder/3) }; lock(FALSE) }, ignoreInit = TRUE)
    observeEvent(input[[id4]], { if(lock()) return(); lock(TRUE); remainder <- 1 - input[[id4]]; sum_others <- isolate(input[[id1]])+isolate(input[[id2]])+isolate(input[[id3]]); if (sum_others > 0) { updateSliderInput(session,id1,value=remainder*(isolate(input[[id1]])/sum_others)); updateSliderInput(session,id2,value=remainder*(isolate(input[[id2]])/sum_others)); updateSliderInput(session,id3,value=remainder*(isolate(input[[id3]])/sum_others)) } else { updateSliderInput(session,id1,value=remainder/3); updateSliderInput(session,id2,value=remainder/3); updateSliderInput(session,id3,value=remainder/3) }; lock(FALSE) }, ignoreInit = TRUE)
  }
  calib_tdm_lock <- reactiveVal(FALSE); make_3_slider_observers("calib_weight_wf", "calib_weight_sm", "calib_weight_martin", calib_tdm_lock)
  single_hydro_lock <- reactiveVal(FALSE); make_2_slider_observers("weight_hydro_2017", "weight_hydro_2020", single_hydro_lock)
  single_tdm_lock <- reactiveVal(FALSE); make_3_slider_observers("weight_wf", "weight_sm", "weight_martin", single_tdm_lock)
  cmp_hydro_lock <- reactiveVal(FALSE); make_2_slider_observers("cmp_weight_hydro_2017", "cmp_weight_hydro_2020", cmp_hydro_lock)
  cmp_tdm_lock <- reactiveVal(FALSE); make_3_slider_observers("cmp_weight_wf", "cmp_weight_sm", "cmp_weight_martin", cmp_tdm_lock)
  objective_lock <- reactiveVal(FALSE); make_4_slider_observers("w_chinook", "w_steelhead", "w_hydro", "w_hatchery", objective_lock)
  
  # --- Simulation Observers ---
  observeEvent(input$run_calib, {
    tdm_weights <- c(exp_WF = input$calib_weight_wf, exp_SM = input$calib_weight_sm, lin_Martin = input$calib_weight_martin); tdm_weights_norm <- normalize_weights(tdm_weights)
    ref_env <- "1" 
    
    # FIX: Explicitly calculate degree days for the full calibration period (2011-2024).
    deg_day_cal <- compute_deg_day_adult(
      env_nm = ref_env,
      sim_years = real_years,
      spawn_dates = spawn_dates_by_alt[[ref_env]][match(real_years, sim_years)],
      env_ext_list = env_ext_list
    )
    
    calib_surv_vec <- (surv_lookup_full[[paste(ref_env, "exp_WF", sep="_")]]*tdm_weights_norm["exp_WF"])+(surv_lookup_full[[paste(ref_env, "exp_SM", sep="_")]]*tdm_weights_norm["exp_SM"])+(surv_lookup_full[[paste(ref_env, "lin_Martin", sep="_")]]*tdm_weights_norm["lin_Martin"])
    P_calib <- base_P
    for(param in c("SAR_mean","rear_surv")){P_calib[[param]]<-(base_P_list$exp_WF[[ref_env]][[param]]*tdm_weights_norm["exp_WF"])+(base_P_list$exp_SM[[ref_env]][[param]]*tdm_weights_norm["exp_SM"])+(base_P_list$lin_Martin[[ref_env]][[param]]*tdm_weights_norm["lin_Martin"])}
    
    out <- simulate_variant(surv_vec=calib_surv_vec[1:n_calib],P=P_calib,years=n_calib,S_init=S_seed_calib,
                            SAR_vec=rep(P_calib$SAR_mean,n_calib),K_spawners_vec=rep(P_calib$K_spawners,n_calib),
                            deg_day_adult=deg_day_cal, sim_years_vec=real_years)
    
    calib_data(tibble(year=real_years,observed=esc_obs$spawners,predicted=out$spawners))
  })
  observeEvent(input$run_sim, {
    shiny::showNotification("Running simulation…",id="sim_busy",duration=NULL)
    tdm_weights<-c(exp_WF=input$weight_wf,exp_SM=input$weight_sm,lin_Martin=input$weight_martin)
    hydro_weights<-c(`2017`=input$weight_hydro_2017,`2020`=input$weight_hydro_2020)
    sim_params<-list(yrs=input$sim_years)
    result<-run_weighted_simulation(as.numeric(input$alternative),hydro_weights,tdm_weights,sim_params)
    sim_data(result%>%mutate(alternative=input$alternative))
    shiny::removeNotification("sim_busy");shiny::showNotification("Simulation complete!",type="message",duration=2)
  })
  observeEvent(input$run_cmp, {
    req(input$alts_cmp);shiny::showNotification("Running comparison…",id="cmp_busy",duration=NULL)
    tdm_weights<-c(exp_WF=input$cmp_weight_wf,exp_SM=input$cmp_weight_sm,lin_Martin=input$cmp_weight_martin)
    hydro_weights<-c(`2017`=input$cmp_weight_hydro_2017,`2020`=input$cmp_weight_hydro_2020)
    sim_params<-list(yrs=input$cmp_years)
    results_list<-purrr::map(input$alts_cmp,function(alt){res<-run_weighted_simulation(as.numeric(alt),hydro_weights,tdm_weights,sim_params);res%>%mutate(alternative=alt)})
    comparison_results<-bind_rows(results_list);cmp_data(comparison_results)
    if(nrow(comparison_results)>0){
      median_summary<-comparison_results%>%group_by(alternative)%>%summarise(Median=median(spawners,na.rm=TRUE),.groups="drop")
      chinook_scores<-normalize_scores(median_summary$Median)
      current_scores<-performance_scores();current_scores$Chinook<-NA_real_
      match_idx<-match(median_summary$alternative,current_scores$alternative)
      current_scores$Chinook[match_idx]<-chinook_scores
      performance_scores(current_scores)
    }
    shiny::removeNotification("cmp_busy");shiny::showNotification("Comparison complete! Check the 'Objectives' tab.",type="message",duration=4)
  })
  
  # --- Dynamic UI for Performance Metrics ---
  output$performance_metric_inputs <- renderUI({
    metrics <- c("Steelhead", "Hydropower", "Hatchery"); alts <- 1:5; scores <- performance_scores()
    fluidPage(
      fluidRow(column(3, strong("Chinook (auto)")), lapply(alts, function(alt) { column(2, p(strong(round(scores$Chinook[scores$alternative == alt], 1)))) })),
      lapply(metrics, function(metric) {
        fluidRow(
          column(3, strong(metric)),
          lapply(alts, function(alt) {
            column(2, numericInput(inputId=paste0("pm_",tolower(metric),"_",alt),label=paste("Alt",alt),value=scores[[metric]][alt],min=0,max=100))
          })
        )
      })
    )
  })
  
  # --- Reactive for Final Weighted Scores ---
  final_scores <- reactive({
    scores <- performance_scores()
    for (alt in 1:5) { 
      if(!is.null(input[[paste0("pm_steelhead_", alt)]])) scores$Steelhead[alt]<-input[[paste0("pm_steelhead_",alt)]]
      if(!is.null(input[[paste0("pm_hydropower_",alt)]])) scores$Hydropower[alt]<-input[[paste0("pm_hydropower_",alt)]]
      if(!is.null(input[[paste0("pm_hatchery_",alt)]])) scores$Hatchery[alt]<-input[[paste0("pm_hatchery_",alt)]]
    }
    weights <- c(Chinook=input$w_chinook,Steelhead=input$w_steelhead,Hydropower=input$w_hydro,Hatchery=input$w_hatchery)
    scores %>% mutate(across(c(Chinook,Steelhead,Hydropower,Hatchery),~replace_na(.x,0)))%>%
      mutate(w_Chinook=Chinook*weights["Chinook"],w_Steelhead=Steelhead*weights["Steelhead"],w_Hydropower=Hydropower*weights["Hydropower"],w_Hatchery=Hatchery*weights["Hatchery"],
             TotalScore=w_Chinook+w_Steelhead+w_Hydropower+w_Hatchery)
  })
  
  # --- Decision Support Plot Outputs ---
  output$plot_total_score <- renderPlot({
    df <- final_scores()
    ggplot(df, aes(x = alternative, y = TotalScore, fill = alternative)) + geom_col() +
      scale_fill_viridis_d(guide = "none") + labs(title = "Overall Performance Score by Alternative", x = "Management Alternative", y = "Total Weighted Score") + theme_minimal(base_size = 16)
  })
  output$plot_score_contribution <- renderPlot({
    df <- final_scores() %>% pivot_longer(cols = starts_with("w_"), names_to = "objective", values_to = "contribution") %>% mutate(objective = str_remove(objective, "w_"))
    ggplot(df, aes(x = alternative, y = contribution, fill = objective)) + geom_col(position = "stack") +
      scale_fill_viridis_d(name = "Objective") + labs(title = "Contribution to Score by Objective", x = "Management Alternative", y = "Weighted Score Contribution") + theme_minimal(base_size = 16) + theme(legend.position = "bottom")
  })
  
  # --- All Other Outputs ---
  output$calib_table<-renderDT({
    req(calib_data())%>%
      mutate(across(where(is.numeric),round,0)) %>%
      datatable(
        rownames=FALSE,
        options=list(
          dom='t', # 't' stands for table, removes search/paging controls
          paging = FALSE # Explicitly disable pagination to show all rows
        )
      )
  })
  output$calib_ts_plot<-renderPlot({df<-req(calib_data());ggplot(df,aes(x=year))+geom_line(aes(y=observed,color="Observed"),linewidth=1)+geom_line(aes(y=predicted,color="Predicted"),linewidth=1)+scale_color_manual(name="Series",values=c(Observed="black",Predicted="steelblue"))+scale_y_continuous(labels=comma)+labs(title="Weighted Calibration: Observed vs Predicted (2011–2024)",x="Year",y="Spawner Abundance")+theme_minimal(base_size=16)+theme(legend.position="bottom")})
  output$results_table<-renderDT({req(sim_data())%>%mutate(spawners=round(spawners))%>%dplyr::select(Year=year,`Forecasted Spawners`=spawners)%>%datatable(rownames=FALSE,extensions='Buttons',options=list(dom='Bfrtip',buttons=c('csv')))})
  output$ts_plot<-renderPlot({df<-req(sim_data());ggplot(df,aes(x=year,y=spawners))+geom_line(color="steelblue",size=1.2)+expand_limits(y=0)+scale_y_continuous(labels=scales::comma)+labs(title=paste0("Weighted Forecast for Alternative ",unique(df$alternative)),x="Year",y="Forecasted Spawner Abundance")+theme_minimal(base_size=16)})
  output$dist_plot<-renderPlot({df<-req(sim_data());ggplot(df,aes(x=spawners))+geom_histogram(bins=30,fill="steelblue",alpha=0.8)+scale_x_continuous(labels=scales::comma)+labs(title=paste0("Distribution of Forecasted Spawners for Alt ",unique(df$alternative)),x="Forecasted Spawner Abundance",y="Frequency")+theme_minimal(base_size=16)})
  output$heatmap_plot<-renderPlot({df<-req(sim_data());ggplot(df,aes(x=year,y=alternative,fill=spawners))+geom_tile()+scale_fill_viridis_c(name="Spawners",labels=comma)+labs(title="Forecast Heatmap",x="Year",y="Alternative")+theme_minimal(base_size=14)})
  output$fry_dd_plot<-renderPlot({df<-req(sim_data());ggplot(df,aes(x=year,y=fry_dd))+geom_line()+expand_limits(y=0)+scale_y_continuous(labels=label_number(scale=1e-6,suffix=" M"))+labs(title=paste0("Fry × DD — Alt ",unique(df$alternative)),x="Year",y="Fry after TDM & Density‐dep")+theme_minimal(base_size=14)})
  output$cmp_table<-renderDT({req(cmp_data())%>%mutate(spawners=round(spawners))%>%dplyr::select(Alternative=alternative,Year=year,`Forecasted Spawners`=spawners)%>%datatable(rownames=FALSE,extensions='Buttons',options=list(dom='Bfrtip',buttons=c('csv')))})
  output$cmp_ts_plot<-renderPlot({df<-req(cmp_data());ggplot(df,aes(x=year,y=spawners,color=factor(alternative),group=alternative))+geom_line(size=1.2)+expand_limits(y=0)+scale_y_continuous(labels=scales::comma)+scale_color_viridis_d(name="Alternative")+labs(title="Comparison of Weighted Scenarios",x="Year",y="Forecasted Spawner Abundance")+theme_minimal(base_size=16)+theme(legend.position="bottom")})
  output$cmp_box_plot<-renderPlot({req(cmp_data());last_yr<-max(cmp_data()$year);df_filt<-cmp_data()%>%filter(year>=(last_yr-input$last_n+1));ggplot(df_filt,aes(x=factor(alternative),y=spawners,fill=factor(alternative)))+geom_boxplot()+scale_fill_viridis_d(guide="none")+scale_y_continuous(labels=scales::comma)+labs(title=paste0("Spawner Distribution: Last ",input$last_n," Years"),x="Alternative",y="Forecasted Spawner Abundance")+theme_minimal(base_size=16)})
  output$cmp_boxplot_stats<-renderTable({req(cmp_data());last_yr<-max(cmp_data()$year);cmp_data()%>%filter(year>=(last_yr-input$last_n+1))%>%group_by(alternative)%>%summarise(Minimum=min(spawners,na.rm=TRUE),`1st Qu.`=quantile(spawners,0.25,na.rm=TRUE),Median=median(spawners,na.rm=TRUE),`3rd Qu.`=quantile(spawners,0.75,na.rm=TRUE),Maximum=max(spawners,na.rm=TRUE),.groups="drop")%>%mutate(across(where(is.numeric),round,0))},rownames=FALSE)
  output$cmp_heatmap<-renderPlot({df<-req(cmp_data());ggplot(df,aes(year,alternative,fill=spawners))+geom_tile()+scale_fill_viridis_c(name="Spawners",labels=comma)+labs(title="Comparison Heatmap",x="Year",y="Alternative")+theme_minimal(base_size=14)+theme(axis.text.x=element_text(angle=90,vjust=0.5))})
  output$cmp_fry_dd_plot<-renderPlot({df<-req(cmp_data());ggplot(df,aes(x=year,y=fry_dd,color=alternative))+geom_line()+expand_limits(y=0)+scale_color_viridis_d(name="Alternative")+scale_y_continuous(labels=label_number(scale=1e-6,suffix=" M"))+labs(title="Fry × DD Comparison",x="Year",y="Fry after TDM & Density‐dep")+theme_minimal(base_size=14)})
  
}

shinyApp(ui, server)


