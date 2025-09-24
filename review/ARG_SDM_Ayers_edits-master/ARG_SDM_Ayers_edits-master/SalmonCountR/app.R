library(shiny)
library(tidyverse)
library(DT)

source("global.R")

ui <- navbarPage("American River Fall-run Chinook Salmon Power Bypass Alternative Simulator",
                 
                 # ---- Info tab -----------------------------------------------------------
                 tabPanel("About",
                          fluidRow(
                            column(12,
                                   h3("Author & Contact Information"),
                                   tags$ul(
                                     tags$li(strong("Author:"), "Alexander Vaisvil"),
                                     tags$li(strong("Email:"), tags$a(href = "mailto:avaisvil@usbr.gov", "avaisvil@usbr.gov")),
                                     tags$li(strong("Institution/Organization:"), "U.S. Bureau of Reclamation"),
                                     tags$li(strong("GitHub Repository:"), tags$a(href = "https://github.com/BDO-Science/ARG_SDM_modeling", "https://github.com/BDO-Science/ARG_SDM_modeling")),
                                     tags$li(strong("Date Last Updated:"), format(Sys.Date(), "%B %d, %Y"))
                                   ),
                                   hr(),
                                   h3("Summary"),
                                   p("This interactive Shiny app lets you explore how different power-bypass flow and temperature alternatives affect fall-run Chinook salmon modeled spawner abundance on the American River. You’ll:"),
                                   tags$ul(
                                     tags$li("Load daily temperature & flow data and redd survey dates (2011–2024)."),
                                     tags$li("Model spawn timing across 10-day bins from October–November temperatures."),
                                     tags$li("Compute egg development and fry emergence times from river temperatures."),
                                     tags$li("Apply temperature-dependent mortality (TDM) to estimate egg→fry survival."),
                                     tags$li("Model pre-spawn survival of adults from cumulative degree-days (°C·days)."),
                                     tags$li("Apply density-dependence to fry (fry→smolt)."),
                                     tags$li("Recruit smolts to modeled spawners via a calibrated smolt-to-adult ratio (SAR) and fixed age lags."),
                                     tags$li("Calibrate SAR_mean and rear_surv against observed escapement (2011–2024)."),
                                     tags$li("Forecast up to 100 years deterministically or with stochastic SAR (Normal/Lognormal/Beta/Gamma; full, block, or pulse timing)."),
                                     tags$li("Explore results with tables and plots: time series, distributions, fry×DD, egg→fry survival, SAR series, heatmaps, and comparisons.")
                                   ),
                                   
                                   
                                   h3("Details"),
                                   tags$ul(
                                     tags$li(strong("Spawn timing (redd dates):"),
                                             " Ordered multinomial (cumulative logit) over 10-day bins with fixed effects ",
                                             HTML("<code>Oct_std</code> and <code>Nov_std</code>"),
                                             " and a brood-year random intercept. Forecasts simulate the year effect with an AR(1) process ",
                                             HTML("(φ, σ)"),
                                             "; flow is not used as a covariate. Simulated spawn dates feed degree-day and TDM windows."),
                                     tags$li(strong("Egg development & fry emergence:"),
                                             HTML("<code>hatch_model(T) = 958/T</code> (days to hatch), "),
                                             HTML("<code>emergence_model(T) = 417/T</code> (days to emergence).")),
                                     tags$li(strong("Temperature-Dependent Mortality (TDM):"),
                                             " evaluated over the incubation window starting on each redd date.",
                                             tags$ul(
                                               tags$li(HTML("<em>Exponential (WF2020 / SALMOD2006):</em> daily hazard "),
                                                       HTML("<code>h_i = α · exp(β·T_i)</code>; cumulative survival "),
                                                       HTML("<code>S = exp(-∑ h_i)</code>. Stage-specific (egg & alevin) parameter sets.")),
                                               tags$li(HTML("<em>Linear-threshold (Martin 2017):</em> daily hazard "),
                                                       HTML("<code>h_i = α · max(T_i − β, 0)</code>; "),
                                                       HTML("<code>S = exp(-∑ h_i)</code>."))
                                             )),
                                     tags$li(strong("Adult pre-spawn survival:"),
                                             HTML("<code>S_pre = plogis(θ₀ + θ₁·deg_day)</code>"),
                                             ", where ",
                                             HTML("<code>deg_day = Σ max(T − base_temp, 0)</code> from Oct 1 (default) to spawn.")),
                                     tags$li(strong("Density dependence (fry):"),
                                             HTML("<code>dd = S0 / (1 + redds / K_spawners)</code>"), 
                                             ". ", HTML("<code>K_spawners</code>"),
                                             " may be set from flow (WUA ÷ redd size) if provided."),
                                     tags$li(strong("Life-cycle recursion:"),
                                             HTML("modeled spawners are updated at ages 3–5: "),
                                             HTML("<code>S[y+a] ← S[y+a] + reared[y]·SAR[y]·P(age=a)</code>"),
                                             ", with age lag probabilities ", HTML("<code>{0.75, 0.249, 0.001}</code>.")),
                                     tags$li(strong("Calibration:"),
                                             " jointly fits ", HTML("<code>SAR_mean</code>"), " and ", HTML("<code>rear_surv</code>"),
                                             " by minimizing SSE on years 4–14 of 2011–2024 (modeled spawners vs. observed escapement)."),
                                     tags$li(strong("Forecasting:"),
                                             " produces modeled spawner abundance trajectories under alternative flow/temperature scenarios.")
                                   ),
                                   
                                   h3("Key Equations"),
                                   p("All equations below generate survival or fecundity terms that contribute to modeled spawner abundance (not direct escapement counts)."),
                                   tags$ul(
                                     tags$li(HTML("1. Egg development days: <code>D<sub>egg</sub> = 958 / T</code>")),
                                     tags$li(HTML("2. Fry emergence days: <code>D<sub>fry</sub> = 417 / T</code>")),
                                     tags$li(HTML("3. Exponential TDM (daily): <code>h_i = α·e^{βT_i}</code>; cumulative: <code>S = exp(-∑ h_i)</code>")),
                                     tags$li(HTML("4. Martin TDM (daily): <code>h_i = α·max(T_i - β, 0)</code>; cumulative: <code>S = exp(-∑ h_i)</code>")),
                                     tags$li(HTML("5. Pre-spawn survival: <code>S<sub>pre</sub> = plogis(θ₀ + θ₁·deg_day)</code>")),
                                     tags$li(HTML("6. Degree-days: <code>deg_day = Σ max(T − base_temp, 0)</code>")),
                                     tags$li(HTML("7. Redds: <code>redds = S · female_fraction · S<sub>pre</sub></code>")),
                                     tags$li(HTML("8. Density dependence: <code>dd = S₀ / (1 + redds / K<sub>spawners</sub>)</code>")),
                                     tags$li(HTML("9. Fry after TDM & DD: <code>fry_dd = eggs · S<sub>egg</sub> · dd</code>")),
                                     tags$li(HTML("10. Effective survival: <code>eff_surv = S<sub>egg</sub> · dd</code>")),
                                     tags$li(HTML("11. Age-structured returns (a = 3,4,5): <code>S[y+a] ← S[y+a] + reared[y]·SAR[y]·P(age=a)</code>"))
                                   )
                                   ,
                                   h4("References"),
                                   tags$ul(
                                     tags$li("Bartholow, J.M. & Heasley, J. 2006. Evaluation of Shasta Dam Scenarios Using a Salmon Production Model. USGS."),
                                     tags$li("Bratovich, P., Neal, M., Ransom, A., et al. 2020. Chinook Salmon Early Lifestage Survival and Folsom Dam Power Bypass Considerations. Water Forum Tech. Memo."),
                                     tags$li("CFS 2010. A Revised Sacramento River Winter Chinook Salmon Juvenile Production Model. Cramer Fish Sciences."),
                                     tags$li("HCI 1996. Chinook Salmon Mortality Model: Development, Evaluation, and Application. Hydrologic Consultants Inc."),
                                     tags$li("\"Heat-stress & pre-spawn survival,\" Wiley 2021. https://doi.org/10.1111/rec.13244"),
                                     tags$li("Martin, B.T., et al. 2017. Phenomenological vs. biophysical models of thermal stress in aquatic eggs. Ecol Lett 20:50–59."),
                                     tags$li("USBR 2008. Biological Assessment on Continued Operations of CVP & SWP. US Bureau of Reclamation."),
                                     tags$li("Zeug, S., Bergman, P., Cavallo, B., & Jones, K. 2012. Application of a Life Cycle Simulation Model. Environ Model Assess.")
                                   )
                            )
                          )
                 ),
                 
                 # ---- Calibration ---------------------------------------------------------
                 tabPanel("Calibration",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("calib_variant", "TDM Variant:", choices = names(base_P_list)),
                              actionButton("run_calib",    "Run Calibration")
                            ),
                            mainPanel(
                              DTOutput("calib_table"),
                              plotOutput("calib_ts_plot")
                            )
                          )
                 ),
                 
                 # ---- Single Alternative ---------------------------------------------------
                 tabPanel("Single Alternative",
                          sidebarLayout(
                            sidebarPanel(
                              h4("Simulation Settings"),
                              sliderInput("sim_years", "Forecast Horizon (years)",
                                          min = 10, max = 100, value = 50, step = 1
                              ),
                              selectInput("mode", "Mode:", choices = c("Deterministic", "Stochastic")),
                              conditionalPanel(
                                condition = "input.mode == 'Stochastic'",
                                selectInput("stoch_type", "Distribution Type:",
                                            choices = c("Normal", "Lognormal", "Beta", "Gamma")
                                ),
                                sliderInput("stoch_sd", "SAR SD:", min = 0, max = 0.01,
                                            value = stoch_SAR_opts$sd, step = 0.0001
                                ),
                                conditionalPanel(
                                  condition = "input.stoch_type == 'Beta' || input.stoch_type == 'Gamma'",
                                  numericInput("shape1", "Shape1:", value = stoch_SAR_opts$shape1, min = 0.1, step = 0.1),
                                  numericInput("shape2", "Shape2:", value = stoch_SAR_opts$shape2, min = 0.1, step = 0.1)
                                ),
                                selectInput("stoch_timing", "Timing:", choices = c("all", "block", "pulse"),
                                            selected = stoch_SAR_opts$timing
                                ),
                                conditionalPanel(
                                  condition = "input.stoch_timing == 'block'",
                                  sliderInput("block_start", "Block Start Year:", min = 1, max = 100,
                                              value = min(stoch_SAR_opts$block_years), step = 1
                                  ),
                                  sliderInput("block_end",   "Block End Year:",   min = 1, max = 100,
                                              value = max(stoch_SAR_opts$block_years),   step = 1
                                  )
                                ),
                                conditionalPanel(
                                  condition = "input.stoch_timing == 'pulse'",
                                  selectizeInput("pulse_years", "Pulse Years:",
                                                 choices = 1:100, multiple = TRUE, selected = stoch_SAR_opts$pulse_years
                                  )
                                )
                              ),
                              selectInput("tdm_variant", "TDM Variant:", choices = names(base_P_list)),
                              selectInput("alternative", "Alternative:", choices = unique(egg_summary$env)),
                              numericInput("flow_cfs", "Flow (cfs):",
                                           value = 1500, min = min(instream$flow_cfs), max = max(instream$flow_cfs), step = 50
                              ),
                              actionButton("run_sim", "Run Simulation"),
                              width = 3
                            ),
                            mainPanel(
                              tabsetPanel(
                                tabPanel("Table", DTOutput("results_table")),
                                tabPanel("Time Series", plotOutput("ts_plot")),
                                tabPanel("Distribution", plotOutput("dist_plot")),
                                tabPanel("Heatmap", plotOutput("heatmap_plot")),
                                tabPanel("Fry × DD", plotOutput("fry_dd_plot")),
                                #tabPanel("SAR Series", plotOutput("sar_plot"))
                              ),
                              width = 9
                            )
                          )
                 ),
                 
                 # ---- Compare Alternatives -------------------------------------------------
                 tabPanel("Compare Alternatives",
                          sidebarLayout(
                            sidebarPanel(
                              h4("Comparison Settings"),
                              sliderInput("cmp_years", "Forecast Horizon (years)",
                                          min = 10, max = 100, value = 50, step = 1
                              ),
                              selectInput("cmp_mode", "Mode:", choices = c("Deterministic", "Stochastic")),
                              conditionalPanel(
                                condition = "input.cmp_mode == 'Stochastic'",
                                selectInput("cmp_stoch_type", "Distribution Type:",
                                            choices = c("Normal", "Lognormal", "Beta", "Gamma")
                                ),
                                sliderInput("cmp_sd", "SAR SD:", min = 0, max = 0.01,
                                            value = stoch_SAR_opts$sd, step = 0.0001
                                ),
                                conditionalPanel(
                                  condition = "input.cmp_stoch_type == 'Beta' || input.cmp_stoch_type == 'Gamma'",
                                  numericInput("cmp_shape1", "Shape1:", value = stoch_SAR_opts$shape1, min = 0.1, step = 0.1),
                                  numericInput("cmp_shape2", "Shape2:", value = stoch_SAR_opts$shape2, min = 0.1, step = 0.1)
                                ),
                                selectInput("cmp_timing", "Timing:", choices = c("all", "block", "pulse"),
                                            selected = stoch_SAR_opts$timing
                                ),
                                conditionalPanel(
                                  condition = "input.cmp_timing == 'block'",
                                  sliderInput("cmp_block_start", "Block Start Year:", min = 1, max = 100,
                                              value = min(stoch_SAR_opts$block_years), step = 1
                                  ),
                                  sliderInput("cmp_block_end",   "Block End Year:",   min = 1, max = 100,
                                              value = max(stoch_SAR_opts$block_years),   step = 1
                                  )
                                ),
                                conditionalPanel(
                                  condition = "input.cmp_timing == 'pulse'",
                                  selectizeInput("cmp_pulse_years", "Pulse Years:",
                                                 choices = 1:100, multiple = TRUE, selected = stoch_SAR_opts$pulse_years
                                  )
                                )
                              ),
                              selectInput("cmp_variant", "TDM Variant:", choices = names(base_P_list)),
                              selectInput("alts_cmp", "Alternatives to compare:",
                                          choices = unique(egg_summary$env), multiple = TRUE,
                                          selected = unique(egg_summary$env)
                              ),
                              numericInput("flow_cfs", "Flow (cfs):",
                                           value = 1500, min = min(instream$flow_cfs), max = max(instream$flow_cfs), step = 50
                              ),
                              actionButton("run_cmp", "Run Comparison"),
                              width = 3
                            ),
                            mainPanel(
                              tabsetPanel(
                                tabPanel("Table",    DTOutput("cmp_table")),
                                tabPanel("Time Series", plotOutput("cmp_ts_plot")),
                                tabPanel("Boxplot (Last N Years)",
                                         sliderInput("last_n", "Last N Years:",
                                                     min = 5, max = 100, value = 10, step = 1
                                         ),
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
                 )
)

server <- function(input, output, session) {
  calib_end  <- length(real_years)   # 14
  sim_data <- reactiveVal(NULL)
  cmp_data <- reactiveVal(NULL)
  
  # Helper to surface errors in the UI and console
  run_safe <- function(expr, tag = "task") {
    tryCatch(
      force(expr),
      error = function(e) {
        msg <- paste0("[", tag, "] ", conditionMessage(e))
        cat(msg, "\n")
        shiny::showNotification(msg, type = "error", duration = NULL)
        stop(e)
      }
    )
  }
  
  options(shiny.fullstacktrace = TRUE)
  
  calib_data <- eventReactive(input$run_calib, {
    v <- req(input$calib_variant)
    df <- calib_pred_by_variant[[v]]
    validate(need(!is.null(df), "No precomputed calibration for this variant."))
    df %>% mutate(variant = v)   # optional, handy for labeling
  })
  
  ## put this ONCE in server(), near the top, after calib_pred_by_variant is loaded:
  ymax_global <- max(
    unlist(lapply(calib_pred_by_variant, function(df) c(df$observed, df$predicted))),
    na.rm = TRUE
  ) * 1.1
  
  ## ---- Calibration table (single definition) ----
  output$calib_table <- renderDT({
    df <- req(calib_data())
    df %>%
      mutate(
        observed  = round(observed,  0),
        predicted = round(predicted, 0),
        SAR_mean  = signif(SAR_mean, 4),
        rear_surv = if ("rear_surv" %in% names(df)) signif(rear_surv, 3) else NA
      ) %>%
      datatable(rownames = FALSE, extensions = "Buttons",
                options = list(dom = "Bfrtip", buttons = c("csv")))
  })
  
  ## ---- Calibration plot (single definition) ----
  output$calib_ts_plot <- renderPlot({
    df <- req(calib_data())         # <- use the reactive you already made
    v  <- unique(df$variant)
    
    ggplot(df, aes(x = year)) +
      geom_line(aes(y = observed,  color = "Observed"),  linewidth = 1) +
      geom_line(aes(y = predicted, color = "Predicted"), linewidth = 1) +
      scale_color_manual(name = "Series",
                         values = c("Observed" = "black", "Predicted" = "steelblue")) +
      labs(
        title = paste0("Calibration: Observed vs Predicted (2011–2024) — ", v),
        x = "Year", y = "Spawner Abundance"
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "bottom") +
      scale_y_continuous(
        expand = expansion(mult = c(0, 0.05)),   # forces baseline at 0
        limits = c(0, ymax_global)               # shared max for all variants
      )
  })
  
  # ---- SINGLE-ALT SIMULATION ----
  observeEvent(input$run_sim, {
    cat("run_sim clicked\n")
    shiny::showNotification("Running simulation…", id = "sim_busy", duration = NULL)
    withProgress(message = "Simulating spawners", value = 0, {
      run_safe({
        yrs     <- input$sim_years
        alt     <- input$alternative
        variant <- input$tdm_variant
        seed_vec <- S_seed_fore_list[[variant]]
        incProgress(0.1)
        
        # If your simulate_variant() uses a global n_sim, set it; if not, ignore this.
        if (exists("sim_years", inherits = TRUE)) {
          n_sim <<- length(sim_years)  # or n_sim <<- yrs if that's what simulate_variant expects
        }
        
        # Stochastic SAR controls
        use_stochastic_SAR <<- (input$mode == "Stochastic")
        P_sim <- base_P_list[[variant]][[alt]]
        stoch_SAR_opts$mean <<- P_sim$SAR_mean
        stoch_SAR_opts$model <<- tolower(input$stoch_type)
        stoch_SAR_opts$sd    <<- input$stoch_sd
        if (stoch_SAR_opts$model %in% c("beta","gamma")) {
          stoch_SAR_opts$shape1 <<- input$shape1
          stoch_SAR_opts$shape2 <<- input$shape2
        }
        stoch_SAR_opts$timing <<- input$stoch_timing
        if (stoch_SAR_opts$timing == "block") {
          stoch_SAR_opts$block_years <<- seq(input$block_start, input$block_end)
        }
        if (stoch_SAR_opts$timing == "pulse") {
          stoch_SAR_opts$pulse_years <<- input$pulse_years
          stoch_SAR_opts$pulse_sd    <<- input$stoch_sd
        }
        incProgress(0.2)
        
        # Build & run forecast (no 'surv_method' arg anymore)
        fc_fn <- sim_forecast_fn(
          var_nm   = variant,
          env_nm   = alt,
          flow_cfs = input$flow_cfs,
          S_seed   = seed_vec,
          spawn_dates_by_env = spawn_dates_by_env  # <-- add
        )
        full_out <- fc_fn()
        incProgress(0.4)
        
        # Post-process to requested horizon
        out <- full_out %>%
          dplyr::filter(year > max(real_years)) %>%
          dplyr::slice_head(n = yrs) %>%
          dplyr::mutate(
            spawners  = round(spawners),
            deg_day   = round(deg_day, 0),
            pre_spawn = round(pre_spawn, 3),
            rear_surv = round(rear_surv, 3),
            SAR_used  = round(SAR_used, 4),
            dd        = round(dd, 2),
            fry_dd    = round(fry_dd, 2),
            egg_surv  = round(egg_surv, 2),
            eff_surv  = round(eff_surv, 2),
            alternative = alt,
            variant     = variant
          )
        
        sim_data(out)
        incProgress(0.3)
      }, tag = "single-sim")
    })
    shiny::removeNotification("sim_busy")
    shiny::showNotification("Simulation complete!", type = "message", duration = 2)
  })
  
  # ---- COMPARE-ALTS SIMULATION ----
  observeEvent(input$run_cmp, {
    cat("run_cmp clicked\n")
    shiny::showNotification("Running comparison…", id = "cmp_busy", duration = NULL)
    withProgress(message = "Comparing alternatives", value = 0, {
      run_safe({
        yrs     <- input$cmp_years
        variant <- input$cmp_variant
        
        # Stochastic SAR controls
        use_stochastic_SAR <<- (input$cmp_mode == "Stochastic")
        stoch_SAR_opts$model  <<- tolower(input$cmp_stoch_type)
        # mean per-alt derived below from P for that alt when we run the loop
        stoch_SAR_opts$sd     <<- input$cmp_sd
        if (stoch_SAR_opts$model %in% c("beta","gamma")) {
          stoch_SAR_opts$shape1 <<- input$cmp_shape1
          stoch_SAR_opts$shape2 <<- input$cmp_shape2
        }
        stoch_SAR_opts$timing <<- input$cmp_timing
        if (stoch_SAR_opts$timing == "block") {
          stoch_SAR_opts$block_years <<- seq(input$cmp_block_start, input$cmp_block_end)
        }
        if (stoch_SAR_opts$timing == "pulse") {
          stoch_SAR_opts$pulse_years <<- input$cmp_pulse_years
          stoch_SAR_opts$pulse_sd    <<- input$cmp_sd
        }
        incProgress(0.2)
        
        result_df <- purrr::map_dfr(input$alts_cmp, function(alt) {
          # set SAR mean from the chosen alt here
          stoch_SAR_opts$mean <<- base_P_list[[variant]][[alt]]$SAR_mean
          
          seed_vec <- S_seed_fore_list[[variant]]
          fc_fn <- sim_forecast_fn(
            var_nm   = variant,
            env_nm   = alt,
            flow_cfs = input$flow_cfs,
            S_seed   = seed_vec,
            spawn_dates_by_env = spawn_dates_by_env   # <-- explicit
          )
          full_out <- fc_fn()
          
          if (is.null(full_out) || !is.data.frame(full_out)) return(tibble())
          
          out <- full_out %>%
            dplyr::filter(year > max(real_years)) %>%
            dplyr::slice_head(n = yrs) %>%
            dplyr::mutate(
              spawners    = round(spawners),
              deg_day     = round(deg_day, 0),
              pre_spawn   = round(pre_spawn, 3),
              dd          = round(dd, 2),
              fry_dd      = round(fry_dd, 2),
              egg_surv    = round(egg_surv, 2),
              eff_surv    = round(eff_surv, 2),
              rear_surv   = round(rear_surv, 2),
              SAR_used    = round(SAR_used, 4),
              K_spawners  = round(K_spawners, 0),
              alternative = alt,
              variant     = variant
            )
          out
        })
        
        cmp_data(result_df)
        incProgress(0.6)
      }, tag = "compare-sim")
    })
    shiny::removeNotification("cmp_busy")  # ← was removing sim_busy before
    shiny::showNotification("Comparison complete!", type = "message", duration = 2)
  })
  
  # ---- OUTPUTS ----
  output$results_table <- renderDT({
    df <- req(sim_data())
    
    # 1) Add an iteration column, select & rename everything in the exact order
    df2 <- df %>%
      mutate(
        # round to sensible precision
        spawners   = round(spawners),
        deg_day    = round(deg_day,0),
        pre_spawn = round(pre_spawn,2),
        dd         = round(dd,2),
        fry_dd     = round(fry_dd,2),
        egg_surv   = round(egg_surv,2),
        eff_surv   = round(eff_surv,2),
        SAR_used   = round(SAR_used,4),
        K_spawners = round(K_spawners,0)
      ) %>%
      # 2) select exactly the cols you want, in the order you want:
      dplyr::select(
        year, alternative, spawners, deg_day, pre_spawn,
        dd, fry_dd, egg_surv, eff_surv, rear_surv, SAR_used, K_spawners
      )
    
    # 3) supply a vector of display names that matches that order
    display_names <- c(
      "Year",
      "Alternative",
      "Forecasted spawner abundance",
      "Degree days",
      "Pre-spawn surv",
      "Density-dep (dd)",
      "Fry × DD",
      "Etf TDM",
      "Etf survival",
      "Fry to smolt",
      "SAR",
      "K (spawners)"
    )
    
    # 4) renderDataTable
    datatable(
      df2,
      rownames  = FALSE,
      extensions = 'Buttons',
      options = list(
        dom     = 'Bfrtip',
        buttons = c('csv')
      ),
      colnames = display_names
    )
  })
  
  output$cmp_table <- renderDT({
    # 1) Pull in and massage the raw cmp_data()
    df2 <- req(cmp_data()) %>%
      mutate(
        # round to sensible precision
        spawners   = round(spawners),
        deg_day    = round(deg_day,0),
        pre_spawn = round(pre_spawn,2),
        dd         = round(dd,2),
        fry_dd     = round(fry_dd,2),
        egg_surv   = round(egg_surv,2),
        eff_surv   = round(eff_surv,2),
        SAR_used   = round(SAR_used,4),
        K_spawners = round(K_spawners,0)
      ) %>%
      # 2) select exactly the cols you want, in the order you want:
      dplyr::select(
        year, alternative, spawners, deg_day, pre_spawn,
        dd, fry_dd, egg_surv, eff_surv, rear_surv, SAR_used, K_spawners
      )
    
    # 3) supply a vector of display names that matches that order
    display_names <- c(
      "Year",
      "Alternative",
      "Forecasted spawner abundance",
      "Degree days",
      "Pre‑spawn surv",
      "Density‑dep (dd)",
      "Fry × DD",
      "Etf TDM",
      "Etf survival",
      "SAR",
      "Fry to smolt",
      "K (spawners)"
    )
    
    # 4) renderDataTable
    datatable(
      df2,
      rownames  = FALSE,
      extensions = 'Buttons',
      options = list(
        dom     = 'Bfrtip',
        buttons = c('csv')
      ),
      colnames = display_names
    )
  })
  
  # Single Alternative – Time Series
  output$ts_plot <- renderPlot({
    df <- req(sim_data())
    ggplot(df, aes(x = year, y = spawners)) +
      geom_line(color = "black", size = 1) +                 # static black
      expand_limits(y = 0) +
      labs(
        title = paste0("Spawner Estimate over Time — Alt ", input$alternative),
        x     = "Year",
        y     = "Forecasted Spawner Abundance"
      ) +
      theme_minimal(base_size = 14)
  })
  
  # Single Alternative – Distribution
  output$dist_plot <- renderPlot({
    df <- req(sim_data())
    ggplot(df, aes(x = spawners)) +
      geom_histogram(bins = 30, fill = "black", alpha = 0.7) +  # static black
      labs(
        title = paste0("Spawner Distribution — Alt ", input$alternative),
        x     = "Forecasted Spawner Abundance",
        y     = "Count"
      ) +
      theme_minimal(base_size = 14)
  })
  
  # Single Alternative – Heatmap
  output$heatmap_plot <- renderPlot({
    df <- req(sim_data()) %>%
      mutate(
        year    = factor(year, levels = as.character(seq_len(max(year)))),
        variant = factor(variant, levels = unique(variant))
      )
    ggplot(df, aes(x = year, y = variant, fill = spawners)) +
      geom_tile() +
      scale_fill_viridis_c(name = "Forecasted Spawner Abundance") +               # Viridis continuous
      labs(
        title = paste0("Spawner Heatmap — Alt ", input$alternative),
        x     = "Year",
        y     = "TDM Variant"
      ) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  })
  
  output$fry_dd_plot <- renderPlot({
    df <- req(sim_data())
    ggplot(df, aes(x=year, y=fry_dd)) +
      geom_line() +
      expand_limits(y=0) +
      labs(
        title = paste0("Fry × DD — Alt ", input$alternative),
        x     = "Year",
        y     = "Fry after TDM & Density‐dep"
      ) +
      theme_minimal(base_size=14)
  })
  
  # Compare Alternatives – Time Series
  output$cmp_ts_plot <- renderPlot({
    df <- req(cmp_data()) %>%
      mutate(
        alternative = factor(alternative, levels = env_levels)
      )
    
    # build the base plot
    p <- ggplot(df, aes(year, spawners, color=alternative, group=alternative)) +
      geom_line(size=1) +
      expand_limits(y=0) +
      scale_color_viridis_d(name="Alternative") +
      labs(
        title = "Comparison: Spawners over Time",
        x     = "Year", y = "Forecasted Spawner Abundance"
      ) +
      theme_minimal(base_size = 14)
    
    # conditionally add points only if TRUE
    p + (if (isTRUE(input$show_points)) geom_point(size=1) else NULL)
  })
  
  
  # Compare Alternatives – Distribution would use faceting or overlaid histograms,
  # but if you had it you could do:
  #    scale_fill_viridis_d() similarly.
  
  # Compare Alternatives – Boxplot (Last N Years)
  output$cmp_box_plot <- renderPlot({
    df <- req(cmp_data()) %>%
      mutate(alternative = factor(alternative, levels = env_levels))
    last_yr <- max(df$year)
    df %>%
      filter(year >= (last_yr - input$last_n + 1)) %>%
      ggplot(aes(alternative, spawners, fill = alternative)) +
      geom_boxplot() +
      scale_fill_viridis_d(name = "Alternative") +          # Viridis discrete
      expand_limits(y = 0) +
      labs(
        title = paste0("Spawner Boxplot: Last ", input$last_n, " Years"),
        x     = "Alternative",
        y     = "Forecasted Spawner Abundance"
      ) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$cmp_boxplot_stats <- renderTable({
    req(cmp_data())
    
    cmp_data() %>%
      mutate(alternative = factor(alternative, levels = env_levels)) %>%
      group_by(alternative) %>%
      summarise(
        Minimum   = min(spawners, na.rm = TRUE),
        `1st Qu.` = quantile(spawners, 0.25, na.rm = TRUE),
        Median    = median(spawners, na.rm = TRUE),
        `3rd Qu.` = quantile(spawners, 0.75, na.rm = TRUE),
        Maximum   = max(spawners, na.rm = TRUE),
        .groups   = "drop"
      ) %>%
      arrange(alternative) %>%
      mutate(alternative = as.character(alternative))
  }, rownames = FALSE, digits = 0)
  
  
  # Compare Alternatives – Heatmap
  output$cmp_heatmap <- renderPlot({
    df <- req(cmp_data()) %>%
      mutate(
        year        = factor(year, levels = as.character(seq_len(max(year)))),
        alternative = factor(alternative, levels = env_levels)
      )
    ggplot(df, aes(year, alternative, fill = spawners)) +
      geom_tile() +
      scale_fill_viridis_c(name = "Forecasted Spawner Abundance") +              # Viridis continuous
      labs(
        title = "Comparison Heatmap",
        x     = "Year",
        y     = "Alternative"
      ) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  })
  
  output$cmp_fry_dd_plot <- renderPlot({
    df <- req(cmp_data()) %>%
      mutate(alternative = factor(alternative, levels=env_levels))
    ggplot(df, aes(x=year, y=fry_dd, color=alternative)) +
      geom_line() +
      expand_limits(y=0) +
      scale_color_viridis_d(name="Alternative") +
      labs(
        title = paste0("Fry × DD Comparison"),
        x     = "Year",
        y     = "Fry after TDM & Density‐dep"
      ) +
      theme_minimal(base_size=14)
  })
}

shinyApp(ui, server)
