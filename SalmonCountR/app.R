library(shiny)
library(tidyverse)
library(DT)

source("global.R")

ui <- navbarPage("American River Fall‑run Chinook Salmon Power Bypass Alternative Simulator",
                 
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
                                   h3("Simulator Overview"),
                                   p("This Shiny app models salmon spawner dynamics under power bypass water temperature alternatives. \
        It couples daily temperature–driven mortality on eggs (TDM) with an age‑structured \
        life‑cycle model, so you can explore deterministic and stochastic forecasts for \
        different management scenarios."),
                                   
                                   h4("Alternatives"),
                                   tags$ul(
                                     tags$li(strong("Baseline:"), 
                                             "Current operating rules and historical temperature/flow patterns."),
                                     tags$li(strong("Pulse Release:"), 
                                             "Periodic high‐flow pulses timed to flush redds and improve incubation conditions."),
                                     tags$li(strong("Constant Flow:"), 
                                             "A fixed year‐round flow setpoint (e.g. 2,000 cfs) to minimize thermal stress."),
                                     tags$li(strong("Adaptive Management:"), 
                                             "Flows adjusted in real time based on observed river temperature forecasts.")
                                   ),
                                   
                                   h4("Model Components"),
                                   tags$ul(
                                     tags$li(strong("Temperature‑Dependent Mortality (TDM):"),
                                             "Egg‑to‑fry survival via either exponential or linear functions of daily °C‑days  \
                (exponential variants from SALMOD 2006; Water Forum 2020; linear from Martin et al. 2017)"),
                                     tags$li(strong("Life Cycle Dynamics:"),
                                             HTML("Redds → eggs → fry (density‐dependence) → rearing → smolt → adult returns.")),
                                     tags$li(strong("Stochastic SAR:"), 
                                             "User‐specified Normal, Lognormal, Beta or Gamma noise on smolt‑to‑adult ratio,  \
                applied annually, in a contiguous block, or in select pulse years."),
                                     tags$li(strong("Calibration:"),
                                             "Rearing survival (S_rear), SAR mean, AND age‑specific lag probabilities (ages 3–5)  \
                are fit to observed escapement counts (2011–2024)."),
                                     tags$li(strong("Carrying Capacity (K_spawners):"),
                                             "Computed dynamically from the user‑selected flow (cfs) via an instream lookup directly from DSMHabitat https://github.com/CVPIA-OSC/DSMhabitat/tree/main; https://cvpia-osc.github.io/DSMhabitat/.  \
                Calculated from the instream spawning habitat area (FR_spawn_wua, m²) divided by an average redd size (9.29 m²).  \
                In the app, users select a flow (cfs), which is converted to FR_spawn_wua via lookup, \
                then K_spawners = FR_spawn_wua / 9.29."),
                                     tags$li(strong("Pre‑spawn mortality:"), 
                                             HTML("Adults incur thermal stress prior to spawning via a logistic model:  
                 <code>S<sub>pre</sub> = plogis(intercept + β × degree‑days)</code>,  
                 where “degree‑days” are the sum of daily river temperatures  
                 from adult river entry (e.g. Sept 15) through their spawn date.  
                 Users can inspect cumulative degree‑days and resulting survival in the output table.")
                                     )
                                   ),
                                   
                                   h4("Key Equations"),
                                   tags$ul(
                                     tags$li(HTML("Egg development days = <code>958 / T</code>; Fry emergence days = <code>417 / T</code>")),
                                     tags$li(HTML("<strong>Exponential TDM:</strong>  S<sub>egg</sub> = exp(-Σ α exp(β T<sub>i</sub>))")),
                                     tags$li(HTML("<strong>Linear TDM:</strong>  S<sub>egg</sub> = exp(-α Σ max(T<sub>i</sub> − β, 0))")),
                                     tags$li(HTML("<strong>Density dependence:</strong>  dd = S₀ / (1 + redds / K<sub>spawners</sub>)")),
                                     tags$li(HTML("<strong>Rearing survival:</strong>  fry → smolt = eggs × S<sub>egg</sub> × dd × S<sub>rear</sub>")),
                                     tags$li(HTML("<strong>Age‑structured returns:</strong>  S[t+age] += reared[t] × SAR[t] × lag<sub>age</sub>  (ages 3–5)"))
                                   ),
                                   
                                   h4("Data & Calibration"),
                                   tags$ul(
                                     tags$li("Daily temperature and flow: pre‑computed in env_ext_list."),
                                     tags$li("Observed redd dates: carcass surveys → sim_redds (2011–2024)."),
                                     tags$li("Egg‑to‑fry survival (egg_summary) from TDM simulations."),
                                     tags$li("Age‑structured calibration using optim() against escapement counts."),
                                     tags$li("Dynamic carrying capacity (K_spawners) derived from flow via instream habitat lookup.")
                                   ),
                                   
                                   h4("How to Use This App"),
                                   tags$ol(
                                     tags$li(strong("About:"), "Review model structure, assumptions, and references."),
                                     tags$li(strong("Single Alternative:"), "Pick one scenario, TDM variant, stochastic options,  \
                     and flow value; then click Run Simulation to see tables and plots (incl. K_spawners)."),
                                     tags$li(strong("Compare Alternatives:"), "Select multiple scenarios & flow to compare via time series,  \
                     boxplots, heatmaps, and tables."),
                                     tags$li(strong("Settings:"), "Toggle gridlines or point overlays on plots.")
                                   ),
                                   
                                   h4("References"),
                                   tags$ul(
                                     tags$li("Bartholow, J.M. & Heasley, J. 2006. Evaluation of Shasta Dam Scenarios Using a Salmon Production Model. USGS Report."),
                                     tags$li("Bratovich, P., Neal, M., Ransom, A., et al. 2020. Chinook Salmon Early Lifestage Survival and Folsom Dam Power Bypass Considerations. Water Forum Technical Memo, Aug 2020."),
                                     tags$li("CFS 2010. A Revised Sacramento River Winter Chinook Salmon Juvenile Production Model. Cramer Fish Sciences."),
                                     tags$li("HCI 1996. Chinook Salmon Mortality Model: Development, Evaluation, and Application. Hydrologic Consultants Inc."),
                                     tags$li("\"Heat‐stress & pre‐spawn survival,\" Wiley 2021. https://onlinelibrary.wiley.com/doi/full/10.1111/rec.13244"),
                                     tags$li("Martin, B.T., et al. 2017. Phenomenological vs. biophysical models of thermal stress in aquatic eggs. Ecol Lett 20:50–59. https://doi.org/10.1111/ele.12705"),
                                     tags$li("USBR 2008. Biological Assessment on Continued Operations of CVP & SWP. US Bureau of Reclamation."),
                                     tags$li("Zeug, S., Bergman, P., Cavallo, B., & Jones, K. 2012. Application of a Life Cycle Simulation Model… Environ Model Assess. https://doi.org/10.1007/s10666-012-9306-6")
                                   )
                            )
                          )
                 )
                 ,
                 
                 # ---- Single-alt tab ------------------------------------------------------
                 tabPanel("Single Alternative",
                          sidebarLayout(
                            sidebarPanel(
                              h4("Simulation Settings"),
                              sliderInput("sim_years", "Simulation Length (years)",
                                          min = 10, max = 100, value = n_sim, step = 1),
                              selectInput("mode", "Mode:", choices = c("Deterministic", "Stochastic")),
                              conditionalPanel(
                                condition = "input.mode == 'Stochastic'",
                                selectInput("stoch_type", "Distribution Type:",
                                            choices = c("Normal", "Lognormal", "Beta", "Gamma")),
                                sliderInput("stoch_sd", "SAR SD:", min = 0, max = 0.01,
                                            value = stoch_SAR_opts$sd, step = 0.0001),
                                conditionalPanel(
                                  condition = "input.stoch_type == 'Beta' || input.stoch_type == 'Gamma'",
                                  numericInput("shape1", "Shape1:", value = stoch_SAR_opts$shape1, min = 0.1, step = 0.1),
                                  numericInput("shape2", "Shape2:", value = stoch_SAR_opts$shape2, min = 0.1, step = 0.1)
                                ),
                                selectInput("stoch_timing", "Timing:", choices = c("all", "block", "pulse"),
                                            selected = stoch_SAR_opts$timing),
                                conditionalPanel(
                                  condition = "input.stoch_timing == 'block'",
                                  sliderInput("block_start", "Block Start Year:", min = 1, max = 100, value = min(stoch_SAR_opts$block_years), step = 1),
                                  sliderInput("block_end",   "Block End Year:",   min = 1, max = 100, value = max(stoch_SAR_opts$block_years),   step = 1)
                                ),
                                conditionalPanel(
                                  condition = "input.stoch_timing == 'pulse'",
                                  selectizeInput("pulse_years", "Pulse Years:",
                                                 choices = 1:100, multiple = TRUE,
                                                 selected = stoch_SAR_opts$pulse_years)
                                )
                              ),
                              selectInput("tdm_variant", "TDM Variant:",
                                          choices = names(base_P_list)),
                              selectInput("alternative", "Alternative:",
                                          choices = unique(egg_summary$env)),
                              numericInput("flow_cfs", "Flow (cfs):",
                                           value = 2000,
                                           min   = min(instream$flow_cfs),
                                           max   = max(instream$flow_cfs),
                                           step  = 50),
                              actionButton("run_sim", "Run Simulation"),
                              width = 3
                            ),
                            mainPanel(
                              tabsetPanel(
                                tabPanel("Table", DTOutput("results_table")),
                                tabPanel("Time Series", plotOutput("ts_plot")),
                                tabPanel("Distribution", plotOutput("dist_plot")),
                                tabPanel("Heatmap", plotOutput("heatmap_plot")),
                                tabPanel("SAR Series", plotOutput("sar_plot"))
                              ),
                              width = 9
                            )
                          )
                 ),
                 
                 # ---- Compare-alts tab ----------------------------------------------------
                 tabPanel("Compare Alternatives",
                          sidebarLayout(
                            sidebarPanel(
                              h4("Comparison Settings"),
                              sliderInput("cmp_years", "Simulation Length (years)",
                                          min = 10, max = 100, value = n_sim, step = 1),
                              selectInput("cmp_mode", "Mode:", choices = c("Deterministic", "Stochastic")),
                              conditionalPanel(
                                condition = "input.cmp_mode == 'Stochastic'",
                                selectInput("cmp_stoch_type", "Distribution Type:",
                                            choices = c("Normal", "Lognormal", "Beta", "Gamma")),
                                sliderInput("cmp_sd", "SAR SD:", min = 0, max = 0.01,
                                            value = stoch_SAR_opts$sd, step = 0.0001),
                                conditionalPanel(
                                  condition = "input.cmp_stoch_type == 'Beta' || input.cmp_stoch_type == 'Gamma'",
                                  numericInput("cmp_shape1", "Shape1:", value = stoch_SAR_opts$shape1, min = 0.1, step = 0.1),
                                  numericInput("cmp_shape2", "Shape2:", value = stoch_SAR_opts$shape2, min = 0.1, step = 0.1)
                                ),
                                selectInput("cmp_timing", "Timing:", choices = c("all", "block", "pulse"),
                                            selected = stoch_SAR_opts$timing),
                                conditionalPanel(
                                  condition = "input.cmp_timing == 'block'",
                                  sliderInput("cmp_block_start", "Block Start Year:", min = 1, max = 100, value = min(stoch_SAR_opts$block_years), step = 1),
                                  sliderInput("cmp_block_end",   "Block End Year:",   min = 1, max = 100, value = max(stoch_SAR_opts$block_years),   step = 1)
                                ),
                                conditionalPanel(
                                  condition = "input.cmp_timing == 'pulse'",
                                  selectizeInput("cmp_pulse_years", "Pulse Years:",
                                                 choices = 1:100, multiple = TRUE,
                                                 selected = stoch_SAR_opts$pulse_years)
                                )
                              ),
                              selectInput("cmp_variant", "TDM Variant:",
                                          choices = names(base_P_list)),
                              selectInput("alts_cmp", "Alternatives to compare:",
                                          choices = unique(egg_summary$env), multiple = TRUE,
                                          selected = unique(egg_summary$env)),
                              numericInput("flow_cfs", "Flow (cfs):",
                                           value = 2000,
                                           min   = min(instream$flow_cfs),
                                           max   = max(instream$flow_cfs),
                                           step  = 50),
                              actionButton("run_cmp", "Run Comparison"),
                              width = 3
                            ),
                            mainPanel(
                              tabsetPanel(
                                tabPanel("Table", DTOutput("cmp_table")),
                                tabPanel("Time Series", plotOutput("cmp_ts_plot")),
                                tabPanel("Boxplot (Last N Years)",
                                         sliderInput("last_n", "Last N Years:",
                                                     min = 5, max = 50, value = 10, step = 1),
                                         plotOutput("cmp_box_plot")
                                ),
                                tabPanel("Heatmap", plotOutput("cmp_heatmap")),
                                
                              ),
                              width = 9
                            )
                          )
                 ),
                 
                 # ---- Optional settings ---------------------------------------------------
                 # ---- Enhanced Settings Tab ---------------------------------------------------
                 tabPanel("Settings",
                          sidebarPanel(
                            h4("Visualization Options"),
                            checkboxInput("show_grid", "Show grid lines on plots", TRUE),
                            checkboxInput("show_points", "Overlay points on time series", FALSE),
                            selectInput("color_scheme", "Color Scheme:",
                                        choices = c("Viridis", "Plasma", "Cividis", "Inferno"),
                                        selected = "Viridis"),
                            selectInput("plot_theme", "Plot Theme:",
                                        choices = c("Minimal", "Classic", "Black & White", "Gray"),
                                        selected = "Minimal"),
                            hr(),
                            
                            h4("Model Options"),
                            checkboxInput("verbose_output", "Enable verbose model output", FALSE),
                            sliderInput("default_sim_years", "Default simulation length:",
                                        min = 10, max = 100, value = n_sim, step = 1),
                            hr(),
                            
                            h4("Export Options"),
                            numericInput("csv_precision", "CSV export precision (digits):",
                                         value = 2, min = 0, max = 8, step = 1),
                            textInput("export_folder", "Export folder (optional):",
                                      value = ""),
                            checkboxInput("auto_export", "Automatically export results on simulation", FALSE)
                          ),
                          mainPanel(
                            h5("Customize visualization, model behavior, and export settings here."),
                            p("Use these settings to tailor the simulator to your preferences. Changes here affect all tabs and outputs."),
                            width = 9
                          )
                 )
)

server <- function(input, output, session) {
  env_levels <- as.character(sort(as.numeric(unique(egg_summary$env))))
  real_years <- 2011:2024
  calib_end  <- length(real_years)        # = 14
  shade      <- data.frame(xmin=0.5, xmax=calib_end + 0.5)
  
  # storage for the latest sim & cmp results
  sim_data <- reactiveVal(NULL)
  cmp_data <- reactiveVal(NULL)
  
  # ---- SINGLE‑ALT SIMULATION ----
  observeEvent(input$run_sim, {
    showNotification("Running simulation…", id = "sim_busy", duration = NULL)
    withProgress(message = "Simulating spawners", value = 0, {
      incProgress(0.1)
      
      # 1) grab inputs
      yrs     <- input$sim_years
      alt     <- input$alternative
      variant <- input$tdm_variant
      P_out   <- base_P_list[[variant]]
      
      # 2) build SAR options
      opts <- stoch_SAR_opts
      opts$model  <- tolower(input$stoch_type)
      opts$mean   <- P_out$SAR_mean
      opts$sd     <- input$stoch_sd
      if (opts$model %in% c("beta","gamma")) {
        opts$shape1 <- input$shape1
        opts$shape2 <- input$shape2
      }
      opts$timing <- input$stoch_timing
      if (opts$timing == "block") {
        opts$block_years <- seq(input$block_start, input$block_end)
      }
      if (opts$timing == "pulse") {
        req(length(input$pulse_years) > 0)
        opts$pulse_years <- input$pulse_years
        opts$pulse_sd    <- input$stoch_sd
      }
      
      incProgress(0.1)
      # 3) generate SAR
      SAR <- if (input$mode=="Stochastic") {
        generate_SAR_vec(yrs, opts)
      } else {
        rep(P_out$SAR_mean, yrs)
      }
      
      message("SAR[1:20] = ", paste0(round(SAR[1:20],4), collapse = ", "))
      incProgress(0.3)
      
      # 4) dynamic K_spawners from user flow
      K_spaw_vec <- rep(get_K_spawners(input$flow_cfs), yrs)
      
      years_vec<- sim_years_full[1:yrs]
      
      deg_days <- compute_deg_day_adult(
        env_nm       = alt,
        sim_years    = years_vec,
        spawn_dates  = spawn_dates_vec,
        env_ext_list = env_ext_list
      )
      
      # 7) run the life‑cycle + pre‑spawn mortality
      out <- simulate_variant(
        surv_vec       = surv_lookup_full[[paste(alt,variant,sep="_")]],
        P              = P_out,
        years          = yrs,
        S_init         = S_seed,
        SAR_vec        = SAR,
        K_spawners_vec = K_spaw_vec,
        deg_day_adult  = deg_days
      ) %>% mutate(
        year        = seq_len(yrs),
        deg_day     = round(deg_days,    0),
        pre_spawn   = round(pre_spawn,   3),
        calibrated  = year <= calib_end,     # so you can grey‑out those rows
        alternative = alt,
        variant     = variant,
        S_rear      = P_out$S_rear
      )
      
      incProgress(1)
      sim_data(out)
    })
    removeNotification("sim_busy")
    showNotification("Simulation complete!", type = "message", duration = 2)
  })
  
  
  # ---- COMPARE-ALTS SIMULATION ----
  observeEvent(input$run_cmp, {
    showNotification("Running comparison…", id = "cmp_busy", duration = NULL)
    withProgress(message = "Comparing alternatives", value = 0, {
      incProgress(0.1)
      
      yrs   <- input$cmp_years
      P_cmp <- base_P_list[[input$cmp_variant]]
      
      # build SAR vector
      if (input$cmp_mode == "Stochastic") {
        # 1) start from your defaults
        opts <- stoch_SAR_opts
        
        # 2) overwrite with cmp inputs
        opts$model  <- tolower(input$cmp_stoch_type)
        opts$mean   <- P_cmp$SAR_mean
        opts$sd     <- input$cmp_sd
        # if beta/gamma, overwrite shapes too
        opts$shape1 <- input$cmp_shape1
        opts$shape2 <- input$cmp_shape2
        
        # timing
        opts$timing <- input$cmp_timing
        if (opts$timing == "block") {
          opts$block_years <- seq(input$cmp_block_start, input$cmp_block_end)
        }
        if (opts$timing == "pulse") {
          opts$pulse_years <- as.integer(input$pulse_years)
          opts$pulse_sd    <- input$cmp_sd
        }
        
        SAR <- generate_SAR_vec(yrs, opts)
      } else {
        # deterministic fallback
        SAR <- rep(P_cmp$SAR_mean, yrs)
      }
      
      incProgress(0.3)
      
      # single flow → capacity
      K_spaw_vec <- rep(get_K_spawners(input$flow_cfs), yrs)
      
      result_df <- map_df(input$alts_cmp, function(alt) {
        key <- paste(alt, input$cmp_variant, sep="_")
        
        years_vec<- sim_years_full[1:yrs]
        
        deg_days <- compute_deg_day_adult(
          env_nm       = alt,
          sim_years    = years_vec,
          spawn_dates  = spawn_dates_vec,
          env_ext_list = env_ext_list
        )
        
        # --- now call the model, passing in deg_days ---
        simulate_variant(
          surv_vec         = surv_lookup_full[[key]],
          P                = P_cmp,
          years            = yrs,
          S_init           = S_seed,
          SAR_vec          = SAR,
          K_spawners_vec   = K_spaw_vec,
          deg_day_adult    = deg_days
        ) %>%
          mutate(
            year        = seq_len(yrs),
            deg_day     = deg_days,    # cumulative °C·days
            pre_spawn   = pre_spawn,   # logistic survival
            calibrated  = year <= calib_end,
            alternative = factor(alt, levels = env_levels),
            S_rear      = P_cmp$S_rear
          )
      }) 
      
      incProgress(1)
      cmp_data(result_df)
    })
    removeNotification("cmp_busy")
    showNotification("Comparison complete!", type = "message", duration = 2)
  })
  
  
  # ---- OUTPUTS ----
  output$results_table <- renderDT({
    df <- req(sim_data())
    
    # 1) Add an iteration column, select & rename everything in the exact order
    df2 <- df %>%
      mutate(
        # round to sensible precision
        calibrated = ifelse(calibrated, "✓", ""),
        spawners   = round(spawners),
        deg_day    = round(deg_day,0),
        pre_spawn = round(pre_spawn,4),
        dd         = round(dd,2),
        fry_dd     = round(fry_dd,2),
        eff_surv   = round(eff_surv,2),
        SAR_used   = round(SAR_used,4),
        S_rear     = round(S_rear,2),
        K_spawners = round(K_spawners,0)
      ) %>%
      # 2) select exactly the cols you want, in the order you want:
      select(
        year,
        alternative,
        spawners,
        deg_day,
        pre_spawn,
        dd,
        fry_dd,
        eff_surv,
        SAR_used,
        S_rear,
        K_spawners,
        calibrated
      )
    
    # 3) supply a vector of display names that matches that order
    display_names <- c(
      "Year",
      "Alternative",
      "Spawners",
      "Degree days",
      "Pre‑spawn surv",
      "Density‑dep (dd)",
      "Fry × DD",
      "Etf survival",
      "SAR",
      "Rearing surv",
      "K (spawners)",
      "Calibrated?"
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
    ) %>%
      # 5) grey out the calibrated rows
      formatStyle(
        'calibrated',
        target = 'row',
        backgroundColor = styleEqual("✓", "lightgray")
      )
  })
  
  output$cmp_table <- renderDT({
    # 1) Pull in and massage the raw cmp_data()
    df2 <- req(cmp_data()) %>%
      mutate(
        # round to sensible precision
        calibrated = ifelse(calibrated, "✓", ""),
        spawners   = round(spawners),
        deg_day    = round(deg_day,0),
        pre_spawn = round(pre_spawn,4),
        dd         = round(dd,2),
        fry_dd     = round(fry_dd,2),
        eff_surv   = round(eff_surv,2),
        SAR_used   = round(SAR_used,4),
        S_rear     = round(S_rear,2),
        K_spawners = round(K_spawners,0)
      ) %>%
      # 2) select exactly the cols you want, in the order you want:
      select(
        year,
        alternative,
        spawners,
        deg_day,
        pre_spawn,
        dd,
        fry_dd,
        eff_surv,
        SAR_used,
        S_rear,
        K_spawners,
        calibrated
      )
    
    # 3) supply a vector of display names that matches that order
    display_names <- c(
      "Year",
      "Alternative",
      "Spawners",
      "Degree days",
      "Pre‑spawn surv",
      "Density‑dep (dd)",
      "Fry × DD",
      "Etf survival",
      "SAR",
      "Rearing surv",
      "K (spawners)",
      "Calibrated?"
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
    ) %>%
      # 5) grey out the calibrated rows
      formatStyle(
        'calibrated',
        target = 'row',
        backgroundColor = styleEqual("✓", "lightgray")
      )
  })
  
  # Single Alternative – Time Series
  output$ts_plot <- renderPlot({
    df <- req(sim_data())
    ggplot(df, aes(x = year, y = spawners)) +
      annotate("rect",
               xmin = shade$xmin, xmax = shade$xmax,
               ymin = -Inf,      ymax = +Inf,
               fill = "grey80",  alpha = 0.3) +
      geom_line(color = "black", size = 1) +                 # static black
      { if (input$show_points) geom_point(color = "black") }+
      expand_limits(y = 0) +
      labs(
        title = paste0("Spawners over Time — Alt ", input$alternative),
        x     = "Year",
        y     = "Spawners"
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
        x     = "Spawners",
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
      scale_fill_viridis_c(name = "Spawners") +               # Viridis continuous
      labs(
        title = paste0("Spawner Heatmap — Alt ", input$alternative),
        x     = "Year",
        y     = "TDM Variant"
      ) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  })
  
  # Single Alternative – SAR Series
  output$sar_plot <- renderPlot({
    df <- req(sim_data())
    ggplot(df, aes(x = year, y = SAR_used)) +
      annotate("rect",
               xmin = shade$xmin, xmax = shade$xmax,
               ymin = -Inf,      ymax = +Inf,
               fill = "grey80",  alpha = 0.3) +
      geom_line(color = "black", size = 1) +                 # static black
      { if (input$show_points) geom_point(color = "black") }+
      expand_limits(y = 0) +
      labs(
        title = paste0("SAR over Time — Alt ", input$alternative),
        x     = "Year",
        y     = "Smolt‑to‑Adult Ratio"
      ) +
      theme_minimal(base_size = 14)
  })
  
  # Compare Alternatives – Time Series
  output$cmp_ts_plot <- renderPlot({
    df <- req(cmp_data()) %>%
      mutate(
        year        = factor(year, levels = as.character(seq_len(max(year)))),
        alternative = factor(alternative, levels = env_levels)
      )
    ggplot(df, aes(year, spawners, color = alternative, group = alternative)) +
      annotate("rect",
               xmin = shade$xmin, xmax = shade$xmax,
               ymin = -Inf,      ymax = +Inf,
               fill = "grey80",  alpha = 0.3) +
      geom_line(size = 1) +
      { if (input$show_points) geom_point(size = 1) } +
      expand_limits(y = 0) +
      scale_color_viridis_d(name = "Alternative") +          # Viridis discrete
      labs(
        title = "Comparison: Spawners over Time",
        x     = "Year",
        y     = "Spawners"
      ) +
      theme_minimal(base_size = 14)
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
        y     = "Spawners"
      ) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # Compare Alternatives – Heatmap
  output$cmp_heatmap <- renderPlot({
    df <- req(cmp_data()) %>%
      mutate(
        year        = factor(year, levels = as.character(seq_len(max(year)))),
        alternative = factor(alternative, levels = env_levels)
      )
    ggplot(df, aes(year, alternative, fill = spawners)) +
      geom_tile() +
      scale_fill_viridis_c(name = "Spawners") +              # Viridis continuous
      labs(
        title = "Comparison Heatmap",
        x     = "Year",
        y     = "Alternative"
      ) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  })
}

shinyApp(ui, server)
