library(shiny)
library(tidyverse)
library(DT)

source("global.R")

# helper for plotting order of alternatives
env_levels <- as.character(sort(as.numeric(unique(egg_summary$env))))

ui <- navbarPage("Salmon Life Cycle Simulator",
                 
                 # ---- Info tab -----------------------------------------------------------
                 tabPanel("About",
                          fluidRow(
                            column(12,
                                   h3("Life Cycle Simulator Overview"),
                                   p("This Shiny app models salmon spawner dynamics under alternative thermal regimes.  
        It couples daily temperature–driven mortality on eggs (TDM) with an age‑structured  
        life‑cycle model, so you can explore deterministic and stochastic forecasts for  
        different management scenarios."),
                                   
                                   h4("Model Components"),
                                   tags$ul(
                                     tags$li(strong("Temperature‑Dependent Mortality (TDM):"),
                                             "Egg‑to‑fry survival via either exponential or linear functions of daily °C‑days  
                (exponential variants from SALMOD 2006; Water Forum 2020; linear from Martin et al. 2017" ,  
                                             tags$em("(Martin et al. 2017)")),  
                                     tags$li(strong("Life Cycle Dynamics:"),
                                             HTML("Redds → eggs → fry (density‐dependence) → rearing → smolt → adult returns.")),
                                     tags$li(strong("Stochastic SAR:"), 
                                             "User‐specified Normal, Lognormal, Beta or Gamma noise on smolt‑to‑adult ratio,  
                applied annually, in a contiguous block, or in select pulse years."),
                                     tags$li(strong("Calibration:"),
                                             "Rearing survival (S_rear), SAR mean, AND age‑specific lag probabilities (ages 3–5)  
                are fit to observed escapement counts (2011–2024)."),
                                     tags$li(strong("Pre‑spawn mortality:"), 
                                             HTML("Not yet implemented; future versions will add adult thermal stress prior to spawning."))
                                   ),
                                   
                                   h4("Key Equations"),
                                   tags$ul(
                                     tags$li(HTML("Egg development days = <code>958 / T</code>; Fry emergence days = <code>417 / T</code>")),
                                     tags$li(HTML("<strong>Exponential TDM:</strong>  
                      S<sub>egg</sub> = exp(-Σ α exp(β T<sub>i</sub>))")),
                                     tags$li(HTML("<strong>Linear TDM:</strong>  
                      S<sub>egg</sub> = exp(-α Σ max(T<sub>i</sub> − β, 0))")),
                                     tags$li(HTML("<strong>Density dependence:</strong>  
                      dd = S₀ / (1 + redds / K<sub>spawners</sub>)")),
                                     tags$li(HTML("<strong>Rearing survival:</strong>  
                      fry → smolt = eggs × S<sub>egg</sub> × dd × S<sub>rear</sub>")),
                                     tags$li(HTML("<strong>Age‑structured returns:</strong>  
                      S[t+age] += reared[t] × SAR[t] × lag<sub>age</sub>  (ages 3–5)"))
                                   ),
                                   
                                   h4("Data & Calibration"),
                                   tags$ul(
                                     tags$li("Daily temperature and flow: pre‑computed in env_ext_list."),
                                     tags$li("Observed redd dates: carcass surveys → sim_redds (2011–2024)."),
                                     tags$li("Egg‑to‑fry survival (egg_summary) from TDM simulations."),
                                     tags$li("Age‑structured calibration using optim() against escapement counts.")
                                   ),
                                   
                                   h4("How to Use This App"),
                                   tags$ol(
                                     tags$li(strong("About:"), "Review model structure, assumptions, and references."),
                                     tags$li(strong("Single Alternative:"), "Pick one scenario, TDM variant, stochastic options,  
                     then click Run Simulation to see tables and plots."),
                                     tags$li(strong("Compare Alternatives:"), "Select multiple scenarios to compare via time series,  
                     boxplots, heatmaps, and tables."),
                                     tags$li(strong("Settings:"), "Toggle gridlines or point overlays on plots.")
                                   ),
                                   
                                   h4("References"),
                                   tags$ul(
                                     tags$li("Martin, B.T., et al. 2017. Phenomenological vs. biophysical models of thermal stress in aquatic eggs. Ecol Lett 20:50–59. https://doi.org/10.1111/ele.12705"),
                                     tags$li("Bratovich, P., Neal, M., Ransom, A., et al. 2020. Chinook Salmon Early Lifestage Survival and Folsom Dam Power Bypass Considerations. Water Forum Technical Memo, Aug 2020."),
                                     tags$li("HCI 1996. Chinook Salmon Mortality Model: Development, Evaluation, and Application. Hydrologic Consultants Inc."),
                                     tags$li("Bartholow, J.M. & Heasley, J. 2006. Evaluation of Shasta Dam Scenarios Using a Salmon Production Model. USGS Report."),
                                     tags$li("USBR 2008. Biological Assessment on Continued Operations of CVP & SWP. US Bureau of Reclamation."),
                                     tags$li("Zeug, S., Bergman, P., Cavallo, B., & Jones, K. 2012. Application of a Life Cycle Simulation Model… Environ Model Assess. https://doi.org/10.1007/s10666-012-9306-6"),
                                     tags$li("CFS 2010. A Revised Sacramento River Winter Chinook Salmon Juvenile Production Model. Cramer Fish Sciences."),
                                     tags$li("“Heat‐stress & pre‐spawn survival,” Wiley 2021. https://onlinelibrary.wiley.com/doi/full/10.1111/rec.13244")
                                   )
                            )
                          )
                 ),
                 
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
                              selectInput("alternative", "Alternative (env):",
                                          choices = unique(egg_summary$env)),
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
                              actionButton("run_cmp", "Run Comparison"),
                              width = 3
                            ),
                            mainPanel(
                              tabsetPanel(
                                tabPanel("Time Series", plotOutput("cmp_ts_plot")),
                                tabPanel("Boxplot (Last N Years)",
                                         sliderInput("last_n", "Last N Years:",
                                                     min = 5, max = 50, value = 10, step = 1),
                                         plotOutput("cmp_box_plot")
                                ),
                                tabPanel("Heatmap", plotOutput("cmp_heatmap")),
                                tabPanel("Data Table", DTOutput("cmp_table"))
                              ),
                              width = 9
                            )
                          )
                 ),
                 
                 # ---- Optional settings ---------------------------------------------------
                 tabPanel("Settings",
                          sidebarPanel(
                            checkboxInput("show_grid", "Show grid lines on plots", TRUE),
                            checkboxInput("show_points", "Overlay points on time series", FALSE)
                          ),
                          mainPanel(
                            h5("Toggle grid and points for enhanced visualization."),
                            width = 9
                          )
                 )
)

server <- function(input, output, session) {
  env_levels <- as.character(sort(as.numeric(unique(egg_summary$env))))
  
  # storage for the latest sim & cmp results
  sim_data <- reactiveVal(NULL)
  cmp_data <- reactiveVal(NULL)
  
  # ---- SINGLE-ALT SIMULATION ----
  observeEvent(input$run_sim, {
    showNotification("Running simulation…", id = "sim_busy", duration = NULL)
    withProgress(message = "Simulating spawners", value = 0, {
      incProgress(0.1)
      
      ## 1) grab inputs
      yrs   <- input$sim_years
      P_out <- base_P_list[[input$tdm_variant]]
      
      ## 2) start from your defaults and overwrite
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
        opts$pulse_years <- input$pulse_years      # still character, but generate_SAR_vec() will coerce
        opts$pulse_sd    <- input$stoch_sd
      }
      
      incProgress(0.1)
      ## 3) generate SAR entirely inside generate_SAR_vec()
      SAR <- if (input$mode == "Stochastic") {
        generate_SAR_vec(yrs, opts)
      } else {
        rep(P_out$SAR_mean, yrs)
      }
      
      message("SAR[1:20] = ", paste0(round(SAR[1:20],4), collapse = ", "))
      incProgress(0.3)
      
      ## 4) run the model
      out <- simulate_variant(
        surv_lookup_full[[paste(input$alternative, input$tdm_variant, sep = "_")]],
        P_out, yrs, S_seed, SAR
      ) %>%
        mutate(
          env     = factor(input$alternative, levels = env_levels),
          variant = input$tdm_variant,
          S_rear  = P_out$S_rear
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
      result_df <- map_df(input$alts_cmp, function(alt) {
        key <- paste(alt, input$cmp_variant, sep = "_")
        simulate_variant(
          surv_lookup_full[[key]], P_cmp, yrs, S_seed, SAR
        ) %>%
          mutate(
            env    = factor(alt, levels = env_levels),
            S_rear = P_cmp$S_rear
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
    df %>%
      mutate(
        spawners = round(spawners),
        dd       = round(dd,2),
        fry_dd   = round(fry_dd,2),
        eff_surv = round(eff_surv,2),
        SAR_used = round(SAR_used,4),
        S_rear   = round(S_rear,2)
      ) %>%
      datatable(options = list(pageLength = 10))
  })
  
  output$cmp_table <- renderDT({
    df <- req(cmp_data())
    df %>%
      mutate(
        spawners = round(spawners),
        dd       = round(dd,2),
        fry_dd   = round(fry_dd,2),
        eff_surv = round(eff_surv,2),
        SAR_used = round(SAR_used,4),
        S_rear   = round(S_rear,2)
      ) %>%
      datatable(options = list(pageLength = 10))
  })
  
  output$ts_plot <- renderPlot({
    df <- req(sim_data())   # wait until we have simulation data
    ggplot(df, aes(year, spawners)) +
      geom_line(color = "black") +
      expand_limits(y = 0) +
      { if (input$show_points) geom_point() } +
      { if (input$show_grid)  theme(panel.grid.minor = element_line()) } +
      labs(title = "Spawners over Time", y = "Spawners") +
      theme_minimal()
  })
  
  output$dist_plot <- renderPlot({
    df <- req(sim_data())
    ggplot(df, aes(spawners)) +
      geom_histogram(bins = 30, fill = "black", alpha = 0.7) +
      labs(title = "Spawner Distribution", x = "Spawners") +
      theme_minimal()
  })
  
  output$heatmap_plot <- renderPlot({
    df <- req(sim_data())
    df %>%
      mutate(year = factor(year)) %>%
      ggplot(aes(year, variant, fill = spawners)) +
      geom_tile() +
      scale_fill_viridis_c() +
      expand_limits(y = 0) +
      labs(title = "Spawner Heatmap") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  })
  
  output$cmp_ts_plot <- renderPlot({
    df <- req(cmp_data())   # wait until we have comparison data
    ggplot(df, aes(year, spawners, color = env)) +
      geom_line() +
      expand_limits(y = 0) +
      { if (input$show_points) geom_point(size = 1) } +
      { if (input$show_grid)  theme(panel.grid.minor = element_line()) } +
      labs(title = "Comparison: Spawners over Time", color = "Env") +
      theme_minimal()
  })
  
  output$cmp_box_plot <- renderPlot({
    df <- req(cmp_data())
    last_yr <- max(df$year)
    df %>%
      filter(year >= (last_yr - input$last_n + 1)) %>%
      ggplot(aes(env, spawners, fill = env)) +
      geom_boxplot() +
      expand_limits(y = 0) +
      labs(title = paste0("Spawner Boxplot: Last ", input$last_n, " Years")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$cmp_heatmap <- renderPlot({
    df <- req(cmp_data())
    df %>%
      mutate(year = factor(year)) %>%
      ggplot(aes(year, env, fill = spawners)) +
      geom_tile() +
      scale_fill_viridis_c() +
      expand_limits(y = 0) +
      labs(title = "Comparison Heatmap") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  })
  
  output$sar_plot <- renderPlot({
    df <- req(sim_data())
    tibble(year = seq_len(nrow(df)), SAR = df$SAR_used) %>%
      ggplot(aes(year, SAR)) +
      geom_point() +
      geom_line() +
      labs(title = "SAR over Time", y = "SAR", x = "Year") +
      theme_minimal()
  })
}

shinyApp(ui, server)
