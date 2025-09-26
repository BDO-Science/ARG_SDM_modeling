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
# Updated to use a fixed, robust range for standardization
normalize_scores <- function(scores, objective) {
  # Define robust min-max ranges based on analysis of model outputs
  bounds <- list(
    salmon = c(15000, 40000),
    hydro = c(0, 450000),
    steelhead = c(50, 61)
  )
  
  min_s <- bounds[[objective]][1]
  max_s <- bounds[[objective]][2]
  
  if (max_s == min_s) return(rep(100, length(scores))) 
  
  # Scale scores to the 0-100 range and cap at 0 and 100
  scaled <- (scores - min_s) / (max_s - min_s) * 100
  pmax(0, pmin(100, scaled)) # Cap values between 0 and 100
}


ui <- navbarPage("Fall-run Chinook Power Bypass Simulator",
                 useShinyjs(), 
                 
                 # ---- About Tab ----
                 tabPanel("About",
                          fluidRow(
                            column(12,
                                   h3("Author & Contact Information"),
                                   tags$ul(
                                     tags$li("Author: [Your Name]"),
                                     tags$li("Contact: [Your Email]")
                                   ),
                                   h3("Application Overview"),
                                   p("This Shiny app is designed to simulate and evaluate the effects of different Folsom Dam power bypass scenarios on the fall-run Chinook salmon population in the lower American River. It integrates various sub-models to forecast long-term spawner abundance and assess performance against key objectives."),
                                   h3("Key Features"),
                                   tags$ul(
                                     tags$li(strong("Spawner Forecasting:"), " Projects future salmon populations under different assumptions."),
                                     tags$li(strong("Performance Metrics:"), " Compares alternatives across biological and economic objectives."),
                                     tags$li(strong("Decision Support:"), " Utilizes swing weighting to rank alternatives based on user preferences."),
                                     tags$li(strong("Data Exploration:"), " Allows for detailed exploration of temperature and flow data.")
                                   )
                            )
                          )
                 ),
                 
                 # ---- Main Analysis Tab ----
                 tabPanel("Spawner Forecast & Performance",
                          sidebarLayout(
                            sidebarPanel(
                              h4("Model Configuration"),
                              selectInput("tdm_variant", "Select TDM Variant:",
                                          choices = unique(df_all$variant),
                                          selected = "v1"),
                              
                              selectInput("sar_option", "Select SAR Assumption:",
                                          choices = unique(df_all$sar_name),
                                          selected = "SAR_med"),
                              
                              sliderTextInput("hydrology_year", "Select Hydrology Year (for Steelhead Metric):",
                                              choices = sort(unique(egg_summary$year)),
                                              selected = 2011,
                                              grid = TRUE),
                              
                              h4("Plot Options"),
                              sliderInput("last_n", "Focus on Last N Years of Simulation:",
                                          min = 5, max = 50, value = 20, step = 5)
                            ),
                            mainPanel(
                              tabsetPanel(
                                tabPanel("Performance Metrics", 
                                         h4("Performance Metrics Summary"),
                                         p("This table summarizes the average performance of each alternative across key metrics. Steelhead days are reactive to the selected hydrology year."),
                                         DT::dataTableOutput("perf_metrics_table")
                                ),
                                tabPanel("Spawner Forecast Plot", 
                                         plotOutput("forecast_plot", height = "600px")
                                ),
                                tabPanel("Boxplot Comparison", 
                                         h4("Spawner Distribution Comparison"),
                                         p("Comparing the distribution of forecasted spawner abundance for the last N years of the simulation period."),
                                         plotOutput("cmp_boxplot"),
                                         h4("Summary Statistics"),
                                         tableOutput("cmp_boxplot_stats")
                                ),
                                tabPanel("Heatmap", 
                                         plotOutput("cmp_heatmap", height = "800px")
                                )
                              )
                            )
                          )
                 ),
                 
                 # ---- Swing Weighting Tab ----
                 tabPanel("Swing Weighting",
                          fluidRow(
                            column(12,
                                   h3("Objective Weighting using Swing Weights"),
                                   p("Use this tool to rank and score hypothetical alternatives. Your preferences will be translated into objective weights used in the 'Decision Support' tab. This method helps clarify priorities by evaluating the 'swing' from the worst to the best outcome for each objective."),
                                   hr()
                            )
                          ),
                          fluidRow(
                            column(3, strong("Objective Description")),
                            column(2, strong("Worst Case")),
                            column(2, strong("Best Case")),
                            column(5, strong("Hypothetical Alternatives"))
                          ),
                          hr(),
                          # Salmon
                          fluidRow(
                            column(3, "Fall-run Chinook Salmon Abundance"),
                            column(2, "15,000"),
                            column(2, "40,000"),
                            column(5, 
                                   fluidRow(
                                     column(4, p(strong("Alt 1: Best in Salmon"), br(), "Gets best salmon outcome, worst on others.")),
                                     column(4, numericInput("rank_salmon", "Rank (1=best):", value = 1, min = 1, max = 3)),
                                     column(4, numericInput("score_salmon", "Score (0-100):", value = 100, min = 0, max = 100))
                                   )
                            )
                          ),
                          hr(),
                          # Hydro
                          fluidRow(
                            column(3, "Hydropower Generation Cost ($ lost)"),
                            column(2, "$450,000"),
                            column(2, "$0"),
                            column(5, 
                                   fluidRow(
                                     column(4, p(strong("Alt 2: Best in Hydro"), br(), "Gets best hydro outcome, worst on others.")),
                                     column(4, numericInput("rank_hydro", "Rank (1=best):", value = 2, min = 1, max = 3)),
                                     column(4, numericInput("score_hydro", "Score (0-100):", value = 80, min = 0, max = 100))
                                   )
                            )
                          ),
                          hr(),
                          # Steelhead
                          fluidRow(
                            column(3, "Steelhead Health (# days < 65°F)"),
                            column(2, "50"),
                            column(2, "61"),
                            column(5, 
                                   fluidRow(
                                     column(4, p(strong("Alt 3: Best in Steelhead"), br(), "Gets best steelhead outcome, worst on others.")),
                                     column(4, numericInput("rank_steelhead", "Rank (1=best):", value = 3, min = 1, max = 3)),
                                     column(4, numericInput("score_steelhead", "Score (0-100):", value = 50, min = 0, max = 100))
                                   )
                            )
                          ),
                          hr(),
                          fluidRow(
                            column(12, 
                                   h4("Calculated Objective Weights"),
                                   p("These weights are automatically applied on the 'Decision Support' page."),
                                   tableOutput("weights_table")
                            )
                          )
                 ),
                 
                 # ---- Decision Support Tab ----
                 tabPanel("Decision Support",
                          fluidRow(
                            column(12, 
                                   h3("Multi-Objective Decision Analysis"),
                                   p("This table ranks the management alternatives based on a weighted score. The scores for each objective are first normalized to a common scale (0-100) using a fixed range, and then combined using the weights from the 'Swing Weighting' tab. Higher weighted scores indicate better overall performance.")
                            )
                          ),
                          hr(),
                          DT::dataTableOutput("decision_table")
                 ),
                 
                 # ---- Temperature Data Exploration Tab ----
                 tabPanel("Temperature Explorer",
                          sidebarLayout(
                            sidebarPanel(
                              h4("Display Options"),
                              selectInput("temp_hydro_select", "Select Hydrology Year:",
                                          choices = hydro_years_placeholder,
                                          selected = 2011),
                              pickerInput("temp_alt_select", "Select Alternatives:",
                                          choices = levels(df_all$alternative),
                                          selected = levels(df_all$alternative)[1:2],
                                          options = list(`actions-box` = TRUE),
                                          multiple = TRUE)
                            ),
                            mainPanel(
                              h4("Daily Temperature Patterns"),
                              p("Explore the simulated daily water temperature for the selected hydrology and management alternatives."),
                              plotOutput("temp_plot", height = "600px")
                            )
                          ))
)


server <- function(input, output, session) {
  
  # ---- Reactive Data for Main Tab ----
  
  # Reactive expression for filtered forecast data
  cmp_data <- reactive({
    req(input$tdm_variant, input$sar_option)
    df_all %>%
      filter(variant == input$tdm_variant, sar_name == input$sar_option)
  })
  
  # Reactive expression for performance metrics
  perf_data <- reactive({
    req(input$tdm_variant, input$sar_option, input$hydrology_year)
    
    # Spawner data is averaged across all simulation years
    spawner_data <- cmp_data() %>%
      group_by(alternative) %>%
      summarise(
        avg_spawners = mean(spawners, na.rm = TRUE),
        .groups = 'drop'
      )
    
    # Steelhead data is now REACTIVE to the selected hydrology year
    steelhead_data <- egg_summary %>%
      filter(variant == input$tdm_variant, year == input$hydrology_year) %>%
      group_by(alternative) %>%
      summarise(
        days_lt_65F = `days_lt_18.3C`[1], # Just take the value for that single year
        .groups = 'drop'
      )
    
    # Join them and add placeholder for cost
    spawner_data %>%
      left_join(steelhead_data, by = "alternative") %>%
      mutate(
        cost = runif(n(), 50000, 400000) # Placeholder for hydropower cost
      )
  })
  
  # ---- Outputs for Main Tab ----
  
  output$perf_metrics_table <- DT::renderDataTable({
    df <- perf_data() %>%
      mutate(across(where(is.numeric), round, 0)) %>%
      rename(
        "Alternative" = alternative,
        "Avg. Spawners" = avg_spawners,
        "Days < 65°F" = days_lt_65F,
        "Hydropower Cost ($)" = cost
      )
    DT::datatable(df, options = list(pageLength = 20, scrollX = TRUE), rownames = FALSE)
  })
  
  output$forecast_plot <- renderPlot({
    df <- req(cmp_data())
    ggplot(df, aes(x = year, y = spawners, color = alternative)) +
      geom_smooth(se = FALSE, span = 0.1) +
      scale_y_continuous(labels = scales::comma) +
      labs(
        title = "Forecasted Spawner Abundance Over Time",
        x = "Year",
        y = "Spawner Abundance",
        color = "Alternative"
      ) +
      theme_minimal(base_size = 16) +
      theme(legend.position = "bottom")
  })
  
  output$cmp_boxplot <- renderPlot({
    req(cmp_data(), input$last_n)
    last_yr <- max(cmp_data()$year)
    df_filt <- cmp_data() %>% filter(year >= (last_yr - input$last_n + 1))
    ggplot(df_filt, aes(x = factor(alternative), y = spawners, fill = factor(alternative))) +
      geom_boxplot() +
      scale_fill_viridis_d(guide = "none") +
      scale_y_continuous(labels = scales::comma) +
      labs(
        title = paste0("Spawner Distribution: Last ", input$last_n, " Years"),
        x = "Alternative",
        y = "Forecasted Spawner Abundance"
      ) +
      theme_minimal(base_size = 16) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$cmp_boxplot_stats <- renderTable({
    req(cmp_data(), input$last_n)
    last_yr <- max(cmp_data()$year)
    cmp_data() %>%
      filter(year >= (last_yr - input$last_n + 1)) %>%
      group_by(alternative) %>%
      summarise(
        Minimum = min(spawners, na.rm = TRUE),
        `1st Qu.` = quantile(spawners, 0.25, na.rm = TRUE),
        Median = median(spawners, na.rm = TRUE),
        `3rd Qu.` = quantile(spawners, 0.75, na.rm = TRUE),
        Maximum = max(spawners, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(across(where(is.numeric), round, 0))
  }, rownames = FALSE)
  
  output$cmp_heatmap <- renderPlot({
    df <- req(cmp_data())
    ggplot(df, aes(year, alternative, fill = spawners)) +
      geom_tile() +
      scale_fill_viridis_c(name = "Spawners", labels = comma) +
      labs(title = "Comparison Heatmap", x = "Year", y = "Alternative") +
      theme_minimal(base_size = 14) +
      theme(axis.text.y = element_text(size = 8))
  })
  
  # ---- Swing Weighting Server Logic ----
  
  objective_weights <- reactive({
    scores <- c(
      salmon = input$score_salmon,
      hydro = input$score_hydro,
      steelhead = input$score_steelhead
    )
    
    # Normalize scores to get weights
    total_score <- sum(scores, na.rm = TRUE)
    if (total_score > 0) {
      weights <- scores / total_score
    } else {
      weights <- setNames(rep(1/3, 3), names(scores))
    }
    weights
  })
  
  output$weights_table <- renderTable({
    req(objective_weights())
    data.frame(
      Objective = c("Fall-run Chinook Salmon", "Hydropower Generation", "Steelhead Health"),
      Weight = scales::percent(objective_weights(), accuracy = 0.1)
    )
  })
  
  # ---- Decision Support Server Logic ----
  
  scores_data <- reactive({
    req(perf_data(), objective_weights())
    weights <- objective_weights()
    
    final_scores <- perf_data() %>%
      filter(!is.na(alternative)) %>% # Remove potentially funky rows
      mutate(
        # Normalize each metric to a 0-100 score using robust ranges
        norm_salmon = normalize_scores(avg_spawners, "salmon"),
        # Negative cost because we want to MINIMIZE cost (higher cost is worse)
        norm_hydro = normalize_scores(-cost, "hydro"), 
        norm_steelhead = normalize_scores(days_lt_65F, "steelhead"),
        
        # Calculate the final weighted score
        weighted_score = (norm_salmon * weights['salmon']) + 
          (norm_hydro * weights['hydro']) + 
          (norm_steelhead * weights['steelhead'])
      ) %>%
      arrange(desc(weighted_score))
    
    final_scores
  })
  
  output$decision_table <- DT::renderDataTable({
    df <- scores_data() %>%
      select(
        alternative,
        avg_spawners, 
        cost, 
        days_lt_65F,
        weighted_score
      ) %>%
      mutate(
        across(where(is.numeric), round, 1),
        rank = row_number()
      ) %>%
      rename(
        "Rank" = rank,
        "Alternative" = alternative,
        "Avg Spawners" = avg_spawners,
        "Hydropower Cost ($)" = cost,
        "Steelhead Days (<65°F)" = days_lt_65F,
        "Weighted Score" = weighted_score
      ) %>% 
      select(Rank, Alternative, `Weighted Score`, everything())
    
    DT::datatable(df, options = list(pageLength = 20, scrollX = TRUE), rownames = FALSE)
  })
  
  # ---- Temperature Explorer Server Logic ----
  output$temp_plot <- renderPlot({
    req(input$temp_hydro_select, input$temp_alt_select)
    
    hydro_year <- input$temp_hydro_select
    alts_to_show <- input$temp_alt_select
    
    df <- env_ext_list[[hydro_year]] %>%
      filter(alternative %in% alts_to_show)
    
    ggplot(df, aes(x = Date, y = temp_f, color = alternative)) +
      geom_line(alpha = 0.8) +
      scale_x_datetime(date_labels = "%b %d") +
      labs(
        title = paste("Simulated Daily Water Temperature in", hydro_year),
        x = "Date",
        y = "Temperature (°F)",
        color = "Alternative"
      ) +
      theme_minimal(base_size = 16) +
      theme(legend.position = "bottom")
  })
  
}

shinyApp(ui, server)
