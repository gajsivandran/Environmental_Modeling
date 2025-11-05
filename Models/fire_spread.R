# app.R — Forest Disturbance Model with Age-Dependent Ignition and Simple Fire Spread
# Adds: (1) ignition prob increases with age, (2) 4-neighbor spread per burn
# Simple, interpretable, and visually engaging

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(shinythemes)

ui <- fluidPage(
  theme = shinytheme("flatly"),
  titlePanel("Forest Disturbance — Age-Dependent Ignition & Simple Fire Spread"),

  sidebarLayout(
    sidebarPanel(
      sliderInput("grid_size", "Grid Size (NxN):",
                  min = 5, max = 50, value = 20),
      sliderInput("r", "Intrinsic Growth Rate (r):",
                  min = 0.01, max = 0.6, value = 0.2, step = 0.01),
      sliderInput("K", "Carrying Capacity (K):",
                  min = 10, max = 500, value = 100, step = 10),

      tags$hr(),
      h4("Disturbance Settings"),
      sliderInput("p_base", "Baseline ignition probability:",
                  min = 0, max = 0.1, value = 0.01, step = 0.005),
      sliderInput("p_age_slope", "Increase in ignition per year of age:",
                  min = 0, max = 0.02, value = 0.005, step = 0.001),
      sliderInput("p_spread", "Neighbor spread probability:",
                  min = 0, max = 1, value = 0.3, step = 0.05),

      tags$hr(),
      sliderInput("sim_years", "Years to Simulate (per run):",
                  min = 1, max = 250, value = 50),
      sliderInput("b0_frac", HTML("B\u2080 as fraction of K (post-disturbance):"),
                  min = 0.01, max = 0.5, value = 0.05, step = 0.01),
      checkboxInput("b0_jitter", "Randomize B\u2080 ±20%", value = TRUE),

      tags$hr(),
      actionButton("run_sim", "Run Simulation", class = "btn-success"),
      actionButton("reset_sim", "Reset Simulation", class = "btn-danger")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Model Overview",
                 h3("Simplified Fire Spread Model"),
                 p("Each patch (cell) grows logistically. Older stands are more likely to ignite,
             representing the buildup of fuel with age. Fires can spread to neighboring patches
             with a fixed probability, creating clustered disturbances."),
                 tags$hr(),
                 h4("Ignition Probability"),
                 withMathJax("
          \\[
          p_{\\text{ignite}} = \\min(1,\\ p_{\\text{base}} + s_{\\text{age}} \\times \\text{Age})
          \\]
          "),
                 h4("Fire Spread"),
                 withMathJax("
          \\[
          P_{\\text{spread}} = p_{\\text{spread}}
          \\]
          "),
                 p("If a cell burns, each of its four immediate neighbors (N, S, E, W)
             may also ignite with probability \\(p_{\\text{spread}}\\).
             Burned cells reset their age and biomass, then regrow after one year.")
        ),

        tabPanel("Grid View", plotOutput("forest_plot", height = "580px")),
        tabPanel("Age Distribution", plotOutput("age_dist_plot", height = "420px")),
        tabPanel("Total Biomass", plotOutput("biomass_plot", height = "420px"))
      )
    )
  )
)

server <- function(input, output, session) {

  # Helper to shift matrices for neighbor spread
  shift_mat <- function(M, dx, dy, fill = FALSE) {
    n1 <- nrow(M); n2 <- ncol(M)
    out <- matrix(fill, n1, n2)
    x_src <- (max(1, 1 - dx)):(min(n1, n1 - dx))
    y_src <- (max(1, 1 - dy)):(min(n2, n2 - dy))
    x_dst <- x_src + dx
    y_dst <- y_src + dy
    out[x_dst, y_dst] <- M[x_src, y_src]
    out
  }

  # Reactive state
  state <- reactiveValues(
    B = NULL, lag = NULL, age = NULL,
    total_biomass = numeric(0), current_year = 0
  )

  # Initialize forest
  initialize_forest <- function(n, K) {
    list(
      B = matrix(0.1 * K, n, n),
      lag = matrix(0L, n, n),
      age = matrix(0L, n, n)
    )
  }

  # Reset
  observeEvent(input$reset_sim, {
    n <- input$grid_size; K <- input$K
    init <- initialize_forest(n, K)
    state$B <- init$B; state$lag <- init$lag; state$age <- init$age
    state$total_biomass <- numeric(0); state$current_year <- 0
  }, ignoreInit = TRUE)

  # Startup init
  observeEvent(TRUE, {
    n <- input$grid_size; K <- input$K
    init <- initialize_forest(n, K)
    state$B <- init$B; state$lag <- init$lag; state$age <- init$age
    state$total_biomass <- numeric(0); state$current_year <- 0
  }, once = TRUE)

  # Run simulation
  observeEvent(input$run_sim, {
    req(state$B)
    n <- input$grid_size
    r <- input$r; K <- input$K
    years <- input$sim_years
    p_base <- input$p_base
    p_slope <- input$p_age_slope
    p_spread <- input$p_spread
    b0_frac <- input$b0_frac; jitter_b0 <- input$b0_jitter
    lag_years <- 1L

    B <- state$B; lag <- state$lag; age <- state$age
    total_biomass <- state$total_biomass

    for (t in 1:years) {
      age <- age + 1L

      # Logistic growth for active cells
      can_grow <- (lag == 0L)
      if (any(can_grow)) {
        B[can_grow] <- B[can_grow] + r * B[can_grow] * (1 - B[can_grow]/K)
        B[can_grow] <- pmin(pmax(B[can_grow], 0), K)
      }

      # 🔥 Step 1: ignition increases with age
      p_ignite <- pmin(1, p_base + p_slope * age)
      ignitions <- matrix(runif(n*n) < p_ignite, n, n)

      # 🔥 Step 2: neighbor spread (4-neighbor)
      spread_mask <- matrix(FALSE, n, n)
      for (off in list(c(1,0), c(-1,0), c(0,1), c(0,-1))) {
        spread_mask <- spread_mask | (shift_mat(ignitions, off[1], off[2]) & (runif(n*n) < p_spread))
      }

      burn_mask <- ignitions | spread_mask

      # Apply burns
      if (any(burn_mask)) {
        B[burn_mask] <- 0
        lag[burn_mask] <- lag_years
        age[burn_mask] <- 0L
      }

      # Decrement lag and reseed
      finished_lag <- (lag > 0L) & ((lag - 1L) == 0L)
      lag[lag > 0L] <- lag[lag > 0L] - 1L
      if (any(finished_lag)) {
        B0 <- b0_frac * K
        if (jitter_b0) {
          J <- runif(sum(finished_lag), 0.8, 1.2)
          B[finished_lag] <- pmax(1e-8, B0 * J)
        } else {
          B[finished_lag] <- pmax(1e-8, B0)
        }
      }

      total_biomass <- c(total_biomass, sum(B))
    }

    state$B <- B; state$lag <- lag; state$age <- age
    state$total_biomass <- total_biomass
    state$current_year <- state$current_year + years
  })

  # --- Plots ---
  output$forest_plot <- renderPlot({
    req(state$B)
    n <- nrow(state$B)
    df <- expand.grid(x = 1:n, y = 1:n)
    df$B <- as.vector(state$B)
    df$age <- as.vector(state$age)

    ggplot(df, aes(x, y, fill = B)) +
      geom_tile(color = "gray90", size = 0.25) +
      geom_text(aes(label = age), size = ifelse(n > 25, 1.8, 2.4), color = "gray20") +
      scale_fill_gradientn(colors = c("red", "orange", "yellow", "darkgreen"),
                           name = "Biomass", limits = c(0, input$K)) +
      coord_fixed() +
      labs(title = paste0("Biomass and Age (Years) — Year ", state$current_year),
           x = NULL, y = NULL) +
      theme_minimal(base_size = 14) +
      theme(axis.text = element_blank(), panel.grid = element_blank())
  })

  output$age_dist_plot <- renderPlot({
    req(state$age)
    df <- data.frame(Age = as.vector(state$age))
    ggplot(df, aes(x = Age)) +
      geom_histogram(aes(y = after_stat(density)), binwidth = 1,
                     fill = "steelblue", color = "white", alpha = 0.6) +
      geom_density(linewidth = 1.2, color = "darkblue") +
      theme_minimal(base_size = 14) +
      labs(title = paste0("Patch Age Distribution — Year ", state$current_year),
           x = "Years Since Last Disturbance", y = "Probability Density")
  })

  output$biomass_plot <- renderPlot({
    req(state$total_biomass)
    df <- data.frame(Year = seq_along(state$total_biomass),
                     TotalBiomass = state$total_biomass)
    ggplot(df, aes(Year, TotalBiomass)) +
      geom_line(linewidth = 1.2, color = "forestgreen") +
      theme_minimal(base_size = 14) +
      labs(title = "Total Biomass Over Time", y = "Sum of Biomass", x = "Year")
  })
}

shinyApp(ui, server)
