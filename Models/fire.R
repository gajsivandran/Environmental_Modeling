# app.R — Forest Disturbance Model with Monte Carlo Description in Overview
# Simulation up to 250 years; Monte Carlo replicates slider inside its tab.

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(shinythemes)

ui <- fluidPage(
  theme = shinytheme("flatly"),
  titlePanel("Forest Disturbance Simulation — Monte Carlo Exploration"),

  sidebarLayout(
    sidebarPanel(
      sliderInput("grid_size", "Grid Size (NxN):",
                  min = 5, max = 50, value = 20),
      sliderInput("r", "Intrinsic Growth Rate (r):",
                  min = 0.01, max = 0.6, value = 0.2, step = 0.01),
      sliderInput("K", "Carrying Capacity (K):",
                  min = 10, max = 500, value = 100, step = 10),
      sliderInput("disturb_prob", "Annual Disturbance Probability:",
                  min = 0, max = 0.1, value = 0.05, step = 0.01),

      # ✅ Allow long-term simulations up to 250 years
      sliderInput("sim_years", "Years to Simulate (per run):",
                  min = 1, max = 250, value = 50),

      tags$hr(),
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
                 h3("Forest Disturbance Model Overview"),
                 p("This model represents a forest as a grid of patches (cells),
                    each growing logistically until a disturbance (e.g., fire) resets its biomass.
                    The model captures both deterministic growth and stochastic disturbance effects."),

                 tags$hr(),
                 h4("Logistic Growth Equation"),
                 withMathJax("
                 \\[
                 B(t) = \\frac{K}{1 + A e^{-r t}}
                 \\]
                 "),
                 p("where \\(B(t)\\) is biomass at time \\(t\\), \\(K\\) is carrying capacity,
                    and \\(r\\) is intrinsic growth rate."),

                 h4("Differential Form"),
                 withMathJax("
                 \\[
                 \\frac{dB}{dt} = rB\\left(1 - \\frac{B}{K}\\right)
                 \\]
                 "),

                 h4("Disturbance Process"),
                 withMathJax("
                 \\[
                 D_{i,j}(t) \\sim \\text{Bernoulli}(p_{disturb})
                 \\]
                 "),
                 p("Each year, each patch \\((i,j)\\) has a probability \\(p_{disturb}\\)
                    of burning. Disturbed cells reset biomass to a small fraction of \\(K\\)
                    and restart growth. The system’s behavior emerges from many such stochastic events."),

                 tags$hr(),
                 h4("Monte Carlo Simulation"),
                 p("Because disturbances occur randomly, every model run produces a slightly different outcome.
                    The Monte Carlo tab repeats the same simulation many times, each with new random events,
                    to estimate the expected trajectory and variability of total forest biomass over time."),

                 withMathJax("
                 \\[
                 E[B(t)] \\approx \\frac{1}{N} \\sum_{i=1}^{N} B_i(t)
                 \\]
                 "),
                 p("where \\(B_i(t)\\) is the biomass from the \\(i^{th}\\) simulation and \\(N\\) is the number of replicates."),
                 p("The result is a mean biomass curve with uncertainty bands that reflect stochastic variability
                    — illustrating how repeated random disturbances can shape long-term forest dynamics.")
        ),

        tabPanel("Grid View", plotOutput("forest_plot", height = "580px")),
        tabPanel("Age Distribution", plotOutput("age_dist_plot", height = "420px")),
        tabPanel("Total Biomass", plotOutput("biomass_plot", height = "420px")),

        # ✅ Monte Carlo Tab with internal slider
        tabPanel("Monte Carlo Simulation",
                 br(),
                 sliderInput("n_reps", "Monte Carlo Replicates:",
                             min = 1, max = 100, value = 20, step = 1),
                 actionButton("run_montecarlo", "Run Monte Carlo Simulation", class = "btn-primary"),
                 br(), br(),
                 plotOutput("montecarlo_plot", height = "450px"))
      )
    )
  )
)

server <- function(input, output, session) {

  # Store ongoing simulation state
  state <- reactiveValues(
    B = NULL,
    lag = NULL,
    age = NULL,
    total_biomass = numeric(0),
    current_year = 0
  )

  # Function to initialize forest
  initialize_forest <- function(n, K) {
    list(
      B   = matrix(0.1 * K, n, n),
      lag = matrix(0L, n, n),
      age = matrix(0L, n, n)
    )
  }

  # Function for one realization
  simulate_forest <- function(n, r, K, p_disturb, years, b0_frac, jitter_b0) {
    B <- matrix(0.1 * K, n, n)
    lag <- matrix(0L, n, n)
    total_biomass <- numeric(0)
    age <- matrix(0L, n, n)
    lag_years <- 1L

    for (t in 1:years) {
      age <- age + 1L

      can_grow <- (lag == 0L)
      if (any(can_grow)) {
        B[can_grow] <- B[can_grow] + r * B[can_grow] * (1 - B[can_grow] / K)
        B[can_grow] <- pmax(0, pmin(B[can_grow], K))
      }

      disturb_mask <- matrix(runif(n * n) < p_disturb, n, n)
      if (any(disturb_mask)) {
        B[disturb_mask] <- 0
        lag[disturb_mask] <- lag_years
        age[disturb_mask] <- 0L
      }

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
    total_biomass
  }

  # Reset simulation
  observeEvent(input$reset_sim, {
    init <- initialize_forest(input$grid_size, input$K)
    state$B <- init$B
    state$lag <- init$lag
    state$age <- init$age
    state$total_biomass <- numeric(0)
    state$current_year <- 0
  }, ignoreInit = TRUE)

  # Initialize on startup
  observeEvent(TRUE, {
    init <- initialize_forest(input$grid_size, input$K)
    state$B <- init$B
    state$lag <- init$lag
    state$age <- init$age
    state$total_biomass <- numeric(0)
    state$current_year <- 0
  }, once = TRUE)

  # Run single simulation (Grid/Age)
  observeEvent(input$run_sim, {
    req(state$B)
    n <- input$grid_size
    r <- input$r
    K <- input$K
    p_disturb <- input$disturb_prob
    years <- input$sim_years
    b0_frac <- input$b0_frac
    jitter_b0 <- input$b0_jitter
    lag_years <- 1L

    B <- state$B
    lag <- state$lag
    age <- state$age
    total_biomass <- state$total_biomass

    for (t in 1:years) {
      age <- age + 1L
      can_grow <- (lag == 0L)
      if (any(can_grow)) {
        B[can_grow] <- B[can_grow] + r * B[can_grow] * (1 - B[can_grow] / K)
        B[can_grow] <- pmax(0, pmin(B[can_grow], K))
      }

      disturb_mask <- matrix(runif(n * n) < p_disturb, n, n)
      if (any(disturb_mask)) {
        B[disturb_mask] <- 0
        lag[disturb_mask] <- lag_years
        age[disturb_mask] <- 0L
      }

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

    state$B <- B
    state$lag <- lag
    state$age <- age
    state$total_biomass <- total_biomass
    state$current_year <- state$current_year + years
  })

  # --- Grid View ---
  output$forest_plot <- renderPlot({
    req(state$B)
    n <- nrow(state$B)
    df <- expand.grid(x = 1:n, y = 1:n)
    df$B <- as.vector(state$B)
    df$age <- as.vector(state$age)

    ggplot(df, aes(x, y, fill = B)) +
      geom_tile(color = "gray90", size = 0.25) +
      geom_text(aes(label = age), size = ifelse(n > 25, 1.8, 2.5), color = "gray20") +
      scale_fill_gradientn(colors = c("red", "orange", "yellow", "darkgreen"),
                           name = "Biomass", limits = c(0, input$K)) +
      coord_fixed() +
      labs(title = paste0("Biomass and Age (Years) — Year ", state$current_year),
           x = NULL, y = NULL) +
      theme_minimal(base_size = 14) +
      theme(axis.text = element_blank(), panel.grid = element_blank())
  })

  # --- Age Distribution ---
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

  # --- Total Biomass ---
  output$biomass_plot <- renderPlot({
    req(state$total_biomass)
    df <- data.frame(Year = seq_along(state$total_biomass),
                     TotalBiomass = state$total_biomass)
    ggplot(df, aes(Year, TotalBiomass)) +
      geom_line(linewidth = 1.2, color = "forestgreen") +
      theme_minimal(base_size = 14) +
      labs(title = "Total Biomass Over Time", y = "Sum of Biomass", x = "Year")
  })

  # --- Monte Carlo Simulation ---
  observeEvent(input$run_montecarlo, {
    n <- input$grid_size
    r <- input$r
    K <- input$K
    p_disturb <- input$disturb_prob
    years <- input$sim_years
    b0_frac <- input$b0_frac
    jitter_b0 <- input$b0_jitter
    n_reps <- input$n_reps

    results <- replicate(
      n_reps,
      simulate_forest(n, r, K, p_disturb, years, b0_frac, jitter_b0)
    )

    df <- as.data.frame(results)
    df$Year <- 1:years
    df_long <- df %>%
      pivot_longer(-Year, names_to = "Run", values_to = "Biomass")

    summary_df <- df_long %>%
      group_by(Year) %>%
      summarise(mean = mean(Biomass),
                sd = sd(Biomass),
                .groups = "drop")

    output$montecarlo_plot <- renderPlot({
      ggplot() +
        geom_line(data = df_long, aes(Year, Biomass, group = Run),
                  color = "gray70", alpha = 0.4) +
        geom_ribbon(data = summary_df,
                    aes(Year, ymin = mean - sd, ymax = mean + sd),
                    fill = "palegreen3", alpha = 0.4) +
        geom_line(data = summary_df, aes(Year, mean),
                  color = "darkgreen", linewidth = 1.2) +
        theme_minimal(base_size = 14) +
        labs(title = paste0("Monte Carlo Biomass Simulation (", n_reps, " runs)"),
             y = "Total Biomass", x = "Year")
    })
  })
}

shinyApp(ui, server)
