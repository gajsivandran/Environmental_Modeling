# app.R — Plant Succession (Markov Chain + Monte Carlo)
# Starts fully bare at Year 0 (state 1 everywhere)

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)
library(shinythemes)

ui <- fluidPage(
  theme = shinytheme("flatly"),
  titlePanel("🌿 Plant Succession — Interactive Markov Chain with Monte Carlo"),

  sidebarLayout(
    sidebarPanel(
      width = 5,
      h4("Edit Transition Matrix (P)"),
      p("Rows = current state, columns = next state.
        Rows renormalize automatically to sum = 1."),
      DTOutput("transition_matrix"),

      tags$hr(),
      sliderInput("grid_size", "Grid Size (NxN):",
                  min = 5, max = 50, value = 20),
      sliderInput("sim_years", "Years per Run:",
                  min = 1, max = 250, value = 50),
      sliderInput("disturb_prob", "Disturbance probability (reset → bare):",
                  min = 0, max = 0.2, value = 0.02, step = 0.01),

      tags$hr(),
      actionButton("run_sim", "Run Simulation", class = "btn-success"),
      actionButton("reset_sim", "Reset Simulation", class = "btn-danger")
    ),

    mainPanel(
      width = 7,
      tabsetPanel(
        tabPanel("Model Overview",
                 h3("1️⃣ Ecological Succession as a Markov Chain"),
                 p("In this model, each grid cell represents a patch of land that can exist
             in one of several successional stages — from bare ground to mature forest.
             Over time, each patch transitions between stages with certain probabilities."),

                 h4("Mathematical Formulation"),
                 p("A **Markov chain** is a stochastic process that moves through a set of discrete states
             according to fixed transition probabilities.
             The key assumption, known as the *Markov property*, is that the future depends
             only on the current state — not on the past history of states."),
                 withMathJax("
          \\[
          P =
          \\begin{bmatrix}
          p_{11} & p_{12} & p_{13} & p_{14} & p_{15} \\\\
          p_{21} & p_{22} & p_{23} & p_{24} & p_{25} \\\\
          p_{31} & p_{32} & p_{33} & p_{34} & p_{35} \\\\
          p_{41} & p_{42} & p_{43} & p_{44} & p_{45} \\\\
          p_{51} & p_{52} & p_{53} & p_{54} & p_{55}
          \\end{bmatrix}
          \\]
          "),
                 p("Each row of the transition matrix \\( P \\) corresponds to the *current* state,
             and each column gives the probability of moving to the *next* state."),

                 h4("State Vector Dynamics"),
                 p("Let \\( \\mathbf{x}_t = [x_{1,t}, x_{2,t}, ..., x_{5,t}] \\) be the proportion of the landscape
             in each state at time \\( t \\).
             The dynamics of succession across the whole landscape are given by matrix multiplication:"),
                 withMathJax("
          \\[
          \\mathbf{x}_{t+1} = \\mathbf{x}_t P
          \\]
          "),
                 p("This equation predicts how the composition of the landscape evolves year by year.
             By iterating this relationship, we can trace the long-term successional trajectory."),

                 h4("Stationary Distribution"),
                 p("If the transition matrix remains constant, the system often converges to a
             **stationary distribution** \\( \\pi \\), where the proportions of each state stop changing:"),
                 withMathJax("
          \\[
          \\pi P = \\pi
          \\]
          "),
                 p("Ecologically, this represents a *dynamic equilibrium*:
             even though individual patches keep changing, the overall composition of the landscape stays stable."),

                 tags$hr(),
                 h3("2️⃣ Disturbance and Stochasticity"),
                 p("Real ecosystems experience disturbances such as fire, flooding, or disease that can reset a patch
             to the early (bare ground) state.
             In this model, each patch has a small probability of disturbance each year,
             defined by the **Disturbance Probability** slider in the sidebar."),

                 p("This addition introduces stochastic behavior that prevents the system
             from reaching a perfectly stable equilibrium,
             instead maintaining a constantly shifting mosaic of successional stages."),

                 tags$hr(),
                 h3("3️⃣ Monte Carlo Simulations"),
                 p("Even with the same transition probabilities, random outcomes differ between runs.
             To capture this uncertainty, the **Monte Carlo** tab repeats the simulation
             many times (replicates) and averages the results."),
                 withMathJax("
          \\[
          \\hat{E}[x_{i,t}] = \\frac{1}{N} \\sum_{r=1}^{N} x_{i,t}^{(r)}
          \\]
          "),
                 p("where \\( x_{i,t}^{(r)} \\) is the proportion of the landscape in state \\( i \\)
             at time \\( t \\) during replicate \\( r \\),
             and \\( N \\) is the total number of Monte Carlo replicates."),
                 p("The shaded region in the Monte Carlo plot represents the range
             (minimum–maximum) of outcomes across replicates,
             while the solid line shows the mean trajectory."),
                 tags$hr(),
                 h4("Summary"),
                 tags$ul(
                   tags$li("Each grid cell evolves independently according to the Markov matrix \\( P \\)."),
                   tags$li("At Year 0, all patches begin as **bare ground** (state 1)."),
                   tags$li("Disturbance events can reset patches back to bare."),
                   tags$li("Monte Carlo replicates reveal uncertainty in the successional outcomes.")
                 )
        ),

        tabPanel("Grid View", plotOutput("succession_plot", height = "580px")),
        tabPanel("State Distribution", plotOutput("state_dist_plot", height = "420px")),
        tabPanel("Monte Carlo",
                 sliderInput("n_reps", "Monte Carlo Replicates:",
                             min = 5, max = 25, value = 10, step = 5),
                 actionButton("run_mc", "Run Monte Carlo", class = "btn-primary"),
                 br(), br(),
                 plotOutput("monte_plot", height = "450px"),
                 p("Each replicate simulates the same landscape size and transition matrix P.
            The ribbon shows the range across replicates, and the line shows the mean proportion.")
        )
      )
    )
  )
)

server <- function(input, output, session) {

  # ---- Setup ----
  state_colors <- c("saddlebrown", "lightgreen", "darkolivegreen", "forestgreen", "darkgreen")

  base_P <- matrix(c(
    0.70, 0.30, 0.00, 0.00, 0.00,
    0.05, 0.70, 0.25, 0.00, 0.00,
    0.00, 0.10, 0.70, 0.20, 0.00,
    0.00, 0.00, 0.15, 0.70, 0.15,
    0.05, 0.00, 0.00, 0.10, 0.85
  ), nrow = 5, byrow = TRUE)
  P_values <- reactiveVal(base_P)

  # ---- Matrix editor ----
  output$transition_matrix <- renderDT({
    datatable(
      as.data.frame(round(P_values(), 2)),
      rownames = paste("State", 1:5),
      colnames = paste("→", 1:5),
      editable = TRUE,
      options = list(dom = "t", paging = FALSE, ordering = FALSE)
    )
  })

  observeEvent(input$transition_matrix_cell_edit, {
    info <- input$transition_matrix_cell_edit
    new_P <- P_values()
    new_P[info$row, info$col] <- as.numeric(info$value)
    new_P <- t(apply(new_P, 1, function(r) if (sum(r) > 0) r / sum(r) else rep(0, length(r))))
    P_values(new_P)
  })

  # ---- Simulation helpers ----
  state <- reactiveValues(grid = NULL, current_year = 0, dist_history = NULL)

  # ✅ Start completely bare (state = 1 everywhere)
  initialize_grid <- function(n) {
    matrix(1L, n, n)
  }

  observeEvent(input$reset_sim, {
    n <- input$grid_size
    state$grid <- initialize_grid(n)
    state$current_year <- 0
    state$dist_history <- data.frame()
  }, ignoreInit = TRUE)

  observeEvent(TRUE, {
    n <- input$grid_size
    state$grid <- initialize_grid(n)
    state$current_year <- 0
    state$dist_history <- data.frame()
  }, once = TRUE)

  run_single_sim <- function(P, disturb_p, n, years, grid0 = NULL) {
    grid <- if (is.null(grid0)) initialize_grid(n) else grid0
    dist_hist <- data.frame()

    for (t in 1:years) {
      new_grid <- matrix(NA, n, n)
      for (i in 1:n) {
        for (j in 1:n) {
          new_grid[i,j] <- sample(1:5, 1, prob = P[grid[i,j],])
        }
      }
      disturb_mask <- matrix(runif(n*n) < disturb_p, n, n)
      new_grid[disturb_mask] <- 1

      freqs <- as.numeric(table(factor(new_grid, levels = 1:5))) / (n*n)
      dist_hist <- rbind(dist_hist,
                         data.frame(Year = t,
                                    Bare = freqs[1], Herbs = freqs[2],
                                    Shrubs = freqs[3], Young = freqs[4],
                                    Mature = freqs[5]))
      grid <- new_grid
    }
    dist_hist
  }

  # ---- Single-run simulation ----
  observeEvent(input$run_sim, {
    req(state$grid)
    n <- input$grid_size
    years <- input$sim_years
    P <- P_values()
    disturb_p <- input$disturb_prob

    # Run sim AND keep track of final grid state
    grid <- state$grid
    for (t in 1:years) {
      new_grid <- matrix(NA, n, n)
      for (i in 1:n) {
        for (j in 1:n) {
          new_grid[i, j] <- sample(1:5, 1, prob = P[grid[i, j], ])
        }
      }
      disturb_mask <- matrix(runif(n * n) < disturb_p, n, n)
      new_grid[disturb_mask] <- 1
      grid <- new_grid
    }

    # compute composition over time
    df <- run_single_sim(P, disturb_p, n, years, state$grid)
    df$Year <- df$Year + state$current_year

    # ✅ keep final grid
    state$grid <- grid
    state$current_year <- state$current_year + years
    state$dist_history <- bind_rows(state$dist_history, df)
  })


  # ---- Monte Carlo simulations ----
  observeEvent(input$run_mc, {
    n_reps <- input$n_reps
    n <- input$grid_size
    years <- input$sim_years
    P <- P_values()
    disturb_p <- input$disturb_prob

    replicate_results <- replicate(n_reps, run_single_sim(P, disturb_p, n, years), simplify = FALSE)
    df_all <- bind_rows(replicate_results, .id = "Replicate")
    df_summary <- df_all %>%
      pivot_longer(cols = -c(Replicate, Year), names_to = "State", values_to = "Prop") %>%
      group_by(Year, State) %>%
      summarise(mean = mean(Prop), min = min(Prop), max = max(Prop), .groups = "drop")

    output$monte_plot <- renderPlot({
      ggplot(df_summary, aes(Year, mean, color = State, fill = State)) +
        geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.2, color = NA) +
        geom_line(linewidth = 1.1) +
        scale_color_manual(values = state_colors) +
        scale_fill_manual(values = state_colors) +
        labs(title = paste0("Monte Carlo (", n_reps, " replicates)"),
             y = "Proportion of Landscape", x = "Year") +
        theme_minimal(base_size = 14)
    })
  })

  # ---- Plots ----
  output$succession_plot <- renderPlot({
    req(state$grid)
    n <- nrow(state$grid)
    df <- expand.grid(x = 1:n, y = 1:n)
    df$state <- as.vector(state$grid)

    ggplot(df, aes(x, y, fill = factor(state))) +
      geom_tile(color = "gray80") +
      scale_fill_manual(values = state_colors,
                        labels = c("Bare","Herbs","Shrubs","Young","Mature"),
                        name = "Stage") +
      coord_fixed() +
      labs(title = paste0("Successional States — Year ", state$current_year),
           x = NULL, y = NULL) +
      theme_minimal(base_size = 14) +
      theme(axis.text = element_blank(), panel.grid = element_blank())
  })

  output$state_dist_plot <- renderPlot({
    req(!is.null(state$dist_history) && nrow(state$dist_history) > 0)
    df_long <- state$dist_history %>%
      pivot_longer(cols = c(Bare, Herbs, Shrubs, Young, Mature),
                   names_to = "State", values_to = "Proportion")
    ggplot(df_long, aes(Year, Proportion, color = State)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = state_colors) +
      theme_minimal(base_size = 14) +
      labs(title = "Successional Stage Composition Over Time",
           y = "Proportion of Landscape", x = "Year")
  })
}

shinyApp(ui, server)
