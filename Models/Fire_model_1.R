# ============================================================
# Forest Dynamics — Part 1: Markov Chains & Monte Carlo
# ============================================================

library(shiny)
library(tidyverse)
library(plotly)

# ------------------------------------------------------------
# GLOBAL SETTINGS
# ------------------------------------------------------------
GRID_SIZE <- 20
MAX_AGE   <- 50

# Color palette by age
age_palette <- colorRampPalette(c("lightyellow", "darkgreen"))
age_to_color <- function(age) {
  n_colors <- MAX_AGE + 1
  pal <- age_palette(n_colors)
  pal[pmin(age, MAX_AGE) + 1]
}

# ------------------------------------------------------------
# UI
# ------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Forest Dynamics — Part 1: Markov Chains & Monte Carlo"),

  tabsetPanel(
    id = "main_tabs",

    # ========================================================
    # TAB 0 — MODEL OVERVIEW
    # ========================================================
    tabPanel(
      title = "Model Overview",
      br(),
      h3("Overview of the Forest Growth and Fire Model"),
      p("This model simulates the long-term dynamics of a forest landscape represented as a grid of individual patches (cells)."),
      p("Each patch follows simple stochastic rules for growth and fire disturbance, demonstrating how a Markov process and Monte Carlo simulation can represent ecological succession."),

      hr(),
      h4("1. State Representation"),
      p("Each cell has an integer-valued state representing its stand age (in years)."),
      tags$ul(
        tags$li("Ages range from 0 (bare or recently burned) to 50 (old growth)."),
        tags$li("Each year, all cells age by 1, unless they burn and reset to age 0.")
      ),

      hr(),
      h4("2. Transition Dynamics"),
      p("Each patch transitions independently each year according to these rules:"),
      tags$ol(
        tags$li("Age increases by one year (up to a maximum of 50)."),
        tags$li("Fire may occur, resetting the age to zero."),
        tags$li("Fire probability depends on both a user-controlled base rate and age-dependent fuel buildup.")
      ),

      hr(),
      h4("3. Fire Probability Equation"),
      withMathJax(),
      p("The probability that a given cell burns in a given year is:"),
      helpText("$$ P(\\text{burn}) = \\min( p_{\\text{base}} + 0.02 \\times \\frac{\\text{age}}{50}, 0.3 ) $$"),
      p("where:"),
      tags$ul(
        tags$li("\\( p_{\\text{base}} \\): user-defined base probability of ignition (slider control)"),
        tags$li("Age term increases flammability as fuel accumulates"),
        tags$li("0.3 is the upper bound to avoid unrealistic full-burn scenarios")
      ),

      hr(),
      h4("4. Transition Matrix (Conceptual Form)"),
      p("If we treat each patch as a Markov chain with discrete age states \\(0, 1, 2, ..., 50\\), the transition probabilities can be expressed as:"),
      helpText("$$ P_{a,a+1} = 1 - P(\\text{burn at age } a) $$"),
      helpText("$$ P_{a,0} = P(\\text{burn at age } a) $$"),
      p("All other transitions are 0. Once at age 50, patches remain there unless burned."),

      hr(),
      h4("5. Model Assumptions"),
      tags$ul(
        tags$li("Each cell acts independently (no spatial fire spread between neighbors)."),
        tags$li("Fire resets age to 0 instantly, with no partial damage or lag recovery."),
        tags$li("Age increases deterministically (1 year per step) unless burned."),
        tags$li("All cells share identical parameters for growth and ignition."),
        tags$li("No climate feedbacks or fuel moisture effects — fire probability is static except for the slider-controlled base rate."),
        tags$li("Landscape-level dynamics emerge statistically from many independent stochastic patches (Monte Carlo principle).")
      ),

      hr(),
      h4("6. Interpretation"),
      p("This model represents the simplest possible stochastic disturbance–succession cycle. Despite its simplicity, it reveals how changing ignition pressure or fuel buildup rules can shift the steady-state age structure of a forest."),
      p("In later parts, we’ll add spatial spread, weather variability, and feedbacks between climate and fire.")
    ),

    # ========================================================
    # TAB 1 — PROCEDURE
    # ========================================================
    tabPanel(
      title = "Procedure",
      br(),
      fluidRow(
        column(
          width = 6,
          h4("Interactive 20×20 Forest Grid"),
          plotOutput("forest_plot", height = "450px"),
          br(),
          fluidRow(
            column(6,
                   sliderInput("fire_base", "Base Fire Probability:",
                               min = 0.0, max = 0.05, value = 0.01, step = 0.005)),
            column(6,
                   actionButton("step1", "Step 1 Year", class = "btn-success"),
                   actionButton("step10", "Step 10 Years", class = "btn-primary"))
          ),
          br(),
          wellPanel(
            h5("What You See:"),
            tags$ul(
              tags$li("Each square = one forest patch (cell)."),
              tags$li("Color scale: yellow → green indicates stand age (0–50 years)."),
              tags$li("Red = patch that burned this year."),
              tags$li("Outlined cells = 5 tracked patches shown in the right-hand plot.")
            )
          )
        ),

        column(
          width = 6,
          h4("Simulation Outputs"),
          plotlyOutput("composition_plot", height = "250px"),
          br(),
          plotlyOutput("trajectory_plot", height = "250px"),
          br(),
          wellPanel(
            h5("Variable Definitions"),
            tags$ul(
              tags$li(tags$b("mean_age:"), " average age of all cells in that year."),
              tags$li(tags$b("prop_young:"), " proportion of patches younger than 10 years."),
              tags$li(tags$b("prop_old:"), " proportion of patches older than 10 years (1 - prop_young).")
            )
          )
        )
      )
    ),

    # ========================================================
    # TAB 2 — SIMULATION
    # ========================================================
    tabPanel(
      title = "Simulation",
      sidebarLayout(
        sidebarPanel(
          h4("Monte Carlo Simulation"),
          numericInput("n_patches", "Number of patches:", 1000, min = 10, max = 10000, step = 100),
          numericInput("n_years", "Number of years:", 100, min = 10, max = 500, step = 10),
          actionButton("run_btn", "Run Simulation", class = "btn-success"),
          br(), br(),
          helpText("This runs many independent forest patches to show average dynamics.")
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Proportion Through Time", plotlyOutput("time_plot")),
            tabPanel("Final-Year Composition", plotlyOutput("pie_plot"))
          )
        )
      )
    )
  )
)

# ------------------------------------------------------------
# SERVER
# ------------------------------------------------------------
server <- function(input, output, session) {

  # ==========================================================
  # TAB 1 — PROCEDURE SIMULATION
  # ==========================================================

  forest <- reactiveVal(matrix(MAX_AGE, nrow = GRID_SIZE, ncol = GRID_SIZE))
  burned <- reactiveVal(matrix(0, nrow = GRID_SIZE, ncol = GRID_SIZE))
  current_year <- reactiveVal(0)

  history_data <- reactiveVal(tibble(year = numeric(), mean_age = numeric(), prop_young = numeric()))

  observe({
    ages <- isolate(forest())
    history_data(
      tibble(
        year = 0,
        mean_age = mean(ages),
        prop_young = mean(ages < 10)
      )
    )
  })

  set.seed(1)
  tracked_cells <- sample(1:(GRID_SIZE^2), 5)
  tracked_history <- reactiveVal(
    tibble(cell = tracked_cells, year = 0, age = rep(MAX_AGE, length(tracked_cells)))
  )

  step_year <- function(n_steps = 1) {
    ages <- forest()
    fires <- burned()
    hist <- history_data()
    yr <- current_year()
    track <- tracked_history()

    for (s in seq_len(n_steps)) {
      fire_prob <- pmin(input$fire_base + 0.02 * (ages / MAX_AGE), 0.3)
      fire <- matrix(rbinom(length(ages), 1, fire_prob), nrow = GRID_SIZE)
      fires <- fire
      ages <- ifelse(fire == 1, 0, pmin(ages + 1, MAX_AGE))
      yr <- yr + 1
      hist <- bind_rows(
        hist,
        tibble(
          year = yr,
          mean_age = mean(ages),
          prop_young = mean(ages < 10)
        )
      )
      track <- bind_rows(
        track,
        tibble(cell = tracked_cells, year = yr, age = ages[tracked_cells])
      )
    }

    forest(ages)
    burned(fires)
    current_year(yr)
    history_data(hist)
    tracked_history(track)
  }

  observeEvent(input$step1, step_year(1))
  observeEvent(input$step10, step_year(10))

  # ---- Forest grid plot ----
  output$forest_plot <- renderPlot({
    ages <- forest()
    fire <- burned()

    col_mat <- matrix(age_to_color(ages), nrow = GRID_SIZE, ncol = GRID_SIZE)
    col_mat[fire == 1] <- "red"

    image(
      1:GRID_SIZE, 1:GRID_SIZE,
      t(matrix(seq_len(GRID_SIZE^2), GRID_SIZE, GRID_SIZE)[, GRID_SIZE:1]),
      col = as.vector(col_mat),
      axes = FALSE,
      xlab = "", ylab = "",
      main = paste("Year:", current_year())
    )

    tracked_xy <- expand.grid(x = 1:GRID_SIZE, y = 1:GRID_SIZE)[tracked_cells, ]
    for (i in 1:nrow(tracked_xy)) {
      rect(
        tracked_xy$x[i] - 0.5, tracked_xy$y[i] - 0.5,
        tracked_xy$x[i] + 0.5, tracked_xy$y[i] + 0.5,
        border = "black", lwd = 2
      )
    }
  })

  output$composition_plot <- renderPlotly({
    hist <- history_data() |>
      mutate(prop_old = 1 - prop_young) |>
      pivot_longer(cols = starts_with("prop"), names_to = "class", values_to = "prop")

    p <- ggplot(hist, aes(x = year, y = prop, color = class)) +
      geom_line(linewidth = 1) +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(x = "Year", y = "Proportion", title = "Forest Composition Over Time") +
      theme_minimal(base_size = 13)
    ggplotly(p)
  })

  output$trajectory_plot <- renderPlotly({
    df <- tracked_history()
    p <- ggplot(df, aes(x = year, y = age, color = factor(cell))) +
      geom_line(linewidth = 1) +
      labs(
        title = "Ages of 5 Tracked Forest Patches",
        x = "Year", y = "Age", color = "Cell"
      ) +
      theme_minimal(base_size = 13)
    ggplotly(p)
  })

  # ==========================================================
  # TAB 2 — MONTE CARLO SIMULATION
  # ==========================================================
  states <- c("Bare", "Early", "Young", "Mature", "Old")
  default_P <- tribble(
    ~From,    ~Bare, ~Early, ~Young, ~Mature, ~Old,
    "Bare",    0.3,   0.7,    0.0,    0.0,    0.0,
    "Early",   0.0,   0.4,    0.6,    0.0,    0.0,
    "Young",   0.0,   0.0,    0.6,    0.4,    0.0,
    "Mature",  0.0,   0.0,    0.1,    0.7,    0.2,
    "Old",     0.0,   0.0,    0.0,    0.1,    0.9
  )

  simulate_forest <- function(P, n_patches = 1000, years = 100, init_probs = NULL) {
    S <- nrow(P)
    if (is.null(init_probs)) init_probs <- rep(1/S, S)
    cur <- sample.int(S, n_patches, replace = TRUE, prob = init_probs)
    out <- matrix(NA_integer_, nrow = years + 1, ncol = n_patches)
    out[1, ] <- cur
    for (t in 1:years) {
      next_states <- vapply(cur, function(s) sample.int(S, 1, prob = P[s, ]), integer(1))
      out[t + 1, ] <- next_states
      cur <- next_states
    }
    tibble(
      year = rep(0:years, each = n_patches),
      patch = rep(1:n_patches, times = years + 1),
      state = factor(states[out], levels = states)
    )
  }

  sim_data <- eventReactive(input$run_btn, {
    P <- as.matrix(default_P[ , -1])
    rownames(P) <- default_P$From
    simulate_forest(P, n_patches = input$n_patches, years = input$n_years)
  })

  output$time_plot <- renderPlotly({
    df <- sim_data() |>
      count(year, state) |>
      group_by(year) |>
      mutate(prop = n / sum(n))
    p <- ggplot(df, aes(x = year, y = prop, color = state)) +
      geom_line(linewidth = 1) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(
        title = "Forest Composition Over Time (Monte Carlo)",
        x = "Year", y = "Proportion", color = "State"
      ) +
      theme_minimal(base_size = 14)
    ggplotly(p)
  })

  output$pie_plot <- renderPlotly({
    df_final <- sim_data() |>
      filter(year == max(year)) |>
      count(state) |>
      mutate(prop = n / sum(n))
    plot_ly(df_final, labels = ~state, values = ~prop, type = "pie",
            textinfo = "label+percent",
            insidetextorientation = "radial") |>
      layout(title = "Final-Year Composition (Monte Carlo Estimate)")
  })
}

# ------------------------------------------------------------
# RUN APP
# ------------------------------------------------------------
shinyApp(ui, server)
