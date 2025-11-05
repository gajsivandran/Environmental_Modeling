# app.R — Mosquito Life Cycle Markov Chain (Fast Multinomial + Full Overview)

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)
library(shinythemes)
library(shinyjs)

ui <- fluidPage(
  useShinyjs(),
  theme = shinytheme("flatly"),
  titlePanel("🦟 Mosquito Life Cycle — Editable Markov Chain"),

  # ---- Custom CSS ----
  tags$head(
    tags$style(HTML("
      table.dataTable tbody td {
        font-size: 13px !important;
        padding: 4px 8px !important;
      }
      table.dataTable thead th {
        font-size: 13px !important;
        padding: 4px 8px !important;
      }
      table.dataTable tbody td input {
        font-size: 14px !important;
        text-align: center !important;
        height: 26px !important;
        padding: 0px 2px !important;
      }
      .dataTables_wrapper { font-size: 13px; }
    "))
  ),

  sidebarLayout(
    sidebarPanel(
      width = 5,
      h4("Edit Transition Matrix (P)"),
      p("Rows = current stage, columns = next stage."),
      p("⚠️ Each row should sum to 1 for a valid Markov chain. If not, a warning will appear below."),
      DTOutput("transition_matrix"),
      uiOutput("row_warning"),

      tags$hr(),
      p("Initial population fixed at 100 individuals."),
      sliderInput("sim_days", "Days per Run:", min = 1, max = 150, value = 40),
      sliderInput("mortality_prob", "Environmental mortality (reset → Dead):",
                  min = 0, max = 0.1, value = 0.01, step = 0.005),
      sliderInput("repro_rate", "Adult Reproduction Rate (Eggs per Adult per Day):",
                  min = 0, max = 2, value = 0.5, step = 0.1),

      radioButtons("view_mode", "Display Mode:",
                   choices = c("Proportion" = "prop", "Count" = "count"),
                   selected = "prop", inline = TRUE),

      tags$hr(),
      actionButton("run_sim", "Run Simulation", class = "btn-success")
    ),

    mainPanel(
      width = 7,
      tabsetPanel(
        tabPanel("Overview",
                 h3("Mosquito Life Cycle as a Markov Chain"),
                 p("This model represents the daily transitions between different life stages of a mosquito population
             — from egg to larva to pupa to adult — using a probabilistic framework called a Markov chain.
             Each day, individuals have certain probabilities of surviving, developing to the next stage, or dying."),

                 h4("1️⃣ Model Structure"),
                 p("The system is composed of five discrete states:"),
                 tags$ul(
                   tags$li(strong("Egg:"), " newly laid eggs that may hatch into larvae."),
                   tags$li(strong("Larva:"), " aquatic feeding stage before pupation."),
                   tags$li(strong("Pupa:"), " transitional, non-feeding stage leading to emergence."),
                   tags$li(strong("Adult:"), " reproductive, flying stage that lays new eggs."),
                   tags$li(strong("Dead:"), " terminal absorbing state — once entered, mosquitoes remain dead.")
                 ),
                 p("Each state represents a group (or cohort) of individuals, and transitions occur each day according to
             a set of probabilities stored in a transition matrix."),

                 h4("2️⃣ The Transition Matrix"),
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
                 p("Each row corresponds to the *current* state, and each column represents the probability of moving to
             a *next* state on the following day. For example, \\( p_{23} \\) is the probability that a larva
             becomes a pupa during one day."),

                 p("Because the model is a Markov chain, the system is assumed to have no memory:
             the probability of transitioning depends only on the current state, not on how the mosquito arrived there."),

                 h4("3️⃣ Population Dynamics"),
                 p("Let \\( \\mathbf{x}_t = [x_{Egg,t}, x_{Larva,t}, x_{Pupa,t}, x_{Adult,t}, x_{Dead,t}] \\)
             be the number of individuals in each state at day \\( t \\). The population evolves through matrix multiplication:"),
                 withMathJax("
            \\[
            \\mathbf{x}_{t+1} = \\mathbf{x}_t P
            \\]
          "),
                 p("This equation predicts the *expected number* of individuals in each stage over time.
             In this app, stochasticity is introduced using multinomial sampling to make the simulation more realistic —
             small random differences accumulate across runs, mimicking natural variation in survival and development."),

                 h4("4️⃣ Mortality and Reproduction"),
                 p("At each time step:"),
                 tags$ul(
                   tags$li("A small fraction of each living stage dies, controlled by the ", strong("Environmental Mortality"), " parameter."),
                   tags$li("Adults lay eggs at a rate set by ", strong("Adult Reproduction Rate"), ", adding new individuals to the Egg stage.")
                 ),
                 p("Unlike succession models, adults here do not transition back to eggs directly — instead, reproduction creates
             new eggs externally, while dead individuals accumulate in the absorbing state."),

                 h4("5️⃣ Monte Carlo Interpretation"),
                 p("Although this app runs a single simulation, each run is equivalent to one realization of a stochastic process.
             Running it multiple times (with the same parameters) will produce slightly different trajectories,
             reflecting uncertainty in biological processes like survival, development time, and reproduction."),

                 h4("6️⃣ Ecological Interpretation"),
                 p("By adjusting the transition probabilities, you can simulate different environmental conditions:
             cooler temperatures slow development (reducing progression probabilities), while high mortality mimics predation
             or pesticide effects. Increasing reproduction rate leads to faster population growth and potentially
             cyclical boom-and-bust patterns if mortality is low."),
                 p("This Markov chain model provides a simplified but powerful framework for understanding how small daily probabilities
             translate into large-scale population trends in mosquito ecology.")
        ),

        tabPanel("Stage Dynamics", plotOutput("stage_plot", height = "480px")),
        tabPanel("Population Size", plotOutput("pop_plot", height = "420px"))
      )
    )
  )
)

server <- function(input, output, session) {

  # ---- Setup ----
  state_labels <- c("Egg", "Larva", "Pupa", "Adult", "Dead")
  base_P <- matrix(c(
    0.10, 0.85, 0.00, 0.00, 0.05,
    0.00, 0.50, 0.40, 0.00, 0.10,
    0.00, 0.00, 0.60, 0.35, 0.05,
    0.00, 0.00, 0.00, 0.90, 0.10,
    0.00, 0.00, 0.00, 0.00, 1.00
  ), nrow = 5, byrow = TRUE)

  P_values <- reactiveVal(base_P)

  # ---- Validation helper ----
  validate_rows <- function(P) {
    row_sums <- rowSums(P)
    deviations <- abs(row_sums - 1)
    bad_rows <- which(deviations > 0.01)
    list(bad_rows = bad_rows, row_sums = row_sums)
  }

  # ---- Editable Matrix ----
  output$transition_matrix <- renderDT({
    df <- as.data.frame(round(P_values(), 2))
    rownames(df) <- state_labels
    colnames(df) <- paste("→", state_labels)
    datatable(
      df,
      editable = TRUE,
      rownames = TRUE,
      options = list(
        dom = "t",
        paging = FALSE,
        ordering = FALSE,
        searching = FALSE,
        autoWidth = TRUE
      )
    )
  }, server = TRUE)

  proxy <- dataTableProxy("transition_matrix")

  observeEvent(input$transition_matrix_cell_edit, {
    info <- input$transition_matrix_cell_edit
    i <- info$row; j <- info$col
    v <- as.numeric(info$value)
    if (!is.na(v)) {
      v <- round(v / 0.05) * 0.05
      new_P <- P_values()
      new_P[i, j] <- v
      P_values(new_P)
      replaceData(proxy, as.data.frame(round(new_P, 2)), resetPaging = FALSE)
    }
  })

  output$row_warning <- renderUI({
    val <- validate_rows(P_values())
    if (length(val$bad_rows) > 0) {
      bad_names <- paste(state_labels[val$bad_rows], collapse = ", ")
      div(style = "color:#b30000; font-weight:bold; margin-top:10px;",
          paste("⚠️ Warning: Row(s)", bad_names, "do not sum to 1."))
    } else {
      div(style = "color:darkgreen; font-weight:bold; margin-top:10px;",
          "✅ All rows sum to 1.")
    }
  })

  # ---- Fast stochastic simulation ----
  run_single_sim <- function(P, mortality_p, repro_rate, days) {
    stages <- c("Egg", "Larva", "Pupa", "Adult", "Dead")
    n_stages <- length(stages)
    pop <- c(100, 0, 0, 0, 0)
    hist <- matrix(0, nrow = days, ncol = n_stages + 2)
    colnames(hist) <- c("Day", stages, "Total")

    for (t in seq_len(days)) {
      next_pop <- numeric(n_stages)
      for (i in seq_len(4)) {
        survivors <- round(pop[i] * (1 - mortality_p))
        if (survivors > 0) {
          trans <- rmultinom(1, size = survivors, prob = P[i, ])
          next_pop <- next_pop + trans[, 1]
        }
      }
      next_pop[1] <- next_pop[1] + round(pop[4] * repro_rate)
      pop <- next_pop
      hist[t, ] <- c(t, pop, sum(pop))
    }
    as.data.frame(hist)
  }

  # ---- Run Simulation ----
  observeEvent(input$run_sim, {
    updateActionButton(session, "run_sim", label = "Running...")
    P <- P_values()
    mortality_p <- input$mortality_prob
    repro_rate <- input$repro_rate
    days <- input$sim_days

    df <- run_single_sim(P, mortality_p, repro_rate, days)
    updateActionButton(session, "run_sim", label = "Run Simulation")

    # ---- Stage plot ----
    output$stage_plot <- renderPlot({
      df_long <- df %>%
        select(-Dead, -Total) %>%
        pivot_longer(cols = c(Egg, Larva, Pupa, Adult),
                     names_to = "Stage", values_to = "Count")
      stage_colors_named <- c(
        "Egg"   = "#56B4E9",
        "Larva" = "#E69F00",
        "Pupa"  = "#D55E00",
        "Adult" = "#7C0000"
      )
      if (input$view_mode == "prop") {
        df_long <- df_long %>%
          group_by(Day) %>%
          mutate(Proportion = Count / sum(Count))
        ggplot(df_long, aes(Day, Proportion, color = Stage)) +
          geom_line(linewidth = 1.5) +
          scale_color_manual(values = stage_colors_named) +
          theme_minimal(base_size = 14) +
          labs(title = "Active Mosquito Life Stages Over Time",
               y = "Proportion of Active Population", x = "Day")
      } else {
        ggplot(df_long, aes(Day, Count, color = Stage)) +
          geom_line(linewidth = 1.5) +
          scale_color_manual(values = stage_colors_named) +
          theme_minimal(base_size = 14) +
          labs(title = "Active Mosquito Life Stages Over Time",
               y = "Number of Individuals", x = "Day")
      }
    })

    # ---- Population plot ----
    output$pop_plot <- renderPlot({
      ggplot(df, aes(Day, Total)) +
        geom_line(linewidth = 1.4, color = "darkgreen") +
        theme_minimal(base_size = 14) +
        labs(title = "Total Population Size Over Time",
             y = "Number of Individuals", x = "Day")
    })
  })
}

shinyApp(ui, server)
