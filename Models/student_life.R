# app.R — Markov Chains: Hourly Energy Simulation (Simplified Compare Tab)

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)

# ----------------------------------------------------------
# Shared setup
# ----------------------------------------------------------
stages <- c("Tired", "Awake", "Focused", "Distracted", "Exhausted")

# Transition matrices
P_base <- matrix(c(
  0.2,0.7,0.0,0.0,0.1,
  0.0,0.3,0.5,0.1,0.1,
  0.0,0.1,0.6,0.2,0.1,
  0.1,0.2,0.2,0.3,0.2,
  0.5,0.0,0.0,0.0,0.5
), nrow = 5, byrow = TRUE, dimnames = list(stages, stages))

P_caff <- matrix(c(
  0.1,0.8,0.1,0.0,0.0,
  0.0,0.2,0.6,0.1,0.1,
  0.0,0.1,0.6,0.2,0.1,
  0.0,0.1,0.3,0.4,0.2,
  0.5,0.1,0.0,0.0,0.4
), nrow = 5, byrow = TRUE, dimnames = list(stages, stages))

P_steady <- matrix(c(
  0.3,0.6,0.0,0.0,0.1,
  0.0,0.5,0.3,0.1,0.1,
  0.0,0.1,0.7,0.1,0.1,
  0.1,0.3,0.2,0.3,0.1,
  0.3,0.2,0.0,0.0,0.5
), nrow = 5, byrow = TRUE, dimnames = list(stages, stages))

P_dreamer <- matrix(c(
  0.4,0.4,0.0,0.1,0.1,
  0.1,0.3,0.1,0.4,0.1,
  0.0,0.1,0.4,0.3,0.2,
  0.1,0.1,0.1,0.5,0.2,
  0.5,0.0,0.0,0.0,0.5
), nrow = 5, byrow = TRUE, dimnames = list(stages, stages))

# Sampling step
step_once <- function(current_state_index, P) {
  sample.int(5, size = 1, prob = P[current_state_index, ])
}

# ----------------------------------------------------------
# Module for single personality
# ----------------------------------------------------------
matrix_tab_ui <- function(id, title) {
  ns <- NS(id)
  tagList(
    h3(title),
    fluidRow(
      column(6,
             tags$div(style="margin-bottom:6px;",
                      strong("Current hour: "), textOutput(ns("hour"), inline = TRUE),
                      HTML("&nbsp;&nbsp;|&nbsp;&nbsp;"),
                      strong("Current state: "), textOutput(ns("state_label"), inline = TRUE)
             ),
             actionButton(ns("sim1"), "Simulate 1 hour", class = "btn-primary"),
             actionButton(ns("sim10"), "Simulate 10 hours", class = "btn-secondary", style="margin-left:8px;"),
             actionButton(ns("reset"), "Reset", class = "btn-outline-danger", style="margin-left:8px;"),
             br(), br(),
             plotOutput(ns("hist_plot"), height = "350px")
      ),
      column(6,
             h4("Transition Matrix (row = current, column = next)"),
             tableOutput(ns("matrix_tbl")),
             tags$small("Each row sums to 1. The next state is sampled from that row.")
      )
    )
  )
}

matrix_tab_server <- function(id, P, color_palette) {
  moduleServer(id, function(input, output, session) {
    rv <- reactiveValues(
      hour = 0,
      current = 1L,
      counts = setNames(rep(0L, 5), stages)
    )

    output$hour <- renderText(rv$hour)
    output$state_label <- renderText(stages[rv$current])

    output$matrix_tbl <- renderTable({
      as.data.frame(round(P, 2)) |>
        `colnames<-`(stages) |>
        tibble::rownames_to_column("From")
    }, bordered = TRUE, striped = TRUE, digits = 2)

    output$hist_plot <- renderPlot({
      df <- tibble(Stage = factor(names(rv$counts), levels = stages),
                   Hours = as.numeric(rv$counts))
      ggplot(df, aes(Stage, Hours, fill = Stage)) +
        geom_col(width = 0.7) +
        scale_fill_manual(values = color_palette, drop = FALSE) +
        labs(title = "Accumulated Time in Each Stage",
             x = NULL, y = "Hours spent") +
        theme_minimal(base_size = 14) +
        theme(legend.position = "none")
    })

    observeEvent(input$sim1, {
      next_state <- step_once(rv$current, P)
      rv$hour <- rv$hour + 1
      rv$current <- next_state
      rv$counts[next_state] <- rv$counts[next_state] + 1L
    })

    observeEvent(input$sim10, {
      for (k in 1:10) {
        next_state <- step_once(rv$current, P)
        rv$hour <- rv$hour + 1
        rv$current <- next_state
        rv$counts[next_state] <- rv$counts[next_state] + 1L
      }
    })

    observeEvent(input$reset, {
      rv$hour <- 0
      rv$current <- 1L
      rv$counts <- setNames(rep(0L, 5), stages)
    })
  })
}

# ----------------------------------------------------------
# Compare Panel (Manual Simulation)
# ----------------------------------------------------------
compare_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Compare Personalities"),
    p("Simulate all three personality types side by side. Each simulation step represents one hour."),
    actionButton(ns("sim1_all"), "Simulate 1 hour", class = "btn-primary"),
    actionButton(ns("sim10_all"), "Simulate 10 hours", class = "btn-secondary", style="margin-left:8px;"),
    actionButton(ns("reset_all"), "Reset All", class = "btn-outline-danger", style="margin-left:8px;"),
    br(), br(),
    plotOutput(ns("compare_plot"), height = "500px")
  )
}

compare_server <- function(id, matrices, colors) {
  moduleServer(id, function(input, output, session) {
    rv <- reactiveValues(
      hour = 0,
      counts = replicate(3, setNames(rep(0, 5), stages), simplify = FALSE),
      current = rep(1L, 3)
    )

    observeEvent(input$sim1_all, {
      for (i in 1:3) {
        next_state <- step_once(rv$current[i], matrices[[i]])
        rv$current[i] <- next_state
        rv$counts[[i]][next_state] <- rv$counts[[i]][next_state] + 1
      }
      rv$hour <- rv$hour + 1
    })

    observeEvent(input$sim10_all, {
      for (k in 1:10) {
        for (i in 1:3) {
          next_state <- step_once(rv$current[i], matrices[[i]])
          rv$current[i] <- next_state
          rv$counts[[i]][next_state] <- rv$counts[[i]][next_state] + 1
        }
        rv$hour <- rv$hour + 1
      }
    })

    observeEvent(input$reset_all, {
      rv$hour <- 0
      rv$current <- rep(1L, 3)
      rv$counts <- replicate(3, setNames(rep(0, 5), stages), simplify = FALSE)
    })

    output$compare_plot <- renderPlot({
      df <- bind_rows(
        tibble(Personality = "Over-Caffeinated", Stage = names(rv$counts[[1]]), Hours = rv$counts[[1]]),
        tibble(Personality = "Steady Strategist", Stage = names(rv$counts[[2]]), Hours = rv$counts[[2]]),
        tibble(Personality = "Procrastinating Dreamer", Stage = names(rv$counts[[3]]), Hours = rv$counts[[3]])
      ) %>%
        mutate(
          Stage = factor(Stage, levels = stages),
          Total = ave(Hours, Personality, FUN = sum),
          Proportion = ifelse(Total > 0, Hours / Total, 0)
        )

      ggplot(df, aes(Stage, Proportion, fill = Stage)) +
        geom_col(width = 0.7) +
        scale_fill_manual(values = colors) +
        facet_wrap(~Personality, nrow = 1) +
        labs(title = paste("Accumulated Proportion of Time per Stage — Hour", rv$hour),
             y = "Proportion of Hours", x = NULL) +
        theme_minimal(base_size = 14) +
        theme(legend.position = "bottom")
    })
  })
}


# ----------------------------------------------------------
# Colors + Main App
# ----------------------------------------------------------
stage_colors <- c(
  "Tired" = "#8da0cb",
  "Awake" = "#66c2a5",
  "Focused" = "#e78ac3",
  "Distracted" = "#fc8d62",
  "Exhausted" = "#a6d854"
)

ui <- fluidPage(
  titlePanel("Markov Chains — Hourly Energy Simulation"),
  p("Simulate how energy and focus fluctuate throughout the day for different personality types."),
  tabsetPanel(
    tabPanel("Base", matrix_tab_ui("base", "Base Transition Matrix")),
    tabPanel("Over-Caffeinated Optimist", matrix_tab_ui("caff", "Over-Caffeinated Optimist")),
    tabPanel("Steady Strategist", matrix_tab_ui("steady", "Steady Strategist")),
    tabPanel("Procrastinating Dreamer", matrix_tab_ui("dream", "Procrastinating Dreamer")),
    tabPanel("Compare Personalities", compare_ui("compare"))
  )
)

server <- function(input, output, session) {
  matrix_tab_server("base",   P = P_base,   color_palette = stage_colors)
  matrix_tab_server("caff",   P = P_caff,   color_palette = stage_colors)
  matrix_tab_server("steady", P = P_steady, color_palette = stage_colors)
  matrix_tab_server("dream",  P = P_dreamer, color_palette = stage_colors)
  compare_server("compare",
                 matrices = list(P_caff, P_steady, P_dreamer),
                 colors = stage_colors)
}

shinyApp(ui, server)
