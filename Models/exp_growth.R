# app_exp.R
# Interactive exponential fit tuner (with best-RMSE reference) for cod population data
# All values and sliders use 3 significant figures

library(shiny)
library(tidyverse)

# ---------------------------
# Data load (with safe fallback)
# ---------------------------
load_cod_data <- function() {
  path <- "data/cod_timeseries.csv"
  if (file.exists(path)) {
    read_csv(path, show_col_types = FALSE) |>
      select(Year, Pop) |>
      arrange(Year)
  } else {
    # Fallback synthetic series
    set.seed(42)
    tibble(
      Year = 1980:2009,
      Pop  = 200 / (1 + exp(-0.18 * (Year - 1990))) *
        exp(-0.03 * pmax(Year - 1995, 0)) + rnorm(30, 0, 3)
    ) |>
      mutate(Pop = pmax(Pop, 0.1)) |>
      arrange(Year)
  }
}

cod <- load_cod_data() |>
  mutate(t = Year - min(Year))

# ---------------------------
# Helpers
# ---------------------------
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
exp_fun <- function(t, a, b) a * exp(b * t)

# Try to get the best-fit exponential model
start_vals <- list(a = min(cod$Pop, na.rm = TRUE), b = 0.05)
fit_best <- tryCatch(
  nls(Pop ~ a * exp(b * t), data = cod, start = start_vals,
      control = nls.control(maxiter = 2000, warnOnly = TRUE)),
  error = function(e) NULL
)

has_best <- !is.null(fit_best)
best_par <- if (has_best) as.list(coef(fit_best)) else start_vals
best_rmse <- if (has_best) rmse(cod$Pop, fitted(fit_best)) else NA_real_

# Define slider ranges
a_rng <- c(0.1, 2 * max(cod$Pop, na.rm = TRUE))
b_rng <- c(-0.2, 0.5)

# ---------------------------
# UI
# ---------------------------
ui <- fluidPage(
  titlePanel("Interactive Exponential Fit (Base nls): Cod Population"),
  sidebarLayout(
    sidebarPanel(
      helpText(HTML("Model: <b>Pop(t) = a e<sup>b t</sup></b>, with t = Year - min(Year)")),

      sliderInput("a", "a (initial value):",
                  min = signif(a_rng[1], 3), max = signif(a_rng[2], 3),
                  value = signif(best_par$a, 3),
                  step = signif(diff(a_rng) / 200, 3)),

      sliderInput("b", "b (growth rate):",
                  min = signif(b_rng[1], 3), max = signif(b_rng[2], 3),
                  value = signif(best_par$b, 3),
                  step = signif(diff(b_rng) / 200, 3)),

      actionButton("snap", "Snap to Best Fit (nls)"),
      hr(),
      strong("Your exponential curve:"),
      verbatimTextOutput("eqn_user", placeholder = TRUE),
      div("RMSE (your curve):"),
      h3(textOutput("rmse_user"), style = "margin-top:-8px;"),
      hr(),
      strong("Best fit (dashed line):"),
      verbatimTextOutput("eqn_best", placeholder = TRUE),
      div("RMSE (best fit):"),
      h4(textOutput("rmse_best"), style = "margin-top:-8px;")
    ),

    mainPanel(
      plotOutput("fitplot", height = 480),
      br(),
      tags$small(em(
        if (file.exists("data/cod_timeseries.csv")) {
          "Loaded data/cod_timeseries.csv"
        } else {
          "No data/cod_timeseries.csv found — using synthetic demo data."
        }
      ))
    )
  )
)

# ---------------------------
# Server
# ---------------------------
server <- function(input, output, session) {

  # Snap sliders to best-fit
  observeEvent(input$snap, {
    updateSliderInput(session, "a", value = signif(best_par$a, 3))
    updateSliderInput(session, "b", value = signif(best_par$b, 3))
  })

  preds_user <- reactive({
    tibble(
      Year = cod$Year,
      t    = cod$t,
      Pred = exp_fun(cod$t, input$a, input$b)
    )
  })

  # Text outputs
  output$eqn_user <- renderText({
    paste0("Pop(t) = ", signif(input$a, 3), " × exp(", signif(input$b, 3), " × t)")
  })

  output$rmse_user <- renderText({
    sprintf("%.3f", rmse(cod$Pop, preds_user()$Pred))
  })

  output$eqn_best <- renderText({
    if (has_best) {
      paste0("Pop(t) = ", signif(best_par$a, 3),
             " × exp(", signif(best_par$b, 3), " × t)")
    } else {
      "Best fit unavailable (nls did not converge)."
    }
  })

  output$rmse_best <- renderText({
    if (has_best) sprintf("%.3f", best_rmse) else "—"
  })

  # Plot
  output$fitplot <- renderPlot({
    grid <- tibble(Year = seq(min(cod$Year), max(cod$Year), length.out = 400)) |>
      mutate(t = Year - min(cod$Year),
             Pred_user = exp_fun(t, input$a, input$b),
             Pred_best = if (has_best)
               exp_fun(t, best_par$a, best_par$b) else NA_real_)

    ggplot(cod, aes(Year, Pop)) +
      geom_point(size = 2, alpha = 0.9) +
      geom_line(alpha = 0.35) +
      geom_line(data = grid, aes(y = Pred_user), linewidth = 1, color = "#0072B2") +
      { if (has_best) geom_line(data = grid, aes(y = Pred_best),
                                linewidth = 1, linetype = "dashed", color = "grey40") } +
      labs(
        title = "Cod Population with Interactive Exponential Fit",
        subtitle = paste0(
          "a = ", signif(input$a, 3),
          ", b = ", signif(input$b, 3),
          " | RMSE (yours) = ", sprintf("%.3f", rmse(cod$Pop, preds_user()$Pred)),
          if (has_best) paste0(" | RMSE (best) = ", sprintf("%.3f", best_rmse)) else ""
        ),
        x = "Year", y = "Population (units)"
      ) +
      theme_classic()
  })
}

shinyApp(ui, server)
