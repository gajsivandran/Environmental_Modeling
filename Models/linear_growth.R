# app.R
# Interactive linear fit tuner (with best-RMSE reference) for cod population data
# All displayed values and sliders use 3 significant figures

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
    # Fallback synthetic series (so the app still runs in class)
    set.seed(42)
    tibble(
      Year = 1980:2009,
      Pop  = 200 / (1 + exp(-0.18*(Year - 1990))) *
        exp(-0.03*pmax(Year-1995,0)) + rnorm(30, 0, 3)
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

# OLS (min-RMSE for linear model)
ols <- lm(Pop ~ t, data = cod)
a_best <- unname(coef(ols)[1])
b_best <- unname(coef(ols)[2])
best_rmse <- rmse(cod$Pop, fitted(ols))

# Slider ranges (reasonable defaults)
t_range <- range(cod$t, na.rm = TRUE)
pop_rng <- range(cod$Pop, na.rm = TRUE)
a_rng   <- c(pop_rng[1] - 0.5 * diff(pop_rng),
             pop_rng[2] + 0.5 * diff(pop_rng))
b_scale <- ifelse(diff(t_range) > 0, diff(pop_rng) / diff(t_range), 1)
b_rng   <- c(-3, 3) * b_scale

# ---------------------------
# UI
# ---------------------------
ui <- fluidPage(
  titlePanel("Interactive Linear Fit (3 sig figs): Cod Population"),
  sidebarLayout(
    sidebarPanel(
      helpText("Model:  Pop(t) = a + b * t   with  t = Year - min(Year)"),

      sliderInput("a", "Intercept (a):",
                  min = signif(a_rng[1], 3), max = signif(a_rng[2], 3),
                  value = signif(a_best, 3),
                  step = signif(diff(a_rng) / 200, 3)),

      sliderInput("b", "Slope (b):",
                  min = signif(b_rng[1], 3), max = signif(b_rng[2], 3),
                  value = signif(b_best, 3),
                  step = signif(diff(b_rng) / 200, 3)),

      actionButton("snap", "Snap to OLS (min RMSE)"),
      hr(),
      strong("Your line:"),
      verbatimTextOutput("eqn_user", placeholder = TRUE),
      div("RMSE (your line):"),
      h3(textOutput("rmse_user"), style = "margin-top:-8px;"),
      hr(),
      strong("OLS reference (dashed line):"),
      verbatimTextOutput("eqn_best", placeholder = TRUE),
      div("RMSE (OLS):"),
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

  # Snap sliders to OLS on click
  observeEvent(input$snap, {
    updateSliderInput(session, "a", value = signif(a_best, 3))
    updateSliderInput(session, "b", value = signif(b_best, 3))
  })

  # Predictions for user's sliders
  preds_user <- reactive({
    tibble(
      Year = cod$Year,
      t    = cod$t,
      Pred = input$a + input$b * cod$t
    )
  })

  # Text outputs
  output$eqn_user <- renderText({
    paste0("Pop(t) = ", signif(input$a, 3),
           " + ", signif(input$b, 3), " * t")
  })

  output$rmse_user <- renderText({
    sprintf("%.3f", rmse(cod$Pop, preds_user()$Pred))
  })

  output$eqn_best <- renderText({
    paste0("Pop(t) = ", signif(a_best, 3),
           " + ", signif(b_best, 3), " * t")
  })

  output$rmse_best <- renderText({
    sprintf("%.3f", best_rmse)
  })

  # Plot
  output$fitplot <- renderPlot({
    # Smooth lines for display
    grid <- tibble(
      Year = seq(min(cod$Year), max(cod$Year), length.out = 400)
    ) |>
      mutate(t = Year - min(cod$Year),
             Pred_user = input$a + input$b * t,
             Pred_best = a_best + b_best * t)

    ggplot(cod, aes(Year, Pop)) +
      geom_point(size = 2, alpha = 0.9) +
      geom_line(alpha = 0.35) +
      # User-selected line
      geom_line(data = grid, aes(y = Pred_user),
                linewidth = 1, color = "#0072B2") +
      # OLS (best fit)
      geom_line(data = grid, aes(y = Pred_best),
                linewidth = 1, linetype = "dashed", color = "grey40") +
      labs(
        title = "Cod Population with Interactive Linear Fit",
        subtitle = paste0(
          "a = ", signif(input$a, 3),
          ", b = ", signif(input$b, 3),
          " | RMSE (yours) = ", sprintf("%.3f", rmse(cod$Pop, preds_user()$Pred)),
          " | RMSE (OLS) = ", sprintf("%.3f", best_rmse)
        ),
        x = "Year", y = "Population (units)"
      ) +
      theme_classic()
  })
}

shinyApp(ui, server)
