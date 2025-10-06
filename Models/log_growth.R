# app_logistic_base.R
# Interactive logistic fit (base R nls only) for cod population data
# Sliders and outputs use 3 significant figures

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
logistic_fun <- function(t, K, A, r) K / (1 + A * exp(-r * t))

# Heuristic starting values for nls and slider defaults
pop_max <- max(cod$Pop, na.rm = TRUE)
pop_min <- min(cod$Pop, na.rm = TRUE)
K_start <- 1.2 * pop_max
r_start <- 0.1
P0 <- cod$Pop[cod$t == 0][1]
if (is.na(P0)) P0 <- pop_min
A_start <- (K_start / (P0 + 1e-6)) - 1
if (!is.finite(A_start) || A_start <= 0) A_start <- 1

# Fit best logistic model (nls)
fit_best <- tryCatch(
  nls(
    Pop ~ K / (1 + A * exp(-r * t)),
    data = cod,
    start = list(K = K_start, A = A_start, r = r_start),
    control = nls.control(maxiter = 2000, warnOnly = TRUE)
  ),
  error = function(e) NULL
)

has_best <- !is.null(fit_best)
best_par <- if (has_best) as.list(coef(fit_best)) else list(K = K_start, A = A_start, r = r_start)
best_rmse <- if (has_best) rmse(cod$Pop, fitted(fit_best)) else NA_real_

# Slider ranges (use rough scale spacing)
K_rng <- c(0.5 * pop_max, 3 * pop_max)
A_rng <- c(1e-3, 20)
r_rng <- c(-1, 1)

# ---------------------------
# UI
# ---------------------------
ui <- fluidPage(
  titlePanel("Interactive Logistic Fit (Base nls): Cod Population"),
  sidebarLayout(
    sidebarPanel(
      helpText(HTML("Model: <b>Pop(t) = K / (1 + A e<sup>-r t</sup>)</b>, with t = Year - min(Year)")),

      sliderInput("K", "K (carrying capacity):",
                  min = signif(K_rng[1], 3), max = signif(K_rng[2], 3),
                  value = signif(best_par$K, 3), step = signif(diff(K_rng)/200, 3)),

      sliderInput("A", "A (initial position):",
                  min = signif(A_rng[1], 3), max = signif(A_rng[2], 3),
                  value = signif(best_par$A, 3), step = 0.01),

      sliderInput("r", "r (growth rate):",
                  min = signif(r_rng[1], 3), max = signif(r_rng[2], 3),
                  value = signif(best_par$r, 3), step = 0.001),

      actionButton("snap", "Snap to Best Fit (nls)"),
      hr(),
      strong("Your logistic curve:"),
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
      plotOutput("fitplot", height = 500),
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

  observeEvent(input$snap, {
    updateSliderInput(session, "K", value = signif(best_par$K, 3))
    updateSliderInput(session, "A", value = signif(best_par$A, 3))
    updateSliderInput(session, "r", value = signif(best_par$r, 3))
  })

  preds_user <- reactive({
    tibble(
      Year = cod$Year,
      t    = cod$t,
      Pred = logistic_fun(cod$t, input$K, input$A, input$r)
    )
  })

  output$eqn_user <- renderText({
    paste0("Pop(t) = ", signif(input$K, 3),
           " / (1 + ", signif(input$A, 3),
           " * exp(-", signif(input$r, 3), " * t))")
  })

  output$rmse_user <- renderText({
    sprintf("%.3f", rmse(cod$Pop, preds_user()$Pred))
  })

  output$eqn_best <- renderText({
    if (has_best) {
      paste0("Pop(t) = ", signif(best_par$K, 3),
             " / (1 + ", signif(best_par$A, 3),
             " * exp(-", signif(best_par$r, 3), " * t))")
    } else {
      "Best fit unavailable (nls did not converge). Try adjusting sliders."
    }
  })

  output$rmse_best <- renderText({
    if (has_best) sprintf("%.3f", best_rmse) else "—"
  })

  output$fitplot <- renderPlot({
    grid <- tibble(Year = seq(min(cod$Year), max(cod$Year), length.out = 400)) |>
      mutate(t = Year - min(cod$Year),
             Pred_user = logistic_fun(t, input$K, input$A, input$r),
             Pred_best = if (has_best)
               logistic_fun(t, best_par$K, best_par$A, best_par$r) else NA_real_)

    ggplot(cod, aes(Year, Pop)) +
      geom_point(size = 2, alpha = 0.9) +
      geom_line(alpha = 0.35) +
      geom_line(data = grid, aes(y = Pred_user), linewidth = 1, color = "#0072B2") +
      { if (has_best) geom_line(data = grid, aes(y = Pred_best),
                                linewidth = 1, linetype = "dashed", color = "grey40") } +
      labs(
        title = "Cod Population with Interactive Logistic Fit (Base nls)",
        subtitle = paste0(
          "K = ", signif(input$K, 3),
          ", A = ", signif(input$A, 3),
          ", r = ", signif(input$r, 3)
        ),
        x = "Year", y = "Population (units)"
      ) +
      theme_classic()
  })
}

shinyApp(ui, server)
