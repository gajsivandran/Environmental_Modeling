# app_harvest_P0_delay.R
# Logistic growth with constant harvesting + delayed start (Adjust H, P₀, and start delay)
# All displayed values use 3 significant figures

library(shiny)
library(tidyverse)

# ---------------------------
# Data load (with safe fallback)
# ---------------------------
load_cod_data <- function() {
  path <- "data/cod_timeseries.csv"
  read_csv(path, show_col_types = FALSE) |>
    select(Year, Pop) |>
    arrange(Year)

}

cod <- load_cod_data() |>
  mutate(t = Year - min(Year))

# ---------------------------
# Helpers
# ---------------------------
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
logistic_fun <- function(t, K, A, r) K / (1 + A * exp(-r * t))
dP_dt <- function(P, r, K, H) r * P * (1 - P / K) - H

# Delayed-harvest simulator: harvesting turns on when t >= H_start
simulate_harvest <- function(P0, r, K, H, H_start = 0,
                             years = 40, dt = 0.1, year0 = 0) {
  n <- max(2, ceiling(years / dt) + 1)
  t <- seq(0, years, length.out = n)
  P <- numeric(n)
  P[1] <- max(0, P0)
  for (i in 1:(n - 1)) {
    H_eff <- ifelse(t[i] >= H_start, H, 0)
    P[i + 1] <- max(0, P[i] + dt * dP_dt(P[i], r, K, H_eff))
  }
  tibble(time = t, Year = year0 + t, Pop = P)
}

# ---------------------------
# Fit logistic with H = 0
# ---------------------------
pop_max <- max(cod$Pop, na.rm = TRUE)
pop_min <- min(cod$Pop, na.rm = TRUE)
K_start <- 1.2 * pop_max
r_start <- 0.1
P0_obs  <- cod$Pop[cod$t == 0][1]
if (is.na(P0_obs)) P0_obs <- pop_min
A_start <- (K_start / (P0_obs + 1e-6)) - 1
if (!is.finite(A_start) || A_start <= 0) A_start <- 1

fit_logis <- tryCatch(
  nls(
    Pop ~ K / (1 + A * exp(-r * t)),
    data = cod,
    start = list(K = K_start, A = A_start, r = r_start),
    control = nls.control(maxiter = 2000, warnOnly = TRUE)
  ),
  error = function(e) NULL
)

if (!is.null(fit_logis)) {
  pars <- as.list(coef(fit_logis))
  K_fix <- pars$K
  A_fix <- pars$A
  r_fix <- pars$r
  P0_fix <- K_fix / (1 + A_fix)
  rmse_best <- rmse(cod$Pop, logistic_fun(cod$t, K_fix, A_fix, r_fix))
} else {
  K_fix <- K_start; A_fix <- A_start; r_fix <- r_start
  P0_fix <- K_fix / (1 + A_fix)
  rmse_best <- NA_real_
}

H_MSY_fix <- r_fix * K_fix / 4

# ---------------------------
# UI
# ---------------------------
ui <- fluidPage(
  titlePanel("Logistic Growth with Constant Harvesting (Adjust H, P₀, and Harvest Start)"),
  sidebarLayout(
    sidebarPanel(
      helpText(HTML(
        "Baseline (H=0) fitted logistic model:<br/>",
        "<b>Pop(t) = K / (1 + A e<sup>-r t</sup>)</b><br/>",
        "Adjust harvest rate <b>H</b>, initial population <b>P<sub>0</sub></b>, ",
        "and when harvesting begins."
      )),
      sliderInput(
        "H", "Harvest rate H (units/year):",
        min = 0,
        max = signif(1.5 * H_MSY_fix, 3),
        value = 0,
        step = signif((1.5 * H_MSY_fix) / 200, 3)
      ),
      sliderInput(
        "H_start", "Harvest begins after (years):",
        min = 0, max = 100, value = 0, step = 1
      ),
      sliderInput(
        "P0", "Initial population P₀:",
        min = 0, max = signif(1.2 * K_fix, 3),
        value = signif(P0_fix, 3),
        step = signif((1.2 * K_fix) / 200, 3)
      ),
      sliderInput(
        "years", "Simulation horizon (years):",
        min = 5, max = 100, value = 40, step = 1
      ),
      hr(),
      strong("Fixed model parameters (from H=0 fit):"),
      tags$div(HTML(
        paste0(
          "K = ", signif(K_fix, 3), "<br/>",
          "r = ", signif(r_fix, 3), "<br/>",
          "A = ", signif(A_fix, 3),
          "  (implies baseline P<sub>0</sub> = ", signif(P0_fix, 3), ")<br/>",
          "H<sub>MSY</sub> = ", signif(H_MSY_fix, 3)
        )
      )),
      hr(),
      strong("Outcome:"),
      textOutput("outcome_text")
    ),
    mainPanel(
      plotOutput("simplot", height = 520),
      br(),
      tags$small(em(
        if (file.exists("data/cod_timeseries.csv")) {
          "Loaded data/cod_timeseries.csv (points)."
        } else {
          "No data/cod_timeseries.csv found — using synthetic demo data (points)."
        }
      ))
    )
  )
)

# ---------------------------
# Server
# ---------------------------
server <- function(input, output, session) {

  msy <- reactive({ r_fix * K_fix / 4 })

  outcome_label <- reactive({
    H <- input$H
    if (H <= msy()) {
      Pstar <- (K_fix / 2) * (1 + sqrt(1 - 4 * H / (r_fix * K_fix)))
      paste0("Sustainable (once harvesting starts): equilibrium ~ ", signif(Pstar, 3))
    } else {
      "Unsustainable: no positive equilibrium (collapse likely once harvesting starts)"
    }
  })
  output$outcome_text <- renderText(outcome_label())

  sim <- reactive({
    simulate_harvest(
      P0 = input$P0, r = r_fix, K = K_fix,
      H = input$H, H_start = input$H_start,
      years = input$years, dt = 0.1,
      year0 = min(cod$Year)
    )
  })

  output$simplot <- renderPlot({
    df_sim <- sim()
    year_min <- min(cod$Year)
    year_max <- year_min + input$years

    grid_fit <- tibble(
      Year = seq(min(cod$Year), max(cod$Year), length.out = 400)
    ) |>
      mutate(t = Year - min(cod$Year),
             Pop_fit = logistic_fun(t, K_fix, A_fix, r_fix))

    # Sustainable equilibrium value (if any) for the chosen H
    eq_y <- NA_real_
    if (input$H <= msy()) {
      eq_y <- (K_fix / 2) * (1 + sqrt(1 - 4 * input$H / (r_fix * K_fix)))
    }

    # Build plot
    p <- ggplot() +
      geom_point(data = cod, aes(Year, Pop), size = 2, alpha = 0.85) +
      geom_line(data = cod, aes(Year, Pop), alpha = 0.25) +
      geom_line(data = grid_fit, aes(Year, Pop_fit),
                linetype = "dashed", color = "grey40", linewidth = 1) +
      geom_line(data = df_sim, aes(Year, Pop),
                linewidth = 1.2, color = "#0072B2") +
      # Harvest start marker
      geom_vline(xintercept = year_min + input$H_start,
                 linetype = "dashed", color = "red", alpha = 0.6) +
      labs(
        title = "Logistic Growth with Constant Harvesting (Delayed Start)",
        subtitle = paste0(
          "r = ", signif(r_fix, 3),
          ", K = ", signif(K_fix, 3),
          " | H = ", signif(input$H, 3),
          ", P₀ = ", signif(input$P0, 3),
          " | H starts at t = ", signif(input$H_start, 3), " yr",
          " | H_MSY = ", signif(msy(), 3)
        ),
        x = "Year", y = "Population (units)"
      ) +
      theme_classic()

    # Add equilibrium segment only after harvesting begins and only if sustainable
    if (is.finite(eq_y)) {
      x0 <- year_min + input$H_start
      x1 <- year_max
      if (x1 > x0) {
        p <- p + geom_segment(x = x0, xend = x1, y = eq_y, yend = eq_y,
                              linetype = "dotdash", color = "grey30")
      }
    }

    p
  })
}

shinyApp(ui, server)
