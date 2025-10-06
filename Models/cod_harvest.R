# app_harvest_P0.R
# Logistic growth with constant harvesting (Adjust H and P₀)
# All displayed values use 3 significant figures

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
dP_dt <- function(P, r, K, H) r * P * (1 - P / K) - H

simulate_harvest <- function(P0, r, K, H, years = 40, dt = 0.1, year0 = 0) {
  n <- max(2, ceiling(years / dt) + 1)
  t <- seq(0, years, length.out = n)
  P <- numeric(n)
  P[1] <- max(0, P0)
  for (i in 1:(n - 1)) {
    P[i + 1] <- max(0, P[i] + dt * dP_dt(P[i], r, K, H))
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
  titlePanel("Logistic Growth with Constant Harvesting (Adjust H and P₀)"),
  sidebarLayout(
    sidebarPanel(
      helpText(HTML(
        "Baseline (H=0) fitted logistic model:<br/>",
        "<b>Pop(t) = K / (1 + A e<sup>-r t</sup>)</b><br/>",
        "You can vary harvest rate <b>H</b> and initial population <b>P<sub>0</sub></b>."
      )),
      sliderInput(
        "H", "Harvest rate H (units/year):",
        min = 0,
        max = signif(1.5 * H_MSY_fix, 3),
        value = 0,
        step = signif((1.5 * H_MSY_fix) / 200, 3)
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
      paste0("Sustainable: equilibrium at ~ ", signif(Pstar, 3))
    } else {
      "Unsustainable: no positive equilibrium (collapse likely)"
    }
  })
  output$outcome_text <- renderText(outcome_label())

  sim <- reactive({
    simulate_harvest(
      P0 = input$P0, r = r_fix, K = K_fix,
      H = input$H, years = input$years,
      dt = 0.1, year0 = min(cod$Year)
    )
  })

  output$simplot <- renderPlot({
    df_sim <- sim()
    grid_fit <- tibble(
      Year = seq(min(cod$Year), max(cod$Year), length.out = 400)
    ) |>
      mutate(t = Year - min(cod$Year),
             Pop_fit = logistic_fun(t, K_fix, A_fix, r_fix))

    eq_y <- NA_real_
    if (input$H <= msy()) {
      eq_y <- (K_fix / 2) * (1 + sqrt(1 - 4 * input$H / (r_fix * K_fix)))
    }

    ggplot() +
      geom_point(data = cod, aes(Year, Pop), size = 2, alpha = 0.85) +
      geom_line(data = cod, aes(Year, Pop), alpha = 0.25) +
      geom_line(data = grid_fit, aes(Year, Pop_fit),
                linetype = "dashed", color = "grey40", linewidth = 1) +
      geom_line(data = df_sim, aes(Year, Pop),
                linewidth = 1.2, color = "#0072B2") +
      { if (is.finite(eq_y)) geom_hline(yintercept = eq_y,
                                        linetype = "dotdash", color = "grey30") } +
      labs(
        title = "Logistic Growth with Constant Harvesting",
        subtitle = paste0(
          "r = ", signif(r_fix, 3),
          ", K = ", signif(K_fix, 3),
          " | H = ", signif(input$H, 3),
          ", P₀ = ", signif(input$P0, 3),
          " | H_MSY = ", signif(msy(), 3)
        ),
        x = "Year", y = "Population (units)"
      ) +
      theme_classic()
  })
}

shinyApp(ui, server)
