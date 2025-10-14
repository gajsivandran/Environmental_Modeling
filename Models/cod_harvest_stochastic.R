# app_logistic_harvest_badyears.R
# Logistic growth with harvesting and stochastic "bad-year" shocks to r
# - P0 is fixed from the H = 0 logistic fit (no slider)
# - Second panel shows r_eff(t) over time
# - All displayed values use 3 significant figures

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
  mutate(t = Year - min(Year))  # centered time for fitting

# ---------------------------
# Helpers
# ---------------------------
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
logistic_fun <- function(t, K, A, r) K / (1 + A * exp(-r * t))

# Euler integration of dP/dt = r_eff(t)*P*(1 - P/K) - H
# Returns full series including r_eff(t) and bad_year flag
simulate_badyears <- function(P0, r, K, H, years, dt, year0,
                              by_p, bad_factor, seed = NULL,
                              clip_nonneg = TRUE) {
  if (!is.null(seed)) set.seed(seed)

  n <- max(2, ceiling(years / dt) + 1)
  t <- seq(0, years, length.out = n)
  P <- numeric(n)
  P[1] <- max(0, P0)

  # Year indices (0,1,2,...) from start; draw bad years once per whole year
  yrs_idx <- floor(t)
  unique_years <- unique(yrs_idx)
  bad_year_flags <- rbinom(length(unique_years), size = 1, prob = by_p) == 1
  names(bad_year_flags) <- as.character(unique_years)

  r_eff_vec <- numeric(n)

  for (i in 1:(n - 1)) {
    yk <- as.character(yrs_idx[i])
    r_eff <- if (!is.null(bad_year_flags[yk]) && bad_year_flags[yk]) r * bad_factor else r
    r_eff_vec[i] <- r_eff
    dP <- r_eff * P[i] * (1 - P[i] / K) - H
    P[i + 1] <- P[i] + dt * dP
    if (clip_nonneg) P[i + 1] <- max(0, P[i + 1])
  }
  # last step r_eff
  yk_last <- as.character(yrs_idx[n])
  r_eff_vec[n] <- if (!is.null(bad_year_flags[yk_last]) && bad_year_flags[yk_last]) r * bad_factor else r

  tibble(
    time     = t,
    Year     = year0 + t,
    Pop      = P,
    r_eff    = r_eff_vec,
    bad_year = bad_year_flags[as.character(floor(t))]
  )
}

# ---------------------------
# Fit baseline (H = 0) logistic for defaults (K, r, A, implied P0)
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
  P0_fix <- K_fix / (1 + A_fix) # implied initial pop at t=0 (H=0 fit)
  rmse_best <- rmse(cod$Pop, logistic_fun(cod$t, K_fix, A_fix, r_fix))
} else {
  K_fix <- K_start; A_fix <- A_start; r_fix <- r_start
  P0_fix <- K_fix / (1 + A_fix)
  rmse_best <- NA_real_
}

H_MSY_baseline <- r_fix * K_fix / 4  # for slider scaling

# ---------------------------
# UI
# ---------------------------
ui <- fluidPage(
  titlePanel("Logistic with Harvesting and Random Bad-Year Shocks to r (P₀ fixed)"),
  sidebarLayout(
    sidebarPanel(
      helpText(HTML(
        "Model: <b>dP/dt = r_eff(t) P (1 - P/K) - H</b><br/>",
        "Baseline <i>r</i>, <i>K</i>, and <i>P<sub>0</sub></i> from the H=0 logistic fit.<br/>",
        "Each year is 'bad' with probability <i>p</i>; then <i>r</i> is scaled by <i>bad_factor</i>."
      )),
      sliderInput("H", "Harvest rate H (units/year):",
                  min = 0,
                  max = signif(1.5 * H_MSY_baseline, 3),
                  value = signif(0.5 * H_MSY_baseline, 3),
                  step = signif((1.5 * H_MSY_baseline) / 200, 3)),
      sliderInput("years", "Simulation horizon (years):",
                  min = 5, max = 100, value = 40, step = 1),
      hr(),
      strong("Bad-year process on r:"),
      sliderInput("by_p", "Bad-year probability p (per year):",
                  min = 0, max = 0.7, value = 0.25, step = 0.01),
      sliderInput("bad_factor", "Bad-year multiplier on r (e.g., 0.5 halves r):",
                  min = 0.05, max = 1.0, value = 0.5, step = 0.01),
      numericInput("seed", "Random seed", value = 123, min = 1, step = 1),
      actionButton("resim", "Resimulate"),
      hr(),
      strong("Fixed parameters (from H = 0 fit):"),
      tags$div(HTML(
        paste0(
          "K = ", signif(K_fix, 3), "<br/>",
          "r = ", signif(r_fix, 3), "<br/>",
          "A = ", signif(A_fix, 3),
          "  ⇒ P<sub>0</sub> (fixed) = ", signif(P0_fix, 3), "<br/>",
          "RMSE(H=0 fit) = ", ifelse(is.na(rmse_best), "—", sprintf('%.3f', rmse_best))
        )
      ))
    ),
    mainPanel(
      plotOutput("popplot", height = 420),
      plotOutput("rplot",  height = 200),
      br(),
      tags$small(em(
        if (file.exists("data/cod_timeseries.csv")) {
          "Loaded data/cod_timeseries.csv (points)."
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

  # Resimulation trigger to redraw with same settings
  bump <- reactiveVal(0)
  observeEvent(input$resim, { bump(isolate(bump()) + 1) })

  # Simulate with bad-year shocks (P0 fixed)
  sim_df <- reactive({
    bump()
    simulate_badyears(
      P0 = P0_fix,            # <-- P0 fixed here
      r  = r_fix,
      K  = K_fix,
      H  = input$H,
      years = input$years,
      dt = 0.1,
      year0 = min(cod$Year),
      by_p = input$by_p,
      bad_factor = input$bad_factor,
      seed = input$seed,
      clip_nonneg = TRUE
    )
  })

  output$popplot <- renderPlot({
    df_sim <- sim_df()

    # Baseline (H=0) best-fit over observed years (for reference)
    grid_fit <- tibble(Year = seq(min(cod$Year), max(cod$Year), length.out = 400)) |>
      mutate(t = Year - min(cod$Year),
             Pop_fit = logistic_fun(t, K_fix, A_fix, r_fix))

    # MSY numbers for subtitle
    Hmsy_base <- r_fix * K_fix / 4
    Hmsy_bad  <- input$bad_factor * r_fix * K_fix / 4

    # Rectangles marking bad years in the simulation
    bad_years <- df_sim |>
      mutate(y0 = floor(time)) |>
      group_by(y0) |>
      summarise(is_bad = any(bad_year), .groups = "drop") |>
      filter(is_bad) |>
      mutate(xmin = min(cod$Year) + y0,
             xmax = xmin + 1)

    ggplot() +
      { if (nrow(bad_years) > 0)
        geom_rect(data = bad_years,
                  aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                  fill = "red", alpha = 0.08) } +
      geom_point(data = cod, aes(Year, Pop), size = 2, alpha = 0.85) +
      geom_line(data = cod, aes(Year, Pop), alpha = 0.25) +
      geom_line(data = grid_fit, aes(Year, Pop_fit),
                linetype = "dashed", color = "grey40", linewidth = 1) +
      geom_line(data = df_sim, aes(Year, Pop),
                linewidth = 1.2, color = "#0072B2") +
      labs(
        title = "Logistic Growth with Harvesting and Bad-Year Shocks to r (P₀ fixed)",
        subtitle = paste0(
          "H = ", signif(input$H, 3),
          " | MSY baseline = ", signif(Hmsy_base, 3),
          ", MSY in bad year = ", signif(Hmsy_bad, 3),
          " | p(bad year) = ", signif(input$by_p, 3),
          ", factor = ", signif(input$bad_factor, 3)
        ),
        x = "Year", y = "Population (units)"
      ) +
      theme_classic()
  })

  output$rplot <- renderPlot({
    df_sim <- sim_df()

    # Build yearly shading again for alignment
    bad_years <- df_sim |>
      mutate(y0 = floor(time)) |>
      group_by(y0) |>
      summarise(is_bad = any(bad_year), .groups = "drop") |>
      filter(is_bad) |>
      mutate(xmin = min(cod$Year) + y0,
             xmax = xmin + 1)

    ggplot() +
      { if (nrow(bad_years) > 0)
        geom_rect(data = bad_years,
                  aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                  fill = "red", alpha = 0.08) } +
      geom_step(data = df_sim, aes(Year, r_eff), linewidth = 1) +
      geom_hline(yintercept = r_fix, linetype = "dashed", color = "grey40") +
      geom_hline(yintercept = r_fix * input$bad_factor,
                 linetype = "dotted", color = "grey50") +
      labs(
        title = "Effective Growth Rate r_eff(t) (bad years scale r)",
        subtitle = paste0(
          "Baseline r = ", signif(r_fix, 3),
          " | Bad-year r = ", signif(r_fix * input$bad_factor, 3)
        ),
        x = "Year", y = "r_eff(t)"
      ) +
      theme_classic()
  })
}

shinyApp(ui, server)
