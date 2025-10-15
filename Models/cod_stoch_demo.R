# app_logistic_harvest_badyears.R
# Logistic growth with harvesting and stochastic "bad-year" shocks to r
# - P0 is fixed from the H = 0 logistic fit (no slider)
# - r is randomly scaled in "bad years"
# - Tracks Pop at target year (2050): scatter of runs + histogram with fitted distributions
# - Buttons: Resimulate (1 run), Run 10 simulations, Run 100 simulations, Clear points
# - All displayed values use 3 significant figures

library(shiny)
library(tidyverse)
library(MASS)   # for fitdistr

# ---------------------------
# Data load (with safe fallback)
# ---------------------------
load_cod_data <- function() {
  path <- "data/cod_timeseries.csv"
  readr::read_csv(path, show_col_types = FALSE) |>
    dplyr::select(Year, Pop) |>
    dplyr::arrange(Year)
}

cod <- load_cod_data() |>
  mutate(t = Year - min(Year))  # centered time for fitting

# ---------------------------
# Helpers
# ---------------------------
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
logistic_fun <- function(t, K, A, r) K / (1 + A * exp(-r * t))

# Euler integration of dP/dt = r_eff(t)*P*(1 - P/K) - H_eff(t)
# H_eff(t) = 0 for Year < H_start_year, and = H otherwise
# Returns full series including r_eff(t) and bad_year flag
simulate_badyears <- function(P0, r, K, H, years, dt, year0,
                              by_p, bad_factor, seed = NULL,
                              H_start_year = -Inf,
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

    current_year <- year0 + t[i]
    H_eff <- if (current_year >= H_start_year) H else 0

    dP <- r_eff * P[i] * (1 - P[i] / K) - H_eff
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

pop_at_year <- function(df, target_year) {
  year_min <- min(df$Year, na.rm = TRUE)
  year_max <- max(df$Year, na.rm = TRUE)
  if (is.finite(year_min) && is.finite(year_max) &&
      target_year >= year_min && target_year <= year_max) {
    as.numeric(approx(x = df$Year, y = df$Pop, xout = target_year)$y)
  } else {
    NA_real_
  }
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
start_year <- min(cod$Year)
TARGET_YEAR <- 2050

# ---------------------------
# UI
# ---------------------------
ui <- fluidPage(
  titlePanel("Logistic with Harvesting, Bad-Year Shocks, and Pop(2050) Tracking"),
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
      sliderInput("H_start", "Harvest starts (calendar year):",
                  min = start_year,
                  max = start_year + 100,
                  value = start_year + 5,
                  step = 1),
      sliderInput("years", "Simulation horizon (years):",
                  min = 5, max = 120, value = 60, step = 1),
      hr(),
      strong("Bad-year process on r:"),
      sliderInput("by_p", "Bad-year probability p (per year):",
                  min = 0, max = 0.7, value = 0.25, step = 0.01),
      sliderInput("bad_factor", "Bad-year multiplier on r (e.g., 0.5 halves r):",
                  min = 0.05, max = 1.0, value = 0.5, step = 0.01),
      numericInput("seed", "Random seed (base)", value = 123, min = 1, step = 1),
      div(
        actionButton("resim", "Resimulate (1 run)"),
        actionButton("run10", "Run 10 simulations"),
        actionButton("run100", "Run 100 simulations"),
        style = "display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 8px; margin-top: 6px;"
      ),
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
      )),
      hr(),
      strong("Per-run capture:"),
      actionButton("clear_runs", "Clear 2050 points")
    ),
    mainPanel(
      plotOutput("popplot", height = 420),
      plotOutput("rplot",  height = 200),
      plotOutput("scatter_runs", height = 240),
      plotOutput("hist_runs", height = 260),
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

  # Store latest simulation for plotting the top two panels
  last_sim <- reactiveVal(NULL)

  # Storage for (run, pop2050) pairs
  runs_tbl <- reactiveVal(tibble(run = integer(), pop2050 = numeric()))

  # Ensure the simulation horizon reaches TARGET_YEAR for capture
  need_years_to_target <- reactive(TARGET_YEAR - start_year)

  run_one_sim <- function(seed_val) {
    years_span <- max(input$years, need_years_to_target())
    simulate_badyears(
      P0 = P0_fix,
      r  = r_fix,
      K  = K_fix,
      H  = input$H,
      years = years_span,
      dt = 0.1,
      year0 = start_year,
      by_p = input$by_p,
      bad_factor = input$bad_factor,
      seed = seed_val,
      H_start_year = input$H_start,
      clip_nonneg = TRUE
    )
  }

  append_point_target <- function(df) {
    p <- pop_at_year(df, TARGET_YEAR)
    cur <- runs_tbl()
    new_row <- tibble(run = nrow(cur) + 1L, pop2050 = p)
    runs_tbl(bind_rows(cur, new_row))
  }

  # Initial run on app load so plots aren't empty
  observeEvent(TRUE, {
    df0 <- run_one_sim(input$seed)
    last_sim(df0)
    append_point_target(df0)
  }, once = TRUE)

  # Single run
  observeEvent(input$resim, {
    seed_val <- input$seed + nrow(runs_tbl()) + 1L
    df <- run_one_sim(seed_val)
    last_sim(df)
    append_point_target(df)
  })

  # Batch: run 10 simulations and append 10 points
  observeEvent(input$run10, {
    base <- nrow(runs_tbl())
    last_df <- NULL
    for (i in 1:10) {
      seed_val <- input$seed + base + i
      df_i <- run_one_sim(seed_val)
      append_point_target(df_i)
      last_df <- df_i
    }
    if (!is.null(last_df)) last_sim(last_df)
  })

  # Batch: run 100 simulations and append 100 points (with progress bar)
  observeEvent(input$run100, {
    base <- nrow(runs_tbl())
    last_df <- NULL
    withProgress(message = "Running 100 simulations...", value = 0, {
      for (i in 1:100) {
        seed_val <- input$seed + base + i
        df_i <- run_one_sim(seed_val)
        append_point_target(df_i)
        last_df <- df_i
        incProgress(1/100)
      }
    })
    if (!is.null(last_df)) last_sim(last_df)
  })

  # Clear the accumulated points (does not clear the latest sim)
  observeEvent(input$clear_runs, {
    runs_tbl(tibble(run = integer(), pop2050 = numeric()))
  })

  # ----- Plots -----
  output$popplot <- renderPlot({
    df_sim <- last_sim()
    req(df_sim)

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
          " | Harvest starts = ", input$H_start,
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
    df_sim <- last_sim()
    req(df_sim)

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

  # Scatter: Simulation # vs Pop in 2050
  output$scatter_runs <- renderPlot({
    pts <- runs_tbl()
    ggplot(pts, aes(x = run, y = pop2050)) +
      geom_point(size = 2) +
      geom_line(alpha = 0.4) +
      labs(
        title = paste0("Population at Year ", TARGET_YEAR, " per Simulation"),
        x = "Simulation number",
        y = paste0("Population in ", TARGET_YEAR)
      ) +
      theme_classic()
  })

  # Histogram of Pop(2050) with Normal & Lognormal fits
  output$hist_runs <- renderPlot({
    pts <- runs_tbl() |> filter(is.finite(pop2050))
    validate(
      need(nrow(pts) >= 2, "Run more simulations to view histogram and fits.")
    )

    vals <- pts$pop2050

    # Fit Normal
    nf <- tryCatch(MASS::fitdistr(vals, "normal"), error = function(e) NULL)

    # Fit Lognormal (requires positive values)
    vals_pos <- vals[vals > 0]
    lf <- if (length(vals_pos) >= 2) {
      tryCatch(MASS::fitdistr(vals_pos, "lognormal"), error = function(e) NULL)
    } else NULL

    # Density grid
    xgrid <- seq(min(vals, na.rm = TRUE), max(vals, na.rm = TRUE), length.out = 400)
    dens_df <- tibble(x = xgrid)

    if (!is.null(nf)) {
      mu  <- nf$estimate["mean"]; sdv <- nf$estimate["sd"]
      dens_df <- dens_df |>
        mutate(norm = dnorm(x, mean = mu, sd = sdv))
    }

    if (!is.null(lf)) {
      mlog <- lf$estimate["meanlog"]; slog <- lf$estimate["sdlog"]
      dens_df <- dens_df |>
        mutate(lognorm = ifelse(x > 0, dlnorm(x, meanlog = mlog, sdlog = slog), NA_real_))
    }

    ggplot() +
      geom_histogram(data = pts, aes(x = pop2050, y = ..density..),
                     bins = 30, fill = "grey85", color = "grey40") +
      { if (!is.null(nf))
        geom_line(data = dens_df, aes(x = x, y = norm), linewidth = 1.2, color = "#0072B2") } +
      { if (!is.null(lf))
        geom_line(data = dens_df, aes(x = x, y = lognorm), linewidth = 1.2, linetype = "dashed", color = "#D55E00") } +
      labs(
        title = paste0("Distribution of Population in ", TARGET_YEAR),
        subtitle = paste(
          c(
            if (!is.null(nf)) sprintf("Normal fit: μ=%.3g, σ=%.3g", nf$estimate["mean"], nf$estimate["sd"]) else NULL,
            if (!is.null(lf)) sprintf("Lognormal fit: meanlog=%.3g, sdlog=%.3g", lf$estimate["meanlog"], lf$estimate["sdlog"]) else NULL
          ),
          collapse = "   |   "
        ),
        x = paste0("Population in ", TARGET_YEAR),
        y = "Density"
      ) +
      theme_classic()
  })
}

shinyApp(ui, server)
