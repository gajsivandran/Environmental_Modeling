# app.R — Rainfall-Proportional Logistic Forest Model (with mortality)
# -------------------------------------------------------------------
# Features:
# ✅ Tabs 1–6 (Gamma fits, Logistic, Rainfall Dependence scenarios)
# ✅ r_t = r0 * min(1 + excess, 1.2) - m * deficit  (fixed historical mean reference)
#    where excess  = max(R/Rbar - 1, 0), deficit = max(1 - R/Rbar, 0)
# ✅ User tunes drought mortality m (die-off when severe drought)
# ✅ Reference curve overlay (no rainfall dependence)
# ✅ Rainfall statistics displayed under plots
# ✅ Persistent multi-run comparison + Reset in Tab 6
# ✅ All rainfall-related and simulated values forced ≥ 0
# -------------------------------------------------------------------

library(shiny)
library(tidyverse)
library(lubridate)
library(readr)
library(scales)

# ---------------------------
# Data load (SeaTac rainfall)
# ---------------------------
DATA_PATH <- "data/seatac_data.csv"
stopifnot(file.exists(DATA_PATH))

raw0 <- read_csv(
  DATA_PATH,
  col_types = cols(.default = col_character()),
  na = c("", "NA", "NaN"),
  show_col_types = FALSE
) |>
  rename_with(toupper)

date_parsed <- parse_date_time(
  raw0$DATE,
  orders = c("Y-m-d", "m/d/Y", "Y/m/d", "d/m/Y", "Ymd"),
  quiet = TRUE
)
DATE <- as.Date(date_parsed)
PRCP <- parse_number(raw0$PRCP) / 10  # tenths of mm → mm

data <- tibble(DATE = DATE, PRCP = PRCP) |>
  filter(!is.na(DATE)) |>
  mutate(PRCP = pmax(PRCP, 0)) |>  # ensure nonnegative
  arrange(DATE)

annual_mm <- data |>
  mutate(Year = year(DATE)) |>
  group_by(Year) |>
  summarise(AnnualPRCP_mm = sum(PRCP, na.rm = TRUE), .groups = "drop") |>
  filter(!is.na(AnnualPRCP_mm), AnnualPRCP_mm >= 0) |>
  arrange(Year)

# Fixed historical mean (used in ALL scenarios)
Rbar_ref <- mean(annual_mm$AnnualPRCP_mm)

# ---------------------------
# Helper functions
# ---------------------------
gamma_ms_to_shape_rate <- function(mu, sd) {
  mu <- max(mu, 0)
  sd <- max(sd, 1e-6)
  shape <- (mu / sd)^2
  rate  <- mu / (sd^2)
  list(shape = shape, rate = rate)
}

logistic_B <- function(t, r, K, B0) {
  K / (1 + ((K - B0) / B0) * exp(-r * t))
}

logistic_dynamic <- function(r_series, K, B0, dt = 1) {
  n <- length(r_series)
  B <- numeric(n)
  B[1] <- B0
  for (t in 2:n) {
    B[t] <- B[t - 1] + r_series[t - 1] * B[t - 1] * (1 - B[t - 1] / K) * dt
    B[t] <- max(B[t], 0) # no negatives
  }
  tibble(Year = seq_along(B), B = B)
}

spell_lengths <- function(prcp_vec, wet_threshold = 0) {
  wet_flag <- prcp_vec > wet_threshold
  r <- rle(wet_flag)
  wet_lengths <- pmax(r$lengths[r$values], 0)
  dry_lengths <- pmax(r$lengths[!r$values], 0)
  list(wet = wet_lengths, dry = dry_lengths)
}

simulate_daily_rain <- function(n_years = 100, mu_vol, sd_vol, mu_wet, sd_wet, mu_dry, sd_dry) {
  mu_vol <- max(mu_vol, 0); sd_vol <- max(sd_vol, 1e-6)
  mu_wet <- max(mu_wet, 0); sd_wet <- max(sd_wet, 1e-6)
  mu_dry <- max(mu_dry, 0); sd_dry <- max(sd_dry, 1e-6)

  pr_vol <- gamma_ms_to_shape_rate(mu_vol, sd_vol)
  pr_wet <- gamma_ms_to_shape_rate(mu_wet, sd_wet)
  pr_dry <- gamma_ms_to_shape_rate(mu_dry, sd_dry)
  n_days <- n_years * 365
  rain <- c()
  while (length(rain) < n_days) {
    dry <- max(rgamma(1, pr_dry$shape, pr_dry$rate), 0)
    wet <- max(rgamma(1, pr_wet$shape, pr_wet$rate), 0)
    dry_days <- rep(0, round(dry))
    wet_days <- pmax(rgamma(round(wet), pr_vol$shape, pr_vol$rate), 0)
    rain <- c(rain, dry_days, wet_days)
  }
  rain <- rain[1:n_days]
  tibble(Day = 1:n_days, Rain = rain, Year = ceiling((1:n_days) / 365))
}

# --- Rainfall → growth coupling (Option 4: proportional + mortality, with cap) ---
# r_t = r0 * min(1 + excess, cap) - m * deficit
# where ratio = R / Rbar_ref, excess = max(ratio - 1, 0), deficit = max(1 - ratio, 0)
rain_adjusted_growth <- function(R, r0 = 0.2, Rbar_ref, m = 0.6, cap = 1.2) {
  R <- pmax(R, 0)
  ratio  <- R / Rbar_ref
  excess <- pmax(ratio - 1, 0)
  deficit <- pmax(1 - ratio, 0)       # fractional shortfall
  grow <- r0 * pmin(1 + excess, cap)  # cap benefit at +20%
  r <- grow - m * deficit             # mortality can push negative
  r
}

# ---------------------------
# Empirical vectors
# ---------------------------
spells <- spell_lengths(data$PRCP, wet_threshold = 0)
wet_day_volumes <- pmax(data$PRCP[data$PRCP > 0], 0)
wet_durations <- pmax(spells$wet, 0)
interstorm_durations <- pmax(spells$dry, 0)

# ---------------------------
# UI
# ---------------------------
ui <- navbarPage(
  title = "Rainfall → Forest Model",
  id = "main_nav",

  # TAB 1 — Distributions
  tabPanel(
    title = "Distributions",
    fluidPage(
      h3("Step 1 — Empirical Distributions & Gamma Fits"),
      plotOutput("plot_volume", height = "250px"),
      plotOutput("plot_wet_duration", height = "250px"),
      plotOutput("plot_dry_duration", height = "250px")
    )
  ),

  # TAB 2 — Logistic Growth
  tabPanel(
    title = "Logistic Growth",
    fluidPage(
      h3("Step 2 — Simple Logistic Forest Growth Model"),
      sidebarLayout(
        sidebarPanel(
          sliderInput("r", "Growth rate (r)", min = 0.01, max = 1, value = 0.2, step = 0.01),
          sliderInput("K", "Carrying capacity (K)", min = 50, max = 1000, value = 500, step = 10),
          sliderInput("B0", "Initial biomass (B0)", min = 1, max = 500, value = 50, step = 5),
          sliderInput("tmax", "Simulation time (years)", min = 10, max = 200, value = 100, step = 10)
        ),
        mainPanel(plotOutput("logistic_plot", height = "400px"))
      )
    )
  ),

  # TAB 3 — Rainfall Dependence
  tabPanel(
    title = "Rainfall Dependence",
    fluidPage(
      h3("Step 3 — Rainfall-Driven Forest Growth Scenarios"),
      withMathJax(),
      wellPanel(
        HTML('
    <h4>Model Background: Logistic Forest Growth</h4>
    <p>The forest biomass \\(B(t)\\) follows a <strong>logistic growth model</strong> driven by rainfall-dependent growth rate \\(r_t\\):</p>
    <p style="text-align:center;">
    $$B(t) = \\frac{K}{1 + \\left(\\frac{K - B_0}{B_0}\\right)e^{-r_t t}}$$
    </p>
    <p>The corresponding <strong>differential form</strong> of the logistic equation is:</p>
    <p style="text-align:center;">
    $$\\frac{dB}{dt} = r_t B \\left(1 - \\frac{B}{K}\\right)$$
    </p>

    <p>Here:</p>
    <ul>
      <li>\\(r_t\\) = intrinsic growth rate (may vary annually with rainfall)</li>
      <li>\\(K\\) = carrying capacity of the forest</li>
      <li>\\(B_0\\) = initial biomass</li>
    </ul>

    <hr/>
    <h4>Rainfall Dependence</h4>
    <p>We link the intrinsic rate \\(r_t\\) to annual rainfall \\(R_t\\) relative to a fixed historical mean \\(\\overline{R}\\):</p>
    <p style="text-align:center;">
    $$r_t = r_0\\,\\min\\!\\big(1 + \\max(\\tfrac{R_t}{\\overline{R}}-1,\\,0),\\,1.2\\big)\\; -\\; m\\,\\max\\!\\big(1-\\tfrac{R_t}{\\overline{R}},\\,0\\big)$$
    </p>
    <p>This formulation caps wet-year benefits at +20% of \\(r_0\\) and allows strong drought-driven mortality (negative \\(r_t\\)), controlled by parameter \\(m\\).</p>
  '),
        sliderInput("m", "Drought mortality strength (m)",
                    min = 0, max = 1.0, value = 0.6, step = 0.05, width = "60%")
      ),

      tabsetPanel(
        tabPanel("1️⃣ No Rainfall Dependence", plotOutput("plot_no_rain", height = "350px")),
        tabPanel("2️⃣ Historical Rainfall", plotOutput("plot_hist_rain", height = "350px"), verbatimTextOutput("stats_hist")),
        tabPanel("3️⃣ Simulated Annual Rainfall", plotOutput("plot_sim_rain", height = "350px"), verbatimTextOutput("stats_sim")),
        tabPanel("4️⃣ Adjust Annual Gamma Parameters",
                 sidebarLayout(
                   sidebarPanel(
                     sliderInput("mu_pct", "Change mean (%)", -50, 100, 0, 5),
                     sliderInput("sd_pct", "Change sd (%)", -50, 100, 0, 5),
                     actionButton("regen", "Regenerate Simulation")
                   ),
                   mainPanel(plotOutput("plot_adjusted", height = "350px"), verbatimTextOutput("stats_adjusted"))
                 )
        ),
        tabPanel("5️⃣ Synthetic Daily Rainfall",
                 fluidPage(
                   p("Generates synthetic daily rainfall from interstorm, wet-spell, and wet-day Gamma distributions, then aggregates to annual rainfall."),
                   plotOutput("plot_synthetic", height = "350px"),
                   verbatimTextOutput("stats_synthetic")
                 )
        ),
        # ----- UI for Tab 6 -----
        tabPanel("6️⃣ Adjust Event Distributions",
                 sidebarLayout(
                   sidebarPanel(
                     h5("Wet-day volume (mm)"),
                     sliderInput("mu_vol", "Mean (%)", -50, 100, 0, 5),
                     sliderInput("sd_vol", "SD (%)", -50, 100, 0, 5),
                     h5("Wet-spell duration (days)"),
                     sliderInput("mu_wet", "Mean (%)", -50, 100, 0, 5),
                     sliderInput("sd_wet", "SD (%)", -50, 100, 0, 5),
                     h5("Interstorm duration (days)"),
                     sliderInput("mu_dry", "Mean (%)", -50, 100, 0, 5),
                     sliderInput("sd_dry", "SD (%)", -50, 100, 0, 5),

                     tags$hr(),
                     actionButton("regen6", "Add 1 Simulation"),
                     actionButton("regen6_10", "Add 10 Simulations"),
                     actionButton("reset6", "Reset Plot", class = "btn-danger"),

                     tags$hr(),
                     sliderInput("year_x", "Year X for population distribution", min = 1, max = 100, value = 50, step = 1)
                   ),
                   mainPanel(
                     plotOutput("plot_adjusted_synth", height = "350px"),
                     verbatimTextOutput("stats_adjusted_synth"),
                     tags$hr(),
                     h4("Distribution of Biomass at Year X (across simulations)"),
                     plotOutput("hist_yearX", height = "300px"),
                     verbatimTextOutput("stats_yearX")
                   )
                 )
        )

      )
    )
  )
)

# ---------------------------
# SERVER
# ---------------------------
server <- function(input, output, session) {
  # --- Tab 1: Gamma Fits ---
  output$plot_volume <- renderPlot({
    x <- pmax(wet_day_volumes, 0)
    mu <- mean(x); sdv <- sd(x); pr <- gamma_ms_to_shape_rate(mu, sdv)
    ggplot(tibble(x), aes(x)) +
      geom_histogram(aes(y = after_stat(density)), bins = 40, alpha = 0.5) +
      stat_function(fun = function(z) dgamma(pmax(z, 0), pr$shape, pr$rate)) +
      theme_minimal() + labs(x = "Wet-day rainfall (mm)", y = "Density")
  })

  output$plot_wet_duration <- renderPlot({
    x <- pmax(wet_durations, 0)
    mu <- mean(x); sdv <- sd(x); pr <- gamma_ms_to_shape_rate(mu, sdv)
    ggplot(tibble(x), aes(x)) +
      geom_histogram(aes(y = after_stat(density)), bins = 40, alpha = 0.5) +
      stat_function(fun = function(z) dgamma(pmax(z, 0), pr$shape, pr$rate)) +
      theme_minimal() + labs(x = "Wet-spell duration (days)", y = "Density")
  })

  output$plot_dry_duration <- renderPlot({
    x <- pmax(interstorm_durations, 0)
    mu <- mean(x); sdv <- sd(x); pr <- gamma_ms_to_shape_rate(mu, sdv)
    ggplot(tibble(x), aes(x)) +
      geom_histogram(aes(y = after_stat(density)), bins = 40, alpha = 0.5) +
      stat_function(fun = function(z) dgamma(pmax(z, 0), pr$shape, pr$rate)) +
      theme_minimal() + labs(x = "Interstorm duration (days)", y = "Density")
  })

  # --- Tab 2: Logistic Growth ---
  output$logistic_plot <- renderPlot({
    t <- seq(0, input$tmax, 1)
    B <- pmax(logistic_B(t, input$r, input$K, input$B0), 0)
    ggplot(tibble(t, B), aes(t, B)) +
      geom_line(color = "steelblue", linewidth = 1.2) +
      geom_hline(yintercept = input$K, linetype = "dashed", color = "gray50") + theme_minimal()
  })

  # --- Reference curve (no rainfall) ---
  ref_curve <- reactive({
    t <- seq(0, 100, 1)
    tibble(t, Bref = logistic_B(t, 0.2, 500, 50))
  })

  # --- Tabs 3–6: Rainfall dependence ---
  output$plot_no_rain <- renderPlot({
    ggplot(ref_curve(), aes(t, Bref)) +
      geom_line(color = "darkgreen", linewidth = 1.2) + theme_minimal() +
      labs(title = "No rainfall dependence (reference)", x = "Year", y = "Biomass (B)")
  })

  # Historical rainfall
  output$plot_hist_rain <- renderPlot({
    rain <- pmax(annual_mm$AnnualPRCP_mm, 0)
    r_t <- rain_adjusted_growth(rain, r0 = 0.2, Rbar_ref = Rbar_ref, m = input$m)
    df <- logistic_dynamic(r_t, 500, 50)
    ggplot() +
      geom_line(data = ref_curve(), aes(t, Bref), color = "gray60", linetype = "dashed") +
      geom_line(data = df, aes(Year, B), color = "dodgerblue", linewidth = 1.2) +
      theme_minimal() + labs(title = "Historical rainfall-driven growth", x = "Year", y = "Biomass (B)")
  })
  output$stats_hist <- renderText({
    paste0("Historical rainfall — mean = ", round(mean(pmax(annual_mm$AnnualPRCP_mm, 0)), 1),
           " mm, sd = ", round(sd(pmax(annual_mm$AnnualPRCP_mm, 0)), 1),
           " (reference mean R̄ fixed at ", round(Rbar_ref, 1), " mm)")
  })

  # Simulated annual rainfall (Gamma)
  output$plot_sim_rain <- renderPlot({
    mu <- mean(pmax(annual_mm$AnnualPRCP_mm, 0))
    sdv <- sd(pmax(annual_mm$AnnualPRCP_mm, 0))
    pr <- gamma_ms_to_shape_rate(mu, sdv)
    sim_rain <- pmax(rgamma(nrow(annual_mm), pr$shape, pr$rate), 0)
    r_t <- rain_adjusted_growth(sim_rain, r0 = 0.2, Rbar_ref = Rbar_ref, m = input$m)
    df <- logistic_dynamic(r_t, 500, 50)
    ggplot() +
      geom_line(data = ref_curve(), aes(t, Bref), color = "gray60", linetype = "dashed") +
      geom_line(data = df, aes(Year, B), color = "orange", linewidth = 1.2) + theme_minimal() +
      labs(title = "Simulated annual rainfall (Gamma)", x = "Year", y = "Biomass (B)")
  })
  output$stats_sim <- renderText({
    mu <- mean(pmax(annual_mm$AnnualPRCP_mm, 0))
    sdv <- sd(pmax(annual_mm$AnnualPRCP_mm, 0))
    pr <- gamma_ms_to_shape_rate(mu, sdv)
    sim_rain <- pmax(rgamma(10000, pr$shape, pr$rate), 0)
    paste0("Observed mean = ", round(mu, 1), " mm (sd ", round(sdv, 1),
           "); Simulated (Gamma) mean ≈ ", round(mean(sim_rain), 1),
           " mm (sd ", round(sd(sim_rain), 1), "). R̄ (reference) = ", round(Rbar_ref, 1), " mm.")
  })

  # Adjusted annual rainfall (Tab 4)
  observeEvent(input$regen, {
    mu <- max(mean(pmax(annual_mm$AnnualPRCP_mm, 0)) * (1 + input$mu_pct / 100), 0)
    sdv <- max(sd(pmax(annual_mm$AnnualPRCP_mm, 0)) * (1 + input$sd_pct / 100), 1e-6)
    pr <- gamma_ms_to_shape_rate(mu, sdv)
    sim_rain <- pmax(rgamma(nrow(annual_mm), pr$shape, pr$rate), 0)
    r_t <- rain_adjusted_growth(sim_rain, r0 = 0.2, Rbar_ref = Rbar_ref, m = input$m)
    df <- logistic_dynamic(r_t, 500, 50)
    output$plot_adjusted <- renderPlot({
      ggplot() +
        geom_line(data = ref_curve(), aes(t, Bref), color = "gray60", linetype = "dashed") +
        geom_line(data = df, aes(Year, B), color = "purple", linewidth = 1.2) + theme_minimal() +
        labs(title = "Adjusted annual rainfall → forest growth", x = "Year", y = "Biomass (B)")
    })
    output$stats_adjusted <- renderText({
      paste0("Adjusted annual rainfall — mean = ", round(mu, 1),
             " mm, sd = ", round(sdv, 1), " mm. R̄ (reference) = ", round(Rbar_ref, 1), " mm.")
    })
  })

  # Synthetic daily rainfall (Tab 5)
  output$plot_synthetic <- renderPlot({
    mu_vol <- max(mean(wet_day_volumes), 0); sd_vol <- max(sd(wet_day_volumes), 1e-6)
    mu_wet <- max(mean(wet_durations), 0);   sd_wet <- max(sd(wet_durations), 1e-6)
    mu_dry <- max(mean(interstorm_durations), 0); sd_dry <- max(sd(interstorm_durations), 1e-6)
    daily_sim <- simulate_daily_rain(100, mu_vol, sd_vol, mu_wet, sd_wet, mu_dry, sd_dry)
    annual_sim <- daily_sim |> group_by(Year) |> summarise(Annual = pmax(sum(Rain), 0))
    r_t <- rain_adjusted_growth(annual_sim$Annual, r0 = 0.2, Rbar_ref = Rbar_ref, m = input$m)
    df <- logistic_dynamic(r_t, 500, 50)
    ggplot() +
      geom_line(data = ref_curve(), aes(t, Bref), color = "gray60", linetype = "dashed") +
      geom_line(data = df, aes(Year, B), color = "darkred", linewidth = 1.2) + theme_minimal() +
      labs(title = "Synthetic daily rainfall → forest growth", x = "Year", y = "Biomass (B)")
  })
  output$stats_synthetic <- renderText({
    paste0("Synthetic driver built from daily spell sampling. Historical annual mean = ",
           round(mean(annual_mm$AnnualPRCP_mm), 1), " mm; R̄ (reference) = ",
           round(Rbar_ref, 1), " mm.")
  })

  # ----- Tab 6: Adjustable event distributions (persistent + 1/10 add + reset + Year-X hist) -----
  rv <- reactiveValues(df = tibble(), stats = tibble())

  simulate_one_run <- function(mu_vol, sd_vol, mu_wet, sd_wet, mu_dry, sd_dry) {
    daily_sim <- simulate_daily_rain(100, mu_vol, sd_vol, mu_wet, sd_wet, mu_dry, sd_dry)
    annual_sim <- daily_sim |> dplyr::group_by(Year) |> dplyr::summarise(Annual = pmax(sum(Rain), 0), .groups = "drop")

    rain_mean <- mean(annual_sim$Annual)
    rain_sd   <- sd(annual_sim$Annual)

    r_t <- rain_adjusted_growth(annual_sim$Annual, r0 = 0.2, Rbar_ref = Rbar_ref, m = input$m)
    df  <- logistic_dynamic(r_t, 500, 50)

    list(df = df, rain_mean = rain_mean, rain_sd = rain_sd)
  }

  # Helper to append a run
  append_run <- function(run_result) {
    run_id <- ifelse(nrow(rv$df) == 0, 1, max(rv$df$run) + 1)
    df_new <- run_result$df
    df_new$run <- run_id
    rv$df <- dplyr::bind_rows(rv$df, df_new)
    rv$stats <- dplyr::bind_rows(
      rv$stats,
      tibble::tibble(run = run_id,
                     rain_mean = run_result$rain_mean,
                     rain_sd   = run_result$rain_sd)
    )
    invisible(run_id)
  }

  # Build means/SDs from adjusted sliders
  get_adjusted_event_params <- reactive({
    mu_vol <- max(mean(wet_day_volumes) * (1 + input$mu_vol / 100), 0)
    sd_vol <- max(sd(wet_day_volumes)   * (1 + input$sd_vol / 100), 1e-6)
    mu_wet <- max(mean(wet_durations)   * (1 + input$mu_wet / 100), 0)
    sd_wet <- max(sd(wet_durations)     * (1 + input$sd_wet / 100), 1e-6)
    mu_dry <- max(mean(interstorm_durations) * (1 + input$mu_dry / 100), 0)
    sd_dry <- max(sd(interstorm_durations)   * (1 + input$sd_dry / 100), 1e-6)
    list(mu_vol=mu_vol, sd_vol=sd_vol, mu_wet=mu_wet, sd_wet=sd_wet, mu_dry=mu_dry, sd_dry=sd_dry)
  })

  # Add 1 simulation
  observeEvent(input$regen6, {
    p <- get_adjusted_event_params()
    run_result <- simulate_one_run(p$mu_vol, p$sd_vol, p$mu_wet, p$sd_wet, p$mu_dry, p$sd_dry)
    append_run(run_result)
  })

  # Add 10 simulations
  observeEvent(input$regen6_10, {
    p <- get_adjusted_event_params()
    for (i in 1:10) {
      run_result <- simulate_one_run(p$mu_vol, p$sd_vol, p$mu_wet, p$sd_wet, p$mu_dry, p$sd_dry)
      append_run(run_result)
    }
  })

  # Reset button: clears previous simulations and plots
  observeEvent(input$reset6, {
    rv$df <- tibble::tibble()
    rv$stats <- tibble::tibble()
  })

  # Plot all time-series with reference
  output$plot_adjusted_synth <- renderPlot({
    ggplot() +
      geom_line(data = ref_curve(), aes(t, Bref), color = "gray60", linetype = "dashed") +
      { if (nrow(rv$df) > 0) geom_line(data = rv$df, aes(Year, B, group = run, color = factor(run)), alpha = 0.8) } +
      theme_minimal() +
      labs(title = "Adjusted Synthetic Rainfall → Forest Growth",
           x = "Year", y = "Biomass (B)") +
      guides(color = "none")
  })

  # Text summary of latest run (or cleared)
  output$stats_adjusted_synth <- renderText({
    if (nrow(rv$stats) == 0) return("No simulations yet. Use 'Add 1 Simulation' or 'Add 10 Simulations'.")
    last <- rv$stats |> dplyr::slice(dplyr::n())
    paste0("Latest simulation — annual rainfall mean = ", round(last$rain_mean, 1),
           " mm, sd = ", round(last$rain_sd, 1), " mm. R̄ (reference) = ",
           round(Rbar_ref, 1), " mm. Total runs: ", nrow(rv$stats), ".")
  })

  # Histogram of biomass at selected Year X across all runs
  output$hist_yearX <- renderPlot({
    req(nrow(rv$df) > 0)
    yr <- input$year_x
    df_year <- rv$df |> dplyr::filter(Year == yr)
    validate(need(nrow(df_year) > 0, "No data at this year yet — add simulations."))

    muB <- mean(df_year$B); sdB <- sd(df_year$B)

    ggplot(df_year, aes(x = pmax(B, 0))) +
      geom_histogram(aes(y = after_stat(density)), bins = 20, alpha = 0.5) +
      geom_density(linewidth = 1.1) +
      stat_function(fun = function(x) dnorm(x, mean = muB, sd = ifelse(sdB > 0, sdB, 1e-6)),
                    linetype = "dashed") +
      theme_minimal() +
      labs(title = paste0("Biomass distribution at Year ", yr, " (across runs)"),
           x = "Biomass (B)", y = "Density")
  })

  # Stats for Year X
  output$stats_yearX <- renderText({
    if (nrow(rv$df) == 0) return("No simulations yet.")
    yr <- input$year_x
    df_year <- rv$df |> dplyr::filter(Year == yr)
    if (nrow(df_year) == 0) return("No data at this year yet — add simulations.")
    paste0("Runs = ", dplyr::n_distinct(df_year$run),
           " | Mean B = ", round(mean(df_year$B), 2),
           " | SD B = ", round(sd(df_year$B), 2),
           " | Min B = ", round(min(df_year$B), 2),
           " | Max B = ", round(max(df_year$B), 2))
  })

}

# ---------------------------
# Run Shiny App
# ---------------------------
shinyApp(ui, server)
