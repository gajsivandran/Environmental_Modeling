# app.R — Sea-Tac Precip Explorer
# -------------------------------------------------------
# Tabs:
#  - Daily: time series with brush-zoom + non-zero daily histogram
#  - Annual: totals by year
#  - PDF: density view of annual totals
#  - Projection: fit parametric dist to annual totals, simulate 20-year path
#  - Ensemble: fit parametric dist, simulate runs × years; overlay distributions
#
# Data requirement:
#  - data/seatac_data.csv with columns DATE (m/d/Y) and PRCP (tenths of mm or mm)
#    We convert DATE with lubridate::mdy and PRCP <- as.numeric(PRCP) / 10 to mm

library(shiny)
library(tidyverse)
library(lubridate)
library(MASS)  # fitdistr

# ---------------- Helpers ----------------
`%||%` <- function(x, y) if (is.null(x)) y else x

safe_read <- function(path) {
  validate <- function(df) {
    req(all(c("DATE","PRCP") %in% names(df)))
    df
  }
  readr::read_csv(path, show_col_types = FALSE) |> validate()
}

clean_precip <- function(df) {
  df |>
    mutate(
      DATE = mdy(DATE),
      PRCP = suppressWarnings(as.numeric(PRCP) / 10)  # convert to mm
    ) |>
    filter(!is.na(DATE), !is.na(PRCP)) |>
    arrange(DATE)
}

aggregate_annual <- function(df) {
  df |>
    group_by(Year = year(DATE)) |>
    summarise(Annual_mm = sum(PRCP, na.rm = TRUE), .groups = "drop") |>
    arrange(Year)
}

# ---------------- UI ----------------
ui <- fluidPage(
  titlePanel("Sea-Tac Precipitation Explorer"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Controls"),
      uiOutput("date_ui"),
      tags$hr(),
      h5("Daily histogram (PRCP > 0)"),
      sliderInput("daily_bin", "Binwidth (mm):", min = 0.2, max = 10, value = 1, step = 0.2),
      tags$hr(),
      h5("Annual density (PDF)"),
      selectInput("kernel", "Kernel:",
                  choices = c("gaussian","epanechnikov","rectangular","triangular","biweight","cosine","optcosine"),
                  selected = "gaussian"),
      sliderInput("bw", "Bandwidth (mm):", min = 5, max = 100, value = 20, step = 1),
      helpText("Tip: Smaller bandwidth shows more detail; larger bandwidth smooths the curve.")
    ),
    mainPanel(
      tabsetPanel(
        id = "tabs",

        # ---------- Daily ----------
        tabPanel(
          "Daily",
          br(),
          fluidRow(
            column(
              width = 12,
              h4("Daily precipitation time series"),
              p("Drag on the plot to zoom. Double-click to reset."),
              plotOutput("ts_plot", height = 320,
                         brush = brushOpts(id = "ts_brush", resetOnNew = TRUE),
                         dblclick = "ts_dbl"),
              verbatimTextOutput("ts_stats")
            )
          ),
          tags$hr(),
          fluidRow(
            column(
              width = 12,
              h4("Histogram of wet-day precipitation (PRCP > 0)"),
              plotOutput("daily_hist", height = 300),
              verbatimTextOutput("daily_hist_stats")
            )
          )
        ),

        # ---------- Annual ----------
        tabPanel(
          "Annual",
          br(),
          uiOutput("year_ui"),
          fluidRow(
            column(
              width = 12,
              h4("Annual rainfall totals (mm)"),
              plotOutput("annual_plot", height = 300),
              verbatimTextOutput("annual_stats")
            )
          )
        ),

        # ---------- PDF of Annual ----------
        tabPanel(
          "PDF",
          br(),
          fluidRow(
            column(
              width = 12,
              h4("PDF of annual rainfall (histogram + density)"),
              plotOutput("annual_pdf", height = 320),
              verbatimTextOutput("annual_pdf_stats")
            )
          )
        ),

        # ---------- Projection (20-year path) ----------
        tabPanel(
          "Projection",
          br(),
          uiOutput("proj_year_ui"),
          fluidRow(
            column(
              width = 6,
              selectInput("fit_dist", "Distribution to fit:",
                          choices = c("Normal", "Lognormal", "Gamma"),
                          selected = "Lognormal")
            ),
            column(
              width = 3,
              numericInput("proj_seed", "Random seed:", value = 2025, min = 0, step = 1)
            ),
            column(
              width = 3,
              actionButton("proj_go", "Resimulate 20 years", class = "btn-primary", width = "100%")
            )
          ),
          p(em("Projection is fit to the currently selected historical year range.")),
          plotOutput("proj_plot", height = 360),
          verbatimTextOutput("proj_params"),
          verbatimTextOutput("proj_stats")
        ),

        # ---------- Ensemble (runs × years) ----------
        tabPanel(
          "Ensemble",
          br(),
          uiOutput("ens_year_ui"),
          fluidRow(
            column(
              width = 4,
              selectInput("ens_fit_dist", "Distribution to fit:",
                          choices = c("Normal", "Lognormal", "Gamma"),
                          selected = "Lognormal")
            ),
            column(
              width = 3,
              numericInput("ens_years", "Years per run:", value = 20, min = 1, step = 1)
            ),
            column(
              width = 3,
              sliderInput("ens_runs", "Number of runs:", min = 100, max = 5000, value = 1000, step = 100)
            ),
            column(
              width = 2,
              numericInput("ens_seed", "Random seed:", value = 2025, min = 0, step = 1)
            )
          ),
          actionButton("ens_go", "Simulate ensemble", class = "btn-primary"),
          p(em("Fits to the selected historical year range, then simulates runs × years annual totals.")),
          plotOutput("ens_overlay", height = 360),
          verbatimTextOutput("ens_params"),
          verbatimTextOutput("ens_stats")
        )

      )
    )
  )
)

# ---------------- Server ----------------
server <- function(input, output, session) {

  # -- Load & clean --
  raw <- reactive({
    path <- "data/seatac_data.csv"
    validate(need(file.exists(path), paste0("File not found: ", path)))
    safe_read(path)
  })

  daily <- reactive({
    clean_precip(raw())
  })

  # -- Date controls based on data extent --
  output$date_ui <- renderUI({
    d <- daily()
    req(nrow(d) > 0)
    dateRangeInput(
      "date_range", "Date range:",
      start = min(d$DATE, na.rm = TRUE),
      end   = max(d$DATE, na.rm = TRUE),
      min   = min(d$DATE, na.rm = TRUE),
      max   = max(d$DATE, na.rm = TRUE)
    )
  })

  # -- Filtered by date --
  daily_filtered <- reactive({
    d <- daily()
    req(input$date_range)
    d |> filter(DATE >= input$date_range[1], DATE <= input$date_range[2])
  })

  # -- Annual aggregation + year filter UI --
  annual_all <- reactive({
    aggregate_annual(daily())
  })

  output$year_ui <- renderUI({
    a <- annual_all()
    req(nrow(a) > 0)
    sliderInput(
      "year_range", "Year range:",
      min = min(a$Year, na.rm = TRUE),
      max = max(a$Year, na.rm = TRUE),
      value = c(min(a$Year, na.rm = TRUE), max(a$Year, na.rm = TRUE)),
      step = 1, sep = ""
    )
  })

  annual_filtered <- reactive({
    a <- annual_all()
    req(input$year_range)
    a |> filter(Year >= input$year_range[1], Year <= input$year_range[2])
  })

  # ---------------- Daily: Time series with brush-zoom ----------------
  ranges <- reactiveValues(x = NULL, y = NULL)

  observeEvent(input$ts_brush, {
    brush <- input$ts_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
    }
  })

  observeEvent(input$ts_dbl, {
    ranges$x <- NULL
    ranges$y <- NULL
  })

  output$ts_plot <- renderPlot({
    d <- daily_filtered()
    req(nrow(d) > 0)
    base <- ggplot(d, aes(DATE, PRCP)) +
      geom_line(linewidth = 0.5) +
      labs(x = "Date", y = "Precipitation (mm)", title = "Daily precipitation") +
      theme_minimal()
    if (!is.null(ranges$x) && !is.null(ranges$y)) {
      base + coord_cartesian(xlim = ranges$x, ylim = ranges$y)
    } else {
      base
    }
  })

  output$ts_stats <- renderText({
    d <- daily_filtered()
    sprintf(
      "Daily stats (selected range): n=%d | mean=%.2f mm | sd=%.2f mm | min=%.2f | max=%.2f",
      nrow(d), mean(d$PRCP), sd(d$PRCP), min(d$PRCP), max(d$PRCP)
    )
  })

  # ---------------- Daily: Histogram of PRCP > 0 ----------------
  output$daily_hist <- renderPlot({
    d <- daily_filtered() |> filter(PRCP > 0)
    req(nrow(d) > 0)
    ggplot(d, aes(PRCP)) +
      geom_histogram(binwidth = input$daily_bin, color = "white") +
      labs(x = "Precipitation (mm)", y = "Count",
           title = "Wet-day precipitation (PRCP > 0)") +
      theme_minimal()
  })

  output$daily_hist_stats <- renderText({
    d <- daily_filtered()
    wet <- d$PRCP[d$PRCP > 0]
    if (length(wet) == 0) return("No wet days in selection.")
    sprintf(
      "Wet-day stats: n=%d | mean=%.2f mm | sd=%.2f mm | min=%.2f | max=%.2f",
      length(wet), mean(wet), sd(wet), min(wet), max(wet)
    )
  })

  # ---------------- Annual: Totals plot ----------------
  output$annual_plot <- renderPlot({
    a <- annual_filtered()
    req(nrow(a) > 0)
    ggplot(a, aes(Year, Annual_mm)) +
      geom_col() +
      geom_line() +
      geom_point() +
      labs(x = "Year", y = "Annual precipitation (mm)",
           title = "Annual rainfall totals") +
      theme_minimal()
  })

  output$annual_stats <- renderText({
    a <- annual_filtered()
    sprintf(
      "Annual totals (selected years): n=%d | mean=%.1f mm | sd=%.1f mm | min=%.1f | max=%.1f",
      nrow(a), mean(a$Annual_mm), sd(a$Annual_mm),
      min(a$Annual_mm), max(a$Annual_mm)
    )
  })

  # ---------------- PDF of Annual totals ----------------
  output$annual_pdf <- renderPlot({
    a <- annual_filtered()
    req(nrow(a) >= 3)  # need a few points for density
    ggplot(a, aes(Annual_mm)) +
      geom_histogram(aes(y = after_stat(density)), bins = 15, alpha = 0.5) +
      geom_density(kernel = input$kernel, bw = input$bw, linewidth = 1) +
      labs(x = "Annual precipitation (mm)",
           y = "Density",
           title = "PDF of annual rainfall (histogram + density)") +
      theme_minimal()
  })

  output$annual_pdf_stats <- renderText({
    a <- annual_filtered()
    sprintf(
      "Density input: kernel=%s | bandwidth=%.1f mm | n(years)=%d | mean=%.1f mm | sd=%.1f mm",
      input$kernel, input$bw, nrow(a), mean(a$Annual_mm), sd(a$Annual_mm)
    )
  })

  # ---------------- Projection: Fit + simulate 20 years ----------------
  output$proj_year_ui <- renderUI({
    a <- annual_all()
    req(nrow(a) > 0)
    yrs <- range(a$Year, na.rm = TRUE)
    tagList(
      helpText(sprintf("Historical data span: %d–%d", yrs[1], yrs[2]))
    )
  })

  fit_params <- reactive({
    a <- annual_filtered()
    x <- a$Annual_mm
    validate(need(length(x) >= 5, "Need at least 5 years in the selected range to fit a distribution."))
    x_pos <- x[x > 0]
    dist <- input$fit_dist

    if (dist == "Normal") {
      fit <- MASS::fitdistr(x, "normal")
      list(dist = "Normal", params = as.list(fit$estimate))
    } else if (dist == "Lognormal") {
      validate(need(length(x_pos) >= 5, "Lognormal fit needs positive data (remove zero/negative years)."))
      fit <- MASS::fitdistr(x_pos, "lognormal")
      list(dist = "Lognormal", params = as.list(fit$estimate))
    } else {
      validate(need(length(x_pos) >= 5, "Gamma fit needs positive data (remove zero/negative years)."))
      fit <- MASS::fitdistr(x_pos, "gamma")
      list(dist = "Gamma", params = as.list(fit$estimate))
    }
  })

  projection <- eventReactive(input$proj_go, {
    set.seed(input$proj_seed %||% 2025)
    a <- annual_filtered()
    req(nrow(a) > 0)
    last_year <- max(a$Year, na.rm = TRUE)
    years_out <- 20
    proj_years <- seq.int(from = last_year + 1, by = 1, length.out = years_out)

    fp <- fit_params()
    dname <- fp$dist
    p <- fp$params

    sim <- switch(
      dname,
      "Normal"    = rnorm(years_out, mean = p$mean, sd = p$sd),
      "Lognormal" = rlnorm(years_out, meanlog = p$meanlog, sdlog = p$sdlog),
      "Gamma"     = rgamma(years_out, shape = p$shape, rate = p$rate)
    )
    tibble(Year = proj_years, Annual_mm = pmax(sim, 0))  # clamp negs for Normal
  }, ignoreInit = FALSE)

  output$proj_plot <- renderPlot({
    hist_a <- annual_filtered()
    req(nrow(hist_a) > 0)
    proj_a <- projection()

    y_max <- max(c(hist_a$Annual_mm, proj_a$Annual_mm), na.rm = TRUE)

    ggplot() +
      # Historical
      geom_col(data = hist_a, aes(Year, Annual_mm), alpha = 0.8) +
      geom_line(data = hist_a, aes(Year, Annual_mm)) +
      geom_point(data = hist_a, aes(Year, Annual_mm)) +
      # Projection background band
      {
        if (nrow(proj_a) > 0)
          annotate("rect",
                   xmin = min(proj_a$Year) - 0.5,
                   xmax = max(proj_a$Year) + 0.5,
                   ymin = 0, ymax = y_max,
                   alpha = 0.08)
      } +
      # Projected
      geom_col(data = proj_a, aes(Year, Annual_mm), alpha = 0.6) +
      geom_line(data = proj_a, aes(Year, Annual_mm), linetype = "dashed") +
      geom_point(data = proj_a, aes(Year, Annual_mm), shape = 21, stroke = 0.6) +
      labs(x = "Year", y = "Annual precipitation (mm)",
           title = "Historical vs. 20-year projected annual rainfall") +
      theme_minimal()
  })

  output$proj_params <- renderText({
    fp <- fit_params()
    if (fp$dist == "Normal") {
      sprintf("Fitted distribution: Normal(mean = %.1f, sd = %.1f)", fp$params$mean, fp$params$sd)
    } else if (fp$dist == "Lognormal") {
      sprintf("Fitted distribution: Lognormal(meanlog = %.3f, sdlog = %.3f)", fp$params$meanlog, fp$params$sdlog)
    } else {
      sprintf("Fitted distribution: Gamma(shape = %.3f, rate = %.6f)", fp$params$shape, fp$params$rate)
    }
  })

  output$proj_stats <- renderText({
    proj_a <- projection()
    sprintf("Projection summary (20 years): mean = %.1f mm | sd = %.1f mm | min = %.1f | max = %.1f",
            mean(proj_a$Annual_mm), sd(proj_a$Annual_mm),
            min(proj_a$Annual_mm), max(proj_a$Annual_mm))
  })

  # ---------------- Ensemble: fit params just for this tab ----------------
  output$ens_year_ui <- renderUI({
    a <- annual_all()
    req(nrow(a) > 0)
    yrs <- range(a$Year, na.rm = TRUE)
    tagList(
      helpText(sprintf("Historical data span: %d–%d", yrs[1], yrs[2]))
    )
  })

  ens_fit_params <- reactive({
    a <- annual_filtered()
    x <- a$Annual_mm
    validate(need(length(x) >= 5, "Need at least 5 years in the selected range to fit a distribution."))
    x_pos <- x[x > 0]
    dist <- input$ens_fit_dist

    if (dist == "Normal") {
      fit <- MASS::fitdistr(x, "normal")
      list(dist = "Normal", params = as.list(fit$estimate))
    } else if (dist == "Lognormal") {
      validate(need(length(x_pos) >= 5, "Lognormal fit needs positive annual totals."))
      fit <- MASS::fitdistr(x_pos, "lognormal")
      list(dist = "Lognormal", params = as.list(fit$estimate))
    } else {
      validate(need(length(x_pos) >= 5, "Gamma fit needs positive annual totals."))
      fit <- MASS::fitdistr(x_pos, "gamma")
      list(dist = "Gamma", params = as.list(fit$estimate))
    }
  })

  ens_sims <- eventReactive(input$ens_go, {
    set.seed(input$ens_seed %||% 2025)
    runs  <- input$ens_runs  %||% 1000
    years <- input$ens_years %||% 20
    validate(need(runs > 0 && years > 0, "Runs and years must be positive."))

    fp <- ens_fit_params()
    dname <- fp$dist
    p <- fp$params
    n <- runs * years

    sim_vals <- switch(
      dname,
      "Normal"    = pmax(rnorm(n, mean = p$mean, sd = p$sd), 0),          # clamp negatives to 0
      "Lognormal" = rlnorm(n, meanlog = p$meanlog, sdlog = p$sdlog),
      "Gamma"     = rgamma(n, shape = p$shape, rate = p$rate)
    )

    tibble(Annual_mm = sim_vals)
  }, ignoreInit = TRUE)

  output$ens_overlay <- renderPlot({
    hist_a <- annual_filtered()
    req(nrow(hist_a) > 0)

    sims <- ens_sims()
    df_hist <- hist_a |> transmute(Annual_mm, Source = "Historical")
    df_sim  <- sims   |> transmute(Annual_mm, Source = "Simulated")
    df_all <- bind_rows(df_hist, df_sim)

    ggplot(df_all, aes(x = Annual_mm, fill = Source)) +
      geom_histogram(aes(y = after_stat(density)), bins = 25, alpha = 0.35, position = "identity") +
      geom_density(linewidth = 1) +
      labs(
        x = "Annual precipitation (mm)",
        y = "Density",
        title = "Distribution of annual rainfall: Historical vs Simulated Ensemble"
      ) +
      theme_minimal()
  })

  output$ens_params <- renderText({
    fp <- ens_fit_params()
    if (fp$dist == "Normal") {
      sprintf("Fitted to historical annual totals: Normal(mean = %.1f, sd = %.1f)",
              fp$params$mean, fp$params$sd)
    } else if (fp$dist == "Lognormal") {
      sprintf("Fitted to historical annual totals: Lognormal(meanlog = %.3f, sdlog = %.3f)",
              fp$params$meanlog, fp$params$sdlog)
    } else {
      sprintf("Fitted to historical annual totals: Gamma(shape = %.3f, rate = %.6f)",
              fp$params$shape, fp$params$rate)
    }
  })

  output$ens_stats <- renderText({
    hist_a <- annual_filtered()
    sims   <- ens_sims()
    if (is.null(sims) || nrow(sims) == 0) {
      return(sprintf(
        "Historical: n=%d | mean=%.1f | sd=%.1f | min=%.1f | max=%.1f",
        nrow(hist_a), mean(hist_a$Annual_mm), sd(hist_a$Annual_mm),
        min(hist_a$Annual_mm), max(hist_a$Annual_mm)
      ))
    }
    sprintf(
      paste0(
        "Historical: n=%d | mean=%.1f | sd=%.1f | min=%.1f | max=%.1f\n",
        "Simulated (runs×years = %d×%d = %d): mean=%.1f | sd=%.1f | ",
        "min=%.1f | 5th=%.1f | 50th=%.1f | 95th=%.1f | max=%.1f"
      ),
      nrow(hist_a), mean(hist_a$Annual_mm), sd(hist_a$Annual_mm),
      min(hist_a$Annual_mm), max(hist_a$Annual_mm),
      input$ens_runs, input$ens_years, input$ens_runs * input$ens_years,
      mean(sims$Annual_mm), sd(sims$Annual_mm),
      min(sims$Annual_mm),
      quantile(sims$Annual_mm, 0.05),
      quantile(sims$Annual_mm, 0.50),
      quantile(sims$Annual_mm, 0.95),
      max(sims$Annual_mm)
    )
  })
}

# ---------------- Run ----------------
shinyApp(ui, server)
