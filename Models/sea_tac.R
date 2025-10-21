# app.R — Sea-Tac PRCP explorer
# - Daily time series with wet-day metrics
# - Annual totals with 5-yr moving average + decadal means
# - Distribution comparison across two periods (with zoom)
# - Monthly tab (choose a month across years), also with zoom on distributions

library(shiny)
library(tidyverse)
library(lubridate)
library(readr)
library(scales)
library(zoo)       # rollmean
library(plotly)    # interactive zoom for distributions

# ---------------------------
# Data load
# ---------------------------
DATA_PATH <- "data/seatac_data.csv"
stopifnot(file.exists(DATA_PATH))

# Parse DATE as m/d/YYYY, PRCP numeric
raw <- read_csv(
  DATA_PATH,
  col_types = cols(
    DATE = col_date(format = "%m/%d/%Y"),
    PRCP = col_double(),
    .default = col_guess()
  )
)

data <- raw |>
  select(DATE, PRCP) |>
  mutate(
     PRCP = suppressWarnings(as.numeric(PRCP))/10
  ) |>
  filter(!is.na(DATE)) |>
  arrange(DATE)

date_min <- min(data$DATE, na.rm = TRUE)
date_max <- max(data$DATE, na.rm = TRUE)

# Helpers
filter_zero <- function(df, exclude_zero) {
  if (isTRUE(exclude_zero)) dplyr::filter(df, PRCP > 0) else df
}
as_num_safely <- function(x, default = 0) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), default, x)
}

# ---------------------------
# UI
# ---------------------------
ui <- fluidPage(
  titlePanel("Sea-Tac Rainfall — PRCP (mm)"),
  tabsetPanel(
    id = "tabs",

    # -------- Daily Time Series --------
    tabPanel(
      "Daily Time Series",
      sidebarLayout(
        sidebarPanel(
          dateRangeInput(
            "dates", "Date range:",
            start = date_min, end = date_max,
            min = date_min, max = date_max
          ),
          numericInput(
            "wet_thr", "Wet-day threshold (mm):",
            value = 0.1, min = 0, max = 1000, step = 0.1
          ),
          actionButton("go", "Update plot", class = "btn btn-primary"),
          helpText("PRCP units: millimeters (GHCND daily).")
        ),
        mainPanel(
          br(),
          plotOutput("ts_plot", height = "380px"),
          br(),
          fluidRow(
            column(2, strong("Days in range:"), textOutput("n_days")),
            column(2, strong("Total PRCP (mm):"), textOutput("sum_mm")),
            column(2, strong("Mean (mm/day):"), textOutput("mean_mm")),
            column(2, strong("Max 1-day (mm):"), textOutput("max_mm")),
            column(2, strong("Wet days:"), textOutput("wet_days")),
            column(2, strong("% days wet:"), textOutput("wet_pct"))
          ),
          br(),
          textOutput("note_na")
        )
      )
    ),

    # -------- Annual Stats & Distributions --------
    tabPanel(
      "Annual Stats & Distributions",
      sidebarLayout(
        sidebarPanel(
          h4("Annual totals"),
          dateRangeInput(
            "annual_range", "Annual stats range:",
            start = date_min, end = date_max,
            min = date_min, max = date_max
          ),
          actionButton("go_annual", "Update annual", class = "btn btn-primary"),
          hr(),
          h4("Distribution comparison (daily PRCP)"),
          dateRangeInput(
            "dist_r1", "Period A:",
            start = date_min, end = date_max, min = date_min, max = date_max
          ),
          checkboxInput("enable_r2", "Overlay Period B", value = FALSE),
          conditionalPanel(
            "input.enable_r2",
            dateRangeInput(
              "dist_r2", "Period B:",
              start = date_min, end = date_max, min = date_min, max = date_max
            )
          ),
          checkboxInput("exclude_zero", "Exclude zero-rain days", value = FALSE),
          actionButton("go_dist", "Update distributions", class = "btn btn-secondary")
        ),
        mainPanel(
          br(),
          h4("Annual Total Rainfall (mm)"),
          plotOutput("annual_bar", height = "380px"),
          br(),
          fluidRow(
            column(4, strong("Years included:"), textOutput("n_years")),
            column(4, strong("Mean annual total (mm):"), textOutput("mean_ann")),
            column(4, strong("SD of annual totals (mm):"), textOutput("sd_ann"))
          ),
          br(),
          fluidRow(
            column(6, strong("Latest 5-year mean (mm):"), textOutput("ma5_latest")),
            column(6, strong("Decadal means (mm):"), textOutput("dec_means"))
          ),
          hr(),
          h4("Distribution of Daily Rainfall (PRCP, mm)"),
          plotlyOutput("dist_plot", height = "360px"),   # <-- interactive
          br(),
          fluidRow(
            column(6, strong("Period A — mean (mm), sd (mm):"), textOutput("stats_r1")),
            column(6, strong("Period B — mean (mm), sd (mm):"), textOutput("stats_r2"))
          )
        )
      )
    ),

    # ---------- Monthly (single-month across years) ----------
    tabPanel(
      "Monthly",
      sidebarLayout(
        sidebarPanel(
          h4("Monthly totals (pick one month)"),
          selectInput(
            "month_sel", "Month:",
            choices = setNames(1:12, month.name), selected = 1
          ),
          dateRangeInput(
            "monthly_range", "Monthly stats range:",
            start = date_min, end = date_max,
            min = date_min, max = date_max
          ),
          actionButton("go_month", "Update monthly", class = "btn btn-primary"),
          hr(),
          h4("Distribution comparison (daily PRCP in chosen month)"),
          dateRangeInput(
            "mdist_r1", "Period A:",
            start = date_min, end = date_max, min = date_min, max = date_max
          ),
          checkboxInput("menable_r2", "Overlay Period B", value = FALSE),
          conditionalPanel(
            "input.menable_r2",
            dateRangeInput(
              "mdist_r2", "Period B:",
              start = date_min, end = date_max, min = date_min, max = date_max
            )
          ),
          checkboxInput("mexclude_zero", "Exclude zero-rain days", value = FALSE),
          actionButton("go_mdist", "Update distributions", class = "btn btn-secondary")
        ),
        mainPanel(
          br(),
          h4(textOutput("monthly_title")),
          plotOutput("monthly_bar", height = "380px"),
          br(),
          fluidRow(
            column(4, strong("Years included:"), textOutput("m_n_years")),
            column(4, strong("Mean monthly total (mm):"), textOutput("m_mean")),
            column(4, strong("SD of monthly totals (mm):"), textOutput("m_sd"))
          ),
          br(),
          fluidRow(
            column(6, strong("Latest 5-year mean (mm):"), textOutput("m_ma5_latest")),
            column(6, strong("Decadal means (mm):"), textOutput("m_dec_means"))
          ),
          hr(),
          h4(textOutput("mdist_title")),
          plotlyOutput("mdist_plot", height = "360px"),  # <-- interactive
          br(),
          fluidRow(
            column(6, strong("Period A — mean (mm), sd (mm):"), textOutput("m_stats_r1")),
            column(6, strong("Period B — mean (mm), sd (mm):"), textOutput("m_stats_r2"))
          )
        )
      )
    )
  )
)

# ---------------------------
# Server
# ---------------------------
server <- function(input, output, session) {

  # ---- Daily tab ----
  filtered <- eventReactive(input$go, {
    req(input$dates)
    df <- data |>
      filter(DATE >= input$dates[1], DATE <= input$dates[2]) |>
      filter(!is.na(PRCP))
    list(df = df, wet_thr = as_num_safely(input$wet_thr, 0))
  }, ignoreInit = FALSE)

  output$ts_plot <- renderPlot({
    f <- filtered(); df <- f$df
    validate(need(nrow(df) > 0, "No data in selected range."))
    ggplot(df, aes(DATE, PRCP)) +
      geom_col(width = 1) +
      labs(
        x = "Date", y = "PRCP (mm)",
        title = "Daily Precipitation (PRCP, mm)",
        subtitle = paste(format(min(df$DATE), "%Y-%m-%d"),
                         "to", format(max(df$DATE), "%Y-%m-%d"))
      ) +
      scale_y_continuous(labels = label_number()) +
      theme_minimal(base_size = 13)
  })

  output$n_days  <- renderText(nrow(filtered()$df))
  output$sum_mm  <- renderText(round(sum(filtered()$df$PRCP, na.rm = TRUE), 3))
  output$mean_mm <- renderText(round(mean(filtered()$df$PRCP, na.rm = TRUE), 3))
  output$max_mm  <- renderText(round(max(filtered()$df$PRCP, na.rm = TRUE), 3))

  output$wet_days <- renderText({
    f <- filtered(); sum(f$df$PRCP >= f$wet_thr, na.rm = TRUE)
  })
  output$wet_pct <- renderText({
    f <- filtered(); n <- nrow(f$df); if (n == 0) return("—")
    paste0(round(100 * sum(f$df$PRCP >= f$wet_thr, na.rm = TRUE) / n, 1), "%")
  })
  output$note_na <- renderText({
    rng <- input$dates; if (is.null(rng)) return("")
    n_all <- data |> filter(DATE >= rng[1], DATE <= rng[2]) |> nrow()
    n_kept <- nrow(filtered()$df)
    n_drop <- n_all - n_kept
    if (n_drop > 0) paste0("Note: ", n_drop, " rows with NA PRCP were omitted.") else ""
  })

  # ---- Annual stats ----
  annual_dat <- eventReactive(input$go_annual, {
    req(input$annual_range)
    data |>
      filter(DATE >= input$annual_range[1], DATE <= input$annual_range[2]) |>
      mutate(Year = year(DATE)) |>
      group_by(Year) |>
      summarise(AnnualTotal = sum(PRCP, na.rm = TRUE), .groups = "drop") |>
      arrange(Year)
  }, ignoreInit = FALSE)

  output$annual_bar <- renderPlot({
    ad <- annual_dat()
    validate(need(nrow(ad) > 0, "No years in selected range."))

    ad <- ad |>
      mutate(
        MA5 = zoo::rollmean(AnnualTotal, 5, align = "right", fill = NA),
        Decade = floor(Year / 10) * 10
      )

    dec <- ad |>
      group_by(Decade) |>
      summarise(
        y = mean(AnnualTotal, na.rm = TRUE),
        x_start = min(Year),
        x_end   = max(Year),
        x_mid   = (x_start + x_end) / 2,
        .groups = "drop"
      )

    ggplot(ad, aes(x = Year, y = AnnualTotal)) +
      geom_col() +
      geom_segment(
        data = dec,
        aes(x = x_start, xend = x_end, y = y, yend = y),
        linewidth = 1
      ) +
      geom_text(
        data = dec,
        aes(x = x_mid, y = y, label = paste0(Decade, "s")),
        vjust = -0.7, size = 3
      ) +
      geom_line(aes(y = MA5), color = "red", linewidth = 1) +
      labs(
        x = "Year",
        y = "Total rainfall (mm)",
        title = "Annual Total Rainfall with 5-year Moving Average (red) and Decadal Means"
      ) +
      scale_x_continuous(breaks = pretty_breaks()) +
      scale_y_continuous(labels = label_number()) +
      theme_minimal(base_size = 13)
  })

  output$n_years <- renderText(nrow(annual_dat()))
  output$mean_ann <- renderText(round(mean(annual_dat()$AnnualTotal, na.rm = TRUE), 1))
  output$sd_ann   <- renderText(round(sd(annual_dat()$AnnualTotal,   na.rm = TRUE), 1))

  output$ma5_latest <- renderText({
    ad <- annual_dat()
    if (nrow(ad) < 5) return("— (need ≥5 years)")
    ad <- ad |> mutate(MA5 = zoo::rollmean(AnnualTotal, 5, align = "right", fill = NA))
    idx <- which(!is.na(ad$MA5))
    if (length(idx) == 0) return("—")
    last_idx <- tail(idx, 1)
    yr_end <- ad$Year[last_idx]
    yr_start <- yr_end - 4
    paste0(round(ad$MA5[last_idx], 1), " (", yr_start, "–", yr_end, ")")
  })

  output$dec_means <- renderText({
    ad <- annual_dat()
    if (nrow(ad) == 0) return("—")
    dec <- ad |>
      mutate(Decade = floor(Year / 10) * 10) |>
      group_by(Decade) |>
      summarise(mu = mean(AnnualTotal, na.rm = TRUE), .groups = "drop") |>
      arrange(Decade)
    paste(
      purrr::map2_chr(dec$Decade, dec$mu, ~ paste0(.x, "s: ", round(.y, 1))),
      collapse = "; "
    )
  })

  # ---- Annual distribution comparison (interactive) ----
  dist_inputs <- eventReactive(input$go_dist, {
    req(input$dist_r1)
    list(
      r1 = input$dist_r1,
      r2 = if (isTRUE(input$enable_r2)) input$dist_r2 else NULL,
      excl0 = isTRUE(input$exclude_zero)
    )
  }, ignoreInit = FALSE)

  output$dist_plot <- renderPlotly({
    di <- dist_inputs()

    make_period <- function(rng, label) {
      df <- data |>
        dplyr::filter(DATE >= rng[1], DATE <= rng[2]) |>
        dplyr::select(PRCP) |>
        dplyr::filter(!is.na(PRCP))
      if (isTRUE(di$excl0)) df <- dplyr::filter(df, PRCP > 0)
      dplyr::mutate(df, Period = label)
    }

    r1 <- make_period(di$r1, "A")
    validate(need(nrow(r1) > 0, "No data in Period A."))

    df_all <- r1
    if (!is.null(di$r2)) {
      r2 <- make_period(di$r2, "B")
      validate(need(nrow(r2) > 0, "No data in Period B."))
      df_all <- dplyr::bind_rows(r1, r2)
    }

    p <- ggplot(df_all, aes(x = PRCP, fill = Period, color = Period)) +
      geom_histogram(aes(y = after_stat(density)),
                     bins = 40, alpha = 0.15, position = "identity", na.rm = TRUE) +
      geom_density(linewidth = 1, alpha = 0.5, na.rm = TRUE) +
      labs(
        x = "Daily PRCP (mm)", y = "Density",
        title = "Distribution of Daily Rainfall",
        subtitle = if (is.null(di$r2)) "Period A only" else "Overlay of Periods A and B"
      ) +
      theme_minimal(base_size = 13)

    ggplotly(p) |> layout(dragmode = "zoom")
  })

  output$stats_r1 <- renderText({
    di <- dist_inputs()
    r1 <- data |>
      dplyr::filter(DATE >= di$r1[1], DATE <= di$r1[2]) |>
      dplyr::select(PRCP) |>
      dplyr::filter(!is.na(PRCP))
    if (isTRUE(di$excl0)) r1 <- dplyr::filter(r1, PRCP > 0)
    if (nrow(r1) < 1) return("—")
    paste(round(mean(r1$PRCP), 3), ", ", round(stats::sd(r1$PRCP), 3))
  })

  output$stats_r2 <- renderText({
    di <- dist_inputs()
    if (is.null(di$r2)) return("—")
    r2 <- data |>
      dplyr::filter(DATE >= di$r2[1], DATE <= di$r2[2]) |>
      dplyr::select(PRCP) |>
      dplyr::filter(!is.na(PRCP))
    if (isTRUE(di$excl0)) r2 <- dplyr::filter(r2, PRCP > 0)
    if (nrow(r2) < 1) return("—")
    paste(round(mean(r2$PRCP), 3), ", ", round(stats::sd(r2$PRCP), 3))
  })

  # ===== Monthly tab logic (single month across years) =====
  output$monthly_title <- renderText({
    paste0("Monthly Total Rainfall — ", month.name[as.integer(input$month_sel)],
           " (bars), with 5-year Moving Average (red) and Decadal Means")
  })

  monthly_dat <- eventReactive(input$go_month, {
    req(input$monthly_range, input$month_sel)
    mnum <- as.integer(input$month_sel)
    data |>
      filter(DATE >= input$monthly_range[1], DATE <= input$monthly_range[2]) |>
      mutate(Year = year(DATE), Mon = month(DATE)) |>
      filter(Mon == mnum) |>
      group_by(Year) |>
      summarise(MonthTotal = sum(PRCP, na.rm = TRUE), .groups = "drop") |>
      arrange(Year)
  }, ignoreInit = FALSE)

  output$monthly_bar <- renderPlot({
    md <- monthly_dat()
    validate(need(nrow(md) > 0, "No years in selected range for this month."))

    md <- md |>
      mutate(
        MA5 = zoo::rollmean(MonthTotal, 5, align = "right", fill = NA),
        Decade = floor(Year / 10) * 10
      )

    decm <- md |>
      group_by(Decade) |>
      summarise(
        y = mean(MonthTotal, na.rm = TRUE),
        x_start = min(Year),
        x_end   = max(Year),
        x_mid   = (x_start + x_end) / 2,
        .groups = "drop"
      )

    ggplot(md, aes(x = Year, y = MonthTotal)) +
      geom_col() +
      geom_segment(data = decm,
                   aes(x = x_start, xend = x_end, y = y, yend = y),
                   linewidth = 1) +
      geom_text(data = decm,
                aes(x = x_mid, y = y, label = paste0(Decade, "s")),
                vjust = -0.7, size = 3) +
      geom_line(aes(y = MA5), color = "red", linewidth = 1) +
      labs(x = "Year", y = "Total rainfall (mm)") +
      scale_x_continuous(breaks = pretty_breaks()) +
      scale_y_continuous(labels = label_number()) +
      theme_minimal(base_size = 13)
  })

  output$m_n_years <- renderText(nrow(monthly_dat()))
  output$m_mean    <- renderText(round(mean(monthly_dat()$MonthTotal, na.rm = TRUE), 1))
  output$m_sd      <- renderText(round(sd(monthly_dat()$MonthTotal,   na.rm = TRUE), 1))

  output$m_ma5_latest <- renderText({
    md <- monthly_dat()
    if (nrow(md) < 5) return("— (need ≥5 years)")
    md <- md |> mutate(MA5 = zoo::rollmean(MonthTotal, 5, align = "right", fill = NA))
    idx <- which(!is.na(md$MA5))
    if (length(idx) == 0) return("—")
    last_idx <- tail(idx, 1)
    yr_end <- md$Year[last_idx]
    yr_start <- yr_end - 4
    paste0(round(md$MA5[last_idx], 1), " (", yr_start, "–", yr_end, ")")
  })

  output$m_dec_means <- renderText({
    md <- monthly_dat()
    if (nrow(md) == 0) return("—")
    dec <- md |>
      mutate(Decade = floor(Year / 10) * 10) |>
      group_by(Decade) |>
      summarise(mu = mean(MonthTotal, na.rm = TRUE), .groups = "drop") |>
      arrange(Decade)
    paste(purrr::map2_chr(dec$Decade, dec$mu, ~ paste0(.x, "s: ", round(.y, 1))),
          collapse = "; ")
  })

  output$mdist_title <- renderText({
    paste0("Distribution of Daily Rainfall (", month.name[as.integer(input$month_sel)], ")")
  })

  mdist_inputs <- eventReactive(input$go_mdist, {
    req(input$mdist_r1, input$month_sel)
    list(
      mon = as.integer(input$month_sel),
      r1 = input$mdist_r1,
      r2 = if (isTRUE(input$menable_r2)) input$mdist_r2 else NULL,
      excl0 = isTRUE(input$mexclude_zero)
    )
  }, ignoreInit = FALSE)

  output$mdist_plot <- renderPlotly({
    di <- mdist_inputs()

    make_period_m <- function(rng, label) {
      df <- data |>
        filter(DATE >= rng[1], DATE <= rng[2]) |>
        mutate(Mon = month(DATE)) |>
        filter(Mon == di$mon) |>
        select(PRCP) |>
        filter(!is.na(PRCP))
      if (isTRUE(di$excl0)) df <- filter(df, PRCP > 0)
      mutate(df, Period = label)
    }

    r1 <- make_period_m(di$r1, "A")
    validate(need(nrow(r1) > 0, "No data in Period A for this month."))

    df_all <- r1
    if (!is.null(di$r2)) {
      r2 <- make_period_m(di$r2, "B")
      validate(need(nrow(r2) > 0, "No data in Period B for this month."))
      df_all <- bind_rows(r1, r2)
    }

    p <- ggplot(df_all, aes(x = PRCP, fill = Period, color = Period)) +
      geom_histogram(aes(y = after_stat(density)),
                     bins = 40, alpha = 0.15, position = "identity", na.rm = TRUE) +
      geom_density(linewidth = 1, alpha = 0.5, na.rm = TRUE) +
      labs(
        x = "Daily PRCP (mm)", y = "Density",
        title = paste0("Distribution of Daily Rainfall — ", month.name[di$mon]),
        subtitle = if (is.null(di$r2)) "Period A only" else "Overlay of Periods A and B"
      ) +
      theme_minimal(base_size = 13)

    ggplotly(p) |> layout(dragmode = "zoom")
  })

  output$m_stats_r1 <- renderText({
    di <- mdist_inputs()
    r1 <- data |>
      filter(DATE >= di$r1[1], DATE <= di$r1[2]) |>
      mutate(Mon = month(DATE)) |>
      filter(Mon == di$mon) |>
      select(PRCP) |>
      filter(!is.na(PRCP))
    if (isTRUE(di$excl0)) r1 <- filter(r1, PRCP > 0)
    if (nrow(r1) < 1) return("—")
    paste(round(mean(r1$PRCP), 3), ", ", round(stats::sd(r1$PRCP), 3))
  })

  output$m_stats_r2 <- renderText({
    di <- mdist_inputs()
    if (is.null(di$r2)) return("—")
    r2 <- data |>
      filter(DATE >= di$r2[1], DATE <= di$r2[2]) |>
      mutate(Mon = month(DATE)) |>
      filter(Mon == di$mon) |>
      select(PRCP) |>
      filter(!is.na(PRCP))
    if (isTRUE(di$excl0)) r2 <- filter(r2, PRCP > 0)
    if (nrow(r2) < 1) return("—")
    paste(round(mean(r2$PRCP), 3), ", ", round(stats::sd(r2$PRCP), 3))
  })
}

shinyApp(ui, server)
