# app.R — Sea-Tac Rainfall Explorer (no sidebar, advanced distribution toggles)
# - Fixed units: mm
# - Annual moving average fixed to 7 years (always shown)
# - Annual distribution tab: toggle multiple fits (Gamma, Normal, Log-Pearson III, Lognormal, Exponential)
#   * Each curve has its own color
#   * Parameters shown only for selected fits
#   * PDF + CDF, with tails shaded below 10% and above 90% for each selected fit
# - Daily plot covers full data range (zoom/pan with Plotly)
# - Simulation tabs unchanged

library(shiny)
library(tidyverse)
library(lubridate)
library(readr)
library(scales)
library(zoo)
library(plotly)
library(MASS)
library(purrr)
library(fitdistrplus)
library(RColorBrewer)
# ---------------------------
# Data load (bullet-proof)
# ---------------------------
DATA_PATH <- "data/seatac_data.csv"
stopifnot(file.exists(DATA_PATH))

raw0 <- readr::read_csv(
  DATA_PATH,
  col_types = readr::cols(.default = readr::col_character()),
  na = c("", "NA", "NaN"),
  show_col_types = FALSE
) |>
  dplyr::rename_with(toupper)

date_parsed <- lubridate::parse_date_time(
  raw0$DATE,
  orders = c("Y-m-d", "m/d/Y", "Y/m/d", "d/m/Y", "Ymd"),
  quiet = TRUE
)
DATE <- as.Date(date_parsed)
PRCP <- readr::parse_number(raw0$PRCP) / 10  # tenths of mm -> mm

data <- tibble::tibble(DATE = DATE, PRCP = PRCP) |>
  dplyr::filter(!is.na(DATE)) |>
  dplyr::arrange(DATE)

if (any(is.na(date_parsed))) {
  bad_n <- sum(is.na(date_parsed))
  message(sprintf("Warning: %d DATE values could not be parsed and were dropped.", bad_n))
  print(utils::head(raw0$DATE[is.na(date_parsed)], 5))
}

stopifnot(nrow(data) > 0)
date_min <- min(data$DATE, na.rm = TRUE)
date_max <- max(data$DATE, na.rm = TRUE)

annual_mm <- data |>
  dplyr::mutate(Year = lubridate::year(DATE)) |>
  dplyr::group_by(Year) |>
  dplyr::summarise(AnnualPRCP_mm = sum(PRCP, na.rm = TRUE), .groups = "drop") |>
  dplyr::filter(!is.na(AnnualPRCP_mm), AnnualPRCP_mm >= 0) |>
  dplyr::arrange(Year)

# ---------------------------
# Helpers (fixed to mm)
# ---------------------------
label_units <- function() "mm"
convert_units <- function(x_mm) x_mm

wet_run_lengths <- function(x_mm) {
  wet <- x_mm > 0
  r <- rle(wet)
  lens <- r$lengths[r$values]
  lens[!is.na(lens)]
}
dry_run_lengths <- function(x_mm) {
  wet <- x_mm > 0
  r <- rle(wet)
  lens <- r$lengths[!r$values]
  lens[!is.na(lens)]
}

get_interstorm_lengths <- function(df) {
  if (nrow(df) == 0) return(numeric(0))
  rle_block <- rle(df$PRCP == 0)
  rle_block$lengths[rle_block$values]
}

get_duration_lengths <- function(df) {
  if (nrow(df) == 0) return(numeric(0))
  rle_block <- rle(df$PRCP > 0)
  rle_block$lengths[rle_block$values]
}


sim_stats <- function(series, units_label) {
  tibble::tibble(
    Metric = c("Total rainfall", "Mean daily", "Wet days (>0)", "Dry days (=0)", "Max daily"),
    Value  = c(
      sprintf("%.1f %s", sum(series, na.rm = TRUE), units_label),
      sprintf("%.2f %s", mean(series, na.rm = TRUE), units_label),
      sum(series > 0, na.rm = TRUE),
      sum(series <= 0, na.rm = TRUE),
      sprintf("%.2f %s", suppressWarnings(max(series, na.rm = TRUE)), units_label)
    )
  )
}

# ---------------------------
# Distribution utilities (annual totals)
# ---------------------------
# Colors (consistent across PDF & CDF)
dist_colors <- c(
  gamma        = "#1b9e77",
  normal       = "#d95f02",
  logpearson3  = "#7570b3",
  lognormal    = "#e7298a", # "another fat-tailed distribution"
  exponential  = "#66a61e"
)

# Fitters return lists with: name, params (named), dens, pfun, qfun, label (for legend)
fit_distributions <- function(x) {
  x_pos <- x[x > 0 & is.finite(x)]
  out <- list()
  # ---------- Gamma(shape, rate) with safe start + rate/scale-agnostic extraction ----------
  if (length(x_pos) > 1) {
    m <- mean(x_pos); v <- var(x_pos)
    shape0 <- if (v > 0) m^2 / v else 1
    rate0  <- if (v > 0) m / v else 1

    fg <- try(suppressWarnings(
      MASS::fitdistr(x_pos, densfun = "gamma", start = list(shape = shape0, rate = rate0))
    ), silent = TRUE)

    if (!inherits(fg, "try-error")) {
      est <- fg$estimate
      shp  <- unname(est[["shape"]])
      rate <- if ("rate" %in% names(est)) {
        unname(est[["rate"]])
      } else if ("scale" %in% names(est)) {
        1 / unname(est[["scale"]])
      } else rate0
      # guard
      if (!is.finite(shp) || shp <= 0 || !is.finite(rate) || rate <= 0) {
        shp  <- max(1e-6, shape0)
        rate <- max(1e-6, rate0)
      }
    } else {
      shp  <- max(1e-6, shape0)
      rate <- max(1e-6, rate0)
    }

    out$gamma <- list(
      name = "gamma",
      params = c(shape = shp, rate = rate,
                 mean = shp / rate, sd = sqrt(shp) / rate),
      dens = function(z) dgamma(z, shape = shp, rate = rate),
      pfun = function(z) pgamma(z, shape = shp, rate = rate),
      qfun = function(p) qgamma(p, shape = shp, rate = rate),
      label = "Gamma"
    )
  }



  # ---------- Normal(mean, sd) ----------
  fn <- try(suppressWarnings(MASS::fitdistr(x, densfun = "normal")), silent = TRUE)
  if (!inherits(fn, "try-error")) {
    mu <- unname(fn$estimate[["mean"]]); sdv <- max(1e-12, unname(fn$estimate[["sd"]]))
    out$normal <- list(
      name = "normal",
      params = c(mean = mu, sd = sdv),
      dens = function(z) dnorm(z, mean = mu, sd = sdv),
      pfun = function(z) pnorm(z, mean = mu, sd = sdv),
      qfun = function(p) qnorm(p, mean = mu, sd = sdv),
      label = "Normal"
    )
  }

  # ---------- Log-Pearson III (Pearson III on Y = ln X) ----------
  # Safe for negative skew: use positive scale, reflect when g < 0
  if (length(x_pos) > 2) {
    y <- log(x_pos)
    m <- mean(y); s <- sd(y)
    if (is.finite(s) && s > 0) {
      g <- mean(((y - m)/s)^3)  # skew of log data
      if (is.finite(g) && abs(g) > 1e-6) {
        k      <- 4 / (g^2)                  # shape > 0
        theta0 <- s * abs(g) / 2             # strictly positive scale
        sg     <- sign(g)                    # +1 or -1
        # location; keep original sign in the shift
        loc    <- m - sg * k * theta0

        dY <- function(yy) {
          if (sg >= 0) {
            ifelse(yy > loc, dgamma(yy - loc, shape = k, scale = theta0), 0)
          } else {
            ifelse(yy < loc, dgamma(loc - yy, shape = k, scale = theta0), 0)
          }
        }
        pY <- function(yy) {
          if (sg >= 0) {
            ifelse(yy > loc, pgamma(yy - loc, shape = k, scale = theta0), 0)
          } else {
            # P(Y <= y) = 1 - P(Gamma < loc - y) for y < loc
            ifelse(yy < loc, 1 - pgamma(loc - yy, shape = k, scale = theta0), 1)
          }
        }

        out$logpearson3 <- list(
          name = "logpearson3",
          params = c(k = k, theta = theta0, loc = loc, skew_log = g),
          dens = function(z) {
            zz <- pmax(z, .Machine$double.eps)
            dY(log(zz)) / zz
          },
          pfun = function(z) {
            zz <- pmax(z, .Machine$double.eps)
            pY(log(zz))
          },
          qfun = function(p) {
            p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
            # bracket in x-space using wide log-range and sort to ensure lower<upper
            lohi <- sort(c(exp(m - 8 * s), exp(m + 8 * s)))
            uniroot(function(xx) out$logpearson3$pfun(xx) - p,
                    lower = lohi[1], upper = lohi[2], tol = 1e-8)$root
          },
          label = "Log-Pearson III"
        )
      }
    }
  }

  # ---------- Lognormal(meanlog, sdlog) ----------
  if (length(x_pos) > 1) {
    ml <- mean(log(x_pos)); sl <- sd(log(x_pos))
    if (is.finite(ml) && is.finite(sl) && sl > 0) {
      out$lognormal <- list(
        name = "lognormal",
        params = c(meanlog = ml, sdlog = sl,
                   mean = exp(ml + 0.5 * sl^2),
                   sd = sqrt((exp(sl^2) - 1) * exp(2 * ml + sl^2))),
        dens = function(z) dlnorm(z, meanlog = ml, sdlog = sl),
        pfun = function(z) plnorm(z, meanlog = ml, sdlog = sl),
        qfun = function(p) qlnorm(p, meanlog = ml, sdlog = sl),
        label = "Lognormal (fat-tailed)"
      )
    }
  }

  # ---------- Exponential(rate) ----------
  if (length(x_pos) > 0) {
    rate <- 1 / mean(x_pos)
    out$exponential <- list(
      name = "exponential",
      params = c(rate = rate, mean = 1 / rate),
      dens = function(z) dexp(z, rate = rate),
      pfun = function(z) pexp(z, rate = rate),
      qfun = function(p) qexp(p, rate = rate),
      label = "Exponential"
    )
  }

  out
}

# ---------------------------
# UI (no sidebar)
# ---------------------------
ui <- fluidPage(
  titlePanel("Sea-Tac Rainfall Explorer"),
  fluidRow(
    column(
      width = 12,
      tabsetPanel(
        id = "parent_tabs",
        type = "tabs",

        # ===== Annual Data =====
        tabPanel(
          title = "Annual Data",
          br(),
          tabsetPanel(
            id = "annual_tabs",
            type = "tabs",

            tabPanel(
              title = "Daily precipitation",
              br(),
              plotlyOutput("daily_plot", height = 420),
              br(),
              fluidRow(
                column(6, wellPanel(h4("Daily summary (full record)"), verbatimTextOutput("daily_summary"))),
                column(6, wellPanel(h4("Notes"),
                                    tags$ul(
                                      tags$li("Drag to zoom, double-click to reset zoom."),
                                      tags$li("Units are fixed to millimeters (mm)."),
                                      tags$li("This view spans the full available record.")
                                    )))
              )
            ),

            tabPanel(
              title = "Annual totals",
              br(),
              plotOutput("annual_plot", height = 420),
              br(),
              fluidRow(
                column(6, wellPanel(h4("Annual statistics"), verbatimTextOutput("annual_stats"))),
                column(6, wellPanel(h4("What to look for"),
                                    tags$ul(
                                      tags$li("Long-term trends or shifts in typical yearly rainfall."),
                                      tags$li("Multi-year wet/dry spells (7-yr moving average shown)."),
                                      tags$li("Outlier years with unusually high or low totals.")
                                    )))
              )
            ),

            tabPanel(
              title = "Distribution of annual totals",
              br(),
              # ---- controls for toggling distributions ----
              fluidRow(
                column(
                  12,
                  checkboxGroupInput(
                    "dist_choices", "Overlay distributions (select one or more):",
                    choices = c(
                      "Gamma" = "gamma",
                      "Normal" = "normal",
                      "Log-Pearson III" = "logpearson3",
                      "Lognormal (fat-tailed)" = "lognormal",
                      "Exponential" = "exponential"
                    ),
                    selected = c("gamma", "normal")  # sensible defaults
                  )
                )
              ),
              # ---- PDF + CDF ----
              fluidRow(
                column(6, plotOutput("annual_pdf", height = 420)),
                column(6, plotlyOutput("annual_cdf", height = 420))
              ),
              br(),
              fluidRow(
                column(12, wellPanel(h4("Fit parameters (selected)"),
                                     htmlOutput("fit_params_html")))
              ),
              tags$hr(),
              div(style="font-size:12px; color:#666;", "Data source: ", code(DATA_PATH))
            )
          )
        ),

        # ===== Daily Data =====
        tabPanel(
          title = "Daily Data",
          br(),
          tabsetPanel(
            id = "daily_tabs",
            type = "tabs",

            tabPanel(
              title = "Volume",
              br(),
              plotOutput("volume_hist", height = 420),
              br(),
              wellPanel(h4("Fit summary (Gamma on wet-day amounts)"), verbatimTextOutput("volume_fit"))
            ),

            tabPanel(
              title = "Duration",
              br(),
              plotOutput("duration_hist", height = 420),
              br(),
              wellPanel(h4("Fit summary (Gamma on wet-spell length)"), verbatimTextOutput("duration_fit"))
            ),

            tabPanel(
              title = "Interstorm",
              br(),
              plotOutput("interstorm_hist", height = 420),
              br(),
              wellPanel(h4("Fit summary (Gamma on dry-spell length)"), verbatimTextOutput("interstorm_fit"))
            )
          )
        ),

        # ===== Simulation =====
        tabPanel(
          title = "Simulation",
          br(),
          tabsetPanel(
            id = "sim_tabs",
            type = "tabs",

            tabPanel(
              title = "Procedure",
              br(),

              fluidRow(
                column(4, plotOutput("sim_plot_interstorm", height = 220)),
                column(4, plotOutput("sim_plot_duration", height = 220)),
                column(4, plotOutput("sim_plot_volume", height = 220))
              ),
              br(),
              fluidRow(
                column(4, actionButton("sim_step", "Sample Step", class = "btn-primary")),
                column(4, actionButton("sim_complete", "Complete 1-Year Simulation", class = "btn-success")),
                column(4, actionButton("sim_reset", "Reset", class = "btn-danger"))
              ),
              br(),
              plotOutput("sim_series", height = 360),
              br(),
              uiOutput("sim_year_summary_ui"),
              br(),
              fluidRow(
                column(4, plotOutput("overlay_volume", height = 220)),
                column(4, plotOutput("overlay_duration", height = 220)),
                column(4, plotOutput("overlay_interstorm", height = 220))
              )
            ),

            tabPanel(
              title = "Bulk sim",
              br(),
              fluidRow(
                column(
                  6,
                  selectInput(
                    "bulk_start_year",
                    "Pick a start year for 3-year historical window",
                    choices = {
                      yrs <- sort(unique(lubridate::year(data$DATE)))
                      yrs[yrs <= max(yrs, na.rm = TRUE) - 2]
                    },
                    selected = {
                      yrs <- sort(unique(lubridate::year(data$DATE)))
                      if (length(yrs) >= 3) yrs[1] else yrs[1]
                    }
                  )
                ),
                column(
                  6,
                  div(style = "margin-top:27px;",
                      actionButton("bulk_run", "Simulate 3 years", class = "btn-primary"),
                      actionButton("bulk_reset", "Clear", class = "btn-warning", style = "margin-left:8px;")
                  )
                )
              ),
              br(),
              fluidRow(
                column(12, plotOutput("bulk_hist_series", height = 260))
              ),
              br(),
              fluidRow(
                column(12, plotOutput("bulk_sim_series", height = 260))
              ),
              br(),
              h4("3-year stats: Historical vs. Simulated"),
              fluidRow(
                column(6, tableOutput("bulk_hist_stats")),
                column(6, tableOutput("bulk_sim_stats"))
              )
            ),
            tabPanel(
              "Monthly Simulation",
              sidebarLayout(
                sidebarPanel(
                  h4("Select months to overlay"),
                  checkboxGroupInput(
                    "monthly_months", NULL,
                    choices = month.name,
                    selected = c("January", "July")
                  )
                ),
                mainPanel(
                  uiOutput("monthly_plot_ui")  # <- NEW (was plotOutput)
                )
              )
            ),
            tabPanel(
              title = "Bulk Month sim",
              br(),
              fluidRow(
                column(
                  6,
                  selectInput(
                    "bulk_month_start_year",
                    "Pick a start year for 3-year historical window",
                    choices = {
                      yrs <- sort(unique(lubridate::year(data$DATE)))
                      yrs[yrs <= max(yrs, na.rm = TRUE) - 2]
                    },
                    selected = {
                      yrs <- sort(unique(lubridate::year(data$DATE)))
                      if (length(yrs) >= 3) yrs[1] else yrs[1]
                    }
                  )
                ),
                column(
                  6,
                  div(style = "margin-top:27px;",
                      actionButton("bulk_month_run", "Simulate 3 years (monthly distributions)", class = "btn-primary"),
                      actionButton("bulk_month_reset", "Clear", class = "btn-warning", style = "margin-left:8px;")
                  )
                )
              ),
              br(),
              fluidRow(
                column(12, plotOutput("bulk_month_hist_series", height = 260))
              ),
              br(),
              fluidRow(
                column(12, plotOutput("bulk_month_sim_series", height = 260))
              ),
              br(),
              h4("3-year stats: Historical vs. Month-based Simulated"),
              fluidRow(
                column(6, tableOutput("bulk_month_hist_stats")),
                column(6, tableOutput("bulk_month_sim_stats"))
              )
            )

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

  # ---- Daily: full record in mm ----
  daily_full <- reactive({
    data |>
      mutate(PRCP_u = convert_units(PRCP))
  })

  output$daily_plot <- renderPlotly({
    d <- daily_full()
    y_lab <- paste0("Daily precipitation (", label_units(), ")")
    p <- ggplot(d, aes(x = DATE, y = PRCP_u)) +
      geom_col(width = 0.9, alpha = 0.85) +
      labs(x = NULL, y = y_lab, title = "Daily precipitation (full record)") +
      scale_y_continuous(labels = label_number(accuracy = 0.1)) +
      theme_minimal(base_size = 13)
    ggplotly(p, dynamicTicks = TRUE) |> layout(hovermode = "x unified")
  })

  output$daily_summary <- renderText({
    d <- daily_full()
    u <- label_units()
    n_days <- nrow(d)
    total <- sum(d$PRCP_u, na.rm = TRUE)
    wetdays <- sum(d$PRCP_u > 0, na.rm = TRUE)
    max_day <- suppressWarnings(max(d$PRCP_u, na.rm = TRUE))
    max_date <- d$DATE[which.max(d$PRCP_u)]
    paste0(
      sprintf("Days: %d\n", n_days),
      sprintf("Total precipitation: %.1f %s\n", total, u),
      sprintf("Wet days (>0): %d\n", wetdays),
      sprintf("Max daily: %.2f %s on %s\n", max_day, u, format(max_date, "%Y-%m-%d"))
    )
  })

  # ---- Annual (fixed MA k = 7, always shown) ----
  annual_reactive <- reactive({
    annual_mm |>
      transmute(Year, Annual_u = convert_units(AnnualPRCP_mm)) |>
      arrange(Year)
  })

  output$annual_plot <- renderPlot({
    a <- annual_reactive()
    u <- label_units()
    a_ma <- a |> mutate(MA = zoo::rollmean(Annual_u, k = 7, fill = NA, align = "center"))
    ggplot(a, aes(x = Year, y = Annual_u)) +
      geom_col(fill = "lightblue", alpha = 0.9) +
      geom_line(data = a_ma, aes(y = MA), linewidth = 1.1, alpha = 0.95) +
      labs(x = NULL, y = paste0("Annual total (", u, ")"),
           title = "Annual precipitation totals (bars) + 7-yr moving average (line)") +
      scale_x_continuous(breaks = pretty_breaks()) +
      scale_y_continuous(labels = label_number(accuracy = 1)) +
      theme_minimal(base_size = 13)
  })

  output$annual_stats <- renderText({
    a <- annual_reactive()
    u <- label_units()
    n <- nrow(a); meanv <- mean(a$Annual_u); sdv <- sd(a$Annual_u)
    minv <- min(a$Annual_u); maxv <- max(a$Annual_u)
    y_min <- a$Year[which.min(a$Annual_u)]; y_max <- a$Year[which.max(a$Annual_u)]
    paste0(
      sprintf("Years: %d\n", n),
      sprintf("Mean annual: %.1f %s\n", meanv, u),
      sprintf("Std. dev.: %.1f %s\n", sdv, u),
      sprintf("Min: %.1f %s (%d)\n", minv, u, y_min),
      sprintf("Max: %.1f %s (%d)\n", maxv, u, y_max)
    )
  })

  # ================= Advanced Annual Distribution Tab =================
  annual_vals <- reactive({
    annual_reactive()$Annual_u
  })

  all_fits <- reactive({
    fit_distributions(annual_vals())
  })

  selected_fits <- reactive({
    req(input$dist_choices)
    fits <- all_fits()
    keep <- intersect(names(fits), input$dist_choices)
    fits[keep]
  })

  # PDF with overlays and shaded tails for each selected distribution
  output$annual_pdf <- renderPlot({
    a <- annual_reactive()
    x <- a$Annual_u
    df_hist <- tibble(x = x)

    # base histogram as density
    p <- ggplot(df_hist, aes(x = x)) +
      geom_histogram(aes(y = after_stat(density)), bins = 25, alpha = 0.9) +
      labs(x = paste0("Annual total (", label_units(), ")"),
           y = "Density", title = "PDF: Distribution of annual totals") +
      theme_minimal(base_size = 13)

    fits <- selected_fits()
    if (length(fits) == 0) return(p)

    # x-range for curves
    xr <- range(x, na.rm = TRUE)
    xs <- seq(from = max(1e-6, xr[1]), to = xr[2], length.out = 400)

    # Add shaded tails + lines per fit
    for (nm in names(fits)) {
      f <- fits[[nm]]
      col <- dist_colors[[nm]]
      q10 <- tryCatch(f$qfun(0.10), error = function(e) NA_real_)
      q90 <- tryCatch(f$qfun(0.90), error = function(e) NA_real_)
      dens_vals <- f$dens(xs)

      # If gamma curve looks numerically flat, draw with moments instead (for display only)
      if (nm == "gamma" && (all(!is.finite(dens_vals)) || max(dens_vals, na.rm = TRUE) < 1e-12)) {
        m <- mean(x, na.rm = TRUE); v <- var(x, na.rm = TRUE)
        shp_m  <- max(1e-6, m^2 / v)
        rate_m <- max(1e-6, m / v)
        dens_vals <- dgamma(xs, shape = shp_m, rate = rate_m)
        q10 <- qgamma(0.10, shape = shp_m, rate = rate_m)
        q90 <- qgamma(0.90, shape = shp_m, rate = rate_m)
      }
      # Tail shading (below 10%)
      if (is.finite(q10)) {
        xs_l <- xs[xs <= q10]
        if (length(xs_l) > 1) {
          p <- p + geom_ribbon(
            data = tibble(xx = xs_l, yy = f$dens(xs_l)),
            aes(x = xx, ymin = 0, ymax = yy),
            fill = col, alpha = 0.12, inherit.aes = FALSE
          )
        }
      }
      # Tail shading (above 90%)
      if (is.finite(q90)) {
        xs_r <- xs[xs >= q90]
        if (length(xs_r) > 1) {
          p <- p + geom_ribbon(
            data = tibble(xx = xs_r, yy = f$dens(xs_r)),
            aes(x = xx, ymin = 0, ymax = yy),
            fill = col, alpha = 0.12, inherit.aes = FALSE
          )
        }
      }

      # Overlay density curve
      p <- p + geom_line(
        data = tibble(xx = xs, yy = dens_vals),
        aes(x = xx, y = yy),
        linewidth = 1.1, color = col, inherit.aes = FALSE
      )
    }

    p
  })

  # CDF with overlays, dotted y=0.1 & 0.9, labeled tails and shaded regions

  output$annual_cdf <- plotly::renderPlotly({
    x <- annual_vals()
    fits <- selected_fits()

    # Empirical CDF
    ec <- ecdf(x)
    xr <- range(x, na.rm = TRUE)
    xs <- seq(from = max(1e-6, xr[1]), to = xr[2], length.out = 400)
    df_ecdf <- tibble(xx = xs, Femp = ec(xs))

    p <- ggplot(df_ecdf, aes(x = xx, y = Femp)) +
      geom_step(linewidth = 1, color = "black", alpha = 0.7) +
      labs(x = paste0("Annual total (", label_units(), ")"),
           y = "Cumulative probability",
           title = "CDF: Annual totals (empirical + selected fits)") +
      theme_minimal(base_size = 13) +
      geom_hline(yintercept = c(0.1, 0.9), linetype = "dotted") +
      annotate("text", x = xr[1], y = 0.05, label = "Lower 10%",
               hjust = 0, vjust = 0.5, size = 3) +
      annotate("text", x = xr[2], y = 0.95, label = "Upper 10%",
               hjust = 1, vjust = 0.5, size = 3)

    if (length(fits) > 0) {
      for (nm in names(fits)) {
        f <- fits[[nm]]
        col <- dist_colors[[nm]]

        # --- Fix: compute Gamma directly to avoid closure issues ---
        if (nm == "gamma") {
          shp  <- as.numeric(f$params[["shape"]])
          rate <- as.numeric(f$params[["rate"]])
          Fvals <- pgamma(xs, shape = shp, rate = rate)
          q10 <- qgamma(0.10, shape = shp, rate = rate)
          q90 <- qgamma(0.90, shape = shp, rate = rate)
        } else {
          Fvals <- f$pfun(xs)
          q10 <- tryCatch(f$qfun(0.10), error = function(e) NA_real_)
          q90 <- tryCatch(f$qfun(0.90), error = function(e) NA_real_)
        }

        # --- Lower tail shading ---
        if (is.finite(q10)) {
          xs_l <- xs[xs <= q10]
          if (length(xs_l) > 1) {
            p <- p + geom_ribbon(
              data = tibble(xx = xs_l,
                            yy = if (nm == "gamma") pgamma(xs_l, shp, rate) else f$pfun(xs_l)),
              aes(x = xx, ymin = 0, ymax = yy),
              fill = col, alpha = 0.10, colour = NA, inherit.aes = FALSE
            )
          }
        }

        # --- Upper tail shading ---
        if (is.finite(q90)) {
          xs_r <- xs[xs >= q90]
          if (length(xs_r) > 1) {
            p <- p + geom_ribbon(
              data = tibble(xx = xs_r,
                            yy = if (nm == "gamma") pgamma(xs_r, shp, rate) else f$pfun(xs_r)),
              aes(x = xx, ymin = yy, ymax = 1),
              fill = col, alpha = 0.10, colour = NA, inherit.aes = FALSE
            )
          }
        }

        # --- Overlay CDF line ---
        p <- p + geom_line(
          data = tibble(xx = xs, yy = Fvals),
          aes(x = xx, y = yy),
          linewidth = 1.1, color = col, inherit.aes = FALSE
        )
      }
    }

    # Make interactive (drag to zoom, double-click to reset)
    plotly::ggplotly(p, dynamicTicks = TRUE) |>
      plotly::layout(hovermode = "x unified", dragmode = "zoom")
  })


  output$fit_params_html <- renderUI({
    fits <- selected_fits()
    if (length(fits) == 0) return(HTML("<em>No fits selected.</em>"))
    # Build colored parameter blocks
    blocks <- purrr::imap_chr(fits, function(f, nm) {
      col <- dist_colors[[nm]]
      par_lines <- paste(
        sprintf("%s = %.4f", names(f$params), as.numeric(f$params)),
        collapse = "<br/>"
      )
      sprintf(
        '<div style="border-left:6px solid %s; padding-left:10px; margin:6px 0;">
           <b>%s</b><br/><span style="font-family:monospace;">%s</span>
         </div>',
        col, f$label, par_lines
      )
    })
    HTML(paste(blocks, collapse = ""))
  })

  # ================= Daily Data distributions & fits (mm) =================
  volume_data <- reactive({
    wd <- data$PRCP
    wd <- wd[wd > 0 & is.finite(wd)]
    validate(need(length(wd) > 1, "Not enough wet days to plot volume distribution."))
    wd
  })
  volume_fit <- reactive({
    wd <- volume_data()
    fit <- try(suppressWarnings(MASS::fitdistr(wd, densfun = "gamma")), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    list(shape = unname(fit$estimate[["shape"]]), rate = unname(fit$estimate[["rate"]]))
  })
  output$volume_hist <- renderPlot({
    wd <- volume_data()
    vf <- volume_fit()
    p <- ggplot(tibble(x = wd), aes(x = x)) +
      geom_histogram(aes(y = after_stat(density)), bins = 40, alpha = 0.9) +
      labs(x = paste0("Wet-day rainfall (", label_units(), ")"), y = "Density",
           title = "Wet-day daily rainfall (Gamma fit)") +
      theme_minimal(base_size = 13)
    if (!is.null(vf)) p <- p + stat_function(fun = function(x) dgamma(x, shape = vf$shape, rate = vf$rate), linewidth = 1.1)
    p
  })
  output$volume_fit <- renderText({
    vf <- volume_fit()
    if (is.null(vf)) return("Gamma fit failed.")
    u <- label_units()
    paste0(
      "Gamma fit on wet-day rainfall:\n",
      sprintf("  shape = %.3f\n", vf$shape),
      sprintf("  rate  = %.6f\n", vf$rate),
      sprintf("  mean  = %.2f %s, SD = %.2f %s", vf$shape / vf$rate, u, sqrt(vf$shape) / vf$rate, u)
    )
  })

  duration_data <- reactive({
    lens <- wet_run_lengths(data$PRCP)
    lens <- lens[is.finite(lens)]
    validate(need(length(lens) > 0, "No wet spells found."))
    lens
  })
  duration_fit <- reactive({
    x <- duration_data()
    if (!all(x > 0)) return(NULL)
    fit <- try(suppressWarnings(MASS::fitdistr(x, densfun = "gamma")), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    list(shape = unname(fit$estimate[["shape"]]), rate = unname(fit$estimate[["rate"]]))
  })
  output$duration_hist <- renderPlot({
    x <- duration_data(); df <- tibble(k = x); f <- duration_fit()
    p <- ggplot(df, aes(x = k)) +
      geom_histogram(aes(y = after_stat(density)), binwidth = 1, boundary = 0.5, closed = "left", alpha = 0.9) +
      scale_x_continuous(breaks = pretty_breaks()) +
      labs(x = "Wet-spell length (days)", y = "Density", title = "Wet-spell length (Gamma fit)") +
      theme_minimal(base_size = 13)
    if (!is.null(f)) p <- p + stat_function(fun = function(x) dgamma(x, shape = f$shape, rate = f$rate), linewidth = 1.1)
    p
  })
  output$duration_fit <- renderText({
    f <- duration_fit()
    if (is.null(f)) return("Gamma fit failed.")
    paste0(
      "Gamma fit on wet-spell length (days):\n",
      sprintf("  shape = %.3f\n", f$shape),
      sprintf("  rate  = %.6f\n", f$rate),
      sprintf("  mean  = %.2f days, SD = %.2f days", f$shape / f$rate, sqrt(f$shape) / f$rate)
    )
  })

  interstorm_data <- reactive({
    lens <- dry_run_lengths(data$PRCP)
    lens <- lens[is.finite(lens)]
    validate(need(length(lens) > 0, "No dry spells found."))
    lens
  })
  interstorm_fit <- reactive({
    x <- interstorm_data()
    if (!all(x > 0)) return(NULL)
    fit <- try(suppressWarnings(MASS::fitdistr(x, densfun = "gamma")), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    list(shape = unname(fit$estimate[["shape"]]), rate = unname(fit$estimate[["rate"]]))
  })
  output$interstorm_hist <- renderPlot({
    x <- interstorm_data(); df <- tibble(k = x); f <- interstorm_fit()
    p <- ggplot(df, aes(x = k)) +
      geom_histogram(aes(y = after_stat(density)), binwidth = 1, boundary = 0.5, closed = "left", alpha = 0.9) +
      scale_x_continuous(breaks = pretty_breaks()) +
      labs(x = "Dry-spell length (days)", y = "Density", title = "Dry-spell length (Gamma fit)") +
      theme_minimal(base_size = 13)
    if (!is.null(f)) p <- p + stat_function(fun = function(x) dgamma(x, shape = f$shape, rate = f$rate), linewidth = 1.1)
    p
  })
  output$interstorm_fit <- renderText({
    f <- interstorm_fit()
    if (is.null(f)) return("Gamma fit failed.")
    paste0(
      "Gamma fit on dry-spell length (days):\n",
      sprintf("  shape = %.3f\n", f$shape),
      sprintf("  rate  = %.6f\n", f$rate),
      sprintf("  mean  = %.2f days, SD = %.2f days", f$shape / f$rate, sqrt(f$shape) / f$rate)
    )
  })
  # ================= Simulation: Procedure (1-year limit + summaries) =================
  MAX_DAYS <- 365

  # consistent colors
  sim_colors <- list(
    interstorm = "#d95f02",  # orange
    duration   = "#1b9e77",  # green
    volume     = "#7570b3"   # purple
  )

  rv <- reactiveValues(
    phase = 1,
    series = numeric(0),
    colors = character(0),      # fill color for each bar
    shade = NULL,               # current shaded range (for interstorm/duration)
    last_inter = NA_real_,
    last_dur = NA_real_,
    last_vols = numeric(0),
    finished = FALSE
  )

  # --- Reset simulation ---
  observeEvent(input$sim_reset, {
    rv$phase <- 1
    rv$series <- numeric(0)
    rv$colors <- character(0)
    rv$shade <- NULL
    rv$last_inter <- NA_real_
    rv$last_dur <- NA_real_
    rv$last_vols <- numeric(0)
    rv$finished <- FALSE
  })

  # --- Complete full 1-year simulation ---
  observeEvent(input$sim_complete, {
    if (rv$finished) return(NULL)

    vf <- volume_fit(); df <- duration_fit(); ift <- interstorm_fit()
    if (is.null(vf) || is.null(df) || is.null(ift)) return(NULL)

    while (!rv$finished) {
      remaining <- MAX_DAYS - length(rv$series)
      if (remaining <= 0) {
        rv$finished <- TRUE
        break
      }

      # --- Phase 1: Interstorm ---
      k <- ceiling(rgamma(1, shape = ift$shape, rate = ift$rate))
      k <- min(max(0, as.integer(k)), remaining)
      rv$series <- c(rv$series, rep(0, k))
      rv$colors <- c(rv$colors, rep("gray70", k))
      if (length(rv$series) >= MAX_DAYS) { rv$finished <- TRUE; break }

      # --- Phase 2: Duration ---
      L <- ceiling(rgamma(1, shape = df$shape, rate = df$rate))
      L <- min(max(1, as.integer(L)), MAX_DAYS - length(rv$series))
      rv$series <- c(rv$series, rep(NA, L))
      rv$colors <- c(rv$colors, rep("gray70", L))
      if (length(rv$series) >= MAX_DAYS) { rv$finished <- TRUE; break }

      # --- Phase 3: Volume ---
      vols <- ceiling(pmax(0, rgamma(L, vf$shape, vf$rate)))
      start <- length(rv$series) - L + 1
      end <- length(rv$series)
      rv$series[start:end] <- vols
      rv$colors[start:end] <- sim_colors$volume

      # Continue loop until year complete
      if (length(rv$series) >= MAX_DAYS) rv$finished <- TRUE
    }

    # Ensure shading cleared and final summary set
    rv$shade <- NULL
    rv$phase <- 1
  })

  # --- Step simulation ---
  observeEvent(input$sim_step, {
    if (rv$finished) return(NULL)
    vf <- volume_fit(); df <- duration_fit(); ift <- interstorm_fit()
    remaining <- MAX_DAYS - length(rv$series)
    if (remaining <= 0) { rv$finished <- TRUE; return(NULL) }

    # --- 1: Interstorm ---
    if (rv$phase == 1) {
      if (is.null(ift)) return(NULL)
      k <- ceiling(rgamma(1, shape = ift$shape, rate = ift$rate))
      k <- min(max(0, as.integer(k)), remaining)
      start <- length(rv$series) + 1
      end <- start + k - 1
      rv$series <- c(rv$series, rep(0, k))
      rv$colors <- c(rv$colors, rep("gray70", k))
      rv$shade <- list(start = start, end = end, col = sim_colors$interstorm)
      rv$last_inter <- k
      rv$phase <- 2
      if (length(rv$series) >= MAX_DAYS) rv$finished <- TRUE

      # --- 2: Duration ---
    } else if (rv$phase == 2) {
      if (is.null(df)) return(NULL)
      L <- ceiling(rgamma(1, shape = df$shape, rate = df$rate))
      L <- min(max(1, as.integer(L)), remaining)
      start <- length(rv$series) + 1
      end <- start + L - 1
      rv$series <- c(rv$series, rep(NA, L))     # placeholder for rainfall
      rv$colors <- c(rv$colors, rep("gray70", L))
      rv$shade <- list(start = start, end = end, col = sim_colors$duration)
      rv$last_dur <- L
      rv$phase <- 3

      # --- 3: Volume ---
    } else if (rv$phase == 3) {
      if (is.null(vf) || is.na(rv$last_dur)) return(NULL)
      L <- min(rv$last_dur, remaining)
      vols <- ceiling(pmax(0, rgamma(L, vf$shape, vf$rate)))
      start <- length(rv$series) - L + 1
      end <- length(rv$series)
      rv$series[start:end] <- vols
      # turn previous purple bars to grey, new volumes purple
      rv$colors[] <- "gray70"
      rv$colors[start:end] <- sim_colors$volume
      rv$shade <- NULL
      rv$last_vols <- vols
      rv$phase <- 1
      if (length(rv$series) >= MAX_DAYS) rv$finished <- TRUE
    }
  })

  # --- Status text ---
  output$sim_status <- renderText({
    rem <- MAX_DAYS - length(rv$series)
    paste0(
      "Phase: ",
      c(`1` = "Interstorm → sample dry days",
        `2` = "Duration → sample wet-spell length",
        `3` = "Volume → sample daily rainfall")[as.character(rv$phase)],
      "\nDays simulated: ", length(rv$series), "/", MAX_DAYS,
      if (rv$finished) "\nStatus: Year complete." else paste0("\nRemaining: ", rem, " day(s)."),
      "\nLast interstorm: ", ifelse(is.na(rv$last_inter), "—", rv$last_inter),
      " day(s)\nLast duration: ", ifelse(is.na(rv$last_dur), "—", rv$last_dur),
      " day(s)\nLast volumes: ",
      ifelse(length(rv$last_vols) == 0, "—", paste(rv$last_vols, collapse = ", "))
    )
  })

  # --- Distribution plots: Interstorm → Duration → Volume ---
  # ---- Interstorm (green) ----
  output$sim_plot_interstorm <- renderPlot({
    x <- interstorm_data(); f <- interstorm_fit()
    p <- ggplot(tibble(k = x), aes(x = k)) +
      geom_histogram(aes(y = after_stat(density)), binwidth = 1,
                     boundary = 0.5, fill = sim_colors$interstorm, alpha = 0.6) +
      {if (!is.null(f)) stat_function(fun = function(z) dgamma(z, f$shape, f$rate),
                                      linewidth = 1.1, color = sim_colors$interstorm)} +
      {if (!is.na(rv$last_inter))
        geom_vline(xintercept = rv$last_inter, color = sim_colors$interstorm, linewidth = 1)} +
      labs(x = "Dry-spell (days)", y = "Density", title = "Interstorm (Gamma)") +
      theme_minimal(base_size = 12)

    # Overlay sampled runs once year is complete
    if (rv$finished) {
      samp_len <- sampled_runs()$dry
      if (length(samp_len) > 1)
        p <- p + geom_density(data = tibble(k = samp_len), aes(x = k),
                              color = "red", linewidth = 1.1)
    }

    p
  })

  # ---- Duration (orange) ----
  output$sim_plot_duration <- renderPlot({
    x <- duration_data(); f <- duration_fit()
    p <- ggplot(tibble(k = x), aes(x = k)) +
      geom_histogram(aes(y = after_stat(density)), binwidth = 1,
                     boundary = 0.5, fill = sim_colors$duration, alpha = 0.6) +
      {if (!is.null(f)) stat_function(fun = function(z) dgamma(z, f$shape, f$rate),
                                      linewidth = 1.1, color = sim_colors$duration)} +
      {if (!is.na(rv$last_dur))
        geom_vline(xintercept = rv$last_dur, color = sim_colors$duration, linewidth = 1)} +
      labs(x = "Wet-spell (days)", y = "Density", title = "Duration (Gamma)") +
      theme_minimal(base_size = 12)

    # Overlay sampled runs once year is complete
    if (rv$finished) {
      samp_len <- sampled_runs()$wet
      if (length(samp_len) > 1)
        p <- p + geom_density(data = tibble(k = samp_len), aes(x = k),
                              color = "red", linewidth = 1.1)
    }

    p
  })

  # ---- Volume (purple) ----
  output$sim_plot_volume <- renderPlot({
    wd <- volume_data(); f <- volume_fit()
    p <- ggplot(tibble(x = wd), aes(x = x)) +
      geom_histogram(aes(y = after_stat(density)), bins = 40,
                     fill = sim_colors$volume, alpha = 0.6) +
      {if (!is.null(f)) stat_function(fun = function(z) dgamma(z, f$shape, f$rate),
                                      linewidth = 1.1, color = sim_colors$volume)} +
      {if (length(rv$last_vols) > 0)
        geom_vline(xintercept = rv$last_vols, color = sim_colors$volume, alpha = 0.5)} +
      labs(x = paste0("Wet-day rainfall (", label_units(), ")"),
           y = "Density", title = "Volume (Gamma)") +
      theme_minimal(base_size = 12)

    # Overlay sampled volume once year is complete
    if (rv$finished) {
      samp_wd <- sampled_volume()
      if (length(samp_wd) > 1)
        p <- p + geom_density(data = tibble(x = samp_wd), aes(x = x),
                              color = "red", linewidth = 1.1)
    }

    p
  })


  # --- Time series plot ---
  output$sim_series <- renderPlot({
    s <- rv$series
    if (length(s) == 0) return(NULL)
    df <- tibble(t = seq_along(s), y = s, col = rv$colors)

    p <- ggplot(df, aes(x = t, y = y)) +
      geom_col(aes(fill = col), color = NA, show.legend = FALSE) +
      scale_fill_identity() +
      labs(x = "Time (days)", y = paste0("Amount (", label_units(), ")"),
           title = "Simulated rainfall (colored by current phase)") +
      theme_minimal(base_size = 13)

    # shaded active window (interstorm or duration)
    if (!is.null(rv$shade)) {
      p <- p + annotate("rect",
                        xmin = rv$shade$start - 0.5,
                        xmax = rv$shade$end + 0.5,
                        ymin = -Inf, ymax = Inf,
                        fill = rv$shade$col, alpha = 0.2)
    }

    p
  })


  # ======= Year summary + overlays (shown after year complete) =======
  output$sim_year_summary_ui <- renderUI({
    if (!rv$finished) return(NULL)
    u <- label_units()
    total_sim <- sum(rv$series, na.rm = TRUE)
    mean_hist_annual_u <- convert_units(mean(annual_mm$AnnualPRCP_mm, na.rm = TRUE))
    zeros <- sum(rv$series <= 0, na.rm = TRUE)
    tagList(
      wellPanel(
        h4("One-year simulation summary"),
        HTML(sprintf(
          "<b>Total simulated rainfall:</b> %.1f %s<br/>
           <b>Historical mean annual rainfall:</b> %.1f %s<br/>
           <b>Difference:</b> %.1f %s (sim − hist mean)<br/>
           <b>Dry days (zeros):</b> %d / %d",
          total_sim, u, mean_hist_annual_u, u, total_sim - mean_hist_annual_u, u, zeros, 365
        ))
      ),
      h4("Sampled vs. Historical Distributions")
    )
  })

  sampled_volume <- reactive({
    if (!rv$finished) return(numeric(0))
    x <- rv$series
    x[x > 0]
  })
  sampled_runs <- reactive({
    if (!rv$finished) return(list(wet = integer(0), dry = integer(0)))
    r <- rle(rv$series > 0)
    list(
      wet = r$lengths[r$values],
      dry = r$lengths[!r$values]
    )
  })

  output$overlay_volume <- renderPlot({
    if (!rv$finished) return(NULL)
    hist_wd <- volume_data()
    samp_wd <- sampled_volume()
    p <- ggplot(tibble(x = hist_wd), aes(x = x)) +
      geom_histogram(aes(y = after_stat(density)), bins = 40, alpha = 0.6) +
      labs(x = paste0("Wet-day rainfall (", label_units(), ")"),
           y = "Density", title = "Volume: Historical (hist) vs. Sampled (line)") +
      theme_minimal(base_size = 12)
    if (length(samp_wd) > 1) {
      p <- p + geom_density(data = tibble(x = samp_wd), aes(x = x), linewidth = 1.2)
    }
    p
  })
  output$overlay_duration <- renderPlot({
    if (!rv$finished) return(NULL)
    hist_len <- duration_data()
    samp_len <- sampled_runs()$wet
    p <- ggplot(tibble(k = hist_len), aes(x = k)) +
      geom_histogram(aes(y = after_stat(density)), binwidth = 1, boundary = 0.5,
                     closed = "left", alpha = 0.6) +
      scale_x_continuous(breaks = pretty_breaks()) +
      labs(x = "Wet-spell length (days)", y = "Density",
           title = "Duration: Historical (hist) vs. Sampled (line)") +
      theme_minimal(base_size = 12)
    if (length(samp_len) > 1) {
      p <- p + geom_density(data = tibble(k = samp_len), aes(x = k), linewidth = 1.2)
    }
    p
  })
  output$overlay_interstorm <- renderPlot({
    if (!rv$finished) return(NULL)
    hist_len <- interstorm_data()
    samp_len <- sampled_runs()$dry
    p <- ggplot(tibble(k = hist_len), aes(x = k)) +
      geom_histogram(aes(y = after_stat(density)), binwidth = 1, boundary = 0.5,
                     closed = "left", alpha = 0.6) +
      scale_x_continuous(breaks = pretty_breaks()) +
      labs(x = "Dry-spell length (days)", y = "Density",
           title = "Interstorm: Historical (hist) vs. Sampled (line)") +
      theme_minimal(base_size = 12)
    if (length(samp_len) > 1) {
      p <- p + geom_density(data = tibble(k = samp_len), aes(x = k), linewidth = 1.2)
    }
    p
  })

  # ===================== BULK SIM (3-year historical vs simulated) =====================
  MAX_DAYS_BULK <- 3 * 365

  bulk_hist_series <- reactive({
    req(input$bulk_start_year)
    start_year <- as.integer(input$bulk_start_year)
    end_year <- start_year + 2
    d <- data |>
      mutate(Year = year(DATE), PRCP_u = convert_units(PRCP)) |>
      filter(Year >= start_year, Year <= end_year) |>
      arrange(DATE)
    vals <- d$PRCP_u
    if (length(vals) >= MAX_DAYS_BULK) vals <- vals[1:MAX_DAYS_BULK]
    vals
  })

  bulk_sim <- reactiveVal(NULL)

  observeEvent(input$bulk_reset, {
    bulk_sim(NULL)
  })

  observeEvent(input$bulk_run, {
    vf <- volume_fit(); df <- duration_fit(); ift <- interstorm_fit()
    if (is.null(vf) || is.null(df) || is.null(ift)) {
      bulk_sim(NULL); return(NULL)
    }
    s <- numeric(0)
    phase <- 1
    L <- 1
    while (length(s) < MAX_DAYS_BULK) {
      remaining <- MAX_DAYS_BULK - length(s)
      if (phase == 1) {
        k <- ceiling(rgamma(1, shape = ift$shape, rate = ift$rate))
        k <- max(0, as.integer(k))
        add_k <- min(k, remaining)
        s <- c(s, rep(0, add_k))
        phase <- 2
      } else if (phase == 2) {
        L <- ceiling(rgamma(1, shape = df$shape, rate = df$rate))
        L <- max(1, as.integer(L))
        L <- min(L, remaining)
        phase <- 3
      } else if (phase == 3) {
        L_curr <- min(L, remaining)
        vols <- rgamma(L_curr, shape = vf$shape, rate = vf$rate)
        vols <- ceiling(pmax(0, vols))
        s <- c(s, vols)
        phase <- 1
      }
    }
    if (length(s) > MAX_DAYS_BULK) s <- s[1:MAX_DAYS_BULK]
    bulk_sim(s)
  })

  output$bulk_hist_series <- renderPlot({
    hist_vals <- bulk_hist_series()
    n <- length(hist_vals)
    df <- tibble(t = seq_len(n), y = hist_vals)
    ggplot(df, aes(x = t, y = y)) +
      geom_col(alpha = 0.9) +
      labs(
        title = "Historical rainfall (3 years)",
        x = "Time (days)", y = paste0("Amount (", label_units(), ")")
      ) +
      theme_minimal(base_size = 13)
  })

  output$bulk_sim_series <- renderPlot({
    s <- bulk_sim()
    validate(need(!is.null(s), "Click 'Simulate 3 years' to generate a simulated series."))
    df <- tibble(t = seq_along(s), y = s)
    ggplot(df, aes(x = t, y = y)) +
      geom_col(alpha = 0.9) +
      labs(
        title = "Simulated rainfall (3 years)",
        x = "Time (days)", y = paste0("Amount (", label_units(), ")")
      ) +
      theme_minimal(base_size = 13)
  })

  output$bulk_hist_stats <- renderTable({
    vals <- bulk_hist_series()
    u <- label_units()
    sim_stats(vals, u)
  }, striped = TRUE, bordered = TRUE, spacing = "s", align = "l")

  output$bulk_sim_stats <- renderTable({
    vals <- bulk_sim()
    validate(need(!is.null(vals), "—"))
    u <- label_units()
    sim_stats(vals, u)
  }, striped = TRUE, bordered = TRUE, spacing = "s", align = "l")

  # ================== Monthly Simulation Tab ==================
  output$monthly_plot_ui <- renderUI({
    req(input$monthly_months)

    month_colors <- setNames(RColorBrewer::brewer.pal(12, "Paired"), month.name)
    df <- data |>
      filter(!is.na(DATE), !is.na(PRCP)) |>
      mutate(month = factor(month(DATE, label = TRUE, abbr = FALSE), levels = month.name))
    selected_months <- intersect(input$monthly_months, month.name)
    validate(need(length(selected_months) > 0, "Select at least one month."))

    get_var_data <- function(var, m) {
      d <- df |> filter(month == m)
      if (nrow(d) == 0) return(numeric(0))
      if (var == "volume") {
        x <- d |> filter(PRCP > 0) |> pull(PRCP)
      } else if (var == "interstorm") {
        x <- get_interstorm_lengths(d)
      } else if (var == "duration") {
        x <- get_duration_lengths(d)
      } else numeric(0)
    }

    make_plotly <- function(var_name, xlab) {
      all_data <- unlist(lapply(month.name, \(m) get_var_data(var_name, m)))
      if (length(all_data) < 3) return(NULL)
      fit_all <- suppressWarnings(fitdistrplus::fitdist(all_data, "gamma"))
      shape <- fit_all$estimate[1]; rate <- fit_all$estimate[2]

      p <- ggplot(tibble(x = all_data), aes(x)) +
        geom_histogram(aes(y = after_stat(density)), bins = 40,
                       fill = "grey80", color = "white") +
        stat_function(fun = \(z) dgamma(z, shape, rate),
                      color = "black", linewidth = 1.1) +
        theme_minimal(base_size = 13) +
        labs(x = xlab, y = "Density",
             title = paste0(str_to_title(var_name), ": Overall (gray) + Monthly overlays"))

      overlay_df <- map_dfr(selected_months, function(m) {
        x <- get_var_data(var_name, m)
        if (length(x) > 3) {
          fit_m <- suppressWarnings(fitdistrplus::fitdist(x, "gamma"))
          tibble(
            z = seq(min(x), max(x), length.out = 200),
            density = dgamma(seq(min(x), max(x), length.out = 200),
                             fit_m$estimate[1], fit_m$estimate[2]),
            Month = m
          )
        }
      })

      if (nrow(overlay_df) > 0) {
        p <- p +
          geom_line(data = overlay_df,
                    aes(x = z, y = density, color = Month),
                    linewidth = 1.1, alpha = 0.9) +
          scale_color_manual(values = month_colors[selected_months]) +
          guides(color = guide_legend(title = "Month"))
      }

      ggplotly(p, dynamicTicks = TRUE) |>
        layout(hovermode = "x unified", dragmode = "zoom")
    }

    fluidRow(
      column(4, make_plotly("volume", "Wet-day rainfall (mm)")),
      column(4, make_plotly("duration", "Wet-spell length (days)")),
      column(4, make_plotly("interstorm", "Dry-spell length (days)"))
    )
  })

  # ===================== BULK MONTH SIM (3-year historical vs seasonal simulated) =====================
  MAX_DAYS_BULK_MONTH <- 3 * 365

  # --- Monthly parameter fits ---
  monthly_gamma_fits <- reactive({
    df <- data |>
      filter(!is.na(DATE), !is.na(PRCP)) |>
      mutate(month = month(DATE, label = TRUE, abbr = FALSE))

    fits <- list(volume = list(), duration = list(), interstorm = list())
    for (m in month.name) {
      d_m <- df |> filter(month == m)
      if (nrow(d_m) < 5) next

      vols <- d_m |> filter(PRCP > 0) |> pull(PRCP)
      dur  <- get_duration_lengths(d_m)
      inter <- get_interstorm_lengths(d_m)

      safe_fit <- function(x) {
        if (length(x) < 3) return(NULL)
        f <- try(suppressWarnings(fitdistrplus::fitdist(x, "gamma")), silent = TRUE)
        if (inherits(f, "try-error")) return(NULL)
        list(shape = f$estimate[1], rate = f$estimate[2])
      }

      fits$volume[[m]] <- safe_fit(vols)
      fits$duration[[m]] <- safe_fit(dur)
      fits$interstorm[[m]] <- safe_fit(inter)
    }
    fits
  })

  bulk_month_hist_series <- reactive({
    req(input$bulk_month_start_year)
    start_year <- as.integer(input$bulk_month_start_year)
    end_year <- start_year + 2
    d <- data |>
      mutate(Year = year(DATE), PRCP_u = convert_units(PRCP)) |>
      filter(Year >= start_year, Year <= end_year) |>
      arrange(DATE)
    vals <- d$PRCP_u
    if (length(vals) >= MAX_DAYS_BULK_MONTH) vals <- vals[1:MAX_DAYS_BULK_MONTH]
    vals
  })

  bulk_month_sim <- reactiveVal(NULL)

  observeEvent(input$bulk_month_reset, {
    bulk_month_sim(NULL)
  })

  observeEvent(input$bulk_month_run, {
    vf_all <- volume_fit(); df_all <- duration_fit(); ift_all <- interstorm_fit()
    monthly_fits <- monthly_gamma_fits()
    if (is.null(vf_all) || is.null(df_all) || is.null(ift_all)) {
      bulk_month_sim(NULL)
      return(NULL)
    }

    s <- numeric(0)
    phase <- 1
    L <- 1

    # Start from January
    day_counter <- 1
    current_month <- 1

    while (length(s) < MAX_DAYS_BULK_MONTH) {
      remaining <- MAX_DAYS_BULK_MONTH - length(s)
      # Determine month from day number
      current_month <- ((day_counter - 1) %/% 30) %% 12 + 1
      mname <- month.name[current_month]

      # Choose month-specific fits (fallback to global if missing)
      vf <- monthly_fits$volume[[mname]] %||% vf_all
      df <- monthly_fits$duration[[mname]] %||% df_all
      ift <- monthly_fits$interstorm[[mname]] %||% ift_all

      if (phase == 1) {
        k <- ceiling(rgamma(1, shape = ift$shape, rate = ift$rate))
        k <- max(0, as.integer(k))
        add_k <- min(k, remaining)
        s <- c(s, rep(0, add_k))
        phase <- 2
        day_counter <- day_counter + add_k
      } else if (phase == 2) {
        L <- ceiling(rgamma(1, shape = df$shape, rate = df$rate))
        L <- max(1, as.integer(L))
        L <- min(L, remaining)
        phase <- 3
      } else if (phase == 3) {
        L_curr <- min(L, remaining)
        vols <- rgamma(L_curr, shape = vf$shape, rate = vf$rate)
        vols <- ceiling(pmax(0, vols))
        s <- c(s, vols)
        phase <- 1
        day_counter <- day_counter + L_curr
      }
    }

    if (length(s) > MAX_DAYS_BULK_MONTH) s <- s[1:MAX_DAYS_BULK_MONTH]
    bulk_month_sim(s)
  })

  output$bulk_month_hist_series <- renderPlot({
    hist_vals <- bulk_month_hist_series()
    n <- length(hist_vals)
    df <- tibble(t = seq_len(n), y = hist_vals)
    ggplot(df, aes(x = t, y = y)) +
      geom_col(alpha = 0.9) +
      labs(
        title = "Historical rainfall (3 years)",
        x = "Time (days)", y = paste0("Amount (", label_units(), ")")
      ) +
      theme_minimal(base_size = 13)
  })

  output$bulk_month_sim_series <- renderPlot({
    s <- bulk_month_sim()
    validate(need(!is.null(s), "Click 'Simulate 3 years (monthly distributions)' to run."))
    df <- tibble(t = seq_along(s), y = s)
    ggplot(df, aes(x = t, y = y)) +
      geom_col(alpha = 0.9, fill = "#377eb8") +
      labs(
        title = "Simulated rainfall (3 years, month-varying distributions)",
        x = "Time (days)", y = paste0("Amount (", label_units(), ")")
      ) +
      theme_minimal(base_size = 13)
  })

  output$bulk_month_hist_stats <- renderTable({
    vals <- bulk_month_hist_series()
    u <- label_units()
    sim_stats(vals, u)
  }, striped = TRUE, bordered = TRUE, spacing = "s", align = "l")

  output$bulk_month_sim_stats <- renderTable({
    vals <- bulk_month_sim()
    validate(need(!is.null(vals), "—"))
    u <- label_units()
    sim_stats(vals, u)
  }, striped = TRUE, bordered = TRUE, spacing = "s", align = "l")

}

shinyApp(ui, server)
