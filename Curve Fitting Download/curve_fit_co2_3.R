# app.R — Fit-by-eye linear model + Range selection + Best Fit button (3 s.f.)
# Expects data at Models/data/co2_data.txt relative to this app.

library(shiny)
library(ggplot2)
library(dplyr)

# ---------- Helpers for 3 significant figures ----------
sf3  <- function(x) signif(x, 3)
fmt3 <- function(x) format(sf3(x), trim = TRUE, scientific = FALSE)
label_sf3 <- function(x) fmt3(x)
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))

# ---------- Load data ----------
data_path <- file.path("Models", "data", "co2_data.txt")
if (!file.exists(data_path)) {
  stop(paste0("Data file not found at: ", normalizePath(data_path, winslash = "/"),
              "\nIf app.R is inside 'Models/', change path to file.path('data','co2_data.txt')."))
}

co2 <- read.table(
  data_path,
  header = FALSE,
  comment.char = "#",
  col.names = c("Year","Month","DecimalDate","CO2_Monthly","CO2_Deseason",
                "NumDays","StdDev","Uncertainty")
) %>%
  as_tibble() %>%
  filter(!is.na(DecimalDate), !is.na(CO2_Monthly)) %>%
  arrange(DecimalDate)

# Center time for stability
t0   <- min(co2$DecimalDate)
tau  <- co2$DecimalDate - t0
y    <- co2$CO2_Monthly
tmin <- min(co2$DecimalDate)
tmax <- max(co2$DecimalDate)

# OLS (for slider ranges)
ols    <- lm(y ~ tau)
a_seed <- unname(coef(ols)[1])
b_seed <- unname(coef(ols)[2])
a_min  <- a_seed - 30
a_max  <- a_seed + 30
b_min  <- b_seed - 0.5
b_max  <- b_seed + 0.5

# ---------- UI ----------
ui <- fluidPage(
  titlePanel("Fit by Eye: Linear Trend for Mauna Loa CO\u2082"),
  sidebarLayout(
    sidebarPanel(
      p("Model:  y = a + b · (t − t₀)"),
      p(em(paste0("t₀ = ", fmt3(t0), " (fixed)"))),
      sliderInput("a", "a (intercept @ t₀, ppm)",
                  min = a_min, max = a_max, value = a_min, step = 0.1),
      sliderInput("b", "b (slope, ppm per year)",
                  min = b_min, max = b_max, value = b_min, step = 0.005),

      # Range selector
      sliderInput("trange", "Data range (years)",
                  min = tmin, max = tmax, value = c(tmin, tmax), step = 0.1),

      # Buttons
      actionButton("reset", "Reset sliders"),
      actionButton("bestfit", "Best Fit (Lowest RMSE)")
    ),
    mainPanel(
      plotOutput("p_fit", height = "460px"),
      br(),
      strong(textOutput("rmse_txt")),
      hr(),
      plotOutput("p_resid", height = "240px")
    )
  )
)

# ---------- Server ----------
server <- function(input, output, session) {
  observeEvent(input$reset, {
    updateSliderInput(session, "a", value = a_min)
    updateSliderInput(session, "b", value = b_min)
  })

  # ---- Auto-fit button (Best Fit) ----
  observeEvent(input$bestfit, {
    # Subset data to selected range
    rmin <- input$trange[1]
    rmax <- input$trange[2]
    sel <- co2$DecimalDate >= rmin & co2$DecimalDate <= rmax
    df_sel <- co2[sel, , drop = FALSE]

    # Fit y = a + b * (t - t0)
    fit <- lm(CO2_Monthly ~ I(DecimalDate - t0), data = df_sel)
    a_fit <- unname(coef(fit)[1])
    b_fit <- unname(coef(fit)[2])

    # Update sliders to fitted values
    updateSliderInput(session, "a", value = a_fit)
    updateSliderInput(session, "b", value = b_fit)
  })

  # Reactive fitted data
  fit_data <- reactive({
    a <- sf3(input$a)
    b <- sf3(input$b)
    yhat_all <- a + b * tau
    df_all <- tibble(
      DecimalDate = co2$DecimalDate,
      Observed = y,
      Fitted = yhat_all
    )

    rmin <- input$trange[1]
    rmax <- input$trange[2]
    sel  <- df_all$DecimalDate >= rmin & df_all$DecimalDate <= rmax
    df_sel <- df_all[sel, , drop = FALSE]
    df_sel <- mutate(df_sel, Residual = Observed - Fitted)

    list(
      a = a, b = b,
      rmin = rmin, rmax = rmax,
      all = df_all,
      sel = df_sel
    )
  })

  # ---- Main fit plot ----
  output$p_fit <- renderPlot({
    fd <- fit_data()

    ggplot() +
      geom_line(data = fd$all, aes(DecimalDate, Observed), alpha = 0.35) +
      geom_line(data = fd$sel, aes(DecimalDate, Observed)) +
      geom_line(data = fd$sel, aes(DecimalDate, Fitted), linewidth = 1) +
      scale_x_continuous(labels = label_sf3) +
      scale_y_continuous(labels = label_sf3) +
      labs(
        title = "Fit by eye or snap to best fit",
        subtitle = paste(
          "t₀ =", fmt3(t0),
          "| a =", fmt3(fd$a),
          "| b =", fmt3(fd$b), "ppm/yr",
          "| Range:", fmt3(fd$rmin), "–", fmt3(fd$rmax)
        ),
        x = "Year",
        y = expression("CO"[2]*" (ppm)")
      ) +
      theme_minimal(base_size = 13)
  })

  # ---- RMSE display ----
  output$rmse_txt <- renderText({
    fd <- fit_data()
    paste0("RMSE on selected range: ",
           fmt3(rmse(fd$sel$Observed, fd$sel$Fitted)), " ppm")
  })

  # ---- Residuals plot ----
  output$p_resid <- renderPlot({
    fd <- fit_data()
    ggplot(fd$sel, aes(DecimalDate, Residual)) +
      geom_hline(yintercept = 0) +
      geom_line() +
      scale_x_continuous(labels = label_sf3) +
      scale_y_continuous(labels = label_sf3) +
      labs(
        title = "Residuals (Observed − Fitted) — Selected Range",
        x = "Year", y = "ppm"
      ) +
      theme_minimal(base_size = 13)
  })
}

shinyApp(ui, server)
