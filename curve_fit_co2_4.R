# app.R — Exponential fit with robust "Best Fit" button (3 significant figures)
# Model: y = A * exp(r * (t - t0)) + B

library(shiny)
library(ggplot2)
library(dplyr)

# ---------- Helpers ----------
sf3  <- function(x) signif(x, 3)
fmt3 <- function(x) format(sf3(x), trim = TRUE, scientific = FALSE)
label_sf3 <- function(x) fmt3(x)
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))

# ---------- Load data ----------
data_path <- file.path("Models", "data", "co2_data.txt")
if (!file.exists(data_path)) {
  stop(paste0("Data file not found at: ", normalizePath(data_path, winslash = "/")))
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

# ---------- Setup ----------
t0   <- min(co2$DecimalDate)
tau  <- co2$DecimalDate - t0
y    <- co2$CO2_Monthly
tmin <- min(co2$DecimalDate)
tmax <- max(co2$DecimalDate)

# Slider ranges
A_guess <- 10
r_guess <- 0.01
B_guess <- min(y)
A_min <- 0
A_max <- 200
r_min <- 0
r_max <- 0.05
B_min <- B_guess - 50
B_max <- B_guess + 100

# ---------- UI ----------
ui <- fluidPage(
  titlePanel("Fit by Eye: Exponential Model for Mauna Loa CO\u2082"),
  sidebarLayout(
    sidebarPanel(
      p("Model:  y = A · exp(r · (t − t₀)) + B"),
      p(em(paste0("t₀ = ", fmt3(t0), " (fixed)"))),

      sliderInput("A", "A (scale)", min = A_min, max = A_max, value = A_guess, step = 1),
      sliderInput("r", "r (rate, 1/year)", min = r_min, max = r_max, value = r_guess, step = 0.001),
      sliderInput("B", "B (baseline, ppm)", min = B_min, max = B_max, value = B_guess, step = 1),

      sliderInput("trange", "Data range (years)",
                  min = tmin, max = tmax, value = c(tmin, tmax), step = 0.1),

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

# ---------- SERVER ----------
server <- function(input, output, session) {
  observeEvent(input$reset, {
    updateSliderInput(session, "A", value = A_guess)
    updateSliderInput(session, "r", value = r_guess)
    updateSliderInput(session, "B", value = B_guess)
  })

  # ---- Auto-fit (robust best fit) ----
  observeEvent(input$bestfit, {
    # Subset to selected range
    rmin <- input$trange[1]
    rmax <- input$trange[2]
    df_sel <- co2 %>% filter(DecimalDate >= rmin, DecimalDate <= rmax)
    t  <- df_sel$DecimalDate
    yy <- df_sel$CO2_Monthly
    tau_sel <- t - t0

    # SSE objective
    sse <- function(par) {
      A <- par[1]; r <- par[2]; B <- par[3]
      yhat <- A * exp(r * tau_sel) + B
      sum((yy - yhat)^2)
    }

    # ---- Start 1: current slider values ----
    par1 <- c(max(0, as.numeric(input$A)),
              max(0, as.numeric(input$r)),
              as.numeric(input$B))

    # ---- Start 2: log-linear seeded values ----
    B0 <- min(yy) - 5
    y_adj <- pmax(yy - B0, 1e-6)
    fit_log <- try(lm(log(y_adj) ~ tau_sel), silent = TRUE)
    if (!inherits(fit_log, "try-error")) {
      r_init <- as.numeric(coef(fit_log)[2])
      A_init <- as.numeric(exp(coef(fit_log)[1]))
    } else {
      r_init <- 0.01; A_init <- 10
    }
    par2 <- c(max(0, A_init), max(0, r_init), B0)

    # ---- Bounds (keep things reasonable) ----
    lower <- c(0,     0,      min(yy) - 100)
    upper <- c(500,   0.10,   max(yy) + 100)

    # ---- Optimize from both starts and keep the best ----
    fit1 <- optim(par1, sse, method = "L-BFGS-B", lower = lower, upper = upper)
    fit2 <- optim(par2, sse, method = "L-BFGS-B", lower = lower, upper = upper)
    best <- if (fit1$value <= fit2$value) fit1$par else fit2$par

    # Update sliders with 3 s.f.
    updateSliderInput(session, "A", value = signif(best[1], 3))
    updateSliderInput(session, "r", value = signif(best[2], 3))
    updateSliderInput(session, "B", value = signif(best[3], 3))
  })

  # ---- Reactive fit data ----
  fit_data <- reactive({
    A <- sf3(input$A)
    r <- sf3(input$r)
    B <- sf3(input$B)
    yhat_all <- A * exp(r * tau) + B
    df_all <- tibble(DecimalDate = co2$DecimalDate,
                     Observed = y,
                     Fitted = yhat_all)

    rmin <- input$trange[1]
    rmax <- input$trange[2]
    df_sel <- df_all %>%
      filter(DecimalDate >= rmin, DecimalDate <= rmax) %>%
      mutate(Residual = Observed - Fitted)

    list(A = A, r = r, B = B,
         rmin = rmin, rmax = rmax,
         all = df_all, sel = df_sel)
  })

  # ---- Main plot ----
  output$p_fit <- renderPlot({
    fd <- fit_data()
    ggplot() +
      geom_line(data = fd$all, aes(DecimalDate, Observed), color = "gray70", alpha = 0.5) +
      geom_line(data = fd$sel, aes(DecimalDate, Observed), color = "black") +
      geom_line(data = fd$sel, aes(DecimalDate, Fitted), color = "blue", linewidth = 1) +
      scale_x_continuous(labels = label_sf3) +
      scale_y_continuous(labels = label_sf3) +
      labs(
        title = "Exponential Fit to CO\u2082 Data",
        subtitle = paste(
          "t₀ =", fmt3(t0),
          "| A =", fmt3(fd$A),
          "| r =", fmt3(fd$r), "/yr",
          "| B =", fmt3(fd$B),
          "| Range:", fmt3(fd$rmin), "–", fmt3(fd$rmax)
        ),
        x = "Year",
        y = expression("CO"[2]*" (ppm)")
      ) +
      theme_minimal(base_size = 13)
  })

  # ---- RMSE text ----
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
