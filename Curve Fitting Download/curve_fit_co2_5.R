# app.R — Exponential + Sinusoid model for Mauna Loa CO2
# Model: y = A * exp(r*(t - t0)) + C * sin(2*pi*(t - phi)) + B
# Includes range selection, RMSE, and robust Best Fit button (3 s.f.)

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
freq <- 1   # fixed annual cycle

# Initial guesses
A_guess <- 10
r_guess <- 0.01
B_guess <- min(y)
C_guess <- 3
phi_guess <- 0.25

# Slider ranges
A_min <- 0;   A_max <- 200
r_min <- 0;   r_max <- 0.05
B_min <- B_guess - 50; B_max <- B_guess + 100
C_min <- 0;   C_max <- 10
phi_min <- 0; phi_max <- 1

# ---------- UI ----------
ui <- fluidPage(
  titlePanel("Fit by Eye: Exponential + Sinusoid Model for Mauna Loa CO\u2082"),
  sidebarLayout(
    sidebarPanel(
      p("Model:  y = A·exp(r·(t − t₀)) + C·sin(2π·(t − φ)) + B"),
      p(em(paste0("t₀ = ", fmt3(t0), " (fixed), frequency = 1 cycle/year"))),

      sliderInput("A", "A (exponential scale)", min = A_min, max = A_max, value = A_guess, step = 1),
      sliderInput("r", "r (growth rate, 1/year)", min = r_min, max = r_max, value = r_guess, step = 0.001),
      sliderInput("B", "B (baseline, ppm)", min = B_min, max = B_max, value = B_guess, step = 1),
      sliderInput("C", "C (sinusoid amplitude, ppm)", min = C_min, max = C_max, value = C_guess, step = 0.1),
      sliderInput("phi", "φ (phase shift, years)", min = phi_min, max = phi_max, value = phi_guess, step = 0.01),

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
    updateSliderInput(session, "C", value = C_guess)
    updateSliderInput(session, "phi", value = phi_guess)
  })

  # ---- Robust Best Fit using optim() ----
  observeEvent(input$bestfit, {
    rmin <- input$trange[1]
    rmax <- input$trange[2]
    df_sel <- co2 %>% filter(DecimalDate >= rmin, DecimalDate <= rmax)
    t <- df_sel$DecimalDate
    yy <- df_sel$CO2_Monthly
    tau_sel <- t - t0

    # Objective: SSE for exponential + sinusoid
    sse <- function(par) {
      A <- par[1]; r <- par[2]; B <- par[3]; C <- par[4]; phi <- par[5]
      yhat <- A * exp(r * tau_sel) + C * sin(2*pi*freq*(t - phi)) + B
      sum((yy - yhat)^2)
    }

    # Starting points: current sliders
    start <- c(as.numeric(input$A),
               as.numeric(input$r),
               as.numeric(input$B),
               as.numeric(input$C),
               as.numeric(input$phi))

    # Bounds
    lower <- c(0, 0, min(yy)-100, 0, 0)
    upper <- c(500, 0.1, max(yy)+100, 15, 1)

    # Optimize
    fit <- optim(start, sse, method = "L-BFGS-B", lower = lower, upper = upper)

    updateSliderInput(session, "A", value = sf3(fit$par[1]))
    updateSliderInput(session, "r", value = sf3(fit$par[2]))
    updateSliderInput(session, "B", value = sf3(fit$par[3]))
    updateSliderInput(session, "C", value = sf3(fit$par[4]))
    updateSliderInput(session, "phi", value = sf3(fit$par[5]))
  })

  # ---- Reactive data ----
  fit_data <- reactive({
    A <- sf3(input$A)
    r <- sf3(input$r)
    B <- sf3(input$B)
    C <- sf3(input$C)
    phi <- sf3(input$phi)

    yhat_all <- A * exp(r * tau) + C * sin(2*pi*freq*(co2$DecimalDate - phi)) + B
    df_all <- tibble(DecimalDate = co2$DecimalDate,
                     Observed = y,
                     Fitted = yhat_all)

    rmin <- input$trange[1]
    rmax <- input$trange[2]
    df_sel <- df_all %>%
      filter(DecimalDate >= rmin, DecimalDate <= rmax) %>%
      mutate(Residual = Observed - Fitted)

    list(A = A, r = r, B = B, C = C, phi = phi,
         rmin = rmin, rmax = rmax,
         all = df_all, sel = df_sel)
  })

  # ---- Main Plot ----
  output$p_fit <- renderPlot({
    fd <- fit_data()
    ggplot() +
      geom_line(data = fd$all, aes(DecimalDate, Observed), color = "gray70", alpha = 0.5) +
      geom_line(data = fd$sel, aes(DecimalDate, Observed), color = "black") +
      geom_line(data = fd$sel, aes(DecimalDate, Fitted), color = "blue", linewidth = 1) +
      scale_x_continuous(labels = label_sf3) +
      scale_y_continuous(labels = label_sf3) +
      labs(
        title = "Exponential + Sinusoid Fit to CO\u2082 Data",
        subtitle = paste(
          "t₀ =", fmt3(t0),
          "| A =", fmt3(fd$A),
          "| r =", fmt3(fd$r), "/yr",
          "| B =", fmt3(fd$B),
          "| C =", fmt3(fd$C),
          "| φ =", fmt3(fd$phi),
          "| Range:", fmt3(fd$rmin), "–", fmt3(fd$rmax)
        ),
        x = "Year",
        y = expression("CO"[2]*" (ppm)")
      ) +
      theme_minimal(base_size = 13)
  })

  # ---- RMSE ----
  output$rmse_txt <- renderText({
    fd <- fit_data()
    paste0("RMSE on selected range: ",
           fmt3(rmse(fd$sel$Observed, fd$sel$Fitted)), " ppm")
  })

  # ---- Residuals ----
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
