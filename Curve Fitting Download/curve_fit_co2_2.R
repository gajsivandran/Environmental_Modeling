# app.R — Fit-by-eye linear model + Residuals + RMSE (3 significant figures)
# Expects data at Models/data/co2_data.txt relative to this app.

library(shiny)
library(ggplot2)
library(dplyr)

# ---------- Helpers for 3 significant figures ----------
sf3  <- function(x) signif(x, 3)                                # numeric rounding
fmt3 <- function(x) format(sf3(x), trim = TRUE, scientific = FALSE)  # pretty text
label_sf3 <- function(x) fmt3(x)                                 # axis labels

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
  filter(!is.na(DecimalDate), !is.na(CO2_Monthly))

# Center time for stable sliders
t0  <- min(co2$DecimalDate, na.rm = TRUE)
tau <- co2$DecimalDate - t0
y   <- co2$CO2_Monthly

# Use OLS only to set sensible slider ranges
ols    <- lm(y ~ tau)
a_seed <- unname(coef(ols)[1])
b_seed <- unname(coef(ols)[2])

# Slider bounds
a_min <- a_seed - 30
a_max <- a_seed + 30
b_min <- b_seed - 0.5
b_max <- b_seed + 0.5

# RMSE helper
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))

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
      actionButton("reset", "Reset sliders")
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

server <- function(input, output, session) {
  observeEvent(input$reset, {
    updateSliderInput(session, "a", value = a_min)
    updateSliderInput(session, "b", value = b_min)
  })

  # Reactive fitted values and residuals
  fit_data <- reactive({
    a <- sf3(input$a)  # snap to 3 s.f.
    b <- sf3(input$b)
    yhat <- a + b * tau
    tibble(
      DecimalDate = co2$DecimalDate,
      Observed    = y,
      Fitted      = yhat,
      Residual    = y - yhat,
      a = a, b = b
    )
  })

  output$p_fit <- renderPlot({
    df <- fit_data()

    ggplot(df, aes(DecimalDate, Observed)) +
      geom_line(alpha = 0.7) +
      geom_line(aes(y = Fitted), linewidth = 1) +
      scale_x_continuous(labels = label_sf3) +
      scale_y_continuous(labels = label_sf3) +
      labs(
        title = "Try to match the data by eye",
        subtitle = paste(
          "t₀ =", fmt3(t0),
          "| a =", fmt3(df$a[1]),
          "| b =", fmt3(df$b[1]), "ppm/yr"
        ),
        x = "Year",
        y = expression("CO"[2]*" (ppm)")
      ) +
      theme_minimal(base_size = 13)
  })

  output$rmse_txt <- renderText({
    df <- fit_data()
    paste0("RMSE: ", fmt3(rmse(df$Observed, df$Fitted)), " ppm")
  })

  output$p_resid <- renderPlot({
    df <- fit_data()
    ggplot(df, aes(DecimalDate, Residual)) +
      geom_hline(yintercept = 0) +
      geom_line() +
      scale_x_continuous(labels = label_sf3) +
      scale_y_continuous(labels = label_sf3) +
      labs(title = "Residuals (Observed − Fitted)", x = "Year", y = "ppm") +
      theme_minimal(base_size = 13)
  })
}

shinyApp(ui, server)
