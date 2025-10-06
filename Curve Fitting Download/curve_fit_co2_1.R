# app.R — Fit-by-eye linear model (sliders start at minimum values)
# Expects data at Models/data/co2_data.txt relative to this app.

library(shiny)
library(ggplot2)
library(dplyr)

# ---------- Helpers for 3 significant figures ----------
sf3  <- function(x) signif(x, 3)               # numeric rounding
fmt3 <- function(x) format(sf3(x), trim = TRUE, scientific = FALSE)
label_sf3 <- function(x) fmt3(x)

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

# Use OLS only to set slider ranges
ols    <- lm(y ~ tau)
a_seed <- unname(coef(ols)[1])
b_seed <- unname(coef(ols)[2])

# Slider bounds
a_min <- a_seed - 30
a_max <- a_seed + 30
b_min <- b_seed - 0.5
b_max <- b_seed + 0.5

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
      plotOutput("p", height = "500px")
    )
  )
)

server <- function(input, output, session) {
  observeEvent(input$reset, {
    updateSliderInput(session, "a", value = a_min)
    updateSliderInput(session, "b", value = b_min)
  })

  output$p <- renderPlot({
    # snap slider values to 3 s.f.
    a <- sf3(input$a)
    b <- sf3(input$b)

    yhat <- a + b * tau
    df   <- tibble(DecimalDate = co2$DecimalDate,
                   Observed = y,
                   Fitted   = yhat)

    ggplot(df, aes(DecimalDate, Observed)) +
      geom_line(alpha = 0.7) +
      geom_line(aes(y = Fitted), linewidth = 1) +
      scale_x_continuous(labels = label_sf3) +
      scale_y_continuous(labels = label_sf3) +
      labs(
        title = "Try to match the data by eye",
        subtitle = paste(
          "t₀ =", fmt3(t0),
          "| a =", fmt3(a),
          "| b =", fmt3(b), "ppm/yr"
        ),
        x = "Year",
        y = expression("CO"[2]*" (ppm)")
      ) +
      theme_minimal(base_size = 13)
  })
}

shinyApp(ui, server)
