# app.R
# -------------------------------------------------------------
# Interactive curve fitting for Mauna Loa CO2
# Models:
# 1) Linear
# 2) Exponential
# 3) Linear + Sinusoid
# 4) Exponential + Sinusoid
# Students can adjust parameters and see RMSE.
# -------------------------------------------------------------

library(shiny)
library(ggplot2)
library(dplyr)

# ---- Load data ----
co2_data <- read.table("Models/data/co2_data.txt",
                       header = FALSE,
                       comment.char = "#",
                       col.names = c("Year","Month","DecimalDate",
                                     "CO2_Monthly","CO2_Deseason",
                                     "NumDays","StdDev","Uncertainty")) |>
  as_tibble() |>
  filter(!is.na(DecimalDate), !is.na(CO2_Monthly))

t0 <- min(co2_data$DecimalDate, na.rm = TRUE)  # origin for stability

# ---- Helper: RMSE ----
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))

# ---- UI ----
ui <- fluidPage(
  titlePanel("Curve Fitting Explorer: Mauna Loa CO₂"),
  sidebarLayout(
    sidebarPanel(
      selectInput("model_type", "Choose model:",
                  choices = c("Line" = "lin",
                              "Exponential" = "exp",
                              "Line + Sinusoid" = "lin_sin",
                              "Exponential + Sinusoid" = "exp_sin"),
                  selected = "lin_sin"),
      tags$hr(),

      # ----- Linear parameters -----
      conditionalPanel(
        condition = "['lin','lin_sin'].includes(input.model_type)",
        h4("Linear trend:  y = a + b·(t - t0)"),
        sliderInput("a_lin", "a (intercept @ t0)", min = 250, max = 400, value = 315, step = 0.5),
        sliderInput("b_lin", "b (ppm per year)",  min = -1,  max = 5,   value = 1.5, step = 0.05)
      ),

      # ----- Exponential parameters -----
      conditionalPanel(
        condition = "['exp','exp_sin'].includes(input.model_type)",
        h4("Exponential trend:  y = B + A·exp(r·(t - t0))"),
        sliderInput("B_exp", "B (baseline)", min = 250, max = 400, value = 300, step = 0.5),
        sliderInput("A_exp", "A (scale)",    min =  1,  max = 100, value = 30,  step = 1),
        sliderInput("r_exp", "r (rate/yr)",  min = 0,   max = 0.05, value = 0.01, step = 0.001)
      ),

      # ----- Sinusoid parameters -----
      conditionalPanel(
        condition = "['lin_sin','exp_sin'].includes(input.model_type)",
        h4("Seasonality:  + Asin · sin(2π·(t - φ))"),
        sliderInput("A_sin", "Asin (amplitude, ppm)", min = 0, max = 10, value = 3, step = 0.1),
        sliderInput("phi",   "φ (phase shift, years)", min = 0, max = 1, value = 0.3, step = 0.01),
        helpText("Frequency fixed at 1 cycle/year to represent annual seasonality.")
      ),

      checkboxInput("show_deseason", "Also show deseasonalized series", value = FALSE),
      checkboxInput("show_resid", "Show residuals panel", value = FALSE)
    ),

    mainPanel(
      plotOutput("fit_plot", height = "450px"),
      tags$br(),
      strong(textOutput("rmse_text")),
      conditionalPanel(
        condition = "input.show_resid",
        tags$hr(),
        plotOutput("resid_plot", height = "250px")
      )
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {

  preds <- reactive({
    t <- co2_data$DecimalDate
    y <- co2_data$CO2_Monthly
    tau <- t - t0

    # Base components
    trend <- switch(input$model_type,
                    "lin"     = input$a_lin + input$b_lin * tau,
                    "lin_sin" = input$a_lin + input$b_lin * tau,
                    "exp"     = input$B_exp + input$A_exp * exp(input$r_exp * tau),
                    "exp_sin" = input$B_exp + input$A_exp * exp(input$r_exp * tau))

    seas <- switch(input$model_type,
                   "lin"     = 0,
                   "exp"     = 0,
                   "lin_sin" = input$A_sin * sin(2*pi*(t - input$phi)),
                   "exp_sin" = input$A_sin * sin(2*pi*(t - input$phi)))

    yhat <- trend + seas
    list(t = t, y = y, yhat = yhat)
  })

  output$fit_plot <- renderPlot({
    pr <- preds()
    p <- ggplot(co2_data, aes(x = DecimalDate, y = CO2_Monthly)) +
      geom_line(alpha = 0.7) +
      geom_line(aes(y = pr$yhat), linewidth = 1) +
      labs(title = "Curve Fitting Explorer",
           subtitle = paste("Model:", switch(input$model_type,
                                             lin = "Line",
                                             exp = "Exponential",
                                             lin_sin = "Line + Sinusoid",
                                             exp_sin = "Exponential + Sinusoid")),
           x = "Year",
           y = expression("CO"[2]*" (ppm)")) +
      theme_minimal()

    if (input$show_deseason) {
      p <- p + geom_line(aes(y = CO2_Deseason), linetype = "dashed")
    }
    p
  })

  output$rmse_text <- renderText({
    pr <- preds()
    paste0("RMSE: ", round(rmse(pr$y, pr$yhat), 3), " ppm")
  })

  output$resid_plot <- renderPlot({
    pr <- preds()
    resid <- pr$y - pr$yhat
    ggplot(data.frame(DecimalDate = pr$t, Residual = resid),
           aes(x = DecimalDate, y = Residual)) +
      geom_hline(yintercept = 0, linewidth = 0.5) +
      geom_line() +
      labs(title = "Residuals (Observed - Predicted)",
           x = "Year", y = "ppm") +
      theme_minimal()
  })
}

shinyApp(ui, server)
