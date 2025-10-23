library(shiny)
library(ggplot2)

ui <- fluidPage(
  titlePanel("Roulette Wheel Simulator with Betting and Net Winnings"),

  sidebarLayout(
    sidebarPanel(
      numericInput("num_spins", "Number of Spins:", value = 100, min = 1),
      selectInput("bet_number", "Pick a Number to Bet On (0–36):", choices = 0:36),
      actionButton("spin_btn", "Spin the Wheel")
    ),

    mainPanel(
      plotOutput("histogram"),
      verbatimTextOutput("summary"),
      verbatimTextOutput("winnings")
    )
  )
)

server <- function(input, output) {
  spins <- reactiveVal(NULL)

  observeEvent(input$spin_btn, {
    wheel <- 0:36
    outcomes <- sample(wheel, input$num_spins, replace = TRUE)
    spins(outcomes)
  })

  output$histogram <- renderPlot({
    req(spins())
    ggplot(data.frame(Number = spins()), aes(x = factor(Number))) +
      geom_bar(fill = "darkgreen", color = "black") +
      labs(title = "Roulette Spin Results", x = "Number", y = "Frequency") +
      theme_minimal()
  })

  output$summary <- renderPrint({
    req(spins())
    table(spins())
  })

  output$winnings <- renderPrint({
    req(spins())
    bet <- as.numeric(input$bet_number)
    hits <- sum(spins() == bet)
    payout <- hits * 35         # Standard payout for single number
    cost <- input$num_spins * 1 # $1 per spin
    net <- payout - cost        # Net winnings
    cat("Your number:", bet, "\n")
    cat("Hits:", hits, "\n")
    cat("Total Winnings: $", payout, "\n")
    cat("Total Cost: $", cost, "\n")
    cat("Net Winnings: $", net, "\n")
  })
}

shinyApp(ui = ui, server = server)
