# app.R — Coin Flip Simulator + CLT Demo
# - Tab 1: Simulator (1, 2, or 3 coins per trial) with Flip and Flip ×10
# - Tab 2: Central Limit Theorem demo for coin flips (with bias slider p),
#          showing histogram of sample proportions with a normal overlay.

library(shiny)
library(ggplot2)

ui <- navbarPage(
  title = "Coin Flips",
  id = "tabs",

  # ----------------------------
  # Tab 1: Flip Simulator
  # ----------------------------
  tabPanel(
    title = "Simulator",
    fluidPage(
      br(),
      tags$small("Choose 1, 2, or 3 coins per trial. Flip once or flip ×10. Histogram updates live."),
      br(),
      fluidRow(
        column(
          width = 5,
          radioButtons(
            "mode",
            "Flip mode (coins per trial):",
            choices = c("1 coin" = "one", "2 coins" = "two", "3 coins" = "three"),
            selected = "one",
            inline = TRUE
          ),
          div(
            actionButton("flip1", "Flip", class = "btn btn-primary btn-lg me-2"),
            actionButton("flip10", "Flip ×10", class = "btn btn-primary btn-lg me-2"),
            actionButton("reset", "Reset", class = "btn btn-outline-secondary")
          ),
          br(), br(),
          tags$h4("Results"),
          uiOutput("summary_boxes"),
          br(),
          tags$div(
            style = "padding: 12px; border: 1px solid #ddd; border-radius: 8px;",
            uiOutput("prop_box")
          )
        ),
        column(
          width = 7,
          plotOutput("hist", height = "420px")
        )
      )
    )
  ),

  # ----------------------------
  # Tab 2: CLT Demo
  # ----------------------------
  tabPanel(
    title = "CLT Demo",
    fluidPage(
      br(),
      tags$small("Central Limit Theorem with coin flips: sample proportions from repeated samples approach a normal curve as n grows."),
      br(),
      fluidRow(
        column(
          width = 4,
          sliderInput(
            "clt_p", "Probability of Heads (p):", min = 0, max = 1, value = 0.5, step = 0.01
          ),
          sliderInput(
            "clt_n", "Sample size per experiment (n):", min = 1, max = 1000, value = 30, step = 1
          ),
          sliderInput(
            "clt_M", "Number of experiments (M):", min = 100, max = 20000, value = 5000, step = 100
          ),
          actionButton("clt_run", "Run Simulation", class = "btn btn-primary me-2"),
          actionButton("clt_clear", "Clear", class = "btn btn-outline-secondary")
        ),
        column(
          width = 8,
          plotOutput("clt_hist", height = "420px"),
          br(),
          uiOutput("clt_summary")
        )
      )
    )
  )
)

server <- function(input, output, session) {

  # ============================
  # Simulator state & helpers
  # ============================
  rv <- reactiveValues(
    outcomes1 = character(0),   # 1 coin: "Heads"/"Tails"
    outcomes2 = character(0),   # 2 coins: "HH"/"HT"/"TT" (HT includes TH)
    outcomes3 = character(0),   # 3 coins: "HHH"/"HHT"/"HTT"/"TTT"
    last1 = NA_character_,
    last2 = NA_character_,
    last3 = NA_character_
  )

  run_trials_one <- function(k) {
    res <- if (k <= 0) character(0) else sample(c("Heads", "Tails"), size = k, replace = TRUE)
    list(all = res, last = if (length(res)) tail(res, 1) else NA_character_)
  }

  run_trials_two <- function(k) {
    if (k <= 0) return(list(all = character(0), last = NA_character_))
    draws <- matrix(sample(c("H", "T"), size = 2 * k, replace = TRUE), nrow = 2)
    hcount <- colSums(draws == "H")
    cats <- c("TT", "HT", "HH")[hcount + 1L]   # 0->TT, 1->HT, 2->HH
    list(all = cats, last = tail(cats, 1))
  }

  run_trials_three <- function(k) {
    if (k <= 0) return(list(all = character(0), last = NA_character_))
    draws <- matrix(sample(c("H", "T"), size = 3 * k, replace = TRUE), nrow = 3)
    hcount <- colSums(draws == "H")
    map <- c("TTT", "HTT", "HHT", "HHH")       # 0->TTT, 1->HTT, 2->HHT, 3->HHH
    cats <- map[hcount + 1L]
    list(all = cats, last = tail(cats, 1))
  }

  append_trials <- function(mode, k) {
    if (mode == "one") {
      out <- run_trials_one(k)
      rv$outcomes1 <- c(rv$outcomes1, out$all)
      rv$last1 <- out$last
    } else if (mode == "two") {
      out <- run_trials_two(k)
      rv$outcomes2 <- c(rv$outcomes2, out$all)
      rv$last2 <- out$last
    } else {
      out <- run_trials_three(k)
      rv$outcomes3 <- c(rv$outcomes3, out$all)
      rv$last3 <- out$last
    }
  }

  observeEvent(input$flip1,  { append_trials(input$mode, 1L) })
  observeEvent(input$flip10, { append_trials(input$mode, 10L) })

  observeEvent(input$reset, {
    rv$outcomes1 <- character(0)
    rv$outcomes2 <- character(0)
    rv$outcomes3 <- character(0)
    rv$last1 <- NA_character_
    rv$last2 <- NA_character_
    rv$last3 <- NA_character_
  })

  # ---- Summary boxes (mode-aware) ----
  output$summary_boxes <- renderUI({
    if (input$mode == "one") {
      heads <- sum(rv$outcomes1 == "Heads")
      tails <- sum(rv$outcomes1 == "Tails")
      total <- heads + tails
      last <- ifelse(is.na(rv$last1), "—", rv$last1)

      fluidRow(
        column(
          width = 6,
          tags$div(style = "padding:12px; border:1px solid #ddd; border-radius:8px; margin-bottom:8px;",
                   tags$strong("Heads:"), span(style="float:right;", heads)),
          tags$div(style = "padding:12px; border:1px solid #ddd; border-radius:8px;",
                   tags$strong("Tails:"), span(style="float:right;", tails))
        ),
        column(
          width = 6,
          tags$div(style = "padding:12px; border:1px solid #ddd; border-radius:8px; margin-bottom:8px;",
                   tags$strong("Total trials:"), span(style="float:right;", total)),
          tags$div(style = "padding:12px; border:1px solid #ddd; border-radius:8px;",
                   tags$strong("Last flip:"), span(style="float:right;", last))
        )
      )

    } else if (input$mode == "two") {
      hh <- sum(rv$outcomes2 == "HH")
      ht <- sum(rv$outcomes2 == "HT")
      tt <- sum(rv$outcomes2 == "TT")
      total <- hh + ht + tt
      last <- ifelse(is.na(rv$last2), "—", rv$last2)

      fluidRow(
        column(
          width = 6,
          tags$div(style = "padding:12px; border:1px solid #ddd; border-radius:8px; margin-bottom:8px;",
                   tags$strong("HH:"), span(style="float:right;", hh)),
          tags$div(style = "padding:12px; border:1px solid #ddd; border-radius:8px;",
                   tags$strong("HT (mixed):"), span(style="float:right;", ht))
        ),
        column(
          width = 6,
          tags$div(style = "padding:12px; border:1px solid #ddd; border-radius:8px; margin-bottom:8px;",
                   tags$strong("TT:"), span(style="float:right;", tt)),
          tags$div(style = "padding:12px; border:1px solid #ddd; border-radius:8px;",
                   tags$strong("Last pair:"), span(style="float:right;", last))
        ),
        column(
          width = 12,
          br(),
          tags$div(style = "padding:12px; border:1px solid #ddd; border-radius:8px;",
                   tags$strong("Total trials:"), span(style="float:right;", total))
        )
      )

    } else {
      hhh <- sum(rv$outcomes3 == "HHH")
      hht <- sum(rv$outcomes3 == "HHT")   # exactly two heads
      htt <- sum(rv$outcomes3 == "HTT")   # exactly one head
      ttt <- sum(rv$outcomes3 == "TTT")
      total <- hhh + hht + htt + ttt
      last <- ifelse(is.na(rv$last3), "—", rv$last3)

      fluidRow(
        column(
          width = 6,
          tags$div(style = "padding:12px; border:1px solid #ddd; border-radius:8px; margin-bottom:8px;",
                   tags$strong("HHH (3 heads):"), span(style="float:right;", hhh)),
          tags$div(style = "padding:12px; border:1px solid #ddd; border-radius:8px;",
                   tags$strong("HHT (2H 1T):"), span(style="float:right;", hht))
        ),
        column(
          width = 6,
          tags$div(style = "padding:12px; border:1px solid #ddd; border-radius:8px; margin-bottom:8px;",
                   tags$strong("HTT (1H 2T):"), span(style="float:right;", htt)),
          tags$div(style = "padding:12px; border:1px solid #ddd; border-radius:8px;",
                   tags$strong("TTT (3 tails):"), span(style="float:right;", ttt))
        ),
        column(
          width = 12,
          br(),
          tags$div(style = "padding:12px; border:1px solid #ddd; border-radius:8px; margin-bottom:8px;",
                   tags$strong("Total trials:"), span(style="float:right;", total)),
          tags$div(style = "padding:12px; border:1px solid #ddd; border-radius:8px;",
                   tags$strong("Last triple:"), span(style="float:right;", last))
        )
      )
    }
  })

  # ---- Proportion box (mode-aware) ----
  output$prop_box <- renderUI({
    if (input$mode == "one") {
      heads <- sum(rv$outcomes1 == "Heads")
      tails <- sum(rv$outcomes1 == "Tails")
      total <- heads + tails
      ptxt <- if (total == 0) "—" else sprintf("Heads: %.3f | Tails: %.3f", heads/total, tails/total)
      tagList(tags$strong("Outcome proportions:"), span(style="float:right;", ptxt))

    } else if (input$mode == "two") {
      hh <- sum(rv$outcomes2 == "HH")
      ht <- sum(rv$outcomes2 == "HT")
      tt <- sum(rv$outcomes2 == "TT")
      total <- hh + ht + tt
      ptxt <- if (total == 0) "—" else sprintf("HH: %.3f | HT: %.3f | TT: %.3f", hh/total, ht/total, tt/total)
      tagList(tags$strong("Outcome proportions:"), span(style="float:right;", ptxt))

    } else {
      hhh <- sum(rv$outcomes3 == "HHH")
      hht <- sum(rv$outcomes3 == "HHT")
      htt <- sum(rv$outcomes3 == "HTT")
      ttt <- sum(rv$outcomes3 == "TTT")
      total <- hhh + hht + htt + ttt
      ptxt <- if (total == 0) "—" else sprintf("HHH: %.3f | HHT: %.3f | HTT: %.3f | TTT: %.3f",
                                               hhh/total, hht/total, htt/total, ttt/total)
      tagList(tags$strong("Outcome proportions:"), span(style="float:right;", ptxt))
    }
  })

  # ---- Histogram (mode-aware)
  output$hist <- renderPlot({
    if (input$mode == "one") {
      levels_fixed <- c("Heads", "Tails")
      counts <- table(factor(rv$outcomes1, levels = levels_fixed))
      df <- data.frame(Outcome = factor(levels_fixed, levels = levels_fixed),
                       Count = as.integer(counts))
      title_txt <- "Histogram of Outcomes (1 coin per trial)"
    } else if (input$mode == "two") {
      levels_fixed <- c("HH", "HT", "TT")
      counts <- table(factor(rv$outcomes2, levels = levels_fixed))
      df <- data.frame(Outcome = factor(levels_fixed, levels = levels_fixed),
                       Count = as.integer(counts))
      title_txt <- "Histogram of Outcomes (2 coins per trial)"
    } else {
      levels_fixed <- c("HHH", "HHT", "HTT", "TTT")
      counts <- table(factor(rv$outcomes3, levels = levels_fixed))
      df <- data.frame(Outcome = factor(levels_fixed, levels = levels_fixed),
                       Count = as.integer(counts))
      title_txt <- "Histogram of Outcomes (3 coins per trial)"
    }

    ymax <- max(1L, df$Count)
    ggplot(df, aes(x = Outcome, y = Count)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = Count), vjust = -0.4, size = 5) +
      coord_cartesian(ylim = c(0, ymax + ceiling(0.12 * ymax))) +
      labs(title = title_txt, x = NULL, y = "Count") +
      theme_minimal(base_size = 14) +
      theme(panel.grid.minor = element_blank())
  })

  # ============================
  # CLT Demo
  # ============================
  clt_state <- reactiveValues(props = numeric(0))

  observeEvent(input$clt_run, {
    n <- as.integer(input$clt_n)
    M <- as.integer(input$clt_M)
    p <- as.numeric(input$clt_p)

    # simulate M experiments of n flips each; compute sample proportion of Heads for each
    # rbinom is fast and exact for Binomial(n, p)
    heads_counts <- rbinom(M, size = n, prob = p)
    props <- heads_counts / n
    clt_state$props <- props
  })

  observeEvent(input$clt_clear, {
    clt_state$props <- numeric(0)
  })

  output$clt_hist <- renderPlot({
    props <- clt_state$props
    n <- as.integer(input$clt_n)
    p <- as.numeric(input$clt_p)

    # theoretical normal approximation for sample proportion
    mu <- p
    sigma <- sqrt(p * (1 - p) / n)

    # choose sensible bins for proportions
    bins <- max(20, ceiling(sqrt(length(props))))
    df <- data.frame(prop = props)

    gg <- ggplot(df, aes(x = prop)) +
      geom_histogram(aes(y = ..density..), bins = bins, fill = "gray70", color = "white") +
      labs(
        title = sprintf("Sampling Distribution of Sample Proportion (n = %d, p = %.2f, M = %s)", n, p, format(length(props), big.mark=",")),
        x = "Sample proportion of Heads (ȳ = #Heads / n)",
        y = "Density"
      ) +
      theme_minimal(base_size = 14) +
      theme(panel.grid.minor = element_blank())

    # add normal overlay only if sigma > 0 and we have data
    if (length(props) > 0 && sigma > 0) {
      gg <- gg +
        stat_function(fun = dnorm, args = list(mean = mu, sd = sigma), linewidth = 1)
    }

    gg
  })

  output$clt_summary <- renderUI({
    props <- clt_state$props
    n <- as.integer(input$clt_n)
    p <- as.numeric(input$clt_p)
    mu <- p
    sigma <- sqrt(p * (1 - p) / n)

    if (length(props) == 0) {
      HTML("<em>No simulation yet. Choose p, n, M and click <strong>Run Simulation</strong>.</em>")
    } else {
      emp_mean <- mean(props)
      emp_sd   <- sd(props)
      tagList(
        tags$div(
          style = "padding:12px; border:1px solid #ddd; border-radius:8px;",
          HTML(sprintf(
            "<strong>Theoretical (Normal Approx.):</strong> mean = %.4f, sd = %.4f<br>
             <strong>Empirical (from simulation):</strong> mean = %.4f, sd = %.4f<br>
             <small>By the CLT, as n grows the sampling distribution of the sample proportion approaches Normal(%.4f, %.4f²).</small>",
            mu, sigma, emp_mean, emp_sd, mu, sigma
          ))
        )
      )
    }
  })
}

shinyApp(ui, server)
