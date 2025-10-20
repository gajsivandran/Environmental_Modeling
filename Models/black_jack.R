# app.R — Blackjack: Monte Carlo + Play (hidden hole card + W/L tracking)
# Tabs:
#   1) Baseline Monte Carlo
#   2) Vs Dealer Upcard
#   3) Play: Deal / Hit / Stand / Shuffle; 1 or 7 decks; S17/H17; counts W/L/P + bankroll
# Notes:
#   - In Play tab, dealer's hole card remains hidden until Stand (or natural resolution).
#   - Unit bet = 1. Payoffs: win +1, loss -1, player blackjack +1.5, push 0.
#   - No splits, doubles, insurance.

library(shiny)
library(ggplot2)
library(dplyr)
library(DT)

# ---------------------------
# Card utilities (shared)
# ---------------------------
rank_of <- function(cards) ((cards - 1) %% 13) + 1               # 1..13 (Ace=1,...,K=13)
suit_of <- function(cards) ((cards - 1) %/% 13) + 1              # 1..4
rank_value <- function(ranks) ifelse(ranks == 1, 11, pmin(ranks, 10))

hand_total <- function(cards) {
  r <- rank_of(cards)
  vals <- rank_value(r)
  total <- sum(vals)
  nA <- sum(r == 1)
  while (total > 21 && nA > 0) { total <- total - 10; nA <- nA - 1 }
  total
}

is_blackjack <- function(cards) {
  if (length(cards) != 2) return(FALSE)
  r <- rank_of(cards)
  (1 %in% r) && any(r %in% c(10,11,12,13))
}

# "Soft" = at least one Ace counted as 11 in best total (<=21)
is_soft <- function(cards) {
  r <- rank_of(cards)
  vals <- rank_value(r)
  total <- sum(vals)
  nA <- sum(r == 1)
  downgraded <- 0
  while (total > 21 && nA > 0) { total <- total - 10; nA <- nA - 1; downgraded <- downgraded + 1 }
  (sum(rank_of(cards) == 1) - downgraded) > 0 && total <= 21
}

# Pretty labels
rank_label <- function(r) c("A","2","3","4","5","6","7","8","9","10","J","Q","K")[r]
suit_label <- function(s) c("\u2660","\u2665","\u2666","\u2663")[s]  # spade, heart, diamond, club
card_label <- function(id_vec) {
  if (length(id_vec) == 0) return(character(0))
  r <- rank_of(id_vec); s <- suit_of(id_vec)
  paste0(rank_label(r), suit_label(s))
}

# ---------------------------
# Single-hand simulator (tabs 1 & 2)
# ---------------------------
play_hand <- function(stand_on = 17, dealer_S17 = TRUE, dealer_upcard_filter = NULL) {
  draw_hand <- function() {
    deck <- sample(1:52)
    draw <- function(n = 1) { out <- deck[seq_len(n)]; deck <<- deck[-seq_len(n)]; out }
    p <- c(draw(), draw())
    d <- c(draw(), draw())  # dealer[1] is upcard
    list(p = p, d = d, draw = draw)
  }

  if (!is.null(dealer_upcard_filter)) {
    repeat {
      h <- draw_hand()
      if (rank_of(h$d[1]) == dealer_upcard_filter) break
    }
  } else {
    h <- draw_hand()
  }

  p <- h$p; d <- h$d; draw <- h$draw

  # Naturals
  p_bj <- is_blackjack(p)
  d_bj <- is_blackjack(d)
  if (p_bj || d_bj) {
    if (p_bj && !d_bj) return(list(result = "player_blackjack", payoff = 1.5))
    if (!p_bj && d_bj) return(list(result = "dealer_blackjack", payoff = -1))
    return(list(result = "push_blackjacks", payoff = 0))
  }

  # Player: hit to threshold
  while (hand_total(p) < stand_on) {
    p <- c(p, draw())
    if (hand_total(p) > 21) return(list(result = "player_bust", payoff = -1))
  }

  # Dealer
  repeat {
    dtot <- hand_total(d)
    if (dtot > 21) return(list(result = "dealer_bust", payoff = +1))
    if (dtot > 17) break
    if (dtot == 17) {
      if (dealer_S17) break
      if (!is_soft(d)) break
    }
    d <- c(d, draw())
  }

  pt <- hand_total(p); dt <- hand_total(d)
  if (pt > dt) return(list(result = "player_win", payoff = +1))
  if (pt < dt) return(list(result = "dealer_win", payoff = -1))
  list(result = "push", payoff = 0)
}

summarize_results <- function(res, n, stand_on, dealer_S17, label_extra = NULL) {
  outcomes <- vapply(res, `[[`, character(1), "result")
  payoffs  <- vapply(res, `[[`, numeric(1),   "payoff")
  tab <- sort(table(outcomes))
  probs <- as.numeric(tab) / n
  tib <- tibble(outcome = names(tab), count = as.integer(tab), prob = round(probs, 4))
  list(
    n = n,
    stand_on = stand_on,
    dealer_rule = if (dealer_S17) "S17 (stand soft 17)" else "H17 (hit soft 17)",
    extra = label_extra,
    table = tib,
    mean_return_per_dollar = mean(payoffs),
    sd_return = sd(payoffs),
    house_edge_pct = -100 * mean(payoffs)
  )
}

simulate_blackjack <- function(n = 10000, stand_on = 17, dealer_S17 = TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  res <- replicate(n, play_hand(stand_on = stand_on, dealer_S17 = dealer_S17), simplify = FALSE)
  summarize_results(res, n, stand_on, dealer_S17, label_extra = NULL)
}

simulate_blackjack_vs_upcard <- function(n = 10000, stand_on = 17, dealer_S17 = TRUE,
                                         upcard_value = 10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  res <- replicate(
    n,
    play_hand(stand_on = stand_on, dealer_S17 = dealer_S17, dealer_upcard_filter = upcard_value),
    simplify = FALSE
  )
  label <- paste0("Dealer upcard = ", ifelse(upcard_value == 11, "Ace", upcard_value))
  summarize_results(res, n, stand_on, dealer_S17, label_extra = label)
}

# ---------------------------
# UI
# ---------------------------
ui <- fluidPage(
  titlePanel("Blackjack — Monte Carlo & Play"),
  tabsetPanel(
    id = "tabs",

    # ---- Tab 1: Baseline
    tabPanel(
      title = "Baseline",
      sidebarLayout(
        sidebarPanel(
          numericInput("n", "Number of simulated hands", value = 20000, min = 1000, step = 1000),
          sliderInput("stand_on", "Player stands at (simple policy):", min = 12, max = 20, value = 17),
          radioButtons("dealer_rule", "Dealer rule:",
                       choices = c("S17 (stand on soft 17)" = "S17", "H17 (hit soft 17)" = "H17"),
                       selected = "S17"),
          numericInput("seed", "Random seed (optional):", value = NA),
          actionButton("run", "Run simulations", class = "btn-primary"),
          br(), br(),
          helpText("Assumptions: fresh deck each hand (no counting), no splits/doubles/insurance.")
        ),
        mainPanel(
          fluidRow(
            column(4, h4("Expected return per $1"), textOutput("mean_return"),
                   tags$small("Positive = player edge; Negative = house edge")),
            column(4, h4("House edge (%)"), textOutput("house_edge")),
            column(4, h4("Std. dev. of payoff"), textOutput("sd_return"))
          ),
          hr(),
          h4("Outcome frequencies"),
          DTOutput("tbl"),
          hr(),
          h4("Outcome distribution"),
          plotOutput("barplot", height = 320)
        )
      )
    ),

    # ---- Tab 2: Vs Dealer Upcard
    tabPanel(
      title = "Vs Dealer Upcard",
      sidebarLayout(
        sidebarPanel(
          numericInput("n2", "Number of simulated hands", value = 20000, min = 1000, step = 1000),
          sliderInput("upcard", "Dealer upcard:", min = 2, max = 11, value = 10, step = 1, ticks = TRUE),
          helpText("11 = Ace"),
          sliderInput("stand_on2", "Player stands at (for this upcard):", min = 12, max = 20, value = 16),
          radioButtons("dealer_rule2", "Dealer rule:",
                       choices = c("S17 (stand on soft 17)" = "S17", "H17 (hit soft 17)" = "H17"),
                       selected = "S17"),
          numericInput("seed2", "Random seed (optional):", value = NA),
          actionButton("run2", "Run conditional simulations", class = "btn-primary"),
          br(), br(),
          helpText("Only simulates hands with the chosen visible dealer upcard.")
        ),
        mainPanel(
          fluidRow(
            column(4, h4("Expected return per $1"), textOutput("mean_return2"),
                   tags$small("Positive = player edge; Negative = house edge")),
            column(4, h4("House edge (%)"), textOutput("house_edge2")),
            column(4, h4("Std. dev. of payoff"), textOutput("sd_return2"))
          ),
          hr(),
          uiOutput("cond_header"),
          DTOutput("tbl2"),
          hr(),
          h4("Outcome distribution"),
          plotOutput("barplot2", height = 320)
        )
      )
    ),

    # ---- Tab 3: Play (interactive blackjack)
    tabPanel(
      title = "Play",
      sidebarLayout(
        sidebarPanel(
          selectInput("num_decks", "Decks in shoe:", choices = c("1 Deck" = 1, "7 Decks" = 7), selected = 1),
          radioButtons("dealer_rule_play", "Dealer rule:",
                       choices = c("S17 (stand on soft 17)" = "S17", "H17 (hit soft 17)" = "H17"),
                       selected = "S17"),
          actionButton("shuffle", "Shuffle / New Shoe", class = "btn-secondary"),
          tags$hr(),
          actionButton("deal", "Deal", class = "btn-primary"),
          actionButton("hit",  "Hit"),
          actionButton("stand","Stand"),
          tags$hr(),
          h5("Shoe status"),
          textOutput("shoe_counts"),
          textOutput("ten_prob"),
          textOutput("ten_or_ace_prob"),
          helpText("10-value = 10/J/Q/K. 10-or-Ace = (10/J/Q/K/A)."),
          tags$hr(),
          h5("Session stats"),
          textOutput("hands_txt"),
          textOutput("wl_txt"),
          textOutput("bankroll_txt")
        ),
        mainPanel(
          fluidRow(
            column(6,
                   h4("Player"),
                   uiOutput("player_cards"),
                   h5(textOutput("player_total"))
            ),
            column(6,
                   h4("Dealer"),
                   uiOutput("dealer_cards"),
                   h5(textOutput("dealer_total"))
            )
          ),
          hr(),
          h4("Result"),
          htmlOutput("result_text"),
          tags$small("No splits, doubles, or insurance in this mode.")
        )
      )
    )
  )
)

# ---------------------------
# Server
# ---------------------------
server <- function(input, output, session) {

  # ======== Tab 1: Baseline ========
  sim <- eventReactive(input$run, {
    dealer_S17 <- (input$dealer_rule == "S17")
    seed_val <- if (is.na(input$seed)) NULL else as.integer(input$seed)
    simulate_blackjack(
      n = as.integer(input$n),
      stand_on = as.integer(input$stand_on),
      dealer_S17 = dealer_S17,
      seed = seed_val
    )
  }, ignoreInit = TRUE)

  output$mean_return <- renderText({ req(sim()); sprintf("%.4f", sim()$mean_return_per_dollar) })
  output$house_edge  <- renderText({ req(sim()); sprintf("%.3f", sim()$house_edge_pct) })
  output$sd_return   <- renderText({ req(sim()); sprintf("%.4f", sim()$sd_return) })

  output$tbl <- renderDT({
    req(sim())
    dat <- sim()$table %>% arrange(desc(count))
    datatable(
      dat, rownames = FALSE,
      options = list(pageLength = 10),
      caption = htmltools::tags$caption(
        style = "caption-side: top; text-align: left;",
        sprintf("n = %s, Player stands on %s, Dealer: %s",
                scales::comma(sim()$n), sim()$stand_on, sim()$dealer_rule)
      )
    )
  })

  output$barplot <- renderPlot({
    req(sim())
    dat <- sim()$table %>% mutate(outcome = factor(outcome, levels = outcome))
    ggplot(dat, aes(x = outcome, y = count)) +
      geom_col() +
      geom_text(aes(label = paste0(count, " (", sprintf("%.2f%%", 100*prob), ")")),
                vjust = -0.3, size = 3) +
      labs(x = NULL, y = "Count") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 20, hjust = 1))
  })

  # ======== Tab 2: Vs Dealer Upcard ========
  sim2 <- eventReactive(input$run2, {
    dealer_S17 <- (input$dealer_rule2 == "S17")
    seed_val <- if (is.na(input$seed2)) NULL else as.integer(input$seed2)
    simulate_blackjack_vs_upcard(
      n = as.integer(input$n2),
      stand_on = as.integer(input$stand_on2),
      dealer_S17 = dealer_S17,
      upcard_value = as.integer(input$upcard),
      seed = seed_val
    )
  }, ignoreInit = TRUE)

  output$mean_return2 <- renderText({ req(sim2()); sprintf("%.4f", sim2()$mean_return_per_dollar) })
  output$house_edge2  <- renderText({ req(sim2()); sprintf("%.3f", sim2()$house_edge_pct) })
  output$sd_return2   <- renderText({ req(sim2()); sprintf("%.4f", sim2()$sd_return) })

  output$cond_header <- renderUI({
    req(sim2())
    tags$h4(sprintf("Outcome frequencies — %s", sim2()$extra %||% ""))
  })

  output$tbl2 <- renderDT({
    req(sim2())
    dat <- sim2()$table %>% arrange(desc(count))
    datatable(
      dat, rownames = FALSE,
      options = list(pageLength = 10),
      caption = htmltools::tags$caption(
        style = "caption-side: top; text-align: left;",
        sprintf("n = %s, Player stands on %s, Dealer: %s, %s",
                scales::comma(sim2()$n), sim2()$stand_on, sim2()$dealer_rule, sim2()$extra %||% "")
      )
    )
  })

  output$barplot2 <- renderPlot({
    req(sim2())
    dat <- sim2()$table %>% mutate(outcome = factor(outcome, levels = outcome))
    ggplot(dat, aes(x = outcome, y = count)) +
      geom_col() +
      geom_text(aes(label = paste0(count, " (", sprintf("%.2f%%", 100*prob), ")")),
                vjust = -0.3, size = 3) +
      labs(x = NULL, y = "Count") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 20, hjust = 1))
  })

  # ======== Tab 3: Play (interactive blackjack) ========
  # Reactive shoe + game state + session stats
  rv <- reactiveValues(
    shoe = integer(0),
    decks = 1L,
    player = integer(0),
    dealer = integer(0),
    in_round = FALSE,
    resolved = FALSE,
    result = "",
    # Stats
    bankroll = 0,   # cumulative units
    hands = 0L,
    wins = 0L,
    losses = 0L,
    pushes = 0L
  )

  # Build a fresh shoe
  build_shoe <- function(n_decks = 1L) {
    shoe <- unlist(rep(list(1:52), n_decks))
    sample(shoe)
  }

  # Draw a card from rv$shoe; auto-reshuffle if empty
  draw_card <- function() {
    if (length(rv$shoe) == 0) rv$shoe <- build_shoe(rv$decks)
    c <- rv$shoe[[1]]
    rv$shoe <- rv$shoe[-1]
    c
  }

  # Probability helpers
  prob_next_10_value <- function() {
    if (length(rv$shoe) == 0) return(NA_real_)
    r <- rank_of(rv$shoe)
    mean(r %in% c(10,11,12,13))
  }
  prob_next_10_or_ace <- function() {
    if (length(rv$shoe) == 0) return(NA_real_)
    r <- rank_of(rv$shoe)
    mean(r %in% c(1,10,11,12,13))
  }

  # Settlement helper: update W/L/P, bankroll, hands, and result text
  settle <- function(payoff, msg) {
    rv$bankroll <- rv$bankroll + payoff
    rv$hands <- rv$hands + 1L
    if (payoff > 0) rv$wins <- rv$wins + 1L
    else if (payoff < 0) rv$losses <- rv$losses + 1L
    else rv$pushes <- rv$pushes + 1L
    rv$result <- msg
    rv$resolved <- TRUE
    rv$in_round <- FALSE
  }

  # Reset shoe (and keep stats) when decks selection changes
  observeEvent(input$num_decks, {
    rv$decks <- as.integer(input$num_decks)
    rv$shoe <- build_shoe(rv$decks)
    rv$player <- integer(0); rv$dealer <- integer(0)
    rv$in_round <- FALSE; rv$resolved <- FALSE; rv$result <- ""
  }, ignoreInit = FALSE)

  # Shuffle button: resets shoe and clears current hand (keeps stats)
  observeEvent(input$shuffle, {
    rv$shoe <- build_shoe(rv$decks)
    rv$player <- integer(0); rv$dealer <- integer(0)
    rv$in_round <- FALSE; rv$resolved <- FALSE; rv$result <- ""
  })

  # Deal button
  observeEvent(input$deal, {
    req(length(rv$shoe) >= 4)
    if (rv$in_round) return(NULL)  # ignore if a round is active
    rv$player <- c(draw_card(), draw_card())
    rv$dealer <- c(draw_card(), draw_card())
    rv$in_round <- TRUE
    rv$resolved <- FALSE
    rv$result <- ""

    # Check naturals; if round ends here, reveal dealer hole card and settle
    p_bj <- is_blackjack(rv$player)
    d_bj <- is_blackjack(rv$dealer)
    if (p_bj || d_bj) {
      if (p_bj && !d_bj)      settle(1.5, "Player Blackjack! +1.5")
      else if (!p_bj && d_bj) settle(-1,  "Dealer Blackjack. -1")
      else                    settle(0,   "Both Blackjack — Push.")
    }
  })

  # Hit button
  observeEvent(input$hit, {
    req(rv$in_round, !rv$resolved)
    rv$player <- c(rv$player, draw_card())
    if (hand_total(rv$player) > 21) {
      settle(-1, "Player busts. -1")
    }
  })

  # Stand button (dealer plays out; then reveal & settle)
  observeEvent(input$stand, {
    req(rv$in_round, !rv$resolved)
    dealer_S17 <- (input$dealer_rule_play == "S17")

    repeat {
      dt <- hand_total(rv$dealer)
      if (dt > 21) { settle(+1, "Dealer busts! +1"); return(NULL) }
      if (dt > 17) break
      if (dt == 17) {
        if (dealer_S17) break
        if (!is_soft(rv$dealer)) break
      }
      rv$dealer <- c(rv$dealer, draw_card())
    }

    pt <- hand_total(rv$player); dt <- hand_total(rv$dealer)
    if (pt > dt)      settle(+1, "Player wins! +1")
    else if (pt < dt) settle(-1, "Dealer wins. -1")
    else              settle(0,  "Push.")
  })

  # --- Outputs for Play tab
  # Player cards are always fully shown
  output$player_cards <- renderUI({
    labs <- card_label(rv$player)
    if (length(labs) == 0) tags$p("(no cards)") else tags$p(paste(labs, collapse = "  "))
  })
  output$player_total <- renderText({
    if (length(rv$player) == 0) return("")
    sprintf("Total: %d", hand_total(rv$player))
  })

  # Dealer cards: hide hole card until round resolves
  output$dealer_cards <- renderUI({
    if (length(rv$dealer) == 0) return(tags$p("(no cards)"))
    labs <- card_label(rv$dealer)
    if (rv$in_round && !rv$resolved) {
      # show only upcard and a hidden placeholder
      tags$p(paste0(labs[1], "  [\u25A0\u25A0]"))  # ■■ as hidden card symbol
    } else {
      # show all after Stand or natural resolution
      tags$p(paste(labs, collapse = "  "))
    }
  })
  output$dealer_total <- renderText({
    if (length(rv$dealer) == 0) return("")
    if (rv$in_round && !rv$resolved) {
      # Only show upcard label while hole card is hidden
      up <- card_label(rv$dealer[1])
      sprintf("Showing: %s", up)
    } else {
      sprintf("Total: %d", hand_total(rv$dealer))
    }
  })

  output$result_text <- renderUI({
    if (rv$result == "") HTML("<em>Click Deal to begin.</em>")
    else HTML(sprintf("<strong>%s</strong>", rv$result))
  })

  output$shoe_counts <- renderText({
    sprintf("Cards remaining in shoe: %d", length(rv$shoe))
  })
  output$ten_prob <- renderText({
    p <- prob_next_10_value()
    if (is.na(p)) "P(10-value next): n/a" else sprintf("P(10-value next): %.2f%%", 100*p)
  })
  output$ten_or_ace_prob <- renderText({
    p <- prob_next_10_or_ace()
    if (is.na(p)) "P(10-or-Ace next): n/a" else sprintf("P(10-or-Ace next): %.2f%%", 100*p)
  })

  # Stats text
  output$hands_txt <- renderText({ sprintf("Hands played: %d", rv$hands) })
  output$wl_txt <- renderText({ sprintf("Wins: %d  |  Losses: %d  |  Pushes: %d", rv$wins, rv$losses, rv$pushes) })
  output$bankroll_txt <- renderText({
    sprintf("Bankroll (units): %+.1f", rv$bankroll)
  })
}

shinyApp(ui, server)
