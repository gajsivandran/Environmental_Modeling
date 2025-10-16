# app.R
# PDFs & CDFs Explorer — with "Sampling from CDF" + Running PDF overlay

library(shiny)
library(ggplot2)
library(scales)
library(dplyr)

ui <- fluidPage(
  titlePanel("PDF & CDF Explorer — Tree Heights (Normal)"),
  sidebarLayout(
    sidebarPanel(
      h4("Dataset: Tree heights (meters)"),
      helpText("Loaded from data/tree_height.csv (expects a column named 'height_m')"),
      sliderInput("a", "Lower bound (a)", min = 0, max = 60, value = 18, step = 0.1),
      sliderInput("b", "Upper bound (b)", min = 0, max = 60, value = 24, step = 0.1),
      checkboxInput("show_kde", "Overlay kernel density (empirical)", value = TRUE),
      checkboxInput("show_theory", "Overlay Normal model (theoretical)", value = TRUE)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Overview",
                 br(),
                 plotOutput("bar_plot", height = 360),
                 br(),
                 verbatimTextOutput("sd_txt"),
                 br(),
                 plotOutput("overview_hist", height = 360)
        ),
        tabPanel("PDF (Density)",
                 br(),
                 plotOutput("pdf_plot", height = 420),
                 br(),
                 verbatimTextOutput("prob_txt")
        ),
        tabPanel("CDF (Cumulative)",
                 br(),
                 plotOutput("cdf_plot", height = 420),
                 br(),
                 verbatimTextOutput("cdf_interval_txt")
        ),

        # ---------- NEW: Sampling from CDF ----------
        tabPanel("Sampling from CDF",
                 br(),
                 fluidRow(
                   column(
                     7,
                     radioButtons(
                       "samp_source", "Sample from:",
                       c("Empirical ECDF (data)" = "emp",
                         "Normal(μ, σ) fit"     = "norm"),
                       inline = TRUE
                     )
                   ),
                   column(
                     5,
                     numericInput("samp_bins", "Histogram bins:", value = 15, min = 5, max = 80, step = 1)
                   )
                 ),
                 fluidRow(
                   column(
                     7,
                     actionButton("draw1", "Draw 1 sample"),
                     actionButton("draw10", "Draw 10 samples"),
                     actionButton("clear_samp", "Clear samples"),
                     tags$small(em("Each draw uses a fresh U ~ Uniform(0,1)"))
                   ),
                   column(5, verbatimTextOutput("last_sample_txt"))
                 ),
                 br(),
                 plotOutput("sample_cdf_plot", height = 360),
                 br(),
                 plotOutput("sample_single_bar", height = 140),
                 br(),
                 plotOutput("sample_hist", height = 300),
                 br(),
                 plotOutput("sample_pdf", height = 320)   # <- NEW running PDF
        ),

        tabPanel("Data (Preview)", tableOutput("head_tbl"))
      )
    )
  )
)

server <- function(input, output, session){

  # ---- Read data safely ----
  tree_data <- reactive({
    df <- read.csv("data/tree_height.csv")
    validate(need("height_m" %in% names(df),
                  "Error: data/tree_height.csv must contain a column named 'height_m'."))
    df
  })

  # ---- Parameters (μ, σ) from the data ----
  params <- reactive({
    x <- tree_data()$height_m
    list(mu = mean(x, na.rm = TRUE), sigma = sd(x, na.rm = TRUE))
  })

  # ======================
  #   OVERVIEW TAB
  # ======================

  output$bar_plot <- renderPlot({
    df <- tree_data()
    x <- df$height_m
    mu <- params()$mu
    df_bar <- data.frame(idx = seq_along(x), height_m = x)

    ggplot(df_bar, aes(idx, height_m)) +
      geom_col(fill = "gray85", color = "gray40") +
      geom_hline(yintercept = mu, color = "darkgreen", linewidth = 1) +
      annotate("text", x = max(df_bar$idx) * 0.98, y = mu,
               label = sprintf("Mean = %.2f m", mu),
               hjust = 1, vjust = -0.6, color = "darkgreen") +
      labs(title = "Tree heights — bar plot with mean overlay",
           x = "Observation index", y = "Height (m)") +
      theme_minimal(base_size = 14)
  })

  output$sd_txt <- renderText({
    mu <- params()$mu
    sigma <- params()$sigma
    sprintf("Mean height (μ): %.3f m\nStandard deviation (σ): %.3f m", mu, sigma)
  })

  output$overview_hist <- renderPlot({
    df <- tree_data()
    x <- df$height_m
    mu <- params()$mu
    sdv <- params()$sigma

    p <- ggplot(df, aes(height_m)) +
      geom_histogram(bins = max(20, floor(sqrt(length(x)))),
                     color = "gray35", fill = "gray85") +
      labs(title = "Histogram with mean and ±1σ / ±2σ markers",
           x = "Height (m)", y = "Count") +
      theme_minimal(base_size = 14)

    p <- p +
      geom_vline(xintercept = mu, color = "black", linewidth = 1) +
      geom_vline(xintercept = c(mu - sdv, mu + sdv),
                 linetype = "dashed", color = "red") +
      geom_vline(xintercept = c(mu - 2*sdv, mu + 2*sdv),
                 linetype = "dotted", color = "blue")

    ymax <- ggplot_build(p)$data[[1]]$y %>% max(na.rm = TRUE)
    p +
      annotate("text", x = mu, y = ymax, label = "μ", vjust = -0.4) +
      annotate("text", x = mu - sdv, y = ymax, label = "μ−1σ", vjust = -0.4, color = "red") +
      annotate("text", x = mu + sdv, y = ymax, label = "μ+1σ", vjust = -0.4, color = "red") +
      annotate("text", x = mu - 2*sdv, y = ymax, label = "μ−2σ", vjust = -0.4, color = "blue") +
      annotate("text", x = mu + 2*sdv, y = ymax, label = "μ+2σ", vjust = -0.4, color = "blue")
  })

  # ======================
  #   PDF / CDF TABS
  # ======================

  prob_info <- reactive({
    x <- tree_data()$height_m
    mu <- params()$mu
    sigma <- params()$sigma
    a <- min(input$a, input$b); b <- max(input$a, input$b)
    p_emp <- mean(x > a & x < b)
    p_th  <- pnorm(b, mu, sigma) - pnorm(a, mu, sigma)
    list(a=a, b=b, emp=p_emp, th=p_th, mu=mu, sigma=sigma)
  })

  output$pdf_plot <- renderPlot({
    df <- tree_data()
    x <- df$height_m
    a <- prob_info()$a; b <- prob_info()$b
    mu <- prob_info()$mu; sigma <- prob_info()$sigma
    xr <- range(x, na.rm = TRUE)

    p <- ggplot(df, aes(x)) +
      geom_histogram(aes(y = after_stat(density)),
                     bins = max(20, floor(sqrt(length(x)))),
                     color = "gray35", fill = "gray85") +
      labs(title = "PDF view — Tree height distribution",
           x = "Height (m)", y = "Density") +
      theme_minimal(base_size = 14)

    grid <- seq(min(xr)-2, max(xr)+2, length.out = 800)
    dens_th <- dnorm(grid, mu, sigma)
    shade_df <- data.frame(grid = grid, dens = dens_th,
                           in_band = grid >= a & grid <= b)

    if (input$show_theory) {
      p <- p +
        geom_area(data = subset(shade_df, in_band),
                  aes(x = grid, y = dens), alpha = 0.3, fill = "skyblue") +
        geom_line(data = data.frame(x = grid, y = dens_th),
                  aes(x, y), linewidth = 1.1, color = "blue")
    }
    if (input$show_kde) {
      p <- p + geom_density(linewidth = 1, linetype = "dashed", color = "black")
    }

    p + geom_vline(xintercept = c(a,b), linetype = "dashed", color = "red") +
      geom_vline(xintercept = mu, color = "darkgreen") +
      coord_cartesian(xlim = c(min(xr)-1, max(xr)+1))
  })

  output$prob_txt <- renderText({
    pi <- prob_info()
    paste0("Empirical P(", pi$a, " < X < ", pi$b, ") = ",
           sprintf("%.3f", pi$emp),
           "\nTheoretical Normal(μ=", round(pi$mu,2), ", σ=", round(pi$sigma,2),
           ") gives P(a < X < b) = ", sprintf("%.3f", pi$th))
  })

  output$cdf_plot <- renderPlot({
    df <- tree_data()
    x <- df$height_m
    mu <- params()$mu; sigma <- params()$sigma
    a <- prob_info()$a; b <- prob_info()$b

    df_ecdf <- data.frame(x = sort(x), F = ecdf(x)(sort(x)))
    grid <- seq(min(x)-2, max(x)+2, length.out = 800)
    F_th <- pnorm(grid, mu, sigma)

    p <- ggplot() +
      geom_step(data = df_ecdf, aes(x, F), linewidth = 1) +
      geom_line(data = data.frame(x = grid, F = F_th),
                aes(x, F), color = "blue", linewidth = 1) +
      labs(title = "CDF view — Empirical (step) vs Theoretical (line)",
           x = "Height (m)", y = "F(x) = P(X ≤ x)") +
      theme_minimal(base_size = 14)

    Fa <- pnorm(a, mu, sigma); Fb <- pnorm(b, mu, sigma)
    p + geom_vline(xintercept = c(a,b), linetype = "dashed", color = "red") +
      geom_hline(yintercept = c(Fa,Fb), linetype = "dotted") +
      annotate("rect", xmin=a, xmax=b, ymin=Fa, ymax=Fb,
               alpha=0.15, fill="skyblue")
  })

  output$cdf_interval_txt <- renderText({
    pi <- prob_info()
    paste0("CDF-based probability:\nP(a < X < b) = F(b) - F(a) = ",
           sprintf("%.3f", pi$th))
  })

  output$head_tbl <- renderTable({ head(tree_data(), 10) })

  # ======================
  #   SAMPLING TAB (state + actions)
  # ======================

  inv_ecdf <- function(u, x) {
    stats::quantile(x, probs = u, type = 1, names = FALSE, na.rm = TRUE)
  }

  samp <- reactiveValues(u = NULL, x = NULL, xs = numeric(0))

  observeEvent(input$draw1, {
    xdat <- tree_data()$height_m
    mu   <- params()$mu; sigma <- params()$sigma
    u <- runif(1)
    x <- if (input$samp_source == "emp") inv_ecdf(u, xdat) else qnorm(u, mu, sigma)
    samp$u <- u; samp$x <- x; samp$xs <- c(samp$xs, x)
  })

  observeEvent(input$draw10, {
    xdat <- tree_data()$height_m
    mu   <- params()$mu; sigma <- params()$sigma
    uu <- runif(10)
    xx <- if (input$samp_source == "emp") sapply(uu, inv_ecdf, x = xdat) else qnorm(uu, mu, sigma)
    samp$u <- tail(uu, 1); samp$x <- tail(xx, 1); samp$xs <- c(samp$xs, xx)
  })

  observeEvent(input$clear_samp, { samp$u <- NULL; samp$x <- NULL; samp$xs <- numeric(0) })

  output$last_sample_txt <- renderText({
    if (is.null(samp$u) || is.null(samp$x)) "No samples yet. Click 'Draw 1 sample' to begin."
    else sprintf("Last draw: U = %.4f  →  X = %.3f (meters)", samp$u, samp$x)
  })

  output$sample_cdf_plot <- renderPlot({
    df <- tree_data()
    x  <- df$height_m
    mu <- params()$mu; sigma <- params()$sigma

    df_ecdf <- data.frame(x = sort(x), F = ecdf(x)(sort(x)))
    grid <- seq(min(x)-2, max(x)+2, length.out = 800)
    F_th <- pnorm(grid, mu, sigma)

    p <- ggplot() +
      geom_step(data = df_ecdf, aes(x, F), linewidth = 1) +
      geom_line(data = data.frame(x = grid, F = F_th),
                aes(x, F), color = "blue", linewidth = 1) +
      labs(title = "Sampling from the CDF — locate U and its mapped X",
           x = "Height (m)", y = "F(x) = P(X ≤ x)") +
      theme_minimal(base_size = 14)

    if (!is.null(samp$u) && !is.null(samp$x)) {
      p <- p +
        geom_hline(yintercept = samp$u, linetype = "dotted") +
        geom_vline(xintercept = samp$x, linetype = "dashed", color = "red") +
        geom_point(aes(x = samp$x, y = samp$u), size = 3)
    }
    p
  })

  output$sample_single_bar <- renderPlot({
    if (is.null(samp$x)) {
      ggplot() + annotate("text", x = 0.5, y = 0.5,
                          label = "Draw a sample to see its bin.", size = 5) +
        theme_void()
    } else {
      xdat <- tree_data()$height_m
      rng  <- range(xdat, na.rm = TRUE)
      brks <- pretty(rng, n = input$samp_bins)
      bn   <- cut(samp$x, breaks = brks, include.lowest = TRUE, right = TRUE)
      dfb  <- data.frame(bin = levels(cut(brks[-1], breaks = brks)),
                         count = as.numeric(levels(cut(brks[-1], breaks = brks)) == as.character(bn)))
      ggplot(dfb, aes(x = bin, y = count)) +
        geom_col(fill = "gray70", color = "gray30") +
        labs(title = "Current draw: which bin did it increment?",
             x = "Height (m) bin", y = "Count (this draw only)") +
        theme_minimal(base_size = 13) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
  })

  output$sample_hist <- renderPlot({
    if (length(samp$xs) == 0) {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No samples yet.", size = 5) +
        theme_void()
    } else {
      df <- data.frame(x = samp$xs)
      ggplot(df, aes(x)) +
        geom_histogram(bins = input$samp_bins, color = "gray35", fill = "gray80") +
        labs(title = sprintf("Running histogram of samples (n = %d)", length(samp$xs)),
             x = "Sampled height (m)", y = "Count") +
        theme_minimal(base_size = 14)
    }
  })

  # ---------- NEW: Running PDF (overlay) ----------
  output$sample_pdf <- renderPlot({
    x_full <- tree_data()$height_m
    mu     <- params()$mu; sigma <- params()$sigma

    # Grid for theoretical normal
    xr   <- range(x_full, na.rm = TRUE)
    grid <- seq(xr[1] - 2, xr[2] + 2, length.out = 800)
    d_th <- dnorm(grid, mean = mu, sd = sigma)

    base <- ggplot() +
      labs(
        title = sprintf("Running PDF (Density) — Samples vs Full Dataset (n_samples = %d)", length(samp$xs)),
        x = "Height (m)", y = "Density"
      ) +
      theme_minimal(base_size = 14)

    # Always plot full dataset KDE and theoretical curve
    base <- base +
      geom_density(data = data.frame(x = x_full),
                   aes(x = x, linetype = "Empirical KDE (full data)"),
                   linewidth = 1) +
      geom_line(data = data.frame(x = grid, y = d_th),
                aes(x, y, linetype = "Normal(μ,σ) fit"),
                linewidth = 1)

    # Add running sample KDE if we have samples
    if (length(samp$xs) > 1) {
      base <- base +
        geom_density(
          data = data.frame(x = samp$xs),
          aes(x = x,
              color    = "Running KDE (samples)",
              linetype = "Running KDE (samples)"),
          linewidth = 1
        )
    }

    base +
      scale_color_manual(
        values = c(
          "Empirical KDE (full data)" = "blue",
          "Normal(μ,σ) fit"           = "green",
          "Running KDE (samples)"     = "red"
        ),
        name = NULL
      ) +
      scale_linetype_manual(
        values = c(
          "Empirical KDE (full data)" = "dashed",
          "Normal(μ,σ) fit"           = "solid",
          "Running KDE (samples)"     = "dotdash"
        ),
        name = NULL
      ) +
      guides(
        color = guide_legend(override.aes = list(size = 1.2)),
        linetype = guide_legend(override.aes = list(size = 1.2))
      ) +
      theme(
        legend.position = "top",
        legend.text = element_text(size = 12)
      )

  })
}

shinyApp(ui, server)
