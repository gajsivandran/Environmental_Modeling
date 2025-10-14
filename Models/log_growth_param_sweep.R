# cod_influence_2025.R
# ------------------------------------------------------------
# Goal: quantify influence of r and K on cod population in year 2025
# Model: P(t) = K / (1 + A * exp(-r * t)), with t = Year - min(Year)
# HOLD A at 5.29   (as requested)
# Outputs:
#   - out/cod_2025_influence.csv (grid with sensitivities & elasticities)
#   - plots/cod_2025_heatmap.png
#   - plots/cod_2025_slices.png
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

# ----------------------------
# Settings (edit if needed)
# ----------------------------
A_fixed <- 5.29         # HOLD A at 5.29
data_path <- "data/cod_timeseries.csv"

# ----------------------------
# Data + time origin (t = Year - min_year)
# ----------------------------
load_cod_data <- function() {
  path <- "data/cod_timeseries.csv"
  read_csv(path, show_col_types = FALSE) |>
    select(Year, Pop) |>
    arrange(Year)

}

cod <- load_cod_data() |>
  mutate(t = Year - min(Year))

  min_year <- min(cod$Year, na.rm = TRUE)
  t_2025   <- 2025 - min_year

  # Suggest K range from data
  pop_max <- max(cod$Pop, na.rm = TRUE)
  K_min_suggest <- 0.5 * pop_max
  K_max_suggest <- 3.0 * pop_max


# ----------------------------
# Logistic & sensitivities
# ----------------------------
logistic_fun <- function(t, K, A, r) {
  K / (1 + A * exp(-r * t))
}

# P = K / E, E = 1 + A e^{-rt}
# dP/dK = 1/E = P/K
# dP/dr = K * (A t e^{-rt}) / E^2
# Elasticities: (K/P)(dP/dK) = 1;  (r/P)(dP/dr) = r * K * (A t e^{-rt}) / (P E^2)
dP_dK_fun <- function(t, K, A, r) {
  E <- 1 + A * exp(-r * t)
  1 / E
}

dP_dr_fun <- function(t, K, A, r) {
  E <- 1 + A * exp(-r * t)
  K * (A * t * exp(-r * t)) / (E^2)
}

elasticity_K <- function(t, K, A, r) {
  1.0
}

elasticity_r <- function(t, K, A, r) {
  P  <- logistic_fun(t, K, A, r)
  E  <- 1 + A * exp(-r * t)
  (r / P) * (K * (A * t * exp(-r * t)) / (E^2))
}

# ----------------------------
# Parameter grids
# ----------------------------
# r range: allow negative to positive to see threshold behavior
r_seq <- seq(0, 0.5, length.out = 201)

# K range: from data-informed suggestion (or manual above)
K_seq <- seq(K_min_suggest, K_max_suggest, length.out = 201)

# ----------------------------
# Compute grid at year 2025
# ----------------------------
grid <- tidyr::crossing(r = r_seq, K = K_seq) |>
  mutate(
    t_year   = t_2025,
    A        = A_fixed,
    P2025    = logistic_fun(t_year, K, A, r),
    dP_dK    = dP_dK_fun(t_year, K, A, r),
    dP_dr    = dP_dr_fun(t_year, K, A, r),
    e_K      = elasticity_K(t_year, K, A, r),
    e_r      = elasticity_r(t_year, K, A, r)
  )

# ----------------------------
# Save outputs
# ----------------------------
dir.create("out",   showWarnings = FALSE, recursive = TRUE)
dir.create("plots", showWarnings = FALSE, recursive = TRUE)

readr::write_csv(grid, file = "out/cod_2025_influence.csv")
message("Wrote grid to out/cod_2025_influence.csv")

# ----------------------------
# Plots
# ----------------------------
# Heatmap of P2025 across r and K
p_heat <- ggplot(grid, aes(r, K, fill = P2025)) +
  geom_raster(interpolate = TRUE) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, alpha = 0.7) +
  scale_fill_viridis_c(name = "P(2025)") +
  labs(
    title = "Cod population in 2025 across r and K",
    subtitle = paste0("A fixed at ", signif(A_fixed, 3),
                      " · t_2025 = ", t_2025, " (Year 2025 − ", min_year, ")"),
    x = "r (growth rate)",
    y = "K (carrying capacity)"
  ) +
  theme_classic()

ggsave("plots/cod_2025_heatmap.png", p_heat, width = 8, height = 5.2, dpi = 300)
message("Saved plots/cod_2025_heatmap.png")

# Slices: P2025 vs K for selected r; and P2025 vs r for selected K
r_vals <- c(0, 0.05, 0.2)
slice_K <- tibble(K = K_seq) |>
  tidyr::crossing(r = r_vals) |>
  mutate(
    A     = A_fixed,
    P2025 = logistic_fun(t_2025, K, A, r),
    r_lab = paste0("r = ", signif(r, 3))
  )

K_vals <- quantile(K_seq, probs = c(0.15, 0.5, 0.85), names = FALSE)
slice_r <- tibble(r = r_seq) |>
  tidyr::crossing(K = as.numeric(K_vals)) |>
  mutate(
    A     = A_fixed,
    P2025 = logistic_fun(t_2025, K, A, r),
    K_lab = paste0("K = ", formatC(K, format = "fg", digits = 3))
  )

p_slices <- (
  ggplot(slice_K, aes(K, P2025, linetype = r_lab)) +
    geom_line() +
    labs(title = "P(2025) vs K (selected r values)",
         x = "K", y = "P(2025)", linetype = "r") +
    theme_classic()
) /
  (
    ggplot(slice_r, aes(r, P2025, linetype = K_lab)) +
      geom_line() +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, alpha = 0.7) +
      labs(title = "P(2025) vs r (selected K values)",
           x = "r", y = "P(2025)", linetype = "K") +
      theme_classic()
  )

# patchwork is optional; if not installed, save two separate plots
save_with_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (save_with_patchwork) {
  p_combined <- patchwork::wrap_plots(p_slices, ncol = 1)
  ggsave("plots/cod_2025_slices.png", p_combined, width = 8, height = 8, dpi = 300)
} else {
  ggsave("plots/cod_2025_slices_top.png",
         ggplot(slice_K, aes(K, P2025, linetype = r_lab)) +
           geom_line() +
           labs(title = "P(2025) vs K (selected r values)",
                x = "K", y = "P(2025)", linetype = "r") +
           theme_classic(),
         width = 8, height = 4, dpi = 300)

  ggsave("plots/cod_2025_slices_bottom.png",
         ggplot(slice_r, aes(r, P2025, linetype = K_lab)) +
           geom_line() +
           geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, alpha = 0.7) +
           labs(title = "P(2025) vs r (selected K values)",
                x = "r", y = "P(2025)", linetype = "K") +
           theme_classic(),
         width = 8, height = 4, dpi = 300)
  message("patchwork not installed; wrote two slice plots instead.")
}

# ----------------------------
# Quick console summary
# ----------------------------
summ <- grid |>
  summarize(
    min_P2025 = min(P2025, na.rm = TRUE),
    med_P2025 = median(P2025, na.rm = TRUE),
    max_P2025 = max(P2025, na.rm = TRUE)
  )

message("Summary P(2025) over grid:\n",
        "  min = ", signif(summ$min_P2025, 3),
        " | median = ", signif(summ$med_P2025, 3),
        " | max = ", signif(summ$max_P2025, 3))

message("Elasticity wrt K is identically 1 (P scales linearly with K). ",
        "Elasticity wrt r varies by (r, K, A, t_2025); see CSV.")
