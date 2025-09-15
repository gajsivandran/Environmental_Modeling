# ----------------------------------------------
# 1-D Advection–Diffusion–Decay (ADE) Template
# Crank–Nicolson (diffusion) + upwind (advection)
# ----------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
})

run_ade <- function(
    L = 10000,            # domain length (m)
    T_hours = 6,          # total time (hours)
    v = 0.2,              # velocity (m/s), + downstream
    D = 10,               # dispersion (m^2/s)
    k = 1e-5,             # decay rate (1/s)
    C_left = 0,           # left boundary concentration (mg/L) - Dirichlet
    C0 = 0,               # background concentration (mg/L)
    spill_center = 1000,  # m
    spill_width  = 200,   # m (total width)
    spill_conc   = 5,     # mg/L (added on top of background)
    dx = 10,              # spatial step (m) -> N ~ L/dx + 1
    CFL_target = 0.8,     # target advective CFL
    stations_m = c(1000, 5000, 9000),     # where to record time series
    snapshot_hours = c(0, 1, 3, 6)        # profiles to plot
) {
  # ---- Derived quantities ----
  T <- T_hours * 3600       # seconds
  x <- seq(0, L, by = dx)
  N <- length(x)

  # Auto time step from advection CFL (diffusion handled implicitly)
  if (abs(v) > 0) {
    dt_adv <- CFL_target * dx / abs(v)
  } else {
    dt_adv <- 60  # fall-back if v == 0
  }
  dt <- min(dt_adv, 120)      # keep dt reasonable if flow is slow
  nt <- ceiling(T / dt)
  dt <- T / nt                # make it divide T exactly
  time <- seq(0, T, by = dt)
  nt <- length(time) - 1

  # Diagnostics
  CFL <- abs(v) * dt / dx
  rD  <- D * dt / dx^2
  Pe  <- if (D > 0) abs(v) * dx / D else Inf
  message(sprintf("Grid: N=%d (dx=%.2f m), steps=%d (dt=%.2f s)", N, dx, nt, dt))
  message(sprintf("CFL (advection)=%.3f, rD (diffusion)=%.3f, Péclet (dx-scale)=%.3f",
                  CFL, rD, Pe))

  # ---- Initial condition: box spill over background ----
  C <- rep(C0, N)
  half_w <- spill_width / 2
  idx_spill <- which(x >= (spill_center - half_w) & x <= (spill_center + half_w))
  C[idx_spill] <- C[idx_spill] + spill_conc

  # ---- Helper: upwind gradient ∂C/∂x (full grid) ----
  upwind_grad <- function(C, v, dx) {
    g <- numeric(length(C))
    if (v >= 0) {
      # backward difference; g[1] unused (left boundary fixed)
      g[2:length(C)] <- (C[2:length(C)] - C[1:(length(C)-1)]) / dx
      g[1] <- (C[2] - C[1]) / dx  # not used but defined
    } else {
      # forward difference; right boundary: zero-gradient (ghost = last)
      g[1:(length(C)-1)] <- (C[2:length(C)] - C[1:(length(C)-1)]) / dx
      g[length(C)] <- 0
    }
    g
  }

  # ---- Build Crank–Nicolson matrices for diffusion on unknowns C[2..N] ----
  # Lm encodes the second derivative stencil WITHOUT the 1/dx^2 factor.
  nU <- N - 1                # unknowns are i = 2..N  (i=1 is fixed boundary)
  main  <- rep(-2, nU)
  upper <- rep( 1, nU - 1)
  lower <- rep( 1, nU - 1)
  # Right boundary Neumann (∂C/∂x=0): last row becomes [-2, +2]
  lower[nU - 1] <- 2

  # Tridiagonal constructor
  tri <- function(lower, main, upper) {
    M <- matrix(0, nU, nU)
    diag(M) <- main
    M[row(M) == col(M) + 1] <- lower  # below diagonal
    M[row(M) + 1 == col(M)] <- upper  # above diagonal
    M
  }
  Lm <- tri(lower, main, upper)  # stencil (no 1/dx^2)

  lam <- D * dt / (2 * dx^2)     # Crank–Nicolson factor
  A <- diag(nU) - lam * Lm       # left-hand matrix
  B <- diag(nU) + lam * Lm       # right-hand matrix

  # Constant boundary contribution for CN at j=1 (depends on C_left)
  b_cn <- numeric(nU); b_cn[1] <- 2 * lam * C_left

  # ---- Storage ----
  C_store <- matrix(NA_real_, nrow = nt + 1, ncol = N)
  C_store[1, ] <- C

  # ---- Time stepping ----
  for (n in 1:nt) {
    # 1) Advection (explicit, upwind)
    grad <- upwind_grad(C, v, dx)
    adv  <- -v * grad

    # 2) Decay (explicit)
    dec  <- -k * C

    # 3) Assemble RHS on unknowns (i = 2..N)
    CU   <- C[2:N]
    rhs  <- as.vector(B %*% CU) + dt * (adv[2:N] + dec[2:N]) + b_cn

    # 4) Solve for new unknowns with diffusion CN
    CUn  <- solve(A, rhs)

    # 5) Recompose solution, enforce boundaries
    Cnew <- C
    Cnew[1] <- C_left                 # Dirichlet at x=0
    Cnew[2:N] <- CUn
    # Neumann at x=L already built into A/B; no extra step needed

    C <- Cnew
    C_store[n + 1, ] <- C
  }

  # ---- Prepare outputs ----
  df_all <- as_tibble(C_store) |>
    mutate(time_s = time) |>
    pivot_longer(cols = -time_s,
                 names_to = "col", values_to = "C") |>
    mutate(i = as.integer(gsub("V", "", col)),
           x = x[i],
           time_h = time_s / 3600) |>
    select(time_h, x, C)

  # Snapshots
  snapshots <- df_all |> filter(time_h %in% snapshot_hours)

  # Time series at stations (nearest grid points)
  pick_idx <- sapply(stations_m, function(s) which.min(abs(x - s)))
  station_labels <- paste0("x = ", round(x[pick_idx]), " m")
  timeseries <- tibble(time_h = time) |>
    mutate(across(everything(), as.numeric))
  for (j in seq_along(pick_idx)) {
    timeseries[[station_labels[j]]] <- C_store[, pick_idx[j]]
  }
  timeseries <- timeseries |> pivot_longer(-time_h,
                                           names_to = "station",
                                           values_to = "C")

  list(x = x,
       time = time,
       Cmat = C_store,
       snapshots = snapshots,
       timeseries = timeseries,
       grid_info = list(dx = dx, dt = dt, CFL = CFL, rD = rD, Pe = Pe))
}

# ---------------------------
# Run with defaults + plots
# ---------------------------
out <- run_ade()

# Snapshot profiles over distance
p1 <- ggplot(out$snapshots, aes(x = x, y = C, group = factor(time_h))) +
  geom_line() +
  facet_wrap(~ paste0("t = ", time_h, " h"), ncol = 2) +
  labs(x = "Distance (m)", y = "Concentration (mg/L)",
       title = "Advection–Diffusion–Decay: Snapshots") +
  theme_minimal(base_size = 12)

# Time series at stations
p2 <- ggplot(out$timeseries, aes(x = time_h, y = C, color = station)) +
  geom_line() +
  labs(x = "Time (hours)", y = "Concentration (mg/L)",
       color = "Station",
       title = "Concentration vs. Time at Selected Stations") +
  theme_minimal(base_size = 12)

print(p1); print(p2)

# ---------------------------------------
# OPTIONAL: Animation (uncomment to use)
# ---------------------------------------
# install.packages(c("gganimate", "gifski"))
# library(gganimate)
# anim <- ggplot(as.data.frame(out$Cmat) |>
#                  mutate(time_h = out$time/3600) |>
#                  pivot_longer(-time_h, names_to="col", values_to="C") |>
#                  mutate(i = as.integ
#                         x = out$x[i]),
#                aes(x = x, y = C)) +
#   geom_line() +
#   labs(x="Distance (m)", y="Concentration (mg/L)",
#        title = "t = {round(frame_time,2)} h") +
#   transition_time(time_h) +
#   ease_aes('linear') +
#   theme_minimal(base_size=12)
# animate(anim, nframes = 200, fps = 20, renderer = gifski_renderer())
