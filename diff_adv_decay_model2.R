# ===============================================================
# 1-D Advection–Diffusion–Decay (ADE) teaching model
#   C_t = D * C_xx - v * C_x - k * C
# Domain: x in [0, L], time in [0, T]
# BCs: Dirichlet at x=0 (fixed C); Neumann at x=L (zero gradient)
# Scheme: upwind advection + FTCS diffusion + explicit decay
# Plots: profiles at t = 0, 1, 3, 6 h; time-series at x = 1, 5, 9 km
# Units: C mg/L, v m/s, D m^2/s, k 1/s
# ===============================================================

# ----- Assumptions (prints 5 bullets + where each breaks) -----
ade_assumptions <- function() {
  cat("\nModeling assumptions (and when they break):\n",
      "  • 1-D flow (no lateral/vertical variation)\n",
      "    – Breaks with lateral stratification or side embayments.\n",
      "  • Constant v, D, k in space and time\n",
      "    – Breaks during storms or strong diurnal temperature swings.\n",
      "  • Well-mixed within each grid cell\n",
      "    – Breaks with hyporheic exchange or strong sub-cell heterogeneity.\n",
      "  • First-order (linear) decay: -k C\n",
      "    – Breaks for saturating/multi-step kinetics (e.g., Michaelis–Menten).\n",
      "  • No internal sources/sinks beyond decay\n",
      "    – Breaks with tributaries, groundwater seepage, outfalls, or sorption.\n", sep = "")
}

# ----- Diagnostics: CFL, r_D, Péclet, Units Table -----
diagnostics <- function(p, dx, dt) {
  CFL   <- abs(p$v_ms) * dt / dx
  r_D   <- p$D_m2s * dt / dx^2
  Pe_dx <- if (p$D_m2s > 0) abs(p$v_ms) * dx / p$D_m2s else Inf
  Pe_L  <- if (p$D_m2s > 0) abs(p$v_ms) * p$L_m / p$D_m2s else Inf

  cat(sprintf("\nStability / scaling:\n  CFL = %.3f   r_D = %.3f   Pe_cell(dx) = %.3f   Pe_domain(L) = %.3f\n",
              CFL, r_D, Pe_dx, Pe_L))
  if (CFL > 0.9) warning(sprintf("CFL %.3f > 0.9: decrease dt or increase dx.", CFL), call. = FALSE)
  if (r_D > 0.5) warning(sprintf("r_D %.3f > 0.5 (explicit FTCS diffusion limit).", r_D), call. = FALSE)

  units_table <- data.frame(
    Quantity = c("Concentration C", "Velocity v", "Diffusivity D", "Decay rate k",
                 "Domain length L", "Grid spacing dx", "Time step dt"),
    Units    = c("mg/L", "m/s", "m^2/s", "1/s", "m", "m", "s"),
    Value    = c(NA, p$v_ms, p$D_m2s, p$k_1s, p$L_m, dx, dt)
  )
  cat("\nUnits / values:\n"); print(units_table, row.names = FALSE)
  invisible(list(CFL = CFL, r_D = r_D, Pe_cell = Pe_dx, Pe_domain = Pe_L, units = units_table))
}

# ----- Main simulator with plotting -----
run_ade <- function(params = list()) {
  # Sensible defaults for the requested figures
  p <- modifyList(list(
    L_m = 10000, nx = 201, dt_s = 5, T_end_h = 6,
    v_ms = 0.20, D_m2s = 5.0, k_1s = 0.0,
    C_left_mgL = 1.0, C_init_mgL = 0.0,
    times_h = c(0, 1, 3, 6), stations_km = c(1, 5, 9),
    clip_nonneg = TRUE, verbose = TRUE, make_plots = TRUE
  ), params)

  # Grid/time
  dx <- p$L_m / (p$nx - 1)
  t_end_s <- p$T_end_h * 3600
  nt <- ceiling(t_end_s / p$dt_s)
  x <- seq(0, p$L_m, length.out = p$nx)
  t <- seq(0, nt) * p$dt_s

  # Initial & boundary
  C <- matrix(p$C_init_mgL, nrow = nt + 1, ncol = p$nx)
  C[1, 1] <- p$C_left_mgL
  C[1, p$nx] <- C[1, p$nx - 1]

  # Diagnostics
  diag <- diagnostics(p, dx, p$dt_s)

  # Time stepping
  for (n in 1:nt) {
    Cn <- C[n, ]
    Cn[1] <- p$C_left_mgL
    Cn[p$nx] <- Cn[p$nx - 1]

    # Upwind advection
    adv <- numeric(p$nx)
    if (p$v_ms >= 0) {
      adv[2:(p$nx - 1)] <- - (p$v_ms * p$dt_s / dx) * (Cn[2:(p$nx - 1)] - Cn[1:(p$nx - 2)])
    } else {
      adv[2:(p$nx - 1)] <- - (abs(p$v_ms) * p$dt_s / dx) * (Cn[3:p$nx] - Cn[2:(p$nx - 1)])
    }

    # FTCS diffusion and explicit decay
    rD <- p$D_m2s * p$dt_s / dx^2
    diff <- numeric(p$nx)
    diff[2:(p$nx - 1)] <- rD * (Cn[3:p$nx] - 2 * Cn[2:(p$nx - 1)] + Cn[1:(p$nx - 2)])
    decay <- - p$k_1s * p$dt_s * Cn

    C_new <- Cn + adv + diff + decay
    C_new[1] <- p$C_left_mgL
    C_new[p$nx] <- C_new[p$nx - 1]
    if (p$clip_nonneg) C_new <- pmax(C_new, 0)
    C[n + 1, ] <- C_new
  }

  # Visualization (≤ 20 lines; correct legends)
  if (isTRUE(p$make_plots)) {
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    if (names(dev.cur()) == "null device" && interactive()) dev.new()
    par(mfrow = c(1, 2), mar = c(4.2, 4.8, 2, 1))

    times_h <- p$times_h[p$times_h <= max(t) / 3600]
    ti <- sapply(times_h, function(h) which.min(abs(t - h * 3600)))
    stations_km <- p$stations_km[p$stations_km <= p$L_m / 1000]
    si <- sapply(stations_km * 1000, function(m) which.min(abs(x - m)))

    cols1 <- seq_along(ti); cols2 <- seq_along(si)

    matplot(x / 1000, t(C[ti, , drop = FALSE]), type = "l", lty = 1, lwd = 2, col = cols1,
            xlab = "x (km)", ylab = "C (mg/L)", main = "Profiles")
    legend("topright", legend = paste0(times_h, " h"), col = cols1, lty = 1, lwd = 2,
           bty = "n", title = "t =")

    matplot(t / 3600, C[, si, drop = FALSE], type = "l", lty = 1, lwd = 2, col = cols2,
            xlab = "time (h)", ylab = "C (mg/L)", main = "Stations")
    legend("right", legend = paste0(stations_km, " km"), col = cols2, lty = 1, lwd = 2,
           bty = "n", title = "x =")
  }

  invisible(list(x_m = x, t_s = t, C_mgL = C, params = p, diagnostics = diag))
}
out <- run_ade()
# ----- Quick start (uncomment to run) -----
# ade_assumptions()
# out <- run_ade()
# str(out)  # explore outputs
