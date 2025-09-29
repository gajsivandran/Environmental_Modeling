# -------- Scatter plot + single polynomial fit --------

# 1) Define your numbers
x <- 1:10
y <- c(6, 1, 7, 2, 3, 3, 9, 3, 3, 0)

# 2) Plot the scatter (base layer)
plot(x, y,
     pch = 19,
     xlab = "Index (1..10)",
     ylab = "Given Value",
     main = "Scatter with Single Polynomial Fit")

# 3) Ask the user for the polynomial order
get_degree <- function() {
  repeat {
    ans <- readline(prompt = "Enter the polynomial order for the best-fit line: ")
    k <- suppressWarnings(as.integer(ans))
    if (!is.na(k) && k >= 1) return(k)
    cat("Please enter a valid integer >= 1.\n")
  }
}

# Helper to format polynomial equation
poly_equation <- function(model) {
  coefs <- coef(model)
  terms <- c()
  for (i in seq_along(coefs)) {
    pow <- i - 1
    term <- if (pow == 0) {
      sprintf("%.3f", coefs[i])
    } else if (pow == 1) {
      sprintf("%.3f*x", coefs[i])
    } else {
      sprintf("%.3f*x^%d", coefs[i], pow)
    }
    terms <- c(terms, term)
  }
  paste("y =", paste(terms, collapse = " + "))
}

# 4) Fit polynomial
k <- get_degree()
fit <- lm(y ~ poly(x, degree = k, raw = TRUE))

# Print polynomial equation to console
cat("Polynomial fit (order", k, "):\n")
cat(poly_equation(fit), "\n\n")

# Smooth curve
x_new <- seq(min(x), max(x), length.out = 200)
y_hat <- predict(fit, newdata = data.frame(x = x_new))

# Draw line
lines(x_new, y_hat, lwd = 2, col = "blue")

legend("topleft",
       legend = paste("Polynomial fit (order =", k, ")"),
       lwd = 2,
       col = "blue",
       bty = "n")


cat("Finished plotting. You can now close the plot window.\n")
# ------------------------------------------------------------------------
