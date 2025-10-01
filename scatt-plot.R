# Create data
x <- 1:10
y <- x^2

# Scatter plot with upside-down triangles
plot(x, y,
     main = "Simple Plot",
     xlab = "x",
     ylab = "y^2",
     pch = 25,        # upside-down triangle
     col = "green",     # border color
     bg = "yellow",      # fill color
     cex = 1.5        # point size scaling
)

# Add grid lines
grid




library(ggplot2)

# Create data
df <- data.frame(x = 1:10, y = (1:10)^2)

ggplot(df, aes(x, y)) +
  geom_point(shape = 25, color = "green", fill = "yellow", size = 4) +
  labs(
    title = "Simple Plot",
    x = "x",
    y = "y^2"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey80", size = 0.5),
    panel.grid.minor = element_line(color = "grey90", size = 0.25)
  )

library(ggplot2)

# Create data
df <- data.frame(x = 1:10, y = (1:10)^2)

ggplot(df, aes(x, y)) +
  geom_point(shape = 25, color = "red", fill = "red", size = 4) +
  labs(
    title = "Simple Plot",
    x = "x",
    y = "y^2"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey80", size = 0.5),
    panel.grid.minor = element_line(color = "grey90", size = 0.25)
  )

