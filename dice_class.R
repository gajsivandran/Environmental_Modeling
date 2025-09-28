# Load required package
library(ggplot2)

# Set seed for reproducibility (optional)
set.seed(123)

# Repeat the experiment 1000 times
results <- replicate(1000000, {
  # Generate 2 random integers between 1 and 6
  dice <- sample(1:6, 5, replace = TRUE)

  # Display each roll (commented out to avoid huge output, uncomment if needed)
  # print(dice)

  # Add them
  sum(dice)
})

# Convert to data frame for plotting
results_df <- data.frame(Sum = results)

# Plot frequency of numbers
ggplot(results_df, aes(x = factor(Sum))) +
  geom_bar(fill = "skyblue", color = "black") +
  labs(
    title = "Frequency of Dice Sums (1000 Trials)",
    x = "Sum of Two Dice",
    y = "Frequency"
  ) +
  theme_minimal()
