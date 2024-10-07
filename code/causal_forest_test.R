library(ggplot2)
library(gridExtra)
library(grf)

# Generate new data
set.seed(1)
n <- 10000
p <- 20
X <- matrix(rnorm(n * p), n, p)

# Modified TAU function to include X1
TAU <- function(x1, x3) {
  1 / (1 + exp(-x3)) + 0.5 * (x1 > 0)
}

TAU_values <- TAU(X[, 1], X[, 3])
W <- rbinom(n, 1, 1 / (1 + exp(-X[, 1] - X[, 2])))
Y <- pmax(X[, 2] + X[, 3], 0) + rowMeans(X[, 4:6]) / 2 + W * TAU_values + rnorm(n)

# Fit models
lm_model <- lm(Y ~ . + W:X1 + W:X3, data = data.frame(X, W = W, Y = Y))
cf_model <- causal_forest(X, Y, W, tune.parameters = "all")

# Create dataframe for plotting
plot_data <- data.frame(
  X1 = X[, 1],
  X3 = X[, 3],
  tau_true = TAU_values,
  tau_lm = predict(lm_model, newdata = data.frame(X, W = 1)) - 
    predict(lm_model, newdata = data.frame(X, W = 0)),
  tau_cf = predict(cf_model)$predictions
)

# Determine global min and max for color scale
global_min <- min(plot_data$tau_true, plot_data$tau_lm, plot_data$tau_cf)
global_max <- max(plot_data$tau_true, plot_data$tau_lm, plot_data$tau_cf)

# Plotting function
create_plot <- function(data, x, y, z, title) {
  ggplot(data, aes(x = .data[[x]], y = .data[[y]], z = .data[[z]])) +
    stat_summary_2d(fun = mean, bins = 30) +
    scale_fill_viridis_c(option = "plasma", limits = c(global_min, global_max)) +
    labs(title = title, x = x, y = y, fill = "TE") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, hjust = 0.5),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 6),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.5, "cm"))
}

# Create plots
p_true <- create_plot(plot_data, "X1", "X3", "tau_true", "True Treatment Effect (TE)")
p_lm <- create_plot(plot_data, "X1", "X3", "tau_lm", "Predicted TE by LM")
p_cf <- create_plot(plot_data, "X1", "X3", "tau_cf", "Predicted TE by CF")

# Combine plots
combined_plot <- grid.arrange(p_true, p_lm, p_cf, nrow = 1)

# Display the plot
print(combined_plot)
