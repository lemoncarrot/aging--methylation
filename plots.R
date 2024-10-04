library(cowplot)
library(tidyverse)
set.seed(123)

generate_data <- function(x_range = 0:20, absolute_change = 0.4, corr_coef = 0.8, base_mean = 0.3, base_sd = 0.01, error_sd = 0.005) {
  x <- x_range
  
  y <- rnorm(length(x), mean = base_mean, sd = base_sd)
  
  linear_component <- seq(from = 0, to = absolute_change, length.out = length(x))
  noise <- rnorm(length(x), mean = 0, sd = sqrt((1 - corr_coef^2) * var(linear_component)))
  
  y <- base_mean + linear_component + noise - (base_mean + linear_component[1])
  
  data <- data.frame(x = x, y = y)
  
  data$color <- ifelse(data$x %in% c(7, 9, 16), "red", "green")
  data$ymin <- data$y - error_sd
  data$ymax <- data$y + error_sd
  
  larger_error_indices <- c(8, 10, 17)
  data$ymin[larger_error_indices] <- data$ymin[larger_error_indices] - 0.15
  data$ymax[larger_error_indices] <- data$ymax[larger_error_indices] + 0.15
  
  return(data)
}

data <- generate_data(x_range = 0:20, absolute_change = 0.4, corr_coef = 0.3)
ggplot(data, aes(x = x, y = y)) +
  #geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, size = 0.8, color = "darkblue") +  # Add error bars
  geom_point(size=6) +  # Plot points
  #scale_color_manual(values = c("black", "green", "red")) +  # Define custom colors
  theme_cowplot(12) +  # Use cowplot theme
  #geom_vline(xintercept = 40, linetype = "dotted", color = "blue", size = 2) +
  #geom_vline(xintercept = 60, linetype = "dotted", color = "blue", size = 2) +
  #annotate("rect", xmin = 40, xmax = 60, ymin = 0, ymax = 1, alpha = 0.2, fill = "darkblue") +
  ylim(0, 1) +  # Set y-axis limits
  geom_smooth(method="lm", se=FALSE, ) +
  theme(legend.position = "none") + 
  labs(x = "Age", y = "Beta Value") +  # Label the axes
  theme(
    legend.position = "none",       
    text = element_text(face = "bold"), 
    axis.title = element_text(face = "bold"), 
    axis.text = element_text(face = "bold"),  
    plot.title = element_text(face = "bold"),  
    plot.subtitle = element_text(face = "bold") 
  ) +
  xlim(0, 20) 
ggsave("badcorr.png", plot=last_plot(), dpi=300, height = 8, width = 5)


# Load necessary libraries
library(ggplot2)
library(cowplot)

# Generate some example data
set.seed(123)
n <- 200
x <- rnorm(n)
y <- rnorm(n)
data <- data.frame(x = x, y = y)

# Apply k-means clustering
set.seed(123)  # For reproducibility
kmeans_result <- kmeans(data, centers = 4)
data$cluster <- as.factor(kmeans_result$cluster)
centers <- as.data.frame(kmeans_result$centers)
centers$cluster <- as.factor(1:4)

# Plot the clustered data
x
ggsave("examplekmeans.png", plot=last_plot(), dpi=300, height = 5, width = 10)



