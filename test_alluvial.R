# Load necessary libraries
library(ggplot2)
library(ggalluvial)
library(dplyr)

# Example dataset
data <- data.frame(
  sample = rep(c("Sample1", "Sample2", "Sample3"), each = 3),
  category = c("A", "B", "C", "A", "B", "D", "A", "C", "D"),
  count = c(10, 20, 15, 5, 25, 10, 8, 12, 7)
)

# Identify shared flows (categories present in multiple samples)
shared_flows <- data %>%
  group_by(category) %>%
  filter(n_distinct(sample) > 1)

# Filter the data to keep only shared flows
filtered_data <- data %>%
  semi_join(shared_flows, by = "category")

# Create the alluvial plot
ggplot(filtered_data,
       aes(axis1 = sample, axis2 = category, y = count)) +
  geom_alluvium(aes(fill = category)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  labs(title = "Alluvial Plot with Shared Flows",
       x = "Sample and Category",
       y = "Count")
