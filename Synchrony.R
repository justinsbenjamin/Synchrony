# Justin Benjamin
# Synchrony Project


# I am investigting hatching synchrony in Pūkeko (Porphyrio melanotus melanotus) 
# as a possible evolutionary route to joint laying behaviour. 

# Load libraries
library(dplyr)
library(ggplot2)
library(effects)
library(glmmTMB)
library(MASS)
library(performance)
library(tidyr)


read_xlsx()


library(ggplot2)

# Example data: Hatching days for individual eggs
hatching_data <- data.frame(
  hatching_day = c(1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 5, 1, 2, 4)  # Replace with your actual data
)

# Count occurrences of each hatching day
hatch_counts <- table(hatching_data$hatching_day)

# Convert to proportions
hatch_proportions <- hatch_counts / sum(hatch_counts)

# Convert to data frame for ggplot
hatch_df <- data.frame(
  hatching_day = as.numeric(names(hatch_counts)),
  proportion = as.numeric(hatch_proportions)
)

# Plot histogram with proportions
ggplot(hatch_df, aes(x = hatching_day, y = proportion)) +
  geom_col(fill = "skyblue", color = "black") +
  labs(title = "Proportion of Eggs Hatching by Day",
       x = "Hatching Day",
       y = "Proportion of Total Eggs") +
  theme_classic()




library(dplyr)
library(ggplot2)

# Example dataset (each row is an egg)
hatching_data <- data.frame(
  nest_id = c(1,1,1,2,2,2,2,3,3,3),  # Nest IDs
  hatching_day = c(1,2,3,1,1,2,3,1,2,2)  # Hatching days
)

# Calculate clutch size per nest
hatch_summary <- hatching_data %>%
  group_by(nest_id) %>%
  mutate(clutch_size = n()) %>%
  group_by(nest_id, hatching_day) %>%
  summarise(n_hatched = n(), clutch_size = max(clutch_size), .groups = "drop") %>%
  mutate(prop_hatched = n_hatched / clutch_size)

# Aggregate proportions across all nests
hatch_distribution <- hatch_summary %>%
  group_by(hatching_day) %>%
  summarise(avg_prop_hatched = mean(prop_hatched), .groups = "drop")

# Plot
ggplot(hatch_distribution, aes(x = hatching_day, y = avg_prop_hatched)) +
  geom_col(fill = "skyblue", color = "black") +
  labs(x = "Hatching Day",
       y = "Mean Clutch Proportion") +
  theme_classic()

