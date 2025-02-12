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




# Example data set
hatching_data <- data.frame(
  year = c(2023, 2023, 2023, 2024, 2024, 2024, 2024),
  nest_ID = c("A", "B", "C", "A", "B", "C", "D"),
  hatching_day = c(1, 2, 3, 1, 1, 2, 3))

# Create a unique nest identifier
hatching_data <- hatching_data %>%
  mutate(unique_nest_ID = paste(year, nest_ID, sep = "_"))

print(hatching_data)




# Nests with 6 or more eggs considered a joint clutch
# Nests with 5 or less eggs considered a single female clutch
# Do one histogram of clutch proportion for each type of nest

# Do mean and median hatch spread for single female and joint clutches 

# Do mean and median hatching success for single and joint clutches 



