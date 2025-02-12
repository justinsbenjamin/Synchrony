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
library(readxl)
library(lubridate)

# Read in data file
data <- read_excel("GRD_data_2025.xlsx")
View(data)

data_filtered <- data %>% filter(is.na(Exclusion) | Exclusion == "") %>%
  mutate(unique_nest_ID = paste(Year, Nest_ID, sep = "_"))
View(data_filtered)

hatch_spread <- data_filtered %>%
  filter(Hatch_begin != "NA" & Hatch_end != "NA") %>%
  mutate(
    Hatch_begin = as.numeric(Hatch_begin),
    Hatch_end = as.numeric(Hatch_end),
    Hatch_spread = Hatch_end - Hatch_begin + 1)

single_female <- data_filtered %>% filter(Clutch_size <= 5)
joint_female <- data_filtered %>% filter(Clutch_size >= 6)

single_hatch_spread <- hatch_spread %>% filter(Clutch_size <= 5)
joint_hatch_spread <- hatch_spread %>% filter(Clutch_size >= 6)

summary(single_hatch_spread$Hatch_spread)
summary(joint_hatch_spread$Hatch_spread)


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
data_filtered <- data_filtered %>%
  mutate(unique_nest_ID = paste(year, nest_ID, sep = "_"))

print(hatching_data)


# Do one histogram of clutch proportion for each type of nest





