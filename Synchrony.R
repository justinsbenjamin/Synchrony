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

# Read in data files
nest_data <- read_excel("GRD_data_2025.xlsx")
chick_data <- read_excel("Chick_GRD_2025.xlsx")
View(nest_data)
View(chick_data)

# Filter nests that aren't excluded
nest_data_filtered <- nest_data %>% filter(is.na(Exclusion) | Exclusion == "") %>%
  mutate(unique_nest_ID = paste(Year, Nest_ID, sep = "_"))
View(nest_data_filtered)

chick_data_filtered <- chick_data %>% filter(is.na(Exclusion) | Exclusion == "") %>%
  mutate(unique_nest_ID = paste(Year, Nest_ID, sep = "_"))
View(chick_data_filtered)

# Calculating hatch spread
hatch_spread <- nest_data_filtered %>%
  filter(Hatch_begin != "NA" & Hatch_end != "NA") %>%
  mutate(
    Hatch_begin = as.numeric(Hatch_begin),
    Hatch_end = as.numeric(Hatch_end),
    Hatch_spread = Hatch_end - Hatch_begin + 1)

# Separating single and joint female nests
single_female <- nest_data_filtered %>% filter(Clutch_size <= 5)
joint_female <- nest_data_filtered %>% filter(Clutch_size >= 6)

single_hatch_spread <- hatch_spread %>% filter(Clutch_size <= 5)
joint_hatch_spread <- hatch_spread %>% filter(Clutch_size >= 6)

summary(single_hatch_spread$Hatch_spread)
summary(joint_hatch_spread$Hatch_spread)






# Assuming df1 has the clutch sizes and df2 has the nest information
chick_data_w_clutch <- chick_data_filtered %>%
  left_join(nest_data_filtered, by = "Nest_ID")

# Remove duplicates in nest_data_filtered
nest_data_filtered_unique <- nest_data_filtered %>%
  distinct(Nest_ID, .keep_all = TRUE)

# Perform the join
chick_data_w_clutch <- chick_data_filtered %>%
  left_join(nest_data_filtered_unique, by = "Nest_ID")




# Example dataset (each row is an egg)
hatching_data <- data.frame(
  nest_id = c(1,1,1,2,2,2,2,3,3,3),  # Nest IDs
  hatching_day = c(1,2,3,1,1,2,3,1,2,2)  # Hatching days
)

# Calculate clutch size per nest
hatch_summary <- chick_data_w_clutch %>%
  group_by(Nest_ID, Hatch_day) %>%
  summarise(n_hatched = n(), clutch_size = max(Clutch_size), .groups = "drop") %>%
  mutate(prop_hatched = n_hatched / clutch_size)

# Aggregate proportions across all nests
hatch_distribution <- hatch_summary %>%
  group_by(Hatch_day) %>%
  summarise(avg_prop_hatched = mean(prop_hatched), .groups = "drop")

# Plot
ggplot(hatch_distribution, aes(x = Hatch_day, y = avg_prop_hatched)) +
  geom_col(fill = "skyblue", color = "black") +
  labs(x = "Hatching Day",
       y = "Mean Clutch Proportion") +
  theme_classic()


ggplot(hatch_distribution, aes(x = Hatch_day, y = avg_prop_hatched)) +
  geom_col(fill = "skyblue", color = "black") +
  labs(x = "Hatching Day",
       y = "Mean Clutch Proportion") +
  scale_x_continuous(limits = c(0, 15))
  theme_classic()



# Example data set
hatching_data <- data.frame(
  year = c(2023, 2023, 2023, 2024, 2024, 2024, 2024),
  nest_ID = c("A", "B", "C", "A", "B", "C", "D"),
  hatching_day = c(1, 2, 3, 1, 1, 2, 3))










