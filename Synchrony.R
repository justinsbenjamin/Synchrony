# Justin Benjamin
# Synchrony Project

# I am investigting hatching synchrony in Pūkeko 
# (Porphyrio melanotus melanotus) as a possible evolutionary route 
# to joint laying behaviour. I will be exploring the hatching spread, 
# hatching rate, and survival of single female and joint female nests 
# and looking at the distribution of hatching. 

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
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>%
  mutate(Females = ifelse(Clutch_size <= 5, "Single Female", "Joint Females"))
View(nest_data_filtered)

chick_data_filtered <- chick_data %>% filter(is.na(Exclusion) | Exclusion == "") %>%
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>%
  left_join(nest_data_filtered %>% select(Nest_ID, Clutch_size, Hatched_eggs), by = "Nest_ID") %>%
  mutate(Females = ifelse(Clutch_size <= 5, "Single Female", "Joint Females")) %>%
  filter(!is.na(Females))
View(chick_data_filtered)

# Calculating hatch spread
hatch_spread <- nest_data_filtered %>%
  filter(Hatch_begin != "NA" & Hatch_end != "NA") %>%
  mutate(
    Hatch_begin = as.numeric(Hatch_begin),
    Hatch_end = as.numeric(Hatch_end),
    Hatch_spread = Hatch_end - Hatch_begin + 1)


data_proportions <- chick_data_filtered %>%
  mutate(Hatch_proportion = Hatched_eggs / Clutch_size) %>%
  filter(!is.na(Hatch_proportion)) %>% # Remove NAs if needed
  group_by(Nest_ID)

data_filtered <- chick_data_filtered %>%
  filter(!(Nest_ID %in% data_proportions))

ggplot(data_proportions, aes(x = Hatch_day, weight = Hatch_proportion)) +
  geom_histogram(fill = "skyblue", colour = "black", binwidth = 1) +
  labs(x = "Relative Hatch Day", y = "Proportion of Eggs Hatched") +
  facet_wrap(~ Females)
  scale_y_continuous(breaks = seq(0, 120, by = 10)) +
  theme_classic()



# Plot
ggplot(chick_data_filtered, aes(x = Hatch_day, weight = Hatch_day)) +
  geom_histogram(fill = "skyblue", colour = "black", binwidth = 1) +
  facet_wrap(~ Females) +
  labs(x = "Hatching Day",
       y = "Chicks hatched") +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  scale_x_continuous(breaks = seq(1, 12, by = 1)) +
  theme_classic() +
  theme(panel.spacing = unit(1, "lines"))





hatch_spread <- nest_data_filtered %>%
  filter(Hatch_begin != "NA" & Hatch_end != "NA") %>%
  mutate(
    Hatch_begin = as.numeric(Hatch_begin),
    Hatch_end = as.numeric(Hatch_end),
    Hatch_spread = Hatch_end - Hatch_begin + 1)

data_proportions <- chick_data_filtered %>%
  mutate(Hatch_proportion = Hatched_eggs / Clutch_size) %>%
  filter(!is.na(Hatch_proportion)) %>% # Remove NAs if needed
  group_by(Nest_ID)

data_filtered <- chick_data_filtered %>%
  filter(!(Nest_ID %in% data_proportions))

ggplot(data_filtered, aes(x = Hatch_day, y = data_proportions), fill = factor(Hatch_day)) +
  geom_violin(alpha = 0.5) +  
  geom_quasirandom(width = 0.2, size = 2, alpha = 0.7, color = "black", method = "smiley") + 
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +  
  labs(x = "Hatching Day",
       y = "Percentage of Clutch Hatched",
       title = "Distribution of % Clutch Hatched Per Day") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14))







library(ggplot2)
library(ggbeeswarm)  # For better point jittering

# Calculate hatch spread
hatch_spread <- nest_data_filtered %>%
  filter(Hatch_begin != "NA" & Hatch_end != "NA") %>%
  mutate(
    Hatch_begin = as.numeric(Hatch_begin),
    Hatch_end = as.numeric(Hatch_end),
    Hatch_spread = Hatch_end - Hatch_begin + 1
  )

# Compute proportion of clutch hatched per day
data_proportions <- chick_data_filtered %>%
  mutate(Hatch_proportion = (1 / Clutch_size) * 100) %>%  # Convert to percentage
  filter(!is.na(Hatch_proportion))  # Remove NAs if needed

# Plot violin graph
ggplot(data_proportions, aes(x = as.factor(Hatch_day), y = Hatch_proportion)) +
  geom_violin(fill = "lightblue", alpha = 0.5) +  
  geom_quasirandom(width = 0.2, size = 2, alpha = 0.7, color = "black", method = "smiley") +  
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 60)) +  
  scale_x_discrete(drop = FALSE) + 
  labs(x = "Hatching Day",
       y = "Percentage of Clutch Hatched",
       title = "Distribution of % Clutch Hatched Per Day") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14))

# Compute proportion of clutch hatched per day
data_proportions <- chick_data_filtered %>%
  mutate(Hatch_proportion = (1 / Hatched_eggs) * 100) %>%  # Convert to percentage
  filter(!is.na(Hatch_proportion))  # Remove NAs if needed


ggplot(data_proportions, aes(x = as.factor(Hatch_day), y = Hatch_proportion)) +
  geom_violin(fill = "lightblue", alpha = 0.5) +  
  geom_quasirandom(width = 0.2, size = 2, alpha = 0.7, color = "black", method = "smiley") +  
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 60)) +  
  scale_x_discrete(drop = FALSE) + 
  labs(x = "Hatching Day",
       y = "Percentage of Clutch Hatched",
       title = "Distribution of % Clutch Hatched Per Day") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14))














library(ggplot2)
library(dplyr)

# Create the example dataset
hatch_data <- data.frame(
  Nest_ID = rep(c("Nest_1", "Nest_2", "Nest_3"), c(4, 4, 6)),  # Assign nests
  Clutch_size = rep(c(8, 5, 10), c(4, 4, 6)),  # Repeat clutch sizes accordingly
  Hatched_eggs = rep(c(4, 4, 6), c(4, 4, 6)),  # Repeat hatched eggs accordingly
  Relative_hatch_day = c(1, 1, 1, 3,  # Nest 1
                         1, 1, 1, 1,  # Nest 2
                         1, 1, 2, 2, 5, 7)  # Nest 3
)

# View dataset
print(hatch_data)

ggplot(hatch_data, aes(x = Nest_ID, y = Relative_hatch_day, fill = Nest_ID)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  # Violin plot to show density
  geom_jitter(width = 0.1, size = 2, alpha = 0.8) +  # Add individual points for better visualization
  scale_y_continuous(breaks = seq(1, 7, by = 1)) +  # Ensure all days appear on the y-axis
  labs(x = "Nest", y = "Relative Hatch Day",
       title = "Distribution of Hatching Days Across Nests") +
  theme_minimal() +
  theme(legend.position = "none",  # Remove legend since Nest_ID is already on x-axis
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14))













library(ggplot2)
library(dplyr)

# Example dataset
hatch_data <- data.frame(
  Nest_ID = rep(c("Nest_1", "Nest_2", "Nest_3"), c(4, 4, 6)),
  Clutch_size = rep(c(8, 5, 10), c(4, 4, 6)),
  Hatched_eggs = rep(c(4, 4, 6), c(4, 4, 6)),
  Relative_hatch_day = c(1, 1, 1, 3,  # Nest 1
                         1, 1, 1, 1,  # Nest 2
                         1, 1, 2, 2, 5, 7)  # Nest 3
)

# Calculate proportions for each nest and hatch day
hatch_summary <- hatch_data %>%
  group_by(Nest_ID, Relative_hatch_day) %>%
  summarise(Hatched_count = n(), .groups = "drop") %>%
  left_join(hatch_data %>% select(Nest_ID, Clutch_size) %>% distinct(), by = "Nest_ID") %>%
  mutate(Proportion = (Hatched_count / Clutch_size) * 100)  # Convert to percentage

# View summary dataset
print(hatch_summary)

ggplot(hatch_summary, aes(x = factor(Relative_hatch_day), y = Proportion, fill = factor(Relative_hatch_day))) +
  geom_violin(alpha = 0.5) +  # Violin plot for density
  geom_jitter(width = 0.1, size = 2, alpha = 0.8) +  # Add points for individual nests
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +  # Format y-axis as percentages
  labs(x = "Hatching Day",
       y = "Percentage of Eggs Hatched",
       title = "Distribution of Hatching Proportions by Day") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14))















library(ggplot2)
library(dplyr)

# Example dataset
hatch_data <- data.frame(
  Nest_ID = rep(c("Nest_1", "Nest_2", "Nest_3"), c(4, 4, 6)),
  Clutch_size = rep(c(8, 5, 10), c(4, 4, 6)),
  Hatched_eggs = rep(c(4, 4, 6), c(4, 4, 6)),
  Relative_hatch_day = c(1, 1, 1, 3,  # Nest 1
                         1, 1, 1, 1,  # Nest 2
                         1, 1, 2, 2, 5, 7)  # Nest 3
)

# Calculate number of eggs hatched per day per nest
hatch_summary <- hatch_data %>%
  group_by(Nest_ID, Relative_hatch_day) %>%
  summarise(Hatched_count = n(), .groups = "drop") %>%
  left_join(hatch_data %>% select(Nest_ID, Clutch_size, Hatched_eggs) %>% distinct(), by = "Nest_ID") %>%
  mutate(
    `%_of_Clutch` = (Hatched_count / Clutch_size) * 100,   # % of total clutch size
    `%_of_Hatched` = (Hatched_count / Hatched_eggs) * 100  # % of total hatched eggs in that nest
  )

# View dataset
print(hatch_summary)


ggplot(hatch_summary, aes(x = factor(Relative_hatch_day), y = `%_of_Clutch`, fill = factor(Relative_hatch_day))) +
  geom_violin(alpha = 0.5) +  
  geom_jitter(width = 0.1, size = 2, alpha = 0.8) +  
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +  
  labs(x = "Hatching Day",
       y = "Percentage of Clutch Hatched",
       title = "Distribution of % Clutch Hatched Per Day") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14))



ggplot(hatch_summary, aes(x = factor(Relative_hatch_day), y = `%_of_Hatched`, fill = factor(Relative_hatch_day))) +
  geom_violin(alpha = 0.5) +  
  geom_jitter(width = 0.1, size = 2, alpha = 0.8) +  
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +  
  labs(x = "Hatching Day",
       y = "Percentage of Hatched Eggs",
       title = "Distribution of % Hatched Eggs Per Day") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14))






library(ggplot2)
library(dplyr)

set.seed(123)  # For reproducibility

# Generate random nest data
nests <- data.frame(
  Nest_ID = paste0("Nest_", 1:10),
  Clutch_size = sample(4:10, 10, replace = TRUE),
  Hatched_eggs = sample(2:8, 10, replace = TRUE)
)

# Simulate relative hatch days per nest
hatch_data <- nests %>%
  rowwise() %>%
  mutate(
    Relative_hatch_days = list(sort(sample(1:7, Hatched_eggs, replace = TRUE))) # Hatch over 1-7 days
  ) %>%
  unnest(cols = c(Relative_hatch_days)) %>%
  rename(Relative_hatch_day = Relative_hatch_days) %>%
  ungroup()

# Count hatched eggs per day per nest
hatch_summary <- hatch_data %>%
  group_by(Nest_ID, Relative_hatch_day) %>%
  summarise(Hatched_count = n(), .groups = "drop") %>%
  left_join(nests, by = "Nest_ID") %>%
  mutate(
    `%_of_Clutch` = (Hatched_count / Clutch_size) * 100,   # % of total clutch size
    `%_of_Hatched` = (Hatched_count / Hatched_eggs) * 100  # % of total hatched eggs in that nest
  )

# View dataset
print(hatch_summary)


ggplot(hatch_summary, aes(x = factor(Relative_hatch_day), y = `%_of_Clutch`, fill = factor(Relative_hatch_day))) +
  geom_violin(alpha = 0.5) +  
  geom_quasirandom(width = 0.2, size = 2, alpha = 0.7, color = "black", method = "smiley") + 
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +  
  labs(x = "Hatching Day",
       y = "Percentage of Clutch Hatched",
       title = "Distribution of % Clutch Hatched Per Day") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14))


ggplot(hatch_data, aes(x = as.factor(Relative_hatch_day), y = `%_of_Clutch`)) +
  geom_violin(fill = "lightblue", alpha = 0.5) +
  geom_quasirandom(width = 0.2, size = 2, alpha = 0.7, color = "black", method = "smiley") +
  labs(x = "Hatch Day", y = "Hatch Proportion (0-100%)") +
  theme_classic()

