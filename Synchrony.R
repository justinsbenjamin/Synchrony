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
library(MASS)
library(performance)
library(tidyr)
library(readxl)
library(lubridate)
library(ggbeeswarm)

# Read in data files
nest_data <- read_excel("GRD_data_2025.xlsx")
chick_data <- read_excel("Chick_GRD_2025.xlsx")
View(nest_data)
View(chick_data)

# Filter nests that aren't excluded
# Create unique nest IDs
# Add single female and joint female categories
# Calculate hatch spread
nest_data_filtered <- nest_data %>% filter(is.na(Exclusion) | Exclusion == "") %>%
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>%
  mutate(Females = ifelse(Clutch_size <= 5, "Single Female", "Joint Females")) %>%
  filter(Hatch_begin != "NA" & Hatch_end != "NA") %>%
  mutate(
    Hatch_begin = as.numeric(Hatch_begin),
    Hatch_end = as.numeric(Hatch_end),
    Hatch_spread = Hatch_end - Hatch_begin + 1)
View(nest_data_filtered)

# Filter nests that aren't excluded
# Create unique nest IDs
# Merge clutch size and hatched eggs from nest data set to chick data set
# Add single female and joint female categories
chick_data_filtered <- chick_data %>% filter(is.na(Exclusion) | Exclusion == "") %>%
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>%
  left_join(nest_data_filtered %>% dplyr::select(Nest_ID, Clutch_size, Hatched_eggs), by = "Nest_ID") %>%
  mutate(Females = ifelse(Clutch_size <= 5, "Single Female", "Joint Females")) %>%
  filter(!is.na(Females)) %>%
  mutate(Hatch_rate = Hatched_eggs / Clutch_size) %>%
  filter(!is.na(Hatch_rate)) %>% # Remove NAs if needed
  group_by(Nest_ID)
View(chick_data_filtered)


# Histogram showing counts of hatching based on relative hatch days
# separated between single female and joint female nests. 
chick_data_filtered <- chick_data_filtered %>% filter(!is.na(Hatch_day))
ggplot(chick_data_filtered, aes(x = Hatch_day, weight = Hatch_day)) +
  geom_histogram(fill = "skyblue", colour = "black", binwidth = 1) +
  facet_wrap(~ Females, scales = "free_x") +
  labs(x = "Relative hatching day",
       y = "Chicks hatched") +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  scale_x_continuous(breaks = seq(1, 11, by = 1)) +
  theme_classic() +
  theme(panel.spacing = unit(1, "lines"))




ggplot(chick_data_filtered, aes(x = Hatch_day, weight = Hatch_rate)) +
  geom_histogram(fill = "skyblue", colour = "black", binwidth = 1) +
  labs(x = "Relative Hatch Day", y = "Proportion of Eggs Hatched") +
  facet_wrap(~ Females)
scale_y_continuous(breaks = seq(0, 120, by = 10)) +
  theme_classic()


data_proportions <- chick_data_filtered %>%
  group_by(Nest_ID) %>%
  mutate(Hatch_proportion = (n() / Hatched_eggs) * 100) %>%  # Compute correct proportion
  filter(!is.na(Hatch_proportion))  # Remove NAs if needed


data_proportions <- chick_data_filtered %>%
  mutate(Hatch_proportion = Hatched_eggs / Clutch_size) %>%
  filter(!is.na(Hatch_proportion)) %>% # Remove NAs if needed
  group_by(Nest_ID)

data_filtered <- chick_data_filtered %>%
  filter(!(Nest_ID %in% data_proportions))

ggplot(data_filtered, aes(x = Hatch_day, y = data_proportions), 
       fill = factor(Hatch_day)) +
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



# Compute proportion of clutch hatched per day

data_proportions <- chick_data_filtered %>%
  group_by(Nest_ID) %>%
  mutate(Hatch_proportion = (n() / Hatched_eggs) * 100) %>%  # Compute correct proportion
  filter(!is.na(Hatch_proportion))  # Remove NAs if needed

# Plot violin graph
ggplot(data_proportions, aes(x = as.factor(Hatch_day), y = Hatch_proportion)) +
  geom_violin(fill = "lightblue", alpha = 0.5) +  
  geom_jitter(width = 0.1, size = 2, alpha = 0.8) + 
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +  
  scale_x_discrete(drop = FALSE) + 
  labs(x = "Relative hatching day",
       y = "Percentage of clutch hatched") +
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
  geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
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



# Create the example dataset
hatch_data <- data.frame(
  Nest_ID = rep(c("Nest_1", "Nest_2", "Nest_3"), c(4, 4, 6)),  # Assign nests
  Clutch_size = rep(c(8, 5, 10), c(4, 4, 6)),  # Repeat clutch sizes accordingly
  Hatched_eggs = rep(c(4, 4, 6), c(4, 4, 6)),  # Repeat hatched eggs accordingly
  Relative_hatch_day = c(1, 1, 1, 3,  # Nest 1
                         1, 1, 1, 1,  # Nest 2
                         1, 1, 2, 2, 5, 7)  # Nest 3
)



###############################
##### This is the good figure!!!!
###############################

data_proportions <- chick_data_filtered %>%
  group_by(Nest_ID, Hatch_day) %>%
  summarise(Chicks_hatched = n(), .groups = "drop") %>%  # Count how many chicks hatched per day
  left_join(nest_data_filtered %>% dplyr::select(Nest_ID, Hatched_eggs, Females), by = "Nest_ID") %>%
  mutate(Hatch_proportion = (Chicks_hatched / Hatched_eggs) * 100)  # Compute correct percentage

# Plot violin graph
ggplot(data_proportions, aes(x = as.factor(Hatch_day), y = Hatch_proportion, group = Nest_ID)) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), 
              size = 1, alpha = 0.8) + 
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100), 
                     breaks = seq(0, 100, by = 10)) +  
  labs(x = "Relative hatching day",
       y = "Percentage of hatched chicks per nest") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14)) +
  facet_wrap(~ Females, scales = "free_x")

scale_x_(drop = FALSE) + 




















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



