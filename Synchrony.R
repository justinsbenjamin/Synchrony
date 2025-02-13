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
















