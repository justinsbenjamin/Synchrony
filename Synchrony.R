# Justin Benjamin
# Synchrony Project

# I am investigating hatching synchrony in Pūkeko (Porphyrio melanotus melanotus) 
# as a possible evolutionary route to joint laying behaviour in this species. 
# I hypothesize that jointly laying synchronously with a close relative increases
# the inclusive fitness of the dominant female compared to hypothetical 
# prolonged laying by a single female. 
# I will be exploring the estimated and true hatching spread, hatching rate, 
# and survival of single female and joint female nests and looking at 
# the distribution of hatching. I am interested to see if factors like 
# clutch size, group size, clutch number, joint/single female clutches, and level 
# of inbreeding have impacts on hatching synchrony and reproductive success. 

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
  left_join(nest_data_filtered %>% dplyr::select(Nest_ID, Clutch_size, Hatched_eggs, Hatch_spread), by = "Nest_ID") %>%
  mutate(Females = ifelse(Clutch_size <= 5, "Single Female", "Joint Females")) %>%
  filter(!is.na(Females)) %>%
  mutate(Hatch_rate = Hatched_eggs / Clutch_size) %>%
  filter(!is.na(Hatch_rate)) %>% # Remove NAs if needed
  group_by(Nest_ID) %>% 
  mutate(Synchrony = case_when(
    Hatch_spread >= 1 & Hatch_spread <= 5 ~ "Synchronous",
    Hatch_spread >= 6 ~ "Asynchronous",
    TRUE ~ NA_character_  # Handles missing or unexpected values
  ))
View(chick_data_filtered)


single_female <- chick_data_filtered %>% 
  filter(Females == "Single Female")
summary(single_female)

joint_female <- chick_data_filtered %>% 
  filter(Females == "Joint Females")
summary(joint_female)


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



###############################
##### This is the good figure!!!! ####
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



experiment_data <- read_excel("Synchrony_experiment_GRD.xlsx")
View(experiment_data)

ggplot(experiment_data, aes(x = Treatment_Type, fill = factor(first_sur))) +
  geom_bar(position = "dodge") +
  labs(x = "Treatment Type", y = "Count", fill = "Survival") +
  scale_fill_manual(values = c("0" = "Black", "1" = "Grey")) +  
  theme_classic()
  

chick_data_filtered %>% 
  filter(Hatched_eggs > 1) 

ggplot(chick_data_filtered, aes(x = Hatch_spread, y = Hatch_rate, fill = Hatch_spread)) +
  geom_jitter(width = 0.2, alpha = 0.5) +  # Adds points to show spread
  labs(x = "Relative Hatch Day", y = "Hatching Success") +
  scale_x_continuous(breaks = seq(1, 11, by = 1)) +
  theme_classic() +
  facet_wrap(~ Females) 



# Make a scatter plot with hatching spread on the x axis and hatching success
# on the y axis. See if there is any trends 



# library(glmmTMB)
# library(MuMIn)  # For model selection and AIC ranking

# model_full <- glmmTMB(Survived ~ (NestType + HatchingSpread + ClutchSize + GroupSize + 
                                    # ObservationPeriod + ClutchNumber)^2 + 
                        # (1 | Year/NestID), 
                      # family = binomial, data = df)

# models <- dredge(model_full, rank = "AICc")

# model_list <- list(
  # full = model_full,
  # no_interactions = glmmTMB(Survived ~ NestType + HatchingSpread + ClutchSize + 
                            # GroupSize + ClutchNumber + ObservationPeriod (1 | Year/NestID), 
                            # family = binomial, data = df),
  # null_model = glmmTMB(Survived ~ (1 | Year/NestID), 
                          # family = binomial, data = df))

# model_selection <- model.sel(model_list)
# print(model_selection)











#### Synchrony experiment stuff ####


# library(glmmTMB)
# library(MuMIn)  # For model selection and AIC ranking

# model_full <- glmmTMB(Survived ~ (Treatment + Estimated_Hatch_Spread + ClutchSize + 
# GroupSize + ObservationPeriod + ClutchNumber + native_foreign eggs)^2 + 
# (1 | Year/NestID), 
# family = binomial, data = df)

# models <- dredge(model_full, rank = "AICc")

# model_list <- list(
# full = model_full,
# no_interactions = glmmTMB(Survived ~ NestType + HatchingSpread + ClutchSize + 
# GroupSize + ClutchNumber + ObservationPeriod (1 | Year/NestID), 
# family = binomial, data = df),
# null_model = glmmTMB(Survived ~ (1 | Year/NestID), 
# family = binomial, data = df))

# model_selection <- model.sel(model_list)
# print(model_selection)



experiment_data <- read_excel("Compiled_experiment_data.xlsx")  %>% 
  mutate(Hatch_success = Hatched/Manipulated_clutch_size) %>%
  mutate(Foreign_percentage = Foreign_eggs/Manipulated_clutch_size) %>%
  filter(Hatch_success != "NA") %>%
  filter(Treatment != "Control") %>% filter(Treatment != "Other")
View(experiment_data)

tapply(experiment_data$Hatch_success, experiment_data$Treatment, summary, na.rm = TRUE)
tapply(experiment_data$Surival_60, experiment_data$Treatment, summary, na.rm = TRUE)






# Boxplot with jittered points

# Hatch success plot
successful_nests <- experiment_data %>% filter(Hatch_success >0)
ggplot(successful_nests, aes(x = Treatment, y = Hatch_success, color = Treatment)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +  
  geom_jitter(position = position_jitter(width = 0.05, height = 0), size = 3, alpha = 0.7) + 
  labs(y = "Hatch rate", x = "Treatment") +
  theme_classic() +
  scale_color_manual(values = c("Asynch" = "darkblue", "Synch" = "firebrick2"))

# Number of swapped (foreign) eggs in nest plot
ggplot(experiment_data, aes(x = Treatment, y = Foreign_percentage, color = Treatment)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +  
  geom_jitter(position = position_jitter(width = 0.05, height = 0), size = 3, alpha = 0.7) +  
  labs(y = "Percentage of foreign eggs in nest", x = "Treatment") +
  theme_classic() +
  scale_color_manual(values = c("Asynch" = "darkblue", "Synch" = "firebrick2"))

# Clutch size plot
ggplot(experiment_data, aes(x = Treatment, y = Manipulated_clutch_size, color = Treatment)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.1, height = 0), size = 3, alpha = 0.7) +
  labs(y = "Clutch size", x = "Treatment") +
  theme_classic() +
  scale_color_manual(values = c("Asynch" = "darkblue", "Synch" = "firebrick2"))



data_expanded <- experiment_data %>%
  mutate(NonSurvivors = Hatched - Surival_60) %>% # Calculate non-survivors
  select(Nest, Treatment, Surival_60, NonSurvivors) %>%
  pivot_longer(cols = c(Surival_60, NonSurvivors), 
                      names_to = "Status", 
                      values_to = "Count") %>%
  mutate(Surival_60 = ifelse(Status == "Surival_60", 1, 0)) %>%
  uncount(Count) %>%
  select(Nest, Treatment, Surival_60)
View(data_expanded)


# Bar plot of survival data
ggplot(data_expanded, aes(x = Treatment, fill = factor(Surival_60))) +
  geom_bar(position = "dodge") +
  labs(x = "Treatment", y = "count", fill = "Surival_60") +
  scale_fill_manual(values = c("0" = "Black", "1" = "Grey")) +  
  theme_classic()
  










(t.test(Hatch_success ~ Treatment, data = experiment_data))



#### POWER ANALYSIS STUFF ####

library(pwr)
library(rstatix)

cohens_d(Hatch_success ~ Treatment, data = experiment_data)
cohens_d(Surival_60 ~ Treatment, data = experiment_data)

# Setting the parameters
sample_size <- 17    # Sample size
alpha <- 0.05        # Significance level
power <- 0.8         # Power
effect_size <- 0.6   # Effect size  

# Solve for effect size needed
effect_size <- pwr.p.test(n = sample_size, sig.level = alpha, power = power, alternative = "two.sided")$h
print(effect_size)

# Solve for sample size needed
sample_size <- pwr.p.test(h = effect_size, sig.level = alpha, power = power, alternative = "two.sided")$n
print(sample_size)

