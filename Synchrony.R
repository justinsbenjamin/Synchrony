# Justin Benjamin
# Synchrony Project

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
library(pwr)
library(rstatix)
library(glmmTMB)
library(effsize)
library(DHARMa)


########### CHAPTER 2 ##########
#### NATURAL SYNCHRONY DATA ####
################################

# I am investigating hatching synchrony in Pūkeko (Porphyrio melanotus melanotus) 
# as a possible evolutionary route to joint laying behaviour in this species. 
# This data has been compiled from field books and spreadsheets to create a 
# master dataset with both summary data for all nests along with specific 
# hatching data for individual nests. I am exploring if factors like 
# clutch size, group size, clutch number, joint/single female clutches, 
# have impacts on hatching synchrony and reproductive success. And also looking
# how hatching synchrony affects hatching success and survival.

# Read in data and filtering/data modifications
masterlist_data <- read_excel("Nests_masterlist.xlsx") %>%
  filter(Exclusion == "GOOD" | Exclusion == "ABCL") %>% # Keeping only the good nests
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>% # Unique nest ID 
  mutate(Clutch_size = as.numeric(Clutch_size)) %>%
  mutate(Females = ifelse(Clutch_size <= 5, "Single Female", "Joint Females")) %>% # assign single/joint female
  filter(Hatch_begin != "NA" & Hatch_end != "NA") %>%
  mutate(Hatch_begin = as.numeric(Hatch_begin),
         Hatch_end = as.numeric(Hatch_end),
         Hatch_spread = Hatch_end - Hatch_begin + 1) %>% # First to last day inclusive
  mutate(Hatch_begin = as.numeric(Hatch_begin), Date_found = as.numeric(Date_found), 
       Observed_incubation_period = as.numeric(Hatch_begin - Date_found)) %>%
  mutate(Hatch_success = as.numeric(Hatched_eggs)/as.numeric(Clutch_size))
View(masterlist_data)

Clutch_size <- (masterlist_data$Clutch_size)
hist(as.numeric(Clutch_size), 
     xlab = "Clutch size",
     xlim = c(0, 20),
     breaks = seq(0,20,1),
     main = "")

Hatch_spread <- (masterlist_data$Hatch_spread)
hist(as.numeric(Hatch_spread), 
     xlab = "Hatch spread (days)",
     xlim = c(0, 15),
     breaks = seq(0,15,1),
     main = "")

ggplot(masterlist_data, aes(x = as.numeric(Clutch_size), y = as.numeric(Hatch_spread), colour = Females)) +
  geom_point() +
  theme_classic()


ggplot(masterlist_data, aes(x = as.numeric(Clutch_size), 
                            y = as.numeric(Hatch_success), colour = Females)) +
  geom_point() +
  theme_classic()



df <- masterlist_data %>%
  filter(!is.na(as.numeric(Survived_2))) %>%
  group_by(Hatched_eggs, Survived_2) %>%
  summarise(count = n(), .groups = "drop")
View(df)

ggplot(df, aes(x = Hatched_eggs, y = Survived_2, size = factor(count))) +
  geom_point(alpha = 0.6) +
  scale_size_manual(values = c("1" = 1, "2" = 2, "3" = 3, "4" = 4, "5" = 6, "7" = 7, "11" = 11)) + 
  theme_classic() +
  labs(y = "Survived", x = "Hatched eggs") 

ggplot(masterlist_data, aes(x = as.numeric(Hatched_eggs), as.numeric(Survived_2))) +
  geom_point() +
  theme_classic()




df <- masterlist_data %>%
  filter(!is.na(as.numeric(Survived_2))) %>%
  group_by(Hatch_spread, Survived_2) %>%
  summarise(count = n(), .groups = "drop")
View(df)

ggplot(df, aes(x = Hatch_spread, y = Survived_2, size = factor(count))) +
  geom_point(alpha = 0.6) +
  scale_size_manual(values = c("1" = 1, "2" = 2, "3" = 3, "4" = 4, "6" = 6, "13" = 10)) + 
  theme_classic() +
  labs(y = "Survived", x = "Hatch spread (days)") 


df <- masterlist_data %>%
  filter(!is.na(as.numeric(Survived_2))) %>%
  group_by(Hatched_eggs, Survived_2) %>%
  summarise(count = n(), .groups = "drop")
View(df)

ggplot(df, aes(x = Hatched_eggs, y = Survived_2, size = factor(count))) +
  geom_point(alpha = 0.6) +
  scale_size_manual(values = c("1" = 1, "2" = 2, "3" = 3, "4" = 4, "5" = 6, "7" = 7, "11" = 11)) + 
  theme_classic() +
  labs(y = "Survived", x = "Hatched eggs") 





ggplot(masterlist_data, aes(x = as.numeric(Hatched_eggs), y = survive_1)) +
  geom_point()


filter(survive_1 != "NA" & Survived_2 != "NA") 



masterlist_data <- read_excel("Nests_masterlist.xlsx")
masterlist_data <- masterlist_data %>% filter(Exclusion == "GOOD" | Exclusion == "ABCL") %>%
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>%
  mutate(Females = ifelse(Clutch_size <= 5, "Single Female", "Joint Females")) %>%
  mutate(
    Hatch_begin = as.numeric(Hatch_begin),
    Hatch_end = as.numeric(Hatch_end),
    Hatch_spread = Hatch_end - Hatch_begin + 1) %>%
  filter(survive_1 != 0 & survive_1 != "NA") %>%
  filter(Survived_2 != "NA" & Survived_2 != "2?")
View(masterlist_data)

ggplot(masterlist_data, aes(x = as.numeric(Hatched_eggs), y = survive_1)) +
  geom_point()

# Survival data conservative estimate with number of banded juveniles as true number

# Survival data with liberal estimate number of chicks that hatched
mutate(Survival_60 = Hatched_eggs)
mutate(Survival_60 = Hatched_eggs)

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


chick_data_filtered %>% 
  filter(Hatched_eggs > 1) 

ggplot(chick_data_filtered, aes(x = Hatch_spread, y = Hatch_rate, fill = Hatch_spread)) +
  geom_jitter(width = 0.2, alpha = 0.5) +  # Adds points to show spread
  labs(x = "Relative Hatch Day", y = "Hatching Success") +
  scale_x_continuous(breaks = seq(1, 11, by = 1)) +
  theme_classic() +
  facet_wrap(~ Females) 



########## CHAPTER 3 #########
#### SYNCHRONY EXPERIMENT ####
##############################


# 
# I hypothesize that jointly laying synchronously with a close relative increases
# the inclusive fitness of the dominant female compared to hypothetical 
# prolonged laying by a single female.

# Read in data and filtering/data modifications
experiment_data <- read_excel("Compiled_synchrony_experiment_data.xlsx")  %>%
  mutate(Nest_ID = paste(Year, Nest, sep = "_")) %>%
  group_by(Nest_ID, Final_manipulated_clutch) %>%
  mutate(Hatch_success = sum(Hatched, na.rm = TRUE) / Manipulated_clutch_size) %>%
  filter(!Treatment %in% c("Control", "Other")) %>%
  filter(! Nest_ID %in% c("2018_E", "2018_AD", "2018_AZ", "2018_BD", "2024_W")) %>%
  mutate(Foreign_percentage = Foreign_eggs/Manipulated_clutch_size) %>%
  mutate(Treatment = recode(Treatment,
                            "Synch" = "Synchronous",
                            "Asynch" = "Asynchronous")) %>%
  ungroup() 

# Hatch success plot
successful_nests <- experiment_data %>% filter(Hatch_success >0) %>%
  mutate(Hatch_begin = as.Date(Hatch_begin), Date_found = as.Date(Date_found), 
         Observed_incubation_period = as.numeric(Hatch_begin - Date_found)) %>%
  mutate(Hatch_begin = as.Date(Hatch_begin), Date_transfer = as.Date(Date_transfer),
    Transfered_time = as.numeric(Hatch_begin - Date_transfer)) %>%
  mutate(Estimated_hatch_spread = as.numeric(Estimated_hatch_spread)) %>%
  mutate(True_hatch_spread = as.numeric(True_hatch_spread))

ggplot(successful_nests, aes(x = Treatment, y = Hatch_success, colour = Treatment)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +  
  geom_jitter(position = position_jitter(width = 0.05, height = 0), size = 3, alpha = 0.7) + 
  labs(y = "Hatch success", x = "Treatment") +
  guides(color = FALSE) +
  theme_classic() +
  scale_color_manual(values = c("Asynchronous" = "darkblue", "Synchronous" = "firebrick2"))

# Number of swapped (foreign) eggs in nest plot
ggplot(successful_nests, aes(x = Treatment, y = Foreign_percentage, colour = Treatment)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +  
  geom_jitter(position = position_jitter(width = 0.05, height = 0), size = 3, alpha = 0.7) +  
  labs(y = "% of foreign eggs in nest", x = "Treatment") +
  guides(color = FALSE) +
  theme_classic() +
  scale_color_manual(values = c("Asynchronous" = "darkblue", "Synchronous" = "firebrick2"))

# Clutch size plot
ggplot(successful_nests, aes(x = Treatment, y = Manipulated_clutch_size, color = Treatment)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +  
  geom_jitter(position = position_jitter(width = 0.05, height = 0), size = 3, alpha = 0.7) +
  labs(y = "Clutch size", x = "Treatment") +
  guides(color = FALSE) +
  theme_classic() +
  scale_color_manual(values = c("Asynchronous" = "darkblue", "Synchronous" = "firebrick2"))

# Estimated hatch spread plot
ggplot(successful_nests, aes(x = Treatment, y = as.numeric(Estimated_hatch_spread), colour = Treatment)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.1, height = 0), size = 3, alpha = 0.7) +
  labs(y = "Estimated hatch spread (days)", x = "Treatment") +
  guides(color = FALSE) +
  theme_classic() +
  scale_color_manual(values = c("Asynchronous" = "darkblue", "Synchronous" = "firebrick2"))

# Survival to 60 days
data_long_format <- successful_nests %>%
  mutate(NonSurvivors = Hatched - Survival_60) %>% # Calculate non-survivors
  pivot_longer(cols = c(Survival_60, NonSurvivors), 
               names_to = "Status", 
               values_to = "Count") %>%
  mutate(Survival_60 = ifelse(Status == "Survival_60", 1, 0)) %>%
  uncount(Count) %>%
  mutate(Treatment_survival = paste(Treatment, Survival_60, sep = "_"))
View(data_long_format)

data_long_format %>%
  group_by(Treatment, Survival_60) %>%
  summarise(Count = n(), .groups = "drop")

ggplot(data_long_format, aes(x = Treatment, fill = factor(Treatment_survival))) +
  geom_bar(position = "dodge") +
  labs(x = "Treatment", y = "Count", fill = "Survived") +
  scale_fill_manual(values = c("Asynchronous_1" = "darkblue", "Asynchronous_0" = "blue", 
                               "Synchronous_1" = "firebrick4", "Synchronous_0" = "firebrick1"), 
                    guide = "none") +
  theme_classic()



# Summary numbers of sample size, min, max, mean, SD, and SE for variables. 
df_long <- successful_nests %>%
  pivot_longer(
    cols = where(is.numeric) & !all_of("Year"),
    names_to = "Variable",
    values_to = "Value")

df_summary <- df_long %>%
  group_by(Treatment, Variable) %>%
  summarise(
    n = n(),
    min = min(Value, na.rm = TRUE),
    max = max(Value, na.rm = TRUE),
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    se = sd / sqrt(n),
    .groups = "drop")
print(df_summary, n = 30)




#### T-TEST STUFF ####
t.test(Hatch_success ~ Treatment, data = successful_nests)
t.test(Foreign_percentage ~ Treatment, data = successful_nests)
t.test(Manipulated_clutch_size ~ Treatment, data = successful_nests)
t.test(Observed_incubation_period ~ Treatment, data = successful_nests)
t.test(Transfered_time ~ Treatment, data = successful_nests)
t.test(True_hatch_spread ~ Treatment, data = successful_nests)


  
#### POWER ANALYSIS STUFF ####

# Cohen's D power analysis stuff
cohens_d(Hatch_success ~ Treatment, data = successful_nests)
cohens_d(Survival_60 ~ Treatment, data = successful_nests)

SynchHatch <- successful_nests$Hatch_success[successful_nests$Treatment == "Synchronous"]
AsynchHatch <- successful_nests$Hatch_success[successful_nests$Treatment == "Asynchronous"]

hedges_g <- cohen.d(AsynchHatch, SynchHatch, hedges.correction = TRUE)
print(hedges_g)

# Setting the parameters
sample_size <- # Sample size
alpha <- 0.05        # Significance level
power <- 0.8         # Power
effect_size <- -0.3943032 # Hedge's G effect size  

# Solve for effect size needed
effect_size <- pwr.p.test(n = sample_size, sig.level = alpha, 
                          power = power, alternative = "two.sided")$h
print(effect_size)

# Solve for sample size needed
g_sample_size <- pwr.p.test(h = effect_size, sig.level = alpha, 
                          power = power, alternative = "less")$n
print(g_sample_size)



#### GLMM STUFF ####

model_1 <- glmmTMB(True_hatch_spread ~ Manipulated_clutch_size + (1|Year),
                       family = poisson, data = successful_nests)
check_model(model_1)

model_2 <- glmmTMB(Survival_60 ~ Treatment + Manipulated_clutch_size +
                      (1|Year), family = binomial, data = data_long_format)
check_model(model_2)

model_3 <- glmmTMB(Survival_60 ~ Treatment + Hatched +
                     (1|Year), family = binomial, data = data_long_format)
check_model(model_3)

model_4 <- glmmTMB(Survival_60 ~ Treatment + Hatched +
                  (1|Year), family = binomial, data = data_long_format)
check_model(model_4)

model_5 <- glmmTMB(Hatch_success ~ Treatment + Manipulated_clutch_size +
                     (1|Year), family = gaussian, data = data_long_format)
check_model(model_5)









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



#### Weighting nests by hatching synchrony using Shannon entropy and by hatching success. 
# Score calculated by a weighted hatching synchrony multiplied by hatching
# success to give new weight for nests combining both factors. Nests with higher
# synchrony and higher success weighted the highest. 

# 1. Entropy-weighting function
# alpha: emphasis on synchrony (entropy). Can choose to weight synchrony greater 
# or less than hatching success. 
# beta: emphasis on hatching success. Can choose to weight hatching success greater 
# or less than synchrony. 
# epsilon: small constant to prevent divide-by-zero
calculate_synchrony_weights <- function(data, alpha = 1, beta = 1, epsilon = 0.1) {
  
  results <- lapply(names(data), function(nest_name) {
    nest <- nest_data[[nest_name]]
    hatch_days <- nest$hatch_days
    eggs_laid <- nest$eggs_laid
    eggs_hatched <- length(hatch_days)
    
    # Hatching success
    success_rate <- eggs_hatched / eggs_laid
    
    # Entropy of hatch days
    p <- table(hatch_days) / eggs_hatched
    entropy <- -sum(p * log2(p))
    
    # Inversing synchrony weight (lower entropy = higher weight)
    synchrony_weight <- 1 / (entropy + epsilon)
    
    # Final combined weight with adjustable emphasis
    combined_weight <- (synchrony_weight ^ alpha) * (success_rate ^ beta)
    
    return(data.frame(
      Nest = nest_name,
      Entropy = round(entropy, 3),
      SuccessRate = round(success_rate, 3),
      RawWeight = combined_weight
    ))
  }) %>% bind_rows()
  
  # Re-scaling weights to between 0 and 1. 
  results$ScaledWeight <- round(results$RawWeight / max(results$RawWeight), 3)
  return(results)
}


# Example data
nests <- list(
  A = list(hatch_days = c(1,1,1,1,1), clutch_size = 5),
  B = list(hatch_days = c(1,1,1,1,2), clutch_size = 8),
  C = list(hatch_days = c(1,1,2,3), clutch_size = 6),
  D = list(hatch_days = c(1,2,4), clutch_size = 5),
  E = list(hatch_days = c(1,2,3,4,5), clutch_size = 5),
  F = list(hatch_days = c(1,1,1,2), clutch_size = 5),
  G = list(hatch_days = c(1,1,1,1,2,2,2), clutch_size = 8),
  H = list(hatch_days = c(1,1,2,2), clutch_size = 6),
  I = list(hatch_days = c(1,1,3), clutch_size = 4),
  J = list(hatch_days = c(1,1,2,4,5), clutch_size = 7)  
)


# Calculating Shannon entropy function
calculate_entropy <- function(hatch_days) {
  if (length(hatch_days) == 0) return(NA_real_) # No chicks hatched
  p <- table(hatch_days) / length(hatch_days)
  if (length(p) == 1) return(0)  # All on same day = perfect synchrony
  return(-sum(p * log2(p)))
}

# Main function 
calculate_synchrony_weights <- function(nest_data, alpha = 1, beta = 1, epsilon = 0.1) {
  results <- lapply(names(nest_data), function(nest_name) {
    nest <- nest_data[[nest_name]]
    hatch_days <- nest$hatch_days
    clutch_size <- nest$clutch_size
    eggs_hatched <- length(hatch_days)
    
    # Bad nests get skipped
    if (clutch_size == 0 || eggs_hatched == 0) return(NULL)
    
    hatch_success <- eggs_hatched/clutch_size
    entropy <- calculate_entropy(hatch_days)
    if (is.na(entropy)) return(NULL)  # just in case
    
    synchrony_weight <- 1 / (entropy + epsilon)
    combined_weight <- (synchrony_weight ^ alpha) * (hatch_success ^ beta)
    
    data.frame(
      Nest = nest_name,
      Entropy = round(entropy, 3),
      Hatch_success = round(hatch_success, 3),
      RawWeight = combined_weight
    )
  }) %>% bind_rows()
  
  # Re-scale weights between 0 and 1
  results$ScaledWeight <- round(results$RawWeight / max(results$RawWeight), 3)
  return(results)
}

# Printed results by nest in alphabetical order
weights <- calculate_synchrony_weights(nests)
print(weights)

# Printed results by nest in order of weights highest to lowest
weights <- calculate_synchrony_weights(nests, alpha = 1, beta = 1, epsilon = 0.1) %>%
  arrange(desc(ScaledWeight))
print(weights)
####

################################
############CHAPTER 2 ##########
#### NATURAL SYNCHRONY DATA ####
################################

masterlist_data <- read_excel("Nests_masterlist.xlsx")
masterlist_data <- masterlist_data %>% filter(Exclusion == "GOOD" | Exclusion == "ABCL") %>%
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>%
  mutate(Females = ifelse(Clutch_size <= 5, "Single Female", "Joint Females")) %>%
  filter(Hatch_begin != "NA" & Hatch_end != "NA") %>%
  mutate(
    Hatch_begin = as.numeric(Hatch_begin),
    Hatch_end = as.numeric(Hatch_end),
    Hatch_spread = Hatch_end - Hatch_begin + 1) %>%
  filter(survive_1 != 0 & survive_1 != "NA") %>%
  filter(Survived_2 != "NA" & Survived_2 != "2?")
View(masterlist_data)

ggplot(masterlist_data, aes(x = as.numeric(Hatched_eggs), y = survive_1)) +
  geom_point()



masterlist_data <- read_excel("Nests_masterlist.xlsx")
masterlist_data <- masterlist_data %>% filter(Exclusion == "GOOD" | Exclusion == "ABCL") %>%
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>%
  mutate(Females = ifelse(Clutch_size <= 5, "Single Female", "Joint Females")) %>%
  mutate(
    Hatch_begin = as.numeric(Hatch_begin),
    Hatch_end = as.numeric(Hatch_end),
    Hatch_spread = Hatch_end - Hatch_begin + 1) %>%
  filter(survive_1 != 0 & survive_1 != "NA") %>%
  filter(Survived_2 != "NA" & Survived_2 != "2?")
View(masterlist_data)

ggplot(masterlist_data, aes(x = as.numeric(Hatched_eggs), y = survive_1)) +
  geom_point()



# Survival data conservative estimate with number of banded juveniles as true number



# Survival data with liberal estimate number of chicks that hatched
mutate(Survival_60 = Hatched_eggs)
mutate(Survival_60 = Hatched_eggs)