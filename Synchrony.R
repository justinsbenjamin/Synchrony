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


# Get residuals: asynchrony relative to expected span for brood size
Hatch_success <- masterlist_data %>%
  filter(Hatch_success != "NA")
model <- lm(Hatch_success ~ Clutch_size, data = Hatch_success)
Hatch_success$ha_resid <- resid(model)

hist(Hatch_success$ha_resid, breaks = 20)




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


model <- lm(Hatch_spread ~ Hatched_eggs, data = masterlist_data)

# Get residuals: asynchrony relative to expected span for brood size
masterlist_data$hatching_asynchrony <- resid(model)
masterlist_data <- masterlist_data %>%
  filter(Clutch_size => 5)

hist(masterlist_data$hatching_asynchrony)

masterlist_data <- masterlist_data %>%
  filter(Hatched_eggs > 1)
ggplot(masterlist_data, aes(x = as.numeric(Hatched_eggs), y = as.numeric(Hatch_spread), colour = Females)) +
  geom_point() +
  geom_jitter(position = position_jitter(width = 0.5, height = 0), size = 2, alpha = 0.7) +
  theme_classic()


ggplot(masterlist_data, aes(x = as.numeric(hatching_asynchrony), y = as.numeric(Survived_2), colour = Females)) +
  geom_point() +
  geom_jitter(position = position_jitter(width = 0.5, height = 0), size = 2, alpha = 0.7) +
  scale_x_continuous(limits = c(-3, 3)) +
  theme_classic()


ggplot(masterlist_data, aes(x = as.numeric(Clutch_size), 
                            y = as.numeric(Hatch_success), colour = Females)) +
  geom_point() +
  theme_classic()

masterlist_data <- masterlist_data %>%
  filter(Hatched_eggs >1) 
  
ggplot(masterlist_data, aes(x = as.numeric(Hatch_spread), 
                            y = as.numeric(Hatch_success), colour = Females)) +
  geom_point() +
  geom_jitter(position = position_jitter(width = 0.5, height = 0), size = 2, alpha = 0.7) +
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


ggplot(masterlist_data, aes(x = as.numeric(Hatched_eggs), y = survive_1)) +
  geom_point()

# Survival data conservative estimate with number of banded juveniles as true number

# Survival data with liberal estimate number of chicks that hatched
mutate(Survival_60 = Hatched_eggs)
mutate(Survival_60 = Hatched_eggs)







########## CHAPTER 3 #########
#### SYNCHRONY EXPERIMENT ####
##############################

# I hypothesize that jointly laying synchronously with a close relative increases
# the inclusive fitness of the dominant female compared to hypothetical 
# prolonged laying by a single female.

# Read in data and filtering/data modifications
# Nests 2018_E, 2018_AD, 2018_AZ, 2018_BD were not manipulated and thus excluded. 
# Nest 2024_W predated during hatching, also excluded. 
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

# Further cleaning up the data set with successful nests
successful_nests <- experiment_data %>% filter(Hatch_success >0) %>%
  mutate(Hatch_begin = as.Date(Hatch_begin), Date_found = as.Date(Date_found), 
         Observed_incubation_period = as.numeric(Hatch_begin - Date_found)) %>%
  mutate(Hatch_begin = as.Date(Hatch_begin), Date_transfer = as.Date(Date_transfer),
    Transfered_time = as.numeric(Hatch_begin - Date_transfer)) %>%
  mutate(Estimated_hatch_spread = as.numeric(Estimated_hatch_spread)) %>%
  mutate(True_hatch_spread = as.numeric(True_hatch_spread))

# Nest data in pivoted longer format for faceted figure and summary stats.
df_long <- successful_nests %>%
  pivot_longer(
    cols = where(is.numeric) & !all_of("Year"),
    names_to = "Variable",
    values_to = "Value")

# Summary numbers of sample size, min, max, mean, SD, and SE for variables.
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

# Faceted figure with all numeric things. Can filter it down in the df_long_filtered
# data set to only do a handful of panels. 
df_long_filtered <- subset(df_long, Variable != "Number_transfers")

ggplot(df_long_filtered, aes(x = Treatment, y = as.numeric(Value), colour = Treatment)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.1, height = 0), size = 3, alpha = 0.7) +
  facet_wrap(~Variable, scales="free_y") +
  guides(color = FALSE) +
  theme_classic() +
  scale_color_manual(values = c("Asynchronous" = "darkblue", "Synchronous" = "firebrick2")) +
  labs(x = NULL, y = NULL) 

# Example template for doing a single panel instead of facet wrapping them together. 
ggplot(successful_nests, aes(x = Treatment, y = as.numeric(Estimated_hatch_spread), colour = Treatment)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.1, height = 0), size = 3, alpha = 0.7) +
  labs(y = "Estimated hatch spread (days)", x = "Treatment") +
  guides(color = FALSE) +
  theme_classic() +
  scale_color_manual(values = c("Asynchronous" = "darkblue", "Synchronous" = "firebrick2"))

# Survival to 60 days in pivoted longer format for binomial stuff and figure.
data_long_format <- successful_nests %>%
  mutate(NonSurvivors = Hatched - Survival_60) %>% # Calculate non-survivors
  pivot_longer(cols = c(Survival_60, NonSurvivors), 
               names_to = "Status", 
               values_to = "Count") %>%
  mutate(Survival_60 = ifelse(Status == "Survival_60", 1, 0)) %>%
  uncount(Count) %>%
  mutate(Treatment_survival = paste(Treatment, Survival_60, sep = "_")) %>%

# Barplot of survival data.
ggplot(data_long_format, aes(x = Treatment, fill = factor(Treatment_survival))) +
  geom_bar(position = "dodge") +
  labs(x = "Treatment", y = "Count", fill = "Survived") +
  scale_fill_manual(values = c("Asynchronous_1" = "darkblue", "Asynchronous_0" = "blue", 
                               "Synchronous_1" = "firebrick4", "Synchronous_0" = "firebrick1"), 
                    guide = "none") +
  theme_classic()

data_long_format %>%
  group_by(Treatment, Survival_60) %>%
  summarise(Count = n(), .groups = "drop")

# T-TEST STUFF 
t.test(Hatch_success ~ Treatment, data = successful_nests)
t.test(Foreign_percentage ~ Treatment, data = successful_nests)
t.test(Manipulated_clutch_size ~ Treatment, data = successful_nests)
t.test(Observed_incubation_period ~ Treatment, data = successful_nests)
t.test(Transfered_time ~ Treatment, data = successful_nests)
t.test(True_hatch_spread ~ Treatment, data = successful_nests)

# POWER ANALYSIS STUFF

# Cohen's D effect size calculation
cohens_d(Hatch_success ~ Treatment, data = successful_nests)
cohens_d(Survival_60 ~ Treatment, data = successful_nests)

# Hedge's G effect size calculation
# This is better for my small sample size
SynchHatch <- successful_nests$Hatch_success[successful_nests$Treatment == "Synchronous"]
AsynchHatch <- successful_nests$Hatch_success[successful_nests$Treatment == "Asynchronous"]

hedges_g <- cohen.d(AsynchHatch, SynchHatch, hedges.correction = TRUE)
print(hedges_g)

# Setting the parameters for power analysis
sample_size <-            # Sample size
alpha <- 0.05             # Significance level
power <- 0.8              # Power
effect_size <- hedges_g # Hedge's G effect size  

# Solve for effect size needed
## Figure out 2 side or one sided analysis ##
effect_size <- pwr.p.test(n = sample_size, sig.level = alpha, 
                          power = power, alternative = "two.sided")$h
print(effect_size)

# Solve for sample size needed
g_sample_size <- pwr.p.test(h = effect_size, sig.level = alpha, 
                          power = power, alternative = "less")$n
print(g_sample_size)


# GLMM STUFF

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




# Figures made for GRD, not sure I'll need these again. 
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







