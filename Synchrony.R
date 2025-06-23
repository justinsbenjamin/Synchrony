########### JUSTIN BENJAMIN ##########
########### MSC THESIS CODE ##########
######################################

# Load libraries
library(dplyr)
library(ggplot2)
library(ggpattern)
library(effects)
library(MASS)
library(performance)
library(tidyr)
library(readxl)
library(lubridate)
library(pwr)
library(rstatix)
library(glmmTMB)
library(effsize)
library(DHARMa)

########### CHAPTER 2 ##########
###### OBSERVATIONAL DATA ######
################################

#### CH.2 INTRO AND HYPOTHESES ####
# I am investigating nesting factors in Pūkeko (Porphyrio melanotus melanotus) 
# to explore the effects of hatching spread (HS) on hatching success and 
# survival. The data for chapter 2 has been compiled from field books and 
# spreadsheets to create a master data set with both summary data for all 
# nests along with specific hatching data for individual nests. Data is from 
# the years 2008, 2010, 2013, 2014-2018, 2022-2024. 

# HYPOTHESIS 1: Dominant females lay synchronously with another female to 
# increase their inclusive fitness. 
# PREDICTION 1a: Shorter HS have better hatching success
# PREDICTION 1b: Shorter HS have better survival

# HYPOTHESIS 2: Females lay together because the optimal clutch size is 
# greater than the size of single female clutches regardless of hatch spread.
# PREDICTION 2a: Larger clutches hatch more eggs 
# PREDICTION 2b: Larger broods have more survivors
# PREDICTION 2c: Larger clutches have greater hatching spread 

# HYPOTHESIS 3: Subordinate females benefit more from synchronous nests
# PREDICTION 3a: Hatch order less important with short HS

#### DATA FILTERING AND ORGANIZATION ####
# Sample sizes of excluded and included groups
nest_counts <- read_excel("Nests_masterlist.xlsx") %>%
count(Exclusion, name = "nest_counts") %>%
arrange(desc("nest_counts"))
print(nest_counts)

masterlist_data <- read_excel("Nests_masterlist.xlsx") %>%
  mutate(Nest_ID                 = paste(Year, Nest_ID, sep = "_"), # Unique Nest ID 
         Hatch_success           = as.numeric(Hatched_eggs)/as.numeric(Clutch_size), 
         Survived_2_cons         = ifelse(is.na(Survived_2), 0, Survived_2), 
         Hatch_begin             = as.numeric(Hatch_begin),
         Hatch_end               = as.numeric(Hatch_end),
         Hatch_spread            = Hatch_end - Hatch_begin + 1, 
         Hatch_begin             = as.numeric(Hatch_begin), 
         Date_found              = as.numeric(Date_found), 
         Observed_nesting_period = as.numeric(Hatch_begin - Date_found))
  mutate(across(c(Clutch_size, Hatched_eggs, Survived_2, Survived_2_cons,
                  Hatch_begin, Hatch_end, Hatch_spread, Year), 
                  ~ as.numeric(as.character(.))))
View(masterlist_data)

#### HYPOTHESIS 1 ####
# Dominant females lay synchronously with another female to increase 
# their inclusive fitness. 

#### PREDICTION 1a: Shorter HS have greater hatching success ####
model_2a <- glmmTMB(Hatched_eggs ~ Hatch_spread + (1|Year) 
                   ,family = poisson, data = masterlist_data)
check_model(model_2a)
res <- simulateResiduals(fittedModel = model_2a, plot = TRUE)
testDispersion(res)

model_2b <- glmmTMB(Hatched_eggs ~ Hatch_spread + I(Hatch_spread^2) + (1|Year) 
                    ,family = poisson, data = masterlist_data)
check_model(model_2b)
res <- simulateResiduals(fittedModel = model_2b, plot = TRUE)
testDispersion(res)

model_2c <- glmmTMB(Hatched_eggs ~ I(Hatch_spread^2) + (1|Year) 
                    ,family = nbinom2, data = masterlist_data)
check_model(model_2c)
res <- simulateResiduals(fittedModel = model_2c, plot = TRUE)
testDispersion(res)

null_model <- glmmTMB(Hatched_eggs ~ 1 + (1|Year), family = nbinom2, 
                      data = masterlist_data)
check_model(null_model)
res_null <- simulateResiduals(fittedModel = null_model, plot = TRUE)
plot(res_null)

summary(masterlist_data$Hatched_eggs)
hist(masterlist_data$Hatched_eggs)


#### PREDICTION 1b: Shorter HS have better survival ####

# Model analyzing survival by HS with conservative survival estimate
model_2b <- glmmTMB(Survived_2_cons ~ Hatch_spread + (1|Year) 
                   ,family = poisson(link = log), data = masterlist_data)
check_model(model_2b)

res <- simulateResiduals(fittedModel = model_2b, plot = TRUE)
testDispersion(res)
testZeroInflation(res) 

summary(model_2b)

# liberal survival estimate
lib_data <- masterlist_data %>%
  filter(!is.na(Survived_2))
model_2c <- glmmTMB(Survived_2 ~ Hatch_spread + (1|Year) 
                   ,family = poisson(link = log), data = lib_data)
check_model(model_2c)

res <- simulateResiduals(fittedModel = model_2c, plot = TRUE)
testDispersion(res)
testZeroInflation(res) 
summary(model_2c)

# Faceted figure of hatching rate, hatched eggs, and conservative and liberal
# estimates of survival by hatch spread

# Pivoting data longer to make faceted figure
longer_data <- masterlist_data %>%
  pivot_longer(cols = c('Hatched_eggs', 'Hatch_success',
                        'Survived_2'), 
               names_to = 'Data', 
               values_to = 'Value') 
View(longer_data)

df_long <- df %>%
  mutate(liberal = ifelse(is.na(survival), NA, survival),
         conservative = ifelse(is.na(survival), 0, survival)) %>%
  pivot_longer(cols = c(liberal, conservative), names_to = "estimate", 
               values_to = "survival")


df_long <- df %>%
  pivot_longer(cols = c('Hatched_eggs', 'Hatch_success', 'Survived_2'),
               names_to = "Data",
               values_to = "Value") %>%
  mutate(converted_from_na = ifelse(Data == "Survived_2", 
                                    masterlist_data$converted_from_na, FALSE))


longer_data <- masterlist_data %>%
  pivot_longer(cols = c(Hatched_eggs, Hatch_success, Survived_2_cons),
               names_to = "Data",
               values_to = "Value") %>%
  mutate(converted_from_na = ifelse(Data == "Survived_2_cons", 
                                    converted_from_na, FALSE))



ggplot(longer_data, aes(x = Hatch_spread, y = Value, color = converted_from_na)) +
  geom_jitter(width = 0.2, height = 0.05, alpha = 0.8) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                     name = "Converted from NA") +
  facet_wrap(~Data, scales = "free_y") +
  theme_classic() +
  labs(x = NULL, y = NULL)


#### HYPOTHESIS 3 ####
#### PREDICTION 3a: Hatch order less important in synchronous nests ####

chick_data <- read_excel("Final_Compiled_Chick_Data.xlsx") %>%
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>%
  filter(!is.na(Survived) & Survived %in% c(0, 1)) %>%
  mutate(Survived = if (!is.numeric(Survived)) as.numeric(Survived) else Survived)
View(chick_data)

chick_data$Clutch_size <- masterlist_data$Clutch_size[match(chick_data$Nest_ID,
                                          masterlist_data$Nest_ID)]
chick_data$Hatch_spread <- masterlist_data$Hatch_spread[match(chick_data$Nest_ID, 
                                           masterlist_data$Nest_ID)]

# Model with hatch order as sole predictor variable. 
# Same as Cody's paper but with Nest ID nested within year. 
model_d <- glmmTMB(Survived ~ Hatch_order +
                     (1|Year), family = binomial, data = chick_data)
check_model(model_d)
summary(model_d)

# Model with interaction between hatch order clutch size
model_d <- glmmTMB(Survived_2 ~ Hatch_order*Clutch_size +
                     (Nest_ID/Year), family = binomial, data = chick_data)
check_model(model_d)
summary(model_d)

# Model with interaction between hatch order and hatching spread
model_d <- glmmTMB(Survived ~ Hatch_order*Hatch_spread +
                     (1|Nest_ID), family = binomial, data = chick_data)
check_model(model_d)
res <- simulateResiduals(fittedModel = model_d, plot = TRUE)
testDispersion(res)
testZeroInflation(res) 


summary(model_d)

ggplot(chick_data, aes(x = Hatch_order, y = Survived)) +
  geom_jitter(width = 0.5, height = 0, size = 2, alpha = 0.7) +
  labs(y = "Survived", x = "Hatch_spread") +
  theme_classic()


# Model of laying order on hatching order 
model_d <- glmmTMB(Hatch_order ~ Lay_order +
                     (Nest_ID/Year), family = gaussian, data = chick_data)
check_model(model_d)
summary(model_d)

# Model of laying spread on hatching spread
model_d <- glmmTMB(Hatch_spread ~ Lay_spread +
                     (Nest_ID/Year), family = poisson, data = chick_data)
check_model(model_d)
summary(model_d)

# Binomial figure of survival by hatch order
ggplot(chick_data, aes(x = Hatch_order, y = Survived)) +
  geom_jitter(width = 0.5, height = 0, size = 2, alpha = 0.7) +
  labs(y = "Survived", x = "Hatch order") +
  theme_classic()

# Filtering out data that has unusually large birds/potential recaptures 
# and pivoting data longer. 
chick_data_longer <- chick_data %>%
  filter(Mass < 33) %>%
  filter(`Shield to Tip` < 24) %>%
  pivot_longer(cols = c('Mass', 'Tarsus',
                        'Shield to Tip'), 
               names_to = 'Morphometrics', 
               values_to = 'Value') %>%
  filter(!is.na(Value))

View(chick_data_longer)

# Figure of size by hatch order
# Can also swap out Hatch_order for Hatch_Day
ggplot(chick_data_longer, aes(x = Hatch_order, y = Value, colour = Year)) +
  geom_jitter(width = 0.25, height = 0, size = 2, alpha = 0.7) +
  facet_wrap(~Morphometrics, scales="free_y", labeller = labeller(
    Morphometrics = c(Mass = "Mass (g)",
                      Tarsus = "Left outer tarsus (mm)", 
                      `Shield to Tip` = "Shield to tip (mm)"))) +
  scale_y_continuous() +
  theme_classic()

# Figure of survival by size at hatching 
ggplot(chick_data_longer, aes(x = Value, y = as.character(Survived), 
                              colour = Year)) +
  geom_point(size = 2, alpha = 0.7) +
  facet_wrap(~Morphometrics, scales="free_x", labeller = labeller(
    Morphometrics = c(Mass = "Mass (g)",
                      Tarsus = "Left outer tarsus (mm)", 
                      `Shield to Tip` = "Shield to tip (mm)"))) +
  labs(x = "Measurement", y = "Survival") +
  scale_x_continuous() +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))


chick_data <- chick_data %>%
  filter(Mass < 33) %>%
  filter(`Shield to Tip` < 24)

# Model of laying spread on hatching spread
model_d <- glmmTMB(Mass ~ Hatch_order +
                     (1|Nest_ID), family = gaussian, data = chick_data)
check_model(model_d)
res <- simulateResiduals(fittedModel = model_d, plot = TRUE)
testDispersion(res)
summary(model_d)


model_d <- glmmTMB(`Shield to Tip` ~ Hatch_order +
                     (1|Nest_ID), family = gaussian, data = chick_data)
check_model(model_d)
res <- simulateResiduals(fittedModel = model_d, plot = TRUE)
testDispersion(res)
summary(model_d)

model_quad <- glmmTMB(Survived ~ Mass + I(Mass^2) + (1 | Nest_ID), 
                      family = binomial, data = chick_data)


chick_data$logMass <- log(chick_data$Mass)
model_log <- glmmTMB(Survived ~ logMass + (1|Nest_ID), family = binomial, data = 
                       chick_data)

model_d <- glmmTMB(Survived ~ Mass + `Shield to Tip` +
                     (1|Nest_ID), family = binomial, data = chick_data)
check_model(model_d)
res <- simulateResiduals(fittedModel = model_d, plot = TRUE)
testDispersion(res)
testZeroInflation(res) 


#### HYPOTHESIS 2 ####
# Females lay together because the optimal clutch size is greater 
# than the size of single female clutches.

# PREDICTION 2a: Larger clutches hatch more eggs

# Model analyzing number of hatched eggs by clutch size
masterlist_data <- masterlist_data%>%
  filter(Hatched_eggs > 0)

hatched_eggs <- glmmTMB(Hatched_eggs ~ I(Clutch_size^2) +
                          (1|Year), family = poisson(link = "log"), 
                        data = masterlist_data)
check_model(hatched_eggs)
res <- simulateResiduals(fittedModel = hatched_eggs, plot = TRUE)
testDispersion(res)


# Model analyzing hatching rate by clutch size
hatch_success <- glmmTMB(Hatch_success ~ Clutch_size + 
                          (1|Year), family = , data = masterlist_data)
check_model(hatch_success)
res <- simulateResiduals(fittedModel = hatched_eggs, plot = TRUE)
testDispersion(res)


hatch_spread1 <- glmmTMB(Hatch_spread ~ Clutch_size + Hatched_eggs +
                        (1|Year), family = nbinom2, data = masterlist_data)
check_model(hatch_spread1)
res <- simulateResiduals(fittedModel = hatch_spread1, plot = TRUE)
testDispersion(res)


hatch_spread <- glmmTMB(Hatch_spread ~ Hatched_eggs + 
                       (1|Year), family = poisson, data = masterlist_data)
check_model(hatch_spread)
res <- simulateResiduals(fittedModel = hatch_spread, plot = TRUE)
testDispersion(res)
summary(hatch_spread)

# PREDICTION 2b: Larger broods have more survivors




# Conservative model analyzing number of survivors by clutch size and 
# te number of hatched eggs
survived_cons_clutch <- glmmTMB(Survived_2_cons ~ I(Clutch_size^2) + Hatched_eggs +
                        (1|Year), family = poisson(link = "log"), 
                        data = masterlist_data)
check_model(survived_cons)
res <- simulateResiduals(fittedModel = survived_cons_clutch, plot = TRUE)
testDispersion(res)
testZeroInflation(res) 

survived_cons_clutch <- glmmTMB(Survived_2_cons ~ Clutch_size + Hatched_eggs +
                               (1|Year), family = poisson(link = "log"), 
                               data = masterlist_data)
check_model(survived_cons)
res <- simulateResiduals(fittedModel = survived_cons_clutch, plot = TRUE)
testDispersion(res)
testZeroInflation(res) 


# Liberal model analyzing number of survivors by clutch size and # of hatched eggs
survived_clutch <- glmmTMB(Survived_2 ~ Clutch_size + Hatched_eggs
                          (1|Year), family = poisson, data = lib_data)
check_model(hatched_eggs)
res <- simulateResiduals(fittedModel = survived_cons_clutch, plot = TRUE)
testDispersion(res)
testZeroInflation(res) 

# Pivoting data longer to make faceted figure
longer_data <- masterlist_data %>%
pivot_longer(cols = c('Hatched_eggs', 'Hatch_success',
                      'Survived_2', 'Survived_2_cons'), 
             names_to = 'Data', 
             values_to = 'Value')

# Faceted figure of hatching success, hatched eggs, and conservative and liberal
# estimates of survival by clutch size.
ggplot(longer_data, aes(x = Clutch_size, y = Value)) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 2, 
              alpha = 0.7) +
  facet_wrap(~Data, scales="free_y", labeller = labeller(
    Data = c(Hatched_eggs = "Number of hatched eggs",
             Hatch_success = "Hatching rate (%)", 
             Survived_2 = "Number of survivors",
             Survived_2_cons = "Conservative number of survivors"))) +
  labs(y= "Counts", x = "Clutch size") +
  geom_smooth(method = "loess") +
  theme_classic()

longer_data_2 <- masterlist_data %>%
  pivot_longer(cols = c('Hatch_success',
                        'Survived_2', 'Survived_2_cons'), 
               names_to = 'Data', 
               values_to = 'Value')

# Figure of number survivors by number of hatched eggs with liberal estimate
ggplot(longer_data_2, aes(x = Hatched_eggs, y = Value)) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 2, 
              alpha = 0.7) +
  facet_wrap(~Data, scales="free_y", labeller = labeller(
    Data = c(Hatch_success = "Hatching rate (%)", 
             Survived_2 = "Number of survivors",
             Survived_2_cons = "Conservative number of survivors"))) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 8, by = 1))







### Supplementary Chapter 2 analyses ###
### Most descriptive model to break it ###
# Model with all fixed effects as random effects to break and have singular fit. 
model_a <- glmmTMB(as.numeric(Hatched_eggs) ~ as.numeric(Clutch_size) + 
                     as.numeric(Group_size) + as.numeric(Hatch_spread) +
                     (1|Year) + (1|as.numeric(Clutch_size)) + 
                     (1|as.numeric(Group_size)) +
                     (1|as.numeric(Hatch_spread)), family = poisson(link = log), 
                     data = masterlist_data)
check_model(model_a)
summary(model_a)


########## CHAPTER 3 #########
#### SYNCHRONY EXPERIMENT ####
##############################

#### CH.3 INTRO AND HYPOTHESIS ####
# Here we experimentally lengthened and shortened the HS of nests 
# in Pūkeko by using known laying dates, candling images, and flotation scores.

# HYPOTHESIS: Jointly laying synchronously with a close relative increases
# the inclusive fitness of the dominant female compared to hypothetical 
# prolonged laying of a similar sized clutch by a single female.

# PREDICTION A: Synchronous nests have greater hatching success than 
# Asynchronous nests as there are fewer eggs deserted in the nest after hatching. 
# PREDICTION B: Synchronous nests have greater survival than Asynchronous nests 
# as there are fewer late hatching chicks with poorer survival chances. 

#### DATA FILTERING AND ORGANIZATION #### 
# Nests 2018_E, 2018_AD, 2018_AZ, 2018_BD were not manipulated and thus excluded. 
# Nest 2024_W predated during hatching, also excluded. 
experiment_data <- read_excel("Compiled_synchrony_experiment_data.xlsx") %>%
  mutate(Nest_ID = paste(Year, Nest, sep = "_")) %>%
  filter(!Treatment %in% c("Control","Other")) %>%
  filter(!Nest_ID %in% c("2018_E","2018_AD","2018_AZ","2018_BD","2024_W")) %>%  
  mutate(Treatment = recode(Treatment, "Synch" = "Synchronous",
                                       "Asynch" = "Asynchronous"), 
         Hatch_success = (Hatched/Manipulated_clutch_size), 
         Foreign_percentage = (Foreign_eggs/Manipulated_clutch_size))

# Further cleaning of the data set with successful nests
successful_nests <- experiment_data %>% 
  filter(Hatch_success >0) %>%
  mutate(Date_found = as.Date(Date_found),
         Hatch_begin = as.Date(Hatch_begin), 
         Observed_nesting_period = as.numeric(Hatch_begin - Date_found)) %>%
  mutate(Date_transfer = as.Date(Date_transfer),
         Transfered_time = as.numeric(Hatch_begin - Date_transfer)) %>%
  mutate(Estimated_hatch_spread = as.numeric(Estimated_hatch_spread),
         True_hatch_spread = as.numeric(True_hatch_spread))

# Function to pivot data longer for both the full experiment data set and the 
# successful nests data set, then give summary statistics. 

# Function to pivot data long
pivot_data_long <- function(df, dataset_name = NA) {
  df %>%
    pivot_longer(
      cols = where(is.numeric) & !all_of("Year"),
      names_to = "Variable",
      values_to = "Value"
    ) %>%
    mutate(Source = dataset_name)}

pivot_data_long(experiment_data, "experiment")
View(pivot_data_long)

# Function to summarize pivoted data
summarize_data <- function(df_long) {
  df_long %>%
    group_by(Treatment, Variable) %>%
    summarise(
      n = n(),
      min = min(Value, na.rm = TRUE),
      max = max(Value, na.rm = TRUE),
      mean = mean(Value, na.rm = TRUE),
      sd = sd(Value, na.rm = TRUE),
      se = sd / sqrt(n),
      .groups = "drop"
    )
}








summarize_experiment_data <- function(df, dataset_name = NA) {
  df_long <- df %>%
    pivot_longer(
      cols = where(is.numeric) & !all_of("Year"),
      names_to = "Variable",
      values_to = "Value"
    ) %>%
    mutate(Source = dataset_name)
  
  df_summary <- df_long %>%
    group_by(Treatment, Variable) %>%
    summarise(
      n = n(),
      min = min(Value, na.rm = TRUE),
      max = max(Value, na.rm = TRUE),
      mean = mean(Value, na.rm = TRUE),
      sd = sd(Value, na.rm = TRUE),
      se = sd / sqrt(n),
      .groups = "drop"
    )
  
  # Return both quietly
  list(long = df_long, summary = df_summary)
}

result <- summarize_experiment_data(experiment_data, "experiment")

# View or print what you want
print(result$long, n = 30)
print(result$summary, n = 30)

# Or interactively:
View(result$long)
View(result$summary)






summarize_experiment_data <- function(df, dataset_name = NA) {
  # Pivot data to long format
  df_long <- df %>%
    pivot_longer(
      cols = where(is.numeric) & !all_of("Year"),
      names_to = "Variable",
      values_to = "Value"
    ) %>%
    mutate(Source = dataset_name)
  
  # Print the pivoted data
  cat("\nPivoted (long) data for:", dataset_name, "\n")
  print(df_long, n = Inf)
  
  # Summarize the data
  df_summary <- df_long %>%
    group_by(Treatment, Variable) %>%
    summarise(
      n = n(),
      min = min(Value, na.rm = TRUE),
      max = max(Value, na.rm = TRUE),
      mean = mean(Value, na.rm = TRUE),
      sd = sd(Value, na.rm = TRUE),
      se = sd / sqrt(n),
      .groups = "drop"
    )
  
  # Print the summary
  cat("\nSummary data for:", dataset_name, "\n")
  print(df_summary, n = Inf)
  
  # Return both
  return(list(long = df_long, summary = df_summary))
}

result1 <- summarize_experiment_data(experiment_data, "experiment")
print(result1$long, n = 30)
print(result1$summary, n = 30)
result <- summarize_experiment_data(experiment_data, "experiment")
result <- summarize_experiment_data(experiment_data, "experiment")

print(result)

summarize_experiment_data



summarize_experiment_data <- function(df, dataset_name = NA) {
  df_long <- df %>%
    pivot_longer(
    cols = where(is.numeric) & !all_of("Year"),
    names_to = "Variable",
    values_to = "Value") %>%
    mutate(Source = dataset_name)

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
  
# Print summary
  cat("\nSummary for:", dataset_name, "\n")
  print(df_summary, n = 30)
  
  return(list(long = df_long, summary = df_summary))
}

summarize_experiment_data(experiment_data, "experiment")





result2 <- summarize_experiment_data(successful_nests, "control")

experiment_summary <- result1$summary

control_data_long <- result2$long
control_summary <- result2$summary










# Nest data in pivoted longer format for faceted figure and summary stats.
experiment_data_long <- experiment_data %>%
  pivot_longer(
    cols = where(is.numeric) & !all_of("Year"),
    names_to = "Variable",
    values_to = "Value")

# Create a reusable summary function
summarize_data <- function(df) {
    df %>%
    group_by(Treatment, Variable) %>%
    summarise(
      n = n(),
      min = min(Value, na.rm = TRUE),
      max = max(Value, na.rm = TRUE),
      mean = mean(Value, na.rm = TRUE),
      sd = sd(Value, na.rm = TRUE),
      se = sd / sqrt(n),
      .groups = "drop")}

# Apply it to two different datasets
summary1 <- summarize_data(df_long)
print(summary1, n = 30)


df_summary2 <- summarize_data(_df)

# Print results
print(df_summary1, n = 30)
print(df_summary2, n = 30)


















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

### T-TESTS testing effectiveness of experimental design ###

# Analysis of difference in proportion of swapped or original eggs between treatments. 
t.test(Foreign_percentage ~ Treatment, data = successful_nests)

# Analysis of difference in final manipulated clutch size between treatments. 
t.test(Manipulated_clutch_size ~ Treatment, data = successful_nests)

# Analysis of difference in observed nesting period between treatments. 
t.test(Observed_nesting_period ~ Treatment, data = successful_nests)

# Hatch spread analysis adjusted to only include nests that hatched more than one egg. 
successful_nests_HS <- successful_nests %>%
  filter(True_hatch_spread >1)
t.test(True_hatch_spread ~ Treatment, data = successful_nests_HS)

# Transferred time analysis adjusted to exclude nest that had egg manipulation
# occur after hatching had begun.
successful_nests_MT <- successful_nests %>%
  filter(Transfered_time > 0)
t.test(Transfered_time ~ Treatment, data = successful_nests_MT)

# Faceted figure with all numeric things. Can filter it down in the df_long_filtered
# data set to only do a handful of panels. 
# Can also do this for a single figure instead of faceted. 
df_long_filtered <- subset(df_long, Variable != "Number_transfers")

ggplot(df_long_filtered, aes(x = Treatment, y = as.numeric(Value), 
                             colour = Treatment)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.1, height = 0), 
              size = 3, alpha = 0.7) +
  facet_wrap(~Variable, scales="free_y") +
  guides(color = "none") +
  theme_classic() +
  scale_color_manual(values = c("Asynchronous" = "darkblue", 
                                "Synchronous" = "firebrick2")) +
  labs(x = NULL, y = NULL) 

### PREDICTION A ###: Synch nests have greater hatching success than Asynch nests
# as there are fewer eggs deserted in the nest after hatching. 

# T-test analyzing difference in hatching success between treatments
t.test(Hatch_success ~ Treatment, data = successful_nests)

# POWER ANALYSES

# Cohen's D effect size calculation
cohens_d(Hatch_success ~ Treatment, data = successful_nests)
cohens_d(Survival_60 ~ Treatment, data = successful_nests)

# Hedge's G effect size calculation. This is better for my small sample size
SynchHatch <- successful_nests$Hatch_success[successful_nests$Treatment == "Synchronous"]
AsynchHatch <- successful_nests$Hatch_success[successful_nests$Treatment == "Asynchronous"]

hedges_g <- cohen.d(AsynchHatch, SynchHatch, hedges.correction = TRUE)
print(hedges_g)

# Setting the parameters for power analysis
sample_size <-            # Sample size
  alpha <- 0.05             # Significance level
power <- 0.8              # Power
effect_size <- hedges_g   # Hedge's G effect size  

# Calculating effect size needed
## Figure out 2 side or one sided analysis ##
effect_size <- pwr.p.test(n = sample_size, sig.level = alpha, 
                          power = power, alternative = "two.sided")$h
print(effect_size)

# Calculating sample size needed
g_sample_size <- pwr.p.test(h = effect_size, sig.level = alpha, 
                            power = power, alternative = "less")$n
print(g_sample_size)

# Model analyzing number of hatched eggs by clutch size
model_3a <- glmmTMB(Hatched ~ Manipulated_clutch_size +
                     (1|Year), family = poisson, data = successful_nests)
check_model(model_3a)
summary(model_3a)

# Model analyzing number of hatched eggs by proportion of foreign eggs
model_3b <- glmmTMB(Hatched ~ Foreign_percentage +
                     (1|Year), family = poisson, data = successful_nests)
check_model(model_3b)
summary(model_3c)

# Model analyzing number of hatched eggs by treatment
model_3c <- glmmTMB(Hatched ~ Treatment +
                     (1|Year), family = poisson, data = successful_nests)
check_model(model_3c)
summary(model_3c)


### PREDICTION B ###: Synchronous nests have greater survival than Asynchronous
# nests as there are fewer late hatching chicks with poorer survival.

# Model analyzing number of survivors by treatment.
model_3d <- glmmTMB(Survival_60 ~ Treatment +
                     (1|Year), family = poisson, data = successful_nests)
check_model(model_3d)
summary(model_3d)

# Model analyzing number of survivors by clutch size.
model_3e <- glmmTMB(Survival_60 ~ Manipulated_clutch_size +
                     (1|Year), family = poisson, data = successful_nests)
check_model(model_3e)
summary(model_3e)

# Model analyzing number of survivors by number of hatched eggs. 
model_3f <- glmmTMB(Survival_60 ~ Hatched +
                     (1|Year), family = poisson, data = successful_nests)
check_model(model_3f)
summary(model_3f)

# Survival to 60 days pivoted longer format for the figure.
data_long_format <- successful_nests %>%
  mutate(NonSurvivors = Hatched - Survival_60) %>% # Calculate non-survivors
  pivot_longer(cols = c(Survival_60, NonSurvivors), 
               names_to = "Status", 
               values_to = "Count") %>%
  mutate(Survival_60 = ifelse(Status == "Survival_60", 1, 0)) %>%
  uncount(Count) %>%
  mutate(Treatment_survival = paste(Treatment, Survival_60, sep = "_")) 

# Bar plot of survival data.
ggplot(data_long_format, aes(x = Treatment, fill = factor(Treatment_survival)))+
  geom_bar(position = "dodge") +
  labs(x = "Treatment", y = "Count", fill = "Survived") +
  scale_fill_manual(values = c("Asynchronous_1" = "darkblue", 
                               "Asynchronous_0" = "blue", 
                               "Synchronous_1" = "firebrick4", 
                               "Synchronous_0" = "firebrick1"), 
                    guide = "none") +
  theme_classic() 

### Supplementary Chapter 3 analyses ###

# Model analyzing the hatch spread by clutch size
model_1 <- glmmTMB(True_hatch_spread ~ Manipulated_clutch_size + (1|Year),
                   family = poisson, data = successful_nests)
check_model(model_1)
summary(model_1)
