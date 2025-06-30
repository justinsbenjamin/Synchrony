#########################
#### JUSTIN BENJAMIN ####
#### MSC THESIS CODE ####
#########################

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
library(broom.mixed)

sample_data <- read_excel("Book1.xlsx") %>%
  mutate(Date_found = as.Date(Date_found),
         Hatch_begin = as.Date(Hatch_begin),
         Hatch_end = as.Date(Hatch_end),
         Last_observed = as.Date(Last_observed)) %>%
  mutate(Hatch_success = (Hatched_eggs/Clutch_size*100), 
         Hatch_spread = as.numeric(Hatch_end - Hatch_begin + 1)) %>%
    mutate(Observed_nesting_period = if_else(
        !is.na(Hatch_begin),
        as.numeric(Hatch_begin - Date_found),
        as.numeric(Last_observed - Date_found)))
View(sample_data)

ggplot(sample_data, aes(x = Observed_nesting_period, y = Hatch_success)) +
  geom_point() +
  theme_classic()



############################
####      CHAPTER 2     ####
#### OBSERVATIONAL DATA ####
############################

#### HYPOTHESES & PREDICTIONS ####

# HYPOTHESIS 1: Dominant females lay synchronously with another female to 
# increase their inclusive fitness. 
# PREDICTION 1a: Shorter HS have better hatching success
# PREDICTION 1b: Shorter HS have better survival

# HYPOTHESIS 2: Females lay together because the optimal clutch size is 
# greater than the size of single female clutches regardless of hatch spread.
# PREDICTION 2a: # of hatched eggs will increase with clutch size 
# then begin to decrease when clutch is too large to support.
# PREDICTION 2b: Larger broods have more survivors. 
# PREDICTION 2c: Larger clutches have greater hatching spread. 

# HYPOTHESIS 3: Subordinate females benefit more from synchronous nests
# PREDICTION 3a: Hatch order less important for survival with short HS
# PREDICTION 3b: Last hatched eggs are smaller as they spend more developmental
# energy catching up to early hatching chicks. 

#### DATA FILTERING AND ORGANIZATION ####
# Sample sizes of excluded and included groups
nest_counts <- read_excel("Nests_masterlist.xlsx") %>%
  count(Exclusion, name = "nest_counts") %>%
  arrange(desc("nest_counts"))
print(nest_counts)

observation_period_data <- read_excel("Nests_masterlist.xlsx") %>%
  filter(Exclusion == "GOOD" | Exclusion == "MISSING_HATCH_ORDER" |
         Exclusion == "FAILED") %>%
  mutate(Date_found = as.Date(Date_found),
         Hatch_begin = as.Date(Hatch_begin),
         Hatch_end = as.Date(Hatch_end),
         Last_observed = as.Date(Last_observed)) %>%
  mutate(Hatch_success = (as.numeric(Hatched_eggs)/as.numeric(Clutch_size)*100), 
         Hatch_spread = as.numeric(Hatch_end - Hatch_begin + 1)) %>%
  mutate(Observed_nesting_period = if_else(
    !is.na(Hatch_begin),
    as.numeric(Hatch_begin - Date_found),
    as.numeric(Last_observed - Date_found)))
View(observation_period_data)

ggplot(observation_period_data, aes(x = Observed_nesting_period, y = Hatch_success)) +
  geom_jitter(width = 0.2, height = 0.5, alpha = 0.8) +
  theme_classic()





masterlist_data <- read_excel("Nests_masterlist.xlsx",
  na = c("", "NO_RECORD", "MISSING")) %>%
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>%  
  filter(Exclusion == "GOOD" | Exclusion == "MISSING_HATCH_ORDER") %>%
  mutate(Hatch_success = (as.numeric(Hatched_eggs)/as.numeric(Clutch_size)*100), 
         Survived_2_cons = ifelse(is.na(Survived_2), 0, Survived_2), 
         Hatch_begin = as.Date(as.numeric(Hatch_begin)),
         Hatch_end = as.Date(Hatch_end),
         Hatch_spread = as.numeric(Hatch_end - Hatch_begin + 1), 
         Date_found = as.Date(Date_found), 
         Last_observed = as.Date(Last_observed), 
         Observed_nesting_period = as.numeric(Hatch_begin - Date_found)) %>%
  mutate(across(c(Clutch_size, Hatched_eggs, Survived_2, Survived_2_cons,
         Hatch_begin, Hatch_end, Hatch_spread, Year), 
         ~ as.numeric(as.character(.)))) %>%
  mutate(Survived_2_cons = ifelse(is.na(Survived_2), 0, Survived_2),
        converted_from_na = is.na(Survived_2))
View(masterlist_data)

#### HYPOTHESES 1 & 2

ggplot(masterlist_data, aes(x = Clutch_size, y = Hatch_spread)) +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.8) +
  theme_classic()

ggplot(masterlist_data, aes(x = Hatched_eggs, y = Hatch_spread)) +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.8) +
  theme_classic()

ggplot(masterlist_data, aes(x = Observed_nesting_period, y = Hatch_success)) +
  geom_point() +
  theme_classic()

# Pivoting data longer with both survival estimates
# Two measures of survival: When survival not known, the NAs were converted
# to zeros as a conservative estimate, and then as a more liberal estimate of 
# survival, the NAs were removed. 
longer_data1 <- masterlist_data %>%
  pivot_longer(cols = c('Hatched_eggs', 'Hatch_success',
                        'Survived_2', 'Survived_2_cons'), 
               names_to = 'Data', 
               values_to = 'Value') %>%
  mutate(converted_from_na = ifelse(Data == "Survived_2_cons", 
                                    converted_from_na, FALSE)) %>%
  pivot_longer(
    cols = c("Hatch_spread", "Clutch_size"), 
    names_to = "Predictor", 
    values_to = "Predictor_value")

# Figure of Hatch rate, Number of hatched eggs and Number of survivors with 
# conservative and liberal estimates by HS and Clutch size. 
figure1 <- ggplot(longer_data1, aes(x = Predictor_value, y = Value)) +
  geom_jitter(width = 0.2, height = 0.05, alpha = 0.8) +
  facet_grid(Data ~ Predictor, scales = "free", labeller = labeller(
    Data = c(Hatched_eggs = "# of hatched eggs",
             Hatch_success = "Hatching rate (%)", 
             Survived_2 = "Liberal # of survivors",
             Survived_2_cons = "Conservative # of survivors"), 
    Predictor = c(Clutch_size = "Clutch size", 
                  Hatch_spread = "Hatch spread (days)")), switch = "both") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  theme(strip.background = element_blank(),
        strip.placement = "outside", 
        strip.text = element_text(size = 8))
print(figure1)

# Added colours to the NAs that were converted to 0s
figure2 <- figure1 + aes(x = Predictor_value, y = Value, color = converted_from_na) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"))
print(figure2)

# Added a smoother function to see the distribution of the data. 
figure3 <- figure1 +
  geom_smooth(method = "loess") 
print(figure3)

# Model analyzing proportion of hatched egg by clutch size and HS
model_1 <- glmmTMB(cbind(Hatched_eggs, Clutch_size - Hatched_eggs) ~ 
                  Hatch_spread + Clutch_size + (1|Year), family = betabinomial,
                  data = masterlist_data)

# Function to view model diagnostics 
diagnostics <- function(model) {
print(check_model(model)) 
res <- simulateResiduals(fittedModel = model, plot = TRUE)
testDispersion(res)
testZeroInflation(res)}

diagnostics(model_1)

# Model analyzing proportion of survivors (conservative estimate) by clutch size and HS
model_2 <- update(model1, formula = cbind(Survived_2_cons, Hatched_eggs - Survived_2_cons) ~ 
                    Hatch_spread + Hatched_eggs + Clutch_size + (1|Year))
diagnostics(model_2)

# Model analyzing proportion of survivors (conservative estimate) by clutch size and HS
lib_data <- masterlist_data %>% filter(!is.na(Survived_2))

model_3 <- update(model_2, formula = cbind(Survived_2, 
                         Hatched_eggs - Survived_2) ~ Hatch_spread + Hatched_eggs
                         + Clutch_size + (1|Year),
                         data = lib_data)
diagnostics(model_3)

#### HYPOTHESIS 3 

# Model of laying order on hatching order 
# Still working on finalizing the data for this. 
# model_d <- glmmTMB(Hatch_order ~ Lay_order + (Nest_ID/Year), 
                    # family = gaussian, data = chick_data)
# diagnostics(model_d)

# Model of laying spread on hatching spread
# Probably won't have enough data for this but will try. 
# model_d <- glmmTMB(Hatch_spread ~ Lay_spread +
          # (Nest_ID/Year), family = poisson, data = chick_data)
# diagnostics(model_d)

#### PREDICTION 3a: Earlier hatched eggs have greater survival (Dey et al. 2014).

chick_data <- read_excel("Final_Compiled_Chick_Data.xlsx") %>%
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_"))

chick_data$Clutch_size <- masterlist_data$Clutch_size[match(chick_data$Nest_ID,
                                                      masterlist_data$Nest_ID)]
chick_data$Hatch_spread <- masterlist_data$Hatch_spread[match(chick_data$Nest_ID, 
                                                       masterlist_data$Nest_ID)]
chick_data$Hatched_eggs <- masterlist_data$Hatched_eggs[match(chick_data$Nest_ID, 
                                                        masterlist_data$Nest_ID)]
chick_data_survival <- chick_data %>%
  filter(!is.na(Survived) & Survived %in% c(0, 1)) %>%
  mutate(Survived = if (!is.numeric(Survived)) as.numeric(Survived) else Survived)

chick_data_survival <- chick_data_survival %>%
  filter(Hatched_eggs > 1)

# Figure of survival by hatching order. 
ggplot(chick_data_survival, aes(x = Hatch_order, y = as.character(Survived), colour = Year)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.7) +
  labs(y = "Survived", x = "Hatch_order") +
  theme_classic()

# Model with hatch order as sole predictor variable on survival same as 
# Cody Dey's paper (Dey et al. 2014). 
model_2_3a <- glmmTMB(as.numeric(Survived) ~ as.numeric(Hatch_order) + Clutch_size +
                  (1|Nest_ID) + (1|Year), family = binomial, data = chick_data_survival)
diagnostics(model_2_3a)

# Model with interaction between hatch order clutch size
model_2_3b <- glmmTMB(as.numeric(Survived) ~ as.numeric(Hatch_order)*Clutch_size +
                     (1|Nest_ID) + (1|Year), family = binomial, data = chick_data_survival)
diagnostics(model_2_3b)

# Model with interaction between hatch order and hatching spread
model_2_3c <- glmmTMB(as.numeric(Survived) ~ as.numeric(Hatch_order)*Hatch_spread + Clutch_size +
                     (1|Nest_ID) + (1|Year), family = binomial, data = chick_data_survival)
diagnostics(model_2_3c)

# Filtering out data that has unusually large birds (likely typos or issues with
# the measurements and potential recaptures) then pivoting data longer. 
chick_data_longer <- chick_data %>%
  filter(Mass < 33) %>%
  filter(`Shield to Tip` < 24) %>%
  pivot_longer(cols = c('Mass', 'Tarsus',
                        'Shield to Tip'), 
               names_to = 'Morphometrics', 
               values_to = 'Value') %>%
  filter(!is.na(Value))

# Figure of size by hatch order
ggplot(chick_data_longer, aes(x = Hatch_order, y = Value, colour = Year)) +
  geom_jitter(width = 0.25, height = 0, size = 2, alpha = 0.7) +
  facet_wrap(~Morphometrics, scales="free_y", labeller = labeller(
    Morphometrics = c(Mass = "Mass (g)",
                      Tarsus = "Left outer tarsus (mm)", 
                      `Shield to Tip` = "Shield to tip (mm)"))) +
  scale_y_continuous() +
  labs(x = "Hatch order", y = "Measurement") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

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

##############################
####      CHAPTER 3       ####
#### SYNCHRONY EXPERIMENT ####
##############################

#### HYPOTHESIS & PREDICTIONS ####

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

# Nests removed if only hatched a single egg
successful_nests_HS <- successful_nests %>%
  filter(Hatched >1)

# Nest removed since egg manipulation after hatching began
successful_nests_MT <- successful_nests %>%
  filter(Transfered_time > 0)

# Function to pivot data long
pivot_data_long <- function(df) {
  df %>%
    pivot_longer(
      cols = where(is.numeric) & !all_of("Year"),
      names_to = "Variable",
      values_to = "Value")}

# Run pivoting long function 
experiment_data_long <- pivot_data_long(experiment_data)
successful_nests_long <- pivot_data_long(successful_nests)
successful_nests_long_HS <- pivot_data_long(successful_nests_HS)
successful_nests_long_MT <- pivot_data_long(successful_nests_MT)

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
      .groups = "drop")}

summary_experiment_data <- summarize_data(experiment_data_long)
summary_successful <- summarize_data(successful_nests_long)
summary_successful_HS <- summarize_data(successful_nests_long_HS)
summary_successful_MT <- summarize_data(successful_nests_long_MT)

lapply(list(summary_experiment_data, summary_successful, 
            summary_successful_HS, summary_successful_MT), print, n = 30)

#### T-TESTS analyzing the effectiveness of the experimental design ####

# Analysis of difference in proportion of swapped eggs between treatments. 
t.test(Foreign_percentage ~ Treatment, data = experiment_data)

# Analysis of difference in final manipulated clutch size between treatments. 
t.test(Manipulated_clutch_size ~ Treatment, data = experiment_data)

# Analysis of difference in observed nesting period between treatments. 
t.test(Observed_nesting_period ~ Treatment, data = successful_nests)

# Hatch spread analysis adjusted to only include nests that hatched >1 egg. 
t.test(True_hatch_spread ~ Treatment, data = successful_nests_HS)

# Transferred time analysis adjusted to exclude nest that had egg manipulation
# occur after hatching had begun.
t.test(Transfered_time ~ Treatment, data = successful_nests_MT)

# T-test analyzing difference in number of hatched eggs between treatments.
t.test(Hatched ~ Treatment, data = successful_nests)

# T-test analyzing difference in hatching success between treatments.
t.test(Hatch_success ~ Treatment, data = successful_nests)

# POWER ANALYSIS

# Cohen's D effect size calculation
cohens_d(Hatch_success ~ Treatment, data = successful_nests)

# Hedge's G effect size calculation. This is better for my small sample size
SynchHatch <- successful_nests$Hatch_success[successful_nests$Treatment == "Synchronous"]
AsynchHatch <- successful_nests$Hatch_success[successful_nests$Treatment == "Asynchronous"]

hedges_g <- cohen.d(AsynchHatch, SynchHatch, hedges.correction = TRUE)
print(hedges_g)

# Setting the parameters for power analysis
alpha <- 0.05             # Significance level
power <- 0.8              # Power
effect_size <- hedges_g   # Hedge's G effect size  

# Calculating sample size needed
g_sample_size <- pwr.p.test(h = effect_size, sig.level = alpha, 
                            power = power, alternative = "less")$n
print(g_sample_size)

# GLMMs on hatching success

# Model analyzing number of hatched eggs by clutch size
model_1 <- glmmTMB(cbind(Hatched, Manipulated_clutch_size - Hatched) ~ 
                     Manipulated_clutch_size + (1|Year), family = betabinomial, 
                     data = successful_nests)
diagnostics(model_1)

# Model analyzing proportion of hatched chicks by number of swapped eggs
model_2 <- glmmTMB(cbind(Hatched, Manipulated_clutch_size - Hatched) ~ 
                     Foreign_percentage + (1|Year), family = betabinomial, 
                     data = successful_nests)
diagnostics(model_2)

# Model analyzing number of hatched eggs by treatment
model_2 <- glmmTMB(cbind(Hatched, Manipulated_clutch_size - Hatched) ~ 
                     Treatment + (1|Year), family = betabinomial, 
                   data = successful_nests)
diagnostics(model_2)

# Model analyzing number of survivors by treatment.
model_1 <- glmmTMB(cbind(Survival_60, Hatched - Survival_60) ~ 
                     Treatment + Hatched + (1|Year), 
                   family = betabinomial, data = successful_nests)
diagnostics(model_1)
summary(model_1)

# Survival to 60 days pivoted longer format for the figure.
survival_data_long <- successful_nests %>%
  mutate(NonSurvivors = Hatched - Survival_60) %>% 
  pivot_longer(cols = c(Survival_60, NonSurvivors), 
               names_to = "Status", 
               values_to = "Count") %>%
  mutate(Survival_60 = ifelse(Status == "Survival_60", 1, 0)) %>%
  uncount(Count) %>%
  mutate(Treatment_survival = paste(Treatment, Survival_60, sep = "_")) 

# Bar plot of survival data.
ggplot(survival_data_long, aes(x = Treatment, fill=factor(Treatment_survival)))+
  geom_bar(position = "dodge") +
  labs(x = "Treatment", y = "Count", fill = "Survived") +
  scale_fill_manual(values = c("Asynchronous_1" = "darkblue", 
                               "Asynchronous_0" = "blue", 
                               "Synchronous_1" = "firebrick4", 
                               "Synchronous_0" = "firebrick1"), guide = "none") +
  theme_classic() 

### Supplementary Chapter 3 analyses ###

# Model analyzing the hatch spread by clutch size
model_1 <- glmmTMB(True_hatch_spread ~ Manipulated_clutch_size + (1|Year),
                   family = poisson, data = successful_nests)
check_model(model_1)
summary(model_1)