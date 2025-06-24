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

############################
####      CHAPTER 2     ####
#### OBSERVATIONAL DATA ####
############################

#### CH.2 INTRO AND HYPOTHESES ####
# I am investigating nesting factors in Pūkeko (Porphyrio melanotus melanotus) 
# to explore the effects of hatching spread (HS) on hatching success and 
# survival. The data for chapter 2 has been compiled from field books and 
# spreadsheets to create master data sets with both summary data for all 
# nests along with specific hatching data for individual nests. This data is
# from the years 2008, 2010, 2013, 2014-2018, 2022-2024. 

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

masterlist_data <- read_excel("Nests_masterlist.xlsx",
                              ## BMB: is it OK to treat these all as NA?
                              na = c("", "NO_RECORD", "MISSING")) %>%
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>%  
  filter(Exclusion == "GOOD" | Exclusion == "MISSING_HATCH_ORDER") %>%
  mutate(Hatch_success = as.numeric(Hatched_eggs)/as.numeric(Clutch_size), 
         Survived_2_cons = ifelse(is.na(Survived_2), 0, Survived_2), 
         Hatch_begin = as.numeric(Hatch_begin),
         Hatch_end = as.numeric(Hatch_end),
         Hatch_spread = (Hatch_end - Hatch_begin + 1), 
         Hatch_begin = as.numeric(Hatch_begin), 
         Date_found = as.numeric(Date_found), 
         Observed_nesting_period = as.numeric(Hatch_begin - Date_found)) %>%
  mutate(across(c(Clutch_size, Hatched_eggs, Survived_2, Survived_2_cons,
         Hatch_begin, Hatch_end, Hatch_spread, Year), 
         ~ as.numeric(as.character(.)))) %>%
  mutate(Survived_2_cons = ifelse(is.na(Survived_2), 0, Survived_2),
        converted_from_na = is.na(Survived_2))
View(masterlist_data)




#### HYPOTHESIS 1 ####
# Dominant females lay synchronously with another female to increase 
# their inclusive fitness. 

longer_data1 <- masterlist_data %>%
  pivot_longer(
    cols = c("Hatched_eggs", "Hatch_success", "Survived_2_cons"),
    names_to = "Data",
    values_to = "Value") %>%
  mutate(converted_from_na = ifelse(Data == "Survived_2_cons", 
                                    converted_from_na, FALSE))

# Figure of Hatch rate, number of hatched eggs and number of survivors by 
# Hatching spread. Two measures of survival: When survival not known converted
# to zeros as a conservative estimate, and then the NAs removed as a more 
# liberal estimate of survival. 
ggplot(longer_data1, aes(x = Hatch_spread, y = Value, color = converted_from_na)) +
  geom_jitter(width = 0.2, height = 0.05, alpha = 0.8) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  facet_wrap(~Data, scales = "free_y", labeller = labeller(
    Data = c(Hatched_eggs = "Number of hatched eggs",
             Hatch_success = "Hatching proportion", 
             Survived_2_cons = "Number of survivors where red 
             = NAs converted to 0s"))) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Hatch Spread (days)", y = NULL)

#### PREDICTION 1a: Shorter HS have greater hatching success ####

## Model analyzing the hatching rate by hatching spread. 
## I don't know how to analyze this. Obviously I could do just the number of 
## hatched eggs but then I think I'll need to add clutch size as a co-variate?

## BMB: it should be a binomial (or beta-binomial) model with clutch size as
## the weight, i.e.
## glmmmTMB(prop_survived ~ ..., weights = clutch_size, ...)

#### PREDICTION 1b: Shorter HS have better survival ####

# Model analyzing survival by HS with conservative survival estimate
model_2_1b_cons <- glmmTMB(Survived_2_cons ~ Hatch_spread + (1|Year),
                    family = poisson(link = log), data = masterlist_data)
check_model(model_2_1b_cons)
res <- simulateResiduals(fittedModel = model_2_1b_cons, plot = TRUE)
testDispersion(res)
testZeroInflation(res)
## BMB: surprising that the heteroscedasticity in check_models homogeneity-of-variance
## (scale-location) plot doesn't get flagged by DHARMa
## I'm a little surprised at treating Survived_2_cons as Poisson rather than as
##  something like a binomial.  Do you have the denominators (i.e., the total
##  possible survivors for each case)?

# Model analyzing survival by HS with liberal survival estimate
lib_data <- masterlist_data %>%
  filter(!is.na(Survived_2))
model_2_1b_lib <- glmmTMB(Survived_2 ~ Hatch_spread + (1|Year) 
                   ,family = poisson(link = log), data = lib_data)
check_model(model_2_1b_lib)
res <- simulateResiduals(fittedModel = model_2_1b_lib, plot = TRUE)
testDispersion(res)


#### HYPOTHESIS 2 ####
# Females lay together because the optimal clutch size is greater 
# than the size of single female clutches.

# Figure of the hatching rate, number of hatched eggs, and number of survivors
# with bith survival estimates by clutch size. 
ggplot(longer_data1, aes(x = Clutch_size, y = Value, color = converted_from_na)) +
  geom_jitter(width = 0.2, height = 0.05, alpha = 0.8) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  facet_wrap(~Data, scales = "free_y", labeller = labeller(
    Data = c(Hatched_eggs = "Number of hatched eggs",
             Hatch_success = "Hatching proportion", 
             Survived_2_cons = "Number of survivors where red 
             = NAs converted to 0s"))) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Clutch size", y = NULL)

# Pivoting data longer with both survival estimates
longer_data2 <- masterlist_data %>%
  pivot_longer(cols = c('Hatched_eggs', 'Hatch_success',
                        'Survived_2', 'Survived_2_cons'), 
               names_to = 'Data', 
               values_to = 'Value')

# Separate figure that has both survival estimates in separate panels instead.
# Faceted figure of hatching success, hatched eggs, and conservative and liberal
# estimates of survival by clutch size.
ggplot(longer_data2, aes(x = Clutch_size, y = Value)) +
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

# PREDICTION 2a: Number of hatched eggs will increase with clutch size 
# then begin to decrease when the clutch is too large to support. 

# Model analyzing number of hatched eggs by clutch size
model_2_2a <- glmmTMB(Hatched_eggs ~ I(Clutch_size^2) +
                          (1|Year), family = poisson(link = "log"), 
                        data = masterlist_data)
check_model(model_2_2a)
res <- simulateResiduals(fittedModel = model_2_2a, plot = TRUE)
testDispersion(res)

## BMB: I would suggest something like this ...
bmb1 <- glmmTMB(Hatched_eggs/Clutch_size ~ Clutch_size + (1|Year),
                weights = Clutch_size,
                family = betabinomial, data = masterlist_data)
check_model(bmb1)
## BMB: the posterior predictive check is weird but I suspect that's something
## wonky with what the package is doing, not a real problem with the model
## It would be interesting/useful to dig in and see what's going on and/or post
## a reproducible issue to whichever easystats package github site is relevant
##
## I think it can't handle the model when specified as proportion + weights --
## yes, that's it ...
performance::check_predictions(bmb1)

## re-fit model with c(successes, failures) instead of prop + weights
bmb2 <- update(bmb1, cbind(Hatched_eggs, Clutch_size -Hatched_eggs) ~ .,
                weights = NULL)
performance::check_predictions(bmb2)

res <- simulateResiduals(bmb1, plot = TRUE)

# PREDICTION 2b: More hatched eggs have more survivors regardless of hatch spread.

# Conservative model analyzing number of survivors by clutch size and 
# the number of hatched eggs
model_2_2b_cons <- glmmTMB(Survived_2_cons ~ Hatched_eggs*Hatch_spread +
                        (1|Year), family = poisson(link = "log"), 
                        data = masterlist_data)
check_model(model_2_2b_cons)
res <- simulateResiduals(fittedModel = model_2_2b_cons, plot = TRUE)
testDispersion(res)
testZeroInflation(res) 

# Liberal model analyzing number of survivors by clutch size and # of hatched eggs
model_2_2b_lib <- glmmTMB(Survived_2 ~ Hatched_eggs +
                          (1|Year), family = poisson, data = lib_data)
check_model(model_2_2b_lib)
res <- simulateResiduals(fittedModel = model_2_2b_lib, plot = TRUE)
testDispersion(res)
testZeroInflation(res) 


#### HYPOTHESIS 3 #### Hatch order less important in synchronous nests

# Model of laying order on hatching order 
# Still working on finalizing the data for this. 
# model_d <- glmmTMB(Hatch_order ~ Lay_order + (Nest_ID/Year), 
                    # family = gaussian, data = chick_data)
# res <- simulateResiduals(fittedModel = model_2_3c, plot = TRUE)
# testDispersion(res)
# check_model(model_d)

# Model of laying spread on hatching spread
# Probably won't have enough data for this but will try. 
# model_d <- glmmTMB(Hatch_spread ~ Lay_spread +
          # (Nest_ID/Year), family = poisson, data = chick_data)
# res <- simulateResiduals(fittedModel = model_2_3c, plot = TRUE)
# testDispersion(res)
# check_model(model_d)


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

# Model with hatch order as sole predictor variable. 
# Same as Cody Dey's paper (Dey et al. 2014). 
# Cody did nest_ID as the random effect. If I do Nest_ID nested within year 
# I think it causes too many issues. 

# I think I should also only keep nests that hatch more than one egg. 
chick_data_survival <- chick_data_survival %>%
  filter(Hatched_eggs > 1)

View(chick_data_survival)
model_2_3a <- glmmTMB(Survived ~ Hatch_order +
                  (1|Year), family = binomial, data = chick_data_survival)
check_model(model_2_3a)
res <- simulateResiduals(fittedModel = model_2_3a, plot = TRUE)
testDispersion(res)

model_2_3a <- glmmTMB(Survived ~ Hatch_order +
                     (1|Nest_ID), family = binomial, data = chick_data_survival)
check_model(model_2_3a)
res <- simulateResiduals(fittedModel = model_2_3a, plot = TRUE)
testDispersion(res)

# Figure of survival by hatching order. 
ggplot(chick_data_survival, aes(x = Hatch_order, y = as.numeric(Survived), colour = Year)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.7) +
  labs(y = "Survived", x = "Hatch_order") +
  theme_classic()
  

# PREDICTION 3b: Hatch order more important in larger clutches
# Model with interaction between hatch order clutch size
model_2_3b <- glmmTMB(Survived ~ Hatch_order*Clutch_size +
                     (Nest_ID/Year), family = binomial, data = chick_data_survival)
check_model(model_2_3b)
res <- simulateResiduals(fittedModel = model_2_3b, plot = TRUE)
testDispersion(res)

# PREDICTION 3c: Hatch order more important for survival in asynchronous clutches.
# Model with interaction between hatch order and hatching spread
model_2_3c <- glmmTMB(Survived ~ Hatch_order*Hatch_spread +
                     (1|Nest_ID), family = binomial, data = chick_data)
check_model(model_2_3c)
res <- simulateResiduals(fittedModel = model_2_3c, plot = TRUE)
testDispersion(res)


# PREDICTION 3d: Late hatching eggs will be smaller than early hatching eggs
# because they spend more developmental energy to catch up. 

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
# Tarsus measurements switched protocols over the years, will deal with that.
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


### Supplementary Chapter 2 analyses ###
### Most descriptive model to break it ###
# Model with all fixed effects as random effects to break and have singular fit. 
# model_a <- glmmTMB(as.numeric(Hatched_eggs) ~ as.numeric(Clutch_size) + 
                  #   as.numeric(Group_size) + as.numeric(Hatch_spread) +
                  #   (1|Year) + (1|as.numeric(Clutch_size)) + 
                  #   (1|as.numeric(Group_size)) +
                  #   (1|as.numeric(Hatch_spread)), family = poisson(link = log), 
                  #   data = masterlist_data)
# check_model(model_a)
# summary(model_a)




##############################
####      CHAPTER 3       ####
#### SYNCHRONY EXPERIMENT ####
##############################

#### CH.3 INTRO AND HYPOTHESIS ####
# We experimentally lengthened and shortened the HS of nests Pūkeko by using 
# known laying dates, candling images, and flotation scores in 2018, 2023-24. 

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


#### PREDICTION A ####: Synchronous nests have greater hatching success than 
# Asynchronous nests as there are fewer eggs deserted in the nest after hatching. 

# T-test analyzing difference in number of hatched eggs between treatments.
t.test(Hatched ~ Treatment, data = successful_nests)

# T-test analyzing difference in hatching success between treatments.
t.test(Hatch_success ~ Treatment, data = successful_nests)

# POWER ANALYSIS
# How large would sample size need to be given observed effect size?

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


#### PREDICTION B ####: Synchronous nests have greater survival than 
# Asynchronous nests as there are fewer late hatching chicks with poorer survival.

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
survival_data_long <- successful_nests %>%
  mutate(NonSurvivors = Hatched - Survival_60) %>% # Calculate non-survivors
  pivot_longer(cols = c(Survival_60, NonSurvivors), 
               names_to = "Status", 
               values_to = "Count") %>%
  mutate(Survival_60 = ifelse(Status == "Survival_60", 1, 0)) %>%
  uncount(Count) %>%
  mutate(Treatment_survival = paste(Treatment, Survival_60, sep = "_")) 

# Bar plot of survival data.
# Legend is annoying to make so I will add it separately. 
ggplot(survival_data_long, aes(x = Treatment, fill=factor(Treatment_survival)))+
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
