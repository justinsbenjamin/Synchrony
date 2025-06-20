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
library(ggbeeswarm)
library(pwr)
library(rstatix)
library(glmmTMB)
library(effsize)
library(DHARMa)

########### CHAPTER 2 ##########
####   OBSERVATIONAL DATA   ####
################################

#### CH.2 INTRO AND HYPOTHESES ####
# I am investigating nesting factors in Pūkeko (Porphyrio melanotus melanotus).
# The data for chapter 2 has been compiled from field books and spreadsheets to create a 
# master dataset with both summary data for all nests along with specific 
# hatching data for individual nests. Data is from the years 2008, 2010, 2013, 
# 2014-2018, 2022-2024. 

# HYPOTHESIS 1: Dominant females lay synchronously with another female to increase 
# their inclusive fitness. 

# PREDICTION 1a: Shorter hatching spread (HS) have better survival
# PREDICTION 1b: Hatch order less important in synchronous nests 
# PREDICTION 1c: 

# HYPOTHESIS 2: Females lay together because the optimal clutch size is greater 
# than the size of single female clutches regardless of hatch spread.

# PREDICTION 2a: Larger clutches hatch more eggs 
# PREDICTION 2b: Larger broods have more survivors
# PREDICTION 2c: Larger clutches have greater hatching spread 

#### DATA FILTERING AND ORGANIZATION ####

nest_counts <- read_excel("Nests_masterlist.xlsx") %>%
count(Exclusion, name = "nest_counts") %>%
arrange(desc("nest_counts"))
print(nest_counts)

masterlist_data <- read_excel("Nests_masterlist.xlsx") %>%
  filter(Exclusion == "GOOD") %>% # Keeping only Good nests that meet criteria
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>% # Unique nest ID 
  select(-`...22`, -`...23`, -`...24`) %>%
  mutate(Hatch_success = as.numeric(Hatched_eggs)/as.numeric(Clutch_size)) %>%
  mutate(Survived_2_cons = ifelse(is.na(Survived_2), 0, Survived_2)) 
# 2 estimates of survival, one where NAs are 0 survived and another of just 
# removing the NAs. 
View(masterlist_data)


# mutate(Females = ifelse(Clutch_size <= 5, "Single Female", "Joint Females"))


#### HYPOTHESIS 1 ####
# Dominant females lay synchronously with another female to increase 
# their inclusive fitness. 

masterlist_data <- read_excel("Nests_masterlist.xlsx") %>%
  filter(Exclusion == "GOOD") %>% # Keeping only Good nests that meet criteria
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>% # Unique nest ID 
  select(-`...22`, -`...23`, -`...24`) %>%
  mutate(Hatch_success = as.numeric(Hatched_eggs)/as.numeric(Clutch_size)) %>%
  mutate(Survived_2_cons = ifelse(is.na(Survived_2), 0, Survived_2)) %>%
  mutate(Hatch_begin = as.numeric(Hatch_begin),
         Hatch_end = as.numeric(Hatch_end),
         Hatch_spread = Hatch_end - Hatch_begin + 1) %>%
mutate(Hatch_begin = as.numeric(Hatch_begin), Date_found = as.numeric(Date_found), 
       Observed_incubation_period = as.numeric(Hatch_begin - Date_found)) %>%
  mutate(across(c(Clutch_size, Hatched_eggs, Survived_2, Survived_2_cons,
                  Hatch_begin, Hatch_end, Hatch_spread, Year), ~ as.numeric(as.character(.))))
View(masterlist_data)


#### PREDICTION 1a: Shorter hatching spread (HS) have better survival ####

### Most descriptive model to break it ###
# Model with all fixed effects as random effects to break and have singular fit. 
model_a <- glmmTMB(as.numeric(Hatched_eggs) ~ as.numeric(Clutch_size) + as.numeric(Group_size) + as.numeric(Hatch_spread) +
                     (1|Year) + (1|as.numeric(Clutch_size)) + (1|as.numeric(Group_size)) +
                     (1|as.numeric(Hatch_spread)), family = poisson, data = masterlist_data)
check_model(model_a)
summary(model_a)

# conservative survival estimate
model_e <- glmmTMB(Survived_2_cons ~ Hatch_spread + (1|Year) 
                   ,family = poisson(link = log), data = masterlist_data)
check_model(model_e)
summary(model_e)

# liberal survival estimate
lib_data <- masterlist_data %>%
  filter(!is.na(Survived_2))
model_a <- glmmTMB(Survived_2 ~ Hatch_spread + (1|Year) 
                   ,family = poisson, data = lib_data)
check_model(model_a)
summary(model_a)

# Facetted figure of hatching success, hatched eggs, and conservative and liberal
# estimates of survival by hatch spread

masterlist_data <- masterlist_data %>%
  mutate(Females = ifelse(Clutch_size <= 5, "Single Female", "Joint Females")) %>%
  mutate(Survived_2 = ifelse(Females == "Joint Females", Survived_2 / 2, Survived_2))
  

ggplot(longer_data, aes(x = Hatch_spread, y = Value, colour = Females)) +
  geom_point(alpha = 0.6) +  
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 2, alpha = 0.7) +
  facet_wrap(~Data, scales="free_y") +
  theme_classic()



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


ggplot(masterlist_data, aes(x = as.numeric(Hatch_spread), y = as.numeric(survival_cons), 
                            colour = Females)) +
  geom_point() +
  geom_jitter(position = position_jitter(width = 0.5, height = 0), size = 2, 
              alpha = 0.7) +
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

# Hatch rate and hatch spread plot
# Figure out how to add linear regressions bound by 0 and 1
# Ian suggested pivoting the data longer and doing raw counts of hatched eggs 
# and total clutch size as that will be more informative in a glmm than just 
# hatch success since there is a lot of variation in clutch size. 

# Clean and prepare data
df <- masterlist_data
df$Hatch_success <- as.numeric(df$Hatch_success)
df$Hatch_success[df$Hatch_success == 1] <- 0.999
df$Hatch_success[df$Hatch_success == 0] <- 0.001
df$Hatch_spread <- as.numeric(df$Hatch_spread)

# Fit beta regression
fit <- betareg(Hatch_success ~ Hatch_spread, data = df)

# Generate prediction data
pred_data <- data.frame(Hatch_spread = seq(min(df$Hatch_spread, na.rm = TRUE),
                                           max(df$Hatch_spread, na.rm = TRUE), 
                                           length.out = 100))

# Design matrix for new data (for SE calculation)
X <- model.matrix(~ Hatch_spread, data = pred_data)

# Get predicted values on link (logit) scale
pred_link <- predict(fit, newdata = pred_data, type = "link")

# Get variance-covariance matrix of model coefficients
vcov_mat <- vcov(fit)

# Compute standard errors on the link scale
se_link <- sqrt(rowSums((X %*% vcov_mat[1:2, 1:2]) * X))  # Only use mean model part of vcov

# Inverse logit function
inv_logit <- function(x) exp(x) / (1 + exp(x))

# Compute predicted fit and CI on response scale
pred_data$fit   <- inv_logit(pred_link)
pred_data$lower <- inv_logit(pred_link - 1.96 * se_link)
pred_data$upper <- inv_logit(pred_link + 1.96 * se_link)

ggplot(df, aes(x = Hatch_spread, y = Hatch_success, colour = Females)) +
  geom_jitter(width = 0.5, height = 0, size = 2, alpha = 0.7) +
  geom_ribbon(data = pred_data, aes(x = Hatch_spread, ymin = lower, ymax = upper), 
              fill = "grey70", alpha = 0.4, inherit.aes = FALSE) +
  geom_line(data = pred_data, aes(x = Hatch_spread, y = fit), 
            color = "black", size = 0.75, inherit.aes = FALSE) +
  labs(y = "Hatch rate", x = "Hatch spread (days)") +
  theme_classic()


# Plot
ggplot(df, aes(x = Hatch_spread, y = Hatch_success, colour = Females)) +
  geom_jitter(width = 0.5, height = 0, size = 2, alpha = 0.7) +
  geom_ribbon(data = pred_data, aes(ymin = lower, ymax = upper), inherit.aes = FALSE,
              fill = "grey70", alpha = 0.4) +
  geom_line(data = pred_data, aes(y = fit), inherit.aes = FALSE, color = "black", size = 1.2) +
  labs(y = "Hatch rate", x = "Hatch spread (days)") +
  theme_classic()

ggplot(masterlist_data, aes(x = Hatch_spread, y = Hatch_success, colour = Females)) +
  geom_jitter(width = 0.5, height = 0, size = 2, alpha = 0.7) +
  labs(y = "Hatch rate", x = "Hatch spread (days)") +
  theme_classic()

ggplot(masterlist_data, aes(x = Clutch_size, y = Hatch_success)) +
  geom_jitter(width = 0.5, height = 0, size = 2, alpha = 0.7) +
  labs(y = "Hatch rate", x = "Clutch size") +
  theme_classic()


# Get residuals: asynchrony relative to expected span for brood size
# Doesn't follow assumptions of all nests incubated before or after 
# egg laying period is done though.Needs some touch ups still. 
model <- lm(Hatch_spread ~ Hatched_eggs, data = masterlist_data)
masterlist_data$hatching_asynchrony <- resid(model)
masterlist_data <- masterlist_data 
hist(masterlist_data$hatching_asynchrony)

glm_fit <- glm(as.numeric(Hatch_spread) ~ as.numeric(Hatched_eggs), family = poisson, data = masterlist_data)
summary(glm_fit)$dispersion
overdispersion <- sum(residuals(glm_fit, type = "pearson")^2) / df.residual(glm_fit)
print(overdispersion)

new_data <- data.frame(Hatched_eggs = seq(min(masterlist_data$Hatched_eggs, na.rm = TRUE),
                                          max(masterlist_data$Hatched_eggs, na.rm = TRUE), length.out = 100))

# Get predictions on response scale (counts)
pred <- predict(glm_fit, newdata = new_data, type = "link", se.fit = TRUE)
new_data$fit <- exp(pred$fit)
new_data$lower <- exp(pred$fit - 1.96 * pred$se.fit)
new_data$upper <- exp(pred$fit + 1.96 * pred$se.fit)


#### HYPOTHESIS 1 ####
#### PREDICTION 1b: Hatch order less important in synchronous nests ####

chick_data <- read_excel("Final_Compiled_Chick_Data.xlsx") %>%
  mutate(ChickID = paste(Year, Nest_ID, sep = "_")) %>%
  filter(!is.na(Survived) & Survived %in% c(0, 1)) %>%
  mutate(Survived = if (!is.numeric(Survived)) as.numeric(Survived) else Survived)
View(chick_data)







chick_data$Survived <- as.numeric(as.character(chick_data$Survived))

# Model with hatch order as sole predictor
# Same as Cody's paper but with Nest ID nested within year. 
model_d <- glmmTMB(Survived ~ Hatch_order +
                     (1|Year), family = binomial, data = chick_data)
check_model(model_d)
summary(model_d)

# Model with interaction between hatch order and hatch spread
model_d <- glmmTMB(Survived ~ Hatch_Day +
                     (1|Year), family = binomial, data = chick_data)
check_model(model_d)
summary(model_d)

# Model with interaction between hatch order clutch size
model_d <- glmmTMB(Survived_2 ~ Hatch_order*Clutch_size +
                     (1|Year), family = binomial, data = )
check_model(model_d)
summary(model_d)

# Model of laying order on hatching order 
model_d <- glmmTMB(Survived_2 ~ Hatch_order*Clutch_size +
                     (1|Year), family = binomial, data = )
check_model(model_d)
summary(model_d)

# Model of laying spread on hatching spread
model_d <- glmmTMB(Survived_2 ~ Hatch_order*Clutch_size +
                     (1|Year), family = binomial, data = )
check_model(model_d)
summary(model_d)


# Binomial figure of survival by hatch order
ggplot(chick_data, aes(x = Hatch_order, y = Survived)) +
  geom_jitter(width = 0.5, height = 0, size = 2, alpha = 0.7) +
  labs(y = "Survived", x = "Hatch order") +
  theme_classic()

chick_data_longer <- chick_data %>%
  filter(Mass < 33) %>%
  filter(`Shield to Tip` < 24) %>%
  pivot_longer(cols = c('Mass', 'Tarsus',
                        'Shield to Tip'), 
               names_to = 'Morphometrics', 
               values_to = 'Value') %>%
  filter(!is.na(Value))

View(chick_data_longer)

ggplot(chick_data_longer, aes(x = Hatch_Day, y = Value, colour = Year)) +
  geom_jitter(width = 0.25, height = 0.25, size = 2, alpha = 0.7) +
  facet_wrap(~Morphometrics, scales="free_y", labeller = labeller(
    Morphometrics = c(Mass = "Mass (g)",
             Tarsus = "Left outer tarsus (mm)", 
             `Shield to Tip` = "Shield to tip (mm)"))) +
  scale_y_continuous() +
  theme_classic()

chick_data <- chick_data %>% 
  filter(Mass < 33) %>%
  filter(`Shield to Tip` < 24) 
ggplot(chick_data, aes(x = as.numeric(`Shield to Tip`), y = Survived, colour = Year)) +
  geom_jitter(width = 0.25, height = 0.25, size = 2, alpha = 0.7) + 
  theme_classic()









# Binomial figure of survival by hatch day 
ggplot(chick_data, aes(x = Hatch_Day, y = Hatch_order, colour = Survived)) +
  geom_jitter(width = 0.5, height = 0.5, size = 2, alpha = 0.7) +
  labs(y = "Hatch Order", x = "Hatch Day") +
  theme_classic()



#### HYPOTHESIS 2 ####
# Females lay together because the optimal clutch size is greater 
# than the size of single female clutches.

# PREDICTION 2a: Larger clutches hatch more eggs

# Model analyzing number of hatched eggs by clutch size
hatched_eggs <- glmmTMB(Hatched_eggs ~ I(Clutch_size^2) + 
                          (1|Year), family = poisson, data = wider_data)
check_model(hatched_eggs)

# Model analyzing hatching rate by clutch size
hatch_success <- glmmTMB(Hatch_success ~ Clutch_size + 
                          (1|Year), family = gaussian, data = masterlist_data)
check_model(hatch_success)


# PREDICTION 2b: Larger broods have more survivors

# Conservative model analyzing number of survivors by clutch size and # of hatched eggs
survived_cons_clutch <- glmmTMB(Survived_2_cons ~ I(Clutch_size^2) + Hatched_eggs +
                          (1|Year), family = poisson, data = masterlist_data)
check_model(survived_cons)

# Liberal model analyzing number of survivors by clutch size and # of hatched eggs
survived_clutch <- glmmTMB(Survived_2 ~ I(Clutch_size^2) + Hatched_eggs
                          (1|Year), family = poisson, data = lib_data)
check_model(hatched_eggs)

# Pivoting data longer to make faceted figure
longer_data <- masterlist_data %>%
pivot_longer(cols = c('Hatched_eggs', 'Hatch_success',
                      'Survived_2', 'Survived_2_cons'), 
             names_to = 'Data', 
             values_to = 'Value')

# Faceted figure of hatching success, hatched eggs, and conservative and liberal
# estimates of survival by hatch spread. Could also use geom_jitter instead of geom_count 
# geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 2, alpha = 0.7).
ggplot(longer_data, aes(x = Clutch_size, y = Value)) +
  geom_count() +  
  facet_wrap(~Data, scales="free_y", labeller = labeller(
    Data = c(Hatched_eggs = "Number of hatched eggs",
             Hatch_success = "Hatching rate", 
             Survived_2 = "Number of survivors",
             Survived_2_cons = "Conservative number of survivors"))) +
  scale_size_area(max_size = 6) + 
  labs(y= "Counts", x = "Clutch size") +
  theme_classic()

# Figure of number survivors by number of hatched eggs with liberal estimate
ggplot(masterlist_data, aes(x = Hatched_eggs, y = Survived_2)) +
  geom_count() +
  scale_size_area(max_size = 6) + 
  theme_classic() +
  labs(y = "Survived", x = "Hatched eggs") +
  scale_y_continuous(c(0, 5)) +
  scale_x_continuous(breaks = seq(0, 8, by = 1))

# Figure of number of survivors by number hatched eggs with conservative estimate
ggplot(masterlist_data, aes(x = Hatched_eggs, y = Survived_2_cons)) +
  geom_count() +
  scale_size_area(max_size = 6) + 
  theme_classic() +
  labs(y = "Survived", x = "Hatched eggs") +
  scale_y_continuous(c(0, 5)) +
  scale_x_continuous(breaks = seq(0, 10, by = 1))









########## CHAPTER 3 #########
#### SYNCHRONY EXPERIMENT ####
##############################

#### CH.3 INTRO AND HYPOTHESIS ####
# Here we experimentally lengthened and shortened the HS of nests 
# in Pūkeko by using known laying dates, candling images, and flotation scores.

# HYPOTHESIS: Jointly laying synchronously with a close relative increases
# the inclusive fitness of the dominant female compared to hypothetical 
# prolonged laying of a similar sized clutch by a single female.

# PREDICTION A: Synch nests have greater hatching success than Asynch nests 
# as there are fewer eggs deserted in the nest after hatching. 
# PREDICTION B: Synch nests have greater survival than Asynch nests as there 
# are fewer late hatching chicks with poorer survival chances. 

#### DATA FILTERING AND ORGANIZATION #### 
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
         Observed_nesting_period = as.numeric(Hatch_begin - Date_found)) %>%
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
  guides(color = "none") +
  theme_classic() +
  scale_color_manual(values = c("Asynchronous" = "darkblue", "Synchronous" = "firebrick2")) +
  labs(x = NULL, y = NULL) 

# Example template for doing a single panel instead of facet wrapping them together. 
ggplot(successful_nests_HS, aes(x = Treatment, y = as.numeric(True_hatch_spread), colour = Treatment)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(position = position_jitter(width = 0.1, height = 0), size = 3, alpha = 0.7) +
  labs(y = "Estimated hatch spread (days)", x = "Treatment") +
  guides(color = FALSE) +
  theme_classic() +
  scale_color_manual(values = c("Asynchronous" = "darkblue", "Synchronous" = "firebrick2"))

# Survival to 60 days in pivoted longer format for binomial models and figure.
data_long_format <- successful_nests %>%
  mutate(NonSurvivors = Hatched - Survival_60) %>% # Calculate non-survivors
  pivot_longer(cols = c(Survival_60, NonSurvivors), 
               names_to = "Status", 
               values_to = "Count") %>%
  mutate(Survival_60 = ifelse(Status == "Survival_60", 1, 0)) %>%
  uncount(Count) %>%
  mutate(Treatment_survival = paste(Treatment, Survival_60, sep = "_")) 

# Bar plot of survival data.
ggplot(data_long_format, aes(x = Treatment, fill = factor(Treatment_survival))) +
  geom_bar(position = "dodge") +
  labs(x = "Treatment", y = "Count", fill = "Survived") +
  scale_fill_manual(values = c("Asynchronous_1" = "darkblue", "Asynchronous_0" = "blue", 
                               "Synchronous_1" = "firebrick4", "Synchronous_0" = "firebrick1"), 
                    guide = "none") +
  theme_classic() 


### T-TESTS ###

# Analysis of difference in proportion of swapped or original eggs between treatments. 
t.test(Foreign_percentage ~ Treatment, data = successful_nests)

# Analysis of difference in final manipulated clutch size between treatments. 
t.test(Manipulated_clutch_size ~ Treatment, data = successful_nests)

# Analysis of difference in observed nesting period between treatments. 
t.test(Observed_nesting_period ~ Treatment, data = successful_nests)

# Hatch spread analysis adjusted to only include nests that hatched more than one egg. 
successful_nests_HS <- successful_nests %>%
  filter(True_hatch_spread >1)
t.test(True_hatch_spread ~ Treatment, data = successful_nests)

# Transferred time analysis adjusted to exclude nest that had egg manipulation
# occur after hatching had begun.
t.test(Transfered_time ~ Treatment, data = successful_nests)

# Analysis of difference in hatching success between treatments
t.test(Hatch_success ~ Treatment, data = successful_nests)

### POWER ANALYSES ###

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


# PREDICTION B: Synch nests have greater survival than Asynch nests as there 
# are fewer late hatching chicks with poorer survival chances.

### GLMMS ###

# Model analyzing the hatch spread by clutch size
model_1 <- glmmTMB(True_hatch_spread ~ Manipulated_clutch_size + (1|Year),
                       family = poisson, data = successful_nests)
check_model(model_1)

# Model analyzing survival by treatment. 
model_2 <- glmmTMB(Survival_60 ~ Treatment  + (1|Year), family = binomial,
                   data = data_long_format)
check_model(model_2)
summary(model_2)

# Model analyzing survival by treatment and clutch size
model_3 <- glmmTMB(Survival_60 ~ Treatment + Manipulated_clutch_size +
                      (1|Year), family = binomial, data = data_long_format)
check_model(model_3)
summary(model_3)
