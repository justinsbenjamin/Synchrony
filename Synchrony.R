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
library(simr)

############################
####      CHAPTER 2     ####
#### OBSERVATIONAL DATA ####
############################

# Sample sizes of excluded and included groups
nest_counts <- read_excel("Nests_masterlist.xlsx") %>%
  count(Exclusion, name = "nest_counts") %>%
  arrange(desc("nest_counts"))
print(nest_counts)

masterlist_data <- read_excel("Nests_masterlist.xlsx",
  na = c("", "NO_RECORD", "MISSING")) %>%
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>%  
  filter(Exclusion == "GOOD" | Exclusion == "MISSING_HATCH_ORDER") %>%
  mutate(Hatch_success = (as.numeric(Hatched_eggs)/as.numeric(Clutch_size)*100), 
         Survived_2_cons = ifelse(is.na(Survived_2), 0, Survived_2), 
         Hatch_begin = as.Date(Hatch_begin),
         Hatch_end = as.Date(Hatch_end),
         Hatch_spread = as.numeric(Hatch_end - Hatch_begin + 1),
         Date_found = as.Date(Date_found), 
         Last_observed = as.Date(Last_observed), 
         Observed_nesting_period = as.numeric(Hatch_begin - Date_found), 
         Females = if_else(Clutch_size <= 5,"Single female","Joint female")) %>%
  mutate(across(c(Clutch_size, Hatched_eggs, Survived_2, Survived_2_cons,
         Hatch_begin, Hatch_end, Hatch_spread, Year), 
         ~ as.numeric(as.character(.)))) %>%
  mutate(Survived_2_cons = ifelse(is.na(Survived_2), 0, Survived_2),
        converted_na = is.na(Survived_2))

# Pivoting data longer with two survival estimates
# 1) Conservative estimate: When survival unknown, NAs converted to 0s.
# 2) Liberal estimate: NAs were removed. 
longer_data1 <- masterlist_data %>%
  pivot_longer(cols = c('Hatched_eggs', 'Hatch_success',
                        'Survived_2', 'Survived_2_cons'), 
               names_to = 'Data', 
               values_to = 'Value') %>%
  mutate(converted_na = ifelse(Data == "Survived_2_cons", converted_na, FALSE)) %>%
  pivot_longer(cols = c("Hatch_spread", "Clutch_size"), 
               names_to = "Predictor", 
               values_to = "Predictor_value")

# Figure of hatch rate, number of hatched eggs and number of survivors with 
# conservative and liberal estimates by HS and clutch size. 
figure_1 <- ggplot(longer_data1, aes(x = Predictor_value, y = Value)) +
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
  theme(strip.background = element_blank(), strip.placement = "outside", 
        strip.text = element_text(size = 10))
print(figure_1)

# Added colours to the NAs that were converted to 0s
figure_2 <- figure1 + aes(x = Predictor_value, y = Value, color = converted_na) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  labs(x = NULL, y = NULL) 
print(figure_2)

# Added a smoother function to see the distribution of the data. 
figure_3 <- figure1 +
  geom_smooth(method = "loess") 
print(figure_3)

ggplot(masterlist_data, aes(x = Hatched_eggs, y = Survived_2)) +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.8) +
  theme_classic()

# PREDICTION 1A: Shorter HS have better hatching success
# PREDICTION 2A: Number of hatched eggs increase with clutch size 
# then decrease when clutch is too large. 
  
# Model analyzing proportion of hatched egg by clutch size and HS
model_1 <- glmmTMB(cbind(Hatched_eggs, Clutch_size - Hatched_eggs) ~ 
                  Hatch_spread + Clutch_size + 
                  (1|Year), family = betabinomial,
                  data = masterlist_data)

# Function to view model diagnostics 
diagnostics <- function(model) {
print(check_model(model)) 
res <- simulateResiduals(fittedModel = model, plot = TRUE)
testDispersion(res)
testZeroInflation(res)}

diagnostics(model_1)
summary(model_1)

# PREDICTION 1B: Shorter HS have better survival
# PREDICTION 2B: Larger broods have more survivors. 

# Model analyzing conservative estimate of survivors by clutch size and HS
model_2 <- glmmTMB(cbind(Survived_2_cons, Hatched_eggs - Survived_2_cons) ~ 
                   Hatch_spread + Hatched_eggs + Clutch_size + (1|Year), 
                   family = betabinomial, data = masterlist_data)
diagnostics(model_2)
summary(model_2)

# Model analyzing liberal survivors (liberal estimate) by clutch size and HS
lib_data <- masterlist_data %>% filter(!is.na(Survived_2))

model_3 <- glmmTMB(cbind(Survived_2, Hatched_eggs - Survived_2) ~ 
                   Hatch_spread + Hatched_eggs + Clutch_size + (1|Year), 
                   family = betabinomial,
                   data = lib_data)
diagnostics(model_3)
summary(model_3)

# PREDICTION 2C: Larger clutches have larger hatch spreads
correlation_a <- cor(as.numeric(masterlist_data$Hatch_spread), 
                   as.numeric(masterlist_data$Clutch_size), method = 'pearson')
correlation_b <- cor(as.numeric(masterlist_data$Hatch_spread), 
                   as.numeric(masterlist_data$Hatched_eggs), method = 'pearson')
correlation_a
correlation_b

ggplot(masterlist_data, aes(x = Clutch_size, y = Hatch_spread)) +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.8) +
  theme_classic()

ggplot(masterlist_data, aes(x = Hatched_eggs, y = Hatch_spread)) +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.8) +
  theme_classic()

# HYPOTHESIS 3: Subordinate females benefit more from synchronous nests

# Prediction 3a: Positive correlation between lay order and hatch order
laying_data <- read_excel("Lay_hatch_data.xlsx") %>%
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) 

laying_data <- laying_data %>%
  filter(Known_lay_order < 10)
ggplot(laying_data, aes(x = Known_lay_order, y = Hatch_order)) +
  geom_jitter(width = 0.05, height = 0.05, size = 2, alpha = 0.7) +
  labs(y = "Hatch order", x = "Lay order") +
  theme_classic()  + 
  scale_x_continuous(breaks = 0:10) + 
  scale_y_continuous(breaks = 0:10) +
  geom_smooth(method = "lm", se = TRUE)
  
correlation <- cor.test(as.numeric(laying_data$Known_lay_order), 
               as.numeric(laying_data$Hatch_order), method = 'pearson')
print(c(correlation$estimate, correlation$conf.int, correlation$p.value))

ggplot(laying_clean, aes(x = Known_lay_order, y = Incubation_period)) +
  geom_jitter(width = 0.05, height = 0.05, size = 2, alpha = 1) +
  labs(y = "Observed nesting period (days)", x = "Laying order") +
  theme_classic() 

laying_clean <- laying_data %>% 
  mutate(Date = as.Date(Date)) %>%
  add_count(Nest_ID, Date, Known_lay_order, name = "n_dup") %>% 
  filter(n_dup == 1) %>% 
  select(-n_dup)

break_pt <- 5        
laying_clean <- within(laying_clean, {
  lay_1to5  <- pmin(Known_lay_order, break_pt)
  lay_over5 <- pmax(0, Known_lay_order - break_pt)})

model_abc <- glmmTMB(Incubation_period ~ lay_1to5 + lay_over5, 
                     family = compois(link = "log"), data = laying_clean)
diagnostics(model_abc)
summary(model_abc)

#### PREDICTION 3b: Earlier hatched eggs have greater survival

chick_data <- read_excel("Final_Compiled_Chick_Data.xlsx") %>%
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>%
  mutate("Shield to Tip" = "S")

cols <- c("Clutch_size", "Hatch_spread", "Hatched_eggs")
idx  <- match(chick_data$Nest_ID, masterlist_data$Nest_ID)
chick_data[cols] <- masterlist_data[idx, cols] 

chick_data_survival <- chick_data %>%
  filter(!is.na(Survived) & Survived %in% c(0, 1)) %>%
  mutate(Survived = if (!is.numeric(Survived)) as.numeric(Survived) else Survived)

chick_data_survival <- chick_data_survival %>%
  filter(Hatched_eggs > 1)

correlation <- cor(as.numeric(chick_data_survival$Hatch_order), 
                   as.numeric(chick_data_survival$Clutch_size), method = 'pearson')
correlation

# Figure of survival by hatching order. 
ggplot(chick_data_survival, aes(x = Hatch_order, y = as.character(Survived), colour = Year)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.7) +
  labs(y = "Survived", x = "Hatch_order") +
  theme_classic()

# Model with hatch order as sole predictor variable on survival.
model_2_3a <- glmmTMB(as.numeric(Survived) ~ as.numeric(Hatch_order) +
                  (1|Nest_ID) + (1|Year), family = binomial, data = chick_data_survival)
diagnostics(model_2_3a)
summary(model_2_3a)

# Model with interaction between hatch order and hatching spread
model_2_3c <- glmmTMB(as.numeric(Survived) ~ as.numeric(Hatch_order)*Hatch_spread +
                     (1|Nest_ID) + (1|Year), family = binomial, data = chick_data_survival)
diagnostics(model_2_3c)
summary(model_2_3c)

# PREDICTION 3c: Last hatched eggs are smaller as they spend more developmental
# energy catching up to early hatching chicks. 

model_2_3d <- glmmTMB(as.numeric(Survived) ~ as.numeric(Hatch_order)*Mass + 
                        (1|Nest_ID) + (1|Year), family = binomial, data = chick_data_survival)
diagnostics(model_2_3d)
summary(model_2_3d)

# Pivoting data longer then filtering out data that has unusually large birds
# (likely typos or issues with the measurements and potential recaptures). 
chick_data_morph <- chick_data %>%
  filter(!Mass > 35) &
      !(`Shield to Tip` > 24)  &
      !(Tarsus > 40) %>%
  filter("Tarsus" %in% Year < 2019)

model_mass <- glmmTMB(Mass ~ as.numeric(Hatch_order) + (1|Nest_ID) + (1|Year), 
                      family = poisson, data = chick_data_survival)
diagnostics(model_mass)
summary(model_mass)

model_StoT <- glmmTMB(StoT ~ as.numeric(Hatch_order) + (1|Nest_ID) + (1|Year), 
                      family = gausssian, data = chick_data_survival)
diagnostics(model_StoT)
summary(model_StoT)

model_Tars <- glmmTMB(Tarsus ~ as.numeric(Hatch_order) + (1|Nest_ID) + (1|Year), 
                      family = gaussian, data = chick_data_survival)
diagnostics(model_Tars)
summary(model_Tars)

# Pivoting data longer then filtering out data that has unusually large birds
# (likely typos or issues with the measurements and potential recaptures). 
chick_data_longer <- chick_data %>%
  pivot_longer(cols = c('Mass', 'Tarsus',
                        'Shield to Tip'), 
               names_to = 'Morphometrics', 
               values_to = 'Value') %>%
  filter(!is.na(Value)) %>%
  filter(!(Morphometrics == "Mass" & (Value > 35)) &
           !(Morphometrics == "Shield to Tip" & (Value > 24))  &
           !(Morphometrics == "Tarsus" & Year < 2019)) %>%
  filter(!(Morphometrics == "Tarsus" & Value > 40))

# Figure of chick size by hatch order
ggplot(chick_data_longer, aes(x = Hatch_order, y = Value)) +
  geom_jitter(width = 0, height = 0, size = 2, alpha = 0.7) +
  facet_wrap(~Morphometrics, scales="free_y", labeller = labeller(
    Morphometrics = c(Mass = "Mass (g)",
                      Tarsus = "Left outer tarsus (mm)", 
                      `Shield to Tip` = "Shield to tip (mm)"))) +
  scale_y_continuous() +
  labs(x = "Hatch order", y = "Measurement") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  theme(strip.background = element_blank(), strip.placement = "outside", 
        strip.text = element_text(size = 10))

# Figure of survival by size at hatching 
chick_data_longer <- chick_data_longer %>%
  filter(!is.na(Survived)) %>%
  filter(! Survived == "NA")
  
ggplot(chick_data_longer, aes(x = Value, y = as.character(Survived))) +
  geom_point(size = 2, alpha = 0.7) +
  facet_wrap(~Morphometrics, scales="free_x", labeller = labeller(
  Morphometrics = c(Mass = "Mass (g)", Tarsus = "Left outer tarsus (mm)", 
                    `Shield to Tip` = "Shield to tip (mm)"))) +
  labs(x = "Measurement", y = "Survived") +
  scale_x_continuous() +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

##############################
####      CHAPTER 3       ####
#### SYNCHRONY EXPERIMENT ####
##############################

# Nests 2018_E, 2018_AD, 2018_AZ, 2018_BD were not manipulated and thus excluded. 
# Nest 2024_W predated during hatching, also excluded. 
experiment_data <- read_excel("Compiled_synchrony_experiment_data.xlsx") %>%
  mutate(Nest_ID = paste(Year, Nest, sep = "_")) %>%
  filter(!Treatment %in% c("Control","Other")) %>%
  filter(!Nest_ID %in% c("2018_E","2018_AD","2018_AZ","2018_BD")) %>%  
  mutate(Treatment = recode(Treatment, "Synch" = "Synchronous",
                                       "Asynch" = "Asynchronous"), 
         Hatch_success = (Hatched/Manipulated_clutch_size), 
         Foreign_percentage = (Foreign_eggs/Manipulated_clutch_size))

# Further cleaning of the data set with successful nests
successful_nests <- experiment_data %>% 
  filter(Hatch_success >0) %>%
  filter(!Nest_ID %in% "2024_W") %>%
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
pivot_data_long <- function(df) { df %>%
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
    summarise(n = n(),
      min = min(Value, na.rm = TRUE),
      max = max(Value, na.rm = TRUE),
      mean = mean(Value, na.rm = TRUE),
      sd = sd(Value, na.rm = TRUE),
      se = sd / sqrt(n),
      .groups = "drop")}

summary_experiment_data <- summarize_data(experiment_data_long)
summary_successful_nests <- summarize_data(successful_nests_long)
summary_successful_HS <- summarize_data(successful_nests_long_HS)
summary_successful_MT <- summarize_data(successful_nests_long_MT)

lapply(list(summary_experiment_data, summary_successful, 
            summary_successful_HS, summary_successful_MT), print, n = 30)

#### T-TESTS analyzing the effectiveness of the experimental design ####
t.test(Foreign_percentage ~ Treatment, data = experiment_data) # proportion swapped eggs
t.test(Manipulated_clutch_size ~ Treatment, data = experiment_data) # manipulated clutch size
t.test(Observed_nesting_period ~ Treatment, data = successful_nests) # observed nesting period
t.test(True_hatch_spread ~ Treatment, data = successful_nests_HS) # Hatch spread
t.test(Transfered_time ~ Treatment, data = successful_nests_MT) # Transferred time

# Model analyzing number of hatched eggs by clutch size
model_1 <- glmmTMB(cbind(Hatched, Manipulated_clutch_size - Hatched) ~ 
                     Treatment + Manipulated_clutch_size,
                     family = betabinomial, data = successful_nests)
diagnostics(model_1)
summary(model_1)

# Model analyzing number of survivors by treatment.
model_2 <- glmmTMB(cbind(Survival_60, Hatched - Survival_60) ~ 
                   Treatment + Hatched + (1|Year), 
                   family = betabinomial, data = successful_nests)
diagnostics(model_2)
summary(model_2)

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
