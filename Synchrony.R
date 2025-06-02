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

# Primary filtering/data modifications of the masterlist data set. 
masterlist_data <- read_excel("Nests_masterlist.xlsx") %>%
  filter(Exclusion == "GOOD" | Exclusion == "ABCL") %>% # Keeping only the good nests
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>% # Unique nest ID 
  mutate(Clutch_size = as.numeric(Clutch_size)) %>%
  filter(Hatch_begin != "NA" & Hatch_end != "NA") %>%
  mutate(Females = ifelse(Clutch_size <= 5, "Single Female", "Joint Females")) %>% # assign single/joint female
  mutate(Hatch_begin = as.numeric(Hatch_begin),
         Hatch_end = as.numeric(Hatch_end),
         Hatch_spread = Hatch_end - Hatch_begin + 1) %>% # First to last day inclusive
  mutate(Hatch_begin = as.numeric(Hatch_begin), Date_found = as.numeric(Date_found), 
       Observed_incubation_period = as.numeric(Hatch_begin - Date_found)) %>%
  mutate(Hatch_success = as.numeric(Hatched_eggs)/as.numeric(Clutch_size))
View(masterlist_data)

# Check differences in joint vs singe female nests
masterlist_data %>% count(Females, sort=TRUE)

# Filtering of data for only failed nests 
masterlist_data <- read_excel("Nests_masterlist.xlsx") %>%
  filter(Exclusion == "GOOD" | Exclusion == "ABCL") %>% # Keeping only the good nests
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>% # Unique nest ID 
  filter(Hatch_begin == "NA" & Hatch_end == "NA") %>%
  mutate(Clutch_size = as.numeric(Clutch_size)) %>%
  mutate(Females = ifelse(Clutch_size <= 5, "Single Female", "Joint Females"))  # assign single/joint female
View(masterlist_data)

# Check difference in joint vs single female nests
masterlist_data %>% count(Females, sort=TRUE)


# Get residuals: asynchrony relative to expected span for brood size
# Doesn't follow assumptions of all nests incubated before or after 
# egg laying period is done though.Needs some touch ups still. 
model <- lm(Hatch_spread ~ Hatched_eggs, data = masterlist_data)
masterlist_data$hatching_asynchrony <- resid(model)
masterlist_data <- masterlist_data 
hist(masterlist_data$hatching_asynchrony)


# Do pivot longer and combine things like hatch success, observed period.
# hatched eggs, survival etc. 
masterlist_data <- masterlist_data %>%
  filter(Hatched_eggs > 1)
ggplot(masterlist_data, aes(x = as.numeric(Hatched_eggs), y = as.numeric(Hatch_spread), colour = Females)) +
  geom_point() +
  geom_jitter(position = position_jitter(width = 0.5, height = 0), size = 2, alpha = 0.7) +
  theme_classic()


masterlist_data <- read_excel("Nests_masterlist.xlsx") %>%
  filter(Exclusion == "GOOD" | Exclusion == "ABCL") %>% # Keeping only the good nests
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>% # Unique nest ID 
  mutate(Clutch_size = as.numeric(Clutch_size)) %>%
  mutate(Females = ifelse(Clutch_size <= 5, "Single Female", "Joint Females")) %>% # assign single/joint female
  mutate(Hatch_begin = as.numeric(Hatch_begin),
         Hatch_end = as.numeric(Hatch_end),
         Hatch_spread = Hatch_end - Hatch_begin + 1) %>% # First to last day inclusive
  mutate(Hatch_begin = as.numeric(Hatch_begin), Date_found = as.numeric(Date_found), 
         Observed_incubation_period = as.numeric(Hatch_begin - Date_found)) %>%
  mutate(Hatch_success = as.numeric(Hatched_eggs)/as.numeric(Clutch_size))

mean()







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


ggplot(masterlist_data, aes(x = as.numeric(Hatched_eggs), y = as.numeric(Hatch_spread), colour = Females)) +
  geom_jitter(width = 0.5, height = 0, size = 2, alpha = 0.7) +
  labs(y = "Hatch spread (days)", x = "Number of hatched eggs") +
  theme_classic()

ggplot(masterlist_data, aes(x = as.numeric(Clutch_size), y = as.numeric(Hatch_spread), colour = Females)) +
  geom_jitter(width = 0.5, height = 0, size = 2, alpha = 0.7) +
  labs(y = "Hatch spread (days)", x = "Clutch size") +
  theme_classic()









# Hatch rate and hatch spread plot
# Figure out how to add linear regressions bound by 0 and 1
# Ian suggested pivoting the data longer and doing raw counts of hatched eggs 
# and total clutch size as that will be more informative in a glmm than just 
# hatch success since there is a lot of variation in clutch size. 

library(betareg)
library(ggplot2)

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








# Make a copy and adjust 0s or 1s if needed
df <- masterlist_data
df$Hatch_success <- as.numeric(df$Hatch_success)
df$Hatch_success[df$Hatch_success == 1] <- 0.999
df$Hatch_success[df$Hatch_success == 0] <- 0.001
df$Hatch_spread <- as.numeric(df$Hatch_spread)

library(betareg)
fit <- betareg(Hatch_success ~ Hatch_spread, data = df)

# Create new data for predictions
pred_data <- data.frame(Hatch_spread = seq(min(df$Hatch_spread, na.rm = TRUE),
                                           max(df$Hatch_spread, na.rm = TRUE),
                                           length.out = 100))

# Predict fitted values (on response scale)
pred_data$Hatch_success <- predict(fit, newdata = pred_data, type = "response")



ggplot(masterlist_data, aes(x = as.numeric(Hatch_spread), 
                            y = as.numeric(Hatch_success), colour = Females)) +
  geom_point() +
  geom_jitter(position = position_jitter(width = 0.5, height = 0), size = 2, alpha = 0.7) +
  geom_line(data = pred_data, aes(x = Hatch_spread, y = Hatch_success), inherit.aes = FALSE,
            color = "black", size = 1.2) +
  labs(y = "Hatch rate", x = "Hatch spread (days)") +
  theme_classic()









ggplot(masterlist_data, aes(x= as.numeric(Clutch_size),
                            y= as.numeric(Hatched_eggs), 
                            colour = Females)) +
  geom_point() +
  geom_jitter(position = position_jitter(width = 0.5, height = 0), size = 2, alpha = 0.7) +
  labs(y = "Hatched eggs", x = "Clutch size") +
  theme_classic()

ggplot(masterlist_data, aes(x= as.numeric(Clutch_size),
                            y= as.numeric(Hatch_success),
                            colour = Females)) +
  geom_point() +
  geom_jitter(position = position_jitter(width = 0.5, height = 0), size = 2, alpha = 0.7) +
  labs(y = "Hatched success", x = "Clutch size") +
  theme_classic()

# Survived by hatched eggs.
df <- masterlist_data %>%
  filter(!is.na(as.numeric(Survived_2))) %>%
  group_by(Hatched_eggs, Survived_2) %>%
  summarise(count = n(), .groups = "drop")

ggplot(df, aes(x = Hatched_eggs, y = Survived_2, size = factor(count))) +
     geom_point(alpha = 0.6) +
     scale_size_manual(values = c("1" = 1, "2" = 2, "3" = 3, "4" = 4, "5" = 6, "7" = 7, "11" = 11)) + 
     theme_classic() +
     labs(y = "Survived", x = "Hatched eggs") 

# Surived by hatched eggs again but jittered points and by female.
ggplot(masterlist_data, aes(x = as.numeric(Hatched_eggs), y = as.numeric(Survived_2), colour = Females)) +
  geom_point(alpha = 0.6) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0.2), size = 2, alpha = 0.7) +
  theme_classic() +
  labs(y = "Survived", x = "Hatched eggs") +
  scale_y_continuous(c(0, 5)) +
  scale_x_continuous(breaks = seq(0, 8, by = 1))


# Hatched eggs by hatch spread with jitterpoints and by female
masterlist_data <- masterlist_data %>%
  filter(Hatched_eggs != "NA") %>%
  filter(Hatch_spread != "NA") 

ggplot(masterlist_data, aes(x = as.numeric(Hatch_spread), y = as.numeric(Hatched_eggs), colour = Females)) +
  geom_point(alpha = 0.6) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0.2), size = 2, alpha = 0.7) +
  theme_classic() +
  labs(y = "Hatched eggs", x = "Hatch spread") 

# Hatch survival by hatch spread with jitter points. 
ggplot(masterlist_data, aes(x = as.numeric(Hatch_spread), y = as.numeric(Survived_2), colour = Females)) +
  geom_point(alpha = 0.6) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0.2), size = 2, alpha = 0.7) +
  theme_classic() +
  labs(y = "Survived", x = "Hatch spread") 

ggplot(masterlist_data, aes(x = Hatch_spread, y = as.numeric(Hatched_eggs))) +
  geom_point(alpha = 0.6) +
  theme_classic() +
  labs(y = "Hatched eggs", x = "Hatch spread") 

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


ggplot(masterlist_data, aes(x = as.numeric(Hatched_eggs), y = as.numeric(Survived_2), color = Females)) +
geom_point(alpha = 0.6) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0.2), size = 2, alpha = 0.7) +
  theme_classic() +
  labs(y = "Survived", x = "Hatched eggs")



# Survival data with conservative and liberal estimates
# Conservative estimation that all NAs for survival to 2 months are actually 0 
# whereas a more liberal estimation is removing the NAs from the analysis. 
masterlist_data <- masterlist_data %>%
  mutate(survival_cons = ifelse(is.na(Survived_2), 0, Survived_2))

ggplot(masterlist_data, aes(x = as.numeric(Hatch_spread), y = as.numeric(survival_cons), 
      colour = Females)) +
      geom_point() +
      geom_jitter(position = position_jitter(width = 0.5, height = 0), size = 2, 
      alpha = 0.7) +
      theme_classic()

ggplot(masterlist_data, aes(x = as.numeric(Hatch_spread), y = as.numeric(Survived_2), 
      colour = Females)) +
      geom_point() +
      geom_jitter(position = position_jitter(width = 0.5, height = 0), size = 2, 
      alpha = 0.7) +
      theme_classic()


# Create unique nest IDs
# Merge clutch size and hatched eggs from nest data set to chick data set

mutate(Synchrony = case_when(
  Hatch_spread >= 1 & Hatch_spread <= 5 ~ "Synchronous",
  Hatch_spread >= 6 ~ "Asynchronous",
  TRUE ~ NA_character_ )) # Handles missing or unexpected values


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

#### GLMM STUFF ####


# Primary filtering/data modifications of the masterlist data set. 
masterlist_data <- read_excel("Nests_masterlist.xlsx") %>%
  filter(Exclusion == "GOOD" | Exclusion == "ABCL") %>% # Keeping only the good nests
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>% # Unique nest ID 
  mutate(Clutch_size = as.numeric(Clutch_size)) %>%
  filter(Clutch_size != "" | Clutch_size != " " | Clutch_size == "NA" | is.na(Clutch_size))
  filter(Hatched_eggs != "4?") %>%
  filter(Hatch_begin != "NA" & Hatch_end != "NA") %>%
  mutate(Females = ifelse(Clutch_size <= 5, "Single Female", "Joint Females")) %>% # assign single/joint female
  mutate(Hatch_begin = as.numeric(Hatch_begin),
         Hatch_end = as.numeric(Hatch_end),
         Hatch_spread = Hatch_end - Hatch_begin + 1) %>% # First to last day inclusive
  mutate(Hatch_begin = as.numeric(Hatch_begin), Date_found = as.numeric(Date_found), 
         Observed_incubation_period = as.numeric(Hatch_begin - Date_found)) %>%
  mutate(Hatch_success = as.numeric(Hatched_eggs)/as.numeric(Clutch_size))
View(masterlist_data)


masterlist_data <- read_excel("Nests_masterlist.xlsx") %>%
  filter(Exclusion == "GOOD" | Exclusion == "ABCL") %>% # Keeping only the good nests
  drop_na(Clutch_size) %>%
  filter(Hatched_eggs != "4?") %>%
  filter(Hatch_begin != "NA" & Hatch_end != "NA") %>%
  mutate(Females = ifelse(Clutch_size <= 5, "Single Female", "Joint Females")) %>% # assign single/joint female
  mutate(Hatch_begin = as.numeric(Hatch_begin),
         Hatch_end = as.numeric(Hatch_end),
         Hatch_spread = Hatch_end - Hatch_begin + 1) %>% # First to last day inclusive
  mutate(Hatch_begin = as.numeric(Hatch_begin), Date_found = as.numeric(Date_found), 
         Observed_incubation_period = as.numeric(Hatch_begin - Date_found)) %>%
  mutate(Hatch_success = as.numeric(Hatched_eggs)/as.numeric(Clutch_size))
View(masterlist_data)

masterlist_data$Clutch_size <- as.numeric(as.character(masterlist_data$Clutch_size))
masterlist_data$Hatched_eggs <- as.numeric(as.character(masterlist_data$Hatched_eggs))




# Model with all fixed effects as random effects to break and have singular fit. 
model_a <- glmmTMB(as.numeric(Hatched_eggs) ~ as.numeric(Clutch_size) + as.numeric(Group_size) + as.numeric(Hatch_spread) +
                  (1|Year) + (1|as.numeric(Clutch_size)) + (1|as.numeric(Group_size)) +
                   (1|as.numeric(Hatch_spread)), family = poisson, data = masterlist_data)
check_model(model_a)















masterlist_data$Clutch_size <- as.numeric(masterlist_data$Clutch_size)


hatched_eggs <- glmmTMB(Hatched_eggs ~ I(Clutch_size^2) + 
                (1|Year), family = poisson, data = masterlist_data, 
                control = glmmTMBControl(rank_check = "adjust"))
check_model(hatched_eggs)






hatched_eggs <- glmmTMB(as.numeric(Hatched_eggs) ~ I(as.numeric(Clutch_size)^2) + as.numeric(Hatch_spread) + 
                 + (1|Year), family = poisson, data = masterlist_data)
check_model(hatched_eggs)



model_c <- glmmTMB(Survived_2 ~ Clutch_size + Manipulated_clutch +
                  (1|Year) (1|Clutch_size) + (1|Hatched_eggs) + (1|Females)
                  (1|Hatch_spread), family = binomial, data = )
check_model(model_c)

model_d <- glmmTMB(Survived_2 ~ Clutch_size + Hatch_order + Hatch_day + Hatch_spread
                  (1|Year) (1|Clutch_size) + (1|Hatched_eggs) + (1|Females)
                  (1|Hatch_spread), family = binomial, data = )
check_model(model_d)

model_e <- glmmTMB(Survived_2 ~ Clutch_size + Manipulated_clutch +
                     (1|Year) (1|Clutch_size) + (1|Hatched_eggs) + (1|Females)
                   (1|Hatch_spread), family = binomial, data = )
check_model(model_e)















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
  guides(color = "none") +
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


