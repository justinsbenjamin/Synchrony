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
library(nlme)
library(patchwork)  
library(ggbeeswarm)

############################
####      CHAPTER 2     ####
#### OBSERVATIONAL DATA ####
############################

# Number of excluded and included groups
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
         Survived_1_cons = ifelse(is.na(Survived_1), 0, Survived_1), 
         Hatch_begin = as.Date(Hatch_begin),
         Hatch_end = as.Date(Hatch_end),
         Hatch_spread = as.numeric(Hatch_end - Hatch_begin + 1),
         Date_found = as.Date(Date_found), 
         Last_observed = as.Date(Last_observed), 
         Observed_nesting_period = as.numeric(Hatch_begin - Date_found), 
         Females = if_else(Clutch_size <= 5,"Single female","Joint female"), 
         Survived_1_cons = ifelse(is.na(Survived_1), 0, Survived_1),
         converted_na_1 = is.na(Survived_1),
         Survived_2_cons = ifelse(is.na(Survived_2), 0, Survived_2),
         converted_na_2 = is.na(Survived_2)) %>%
  mutate(across(c(Clutch_size, Hatched_eggs, Survived_1, Survived_1_cons,
                  Survived_2, Survived_2_cons, Hatch_begin, Hatch_end, 
                  Hatch_spread, Year), ~ as.numeric(as.character(.)))) %>%
  mutate(Relative_survive_success_2 = (Survived_2/(Hatched_eggs)*100)) %>%
  mutate(Survive_success_2 = (Survived_2/(Clutch_size)*100)) %>%
  mutate(Relative_survive_success_1 = (Survived_1/(Hatched_eggs)*100)) %>%
  mutate(Survive_success_1 = (Survived_1/(Clutch_size)*100)) %>%
  mutate(Converted_na_1_prop = is.na(Relative_survive_success_1), 
         Converted_na_2_prop = is.na(Relative_survive_success_2))

# Making Figure 2.1
masterlist_data <- within(masterlist_data, {
  Clutch_size   <- as.numeric(Clutch_size)
  Hatched_eggs  <- as.numeric(Hatched_eggs)
  Hatch_spread  <- as.numeric(Hatch_spread)
  Survived_1   <- as.numeric(Survived_1)
  Survived_2   <- as.numeric(Survived_2)})

x_vars <- c("Clutch_size", "Hatched_eggs", "Hatch_spread", "Survived_1")
y_vars <- c("Hatched_eggs", "Hatch_spread", "Survived_1", "Survived_2")
var_labels <- c(Clutch_size = "Clutch size",
                Hatched_eggs = "# Hatched eggs",
                Hatch_spread = "Hatch spread (days)",
                Survived_1 = "# Survived to 30 days",
                Survived_2 = "# Survived to 60 days")

invalid_pairs <- list(c("Hatched_eggs", "Hatched_eggs"),
                      c("Hatched_eggs", "Hatch_spread"),
                      c("Hatched_eggs", "Survived_1"),
                      c("Hatch_spread", "Hatch_spread"),
                      c("Hatch_spread", "Survived_1"),
                      c("Survived_1", "Survived_1"))

is_invalid <- function(y, x) {
  any(sapply(invalid_pairs, function(p) all(p == c(y, x))))}

# Find bottom-most *valid* plot for each column
valid_plot_positions <- expand.grid(
  x = x_vars,
  y = rev(y_vars),  # bottom-to-top
  stringsAsFactors = FALSE)

bottom_labels <- list()
for (x in x_vars) {
  for (y in valid_plot_positions$y[valid_plot_positions$x == x]) {
    if (!is_invalid(y, x)) {
      bottom_labels[[paste0(y, "_", x)]] <- TRUE
      break}}}

plot_cell <- function(xvar, yvar, show_x = FALSE, show_y = FALSE) {
  x_data <- masterlist_data[[xvar]]
  y_data <- masterlist_data[[yvar]]
  
  # Convert to numeric if needed
  if (!is.numeric(x_data)) {
    x_data <- suppressWarnings(as.numeric(as.character(x_data)))
  }
  if (!is.numeric(y_data)) {
    y_data <- suppressWarnings(as.numeric(as.character(y_data)))
  }
  
  # Remove NAs before calculating limits and breaks
  x_data_noNA <- x_data[!is.na(x_data)]
  y_data_noNA <- y_data[!is.na(y_data)]
  
  # Set default limits and breaks if no data remains
  if (length(x_data_noNA) == 0) {
    x_min <- 0; x_max <- 1; x_breaks <- 0:1
  } else {
    x_min <- floor(min(x_data_noNA))
    x_max <- ceiling(max(x_data_noNA))
    x_breaks <- seq(x_min, x_max, by = 1)
  }
  
  if (length(y_data_noNA) == 0) {
    y_min <- 0; y_max <- 1; y_breaks <- 0:1
  } else {
    y_min <- floor(min(y_data_noNA))
    y_max <- ceiling(max(y_data_noNA))
    y_breaks <- seq(y_min, y_max, by = 1)
  }
  ggplot(masterlist_data, aes_string(x = xvar, y = yvar)) +
    geom_jitter(width = 0.2, height = 0.1, alpha = 0.8) +
    xlab(if (show_x) var_labels[[xvar]] else NULL) +
    ylab(if (show_y) var_labels[[yvar]] else NULL) +
    scale_x_continuous(breaks = x_breaks, limits = c(x_min, x_max)) +
    scale_y_continuous(breaks = y_breaks, limits = c(y_min, y_max)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x  = if (show_x) element_text(size = 8) else element_blank(),
      axis.text.y  = if (show_y) element_text(size = 8) else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank(),
      axis.ticks.y = if (show_y) element_line() else element_blank(),
      axis.title.x = if (show_x) element_text(size = 10) else element_blank(),
      axis.title.y = if (show_y) element_text(size = 10) else element_blank(),
      plot.margin = margin(2, 2, 2, 2))}

plots <- list()
for (y in rev(y_vars)) {
  for (x in x_vars) {
    if (is_invalid(y, x)) {
      plots[[length(plots) + 1]] <- plot_spacer()
    } else {
      show_x <- (y == bottom_x_labels[[x]])
      show_y <- (x == x_vars[1])
      plots[[length(plots) + 1]] <- plot_cell(x, y, show_x, show_y)}}}

final_plot <- wrap_plots(plots, ncol = 4) 
final_plot


# Analysing difference in survival between 30 days and 60 days
survival_dif <- masterlist_data %>%
  mutate(Has_1 = !is.na(Survived_1),
         Has_2 = !is.na(Survived_2)) %>%
  summarize(only_1 = sum(Has_1 & !Has_2),
            only_2 = sum(!Has_1 & Has_2),
            both = sum(Has_1 & Has_2))

survival_dif_results <- masterlist_data %>%
  filter(!is.na(Survived_1), !is.na(Survived_2), Survived_1 > 0) %>%
  mutate(diff = Survived_1 - Survived_2) %>%
  count(diff) %>%
  arrange(diff)
survival_dif_results

nests_with_loss <- masterlist_data %>%
  filter(!is.na(Survived_1), !is.na(Survived_2), Survived_1 > 0) %>%
  mutate(diff = Survived_1 - Survived_2) %>%
  filter(diff > 0) %>%  
  select(Nest_ID, Survived_1, Survived_2, diff)
print(nests_with_loss)

# Hatch order of those that survived and died between 30 and 60 days
survived <- c(1,6,4,1,1,3,5,1.5,1.5,3)
died <- c(2,3,2,2,4,3,5,1,2,7,3.5)

mean(survived)
sd(survived)

mean(died)
sd(died)

# Pivoting data longer
longer_data1 <- masterlist_data %>%
  pivot_longer(cols = c('Hatched_eggs','Hatch_success'), names_to = 'Data', 
               values_to = 'Value') %>%
  pivot_longer(cols = c("Hatch_spread", "Clutch_size"), names_to = "Predictor", 
               values_to = "Predictor_value")

# Pivoting data longer with two survival estimates. 1) Conservative estimate: 
# When survival unknown, NAs converted to 0s. 2) Liberal estimate: NAs removed.
longer_data2 <- masterlist_data %>%
  pivot_longer(cols = c(Survived_1_cons, Survived_2_cons,
               Relative_survive_success_1, Relative_survive_success_2),
               names_to = "Data", values_to = "Value") %>%
  mutate(converted_na = case_when(
      Data == "Survived_1_cons" ~ converted_na_1,
      Data == "Survived_2_cons" ~ converted_na_2,
      Data == "Relative_survive_success_1" ~ Converted_na_1_prop,
      Data == "Relative_survive_success_2" ~ Converted_na_2_prop,
      Value = ifelse(converted_na & is.na(Value), 0, Value),
      Data = factor(Data, levels = c("Survived_1_cons", "Survived_2_cons",
      "Relative_survive_success_1", "Relative_survive_success_2")))) %>%
  pivot_longer(cols = c(Hatch_spread, Clutch_size, Hatched_eggs),
               names_to = "Predictor", values_to = "Predictor_value")

# Plot of number and proportion of clutch surviving to 30 and 60 days. 


longer_data2 <- masterlist_data %>%
  pivot_longer(cols = c(Hatch_success,
                        Relative_survive_success_1, Relative_survive_success_2),
               names_to = "Data", values_to = "Value") %>%
  mutate(Data = factor(Data, levels = c("Hatch_success",
                                 "Relative_survive_success_1", "Relative_survive_success_2"))) %>%
  pivot_longer(cols = c(Hatch_spread, Clutch_size, Hatched_eggs),
               names_to = "Predictor", values_to = "Predictor_value")

# Plot of number and proportion of clutch surviving to 30 and 60 days. 
figure_2 <- ggplot(longer_data2, aes(x = Predictor_value, y = Value)) +
  geom_jitter(width = 0.2, height = 0.05, alpha = 0.8) +
  facet_grid(Data ~ Predictor, scales = "free", labeller = labeller(
    Data = c(Relative_survive_success_2 = "Proportion of brood\nsurvived to 60 days", 
             Hatch_success = "Proportion of clutch\n hatched (%) ", 
             Relative_survive_success_1 = "Proportion of brood\nsurvived to 30 days"), 
    Predictor = c(Clutch_size = "Clutch size", Hatch_spread = "Hatch spread (days)", 
                  Hatched_eggs = "Hatched eggs")), switch = "both") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL)+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  theme(strip.background = element_blank(), strip.placement = "outside", 
        strip.text = element_text(size = 8))
print(figure_2) 

# Function to view model diagnostics 
diagnostics <- function(model) {
  print(check_model(model)) 
  res <- simulateResiduals(fittedModel = model, plot = TRUE)
  testDispersion(res)
  testZeroInflation(res)}
  
# Model analyzing proportion hatched eggs by clutch size, HS, single/joint.
Hatch_model <- glmmTMB(cbind(Hatched_eggs, Clutch_size - Hatched_eggs) ~ 
                  Hatch_spread + Clutch_size + Females +
                  (1|Year), family = betabinomial,
                  data = masterlist_data)

CI_95 <- function(model, odds_ratio = FALSE) {
  coefs <- summary(model)$coefficients$cond
  z <- 1.96 
  estimate <- coefs[, "Estimate"]
  se <- coefs[, "Std. Error"]
  lower <- estimate - z * se
  upper <- estimate + z * se
  result <- data.frame(
    Term = rownames(coefs),
    Estimate = round(estimate, 4),
    SE = round(se, 4), 
    CI_lower = round(lower, 4),
    CI_upper = round(upper,4))}

# Model analyzing conservative estimate of survivors to 30 days by clutch size, 
# females and an interaction between hatch spread and number of hatched eggs.
Cons_30_model <- glmmTMB(cbind(Survived_1_cons, Hatched_eggs - Survived_1_cons) ~ 
                     Hatch_spread*Hatched_eggs + Clutch_size + Females + (1|Year), 
                   family = betabinomial, data = masterlist_data)

# Model analyzing conservative estimate of survivors to 60 days by clutch size, 
# females and an interaction between hatch spread and number of hatched eggs.
Cons_60_model <- glmmTMB(cbind(Survived_2_cons, Hatched_eggs - Survived_2_cons) ~ 
                   Hatch_spread*Hatched_eggs + Clutch_size + Females + (1|Year), 
                   family = betabinomial, data = masterlist_data)

# Model analyzing liberal survivors (liberal estimate) to 30 days by clutch size, 
# females and an interaction between hatch spread and number of hatched eggs.
lib_data_1 <- masterlist_data %>% filter(!is.na(Survived_1))
Lib_30_model <- glmmTMB(cbind(Survived_1, Hatched_eggs - Survived_1) ~ 
                   Hatch_spread*Hatched_eggs + Clutch_size + Females + (1|Year), 
                   family = betabinomial,
                   data = lib_data_1)

lib_data_2 <- masterlist_data %>% filter(!is.na(Survived_2))
Lib_60_model <- glmmTMB(cbind(Survived_2, Hatched_eggs - Survived_2) ~ 
                   Hatch_spread*Hatched_eggs + Clutch_size + Females + (1|Year), 
                   family = betabinomial,
                   data = lib_data_2)

# PREDICTION 2C: Larger clutches have larger hatch spreads
longer_data3 <- masterlist_data %>%
  pivot_longer(cols = c(Clutch_size, Hatched_eggs),
               names_to = "Predictor", values_to = "Predictor_value")

correlation_a <- cor.test(as.numeric(masterlist_data$Hatch_spread), 
                          as.numeric(masterlist_data$Clutch_size), method = "spearman")
correlation_b <- cor.test(as.numeric(masterlist_data$Hatch_spread), 
                          as.numeric(masterlist_data$Hatched_eggs), method = "spearman")
correlation_a
correlation_b

# Adding laying and hatching data
laying_data <- read_excel("Lay_hatch_data.xlsx") %>%
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>%
  filter(Hatch_success > 50) 
  
correlation_lay_hatch <- cor.test(as.numeric(laying_data$Hatch_order), 
                        as.numeric(laying_data$Known_lay_order), method = 'pearson')
print(c(correlation_lay_hatch$estimate, correlation_lay_hatch$conf.int, correlation_lay_hatch$p.value))

# Piece wise model for time in the nest from egg order 1-6 and 6 onwards
laying_data <- laying_data %>%
  mutate(lay_order_c1 = pmin(Known_lay_order, 6),
         lay_order_c2 = pmax(Known_lay_order - 6, 0)) 
         
model_time_in_nest <- glmmTMB(Incubation_period ~ lay_order_c1 + lay_order_c2 + 
                             (1|Nest_ID) + (1|Year), data = laying_data, family = gaussian)

# Figure with hatch order and time in nest by lay order.
laying_data_longer <- laying_data %>%
  pivot_longer(cols = c(Hatch_order, Incubation_period),
               names_to = "Data", values_to = "Value")
  
ggplot(laying_data_longer, aes(x = Known_lay_order, y = Value)) +
  geom_jitter(width = 0.05, height = 0.05, size = 2, alpha = 0.7) +
  labs(y = NULL, x = "Lay order") +
  theme_classic() +
  facet_wrap(~ Data, scales = "free_y", ncol = 1, strip.position = "left", 
            labeller = labeller(Data = c(Hatch_order = "Hatch order", 
            Incubation_period = "Time in nest (days)"))) +  
  scale_x_continuous(breaks = seq(0, 12, by = 0.5),
  labels = ifelse(seq(0, 12, by = 0.5) %% 1 == 0, seq(0, 12, by = 0.5), "")) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  theme(strip.background = element_blank(), strip.placement = "outside", 
        strip.text = element_text(size = 12))

# Model on egg size by known lay order.
egg_size_model <- glmmTMB(Volume_cm_cubed ~ Known_lay_order + (1|Nest_ID) + (1|Year), 
                          family = gaussian, data = laying_data)
summary(egg_size_model)

# Adding chick data 
chick_data <- read_excel("Final_Compiled_Chick_Data.xlsx") %>%
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>%
  filter(!(Mass > 35 | Tarsus > 40 | Shield_to_tip > 24))
# Filtering out data that has unusually large birds
# (likely typos, issues with the measurements or recaptures). 

# Matching data from masterlist_data by nest ID
cols <- c("Clutch_size", "Hatch_spread", "Hatched_eggs")
idx  <- match(chick_data$Nest_ID, masterlist_data$Nest_ID)
chick_data[cols] <- masterlist_data[idx, cols] 

# Pivoting data longer then filtering out Tarsus measurements before 2019
chick_data_longer <- chick_data %>%
  pivot_longer(cols = c('Mass', 'Tarsus', 'Shield_to_tip'), 
               names_to = 'Morphometrics', 
               values_to = 'Value') %>%
  filter(!is.na(Value)) %>%
  filter(!(Morphometrics == "Tarsus" & Year < 2019))

chick_data_survival <- chick_data %>%
  filter(!is.na(Survived) & Survived %in% c(0, 1)) %>%
  mutate(Survived = if (!is.numeric(Survived)) as.numeric(Survived) else Survived)

# Filtering out data for nests that hatched more than one egg.
chick_data_survival <- chick_data_survival %>%
  filter(Hatched_eggs > 1)

# Model with hatch order as sole predictor variable on survival.
model_2_3a <- glmmTMB(as.numeric(Survived) ~ as.numeric(Hatch_order) +
                  (1|Nest_ID) + (1|Year), family = binomial, data = chick_data_survival)

model_2_3 <- glmmTMB(as.numeric(Survived) ~ as.numeric(Hatch_Day) + 
                    (1|Nest_ID) + (1|Year), family = binomial, data = chick_data_survival)

# Model with interaction between hatch order and hatching spread
model_2_3b <- glmmTMB(as.numeric(Survived) ~ as.numeric(Hatch_order)*Hatch_spread +
                     (1|Nest_ID) + (1|Year), family = binomial, data = chick_data_survival)

# Model with interaction between hatch order and clutch size
model_2_3c <- glmmTMB(as.numeric(Survived) ~ as.numeric(Hatch_order)*Clutch_size +
                        (1|Nest_ID) + (1|Year), family = binomial, data = chick_data_survival)

# Model of mass by hatch order
model_mass <- glmmTMB(Mass ~ Hatch_order + (1|Nest_ID) + (1|Year), 
                      family = gaussian, data = chick_data)

# Model of shield to tip length by hatch order
model_StoT <- glmmTMB(Shield_to_tip ~ as.numeric(Hatch_order) + (1|Nest_ID) +
                      (1|Year), family = gaussian, data = chick_data)

# Filtering out data before 2019 because they used a different measurement
tars_data <- chick_data %>%
  filter(!(Year < 2019))

# Model 
model_Tars <- glmmTMB(Tarsus ~ as.numeric(Hatch_order) + (1|Nest_ID) + (1|Year), 
                      family = gaussian, data = tars_data)

# Simulating predicted data from models
# Created a sequence of hatch orders for prediction
pred_data <- data.frame(Hatch_order = seq(1, 10, length.out = 200))

# Get predictions for each model
pred_mass <- pred_data %>%
  mutate(pred = predict(model_mass, newdata = ., re.form = NA, se.fit = TRUE)$fit,
         se = predict(model_mass, newdata = ., re.form = NA, se.fit = TRUE)$se.fit,
         lower = pred - 1.96 * se,
         upper = pred + 1.96 * se,
         Morphometrics = "Mass")

pred_shield <- pred_data %>%
  mutate(pred = predict(model_StoT, newdata = ., re.form = NA, se.fit = TRUE)$fit,
         se = predict(model_StoT, newdata = ., re.form = NA, se.fit = TRUE)$se.fit,
         lower = pred - 1.96 * se,
         upper = pred + 1.96 * se,
         Morphometrics = "Shield_to_tip")

pred_tarsus <- pred_data %>%
  mutate(pred = predict(model_Tars, newdata = ., re.form = NA, se.fit = TRUE)$fit,
         se = predict(model_Tars, newdata = ., re.form = NA, se.fit = TRUE)$se.fit,
         lower = pred - 1.96 * se,
         upper = pred + 1.96 * se,
         Morphometrics = "Tarsus")

pred_all <- bind_rows(pred_mass, pred_shield, pred_tarsus)

# Figure of chick size by hatch order with linear regressions and 95% CI
ggplot(chick_data_longer, aes(x = Hatch_order, y = Value)) +
  geom_jitter(width = 0.04, height = 0.1, size = 1, alpha = 0.7) +
  facet_wrap(~Morphometrics, scales="free_y", ncol = 1, strip.position = "left",
             labeller = labeller(Morphometrics = c(Mass = "Mass (g)",
                                                   Tarsus = "Left outer\ntarsus (mm)", 
                                                   Shield_to_tip = "Shield to\ntip (mm)"))) +
  labs(x = "Hatch order", y = NULL) +
  scale_x_continuous(breaks = seq(0, 10, 0.5),
  labels = ifelse(seq(0, 10, by = 0.5) %% 1 == 0, seq(0, 10, by = 0.5), "")) +
  scale_y_continuous(breaks = int_breaks_5, labels = scales::number_format(accuracy = 1)) +
  theme_classic() +
  geom_ribbon(data = pred_all,
              aes(x = Hatch_order, ymin = lower, ymax = upper),
              fill = "lightblue", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = pred_all,
            aes(x = Hatch_order, y = pred),
            color = "blue", linewidth = 0.75, inherit.aes = FALSE) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  theme(strip.background = element_blank(), strip.placement = "outside", 
        strip.text = element_text(size = 11))

# Model of survival by mass and hatch order interaction
model_2_3d <- glmmTMB(as.numeric(Survived) ~ Mass*Hatch_order +
                     (1|Nest_ID) + (1|Year), family = binomial, 
                      data = chick_data_survival)

# Figure of survival by size at hatching 
chick_data_longer <- chick_data_longer %>%
  filter(!is.na(Survived)) %>%
  filter(! Survived == "NA") %>%
  group_by(Morphometrics, Survived) %>%
  summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop")

chick_data_means <- chick_data_longer %>%
  group_by(Morphometrics, Survived) %>%
  summarise(mean_value = mean(Value, na.rm = TRUE),
            sd = sd(Value, na.rm = TRUE),
            min = min(Value, na.rm = TRUE),
            max = max(Value, na.rm = TRUE),
            n = sum(!is.na(Value)), .groups = "drop")
  
ggplot(chick_data_longer, aes(x = Survived, y = Value)) +
  geom_point(size = 2, alpha = 0.7) +
  facet_wrap(~Morphometrics, scales="free_y", ncol = 1, strip.position = "left",
             labeller = labeller(
  Morphometrics = c(Mass = "Mass (g)", Tarsus = "Left outer tarsus (mm)", 
                    Shield_to_tip = "Shield to tip (mm)"))) +
  geom_line(data = chick_data_means, aes(x = Survived, y = mean_value, group = Morphometrics),
            color = "red", linewidth = 0.6) +
  labs(x = NULL, y = NULL) +
  scale_y_continuous() +
  scale_x_discrete(labels = c("0" = "Died", "1" = "Survived")) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  theme(strip.background = element_blank(), strip.placement = "outside", 
        strip.text = element_text(size = 11))

#Supplementary Figure A1
all_data <- read_excel("Nests_masterlist.xlsx",
                       na = c("", "NO_RECORD", "MISSING")) %>%
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>%  
  filter(Exclusion == "GOOD" | Exclusion == "MISSING_HATCH_ORDER" | Exclusion == "FAILED") %>%
  mutate(Hatch_success = (as.numeric(Hatched_eggs)/as.numeric(Clutch_size)*100)) %>%
  mutate(end_date = if_else(Exclusion == "FAILED", Last_observed, Hatch_begin),
         obs_nesting_period = as.numeric(end_date - Date_found, units = "days"))

ggplot(all_data, aes(x = obs_nesting_period, y = Hatch_success)) +
  geom_jitter(width = 1, height = 0.5, alpha = 0.8) +
  theme_classic() +
  labs(x = "Observed nesting period (days)", y = "Hatch success (%)")


##############################
####      CHAPTER 3       ####
#### SYNCHRONY EXPERIMENT ####
##############################

# Nests 2018_E, 2018_AD, 2018_AZ, 2018_BD were not manipulated and thus excluded. 
# Nest 2024_W predated during hatching, also excluded. 
experiment_data <- read_excel("Compiled_synchrony_experiment_data.xlsx") %>%
  mutate(Nest_ID = paste(Year, Nest, sep = "_")) %>%
  filter(!Treatment %in% c("Control","Other")) %>%
  filter(!Nest_ID %in% c("2018_E","2018_AD","2018_AZ","2018_BD", "2023_C")) %>%  
  mutate(Treatment = recode(Treatment, "Synch" = "Synchronous",
                                       "Asynch" = "Asynchronous"), 
         Hatch_success = (Hatched/Manipulated_clutch_size), 
         Foreign_percentage = (Foreign_eggs/Manipulated_clutch_size), 
         Survive_success_brood = (Survival_60/Hatched), 
         Survive_success_clutch = (Survival_60/Manipulated_clutch_size))

# Further cleaning of the data set with successful nests
successful_nests <- experiment_data %>% 
  filter(Hatch_success > 0) %>%
  filter(!Nest_ID %in% "2024_W") %>%
  mutate(Date_found = as.Date(Date_found),
         Hatch_begin = as.Date(Hatch_begin), 
         Observed_nesting_period = as.numeric(Hatch_begin - Date_found)) %>%
  mutate(Date_transfer = as.Date(Date_transfer),
         Transfered_time = as.numeric(Hatch_begin - Date_transfer)) %>%
  mutate(Estimated_hatch_spread = as.numeric(Estimated_hatch_spread),
         True_hatch_spread = as.numeric(True_hatch_spread))

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
summary_successful_MT <- summarize_data(successful_nests_long_MT)

lapply(list(summary_experiment_data, summary_successful, 
            summary_successful_HS, summary_successful_MT), print, n = 30)

# T-TESTS analyzing the effectiveness of the experimental design
t.test(Foreign_percentage~Treatment, data = experiment_data) # proportion swapped
t.test(Manipulated_clutch_size~Treatment, data = experiment_data) # Clutch size
t.test(Observed_nesting_period~Treatment, data = successful_nests) # Nesting period
t.test(True_hatch_spread~Treatment, data = successful_nests_HS) # HS
t.test(Transfered_time~Treatment, data = successful_nests_MT) # Transferred time

# Model analyzing number of hatched eggs by clutch size
model_1 <- glmmTMB(cbind(Hatched, Manipulated_clutch_size - Hatched) ~ 
                     Treatment + Manipulated_clutch_size + Observed_nesting_period
                     + (1|Year), family = betabinomial, data = successful_nests)

# Model analyzing number of survivors by treatment.
model_2 <- glmmTMB(cbind(Survival_60, Hatched - Survival_60) ~ 
                   Treatment + Hatched + Manipulated_clutch_size + (1|Year), 
                   family = betabinomial, data = successful_nests)

# Survival to 60 days pivoted longer format for the figure.
survival_data_long <- successful_nests %>%
  mutate(NonSurvivors = Hatched - Survival_60) %>% 
  pivot_longer(cols = c(Survival_60, NonSurvivors), 
               names_to = "Status", values_to = "Count") %>%
  mutate(Survival_60 = ifelse(Status == "Survival_60", 1, 0)) %>%
  uncount(Count) %>%
  mutate(Treatment_survival = paste(Treatment, Survival_60, sep = "_")) 

# Bar plot of survival data.
survival_data_long$Survival_60 <- factor(survival_data_long$Survival_60, levels = c(0, 1), 
                                      labels = c("Died", "Survived"))
ggplot(survival_data_long, 
       aes(x = Treatment, fill = Treatment, pattern = Survival_60)) +
  geom_bar_pattern(position = position_dodge(width = 0.9),
                   stat = "count",
                   pattern_color = NA,
                   pattern_fill = "grey",
                   pattern_density = 0.4,
                   pattern_spacing = 0.015,
                   pattern_size = 1.0,
                   pattern_angle = 45) +
  theme_classic() +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  theme(axis.text.x = element_text(size = 11, colour = "black")) +
  geom_text(stat = "count", aes(label = ..count..), 
            position = position_dodge(width = 0.9), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("Asynchronous" = "blue", "Synchronous" = "firebrick1")) +
  guides(fill = "none") + 
  scale_pattern_manual(values = c("Survived" = "stripe", "Died" = "none", 
                                  labels = c("Survived" = "Survived", "Died" = "Died"),
                                  name = "Outcome")) +
  labs(x = "", y = "Count", fill = "Treatment", pattern = "Fate") +
  theme(legend.position = "right", legend.box = "vertical")

# Scatterplot of number of hatched eggs and survivors by treatment.
successful_nests_long <- successful_nests %>%
  pivot_longer(cols = c(Hatched, Survival_60), 
               names_to = "Success", values_to = "Count")

ggplot(successful_nests_long, aes(x = Treatment, y = Count, colour = Treatment)) +
  geom_beeswarm(cex = 4.5) +
  theme_classic() +
  scale_colour_manual(values = c("Asynchronous" = "blue", "Synchronous" = "firebrick1")) +
  facet_wrap(~Success, ncol = 1, strip.position = "left",
  labeller = labeller(
  Success = c(Hatched = "Hatched eggs", Survival_60 = "Survivors"))) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  theme(strip.background = element_blank(), strip.placement = "outside", 
        strip.text = element_text(size = 12)) +
  labs(x = "", y = "") +
  scale_y_continuous(limits = c(0, 7), breaks = seq(0, 7, by = 1)) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 12, colour = "black"))

# Data on if they were fostered or non-fostered. 
swapping_data <- read_excel("Egg_swapping.xlsx") %>%
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>%
  filter(Number_nest_transfers == 1) %>%
  rename(Hatched_yes = Hatched) %>%
  left_join(successful_nests %>% select(Nest_ID, Manipulated_clutch_size, Hatched),
            by = "Nest_ID")

# Figure of hatch order by egg status and nest treatment.
ggplot(swapping_data, aes(x = Status, y = Hatch_order, colour = Treatment)) +
  geom_beeswarm(cex = 4.5) +
  theme_classic() +
  labs(y = "Hatch Order", x = "") +
  theme(axis.text.x = element_text(size = 11, colour = "black")) +
  scale_y_continuous(breaks = seq(0, 8, by = 0.5),
                     labels = ifelse(seq(0, 8, by = 0.5) %% 
                                       1 == 0, seq(0, 8, by = 0.5), "")) +
  scale_x_discrete(expand = expansion(mult = c(1, 1)), labels = c("Original" = "Non-fostered", "Swapped" = "Fostered")) +
  guides(color = guide_legend(title.position = "top", direction = "vertical")) +
  theme(legend.position = "right", legend.box = "vertical") +
  scale_color_manual(values=c("blue", "red"))

ggplot(swapping_data, aes(x = Treatment, y = Hatch_order)) +
  geom_beeswarm(cex = 2) +
  theme_classic() +
  labs(y = "Hatch Order", x = "") +
  theme(axis.text.x = element_text(size = 11, colour = "black")) +
  scale_y_continuous(breaks = seq(0, 8, by = 0.5),
                     labels = ifelse(seq(0, 8, by = 0.5) %% 
                                       1 == 0, seq(0, 8, by = 0.5), ""))

### Add clutch size to swapping data
Swap_hatch_model <- glmmTMB(Hatched_yes ~ Status*Treatment + (1|Nest_ID) +
                     + (1|Year), family = binomial, data = swapping_data)

Swap_survive_model <- glmmTMB(Survived ~ Status*Treatment + Hatch_order*Treatment + (1|Nest_ID) +
                     + (1|Year), family = binomial, data = swapping_data)

model <- glmmTMB(Hatch_order ~ Status*Treatment + (1|Nest_ID) + (1|Year), 
                 family = gaussian, data = swapping_data)
