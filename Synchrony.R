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
library(cowplot)
library(emmeans)

############################
####      CHAPTER 2     ####
#### OBSERVATIONAL DATA ####
############################

# Number of excluded and included groups
nest_counts <- read_excel("Nests_masterlist.xlsx") %>%
  count(Exclusion, name = "nest_counts") %>%
  arrange(desc("nest_counts"))
print(nest_counts)


mayfield_data <- read_excel("Nests_masterlist.xlsx",
                 na = c("", "NO_RECORD", "MISSING")) %>%
         mutate(Hatch_begin = as.Date(Hatch_begin),
                Hatch_end = as.Date(Hatch_end),
                Hatch_spread = as.numeric(Hatch_end - Hatch_begin + 1),
                Date_found = as.Date(Date_found), 
                Last_observed = as.Date(Last_observed), 
                Observed_nesting_period = ifelse(
                Hatched_eggs > 0, 
                as.numeric(Hatch_begin - Date_found), 
                as.numeric(Last_observed - Date_found)),
                Hatch_success = as.numeric(Hatched_eggs) / as.numeric(Clutch_size)) %>%
         filter(!is.na(Observed_nesting_period),
                !is.na(Hatched_eggs))

# Extract observation days (exposure)
Obs_period <- mayfield_data$Observed_nesting_period

# Define failures: nests with 0 hatched eggs
Failures <- ifelse(mayfield_data$Hatched_eggs == 0, 1, 0)

# Total nest-days at risk
total_days <- sum(Obs_period, na.rm = TRUE)

# Total failures
total_failures <- sum(Failures, na.rm = TRUE)

# Daily failure and survival rates
daily_failure <- total_failures / total_days
daily_survival <- 1 - daily_failure

daily_failure
daily_survival
prob_nest_survival <- daily_survival^28.6




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

# Function to view model diagnostics 
diagnostics <- function(model) {
  print(check_model(model)) 
  res <- simulateResiduals(fittedModel = model, plot = TRUE)
  testDispersion(res)
  testZeroInflation(res)}

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

##### FIGURE 2.1 #####
##### FIGURE 2.1A Model Predictions #####
Survive_60_clutch_model <- glmmTMB(cbind(as.numeric(Survived_2), as.numeric(Hatched_eggs) - as.numeric(Survived_2)) ~ 
                                     as.numeric(Clutch_size) + (1|Year), family = betabinomial,
                                   data = masterlist_data)
Clutch_size_mod <- 2:18

pred <- emmeans(Survive_60_clutch_model, ~ Clutch_size, at=list(Clutch_size=Clutch_size_mod), 
                type="response", mode="asymptotic")

pred_df <- as.data.frame(pred)
pred_df$Predicted <- pmin(pmax(pred_df$response * pred_df$Clutch_size, 0), pred_df$Clutch_size)
pred_df$Lower <- pmin(pmax(pred_df$asymp.LCL * pred_df$Clutch_size, 0), pred_df$Clutch_size)
pred_df$Upper <- pmin(pmax(pred_df$asymp.UCL * pred_df$Clutch_size, 0), pred_df$Clutch_size)
pred_df

round_preds <- function(pred_df, cols = c("Predicted", "Lower", "Upper"), digits = 2) {
  pred_df[cols] <- lapply(pred_df[cols], function(x) round(x, digits))
  return(pred_df)}

round_preds(pred_df)

predicted_21A <- list(geom_ribbon(data = pred_df,
                                aes(x = Clutch_size, ymin = Lower, ymax = Upper),
                                fill = "lightblue", alpha = 0.3, inherit.aes = FALSE),
                    geom_line(data = pred_df,
                              aes(x = Clutch_size, y = Predicted),
                              colour = "blue", size = 1, inherit.aes = FALSE))

##### FIGURE 2.1B Model Predictions #####
Survive_60_egg_model <- glmmTMB(cbind(as.numeric(Survived_2), as.numeric(Hatched_eggs) - as.numeric(Survived_2)) ~ 
                                  as.numeric(Hatched_eggs) + (1|Year), family = betabinomial,
                                data = masterlist_data)
Hatched_eggs_mod <- 1:10

pred <- emmeans(Survive_60_egg_model, ~ Hatched_eggs, at=list(Hatched_eggs=Hatched_eggs_mod), 
                type="response", mode="asymptotic")

pred_df <- as.data.frame(pred)
pred_df$Predicted <- pmin(pmax(pred_df$response * pred_df$Hatched_eggs, 0), pred_df$Hatched_eggs)
pred_df$Lower <- pmin(pmax(pred_df$asymp.LCL * pred_df$Hatched_eggs, 0), pred_df$Hatched_eggs)
pred_df$Upper <- pmin(pmax(pred_df$asymp.UCL * pred_df$Hatched_eggs, 0), pred_df$Hatched_eggs)
pred_df

predicted_21B <- list(geom_ribbon(data = pred_df,
                                aes(x = Hatched_eggs, ymin = Lower, ymax = Upper),
                                fill = "lightblue", alpha = 0.3, inherit.aes = FALSE),
                    geom_line(data = pred_df,
                              aes(x = Hatched_eggs, y = Predicted),
                              colour = "blue", size = 1, inherit.aes = FALSE))

##### FIGURE 2.1C Model Predictions #####
Survive_60_spread_model <- glmmTMB(cbind(as.numeric(Survived_2), as.numeric(Hatched_eggs) - as.numeric(Survived_2)) ~ 
                                     as.numeric(Hatch_spread) + (1|Year), family = betabinomial,
                                   data = masterlist_data)
Hatch_spread_mod <- 1:11

pred <- emmeans(Survive_60_spread_model, ~ Hatch_spread, at=list(Hatch_spread=Hatch_spread_mod), 
                type="response", mode="asymptotic")

pred_df <- as.data.frame(pred)
pred_df$Predicted <- pmin(pmax(pred_df$response * pred_df$Hatch_spread, 0), pred_df$Hatch_spread)
pred_df$Lower <- pmin(pmax(pred_df$asymp.LCL * pred_df$Hatch_spread, 0), pred_df$Hatch_spread)
pred_df$Upper <- pmin(pmax(pred_df$asymp.UCL * pred_df$Hatch_spread, 0), pred_df$Hatch_spread)
pred_df

predicted_21C <- list(geom_ribbon(data = pred_df,
                                aes(x = Hatch_spread, ymin = Lower, ymax = Upper),
                                fill = "lightblue", alpha = 0.3, inherit.aes = FALSE),
                    geom_line(data = pred_df,
                              aes(x = Hatch_spread, y = Predicted),
                              colour = "blue", size = 1, inherit.aes = FALSE))

##### FIGURE 2.1D Model Predictions #####
Survive_60_30_model <- glmmTMB(cbind(as.numeric(Survived_2), as.numeric(Survived_1) - as.numeric(Survived_2)) ~ 
                                     as.numeric(Survived_1) + (1|Year), family = betabinomial,
                                   data = masterlist_data)
Survived_1_mod <- 0:6

pred <- emmeans(Survive_60_30_model, ~ Survived_1, at=list(Survived_1=Survived_1_mod), 
                type="response", mode = "latent")
pred_df <- as.data.frame(pred)
pred_df$Predicted <- pmin(pmax(pred_df$response * pred_df$Survived_1, 0), pred_df$Survived_1)
pred_df$Lower <- pmin(pmax(pred_df$asymp.LCL * pred_df$Survived_1, 0), pred_df$Survived_1)
pred_df$Upper <- pmin(pmax(pred_df$asymp.UCL * pred_df$Survived_1, 0), pred_df$Survived_1)
pred_df

predicted_21D <- list(geom_ribbon(data = pred_df,
                                aes(x = Survived_1, ymin = Lower, ymax = Upper),
                                fill = "lightblue", alpha = 0.3, inherit.aes = FALSE),
                    geom_line(data = pred_df,
                              aes(x = Survived_1, y = Predicted),
                              colour = "blue", size = 1, inherit.aes = FALSE))

##### FIGURE 2.1E Model Predictions #####
Survive_30_clutch_model <- glmmTMB(cbind(as.numeric(Survived_1), as.numeric(Hatched_eggs) - as.numeric(Survived_1)) ~ 
                                     as.numeric(Clutch_size) + (1|Year), family = betabinomial,
                                   data = masterlist_data)
Clutch_size_mod <- 2:18

pred <- emmeans(Survive_30_clutch_model, ~ Clutch_size, at=list(Clutch_size=Clutch_size_mod), 
                type="response", mode="asymptotic")
pred_df <- as.data.frame(pred)

pred_df$Predicted <- pmin(pmax(pred_df$response * pred_df$Clutch_size, 0), pred_df$Clutch_size)
pred_df$Lower <- pmin(pmax(pred_df$asymp.LCL * pred_df$Clutch_size, 0), pred_df$Clutch_size)
pred_df$Upper <- pmin(pmax(pred_df$asymp.UCL * pred_df$Clutch_size, 0), pred_df$Clutch_size)
pred_df

predicted_21E <- list(geom_ribbon(data = pred_df,
                                aes(x = Clutch_size, ymin = Lower, ymax = Upper),
                                fill = "lightblue", alpha = 0.3, inherit.aes = FALSE),
                    geom_line(data = pred_df,
                              aes(x = Clutch_size, y = Predicted),
                              colour = "blue", size = 1, inherit.aes = FALSE))

##### FIGURE 2.1F Model Predictions #####
Survive_30_egg_model <- glmmTMB(cbind(as.numeric(Survived_1), as.numeric(Hatched_eggs) - as.numeric(Survived_1)) ~ 
                                     as.numeric(Hatched_eggs) + (1|Year), family = betabinomial,
                                   data = masterlist_data)
Hatched_eggs_mod <- 1:11

pred <- emmeans(Survive_30_egg_model, ~ Hatched_eggs, at=list(Hatched_eggs=Hatched_eggs_mod), 
                type="response", mode="asymptotic")
pred_df <- as.data.frame(pred)
pred_df$Predicted <- pmin(pmax(pred_df$response * pred_df$Hatched_eggs, 0), pred_df$Hatched_eggs)
pred_df$Lower <- pmin(pmax(pred_df$asymp.LCL * pred_df$Hatched_eggs, 0), pred_df$Hatched_eggs)
pred_df$Upper <- pmin(pmax(pred_df$asymp.UCL * pred_df$Hatched_eggs, 0), pred_df$Hatched_eggs)
pred_df

predicted_21F <- list(geom_ribbon(data = pred_df,
                                aes(x = Hatched_eggs, ymin = Lower, ymax = Upper),
                                fill = "lightblue", alpha = 0.3, inherit.aes = FALSE),
                    geom_line(data = pred_df,
                              aes(x = Hatched_eggs, y = Predicted),
                              colour = "blue", size = 1, inherit.aes = FALSE))

##### FIGURE 2.1G Model Predictions #####
Survive_30_spread_model <- glmmTMB(cbind(as.numeric(Survived_1), as.numeric(Hatched_eggs) - as.numeric(Survived_1)) ~ 
                                     as.numeric(Hatch_spread) + (1|Year), family = betabinomial,
                                   data = masterlist_data)
Hatch_spread_mod <- 1:11

pred <- emmeans(Survive_30_spread_model, ~ Hatch_spread, at=list(Hatch_spread=Hatch_spread_mod), 
                type="response", mode="asymptotic")
pred_df <- as.data.frame(pred)
pred_df$Predicted <- pmin(pmax(pred_df$response * pred_df$Hatch_spread, 0), pred_df$Hatch_spread)
pred_df$Lower <- pmin(pmax(pred_df$asymp.LCL * pred_df$Hatch_spread, 0), pred_df$Hatch_spread)
pred_df$Upper <- pmin(pmax(pred_df$asymp.UCL * pred_df$Hatch_spread, 0), pred_df$Hatch_spread)
pred_df

predicted_21G <- list(geom_ribbon(data = pred_df,
                                aes(x = Hatch_spread, ymin = Lower, ymax = Upper),
                                fill = "lightblue", alpha = 0.3, inherit.aes = FALSE),
                    geom_line(data = pred_df,
                              aes(x = Hatch_spread, y = Predicted),
                              colour = "blue", size = 1, inherit.aes = FALSE))

##### FIGURE 2.1H Model Predictions #####
Spread_clutch_model <- glmmTMB(as.numeric(Hatch_spread) ~ as.numeric(Clutch_size) 
                               + (1|Year), family = poisson,
                               data = masterlist_data)
Clutch_size_mod <- 2:18

pred <- emmeans(Spread_clutch_model, ~ Clutch_size, at=list(Clutch_size=Clutch_size_mod), 
                type="response", mode="asymptotic")
pred_df <- as.data.frame(pred)
pred_df

predicted_21H <- list(geom_ribbon(data = pred_df,
                                aes(x = Clutch_size, ymin = asymp.LCL, ymax = asymp.UCL),
                                fill = "lightblue", alpha = 0.3, inherit.aes = FALSE),
                    geom_line(data = pred_df,
                              aes(x = Clutch_size, y = response),
                              colour = "blue", size = 1, inherit.aes = FALSE))


##### FIGURE 2.1I Model Predictions #####
Spread_egg_model <- glmmTMB(as.numeric(Hatch_spread) ~ as.numeric(Hatched_eggs) 
                            + (1|Year), family = poisson,
                            data = masterlist_data)
Hatched_eggs_mod <- 1:10

pred <- emmeans(Spread_egg_model, ~ Hatched_eggs, at=list(Hatched_eggs=Hatched_eggs_mod), 
                type="response", mode="asymptotic")
pred_df <- as.data.frame(pred)
pred_df

predicted_21I <- list(geom_ribbon(data = pred_df,
                                aes(x = Hatched_eggs, ymin = asymp.LCL, ymax = asymp.UCL),
                                fill = "lightblue", alpha = 0.3, inherit.aes = FALSE),
                    geom_line(data = pred_df,
                              aes(x = Hatched_eggs, y = response),
                              colour = "blue", size = 1, inherit.aes = FALSE))


##### FIGURE 2.1J Model Predictions #####
Hatch_clutch_model <- glmmTMB(cbind(as.numeric(Hatched_eggs), as.numeric(Clutch_size) - as.numeric(Hatched_eggs)) ~ 
                                as.numeric(Clutch_size) + (1|Year), family = betabinomial,
                              data = masterlist_data)
Clutch_size_mod <- 2:18

pred <- emmeans(Hatch_clutch_model, ~ Clutch_size, at=list(Clutch_size=Clutch_size_mod), 
                type="response", mode="asymptotic")
pred_df <- as.data.frame(pred)
pred_df$Predicted <- pmin(pmax(pred_df$response * pred_df$Clutch_size, 0), pred_df$Clutch_size)
pred_df$Lower <- pmin(pmax(pred_df$asymp.LCL * pred_df$Clutch_size, 0), pred_df$Clutch_size)
pred_df$Upper <- pmin(pmax(pred_df$asymp.UCL * pred_df$Clutch_size, 0), pred_df$Clutch_size)
pred_df

predicted_21J <- 
  list(geom_ribbon(data = pred_df,
                aes(x = Clutch_size, ymin = Lower, ymax = Upper),
                fill = "lightblue", alpha = 0.3, inherit.aes = FALSE),
  geom_line(data = pred_df,
              aes(x = Clutch_size, y = Predicted),
              colour = "blue", size = 1, inherit.aes = FALSE))


##### Figure 2.1 Putting it all together #####
var_labels <- c(Clutch_size = "Clutch size",
                Hatched_eggs = "Hatched eggs",
                Hatch_spread = "Hatch spread (days)",
                Survived_1 = "Survived 30 days",
                Survived_2 = "Survived 60 days")

axis_breaks <- list(Clutch_size = seq(2, 18, 4),
                    Hatched_eggs = seq(1, 11, 2),
                    Hatch_spread = seq(1, 11, 2),
                    Survived_1 = seq(0, 13, 3),
                    Survived_2 = seq(0, 10, 2))

axis_limits <- list(Clutch_size = c(1, 19),
                    Hatched_eggs  = c(0, 12),
                    Hatch_spread  = c(0, 12),
                    Survived_1    = c(0, 13),
                    Survived_2    = c(0, 10))

make_plot <- function(data, xvar, yvar, label = NULL) {
  ggplot(data, aes_string(xvar, yvar)) +
    geom_jitter(width = 0.2, height = 0.1, alpha = 0.8) +
    labs(x = var_labels[xvar], y = var_labels[yvar], title = label) +
    theme_classic(base_size = 9) +  
    scale_x_continuous(breaks = axis_breaks[[xvar]], limits = axis_limits[[xvar]]) +
    scale_y_continuous(breaks = axis_breaks[[yvar]], limits = axis_limits[[yvar]], labels = function(x) sprintf("%2d", x)) +
    theme(plot.title = element_text(hjust = 0, face = "bold", color = "black", size = 12, margin = margin(b=1)),
          axis.text  = element_text(size = 8, color = "black"),
          axis.title = element_text(size = 10), 
          panel.spacing = unit(0.01, "lines"), 
          plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5))}

invisible_x <- theme(axis.title.x = element_text(color = "white"), axis.ticks.x = element_blank(),
                     axis.text.x  = element_text(color = "white"))
invisible_y <- theme(axis.title.y = element_text(color = "white"), axis.ticks.y = element_blank(),
                     axis.text.y  = element_text(color = "white"))
only_left <- theme(axis.title.y = element_text(), axis.text.y  = element_text(),
                   axis.ticks.y = element_line())
only_bottom <- theme(axis.title.x = element_text(), axis.text.x  = element_text(),
                     axis.ticks.x = element_line())

pA21 <- make_plot(masterlist_data, "Clutch_size", "Survived_2") + only_left + 
  invisible_x + ggtitle("A") + predicted_A21
pB21 <- make_plot(masterlist_data, "Hatched_eggs", "Survived_2") + invisible_y + 
  invisible_x + ggtitle("B") + predicted_B21
pC21 <- make_plot(masterlist_data, "Hatch_spread", "Survived_2") + invisible_y + 
  invisible_x + ggtitle("C") + predicted_C21
pD21 <- make_plot(masterlist_data, "Survived_1", "Survived_2") + invisible_y + 
  only_bottom + ggtitle("D") + predicted_D21 + scale_x_continuous(limits = c(0,6))
pE21 <- make_plot(masterlist_data, "Clutch_size", "Survived_1") + only_left + 
  invisible_x + ggtitle("E") + predicted_E21
pF21 <- make_plot(masterlist_data, "Hatched_eggs", "Survived_1") + invisible_y + 
  invisible_x + ggtitle("F") + predicted_F21
pG21 <- make_plot(masterlist_data, "Hatch_spread", "Survived_1") + invisible_y + 
  only_bottom + ggtitle("G") + predicted_G21
pH21 <- make_plot(masterlist_data, "Clutch_size", "Hatch_spread") + only_left + 
  invisible_x + ggtitle("H") + predicted_H21
pI21 <- make_plot(masterlist_data, "Hatched_eggs", "Hatch_spread") + invisible_y + 
  only_bottom + ggtitle("I") + predicted_I21
pJ21 <- make_plot(masterlist_data, "Clutch_size", "Hatched_eggs") + only_left + 
  only_bottom + ggtitle("J") + predicted_J21

row1 <- plot_grid(pA21, pB21, pC21, pD21, ncol = 4)
row2 <- plot_grid(pE21, pF21, pG21, NULL, ncol = 4)
row3 <- plot_grid(pH21, pI21, NULL, NULL, ncol = 4)
row4 <- plot_grid(pJ21, NULL, NULL, NULL, ncol = 4)

figure2_1 <- plot_grid(row1, row2, row3, row4, ncol = 1, align = "v") 
figure2_1


##### Figure 2.2 #####

# Panel 2.2A
Survive_60_clutch_pred_data <- data.frame(Clutch_size = seq(1, 20, length.out = 200))

pred_survive_60_clutch_model <- Survive_60_clutch_pred_data %>%
  mutate(pred = predict(Survive_60_clutch_model, type = "response", newdata = ., re.form = NA, se.fit = TRUE)$fit,
         se = predict(Survive_60_clutch_model, type = "response", newdata = ., re.form = NA, se.fit = TRUE)$se.fit,
         lower  = pmax(0, pred - 1.96 * se) * 100,
         upper  = pmin(1, pred + 1.96 * se) * 100,
         pred   = pred * 100)

# Panel 2.2B
Survive_60_egg_pred_data <- data.frame(Hatched_eggs = seq(1, 11, length.out = 200))

pred_survive_60_egg_model <- Survive_60_egg_pred_data %>%
  mutate(pred = predict(Survive_60_egg_model, type = "response", newdata = ., re.form = NA, se.fit = TRUE)$fit,
         se = predict(Survive_60_egg_model, type = "response", newdata = ., re.form = NA, se.fit = TRUE)$se.fit,
         lower  = pmax(0, pred - 1.96 * se) * 100,
         upper  = pmin(1, pred + 1.96 * se) * 100,
         pred   = pred * 100)

# Panel 2.2C
Survive_60_spread_pred_data <- data.frame(Hatch_spread = seq(1, 11, length.out = 200))

pred_survive_60_spread_model <- Survive_60_spread_pred_data %>%
  mutate(pred = predict(Survive_60_spread_model, type = "response", newdata = ., re.form = NA, se.fit = TRUE)$fit,
         se = predict(Survive_60_spread_model, type = "response", newdata = ., re.form = NA, se.fit = TRUE)$se.fit,
         lower  = pmax(0, pred - 1.96 * se) * 100,
         upper  = pmin(1, pred + 1.96 * se) * 100,
         pred   = pred * 100)

# Panel 2.2D
Survive_30_clutch_model <- glmmTMB(cbind(as.numeric(Survived_1), as.numeric(Hatched_eggs) - as.numeric(Survived_1)) ~ 
                                     as.numeric(Clutch_size) + (1|Year), family = betabinomial,
                                   data = masterlist_data)

Survive_30_clutch_pred_data <- data.frame(Clutch_size = seq(1, 20, length.out = 200))

pred_survive_30_clutch_model <- Survive_30_clutch_pred_data %>%
  mutate(pred = predict(Survive_30_clutch_model, type = "response", newdata = ., re.form = NA, se.fit = TRUE)$fit,
         se = predict(Survive_30_clutch_model, type = "response", newdata = ., re.form = NA, se.fit = TRUE)$se.fit,
         lower  = pmax(0, pred - 1.96 * se) * 100,
         upper  = pmin(1, pred + 1.96 * se) * 100,
         pred   = pred * 100)

#Panel 2.2E
Survive_30_egg_pred_data <- data.frame(Hatched_eggs = seq(1, 11, length.out = 200))

pred_survive_30_egg_model <- Survive_30_egg_pred_data %>%
  mutate(pred = predict(Survive_30_egg_model, type = "response", newdata = ., re.form = NA, se.fit = TRUE)$fit,
         se = predict(Survive_30_egg_model, type = "response", newdata = ., re.form = NA, se.fit = TRUE)$se.fit,
         lower  = pmax(0, pred - 1.96 * se) * 100,
         upper  = pmin(1, pred + 1.96 * se) * 100,
         pred   = pred * 100)

# Panel 2.2F
Survive_30_spread_pred_data <- data.frame(Hatch_spread = seq(1, 11, length.out = 200))

pred_survive_30_spread_model <- Survive_30_spread_pred_data %>%
  mutate(pred = predict(Survive_30_spread_model, type = "response", newdata = ., re.form = NA, se.fit = TRUE)$fit,
         se = predict(Survive_30_spread_model, type = "response", newdata = ., re.form = NA, se.fit = TRUE)$se.fit,
         lower  = pmax(0, pred - 1.96 * se) * 100,
         upper  = pmin(1, pred + 1.96 * se) * 100,
         pred   = pred * 100)

# Panel 2.2G
Hatch_clutch_pred_data <- data.frame(Clutch_size = seq(2, 18, length.out = 200))

pred_hatch_clutch_model <- Hatch_clutch_pred_data %>%
  mutate(pred = predict(Hatch_clutch_model, type = "response", newdata = ., re.form = NA, se.fit = TRUE)$fit,
         se = predict(Hatch_clutch_model, type = "response", newdata = ., re.form = NA, se.fit = TRUE)$se.fit,
         lower  = pmax(0, pred - 1.96 * se) * 100,
         upper  = pmin(1, pred + 1.96 * se) * 100,
         pred   = pred * 100)

# Panel 2.2H
Hatch_egg_model <- glmmTMB(cbind(as.numeric(Hatched_eggs), as.numeric(Clutch_size) - as.numeric(Hatched_eggs)) ~ 
                                as.numeric(Hatched_eggs) + (1|Year), family = betabinomial,
                              data = masterlist_data)

hatch_egg_pred_data <- data.frame(Hatched_eggs = seq(1, 11, length.out = 200))

pred_hatch_egg_model <- hatch_egg_pred_data %>%
  mutate(pred = predict(Hatch_egg_model, type = "response", newdata = ., re.form = NA, se.fit = TRUE)$fit,
         se = predict(Hatch_egg_model, type = "response", newdata = ., re.form = NA, se.fit = TRUE)$se.fit,
         lower  = pmax(0, pred - 1.96 * se) * 100,
         upper  = pmin(1, pred + 1.96 * se) * 100,
         pred   = pred * 100)

#Panel 2.2I
Hatch_spread_model <- glmmTMB(cbind(as.numeric(Hatched_eggs), as.numeric(Clutch_size) - as.numeric(Hatched_eggs)) ~ 
                             as.numeric(Hatch_spread) + (1|Year), family = betabinomial,
                             data = masterlist_data)

hatch_spread_pred_data <- data.frame(Hatch_spread = seq(1, 11, length.out = 200))

pred_hatch_spread_model <- hatch_spread_pred_data %>%
  mutate(pred = predict(Hatch_spread_model, type = "response", newdata = ., re.form = NA, se.fit = TRUE)$fit,
         se = predict(Hatch_spread_model, type = "response", newdata = ., re.form = NA, se.fit = TRUE)$se.fit,
         lower  = pmax(0, pred - 1.96 * se) * 100,
         upper  = pmin(1, pred + 1.96 * se) * 100,
         pred   = pred * 100)

##### Figure 2.2 Putting it all together #####
var_labels <- c(Clutch_size = "Clutch size",
                Hatched_eggs = "Hatched eggs",
                Hatch_spread = "Hatch spread (days)",
                Hatch_success = "Proportion of clutch\n hatched (%)",
                Relative_survive_success_1 = "Proportion of brood\nsurvived to 30 days",
                Relative_survive_success_2 = "Proportion of brood\nsurvived to 60 days")

axis_breaks <- list(Clutch_size = seq(2, 18, 4),
                    Hatched_eggs = seq(1, 11, 2),
                    Hatch_spread = seq(1, 11, 2), 
                    Hatch_success = seq(0, 100, 25),
                    Relative_survive_success_1 = seq(0, 100, 25),
                    Relative_survive_success_2 = seq(0, 100, 25))

axis_limits <- list(Clutch_size = c(1, 19),
                    Hatched_eggs  = c(0, 12),
                    Hatch_spread  = c(0, 12), 
                    Hatch_success = c(0, 100),
                    Relative_survive_success_1 = c(0, 100),
                    Relative_survive_success_2 = c(0, 100))

make_plot <- function(data, xvar, yvar, label = NULL) {
  ggplot(data, aes_string(xvar, yvar)) +
    geom_jitter(width = 0.2, height = 0.1, alpha = 0.8) +
    labs(x = var_labels[xvar], y = var_labels[yvar], title = label) +
    theme_classic(base_size = 9) +  
    scale_x_continuous(breaks = axis_breaks[[xvar]], limits = axis_limits[[xvar]]) +
    scale_y_continuous(breaks = axis_breaks[[yvar]], limits = axis_limits[[yvar]]) +
    theme(plot.title = element_text(hjust = 0, face = "bold", color = "black", size = 12, margin = margin(b=1)),
          axis.text  = element_text(size = 8, color = "black"),
          axis.title = element_text(size = 10), 
          panel.spacing = unit(0.01, "lines"), 
          plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5))}

add_prediction <- function(pred_df, xvar) {
    list(geom_line(
         data = pred_df,
         aes(x = !!sym(xvar), y = pred),
         color = "blue", size = 1, inherit.aes = FALSE),
    geom_ribbon(
         data = pred_df,
         aes(x = !!sym(xvar), ymin = lower, ymax = upper),
         fill = "lightblue", alpha = 0.3, inherit.aes = FALSE))}

pA22 <- make_plot(masterlist_data, "Clutch_size", "Relative_survive_success_2") + only_left + invisible_x + ggtitle("A") + add_prediction(pred_survive_60_clutch_model, "Clutch_size")
pB22 <- make_plot(masterlist_data, "Hatched_eggs", "Relative_survive_success_2") + invisible_y + invisible_x + ggtitle("B") + add_prediction(pred_survive_60_egg_model, "Hatched_eggs")
pC22 <- make_plot(masterlist_data, "Hatch_spread", "Relative_survive_success_2") + invisible_y + invisible_x + ggtitle("C") + add_prediction(pred_survive_60_spread_model, "Hatch_spread")
pD22 <- make_plot(masterlist_data, "Clutch_size", "Relative_survive_success_1") + only_left + invisible_x + ggtitle("D") + add_prediction(pred_survive_30_clutch_model, "Clutch_size")
pE22 <- make_plot(masterlist_data, "Hatched_eggs", "Relative_survive_success_1") + invisible_y + invisible_x + ggtitle("E") + add_prediction(pred_survive_30_egg_model, "Hatched_eggs")
pF22 <- make_plot(masterlist_data, "Hatch_spread", "Relative_survive_success_1") + invisible_y + invisible_x + ggtitle("F") + add_prediction(pred_survive_30_spread_model, "Hatch_spread")
pG22 <- make_plot(masterlist_data, "Clutch_size", "Hatch_success") + only_left + only_bottom + ggtitle("G") + add_prediction(pred_hatch_clutch_model, "Clutch_size")
pH22 <- make_plot(masterlist_data, "Hatched_eggs", "Hatch_success") + invisible_y + only_bottom + ggtitle("H") + add_prediction(pred_hatch_egg_model, "Hatched_eggs")
pI22 <- make_plot(masterlist_data, "Hatch_spread", "Hatch_success") + invisible_y + only_bottom + ggtitle("I") + add_prediction(pred_hatch_spread_model, "Hatch_spread")

row1 <- plot_grid(pA22, pB22, pC22, ncol = 3)
row2 <- plot_grid(pD22, pE22, pF22, ncol = 3)
row3 <- plot_grid(pG22, pH22, pI22, ncol = 3)

figure_2_2 <- plot_grid(row1, row2, row3, ncol = 1, align = "v") 
figure_2_2


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

# Proportional hatching success model by clutch size, HS, single/joint.
Hatch_model <- glmmTMB(cbind(Hatched_eggs, Clutch_size - Hatched_eggs) ~ 
                         Hatch_spread + Clutch_size + Females +
                         (1|Year), family = betabinomial,
                         data = masterlist_data)
Hatch_model_noyear <- glmmTMB(cbind(Hatched_eggs, Clutch_size - Hatched_eggs) ~ 
                         Hatch_spread + Clutch_size + Females, family = betabinomial,
                       data = masterlist_data)
diagnostics(Hatch_model)
summary(Hatch_model)
print(CI_95(Hatch_model))
anova(Hatch_model, Hatch_model_noyear)

# Model prediction from Hatch_model for number of hatched eggs
Clutch_size_mod <- 2:18
Hatch_spread_mod <- 1:11
Females_mod <- c("Joint female", "Single female")

pred <- emmeans(Hatch_model, ~ Clutch_size + Hatch_spread + Females, 
                at=list(Clutch_size=Clutch_size_mod, Hatch_spread = Hatch_spread_mod,
                        Females = Females_mod),
                type="response")

pred_df <- as.data.frame(pred)

pred_df$Predicted <- pmin(pmax(pred_df$prob * pred_df$Clutch_size, 0), pred_df$Clutch_size)
pred_df$Lower <- pmin(pmax(pred_df$asymp.LCL * pred_df$Clutch_size, 0), pred_df$Clutch_size)
pred_df$Upper <- pmin(pmax(pred_df$asymp.UCL * pred_df$Clutch_size, 0), pred_df$Clutch_size)
pred_df

get_subset_table <- function(pred_df,
                             Clutch_sizes = 2:18,
                             Hatch_spreads = c(1,5,9), 
                             Hatched_eggs = c(1,5,9)) {
  
# Assign female type based on clutch size
pred_df <- pred_df %>%
mutate(Females2 = ifelse(Clutch_size < 6, "Single female", "Joint female"))
  
# If Hatched_eggs column doesn't exist, create it as NA
if(!"Hatched_eggs" %in% names(pred_df)) {
pred_df$Hatched_eggs <- NA}
  
# Filter for selected clutch sizes and hatch spreads
pred_df_sub <- pred_df %>%
    filter(Clutch_size %in% Clutch_sizes,
           Hatch_spread %in% Hatch_spreads) %>%
    filter((Clutch_size < 6 & Females == "Single female") |
             (Clutch_size >= 6 & Females == "Joint female")) %>% 
    select(Clutch_size, Hatch_spread, Females = Females2, Hatched_eggs, Predicted, Lower, Upper)
  
# Round nicely
  pred_df_sub <- pred_df_sub %>%
    mutate(
      Predicted = Predicted*0.5456363,
      Lower = Lower*0.5456363,
      Upper = Upper*0.5456363)
  return(pred_df_sub)}

subset_table_hatch <- get_subset_table(pred_df)
subset_table_hatch
write.table(subset_table_hatch, pipe("pbcopy"), sep = "\t", row.names = FALSE, quote = FALSE)



# Figure 2.3
predicted_hatched_eggs_fig <-  ggplot(subset_table_hatch, aes(x = Clutch_size, y = Predicted,
                               color = factor(Hatch_spread),
                               fill  = factor(Hatch_spread))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, color = NA) +
  scale_color_manual(values = c("1" = "red","5" = "grey35","9" = "blue"),name = "Hatch spread") +
  scale_fill_manual(values = c("1" = "lightpink","5" = "grey70","9" = "lightblue"),name = "Hatch spread") +
  scale_x_continuous(expand = expansion(mult = c(0, .05)), breaks = seq(2,18,2), limits = c(2,18)) +
  scale_y_continuous(expand = expansion(mult = c(0, .05)), breaks = seq(0,8,2), limits = c(0,8)) +
  labs(x = "Clutch size",
       y = "Predicted number of hatched eggs") +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")

predicted_hatched_eggs_fig


# Conservative proportional survival to 60 days model
Cons_60_model <- glmmTMB(cbind(Survived_2_cons, Hatched_eggs - Survived_2_cons) ~ 
                           Hatch_spread*Hatched_eggs + Clutch_size + Females + (1|Year), 
                         family = betabinomial, data = masterlist_data)
Cons_60_model_noyear <- glmmTMB(cbind(Survived_2_cons, Hatched_eggs - Survived_2_cons) ~ 
                           Hatch_spread*Hatched_eggs + Clutch_size + Females, 
                         family = betabinomial, data = masterlist_data)
diagnostics(Cons_60_model)
summary(Cons_60_model)
print(CI_95(Cons_60_model))
anova(Cons_60_model,Cons_60_model_noyear)

Clutch_size_mod <- 2:18
Hatch_spread_mod <- 1:11
Hatched_eggs_mod <- 1:10
Females_mod <- c("Joint female", "Single female")

pred <- emmeans(Cons_60_model, ~ Clutch_size + Hatch_spread + Hatched_eggs + Females, 
                at=list(Clutch_size=Clutch_size_mod, Hatch_spread = Hatch_spread_mod,
                        Hatched_eggs = Hatched_eggs_mod, Females = Females_mod),
                type="response")

# Predicted number of survivors
pred_df <- as.data.frame(pred)
pred_df$Predicted <- pmin(pmax(pred_df$prob * pred_df$Hatched_eggs, 0), pred_df$Hatched_eggs)
pred_df$Lower <- pmin(pmax(pred_df$asymp.LCL * pred_df$Hatched_eggs, 0), pred_df$Hatched_eggs)
pred_df$Upper <- pmin(pmax(pred_df$asymp.UCL * pred_df$Hatched_eggs, 0), pred_df$Hatched_eggs)
pred_df

subset_table_cons_60_model <- get_subset_table(pred_df)
subset_table_cons_60_model
write.table(subset_table_hatch, pipe("pbcopy"), sep = "\t", row.names = FALSE, quote = FALSE)

# Conservative proportional survival to 60 days model
Cons_30_model <- glmmTMB(cbind(Survived_1_cons, Hatched_eggs - Survived_1_cons) ~ 
                           Hatch_spread*Hatched_eggs + Clutch_size + Females + (1|Year), 
                         family = betabinomial, data = masterlist_data)
Cons_30_model_noyear <- glmmTMB(cbind(Survived_1_cons, Hatched_eggs - Survived_1_cons) ~ 
                                  Hatch_spread*Hatched_eggs + Clutch_size + Females, 
                                family = betabinomial, data = masterlist_data)
anova(Cons_30_model, Cons_30_model_noyear)
diagnostics(Cons_30_model)
summary(Cons_30_model)
print(CI_95(Cons_30_model))

pred <- emmeans(Cons_30_model, ~ Clutch_size + Hatch_spread + Hatched_eggs + Females, 
                at=list(Clutch_size=Clutch_size_mod, Hatch_spread = Hatch_spread_mod,
                        Hatched_eggs = Hatched_eggs_mod, Females = Females_mod),
                type="response")

# Predicted number of survivors
pred_df <- as.data.frame(pred)
pred_df$Predicted <- pmin(pmax(pred_df$prob * pred_df$Hatched_eggs, 0), pred_df$Hatched_eggs)
pred_df$Lower <- pmin(pmax(pred_df$asymp.LCL * pred_df$Hatched_eggs, 0), pred_df$Hatched_eggs)
pred_df$Upper <- pmin(pmax(pred_df$asymp.UCL * pred_df$Hatched_eggs, 0), pred_df$Hatched_eggs)
pred_df

subset_table_hatch <- get_subset_table(pred_df)
subset_table_hatch
write.table(subset_table_hatch, pipe("pbcopy"), sep = "\t", row.names = FALSE, quote = FALSE)


# Filtered data converting non-observations to 60 days to 0s
lib_data_2 <- masterlist_data %>% filter(!is.na(Survived_2))

# Liberal proportional survival to 60 days model
Lib_60_model <- glmmTMB(cbind(Survived_2, Hatched_eggs - Survived_2) ~ 
                          Hatch_spread*Hatched_eggs + Clutch_size + Females + (1|Year), 
                        family = betabinomial, data = lib_data_2)
Lib_60_model_noyear <- glmmTMB(cbind(Survived_2, Hatched_eggs - Survived_2) ~ 
                          Hatch_spread*Hatched_eggs + Clutch_size + Females, 
                        family = betabinomial, data = lib_data_2)
diagnostics(Lib_60_model)
summary(Lib_60_model)
print(CI_95(Lib_60_model))
anova(Lib_60_model,Lib_60_model_noyear)

# Predicted number of survivors
pred <- emmeans(Lib_60_model, ~ Clutch_size + Hatch_spread + Females + Hatched_eggs, 
                at=list(Clutch_size=Clutch_size_mod, Hatch_spread = Hatch_spread_mod,
                        Females = Females_mod, Hatched_eggs = Hatched_eggs_mod),
                type="response")

pred_df <- as.data.frame(pred)
pred_df$Predicted <- pmin(pmax(pred_df$prob * pred_df$Hatched_eggs, 0), pred_df$Hatched_eggs)
pred_df$Lower <- pmin(pmax(pred_df$asymp.LCL * pred_df$Hatched_eggs, 0), pred_df$Hatched_eggs)
pred_df$Upper <- pmin(pmax(pred_df$asymp.UCL * pred_df$Hatched_eggs, 0), pred_df$Hatched_eggs)
pred_df

subset_table_hatch <- get_subset_table(pred_df)
subset_table_hatch
write.table(subset_table_hatch, pipe("pbcopy"), sep = "\t", row.names = FALSE, quote = FALSE)


# Filtered data converting non-observations to 60 days to 0s
lib_data_1 <- masterlist_data %>% filter(!is.na(Survived_1))

# Liberal proportional survival to 30 days model
Lib_30_model <- glmmTMB(cbind(Survived_1, Hatched_eggs - Survived_1) ~ 
                          Hatch_spread*Hatched_eggs + Clutch_size + Females + (1|Year), 
                        family = betabinomial,
                        data = lib_data_1)
Lib_30_model_noyear <- glmmTMB(cbind(Survived_1, Hatched_eggs - Survived_1) ~ 
                          Hatch_spread*Hatched_eggs + Clutch_size + Females, 
                        family = betabinomial,
                        data = lib_data_1)
diagnostics(Lib_30_model)
summary(Lib_30_model)
print(CI_95(Lib_30_model))
anova(Lib_30_model, Lib_30_model_noyear)

# Predicted number of survivors
pred <- emmeans(Lib_30_model, ~ Clutch_size + Hatch_spread + Females + Hatched_eggs, 
                at=list(Clutch_size=Clutch_size_mod, Hatch_spread = Hatch_spread_mod,
                        Females = Females_mod, Hatched_eggs = Hatched_eggs_mod),
                type="response")

pred_df <- as.data.frame(pred)
pred_df$Predicted <- pmin(pmax(pred_df$prob * pred_df$Hatched_eggs, 0), pred_df$Hatched_eggs)
pred_df$Lower <- pmin(pmax(pred_df$asymp.LCL * pred_df$Hatched_eggs, 0), pred_df$Hatched_eggs)
pred_df$Upper <- pmin(pmax(pred_df$asymp.UCL * pred_df$Hatched_eggs, 0), pred_df$Hatched_eggs)
pred_df

subset_table_hatch <- get_subset_table(pred_df)
subset_table_hatch
write.table(subset_table_hatch, pipe("pbcopy"), sep = "\t", row.names = FALSE, quote = FALSE)


# Correlations of clutch size with number of hatched eggs
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

laying_data_success <- read_excel("Lay_hatch_data.xlsx") %>%
  filter(Known_lay_order == 1)
mean_incubation <- mean(laying_data_success$Incubation_period, na.rm = TRUE)


# Correlation of lay order and hatch order
correlation_lay_hatch <- cor.test(as.numeric(laying_data$Hatch_order), 
                        as.numeric(laying_data$Known_lay_order), method = 'pearson')
print(c(correlation_lay_hatch$estimate, correlation_lay_hatch$conf.int, correlation_lay_hatch$p.value))

# Piece wise model for time in the nest from egg order 1-6 and 6 onwards
laying_data <- laying_data %>%
  mutate(lay_order_c1 = pmin(Known_lay_order, 6),
         lay_order_c2 = pmax(Known_lay_order - 6, 0)) 
         
model_time_in_nest <- glmmTMB(Incubation_period ~ lay_order_c1 + lay_order_c2 + 
                             (1|Nest_ID) + (1|Year), data = laying_data, family = gaussian)
model_time_in_nest_noyear <- glmmTMB(Incubation_period ~ lay_order_c1 + lay_order_c2 + 
                                (1|Nest_ID), data = laying_data, family = gaussian)
model_time_in_nest_nonest <- glmmTMB(Incubation_period ~ lay_order_c1 + lay_order_c2 + 
                                       (1|Year), data = laying_data, family = gaussian)
diagnostics(model_time_in_nest)
summary(model_time_in_nest)
print(CI_95(model_time_in_nest))
anova(model_time_in_nest, model_time_in_nest_nonest)
anova(model_time_in_nest, model_time_in_nest_noyear)


# Model predictions for linear regressions
Lay_hatch_order_model <- glmmTMB(Hatch_order ~ Known_lay_order + (1|Year), family = gaussian,
                                 data = laying_data)
pred_data <- data.frame(Known_lay_order = seq(1, 11, length.out = 200))

pred_df <- predict(Lay_hatch_order_model, 
                   newdata = pred_data, 
                   type = "response", 
                   re.form = NA, 
                   se.fit = TRUE)

add_prediction <- function(pred_df, xvar) {
  list(geom_line(
    data = pred_df,
    aes(x = !!sym(xvar), y = fit),
    color = "blue", size = 1, inherit.aes = FALSE),
    geom_ribbon(
      data = pred_df,
      aes(x = !!sym(xvar), ymin = lower, ymax = upper),
      fill = "lightblue", alpha = 0.3, inherit.aes = FALSE))}

Lay_hatch_order_model <- glmmTMB(
  Hatch_order ~ Known_lay_order + (1|Year),
  family = gaussian, data = laying_data)

pred_data <- data.frame(Known_lay_order = seq(1, 11, length.out = 200))

preds <- predict(Lay_hatch_order_model,
                 newdata = pred_data,
                 type = "response",
                 re.form = NA,
                 se.fit = TRUE)

pred_lay_hatch <- pred_data %>%
  mutate(fit   = preds$fit,
         se    = preds$se.fit,
         lower = fit - 1.96 * se,
         upper = fit + 1.96 * se)

pred_data2 <- data.frame(Known_lay_order = 1:11) %>%
  mutate(lay_order_c1 = pmin(Known_lay_order, 6),
         lay_order_c2 = pmax(Known_lay_order - 6, 0))

preds2 <- predict(model_time_in_nest,
                  newdata = pred_data2,
                  type = "response",
                  re.form = NA,
                  se.fit = TRUE)

incubation_df <- pred_data2 %>%
  mutate(fit   = preds2$fit,
         se    = preds2$se.fit,
         lower = fit - 1.96 * se,
         upper = fit + 1.96 * se)

var_labels <- c(
  Known_lay_order = "Lay order",
  Hatch_order = "Hatch order",
  Incubation_period = "Time in nest (days)")

axis_breaks <- list(
  Known_lay_order = seq(1, 11, 1),
  Hatch_order = seq(1, 10, 1),
  Incubation_period = seq(21, 33, 2))

axis_limits <- list(
  Known_lay_order = c(1, 11),
  Hatch_order = c(1, 10),
  Incubation_period = c(20.5, 33))


make_plot <- function(data, xvar, yvar, label = NULL) {
  ggplot(data, aes_string(xvar, yvar)) +
    geom_jitter(width = 0.1, height = 0.1, alpha = 0.8) +
    labs(x = var_labels[xvar], y = var_labels[yvar], title = label) +
    theme_classic(base_size = 9) +
    scale_x_continuous(breaks = seq(0, 11, 0.5),
                       labels = ifelse(seq(0, 11, by = 0.5) %% 1 == 0,
                                       seq(0, 11, by = 0.5), "")) +
    scale_y_continuous(breaks = axis_breaks[[yvar]],
                       limits = axis_limits[[yvar]],
                       labels = scales::label_number(accuracy = 1, pad = TRUE)) +
    theme(plot.title = element_text(hjust = 0, face = "bold", color = "black", size = 12, margin = margin(b=1)),
          axis.text  = element_text(size = 8, color = "black"),
          axis.title = element_text(size = 10), 
          panel.spacing = unit(0.01, "lines"), 
          plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5))}

pA24 <- make_plot(laying_data, "Known_lay_order", "Hatch_order") +
  ggtitle("A") +
  add_prediction(pred_lay_hatch, "Known_lay_order") +
  scale_y_continuous(breaks = seq(0, 10, 0.5),
                     labels = ifelse(seq(0, 10, by = 0.5) %% 1 == 0,
                                     seq(0, 10, by = 0.5), ""))

pB24 <- make_plot(laying_data, "Known_lay_order", "Incubation_period") +
  only_left + only_bottom + ggtitle("B") +
  geom_line(data = incubation_df,
    aes(x = Known_lay_order, y = fit_smooth),
    colour = "blue", size = 1,
    inherit.aes = FALSE) +
  geom_ribbon(
    data = incubation_df,
    aes(x = Known_lay_order, ymin = lower, ymax = upper),
    fill = "lightblue", alpha = 0.3,
    inherit.aes = FALSE)

figure24 <- plot_grid(pA24, pB24, ncol = 1, align = "v")
figure24


# Model on egg size by known lay order.
egg_size_model <- glmmTMB(Volume_cm_cubed ~ Known_lay_order + (1|Nest_ID) + (1|Year), 
                          family = gaussian, data = laying_data)
egg_size_model_noyear <- glmmTMB(Volume_cm_cubed ~ Known_lay_order + (1|Nest_ID), 
                                 family = gaussian, data = laying_data)
egg_size_model_nonest <- glmmTMB(Volume_cm_cubed ~ Known_lay_order + (1|Year), 
                                 family = gaussian, data = laying_data)
diagnostics(egg_size_model)
summary(egg_size_model)
print(CI_95(egg_size_model))
anova(egg_size_model, egg_size_model_nonest)
anova(egg_size_model, egg_size_model_noyear)


# Adding chick data 
chick_data <- read_excel("Final_Compiled_Chick_Data.xlsx") %>%
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>%
  filter(!(Mass > 35 | Tarsus > 40 | Shield_to_tip > 24)) %>%
  filter(!(Tarsus & Year < 2019))
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

# Filtering out data for nests excluding nests that hatched less than 2 eggs
chick_data_survival <- chick_data_survival %>%
  filter(Hatched_eggs > 1)

# Model with hatch order as sole predictor variable on survival.
model_hatch_order <- glmmTMB(as.numeric(Survived) ~ as.numeric(Hatch_order) +
                  (1|Nest_ID) + (1|Year), family = binomial, data = chick_data_survival)
model_hatch_order_nonest <- glmmTMB(as.numeric(Survived) ~ as.numeric(Hatch_Day) 
                                  + (1|Year), family = binomial, data = chick_data_survival)
model_hatch_order_noyear <- glmmTMB(as.numeric(Survived) ~ as.numeric(Hatch_Day) + 
                                    (1|Nest_ID) , family = binomial, data = chick_data_survival)
diagnostics(model_hatch_order)
summary(model_hatch_order)
print(CI_95(model_hatch_order))
anova(model_order_day, model_hatch_order_noyear)
anova(model_order_day, model_hatch_order_nonest)

# Model with day of hatching in the sequence as sole predictor variable on survival.
model_hatch_day <- glmmTMB(as.numeric(Survived) ~ as.numeric(Hatch_Day) + 
                    (1|Nest_ID) + (1|Year), family = binomial, data = chick_data_survival)
model_hatch_day_nonest <- glmmTMB(as.numeric(Survived) ~ as.numeric(Hatch_Day) 
                        + (1|Year), family = binomial, data = chick_data_survival)
model_hatch_day_noyear <- glmmTMB(as.numeric(Survived) ~ as.numeric(Hatch_Day) + 
                       (1|Nest_ID) , family = binomial, data = chick_data_survival)
diagnostics(model_hatch_day)
summary(model_hatch_day)
print(CI_95(model_hatch_day))
anova(model_hatch_day, model_hatch_day_noyear)
anova(model_hatch_day, model_hatch_day_nonest)

# Model with interaction between hatch order and hatching spread
model_order_spread_int <- glmmTMB(as.numeric(Survived) ~ as.numeric(Hatch_order)*Hatch_spread +
                     (1|Nest_ID) + (1|Year), family = binomial, data = chick_data_survival)
model_order_spread_int_nonest <- glmmTMB(as.numeric(Survived) ~ as.numeric(Hatch_order)*Hatch_spread
                                    + (1|Year), family = binomial, data = chick_data_survival)
model_order_spread_int_noyear <- glmmTMB(as.numeric(Survived) ~ as.numeric(Hatch_order)*Hatch_spread +
                                 (1|Nest_ID), family = binomial, data = chick_data_survival)
diagnostics(model_order_spread_int)
summary(model_order_spread_int)
print(CI_95(model_order_spread_int))
anova(model_order_spread_int, model_order_spread_int_noyear)
anova(model_order_spread_int, model_order_spread_int_nonest)

# Model with interaction between hatch order and clutch size
model_order_clutch_int <- glmmTMB(as.numeric(Survived) ~ as.numeric(Hatch_order)*Clutch_size +
                        (1|Nest_ID) + (1|Year), family = binomial, data = chick_data_survival)
model_order_clutch_int_noyear <- glmmTMB(as.numeric(Survived) ~ as.numeric(Hatch_order)*Clutch_size +
                          (1|Nest_ID), family = binomial, data = chick_data_survival)
model_order_clutch_int_nonest <- glmmTMB(as.numeric(Survived) ~ as.numeric(Hatch_order)*Clutch_size +
                                 (1|Year), family = binomial, data = chick_data_survival)
diagnostics(model_order_clutch_int)
summary(model_order_clutch_int)
print(CI_95(model_order_clutch_int))
anova(model_order_clutch_int, model_order_clutch_int_nonest)
anova(model_order_clutch_int, model_order_clutch_int_year)


# Model of mass by hatch order
model_mass <- glmmTMB(Mass ~ Hatch_order + (1|Nest_ID) + (1|Year), 
                      family = gaussian, data = chick_data)
model_mass_noyear <- glmmTMB(Mass ~ Hatch_order + (1|Nest_ID), 
                      family = gaussian, data = chick_data)
model_mass_nonest <- glmmTMB(Mass ~ Hatch_order + (1|Year), 
                      family = gaussian, data = chick_data)
diagnostics(model_mass)
summary(model_mass)
print(CI_95(model_mass))
anova(model_mass, model_mass_nonest)
anova(model_mass, model_mass_noyear)

# Model of shield to tip length by hatch order
model_StoT <- glmmTMB(Shield_to_tip ~ as.numeric(Hatch_order) + (1|Nest_ID) +
                      (1|Year), family = gaussian, data = chick_data)
model_StoT_nonest <- glmmTMB(Shield_to_tip ~ as.numeric(Hatch_order) +
                        (1|Year), family = gaussian, data = chick_data)
model_StoT_noyear <- glmmTMB(Shield_to_tip ~ as.numeric(Hatch_order) + (1|Nest_ID)
                             , family = gaussian, data = chick_data)
diagnostics(model_StoT)
summary(model_StoT)
print(CI_95(model_StoT))
anova(model_StoT, model_StoT_nonest)
anova(model_StoT, model_StoT_noyear)

# Filtering out data before 2019 because they used a different measurement
tars_data <- chick_data %>%
  filter(!(Year < 2019))

# Model 
model_Tars <- glmmTMB(Tarsus ~ as.numeric(Hatch_order) + (1|Nest_ID) + (1|Year), 
                      family = gaussian, data = tars_data)
model_Tars_noyear <- glmmTMB(Tarsus ~ as.numeric(Hatch_order) + (1|Nest_ID), 
                      family = gaussian, data = tars_data)
model_Tars_nonest <- glmmTMB(Tarsus ~ as.numeric(Hatch_order)+ (1|Year), 
                      family = gaussian, data = tars_data)
diagnostics(model_Tars)
summary(model_Tars)
print(CI_95(model_Tars))
anova(model_Tars, model_Tars_nonest)
anova(model_Tars, model_Tars_noyear)


# Figure 2.5
# Simulating predicted data from morphometrics models

pred_data <- data.frame(Hatch_order = seq(1, 10, length.out = 200))

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

var_labels <- c(Hatch_order = "Hatch order",
                Mass = "Mass (g)",
                Shield_to_tip = "Shield to tip (mm)",
                Tarsus = "Left outer tarsus (mm)")

axis_breaks <- list(Hatch_order = seq(1, 10, 1),
                    Mass = seq(15, 35, 5),
                    Shield_to_tip = seq(18, 23, 1), 
                    Tarsus = seq(20, 32, 2))

axis_limits <- list(Hatch_order = c(1, 10),
                    Mass = c(15, 35),
                    Shield_to_tip  = c(18, 23), 
                    Tarsus = c(20, 32))

make_plot <- function(data, xvar, yvar, label = NULL) {
  ggplot(data, aes_string(xvar, yvar)) +
    geom_jitter(width = 0.075, height = 0.1, alpha = 0.8) +
    labs(x = var_labels[xvar], y = var_labels[yvar], title = label) +
    theme_classic(base_size = 9) +  
    scale_x_continuous(breaks = seq(0, 10, 0.5),
                       labels = ifelse(seq(0, 10, by = 0.5) %% 1 == 0, seq(0, 10, by = 0.5), "")) +
    scale_y_continuous(breaks = axis_breaks[[yvar]], limits = axis_limits[[yvar]]) +
    theme(plot.title = element_text(hjust = 0, face = "bold", color = "black", size = 12, margin = margin(b=1)),
          axis.text  = element_text(size = 8, color = "black"),
          axis.title = element_text(size = 10), 
          panel.spacing = unit(0.01, "lines"), 
          plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5))}

pA25 <- make_plot(chick_data, "Hatch_order", "Mass") + only_left + invisible_x + ggtitle("A") + add_prediction(pred_mass, "Hatch_order")
pB25 <- make_plot(chick_data, "Hatch_order", "Shield_to_tip") + only_left + invisible_x + ggtitle("B") + add_prediction(pred_shield, "Hatch_order")
pC25 <- make_plot(tars_data, "Hatch_order", "Tarsus") + only_left + only_bottom + ggtitle("C") + add_prediction(pred_tarsus, "Hatch_order")

row1 <- plot_grid(pA25, ncol = 1)
row2 <- plot_grid(pB25, ncol = 1)
row3 <- plot_grid(pC25, ncol = 1)

figure25 <- plot_grid(row1, row2, row3, ncol = 1, align = "v") 
figure25


# Model of survival by mass and hatch order interaction
Survived_mass_hatch_order_model <- glmmTMB(as.numeric(Survived) ~ Mass*Hatch_order +
                     (1|Nest_ID) + (1|Year), family = binomial, 
                      data = chick_data_survival)
diagnostics(Survived_mass_hatch_order_model)
summary(Survived_mass_hatch_order_model)
print(CI_95(Survived_mass_hatch_order_model))

Survived_mass_hatch_order_model_nonest <- glmmTMB(as.numeric(Survived) ~ Mass*Hatch_order + (1|Year),
                     data = chick_data_survival, family = binomial)
Survived_mass_hatch_order_model_noyear <- glmmTMB(as.numeric(Survived) ~ Mass*Hatch_order + (1|Nest_ID),
                       data = chick_data_survival, family = binomial)
anova(Survived_mass_hatch_order_model, Survived_mass_hatch_order_model_nonest)
anova(mSurvived_mass_hatch_order_model, Survived_mass_hatch_order_model_noyear)

# Supplemenatary figure of survival by size at hatching. Decided to go with 
# a table instead of figure but here's the code for the plot. 
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

# Supplementary Figure A1
all_data <- read_excel("Nests_masterlist.xlsx",
                       na = c("", "NO_RECORD", "MISSING")) %>%
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>%  
  filter(Exclusion == "GOOD" | Exclusion == "MISSING_HATCH_ORDER" | Exclusion == "FAILED") %>%
  mutate(Hatched_eggs = as.numeric(Hatched_eggs)) %>%
  mutate(Clutch_size = as.numeric(Clutch_size)) %>%
  mutate(Hatch_success = (as.numeric(Hatched_eggs)/as.numeric(Clutch_size)*100)) %>%
  mutate(end_date = if_else(Exclusion == "FAILED", Last_observed, Hatch_begin),
         obs_nesting_period = as.numeric(end_date - Date_found, units = "days"))

# Observation period model
Observation_period_model <- glmmTMB(cbind(Hatched_eggs, Clutch_size - Hatched_eggs) ~ 
                            obs_nesting_period +
                            (1|Year), family = betabinomial,
                            data = all_data)
diagnostics(Observation_period_model)
summary(Observation_period_model)
print(CI_95(Observation_period_model))

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
summarize_data(successful_nests)

lapply(list(summary_experiment_data, summary_successful, 
            summary_successful_HS, summary_successful_MT), print, n = 30)

# T-TESTS analyzing the effectiveness of the experimental design
t.test(Foreign_percentage~Treatment, data = experiment_data) # proportion swapped
t.test(Manipulated_clutch_size~Treatment, data = experiment_data) # Clutch size
t.test(Observed_nesting_period~Treatment, data = successful_nests) # Nesting period
t.test(True_hatch_spread~Treatment, data = successful_nests) # HS
t.test(Transfered_time~Treatment, data = successful_nests_MT) # Transferred time


# Model analyzing number of hatched eggs by clutch size
model_hatch_exp <- glmmTMB(cbind(Hatched, Manipulated_clutch_size - Hatched) ~ 
                     Treatment + Manipulated_clutch_size + Observed_nesting_period
                     + (1|Year), family = betabinomial, data = successful_nests)
model_hatch_exp_noyear <- glmmTMB(cbind(Hatched, Manipulated_clutch_size - Hatched) ~ 
                          Treatment + Manipulated_clutch_size + Observed_nesting_period,
                          family = betabinomial, data = successful_nests)
diagnostics(model_hatch_exp)
summary(model_hatch_exp)
print(CI_95(model_hatch_exp))
anova(model_hatch_exp, model_hatch_exp_noyear)

# Get predictions at a few clutch sizes (and average hatching_spread)
emm_model_hatch_exp <- emmeans(model_hatch_exp, ~ 1, 
               at = list(Manipulated_clutch_size = c(5,7,9,12)),
               Observed_nesting_period = mean(successful_nests$Observed_nesting_period, na.rm=TRUE),
               Treatment = "Control",                # choose values you care about
               type = "response")                    # gives predicted probability (0-1)

df <- as.data.frame(model_hatch_exp)
df$expected_hatched <- df$response * df$Manipulated_clutch_size
df$expected_hatched_lower <- df$lower.CL * df$Manipulated_clutch_size
df$expected_hatched_upper <- df$upper.CL * df$Manipulated_clutch_size
df

# Model analyzing number of survivors by treatment.
successful_nests$Treatment <- factor(successful_nests$Treatment)

model_survive_exp <- glmmTMB(cbind(Survival_60, Hatched - Survival_60) ~ 
                     Treatment + Hatched + Manipulated_clutch_size + (1|Year), 
                   family = betabinomial, data = successful_nests)
model_survive_exp_noyear <- glmmTMB(cbind(Survival_60, Hatched - Survival_60) ~ 
                            Treatment + Hatched + Manipulated_clutch_size, 
                          family = betabinomial, data = successful_nests)
diagnostics(model_survive_exp)
summary(model_survive_exp)
print(CI_95(model_survive_exp))
anova(model_2, model_2_noyear)

# Model predictions of numbers of survivors
Hatched_mod <- 2:10
Clutch_size_mod <- 2:15
Treatment_mod <- c("Synchronous", "Asynchronous")

pred <- emmeans(model_survive_exp, ~ Hatched + Treatment + Manipulated_clutch_size, 
                at=list(Hatched=Hatch_spread_mod, Manipulated_clutch_size = Clutch_size_mod,
                        Treatment_mod = Treatment), 
                type="response", mode="asymptotic")

pred_df <- as.data.frame(pred)
pred_df$Predicted <- pmin(pmax(pred_df$prob * pred_df$Manipulated_clutch_size, 0), pred_df$Manipulated_clutch_size)
pred_df$Lower <- pmin(pmax(pred_df$asymp.LCL * pred_df$Manipulated_clutch_size, 0), pred_df$Manipulated_clutch_size)
pred_df$Upper <- pmin(pmax(pred_df$asymp.UCL * pred_df$Manipulated_clutch_size, 0), pred_df$Manipulated_clutch_size)
pred_df

predicted <- list(geom_ribbon(data = pred_df,
                                aes(x = Manipulated_clutch_size, ymin = Lower, ymax = Upper),
                                fill = "lightblue", alpha = 0.3, inherit.aes = FALSE),
                    geom_line(data = pred_df,
                              aes(x = Manipulated_clutch_size, y = Predicted),
                              colour = "blue", size = 1, inherit.aes = FALSE))

# Figure 3.3
pred_survived_clutch_size <- ggplot(pred_df, aes(x = Manipulated_clutch_size, 
                                                    y = Predicted, 
                                                    color = Treatment, 
                                                    fill = Treatment)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL), alpha = 0.3, color = NA) +
  labs(x = "Number of Hatched eggs",
       y = "Predicted number of survivors") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0, face = "bold", color = "black", size = 12, margin = margin(b=1))) +
  scale_y_continuous(lim = c(0,8), breaks = seq(0, 10, by = 1), expand = expansion(mult = c(0, .05))) +
  scale_x_continuous(lim = c(1,15), breaks = seq(1, 15, by = 2), expand = expansion(mult = c(0, .05))) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Asynchronous" = "blue", "Synchronous" = "firebrick1")) +
  scale_fill_manual(values = c("Asynchronous" = "blue", "Synchronous" = "firebrick1")) +
  theme(plot.title = element_text(hjust = 0, face = "bold", color = "black", size = 12, margin = margin(b=1)))

# Grid for hatched eggs
Hatched_mod <- seq(2, 10, by = 0.1)
Clutch_mean <- mean(successful_nests$Manipulated_clutch_size, na.rm = TRUE)

pred_hatched <- emmeans(model_2,
                        ~ Hatched + Treatment,
                        at = list(Hatched = Hatched_mod,
                                  Manipulated_clutch_size = Clutch_mean,
                                  Treatment = Treatment_mod),
                        type = "response", mode = "asymptotic") %>%
  as.data.frame()

pred_hatched$Predicted <- pmin(pmax(pred_hatched$prob * pred_hatched$Hatched, 0),
                               pred_hatched$Hatched)
pred_hatched$Lower <- pmin(pmax(pred_hatched$asymp.LCL * pred_hatched$Hatched, 0),
                           pred_hatched$Hatched)
pred_hatched$Upper <- pmin(pmax(pred_hatched$asymp.UCL * pred_hatched$Hatched, 0),
                           pred_hatched$Hatched)

p_hatched <- ggplot(pred_hatched, aes(x = Hatched, y = Predicted,
                                      color = Treatment, fill = Treatment)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, color = NA) +
  labs(x = "Number of hatched eggs",
       y = "Predicted number of survivors") +
  theme_classic(base_size = 12) +
  theme(axis.title.y = element_text(colour = "white"), axis.ticks.y = element_line(),
                       axis.text.y  = element_text()) +
  scale_y_continuous(lim = c(0,9), breaks = seq(0, 9, by = 1), expand = expansion(mult = c(0, .05))) +
  scale_x_continuous(lim = c(2,10), breaks = seq(2, 10, by = 2), expand = expansion(mult = c(0, .05))) +
  scale_color_manual(values = c("Asynchronous" = "blue", "Synchronous" = "firebrick1")) +
  theme(plot.title = element_text(hjust = 0, face = "bold", color = "black", size = 12, margin = margin(b=1))) +
  scale_fill_manual(values = c("Asynchronous" = "blue", "Synchronous" = "firebrick1"))

p_clutch <- p_clutch + ggtitle("A")
p_hatched <- p_hatched + ggtitle("B")

final<- (p_clutch + p_hatched) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom") +
  labs(color = NULL, fill = NULL)
final


# Predicted probability on response scale (0-1)
pred <- predict(model_2, newdata = newdata, type = "response", re.form = NA)

# Add predictions to data frame
newdata$pred_prob <- pred
newdata$expected_survivors <- newdata$pred_prob * newdata$Hatched

set.seed(2025)

# 1. Extract fixed effects and variance-covariance matrix
beta <- fixef(model_2)$cond
V <- vcov(model_2)$cond

# 2. Number of simulations
nsim <- 1000
sim_beta <- mvrnorm(nsim, mu = beta, Sigma = V)

# 3. Create model matrix for newdata
newdata$Treatment <- factor(newdata$Treatment, 
                            levels = levels(successful_nests$Treatment))

X <- model.matrix(~ Treatment + Hatched + Manipulated_clutch_size, data = newdata)

# 4. Calculate predicted probabilities for each simulation
eta <- sim_beta %*% t(X)           # linear predictor
prob <- 1 / (1 + exp(-eta))        # inverse logit
expected_survivors_sim <- prob * newdata$Hatched

newdata$expected_survived_lower <- pmax(0, apply(expected_survivors_sim, 2, quantile, 0.025))
newdata$expected_survived_upper <- pmin(newdata$Hatched, 
                                        apply(expected_survivors_sim, 2, quantile, 0.975))
newdata

plotdata_summary <- newdata %>%
  group_by(Hatched, Treatment) %>%
  summarize(
    expected_survivors = mean(expected_survivors),
    lower = mean(expected_survived_lower),
    upper = mean(expected_survived_upper))

invisible_y <- theme(axis.title.y = element_text(color = "white"), axis.ticks.y = element_blank(),
                     axis.text.y  = element_text(color = "white"))

pred_survived_hatch <- ggplot(plotdata_summary, aes(x = Hatched, 
                             y = expected_survivors, 
                             color = Treatment, 
                             fill = Treatment)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  labs(x = "Number of Hatched eggs",
       y = "Predicted number of survivors") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0, face = "bold", color = "black", size = 12, margin = margin(b=1))) +
  scale_y_continuous(lim = c(0,8), breaks = seq(0, 10, by = 1), expand = expansion(mult = c(0, .05))) +
  scale_x_continuous(lim = c(1,11), breaks = seq(1, 11, by = 2), expand = expansion(mult = c(0, .05))) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Asynchronous" = "blue", "Synchronous" = "firebrick1")) +
  scale_fill_manual(values = c("Asynchronous" = "blue", "Synchronous" = "firebrick1"))


# Survival to 60 days pivoted longer format for the figure.
survival_data_long <- successful_nests %>%
  mutate(NonSurvivors = Hatched - Survival_60) %>% 
  pivot_longer(cols = c(Survival_60, NonSurvivors), 
               names_to = "Status", values_to = "Count") %>%
  mutate(Survival_60 = ifelse(Status == "Survival_60", 1, 0)) %>%
  uncount(Count) %>%
  mutate(Treatment_survival = paste(Treatment, Survival_60, sep = "_")) 

# Additional plot of survival combining all nests together.
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

# Figure 3.2
# Boxplot of number of hatched eggs and survivors by treatment.

pA32 <- ggplot(successful_nests, aes(x = Treatment, y = Hatched, colour = Treatment)) +
  geom_boxplot(aes(fill = Treatment), outlier.shape = NA, width = 0.5, colour = "black", alpha = 0.3) +
  geom_beeswarm(cex = 4.0, size = 2) +
  theme_classic() +
  scale_colour_manual(values = c("Asynchronous" = "blue", "Synchronous" = "firebrick1")) +
  scale_fill_manual(values = c("Asynchronous" = "blue", "Synchronous" = "firebrick1")) +
  labs(x = "", y = "Hatched eggs") +
  scale_y_continuous(limits = c(0, 7.5), breaks = seq(0, 8, by = 1)) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.text.x = element_text(size = 12, colour = "black")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0, face = "bold", color = "black", size = 12, margin = margin(b=1)))

# Plot 2: Survivors
pB32 <- ggplot(successful_nests, aes(x = Treatment, y = Survival_60, colour = Treatment)) +
  geom_boxplot(aes(fill = Treatment), outlier.shape = NA, width = 0.5, colour = "black", alpha = 0.3) +
  geom_beeswarm(cex = 4.0, size = 2) +
  theme_classic() +
  scale_colour_manual(values = c("Asynchronous" = "blue", "Synchronous" = "firebrick1")) +
  scale_fill_manual(values = c("Asynchronous" = "blue", "Synchronous" = "firebrick1")) +
  labs(x = "", y = "Survivors") +
  scale_y_continuous(limits = c(0, 7.5), breaks = seq(0, 8, by = 1)) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.text.x = element_text(size = 12, colour = "black")) +
  theme(plot.title = element_text(hjust = 0, face = "bold", color = "black", size = 12, margin = margin(b=1))) +
  theme(legend.position = "none")

pA32<- pA32 + ggtitle("A")
pB32<- pB32 + ggtitle("B") 

combined_plot <- pA32 / pB32 
combined_plot

# Data on if they were fostered or non-fostered. 
swapping_data <- read_excel("Egg_swapping.xlsx") %>%
  mutate(Nest_ID = paste(Year, Nest_ID, sep = "_")) %>%
  filter(Number_nest_transfers == 1) %>%
  rename(Hatched_yes = Hatched) %>%
  left_join(successful_nests %>% select(Nest_ID, Manipulated_clutch_size, Hatched),
            by = "Nest_ID")

# Figure 3.4 of hatch order by egg status and nest treatment.
ggplot(swapping_data, aes(x = Status, y = Hatch_order, colour = Treatment)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, colour = "black", fill = NA) +
  geom_beeswarm(aes(colour = Treatment), cex = 4.0, size = 2, width = 0.3) +
  theme_classic() +
  labs(y = "Hatch Order", x = "") +
  theme(axis.text.x = element_text(size = 12, colour = "black")) +
  scale_y_continuous(breaks = seq(0, 8, by = 0.5),
                     labels = ifelse(seq(0, 8, by = 0.5) %% 
                                       1 == 0, seq(0, 8, by = 0.5), "")) +
  scale_x_discrete(expand = expansion(add = c(0.6, 0.6)), 
    labels = c("Original" = "Non-fostered", "Swapped" = "Fostered")) +
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

### Hatching success by swapping status
Swap_hatch_model <- glmmTMB(Hatched_yes ~ Status*Treatment + (1|Nest_ID) +
                              (1|Year), family = binomial, data = swapping_data) 
Swap_hatch_model_nonest <- glmmTMB(Hatched_yes ~ Status*Treatment + (1|Nest_ID) +
                              (1|Year), family = binomial, data = swapping_data) 
Swap_hatch_model_noyear <- glmmTMB(Hatched_yes ~ Status*Treatment + (1|Nest_ID) +
                              (1|Year), family = binomial, data = swapping_data) 
diagnostics(Swap_hatch_model)
summary(Swap_hatch_model)
print(CI_95(Swap_hatch_model))
anova(Swap_hatch_model, Swap_hatch_model_nonest)
anova(Swap_hatch_model, Swap_hatch_model_noyear)

# Survival model by swapping status
Swap_survive_model <- glmmTMB(Survived ~ Status*Treatment + Hatch_order*Treatment + (1|Nest_ID) +
                     + (1|Year), family = binomial, data = swapping_data)
Swap_survive_model_noyear <- glmmTMB(Survived ~ Status*Treatment + Hatch_order*Treatment + (1|Nest_ID) +
                                + (1|Year), family = binomial, data = swapping_data)
Swap_survive_model_nonest <- glmmTMB(Survived ~ Status*Treatment + Hatch_order*Treatment + (1|Nest_ID) +
                                + (1|Year), family = binomial, data = swapping_data)
diagnostics(Swap_survive_model)
summary(Swap_survive_model)
print(CI_95(Swap_survive_model))
anova(Swap_survive_model, Swap_survive_model_noyear)
anova(Swap_survive_model, Swap_survive_model_nonest)

# Hatch order model based on swap status and treatment
Hatch_order_model <- glmmTMB(Hatch_order ~ Status*Treatment + (1|Nest_ID) + (1|Year), 
                     family = gaussian, data = swapping_data)
Hatch_order_model_noyear <- glmmTMB(Hatch_order ~ Status*Treatment + (1|Nest_ID), 
                             family = gaussian, data = swapping_data)
Hatch_order_model_nonest <- glmmTMB(Hatch_order ~ Status*Treatment  + (1|Year), 
                             family = gaussian, data = swapping_data)
summary(Hatch_order_model)
summary(Hatch_order_model_noyear)
summary(Hatch_order_model_nonest)


# Power analysis at different effect sizes 
brood_sizes <- c(2,4,6,8)
betas <- c(0.2,0.4,0.6,0.8)
p0 <- 0.75
alpha <- 0.05
power <- 0.80

# helper
inv_logit <- function(x) plogis(x)
logit_p0 <- qlogis(p0)

out <- data.frame()
for (B in brood_sizes) {
       for (b in betas) {
             p1 <- inv_logit(logit_p0 + b)
             mu0 <- B * p0
             mu1 <- B * p1
             var0 <- B * p0 * (1 - p0)
             var1 <- B * p1 * (1 - p1)
             pooled_sd <- sqrt((var0 + var1)/2)
             d <- (mu1 - mu0) / pooled_sd    # Cohen's d approximation
             # pwr.t.test for required n per group
               if (d > 0) {pw <- pwr.t.test(d = d, sig.level = alpha, power = power,
                           type = "two.sample", alternative = "two.sided")
                     n_per_group <- ceiling(pw$n)} else {
                       n_per_group <- NA}
             sample_sizes_results <- rbind(out, data.frame(Brood = B, beta = b,
                                          p1 = round(p1,4),
                                          mu0 = round(mu0,3),
                                          mu1 = round(mu1,3),
                                          pooled_sd = round(pooled_sd,3),
                                          cohend = round(d,3),
                                          n_per_group = n_per_group))}}
print(sample_sizes_results)
