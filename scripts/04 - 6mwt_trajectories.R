# ============================================================================ #
# Prediction of 6MWT Change and Survival Analysis
#
# Overview:
#   This script predicts the 6-minute walk test (6MWT) distance over time 
#   using a linear mixed-effects model and evaluates the prognostic value of 
#   the predicted change on mortality using Kaplan-Meier survival analysis.
#
# Key Components:
#   1. Load required libraries and data.
#   2. Predict 6MWT trajectories using mixed-effects models.
#   3. Group patients into quartiles based on predicted change.
#   4. Plot predicted 6MWT values over time by quartile.
#   5. Evaluate the association between predicted change and survival.
#
# Required Packages:
#   - survival, survminer: For survival analysis.
#   - dplyr, tidyr, ggplot2: For data manipulation and visualization.
# ============================================================================ #

# ============================================================================ #
# Load Libraries and Data ---------------------------------------------------- #
# ============================================================================ #
source("scripts/00 - funcs.R")            # Load custom helper functions
load("data/CHF_models.RData")             # Load saved survival and joint models

library(survival)                         # Cox models and survival analysis
library(survminer)                        # Kaplan-Meier plotting
library(dplyr)                            # Data manipulation
library(tidyr)                            # Data tidying
library(ggplot2)                          # Data visualization

# ============================================================================ #
# Create quartile groups ------------------------------------------------------ 
# ============================================================================ #
# -------------------------------------------------------------------------- #
## Create Prediction Grid for Each Patient ------------------------------------
# -------------------------------------------------------------------------- #
# Create a grid for predicted 6MWT over 12 months for each patient
time_id <- long_df_3 %>%
  filter(time <= tests_time, fup_time >= tests_time) %>%
  distinct(id) %>%
  expand(id, time = 0:12)

# -------------------------------------------------------------------------- #
## Predict 6MWT Trajectories -------------------------------------------------- 
# -------------------------------------------------------------------------- #
# Predict 6MWT values using the mixed-effects model
new_df <- long_df_3 %>%
  filter(time <= tests_time, fup_time >= tests_time) %>%
  select(-c(time, x6mw_dist_meter)) %>%                # Remove redundant columns
  distinct(id, .keep_all = TRUE) %>%
  left_join(time_id, by = "id") %>%
  mutate(pred = predict(lme_6min_3, newdata = .)) %>%  # Predict values
  arrange(id, time)

# -------------------------------------------------------------------------- #
## Calculate Predicted Change and Standardize --------------------------------- 
# -------------------------------------------------------------------------- #
# Compute change in predicted 6MWT values over time
quan <- new_df %>%
  group_by(id) %>%
  summarise(change = last(pred) - first(pred)) %>%
  ungroup() %>%
  mutate(
    sd = sd(change),                           # Compute standard deviation
    sds = change / sd,                         # Standardize the change
    quartile = factor(case_when(
      sds < -0.5 ~ "Low",                      # Lower quartile threshold
      sds > 0.5 ~ "High",                      # Upper quartile threshold
      TRUE ~ "Medium"                          # Middle quartile
    ), levels = c("Low", "Medium", "High"))
  )

quan %>%
  group_by(quartile) %>%
  summarise(
    mean_change = mean(change),
    sd_change = sd(change),
    median_change = median(change),
    lower_ci = quantile(change, 0.025),
    upper_ci = quantile(change, 0.975)
  )
# -------------------------------------------------------------------------- #
## Plot Predicted 6MWT Values by Quartile ------------------------------------- 
# -------------------------------------------------------------------------- #
new_df %>%
  left_join(quan %>% select(id, quartile), by = "id") %>%
  ggplot(aes(x = time, y = pred, color = quartile)) +
  geom_line(aes(group = id), alpha = 0.2) +              # Individual trajectories
  geom_smooth(method = "loess", se = FALSE) +             # Smoothed trend line
  theme_minimal() +
  labs(
    title = "Predicted 6-Minute Walk Distance Over Time",
    x = "Time (months)",
    y = "6-Minute Walk Distance (meters)",
    color = "Quartile Group"
  ) +
  theme(legend.position = "bottom")

# ============================================================================ #
# Survival Analysis -----------------------------------------------------------
# ============================================================================ #
# -------------------------------------------------------------------------- #
## Prepare Data for Survival Analysis ----------------------------------------- 
# -------------------------------------------------------------------------- #
# Create survival data frame with updated follow-up time
surv_df <- long_df_3 %>%
  filter(time <= tests_time, fup_time >= tests_time) %>%
  group_by(id) %>%
  arrange(id, time) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(new_fup = fup_time - tests_time) %>%
  left_join(quan %>% select(id, quartile), by = "id")

# Create survival object using follow-up time and mortality status
surv_object <- Surv(time = surv_df$new_fup, event = surv_df$mortality_status)

# -------------------------------------------------------------------------- #
## Kaplan-Meier Curves -------------------------------------------------------- 
# -------------------------------------------------------------------------- #
# Fit Kaplan-Meier model by quartile group
fit <- survfit(surv_object ~ quartile, data = surv_df)

# Plot Kaplan-Meier Curves
ggsurvplot(
  fit,
  data = surv_df,
  conf.int = TRUE,               # Show confidence intervals
  pval = TRUE,                   # Display p-value from log-rank test
  risk.table = TRUE,             # Add risk table below the plot
  risk.table.col = "strata",     # Color the risk table by group
  legend.title = "Quartile Group",
  legend.labs = c("Low", "Medium", "High"),
  xlab = "Time (months)",
  ylab = "Survival Probability",
  palette = "Dark2",             # Use a clean color palette
  ggtheme = theme_minimal()
)