# ============================================================================ #
# Joint Modeling and Dynamic Predictions of 6-Minute Walk Test and Survival
#
# Overview:
#   This script implements joint modeling of longitudinal 6MWT measurements 
#   and time-to-event survival data to assess dynamic predictions in heart 
#   failure patients. The workflow includes:
#
#   1. Loading required functions, libraries, and data.
#   2. Pre-processing the dataset (handling missing values and defining 
#      survival and longitudinal data structures).
#   3. Developing linear mixed-effects models (LME) for longitudinal 6MWT data.
#   4. Fitting Cox proportional hazards models for survival outcomes.
#   5. Constructing joint models that integrate the LME and survival models.
#
# Dependencies:
#   - JMbayes2: For joint modeling of longitudinal and survival data.
#   - lmtest, MuMIn: Model evaluation tools.
#   - rsample: Resampling functions for validation.
# ============================================================================ #

# ============================================================================ #
# Load Functions, Data, and Libraries
# ============================================================================ #
source("scripts/00 - funcs.R")   # Load custom helper functions
load("data/CHF_data.RData")      # Load the cleaned heart failure dataset

library(JMbayes2)                # Joint modeling for longitudinal & survival data
library(lmtest)                  # Likelihood ratio tests and model comparisons
library(MuMIn)                   # R-squared for mixed models

# ============================================================================ #
# Data Pre-processing ---------------------------------------------------------
# ============================================================================ #

# -------------------------------------------------------------------------- #
## Identify and Handle Missing Values ####
# -------------------------------------------------------------------------- #
# Count the total number of unique patients in the dataset.
total_ids <- length(unique(clean_df$id))

# Identify patients with missing values in the 'nyha' variable.
miss_id_vars <- clean_df %>%
  filter(if_any(nyha, ~ is.na(.))) %>%
  distinct(id) %>%
  pull(id)

# Print the number and percentage of patients with missing values.
print(paste0("Number of patients with missing values: ", length(miss_id_vars), 
             " (", round(length(miss_id_vars) / total_ids * 100, 2), "%) out of ",
             total_ids, " patients"))

# -------------------------------------------------------------------------- #
## Create Models Data Frames ####
# -------------------------------------------------------------------------- #
# Create a dataset for the Cox (survival) model.
# We remove patients with missing variables and then retain only the first visit per patient.
cox_df <- clean_df %>%
  filter(!(id %in% miss_id_vars)) %>%   # Exclude patients with missing variables
  arrange(id, time) %>%                 # Sort by patient ID and time
  group_by(id) %>%
  mutate(n_visits = n()) %>%            # Count the number of visits per patient
  filter(row_number() == 1) %>%         # Keep only the first visit (baseline)
  ungroup()

# Create the longitudinal dataset for the repeated measures analysis.
long_df <- clean_df %>%
  filter(id %in% cox_df$id) %>%         # Only include patients present in the Cox dataset
  arrange(id, time) %>% 
  group_by(id) %>%
  # For each patient, retain the baseline age and NYHA classification,
  # and flag if the patient ever had an ischemic etiology.
  mutate(age = first(age),
         nyha = first(nyha)) %>%
  # Group by key variables and summarize the 6MWT distance (averaging if needed)
  group_by(id, fup_time, mortality_status, time, age, gender, nyha, ischemic_etiology) %>%
  summarise(x6mw_dist_meter = mean(x6mw_dist_meter, na.rm = TRUE),
            .groups = "drop")

# Create a subset of the longitudinal data with only patients having 2 or more visits.
long_df_2 <- long_df %>%
  group_by(id) %>%
  filter(n() >= 2) %>%   # Retain only patients with at least 2 visits
  ungroup()

long_df_3 <- long_df %>%
  group_by(id) %>%
  filter(n() >= 3) %>%   # Retain only patients with at least 3 visits
  ungroup()

long_df_4 <- long_df %>%
  group_by(id) %>%
  filter(n() >= 4) %>%   # Retain only patients with at least 4 visits
  ungroup()

# Update the Cox dataset to include only patients with more than 2 visit.
cox_df_3 <- cox_df %>%
  filter(id %in% long_df_3$id)

# (Optional) Count the number of visits per patient.
long_df %>%
  group_by(id) %>%
  summarise(n = n()) %>%
  group_by(n) %>%
  summarise(t = n())

# Identify patient IDs with high within-subject variability in 6MWT (standard deviation > 25)
long_df %>%
  group_by(id) %>%
  mutate(six_sd = sd(x6mw_dist_meter, na.rm = TRUE)) %>%
  filter(six_sd > 25) %>%
  pull(id) -> high_sd_ids

# ============================================================================ #
# Model Development -----------------------------------------------------------
# ============================================================================ #
# -------------------------------------------------------------------------- #
## Linear Mixed Effects Models (LME) for 6-Minute Walk Test Data ####
# -------------------------------------------------------------------------- #
### Finding the Best Fitting Model -------------------------------------------
# Base model with natural splines on time (3 degrees of freedom) plus age and gender.
lme_6min_0 <- lme(x6mw_dist_meter ~ ns(time, df = 3) + age + gender, 
                  data = long_df, 
                  random = ~ ns(time, df = 3) | id, 
                  control = lmeControl(opt = 'optim'))
r.squaredGLMM(lme_6min_0)

# Extend the model by adding the NYHA variable.
lme_6min_n <- lme(x6mw_dist_meter ~ ns(time, df = 3) + age + gender + nyha, 
                  data = long_df, 
                  random = ~ ns(time, df = 3) | id, 
                  control = lmeControl(opt = 'optim'))
r.squaredGLMM(lme_6min_0)

# Compare the two models using ANOVA.
anova(lme_6min_0, lme_6min_n)

# Test a model with an interaction between time (using splines) and age (sparse interaction).
lme_6min_int_sparse <- lme(x6mw_dist_meter ~ ns(time, df = 3) * (age) + gender + nyha, 
                           data = long_df, 
                           random = ~ ns(time, df = 3) | id, 
                           control = lmeControl(opt = 'optim'))
r.squaredGLMM(lme_6min_int_sparse)

# Test a model with an interaction between time and both gender and age.
lme_6min_int <- lme(x6mw_dist_meter ~ ns(time, df = 3) * (gender + age) + nyha, 
                    data = long_df, 
                    random = ~ ns(time, df = 3) | id, 
                    control = lmeControl(opt = 'optim'))
r.squaredGLMM(lme_6min_int)

# Compare the models to assess improvements in fit.
anova(lme_6min_n, lme_6min_int_sparse, lme_6min_int)

# Final mixed model using the best predictors and interactions.
lme_6min <- lme(x6mw_dist_meter ~ ns(time, df = 3) * (age + gender) + nyha, 
                data = long_df, 
                random = ~ ns(time, df = 3) | id, 
                control = lmeControl(opt = 'optim'))

summary(lme_6min)
r.squaredGLMM(lme_6min)

# Fit the same model on the subset of patients with more than one visit.
lme_6min_3 <- lme(x6mw_dist_meter ~ ns(time, df = 3) * (age + gender) + nyha, 
                  data = long_df_3, 
                  random = ~ ns(time, df = 3) | id, 
                  control = lmeControl(opt = 'optim'))

summary(lme_6min_3)
r.squaredGLMM(lme_6min_3)

# (Optional) Models using only the last measurement have been commented out.
# Uncomment if needed for further analyses.

### Plot: Predicted 6MWT Over Time by Age and Gender ---------------------------
# Create a prediction grid spanning time (0 to 72 months) for various ages.
pred_6mwt_data <- data.frame(
  time   = rep(round(seq(0, 72, length.out = 50)), 4),
  age    = rep(c(rep(40, 50), rep(50, 50), rep(60, 50), rep(70, 50)), 2),
  gender = c(rep("Female", 200), rep("Male", 200)),
  nyha   = "II"   # Fixed NYHA class for prediction purposes
)

# Generate predictions from both the full dataset model and the subset model.
bind_rows(
  pred_6mwt_data %>% 
    mutate(pred = predict(lme_6min, newdata = pred_6mwt_data, level = 0),
           model = "model_1"),
  pred_6mwt_data %>% 
    mutate(pred = predict(lme_6min_3, newdata = pred_6mwt_data, level = 0),
           model = "model_3")
) %>%
  ggplot(aes(x = time, y = pred, colour = age)) +
  geom_point() +
  facet_grid(gender ~ model) +
  labs(title = "Predicted 6MWT Over Time by Age and Gender",
       x = "Time (months)",
       y = "Predicted 6MWT Distance") +
  theme_minimal()

# -------------------------------------------------------------------------- #
## Cox Proportional Hazards Models ####
# -------------------------------------------------------------------------- #
# Fit a baseline Cox model using gender, age, and NYHA.
CoxFit <- coxph(Surv(fup_time, mortality_status) ~ gender + age + nyha, 
                data = cox_df, 
                model = TRUE, x = TRUE, y = TRUE)

CoxFit_3 <- coxph(Surv(fup_time, mortality_status) ~ gender + age + nyha, 
                data = cox_df_3, 
                model = TRUE, x = TRUE, y = TRUE)

# Fit a Cox model that additionally includes the 6MWT distance.
CoxFit_six <- coxph(Surv(fup_time, mortality_status) ~ gender + age + nyha + x6mw_dist_meter, 
                    data = cox_df, 
                    model = TRUE, x = TRUE, y = TRUE)

# Fit a Cox model that additionally includes the 6MWT distance.
CoxFit_six_3 <- coxph(Surv(fup_time, mortality_status) ~ gender + age + nyha + x6mw_dist_meter, 
                    data = cox_df_3, 
                    model = TRUE, x = TRUE, y = TRUE)

# -------------------------------------------------------------------------- #
## Joint Models  ####
# -------------------------------------------------------------------------- #
# Joint model 0: Combine the Cox model (CoxFit) with the longitudinal mixed model (lme_6min)
jointFit_0 <- jm(CoxFit, lme_6min, time_var = "time")
summary(jointFit_0)
exp(coef(jointFit_0)$association)

# Joint model 1: Use functional forms including both the current value and slope of x6mw_dist_meter.
jointFit_1 <- jm(CoxFit, lme_6min, time_var = "time",
                 functional_forms = ~ value(x6mw_dist_meter) + slope(x6mw_dist_meter))
summary(jointFit_1)
exp(coef(jointFit_1)$association)

# Joint model 2: Test an interaction between the current value and slope.
jointFit_1_2 <- jm(CoxFit, lme_6min, time_var = "time",
                 functional_forms = ~ value(x6mw_dist_meter) * slope(x6mw_dist_meter))
summary(jointFit_1_2)
exp(coef(jointFit_1_2)$association)

# Joint model 3: Incorporate the 6MWT measure from the Cox model that includes x6mw_dist_meter.
jointFit_3 <- jm(CoxFit_six, lme_6min, time_var = "time",
                 functional_forms = ~ value(x6mw_dist_meter) + slope(x6mw_dist_meter))
summary(jointFit_3)
exp(coef(jointFit_3)$association)

# Joint model 4: Use splines for the time component in the functional form.
form_splines <- ~ value(x6mw_dist_meter) * ns(time, k = c(10, 20, 40), B = c(0, 72)) + slope(x6mw_dist_meter)
jointFit_4 <- jm(CoxFit_six, lme_6min, time_var = "time",
                 functional_forms = form_splines)
summary(jointFit_4)
exp(coef(jointFit_4)$association)

# Compare all joint models to assess their relative performance.
compare_jm(jointFit_0, jointFit_1, jointFit_1_2, jointFit_3, jointFit_4)

# -------------------------------------------------------------------------- #
## Final models  ####
# -------------------------------------------------------------------------- #
# Final joint model (using CoxFit) with value + slope specification.
jointFit <- jm(CoxFit, lme_6min, time_var = "time",
               functional_forms = ~ value(x6mw_dist_meter) + slope(x6mw_dist_meter))
summary(jointFit)
exp(coef(jointFit)$gammas)
exp(coef(jointFit)$association)

# Joint model on the subset of patients with >1 visit.
jointFit_3 <- jm(CoxFit_3, lme_6min_3, time_var = "time",
                 functional_forms = ~ value(x6mw_dist_meter) + slope(x6mw_dist_meter))
summary(jointFit_3)
exp(coef(jointFit_3)$gammas)
exp(coef(jointFit_3)$association)

# (Optional) Save the final joint model for later use.
save(jointFit, file = "data/jointFit.RData")

# ============================================================================ #
# Full Model Evaluation --------------------------------------------------------
# ============================================================================ #
# -------------------------------------------------------------------------- #
## Follow-up Time Selection ####
# -------------------------------------------------------------------------- #
# Visualize the distribution of follow-up times to help choose the initial
# follow-up time (fu time) for model evaluation.
long_df %>%
  ggplot(aes(x = time)) +
  geom_histogram(binwidth = 1, fill = "gray", color = "black") +
  geom_vline(xintercept = 12, color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = 24, color = "blue", linetype = "dashed", linewidth = 1) +
  labs(title = "Distribution of Follow-up Times",
       x = "Time (months)",
       y = "Frequency") +
  theme_minimal()

# -------------------------------------------------------------------------- #
## Discrimination Analysis (ROC) ####
# -------------------------------------------------------------------------- #
fup_times <- c(24, 48)
# Generate time-dependent ROC data for different joint models and data subsets.
# ROC data using the primary joint model (jointFit) and the full longitudinal data.
roc_data <- create_roc_data(
  model = jointFit,
  long_data = long_df,
  Tstart = tests_time,
  follow_up_times = fup_times
)

# ROC data using the alternative joint model (jointFit_2) and patients with >1 visit.
roc_data_3 <- create_roc_data(
  model = jointFit_3,
  long_data = long_df_3,
  Tstart = tests_time,
  follow_up_times = fup_times
)

roc_data_cox <- create_roc_data(
  model = CoxFit,
  long_data = long_df,
  Tstart = tests_time,
  follow_up_times = fup_times
)

# ROC data using the alternative joint model (jointFit_2) and patients with >1 visit.
roc_data_cox_3 <- create_roc_data(
  model = CoxFit_3,
  long_data = long_df_3,
  Tstart = tests_time,
  follow_up_times = fup_times
)

# -------------------------------------------------------------------------- #
## Calibration Analysis ####
# -------------------------------------------------------------------------- #
# Compute calibration metrics for different joint model specifications and data subsets.
cal_data <- calc_cal_metrics(
  model = jointFit,
  newdata = long_df,
  Tstart = tests_time,
  follow_up_times = fup_times
)

cal_data_3 <- calc_cal_metrics(
  model = jointFit_3,
  newdata = long_df_3,
  Tstart = tests_time,
  follow_up_times = fup_times
)

cal_data_cox <- calc_cal_metrics(
  model = CoxFit,
  newdata = long_df,
  Tstart = tests_time,
  follow_up_times = fup_times
)

cal_data_cox_3 <- calc_cal_metrics(
  model = CoxFit_3,
  newdata = long_df_3,
  Tstart = tests_time,
  follow_up_times = fup_times
)

# -------------------------------------------------------------------------- #
## Plotting ROC and Calibration Curves ####
# -------------------------------------------------------------------------- #
# Arrange and display the ROC and calibration plots.
roc_cal_plot <- ggarrange(
  plot_tvROC(roc_data),
  plot_tvROC(roc_data_cox),
  plot_cal(cal_data),
  plot_cal(cal_data_cox),
  nrow = 2, ncol = 2
)

roc_cal_3_plot <- ggarrange(
  plot_tvROC(roc_data_3),
  plot_tvROC(roc_data_cox_3),
  plot_cal(cal_data_3),
  plot_cal(cal_data_cox_3),
  nrow = 2, ncol = 2
)

# Display the plot
print(roc_cal_plot)
print(roc_cal_3_plot)

# ============================================================================ #
# Save the Final Cleaned Data -------------------------------------------------
# ============================================================================ #
save(cox_df, long_df, CoxFit, lme_6min, jointFit, CoxFit_six,
     cox_df_3, long_df_3, CoxFit_3, lme_6min_3, jointFit_3, CoxFit_six_3,
     long_df_2, long_df_4,
     roc_data, roc_data_3, roc_data_cox, roc_data_cox_3,
     cal_data, cal_data_3, cal_data_cox, cal_data_cox_3,
     file = "data/CHF_models.RData")

save(lme_6min_0,lme_6min_n,lme_6min_int_sparse,lme_6min_int,CoxFit_six,
     jointFit_0, jointFit_1, jointFit_1_2, jointFit_3, jointFit_4,
     file = "data/CHF_models_comparisons.RData")
