# ============================================================================ #
# Joint Modeling and Dynamic Predictions Script
#
# This script performs the following:
#   1. Loads required functions and data.
#   2. Pre-processes the data, including identifying missing values,
#      creating datasets for survival (Cox) and longitudinal analyses.
#   3. Develops linear mixed effects models (LME) for 6-minute walk test data.
#   4. Fits Cox proportional hazards models for survival outcomes.
#   5. Constructs joint models combining the LME and survival models.
#   6. Generates dynamic prediction plots for selected subjects.
#
# Libraries used:
#   - JMbayes2: Joint modeling of longitudinal and time-to-event data.
#   - lmtest, MuMIn: Model testing and performance evaluation.
#   - rsample: Resampling tools.
# ============================================================================ #

# ============================================================================ #
# Load Functions, Data, and Libraries
# ============================================================================ #
source("scripts/00 - funcs.R")   # Load custom helper functions (e.g., plot_dyn_pred)
load("data/CHF_data.RData")      # Load the cleaned CHF dataset
library(JMbayes2)                # For joint modeling of longitudinal and survival data
library(lmtest)                  # For likelihood ratio tests and other model testing functions
library(MuMIn)                   # For calculating R-squared for mixed models
library(rsample)                 # For resampling procedures (if needed)

# ============================================================================ #
# Data Pre-processing ---------------------------------------------------------
# ============================================================================ #
# -------------------------------------------------------------------------- #
## Identify IDs with missing variables ####
# -------------------------------------------------------------------------- #
# Count total unique patient IDs in the cleaned dataset.
total_ids <- length(unique(clean_df$id))

# Identify patient IDs with missing values in the 'nyha' variable.
# (Note: you can also use all variables of interest by using if_any() on a vector.)
clean_df %>%
  # Uncomment and adjust the filter below to use additional variables:
  # filter(if_any(c(all_of(seattle_vars), x6mw_dist_meter), ~ is.na(.))) %>%
  filter(if_any(nyha, ~ is.na(.))) %>%   # Check for missing 'nyha' values
  distinct(id) %>%
  pull(id) -> miss_id_vars

# Print the number and percentage of patients with missing variables.
print(paste0("Number of patients with missing variables: ", length(miss_id_vars), 
             " (", round(length(miss_id_vars) / total_ids * 100, 2), " %) out of ",
             total_ids, " total patients"))

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
         nyha = first(nyha),
         ischemic_etiology = any(ischemic_etiology, na.rm = TRUE)) %>%
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
exp(coef(jointFit)$association)

# Joint model on the subset of patients with >1 visit.
jointFit_3 <- jm(CoxFit_3, lme_6min_3, time_var = "time",
                 functional_forms = ~ value(x6mw_dist_meter) + slope(x6mw_dist_meter))
summary(jointFit_3)
exp(coef(jointFit_3)$association)

# (Optional) Save the final joint model for later use.
save(jointFit, file = "data/jointFit.RData")

# ============================================================================ #
# Dynamic Predictions ----------------------------------------------------------
# ============================================================================ #
# Generate dynamic predictions for a selected subject (ID: 152227420) at various landmark times.
ggarrange(
  plot_dyn_pred(jointFit, long_df, id = 152227420, t0 = 2),
  plot_dyn_pred(jointFit, long_df, id = 152227420, t0 = 12),
  plot_dyn_pred(jointFit, long_df, id = 152227420, t0 = 24),
  plot_dyn_pred(jointFit, long_df, id = 152227420, t0 = 36),
  nrow = 4
)

# Generate dynamic predictions for another subject (ID: 274751) at multiple landmark times.
ggarrange(
  plot_dyn_pred(jointFit, long_df, id = 274751, t0 = 12),
  plot_dyn_pred(jointFit, long_df, id = 274751, t0 = 24),
  plot_dyn_pred(jointFit, long_df, id = 274751, t0 = 36),
  nrow = 3
)

# ============================================================================ #
# Save the Final Cleaned Data -------------------------------------------------
# ============================================================================ #
save(cox_df, long_df, CoxFit, lme_6min, jointFit, CoxFit_six,
     cox_df_3, long_df_3, CoxFit_3, lme_6min_3, jointFit_3, CoxFit_six_3,
     long_df_2, long_df_4,
     file = "data/CHF_models.RData")

save(lme_6min_0,lme_6min_n,lme_6min_int_sparse,lme_6min_int,CoxFit_six,
     jointFit_0, jointFit_1, jointFit_1_2, jointFit_3, jointFit_4,
     file = "data/CHF_models_comparisons.RData")
