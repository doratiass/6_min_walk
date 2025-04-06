# ============================================================================ #
# Joint Modeling of 6-Minute Walk Test and Survival
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
#   - lmtest, MuMIn: Tools for model comparison and evaluation.
#   - rsample: For resampling and validation.
# ============================================================================ #

# ============================================================================ #
# Load Functions, Data, and Libraries
# ============================================================================ #
source("scripts/00 - funcs.R") # Load custom helper functions
load("data/CHF_data.RData") # Load the cleaned heart failure dataset

library(JMbayes2) # Joint modeling of longitudinal & survival data
library(lmtest) # Likelihood ratio tests and model comparison
library(MuMIn) # R-squared for mixed models
library(flextable) # Create and format tables
library(officer) # Word document generation
library(future) # Parallel processing framework
library(furrr) # Parallel execution of functions
library(ggbreak) # Custom breaks for ggplot axes

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
print(paste0(
  "Number of patients with missing values: ",
  length(miss_id_vars),
  " (",
  round(length(miss_id_vars) / total_ids * 100, 2),
  "%) out of ",
  total_ids,
  " patients"
))

# -------------------------------------------------------------------------- #
## Create Models Data Frames ####
# -------------------------------------------------------------------------- #
# Create the longitudinal dataset for the repeated measures analysis.
long_df <- clean_df %>%
  filter(!(id %in% miss_id_vars)) %>% # Exclude patients with missing variables
  arrange(id, time) %>%
  group_by(id) %>%
  # For each patient, retain the baseline age and NYHA classification,
  # and flag if the patient ever had an ischemic etiology.
  mutate(age = first(age), nyha = first(nyha)) %>%
  # Group by key variables and summarize the 6MWT distance (averaging if needed)
  group_by(id, fup_time, mortality_status, time, age, gender, nyha) %>%
  summarise(
    x6mw_dist_meter = mean(x6mw_dist_meter, na.rm = TRUE),
    .groups = "drop"
  )

long_df_3 <- long_df %>%
  group_by(id) %>%
  filter(n() >= 3) %>% # Retain only patients with at least 3 visits
  ungroup()

# Create a dataset for the Cox (survival) model.
# We remove patients with missing variables and then retain only the first visit per patient.
cox_df <- long_df %>% # Only include patients present in the Cox dataset
  arrange(id, time) %>% # Sort by patient ID and time
  group_by(id) %>%
  mutate(n_visits = n()) %>% # Count the number of visits per patient
  filter(row_number() == 1) %>% # Keep only the first visit (baseline)
  ungroup()

cox_df_3 <- cox_df %>%
  filter(id %in% long_df_3$id)

# -------------------------------------------------------------------------- #
## Generate Baseline Summary Table (Table 1) ####
# -------------------------------------------------------------------------- #
# bind_rows(
#   clean_df %>% filter(id %in% long_df$id) %>% mutate(group = "Cohort 1"),
#   clean_df %>% filter(id %in% long_df_3$id) %>% mutate(group = "Cohort 3")
# ) %>%
clean_df %>%
  filter(id %in% long_df$id) %>%
  mutate(group = ifelse(id %in% long_df_3$id, "Cohort 1+3", "Cohort 1")) %>%
  filter(time == 0) %>%
  distinct(id, group, .keep_all = TRUE) %>%
  group_by(id, group) %>%
  # Compute the number of 6MWT tests per subject and the overall change from first to last measurement
  mutate(
    n_6mw = length(unique(x6mw_dist_meter)),
    x6mw_change = x6mw_dist_meter[which.max(time)] -
      x6mw_dist_meter[which.min(time)],
  ) %>%
  ungroup() %>%
  # Remove identifiers and time variable from the summary table
  select(-c(id, time, x6mw_change)) %>%
  # Rename variables using the custom label_get function for presentation
  rename_all(function(x) sapply(x, label_get, USE.NAMES = FALSE)) %>%
  # Create a summary table stratified by 6MWT distance quantiles using gtsummary
  tbl_summary(
    by = group,
    type = list(label_get("n_6mw") ~ "continuous"),
    statistic = list(all_continuous() ~ "{mean} ({sd})"),
    value = list(label_get("gender") ~ "Male"),
    missing = "no"
  ) %>%
  add_overall() %>%
  add_n(statistic = "{N_miss} ({p_miss})") %>%
  add_difference(everything() ~ "smd") %>%
  modify_column_hide(conf.low) %>%
  modify_header(estimate = "**SMD**", n = "**Missing**")

# Save the summary table as an HTML file
gt::gtsave(as_gt(tbl_1), file = "export/tbl_1.html")
# ============================================================================ #
# Model Development -----------------------------------------------------------
# ============================================================================ #
# -------------------------------------------------------------------------- #
## Linear Mixed Effects Models (LME) for 6-Minute Walk Test Data ####
# -------------------------------------------------------------------------- #
### Finding the Best Fitting Model -------------------------------------------
# Base model with natural splines on time (3 degrees of freedom) plus age and gender.
lme_6min_0 <- lme(
  x6mw_dist_meter ~ ns(time, df = 3) + age + gender,
  data = long_df,
  random = ~ ns(time, df = 3) | id,
  control = lmeControl(opt = 'optim')
)
r.squaredGLMM(lme_6min_0)

# Extend the model by adding the NYHA variable.
lme_6min_n <- lme(
  x6mw_dist_meter ~ ns(time, df = 3) + age + gender + nyha,
  data = long_df,
  random = ~ ns(time, df = 3) | id,
  control = lmeControl(opt = 'optim')
)
r.squaredGLMM(lme_6min_0)

# Compare the two models using ANOVA.
anova(lme_6min_0, lme_6min_n)

# Test a model with an interaction between time (using splines) and age (sparse interaction).
lme_6min_int_sparse <- lme(
  x6mw_dist_meter ~ ns(time, df = 3) * (age) + gender + nyha,
  data = long_df,
  random = ~ ns(time, df = 3) | id,
  control = lmeControl(opt = 'optim')
)
r.squaredGLMM(lme_6min_int_sparse)

# Test a model with an interaction between time and both gender and age.
lme_6min_int <- lme(
  x6mw_dist_meter ~ ns(time, df = 3) * (gender + age) + nyha,
  data = long_df,
  random = ~ ns(time, df = 3) | id,
  control = lmeControl(opt = 'optim')
)
r.squaredGLMM(lme_6min_int)

# Compare the models to assess improvements in fit.
anova(lme_6min_n, lme_6min_int_sparse, lme_6min_int)

# Final mixed model using the best predictors and interactions.
lme_6min <- lme(
  x6mw_dist_meter ~ ns(time, df = 3) * (age + gender) + nyha,
  data = long_df,
  random = ~ ns(time, df = 3) | id,
  control = lmeControl(opt = 'optim')
)

summary(lme_6min)
r.squaredGLMM(lme_6min)

# Fit the same model on the subset of patients with more than one visit.
lme_6min_3 <- lme(
  x6mw_dist_meter ~ ns(time, df = 3) * (age + gender) + nyha,
  data = long_df_3,
  random = ~ ns(time, df = 3) | id,
  control = lmeControl(opt = 'optim')
)

summary(lme_6min_3)
r.squaredGLMM(lme_6min_3)

# (Optional) Models using only the last measurement have been commented out.
# Uncomment if needed for further analyses.

### Plot: Predicted 6MWT Over Time by Age and Gender ---------------------------
# Create a prediction grid spanning time (0 to 72 months) for various ages.
pred_6mwt_data <- data.frame(
  time = rep(round(seq(0, 72, length.out = 50)), 4),
  age = rep(c(rep(40, 50), rep(50, 50), rep(60, 50), rep(70, 50)), 2),
  gender = c(rep("Female", 200), rep("Male", 200)),
  nyha = "II" # Fixed NYHA class for prediction purposes
)

# Generate predictions from both the full dataset model and the subset model.
bind_rows(
  pred_6mwt_data %>%
    mutate(
      pred = predict(lme_6min, newdata = pred_6mwt_data, level = 0),
      model = "model_1"
    ),
  pred_6mwt_data %>%
    mutate(
      pred = predict(lme_6min_3, newdata = pred_6mwt_data, level = 0),
      model = "model_3"
    )
) %>%
  ggplot(aes(x = time, y = pred, colour = age)) +
  geom_point() +
  facet_grid(gender ~ model) +
  labs(
    title = "Predicted 6MWT Over Time by Age and Gender",
    x = "Time (months)",
    y = "Predicted 6MWT Distance"
  ) +
  theme_minimal()

# -------------------------------------------------------------------------- #
## Cox Proportional Hazards Models ####
# -------------------------------------------------------------------------- #
# Fit a baseline Cox model using gender, age, and NYHA.
CoxFit <- coxph(
  Surv(fup_time, mortality_status) ~ gender + age + nyha,
  data = cox_df,
  model = TRUE,
  x = TRUE,
  y = TRUE
)

CoxFit_3 <- coxph(
  Surv(fup_time, mortality_status) ~ gender + age + nyha,
  data = cox_df_3,
  model = TRUE,
  x = TRUE,
  y = TRUE
)

# Fit a Cox model that additionally includes the 6MWT distance.
CoxFit_six <- coxph(
  Surv(fup_time, mortality_status) ~ gender + age + nyha + x6mw_dist_meter,
  data = cox_df,
  model = TRUE,
  x = TRUE,
  y = TRUE
)

# Fit a Cox model that additionally includes the 6MWT distance.
CoxFit_six_3 <- coxph(
  Surv(fup_time, mortality_status) ~ gender + age + nyha + x6mw_dist_meter,
  data = cox_df_3,
  model = TRUE,
  x = TRUE,
  y = TRUE
)

# -------------------------------------------------------------------------- #
## Joint Models  ####
# -------------------------------------------------------------------------- #
# Joint model 0: Combine the Cox model (CoxFit) with the longitudinal mixed model (lme_6min)
jointFit_int_0 <- jm(CoxFit, lme_6min, time_var = "time")
summary(jointFit_int_0)
exp(coef(jointFit_int_0)$association)

# Joint model 1: Use functional forms including both the current value and slope of x6mw_dist_meter.
jointFit_int_1 <- jm(
  CoxFit,
  lme_6min,
  time_var = "time",
  functional_forms = ~ value(x6mw_dist_meter) + slope(x6mw_dist_meter)
)
summary(jointFit_int_1)
exp(coef(jointFit_int_1)$association)

# Joint model 2: Test an interaction between the current value and slope.
jointFit_int_2 <- jm(
  CoxFit,
  lme_6min,
  time_var = "time",
  functional_forms = ~ value(x6mw_dist_meter) * slope(x6mw_dist_meter)
)
summary(jointFit_int_2)
exp(coef(jointFit_int_2)$association)

# Joint model 3: Incorporate the 6MWT measure from the Cox model that includes x6mw_dist_meter.
jointFit_int_3 <- jm(
  CoxFit_six,
  lme_6min,
  time_var = "time",
  functional_forms = ~ value(x6mw_dist_meter) + slope(x6mw_dist_meter)
)
summary(jointFit_int_3)
exp(coef(jointFit_int_3)$association)

# Joint model 4: Use splines for the time component in the functional form.
form_splines <- ~ value(x6mw_dist_meter) *
  ns(time, k = c(10, 20, 40), B = c(0, 72)) +
  slope(x6mw_dist_meter)
jointFit_int_4 <- jm(
  CoxFit_six,
  lme_6min,
  time_var = "time",
  functional_forms = form_splines
)
summary(jointFit_int_4)
exp(coef(jointFit_int_4)$association)

# Compare all joint models to assess their relative performance.
compare_jm(
  jointFit_int_0,
  jointFit_int_1,
  jointFit_int_2,
  jointFit_int_3,
  jointFit_int_4
)

# -------------------------------------------------------------------------- #
## Final models  ####
# -------------------------------------------------------------------------- #
# Final joint model (using CoxFit) with value + slope specification.
jointFit <- jm(
  CoxFit,
  lme_6min,
  time_var = "time",
  functional_forms = ~ value(x6mw_dist_meter) + slope(x6mw_dist_meter)
)
summary(jointFit)
exp(coef(jointFit)$gammas)
exp(coef(jointFit)$association)

# Joint model on the subset of patients with >1 visit.
jointFit_3 <- jm(
  CoxFit_3,
  lme_6min_3,
  time_var = "time",
  functional_forms = ~ value(x6mw_dist_meter) + slope(x6mw_dist_meter)
)
summary(jointFit_3)
exp(coef(jointFit_3)$gammas)
exp(coef(jointFit_3)$association)

# (Optional) Save the final joint model for later use.
save(jointFit, file = "data/jointFit.RData")

# -------------------------------------------------------------------------- #
## Table 2 ####
# -------------------------------------------------------------------------- #
tbl_2 <- jointfit_tbl(jointFit, jointFit_3)

doc <- read_docx() %>%
  body_add_flextable(tbl_2) %>%
  body_add_par("") # Add space after table

# Save the document
print(doc, target = "export/tbl_2.docx")

# ============================================================================ #
# Save the Final Cleaned Data -------------------------------------------------
# ============================================================================ #
save(
  cox_df,
  long_df,
  CoxFit,
  lme_6min,
  jointFit,
  CoxFit_six,
  cox_df_3,
  long_df_3,
  CoxFit_3,
  lme_6min_3,
  jointFit_3,
  CoxFit_six_3,
  #  boot_eval_summary,
  file = "data/CHF_models.RData"
)

save(
  lme_6min_0,
  lme_6min_n,
  lme_6min_int_sparse,
  lme_6min_int,
  CoxFit_six,
  jointFit_int_0,
  jointFit_int_1,
  jointFit_int_2,
  jointFit_int_3,
  jointFit_int_4,
  file = "data/CHF_models_comparisons.RData"
)
