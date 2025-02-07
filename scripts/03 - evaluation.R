# ============================================================================ #
# Joint Model Evaluation and Bootstrap Analysis
#
# This script evaluates the discrimination and calibration of joint models
# using time-dependent ROC and calibration curves. In addition, it performs a
# bootstrap evaluation (via cross-validation folds) to assess the stability of
# AUC measures.
#
# Required packages:
#   - future, furrr: For parallel processing
#   - tidyverse: For data manipulation and plotting
#   - JMbayes2: For joint modeling and tvAUC computations
# ============================================================================ #

# ============================================================================ #
# Load Functions, Data, and Libraries
# ============================================================================ #
source("scripts/00 - funcs.R")         # Load custom helper functions 
load("data/CHF_models.RData")          # Load the models
library(future)                        # For parallel processing
library(furrr)                         # For parallel map functions (future_map)
library(tidyverse)                     # For data manipulation and visualization
library(JMbayes2)                      # For joint modeling and tvAUC computations

# ============================================================================ #
# Joint models vs Coxph comparison ---------------------------------------------
# ============================================================================ #
# -------------------------------------------------------------------------- #
## Train-Test creation ####
# -------------------------------------------------------------------------- #
# Create cross-validation folds from the longitudinal data, grouping by subject ID.
CVdats <- create_folds(long_df, V = 5, id_var = "id", seed = 229)
data_train <- CVdats$training[[1]]
data_test <- CVdats$testing[[1]]
# --- Fit Models on the Training Data (First Fold) ---

# Fit the linear mixed model on the training fold.
lme_6min_train <- lme(
  x6mw_dist_meter ~ ns(time, df = 3) * (age + gender) + nyha, 
  data = data_train, 
  random = ~ ns(time, df = 3) | id, 
  control = lmeControl(opt = 'optim')
)

# Fit the Cox proportional hazards model on the training fold.
CoxFit_train <- coxph(
  Surv(fup_time, mortality_status) ~ gender + age + nyha, 
  data = data_train, 
  model = TRUE, x = TRUE, y = TRUE
)

# Combine the models into a joint model.
jointFit_train <- jm(
  CoxFit_train, lme_6min_train, time_var = "time",
  functional_forms = list("x6mw_dist_meter" = ~ value(x6mw_dist_meter) + slope(x6mw_dist_meter))
)

# create coxph model with 6min walk distance for comparison
CoxFit_six_train <- coxph(
  Surv(fup_time, mortality_status) ~ gender + age + nyha + x6mw_dist_meter, 
  data = data_train %>%
    arrange(id, time) %>%                 # Sort by patient ID and time
    group_by(id) %>%
    filter(row_number() == 1) %>%         # Keep only the first visit (baseline)
    ungroup(), 
  model = TRUE, x = TRUE, y = TRUE
)

# -------------------------------------------------------------------------- #
## Bootstrap Evaluation ####
# -------------------------------------------------------------------------- #
# This section performs bootstrap evaluation to assess the stability of AUC measures.
unique_ids <- unique(data_test$id)
B <- 100  # Number of bootstrap iterations (increase as needed)

# Set up parallel processing using multisession.
plan(multisession, workers = parallel::detectCores() - 1)

# Run bootstrap iterations in parallel.
bootROC_results <- future_map(1:B, ~ bootstrap_iteration(.x), 
                              .options = furrr_options(seed = TRUE))

# Combine all bootstrap results into a single data frame.
bootROC_results_df <- bind_rows(bootROC_results)
print(bootROC_results_df)

# Switch back to sequential processing.
plan(sequential)

# Summarize the bootstrap results: compute mean, standard deviation, median, and 95% CI.
bootROC_summary <- bootROC_results_df %>%
  pivot_longer(cols = everything(), names_to = "AUC", values_to = "value") %>%
  group_by(AUC) %>%
  summarise(
    mean     = mean(value),
    sd       = sd(value),
    median   = median(value),
    lower_ci = quantile(value, 0.025),
    upper_ci = quantile(value, 0.975)
  )
print(bootROC_summary)

bootROC_summary %>%
  ggplot(aes(x = AUC, y = mean, ymin = lower_ci, ymax = upper_ci)) +
  geom_pointrange() +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(title = "Bootstrap AUC estimates",
       x = "AUC",
       y = "Value") +
  theme_minimal()


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

# Set the follow-up time (in months) for later filtering
follow_up_time <- 12

# Create datasets with only the last measurement up to the defined follow-up time.
long_df_last <- long_df %>% 
  filter(time <= follow_up_time) %>%
  group_by(id) %>%
  mutate(x6mw_dist_meter = last(x6mw_dist_meter)) %>%
  ungroup()

long_df_2_last <- long_df_2 %>% 
  filter(time <= follow_up_time) %>%
  group_by(id) %>%
  mutate(x6mw_dist_meter = last(x6mw_dist_meter)) %>%
  ungroup()

# -------------------------------------------------------------------------- #
## Discrimination Analysis (ROC) ####
# -------------------------------------------------------------------------- #
# Generate time-dependent ROC data for different joint models and data subsets.

# ROC data using the primary joint model (jointFit) and the full longitudinal data.
roc_data <- create_roc_data(
  joint_model = jointFit,
  long_data = long_df,
  Tstart = follow_up_time,
  follow_up_times = c(12, 36, 60)
)

# ROC data using the primary joint model with only the last measurements.
roc_data_last <- create_roc_data(
  joint_model = jointFit,
  long_data = long_df_last,
  Tstart = follow_up_time,
  follow_up_times = c(12, 36, 60)
)

# ROC data using the alternative joint model (jointFit_2) and patients with >1 visit.
roc_data_2 <- create_roc_data(
  joint_model = jointFit_2,
  long_data = long_df_2,
  Tstart = follow_up_time,
  follow_up_times = c(12, 36, 60)
)

# ROC data for the alternative joint model using only the last measurements.
roc_data_2_last <- create_roc_data(
  joint_model = jointFit_2,
  long_data = long_df_2_last,
  Tstart = follow_up_time,
  follow_up_times = c(12, 36, 60)
)

# -------------------------------------------------------------------------- #
## Calibration Analysis ####
# -------------------------------------------------------------------------- #
# Compute calibration metrics for different joint model specifications and data subsets.
cal_data <- calc_cal_metrics(
  joint_model = jointFit,
  newdata = long_df,
  Tstart = follow_up_time,
  follow_up_times = c(12, 36, 60)
)

cal_data_last <- calc_cal_metrics(
  joint_model = jointFit,
  newdata = long_df_last,
  Tstart = follow_up_time,
  follow_up_times = c(12, 36, 60)
)

cal_data_2 <- calc_cal_metrics(
  joint_model = jointFit_2,
  newdata = long_df_2,
  Tstart = follow_up_time,
  follow_up_times = c(12, 36, 60)
)

cal_data_2_last <- calc_cal_metrics(
  joint_model = jointFit_2,
  newdata = long_df_2_last,
  Tstart = follow_up_time,
  follow_up_times = c(12, 36, 60)
)

# -------------------------------------------------------------------------- #
## Plotting ROC and Calibration Curves ####
# -------------------------------------------------------------------------- #
# Arrange and display the ROC and calibration plots.
# First set of plots using the primary joint model (jointFit)
plot1 <- ggarrange(
  plot_tvROC(roc_data),
  plot_tvROC(roc_data_last),
  plot_cal(cal_data),
  plot_cal(cal_data_last),
  nrow = 2, ncol = 2
)

# Second set of plots using the alternative joint model (jointFit_2)
plot2 <- ggarrange(
  plot_tvROC(roc_data_2),
  plot_tvROC(roc_data_2_last),
  plot_cal(cal_data_2),
  plot_cal(cal_data_2_last),
  nrow = 2, ncol = 2
)

# Display the plots
print(plot1)
print(plot2)
