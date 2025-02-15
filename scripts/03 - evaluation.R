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

eval_df_1 <- create_train_test(long_df, follow_up = tests_time)
eval_df_2 <- create_train_test(long_df_2, follow_up = tests_time)
eval_df_3 <- create_train_test(long_df_3, follow_up = tests_time)
eval_df_4 <- create_train_test(long_df_4, follow_up = tests_time)

# --- Fit Models on the Training Data (First Fold) ---
models_1 <- fit_models(eval_df_1)
models_2 <- fit_models(eval_df_2)
models_3 <- fit_models(eval_df_3)
models_4 <- fit_models(eval_df_4)

# -------------------------------------------------------------------------- #
## Bootstrap Evaluation ####
# -------------------------------------------------------------------------- #
B <- 100  # Number of bootstrap iterations (increase as needed)

# Set up parallel processing using multisession.
plan(multisession, workers = parallel::detectCores() - 1)

# Run bootstrap iterations in parallel.
bootROC_results_1 <- future_map(1:B, ~ bootstrap_iteration(.x, eval_df_1$test, models_1),
                                .options = furrr_options(seed = TRUE))
bootROC_results_2 <- future_map(1:B, ~ bootstrap_iteration(.x, eval_df_2$test, models_2),
                                .options = furrr_options(seed = TRUE))
bootROC_results_3 <- future_map(1:B, ~ bootstrap_iteration(.x, eval_df_3$test, models_3),
                                .options = furrr_options(seed = TRUE))
bootROC_results_4 <- future_map(1:B, ~ bootstrap_iteration(.x, eval_df_4$test, models_4),
                                .options = furrr_options(seed = TRUE))

# Combine all bootstrap results into a single data frame.
bootROC_results_df_1 <- bind_rows(bootROC_results_1)
bootROC_results_df_2 <- bind_rows(bootROC_results_2)
bootROC_results_df_3 <- bind_rows(bootROC_results_3)
bootROC_results_df_4 <- bind_rows(bootROC_results_4)

# Switch back to sequential processing.
plan(sequential)

# Summarize the bootstrap results: compute mean, standard deviation, median, and 95% CI.
bootROC_summary <- bind_rows(
  bootROC_results_df_1 %>% mutate(model = "Model 1"),
  bootROC_results_df_2 %>% mutate(model = "Model 2"),
  bootROC_results_df_3 %>% mutate(model = "Model 3"),
  bootROC_results_df_4 %>% mutate(model = "Model 4")
) %>%
  pivot_longer(cols = -model, names_to = "AUC", values_to = "value") %>%
  group_by(model, AUC) %>%
  summarise(
    mean     = mean(value),
    sd       = sd(value),
    median   = median(value),
    lower_ci = quantile(value, 0.025),
    upper_ci = quantile(value, 0.975)
  )

bootROC_summary %>% 
  ggplot(aes(x = AUC, y = mean, color = AUC)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray40") +
 # scale_color_brewer(palette = "Set1") +
  labs(
    title    = "Bootstrap AUC Estimates",
    subtitle = "Mean with 95% Confidence Intervals",
    x        = "AUC Metric",
    y        = "Value",
    color    = "AUC"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "bottom",
    # panel.grid.major = element_line(color = "gray85"),
    # panel.grid.minor = element_blank(),
    # strip.background = element_rect(fill = "gray90", color = "gray70"),
    # plot.title       = element_text(face = "bold"),
    # plot.subtitle    = element_text(face = "italic"),
    axis.text.x = element_blank()
  ) +
  facet_wrap(~ model)#, scales = "free_y")

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
  model = models_1$joint_model,
  long_data = eval_df_1$test,
  Tstart = tests_time,
  follow_up_times = fup_times
)

# ROC data using the alternative joint model (jointFit_2) and patients with >1 visit.
roc_data_3 <- create_roc_data(
  model = models_3$joint_model,
  long_data = eval_df_3$test,
  Tstart = tests_time,
  follow_up_times = fup_times
)

roc_data_cox <- create_roc_data(
  model = models_1$cox_model,
  long_data = eval_df_1$test,
  Tstart = tests_time,
  follow_up_times = fup_times
)

# ROC data using the alternative joint model (jointFit_2) and patients with >1 visit.
roc_data_cox_3 <- create_roc_data(
  model = models_3$cox_model,
  long_data = eval_df_3$test,
  Tstart = tests_time,
  follow_up_times = fup_times
)

# -------------------------------------------------------------------------- #
## Calibration Analysis ####
# -------------------------------------------------------------------------- #
# Compute calibration metrics for different joint model specifications and data subsets.
cal_data <- calc_cal_metrics(
  model = models_1$joint_model,
  newdata = eval_df_1$test,
  Tstart = tests_time,
  follow_up_times = fup_times
)

cal_data_3 <- calc_cal_metrics(
  model = models_3$joint_model,
  newdata = eval_df_3$test,
  Tstart = tests_time,
  follow_up_times = fup_times
)

cal_data_cox <- calc_cal_metrics(
  model = models_1$cox_model,
  newdata = eval_df_1$test,
  Tstart = tests_time,
  follow_up_times = fup_times
)

cal_data_cox_3 <- calc_cal_metrics(
  model = models_3$cox_model,
  newdata = eval_df_3$test,
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
