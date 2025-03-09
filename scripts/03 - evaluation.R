# ============================================================================ #
# Joint Model Evaluation and Bootstrap Analysis
#
# Overview:
#   This script evaluates the discrimination and calibration of joint models
#   using time-dependent ROC curves and calibration metrics. Additionally, it 
#   performs a bootstrap-based evaluation to assess the stability of AUC 
#   estimates across multiple resampling iterations.
#
# Key Components:
#   1. Load required libraries, functions, and pre-trained models.
#   2. Create train-test splits for cross-validation using longitudinal data.
#   3. Fit models on training data across different data subsets.
#   4. Perform bootstrap resampling to estimate AUC and Brier score variability.
#   5. Summarize and visualize the results using time-dependent ROC and calibration plots.
#
# Required Packages:
#   - future, furrr: For parallelized bootstrap execution.
#   - tidyverse: For data manipulation and visualization.
#   - JMbayes2: For joint modeling and time-dependent AUC computations.
# ============================================================================ #

# ============================================================================ #
# Load Functions, Data, and Libraries
# ============================================================================ #
source("scripts/00 - funcs.R")         # Load custom helper functions
load("data/CHF_models.RData")          # Load saved survival and joint models

library(future)                        # Parallel processing
library(furrr)                         # Parallel execution of functions
library(tidyverse)                     # Data wrangling and visualization
library(JMbayes2)                      # Joint modeling and time-dependent AUC

# ============================================================================ #
# Joint Model vs Cox Model: Discrimination and Calibration Analysis -----------
# ============================================================================ #

# -------------------------------------------------------------------------- #
## Train-Test Split for Model Evaluation --------------------------------------
# -------------------------------------------------------------------------- #
# Generate cross-validation folds from the longitudinal dataset.
# The function `create_train_test` stratifies patients into training and test sets.
eval_df_1 <- create_train_test(long_df, follow_up = tests_time)
eval_df_2 <- create_train_test(long_df_2, follow_up = tests_time)
eval_df_3 <- create_train_test(long_df_3, follow_up = tests_time)
# eval_df_4 <- create_train_test(long_df_4, follow_up = tests_time)

# Train models on the first cross-validation fold.
models_1 <- fit_models(eval_df_1)
models_2 <- fit_models(eval_df_2)
models_3 <- fit_models(eval_df_3)
# models_4 <- fit_models(eval_df_4)

# -------------------------------------------------------------------------- #
## Bootstrap-Based Model Evaluation -------------------------------------------
# -------------------------------------------------------------------------- #
B <- 200  # Number of bootstrap iterations

# Enable parallel processing for faster computation.
plan(multisession, workers = parallel::detectCores() - 1)

# Perform bootstrap resampling in parallel for each dataset.
bootEVAL_results_1 <- future_map(1:B, ~ bootstrap_iteration(.x, eval_df_1$test, models_1),
                                 .options = furrr_options(seed = TRUE))
bootEVAL_results_2 <- future_map(1:B, ~ bootstrap_iteration(.x, eval_df_2$test, models_2),
                                 .options = furrr_options(seed = TRUE))
bootEVAL_results_3 <- future_map(1:B, ~ bootstrap_iteration(.x, eval_df_3$test, models_3),
                                 .options = furrr_options(seed = TRUE))
# bootEVAL_results_4 <- future_map(1:B, ~ bootstrap_iteration(.x, eval_df_4$test, models_4),
#                                 .options = furrr_options(seed = TRUE))

# Combine all bootstrap results into a single dataset.
boot_eval_df_1 <- bind_rows(bootEVAL_results_1)
boot_eval_df_2 <- bind_rows(bootEVAL_results_2)
boot_eval_df_3 <- bind_rows(bootEVAL_results_3)
# boot_eval_df_4 <- bind_rows(bootEVAL_results_4)

# Disable parallel processing.
plan(sequential)

# Save the results.
save(boot_eval_df_1, boot_eval_df_2, boot_eval_df_3,
     file = "data/CHF_models_eval.RData")

# -------------------------------------------------------------------------- #
## Summarizing Bootstrap Results ----------------------------------------------
# -------------------------------------------------------------------------- #
# Compute mean, standard deviation, median, and 95% confidence intervals for AUC and Brier scores.
boot_eval_summary <- bind_rows(
  boot_eval_df_1 %>% mutate(cohort = "Cohort 1"),
  boot_eval_df_2 %>% mutate(cohort = "Cohort 2"),
  boot_eval_df_3 %>% mutate(cohort = "Cohort 3")
) %>%
  mutate(
    auc_diff_36      = auc_joint_36 - auc_cox_36,
    auc_diff_60      = auc_joint_60 - auc_cox_60,
    auc_diff_1_36    = auc_joint_1_36 - auc_cox_1_36,
    auc_diff_1_60    = auc_joint_1_60 - auc_cox_1_60,
    brier_diff_36    = brier_joint_36 - brier_cox_36,
    brier_diff_60    = brier_joint_60 - brier_cox_60,
    brier_diff_1_36  = brier_joint_1_36 - brier_cox_1_36,
    brier_diff_1_60  = brier_joint_1_60 - brier_cox_1_60
  ) %>%
  pivot_longer(cols = -cohort, names_to = "res", values_to = "value") %>%
  group_by(cohort, res) %>%
  summarise(
    mean     = mean(value),
    sd       = sd(value),
    median   = median(value),
    lower_ci = quantile(value, 0.025),
    upper_ci = quantile(value, 0.975)
  ) %>%
  ungroup() %>%
  mutate(model_type = case_when(str_detect(res, "joint") ~ "Joint Model",
                                str_detect(res, "cox") ~ "Cox Model",
                                str_detect(res, "diff") ~ "Difference"),
         fup_start = ifelse(str_detect(res, "_1_"), "1-year follow-up", "Baseline"),
         fup_time = ifelse(str_detect(res, "36"), "36 months", "60 months"),
         metric = ifelse(str_detect(res, "auc"), "AUC", "Brier Score"))

# -------------------------------------------------------------------------- #
## Visualization of Bootstrap Results -----------------------------------------
# -------------------------------------------------------------------------- #

# Plot ROC-AUC for different models and time points.
plot_boot_df(boot_eval_summary, 
             metric_filter = "AUC", 
             y_lab = "AUC",
             # plot_title = "Models evaluation Estimates",
             plot_title = "",
             plot_subtitle = "",
             hline_position = 0.5) -> auc_plot

# Save the AUC plot.
ggsave(file.path("export", "boot_auc_graph.jpeg"), last_plot(), 
       width = 30, height = 20, dpi = 300, background = "white", units = "cm")

# Plot ROC-AUC differencec for different models and time points.
plot_boot_df(boot_eval_summary, 
             model_filter = "Difference",
             metric_filter = "AUC", 
             plot_title = "AUC differences between joint and Cox models",
             hline_position = 0)

# Save the AUC plot.
ggsave(file.path("export", "boot_auc_diff_graph.jpeg"), last_plot(), 
       width = 30, height = 20, dpi = 300, background = "white", units = "cm")

# Plot Brier scores.
plot_boot_df(boot_eval_summary, 
             metric_filter = "Brier Score", 
             y_lab = "Brier Score", 
             plot_title = "",
             plot_subtitle = "",
             hline_position = 0) -> brier_plot

# Save the Brier score plot.
ggsave(file.path("export", "boot_brier_graph.jpeg"), last_plot(), 
       width = 30, height = 20, dpi = 300, background = "white", units = "cm")

# Plot Brier scores differencec.
plot_boot_df(boot_eval_summary, 
             model_filter = "Difference",
             metric_filter = "Brier Score", 
             plot_title = "Brier score differences between joint and Cox models",
             hline_position = 0)

# Save the Brier score diff plot.
ggsave(file.path("export", "boot_brier_diff_graph.jpeg"), last_plot(), 
       width = 30, height = 20, dpi = 300, background = "white", units = "cm")


ggarrange(auc_plot +
            rremove("xlab") +
            rremove("x.text"),
          brier_plot,
          nrow = 2,
          common.legend = TRUE,
          legend = "bottom")

ggsave(file.path("export", "fig_1.jpeg"), last_plot(), 
       width = 900, height = 600, dpi = 130, background = "white", units = "px")

# ============================================================================ #
# Dynamic Predictions for Selected Patients -----------------------------------
# ============================================================================ #
# Generate dynamic predictions for different patients at multiple time points.

ggarrange(
  plot_dyn_pred(jointFit, long_df, id = 654, t0 = 0.1), 
  plot_dyn_pred(jointFit, long_df, id = 654, t0 = 6.1),
  plot_dyn_pred(jointFit, long_df, id = 654, t0 = 12.1),
  plot_dyn_pred(jointFit, long_df, id = 654, t0 = 24.1),
  nrow = 4
)

ggarrange(
  plot_dyn_pred(jointFit, long_df, id = 1009, t0 = 2), 
  plot_dyn_pred(jointFit, long_df, id = 1009, t0 = 12),
  plot_dyn_pred(jointFit, long_df, id = 1009, t0 = 24),
  plot_dyn_pred(jointFit, long_df, id = 1009, t0 = 36),
  nrow = 4
)

# Create plots for Patient 282
pat_1 <- ggarrange(
  plot_dyn_pred(jointFit, long_df, id = 282, t0 = 0.1,
                add_labs = FALSE, add_sec_y = TRUE, sec_y_lab = "") + rremove("x.text"),
  plot_dyn_pred(jointFit, long_df, id = 282, t0 = 6.1,
                add_labs = FALSE, add_sec_y = TRUE, sec_y_lab = "") + rremove("x.text"),
  plot_dyn_pred(jointFit, long_df, id = 282, t0 = 12.1,
                add_labs = FALSE, add_sec_y = TRUE, sec_y_lab = "") + rremove("x.text"),
  plot_dyn_pred(jointFit, long_df, id = 282, t0 = 24.1,
                add_labs = FALSE, add_sec_y = TRUE, sec_y_lab = ""),
  nrow = 4, labels = "AUTO"#, align = "hv"
)

# Create plots for Patient 51
pat_2 <- ggarrange(
  plot_dyn_pred(jointFit, long_df, id = 51, t0 = 0.1,
                add_labs = FALSE, sec_y_lab = "") + rremove("x.text"),
  plot_dyn_pred(jointFit, long_df, id = 51, t0 = 6.1,
                add_labs = FALSE, sec_y_lab = "") + rremove("x.text"),
  plot_dyn_pred(jointFit, long_df, id = 51, t0 = 12.1,
                add_labs = FALSE, sec_y_lab = "") + rremove("x.text"),
  plot_dyn_pred(jointFit, long_df, id = 51, t0 = 24.1,
                add_labs = FALSE, sec_y_lab = ""),
  nrow = 4#, align = "hv"
)

# Combine both sets into a single figure
final_plot <- ggarrange(
  pat_1, pat_2,
  ncol = 2,
#  labels = c("Patient 282", "Patient 51")#,
    heights = c(1, 1)
)

# Annotate figure with a main title
annotate_figure(final_plot, 
                top = text_grob("Dynamic Prediction of Mortality Risk Over Time",
                                face = "bold", size = 16),
                left = text_grob("6 minute walk distance (meters)", 
                                 rot = 90, size = 12),
                right = text_grob("Predicted Probability of Mortality", 
                                  rot = 270, size = 12),
                bottom = text_grob("Time (months)", 
                                   rot = 0, vjust = -0.5, size = 12))

ggsave(file.path("export", "fig_2.jpeg"), last_plot(), 
       width = 25, height = 20, dpi = 300, background = "white", units = "cm")
