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
source("scripts/00 - funcs.R") # Load custom helper functions
load("data/CHF_models.RData") # Load saved survival and joint models

library(future) # Parallel processing
library(furrr) # Parallel execution of functions
library(tidyverse) # Data wrangling and visualization
library(JMbayes2) # Joint modeling and time-dependent AUC
library(ggbreak) # Custom breaks for ggplot axes

# ============================================================================ #
# Joint Model vs Cox Model: Discrimination and Calibration Analysis -----------
# ============================================================================ #

# -------------------------------------------------------------------------- #
## Train-Test Split for Model Evaluation --------------------------------------
# -------------------------------------------------------------------------- #
# Generate cross-validation folds from the longitudinal dataset.
# The function `create_train_test` stratifies patients into training and test sets.
eval_df_1 <- create_train_test(long_df, follow_up = tests_time)
eval_df_3 <- create_train_test(long_df_3, follow_up = tests_time)

# Train models on the first cross-validation fold.
models_1 <- fit_models(eval_df_1)
models_3 <- fit_models(eval_df_3)

# -------------------------------------------------------------------------- #
## Bootstrap-Based Model Evaluation -------------------------------------------
# -------------------------------------------------------------------------- #
B <- 200 # Number of bootstrap iterations

# Enable parallel processing for faster computation.
n_cores <- parallel::detectCores() - 1
cl <- parallel::makeCluster(n_cores)
plan(cluster, workers = cl)

# Perform bootstrap resampling in parallel for each dataset.
bootEVAL_results_1 <- future_map_dfr(
  seq_len(B),
  ~ bootstrap_iteration(.x, eval_df_1$test, models_1),
  .options = furrr_options(seed = TRUE)
)
bootEVAL_results_3 <- future_map_dfr(
  seq_len(B),
  ~ bootstrap_iteration(.x, eval_df_3$test, models_3),
  .options = furrr_options(seed = TRUE)
)

# Combine all bootstrap results into a single dataset.
boot_eval_df_1 <- bind_rows(bootEVAL_results_1)
boot_eval_df_3 <- bind_rows(bootEVAL_results_3)

# Disable parallel processing.
plan(sequential)

# -------------------------------------------------------------------------- #
## Summarizing Bootstrap Results ----------------------------------------------
# -------------------------------------------------------------------------- #
# Compute mean, standard deviation, median, and 95% confidence intervals for AUC and Brier scores.
boot_eval_summary <- bind_rows(
  boot_eval_df_1 %>% mutate(cohort = "Cohort 1"),
  boot_eval_df_3 %>% mutate(cohort = "Cohort 3")
) %>%
  mutate(
    auc_diff_36 = auc_joint_36 - auc_cox_36,
    auc_diff_60 = auc_joint_60 - auc_cox_60,
    auc_diff_1_36 = auc_joint_1_36 - auc_cox_1_36,
    auc_diff_1_60 = auc_joint_1_60 - auc_cox_1_60,
    brier_diff_36 = brier_joint_36 - brier_cox_36,
    brier_diff_60 = brier_joint_60 - brier_cox_60,
    brier_diff_1_36 = brier_joint_1_36 - brier_cox_1_36,
    brier_diff_1_60 = brier_joint_1_60 - brier_cox_1_60
  ) %>%
  pivot_longer(cols = -cohort, names_to = "res", values_to = "value") %>%
  group_by(cohort, res) %>%
  summarise(
    mean = mean(value),
    sd = sd(value),
    median = median(value),
    lower_ci = quantile(value, 0.025),
    upper_ci = quantile(value, 0.975)
  ) %>%
  ungroup() %>%
  mutate(
    model_type = case_when(
      str_detect(res, "joint") ~ "Serial 6MWT",
      str_detect(res, "cox") ~ "Single 6MWT",
      str_detect(res, "diff") ~ "Difference"
    ),
    fup_start = ifelse(
      str_detect(res, "_1_"),
      "1-year follow-up",
      "Baseline"
    ),
    fup_time = ifelse(str_detect(res, "36"), "36 (m)", "60 (m)"),
    metric = ifelse(str_detect(res, "auc"), "AUC", "Brier Score")
  )

# Save the results.
save(
  eval_df_1,
  eval_df_3,
  models_1,
  models_3,
  boot_eval_df_1,
  boot_eval_df_3,
  boot_eval_summary,
  file = "data/CHF_models_eval.RData"
)

# -------------------------------------------------------------------------- #
## Visualization of Bootstrap Results -----------------------------------------
# -------------------------------------------------------------------------- #
# Plot ROC-AUC for different models and time points.
plot_boot_df(
  boot_eval_summary,
  metric_filter = "AUC",
  y_lab = "AUC",
  plot_title = "",
  plot_subtitle = "",
  hline_position = 0.5
) -> auc_plot

auc_plot
# Save the AUC plot.
ggsave(
  file.path("export", "boot_auc_graph.tiff"),
  last_plot(),
  width = 30,
  height = 20,
  dpi = 300,
  background = "white",
  units = "cm"
)

# Plot ROC-AUC differencec for different models and time points.
plot_boot_df(
  boot_eval_summary,
  model_filter = "Difference",
  metric_filter = "AUC",
  plot_title = "AUC differences between joint and Cox models",
  hline_position = 0
)

# Save the AUC plot.
ggsave(
  file.path("export", "boot_auc_diff_graph.tiff"),
  last_plot(),
  width = 30,
  height = 20,
  dpi = 300,
  background = "white",
  units = "cm"
)

# Plot Brier scores.
plot_boot_df(
  boot_eval_summary,
  metric_filter = "Brier Score",
  y_lab = "Brier Score",
  ast_y = 0.275,
  plot_title = "",
  plot_subtitle = "",
  hline_position = 0
) -> brier_plot

brier_plot
# Save the Brier score plot.
ggsave(
  file.path("export", "boot_brier_graph.tiff"),
  last_plot(),
  width = 30,
  height = 20,
  dpi = 300,
  background = "white",
  units = "cm"
)

# Plot Brier scores differencec.
plot_boot_df(
  boot_eval_summary,
  model_filter = "Difference",
  metric_filter = "Brier Score",
  plot_title = "Brier score differences between joint and Cox models",
  hline_position = 0
)

# Save the Brier score diff plot.
ggsave(
  file.path("export", "boot_brier_diff_graph.tiff"),
  last_plot(),
  width = 30,
  height = 20,
  dpi = 300,
  background = "white",
  units = "cm"
)

# -------------------------------------------------------------------------- #
## Figure 1 ------------------------------------------------------------------
# -------------------------------------------------------------------------- #
ggarrange(
  auc_plot +
    ggtitle("AUC ↑") +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title.y = element_blank(),
      strip.text.y.right = element_blank(),
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
      panel.spacing = unit(0.01, "cm")
    ) + #,
    rremove("xlab"),
  brier_plot +
    ggtitle("Brier Score ↓") +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title.y = element_blank(),
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
      panel.spacing = unit(0.05, "cm")
    ) + #,
    rremove("xlab"),
  ncol = 2,
  common.legend = TRUE,
  legend = "bottom"
)

ggsave(
  file.path("export", "fig_1.tiff"),
  last_plot(),
  width = 20,
  height = 12,
  dpi = 300,
  background = "white",
  units = "cm"
)
