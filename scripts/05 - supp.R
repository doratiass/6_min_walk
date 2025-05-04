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
missing_ids_full <- clean_df %>%
  filter(if_any(all_of(small_seattle_vars), is.na)) %>%
  distinct(id) %>%
  pull(id)

# Print the number and percentage of patients with missing values.
print(paste0(
  "Number of patients with missing values: ",
  length(missing_ids_full),
  " (",
  round(length(missing_ids_full) / total_ids * 100, 2),
  "%) out of ",
  total_ids,
  " patients"
))

# -------------------------------------------------------------------------- #
## Create Models Data Frames ####
# -------------------------------------------------------------------------- #
# Create the longitudinal dataset for the repeated measures analysis.
long_df_full <- clean_df %>%
  filter(!(id %in% missing_ids_full)) %>% # Exclude patients with missing variables
  arrange(id, time) %>%
  group_by(id) %>%
  # For each patient, retain the baseline age and NYHA classification,
  # and flag if the patient ever had an ischemic etiology.
  mutate(
    age = first(age),
    nyha = first(nyha),
    sbp = first(sbp),
    bmi = first(bmi),
    na = first(na),
    hgb = first(hgb),
    ischemic_etiology = first(ischemic_etiology),
    fusid = first(fusid),
    ,
    k_spare = first(k_spare),
    ,
    allopurinol = first(allopurinol),
    statin = first(statin),
    ace_arb = first(ace_arb),
    beta_blocker = first(beta_blocker),
  ) %>%
  # Group by key variables and summarize the 6MWT distance (averaging if needed)
  group_by(
    id,
    fup_time,
    mortality_status,
    time,
    !!!syms(small_seattle_vars)
  ) %>%
  summarise(
    x6mw_dist_meter = mean(x6mw_dist_meter, na.rm = TRUE),
    .groups = "drop"
  )

long_df_full_3 <- long_df_full %>%
  group_by(id) %>%
  filter(n() >= 3) %>% # Retain only patients with at least 3 visits
  ungroup()

# Create a dataset for the Cox (survival) model.
# We remove patients with missing variables and then retain only the first visit per patient.
cox_df_full <- long_df_full %>% # Exclude patients with missing variables
  arrange(id, time) %>% # Sort by patient ID and time
  group_by(id) %>%
  mutate(n_visits = n()) %>% # Count the number of visits per patient
  filter(row_number() == 1) %>% # Keep only the first visit (baseline)
  ungroup()

# Update the Cox dataset to include only patients with more than 2 visit.
cox_df_full_3 <- cox_df_full %>%
  filter(id %in% long_df_full_3$id)

# ============================================================================ #
# Model Development -----------------------------------------------------------
# ============================================================================ #
# -------------------------------------------------------------------------- #
## Linear Mixed Effects Models (LME) for 6-Minute Walk Test Data ####
# -------------------------------------------------------------------------- #
# Final mixed model using the best predictors and interactions.
lme_6min_full <- lme(
  x6mw_dist_meter ~ ns(time, df = 3) * (age + gender) + nyha,
  data = long_df_full,
  random = ~ ns(time, df = 3) | id,
  control = lmeControl(opt = 'optim')
)

# Fit the same model on the subset of patients with more than one visit.
lme_6min_full_3 <- lme(
  x6mw_dist_meter ~ ns(time, df = 3) * (age + gender) + nyha,
  data = long_df_full_3,
  random = ~ ns(time, df = 3) | id,
  control = lmeControl(opt = 'optim')
)

# -------------------------------------------------------------------------- #
## Cox Proportional Hazards Models ####
# -------------------------------------------------------------------------- #
# Fit a baseline Cox model using gender, age, and NYHA.
CoxFit_small <- coxph(
  Surv(fup_time, mortality_status) ~ gender + age + nyha,
  data = cox_df_full,
  model = TRUE,
  x = TRUE,
  y = TRUE
)

CoxFit_small_3 <- coxph(
  Surv(fup_time, mortality_status) ~ gender + age + nyha,
  data = cox_df_full_3,
  model = TRUE,
  x = TRUE,
  y = TRUE
)

# Fit a full baseline Cox model using gender, age, and NYHA.
CoxFit_full <- coxph(
  as.formula(paste(
    "Surv(fup_time, mortality_status) ~",
    paste(small_seattle_vars, collapse = " + ")
  )),
  data = cox_df_full,
  model = TRUE,
  x = TRUE,
  y = TRUE
)

CoxFit_full_3 <- coxph(
  as.formula(paste(
    "Surv(fup_time, mortality_status) ~",
    paste(small_seattle_vars, collapse = " + ")
  )),
  data = cox_df_full_3,
  model = TRUE,
  x = TRUE,
  y = TRUE
)

# -------------------------------------------------------------------------- #
## Joint Models  ####
# -------------------------------------------------------------------------- #
# Final joint model (using CoxFit) with value + slope specification.
jointFit_small <- jm(
  CoxFit_small,
  lme_6min_full,
  time_var = "time",
  functional_forms = ~ value(x6mw_dist_meter) + slope(x6mw_dist_meter)
)

# Joint model on the subset of patients with >1 visit.
jointFit_small_3 <- jm(
  CoxFit_small_3,
  lme_6min_full_3,
  time_var = "time",
  functional_forms = ~ value(x6mw_dist_meter) + slope(x6mw_dist_meter)
)

# Final joint model (using CoxFit) with value + slope specification.
jointFit_full <- jm(
  CoxFit_full,
  lme_6min_full,
  time_var = "time",
  functional_forms = ~ value(x6mw_dist_meter) + slope(x6mw_dist_meter)
)

# Joint model on the subset of patients with >1 visit.
jointFit_full_3 <- jm(
  CoxFit_full_3,
  lme_6min_full_3,
  time_var = "time",
  functional_forms = ~ value(x6mw_dist_meter) + slope(x6mw_dist_meter)
)

# -------------------------------------------------------------------------- #
## Table 2 supp ####
# -------------------------------------------------------------------------- #
tbl_small <- jointfit_tbl(jointFit_small, jointFit_small_3)

doc <- read_docx() %>%
  body_add_flextable(tbl_small) %>%
  body_add_par("") # Add space after table

# Save the document
print(doc, target = "export/supp_tbl_small.docx")

tbl_full <- jointfit_tbl(jointFit_full, jointFit_full_3)

doc <- read_docx() %>%
  body_add_flextable(tbl_full) %>%
  body_add_par("") # Add space after table

# Save the document
print(doc, target = "export/supp_tbl_full.docx")
# ============================================================================ #
# Joint Model small vs full: Discrimination and Calibration Analysis -----------
# ============================================================================ #
# -------------------------------------------------------------------------- #
## Train-Test Split for Model Evaluation --------------------------------------
# -------------------------------------------------------------------------- #
# Generate cross-validation folds from the longitudinal dataset.
# The function `create_train_test` stratifies patients into training and test sets.
eval_df_full_1 <- create_train_test(long_df_full, follow_up = tests_time)
eval_df_full_3 <- create_train_test(long_df_full_3, follow_up = tests_time)

# Train models on the first cross-validation fold.
models_full_1 <- fit_models(eval_df_full_1, full = TRUE)
models_full_3 <- fit_models(eval_df_full_3, full = TRUE)

# -------------------------------------------------------------------------- #
## Bootstrap-Based Model Evaluation -------------------------------------------
# -------------------------------------------------------------------------- #
B <- 200 # Number of bootstrap iterations
options(future.globals.maxSize = 891289600)

# Enable parallel processing for faster computation.
plan(multisession, workers = parallel::detectCores() - 1)

# Perform bootstrap resampling in parallel for each dataset.
bootEVAL_results_full_1 <- future_map(
  1:B,
  ~ bootstrap_iteration(.x, eval_df_full_1$test, models_full_1, full = TRUE),
  .options = furrr_options(seed = TRUE)
)
bootEVAL_results_full_3 <- future_map(
  1:B,
  ~ bootstrap_iteration(.x, eval_df_full_3$test, models_full_3, full = TRUE),
  .options = furrr_options(seed = TRUE)
)

# Combine all bootstrap results into a single dataset.
boot_eval_full_df_1 <- bind_rows(bootEVAL_results_full_1)
boot_eval_full_df_3 <- bind_rows(bootEVAL_results_full_3)

# Disable parallel processing.
plan(sequential)

# -------------------------------------------------------------------------- #
## Summarizing Bootstrap Results ----------------------------------------------
# -------------------------------------------------------------------------- #
# Compute mean, standard deviation, median, and 95% confidence intervals for AUC and Brier scores.
boot_eval_full_summary <- bind_rows(
  boot_eval_full_df_1 %>% mutate(cohort = "Cohort 1"),
  boot_eval_full_df_3 %>% mutate(cohort = "Cohort 3")
) %>%
  mutate(
    auc_diff_36 = auc_joint_36 - auc_joint_full_36,
    auc_diff_60 = auc_joint_60 - auc_joint_full_60,
    auc_diff_1_36 = auc_joint_1_36 - auc_joint_full_1_36,
    auc_diff_1_60 = auc_joint_1_60 - auc_joint_full_1_60,
    brier_diff_36 = brier_joint_36 - brier_joint_full_36,
    brier_diff_60 = brier_joint_60 - brier_joint_full_60,
    brier_diff_1_36 = brier_joint_1_36 - brier_joint_full_1_36,
    brier_diff_1_60 = brier_joint_1_60 - brier_joint_full_1_60
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
      str_detect(res, "full") ~ "Full Joint Model",
      str_detect(res, "joint") ~ "Joint Model",
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

# -------------------------------------------------------------------------- #
## Visualization of Bootstrap Results -----------------------------------------
# -------------------------------------------------------------------------- #
# Plot ROC-AUC for different models and time points.
plot_boot_df(
  boot_eval_full_summary,
  model_filter = c("Full Joint Model", "Joint Model"),
  metric_filter = "AUC",
  y_lab = "AUC",
  plot_title = "",
  plot_subtitle = "",
  hline_position = 0.5
) -> auc_plot_full

# Plot Brier scores.
plot_boot_df(
  boot_eval_full_summary,
  model_filter = c("Full Joint Model", "Joint Model"),
  metric_filter = "Brier Score",
  y_lab = "Brier Score",
  ast_y = 0.275,
  plot_title = "",
  plot_subtitle = "",
  hline_position = 0
) -> brier_plot_full

ggarrange(
  auc_plot_full +
    ggtitle("AUC") +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title.y = element_blank(),
      strip.text.y.right = element_blank(),
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
      panel.spacing = unit(0.01, "cm")
    ),
  #      rremove("xlab"),
  brier_plot_full +
    ggtitle("Brier Score") +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title.y = element_blank(),
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
      panel.spacing = unit(0.05, "cm")
    ),
  ncol = 2,
  common.legend = TRUE,
  legend = "bottom"
)

ggsave(
  file.path("export", "supp_fig_1.jpeg"),
  last_plot(),
  width = 30,
  height = 20,
  dpi = 300,
  background = "white",
  units = "cm"
)

save(
  jointFit_small,
  jointFit_small_3,
  jointFit_full,
  jointFit_full_3,
  eval_df_full_1,
  eval_df_full_3,
  models_full_1,
  models_full_3,
  bootEVAL_results_full_1,
  bootEVAL_results_full_3,
  file = "data/supp_data.RData"
)
