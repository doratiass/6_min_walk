# ============================================================================ #
# Script for helper functions for Joint Modeling of Repeated 6-Minute Walk Test
# Measurements and Mortality Prediction
#
# This script includes functions for:
#   - Data cleaning and variable labeling
#   - Evaluating model performance with time-dependent ROC curves
#   - Assessing calibration of joint models
#   - Plotting dynamic predictions for longitudinal and event processes
#
# Libraries used:
#   - tidyverse: Data manipulation and visualization
#   - gtsummary: Summary statistics and table generation
#   - ggpubr: Publication-ready plots
#   - gridExtra: Arranging multiple plots in a grid
#   - flextable: Formatting tables for reporting
#   - survminer: Survival analysis visualizations
# ============================================================================ #

# Load required libraries
library(tidyverse) # For data manipulation and visualization
library(gtsummary) # For summary statistics and table generation
library(ggpubr) # For publication ready plots
library(gridExtra) # For arranging multiple plots
library(flextable) # For creating formatted tables
library(survminer) # For survival analysis plots
library(survival) # For survival analysis
library(JMbayes2) # For joint modeling and tvAUC computations
library(pmcalibration) # For calibration plots
library(survivalROC) # For time-dependent ROC curves
library(compiler) # For compiling R functions for performance
# ============================================================================ #
# Definitions  ----------------------------------------------------------------
# ============================================================================ #
# Define scaling factors and follow-up time limits
age_scale <- 1 # Scaling for age (set to 1 for no scaling)
six_min_scale <- 1 # Scaling for 6-minute walk distance (set to 1 for no scaling)
min_fup_time <- 6 # Minimum follow-up time in months
max_fup_time <- 120 # Maximum follow-up time in months
tests_time <- 12 # Time limit for 6-minute walk test measurements

# ============================================================================ #
# Data Cleaning Functions -----------------------------------------------------
# ============================================================================ #

# The merge_rows function is used to combine multiple rows of data.
# If the input is numeric, it calculates the mean (ignoring NA values);
# otherwise, it returns the first non-missing value.
merge_rows <- function(x) {
  if (is.numeric(x)) {
    return(mean(x, na.rm = TRUE))
  } else {
    return(na.omit(x)[1])
  }
}

# Variables to be included in the Seattle model (or other specific analyses)
seattle_vars <- c(
  "age",
  "gender",
  "ischemic_etiology",
  "ef",
  "sbp",
  "nyha",
  "fusid",
  "k_spare",
  "allopurinol",
  "dm",
  "htn",
  "lipid",
  "bmi",
  "statin",
  "ace_arb_arni",
  "beta_blocker",
  "sglt2",
  "na",
  "hgb",
  "smoke"
)

small_seattle_vars <- c(
  "age",
  "gender",
  "ischemic_etiology",
  "sbp",
  "nyha",
  "ef",
  "fusid",
  "k_spare",
  "allopurinol",
  "bmi",
  "statin",
  "ace_arb_arni",
  "beta_blocker",
  "sglt2",
  "na",
  "hgb",
  "smoke"
)
# ============================================================================ #
# Variable Labeling -----------------------------------------------------------
# ============================================================================ #

# Create a dictionary (tibble) mapping variable names to descriptive labels.
vars_dict <- tibble(
  "index_year" = "Year of index date",
  "age" = "Age",
  "gender" = "Sex (M)",
  "kupa" = "Healthcare provider",
  "mortality_status" = "Mortality status",
  "fup_time" = "Follow-up time (months)",
  "nyha" = "NYHA class",
  "ef" = "Ejection fraction",
  "ischemic_etiology" = "Ischemic etiology",
  "sbp" = "Systolic blood pressure",
  "fusid" = "Furosemide treatment",
  "torsemide" = "Torsemide treatment",
  "bumetanide" = "Bumetanide treatment",
  "metolazone" = "Metolazone treatment",
  "ht" = "Height",
  "wt" = "Weight",
  "statin" = "Statins treatment",
  "ace_inhibitor" = "ACE inhibitors treatment",
  "ace_arb_arni" = "ACEi, ARB or ARNI treatment",
  "beta_blocker" = "Beta blockers treatment",
  "arbs" = "ARBs treatment",
  "diuretic" = "Diuretics treatment",
  "k_spare" = "Mineralocorticoid receptor antagonist treatment",
  "sglt2" = "Sodium-glucose cotransporter-2 inhibitors treatment",
  "arni" = "Angiotensin receptor-neprilysin inhibitors treatment",
  "dm" = "Diabetes mellitus",
  "smoke" = "Smoking status",
  "lipid" = "Hyperlipidemia",
  "htn" = "Hypertension",
  "thiazide" = "Thiazide diuretics treatment",
  "allopurinol" = "Allopurinol treatment",
  "na" = "Sodium",
  "hgb" = "Hemoglobin",
  "ua" = "Uric acid",
  "x6mw_dist_meter" = "6MWT distance (m)",
  "x6mw_dist_meter_diff" = "6MWT distance change (m)",
  "x6mw_dist_meter_diff_percent" = "6MWT distance change (%)",
  "x6mw_dist_meter_derivative" = "6MWT distance derivative",
  "x6mw_dist_meter_diff_from_max" = "6MWT distance change from max",
  "x6mw_dist_meter_diff_percent_from_max" = "6MWT distance change from max (%)",
  "x6mw_dist_meter_derivative_from_max" = "6MWT distance derivative from max",
  "x6mw_hr_pre" = "6MWT HR pre",
  "x6mw_sat_pre" = "6MWT Sat pre",
  "x6mw_hr_post" = "6MWT HR post",
  "x6mw_sat_post" = "6MWT Sat post",
  "x6mw_num" = "6MWT number",
  "n_6mw" = "Number of 6MWTs",
  "x6mw_change" = "6MWT change",
  "bmi" = "BMI"
) %>%
  pivot_longer(
    cols = 1:ncol(.),
    names_to = "var",
    values_to = "name"
  )

# Function to retrieve a human-readable label given a variable name.
label_get <- function(x) {
  ifelse(
    x %in% vars_dict$var,
    vars_dict[vars_dict$var == x, "name", drop = TRUE],
    x
  )
}

# Apply the label_get function to a vector of variable names.
vars_label <- function(x) {
  sapply(x, label_get, USE.NAMES = FALSE)
}

# Function to retrieve the original variable name from a given label.
var_get <- function(x) {
  ifelse(
    x %in% vars_dict$name,
    vars_dict[vars_dict$name == x, "var", drop = TRUE],
    x
  )
}

# ============================================================================ #
# Model Evaluation ------------------------------------------------------------
# ============================================================================ #

jointfit_tbl <- function(
  jointFit1,
  jointFit3,
  round_digits_hr = 2,
  round_digits_p = 3,
  heads = c("Cohort 1", "Cohort 3")
) {
  # Extract survival summary objects
  surv1 <- summary(jointFit1)$Survival
  surv3 <- summary(jointFit3)$Survival

  # Build summary tibble
  tbl <- tibble(
    variable = rownames(surv1),
    cohort1_hr = round(exp(surv1[, 1]), round_digits_hr),
    cohort1_lower = round(exp(surv1[, 3]), round_digits_hr),
    cohort1_upper = round(exp(surv1[, 4]), round_digits_hr),
    cohort1_p = round(surv1[, 5], round_digits_p),
    cohort3_hr = round(exp(surv3[, 1]), round_digits_hr),
    cohort3_lower = round(exp(surv3[, 3]), round_digits_hr),
    cohort3_upper = round(exp(surv3[, 4]), round_digits_hr),
    cohort3_p = round(surv3[, 5], round_digits_p)
  ) %>%
    # Rename variables for clarity
    mutate(
      variable = case_when(
        variable == "genderMale" ~ "Male",
        variable == "age" ~ "Age (years)",
        variable == "nyhaI" ~ "NYHA I (ref: II)",
        variable == "nyhaIII" ~ "NYHA III (ref: II)",
        variable == "nyhaIV" ~ "NYHA IV (ref: II)",
        variable == "nyhaUndetermined" ~ "NYHA Undetermined",
        variable == "value(x6mw_dist_meter)" ~ "6MWT Distance (m)",
        variable == "slope(x6mw_dist_meter)" ~ "6MWT Distance Slope",
        variable == "ischemic_etiologyTRUE" ~ "Ischemic Etiology",
        variable == "sbp" ~ "Systolic Blood Pressure (mmHg)",
        variable == "efHFmrEF" ~ "HFmrEF (ref: HFrEF)",
        variable == "efHFpEF" ~ "HFpEF (ref: HFrEF)",
        variable == "fusidTRUE" ~ "Furosemide Treatment",
        variable == "k_spareTRUE" ~
          "Mineralocorticoid Receptor Antagonist Treatment",
        variable == "allopurinolTRUE" ~ "Allopurinol Treatment",
        variable == "statinTRUE" ~ "Statin Treatment",
        variable == "ace_arb_arniTRUE" ~ "ACEi, ARB or ARNI Treatment",
        variable == "beta_blockerTRUE" ~ "Beta Blockers Treatment",
        variable == "sglt2TRUE" ~
          "Sodium-Glucose Cotransporter-2 Inhibitors Treatment",
        variable == "na" ~ "Sodium (mmol/L)",
        variable == "hgb" ~ "Hemoglobin (g/dL)",
        variable == "bmi" ~ "BMI (kg/m²)",
        variable == "smokepast smoker" ~ "Past Smoker",
        variable == "smoke1-10" ~ "1-10 Cigarettes/Day",
        variable == "smoke>10" ~ ">10 Cigarettes/Day",
        TRUE ~ variable
      )
    )

  # Create formatted flextable
  ft <- tbl %>%
    flextable() %>%
    set_header_labels(
      variable = "Variable",
      cohort1_hr = "HR",
      cohort1_lower = "Lower CI",
      cohort1_upper = "Upper CI",
      cohort1_p = "P-value",
      cohort3_hr = "HR",
      cohort3_lower = "Lower CI",
      cohort3_upper = "Upper CI",
      cohort3_p = "P-value"
    ) %>%
    add_header_row(
      values = c("", heads),
      colwidths = c(1, 4, 4)
    ) %>%
    align(align = "center", part = "all") %>%
    bold(part = "header") %>%
    color(i = ~ cohort1_p < 0.05, j = "cohort1_p", color = "red") %>%
    color(i = ~ cohort3_p < 0.05, j = "cohort3_p", color = "red") %>%
    width(j = 1:9, width = 1.5) %>%
    padding(padding = 5, part = "all") %>%
    autofit() %>%
    theme_vanilla()

  return(ft)
}

create_train_test <- function(df, follow_up, V = 5, id_var = "id", seed = 229) {
  # Create the cross-validation folds using create_folds
  CVdats <- create_folds(df, V = V, id_var = id_var, seed = seed)

  # Combine the training fold with the portion of the testing fold that meets the time condition
  # df_train <- bind_rows(
  #   CVdats$training[[1]],
  #   CVdats$testing[[1]] %>% filter(time <= follow_up)
  # )
  df_train <- CVdats$training[[1]]
  # Use the full testing fold as the test set
  df_test <- CVdats$testing[[1]]

  # Return a list containing both training and testing datasets
  list(train = df_train, test = df_test)
}

create_cox_df <- function(data, Tstart) {
  if (Tstart > 0) {
    data <- data %>%
      filter(time <= Tstart, fup_time > Tstart) %>%
      arrange(id, time) %>%
      group_by(id) %>%
      filter(time == max(time)) %>%
      ungroup() %>%
      mutate(fup_time = fup_time - Tstart, time = 0)
  } else {
    data <- data %>%
      arrange(id, time) %>%
      group_by(id) %>%
      slice(1) %>%
      ungroup()
  }
  return(data)
}


fit_models <- function(data_list, full = FALSE) {
  # Extract the training dataset from the input list.
  # Here we assume data_list was created by your previous function and has an element named "train".
  long_df_train <- data_list$train

  # Derive a dataset for the Cox model.
  # We assume that the survival covariates (e.g., gender, age, nyha, fup_time, mortality_status)
  # are the same across repeated measurements for each subject. Hence, we take the first observation per id.
  cox_df_train <- long_df_train %>%
    arrange(id, time) %>%
    group_by(id) %>%
    slice(1) %>%
    ungroup()

  # 1. Fit the linear mixed-effects model on the training data.
  lme_6min_train <- lme(
    x6mw_dist_meter ~ ns(time, df = 3) * (age + gender) + nyha,
    data = long_df_train,
    random = ~ ns(time, df = 3) | id,
    control = lmeControl(opt = "optim")
  )

  # 2. Fit the Cox proportional hazards model on the survival training data.
  CoxFit_train <- coxph(
    Surv(fup_time, mortality_status) ~ gender + age + nyha,
    data = cox_df_train,
    model = TRUE,
    x = TRUE,
    y = TRUE
  )

  # 3. Combine the two models into a joint model.
  jointFit_train <- jm(
    CoxFit_train,
    lme_6min_train,
    time_var = "time",
    functional_forms = ~ value(x6mw_dist_meter) + slope(x6mw_dist_meter)
  )

  # 4. Fit a Cox model that also includes the 6-min walk distance for comparison.
  if (full) {
    CoxFit_train_full <- coxph(
      as.formula(paste(
        "Surv(fup_time, mortality_status) ~",
        paste(small_seattle_vars, collapse = " + ")
      )),
      data = cox_df_train,
      model = TRUE,
      x = TRUE,
      y = TRUE
    )

    jointFit_train_full <- jm(
      CoxFit_train_full,
      lme_6min_train,
      time_var = "time",
      functional_forms = ~ value(x6mw_dist_meter) + slope(x6mw_dist_meter)
    )
    models <- list(
      joint_model = jointFit_train,
      joint_model_full = jointFit_train_full
    )
  } else {
    CoxFit_six_train <- coxph(
      Surv(fup_time, mortality_status) ~ gender + age + nyha + x6mw_dist_meter,
      data = cox_df_train,
      model = TRUE,
      x = TRUE,
      y = TRUE
    )

    models <- list(
      joint_model = jointFit_train,
      cox_model = CoxFit_six_train
    )
  }

  # Return a list containing all the models.
  return(models)
}
# -------------------------------------------------------------------------- #
## Bootstrap Evaluation -------------------------------------------------------
# -------------------------------------------------------------------------- #

# -----------------------------------------------------------------------------#
# Time-Dependent Metric Computation for Dynamic Survival Predictions via Cox Model
#
# Overview:
#   This function computes time-dependent performance metrics for a Cox proportional
#   hazards model, specifically the time-dependent Area Under the Curve (AUC),
#   Receiver Operating Characteristic (ROC) metrics, and Brier score.
#   It evaluates the discriminative and predictive performance of the model
#   at a specified prediction horizon by comparing the model-generated risk scores
#   with the observed survival outcomes.
#
# Parameters:
#   object:
#     A fitted Cox proportional hazards model (an object of class "coxph"),
#     obtained using the 'coxph' function.
#
#   newdata:
#     A data frame containing the prediction dataset. It must include the following columns:
#       - id: Unique identifier for each subject.
#       - time: The time at which covariates were recorded.
#       - fup_time: The follow-up time at which the event or censoring occurred.
#       - mortality_status: Event indicator (1 if event occurred, 0 if censored).
#
#   Tstart:
#     The starting time for prediction. Only subjects at risk at Tstart are included
#     in the analysis. Default is 0.
#
#   Thoriz:
#     The prediction horizon time point. It must be greater than Tstart.
#     If Thoriz is not provided, it is computed as Tstart + Dt.
#
#   Dt:
#     The time interval for prediction, used to compute Thoriz if not explicitly provided.
#
#   id_var, time_var, Time_var, event_var:
#     Column names in 'newdata' corresponding to the subject ID, covariate measurement time,
#     follow-up time, and event indicator, respectively. Defaults: "id", "time", "fup_time", "mortality_status".
#
#   type:
#     A character value specifying the output type:
#       - "auc": Computes time-dependent AUC.
#       - "roc": Returns ROC-related metrics (e.g., TP, FP, thresholds, F1 score, Youden index).
#       - "brier": Computes the Brier score for calibration assessment.
#
#   ...:
#     Additional arguments passed to the underlying prediction functions.
#
# Returns:
#   A list containing the following elements:
#     - auc: Computed time-dependent AUC (if type = "auc" or "roc").
#     - brier: Computed Brier score (if type = "brier").
#     - Tstart: Adjusted starting time for prediction.
#     - Thoriz: Adjusted prediction horizon time.
#     - nr: Number of unique subjects included in computation.
#     - classObject: Class of the input Cox model.
#     - nameObject: Name of the input Cox model object.
#
#   For type = "roc", additional elements include:
#     - TP: True positive rates.
#     - FP: False positive rates.
#     - thrs: Threshold values used in ROC analysis.
#     - F1score: Threshold maximizing the F1 score.
#     - Youden: Threshold maximizing Youden's index.
#
# Notes:
#   - A small perturbation (1e-06) is added to Tstart and Thoriz to mitigate numerical issues.
#   - Only the most recent covariate measurements for subjects at risk at Tstart are used.
#   - The function internally filters the dataset to ensure only subjects at risk at Tstart are included.
#
# -----------------------------------------------------------------------------#
tvEVAL <- function(
  object,
  newdata,
  Tstart = 0,
  Thoriz = NULL,
  Dt = NULL,
  id_var = "id",
  time_var = "time",
  Time_var = "fup_time",
  event_var = "mortality_status",
  type = c("auc", "roc", "brier"),
  ...
) {
  # Ensure the provided model is a Cox proportional hazards model.
  if (!inherits(object, "coxph")) {
    stop("This function can only be used with 'coxph' objects.\n")
  }

  # Validate that 'newdata' is a non-empty data frame.
  if (!is.data.frame(newdata) || nrow(newdata) == 0) {
    stop("'newdata' must be a non-empty data frame.\n")
  }

  # Ensure at least one of 'Thoriz' (prediction horizon) or 'Dt' (time interval) is provided.
  if (is.null(Thoriz) && is.null(Dt)) {
    stop("Either 'Thoriz' or 'Dt' must be specified.\n")
  }

  # If 'Thoriz' is provided, ensure it is greater than 'Tstart'.
  if (!is.null(Thoriz) && Thoriz <= Tstart) {
    stop("'Thoriz' must be larger than 'Tstart'.\n")
  }

  # Compute 'Thoriz' if it is not explicitly provided.
  if (is.null(Thoriz)) {
    Thoriz <- Tstart + Dt
  }

  # Compute 'Dt' if it is not explicitly provided.
  if (is.null(Dt)) {
    Dt <- Thoriz - Tstart
  }

  # Verify that 'newdata' contains the required columns.
  if (is.null(newdata[[id_var]])) {
    stop("Column '", id_var, "' not found in newdata.")
  }
  if (is.null(newdata[[time_var]])) {
    stop("Column '", time_var, "' not found in newdata.")
  }
  if (is.null(newdata[[Time_var]])) {
    stop("Column '", Time_var, "' not found in newdata.")
  }
  if (is.null(newdata[[event_var]])) {
    stop("Column '", event_var, "' not found in newdata.")
  }

  # Apply a small numerical adjustment to Tstart and Thoriz to avoid precision issues.
  Tstart <- Tstart + 1e-06
  Thoriz <- Thoriz + 1e-06

  # Filter dataset to retain only subjects at risk at Tstart.
  # If a subject has multiple records, only the most recent record before Tstart is kept.
  newdata <- create_cox_df(newdata, Tstart)

  if (type %in% c("auc", "roc")) {
    # Compute the linear predictor (risk score) from the fitted Cox model.
    lp <- predict(object, newdata = newdata, type = "lp")

    # Compute the time-dependent AUC using survivalROC.
    auc <- survivalROC(
      Stime = newdata[[Time_var]],
      status = newdata[[event_var]],
      marker = lp,
      predict.time = Thoriz,
      span = 0.25 * nrow(newdata)^(-0.20)
    )

    if (type == "auc") {
      # Output AUC result with additional model metadata.
      out <- list(
        auc = auc$AUC,
        Tstart = Tstart,
        Thoriz = Thoriz,
        nr = length(unique(newdata[[id_var]])),
        classObject = class(object),
        nameObject = deparse(substitute(object))
      )
      class(out) <- "tvAUC_coxph"
    } else {
      # Compute ROC-related metrics, including thresholds, F1 score, and Youden's index.
      metric_df <- tibble(
        TP = auc$TP,
        FP = auc$FP,
        thres = auc$cut.values
      ) %>%
        mutate(
          precision = TP / (TP + FP),
          recall = TP, # Sensitivity is equivalent to recall.
          F1 = 2 * (precision * recall) / (precision + recall),
          Youden = TP - FP # Youden's index = sensitivity + specificity - 1
        )

      # Extract the median threshold corresponding to the highest F1 score.
      F1score <- median(
        metric_df[metric_df$F1 == max(metric_df$F1, na.rm = TRUE), ]$thres,
        na.rm = TRUE
      )

      # Extract the median threshold corresponding to the highest Youden's index.
      Youden <- median(
        metric_df[
          metric_df$Youden == max(metric_df$Youden, na.rm = TRUE),
        ]$thres,
        na.rm = TRUE
      )

      out <- list(
        TP = metric_df$TP,
        FP = metric_df$FP,
        auc = auc$AUC,
        thrs = metric_df$thres,
        F1score = F1score,
        Youden = Youden,
        Tstart = Tstart,
        Thoriz = Thoriz,
        nr = length(unique(newdata[[id_var]])),
        classObject = class(object),
        nameObject = deparse(substitute(object))
      )
      class(out) <- "tvROC_coxph"
    }
  } else {
    # Compute the Brier score for the specified prediction horizon.
    brier <- brier(object, newdata = newdata, times = Thoriz)$brier

    # Output Brier score result with model metadata.
    out <- list(
      brier = brier,
      Tstart = Tstart,
      Thoriz = Thoriz,
      nr = length(unique(newdata[[id_var]])),
      classObject = class(object),
      nameObject = deparse(substitute(object))
    )
    class(out) <- "tvBrier_coxph"
  }

  return(out)
}

# -----------------------------------------------------------------------------#
# Print Methods for Time-Dependent Performance Metrics in Survival Analysis
#
# Overview:
#   These functions provide print methods for objects of class "tvAUC_coxph",
#   "tvROC_coxph", and "tvBrier_coxph". Each function prints relevant
#   time-dependent performance metrics for Cox proportional hazards models or
#   joint models.
#
# Functions:
#   - print.tvAUC_coxph: Prints estimated time-dependent AUC.
#   - print.tvROC_coxph: Prints ROC-related metrics, including key thresholds.
#   - print.tvBrier_coxph: Prints the computed Brier score for calibration assessment.
#
# Notes:
#   - Each function ensures compatibility with its respective class.
#   - If the model class is "jm" (joint model), it is indicated in the output.
#   - The number of subjects at risk at Tstart is reported.
#
# -----------------------------------------------------------------------------#

# Print method for "tvAUC_coxph" objects (Time-Dependent AUC)
print.tvAUC_coxph <- function(x, digits = 4, ...) {
  # Ensure correct object class
  if (!inherits(x, "tvAUC_coxph")) {
    stop("Use only with 'tvAUC_coxph' objects.\n")
  }

  # Print model type (Cox or Joint Model)
  if (x$class == "jm") {
    cat("\n\tTime-dependent AUC for the Joint Model:", x$nameObject)
  } else {
    cat("\n\tTime-dependent AUC for the Cox Model:", x$nameObject)
  }

  # Print estimated AUC
  cat("\n\nEstimated AUC: ", round(x$auc, digits))

  # Optional: Print confidence intervals for AUC if available (currently commented out)
  # cat("\n\nEstimated AUC: ", round(x$auc, digits),
  #     " (95% CI: ", round(x$low_auc, digits), "-", round(x$upp_auc, digits),
  #     ")", sep = "")

  # Print time horizon at which AUC was estimated
  cat("\nAt time:", round(x$Thoriz, digits))

  # Print number of subjects at risk at prediction start time
  cat(
    "\nUsing information up to time:",
    round(x$Tstart, digits),
    " (",
    x$nr,
    " subjects still at risk)",
    sep = ""
  )

  cat("\n\n")

  invisible(x)
}

# Print method for "tvROC_coxph" objects (Time-Dependent ROC Metrics)
print.tvROC_coxph <- function(x, digits = 4, ...) {
  # Ensure correct object class
  if (!inherits(x, "tvROC_coxph")) {
    stop("Use only with 'tvROC_coxph' objects.\n")
  }

  # Print model type (Cox or Joint Model)
  if (x$class == "jm") {
    cat(
      "\n\tTime-dependent ROC Metrics for the Joint Model:",
      x$nameObject
    )
  } else {
    cat("\n\tTime-dependent ROC Metrics for the Cox Model:", x$nameObject)
  }

  # Print estimated AUC
  cat("\n\nEstimated AUC: ", round(x$auc, digits))

  # Print time horizon at which ROC metrics were computed
  cat("\nAt time:", round(x$Thoriz, digits))

  # Print number of subjects at risk at prediction start time
  cat(
    "\nUsing information up to time:",
    round(x$Tstart, digits),
    " (",
    x$nr,
    " subjects still at risk)",
    sep = ""
  )

  # Print key ROC thresholds
  cat("\n\nOptimal Thresholds:")
  cat("\n- F1-score Optimal Threshold:", round(x$F1score, digits))
  cat("\n- Youden’s Index Optimal Threshold:", round(x$Youden, digits))

  cat("\n\n")

  invisible(x)
}

# Print method for "tvBrier_coxph" objects (Time-Dependent Brier Score)
print.tvBrier_coxph <- function(x, digits = 4, ...) {
  # Ensure correct object class
  if (!inherits(x, "tvBrier_coxph")) {
    stop("Use only with 'tvBrier_coxph' objects.\n")
  }

  # Print model type (Cox or Joint Model)
  if (x$class == "jm") {
    cat(
      "\n\tTime-dependent Brier Score for the Joint Model:",
      x$nameObject
    )
  } else {
    cat("\n\tTime-dependent Brier Score for the Cox Model:", x$nameObject)
  }

  # Print estimated Brier score
  cat("\n\nEstimated Brier Score: ", round(x$brier, digits))

  # Print time horizon at which the Brier score was computed
  cat("\nAt time:", round(x$Thoriz, digits))

  # Print number of subjects at risk at prediction start time
  cat(
    "\nUsing information up to time:",
    round(x$Tstart, digits),
    " (",
    x$nr,
    " subjects still at risk)",
    sep = ""
  )

  cat("\n\n")

  invisible(x)
}

# -----------------------------------------------------------------------------#
# Bootstrap Iteration for Assessing Dynamic AUC and Brier Score Stability
#
# Overview:
#   This function performs a single bootstrap iteration to evaluate the stability
#   and variability of time-dependent AUC and Brier score estimates for survival models.
#   It resamples subjects (with replacement) from the test dataset, reassigns
#   sequential IDs, and constructs bootstrap datasets. The function then computes
#   dynamic AUC and Brier scores for both joint and Cox models across multiple
#   prediction intervals.
#
# Parameters:
#   iteration:
#     The current bootstrap iteration number (used for tracking).
#
#   test_df:
#     A data frame containing the test dataset. Must include an 'id' column
#     uniquely identifying each subject.
#
#   models:
#     A list containing the pre-fitted models required for metric computation:
#       - joint_model: A joint model (from JMbayes2) for dynamic survival prediction.
#       - cox_model: A Cox proportional hazards model for dynamic survival prediction.
#
# Returns:
#   A tibble containing the computed AUC and Brier scores for multiple prediction intervals:
#     - auc_joint_X: AUC from the joint model for different time horizons.
#     - auc_cox_X: AUC from the Cox model for different time horizons.
#     - brier_joint_X: Brier score from the joint model for different time horizons.
#     - brier_cox_X: Brier score from the Cox model for different time horizons.
#   The X values indicate the time horizon (e.g., 36, 60 months).
#
# Dependencies:
#   This function requires the 'JMbayes2', 'survivalROC', and 'riskRegression' packages.
#
# -----------------------------------------------------------------------------#

bootstrap_iteration <- cmpfun(function(
  iteration,
  test_df,
  models,
  full = FALSE
) {
  library(JMbayes2)
  library(survivalROC)
  library(riskRegression)

  # -------------------------------#
  # Resample Subject IDs with Replacement
  # -------------------------------#
  # This step creates a bootstrap sample by resampling subject IDs with replacement.
  # The corresponding data for each sampled subject is extracted, and a new sequential
  # ID is assigned to ensure proper structure in the bootstrap dataset.
  unique_ids <- unique(test_df$id)
  boot_ids <- sample(unique_ids, size = length(unique_ids), replace = TRUE)

  boot_data <- bind_rows(lapply(seq_along(boot_ids), function(j) {
    test_df %>%
      filter(id == boot_ids[j]) %>%
      mutate(id = j) # Assign a new sequential ID
  }))

  # -------------------------------#
  # Compute Dynamic AUC Measures
  # -------------------------------#
  # Computes the time-dependent AUC for both joint and Cox models at different
  # prediction horizons (e.g., 36 and 60 months from baseline, as well as from 12 months onward).

  auc_joint_36 <- tvAUC(
    object = models$joint_model,
    newdata = boot_data,
    Tstart = 0,
    Dt = 36,
    type_weights = "IPCW"
  )
  auc_joint_60 <- tvAUC(
    object = models$joint_model,
    newdata = boot_data,
    Tstart = 0,
    Dt = 60,
    type_weights = "IPCW"
  )
  if (full) {
    auc_joint_full_36 <- tvAUC(
      object = models$joint_model_full,
      newdata = boot_data,
      Tstart = 0,
      Dt = 36,
      type_weights = "IPCW"
    )
    auc_joint_full_60 <- tvAUC(
      object = models$joint_model_full,
      newdata = boot_data,
      Tstart = 0,
      Dt = 60,
      type_weights = "IPCW"
    )
  } else {
    auc_cox_36 <- tvEVAL(
      object = models$cox_model,
      newdata = boot_data,
      Tstart = 0,
      Dt = 36,
      type = "auc"
    )
    auc_cox_60 <- tvEVAL(
      object = models$cox_model,
      newdata = boot_data,
      Tstart = 0,
      Dt = 60,
      type = "auc"
    )
  }
  auc_joint_1_36 <- tvAUC(
    object = models$joint_model,
    newdata = boot_data,
    Tstart = 12,
    Dt = 24,
    type_weights = "IPCW"
  )
  auc_joint_1_60 <- tvAUC(
    object = models$joint_model,
    newdata = boot_data,
    Tstart = 12,
    Dt = 48,
    type_weights = "IPCW"
  )
  if (full) {
    auc_joint_full_1_36 <- tvAUC(
      object = models$joint_model_full,
      newdata = boot_data,
      Tstart = 12,
      Dt = 24,
      type_weights = "IPCW"
    )
    auc_joint_full_1_60 <- tvAUC(
      object = models$joint_model_full,
      newdata = boot_data,
      Tstart = 12,
      Dt = 48,
      type_weights = "IPCW"
    )
  } else {
    auc_cox_1_36 <- tvEVAL(
      object = models$cox_model,
      newdata = boot_data,
      Tstart = 12,
      Dt = 24,
      type = "auc"
    )
    auc_cox_1_60 <- tvEVAL(
      object = models$cox_model,
      newdata = boot_data,
      Tstart = 12,
      Dt = 48,
      type = "auc"
    )
  }
  # -------------------------------#
  # Compute Dynamic Brier Measures
  # -------------------------------#
  # The Brier score is a calibration metric that quantifies the accuracy of predicted
  # survival probabilities. It is computed for different time horizons using both
  # the joint and Cox models.

  brier_joint_36 <- tvBrier(
    object = models$joint_model,
    newdata = boot_data,
    Tstart = 0,
    Dt = 36,
    type_weights = "IPCW"
  )
  brier_joint_60 <- tvBrier(
    object = models$joint_model,
    newdata = boot_data,
    Tstart = 0,
    Dt = 60,
    type_weights = "IPCW"
  )
  if (full) {
    brier_joint_full_36 <- tvBrier(
      object = models$joint_model_full,
      newdata = boot_data,
      Tstart = 0,
      Dt = 36,
      type_weights = "IPCW"
    )
    brier_joint_full_60 <- tvBrier(
      object = models$joint_model_full,
      newdata = boot_data,
      Tstart = 0,
      Dt = 60,
      type_weights = "IPCW"
    )
  } else {
    brier_cox_36 <- tvEVAL(
      object = models$cox_model,
      newdata = boot_data,
      Tstart = 0,
      Dt = 36,
      type = "brier"
    )
    brier_cox_60 <- tvEVAL(
      object = models$cox_model,
      newdata = boot_data,
      Tstart = 0,
      Dt = 60,
      type = "brier"
    )
  }

  brier_joint_1_36 <- tvBrier(
    object = models$joint_model,
    newdata = boot_data,
    Tstart = 12,
    Dt = 24,
    type_weights = "IPCW"
  )
  brier_joint_1_60 <- tvBrier(
    object = models$joint_model,
    newdata = boot_data,
    Tstart = 12,
    Dt = 48,
    type_weights = "IPCW"
  )
  if (full) {
    brier_joint_full_1_36 <- tvBrier(
      object = models$joint_model_full,
      newdata = boot_data,
      Tstart = 12,
      Dt = 24,
      type_weights = "IPCW"
    )
    brier_joint_full_1_60 <- tvBrier(
      object = models$joint_model_full,
      newdata = boot_data,
      Tstart = 12,
      Dt = 48,
      type_weights = "IPCW"
    )
  } else {
    brier_cox_1_36 <- tvEVAL(
      object = models$cox_model,
      newdata = boot_data,
      Tstart = 12,
      Dt = 24,
      type = "brier"
    )
    brier_cox_1_60 <- tvEVAL(
      object = models$cox_model,
      newdata = boot_data,
      Tstart = 12,
      Dt = 48,
      type = "brier"
    )
  }
  # -------------------------------#
  # Return AUC and Brier Results as a Tibble
  # -------------------------------#
  # Combines all computed performance metrics into a tibble for further statistical analysis.
  if (full) {
    tibble(
      auc_joint_36 = auc_joint_36$auc,
      auc_joint_full_36 = auc_joint_full_36$auc,
      auc_joint_60 = auc_joint_60$auc,
      auc_joint_full_60 = auc_joint_full_60$auc,
      auc_joint_1_36 = auc_joint_1_36$auc,
      auc_joint_full_1_36 = auc_joint_full_1_36$auc,
      auc_joint_1_60 = auc_joint_1_60$auc,
      auc_joint_full_1_60 = auc_joint_full_1_60$auc,
      brier_joint_36 = brier_joint_36$Brier,
      brier_joint_full_36 = brier_joint_full_36$Brier,
      brier_joint_60 = brier_joint_60$Brier,
      brier_joint_full_60 = brier_joint_full_60$Brier,
      brier_joint_1_36 = brier_joint_1_36$Brier,
      brier_joint_full_1_36 = brier_joint_full_1_36$Brier,
      brier_joint_1_60 = brier_joint_1_60$Brier,
      brier_joint_full_1_60 = brier_joint_full_1_60$Brier
    )
  } else {
    tibble(
      auc_joint_36 = auc_joint_36$auc,
      auc_cox_36 = auc_cox_36$auc,
      auc_joint_60 = auc_joint_60$auc,
      auc_cox_60 = auc_cox_60$auc,
      auc_joint_1_36 = auc_joint_1_36$auc,
      auc_cox_1_36 = auc_cox_1_36$auc,
      auc_joint_1_60 = auc_joint_1_60$auc,
      auc_cox_1_60 = auc_cox_1_60$auc,
      brier_joint_36 = brier_joint_36$Brier,
      brier_cox_36 = brier_cox_36$brier,
      brier_joint_60 = brier_joint_60$Brier,
      brier_cox_60 = brier_cox_60$brier,
      brier_joint_1_36 = brier_joint_1_36$Brier,
      brier_cox_1_36 = brier_cox_1_36$brier,
      brier_joint_1_60 = brier_joint_1_60$Brier,
      brier_cox_1_60 = brier_cox_1_60$brier
    )
  }
})

# -----------------------------------------------------------------------------#
# Plot Bootstrap Estimates for AUC and Brier Score
#
# Overview:
#   This function generates a visualization of bootstrap-estimated performance metrics
#   (AUC and Brier score) for different survival models, such as the Cox model
#   and Joint Model. The function allows filtering, faceting, and customization of
#   aesthetics, making it useful for comparing model performance across follow-up
#   times and cohorts.
#
# Parameters:
#   data:
#     A data frame containing bootstrap results, including columns for model type,
#     metric type, follow-up time, confidence intervals, and estimated values.
#
#   model_filter:
#     A character vector specifying which models to include in the plot.
#     Default: c("Cox model", "Joint Model").
#
#   metric_filter:
#     A character vector specifying which performance metrics to plot.
#     Default: c("AUC", "Brier score").
#
#   x_var, y_var:
#     The column names in `data` that represent the x-axis and y-axis values.
#     Default: x_var = "fup_time" (follow-up time), y_var = "mean" (estimated metric).
#
#   color_var:
#     The column name in `data` that defines color grouping (typically model type).
#     Default: "model_type".
#
#   ymin_var, ymax_var:
#     Column names specifying the lower and upper confidence interval bounds for error bars.
#     Default: "lower_ci" and "upper_ci".
#
#   color_palette:
#     The color palette for the plot. Uses `scale_color_brewer()`.
#     Default: "Set1".
#
#   hline_position:
#     A numeric value specifying a horizontal reference line position (optional).
#
#   hline_color:
#     The color of the horizontal reference line (if `hline_position` is provided).
#     Default: "gray40".
#
#   facet_rows, facet_cols:
#     Column names used for faceting rows and columns.
#     Default: facet_rows = "fup_start", facet_cols = "cohort".
#
#   facet_scales:
#     Scaling behavior for facets ("fixed", "free", "free_x", or "free_y").
#     Default: "fixed".
#
#   facet_type:
#     Type of faceting: "grid" (default) or "wrap".
#
#   plot_title:
#     Title of the generated plot.
#
# Returns:
#   A ggplot2 object displaying mean bootstrap estimates with 95% confidence intervals
#   for the specified models and performance metrics.
#
# Notes:
#   - The function allows flexible filtering of models and metrics.
#   - Faceting enables comparison across different cohorts and follow-up periods.
#   - Error bars represent 95% confidence intervals.
#   - If `hline_position` is provided, a dashed reference line is added to the plot.
#
# -----------------------------------------------------------------------------#

plot_boot_df <- function(
  data,
  model_filter = c("Single 6MWT", "Serial 6MWT"),
  metric_filter = c("AUC", "Brier Score"),
  x_var = "fup_time",
  y_var = "mean",
  x_lab = "Follow-up Time",
  y_lab = "Metric Value",
  color_var = "model_type",
  ast_y = 0.85,
  mirror = FALSE,
  ymin_var = "lower_ci",
  ymax_var = "upper_ci",
  color_palette = "Set1",
  hline_position = NULL,
  hline_color = "gray40",
  facet_rows = "fup_start",
  facet_cols = "cohort",
  facet_scales = "fixed",
  facet_type = "grid",
  plot_title = "Bootstrap Estimates",
  plot_subtitle = "Mean Bootstrap Estimates with 95% Confidence Intervals"
) {
  # -------------------------------#
  # Filter Data Based on User Input
  # -------------------------------#
  # Retains only the selected models and performance metrics.
  if (mirror) {
    p <- data %>%
      mutate_at(
        vars(mean, lower_ci, upper_ci),
        ~ ifelse(fup_start == "Baseline", ., -.)
      ) %>%
      filter(model_type %in% model_filter, metric %in% metric_filter) %>%
      ggplot(aes(
        x = .data[[x_var]],
        y = .data[[y_var]],
        color = .data[[color_var]]
      )) +
      # Add point estimates with 95% confidence interval error bars.
      geom_pointrange(
        aes(xmin = .data[[ymin_var]], xmax = .data[[ymax_var]]),
        size = 0.5,
        linewidth = 1,
        position = position_dodge(width = 0.5)
      ) +
      scale_x_continuous(labels = function(x) abs(x))
  } else {
    p <- data %>%
      filter(model_type %in% model_filter, metric %in% metric_filter) %>%
      ggplot(aes(
        x = .data[[x_var]],
        y = .data[[y_var]],
        color = .data[[color_var]]
      )) +
      # Add point estimates with 95% confidence interval error bars.
      geom_pointrange(
        aes(ymin = .data[[ymin_var]], ymax = .data[[ymax_var]]),
        size = 0.5,
        linewidth = 1,
        position = position_dodge(width = 0.5)
      )
  }

  # -------------------------------#
  # Add Significance Asterisks
  # -------------------------------#
  sig_data <- data %>%
    filter(
      model_type == "Difference",
      metric %in% metric_filter,
      sign(lower_ci) == sign(upper_ci)
    ) # Subset significant results

  if (nrow(sig_data) > 0) {
    p <- p +
      geom_text(
        data = sig_data,
        aes(y = ast_y, label = "*"),
        size = 6,
        color = "black"
      )
  }

  # -------------------------------#
  # Customize Labels and Theme
  # -------------------------------#
  p <- p +
    # Apply color scheme from RColorBrewer.
    scale_color_brewer(palette = color_palette) +
    labs(
      # title    = plot_title,
      # subtitle = plot_subtitle,
      x = x_lab,
      y = y_lab,
      color = "Model Type"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom"
    )

  # -------------------------------#
  # Apply Faceting Based on User Input
  # -------------------------------#
  if (facet_type == "grid") {
    p <- p +
      facet_grid(
        as.formula(paste(facet_rows, "~", facet_cols)),
        scales = facet_scales
      )
  } else if (facet_type == "wrap") {
    p <- p +
      facet_wrap(
        vars(.data[[facet_rows]], .data[[facet_cols]]),
        scales = facet_scales
      )
  }

  # -------------------------------#
  # Add Horizontal Reference Line (Optional)
  # -------------------------------#
  if (!is.null(hline_position)) {
    if (mirror) {
      p <- p +
        geom_vline(
          xintercept = hline_position,
          linetype = "dashed",
          color = hline_color
        )
    } else {
      p <- p +
        geom_hline(
          yintercept = hline_position,
          linetype = "dashed",
          color = hline_color
        )
    }
  }

  return(p)
}

# -------------------------------------------------------------------------- #
## Discrimination Evaluation --------------------------------------------------
# -------------------------------------------------------------------------- #

#------------------------------------------------------------------------------#
# Function: create_roc_data
#
# Description:
#   Generates time-dependent ROC curve data for a given model at multiple
#   follow-up intervals. The function supports only two types of models:
#     - Joint models (class "jm")
#     - Cox models (class "coxph")
#
#   For each specified follow-up time, it computes the false positive rate (FP),
#   true positive rate (TP), and the Area Under the Curve (AUC), then aggregates
#   these metrics into a tidy tibble.
#
# Arguments:
#   model           - A fitted model object. Must be either:
#                      - a joint model (class "jm"), or
#                      - a Cox model (class "coxph").
#   long_data       - A data frame containing the longitudinal data for predictions.
#   Tstart          - The starting time point for predictions.
#   follow_up_times - A numeric vector of follow-up times (Dt) at which to compute ROC metrics.
#
# Returns:
#   A tibble with the following columns:
#     FP       : False Positive Rate (1 - specificity)
#     TP       : True Positive Rate (sensitivity)
#     FollowUp : The follow-up time interval (Dt)
#     AUC      : Area Under the Curve (rounded to 3 decimal places)
#
# Notes:
#   The function will throw an error if the provided model is not of class "jm" or "coxph".
#------------------------------------------------------------------------------#
create_roc_data <- function(model, long_data, Tstart, follow_up_times) {
  # Validate that the model is either a joint model (jm) or a Cox model (coxph)
  if (!inherits(model, "jm") && !inherits(model, "coxph")) {
    stop("Model must be of class 'jm' or 'coxph'.")
  }

  # Compute ROC metrics for each follow-up time using lapply.
  roc_list <- lapply(follow_up_times, function(Dt) {
    # Select the appropriate ROC function based on the model type.
    if (inherits(model, "jm")) {
      cox_data <- create_cox_df(long_data, Tstart)
      roc_result <- tvROC(model, newdata = cox_data, Tstart = Tstart, Dt = Dt)
      # Extract and round the AUC value from the ROC result.
      auc_value <- round(tvAUC(roc_result)$auc, 3)
    } else {
      # Must be a Cox model (coxph)

      roc_result <- tvEVAL(
        model,
        newdata = long_data,
        Tstart = Tstart,
        Dt = Dt,
        type = "roc"
      )
      # Extract and round the AUC value from the ROC result.
      auc_value <- round(roc_result$auc, 3)
    }
    # Return a data frame with the ROC metrics for the current follow-up time.
    data.frame(
      FP = roc_result$FP, # False Positive Rate
      TP = roc_result$TP, # True Positive Rate
      FollowUp = Dt, # Follow-up time interval
      AUC = auc_value # Area Under the Curve
    )
  })

  # Combine the individual ROC data frames into a single tibble.
  roc_data <- dplyr::bind_rows(roc_list)

  return(roc_data)
}

#------------------------------------------------------------------------------#
# Function: plot_tvROC
#
# Description:
#   Generates a time-dependent ROC curve plot using ggplot2. The function takes
#   a data frame containing false positive rates (FP), true positive rates (TP),
#   follow-up intervals, and corresponding AUC values, then returns a ggplot object
#   that displays the ROC curves for each follow-up time.
#
# Parameters:
#   roc_data   : A data frame with the following columns:
#                  - FP       : False Positive Rate (1 - specificity)
#                  - TP       : True Positive Rate (sensitivity)
#                  - FollowUp : Follow-up time interval (e.g., in months)
#                  - AUC      : Area Under the Curve for the respective follow-up time.
#   model_name : (Optional) A character string labeling the model. If provided,
#                the label is appended to the plot title.
#
# Returns:
#   A ggplot object representing the time-dependent ROC curves.
#
# Notes:
#   - The function uses a colorblind-friendly palette ("Set1") for the curve colors.
#   - The plot is fixed to the coordinate range [0, 1] on both axes.
#------------------------------------------------------------------------------#
plot_tvROC <- function(roc_data, model_name = NULL, col_set = "Set1") {
  # Append a prefix to the model name if provided.
  if (!is.null(model_name)) {
    model_name <- paste("for", model_name)
  }

  # Prepare the ROC data by creating a label that combines follow-up time and AUC.
  # Then, plot TP vs. FP with each follow-up interval displayed in a different color.
  roc_plot <- roc_data %>%
    mutate(roc_label = paste("Follow-up:", FollowUp, "months\nAUC:", AUC)) %>%
    ggplot(aes(x = FP, y = TP, color = roc_label)) +
    geom_line(linewidth = 1) +
    scale_color_brewer(palette = col_set) +
    labs(
      title = paste("Time-Dependent ROC Curves", model_name),
      x = "1 - Specificity",
      y = "Sensitivity",
      color = ""
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dotted",
      color = "black",
      linewidth = 1
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal(base_size = 14) +
    theme(
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.position = "bottom",
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      plot.title = element_text(size = 16, hjust = 0.5)
    )

  return(roc_plot)
}

# -------------------------------------------------------------------------- #
## Calibration Evaluation -----------------------------------------------------
# -------------------------------------------------------------------------- #

# Function to calculate calibration metrics for a joint model.
# It computes calibration plots and metrics (such as ICI) for each specified
# follow-up time and returns combined data frames for plotting.
calc_cal_metrics <- function(model, newdata, Tstart, follow_up_times) {
  # Validate that the model is either a joint model (jm) or a Cox model (coxph)
  if (!inherits(model, "jm") && !inherits(model, "coxph")) {
    stop("Model must be of class 'jm' or 'coxph'.")
  }
  calibration_data <- list() # Initialize list to store calibration data

  # Loop through each follow-up time (Dt)
  for (Dt in follow_up_times) {
    if (inherits(model, "jm")) {
      # Get calibration plot data without displaying the plot
      cal <- JMbayes2::calibration_plot(
        model,
        newdata = newdata,
        Tstart = Tstart,
        Dt = Dt,
        plot = FALSE
      )
      # Calculate calibration metrics (e.g., Integrated Calibration Index)
      ici <- calibration_metrics(model, newdata, Tstart = Tstart, Dt = Dt)
      # Store observed vs. predicted data along with the calibration metric
      cal_df <- tibble(
        observed = cal$observed,
        predicted = cal$predicted,
        follow_up = Dt,
        ici = ici[1],
        model_label = paste("Follow-up:", Dt, "months\nICI:", round(ici, 2))
      )

      density_data <- tibble(
        density = density(cal$pi_u_t)$y / max(density(cal$pi_u_t)$y),
        preds = density(cal$pi_u_t)$x,
        ici = ici[1],
        model_label = paste("Follow-up:", Dt, "months\nICI:", round(ici, 2))
      )
    } else {
      # Must be a Cox model (coxph)
      cox_df <- create_cox_df(newdata, Tstart)
      newd <- cox_df
      newd$fup_time <- Dt
      newd$mortality_status <- 1
      p <- 1 - predict(model, type = "survival", newdata = newd)
      y <- with(cox_df, Surv(fup_time, mortality_status))
      cal <- pmcalibration(
        y = y,
        p = p,
        smooth = "rcs",
        nk = 5,
        ci = "pw",
        time = Dt
      )
      # Store observed vs. predicted data along with the calibration metric
      cal_df <- tibble(
        observed = get_cc(cal)$p_c,
        predicted = get_cc(cal)$p,
        follow_up = Dt,
        ici = cal$metrics[1],
        model_label = paste("Follow-up:", Dt, "months\nEavg:", round(ici, 2))
      )

      density_data <- tibble(
        density = density(p)$y / max(density(p)$y),
        preds = density(p)$x,
        ici = cal$metrics[1],
        model_label = paste("Follow-up:", Dt, "months\nEavg:", round(ici, 2))
      )
    }

    # Save both calibration and density data for this follow-up time
    calibration_data[[Dt]] <- list(cal_df = cal_df, density_data = density_data)
  }

  # Combine calibration data for all follow-up times into single data frames
  combined_cal_df <- bind_rows(lapply(calibration_data, `[[`, "cal_df"))
  combined_density_df <- bind_rows(lapply(
    calibration_data,
    `[[`,
    "density_data"
  ))

  # Return a list containing the combined calibration and density data frames
  list(cal_df = combined_cal_df, density_df = combined_density_df)
}

# Function to plot the calibration curves and density curves.
# Inputs:
#   - calibration_data: A list containing calibration and density data frames.
#   - model_name: (Optional) A label for the model.
plot_cal <- function(calibration_data, model_name = NULL) {
  # Append model name if provided
  if (!is.null(model_name)) {
    model_name <- paste("for", model_name)
  }

  # Extract the combined calibration and density data frames
  combined_cal_df <- calibration_data$cal_df
  combined_density_df <- calibration_data$density_df

  # Generate the calibration plot using ggplot2
  plot <- ggplot() +
    # Plot calibration curves (observed vs predicted)
    geom_line(
      data = combined_cal_df,
      aes(x = predicted, y = observed, color = model_label),
      linewidth = 1.2
    ) +
    # Add a diagonal reference line representing perfect calibration
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dotted",
      color = "black",
      linewidth = 1
    ) +
    # Overlay density curves (dashed lines) of predicted probabilities
    geom_line(
      data = combined_density_df,
      aes(x = preds, y = density, color = model_label),
      linetype = "dashed",
      linewidth = 1
    ) +
    # Configure the primary y-axis and a secondary axis for density
    scale_y_continuous(
      name = "Observed Probability",
      sec.axis = sec_axis(~., name = "Density")
    ) +
    scale_color_brewer(palette = "Set1") + # Use a colorblind-friendly palette
    labs(
      title = paste("Calibration Plot", model_name),
      x = "Predicted Probability",
      color = ""
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal(base_size = 14) +
    theme(
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.position = "bottom",
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      plot.title = element_text(size = 16, hjust = 0.5)
    )

  return(plot)
}

# ============================================================================ #
# Dynamic Predictions Plot ----------------------------------------------------
# ============================================================================ #

# Function to generate a dynamic predictions plot for a specific subject.
# This plot combines longitudinal predictions (e.g., 6MWT distance) and
# event predictions (e.g., mortality risk) over time.
#
# Inputs:
#   - joint_model: The fitted joint model.
#   - data: The dataset containing longitudinal measurements.
#   - id: The subject identifier for whom to generate predictions.
#   - t0: The landmark time for predictions.
#   - fu_t: Follow-up time (default is 60 months).
#   - Various graphical parameters for customizing the appearance.
#
# The function returns a ggplot object displaying:
#   - The predicted longitudinal trajectory with confidence intervals.
#   - Observed data points.
#   - The predicted event probability scaled to match the longitudinal outcome.
plot_dyn_pred <- function(
  joint_model,
  data,
  id,
  t0,
  fu_t = 60,
  fill_CI_long = "#0000FF4D",
  fill_CI_event = "#FF00004D",
  col_line_long = "#0000FF",
  col_line_event = "#FF0000",
  col_points = "blue",
  cex_points = 2,
  lwd_long = 1,
  lwd_event = 1,
  main = NULL,
  xlim = NULL,
  add_labs = TRUE,
  x_months = 6,
  x_lab = "Time (months)",
  y_lab = "6 minute walk distance (meters)",
  add_sec_y = TRUE,
  sec_y_lab = "Event Probability"
) {
  # -------------------------------------------------------------------------- #
  # Prepare the data for the specific subject and landmark time (t0)
  # -------------------------------------------------------------------------- #
  ND <- data %>%
    filter(id == !!id, time < t0) %>% # Select records for the given subject before t0
    mutate(
      mortality_status = 0, # Set mortality status to 0 for prediction purposes
      fup_time = t0 # Define the follow-up time as the landmark time t0
    )

  # -------------------------------------------------------------------------- #
  # Generate Longitudinal Predictions -----------------------------------------
  # -------------------------------------------------------------------------- #
  predLong <- predict(
    joint_model,
    newdata = ND,
    times = seq(t0, fu_t, length.out = 51), # Predict from t0 up to fu_t
    type = "subject_specific",
    return_newdata = TRUE
  )
  # Combine the longitudinal predictions from different responses into one dataframe
  pred_Long <- bind_rows(predLong[[1L]], predLong[[2L]])

  # -------------------------------------------------------------------------- #
  # Generate Event (Survival) Predictions -------------------------------------
  # -------------------------------------------------------------------------- #
  pred_Event <- predict(
    joint_model,
    newdata = ND,
    process = "event",
    times = seq(t0, fu_t, length.out = 51), # Predict event probabilities over time
    return_newdata = TRUE
  )

  # -------------------------------------------------------------------------- #
  # Extract Prediction Attributes for Processing -----------------------------
  # -------------------------------------------------------------------------- #
  id_var <- attr(predLong, "id_var")
  time_var <- attr(predLong, "time_var")
  Time_var <- attr(predLong, "Time_var")
  resp_vars <- attr(predLong, "resp_vars")
  ranges <- attr(predLong, "ranges")
  last_times <- attr(predLong, "last_times")
  y <- attr(predLong, "y")
  times_y <- attr(predLong, "times_y")

  # -------------------------------------------------------------------------- #
  # Process Longitudinal Predictions ------------------------------------------
  # -------------------------------------------------------------------------- #
  # Filter predictions to include only time points up to the last available measurement
  pred_Long <- pred_Long %>%
    filter(!!sym(time_var) <= as.numeric(last_times)) %>%
    rename(times_long = !!sym(time_var))

  # Identify the columns that contain the predicted values and their confidence intervals
  outcome_inx <- grep("pred_", names(pred_Long), fixed = TRUE)

  # Create a dataframe with the predicted values and their lower and upper confidence bounds
  long_data <- pred_Long %>%
    transmute(
      times_long,
      preds_long = pred_Long[[outcome_inx]],
      low_long = pred_Long[[outcome_inx + 1]],
      upp_long = pred_Long[[outcome_inx + 2]]
    ) %>%
    filter(!is.na(preds_long))

  # Prepare the observed longitudinal data for plotting
  observed_data <- data.frame(
    x = unlist(times_y),
    y = unlist(y)
  ) %>%
    filter(!is.na(y))

  # -------------------------------------------------------------------------- #
  # Process Event Predictions -------------------------------------------------
  # -------------------------------------------------------------------------- #
  # A scaling factor is applied to event predictions to align their range
  # with the longitudinal outcome (e.g., 6MWT distance).
  y_axis_factor <- 800 # Factor to scale event probabilities for plotting

  event_data <- data.frame(
    times_event = pred_Event[[Time_var]],
    preds_event = pred_Event[[grep("pred_", names(pred_Event))]],
    low_event = pred_Event[[grep("pred_", names(pred_Event)) + 1]],
    upp_event = pred_Event[[grep("pred_", names(pred_Event)) + 2]]
  ) %>%
    filter(!is.na(preds_event)) %>%
    mutate_at(vars(preds_event, low_event, upp_event), ~ . * y_axis_factor)

  # Define the y-axis limits based on the scaling factor
  y_lim <- c(0, 1.01) * y_axis_factor

  # Set default plot title if not provided
  if (is.null(main)) {
    main <- paste("for", round(t0), "months follow-up")
  }

  # -------------------------------------------------------------------------- #
  # Generate the Dynamic Predictions Plot -------------------------------------
  # -------------------------------------------------------------------------- #
  p <- ggplot() +
    # Plot longitudinal predictions with a confidence interval ribbon
    geom_ribbon(
      data = long_data,
      aes(x = times_long, ymin = low_long, ymax = upp_long),
      fill = fill_CI_long,
      alpha = 0.3
    ) +
    geom_line(
      data = long_data,
      aes(x = times_long, y = preds_long, color = "6MWT distance"),
      linewidth = lwd_long
    ) +
    # Plot the observed longitudinal measurements
    geom_point(
      data = observed_data,
      aes(x = x, y = y),
      color = col_points,
      size = cex_points
    ) +
    # Draw a vertical dashed line at the last observed time point
    geom_vline(
      xintercept = last_times + 0.01,
      linetype = "dashed",
      color = "gray"
    ) +
    # Plot event predictions with a confidence interval ribbon
    geom_ribbon(
      data = event_data,
      aes(x = times_event, ymin = low_event, ymax = upp_event),
      fill = fill_CI_event,
      alpha = 0.3
    ) +
    geom_line(
      data = event_data,
      aes(x = times_event, y = preds_event, color = "Event probability"),
      linewidth = lwd_event
    ) +
    # Configure x-axis breaks
    scale_x_continuous(
      breaks = seq(0, max(event_data$times_event), by = x_months)
    ) +
    scale_color_manual(
      values = c(
        "6MWT distance" = col_line_long,
        "Event probability" = col_line_event
      )
    ) +
    coord_cartesian(ylim = y_lim) +
    # Add labels and title to the plot
    labs(
      x = "",
      y = "",
      color = ""
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.position = "bottom",
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      plot.title = element_text(size = 16, hjust = 0.5),
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )

  if (add_labs) {
    p <- p +
      labs(
        title = paste("Dynamic Predictions Plot", main),
        x = x_lab,
        y = y_lab
      )
  }

  if (add_sec_y) {
    p <- p +
      # Configure the y-axis to include a secondary axis for event probability
      scale_y_continuous(
        sec.axis = sec_axis(
          ~ . / y_axis_factor,
          name = sec_y_lab,
          breaks = seq(0, 1, by = 0.2)
        )
      )
  }

  return(p)
}


# Others ####

bootstrap_iteration_cox <- function(iteration) {
  library(JMbayes2)
  # Resample subject IDs with replacement.
  boot_ids <- sample(unique_ids, size = length(unique_ids), replace = TRUE)

  # For each sampled subject, extract the corresponding rows and assign a new sequential ID.
  boot_data <- bind_rows(lapply(seq_along(boot_ids), function(j) {
    long_df_test %>%
      filter(id == boot_ids[j]) %>%
      mutate(id = j)
  }))

  boot_cox <- boot_data %>%
    arrange(id, time) %>% # Sort by patient ID and time
    group_by(id) %>%
    filter(row_number() == 1) %>% # Keep only the first visit (baseline)
    ungroup()

  boot_joint <- boot_data %>%
    filter(time < 12) %>% # Select records for the given subject before t0
    mutate(
      mortality_status = 0, # Set mortality status to 0 for prediction purposes
      fup_time = 12 # Define the follow-up time as the landmark time t0
    )

  # Compute AUC measures using tvAUC.
  auc_36 <- tvAUC(
    object = jointFit_train,
    newdata = boot_data,
    Tstart = 0,
    Dt = 36
  )
  auc_60 <- tvAUC(
    object = jointFit_train,
    newdata = boot_data,
    Tstart = 0,
    Dt = 60
  )
  auc_1_36 <- tvAUC(
    object = jointFit_train,
    newdata = boot_data,
    Tstart = 12,
    Dt = 24
  )
  auc_1_60 <- tvAUC(
    object = jointFit_train,
    newdata = boot_data,
    Tstart = 12,
    Dt = 48
  )

  cox_preds <- survfit(CoxFit_six_train, newdata = boot_cox)
  joint_preds <- predict(
    jointFit_train,
    newdata = boot_joint,
    process = "event",
    times = 36,
    return_newdata = TRUE
  )

  cox_pred <- cox_preds$cumhaz[cox_preds$time == 36]
  joint_pred <- joint_preds[joint_preds$time == 36, "pred_CIF"]

  timeROC_cox <- timeROC::timeROC(
    T = boot_cox$fup_time,
    delta = boot_cox$mortality_status,
    marker = cox_pred,
    cause = 1,
    weighting = "cox",
    times = 36
  )

  timeROC_joint <- timeROC::timeROC(
    T = boot_cox$fup_time,
    delta = boot_cox$mortality_status,
    marker = joint_pred,
    cause = 1,
    weighting = "cox",
    times = 36
  )

  survivalROC_cox <- survivalROC::survivalROC(
    Stime = boot_cox$fup_time,
    status = boot_cox$mortality_status,
    marker = cox_pred,
    predict.time = 36,
    method = "KM"
  )

  survivalROC_joint <- survivalROC::survivalROC(
    Stime = boot_cox$fup_time,
    status = boot_cox$mortality_status,
    marker = joint_pred,
    predict.time = 36,
    method = "KM"
  )

  risksetROC_cox <- risksetROC::risksetAUC(
    Stime = boot_cox$fup_time,
    status = boot_cox$mortality_status,
    marker = cox_pred,
    tmax = 36
  )

  risksetROC_joint <- risksetROC::risksetAUC(
    Stime = boot_cox$fup_time,
    status = boot_cox$mortality_status,
    marker = joint_pred,
    tmax = 36
  )

  # Return the AUC results as a tibble.
  tibble(
    auc_36 = auc_36$auc,
    auc_60 = auc_60$auc,
    auc_1_36 = auc_1_36$auc,
    auc_1_60 = auc_1_60$auc,
    timeROC_cox = timeROC_cox$AUC[2],
    timeROC_joint = timeROC_joint$AUC[2],
    survivalROC_cox = survivalROC_cox$AUC,
    survivalROC_joint = survivalROC_joint$AUC,
    risksetROC_cox = risksetROC_cox$Cindex,
    risksetROC_joint = risksetROC_joint$Cindex
  )
}
