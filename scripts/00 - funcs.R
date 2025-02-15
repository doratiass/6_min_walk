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
library(tidyverse)      # For data manipulation and visualization
library(gtsummary)      # For summary statistics and table generation
library(ggpubr)         # For publication ready plots
library(gridExtra)      # For arranging multiple plots
library(flextable)      # For creating formatted tables
library(survminer)      # For survival analysis plots
library(survival)       # For survival analysis 
library(JMbayes2)       # For joint modeling and tvAUC computations
library(pmcalibration)  # For calibration plots
library(survivalROC)    # For time-dependent ROC curves
# ============================================================================ #
# Definitions  ----------------------------------------------------------------
# ============================================================================ #
# Define scaling factors and follow-up time limits
age_scale     <- 1    # Scaling for age (set to 1 for no scaling)
six_min_scale <- 1    # Scaling for 6-minute walk distance (set to 1 for no scaling)
min_fup_time  <- 6    # Minimum follow-up time in months
max_fup_time  <- 120  # Maximum follow-up time in months
tests_time <- 12      # Time limit for 6-minute walk test measurements

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
  "fup_time" = "Follow-up time",
  "nyha" = "NYHA class",
  "ischemic_etiology" = "Ischemic etiology",
  "sbp" = "Systolic blood pressure",
  "fusid" = "Use furosemide",
  "torsemide" = "Use torsemide",
  "bumetanide" = "Use bumetanide",
  "metolazone" = "Use metolazone",
  "ht" = "Height",
  "wt" = "Weight",
  "statin" = "Use statins",
  "ace_inhibitor" = "Use ACE inhibitors",
  "beta_blocker" = "Use beta blockers",
  "arbs" = "Use ARBs",
  "diuretic" = "Use diuretics",
  "k_spare" = "Use potassium-sparing diuretics",
  "thiazide" = "Use thiazide diuretics",
  "allopurinol" = "Use allopurinol",
  "na" = "Sodium",
  "hgb" = "Hemoglobin",
  "ua" = "Uric acid",
  "x6mw_dist_meter" = "6MWT distance (50m)",
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
  ifelse(x %in% vars_dict$var,
         vars_dict[vars_dict$var == x, "name", drop = TRUE],
         x)
}

# Apply the label_get function to a vector of variable names.
vars_label <- function(x) {
  sapply(x, label_get, USE.NAMES = FALSE)
}

# Function to retrieve the original variable name from a given label.
var_get <- function(x) {
  ifelse(x %in% vars_dict$name,
         vars_dict[vars_dict$name == x, "var", drop = TRUE],
         x)
}

# ============================================================================ #
# Model Evaluation ------------------------------------------------------------
# ============================================================================ #

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
      filter(time <= Tstart,
             fup_time > Tstart) %>%
      arrange(id, time) %>%
      group_by(id) %>%
      filter(time == max(time)) %>% 
      ungroup() %>%
      mutate(fup_time = fup_time - Tstart,
             time = 0)
  } else {
    data <- data %>%
      arrange(id, time) %>%
      group_by(id) %>%
      slice(1) %>%
      ungroup()
  }
  return(data)
}


fit_models <- function(data_list) {
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
    model = TRUE, x = TRUE, y = TRUE
  )
  
  # 3. Combine the two models into a joint model.
  jointFit_train <- jm(
    CoxFit_train, lme_6min_train, time_var = "time",
    functional_forms = ~ value(x6mw_dist_meter) + slope(x6mw_dist_meter)
  )
  
  # 4. Fit a Cox model that also includes the 6-min walk distance for comparison.
  CoxFit_six_train <- coxph(
    Surv(fup_time, mortality_status) ~ gender + age + nyha + x6mw_dist_meter, 
    data = cox_df_train, 
    model = TRUE, x = TRUE, y = TRUE
  )
  
  # Return a list containing all the models.
  return(list(
    joint_model    = jointFit_train,
    cox_model  = CoxFit_six_train
  ))
}
# -------------------------------------------------------------------------- #
## Discrimination Evaluation --------------------------------------------------
# -------------------------------------------------------------------------- #

# Function to create ROC data for different follow-up times using a joint model.
# Inputs:
#   - model: The fitted model. can be either joint model or cox model
#   - long_data: Longitudinal data for prediction.
#   - Tstart: The starting time for prediction.
#   - follow_up_times: Vector of follow-up times (Dt) at which to calculate ROC.
# Returns a tibble containing false positive rates, true positive rates, follow-up time, and AUC.
create_roc_data_2 <- function(model, long_data, Tstart, follow_up_times) {
  roc_list <- list()  # Initialize list to store ROC results
  
  # Loop through each specified follow-up time (Dt)
  for (Dt in follow_up_times) {
    # Compute time-dependent ROC using tvROC function
    if (inherits(model, "jm")) {
      roc_result <- tvROC(model, newdata = long_data, Tstart = Tstart, Dt = Dt, 
                          type_weights = "IPCW")
    } else {
      roc_result <- tvAUC_coxph(model, newdata = long_data, Tstart = Tstart, Dt = Dt)
    }
    
    # Calculate and round the Area Under the Curve (AUC)
    auc <- round(tvAUC(roc_result)$auc, 3)
    # Store both ROC result and AUC in the list (using Dt as the key)
    roc_list[[as.character(Dt)]] <- list(roc = roc_result, auc = auc)
  }
  
  # Combine all ROC data into a single tibble for plotting
  roc_data <- do.call(bind_rows, lapply(names(roc_list), function(Dt) {
    roc <- roc_list[[Dt]]$roc
    auc <- roc_list[[Dt]]$auc
    data.frame(
      FP = roc$FP,             # False positive rate (1 - specificity)
      TP = roc$TP,             # True positive rate (sensitivity)
      FollowUp = as.numeric(Dt),
      AUC = auc                # Area Under the Curve for this follow-up time
    )
  }))
  
  return(roc_data)
}

# -----------------------------------------------------------------------------#
# Time-Dependent AUC Computation for Dynamic Survival Predictions via Cox Model
#
# Overview:
#   This function computes the time-dependent Area Under the Curve (AUC) for dynamic
#   survival predictions based on a Cox proportional hazards model. It evaluates the 
#   discriminative performance of the model at a specified prediction horizon by 
#   comparing the model-generated risk scores with the observed survival outcomes.
#
# Parameters:
#   object:
#     A fitted Cox proportional hazards model (an object of class "coxph"). This model 
#     should be obtained using the 'coxph' function.
#
#   newdata:
#     A data frame containing the prediction data. It must include the following columns:
#       - id: Unique identifier for each subject.
#       - time: The time at which the covariates were recorded (typically at baseline).
#       - fup_time: The follow-up time corresponding to the occurrence of the event or censoring.
#       - mortality_status: The event indicator (1 if the event occurred, 0 if censored).
#
#   Tstart:
#     The starting time for prediction. Only subjects who are at risk at Tstart are included 
#     in the analysis. (Default is 0.)
#
#   Thoriz:
#     The prediction horizon time point. It must be greater than Tstart. If Thoriz is not provided,
#     then the time interval Dt must be specified so that Thoriz can be computed as Tstart + Dt.
#
#   Dt:
#     The duration (time interval) for prediction. This parameter is used to compute Thoriz when 
#     it is not explicitly provided.
#
#   id_var, time_var, Time_var, event_var:
#     Optional arguments to specify the column names in 'newdata' for the subject identifier,
#     time of covariate measurement, follow-up time, and event indicator, respectively.
#     Default values are "id", "time", "fup_time", and "mortality_status".
#
#   type:
#     A character value specifying the type of output. If "auc", the function returns the computed 
#     AUC along with key prediction time points and model metadata. If "roc", the output includes 
#     additional ROC-related metrics such as true positive and false positive rates, threshold values, 
#     F1 score, and Youden's index.
#
#   ...:
#     Additional arguments that are passed to the underlying prediction functions.
#
# Returns:
#   A list containing the following elements:
#     - auc: The computed time-dependent AUC (numeric).
#     - Tstart: The adjusted starting time used for prediction (numeric).
#     - Thoriz: The adjusted prediction horizon time (numeric).
#     - nr: The number of unique subjects included in the computation (integer).
#     - classObject: The class of the input Cox model.
#     - nameObject: The name of the input Cox model object.
#
#   For type "roc", additional elements include:
#     - TP: The true positive rates.
#     - FP: The false positive rates.
#     - thrs: The threshold values used in the ROC analysis.
#     - F1score: The median threshold corresponding to the maximum F1 score.
#     - Youden: The median threshold corresponding to the maximum Youden's index.
#
# Notes:
#   - A small perturbation (1e-06) is added to Tstart and Thoriz to mitigate potential numerical issues.
#   - For dynamic prediction, only the most recent covariate measurements for subjects at risk at Tstart are used.
#
# -----------------------------------------------------------------------------#
tvAUC_coxph <- function(object, newdata, Tstart = 0, Thoriz = NULL, Dt = NULL, 
                        id_var = "id", time_var = "time", Time_var = "fup_time",
                        event_var = "mortality_status", type = c("auc", "roc"), ...) {
  # Confirm that 'object' is of class "coxph".
  if (!inherits(object, "coxph"))
    stop("Use only with 'coxph' objects.\n")
  
  # Ensure that 'newdata' is a non-empty data frame containing the necessary prediction data.
  if (!is.data.frame(newdata) || nrow(newdata) == 0)
    stop("'newdata' must be a data.frame with at least one row.\n")
  
  # Verify that at least one of 'Thoriz' (prediction horizon) or 'Dt' (time interval) is provided.
  if (is.null(Thoriz) && is.null(Dt))
    stop("Either 'Thoriz' or 'Dt' must be non-null.\n")
  
  # If provided, 'Thoriz' must be greater than 'Tstart'; otherwise, the function stops.
  if (!is.null(Thoriz) && Thoriz <= Tstart)
    stop("'Thoriz' must be larger than 'Tstart'.\n")
  
  # Compute 'Thoriz' using 'Dt' if Thoriz is not provided.
  if (is.null(Thoriz))
    Thoriz <- Tstart + Dt
  
  # Calculate 'Dt' as the difference between Thoriz and Tstart if Dt is not provided.
  if (is.null(Dt))
    Dt <- Thoriz - Tstart
  
  # Confirm that 'newdata' contains the necessary columns as specified by the parameter names.
  if (is.null(newdata[[id_var]]))
    stop("cannot find the '", id_var, "' variable in newdata.")
  if (is.null(newdata[[time_var]]))
    stop("cannot find the '", time_var, "' variable in newdata.")
  if (is.null(newdata[[Time_var]]))
    stop("cannot find the '", Time_var, "' variable in newdata.")
  if (is.null(newdata[[event_var]]))
    stop("cannot find the '", event_var, "' variable in newdata.")
  
  # Adjust Tstart and Thoriz slightly to avoid potential numerical precision issues.
  Tstart <- Tstart + 1e-06
  Thoriz <- Thoriz + 1e-06
  
  # Filter the dataset to include only subjects who are at risk at Tstart.
  # For subjects with multiple records, retain the most recent measurement prior to Tstart.
  newdata <- create_cox_df(newdata, Tstart)

  # Compute the linear predictor (risk score) for each subject using the fitted Cox model.
  lp <- predict(object, newdata = newdata, type = "lp")
  
  # Evaluate the time-dependent AUC using the survivalROC function based on the computed risk scores.
  auc <- survivalROC(Stime = newdata[[Time_var]],
                     status = newdata[[event_var]],
                     marker = lp,
                     predict.time = Thoriz,
                     method = "KM")
  
  # Prepare the output list, which varies based on the selected 'type' parameter.
  if (type == "auc") {
    # For 'auc' type: return the computed AUC along with key time points, subject count, and model metadata.
    out <- list(auc = auc$AUC, Tstart = Tstart, Thoriz = Thoriz,
                nr = length(unique(newdata[[id_var]])),
                classObject = class(object),
                nameObject = deparse(substitute(object)))
    class(out) <- "tvAUC_coxph"
  } else {
    # For 'roc' type: generate ROC metrics and derive additional performance measures.
    metric_df <- tibble(
      TP = auc$TP,
      FP = auc$FP,
      thres = auc$cut.values
    ) %>%
      mutate(
        precision = TP / (TP + FP),
        # Here, TP is assumed to represent sensitivity (i.e., recall)
        recall = TP,
        F1 = 2 * (precision * recall) / (precision + recall),
        Youden = TP - FP  # As specificity = 1 - FP, Youden = TP + (1-FP) - 1
      )
    # Determine the median threshold corresponding to the maximum F1 score.
    F1score <- median(metric_df[metric_df$F1 == max(metric_df$F1, na.rm = TRUE),]$thres, na.rm = TRUE)
    # Determine the median threshold corresponding to the maximum Youden's index.
    Youden <- median(metric_df[metric_df$Youden == max(metric_df$Youden, na.rm = TRUE),]$thres, na.rm = TRUE)
    out <- list(TP = metric_df$TP, FP = metric_df$FP, auc = auc$AUC, 
                thrs = metric_df$thres, F1score = F1score, Youden = Youden,
                Tstart = Tstart, Thoriz = Thoriz, nr = length(unique(newdata[[id_var]])),
                classObject = class(object),
                nameObject = deparse(substitute(object)))
    class(out) <- "tvROC_coxph"
  }
  return(out)
}

tvAUC_coxph_old <- function(object, newdata, Tstart = 0, Thoriz = NULL, Dt = NULL, ...) {
  # Check that a Cox model was supplied
  if (!inherits(object, "coxph"))
    stop("Use only with 'coxph' objects.\n")
  
  if (!is.data.frame(newdata) || nrow(newdata) == 0)
    stop("'newdata' must be a data.frame with at least one row.\n")
  
  if (is.null(Thoriz) && is.null(Dt))
    stop("Either 'Thoriz' or 'Dt' must be non-null.\n")
  
  if (!is.null(Thoriz) && Thoriz <= Tstart)
    stop("'Thoriz' must be larger than 'Tstart'.\n")
  
  if (is.null(Thoriz))
    Thoriz <- Tstart + Dt
  
  # For the Cox model we assume that newdata contains the following variables:
  #   id      : subject identifier
  #   time    : time at which covariates are recorded (often 0 for baseline)
  #   time    : the eventual event (or censoring) time (here the same name)
  #   event   : event indicator (1 = event, 0 = censoring)
  # (If your data use different names, you may change the defaults below.)
  id_var    <- "id"
  time_var  <- "time"    # covariate update (or baseline) time
  Time_var  <- "fup_time"    # event time
  event_var <- "mortality_status"
  
  # Check that the required columns exist in newdata:
  if (is.null(newdata[[id_var]]))
    stop("cannot find the '", id_var, "' variable in newdata.")
  if (is.null(newdata[[time_var]]))
    stop("cannot find the '", time_var, "' variable in newdata.")
  if (is.null(newdata[[Time_var]]))
    stop("cannot find the '", Time_var, "' variable in newdata.")
  if (is.null(newdata[[event_var]]))
    stop("cannot find the '", event_var, "' variable in newdata.")
  
  # Slightly perturb the time points (as in the original code)
  Tstart <- Tstart + 1e-06
  Thoriz <- Thoriz + 1e-06
  
  # For dynamic prediction we keep only subjects who are at risk at Tstart.
  # (Assuming one row per subject, newdata[[time_var]] is when the covariates were measured.)
  # newdata <- newdata[order(newdata[[Time_var]]), ]
  # newdata <- newdata[newdata[[Time_var]] > Tstart, ]
  # newdata <- newdata[newdata[[time_var]] <= Tstart, ]
  newdata <- newdata %>%
    filter(time <= Tstart,
           fup_time > Tstart) %>%
    arrange(id, time) %>%
    group_by(id) %>%
    slice(1) %>%
    ungroup()
  
  # Coerce id to a factor (so later tapply() calls work as expected)
  newdata[[id_var]] <- as.factor(newdata[[id_var]])
  
  # Check that at least one subject experienced an event in [Tstart, Thoriz)
  test1 <- newdata[[Time_var]] < Thoriz & newdata[[event_var]] == 1
  if (!any(test1))
    stop("It appears that there are no events in the interval [Tstart, Thoriz).")
  
  # Create a copy in which we “update” the covariate time to Tstart.
  newdata2 <- newdata
  newdata2[[Time_var]] <- Thoriz
  newdata2[[event_var]] <- 0  # event indicator is not used in prediction
  
  # Define a helper function to compute the predicted event probability
  # from Tstart to Thoriz for each subject.
  predict_event_cox <- function(model, data, Thoriz, ...) {
    data[[Time_var]] <- Thoriz
    pred_event <- 1-predict(model, data, type = "survival")
    return(pred_event)
  }
  
  # Compute the predicted event probability for each subject using newdata2.
  # Then the “dynamic survival probability” is one minus that.
  pred_event <- predict_event_cox(object, newdata2, Thoriz, ...)
  si_u_t <- 1-pred_event
  names(si_u_t) <- newdata2[[id_var]]
  
  # For each subject we extract the “last” recorded event time and status.
  id    <- newdata[[id_var]]
  Time  <- newdata[[Time_var]]
  event <- newdata[[event_var]]
  f <- factor(id, levels = unique(id))
  Time  <- tapply(Time,  f, tail, 1L)
  event <- tapply(event, f, tail, 1L)
  names(Time) <- names(event) <- as.character(unique(id))
  if (any(dupl <- duplicated(Time))) {
    Time[dupl] <- Time[dupl] + runif(sum(dupl), 1e-07, 1e-06)
  }
  if (!all(names(si_u_t) == names(Time)))
    stop("Mismatch between subject names in predictions and in newdata.")
  
  # Now compute the time-dependent AUC by comparing all pairs of subjects.
  # (The following block follows the structure of the original code.)
  auc <- if (length(Time) > 1L) {
    pairs <- combn(as.character(unique(id)), 2)
    Ti <- Time[pairs[1, ]]
    Tj <- Time[pairs[2, ]]
    di <- event[pairs[1, ]]
    dj <- event[pairs[2, ]]
    si_u_t_i <- si_u_t[pairs[1, ]]
    si_u_t_j <- si_u_t[pairs[2, ]]
    
    # Define indicators for the various comparisons.
    ind1 <- (Ti <= Thoriz & di == 1) & (Tj > Thoriz)
    ind2 <- (Ti <= Thoriz & di == 0) & (Tj > Thoriz)
    ind3 <- (Ti <= Thoriz & di == 1) & (Tj <= Thoriz & dj == 0)
    ind4 <- (Ti <= Thoriz & di == 0) & (Tj <= Thoriz & dj == 0)
    names(ind1) <- names(ind2) <- names(ind3) <- names(ind4) <-
      paste(names(Ti), names(Tj), sep = "_")
    ind <- ind1 | ind2 | ind3 | ind4
    
    # For pairs where one of the subjects is censored we use inverse‐probability weighting.
    # if (any(ind2)) {
    #   nams <- strsplit(names(ind2[ind2]), "_")
    #   nams_i <- sapply(nams, "[", 1)
    #   unq_nams_i <- unique(nams_i)
    #   pi_u_t <- predict_event_cox(object,
    #                               newdata[newdata[[id_var]] %in% unq_nams_i, ],
    #                               Thoriz, ...)
    #   f_pred <- factor(newdata[newdata[[id_var]] %in% unq_nams_i, ][[id_var]],
    #                    levels = unique(newdata[[id_var]]))
    #   names(pi_u_t) <- f_pred
    #   pi_u_t <- tapply(pi_u_t, f_pred, tail, 1)
    #   ind[ind2] <- ind[ind2] * pi_u_t[nams_i]
    # }
    # if (any(ind3)) {
    #   nams <- strsplit(names(ind3[ind3]), "_")
    #   nams_j <- sapply(nams, "[", 2)
    #   unq_nams_j <- unique(nams_j)
    #   qi_u_t <- predict_event_cox(object,
    #                               newdata[newdata[[id_var]] %in% unq_nams_j, ],
    #                               Thoriz, ...)
    #   f_pred <- factor(newdata[newdata[[id_var]] %in% unq_nams_j, ][[id_var]],
    #                    levels = unique(newdata[[id_var]]))
    #   names(qi_u_t) <- f_pred
    #   qi_u_t <- 1 - tapply(qi_u_t, f_pred, tail, 1)
    #   ind[ind3] <- ind[ind3] * qi_u_t[nams_j]
    # }
    # if (any(ind4)) {
    #   nams <- strsplit(names(ind4[ind4]), "_")
    #   nams_i <- sapply(nams, "[", 1)
    #   nams_j <- sapply(nams, "[", 2)
    #   unq_nams_i <- unique(nams_i)
    #   unq_nams_j <- unique(nams_j)
    #   pi_u_t <- predict_event_cox(object,
    #                               newdata[newdata[[id_var]] %in% unq_nams_i, ],
    #                               Thoriz, ...)
    #   f_pred_i <- factor(newdata[newdata[[id_var]] %in% unq_nams_i, ][[id_var]],
    #                      levels = unique(newdata[[id_var]]))
    #   names(pi_u_t) <- f_pred_i
    #   pi_u_t <- tapply(pi_u_t, f_pred_i, tail, 1)
    #   qi_u_t <- predict_event_cox(object,
    #                               newdata[newdata[[id_var]] %in% unq_nams_j, ],
    #                               Thoriz, ...)
    #   f_pred_j <- factor(newdata[newdata[[id_var]] %in% unq_nams_j, ][[id_var]],
    #                      levels = unique(newdata[[id_var]]))
    #   names(qi_u_t) <- f_pred_j
    #   qi_u_t <- 1 - tapply(qi_u_t, f_pred_j, tail, 1)
    #   ind[ind4] <- ind[ind4] * pi_u_t[nams_i] * qi_u_t[nams_j]
    # }
    
    sum((si_u_t_i < si_u_t_j) * as.numeric(ind), na.rm = TRUE) / sum(ind, na.rm = TRUE)
  } else {
    NA
  }
  
  out <- list(auc = auc, Tstart = Tstart, Thoriz = Thoriz,
              nr = length(unique(id)),
              classObject = class(object),
              nameObject = deparse(substitute(object)))
  class(out) <- "tvAUC_coxph"
  return(out)
}

print.tvAUC_coxph <- function (x, digits = 4, ...) {
  if (!inherits(x, "tvAUC_coxph"))
    stop("Use only with 'tvAUC' objects.\n")
  if (x$class == "jm")
    cat("\n\tTime-dependent AUC for the Joint Model",  x$nameObject)
  else
    cat("\n\tTime-dependent AUC for the Cox Model",  x$nameObject)
  cat("\n\nEstimated AUC: ", round(x$auc, digits))
  #cat("\n\nEstimated AUC: ", round(x$auc, digits),
  #    " (95% CI: ", round(x$low_auc, digits), "-", round(x$upp_auc, digits),
  #    ")", sep = "")
  cat("\nAt time:", round(x$Thoriz, digits))
  cat("\nUsing information up to time: ", round(x$Tstart, digits),
      " (", x$nr, " subjects still at risk)", sep = "")
  cat("\n\n")
  invisible(x)
}

# -----------------------------------------------------------------------------#
# Bootstrap Iteration for Dynamic AUC Evaluation
#
# This function performs a single bootstrap iteration to assess the stability and 
# variability of dynamic AUC measures for survival models. It resamples subjects 
# (with replacement) from the test dataset, reassigns new sequential IDs, and 
# constructs bootstrap datasets. Dynamic AUC metrics are then computed for both 
# joint and Cox models over various prediction intervals.
#
# Parameters:
#   iteration:
#     The current bootstrap iteration number (primarily for tracking purposes).
#
#   test_df:
#     A data frame containing the test dataset. It must include an 'id' column that 
#     uniquely identifies each subject.
#
#   models:
#     A list containing the pre-fitted models required for AUC computation:
#       - joint_model: A joint model (from JMbayes2) for dynamic survival prediction.
#       - cox_model: A Cox proportional hazards model for dynamic survival prediction.
#
# Returns:
#   A tibble with computed AUC metrics for multiple prediction intervals:
#     - auc_joint_36: AUC from the joint model for a 36-month prediction interval (Tstart = 0, Dt = 36).
#     - auc_joint_60: AUC from the joint model for a 60-month prediction interval (Tstart = 0, Dt = 60).
#     - auc_cox_36:   AUC from the Cox model for a 36-month prediction interval (Tstart = 0, Dt = 36).
#     - auc_cox_60:   AUC from the Cox model for a 60-month prediction interval (Tstart = 0, Dt = 60).
#     - auc_joint_1_36: AUC from the joint model starting at 12 months with a 24-month interval (Tstart = 12, Dt = 24).
#     - auc_joint_1_60: AUC from the joint model starting at 12 months with a 48-month interval (Tstart = 12, Dt = 48).
#     - auc_cox_1_36:   AUC from the Cox model starting at 12 months with a 24-month interval (Tstart = 12, Dt = 24).
#     - auc_cox_1_60:   AUC from the Cox model starting at 12 months with a 48-month interval (Tstart = 12, Dt = 48).
#
# Dependencies:
#   This function requires the 'JMbayes2' and 'survivalROC' packages.
#
# -----------------------------------------------------------------------------#
bootstrap_iteration <- function(iteration, test_df, models) {
  library(JMbayes2)
  library(survivalROC)
  library(riskRegression)
  
  # -------------------------------#
  # Resample Subject IDs with Replacement
  # -------------------------------#
  # This section performs bootstrap resampling to generate a new dataset. Subject IDs are
  # sampled with replacement, and for each sampled subject, the corresponding rows are extracted.
  # New sequential IDs are assigned to ensure consistency in the bootstrap dataset.
  unique_ids <- unique(test_df$id)
  boot_ids <- sample(unique_ids, size = length(unique_ids), replace = TRUE)
  
  boot_data <- bind_rows(lapply(seq_along(boot_ids), function(j) {
    test_df %>%
      filter(id == boot_ids[j]) %>%
      mutate(id = j)
  })) 
  
  # -------------------------------#
  # Compute Dynamic AUC Measures
  # -------------------------------#
  # The dynamic AUC metrics are computed for different time intervals and starting points using 
  # both the joint and Cox models. The tvAUC and tvAUC_coxph functions are employed for these calculations.
  auc_joint_36   <- tvAUC(object = models$joint_model, newdata = boot_data, Tstart = 0, Dt = 36, 
                          type_weights = "IPCW")
  auc_joint_60   <- tvAUC(object = models$joint_model, newdata = boot_data, Tstart = 0, Dt = 60, 
                          type_weights = "IPCW")
  auc_cox_36 <- tvAUC_coxph(object = models$cox_model, newdata = boot_data, Tstart = 0, Dt = 36,
                            type = "auc")
  auc_cox_60 <- tvAUC_coxph(object = models$cox_model, newdata = boot_data, Tstart = 0, Dt = 60,
                            type = "auc")
  
  auc_joint_1_36 <- tvAUC(object = models$joint_model, newdata = boot_data, Tstart = 12, Dt = 24, 
                          type_weights = "IPCW")
  auc_joint_1_60 <- tvAUC(object = models$joint_model, newdata = boot_data, Tstart = 12, Dt = 48, 
                          type_weights = "IPCW")
  auc_cox_1_36 <- tvAUC_coxph(object = models$cox_model, newdata = boot_data, Tstart = 12, Dt = 24,
                              type = "auc")
  auc_cox_1_60 <- tvAUC_coxph(object = models$cox_model, newdata = boot_data, Tstart = 12, Dt = 48,
                              type = "auc")
  auc_score <- Score(object = list("Cox" = models$cox_model),  # Model(s) to evaluate
                     formula = Surv(fup_time, mortality_status) ~ gender + age + nyha + x6mw_dist_meter,  # Survival outcome
                     data = boot_data %>%
                       arrange(id, time) %>%
                       group_by(id) %>%
                       slice(1) %>%
                       ungroup(),  # Data
                     times = c(36,60),  # Time points for AUC calculation
                     cens.method = "ipcw",
                     cens.model = "km",  # IPCW using Kaplan-Meier estimator
                     metrics = "AUC")
  # -------------------------------#
  # Return AUC Results as a Tibble
  # -------------------------------#
  # The computed AUC values for the different models and prediction intervals are combined 
  # into a tibble for further analysis.
  tibble(
    auc_joint_36   = auc_joint_36$auc,
    auc_joint_60   = auc_joint_60$auc,
    auc_cox_36     = auc_cox_36$auc,
    auc_cox_60     = auc_cox_60$auc,
    auc_joint_1_36 = auc_joint_1_36$auc,
    auc_joint_1_60 = auc_joint_1_60$auc,
    auc_cox_1_36   = auc_cox_1_36$auc,
    auc_cox_1_60   = auc_cox_1_60$auc,
    auc_score_36 = auc_score$AUC$score$AUC[1],
    auc_score_60 = auc_score$AUC$score$AUC[2]
  )
}

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
    } else {  # Must be a Cox model (coxph)
      
      roc_result <- tvAUC_coxph(model, newdata = long_data, Tstart = Tstart, Dt = Dt, type = "roc")
      # Extract and round the AUC value from the ROC result.
      auc_value <- round(roc_result$auc, 3)
    }
    # Return a data frame with the ROC metrics for the current follow-up time.
    data.frame(
      FP       = roc_result$FP,  # False Positive Rate
      TP       = roc_result$TP,  # True Positive Rate
      FollowUp = Dt,             # Follow-up time interval
      AUC      = auc_value       # Area Under the Curve
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
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", 
                color = "black", linewidth = 1) +
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
  calibration_data <- list()  # Initialize list to store calibration data
  
  # Loop through each follow-up time (Dt)
  for (Dt in follow_up_times) {
    if (inherits(model, "jm")) {
      # Get calibration plot data without displaying the plot
      cal <- JMbayes2::calibration_plot(model, newdata = newdata, Tstart = Tstart, Dt = Dt, plot = FALSE)
      # Calculate calibration metrics (e.g., Integrated Calibration Index)
      ici <- calibration_metrics(model, newdata, Tstart = Tstart, Dt = Dt)
      # Store observed vs. predicted data along with the calibration metric
      cal_df <- tibble(
        observed = cal$observed,
        predicted = cal$predicted,
        follow_up = Dt, 
        ici = ici[1],
        model_label = paste("Follow-up:", Dt, "months\nICI:", round(ici, 2)))
      
      density_data <- tibble(
        density = density(cal$pi_u_t)$y / max(density(cal$pi_u_t)$y),
        preds = density(cal$pi_u_t)$x,
        ici = ici[1],
        model_label = paste("Follow-up:", Dt, "months\nICI:", round(ici, 2))
      )
    } else {  # Must be a Cox model (coxph)
      cox_df <- create_cox_df(newdata, Tstart)
      newd <- cox_df; newd$fup_time <- Dt; newd$mortality_status <- 1
      p <- 1-predict(model, type = "survival", newdata=newd)
      y <- with(cox_df, Surv(fup_time, mortality_status))
      cal <- pmcalibration(y = y, p = p, smooth = "rcs", nk=5, ci = "pw", time = Dt)
      # Store observed vs. predicted data along with the calibration metric
      cal_df <- tibble(
        observed = get_cc(cal)$p_c,
        predicted = get_cc(cal)$p,
        follow_up = Dt, 
        ici = cal$metrics[1],
        model_label = paste("Follow-up:", Dt, "months\nEavg:", round(ici, 2)))
      
      density_data <- tibble(
        density = density(p)$y / max(density(p)$y),
        preds = density(p)$x,
        ici = cal$metrics[1],
        model_label = paste("Follow-up:", Dt, "months\nEavg:", round(ici, 2)))
    }

    # Save both calibration and density data for this follow-up time
    calibration_data[[Dt]] <- list(cal_df = cal_df, density_data = density_data)
  }
  
  # Combine calibration data for all follow-up times into single data frames
  combined_cal_df <- bind_rows(lapply(calibration_data, `[[`, "cal_df"))
  combined_density_df <- bind_rows(lapply(calibration_data, `[[`, "density_data"))
  
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
    geom_line(data = combined_cal_df, aes(x = predicted, y = observed, color = model_label), linewidth = 1.2) +
    # Add a diagonal reference line representing perfect calibration
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", 
                color = "black", linewidth = 1) +
    # Overlay density curves (dashed lines) of predicted probabilities
    geom_line(data = combined_density_df, aes(x = preds, y = density, color = model_label),
              linetype = "dashed", linewidth = 1) +
    # Configure the primary y-axis and a secondary axis for density
    scale_y_continuous(
      name = "Observed Probability",
      sec.axis = sec_axis(~., name = "Density")
    ) +
    scale_color_brewer(palette = "Set1") +  # Use a colorblind-friendly palette
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
plot_dyn_pred <- function(joint_model, data, id, t0, fu_t = 60,
                          fill_CI_long = "#0000FF4D", fill_CI_event = "#FF00004D",
                          col_line_long = "#0000FF", col_line_event = "#FF0000",
                          col_points = "blue", cex_points = 2,
                          lwd_long = 1, lwd_event = 1, main = NULL,
                          xlim = NULL) {
  
  # -------------------------------------------------------------------------- #
  # Prepare the data for the specific subject and landmark time (t0)
  # -------------------------------------------------------------------------- #
  ND <- data %>%
    filter(id == !!id, time < t0) %>%  # Select records for the given subject before t0
    mutate(
      mortality_status = 0,   # Set mortality status to 0 for prediction purposes
      fup_time = t0           # Define the follow-up time as the landmark time t0
    )
  
  # -------------------------------------------------------------------------- #
  # Generate Longitudinal Predictions -----------------------------------------
  # -------------------------------------------------------------------------- #
  predLong <- predict(
    joint_model,
    newdata = ND,
    times = seq(t0, fu_t, length.out = 51),  # Predict from t0 up to 60 months
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
    times = seq(t0, fu_t, length.out = 51),  # Predict event probabilities over time
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
  y_axis_factor <- 800  # Factor to scale event probabilities for plotting
  
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
  ggplot() +
    # Plot longitudinal predictions with a confidence interval ribbon
    geom_ribbon(data = long_data, aes(x = times_long, ymin = low_long, ymax = upp_long),
                fill = fill_CI_long, alpha = 0.3) +
    geom_line(data = long_data, aes(x = times_long, y = preds_long),
              color = col_line_long, linewidth = lwd_long) +
    # Plot the observed longitudinal measurements
    geom_point(data = observed_data, aes(x = x, y = y),
               color = col_points, size = cex_points) +
    # Draw a vertical dashed line at the last observed time point
    geom_vline(xintercept = last_times + 0.01, linetype = "dashed", color = "gray") +
    # Plot event predictions with a confidence interval ribbon
    geom_ribbon(data = event_data, aes(x = times_event, ymin = low_event, ymax = upp_event),
                fill = fill_CI_event, alpha = 0.3) +
    geom_line(data = event_data, aes(x = times_event, y = preds_event),
              color = col_line_event, linewidth = lwd_event) +
    # Configure the y-axis to include a secondary axis for event probability
    scale_y_continuous(
      sec.axis = sec_axis(~ . / y_axis_factor, name = "Event Probability",
                          breaks = seq(0, 1, by = 0.2))
    ) +
    # Configure x-axis breaks
    scale_x_continuous(breaks = seq(0, max(event_data$times_event), by = 6)) +
    coord_cartesian(ylim = y_lim) +
    # Add labels and title to the plot
    labs(
      title = paste("Dynamic Predictions Plot", main),
      x = "Time (months)",
      y = "6 minute walk distance (meters)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.position = "bottom",
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      plot.title = element_text(size = 16, hjust = 0.5)
    )
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
    arrange(id, time) %>%                 # Sort by patient ID and time
    group_by(id) %>%
    filter(row_number() == 1) %>%         # Keep only the first visit (baseline)
    ungroup()
  
  boot_joint <- boot_data %>%
    filter(time < 12) %>%  # Select records for the given subject before t0
    mutate(
      mortality_status = 0,   # Set mortality status to 0 for prediction purposes
      fup_time = 12           # Define the follow-up time as the landmark time t0
    )
  
  # Compute AUC measures using tvAUC.
  auc_36   <- tvAUC(object = jointFit_train, newdata = boot_data, Tstart = 0, Dt = 36)
  auc_60   <- tvAUC(object = jointFit_train, newdata = boot_data, Tstart = 0, Dt = 60)
  auc_1_36 <- tvAUC(object = jointFit_train, newdata = boot_data, Tstart = 12, Dt = 24)
  auc_1_60 <- tvAUC(object = jointFit_train, newdata = boot_data, Tstart = 12, Dt = 48)
  
  cox_preds <- survfit(CoxFit_six_train, newdata = boot_cox)
  joint_preds <- predict(jointFit_train, newdata = boot_joint, 
                         process = "event",
                         times = 36,
                         return_newdata = TRUE)
  
  cox_pred <- cox_preds$cumhaz[cox_preds$time == 36]
  joint_pred <- joint_preds[joint_preds$time == 36, "pred_CIF"]
  
  timeROC_cox <- timeROC::timeROC(T=boot_cox$fup_time,
                                  delta=boot_cox$mortality_status,
                                  marker=cox_pred,
                                  cause=1,weighting="cox",
                                  times= 36)
  
  timeROC_joint <- timeROC::timeROC(T=boot_cox$fup_time,
                                    delta=boot_cox$mortality_status,
                                    marker=joint_pred,
                                    cause=1,weighting="cox",
                                    times= 36)
  
  survivalROC_cox <- survivalROC::survivalROC(Stime=boot_cox$fup_time,
                                              status=boot_cox$mortality_status,      
                                              marker = cox_pred,     
                                              predict.time = 36, method="KM")
  
  survivalROC_joint <- survivalROC::survivalROC(Stime=boot_cox$fup_time,
                                                status=boot_cox$mortality_status,      
                                                marker = joint_pred,     
                                                predict.time = 36, method="KM")
  
  risksetROC_cox <- risksetROC::risksetAUC(Stime=boot_cox$fup_time,
                                           status=boot_cox$mortality_status,      
                                           marker = cox_pred,     
                                           tmax =  36)
  
  risksetROC_joint <- risksetROC::risksetAUC(Stime=boot_cox$fup_time,
                                             status=boot_cox$mortality_status,      
                                             marker = joint_pred,     
                                             tmax =  36)
  
  # Return the AUC results as a tibble.
  tibble(
    auc_36   = auc_36$auc,
    auc_60   = auc_60$auc,
    auc_1_36 = auc_1_36$auc,
    auc_1_60 = auc_1_60$auc,
    timeROC_cox   = timeROC_cox$AUC[2],
    timeROC_joint   = timeROC_joint$AUC[2],
    survivalROC_cox   = survivalROC_cox$AUC,    
    survivalROC_joint   = survivalROC_joint$AUC,
    risksetROC_cox = risksetROC_cox$Cindex,
    risksetROC_joint = risksetROC_joint$Cindex
  )
}
