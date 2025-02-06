library(tidyverse)
library(gtsummary)
library(ggpubr)
library(gridExtra)
library(flextable)
library(survminer)

# Data cleaning functions -----------------------------------------------------
merge_rows <- function(x) {
  if(is.numeric(x)) {
    return(mean(x, na.rm = TRUE))  
  } else {
    return(na.omit(x)[1]) 
  }
}

# Labeling variables ----------------------------------------------------------
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

label_get <- function(x) {ifelse(x %in% vars_dict$var,
                                 vars_dict[vars_dict$var == x,"name",drop = TRUE],
                                 x)}

vars_label <- function(x) {sapply(x, label_get, USE.NAMES = FALSE)}

var_get <- function(x) {ifelse(x %in% vars_dict$name,
                               vars_dict[vars_dict$name == x,"var",drop = TRUE],
                               x)}

# Model Evaluation ------------------------------------------------------------
create_roc_data <- function(joint_model, long_data, Tstart, follow_up_times) {
  roc_list <- list()
  
  # Loop through follow-up times and calculate the ROC
  for (Dt in follow_up_times) {
    roc_result <- tvROC(joint_model, newdata = long_data, Tstart = Tstart, Dt = Dt)
    auc <- round(tvAUC(roc_result)$auc, 3)
    roc_list[[as.character(Dt)]] <- list(roc = roc_result, auc = auc)
  }
  
  # Create a tibble with the required data for plotting
  roc_data <- do.call(bind_rows, lapply(names(roc_list), function(Dt) {
    roc <- roc_list[[Dt]]$roc
    auc <- roc_list[[Dt]]$auc
    data.frame(
      FP = roc$FP,
      TP = roc$TP,
      FollowUp = as.numeric(Dt),
      AUC = auc
    )
  }))
  
  return(roc_data)
}

plot_tvROC <- function(roc_data, model_name = NULL) {
  if (!is.null(model_name)) {
    model_name <- paste("for", model_name)
  }
  roc_data %>%
    mutate(model_name =  paste("Follow-up:", FollowUp, "months\nAUC:", AUC)) %>%
    ggplot(aes(x = FP, y = TP, 
               color = model_name)) +
    geom_line(linewidth = 1) +
    scale_color_brewer(palette = "Set1") +  # Use a colorblind-friendly palette
    labs(
      title = paste("Time-Dependent ROC Curves", model_name),
      x = "1 - Specificity",
      y = "Sensitivity",
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
}

calc_cal_metrics <- function(joint_model, newdata, Tstart, follow_up_times) {
  # Initialize list for storing data
  calibration_data <- list()
  
  # Loop through follow-up times and compute calibration data
  for (Dt in follow_up_times) {
    cal <- JMbayes2::calibration_plot(joint_model, newdata = newdata, Tstart = Tstart, Dt = Dt, plot = FALSE)
    ici <- calibration_metrics(joint_model, newdata, Tstart = Tstart, Dt = Dt)
    
    # Store observed vs predicted and density data
    cal_df <- tibble(
      observed = cal$observed,
      predicted = cal$predicted,
      follow_up = Dt, 
      ici = ici[1],
      model_label = paste("Follow-up:", Dt, "months\nICI:", round(ici, 2))
    )
    density_data <- tibble(
      density = density(cal$pi_u_t)$y / max(density(cal$pi_u_t)$y),  # Normalize density
      preds = density(cal$pi_u_t)$x,
      ici = ici[1],
      model_label = paste("Follow-up:", Dt, "months\nICI:", round(ici, 2))
    )
    calibration_data[[Dt]] <- list(cal_df = cal_df, density_data = density_data)
  }
  
  # Combine data for all follow-up times
  combined_cal_df <- bind_rows(lapply(calibration_data, `[[`, "cal_df"))
  combined_density_df <- bind_rows(lapply(calibration_data, `[[`, "density_data"))
  
  # Return combined data
  list(cal_df = combined_cal_df, density_df = combined_density_df)
}

# Function to plot calibration
plot_cal <- function(calibration_data, model_name = NULL) {
  if (!is.null(model_name)) {
    model_name <- paste("for", model_name)
  }
  # Extract data
  combined_cal_df <- calibration_data$cal_df
  combined_density_df <- calibration_data$density_df
  
  # Generate the calibration plot
  plot <- ggplot() +
    # Calibration curves
    geom_line(data = combined_cal_df, aes(x = predicted, y = observed, color = model_label), linewidth = 1.2) +
    # Diagonal reference line
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", 
                color = "black", linewidth = 1) +
    # Density curves
    geom_line(data = combined_density_df, aes(x = preds, y = density, color = model_label),
              linetype = "dashed", linewidth = 1) +
    # Improved scales and labels
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
    # Enhanced theme for better aesthetics
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

plot_dyn_pred <- function(joint_model, data, id, t0, fu_t = 60,
                          fill_CI_long = "#0000FF4D", fill_CI_event = "#FF00004D",
                          col_line_long = "#0000FF", col_line_event = "#FF0000",
                          col_points = "blue", cex_points = 2,
                          lwd_long = 1, lwd_event = 1, main = NULL,
                          xlim = NULL) {
  
  # Prepare longitudinal data
  ND <- data %>%
    filter(id == !!id, time < t0) %>%
    mutate(
      mortality_status = 0,   # Add mortality status column
      fup_time = t0          # Add follow-up time column
    )
  
  # Generate longitudinal predictions
  predLong <- predict(
    joint_model,
    newdata = ND,
    times = seq(t0, 60, length.out = 51),
   # times = seq(t0, fu_t+t0, length.out = 51),
    type = "subject_specific",
    return_newdata = TRUE
  )
  pred_Long <- bind_rows(predLong[[1L]], predLong[[2L]]) 
  
  # Generate survival predictions
  pred_Event <- predict(
    joint_model,
    newdata = ND,
    process = "event",
    times = seq(t0, 60, length.out = 51),
    #times = seq(t0, fu_t+t0, length.out = 51),
    return_newdata = TRUE
  )
  
  # Extract attributes
  id_var <- attr(predLong, "id_var")
  time_var <- attr(predLong, "time_var")
  Time_var <- attr(predLong, "Time_var")
  resp_vars <- attr(predLong, "resp_vars")
  ranges <- attr(predLong, "ranges")
  last_times <- attr(predLong, "last_times")
  y <- attr(predLong, "y")
  times_y <- attr(predLong, "times_y")
  
  # Process longitudinal predictions
  pred_Long <- pred_Long %>%
    filter(!!sym(time_var) <= as.numeric(last_times)) %>%
    rename(times_long = !!sym(time_var))
  
  outcome_inx <- grep("pred_", names(pred_Long), fixed = TRUE)
  
  long_data <- pred_Long %>%
    transmute(times_long, 
              preds_long = pred_Long[[outcome_inx]],
              low_long = pred_Long[[outcome_inx+1]],
              upp_long = pred_Long[[outcome_inx+2]]) %>%
    filter(!is.na(preds_long))
  
  observed_data <- data.frame(
    x = unlist(times_y),
    y = unlist(y)
  ) %>%
    filter(!is.na(y))
  
  # Process event predictions
  #y_axis_factor <- max(long_data$upp_long) / max(pred_Event[[grep("pred_", names(pred_Event)) + 2]])
  y_axis_factor <-800 
  
  event_data <- data.frame(
    times_event = pred_Event[[Time_var]],
    preds_event = pred_Event[[grep("pred_", names(pred_Event))]],
    low_event = pred_Event[[grep("pred_", names(pred_Event)) + 1]],
    upp_event = pred_Event[[grep("pred_", names(pred_Event)) + 2]]
  ) %>%
    filter(!is.na(preds_event)) %>%
    mutate_at(vars(preds_event, low_event, upp_event), ~ . * y_axis_factor)
  
  # Set limits
  y_lim <- c(0, 1.01) * y_axis_factor
  
  if (is.null(main)) {
    main <- paste("for", round(t0), "months follow-up")
  }
  
  # Plot
  ggplot() +
    # Longitudinal predictions
    geom_ribbon(data = long_data, aes(x = times_long, ymin = low_long, ymax = upp_long),
                fill = fill_CI_long, alpha = 0.3) +
    geom_line(data = long_data, aes(x = times_long, y = preds_long),
              color = col_line_long, linewidth = lwd_long) +
    geom_point(data = observed_data, aes(x = x, y = y),
               color = col_points, size = cex_points) +
    geom_vline(xintercept = last_times + 0.01, linetype = "dashed", color = "gray") +
    # Event predictions
    geom_ribbon(data = event_data, aes(x = times_event, ymin = low_event, ymax = upp_event),
                fill = fill_CI_event, alpha = 0.3) +
    geom_line(data = event_data, aes(x = times_event, y = preds_event),
              color = col_line_event, linewidth = lwd_event) +
    # Scales
    scale_y_continuous(
      sec.axis = sec_axis(~ . / y_axis_factor, name = "Event Probability",
                          breaks = seq(0, 1, by = 0.2))
    ) +
    scale_x_continuous(breaks = seq(0, max(event_data$times_event), by = 6)) +
    coord_cartesian(ylim = y_lim) +
    # Aesthetics
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
