# ============================================================================ #
# Prediction of 6MWT Change and Survival Analysis
#
# Overview:
#   This script predicts the 6-minute walk test (6MWT) distance over time 
#   using a linear mixed-effects model and evaluates the prognostic value of 
#   the predicted change on mortality using Kaplan-Meier survival analysis.
#
# Key Components:
#   1. Load required libraries and data.
#   2. Predict 6MWT trajectories using mixed-effects models.
#   3. Group patients into quartiles based on predicted change.
#   4. Plot predicted 6MWT values over time by quartile.
#   5. Evaluate the association between predicted change and survival.
#
# Required Packages:
#   - survival, survminer: For survival analysis.
#   - dplyr, tidyr, ggplot2: For data manipulation and visualization.
# ============================================================================ #

# ============================================================================ #
# Load Libraries and Data ---------------------------------------------------- #
# ============================================================================ #
source("scripts/00 - funcs.R")            # Load custom helper functions
load("data/CHF_models.RData")             # Load saved survival and joint models

library(survival)                         # Cox models and survival analysis
library(survminer)                        # Kaplan-Meier plotting
library(dplyr)                            # Data manipulation
library(tidyr)                            # Data tidying
library(ggplot2)                          # Data visualization
library(cowplot)
# ============================================================================ #
# Create quartile groups ------------------------------------------------------ 
# ============================================================================ #
# -------------------------------------------------------------------------- #
## Create Prediction Grid for Each Patient ------------------------------------
# -------------------------------------------------------------------------- #
# Create a grid for predicted 6MWT over 12 months for each patient
time_id <- long_df %>%
  filter(time <= tests_time,
         fup_time >= tests_time) %>%  # Select records for the given subject before t0
  mutate(
    mortality_status = 0,   # Set mortality status to 0 for prediction purposes
    fup_time = 12           # Define the follow-up time as the landmark time t0
  )

# -------------------------------------------------------------------------- #
## Predict 6MWT Trajectories -------------------------------------------------- 
# -------------------------------------------------------------------------- #
# Predict 6MWT values using the mixed-effects model
predLong <- predict(
  jointFit,
  newdata = time_id,
  times = c(0:12),  # Predict from t0 up to fu_t
  type = "subject_specific",
  return_newdata = TRUE,
  control = list(all_times = TRUE)
)

new_df <- predLong[[2L]] %>%
  arrange(id, time)
# -------------------------------------------------------------------------- #
## Calculate Predicted Change and Standardize --------------------------------- 
# -------------------------------------------------------------------------- #
# Compute change in predicted 6MWT values over time
quan <- new_df %>%
  group_by(id,gender) %>%
  summarise(change = last(pred_x6mw_dist_meter) - first(pred_x6mw_dist_meter)) %>%
  ungroup() %>%
  group_by(gender) %>%
  mutate(
    sd = sd(change),                           # Compute standard deviation
    sds = change / sd,                         # Standardize the change
    quartile = factor(case_when(
      sds < 0 ~ "Low",                      # Lower quartile threshold
      sds > 1 ~ "High",                      # Upper quartile threshold
      TRUE ~ "Medium"                          # Middle quartile
    ), levels = c("Low", "Medium", "High"))
  ) %>%
  ungroup()

quan %>%
  group_by(quartile) %>%
  summarise(
    mean_change = mean(change),
    sd_change = sd(change),
    median_change = median(change),
    lower_ci = quantile(change, 0.025),
    upper_ci = quantile(change, 0.975)
  )
# -------------------------------------------------------------------------- #
## Plot Predicted 6MWT Values by Quartile ------------------------------------- 
# -------------------------------------------------------------------------- #
new_df %>%
  left_join(quan %>% select(id, quartile), by = "id") %>%
  ggplot(aes(x = time, y = pred_x6mw_dist_meter, color = quartile)) +
  geom_line(aes(group = id), alpha = 0.05) +              # Individual trajectories
  geom_smooth(method = "gam", se = FALSE) +             # Smoothed trend line
  scale_x_continuous(breaks = 0:12, expand = c(0,0)) +                    # X-axis breaks
  scale_y_continuous(breaks = seq(0, 800, 100), limits = c(0, 800), expand = c(0,0)) +               # Y-axis limits
  theme_minimal() +
  labs(
    title = "", #Predicted 6-Minute Walk Distance Over Time
    x = "Time (months)",
    y = "6-Minute Walk Distance (meters)",
    color = "Group"
  ) +
  theme(legend.position = "bottom") +
 # facet_wrap(~gender) +
  scale_color_brewer(palette="Set1") -> pred_6mwt_plot

# ============================================================================ #
# Survival Analysis -----------------------------------------------------------
# ============================================================================ #
# -------------------------------------------------------------------------- #
## Prepare Data for Survival Analysis ----------------------------------------- 
# -------------------------------------------------------------------------- #
# Create survival data frame with updated follow-up time
surv_df <- long_df %>%
  filter(time <= tests_time, fup_time >= tests_time) %>%
  group_by(id) %>%
  arrange(id, time) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(new_fup = fup_time - tests_time) %>%
  left_join(quan %>% select(id, quartile), by = "id")

# Create survival object using follow-up time and mortality status
surv_object <- Surv(time = surv_df$new_fup, event = surv_df$mortality_status)

# -------------------------------------------------------------------------- #
## Kaplan-Meier Curves -------------------------------------------------------- 
# -------------------------------------------------------------------------- #
# Fit Kaplan-Meier model by quartile group
fit <- survfit(surv_object ~ quartile, data = surv_df)

# Plot Kaplan-Meier Curves
ggsurvplot(
  fit,
  data = surv_df,
  xlim = c(0, 60),
  conf.int = TRUE,               # Show confidence intervals
 # pval = TRUE,                   # Display p-value from log-rank test
  risk.table = TRUE,             # Add risk table below the plot
  risk.table.col = "strata",     # Color the risk table by group
  legend.title = "Group",
  legend.labs = c("Low", "Medium", "High"),
  xlab = "Time (months)",
  ylab = "Survival Probability",
  palette = "Set1",             # Use a clean color palette
  ggtheme = theme_minimal()
) -> km_plot 

# -------------------------------------------------------------------------- #
## Figure 2 ------------------------------------------------------------------
# -------------------------------------------------------------------------- #
ggarrange(pred_6mwt_plot,
          km_plot$plot + theme(legend.position = "none"),
          ncol = 2,
          labels = c("A", "B"),
          common.legend = TRUE,
          legend = "bottom")

ggsave(file.path("export", "fig_2.jpeg"), last_plot(), 
       width = 30, height = 15, dpi = 300, background = "white", units = "cm")

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
  plot_dyn_pred(jointFit, long_df, id = 756, t0 = 0.1), 
  plot_dyn_pred(jointFit, long_df, id = 756, t0 = 3.1),
  plot_dyn_pred(jointFit, long_df, id = 756, t0 = 6.1),
  plot_dyn_pred(jointFit, long_df, id = 756, t0 = 12.1),
  nrow = 4
)

# HIGH - 1003, 444, 172, 657, 654
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
