# Load necessary packages
library(future)
library(furrr)
library(tidyverse)   # For data manipulation
library(JMbayes2)    # For tvAUC

# Set seed for reproducibility
set.seed(123)

# Your data and unique subject IDs
data <- CVdats$testing[[1]]
unique_ids <- unique(data$id)
B <- 100  # Number of bootstrap iterations (increase as needed)

# Define a function for a single bootstrap iteration
bootstrap_iteration <- function(iteration) {
  # Resample subject IDs with replacement
  boot_ids <- sample(unique_ids, size = length(unique_ids), replace = TRUE)
  
  # For each sampled subject, extract rows and assign a new unique sequential ID
  boot_data <- bind_rows(lapply(seq_along(boot_ids), function(j) {
    data %>%
      filter(id == boot_ids[j]) %>%
      mutate(id = j)
  }))
  
  # Compute the AUC measures using tvAUC (modify parameters as needed)
  auc_36  <- tvAUC(object = jointFit_train, newdata = boot_data, Tstart = 0, Dt = 36)
  auc_60  <- tvAUC(object = jointFit_train, newdata = boot_data, Tstart = 0, Dt = 60)
  # auc_1_36 <- tvAUC(object = jointFit_train, newdata = boot_data, Tstart = 12, Dt = 24)
  # auc_1_60 <- tvAUC(object = jointFit_train, newdata = boot_data, Tstart = 12, Dt = 48)
  auc_1_36 <- tvAUC(object = jointFit_train, newdata = boot_data, Tstart = 12, Dt = 24)
  auc_1_60 <- tvAUC(object = jointFit_train, newdata = boot_data, Tstart = 12, Dt = 48)
  
  # Return results as a tibble
  tibble(
    auc_36  = auc_36$auc,
    auc_60  = auc_60$auc,
    auc_1_36 = auc_1_36$auc,
    auc_1_60 = auc_1_60$auc
  )
}

# Set up the future plan:
# Use 'multisession' for Windows or cross-platform; use 'multicore' on Unix-like systems if preferred.
plan(multisession, workers = parallel::detectCores() - 1)

# Run bootstrap iterations in parallel using future_map (from furrr)
bootROC_results <- future_map(1:B, ~ bootstrap_iteration(.x), seed = TRUE)

# Combine the results into a single data frame
bootROC_results_df <- bind_rows(bootROC_results)
print(bootROC_results_df)

# Optionally, switch back to sequential processing:
plan(sequential)

bootROC_results_df %>%
  pivot_longer(cols = everything(), names_to = "AUC", values_to = "value") %>%
  group_by(AUC) %>%
  summarise(
    mean = mean(value),
    sd = sd(value),
    median = median(value),
    lower_ci = quantile(value, 0.025),
    upper_ci = quantile(value, 0.975)
  )
