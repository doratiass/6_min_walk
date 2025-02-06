source("scripts/00 - funcs.R")
load("data/CHF_data.RData")
library(JMbayes2)
library(lmtest)
library(MuMIn)
library(rsample)

# Data Pre-processing ----------------------------------------------------------
follow_up_time <- 12

## Identify IDs with missing variables ----------------------------------------
total_ids <- length(unique(clean_df$id))

clean_df %>%
  # filter(if_any(c(all_of(seattle_vars), x6mw_dist_meter), ~ is.na(.))) %>%
  filter(if_any(nyha, ~ is.na(.))) %>%
  distinct(id) %>%
  pull(id) -> miss_id_vars

print(paste0("Number of patients with missing variables: ", length(miss_id_vars), 
             " (", round(length(miss_id_vars)/total_ids*100, 2), " %) out of ",
             total_ids, " total patients"))

## Create models df ------------------------------------------------------------
# Create cox df
cox_df <- clean_df %>%
  filter(!(id %in% miss_id_vars)) %>%
  arrange(id, time) %>%
  group_by(id) %>%
  mutate(n_visits = n()) %>%
  filter(row_number() == 1) %>%
  ungroup()

# Create long df
long_df <- clean_df %>%
  # filter(time <= follow_up_time) %>%
  filter(id %in% cox_df$id) %>%
  arrange(id, time) %>% 
  group_by(id) %>%
  mutate(age = first(age),
         nyha = first(nyha),
         ischemic_etiology = any(ischemic_etiology, na.rm = TRUE)) %>%
  group_by(id, fup_time, mortality_status, time, age, gender, nyha, ischemic_etiology) %>%
  summarise(x6mw_dist_meter = mean(x6mw_dist_meter, na.rm = TRUE),
            .groups = "drop")

### Create long df with only patients with more than 1 visit -------------------
long_df_2 <- long_df %>%
  group_by(id) %>%
  filter(n() >= 2) %>%
  ungroup()

cox_df_2 <- cox_df %>%
  filter(id %in% long_df_2$id)

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

# measure the number of visits per patient
long_df %>%
  # filter(time <= follow_up_time) %>%
  group_by(id) %>%
  summarise(n = n()) %>%
  group_by(n) %>%
  summarise(t = n())

long_df %>%
  group_by(id) %>%
  mutate(six_sd = sd(x6mw_dist_meter, na.rm = TRUE)) %>%
  filter(six_sd > 25) %>%
  pull(id) -> high_sd_ids

# Model Development -----------------------------------------------------------
## Linear mixed models --------------------------------------------------------
###  Finding the best fitting model -------------------------------------------
lme_6min_0 <- lme(x6mw_dist_meter ~ ns(time, df = 3) + age + gender, 
                  data = long_df, random = ~ ns(time, df = 3) | id, 
                  control = lmeControl(opt = 'optim'))
r.squaredGLMM(lme_6min_0)

lme_6min_n <- lme(x6mw_dist_meter ~ ns(time, df = 3) + age + gender + nyha, 
                  data = long_df, random = ~ ns(time, df = 3) | id, 
                  control = lmeControl(opt = 'optim'))
r.squaredGLMM(lme_6min_0)

anova(lme_6min_0, lme_6min_n)

lme_6min_int_sparse <- lme(x6mw_dist_meter ~ ns(time, df = 3)*(age) + gender + nyha, 
                           data = long_df, random = ~ ns(time, df = 3) | id, 
                           control = lmeControl(opt = 'optim'))
r.squaredGLMM(lme_6min_int_sparse)

lme_6min_int <- lme(x6mw_dist_meter ~ ns(time, df = 3)*(gender + age) + nyha, 
                    data = long_df, random = ~ ns(time, df = 3) | id, 
                    control = lmeControl(opt = 'optim'))
r.squaredGLMM(lme_6min_int)

anova(lme_6min_n, lme_6min_int_sparse, lme_6min_int)

lme_6min <- lme(x6mw_dist_meter ~ ns(time, df = 3) * (age + gender) + nyha, 
                data = long_df, random = ~ ns(time, df = 3) | id, 
                control = lmeControl(opt = 'optim'))

summary(lme_6min)
r.squaredGLMM(lme_6min)

lme_6min_2 <- lme(x6mw_dist_meter ~ ns(time, df = 3) * (age + gender) + nyha, 
                data = long_df_2, random = ~ ns(time, df = 3) | id, 
                control = lmeControl(opt = 'optim'))

summary(lme_6min_2)
r.squaredGLMM(lme_6min_2)

# lme_6min_last <- lme(x6mw_dist_meter ~ time * (age + gender) + nyha, 
#                   data = long_df_last, random = ~ time | id, 
#                   control = lmeControl(opt = 'optim'))
# 
# summary(lme_6min_last)
# r.squaredGLMM(lme_6min_last)
# 
# lme_6min_2_last <- lme(x6mw_dist_meter ~ ns(time, df = 3) * (age + gender) + nyha, 
#                      data = long_df_2_last, random = ~ ns(time, df = 3) | id, 
#                      control = lmeControl(opt = 'optim'))
# 
# summary(lme_6min_2_last)
# r.squaredGLMM(lme_6min_2_last)

### plot time vs 6min walk distance by age -------------------------------------
pred_6mwt_data <- data.frame(
  time = rep(round(seq(0, 72, length.out = 50)),4),
  age = rep(c(rep(40,50),rep(50,50),rep(60,50), rep(70,50)),2),
  gender = c(rep("Female",200),rep("Male",200)),
  nyha = "II"
)

bind_rows(
  pred_6mwt_data %>% mutate(pred = predict(lme_6min, newdata = pred_6mwt_data, level = 0),
                       model = "model_1"),
  pred_6mwt_data %>% mutate(pred = predict(lme_6min_2, newdata = pred_6mwt_data, level = 0),
                       model = "model_2")) %>%
  ggplot(aes(x = time, y = pred, colour = age)) +
  geom_point() +
  facet_grid(gender~model)

## Cox proportional hazards models --------------------------------------------
CoxFit <- coxph(Surv(fup_time, mortality_status) ~ gender + age + nyha, 
                data = cox_df, model = TRUE,
                ,x=TRUE,y=TRUE)

CoxFit_six <- coxph(Surv(fup_time, mortality_status) ~ gender + age + nyha + x6mw_dist_meter, 
                    data = cox_df, model = TRUE,
                    ,x=TRUE,y=TRUE)

CoxFit_2 <- coxph(Surv(fup_time, mortality_status) ~ gender + age + nyha, 
                  data = cox_df_2, model = TRUE,
                  x=TRUE,y=TRUE)

## Joint models ---------------------------------------------------------------
jointFit_0 <- jm(CoxFit, lme_6min, time_var = "time")
summary(jointFit_0)

jointFit_1 <- jm(CoxFit, lme_6min, time_var = "time",
                 functional_forms = list("x6mw_dist_meter" = 
                                           ~ value(x6mw_dist_meter) + slope(x6mw_dist_meter)))

summary(jointFit_1)
exp(coef(jointFit_1)$association)


jointFit_2 <- jm(CoxFit, lme_6min, time_var = "time",
                 functional_forms = list("x6mw_dist_meter" = 
                                           ~ value(x6mw_dist_meter)*slope(x6mw_dist_meter)))
summary(jointFit_2)
exp(coef(jointFit_2)$association)

jointFit_3 <- jm(CoxFit_six, lme_6min, time_var = "time",
                 functional_forms = list("x6mw_dist_meter" = 
                                           ~ value(x6mw_dist_meter) + slope(x6mw_dist_meter)))
summary(jointFit_3)
exp(coef(jointFit_3)$association)


form_splines <- ~ value(x6mw_dist_meter) * ns(time, k = c(10, 20, 40), B = c(0, 72)) + slope(x6mw_dist_meter)
jointFit_4 <- jm(CoxFit_six, lme_6min, time_var = "time",
                 functional_forms = form_splines)
summary(jointFit_4)
exp(coef(jointFit_4)$association)

compare_jm(jointFit_0, jointFit_1, jointFit_2, jointFit_3,jointFit_4)

jointFit <- jm(CoxFit, lme_6min, time_var = "time",
               functional_forms = list("x6mw_dist_meter" = 
                                         ~ value(x6mw_dist_meter) + slope(x6mw_dist_meter)))
summary(jointFit)
exp(coef(jointFit)$association)

# jointFit_last <- jm(CoxFit, lme_6min_last, time_var = "time",
#                  functional_forms = list("x6mw_dist_meter" = 
#                                            ~ value(x6mw_dist_meter)+slope(x6mw_dist_meter)))
# summary(jointFit_last)
# exp(coef(jointFit_last)$association)

jointFit_2 <- jm(CoxFit_2, lme_6min_2, time_var = "time",
                 functional_forms = list("x6mw_dist_meter" = 
                                           ~ value(x6mw_dist_meter)+slope(x6mw_dist_meter)))
summary(jointFit_2)
exp(coef(jointFit_2)$association)

# jointFit_2_last <- jm(CoxFit_2, lme_6min_2_last, time_var = "time",
#                  functional_forms = list("x6mw_dist_meter" = 
#                                            ~ value(x6mw_dist_meter)+slope(x6mw_dist_meter)))
# summary(jointFit_2_last)
# exp(coef(jointFit_2_last)$association)

#save(jointFit, file = "data/jointFit.RData")

# Model Evaluation ------------------------------------------------------------
## Choose initial fu time for the model evaluation ----------------------------
# TODO - what is the best time to start the model evaluation
long_df %>%
  ggplot(aes(x = time)) +
  geom_histogram(binwidth = 1) +
  geom_vline(xintercept = 12, color = "red") +
  geom_vline(xintercept = 24, color = "blue")

## Discrimination -------------------------------------------------------------
roc_data <- create_roc_data(
  joint_model = jointFit,
  long_data = long_df,
  Tstart = follow_up_time,
  follow_up_times = c(12, 36, 60)
)

roc_data_last <- create_roc_data(
  joint_model = jointFit,
  long_data = long_df_last,
  Tstart = follow_up_time,
  follow_up_times = c(12, 36, 60)
)

roc_data_2 <- create_roc_data(
  joint_model = jointFit_2,
  long_data = long_df_2,
  Tstart = follow_up_time,
  follow_up_times = c(12, 36, 60)
)

roc_data_2_last <- create_roc_data(
  joint_model = jointFit_2,
  long_data = long_df_2_last,
  Tstart = follow_up_time,
  follow_up_times = c(12, 36, 60)
)


## Calibration ---------------------------------------------------------------
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
## Plots ---------------------------------------------------------------------
# TODO compare with base cox models
ggarrange(
  plot_tvROC(roc_data),
  plot_tvROC(roc_data_last),
  plot_cal(cal_data),
  plot_cal(cal_data_last))

ggarrange(
  plot_tvROC(roc_data_2),
  plot_tvROC(roc_data_2_last),
  plot_cal(cal_data_2),
  plot_cal(cal_data_2_last))

# Dynamic Predictions --------------------------------------------------------
ggarrange(
  plot_dyn_pred(jointFit, long_df, id = 152227420, t0 = 2),
  plot_dyn_pred(jointFit, long_df, id = 152227420, t0 = 12),
  plot_dyn_pred(jointFit, long_df, id = 152227420, t0 = 24),
  plot_dyn_pred(jointFit, long_df, id = 152227420, t0 = 36),
  nrow = 4)

ggarrange(
  plot_dyn_pred(jointFit, long_df, id = 274751, t0 = 12),
  plot_dyn_pred(jointFit, long_df, id = 274751, t0 = 24),
  plot_dyn_pred(jointFit, long_df, id = 274751, t0 = 36),
  nrow = 3)