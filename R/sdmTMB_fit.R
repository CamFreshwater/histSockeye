# Fit sdmTMB model to century of Nass data 
# Faster alternative to brms_fit while still accounting for variable dispersion
# 1) Import Skip's first century of data and clean as necessary
# 2) Estimate temporal trend using GAM and brms to account for changes in
# variability with sampling regime
# Created by C Freshwater Dec 23, 2021
# Update w/ edited data
# -------------------------------------------------

library(tidyverse)
library(dplyr)
library(forcats)
library(sdmTMB)


# import nass data
dat_in <- read.table(here::here("data", "nasscenturyv3.txt"))
colnames(dat_in) <- c("year", "month", "day", "hart_species", "sex", "age", 
                      "fl")

set.seed(123)
dat <- dat_in %>%
  mutate(
    sex = ifelse(sex == 1, "male", "female") %>% as.factor(),
    age = factor(as.character(age)),
    year_f = as.factor(year),
    age_f = age,
    date = as.POSIXct(paste(day, month, year, sep = "-"),
                      format = "%d-%m-%Y"),
    yday = lubridate::yday(date),
    yday_c = yday - mean(yday, na.rm = TRUE)
  ) %>% 
  #drop secondary age classes and years missing capture data
  filter(
    age %in% c("42", "52", "53", "63"),
    !is.na(yday),
    !is.na(fl)
  )  %>%
  # Add a factor representing sample collection to visualize effects 
  mutate(
    period = case_when(
      year < 1957 ~ "Gilbert-Clemens",
      year > 1972 & year < 1994 ~ "Monkley Dump",
      year > 1993 ~ "Nisga'a",
      TRUE ~ "Bilton"
    ),
    period = fct_relevel(as.factor(period), 
                         "Gilbert-Clemens",
                         "Bilton",
                         "Monkley Dump",
                         "Nisga'a"),
    fl_cm = fl / 10,
    # add random identifier to split datasets for subsetting,
    data_group = sample.int(4, nrow(.), replace = T)
  ) %>% 
  droplevels()
levels(dat$age_f) <- c("4[2]", "5[2]", "5[3]", "6[3]")


# palettes for raw data and model predictions
fill_pal <- c("black", "white")
names(fill_pal) <- c("annual average", "individual measurements")
shape_pal <- c(21, 22, 23, 24)
names(shape_pal) <- unique(dat$period)


# RAW DATA PLOTS ---------------------------------------------------------------

# add dummy year levels for missing data
dat2 <- expand.grid(
  year_f = seq(min(dat$year), max(dat$year), by = 1) %>% as.factor(),
  age = unique(dat$age),
  sex = unique(dat$sex)
) %>% 
  left_join(., dat, by = c("year_f", "age", "sex")) 

mean_dat <- dat2 %>% 
  group_by(year_f, age, age_f, period, sex) %>% 
  summarize(
    n = n(),
    mean_yday = mean(yday, na.rm = T),
    sd_yday = sd(yday_c, na.rm = T),
    mean_fl = mean(fl, na.rm = T),
    sd_fl = sd(fl, na.rm = T),
    se_fl = sd_fl / sqrt(n),
    up = mean_fl + (1.96 * se_fl),
    lo = mean_fl - (1.96 * se_fl),
    .groups = "drop"
  ) %>% 
  distinct() 
empty_yrs <- mean_dat %>% filter(is.na(mean_fl)) %>% pull(year_f) %>% unique()

# write.csv(mean_dat, here::here("outputs", "data", "summary_stats.csv"),
#           row.names = FALSE)


# add in averages for years that are missing individual samples
dat_avg <- read.table(here::here("data", "nasscenturyavgv3FLAT.txt"),
                      header = TRUE) %>%
  pivot_longer(cols = -(iy.)) %>%
  mutate(
    year_f = as.factor(iy.),
    sex = ifelse(grepl("m", name), "male", "female") %>% as.factor,
    age = case_when(
      grepl("42", name) ~ "42",
      grepl("52", name) ~ "52",
      grepl("53", name) ~ "53",
      grepl("63", name) ~ "63"
    ) %>%
      as.factor(),
    period = case_when(
      iy. < 1957 ~ "Gilbert-Clemens",
      iy. > 1972 & iy. < 1994 ~ "Monkley Dump",
      iy. > 1993 ~ "Nisga'a",
      TRUE ~ "Bilton"
    ),
    age_f = age,
    period = fct_relevel(as.factor(period),
                         "Gilbert-Clemens",
                         "Bilton",
                         "Monkley Dump",
                         "Nisga'a"),
    lo = NaN,
    up = NaN,
    data = "annual average"
  ) %>%
  filter(year_f %in% empty_yrs) %>% 
  select(year_f, mean_fl = value, period, lo, up, age, age_f, sex, data)
levels(dat_avg$age_f) <- c("4[2]", "5[2]", "5[3]", "6[3]")

mean_dat_plotting <- mean_dat %>% 
  mutate(data = "individual measurements") %>% 
  filter(sex == "female",
         !is.na(period)) %>% 
  select(year_f, mean_fl, period, lo, up, age_f, data) %>%
  rbind(., dat_avg %>% filter(sex == "female") %>% select(-age, -sex)) %>% 
  group_by(age_f) %>% 
  mutate(
    mean_fl = mean_fl / 10,
    lo = lo / 10,
    up = up / 10,
    overall_mean = mean(mean_fl, na.rm = T)) %>% 
  ungroup() %>% 
  droplevels() 

annual_dot <- ggplot(
  mean_dat_plotting, 
  aes(x = as.numeric(as.character(year_f)))
) +
  geom_pointrange(aes(y = mean_fl, fill = data, ymin = lo, ymax = up, 
                      shape = period), size = 0.3) +
  facet_wrap(~age_f, labeller = label_parsed) +
  ggsidekick::theme_sleek() +
  labs(x = "Year", y = "Fork Length") +
  scale_x_continuous(
    breaks = seq(1915, 2015, by = 20),
    expand = c(0.02, 0.02)
  ) +
  geom_hline(aes(yintercept = overall_mean), lty = 2) +
  scale_shape_manual(values = shape_pal, name = "Sampling\nPeriod") +
  scale_fill_manual(values = fill_pal, name = "Dataset") +
  guides(fill = guide_legend(override.aes = list(shape = 21)))


# annual_sd_dot <- ggplot(mean_dat, aes(x = year_f)) +
#   geom_point(aes(y = sd_fl, fill = period), shape = 21) +
#   facet_grid(~age) +
#   ggsidekick::theme_sleek() +
#   labs(x = "Year", y = "SD Fork Length") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# ggplot(mean_dat, aes(x = year_f)) +
#   geom_point(aes(y = mean_yday, fill = period), shape = 21) +
#   facet_wrap(~age) +
#   ggsidekick::theme_sleek() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# ggplot(mean_dat, aes(x = year_f)) +
#   geom_point(aes(y = sd_yday, fill = period), shape = 21) +
#   facet_wrap(~age) +
#   ggsidekick::theme_sleek() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# export
png(here::here("outputs", "figs", "annual_dot.png"), width = 8, height = 5,
    res = 250, units = "in")
annual_dot
dev.off()

# png(here::here("outputs", "figs", "annual_box.png"), width = 11, height = 7,
#     res = 250, units = "in")
# annual_box
# dev.off()


# FIT MODEL  -------------------------------------------------------------------

# dat_trim <- dat %>% sample_n(30000)

# make fake mesh
dat$x <- runif(nrow(dat))
dat$y <- runif(nrow(dat))
dum_mesh <- make_mesh(dat, c("x", "y"), cutoff = 1000)

fit <- sdmTMB(fl_cm ~ 0 + s(yday_c, by = age, m = 2) +
                s(year, by = age, m = 2) +
                period + age + sex,
              dispformula = ~ 0 + period,
              data = dat,
              mesh = dum_mesh,
              spatial = "off",
              spatiotemporal = "off",
              control = sdmTMBcontrol(
                nlminb_loops = 2,
                newton_loops = 2
              ))



# check residuals
sims <- simulate(fit, nsim = 500)
sims %>% 
  dharma_residuals(fit)
# looks good


dat$resids <- resid(fit, type = "response")
dat$resids2 <- resid(fit_gam, type = "response")
trim_dat <- dat %>% 
  sample_n(30000) %>% 
  filter(sex == "male") %>% 
  pivot_longer(cols = c(resids, resids2), names_to = "resids_type")
ggplot(trim_dat) +
  geom_point(aes(x = year_f, y = value, fill = period), shape = 21) +
  facet_grid(resids_type ~ age_f) +
  ggsidekick::theme_sleek()


# COMPARE OUT OF SAMPLE PERFORMANCE --------------------------------------------

# NOTE this comparison was made with identical constraints (k cannot be greater 
# than 5 without leading to convergence issues for global smooth model)

# fit <- sdmTMB(fl_cm ~ 0 + s(yday_c, by = age, m = 2, k = 5) +
#                 s(year, by = age, m = 2, k = 5) +
#                 period + age + sex,
#               dispformula = ~ 0 + period,
#               data = dat,
#               mesh = dum_mesh,
#               spatial = "off",
#               spatiotemporal = "off",
#               control = sdmTMBcontrol(
#                 nlminb_loops = 2,
#                 newton_loops = 2
#               ))

# fit2 <- sdmTMB(fl_cm ~ 0 + s(yday_c, m = 2, k = 5) + 
#                 s(yday_c, by = age, m = 1, k = 5) +
#                 s(year, m = 2, k = 5) +
#                 s(year, by = age, m = 1, k = 5) +
#                 period + age + sex,
#               dispformula = ~ 0 + period,
#               data = dat,
#               mesh = dum_mesh,
#               spatial = "off",
#               spatiotemporal = "off",
#               control = sdmTMBcontrol(
#                 nlminb_loops = 2,
#                 newton_loops = 2
#               ))


# use average measures not used to fit model to determine whether global smooths 
# should be used

# adjust average data for missing years to match predictions format
# dat_avg$fl_cm <- dat_avg$mean_fl / 10
# dat_avg$year <- as.numeric(as.character(dat_avg$year_f))
# dat_avg$yday_c <- 0
# 
# group_pred <- predict(fit, newdata = dat_avg, re_form = NA)
# global_pred <- predict(fit2, newdata = dat_avg, re_form = NA)
# 
# Metrics::rmse(global_pred$est, dat_avg$fl_cm)
# Metrics::rmse(group_pred$est, dat_avg$fl_cm)
# group specific smooths have better predictions


# CATEGORICAL PREDICTIONS  -----------------------------------------------------

new_dat <- expand.grid(
  age = unique(dat$age),
  sex = unique(dat$sex),
  period = unique(dat$period),
  yday_c = 0,
  year = 1969,#mean(dat$year),
  # dummy spatial variables required 
  x = runif(1),
  y = runif(1)
) %>% 
  mutate(
    age_f = as.factor(age)
  )
levels(new_dat$age_f) <- c("4[2]", "5[2]", "5[3]", "6[3]")


fe_preds <- predict(fit, newdata = new_dat, re_form = NA, se_fit = TRUE)
fe_preds$low <- fe_preds$est + (qnorm(0.025) * fe_preds$est_se)
fe_preds$up <- fe_preds$est + (qnorm(0.975) * fe_preds$est_se)

period_eff_dot <- fe_preds %>% 
  filter(age == "42",
         sex == "male") %>% 
  ggplot(.) +
  geom_pointrange(aes(x = period, y = est, shape = period, 
                      ymin = low, ymax = up),
                  fill = "white") +
  scale_shape_manual(values = shape_pal) +
  ggsidekick::theme_sleek() +
  labs(x = "Sampling Period", y = "Estimated Mean") +
  theme(legend.position = "none")


age_eff_dot <- fe_preds %>% 
  filter(period == "Gilbert-Clemens",
         sex == "male") %>% 
  ggplot(.) +
  geom_pointrange(aes(x = age, y = est, ymin = low, ymax = up, fill = age),
                  shape = 21
  ) +
  ggsidekick::theme_sleek() +
  scale_x_discrete("Age Class", 
                   labels = parse(text = as.character(levels(fe_preds$age_f)))) +
  labs(x = "Age Class", y = "Estimated Mean") +
  scale_fill_brewer(palette = 1) +
  theme(legend.position = "none") 


sex_eff_dot <- fe_preds %>% 
  filter(period == "Gilbert-Clemens",
         age == "42") %>% 
  ggplot(.) +
  geom_pointrange(aes(x = sex, y = est, fill = sex, ymin = low, ymax = up),
                  shape = 21) +
  ggsidekick::theme_sleek() +
  labs(x = "Sex", y = "Estimated Mean") +
  scale_fill_brewer(palette = 5, type = "seq") +
  theme(legend.position = "none")


# sigma estimates
period_sig <- data.frame(
  sig_est = fit$sd_report$par.fixed[names(fit$sd_report$par.fixed) == "b_disp_k"],
  sig_se = sqrt(diag(fit$sd_report$cov.fixed[rownames(fit$sd_report$cov.fixed) == "b_disp_k", 
                                              colnames(fit$sd_report$cov.fixed) == "b_disp_k"])),
  period = unique(dat$period)
) %>% 
  mutate(
    low = sig_est + (qnorm(0.025) * sig_se),
    up = sig_est + (qnorm(0.975) * sig_se)
  )

period_sig_dot <- ggplot(period_sig) +
  geom_pointrange(aes(x = period, y = sig_est, shape = period, 
                      ymin = low, ymax = up),
                  fill = "white") +
  scale_shape_manual(values = shape_pal) +
  ggsidekick::theme_sleek() +
  labs(x = "Sampling Period", y = "Estimated Variance") +
  theme(legend.position = "none")


# combine all fixed effects figs
png(here::here("outputs", "figs", "main_effect.png"), 
    height = 5, width = 8,
    units = "in", res = 250)
cowplot::plot_grid(period_eff_dot,
                   age_eff_dot,
                   sex_eff_dot,
                   period_sig_dot, ncol = 2)
dev.off()


# SMOOTH PREDICTIONS  ----------------------------------------------------------

# year effects (assuming fixed period)
new_dat2 <- expand.grid(
  age = unique(dat$age),
  sex = "female",
  period = "Bilton",
  yday_c = 0,
  year = seq(min(dat$year), max(dat$year), length = 100),
  # dummy spatial variables required 
  x = runif(1),
  y = runif(1)
) %>% 
  mutate(
    age_f = as.factor(age)
  )
levels(new_dat2$age_f) <- c("4[2]", "5[2]", "5[3]", "6[3]")

# smoothed predictions
smooth_preds <- predict(fit, newdata = new_dat2, re_form = NA, se_fit = TRUE)
smooth_preds$low <- smooth_preds$est + (qnorm(0.025) * smooth_preds$est_se)
smooth_preds$up <- smooth_preds$est + (qnorm(0.975) * smooth_preds$est_se)


smooth_year <- ggplot(smooth_preds, aes(x = year, y = est)) +
  geom_line() +
  geom_ribbon(aes(ymin = low, ymax = up), alpha = 0.4) +
  ggsidekick::theme_sleek() +
  labs(x = "Year", y = "Fork Length (cm)") +
  facet_wrap(~age_f, labeller = label_parsed) +
  scale_x_continuous(
    breaks = seq(1915, 2015, by = 20),
    expand = c(0.02, 0.02)
  )

png(here::here("outputs", "figs", "smooth_means_no_global_unconstrainedk.png"), 
    height = 5, width = 8.5,
    units = "in", res = 250)
smooth_year
dev.off()


# smoothed predictions with GAM (just to check)
# fit_gam <- mgcv::gam(fl_cm ~ 0  + s(yday_c, by = age, m = 2, k = 6) +
#                        s(year, by = age, m = 2, k = 6) +
#                        period + age + sex,
#                      data = dat)
#
# smooth_preds_gam_vec <- predict(fit_gam, newdata = new_dat2, se.fit = TRUE)
# smooth_preds_gam <- new_dat2
# smooth_preds_gam$est <- smooth_preds_gam_vec$fit
# smooth_preds_gam$low <- smooth_preds_gam_vec$fit + (qnorm(0.025) * smooth_preds_gam_vec$se.fit)
# smooth_preds_gam$up <- smooth_preds_gam_vec$fit + (qnorm(0.975) * smooth_preds_gam_vec$se.fit)
# 
# ggplot(smooth_preds_gam, aes(x = year, y = est)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = low, ymax = up), alpha = 0.4) +
#   ggsidekick::theme_sleek() +
#   labs(x = "Year", y = "Fork Length (cm)") +
#   facet_wrap(~age_f, labeller = label_parsed) +
#   scale_x_continuous(
#     breaks = seq(1915, 2015, by = 20),
#     expand = c(0.02, 0.02)
#   )


# year day effects
new_dat3 <- expand.grid(
  age = unique(dat$age),
  sex = "female",
  period = "Gilbert-Clemens",
  yday_c = seq(-47, 65, length = 100),
  year = 1969,
  # dummy spatial variables required 
  x = runif(1),
  y = runif(1)
) %>% 
  mutate(
    age_f = as.factor(age)
  )
levels(new_dat2$age_f) <- c("4[2]", "5[2]", "5[3]", "6[3]")


smooth_preds2 <- predict(fit, newdata = new_dat3, re_form = NA, se_fit = TRUE)
smooth_preds2$low <- smooth_preds2$est + (qnorm(0.025) * smooth_preds2$est_se)
smooth_preds2$up <- smooth_preds2$est + (qnorm(0.975) * smooth_preds2$est_se)

smooth_yday <- ggplot(smooth_preds2, 
       aes(x = yday_c, est)) +
  geom_line() +
  geom_ribbon(aes(ymin = low, ymax = up), alpha = 0.4) +
  ggsidekick::theme_sleek() +
  labs(x = "Year", y = "Fork Length (cm)") +
  facet_wrap(~age_f, labeller = label_parsed) +
  scale_x_continuous(
    breaks = seq(1915, 2015, by = 20),
    expand = c(0.02, 0.02)
  )

png(here::here("outputs", "figs", "yday_smooth_no_global.png"), 
    height = 5, width = 8.5,
    units = "in", res = 250)
smooth_yday
dev.off()


# year effects (assuming variable period)
period_sig <- data.frame(
  sig_est = fit$sd_report$par.fixed[names(fit$sd_report$par.fixed) == "b_disp_k"],
  period = unique(dat$period)
)

# period specific sampling day (varies by 3 weeks)
# samp_day <- dat %>% 
#   group_by(period, age_f) %>% 
#   summarize(yday_c = mean(yday_c))

new_dat4 <- expand.grid(
  age = unique(dat$age),
  sex = unique(dat$sex),
  yday_c = 0,
  year = seq(min(dat$year), max(dat$year), by = 1),
  # dummy spatial variables required 
  x = runif(1),
  y = runif(1)
) %>% 
  mutate(
    period = case_when(
      year < 1957 ~ "Gilbert-Clemens",
      year > 1972 & year < 1994 ~ "Monkley Dump",
      year > 1993 ~ "Nisga'a",
      TRUE ~ "Bilton"
    ),
    period = fct_relevel(as.factor(period), 
                         "Gilbert-Clemens",
                         "Bilton",
                         "Monkley Dump",
                         "Nisga'a"),
    age_f = as.factor(age)
    ) #%>% 
  # left_join(., period_sig, by = "period") 
levels(new_dat4$age_f) <- c("4[2]", "5[2]", "5[3]", "6[3]")
# new_dat4 <- left_join(new_dat4, samp_day, by = c("period", "age_f"))


## ORIGINAL PREDICTION INTERVALS (DEFUNCT)

# smooth_fit_per <- predict(fit, newdata = new_dat4, re_form = NA)

# calculate prediction interval
# (from https://rpubs.com/aaronsc32/regression-confidence-prediction-intervals)
# y.fit <- predict(fit_gam, re_form = NA)
# n <- length(dat$fl_cm)
# 
# # sse <- sum((dat$fl_cm - y.fit$est)^2)
# sse <- sum((dat$fl_cm - y.fit)^2)
# mse <- sse / (n - 2)
# 
# t.val <- qt(0.975, n - 2) # Critical value of t
# 
# pred.x <- smooth_fit_per$year
# pred.y <- smooth_fit_per$est
# 
# pred.se.fit <- (1 + (1 / n) + (pred.x - mean(dat$year))^2 / (sum((dat$year - mean(dat$year))^2))) # Standard error of the prediction
# 
# smooth_fit_per$pred_up <- pred.y + t.val * sqrt(mse * pred.se.fit)
# smooth_fit_per$pred_low <- pred.y - t.val * sqrt(mse * pred.se.fit)


## Simulate predictions based on estimated parameters
# 1) generate model matrix based on predictors and formula
# 2) draw parameter values based on estimated values (and MVN dist)
# 3) generate new predictions (excluding residual variability)
# 4) add random deviates based on sigma
# https://fromthebottomoftheheap.net/2014/06/16/simultaneous-confidence-intervals-for-derivatives/


source(here::here("R", "utils.R"))

## generate inputs for simulation (functions saved in utils.R; from sdmTMB)
formula_no_sm <- remove_s_and_t2(fit$formula[[1]])
X_ij <- model.matrix(formula_no_sm, data = dat)
sm <- parse_smoothers(fit$formula[[1]], data = dat)
pred_X_ij <- predict(mgcv::gam(formula_no_sm, data = dat), 
                     new_dat4, type = "lpmatrix")
sm_pred <- parse_smoothers(fit$formula[[1]], data = dat, 
                           newdata = new_dat4)
sm_start <- sm$b_smooth_start

# model matrices
pred_X1_ij = pred_X_ij
pred_Zs = sm_pred$Zs
pred_Xs = sm_pred$Xs

# cov matrix and parameter estimates concatenated to same order as matrix rownames
cov <- fit$sd_report$cov.fixed
coef_list <- as.list(fit$sd_report, "Estimate")
coef_vec <- c(coef_list$b_j, coef_list$b_disp_k, coef_list$bs, 
              coef_list$ln_smooth_sigma)

set.seed(1234)
n_sims <- 1000
cov_sim <- MASS::mvrnorm(n_sims, coef_vec, cov)
b1_j <- cov_sim[ , grep("b_j", colnames(cov))] #fixed effect pars
bs <- cov_sim[ , grep("bs", colnames(cov))] # fixed effect smooths
ssdr <- summary(fit$sd_report) 
b_smooth <- ssdr[rownames(ssdr) == "b_smooth", "Estimate"] %>% as.numeric() #random smooths

# save sigmas as separate dataframe that can be joined 
sigmas <- cov_sim[ , grep("b_disp_k", colnames(cov))] # residual variance
colnames(sigmas) <- paste("sigma_", seq(1, ncol(sigmas)), sep = "")
sigma_dat <- sigmas %>% 
  as.data.frame() %>% 
  mutate(
    iter = seq(1, nrow(sigmas), by = 1)
  ) %>% 
  pivot_longer(
    cols = -iter,
    names_to = "sigma_p",
    values_to = "sigma_est"
  ) %>% 
  mutate(
    period = factor(sigma_p, labels = levels(dat$period))
  ) 


# generate predictions (based on code adapted from sdmTMB.cpp)
sim_list <- vector(mode = "list", length = n_sims)

for (i in seq_len(n_sims)) {
  pred_mu1 <- pred_X1_ij %*% b1_j[i, ]
  
  pred_smooth <- matrix(NA, nrow(pred_X1_ij), length(sm_start))
  
  for (ss in seq_along(sm_start)) {
    beta_ss <- rep(NA, times = ncol(pred_Zs[[ss]]))
    for (j in seq_along(beta_ss)) {
      beta_ss[j] <- b_smooth[sm_start[ss] + j]
    }
    pred_smooth[ , ss] <- pred_Zs[[ss]] %*% beta_ss
  }
  
  pred_smooth1_i <- apply(pred_smooth, 1, sum)
  pred_smooth2_i <- pred_Xs %*% bs[i, ]
  pred_smooth_i <- pred_smooth1_i + pred_smooth2_i
  pred_mu_i <- pred_mu1 + pred_smooth_i
  
  sim_list[[i]] <- new_dat4 %>% 
    mutate(
      iter = i,
      pred_mu = as.numeric(pred_mu_i)
    ) %>% 
    left_join(
      ., sigma_dat, by = c("iter", "period")
    ) %>% 
    mutate(
      sim_obs = rnorm(pred_mu, pred_mu, sqrt(sigma_est))
    )
}

sim_dat <- sim_list %>% 
  bind_rows() %>% 
  group_by(age, sex, year, period, age_f) %>% 
  summarize(
    up_pred_mu = quantile(pred_mu, 0.975),
    low_pred_mu = quantile(pred_mu, 0.025),
    pred_mu = mean(pred_mu),
    up_sim_obs = quantile(sim_obs, 0.975),
    low_sim_obs = quantile(sim_obs, 0.025),
    sim_obs = mean(sim_obs),
    .groups = "drop"
  ) 

saveRDS(sim_dat, here::here("outputs", "data", "tmb_post_preds.rds"))

# sim_list %>% 
#   bind_rows %>% 
#   filter(
#     sex == "female",
#     iter %in% seq(1, 10, by = 1)
#   ) %>% 
#   ggplot(., 
#          aes(x = year)) +
#   geom_line(aes(y = sim_obs, colour = as.factor(iter))) +
#   facet_wrap(~age_f)
#   glimpse()

# include average predictions for context
dat_avg_trim <- dat_avg %>% 
  mutate(
    fl_cm = mean_fl / 10,
    year = as.numeric(as.character(year_f)),
    yday_c = 0
  ) %>% 
  filter(sex == "female") %>% 
  select(
    year, fl_cm, dataset = data, period, age_f
  )
  
# join with random subsample of predictions
# point_dat <- dat %>% 
#   mutate(dataset = "individual measurements") %>% 
#   filter(sex == "female") %>% 
#   select(colnames(dat_avg_trim)) %>% 
#   sample_n(1000) %>% 
#   rbind(., dat_avg_trim) 

smooth_year_per <- ggplot(sim_dat %>% filter(sex == "female"), 
       aes(x = year)) +
  geom_line(aes(y = sim_obs)) +
  geom_ribbon(aes(ymin = low_sim_obs, ymax = up_sim_obs), alpha = 0.3) +
  geom_point(data = dat_avg_trim,
             aes(y = fl_cm, shape = period), fill = "red") +
  scale_shape_manual(values = shape_pal, name = "Sampling\nPeriod") +
  # scale_fill_manual(values = fill_pal, name = "Dataset") +
  ggsidekick::theme_sleek() +
  labs(x = "Year", y = "Fork Length (cm)") +
  facet_wrap(~age_f, labeller = label_parsed) +
  scale_x_continuous(
    breaks = seq(1915, 2015, by = 20),
    expand = c(0.02, 0.02)
  )

png(here::here("outputs", "figs", "year_smooth_period_no_global.png"), 
    height = 5, width = 8.5,
    units = "in", res = 250)
smooth_year_per
dev.off()


# SUMMARY STATISTICS -----------------------------------------------------------

# differences between final year and a) time series average or b) beginning of 
# time series for each age; control and don't control for period effects

smooth_list <- split(smooth_preds, smooth_preds$age_f) 
smooth_period_list <- split(smooth_fit_per, smooth_fit_per$age_f) 
age_names <- c("53", "52", "42", "63")

purrr::map2(smooth_list, age_names, function (x, y) {
  ts <- x$est
  #calc time series mean for each iteration
  ts_mean = mean(ts)
  # #difference from mean
  # diff <- ts[length(ts)] - ts_mean
  # #difference from first obs
  # diff2 <- ts[length(ts)] - ts[1]
  data.frame(
    mean_diff = ts[length(ts)] - ts_mean,
    first_diff = ts[length(ts)] - ts[1],
    mean_rel_diff = (ts[length(ts)] - ts_mean) / ts_mean,
    first_rel_diff = (ts[length(ts)] - ts[1]) / ts[1],
    # low = quantile(diff2, probs = 0.05),
    # high = quantile(diff2, probs = 0.95),
    age = y
  )
}) %>% 
  bind_rows

purrr::map2(smooth_period_list, age_names, function (x, y) {
  ts <- x$est
  #calc time series mean for each iteration
  ts_mean = mean(ts)
  data.frame(
    mean_diff = ts[length(ts)] - ts_mean,
    first_diff = ts[length(ts)] - ts[1],
    age = y
  )
}) %>% 
  bind_rows


## DIFFERENCE IN FECUNDITY -----------------------------------------------------

fec_dat <- expand.grid(
  age = unique(dat$age),
  sex = "female",
  period = "Nisga'a",
  yday_c = 0,
  year = seq(min(dat$year), max(dat$year), by = 1),
  # dummy spatial variables required 
  x = runif(1),
  y = runif(1)
) %>% 
  mutate(
    age_f = as.factor(age)
  )
levels(fec_dat$age_f) <- c("4[2]", "5[2]", "5[3]", "6[3]")


# simulate differences in fecundity assuming sigma equal to estimate for Nisga'a
# sampling period
fec_preds <- predict(fit, newdata = fec_dat, re_form = NA)
nsims <- 1000
ndfw_sigma <- 1.335
# Mason and West regression for Fulton River
b1 = 11.52 #11.10 (Pinkut Creek) 
a1 = -2152 #2015.8 (Pinkut Creek)
# sigma = 14.9 #19.5 labelled as s_yx in text, not sure if this is SD or something else...

sim_list <- vector(mode = "list", length = nsims)


set.seed(1234)

for (i in 1:nsims) {
  sim_fl <- fec_preds$est + rnorm(nrow(fec_preds), 0, ndfw_sigma)
  # convert to POH length in mm using inverse of equation in main manuscript
  sim_hl <-  (0.833 * sim_fl * 10) - 3.508
  # convert to fecundity (1000s) using west and mason regression pars
  sim_list[[i]] <- data.frame( 
    iter = i,
    age  = fec_preds$age,
    year = fec_preds$year,
    sim_fec = (b1 * sim_hl + a1) / 1000
    ) %>% 
    group_by(age) %>% 
    mutate(
      mean_sim_fec = mean(sim_fec)
    ) %>% 
    filter(
      year %in% c(min(dat$year), max(dat$year))
    ) %>% 
    pivot_wider(names_from = "year", values_from = "sim_fec") %>% 
    mutate(diff = `2019` - `1914`,
           rel_diff = diff / mean_sim_fec)
}
sim_dat <- data.table::rbindlist(sim_list) %>% 
  mutate(age_f = as.factor(age)) %>% 
  as.data.frame()
levels(sim_dat$age_f) <- c("4[2]", "5[2]", "5[3]", "6[3]")


diff_fecundity <- sim_dat %>% 
  group_by(age_f) %>% 
  summarize(
    median = median(diff),
    low = quantile(diff, probs = 0.05),
    up = quantile(diff, probs = 0.95),
    rel_median = median(rel_diff),
    rel_low = quantile(rel_diff, probs = 0.05),
    rel_up = quantile(rel_diff, probs = 0.95),
    .groups = "drop"
  ) 


png(here::here("outputs", "figs", "fecundity_diff.png"), 
    height = 4, width = 5, units = "in", res = 250)
ggplot(diff_fecundity) +
  geom_pointrange(aes(x = age_f, y = rel_median, ymin = rel_low, 
                      ymax = rel_up)) +
  geom_hline(aes(yintercept = 0), lty = 2) +
  ggsidekick::theme_sleek() +
  # labs(y = "Percent Decline in Fecundity (2019-1914)") +
  scale_x_discrete(
    "Age",
    labels = c(expression("4"["2"]), expression("5"["2"]), 
               expression("5"["3"]), expression("6"["3"]))
  ) +
  scale_y_continuous(
    "Percent Decline in Fecundity (2019-1914)",
    labels = c("-30%", "-20%", "-10%", "0%")
  ) +
  theme(legend.position = "none")
dev.off()



