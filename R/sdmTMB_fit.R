# Fit sdmTMB model to century of Nass data 
# Faster alternative to brms_fit while still accounting for variable dispersion
# 1) Import Skip's first century of data and clean as necessary
# 2) Estimate temporal trend using GAM and sdmTMB to account for changes in
# variability with sampling regime
# Created by C Freshwater Dec 23, 2021
# Update w/ edited data
# -------------------------------------------------

# install specific branch of sdmTMB that includes dispersion formula
devtools::install_github("https://github.com/pbs-assess/sdmTMB",
                         ref = "dispformula2")


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
      year > 1993 ~ "NFWD",
      TRUE ~ "Bilton"
    ),
    period = fct_relevel(as.factor(period), 
                         "Gilbert-Clemens",
                         "Bilton",
                         "Monkley Dump",
                         "NFWD"),
    age_sex = paste(sex, age, sep = "_") %>% as.factor()
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
      iy. > 1993 ~ "NFWD",
      TRUE ~ "Bilton"
    ),
    age_f = age,
    period = fct_relevel(as.factor(period),
                         "Gilbert-Clemens",
                         "Bilton",
                         "Monkley Dump",
                         "NFWD"),
    lo = NaN,
    up = NaN,
    data = "annual average"
  ) %>%
  filter(year_f %in% empty_yrs) %>% 
  select(year_f, mean_fl = value, period, lo, up, age, age_f, sex, data)
levels(dat_avg$age_f) <- c("4[2]", "5[2]", "5[3]", "6[3]")

mean_dat_plotting <- mean_dat %>% 
  mutate(data = "individual measurements") %>% 
  filter(!is.na(period)) %>%
  select(year_f, mean_fl, period, lo, up, age, age_f, sex, data) %>%
  rbind(., dat_avg) %>% 
  group_by(age, age_f, sex) %>% 
  mutate(
    overall_mean = mean(mean_fl, na.rm = T)) %>% 
  ungroup() %>% 
  droplevels() 

annual_dot <- ggplot(
  mean_dat_plotting %>% filter(sex == "female"), 
  aes(x = as.numeric(as.character(year_f)))
) +
  geom_pointrange(aes(y = mean_fl, fill = data, ymin = lo, ymax = up, 
                      shape = period), size = 0.3) +
  facet_wrap(~age_f, labeller = label_parsed) +
  ggsidekick::theme_sleek() +
  labs(x = "Year", y = "Fork Length (mm)") +
  scale_x_continuous(
    breaks = seq(1915, 2015, by = 20),
    expand = c(0.02, 0.02)
  ) +
  geom_hline(aes(yintercept = overall_mean), lty = 2) +
  scale_shape_manual(values = shape_pal, name = "Sampling\nPeriod") +
  scale_fill_manual(values = fill_pal, name = "Dataset") +
  guides(fill = guide_legend(override.aes = list(shape = 21)))


# export
png(here::here("outputs", "figs", "annual_dot.png"), width = 8, height = 5,
    res = 250, units = "in")
annual_dot
dev.off()



# AGE COMPOSITION --------------------------------------------------------------

comp_dat <-  dat_in %>%
  mutate(
    year_f = as.factor(year),
    age_f = age,
    period = case_when(
      year < 1957 ~ "Gilbert-Clemens",
      year > 1972 & year < 1994 ~ "Monkley Dump",
      year > 1993 ~ "NFWD",
      TRUE ~ "Bilton"
    ),
    period = fct_relevel(as.factor(period), 
                         "Gilbert-Clemens",
                         "Bilton",
                         "Monkley Dump",
                         "NFWD"),
    dominant_age = ifelse(
      age %in% c("42", "52", "53", "63"),
      TRUE,
      FALSE
    )
  ) %>%  
  filter(dominant_age == TRUE) %>% 
  group_by(year) %>% 
  mutate(
    nn = n()
  ) %>%
  group_by(dominant_age, age_f, year, period) %>% 
  summarize(
    ppn = n() / nn,
    .groups = "drop"
  ) %>% 
  ungroup %>%
  distinct() %>% 
  arrange(year) 

comp_dot <- ggplot(comp_dat) +
  geom_point(aes(x = year, y = ppn, shape = period)) +
  scale_shape_manual(values = shape_pal) +
  facet_wrap(~age_f, labeller = label_parsed) +
  ggsidekick::theme_sleek() +
  labs(x = "Year", y = "Proportion of Sample") 

png(here::here("outputs", "figs", "annual_dot_comp.png"), width = 8, height = 5,
    res = 250, units = "in")
comp_dot
dev.off()



# FIT MODEL  -------------------------------------------------------------------

fit <- sdmTMB(fl ~ s(yday_c, by = age_sex, m = 2) +
                s(year, by = age_sex, m = 2) +
                age + sex,
              dispformula = ~ 0 + period,
              data = dat,
              spatial = "off",
              spatiotemporal = "off",
              control = sdmTMBcontrol(
                # nlminb_loops = 2,
                newton_loops = 1
              ),
              silent = FALSE)
sanity(fit)

# check residuals
sims <- simulate(fit, nsim = 250)
dharma_sims <- sims %>% 
  dharma_residuals(fit)
# looks good

resid_mat <- sims - dat$fl

resid_long <- resid_mat %>% 
  as.data.frame() %>% 
  cbind(dat %>% select(year, period), .)  %>% 
  pivot_longer(-c(year, period), names_to = "iter", values_to = "resid_est") 

resid_dat <- resid_long %>% 
  group_by(year) %>% 
  summarize(mean_resid = mean(resid_est),
            up_resid = quantile(resid_est, 0.975),
            low_resid = quantile(resid_est, 0.025))

transition_years <- mean_dat_plotting %>% 
  filter(!period == "Gilbert-Clemens") %>% 
  mutate(year = as.numeric(as.character(year_f))) %>% 
  group_by(period) %>% 
  summarize(first_year = min(year)) %>% 
  pull(first_year)

png(here::here("outputs", "figs", "gam_resid_ts.png"), width = 8, height = 5,
    res = 250, units = "in")
ggplot(resid_dat) +
  geom_pointrange(aes(x = year, y = mean_resid, ymin = low_resid, ymax = up_resid)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = c(transition_years), color = "red") +
  labs(x = "Year", y = "GAM Residuals") +
  ggsidekick::theme_sleek()
dev.off()

png(here::here("outputs", "figs", "gam_resid_bp.png"), width = 8, height = 5,
    res = 250, units = "in")
ggplot(resid_long) +
  geom_boxplot(aes(x = period, y = resid_est)) +
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Period", y = "State Residuals") +
  ggsidekick::theme_sleek()
dev.off()



# CATEGORICAL PREDICTIONS  -----------------------------------------------------

# replaced with effect sizes 
new_dat <- expand.grid(
  age = unique(dat$age),
  sex = unique(dat$sex),
  yday_c = 0,
  year = 1969
  ) %>%
  mutate(
    age_f = as.factor(age),
    age_sex = paste(sex, age, sep = "_") %>% as.factor()
  )
levels(new_dat$age_f) <- c("4[2]", "5[2]", "5[3]", "6[3]")


fe_preds <- predict(fit, newdata = new_dat, re_form = NA, se_fit = TRUE)
fe_preds$low <- fe_preds$est + (qnorm(0.025) * fe_preds$est_se)
fe_preds$up <- fe_preds$est + (qnorm(0.975) * fe_preds$est_se)

age_eff_dot <- fe_preds %>%
  filter(sex == "male") %>%
  ggplot(.) +
  geom_pointrange(aes(x = age, y = est, ymin = low, ymax = up, fill = age),
                  shape = 21
  ) +
  ggsidekick::theme_sleek() +
  scale_x_discrete("Age Class",
                   labels = parse(text = as.character(levels(fe_preds$age_f)))) +
  labs(x = "Age Class", y = "Estimated Mean Fork Length (mm)") +
  scale_fill_brewer(palette = 1) +
  theme(legend.position = "none")


sex_eff_dot <- fe_preds %>%
  filter(age == "42") %>%
  ggplot(.) +
  geom_pointrange(aes(x = sex, y = est, fill = sex, ymin = low, ymax = up),
                  shape = 21) +
  ggsidekick::theme_sleek() +
  labs(x = "Sex", y = "Estimated Mean Fork Length (mm)") +
  scale_fill_brewer(palette = 5, type = "seq") +
  theme(legend.position = "none")

# sigma estimates
# note that model estimates are sd (not variance = sd^2)
period_sig <- data.frame(
  sig_est = as.list(fit$sd_report, "Estimate")$b_disp_k,
  sig_se = as.list(fit$sd_report, "Std. Error")$b_disp_k,
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
  labs(x = "Sampling Period", y = "Estimated Residual SD") +
  theme(legend.position = "none")

png(here::here("outputs", "figs", "sigma_ests.png"), 
    height = 4, width = 5,
    units = "in", res = 250)
period_sig_dot
dev.off()


# combine all fixed effects figs
png(here::here("outputs", "figs", "main_effect.png"),
    height = 3, width = 6.5,
    units = "in", res = 250)
cowplot::plot_grid(age_eff_dot,
                   sex_eff_dot,
                   ncol = 2)
dev.off()


# SMOOTH PREDICTIONS  ----------------------------------------------------------

# year effects (assuming fixed period)
new_dat2 <- expand.grid(
  age = unique(dat$age),
  sex = unique(dat$sex),#"female",
  yday_c = 0,
  year = seq(min(dat$year), max(dat$year), length = 100),
  # dummy spatial variables required 
  x = runif(1),
  y = runif(1)
) %>% 
  mutate(
    age_f = as.factor(age),
    age_sex = paste(sex, age, sep = "_") %>% as.factor()
  )
levels(new_dat2$age_f) <- c("4[2]", "5[2]", "5[3]", "6[3]")

# smoothed predictions
smooth_preds <- predict(fit, newdata = new_dat2, re_form = NA, se_fit = TRUE)
smooth_preds$low <- smooth_preds$est + (qnorm(0.025) * smooth_preds$est_se)
smooth_preds$up <- smooth_preds$est + (qnorm(0.975) * smooth_preds$est_se)


smooth_year <- ggplot(smooth_preds, aes(x = year, y = est)) +
  geom_line(aes(colour = sex)) +
  geom_ribbon(aes(ymin = low, ymax = up, fill = sex), alpha = 0.4) +
  ggsidekick::theme_sleek() +
  scale_color_brewer(type = "qual") +
  scale_fill_brewer(type = "qual") +
  labs(x = "Year", y = "Fork Length (mm)") +
  facet_wrap(~age_f, labeller = label_parsed) +
  scale_x_continuous(
    breaks = seq(1915, 2015, by = 20),
    expand = c(0.02, 0.02)
  ) 


png(here::here("outputs", "figs", "year_smooth.png"), 
    height = 5, width = 8.5,
    units = "in", res = 250)
smooth_year
dev.off()


# year day effects
new_dat3 <- expand.grid(
  age = unique(dat$age),
  sex = unique(dat$sex),#"female",
  period_b_cent = 0,
  period_m_cent = 0,
  period_n_cent = 0,
  yday_c = seq(-47, 65, length = 100),
  year = 1969,
  # dummy spatial variables required 
  x = runif(1),
  y = runif(1)
) %>% 
  mutate(
    age_f = as.factor(age),
    age_sex = paste(sex, age, sep = "_") %>% as.factor()
  )
levels(new_dat3$age_f) <- c("4[2]", "5[2]", "5[3]", "6[3]")

smooth_preds2 <- predict(fit, newdata = new_dat3, re_form = NA, se_fit = TRUE)
smooth_preds2$low <- smooth_preds2$est + (qnorm(0.025) * smooth_preds2$est_se)
smooth_preds2$up <- smooth_preds2$est + (qnorm(0.975) * smooth_preds2$est_se)

smooth_yday <- ggplot(smooth_preds2, 
       aes(x = yday_c, est)) +
  geom_line(aes(colour = sex)) +
  geom_ribbon(aes(ymin = low, ymax = up, fill = sex), alpha = 0.4) +
  ggsidekick::theme_sleek() +
  scale_color_brewer(type = "qual") +
  scale_fill_brewer(type = "qual") +
  labs(x = "Year Day", y = "Fork Length (mm)") +
  facet_wrap(~age_f, labeller = label_parsed) +
  scale_x_continuous(
    breaks = c(-23.5, 0, 26.4, 51.4),
    labels = c(175, 200, 225, 250),
    expand = c(0.02, 0.02)
  )

png(here::here("outputs", "figs", "yday_smooth.png"), 
    height = 5, width = 8.5,
    units = "in", res = 250)
smooth_yday
dev.off()


# OUT OF SAMPLE PREDICTIONS ----------------------------------------------------

dat_avg_trim <- dat_avg %>%
  mutate(
    year = as.numeric(as.character(year_f)),
    yday_c = 0) %>%
  select(
    year,
    mean_fl,
    age_f, sex#, age_sex
  ) %>%
  droplevels()

# year effects (assuming variable period)
period_sig <- data.frame(
  sig_est = fit$sd_report$par.fixed[names(fit$sd_report$par.fixed) == "b_disp_k"],
  period = unique(dat$period)
)

# period specific sampling day (varies by 3 weeks)
samp_day <- dat %>%
  group_by(period, age_f) %>%
  summarize(yday_c = mean(yday_c),
            .groups = "drop")


new_dat4 <- expand.grid(
  age = unique(dat$age),
  sex = unique(dat$sex),
  year = unique(dat_avg_trim$year)
) %>%
  mutate(
    period = case_when(
      year < 1957 ~ "Gilbert-Clemens",
      TRUE ~ "Bilton"
    ),
    period = fct_relevel(as.factor(period),
                         "Gilbert-Clemens",
                         "Bilton"),
    age_f = as.factor(age),
    age_sex = paste(sex, age, sep = "_") %>% as.factor()
    ) 
levels(new_dat4$age_f) <- c("4[2]", "5[2]", "5[3]", "6[3]")
new_dat4 <- new_dat4 %>% 
  left_join(., dat_avg_trim, by = c("sex", "age_f", "year")) %>% 
  left_join(., period_sig, by = c("period")) %>% 
  left_join(., samp_day, by = c("period", "age_f"))


# simulate predictions
# components for calculating smooths
source(here::here("R", "utils.R"))
formula_no_sm <- remove_s_and_t2(fit$formula[[1]])
X_ij <- model.matrix(formula_no_sm, data = dat)
sm <- parse_smoothers(fit$formula[[1]], data = dat)
sm_start <- sm$b_smooth_start
pred_X_ij_oos <- predict(mgcv::gam(formula_no_sm, data = dat), 
                         new_dat4, type = "lpmatrix")
sm_pred <- parse_smoothers(fit$formula[[1]], data = dat, 
                           newdata = new_dat4)
pred_X1_ij_oos = pred_X_ij_oos
pred_Zs_oos = sm_pred$Zs
pred_Xs_oos = sm_pred$Xs

# cov matrix and parameter estimates concatenated to same order as matrix rownames
cov <- fit$sd_report$cov.fixed
coef_list <- as.list(fit$sd_report, "Estimate")
coef_vec <- c(coef_list$b_j, coef_list$b_disp_k, coef_list$bs,
              coef_list$ln_smooth_sigma)


n_sims <- 1000
cov_sim <- MASS::mvrnorm(n_sims, coef_vec, cov)
b1_j <- cov_sim[ , grep("b_j", colnames(cov))] #fixed effect pars
bs <- cov_sim[ , grep("bs", colnames(cov))] # fixed effect smooths
ssdr <- summary(fit$sd_report)
b_smooth <- ssdr[rownames(ssdr) == "b_smooth", "Estimate"] %>% as.numeric() #random smooths

sim_list_oos <- vector(mode = "list", length = n_sims)

for (i in seq_len(n_sims)) {
  pred_mu1 <- pred_X1_ij_oos %*% b1_j[i, ]
  
  pred_smooth <- matrix(NA, nrow(pred_X1_ij_oos), length(sm_start))
  
  for (ss in seq_along(sm_start)) {
    beta_ss <- rep(NA, times = ncol(pred_Zs_oos[[ss]]))
    for (j in seq_along(beta_ss)) {
      beta_ss[j] <- b_smooth[sm_start[ss] + j]
    }
    pred_smooth[ , ss] <- pred_Zs_oos[[ss]] %*% beta_ss
  }
  
  pred_smooth1_i <- apply(pred_smooth, 1, sum)
  pred_smooth2_i <- pred_Xs_oos %*% bs[i, ]
  pred_smooth_i <- pred_smooth1_i + pred_smooth2_i
  pred_mu_i <- pred_mu1 + pred_smooth_i + 
    # add residual variance
    rnorm(nrow(new_dat4), mean = 0, new_dat4$sig_est)
  
  sim_list_oos[[i]] <- new_dat4 %>% 
    mutate(
      iter = i,
      pred_mu = as.numeric(pred_mu_i)
    )
}

sim_oos <- sim_list_oos %>%
  bind_rows() %>%
  rename(obs = mean_fl,
         # obs = fl_cm,
         est = pred_mu) %>% 
  group_by(age_f, sex, year, period, obs) %>% 
  summarize(low = quantile(est, 0.025),
            up = quantile(est, 0.975),
            est = mean(est),
            .groups = "drop") %>% 
  pivot_longer(cols = c(est, obs), names_to = "mean_size") %>% 
  mutate(
    low = ifelse(mean_size == "obs", NA, low),
    up = ifelse(mean_size == "obs", NA, up)
  ) 

oos_points <- ggplot(sim_oos) +
  geom_pointrange(aes(x = year, y = value, ymin = low, ymax = up,
                      fill = mean_size), shape = 21) +
  facet_grid(age_f~sex, scales = "free_y", labeller = label_parsed) +
  labs(x = "Year", y = "Fork Length") +
  ggsidekick::theme_sleek() +
  scale_x_continuous() +
  theme(legend.position = "none")

png(here::here("outputs", "figs", "oos_preds.png"), 
    height = 8.5, width = 7.5,
    units = "in", res = 250)
oos_points
dev.off()


# SUMMARY STATISTICS -----------------------------------------------------------

# differences between final year and a) time series average or b) beginning of 
# time series for each age; control and don't control for period effects

smooth_list <- split(smooth_preds, smooth_preds$age_sex) 
age_sex_names <- levels(smooth_preds$age_sex)

purrr::map2(smooth_list, age_sex_names, function (x, y) {
  ts <- x$est
  #calc time series mean for each iteration
  ts_mean = mean(ts)
  data.frame(
    mean_diff = ts[length(ts)] - ts_mean,
    first_diff = ts[length(ts)] - ts[1],
    mean_rel_diff = (ts[length(ts)] - ts_mean) / ts_mean,
    first_rel_diff = (ts[length(ts)] - ts[1]) / ts[1],
    age = y
  )
}) %>% 
  bind_rows


## DIFFERENCE IN FECUNDITY -----------------------------------------------------

fec_dat <- expand.grid(
  age = unique(dat$age),
  sex = "female",
  period = "NFWD",
  yday_c = 0,
  year = seq(min(dat$year), max(dat$year), by = 1)
) %>% 
  mutate(
    age_f = as.factor(age),
    age_sex = paste(sex, age, sep = "_") %>% as.factor()
  )
levels(fec_dat$age_f) <- c("4[2]", "5[2]", "5[3]", "6[3]")

# components for calculating smooths
source(here::here("R", "utils.R"))
formula_no_sm <- remove_s_and_t2(fit$formula[[1]])
X_ij <- model.matrix(formula_no_sm, data = dat)
sm <- parse_smoothers(fit$formula[[1]], data = dat)
sm_start <- sm$b_smooth_start
pred_X_ij_fec <- predict(mgcv::gam(formula_no_sm, data = dat), 
                     fec_dat, type = "lpmatrix")
sm_pred <- parse_smoothers(fit$formula[[1]], data = dat, 
                           newdata = fec_dat)
pred_X1_ij_fec = pred_X_ij_fec
pred_Zs_fec = sm_pred$Zs
pred_Xs_fec = sm_pred$Xs

# cov matrix and parameter estimates concatenated to same order as matrix rownames
cov <- fit$sd_report$cov.fixed
coef_list <- as.list(fit$sd_report, "Estimate")
coef_vec <- c(coef_list$b_j, coef_list$b_disp_k, coef_list$bs,
              coef_list$ln_smooth_sigma)


n_sims <- 1000
cov_sim <- MASS::mvrnorm(n_sims, coef_vec, cov)
b1_j <- cov_sim[ , grep("b_j", colnames(cov))] #fixed effect pars
bs <- cov_sim[ , grep("bs", colnames(cov))] # fixed effect smooths
ssdr <- summary(fit$sd_report)
b_smooth <- ssdr[rownames(ssdr) == "b_smooth", "Estimate"] %>% as.numeric() #random smooths

# simulate as above but assuming sigma in both years equal to NFWD value and
# include calculations; uses same coefficient samples as prediction intervals
# above
sim_list_fec <- vector(mode = "list", length = n_sims)

for (i in seq_len(n_sims)) {
  pred_mu1 <- pred_X1_ij_fec %*% b1_j[i, ]
  
  pred_smooth <- matrix(NA, nrow(pred_X1_ij_fec), length(sm_start))
  
  for (ss in seq_along(sm_start)) {
    beta_ss <- rep(NA, times = ncol(pred_Zs_fec[[ss]]))
    for (j in seq_along(beta_ss)) {
      beta_ss[j] <- b_smooth[sm_start[ss] + j]
    }
    pred_smooth[ , ss] <- pred_Zs_fec[[ss]] %*% beta_ss
  }
  
  pred_smooth1_i <- apply(pred_smooth, 1, sum)
  pred_smooth2_i <- pred_Xs_fec %*% bs[i, ]
  pred_smooth_i <- pred_smooth1_i + pred_smooth2_i
  pred_mu_i <- pred_mu1 + pred_smooth_i
  
  sim_list_fec[[i]] <- fec_dat %>% 
    mutate(
      iter = i,
      pred_mu = as.numeric(pred_mu_i)
    ) %>% 
    mutate(
      # convert to POH length in mm using inverse of equation in main manuscript
      sim_hl =  (0.833 * pred_mu) - 3.508,
      fec_beta = ifelse(age %in% c("42", "53"), 10.67, 9.77),
      fec_int = ifelse(age %in% c("42", "53"), 1811.9, 1244.7),
      sim_fec = (fec_beta * sim_hl - fec_int) / 1000
    ) %>% 
    group_by(age) %>%
    mutate(
      mean_sim_fec = mean(sim_fec)
    ) %>% 
    filter(
      year %in% c(min(dat$year), max(dat$year))
    ) %>% 
    select(year, age, age_f, iter, sim_fec, mean_sim_fec) %>% 
    pivot_wider(names_from = "year", values_from = "sim_fec") %>% 
    mutate(diff = `2019` - `1914`,
           first_diff = diff / `1914`,
           rel_diff = diff / mean_sim_fec) %>% 
    ungroup()
}

sim_dat_fec <- data.table::rbindlist(sim_list_fec) %>% 
  as.data.frame()


diff_fecundity <- sim_dat_fec %>% 
  group_by(age_f) %>% 
  summarize(
    mean = mean(diff),
    low = quantile(diff, probs = 0.025),
    up = quantile(diff, probs = 0.975),
    rel_mean = mean(first_diff),
    rel_low = quantile(first_diff, probs = 0.025),
    rel_up = quantile(first_diff, probs = 0.975),
    .groups = "drop"
  ) 


png(here::here("outputs", "figs", "fecundity_diff.png"), 
    height = 4, width = 5, units = "in", res = 250)
ggplot(diff_fecundity) +
  geom_pointrange(aes(x = age_f, y = rel_mean, ymin = rel_low, 
                      ymax = rel_up)) +
  geom_hline(aes(yintercept = 0), lty = 2) +
  ggsidekick::theme_sleek() +
  scale_x_discrete(
    "Age",
    labels = c(expression("4"["2"]), expression("5"["2"]), 
               expression("5"["3"]), expression("6"["3"]))
  ) +
  scale_y_continuous(
    "Percent Decline in Fecundity (2019-1914)",
    breaks = c(-0.3, -0.2, -0.1, 0.0),
    labels = c("-30%", "-20%", "-10%", "0%")
  ) +
  theme(legend.position = "none")
dev.off()


# STATE-SPACE RESIDUALS --------------------------------------------------------

library(MARSS)

# model assumes each age-sex is unique observation with shared states
marss_dat <- mean_dat_plotting %>%
  mutate(year = as.numeric(as.character(year_f)),
         age_sex = paste(sex, age, sep = "_")) %>% 
  arrange(year) %>% 
  select(period, year, age_sex, mean_fl) %>% 
  pivot_wider(names_from = age_sex, values_from = mean_fl)

marss_mat <- marss_dat %>% 
  select(-year, -period) %>% 
  as.matrix() %>% 
  t()
yrs <- marss_dat$year

# estimate parameters
Z.model <- rep(1, length = nrow(marss_mat)) %>% as.factor()
R.model <- "diagonal and equal"
kem1 <- MARSS(marss_mat, model = list(Z = Z.model, R = R.model))
kem2 <- MARSS(marss_mat, model = list(Z = Z.model, R = R.model, U = matrix(0)))
AIC(kem1, kem2)

resids_1 <- MARSSresiduals(kem1, type = "tT")
resids_2 <- MARSSresiduals(kem2, type = "tT")

marss_dat$mod1_resid <- resids_1$state.residuals[1, ]
marss_dat$mod2_resid <- resids_2$state.residuals[1, ]
marss_dat$state2 <- kem2$states[1, ]

# period values 
transition_years <- marss_dat %>% 
  filter(!period == "Gilbert-Clemens") %>% 
  group_by(period) %>% 
  summarize(first_year = min(year)) %>% 
  pull(first_year)

state_timeseries <- ggplot(marss_dat %>% filter(!is.na(state2))) +
  geom_point(aes(x = year, y = state2, shape = period)) +
  scale_shape_manual(values = shape_pal) +
  geom_vline(xintercept = c(transition_years), color = "red") +
  labs(x = "Year", y = "Estimated State") +
  ggsidekick::theme_sleek() +
  geom_text(aes(x = -Inf, y = Inf, label = "a)"), vjust = 1, hjust = -0.75) +
  theme(legend.position = c(0.92, 0.8),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank())

resid_timeseries <- ggplot(marss_dat %>% filter(!is.na(mod2_resid))) +
  geom_line(aes(x = year, y = mod2_resid)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = c(transition_years), color = "red") +
  labs(x = "Year", y = "State Residuals") +
  ggsidekick::theme_sleek() +
  geom_text(aes(x = -Inf, y = Inf, label = "b)"), vjust = 1, hjust = -0.75)

resid_boxplot <- ggplot(marss_dat %>% filter(!is.na(mod2_resid))) +
  geom_boxplot(aes(x = period, y = mod2_resid)) +
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Period", y = "State Residuals") +
  ggsidekick::theme_sleek() +
  geom_text(aes(x = -Inf, y = Inf, label = "c)"), vjust = 1, hjust = -0.75)

png(here::here("outputs", "figs", "state_space.png"), width = 8, height = 9,
    res = 250, units = "in")
cowplot::plot_grid(state_timeseries,
                   resid_timeseries,
                   resid_boxplot, nrow = 3)
dev.off()


# CHECK 1985 INFLUENCE ---------------------------------------------------------

dum_mesh_85 <- make_mesh(dat %>% filter(!year == "1985"), c("x", "y"), 
                         cutoff = 1000)

fit_85 <- sdmTMB(fl_cm ~ s(yday_c, by = age, m = 2) +
                   s(year, by = age, m = 2) +
                   age + sex,
                 dispformula = ~ 0 + period,
                 data = dat %>% filter(!year == "1985"),
                 mesh = dum_mesh_85,
                 spatial = "off",
                 spatiotemporal = "off",
                 control = sdmTMBcontrol(
                   nlminb_loops = 2,
                   newton_loops = 2
                 ))


# smoothed predictions
smooth_preds <- predict(fit_85, newdata = new_dat2, re_form = NA, se_fit = TRUE)
smooth_preds$low <- smooth_preds$est + (qnorm(0.025) * smooth_preds$est_se)
smooth_preds$up <- smooth_preds$est + (qnorm(0.975) * smooth_preds$est_se)


smooth_year_85 <- ggplot(smooth_preds, aes(x = year, y = est)) +
  geom_line() +
  geom_ribbon(aes(ymin = low, ymax = up), alpha = 0.4) +
  # geom_point(data = dat_avg_trim, aes(x = year, y = fl_cm)) +
  ggsidekick::theme_sleek() +
  labs(x = "Year", y = "Fork Length (cm)") +
  facet_wrap(~age_f, labeller = label_parsed) +
  scale_x_continuous(
    breaks = seq(1915, 2015, by = 20),
    expand = c(0.02, 0.02)
  ) 

png(here::here("outputs", "figs", "smooth_85_removed.png"), 
    height = 4, width = 5, units = "in", res = 250)
smooth_year_85
dev.off()


# FIT SUPP MODEL 1 -------------------------------------------------------------

# fit alternative model with intercepts for sampling period to include in supp
# analysis (improved AIC, but parameter estimates for intercepts are unreliable)

fit_p <- sdmTMB(fl ~ s(yday_c, by = age, m = 2) +
                s(year, by = age, m = 2) +
                period +
                age + sex,
              dispformula = ~ period,
              data = dat,
              spatial = "off",
              spatiotemporal = "off",
              control = sdmTMBcontrol(
                nlminb_loops = 2,
                newton_loops = 2
              )
              )
sanity(fit_p)


# generate annual predictions
new_dat_p <- expand.grid(
  age = unique(dat$age),
  sex = "female",
  yday_c = 0,
  year = seq(min(dat$year), max(dat$year), length = 100)
) %>% 
  mutate(
    age_f = as.factor(age),
    period = "Gilbert-Clemens"
  )

# insufficient vector memory to generate confidence intervals
smooth_preds_p <- predict(fit_p, 
                          newdata = new_dat_p)

smooth_year_p <- ggplot(smooth_preds_p, aes(x = year, y = est)) +
  geom_line() +
  ggsidekick::theme_sleek() +
  labs(x = "Year", y = "Fork Length (mm)") +
  facet_wrap(~age_f, labeller = label_parsed) +
  scale_x_continuous(
    breaks = seq(1915, 2015, by = 20),
    expand = c(0.02, 0.02)
  ) 


png(here::here("outputs", "figs", "year_smooth_period_ints.png"), 
    height = 5, width = 8.5,
    units = "in", res = 250)
smooth_year_p
dev.off()


## calculate difference in size (as above)
smooth_list_p <- split(smooth_preds_p, smooth_preds_p$age_f) 
age_names <- c("42", "52", "53", "63")

purrr::map2(smooth_list_p, age_names, function (x, y) {
  ts <- x$est
  #calc time series mean for each iteration
  ts_mean = mean(ts)
  data.frame(
    mean_diff = ts[length(ts)] - ts_mean,
    first_diff = ts[length(ts)] - ts[1],
    mean_rel_diff = (ts[length(ts)] - ts_mean) / ts_mean,
    first_rel_diff = (ts[length(ts)] - ts[1]) / ts[1],
    age = y
  )
}) %>% 
  bind_rows()


## simulate and attempt to recover pars

sims_p <- simulate(fit_p, 50)
fix_pars_est <- tidy(fit_p, effects = "fixed") %>% 
  mutate(model = "period")
fix_pars_est$iter <- NA_character_
fix_pars_est$up <- fix_pars_est$estimate + (1.96 * fix_pars_est$std.error)
fix_pars_est$lo <- fix_pars_est$estimate - (1.96 * fix_pars_est$std.error)

out_list <- vector(mode = "list", length = 50)
for (i in seq_len(ncol(sims_p))) {
  dumm <- dat %>% mutate(sim_fl = sims_p[ , i])
  fit_dumm <- sdmTMB(sim_fl ~ s(yday_c, by = age, m = 2) +
                    s(year, by = age, m = 2) +
                    period +
                    age + sex,
                  dispformula = ~ 0 + period,
                  data = dumm,
                  spatial = "off",
                  spatiotemporal = "off",
                  control = sdmTMBcontrol(
                    nlminb_loops = 2,
                    newton_loops = 2
                  )
  )
  out_list[[i]] <- tidy(fit_dumm, effects = "fixed") %>% 
    mutate(model = paste("sim", i, sep = "_"))
  
  rm("fit_dumm")
  gc()
  
  saveRDS(out_list, 
          here::here("outputs", "model_fits", "sim_fix_eff.rds"))
}

sim_dat <- out_list %>% 
  bind_rows() %>% 
  mutate(
    iter = model,
    model = "sim"
  )
sim_dat$model <- "period"


## as above but with original model
sims <- simulate(fit, 50)
fix_pars_est1 <- tidy(fit, effects = "fixed") %>% 
  mutate(model = "no_period",
         iter = NA_character_,
         up = estimate + (1.96 * std.error),
         lo = estimate - (1.96 * std.error))

out_list2 <- vector(mode = "list", length = 50)
for (i in seq_len(ncol(sims))) {
  dumm <- dat %>% mutate(sim_fl = sims[ , i])
  fit_dumm <- sdmTMB(sim_fl ~ s(yday_c, by = age, m = 2) +
                       s(year, by = age, m = 2) +
                       # period +
                       age + sex,
                     dispformula = ~ 0 + period,
                     data = dumm,
                     spatial = "off",
                     spatiotemporal = "off",
                     control = sdmTMBcontrol(
                       nlminb_loops = 2,
                       newton_loops = 2
                     )
  )
  out_list2[[i]] <- tidy(fit_dumm, effects = "fixed") %>% 
    mutate(model = paste("sim_no_period", i, sep = "_"))
  
  rm("fit_dumm")
  gc()
  
  saveRDS(out_list2, 
          here::here("outputs", "model_fits", "sim_fix_eff_no_period.rds"))
}

sim_dat2 <- out_list2 %>% 
  bind_rows() %>% 
  mutate(
    iter = model,
    model = "no_period"
  )

# combine
sim_dat_all <- rbind(sim_dat, sim_dat2)
fix_pars_all <- rbind(fix_pars_est,
                      fix_pars_est1)
saveRDS(list(sim_dat = sim_dat_all,
             fix_pars = fix_pars_all),
        here::here("outputs", "model_fits", "sim_recovery_supp_dat.rds"))

png(here::here("outputs", "figs", "par_recovery_period_ints.png"),
    height = 5, width = 5, units = "in", res = 200)
ggplot() +
  geom_boxplot(data = sim_dat_all,
               aes(x = term, y = estimate, fill = model), alpha = 0.4) +
  geom_pointrange(data = fix_pars_all,
                  aes(x = term, y = estimate, ymin = lo, ymax = up, 
                      fill = model, colour = model),
                  position = position_jitterdodge(),
                  # position = position_dodge(width = 0.85),
                  shape = 21) +
  facet_wrap(~term, scales = "free") +
  ggsidekick::theme_sleek()
dev.off()



# FIT SUPP MODEL 2 -------------------------------------------------------------

# fit alternative model with unique age-sex effects
dat$age_sex <- as.factor(paste(dat$age, dat$sex, sep = "_"))

fit_as <- sdmTMB(fl ~ s(yday_c, by = age, m = 2) +
                  s(year, by = age_sex, m = 2) +
                  period +
                  age_sex,
                dispformula = ~ 0 + period,
                data = dat,
                spatial = "off",
                spatiotemporal = "off",
                control = sdmTMBcontrol(
                  nlminb_loops = 2,
                  newton_loops = 2
                )
)
sanity(fit_as)


# generate annual predictions
new_dat_as <- expand.grid(
  age = unique(dat$age),
  sex = unique(dat$sex),
  yday_c = 0,
  year = seq(min(dat$year), max(dat$year), length = 100)
) %>% 
  mutate(
    age_f = as.factor(age),
    period = "Gilbert-Clemens",
    age_sex = as.factor(paste(age, sex, sep = "_"))
  ) 

smooth_preds_as <- predict(fit_as, 
                          newdata = new_dat_as)

png(here::here("outputs", "figs", "year_smooth_age_sex_int.png"),
    height = 5, width = 5, units = "in", res = 200)
ggplot(smooth_preds_as, aes(x = year, y = est, colour = sex)) +
  geom_line() +
  ggsidekick::theme_sleek() +
  labs(x = "Year", y = "Fork Length (mm)") +
  facet_wrap(~age_f, labeller = label_parsed) +
  scale_x_continuous(
    breaks = seq(1915, 2015, by = 20),
    expand = c(0.02, 0.02)
  ) 
dev.off()