# Fit brms model to century of Nass data 
# Extension of gam_fit 
# 1) Import Skip's first century of data and clean as necessary
# 2) Estimate temporal trend using GAM and brms to account for changes in
# variability with sampling regime
# Created by C Freshwater Dec 23, 2021
# -------------------------------------------------

library(tidyverse)
library(brms)
library(tidybayes)
library(posterior)

# import nass data
dat_in <- read.table(here::here("data", "nasscentury.txt"))
colnames(dat_in) <- c("year", "month", "day", "hart_species", "sex", "age", 
                      "fl")
dat <- dat_in %>% 
  mutate(
    year_f = as.factor(year),
    sex = ifelse(sex == 1, "male", "female"),
    age = factor(as.character(age)),
    date = as.POSIXct(paste(day, month, year, sep = "-"),
                      format = "%d-%m-%Y"),
    yday = lubridate::yday(date),
    yday_c = yday - mean(yday, na.rm = TRUE)
  ) %>% 
  #drop secondary age classes and years missing capture data
  filter(
    age %in% c("42", "52", "53", "63"),
    !is.na(yday)
  )  %>%
  # Add a factor representing sample collection to visualize effects 
  # NOTE: period is not modeled, we're simply assigning a color in the plot to 
  # each data set
  mutate(
    period = case_when(
      year < 1947 ~ "Gilbert-Clemens",
      year > 1993 ~ "Nisga'a",
      TRUE ~ "Bilton"
    ),
    period = fct_relevel(as.factor(period), 
                         "Gilbert-Clemens",
                         "Bilton",
                         "Nisga'a"),
    fl_cm = fl / 10,
    # add random identifier to split datasets for subsetting,
    data_group = sample.int(5, nrow(.), replace = T)
  )


# RAW DATA PLOTS ---------------------------------------------------------------

# add dummy year levels for missing data
dat2 <- expand.grid(
  year_f = seq(min(dat$year), max(dat$year), by = 1) %>% as.factor(),
  age = unique(dat$age),
  sex = unique(dat$sex)
) %>% 
  left_join(., dat, by = c("year_f", "age", "sex")) 

ggplot(dat2) +
  geom_point(aes(x = year_f, y = yday_c, fill = period)) +
  # facet_wrap(~age) +
  ggsidekick::theme_sleek() +
  labs(x = "Year", y = "Fork Length") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# box plots
annual_box <- ggplot(dat2 %>% 
                       filter(sex == "female")) +
  geom_boxplot(aes(x = year_f, y = fl, fill = period)) +
  facet_wrap(~age) +
  ggsidekick::theme_sleek() +
  labs(x = "Year", y = "Fork Length") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


mean_dat <- dat2 %>% 
  group_by(year_f, age, period, sex) %>% 
  summarize(
    n = n(),
    mean_yday = mean(yday, na.rm = T),
    sd_yday = sd(yday_c, na.rm = T),
    mean_fl = mean(fl, na.rm = T),
    sd_fl = sd(fl, na.rm = T),
    se_fl = sd_fl / sqrt(n),
    up = mean_fl + sd_fl,#(1.96 * se_fl),
    lo = mean_fl - sd_fl,#(1.96 * se_fl),
    .groups = "drop"
  ) %>% 
  distinct() %>% 
  ungroup() 

annual_dot <- ggplot(mean_dat %>% 
                       filter(sex == "female"), 
                     aes(x = year_f)) +
  geom_pointrange(aes(y = mean_fl, fill = period, ymin = lo, ymax = up),
                  shape = 21) +
  facet_wrap(~age) +
  ggsidekick::theme_sleek() +
  labs(x = "Year", y = "Fork Length") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

annual_sd_dot <- ggplot(mean_dat, aes(x = year_f)) +
  geom_point(aes(y = sd_fl, fill = period), shape = 21) +
  facet_grid(~age) +
  ggsidekick::theme_sleek() +
  labs(x = "Year", y = "SD Fork Length") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(mean_dat, aes(x = year_f)) +
  geom_point(aes(y = mean_yday, fill = period), shape = 21) +
  facet_wrap(~age) +
  ggsidekick::theme_sleek() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# export
png(here::here("outputs", "figs", "annual_dot.png"), width = 11, height = 7,
    res = 250, units = "in")
annual_dot
dev.off()

png(here::here("outputs", "figs", "annual_box.png"), width = 11, height = 7,
    res = 250, units = "in")
annual_box
dev.off()


# FIT MODEL  -------------------------------------------------------------------

#check priors 
brms::get_prior(fl ~ s(yday_c) + s(year, m = 2, k = 5) + 
                  s(year, by = age, m = 1, k = 5) +
                  age + sex + period,
                data = dat)

library(future)
plan(multiprocess)

split_dat <- split(dat, dat$data_group)


# trim_dat <- dat %>% 
#   sample_n(size = 20000)
# 
# 
# ## experiment with multiple; ~4x faster with 20000 samples than running one
# 
# tictoc::tic()
# brm_test <- brm_multiple(
#   bf(
#     fl_cm ~ #s(yday_c, k = 3) + 
#       s(year, m = 2, k = 5) +
#       s(year, by = age, m = 1, k = 5) +
#       age + sex + period,
#     sigma ~ period
#   ),
#   data = split_dat, seed = 17,
#   iter = 1000, thin = 10, chains = 2, refresh = 0,# warmup = 750, 
#   control = list(adapt_delta = 0.94, max_treedepth = 12)
# )
# tictoc::toc()
# 
# post <- as.array(brm_test)
# bayesplot::mcmc_trace(brm_test)
# round(brm_test$rhats, 3)
# conditional_effects(brm_test, "year")


tictoc::tic()
brm1 <- brm_multiple(
  bf(
    fl_cm ~ s(yday_c, k = 3) + s(year, m = 2, k = 5) + 
      s(year, by = age, m = 1, k = 5) +
      age + sex + period,
    sigma ~ period
  ),
  data = split_dat, seed = 17,
  iter = 2000, thin = 10, chains = 3, warmup = 750, #refresh = 0,  
  control = list(adapt_delta = 0.95, max_treedepth = 12
                 ),
  prior=c(prior(normal(60, 15), class="Intercept"),
          prior()
          prior(normal(0, 10), class="b", coef = "age52"),
          prior(normal(0, 10), class="b", coef = "age53"),
          prior(normal(0, 10), class="b", coef = "age63"),
          prior(normal(0, 10), class="b", coef = "sexmale"),
          prior(normal(0, 10), class="b", coef = "periodGC"),
          prior(normal(0, 10), class="b", coef = "periodmod"),
          prior(normal(0, 50), class="b", coef = "syear_1"),
          prior(normal(0, 10), class="b", coef = "syday_c_1"),
          prior(normal(0, 2), class="b", dpar = "sigma")#,
          # prior(normal(0, 2), class="b", coef = "periodmod", dpar = "sigma"),
          # prior(normal(0, 2), class="b", coef = "periodGC", dpar = "sigma")
  )
)
tictoc::toc()

saveRDS(brm1, here::here("outputs", "data", "brms_fits", "ind_ls.rds"))
brm1 <- readRDS(here::here("outputs", "data", "brms_fits", "ind_ls.rds"))

post <- as.array(brm1)
bayesplot::mcmc_trace(post)
round(brm1$rhats, 3)
conditional_effects(brm1, "year")

dat %>%
  add_residual_draws(brm1) %>%
  median_qi() %>%
  ggplot(aes(sample = .residual)) +
  geom_qq() +
  geom_qq_line()

# posterior predictive checks (student-t similar)
pp_check(brm1)
pp_check(brm1, type='stat', stat='mean')
pp_check(brm1, type='error_scatter_avg')
pp_check(brm1, type='error_hist')
pp_check(brm1, type='intervals')
pp_check(brm1, type='scatter_avg')
pp_check(brm1, x = 'yday_c', type='error_scatter_avg_vs_x')

mcmc_plot(brm_test, type = "neff_hist")
mcmc_plot(brm_test, type = "neff")
mcmc_plot(brm_test, type = "hist")

msms <- conditional_smooths(brm1)
tt <- plot(msms)

pdf(here::here("outputs", "figs", "brm_fit_summary.pdf"))
purrr::map(tt, print)
dev.off()


## POSTERIOR ESTIMATES ---------------------------------------------------------

draws <- as_draws_array(brm1)
sum_draws <- summarise_draws(draws, default_summary_measures())

## coefficient estimates for fixed effects for period of mu 
mu_post_bilton <- draws[ , ,"b_Intercept"]
mu_post_gc <- draws[ , , "b_Intercept"] + draws[ , , "b_periodGC"]
mu_post_mod <- draws[ , , "b_Intercept"] + draws[ , , "b_periodmod"] 
mu_list <- list(mu_post_gc,
                   mu_post_bilton,
                   mu_post_mod)
mu_df <- data.frame(
  period = c("Gilbert-Clemens", "Bilton", "Nisga'a"),
  median = purrr::map(mu_list, median) %>% unlist(),
  low = purrr::map(mu_list, quantile, 0.05) %>% unlist(),
  high = purrr::map(mu_list, quantile, 0.95) %>% unlist(),
  parameter = "mu"
)


## as above for sigma
sigma_post_bilton <- draws[ , ,"b_sigma_Intercept"]
sigma_post_gc <- draws[ , , "b_sigma_Intercept"] + 
  draws[ , , "b_sigma_periodGC"]
sigma_post_mod <- draws[ , , "b_sigma_Intercept"] + 
  draws[ , , "b_sigma_periodmod"] 
sigma_list <- list(sigma_post_gc,
                   sigma_post_bilton,
                   sigma_post_mod)
sigma_df <- data.frame(
  period = c("Gilbert-Clemens", "Bilton", "Nisga'a"),
  median = purrr::map(sigma_list, median) %>% unlist(),
  low = purrr::map(sigma_list, quantile, 0.05) %>% unlist(),
  high = purrr::map(sigma_list, quantile, 0.95) %>% unlist(),
  parameter = "sigma"
) 

period_est <- rbind(mu_df, sigma_df) %>% 
  mutate(period = fct_relevel(factor(period),
                              "Gilbert-Clemens", "Bilton", "Nisga'a")) %>% 
  ggplot(.) +
  geom_pointrange(aes(x = period, y = median, fill = period, 
                      ymin = low, ymax = high),
                  shape = 21) +
  facet_wrap(~parameter, scales = "free_y") +
  ggsidekick::theme_sleek() +
  labs(x = "Sampling Regime", y = "Posterior Estimate") 

png(here::here("outputs", "figs", "period_est.png"), height = 4, width = 7.5,
    units = "in", res = 250)
period_est
dev.off()


# first predictive dataset only has years observed in dataset; these are paired
# with appropriate periods to estimate overall effects
pred_dat <- expand.grid(
  year = seq(min(dat$year), max(dat$year), by = 0.1),
  age = unique(dat$age),
  yday_c = 0,
  sex = unique(dat$sex),
  period = "GC"
) %>% 
 droplevels() 

post_pred <- posterior_predict(brm1, pred_dat, allow_new_levels = TRUE)
pred_dat_bayes <- pred_dat %>% 
  mutate(
    median = apply(post_pred, 2, median),
    low = apply(post_pred, 2, function (x) quantile(x, probs = 0.05)),
    up = apply(post_pred, 2, function (x) quantile(x, probs = 0.95))
  ) 

smooth_pred1 <- ggplot(pred_dat_bayes, aes(x = year, y = median)) +
  geom_line() +
  geom_ribbon(aes(ymin = low, ymax = up), alpha = 0.4) +
  facet_grid(sex~age) + 
  ggsidekick::theme_sleek() +
  labs(x = "Year", y = "Fork Length (cm)")

png(here::here("outputs", "figs", "smooth_1period.png"), 
    height = 4, width = 7.5, units = "in", res = 250)
smooth_pred1
dev.off()



# second predictive dataset associates years with periods (approximate)
pred_dat2 <- pred_dat %>% 
  mutate(
    period = case_when(
      year < 1947 ~ "GC",
      year > 1993 ~ "mod",
      TRUE ~ "bilton"
    )
  )
mean_dat$year <- as.numeric(as.character(mean_dat$year_f))
mean_dat$mean_fl_cm <- mean_dat$mean_fl / 10


post_pred2 <- posterior_predict(brm1, pred_dat2, allow_new_levels = TRUE)
post_ribbon2 <- pred_dat2 %>% 
  mutate(
    mean = apply(post_pred2, 2, mean),
    low = apply(post_pred2, 2, function (x) quantile(x, probs = 0.05)),
    up = apply(post_pred2, 2, function (x) quantile(x, probs = 0.95))
  ) %>% 
  ggplot(., aes(x = year, y = mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = low, ymax = up), alpha = 0.4) +
  facet_grid(sex~age) +
  ggsidekick::theme_sleek() 

smooth_pred2 <- post_ribbon2 +
  geom_point(data = mean_dat %>% filter(!is.na(period)), 
             aes(x = year, y = mean_fl_cm, fill = period), 
             shape = 21, alpha = 0.4) +
  labs(x = "Year", y = "Fork Length (cm)")

png(here::here("outputs", "figs", "smooth_3periods.png"), 
    height = 4, width = 7.5,
    units = "in", res = 250)
smooth_pred2
dev.off()


# visualize year-day effect
pred_dat3 <- expand.grid(
  year = mean(dat$year),
  yday_c = seq(-46, 64, by = 0.1),
  age = "42",
  sex = "male",
  period = "GC"
) 

post_pred3 <- posterior_predict(brm1, pred_dat3, allow_new_levels = TRUE)
smooth_yday <- pred_dat3 %>% 
  mutate(
    median = apply(post_pred3, 2, median),
    low = apply(post_pred3, 2, function (x) quantile(x, probs = 0.05)),
    up = apply(post_pred3, 2, function (x) quantile(x, probs = 0.95))
  ) %>% 
  ggplot(., aes(x = yday_c, y = median)) +
  geom_line() +
  geom_ribbon(aes(ymin = low, ymax = up), alpha = 0.4) +
  ggsidekick::theme_sleek() +
  labs(x = "Centered Year Day", y = "Fork Length (cm)")


png(here::here("outputs", "figs", "smooth_yday.png"), 
    height = 4.5, width = 4.5, units = "in", res = 250)
smooth_yday
dev.off()


# predict differences between final year and time series average for each age
pred_dat4 <- expand.grid(
  year = seq(min(dat$year), max(dat$year), by = 1),
  age = unique(dat$age),
  yday_c = 0,
  sex = "male"
  ) %>% 
  mutate(
    period = case_when(
      year < 1947 ~ "GC",
      year > 1993 ~ "mod",
      TRUE ~ "bilton"
    )
  ) %>% 
  droplevels()
post_pred4 <- posterior_predict(brm1, pred_dat4, allow_new_levels = TRUE)

# time series mean by age
age53 <- post_pred4[ , which(pred_dat4$age == "53")]
age52 <- post_pred4[ , which(pred_dat4$age == "52")]
age42 <- post_pred4[ , which(pred_dat4$age == "42")]
age63 <- post_pred4[ , which(pred_dat4$age == "63")]
age_list <- list(age53, age52, age42, age63)
age_names <- c("53", "52", "42", "63")

purrr::map2(age_list, age_names, function (x, y) {
  #calc time series mean for each iteration
  mu <- apply(x, 1, mean)
  diff <- x[ , ncol(x)] - mu
  diff2 <- x[ , ncol(x)] - x[ , 1]
  data.frame(
    mean_diff1 = mean(diff),
    mean_diff2 = mean(diff2),
    low = quantile(diff2, probs = 0.05),
    high = quantile(diff2, probs = 0.95),
    age = y
  )
}) %>% 
  bind_rows


### MGCV COMPARE ---------------------------------------------------------------

library(mgcv)
library(gratia)

gam1 <- gam(
  fl_cm ~ s(yday_c, k = 3) + s(year, m = 2, k = 5) + 
    s(year, by = age, m = 1, k = 5) +
    age + sex + period, 
  data = dat,
  method = "REML"
)

mgcv_splines <- draw(gam1)

pdf(here::here("outputs", "figs", "gam_fit_summary.pdf"))
mgcv_splines
dev.off()