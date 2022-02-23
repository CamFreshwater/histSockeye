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
    sex = ifelse(sex == 1, "male", "female") %>% as.factor(),
    age = factor(as.character(age)),
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
  # NOTE: period is not modeled, we're simply assigning a color in the plot to 
  # each data set
  mutate(
    period = case_when(
      year < 1947 ~ "Gilbert-Clemens",
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
    data_group = sample.int(5, nrow(.), replace = T)
  ) %>% 
  droplevels()


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
    up = mean_fl + (1.96 * se_fl),
    lo = mean_fl - (1.96 * se_fl),
    .groups = "drop"
  ) %>% 
  distinct() %>% 
  group_by(age, sex) %>% 
  mutate(
    overall_mean = mean(mean_fl, na.rm = T)) %>% 
  ungroup() %>% 
  droplevels() 
levels(mean_dat$age) <- c("4[2]", "5[2]", "5[3]", "6[3]")


annual_dot <- ggplot(
  mean_dat %>% 
    filter(sex == "female",
           !is.na(period)), 
  aes(x = as.numeric(as.character(year_f)))
) +
  geom_pointrange(aes(y = mean_fl, fill = period, ymin = lo, ymax = up),
                  shape = 21) +
  facet_wrap(~age, labeller = label_parsed) +
  ggsidekick::theme_sleek() +
  labs(x = "Year", y = "Fork Length") +
  scale_x_continuous(
    breaks = seq(1915, 2015, by = 20),
    expand = c(0.02, 0.02)
  ) +
  geom_hline(aes(yintercept = overall_mean), lty = 2) +
  scale_fill_discrete(name = "Sampling\nPeriod")

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

# parallelize based on operating system (should speed up some of the spatial
# processing calculations)
library(parallel)
ncores <- detectCores() - 2
if (Sys.info()['sysname'] == "Windows") {
  library(doParallel)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
} else {
  doMC::registerDoMC(ncores)
}


#check priors 
brms::get_prior(bf(
  fl_cm ~ s(yday_c, m = 2, k = 3) + s(yday_c, by = age, m = 1, k = 3) +
    s(year, m = 2, k = 5) + s(year, by = age, m = 1, k = 5) +
    age + sex + period,
  sigma ~ period
),
data = dat)

# library(future)
# plan(multiprocess)

split_dat <- split(dat, dat$data_group)


# trim_dat <- dat %>%
#   sample_n(size = 20000) %>% 
#   split(., .$data_group)


tictoc::tic()
brm1 <- brm_multiple(
  bf(
    fl_cm ~ s(yday_c, m = 2, k = 3) + s(yday_c, by = age, m = 1, k = 3) +
      s(year, m = 2, k = 5) + s(year, by = age, m = 1, k = 5) +
      age + sex + period,
    sigma ~ period
  ),
  data = split_dat, seed = 17,
  iter = 2250, thin = 10, chains = 3, warmup = 750, #refresh = 0,  
  control = list(adapt_delta = 0.97, max_treedepth = 14
                 ),
  prior=c(prior(normal(60, 15), class="Intercept"),
          prior(normal(5, 10), class="b", coef = "age52"),
          prior(normal(5, 10), class="b", coef = "age53"),
          prior(normal(5, 10), class="b", coef = "age63"),
          prior(normal(5, 10), class="b", coef = "sexmale"),
          prior(normal(0, 10), class="b", coef = "periodBilton"),
          prior(normal(0, 10), class="b", coef = "periodMonkleyDump"),
          prior(normal(0, 10), class="b", coef = "periodNisgaa"),
          prior(normal(0, 30), class="b", coef = "syear_1"),
          prior(normal(0, 30), class="b", coef = "syday_c_1"),
          prior(normal(0, 10), class="b", dpar = "sigma"),
          prior(exponential(0.5), class = "Intercept", dpar = "sigma")
  )
)
tictoc::toc()

# saveRDS(brm1, here::here("outputs", "data", "brms_fits", "ind_ls.rds"))
brm1 <- readRDS(here::here("outputs", "data", "brms_fits", "ind_ls.rds"))

post <- as.array(brm1)
bayesplot::mcmc_trace(post)
round(brm1$rhats, 3)
conditional_effects(brm1, "year")

dat %>%
  droplevels() %>% 
  add_residual_draws(brm1) %>%
  median_qi() %>%
  ggplot(aes(sample = .residual)) +
  geom_qq() +
  geom_qq_line()

# posterior predictive checks 
pp_check(brm1)
pp_check(brm1, type='stat', stat='mean')
pp_check(brm1, type='error_scatter_avg')
pp_check(brm1, type='error_hist')
pp_check(brm1, type='intervals')
pp_check(brm1, type='scatter_avg')
pp_check(brm1, x = 'yday_c', type='error_scatter_avg_vs_x')

mcmc_plot(brm1, type = "neff_hist")
mcmc_plot(brm1, type = "neff")
mcmc_plot(brm1, type = "hist")

msms <- conditional_smooths(brm1)
tt <- plot(msms)

pdf(here::here("outputs", "figs", "brm_fit_summary.pdf"))
purrr::map(tt, print)
dev.off()


## POSTERIOR ESTIMATES ---------------------------------------------------------

draws <- as_draws_array(brm1)
sum_draws <- summarise_draws(draws, default_summary_measures())

## coefficient estimates for fixed effects for period  
mu_post_bilton <- draws[ , ,"b_Intercept"] + draws[ , , "b_periodBilton"]
mu_post_gc <- draws[ , , "b_Intercept"] 
mu_post_mod <- draws[ , , "b_Intercept"] + draws[ , , "b_periodNisgaa"]
mu_post_dfo <- draws[ , , "b_Intercept"] + draws[ , , "b_periodMonkleyDump"] 
mu_list <- list(mu_post_gc,
                mu_post_bilton,
                mu_post_dfo,
                mu_post_mod)
mu_df <- data.frame(
  period = c("Gilbert-Clemens", "Bilton", "Monkley Dump", "Nisga'a"),
  median = purrr::map(mu_list, median) %>% unlist(),
  low = purrr::map(mu_list, quantile, 0.05) %>% unlist(),
  high = purrr::map(mu_list, quantile, 0.95) %>% unlist()
) 
mu_df$period = fct_relevel(as.factor(mu_df$period),
                           "Gilbert-Clemens", "Bilton", "Monkley Dump",
                           "Nisga'a")

area_eff_dot <- ggplot(mu_df) +
  geom_pointrange(aes(x = period, y = median, fill = period, 
                      ymin = low, ymax = high),
                  shape = 21) +
  ggsidekick::theme_sleek() +
  labs(x = "Sampling Period", y = "Posterior Estimate Mean") +
  theme(legend.position = "none")


## coefficient estimates for fixed effects for age (assuming female and bilton)
mu_post_42 <- mu_post_gc
mu_post_52 <- draws[ , , "b_Intercept"] + draws[ , , "b_age52"]
mu_post_53 <- draws[ , , "b_Intercept"] + draws[ , , "b_age53"]
mu_post_63 <- draws[ , , "b_Intercept"] + draws[ , , "b_age63"]
mu_age_list <- list(mu_post_42, mu_post_52, mu_post_53, mu_post_63)

age_eff <- data.frame(
  age = c("4[2]", "5[2]", "5[3]", "6[3]"),
  median = purrr::map(mu_age_list, median) %>% unlist(),
  low = purrr::map(mu_age_list, quantile, 0.05) %>% unlist(),
  high = purrr::map(mu_age_list, quantile, 0.95) %>% unlist()
) 

age_eff_dot <- ggplot(age_eff) +
  geom_pointrange(aes(x = age, y = median, ymin = low, ymax = high, fill = age),
                  shape = 21
                  ) +
  ggsidekick::theme_sleek() +
  scale_x_discrete("Age Class", labels = parse(text = unique(age_eff$age))) +
  labs(x = "Age Class", y = "Posterior Estimate Mean") +
  scale_fill_brewer(palette = 1) +
  theme(legend.position = "none") 

  

## coefficient estimates for fixed effects for sex (assuming 42 and bilton)
mu_post_fem <- mu_post_gc
mu_post_male <- draws[ , , "b_Intercept"] + draws[ , , "b_sexmale"]
mu_sex_list <- list(mu_post_fem, mu_post_male)

sex_eff_dot <- data.frame(
  sex = c("female", "male"),
  median = purrr::map(mu_sex_list, median) %>% unlist(),
  low = purrr::map(mu_sex_list, quantile, 0.05) %>% unlist(),
  high = purrr::map(mu_sex_list, quantile, 0.95) %>% unlist()
) %>% 
  ggplot(.) +
  geom_pointrange(aes(x = sex, y = median, fill = sex, 
                      ymin = low, ymax = high),
                  shape = 21) +
  ggsidekick::theme_sleek() +
  labs(x = "Sex", y = "Posterior Estimate Mean") +
  scale_fill_brewer(palette = 5, type = "seq") +
  theme(legend.position = "none")


## as above for sigma
sigma_post_bilton <- draws[ , ,"b_sigma_Intercept"] + 
  draws[ , , "b_sigma_periodBilton"]
sigma_post_gc <- draws[ , , "b_sigma_Intercept"] 
sigma_post_mod <- draws[ , , "b_sigma_Intercept"] + 
  draws[ , , "b_sigma_periodNisgaa"]
sigma_post_dfo <- draws[ , , "b_sigma_Intercept"] + 
  draws[ , , "b_sigma_periodMonkleyDump"] 
sigma_list <- list(sigma_post_gc,
                   sigma_post_bilton,
                   sigma_post_dfo,
                   sigma_post_mod)


sigma_df <- data.frame(
  period = c("Gilbert-Clemens", "Bilton", "Monkley Dump", "Nisga'a"),
  median = purrr::map(sigma_list, median) %>% unlist(),
  low = purrr::map(sigma_list, quantile, 0.05) %>% unlist(),
  high = purrr::map(sigma_list, quantile, 0.95) %>% unlist()
) 
sigma_df$period = fct_relevel(as.factor(sigma_df$period),
                           "Gilbert-Clemens", "Bilton", "Monkley Dump",
                           "Nisga'a")

area_sigma_dot <- ggplot(sigma_df) +
  geom_pointrange(aes(x = period, y = median, fill = period, 
                      ymin = low, ymax = high),
                  shape = 21) +
  ggsidekick::theme_sleek() +
  # geom_text(aes(x = -Inf, y = Inf, label = "d)")) +
  labs(x = "Sampling Period", y = "Posterior Estimate SD") +
  theme(legend.position = "none")


png(here::here("outputs", "figs", "main_effect.png"), 
    height = 5, width = 8.5,
    units = "in", res = 250)
cowplot::plot_grid(area_eff_dot,
                   age_eff_dot,
                   sex_eff_dot,
                   area_sigma_dot, ncol = 2)
dev.off()


# first predictive dataset only has years observed in dataset; these are paired
# with appropriate periods to estimate overall effects
pred_dat <- expand.grid(
  year = seq(min(dat$year), max(dat$year), by = 0.1),
  age = unique(dat$age),
  yday_c = 0,
  sex = unique(dat$sex),
  period = "Gilbert-Clemens"
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


# as above but with fitted draws, not posterior predictions
fit_year <- pred_dat %>% 
  add_fitted_draws(model = brm1, re_formula = NULL,
                   allow_new_levels = TRUE) %>%
  ggplot(., aes(x = year)) +
  stat_lineribbon(aes(y = .value), 
                  .width = c(.9, .5), alpha = 0.25) +
  ggsidekick::theme_sleek() +
  labs(x = "Year", y = "Fork Length (cm)") +
  facet_grid(sex~age) 

png(here::here("outputs", "figs", "smooth_1period_fit.png"), 
    height = 4, width = 7.5, units = "in", res = 250)
fit_year
dev.off()



# second predictive dataset associates years with periods (approximate)
pred_dat2 <- pred_dat %>% 
  mutate(
    period = case_when(
      year < 1947 ~ "Gilbert-Clemens",
      year > 1972 & year < 1994 ~ "Monkley Dump",
      year > 1993 ~ "Nisga'a",
      TRUE ~ "Bilton"
    ),
    period = fct_relevel(as.factor(period), 
                         "Gilbert-Clemens",
                         "Bilton",
                         "Monkley Dump",
                         "Nisga'a")
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

png(here::here("outputs", "figs", "smooth_4periods.png"), 
    height = 4, width = 7.5,
    units = "in", res = 250)
smooth_pred2
dev.off()


# visualize year-day effect
pred_dat3 <- expand.grid(
  year = mean(dat$year),
  yday_c = seq(-46, 64, by = 0.1),
  age = unique(dat$age),
  sex = "male",
  period = "Gilbert-Clemens"
) 

pp_yday <- pred_dat3 %>% 
  add_predicted_draws(model = brm1, re_formula = NULL, 
                      allow_new_levels = TRUE) %>% 
  ggplot(., aes(x = yday_c)) +
  stat_lineribbon(aes(y = .prediction), 
                  .width = c(.9, .5), alpha = 0.25) +
  ggsidekick::theme_sleek() +
  labs(x = "Centered Year Day", y = "Fork Length (cm)") +
  facet_wrap(~age)

fit_yday <- pred_dat3 %>% 
  add_fitted_draws(model = brm1, re_formula = NULL,
                        allow_new_levels = TRUE) %>%
  ggplot(., aes(x = yday_c)) +
  stat_lineribbon(aes(y = .value), 
                  .width = c(.9, .5), alpha = 0.25) +
  ggsidekick::theme_sleek() +
  labs(x = "Centered Year Day", y = "Fork Length (cm)") +
  facet_wrap(~age)

png(here::here("outputs", "figs", "smooth_yday.png"), 
    height = 4.5, width = 4.5, units = "in", res = 250)
pp_yday
dev.off()

png(here::here("outputs", "figs", "fit_smooth_yday.png"), 
    height = 4.5, width = 4.5, units = "in", res = 250)
fit_yday
dev.off()


# predict differences between final year and time series average for each age
pred_dat4 <- expand.grid(
  year = seq(min(dat$year), max(dat$year), by = 1),
  age = unique(dat$age),
  yday_c = 0,
  sex = "male"
  ) %>% 
  mutate(
    period = "Gilbert-Clemens"
    # period = case_when(
    #   year < 1947 ~ "Gilbert-Clemens",
    #   year > 1972 & year < 1994 ~ "Monkley Dump",
    #   year > 1993 ~ "Nisga'a",
    #   TRUE ~ "Bilton"
    # ),
    # period = fct_relevel(as.factor(period), 
    #                      "Gilbert-Clemens",
    #                      "Bilton",
    #                      "Monkley Dump",
    #                      "Nisga'a")
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
  #difference from mean
  diff <- x[ , ncol(x)] - mu
  #difference from first obs
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
  fl_cm ~ s(yday_c, k = 3, by = age) + s(year, m = 2, k = 5) + 
    s(year, by = age, m = 1, k = 5) +
    age + sex + period, 
  data = dat,
  method = "REML"
)
gam2 <- gam(
  fl_cm ~ s(yday_c, m = 2, k = 3) + s(yday_c, by = age, m = 1, k = 3) +
    s(year, m = 2, k = 5) + s(year, by = age, m = 1, k = 5) +
    age + sex + period, 
  data = dat,
  method = "REML"
)
gam3 <- gam(
  fl_cm ~ s(yday_c, by = age, m = 2, k = 3) +  s(year, by = age, m = 2, k = 5) +
    age + sex + period, 
  data = dat,
  method = "REML"
)



lm1 <- lm(fl_cm ~ yday_c + year_f + age + sex, 
          data = dat)

AIC(gam1, gam2, gam3)

mgcv_splines <- draw(gam1)

pdf(here::here("outputs", "figs", "gam_fit_summary.pdf"))
mgcv_splines
dev.off()


## check validity of different age-specific smoothing structures
gam2 <- gam(
  fl_cm ~ s(yday_c, k = 3) + #s(year, m = 2, k = 5) + 
    s(year, by = age, m = 2, k = 5) +
    age + sex + period, 
  data = dat,
  method = "REML"
)
gam2b <- gam(
  fl_cm ~ s(yday_c, k = 3, by = age) + s(year, m = 2, k = 5, by = age) + 
    # s(year, by = age, m = , k = 5) +
    age + sex + period, 
  data = dat,
  method = "REML"
)
gam3 <- gam(
  fl_cm ~ s(yday_c, k = 3) + s(year, m = 2, k = 5) + 
    age + sex + period, 
  data = dat,
  method = "REML"
)
# aic supports global smooth
AIC(gam1, gam2, gam3)
AIC(gam2, gam2b)

plot(gam2)
plot(gam2b)
plot(gam3)


# check centering categorical variables' effect
c_dat <- fastDummies::dummy_cols(dat %>% droplevels(), 
                                 select_columns = c("sex", "age")) %>% 
  mutate(
    sex_female = as.numeric(scale(sex_female, center = T, scale = F)),
    sex_male = as.numeric(scale(sex_male, center = T, scale = F)),
    age_42 = as.numeric(scale(age_42, center = T, scale = F)),
    age_52 = as.numeric(scale(age_52, center = T, scale = F)),
    age_53 = as.numeric(scale(age_53, center = T, scale = F)),
    age_63 = as.numeric(scale(age_63, center = T, scale = F))
  ) %>% 
  glimpse()
gam2b <- gam(
  fl_cm ~ 0 + s(yday_c, k = 3) + #s(year, m = 2, k = 5) + 
    s(year, by = age, m = 2, k = 5) +
    sex_female + sex_male + age_42 + age_52 + age_53 + age_63, 
  data = c_dat,
  method = "REML"
)
gam2 <- gam(
  fl_cm ~ s(yday_c, k = 3) + #s(year, m = 2, k = 5) + 
    s(year, by = age, m = 2, k = 5) +
    age + sex, 
  data = c_dat,
  method = "REML"
)
