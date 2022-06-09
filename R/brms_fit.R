# Fit brms model to century of Nass data 
# Extension of gam_fit 
# 1) Import Skip's first century of data and clean as necessary
# 2) Estimate temporal trend using GAM and brms to account for changes in
# variability with sampling regime
# Created by C Freshwater Dec 23, 2021
# Update w/ edited data
# -------------------------------------------------

library(tidyverse)
library(brms)
library(tidybayes)
library(posterior)


# import nass data
# dat_old <- read.table(here::here("data", "nasscentury.txt"))
dat_in <- read.table(here::here("data", "nasscenturyv3.txt"))
colnames(dat_in) <- c("year", "month", "day", "hart_species", "sex", "age", 
                      "fl")
# empty_yrs <- dat_n %>% filter(is.na(mean_fl)) %>% pull(year_f) %>% unique()


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
    data_group = sample.int(8, nrow(.), replace = T)
  ) %>% 
  droplevels()
levels(dat$age_f) <- c("4[2]", "5[2]", "5[3]", "6[3]")


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
  filter(year_f %in% empty_yrs,
         sex == "female") %>% 
  select(year_f, mean_fl = value, period, lo, up, age_f, data)
levels(dat_avg$age_f) <- c("4[2]", "5[2]", "5[3]", "6[3]")

mean_dat_plotting <- mean_dat %>% 
  mutate(data = "individual measurements") %>% 
  filter(sex == "female",
         !is.na(period)) %>% 
  select(year_f, mean_fl, period, lo, up, age_f, data) %>%
  rbind(., dat_avg) %>% 
  group_by(age_f) %>% 
  mutate(
    mean_fl = mean_fl / 10,
    lo = lo / 10,
    up = up / 10,
    overall_mean = mean(mean_fl, na.rm = T)) %>% 
  ungroup() %>% 
  droplevels() 

shape_pal <- c(22, 21)
names(shape_pal) <- c("annual average", "individual measurements")

annual_dot <- ggplot(
  mean_dat_plotting, 
  aes(x = as.numeric(as.character(year_f)))
) +
  geom_pointrange(aes(y = mean_fl, fill = period, ymin = lo, ymax = up, 
                      shape = data)) +
  facet_wrap(~age_f, labeller = label_parsed) +
  ggsidekick::theme_sleek() +
  labs(x = "Year", y = "Fork Length") +
  scale_x_continuous(
    breaks = seq(1915, 2015, by = 20),
    expand = c(0.02, 0.02)
  ) +
  # geom_vline(xintercept = 1966, col = "red", lty = 2) +
  geom_hline(aes(yintercept = overall_mean), lty = 2) +
  scale_fill_discrete(name = "Sampling\nPeriod") +
  scale_shape_manual(values = shape_pal, name = "Dataset") +
  guides(fill = guide_legend(override.aes = list(shape = 21)))


# annual_sd_dot <- ggplot(mean_dat, aes(x = year_f)) +
#   geom_point(aes(y = sd_fl, fill = period), shape = 21) +
#   facet_grid(~age) +
#   ggsidekick::theme_sleek() +
#   labs(x = "Year", y = "SD Fork Length") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(mean_dat, aes(x = year_f)) +
  geom_point(aes(y = mean_yday, fill = period), shape = 21) +
  facet_wrap(~age) +
  ggsidekick::theme_sleek() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(mean_dat, aes(x = year_f)) +
  geom_point(aes(y = sd_yday, fill = period), shape = 21) +
  facet_wrap(~age) +
  ggsidekick::theme_sleek() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

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


split_dat <- split(dat, dat$data_group)


tictoc::tic()
brm1 <- brm_multiple(
  bf(
    fl_cm ~ s(yday_c, m = 2, k = 3) + s(yday_c, by = age, m = 1, k = 3) +
      s(year, m = 2, k = 5) + s(year, by = age, m = 1, k = 5) +
      age + sex + period,
    sigma ~ period
  ),
  data = split_dat, seed = 17,
  iter = 1800, thin = 10, chains = 4, warmup = 750, #refresh = 0,  
  control = list(adapt_delta = 0.97, max_treedepth = 14
                 ),
  prior=c(prior(normal(60, 10), class="Intercept"),
          prior(normal(5, 10), class="b", coef = "age52"),
          prior(normal(5, 10), class="b", coef = "age53"),
          prior(normal(5, 10), class="b", coef = "age63"),
          prior(normal(5, 10), class="b", coef = "sexmale"),
          prior(normal(0, 10), class="b", coef = "periodBilton"),
          prior(normal(0, 10), class="b", coef = "periodMonkleyDump"),
          prior(normal(0, 10), class="b", coef = "periodNisgaa"),
          prior(normal(0, 5), class="b", coef = "syear_1"),
          prior(normal(0, 5), class="b", coef = "syday_c_1"),
          prior(normal(0, 2.5), class="b", dpar = "sigma"),
          prior(exponential(0.5), class = "Intercept", dpar = "sigma")
  )
)
tictoc::toc()

saveRDS(brm1, here::here("outputs", "data", "brms_fits", "ind_ls_v3.rds"))

post <- as.array(brm1)
bayesplot::mcmc_trace(post)
round(brm1$rhats, 3)
conditional_effects(brm1)

resids <- dat %>%
  droplevels() %>% 
  add_residual_draws(brm1) 

# time series of residuals
mean_resids <- resids %>% 
  group_by(year, month, day, sex, age, fl, year_f, age_f, yday, period, .row) %>% 
  summarize(mean_resid = mean(.residual))

resid_point <- mean_resids %>% 
  ungroup() %>% 
  sample_n(size = 30000, replace = FALSE) %>%
  ggplot(.) +
  geom_point(aes(x = year_f, y = mean_resid, fill = period), 
             alpha = 0.4, shape = 21) +
  facet_grid(age_f ~ sex) +
  ggsidekick::theme_sleek()

pdf(here::here("outputs", "figs", "resids_50k_samps.pdf"))
resid_point
dev.off()


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

brm1 <- readRDS(here::here("outputs", "data", "brms_fits", "ind_ls.rds"))

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


## ANNUAL TRENDS ---------------------------------------------------------------

# first predictive dataset only has years observed in dataset; these are paired
# with appropriate periods to estimate overall effects
pred_dat <- expand.grid(
  year = seq(min(dat$year), max(dat$year), by = 0.1),
  age = unique(dat$age),
  yday_c = 0,
  sex = unique(dat$sex),
  period = "Gilbert-Clemens"
) %>%
  left_join(., dat %>% dplyr::select(age, age_f) %>% distinct(), by = "age") %>% 
  droplevels()

# post_pred <- posterior_predict(brm1, pred_dat, allow_new_levels = TRUE)
# pred_dat_bayes <- pred_dat %>% 
#   mutate(
#     median = apply(post_pred, 2, median),
#     low = apply(post_pred, 2, function (x) quantile(x, probs = 0.05)),
#     up = apply(post_pred, 2, function (x) quantile(x, probs = 0.95))
#   ) 
# 
# smooth_pred1 <- ggplot(pred_dat_bayes, aes(x = year, y = median)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = low, ymax = up), alpha = 0.4) +
#   facet_grid(sex~age) + 
#   ggsidekick::theme_sleek() +
#   labs(x = "Year", y = "Fork Length (cm)")
# 
# png(here::here("outputs", "figs", "smooth_1period.png"), 
#     height = 4, width = 7.5, units = "in", res = 250)
# smooth_pred1
# dev.off()


# as above but with fitted draws, not posterior predictions
fit_draws <- add_epred_draws(pred_dat, brm1, re_formula = NULL,
                   allow_new_levels = TRUE) 
fit_draws2 <- fit_draws %>% 
  group_by(year, age, age_f, yday_c, sex, period) %>% 
  summarize(
    median = median(.value),
    low = quantile(.value, probs = 0.05),
    up = quantile(.value, probs = 0.95),
    .groups = "drop"
  ) 
# saveRDS(fit_draws2,
#         here::here("outputs", "data", "brms_fits", "brms_post_preds.rds"))


fit_pred <- ggplot(fit_draws2 %>% filter(sex == "female"), 
                   aes(x = year, median)) +
  geom_line() +
  geom_ribbon(aes(ymin = low, ymax = up), alpha = 0.4) +
  geom_vline(xintercept = 1966, col = "red", lty = 2) +
  ggsidekick::theme_sleek() +
  labs(x = "Year", y = "Fork Length (cm)") +
  facet_wrap(~age_f, labeller = label_parsed) +
  scale_x_continuous(
    breaks = seq(1915, 2015, by = 20),
    expand = c(0.02, 0.02)
  )
  

png(here::here("outputs", "figs", "smooth_1period_fit.png"), 
    height = 4, width = 5, units = "in", res = 250)
fit_pred
dev.off()


# second predictive dataset associates years with periods (approximate) and 
# includes residual variability
pred_dat2 <- pred_dat %>% 
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
                         "Nisga'a")
  )


preds <- posterior_predict(brm1, pred_dat2, allow_new_levels = TRUE) 
post_pred2 <- pred_dat2 %>% 
  mutate(
    mean = apply(preds, 2, mean),
    low = apply(preds, 2, function (x) quantile(x, probs = 0.05)),
    up = apply(preds, 2, function (x) quantile(x, probs = 0.95))
  ) %>% 
  filter(sex == "female")
post_ribbon2 <- ggplot(post_pred2, aes(x = year, y = mean)) +
  geom_line() +
  geom_vline(xintercept = 1966, col = "red", lty = 2) +
  geom_ribbon(aes(ymin = low, ymax = up), alpha = 0.4) +
  facet_wrap(~age_f, labeller = label_parsed) +
  ggsidekick::theme_sleek() +
  scale_x_continuous(
    breaks = seq(1915, 2015, by = 20),
    expand = c(0.02, 0.02)
  )

smooth_pred2 <- post_ribbon2 +
  geom_point(data = dat %>% 
               filter(!is.na(period),
                      sex == "female") %>% 
               sample_n(1000), 
             aes(x = year, y = fl_cm, fill = period), 
             shape = 21, alpha = 0.4) +
  scale_fill_discrete(name = "Sampling\nPeriod") +
  labs(x = "Year", y = "Fork Length (cm)")

png(here::here("outputs", "figs", "smooth_4periods.png"), 
    height = 4, width = 5.5,
    units = "in", res = 250)
smooth_pred2
dev.off()


#duplicate posterior predictions plot with out of sample means from Foskett
fosk_dat <- read.csv(here::here("data", "foskett_mean_size.csv"),
                     fileEncoding = 'UTF-8-BOM') %>% 
  janitor::clean_names() %>%
  mutate(mean_fl = len_mm / 10,
         age = as.character(age)) %>% 
  left_join(., dat %>% dplyr::select(age, age_f) %>% distinct(), by = "age") 
  

png(here::here("outputs", "figs", "smooth_4periods_foskett.png"), 
    height = 4, width = 5.5,
    units = "in", res = 250)
smooth_pred2 +
  geom_point(data = fosk_dat %>% filter(sex == "F"), 
             aes(x = year, y = mean_fl), fill = "red", shape = 21) +
  scale_fill_discrete(name = "Sampling\nPeriod") +
  labs(x = "Year", y = "Fork Length (cm)")
dev.off()


# visualize year-day effect
pred_dat3 <- expand.grid(
  year = mean(dat$year),
  yday_c = seq(-46, 64, by = 0.1),
  age = unique(dat$age),
  sex = "female",
  period = "Gilbert-Clemens"
) %>% 
  mutate(
    yday = yday_c + mean(dat$yday)
  ) %>% 
  left_join(., dat %>% select(age, age_f) %>% distinct(), by = "age") %>% 
  droplevels()

# pp_yday <- pred_dat3 %>%
#   add_predicted_draws(model = brm1, re_formula = NULL,
#                       allow_new_levels = TRUE) 
# 
# pp_yday_plot <- ggplot(pp_yday, aes(x = yday_c)) +
#   stat_lineribbon(aes(y = .prediction),
#                   .width = c(.9, .5), alpha = 0.25) +
#   ggsidekick::theme_sleek() +
#   labs(x = "Centered Year Day", y = "Fork Length (cm)") +
#   facet_wrap(~age)
# 
# pp_yday_plot +
#   geom_point(data = dat %>% 
#                filter(#period == "Gilbert-Clemens",
#                  sex == "female") %>% 
#                sample_n(1000), 
#              aes(x = yday_c, y = fl_cm, fill = period), 
#              shape = 21, alpha = 0.4)

# png(here::here("outputs", "figs", "smooth_yday.png"), 
#     height = 4.5, width = 4.5, units = "in", res = 250)
# pp_yday
# dev.off()


fit_draws_day <- pred_dat3 %>% 
  add_fitted_draws(model = brm1, re_formula = NULL,
                   allow_new_levels = TRUE)  %>% 
  group_by(year, age, age_f, yday, yday_c, sex, period) %>% 
  summarize(
    median = median(.value),
    low = quantile(.value, probs = 0.05),
    up = quantile(.value, probs = 0.95),
    .groups = "drop"
  ) 

fit_yday <- ggplot(fit_draws_day, aes(x = yday, y = median)) +
  geom_line() +
  geom_ribbon(aes(ymin = low, ymax = up), alpha = 0.4) +
  ggsidekick::theme_sleek() +
  labs(x = "Sampling Date", y = "Fork Length (cm)") +
  facet_wrap(~age_f, labeller = label_parsed) +
  scale_x_continuous(
    breaks = c(166, 213, 257),
    expand = c(0.01, 0.01),
    labels = c("Jun 15", "Aug 1", "Sep 15")
  )

png(here::here("outputs", "figs", "fit_smooth_yday.png"), 
    height = 4, width = 5, units = "in", res = 250)
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


### Residuals ------------------------------------------------------------------

resids <- resid(brm1)

dum$resids <- resids[, "Estimate"]
dum$cohort <- ifelse((dum$year %% 2) == 0, "even", "odd")

pdf(here::here("outputs", "figs", "even_odd_resid.pdf"))
ggplot(dum) +
  geom_boxplot(aes(x = cohort, y = resids)) +
  ggsidekick::theme_sleek() +
  facet_grid(sex~age)
dev.off()


### Fecundity effects ----------------------------------------------------------

## Mason and West regression 
b1 = 11.52 #11.10 (other river) 
a1 = -2152 #2015.8 (other river)

pred_dat5 <- pred_dat %>% filter(sex == "female")

post_fem <- posterior_predict(
  brm1, 
  pred_dat5,
  allow_new_levels = TRUE
) 
  
# calculate fecundity; includes conversion to mm and POH
post_fec <- ((b1 * (post_fem * 10 * 0.84)) + a1) / 1000

age53 <- post_fec[ , which(pred_dat5$age == "53")]
age52 <- post_fec[ , which(pred_dat5$age == "52")]
age42 <- post_fec[ , which(pred_dat5$age == "42")]
age63 <- post_fec[ , which(pred_dat5$age == "63")]
age_list <- list(age53, age52, age42, age63)
age_names <- c("53", "52", "42", "63")


# differences in fecundity through time
purrr::map2(age_list, age_names, function (x, y) {
  #calc time series mean for each iteration
  mu <- apply(x, 1, mean)
  #difference from mean
  diff <- x[ , ncol(x)] - mu
  #difference from first obs
  diff2 <- x[ , ncol(x)] - x[ , 1]
  rel_dif <- diff2 / mu
  data.frame(
    mean_diff_mean = mean(diff),
    mean_rel_diff = mean(rel_dif),
    mean_diff_start = mean(diff2),
    low = quantile(diff2, probs = 0.05),
    high = quantile(diff2, probs = 0.95),
    age = y
  )
}) %>% 
  bind_rows


pred_dat5 %>% 
  mutate(
    median = apply(post_fec, 2, median),
    low = apply(post_fec, 2, function (x) quantile(x, probs = 0.05)),
    up = apply(post_fec, 2, function (x) quantile(x, probs = 0.95))
  ) %>% 
  ggplot(., aes(x = year, y = median)) +
  geom_line() +
  geom_ribbon(aes(ymin = low, ymax = up), alpha = 0.4) +
  ggsidekick::theme_sleek() +
  labs(x = "Year", y = "Fecundity (1000s eggs)") +
  facet_wrap(~age_f, labeller = label_parsed) +
  scale_x_continuous(
    expand = c(0.01, 0.01)
  )


# as above but excluding residual variation
fit_fem <- pred_dat5 %>% 
  add_fitted_draws(model = brm1, re_formula = NULL, allow_new_levels = TRUE)

tt <- fit_fem %>% 
  mutate(egg_count = ((b1 * (.value * 10 * 0.84)) + a1) / 1000) %>% 
  group_by(age) %>% 
  mutate(mean_egg_count = mean(egg_count)) %>% 
  filter(year == max(pred_dat5$year) | year == min(pred_dat5$year),
         period == "Gilbert-Clemens") %>%
  ungroup() %>% 
  select(-c(.chain, .iteration, .row, .value)) %>%
  pivot_wider(names_from = "year", values_from = "egg_count") %>% 
  mutate(diff = `2019` - `1914`,
         rel_diff = diff / mean_egg_count) 


diff_fecundity <- tt %>% 
  group_by(age_f, sex) %>% 
  summarize(
    median = median(diff),
    low = quantile(diff, probs = 0.05),
    up = quantile(diff, probs = 0.95),
    rel_median = median(rel_diff),
    rel_low = quantile(rel_diff, probs = 0.05),
    rel_up = quantile(rel_diff, probs = 0.95),
    .groups = "drop"
  ) %>% 
  mutate(
    age_f = factor(
      age_f, 
      labels = c(expression("4"["2"]), expression("5"["2"]), 
                 expression("5"["3"]), expression("6"["3"]))
    )
  )


png(here::here("outputs", "figs", "fecundity_diff.png"), 
    height = 4, width = 5, units = "in", res = 250)
ggplot(diff_fecundity) +
  geom_pointrange(aes(x = age_f, y = rel_median, ymin = rel_low, ymax = rel_up)) +
  geom_hline(aes(yintercept = 0), lty = 2) +
  ggsidekick::theme_sleek() +
  labs(y = "Difference in Relative Fecundity (2019-1914)") +
  scale_x_discrete(
    "Age",
    labels = c(expression("4"["2"]), expression("5"["2"]), 
               expression("5"["3"]), expression("6"["3"]))
    ) +
  theme(legend.position = "none")
dev.off()



fit_fem %>% 
  mutate(egg_count = ((b1 * (.value * 10 * 0.84)) + a1) / 1000) %>% 
  group_by(year, age, age_f) %>% 
  summarize(
    median = median(egg_count),
    low = quantile(egg_count, probs = 0.05),
    up = quantile(egg_count, probs = 0.95),
    .groups = "drop"
  )  %>% 
  ggplot(., aes(x = year, y = median)) +
  geom_line() +
  geom_ribbon(aes(ymin = low, ymax = up), alpha = 0.4) +
  ggsidekick::theme_sleek() +
  labs(x = "Year", y = "Fecundity (1000s eggs)") +
  facet_wrap(~age_f, labeller = label_parsed) +
  scale_x_continuous(
    expand = c(0.01, 0.01)
  )


### MGCV COMPARE ---------------------------------------------------------------

library(mgcv)
library(gratia)

gam1 <- gam(
  fl_cm ~ s(yday_c, m = 2, k = 3) + s(yday_c, by = age, m = 1, k = 3) + 
    s(year, m = 2, k = 5) + s(year, by = age, m = 1, k = 5) +
    age + sex + period, 
  data = dat %>% sample_n(., 10000),
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
