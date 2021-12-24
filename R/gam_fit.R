# Data clean
# 1) Import Skip's first century of data and clean as necessary
# 2) Generate age-/sex-specific annual means using mixed effects model (year as 
# intercept) 
# 3) Estiamte age-/sex-specific trends using GAM
# Created by C Freshwater Dec 23, 2021
# -------------------------------------------------

library(tidyverse)
library(lme4)
library(mgcv)
library(gratia)


# import nass data
dat_in <- read.table(here::here("data", "nasscentury.txt"))
colnames(dat_in) <- c("year", "month", "day", "hart_species", "sex", "age", 
                      "fl")
dat <- dat_in %>% 
  mutate(
    year_f = as.factor(year),
    sex = ifelse(sex == 1, "male", "female"),
    age = as.character(age),
    date = as.POSIXct(paste(day, month, year, sep = "-"),
                      format = "%d-%m-%Y"),
    yday = lubridate::yday(date),
    yday_c = yday - mean(yday, na.rm = TRUE)
  ) %>% 
  #drop secondary age classes and years missing capture data
  filter(
    age %in% c("42", "52", "53", "63"),
    !is.na(yday)
  ) %>% 
  glimpse()


# visualize raw data
ggplot(dat) +
  geom_boxplot(aes(x = year_f, y = fl, fill = sex)) +
  facet_wrap(~age)

trim_dat <- dat %>% 
  sample_n(size = 10000)
ggplot(trim_dat) +
  geom_point(aes(x = yday_c, y = fl), alpha = 0.3) +
  facet_grid(sex~age)
ggplot(trim_dat) +
  geom_histogram(aes(fl)) +
  facet_grid(sex~age)

  
# ANNUAL ESTIMATES FREQUENTIST -------------------------------------------------

# compare fixed effects
mod1a <- lmer(fl ~ yday_c + sex + age + (1 | year_f:age), data = dat,
             REML = FALSE)
mod1b <- lmer(fl ~ yday_c + age + (1 | year_f:age), data = dat,
              REML = FALSE)
AIC(mod1a, mod1b)
# reasonable support for including sex

# fit
mod1 <- lmer(fl ~ yday_c + sex + age + (1 | year_f:age), data = dat,
              REML = TRUE)


#check model diagnostics
plot(mod1, type=c("p","smooth"), col.line=1)
plot(mod1,
     sqrt(abs(resid(.)))~fitted(.),
     type=c("p","smooth"), col.line=1)
lattice::qqmath(mod1)
plot(mod1, rstudent(.) ~ hatvalues(.))
itsadug::acf_resid(mod1)
## evidence of heavy tails; try fitting equivalent model in brms with student-t

# generate predictions using frequentist model
pred_dat <- expand.grid(age = unique(dat$age),
                        yday_c = 0,
                        sex = unique(dat$sex),
                        year_f = unique(dat$year_f)) 
# Generate confidences intervals for predictions based on model
preds <- merTools::predictInterval(mod1, newdata = pred_dat, n.sims = 999)


# Bind confidence intervals to the predictive data for plotting
pred_dat_freq <- cbind(pred_dat, preds) %>%
  # Add a factor representing sample collection to visualize effects 
  # NOTE: period is not modeled, we're simply assigning a color in the plot to 
  # each data set
  mutate(year = as.numeric(as.character(year_f)),
         period = case_when(
           year < 1947 ~ "GC",
           year > 1993 ~ "mod",
           TRUE ~ "bilton"
         ))

# visualize predictions
ggplot(pred_dat_freq, aes(x = year)) +
  # geom_point(aes(y = fit, colour = period)) +
  geom_pointrange(aes(y = fit, ymin = lwr, ymax = upr, colour = period)) +
  facet_grid(sex~age)



# ANNUAL ESTIMATES BAYESIAN ----------------------------------------------------

library(brms)
library(tidybayes)

brms_n <- brm(fl ~ yday_c + sex + age + (1 | year_f:age), data = trim_dat,
              cores = 4, seed = 17,
              iter = 2500, thin = 10, refresh = 0,
              control = list(adapt_delta = 0.90, max_treedepth = 12))
saveRDS(brms_n, here::here("outputs", "data", "brms_fits", "trim_normal.RDS"))
brms_stud <- brm(fl ~ yday_c + sex + age + (1 | year_f:age), data = trim_dat,
                 family = student,
                 cores = 4, seed = 17,
                 iter = 2500, thin = 10, refresh = 0,
                 control = list(adapt_delta = 0.90, max_treedepth = 12))  
saveRDS(brms_stud, here::here("outputs", "data", "brms_fits", "trim_studt.RDS"))

# rhat and neff look good
summary(brms_stud)

trim_dat %>%
  add_residual_draws(brms_n) %>%
  median_qi() %>%
  ggplot(aes(sample = .residual)) +
  geom_qq() +
  geom_qq_line()
trim_dat %>%
  add_residual_draws(brms_stud) %>%
  median_qi() %>%
  ggplot(aes(sample = .residual)) +
  geom_qq() +
  geom_qq_line()

# posterior predictive checks (student-t similar)
pp_check(brms_n)
pp_check(brms_n, type='stat', stat='mean')
pp_check(brms_n, type='error_scatter_avg')
pp_check(brms_n, type='intervals')
pp_check(brms_n, x = 'yday_c', type='error_scatter_avg_vs_x')

# student-t and normal have similar performance to each other and to lmer4;
# pp check looks good but qqplot still wonky; for now use frequentist for 
# simplicity's sake


post_pred <- posterior_predict(brms_n, pred_dat, allow_new_levels = TRUE)
pred_dat_bayes <- pred_dat %>% 
  mutate(
    mean = apply(post_pred, 2, mean),
    low = apply(post_pred, 2, function (x) quantile(x, probs = 0.05)),
    up = apply(post_pred, 2, function (x) quantile(x, probs = 0.95)),
    year = as.numeric(as.character(year_f))
  ) 

ggplot(pred_dat_bayes, aes(x = year)) +
  # geom_point(aes(y = fit, colour = period)) +
  geom_pointrange(aes(y = mean, ymin = low, ymax = up)) +
  facet_grid(sex~age)


# ANNUAL TRENDS ----------------------------------------------------------------

pred_dat_freq$age_sex <- factor(paste(pred_dat_freq$age, pred_dat_freq$sex, sep = "_"))

gam1 <- gam(fit ~ s(year, m = 2) + s(year, by = age, m = 1) + age + sex,
            data = pred_dat_freq)

draw(gam1, residuals = TRUE)

