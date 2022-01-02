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
  ) 

# predictive dataset
pred_dat <- expand.grid(age = unique(dat$age),
                        yday_c = 0,
                        sex = unique(dat$sex),
                        year_f = unique(dat$year_f)) %>% 
  mutate(
    year = as.numeric(as.character(year_f))
  )


# visualize raw data
ggplot(dat) +
  geom_boxplot(aes(x = year_f, y = fl, fill = sex)) +
  facet_wrap(~age)

trim_dat <- dat %>% 
  sample_n(size = 20000)
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
preds <- merTools::predictInterval(mod1, newdata = pred_dat, n.sims = 999)


# Bind confidence intervals to the predictive data for plotting
pred_dat_freq <- cbind(pred_dat, preds) %>%
  # Add a factor representing sample collection to visualize effects 
  # NOTE: period is not modeled, we're simply assigning a color in the plot to 
  # each data set
  mutate(period = case_when(
           year < 1947 ~ "GC",
           year > 1993 ~ "mod",
           TRUE ~ "bilton"
         ))

# visualize predictions
ggplot(pred_dat_freq, aes(x = year)) +
  geom_point(aes(y = fit, colour = period)) +
  # geom_pointrange(aes(y = fit, ymin = lwr, ymax = upr, colour = period)) +
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
    up = apply(post_pred, 2, function (x) quantile(x, probs = 0.95))
    ) 

ggplot(pred_dat_bayes, aes(x = year)) +
  # geom_point(aes(y = fit, colour = period)) +
  geom_pointrange(aes(y = mean, ymin = low, ymax = up)) +
  facet_grid(sex~age)


# ANNUAL TRENDS ----------------------------------------------------------------

# fit gam to annual estimates to identify age specific trends
pred_dat_freq$age_sex <- factor(paste(pred_dat_freq$age, pred_dat_freq$sex,
                                      sep = "_"))
pred_dat$age_sex <- factor(paste(pred_dat$age, pred_dat$sex,
                                      sep = "_"))

gam1 <- gam(fit ~ s(year, m = 2) + s(year, by = age_sex, m = 1) + age + sex,
            data = pred_dat_freq)
gam2 <- gam(fit ~ s(year, by = age_sex, m = 2) + age + sex,
            data = pred_dat_freq)
gam3 <- gam(fit ~ s(year, m = 2) + s(year, by = age, m = 1) + age + sex,
            data = pred_dat_freq)

appraise(gam1)
qq_plot(gam1, method = "simulate")
draw(gam1, residuals = TRUE)


# calculate derivatives representing significant increases
der <- derivatives(gam1, type = "central") 
sig_der_years <- der %>% 
  filter(smooth == "s(year)") %>% 
  mutate(
    year_f = as.factor(floor(data + 0.00001)),
    sig_change = case_when(
      lower > 0 ~ "sig_increase",
      upper < 0 ~ "sig_decrease", 
      TRUE ~ "nonsig")
  ) 

pred_gam_der <- ggplot(der, aes(x = data, y = derivative)) + 
  geom_line(size = 1.2) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
  facet_wrap(~smooth) +
  ggsidekick::theme_sleek()


pred_gam <- predict.gam(gam1, pred_dat, se.fit = TRUE)
pred_dat_gam <- pred_dat %>% 
  mutate(pred_fl = as.numeric(pred_gam$fit),
         pred_se = as.numeric(pred_gam$se.fit),
         up = pred_fl + (1.96 * pred_se),
         low = pred_fl - (1.96 * pred_se))# %>% 
  # left_join(., sig_der_years %>% select(year_f, sig_change), by = "year_f") %>% 
  # mutate(sig_change2 = ifelse(i))

col_pal <- c("grey60", "blue", "red")
names(col_pal) <- unique(pred_dat_gam$sig_change)


# predictions
pred_gam_plot <- ggplot(pred_dat_gam, aes(x = year, y = pred_fl)) +
  geom_line(aes(col = sex)) +
  geom_ribbon(aes(ymin = low, ymax = up, fill = sex), alpha = 0.4) +
  geom_point(data = pred_dat_freq, aes(x = year, y = fit, col = sex),
             alpha = 0.4) +
  facet_wrap(~age) +
  ggsidekick::theme_sleek() +
  labs(y = "Predicted Length", x = "Year")


pdf(here::here("outputs", "figs", "annual_gam_preds.pdf"))
pred_gam_plot
pred_gam_der
dev.off()


# calculate decline relative to long-term (beginning to 2000) mean
final_mean <- pred_dat_gam %>% 
  filter(year > 2014) %>% 
  group_by(age_sex) %>% 
  summarize(new_mean_fl = mean(pred_fl)) 
  
pred_dat_gam %>% 
  filter(year < 2001) %>% 
  group_by(age_sex) %>% 
  summarize(mean_fl = mean(pred_fl)) %>% 
  left_join(., final_mean, by = "age_sex") %>% 
  mutate(diff = mean_fl - new_mean_fl)


# INDIVIDUAL GAM ---------------------------------------------------------------

# compare GAM fit to individual level data relative to annual estimates (includes
# covariate for run timing)
dat$age_sex <- factor(paste(dat$age, dat$sex, sep = "_"))
gam_ind <- gam(fl ~ s(yday_c) + s(year, m = 2) + s(year, by = age_sex, m = 1) +
                 age + sex, data = dat)
gam_ind2 <- gam(fl ~ s(yday_c) + s(year, m = 2) + s(year, by = age, m = 1) +
                 age + sex, data = dat)
gam_ind3 <- gam(fl ~ s(yday_c) + s(year, m = 2) + age + sex, data = dat)
AIC(gam_ind, gam_ind2, gam_ind3)


#check model diagnostics
appraise(gam_ind)
qq_plot(gam_ind, method = "simulate")
draw(gam_ind, residuals = TRUE)


pred_gam_ind <- predict.gam(gam_ind, pred_dat, se.fit = TRUE)
pred_dat_gam_ind <- pred_dat %>% 
  mutate(pred_fl = as.numeric(pred_gam_ind$fit),
         pred_se = as.numeric(pred_gam_ind$se.fit),
         up = pred_fl + (1.96 * pred_se),
         low = pred_fl - (1.96 * pred_se)) 


pred_gam_ind_fig <- ggplot(pred_dat_gam_ind, aes(x = year, y = pred_fl)) +
  geom_line(aes(col = sex), size = 1) +
  geom_ribbon(aes(ymin = low, ymax = up, fill = sex), alpha = 0.4) +
  geom_jitter(data = trim_dat, aes(x = year, y = fl, col = sex),
              alpha = 0.1) +
  facet_wrap(~age) +
  ggsidekick::theme_sleek() +
  labs(y = "Predicted Length", x = "Year")


# individual based GAM's derivatives
der_ind <- derivatives(gam_ind, type = "central") 
ggplot(der_ind %>% filter(var == "year"), aes(x = data, y = derivative)) + 
  geom_line(size = 1.2) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
  facet_wrap(~smooth)


pdf(here::here("outputs", "figs", "ind_gam_preds.pdf"))
pred_gam_ind_fig
dev.off()


# INDIVIDUAL BAYESIAN GAM ------------------------------------------------------

#check priors (WAY TOO LARGE TO FIT)
brms::get_prior(fl ~ s(yday_c) + s(year, m = 2, k = 4) + 
                  s(year, by = age, m = 1, k = 4) +
                  age + sex,
                data = trim_dat)
brm1 <- brm(
  fl ~ s(yday_c) + s(year, m = 2) + s(year, by = age, m = 1) +
       age + sex,
  data = trim_dat, seed = 17,
  iter = 2000, warmup = 750, thin = 10, cores = 4, refresh = 0,
  # iter = 500, warmup = 1000, thin = 10, cores = 4, refresh = 0,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  prior=c(prior(normal(613, 20), class="Intercept"),
          prior(normal(0, 2.5), class="b"))
)

saveRDS(brm1, here::here("outputs", "data", "brms_fits", "trim_ind.rds"))
