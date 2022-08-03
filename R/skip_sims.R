## Skip Questions
# Fit and run alternative models to address Skip's questions
# July 28, 2022


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
    # add dummy variables for sampling period (necessary to make "average"
    # predictions)
    period_b = ifelse(period == "Bilton", 1, 0),
    period_m = ifelse(period == "Monkley Dump", 1, 0),
    period_n = ifelse(period == "Nisga'a", 1, 0),
    age_sex = paste(age, sex, sep = "_") %>% 
      as.factor()
  ) %>% 
  droplevels()
levels(dat$age_f) <- c("4[2]", "5[2]", "5[3]", "6[3]")
dat$period_b_cent <- dat$period_b - mean(dat$period_b)
dat$period_m_cent <- dat$period_m - mean(dat$period_m)
dat$period_n_cent <- dat$period_n - mean(dat$period_n)


## PERIOD EFFECT ---------------------------------------------------------------

fit1 <- mgcv::gam(fl_cm ~ s(yday_c, by = age, m = 2) +
                s(year, by = age, m = 2) +
                period_b_cent + period_m_cent + period_n_cent + age + sex,
              data = dat
              )
fit2 <- mgcv::gam(fl_cm ~ s(yday_c, by = age, m = 2) +
                    s(year, by = age, m = 2) +
                    age + sex,
                  data = dat
)
AIC(fit1, fit2)


## PERTURB ---------------------------------------------------------------------

library(perturb)
library(car)


fit_gam3 <- mgcv::gam(fl_cm ~ s(yday_c, by = age, m = 2) +
                       s(year, by = age, m = 2) +
                       period + #period_b_cent + period_m_cent + period_n_cent +
                       age + sex,
                     data = dat
)

concur

p_gam <- perturb(fit_gam3, pvars = c("year", "yday_c"), prange = c(1, 1),
                 pfac = list("period", pcnt = 95))

fit1 <- lm(fl_cm ~ year + age + sex + period + yday_c,
           data = dat)
p_lm <- perturb(fit1, pvars = c("year"), prange = c(1),
                 pfac = list("period", pcnt = 95))
# looks reasonable


## EXPLORE ALTERNATIVE STRUCTURES ----------------------------------------------

trim_dat <- dat %>% 
  filter(age == "42", sex == "female") %>% 
  mutate(year_c = year - mean(year))

fit1 <- lm(fl_cm ~ year_c + period_b_cent + period_m_cent + period_n_cent,
           data = trim_dat)
fit2 <- lm(fl_cm ~ year_c,
           data = trim_dat)
fit3 <- lm(fl_cm ~ period_b_cent + period_m_cent + period_n_cent,
           data = trim_dat)

summary(fit1)
summary(fit2)
summary(fit3)


fit1 <- lm(fl_cm ~ year + period + yday_c + age + sex,
           data = dat)
car::vif(fit1)



## FIT FOR CIs ----------------------------------------------------------------


# make fake mesh (necessary in current branch)
dat$x <- runif(nrow(dat))
dat$y <- runif(nrow(dat))
dum_mesh <- make_mesh(dat, c("x", "y"), cutoff = 1000)

fit_year <- sdmTMB(fl_cm ~ s(yday_c, by = age, m = 2) +
                s(year, by = age, m = 2) +
                period_b_cent + period_m_cent + period_n_cent + age + sex,
              data = dat,
              mesh = dum_mesh,
              spatial = "off",
              spatiotemporal = "off",
              control = sdmTMBcontrol(
                nlminb_loops = 2,
                newton_loops = 2
              ))

fit_no_year <- sdmTMB(fl_cm ~ s(yday_c, by = age, m = 2) +
                        period_b_cent + period_m_cent + period_n_cent + age + sex,
                      data = dat,
                      mesh = dum_mesh,
                      spatial = "off",
                      spatiotemporal = "off",
                      control = sdmTMBcontrol(
                        nlminb_loops = 2,
                        newton_loops = 2
                      ))
fit_no_year2 <- sdmTMB(fl_cm ~ period_b_cent + period_m_cent + period_n_cent + 
                         age + sex,
                       data = dat,
                       mesh = dum_mesh,
                       spatial = "off",
                       spatiotemporal = "off",
                       control = sdmTMBcontrol(
                         nlminb_loops = 2,
                         newton_loops = 2
                       ))
fit_no_year3 <- sdmTMB(fl_cm ~ period_b_cent + period_m_cent + period_n_cent,
                       data = dat,
                       mesh = dum_mesh,
                       spatial = "off",
                       spatiotemporal = "off",
                       control = sdmTMBcontrol(
                         nlminb_loops = 2,
                         newton_loops = 2
                       ))

tidy(fit_no_year3, effects = "fixed") %>% 
  filter(!term == "(Intercept)") %>% 
  mutate(
    low = estimate + (qnorm(0.025) * std.error),
    up = estimate + (qnorm(0.975) * std.error),
    group = case_when(
      grepl("period", term) ~ "period",
      grepl("age", term) ~ "age",
      grepl("sex", term) ~ "sex",
      TRUE ~ term
    ),
    group_f = factor(group, levels = c("sex", "age", "period"), 
                     labels = c("sex", "age", "sampling\nperiod")),
    term = fct_reorder(as.factor(term), as.numeric(group_f))
  ) %>% 
  glimpse()

dat %>% 
  filter(age == "42",
         sex == "female") %>% 
  group_by(period) %>%
  summarize(n = n(),
            sd = sd(fl_cm),
            se = sd / sqrt(n))
